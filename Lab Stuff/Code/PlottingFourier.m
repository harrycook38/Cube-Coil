clear all
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Preamble
cd('/Users/Harry/Documents/GitHub/Cube-Coil/Lab Stuff/Data/11-5-22/Origin/')
fname = 'QZFM_1.lvm';   %Filename
Fs = 1213;              %Sampling Frequency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Main
%Load in time stamps and data
[rtime,rdata] = eb_read_lvm(fname);
fgmout = rdata(:,1:3);      %FGM Data

Bout = (1e-6*fgmout)./0.1;       %Field (Conversion 0.1 V/uT)
Bcor = Bout - mean(Bout,1); %Correct for Vertical Offset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fourier
N = length(Bcor);
fti = fft(Bcor);
ft = fti(1:N/2+1,:);
freq = 0:Fs/N:Fs/2;

%Find dominant 3 frequencies.
in = find(freq == 100);
ind = max(abs(ft(1:in,:)),[],1);

fx = find(abs(ft(:,1)) == ind(1));
fy = find(abs(ft(:,2)) == ind(2));
fz = find(abs(ft(:,3)) == ind(3));

mu = [freq(1,fx); freq(1,fy); freq(1,fz);];
for i = 1:3
    if mu(i) == 50 %If powerline found, set to zero
        mu(i) = 0;
    end
end 
disp(['Determined Frequencies: ' num2str(mu')])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Filter & Inverse ft
sd = 0.5;                               %Width of filter
a = (1./sqrt(2.*pi.*sd));               %Normalised Amplitude
gau = a.*exp(-(freq-mu).^8/(2*sd^2));   %Modified Gaussian

filt = gau'.*ft;                        %Filtering
inv = ifft(filt,'symmetric');           %Inverse Transform

redt = rtime(1:2:end);
Bshort = Bcor(1:2:end,:);

%Depending on datasize (datalength missmatch correction)
if length(redt) ~= length(inv)
    inv = inv(1:end-1,:);
elseif length(redt)  ~= length(Bcor)/2
    Bcor = Bcor(1:end-1,:);
end

psd = (1/(Fs*N)).*ft.*conj(ft);
psd(2:end-1,:) = 2*psd(2:end-1,:);


cd('/Users/Harry/Documents/GitHub/Cube-Coil/Lab Stuff/PlotData/Fourier Plots')

% subplot(2,2,1)
f1 = figure;
plot(redt(1015:1200),Bshort(1015:1200,2),'LineWidth',2)
xlabel('Time','FontSize',15); ylabel('Uncorrected B Field','FontSize',15);
pbaspect([2 1 1])
saveas(f1,'1','png')

% subplot(2,2,2)
f2 = figure;
plot(freq(1,1:in/2),10.*log10(psd(1:in/2,2)))
ylim([-300 -100])
xlabel('Frequency','FontSize',15); ylabel('Power Spectrum','FontSize',15);
pbaspect([2 1 1])
saveas(f2,'2','png')


% subplot(2,2,3)
f3 = figure;
plot(freq(1,1:in/2),gau(2,1:in/2),'LineWidth',2)
xlabel('Frequency','FontSize',15); ylabel('Gaussian Amplitude','FontSize',15);
pbaspect([2 1 1])
saveas(f3,'3','png')

% subplot(2,2,4)
f4 = figure;
plot(redt(1015:1200),inv(1015:1200,2),'r','LineWidth',2)
xlabel('Time','FontSize',15); ylabel('Filtered B Field','FontSize',15);
pbaspect([2 1 1])
saveas(f4,'4','png')


