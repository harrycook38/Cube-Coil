clear all
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Preamble
cd('/Users/Harry/Documents/GitHub/Cube-Coil/Lab Stuff/Data/10-5-22/zy-uniformity')
fname = 'QZFM_22.lvm';   %Filename
Fs = 1213;              %Sampling Frequency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Main
%Load in time stamps and data
[rtime,rdata] = eb_read_lvm(fname);
fgmout = rdata(:,1:3);      %FGM Data
%Iout = rdata(:,4:6)./10;   %Current out from loops
% figure(1)
% plot(rtime,Iout)
% title("Current outputs from coils")

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

if mu(2) == 50 %If powerline found, set to zero
    mu(2) = 0;
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

for d = 1 %Plotting... 
figure(2)
subplot(3,1,1)
plot(redt,inv(:,1))
hold on
plot(redt,Bshort(:,1))

subplot(3,1,2)
plot(redt,inv(:,2))
hold on
plot(redt,Bshort(:,2))

subplot(3,1,3)
plot(redt,inv(:,3))
hold on
plot(redt,Bshort(:,3))
end
%^^^^^Y axis still out of phase, how to remove phase?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RMS Output
xrms = rms(inv(:,1));
yrms = rms(inv(:,2));
zrms = rms(inv(:,3));
disp(['RMS values: x: ', num2str(xrms),' y: ',num2str(yrms),' z: ',num2str(zrms)])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Creating Magnetic field vector
Bvec = inv; 
or = [0 0 0]; %Origin
top = max(max(max(Bvec)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Angles
%y is the vertical
%Not sure about these yet:
theta = acosd(Bvec(:,2));
phi = atand(Bvec(:,1)./Bvec(:,3));

figure(3)
subplot(2,1,1)
plot(redt,phi)
subplot(2,1,2)
plot(redt,theta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(4)
pl1 = quiver3(or(1),or(2),or(3),Bvec(1,1),Bvec(1,2),Bvec(1,3),"LineWidth",3);
hold on
pl2 = plot3(Bvec(1,1),Bvec(1,2),Bvec(1,3),'r-');
xlim([-1.2*top 1.2*top]); ylim([-1.2*top 1.2*top]);zlim([-1.2*top 1.2*top]);
axis square

for j = 1:length(rtime(1:2:end))
    pl1.UData = Bvec(j,1);
    pl1.VData = Bvec(j,2);
    pl1.WData = Bvec(j,3);

    pl2.XData = Bvec(1:j,1);
    pl2.YData = Bvec(1:j,2);
    pl2.ZData = Bvec(1:j,3);
    drawnow
end

%% Power Spectrum check
% psd = (1/(Fs*N)).*ft.*conj(ft);
% psd(2:end-1,:) = 2*psd(2:end-1,:);
% 
% for d = 1 %Plotting...
% figure(4)
% subplot(3,1,1)
% plot(freq,10*log10(psd(:,1)))
% subplot(3,1,2)
% plot(freq,10*log10(psd(:,2)))
% subplot(3,1,3)
% plot(freq,10.*log10(psd(:,3)))
% end