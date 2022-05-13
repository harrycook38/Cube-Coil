function [rmsvals] = extractrms(fname,Fs)

addpath('/Users/Harry/Documents/GitHub/Cube-Coil/Lab Stuff/Code')
[~,rdata] = eb_read_lvm(fname);

fgmout = rdata(:,1:3);
Bout =(1e-6*fgmout)./0.1;       
Bcor = Bout - mean(Bout,1);
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Filter & Inverse ft
sd = 0.5;                               %Width of filter
a = (1./sqrt(2.*pi.*sd));               %Normalised Amplitude
gau = a.*exp(-(freq-mu).^8/(2*sd^2));   %Modified Gaussian

filt = gau'.*ft;                        %Filtering

inv = ifft(filt,'symmetric');           %Inverse Transform

xrms = rms(inv(:,1));
yrms = rms(inv(:,2));
zrms = rms(inv(:,3));
rmsvals = [xrms,yrms,zrms];
