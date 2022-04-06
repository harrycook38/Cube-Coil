clear all
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Preamble
addpath('Code/')
addpath('Data/5-4-22')

fname = 'QZFM_6.lvm';   %Filename
Fs = 1213;              %Sampling Frequency
mu = 10;                %Target frequency to filter around
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Main
%Load in time stamps and data
[rtime,rdata] = eb_read_lvm(fname);
fgmout = rdata(:,1:3);      %FGM Data
Iout = rdata(:,4:6)./10;    %Current out from loops
figure(1)
plot(rtime,Iout)

Bout = 0.1e-6*fgmout;       %Field (Conversion 0.1 uT/V)
Bcor = Bout - mean(Bout,1); %Correct for Vertical Offset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fourier
N = length(Bcor);
fti = fft(Bcor);
ft = fti(1:N/2+1,:);
freq = 0:Fs/N:Fs/2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Filter & Inverse ft
sd = 0.5;                               %Width of filter
a = (1./sqrt(2.*pi.*sd));               %Normalised Amplitude
gau = a.*exp(-(freq-mu).^8/(2*sd^2));   %Modified Gaussian
filt = gau'.*ft;                        %Filtering
inv = ifft(filt,'symmetric');           %Inverse Transform

for d = 1 %Plotting...
figure(2)
subplot(3,1,1)
plot(rtime(1:2:end),inv(:,1))
subplot(3,1,2)
plot(rtime(1:2:end),inv(:,2))
subplot(3,1,3)
plot(rtime(1:2:end),inv(:,3))
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RMS Output
xrms = rms(inv(:,1));
yrms = rms(inv(:,2));
zrms = rms(inv(:,3));
%% Power Spectrum check
psd = (1/(Fs*N)).*ft.*conj(ft);
psd(2:end-1,:) = 2*psd(2:end-1,:);

for d = 1 %Plotting...
figure(3)
subplot(3,1,1)
plot(freq,10*log10(psd(:,1)))
subplot(3,1,2)
plot(freq,10*log10(psd(:,2)))
subplot(3,1,3)
plot(freq,10.*log10(psd(:,3)))
end