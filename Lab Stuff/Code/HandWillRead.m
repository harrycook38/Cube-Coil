clear all
close all
clc

addpath('Code/')
addpath('Data/5-4-22')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
fname = 'QZFM_6.lvm';

%Load in time stamps and data
[rtime,rdata] = eb_read_lvm(fname);

fgmout = rdata(:,1:3); %FGM Data
Iout = rdata(:,4:6)./10; %Current out from loops

figure(1)
subplot(2,1,1)
Bout = 0.1e-6*fgmout; %Field (Conversion 0.1 uT/V)
Bcor = Bout - mean(Bout,1); %Correct for Vertical Offset
plot(rtime,Bcor)
title("offset-corrected Magnetic Field")

subplot(2,1,2)
plot(rtime,Iout)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
Fs = 1213;
N = length(Bcor);

subplot(3,1,1)
xfti = fft(Bcor(:,1));
xft = xfti(1:N/2+1);
psdx = (1/(Fs*N)).*xft.*conj(xft);
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:Fs/N:Fs/2;
plot(freq,10*log10(psdx))

subplot(3,1,2)
yfti = fft(Bcor(:,2));
yft = yfti(1:N/2+1);
psdy = (1/(Fs*N)).*yft.*conj(yft);
psdy(2:end-1) = 2*psdy(2:end-1);
plot(freq,10*log10(psdy))

subplot(3,1,3)
zfti = fft(Bcor(:,3));
zft = zfti(1:N/2+1);
psdz = (1/(Fs*N)).*zft.*conj(zft);
psdz(2:end-1) = 2*psdz(2:end-1);
plot(freq,10.*log10(psdz))

%Filter 
mu = 10; %Target frequency
sd = 0.5;

a = [max(psdx); max(psdy); max(psdz);];
b = (1./sqrt(2.*pi.*sd));

gau = b.*exp(-(freq-mu).^8/(2*sd^2));

filtx = gau'.*xft;
invx = ifft(filtx,'symmetric');

filty = gau'.*yft;
invy = ifft(filty,'symmetric');

filtz = gau'.*zft;
invz = ifft(filtz,'symmetric');

figure(3)
subplot(3,1,1)
plot(rtime(1:2:end),invx)

subplot(3,1,2)
plot(rtime(1:2:end),invy)

subplot(3,1,3)
plot(rtime(1:2:end),invz)












