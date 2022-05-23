clear all
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Preamble
cd('/Users/Harry/Documents/GitHub/Cube-Coil/Lab Stuff/Data/10-5-22/zy-uniformity')
%cd('/Users/Harry/Documents/GitHub/Cube-Coil/Lab Stuff/Data/6-5-22')
addpath('/Users/Harry/Documents/GitHub/Cube-Coil/Lab Stuff/Code')
Fs = 1213;              %Sampling Frequency
nfiles= 56;             %Number of files/points measured
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Main

%8x7x3 Array (Y and Z) (Y is in channel 1)
fill = zeros(nfiles,3);
for n = 1:nfiles
        [rmsvals] = extractrms(['QZFM_',num2str(n-1)],Fs);
        fill(n,:) = rmsvals;
end

uniform1 = reshape(fill,[8,7,3]);

%7x7x3 Array (X)
cd('/Users/Harry/Documents/GitHub/Cube-Coil/Lab Stuff/Data/9-5-22/x_uniformity/')
Fs = 1213;              %Sampling Frequency
nfiles= 49;   
fill2 = zeros(nfiles,3);
for n = 1:nfiles
        [rmsvals] = extractrms(['QZFM_',num2str(n-1)],Fs);
        fill2(n,:) = rmsvals;
end

uniform2 = reshape(fill2,[7,7,3]);

surfaces = zeros(7,7,3);
surfaces(:,:,1) = uniform2(:,:,1);
surfaces(:,:,2) = uniform1(1:7,:,1);
surfaces(:,:,3) = uniform1(2:8,:,3);


coords = [-6 -4 -2 0 2 4 6];
figure(1)
surf(coords,coords,surfaces(:,:,1))
hold on
surf(coords,coords,surfaces(:,:,2))
surf(coords,coords,surfaces(:,:,3))
colorbar
xlabel('xgrid'); ylabel('ygrid'); zlabel('Magnetic Field');


%% Errors

avg = squeeze(mean(mean(surfaces)));
c = zeros(7,7,3);

figure(2)
hold on
for i = 1:3
c(:,:,i) = ((surfaces(:,:,i)./avg(i))-1).*100;
posc = abs(c);
surf(coords,coords,posc(:,:,i))
end 
xlabel('xgrid'); ylabel('ygrid'); zlabel('Percent Change');

%Avg Percent Change

perc = squeeze(mean(mean(posc)));

xyperc = mean(perc(1:2));

%clearvars -except surfaces posc

xs = surfaces(:,:,1);
ys = surfaces(:,:,2);
zs = surfaces(:,:,3);



