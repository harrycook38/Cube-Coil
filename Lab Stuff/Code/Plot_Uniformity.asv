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

uniform = reshape(fill,[8,7,3]);

coords1 = [-6 -4 -2 0 2 4 6];
coords2 = [-8 -6 -4 -2 0 2 4 6];
coords3 = [-6 -4 -2 0 2 4 6 8];
figure(1)
surf(coords1,coords2,uniform(:,:,1));
hold on
%surf(coords1,coords2,uniform(:,:,2));
surf(coords1,coords3,uniform(:,:,3));
zlim([0,2.5e-7])
xlabel('xgrid'); ylabel('ygrid'); zlabel('Magnetic Field');


%7x7x3 Array (X)
cd('/Users/Harry/Documents/GitHub/Cube-Coil/Lab Stuff/Data/9-5-22/x_uniformity/')
Fs = 1213;              %Sampling Frequency
nfiles= 49;   
fill2 = zeros(nfiles,3);
for n = 1:nfiles
        [rmsvals] = extractrms(['QZFM_',num2str(n-1)],Fs);
        fill2(n,:) = rmsvals;
end

uniform1 = reshape(fill2,[7,7,3]);

coords1 = [-6 -4 -2 0 2 4 6];
surf(coords1,coords1,uniform1(:,:,1))

%% Errors
averg1 = mean(fill,1);
c1 = ((fill./averg1)-1)*100;
error1 = reshape(c1,[8,7,3]);

averg2 = mean(fill,1);
c2 = ((fill2./averg2)-1)*100;
error2 = reshape(c2,[7,7,3]);

figure(2)
surf(coords1,coords2,abs(error1(:,:,1)));
xlabel('xgrid'); ylabel('ygrid'); zlabel('Percent change from mean');

%Percentage Change
pc1 = mean(c1);
pc2 = mean(c2);




