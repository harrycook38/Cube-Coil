clear all
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Checking to see if the Z sensor is bad, Turns out that the Z coils are the
%problem, as expected. 
%% Preamble
cd('/Users/Harry/Documents/GitHub/Cube-Coil/Lab Stuff/Data/10-5-22/zy-check')
addpath('/Users/Harry/Documents/GitHub/Cube-Coil/Lab Stuff/Code')
Fs = 1213;              %Sampling Frequency
nfiles= 9;              %Number of files/points measured
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Main

%3x3x3 Array (Y and Z) (Y is in channel 3, Z is in channel 1)
fill = zeros(nfiles,3);
for n = 1:nfiles
        [rmsvals] = extractrms(['QZFM_',num2str(n-1)],Fs);
        fill(n,:) = rmsvals;
end

uniform = reshape(fill,[3,3,3]);

coords1 = [-6 0 6];
coords2 = [-7,-1,5];
coords3 = [-5,1,7];

figure(1)
surf(coords1,coords2,uniform(:,:,3));
hold on
surf(coords1,coords3,uniform(:,:,1));
zlim([0,2.5e-7])
xlabel('xgrid'); ylabel('ygrid'); zlabel('Magnetic Field');


%% Errors
% averg = mean(fill,1);
% c1 = ((fill./averg)-1)*100;
% error = reshape(c1,[8,7,3]);
% 
% figure(2)
% surf(coords1,coords2,abs(error(:,:,1)));
% xlabel('xgrid'); ylabel('ygrid'); zlabel('Percent change from mean');

