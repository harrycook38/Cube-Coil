clear all
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Preamble
addpath('/Users/Harry/Documents/GitHub/Cube-Coil/Lab Stuff/Code')
cd('/Users/Harry/Documents/GitHub/Cube-Coil/Lab Stuff/Data/11-5-22/Noise')
%Requires 3 datasets, each with the sensor at the origin.

for i = 1:3
[~,rdata] = eb_read_lvm(['QZFM_', num2str(i-1)]);
fgmout = rdata(:,1:3);
Bout =(1e-6*fgmout)./0.1;       
dat(:,:,i) = Bout;
end 

means = zeros(3,3);
for j = 1:3
means(:,j) = mean(dat(:,:,j));
end

