clear all
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Preamble
cd('/Users/Harry/Documents/GitHub/Cube-Coil/Lab Stuff/Data/z_test_dc/')
addpath('/Users/Harry/Documents/GitHub/Cube-Coil/Lab Stuff/Code')
Fs = 1213;              %Sampling Frequency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Main

nfiles = 9;
fill = zeros(nfiles,3);
for n = 1:nfiles
    [~,rdata] = eb_read_lvm(['QZFM_',num2str(n-1)]);
    fgmout = rdata(:,1:3);
    Bout = 0.1e-6*fgmout;
    meanB = mean(Bout);
    fill(n,:) = meanB;
end

uniform = reshape(fill,[sqrt(nfiles),sqrt(nfiles),3]);
coords = [-6 0 6];

figure
surf(coords,coords,uniform(:,:,1));
xlabel('x')