clear all
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Preamble
addpath('/Users/Harry/Documents/GitHub/Cube-Coil/Lab Stuff/Code')
addpath('/Users/Harry/Documents/GitHub/Cube-Coil/Lab Stuff/Data/11-5-22/Noise/')
%Requires 3 datasets, each with the sensor at the origin.

%% Acquire Data

[~,rdata1] = eb_read_lvm('QZFM_0');
Vout = zeros(length(rdata1),3);
Vout(:,1) = rdata1(:,1);   

[~,rdata2] = eb_read_lvm('QZFM_1');
Vout(:,2) = rdata2(:,2);

[~,rdata3] = eb_read_lvm('QZFM_2');
Vout(:,3) = rdata3(:,3);

bl_V = mean(Vout);

Bfield = 1e-6*(bl_V)./0.1;

