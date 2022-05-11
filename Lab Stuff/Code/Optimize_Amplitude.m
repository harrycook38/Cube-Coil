clear all
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Preamble
addpath('/Users/Harry/Documents/GitHub/Cube-Coil/Lab Stuff/Code')
addpath('/Users/Harry/Documents/GitHub/Cube-Coil/Lab Stuff/Data/11-5-22/Origin_Test/')
%Requires 3 datasets, each with the sensor at the origin.

%% Processing
%Signal from non-corrected pulses
rms_I = [0,0,0];
for i = 1:3
    [~,rdata] = eb_read_lvm(['QZFM_',num2str(i-1)]);
    Bout =(1e-6*rdata)./0.1;       
    if i == 1
        rms_I(:,1) = rms(Bout(:,1));
    elseif i == 2
        rms_I(:,2) = rms(Bout(:,2));
    elseif i == 3
        rms_I(:,3) = rms(Bout(:,3));
    end
end

%User input Amplitude in voltage
in = [0,0,0];
in(:,1) = input('Voltage amplitude in X: ');
in(:,2) = input('Voltage amplitude in Y: ');
in(:,3) = input('Voltage amplitude in Z: ');

target_B = 2e-9;

%% Find new Voltages

rel = rms_I./in;





