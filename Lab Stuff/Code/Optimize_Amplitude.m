clear all
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Preamble
addpath('/Users/Harry/Documents/GitHub/Cube-Coil/Lab Stuff/Code')
cd('/Users/Harry/Documents/GitHub/Cube-Coil/Lab Stuff/Data/11-5-22/Noise/')
%Requires 3 datasets, each with the sensor at the origin.

%% Processing
%Signal from non-corrected pulses

rmsx = extractrms('QZFM_0',1213);
rmsy = extractrms('QZFM_1',1213);
rmsz = extractrms('QZFM_2',1213);

rms_I = [rmsx(1),rmsy(2),rmsz(3)];

%User input Amplitude in voltage
in = [5,5,5];
% in(:,1) = input('Voltage amplitude in X: ');
% in(:,2) = input('Voltage amplitude in Y: ');
% in(:,3) = input('Voltage amplitude in Z: ');


%% Find new Voltages

rel = rms_I./(rms_I(1)); %Equalisation to x field strength





