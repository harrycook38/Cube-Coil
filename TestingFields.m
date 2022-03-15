tic
clear all;close all;clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% In 3D
mu = 4*pi*10^(-7); 
I = 1; %Amp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Creating environment
res = 20;
u = linspace(-0.5,0.5,res); %100cm
[x,y,z] = meshgrid(u,u,u); 
plot3(x(:),y(:),z(:),'.k','LineWidth',0.1)
hold on

%% Coil
do2 = 0.2;
a = linspace(-do2,do2,res.*2*do2);


%% Biot-Savart Law
function [Bx,By,Bz] = biot_savart(rs,rc,I,dl)

    sn = length(rs);
    ln = length(dl);
    B = zeros(sn,3);
    
    for i = 1:sn
        dB = [0 0 0];
        for j = 1:ln
            rprim = rs(i,:)-rc(j,:);
            db = (I * mu0/(4*pi)).*(cross(dl(j,:),rprim))./(abs(rprim).^2);
            dB = dB + db;
        end
        B(i,:) = dB;
    end

    bx = B(:, 1);
    by = B(:, 2);
    bz = B(:, 3);

    lb = round(ns.^(1/3));

    Bx = reshape(bx, lb, lb, lb);
    By = reshape(by, lb, lb, lb);
    Bz = reshape(bz, lb, lb, lb);

end