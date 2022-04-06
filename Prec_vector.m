clear all
close all 
clc

%% Preamble

time = 10;
t = linspace(0,time,1000);

numrotxy = 5;  %Number of circles drawn in xy plane per time
numrotz = 0.5; %Number of vertical gain periods

%% xy plane

T1 = max(t)./numrotxy;
f1 = (2*pi)/T1; 

x = sin(f1.*t);
y = cos(f1.*t);

%% z direction

T2 = max(t)./numrotz; %period = double the timespan to go from -1 to 1 
f2 = (2*pi)/T2; 

z = -cos(f2.*t);
figure(1)
plot(t,z)

or = [0 0 0]; %Origin
r2 = [x;y;z;];
%% Plotting

%Angles
theta = acosd(z);
phi = atand(y./x);

figure(10)

subplot(2,2,[1 3])
pl1 = quiver3(or(1),or(2),or(3),r2(1,1),r2(2,1),r2(3,1),"LineWidth",3);
hold on
pl2 = plot3(r2(1,1),r2(2,1),r2(3,1),'r-');
xlim([-1.5 1.5]); ylim([-1.5 1.5]); zlim([-1.5 1.5])
axis square

subplot(2,2,2)
pl3 = plot(t(1),phi(1));
xlim([0 max(t)]); ylim([-95 95]);

subplot(2,2,4)
pl4 = plot(t(1),theta(1));
xlim([0 max(t)]); ylim([0 190]); 

for j = 1:length(t)
    pl1.UData = r2(1,j);
    pl1.VData = r2(2,j);
    pl1.WData = r2(3,j);

    pl2.XData = r2(1,1:j);
    pl2.YData = r2(2,1:j);
    pl2.ZData = r2(3,1:j);

    pl3.XData = t(1:j);
    pl3.YData = phi(1:j);

    pl4.XData = t(1:j);
    pl4.YData = theta(1:j);

    drawnow
end











