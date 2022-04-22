clear all
close all 
clc

%% Preamble

time = 10;
t = linspace(0,time,1000);

numrotxy = 50;  %Number of circles drawn in xy plane per time period
numrotz = 0.5; %Number of vertical gain periods
% numamp = 1; %Period of xy amplitude change

%% xy plane
T1 = max(t)./numrotxy;
f1 = (2*pi)/T1; 

% z direction
T2 = max(t)./numrotz;
f2 = (2*pi)/T2; 

% Tr = max(t)./numamp;
% fr = (2*pi)./Tr;
% r = (sin(fr.*t + ((3*pi)./2))+1)./2;
% 
% figure(1)
% plot(t,r)
% 
% x = r.*sin(f1.*t);
% y = r.*cos(f1.*t);
% z = -cos(f2.*t)

u = f1.*t;
v = f2.*t;
%TEST
r = 1;
x1 = cos(u).*sin(v);
y1 = sin(u).*sin(v);
z1 = cos(v);

%Rewriting as sums of sins
x = 0.5.*(sin(v+u)+sin(v-u));
y = 0.5.*(sin((pi/2) -(u-v))-sin((pi/2)-(u+v)));
z = sin((pi/2)-v);
or = [0 0 0]; %Origin
r2 = [x;y;z;];

figure(1)
subplot(3,1,1)
plot(t,x); xlabel('Time'); ylabel('X Amplitude');
subplot(3,1,2)
plot(t,y); xlabel('Time'); ylabel('Y Amplitude');
subplot(3,1,3)
plot(t,z); xlabel('Time'); ylabel('Z Amplitude');

%% Plotting
%Angles
theta = acosd(z);
phi = atand(y./x);
figure(2)
plot3(r2(1,:),r2(2,:),r2(3,:),'r-');
xlim([-1.5 1.5]); ylim([-1.5 1.5]); zlim([-1.5 1.5])
xlabel("X"); ylabel("Y"); zlabel("Z");
axis square; grid on;

figure(10) %Drawing Vectors and show angles through time

subplot(2,2,[1 3])
pl1 = quiver3(or(1),or(2),or(3),r2(1,1),r2(2,1),r2(3,1),"LineWidth",3);
hold on
pl2 = plot3(r2(1,1),r2(2,1),r2(3,1),'r-');
xlim([-1.5 1.5]); ylim([-1.5 1.5]); zlim([-1.5 1.5])
xlabel("X",FontSize=15); ylabel("Y",FontSize=15); zlabel("Z",FontSize=15);
axis square
title("Field Output",FontSize=20)

subplot(2,2,2)
pl3 = plot(t(1),phi(1));
xlim([0 max(t)]); ylim([-95 95]);
title("Azimuthal Angle",FontSize=20)
xlabel("Time (s)",FontSize=15); ylabel("Phi",FontSize=15)

subplot(2,2,4)
pl4 = plot(t(1),theta(1));
xlim([0 max(t)]); ylim([0 190]); 
title("Polar Angle",FontSize=20)
xlabel("Time (s)",FontSize=15); ylabel("Theta",FontSize=15)

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











