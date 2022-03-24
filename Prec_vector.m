clear all
close all 
clc

%% 2D

theta = linspace(0,10*pi,1000);

x = sin(theta);
y = cos(theta);

% figure(1)
% plot(theta,x)
% hold on
% plot(theta,y)

% figure(2)

% r = [x;y];
% plot(r(1,:),r(2,:),'k*','LineWidth',10);
% hold on
% plot(r(1,2),r(2,2),'k*','LineWidth',10);
% axis square

% figure(3)
% pl = quiver(0,0,r(1,1),r(2,1));
% xlim([-1.5 1.5]); ylim([-1.5 1.5]);
% axis square
% 
% for i = 1:length(theta)
%     pl.UData = r(1,i);
%     pl.VData = r(2,i);
%     drawnow
% end

%% 3D

z = sin(theta./2); % Height of vector

or = [0 0 0]; %Origin

r2 = [x;y;z;];
figure(10)
pl2 = quiver3(or(1),or(2),or(3),r2(1,1),r2(2,1),r2(3,1),"LineWidth",3);
hold on
pl3 = plot3(r2(1,1),r2(2,1),r2(3,1),'r-');
xlim([-1.5 1.5]); ylim([-1.5 1.5]); zlim([-1.5 1.5])
axis square

for j = 1:length(theta)
    pl2.UData = r2(1,j);
    pl2.VData = r2(2,j);
    pl2.WData = r2(3,j);
    pl3.XData = r2(1,1:j);
    pl3.YData = r2(2,1:j);
    pl3.ZData = r2(3,1:j);
    drawnow
end














