clear all
close all
clc


%% 2D
di = 1; 

d = linspace(-0.5*di,0.5*di,20);

[x,y] = meshgrid(d);

l = linspace(-0.2,0.2,length(d).*1.2);

L1 = [l;l(1).*ones(1,length(l))];
L2 = [l(1).*ones(1,length(l));l];

L3 = [l(end).*ones(1,length(l));l];
L4 = [l;l(end).*ones(1,length(l))];

S = [L1 L2 L3 L4];

plot(x,y,'.k')
hold on
plot(S(1,:),S(2,:),'or')

%% 3D
di = 1; 
d = linspace(-0.5*di,0.5*di,10);
[x,y,z] = meshgrid(d);

l = linspace(-0.2,0.2,length(d).*2);
orgsep = (0.5445./2)*0.4;
loc = [-orgsep orgsep];

for i = 1:2
    pos = loc(i).*ones(1,length(l));
    L1(:,:,i) = [l;pos;l(1).*ones(1,length(l))];
    L2(:,:,i) = [l(1).*ones(1,length(l));pos;l];
    L3(:,:,i) = [l(end).*ones(1,length(l));pos;l];
    L4(:,:,i) = [l;pos;l(end).*ones(1,length(l))];
    S(:,:,i) = [L1(:,:,i) L2(:,:,i) L3(:,:,i) L4(:,:,i)];
end

plot3(x(:),y(:),z(:),'.k')
hold on
plot3(S(1,:,1),S(2,:,1),S(3,:,1),'or')
plot3(S(1,:,2),S(2,:,2),S(3,:,2),'or')
axis square

