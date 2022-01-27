clear all
close all
clc


%% 2D
di = 1; 

d = linspace(-0.5*di,0.5*di,20);

[x,y] = meshgrid(d);

l1 = linspace(-0.2,0.2,length(d).*1.2);

X1 = [l1;l1(1).*ones(1,length(l1))];
X2 = [l1(1).*ones(1,length(l1));l1];

X3 = [l1(end).*ones(1,length(l1));l1];
X4 = [l1;l1(end).*ones(1,length(l1))];

X = [X1 X2 X3 X4];

plot(x,y,'.k')
hold on
plot(X(1,:),X(2,:),'or')

%% 3D
di = 1; 
d = linspace(-0.5*di,0.5*di,10);
[x,y,z] = meshgrid(d);

%Pair 1: Fixed along x with 0.4m diameter
l1 = linspace(-0.2,0.2,length(d).*2);
orgsep = (0.5445./2)*0.4;
loc = [-orgsep orgsep];

% Fixed Along X
for i = 1:2
    pos = loc(i).*ones(1,length(l1));
    X1(:,:,i) = [l1;pos;l1(1).*ones(1,length(l1))];
    X2(:,:,i) = [l1(1).*ones(1,length(l1));pos;l1];
    X3(:,:,i) = [l1(end).*ones(1,length(l1));pos;l1];
    X4(:,:,i) = [l1;pos;l1(end).*ones(1,length(l1))];
    X(:,:,i) = [X1(:,:,i) X2(:,:,i) X3(:,:,i) X4(:,:,i)];
end

plot3(x(:),y(:),z(:),'.k')
hold on
plot3(X(1,:,1),X(2,:,1),X(3,:,1),'or')
plot3(X(1,:,2),X(2,:,2),X(3,:,2),'or')
axis square

%Half as big as previous squares
l2 = linspace(-0.1,0.1,length(d).*2);
orgsep = (0.5445./2)*0.4;
loc = [-orgsep orgsep];
% Fixed Along Y
for i = 1:2
    pos = loc(i).*ones(1,length(l2));
    Y1(:,:,i) = [pos;l2;l2(1).*ones(1,length(l2))];
    Y2(:,:,i) = [pos;l2(1).*ones(1,length(l2));l2];
    Y3(:,:,i) = [pos;l2(end).*ones(1,length(l2));l2];
    Y4(:,:,i) = [pos;l2;l2(end).*ones(1,length(l2))];
    Y(:,:,i) = [Y1(:,:,i) Y2(:,:,i) Y3(:,:,i) Y4(:,:,i)];
end

plot3(Y(1,:,1),Y(2,:,1),Y(3,:,1),'or')
plot3(Y(1,:,2),Y(2,:,2),Y(3,:,2),'or')





