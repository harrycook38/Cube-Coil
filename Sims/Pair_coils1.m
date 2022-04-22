clear all
close all
clc

%Environment to calculate B field on
di = 1;
np = 31;
d = linspace(-0.5*di,0.5*di,np);
[x,y,z] = meshgrid(d);
pp = [x(:),y(:),z(:)];

%Coil diameters
d1sep = 0.3; % Different coil sizes to account for physical standoff

%Pair 1
l1 = linspace(-0.5*d1sep,0.5*d1sep,length(d).*5);
orgsep = 0.5.*d1sep; %.*(0.5445./2) %Optimal square coil separation. 
loc = [-orgsep orgsep];

for k = 1:2
    pos = loc(k).*ones(1,length(l1));
    X1(:,:,k) = [l1;pos;l1(1).*ones(1,length(l1))];
    X2(:,:,k) = [l1(end).*ones(1,length(l1));pos;l1];
    X3(:,:,k) = [-1.*l1;pos;l1(end).*ones(1,length(l1))];
    X4(:,:,k) = [l1(1).*ones(1,length(l1));pos;-1.*l1];
    X(:,:,k) = [X1(:,:,k) X2(:,:,k) X3(:,:,k) X4(:,:,k)];
end

X = reshape(X,3,[])';
XYZ = X;

%Line element vectors
dlX = diff(X);
dlX(end+1,:) = dlX(end,:);
dlX(end./2,:) = dlX((end./2)-1,:);

dl = dlX;

%TESTING
figure(500)
plot3(X(:,1),X(:,2),X(:,3),'bo');
hold on
xlabel('x')
ylabel('y')

quiver3(XYZ(:,1),XYZ(:,2),XYZ(:,3),dl(:,1),dl(:,2),dl(:,3),3)
axis equal

%% Clearing unessential variables for next bit
clearvars -except XYZ dl pp

%% Biot-Savart
mu = pi*4e-7; 
I = 1; %Amps
const = I * mu/(4*pi);
Bout = zeros(length(pp),3); %Magnetic field map

tic
for m = 1:length(pp) %Number of points in space
    for k = 1:length(XYZ) %Number of wire elements
            rprim = pp(m,:) - XYZ(k,:);
            dB = const.*(cross(dl(k,:),rprim)./(abs(rprim).^3));
            Bout(m,:) = Bout(m,:) + dB;
    end

            if mod(m,1000) == 0
               fprintf('Iteration %d...\n',m);
            end
end
toc
Bout(isnan(Bout)) = 0; %Setting NaNs to 0
figure(1)
quiver3(pp(:,1),pp(:,2),pp(:,3),Bout(:,1),Bout(:,2),Bout(:,3),4,'c');
hold on
plot3(XYZ(:,1),XYZ(:,2),XYZ(:,3),'or')
axis equal

xlabel('x')
ylabel('y')
zlabel('z')


%% Inside all coils 

q = linspace(-0.15,0.15,31);
[u,v,w] = meshgrid(q);

in = [u(:),v(:),w(:)]; clearvars u v w; 

mu = pi*4e-7; 
I = 10; %Amps
Bin = zeros(length(in),3); %Magnetic field map

tic
%FIELD WITHIN ALL COILS
for m = 1:length(in) %Number of points in space
    for k = 1:length(XYZ) %Number of wire elements
 
            rprim = in(m,:) - XYZ(k,:); 
            dB = const.*(cross(dl(k,:),rprim)./(abs(rprim).^3)); 
            Bin(m,:) = Bin(m,:) + dB;

    end
            if mod(m,1000) == 0
               fprintf('Iteration %d...\n',m);
            end
end
toc

Bin(isnan(Bin))=0;
figure(1)
quiver3(in(:,1),in(:,2),in(:,3),Bin(:,1),Bin(:,2),Bin(:,3),10,'r');
hold on
plot3(XYZ(:,1),XYZ(:,2),XYZ(:,3),'or')
axis equal

xlabel('x')
ylabel('y')
zlabel('z')

%% Taking cross-section

zerox = find(in(:,1)==0);
zeroz = find(in(:,3)==0);

onaxis = intersect(zerox,zeroz);

axisp = in(onaxis,:);
Baxis = Bin(onaxis,:);

figure(3)

plot3(axisp(:,1),axisp(:,2),axisp(:,3),'or')


quiver3(axisp(:,1),axisp(:,2),axisp(:,3),Baxis(:,1),Baxis(:,2),Baxis(:,3),'r');

plot(axisp,Baxis)




















