clear all
close all
clc


%% 3D

%Environment to calculate B field on
di = 1; 
d = linspace(-0.5*di,0.5*di,15);
isout= abs(d) >= 0.1;
d = d.*isout;
[x,y,z] = meshgrid(d);
points = [x(:),y(:),z(:)];

%Coil diameters
d1sep = 0.3;
d2sep = 0.4;
d3sep = 0.5;

%Pair 1
l1 = linspace(-0.5*d1sep,0.5*d1sep,length(d).*2);
orgsep = (0.5445./2)*d1sep;
loc = [-orgsep orgsep];

for k = 1:2
    pos = loc(k).*ones(1,length(l1));
    X1(:,:,k) = [l1;pos;l1(1).*ones(1,length(l1))];
    X2(:,:,k) = [l1(1).*ones(1,length(l1));pos;l1];
    X3(:,:,k) = [l1(end).*ones(1,length(l1));pos;l1];
    X4(:,:,k) = [l1;pos;l1(end).*ones(1,length(l1))];
    X(:,:,k) = [X1(:,:,k) X2(:,:,k) X3(:,:,k) X4(:,:,k)];
end

%Pair 2
l2 = linspace(-0.5*d2sep,0.5*d2sep,length(d).*2);
orgsep = (0.5445./2)*d2sep;
loc = [-orgsep orgsep];

for k = 1:2
    pos = loc(k).*ones(1,length(l2));
    Y1(:,:,k) = [pos;l2;l2(1).*ones(1,length(l2))];
    Y2(:,:,k) = [pos;l2(1).*ones(1,length(l2));l2];
    Y3(:,:,k) = [pos;l2(end).*ones(1,length(l2));l2];
    Y4(:,:,k) = [pos;l2;l2(end).*ones(1,length(l2))];
    Y(:,:,k) = [Y1(:,:,k) Y2(:,:,k) Y3(:,:,k) Y4(:,:,k)];
end

%Pair 3
l3 = linspace(-0.5*d3sep,0.5*d3sep,length(d).*2);
orgsep = (0.5445./2)*d3sep;
loc = [-orgsep orgsep];

for k = 1:2
    pos = loc(k).*ones(1,length(l3));
    Z1(:,:,k) = [l3;l3(1).*ones(1,length(l3));pos];
    Z2(:,:,k) = [l3(1).*ones(1,length(l3));l3;pos];
    Z3(:,:,k) = [l3(end).*ones(1,length(l3));l3;pos];
    Z4(:,:,k) = [l3;l3(end).*ones(1,length(l3));pos];
    Z(:,:,k) = [Z1(:,:,k) Z2(:,:,k) Z3(:,:,k) Z4(:,:,k)];
end

%Combine & Plot
XYZ = [X,Y,Z];
XYZ = reshape(XYZ,3,[])';

%plot3(XYZ(:,1),XYZ(:,2),XYZ(:,3),'or')
%axis equal

%% Clearing for next bit
clearvars -except XYZ points

%% Biot Savart

%Inner environment
q = linspace(-0.1,0.1,30);
[u,v,w] = meshgrid(q);

inner = [u(:),v(:),w(:)]; clearvars u v w;

mu = pi*4e-7; 
I = 1; %Amps
Bout = zeros(length(points),3); %Magnetic field map
Bin = zeros(length(inner),3); 

tic

% FIELD OUTSIDE THE UNIFORM AREA
for m = 1:length(points(:,1))

            p = points(m,:); % Field points 
            rprim = p - XYZ;
            
            dB = [0 0 0];
            for k = 1:length(XYZ) %Field from source
                db = (I * mu/(4*pi)).*(cross(XYZ(k,:),rprim(k,:)))./(abs(rprim(k,:)).^2);%Wrong
                dB = dB + db;
            end

            Bout(m,:) = dB;
  
            if mod(m,1000) == 0
              fprintf('At iteration %d...\n',m);
            end
end

%FIELD WITHIN ALL COILS
for m = 1:length(inner(:,1))

            p = inner(m,:); % Field points 
            rprim = p - XYZ;

            
            dB = [0 0 0];
            for k = 1:length(XYZ) %Field from source
                db = (I * mu/(4*pi)).*(cross(XYZ(k,:),rprim(k,:)))./(abs(rprim(k,:)).^2);
                dB = dB + db;
            end

            Bin(m,:) = dB;
  
            if mod(m,1000) == 0
              fprintf('At iteration %d...\n',m);
            end
end
toc

figure(2)
quiver3(inner(:,1),inner(:,2),inner(:,3),Bin(:,1),Bin(:,2),Bin(:,3),10,'b');
hold on
quiver3(points(:,1),points(:,2),points(:,3),Bout(:,1),Bout(:,2),Bout(:,3),10,'c');
plot3(XYZ(:,1),XYZ(:,2),XYZ(:,3),'or')
axis equal


%Look at centre of the coils (5cm^2)



