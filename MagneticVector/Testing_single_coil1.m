clear all
close all
clc

%Environment to calculate B field on
di = 0.5; 
np = 15;
d = linspace(-0.5*di,0.5*di,np);
[x,y,z] = meshgrid(d);
pp = [x(:),y(:),z(:)];

%Coil diameters
d1sep = 0.3; % Different coil sizes to account for physical standoff


%Writing Pair 1
l1 = linspace(-0.5*d1sep,0.5*d1sep,length(d).*2);
orgsep = (0.5445./2); %Optimal square coil separation for now %%% increased to put more field points between coils 
loc = [-orgsep orgsep];

for k = 1:2
    pos = loc(k).*ones(1,length(l1));
    X1(:,:,k) = [l1;pos;l1(1).*ones(1,length(l1))];
    X2(:,:,k) = [l1(end).*ones(1,length(l1));pos;l1];
    X3(:,:,k) = [-1.*l1;pos;l1(end).*ones(1,length(l1))];
    X4(:,:,k) = [l1(1).*ones(1,length(l1));pos;-1.*l1];
    X(:,:,k) =  [X1(:,:,k) X2(:,:,k) X3(:,:,k) X4(:,:,k)];
end
X = reshape(X,3,[])';

dlX = diff(X); %Line element vectors
dlX(end+1,:) = dlX(end,:);
dlX(end./2,:) = dlX((end./2)-1,:);  

%TESTING

figure(500)
plot3(X(:,1),X(:,2),X(:,3),'bo');
hold on
xlabel('x')
ylabel('y')

quiver3(X(:,1),X(:,2),X(:,3),dlX(:,1),dlX(:,2),dlX(:,3),3) %Check dls are correct
axis equal


%% Clearing unessential variables for next bit
clearvars -except X dlX x y z pp

%% Biot-Savart

mu = pi*4e-7; 
I = 1; %Amps
const = I * mu/(4*pi);
Bout = zeros(length(pp),3); %Magnetic field map

tic
for m = 1:length(pp) %Number of points in space
    for k = 1:length(X) %Number of wire elements
 
            rprim = pp(m,:) - X(k,:);  %%% the loop does this element by element 
         
            dB = const.*(cross(dlX(k,:),rprim)./(abs(rprim).^3)); %% rprim is a single vector

            Bout(m,:) = Bout(m,:)+dB; %%%% need to sum the contributions from all the elements
     end   

            if mod(m,1000) == 0
                fprintf('Iteration %d...\n',m);
            end
end

toc

figure(1)
quiver3(pp(:,1),pp(:,2),pp(:,3),Bout(:,1),Bout(:,2),Bout(:,3),4,'c');
hold on
plot3(X(:,1),X(:,2),X(:,3),'or')
axis equal

xlabel('x')
ylabel('y')
zlabel('z')
