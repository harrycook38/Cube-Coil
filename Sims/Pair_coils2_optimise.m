clear all
close all
clc

%COIL PAIR 1
%Environment to calculate B field on
di1 = 1;
np = 31;
d1 = linspace(-0.5*di1,0.5*di1,np);
[x1,y1,z1] = meshgrid(d1);

L = 0.3; % Separation: constant for all coils
loc = [-0.5*L 0.5*L];

%Coil diameters
diam1 = 0.3;

[XYZ1, dl1] = writecoily(diam1,loc,np,2,1);

%TESTING
figure(1)
plot3(XYZ1(:,1),XYZ1(:,2),XYZ1(:,3),'bo');
hold on
xlabel('x')
ylabel('y')

quiver3(XYZ1(:,1),XYZ1(:,2),XYZ1(:,3),dl1(:,1),dl1(:,2),dl1(:,3),3)
axis equal

%% Clearing unessential variables for next bit
clearvars -except XYZ1 dl1 loc np

%% Biot-Savart along central coil axis
mu = pi*4e-7; 
I = 1; %Amps
const = I * mu/(4*pi);

%Consider inner 2/3 of the coils as the area to optimise
q = linspace((2/3).*loc(1),(2/3).*loc(2),100); 

v = [zeros(size(q)); q; zeros(size(q))]';
Bin = zeros(length(v),3); %Magnetic field map

tic
%Field along symmetric axis
for m = 1:length(v) %Number of points in space
    for k = 1:length(XYZ1) %Number of wire elements
 
            rprim = v(m,:) - XYZ1(k,:); 
            dB = const.*(cross(dl1(k,:),rprim)./(abs(rprim).^3));
            Bin(m,:) = Bin(m,:) + dB;

    end
end

Bin(isnan(Bin))=0; toc

y = v(:,2);
By = Bin(:,2);

%% Calculating Initial polynomial

fits = polyfit(y,By,4);
fits([2,4]) = 0; 

val = polyval(fits,y);

figure(2)
plot(y,By)
hold on
plot(y,val,'LineWidth',3)
xlabel('y axis')
ylabel('Field Strength')

%% Optimising inner loop

dinit = 0.25; %2/3 original loop size
[Y2,dl2] = writecoily(dinit,loc,np,2,2);

figure(400)
plot3([XYZ1(:,1); Y2(:,1)],[XYZ1(:,2); Y2(:,2)],[XYZ1(:,3); Y2(:,3)],'bo');
hold on
xlabel('x')
ylabel('y')

quiver3([XYZ1(:,1); Y2(:,1)],[XYZ1(:,2); Y2(:,2)],[XYZ1(:,3); Y2(:,3)],...
    [dl1(:,1); dl2(:,1)],[dl1(:,2); dl2(:,2)],[dl1(:,3); dl2(:,3)],3)
axis equal

%Finding optimum diameter for fixed current
diam = dinit; 
loops = 100;
citer = 100;
Cmat = zeros(loops,3);
I = linspace(1,0.01,100);
const = I * mu/(4*pi);

Bi = zeros(length(v),loops);

Cmat = zeros(loops,citer,3);

tic
for iter1 = 1:citer

consti = const(citer);
    tic
    for iter2 = 1:loops
    
        [XYZi,dli] = writecoily(diam,loc,np,2,2);
        
        coords = [XYZ1; XYZi];
        dls = [dl1;dli];
    
        Bin = zeros(length(v),3);
    
        for m = 1:length(v) %Number of points in space
            for k = 1:length(coords) %Number of wire elements
        
                    rprim = v(m,:) - coords(k,:); 
                    dB = consti.*(cross(dls(k,:),rprim)./(abs(rprim).^3));
                    Bin(m,:) = Bin(m,:) + dB;
        
            end
        end
        
        Bin(isnan(Bin))=0;
        Bi(:,iter2) = Bin(:,2);
    
        ifit = polyfit(y,Bi(:,iter2),2);
        ifit(2) = 0;
       
        C(iter2,:) = ifit;
        diam = diam-iter2.*0.001;
    
        if mod(iter2,10) == 0
                   fprintf('Length Iteration %d...\n',iter2);
        end
    
    end
    toc
    if mod(iter1,10) == 0
               fprintf('Current Iteration %d...\n',iter1);
    end

    Cmat(iter1,:,:) = C;

end
toc

% figure(500)
% plot(1:length(C),C(:,3),'r')
% hold on 
% plot(1:length(C),C(:,5),'b')
% 
% figure(600)
% plot(y,By,'r')
% hold on
% plot(y,Bi(:,200),'b')







