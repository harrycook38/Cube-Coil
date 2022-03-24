clear all
close all
clc
%% Initial
%Environment to calculate B field on
di1 = 1;
np = 31;
d1 = linspace(-0.5*di1,0.5*di1,np);
[x1,y1,z1] = meshgrid(d1);

L = 0.4; % Separation: constant for all coils
loc = [-0.5*L 0.5*L]; %Coil locations

%Coil diameters
diam1 = 0.4; 

[XYZ1, dl1] = writecoily(diam1,loc,np,2,1);

%% Clearing unessential variables for next bit 
clearvars -except XYZ1 dl1 loc np 

%% Biot-Savart along central coil axis 
mu = pi*4e-7;
I = 1e-3; %Amps
const = I * mu/(4*pi); 

%Consider inner part of the coils as the area to optimise
q = linspace((0.1).*loc(1),(0.1).*loc(2),100);

v = [zeros(size(q)); q; zeros(size(q))]';
Bin = zeros(length(v),3); %Magnetic field map

%Field along symmetric axis
for m = 1:length(v) %Number of points in space
    for k = 1:length(XYZ1) %Number of wire elements
 
            rprim = v(m,:) - XYZ1(k,:); 
            dB = const.*(cross(dl1(k,:),rprim)./(abs(rprim).^3));
            Bin(m,:) = Bin(m,:) + dB;

    end
end

Bin(isnan(Bin)) = 0; 

y = v(:,2);
By = Bin(:,2);

%% Calculating Initial polynomial

fit = [y.^2 ones(size(y))]\By; % B1x^2 and B0 constants
fit = [fit(1),0,fit(2)];
targ = fit(1); 

%% Optimising inner loop

%Setting ranging values 
Is = 0.5.*I; %Current scales linearly w field (Picked current Will found)
const = Is .* mu/(4*pi);
Ls = linspace(0.1,0.35,26); %Lengths

for i = 1:length(Is)
    for l = 1:length(Ls)

        [XYZi, dli] = writecoily(Ls(l),loc,np,2,2); %Reverse direction
        Bin = zeros(length(v),3);

        for m = 1:length(v)
            for k = 1:length(XYZi)
     
                rprim = v(m,:) - XYZi(k,:);
                dB = const(i).*(cross(dli(k,:),rprim)./(abs(rprim).^3));
                Bin(m,:) = Bin(m,:) + dB;

            end
        end

        Bin(isnan(Bin))=0;
        B(:,l,i) = Bin(:,2);

        fits = [y.^2 ones(size(y))]\B(:,l,i); %Solve for constants
        itarg(l,i) = fits(1);

    end
        if mod(i,5) == 0
                   fprintf('Iteration %d...\n',i);
        end

end

diffs = targ - itarg; %Taking the difference
idx = find(diffs == min(diffs)); %Finding minimum distance between x^2 constants

%disp(['Current: ' num2str(Is(row)) ', Length: ' num2str(Ls(col))])

[XYZdet,dldet] = writecoily(Ls(idx),loc,np,2,2);

figure(600)
plot(y,By,'r')
hold on
plot(y,B(:,idx),'b') %Far too big


%% Troubleshooting

%Look at coils
figure(1000)
quiver3(XYZ1(:,1),XYZ1(:,2),XYZ1(:,3),dl1(:,1),dl1(:,2),dl1(:,3),1)
axis equal
hold on
quiver3(XYZdet(:,1),XYZdet(:,2),XYZdet(:,3),dldet(:,1),dldet(:,2),dldet(:,3),1)


%Look at initial fit

% val = polyval(fit,y);
% 
% figure(2000)
% plot(y,val)

%Look at x2 constants through loops

figure(3000)
plot(diffs) %Looking at where the difference between x2 constants is lowest










