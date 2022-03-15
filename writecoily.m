function [XYZ, dl] = writecoily(diam,loc,n,numcoil,dir)
%WRITE COIL WITH SYMMETRY IN THE Y AXIS

%sep - Coil diameter
%loc - location(s) of the coil
%N - number of wire elements on the side of a square
%numcoil - number of coils
%dir - 1 or 2, direction of current flow

%XYZ - Coordinates of wire
%dl - current elements

%Pair 1
l = linspace(-0.5*diam,0.5*diam,n); 

for k = 1:numcoil
    pos = loc(k).*ones(1,length(l));
    Y1(:,:,k) = [l;pos;l(1).*ones(1,length(l))];
    Y2(:,:,k) = [l(end).*ones(1,length(l));pos;l];
    Y3(:,:,k) = [-1.*l;pos;l(end).*ones(1,length(l))];
    Y4(:,:,k) = [l(1).*ones(1,length(l));pos;-1.*l];
    Y(:,:,k) = [Y1(:,:,k) Y2(:,:,k) Y3(:,:,k) Y4(:,:,k)];
end

Y = reshape(Y,3,[])';
XYZ = Y;

%Line element vectors
if dir == 1
    dlY = diff(Y);
elseif dir == 2
    dlY = diff(1-Y);
end
dlY(end+1,:) = dlY(end,:);
dlY(end./2,:) = dlY((end./2)-1,:);

dl = dlY;

end