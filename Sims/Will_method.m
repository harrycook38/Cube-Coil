clear all
close all
clc


%% 3D

wire = linspace(-0.25, 0.25, 50);
wire_pos = 0.25*ones(size(wire));
wire_neg = -0.25*ones(size(wire));

d = wire(1,2) - wire(1,1);

%Field point environment
x = linspace(-0.3, 0.3, 20);
y = linspace(-0.4, 0.4, 20); % coil sizes are linked to field calculations
z = linspace(-0.5, 0.5, 20);

[X,Y,Z] = meshgrid(x, y, z);

I = 1;
mu = pi*4e-7;
C = (mu*I)/(4*pi);

%zx Coils

[Bxz_x_ur, Bxz_y_ur, Bxz_z_ur] = B_xwire(X, Y, Z, wire, wire_pos, wire_pos, d, C);
[Bxz_x_ul, Bxz_y_ul, Bxz_z_ul] = B_xwire(X, Y, Z, wire, wire_neg, wire_pos, d, C);
[Bxz_x_lr, Bxz_y_lr, Bxz_z_lr] = B_xwire(X, Y, Z, wire, wire_pos, wire_neg, d, -C);
[Bxz_x_ll, Bxz_y_ll, Bxz_z_ll] = B_xwire(X, Y, Z, wire, wire_neg, wire_neg, d, -C);

Bxz_x = Bxz_x_ur + Bxz_x_ul + Bxz_x_lr + Bxz_x_ll;
Bxz_y = Bxz_y_ur + Bxz_y_ul + Bxz_y_lr + Bxz_y_ll;
Bxz_z = Bxz_z_ur + Bxz_z_ul + Bxz_z_lr + Bxz_z_ll;

[Bzx_x_ur, Bzx_y_ur, Bzx_z_ur] = B_zwire(X, Y, Z, wire_pos, wire_pos, wire, d, -C);
[Bzx_x_ul, Bzx_y_ul, Bzx_z_ul] = B_zwire(X, Y, Z, wire_pos, wire_neg, wire, d, -C);
[Bzx_x_lr, Bzx_y_lr, Bzx_z_lr] = B_zwire(X, Y, Z, wire_neg, wire_pos, wire, d, C);
[Bzx_x_ll, Bzx_y_ll, Bzx_z_ll] = B_zwire(X, Y, Z, wire_neg, wire_neg, wire, d, C);

Bzx_x = Bzx_x_ur + Bzx_x_ul + Bzx_x_lr + Bzx_x_ll;
Bzx_y = Bzx_y_ur + Bzx_y_ul + Bzx_y_lr + Bzx_y_ll;
Bzx_z = Bzx_z_ur + Bzx_z_ul + Bzx_z_lr + Bzx_z_ll;

BXZ_X = Bzx_x + Bxz_x;
BXZ_Y = Bzx_y + Bxz_y;
BXZ_Z = Bzx_z + Bxz_z;

%yz coils 

[Bzy_x_ur, Bzy_y_ur, Bzy_z_ur] = B_zwire(X, Y, Z, wire_pos, wire_pos, wire, d, -C);
[Bzy_x_ul, Bzy_y_ul, Bzy_z_ul] = B_zwire(X, Y, Z, wire_pos, wire_neg, wire, d, C);
[Bzy_x_lr, Bzy_y_lr, Bzy_z_lr] = B_zwire(X, Y, Z, wire_neg, wire_pos, wire, d, -C);
[Bzy_x_ll, Bzy_y_ll, Bzy_z_ll] = B_zwire(X, Y, Z, wire_neg, wire_neg, wire, d, C);

Bzy_x = Bzy_x_ur + Bzy_x_ul + Bzy_x_lr + Bzy_x_ll;
Bzy_y = Bzy_y_ur + Bzy_y_ul + Bzy_y_lr + Bzy_y_ll;
Bzy_z = Bzy_z_ur + Bzy_z_ul + Bzy_z_lr + Bzy_z_ll;

[Byz_x_ur, Byz_y_ur, Byz_z_ur] = B_ywire(X, Y, Z, wire_pos, wire, wire_pos, d, C);
[Byz_x_ul, Byz_y_ul, Byz_z_ul] = B_ywire(X, Y, Z, wire_pos, wire, wire_neg, d, -C);
[Byz_x_lr, Byz_y_lr, Byz_z_lr] = B_ywire(X, Y, Z, wire_neg, wire, wire_pos, d, C);
[Byz_x_ll, Byz_y_ll, Byz_z_ll] = B_ywire(X, Y, Z, wire_neg, wire, wire_neg, d, -C);

Byz_x = Byz_x_ur + Byz_x_ul + Byz_x_lr + Byz_x_ll;
Byz_y = Byz_y_ur + Byz_y_ul + Byz_y_lr + Byz_y_ll;
Byz_z = Byz_z_ur + Byz_z_ul + Byz_z_lr + Byz_z_ll;

BYZ_X = Bzy_x + Byz_x;
BYZ_Y = Bzy_y + Byz_y;
BYZ_Z = Bzy_z + Byz_z;

%xy coil

[Bxy_x_ur, Bxy_y_ur, Bxy_z_ur] = B_xwire(X, Y, Z, wire, wire_pos, wire_pos, d, -C);
[Bxy_x_ul, Bxy_y_ul, Bxy_z_ul] = B_xwire(X, Y, Z, wire, wire_neg, wire_pos, d, C);
[Bxy_x_lr, Bxy_y_lr, Bxy_z_lr] = B_xwire(X, Y, Z, wire, wire_pos, wire_neg, d, -C);
[Bxy_x_ll, Bxy_y_ll, Bxy_z_ll] = B_xwire(X, Y, Z, wire, wire_neg, wire_neg, d, C);

Bxy_x = Bxy_x_ur + Bxy_x_ul + Bxy_x_lr + Bxy_x_ll;
Bxy_y = Bxy_y_ur + Bxy_y_ul + Bxy_y_lr + Bxy_y_ll;
Bxy_z = Bxy_z_ur + Bxy_z_ul + Bxy_z_lr + Bxy_z_ll;

[Byx_x_ur, Byx_y_ur, Byx_z_ur] = B_ywire(X, Y, Z, wire_pos, wire, wire_pos, d, C);
[Byx_x_ul, Byx_y_ul, Byx_z_ul] = B_ywire(X, Y, Z, wire_pos, wire, wire_neg, d, C);
[Byx_x_lr, Byx_y_lr, Byx_z_lr] = B_ywire(X, Y, Z, wire_neg, wire, wire_pos, d, -C);
[Byx_x_ll, Byx_y_ll, Byx_z_ll] = B_ywire(X, Y, Z, wire_neg, wire, wire_neg, d, -C);

Byx_x = Byx_x_ur + Byx_x_ul + Byx_x_lr + Byx_x_ll;
Byx_y = Byx_y_ur + Byx_y_ul + Byx_y_lr + Byx_y_ll;
Byx_z = Byx_z_ur + Byx_z_ul + Byx_z_lr + Byx_z_ll;

BXY_X = Bxy_x + Byx_x;
BXY_Y = Bxy_y + Byx_y;
BXY_Z = Bxy_z + Byx_z;


Bx = BXZ_X + BYZ_X + BXY_X;
By = BXZ_Y + BYZ_Y + BXY_Y;
Bz = BXZ_Z + BYZ_Z + BXY_Z;

clearvars -except Bx By Bz X Y Z wire wire_pos wire_neg

figure(1)
plot3(wire, wire_pos, wire_pos)
hold on
plot3(wire, wire_pos, wire_neg)
plot3(wire, wire_neg, wire_pos)
plot3(wire, wire_neg, wire_neg)

plot3(wire_pos, wire, wire_pos)
plot3(wire_pos, wire, wire_neg)
plot3(wire_neg, wire, wire_pos)
plot3(wire_neg, wire, wire_neg)

plot3(wire_pos, wire_pos, wire)
plot3(wire_neg, wire_pos, wire)
plot3(wire_pos, wire_neg, wire)
plot3(wire_neg, wire_neg, wire)

quiver3(X,Y,Z,Bx,By,Bz)


%% Writing functions to calcuate field


function [Bxx,Byy,Bzz] = B_xwire(x,y,z,xp,yp,zp,dx,C)
%Field from x wire

Bxx = zeros(size(x));
Byy = zeros(size(y));
Bzz = zeros(size(z));

for i = 1:length(xp)


    xd = x - xp(i);
    yd = y - yp(i);
    zd = z - zp(i);

    r = sqrt(xd.^2+yd.^2+zd.^2);
    theta = acos(zd./r);
    psi = atan2(yd,xd);

    By = (C/(r.^2)).*dx.*(-1).*cos(theta);
    Bz = (C/(r.^2)).*dx.*sin(theta).*sin(psi);
    
    Byy = Byy +By;
    Bzz = Bzz +Bz;

end
end

function [Bxx,Byy,Bzz] = B_ywire(x,y,z,xp,yp,zp,dy,C)
%Field from y wire 

Bxx = zeros(size(x));
Byy = zeros(size(y));
Bzz = zeros(size(z));

for i = 1:length(xp)

    xd = x - xp(i);
    yd = y - yp(i);
    zd = z - zp(i);

    r = sqrt(xd.^2+yd.^2+zd.^2);
    theta = acos(zd./r);
    psi = atan2(yd,xd);

    Bx = (C/(r.^2)).*dy.*cos(theta);
    Bz = (C/(r.^2)).*dy.*(-1).*sin(theta).*cos(psi);
    
    Bxx = Bxx + Bx;
    Bzz = Bzz + Bz;

end
end

function [Bxx,Byy,Bzz] = B_zwire(x,y,z,xp,yp,zp,dz,C)
%Field from y wire 

Bxx = zeros(size(x));
Byy = zeros(size(y));
Bzz = zeros(size(z));

for i = 1:length(xp)

    xd = x - xp(i);
    yd = y - yp(i);
    zd = z - zp(i);

    r = sqrt(xd.^2+yd.^2+zd.^2);
    theta = acos(zd./r);
    psi = atan2(yd,xd);

    Bx = (C/(r.^2)).*dz.*(-1).*sin(theta).*sin(psi);
    By = (C/(r.^2)).*dz.*sin(theta).*cos(psi);
    
    Bxx = Bxx + Bx;
    Byy = Byy + By;

end
end



