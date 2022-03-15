clear all
close all
clc

%Environment to calculate B field on
di = 1; 
np = 30;
d = linspace(-0.5*di,0.5*di,np);

[X,Y,Z] = meshgrid(d);
pp = [X(:),Y(:),Z(:)];

%Coil diameters
d1sep = 0.3; % Different coil sizes to account for physical standoff

%Pair 1
l1 = linspace(-0.5*d1sep,0.5*d1sep,2.*length(d));
orgsep = (0.5445./2)*d1sep; %Optimal square coil separation. 
loc = [-orgsep orgsep];

for k = 1:2
    pos = loc(k).*ones(1,length(l1));
    X1(:,:,k) = [l1;pos;l1(1).*ones(1,length(l1))];
    X2(:,:,k) = [l1(end).*ones(1,length(l1));pos;l1];
    X3(:,:,k) = [-1.*l1;pos;l1(end).*ones(1,length(l1))];
    X4(:,:,k) = [l1(1).*ones(1,length(l1));pos;-1.*l1];
    W(:,:,k) = [X1(:,:,k) X2(:,:,k) X3(:,:,k) X4(:,:,k)];
end

W = reshape(W,3,[])';

dlX = diff(W); %Line element vectors
dlX(end+1,:) = dlX(end,:);
dlX(end./2,:) = dlX((end./2)-1,:);

%TESTING

figure(500)
plot3(W(:,1),W(:,2),W(:,3),'bo');
hold on
xlabel('x')
ylabel('y')

quiver3(W(:,1),W(:,2),W(:,3),dlX(:,1),dlX(:,2),dlX(:,3),3) %Check dls are correct
axis equal


%% Clearing unessential variables for next bit
clearvars -except W dlX X Y Z np pp

%% Biot-Savart

mu = pi*4e-7; 
I = 1; %Amps
const = I * mu/(4*pi);
N = 60;


Xw = W(:,1); dlx = dlX(:,1);
Yw = W(:,2); dly = dlX(:,2);
Zw = W(:,3); dlz = dlX(:,3);

% Bout = zeros(s)

for m = 1:length(pp) %Number of points in space

            R = pp(m,:) - W;
            Rx = R(:,1);
            Ry = R(:,2);
            Rz = R(:,3);

      dBx = 0; dBy = 0; dBz = 0;
     for k = 1:length(Xw) %Number of wire elements
         
            Xc=(dly(k).*Rz(k))-(dlz(k).*Ry(k));
            Yc=(dlz(k).*Rx(k))-(dlx(k).*Rz(k));
            Zc=(dlx(k).*Ry(k))-(dly(k).*Rx(k));

            Bx1=(const./(abs(R(k)).^3)).*Xc;
            By1=(const./(abs(R(k)).^3)).*Yc;
            Bz1=(const./(abs(R(k)).^3)).*Zc;

            dBx = dBx + Bx1;
            dBy = dBy + By1;
            dBz = dBz + Bz1;

     end   
    dB = [dBx dBy dBz];
    Bout(m,:) = dB; 

            if mod(m,1000) == 0
                fprintf('Iteration %d...\n',m);
            end
end

figure(1)
quiver3(pp(:,1),pp(:,2),pp(:,3),Bout(:,1),Bout(:,2),Bout(:,3),4,'c');
hold on
plot3(Xw,Yw,Zw,'or')
axis equal
  
xlabel('x')
ylabel('y')
zlabel('z')


% for a = 1:np
% for b = 1:np
% for c = 1:np
% for i = 1:length(Xw)
% 
%     Rx(i) = X(a,b,c) - Xw(i);
%     Ry(i) = Y(a,b,c) - Yw(i);
%     Rz(i) = Z(a,b,c) - Zw(i);
%     R(i) = sqrt(Rx(i).^2+Ry(i).^2+Rz(i).^2);
%     
% end
% end
% end
% end
% 


% Rx(N)=(X(a,b,c)-0.5*(Xw(N)+1));
% Ry(N)=(Y(a,b,c)-(0.5*(Yw(N)+1)));
% Rz(N)=(Z(a,b,c)-(0.5*(Zw(N)+1)));
% dlx(N)=-Xw(N)+1;
% dly(N)=-Yw(N)+1;
% dlz(N)=-Zw(N)+1;
% for i=1:N
% Xcross(i)=(dly(i).*Rz(i))-(dlz(i).*Ry(i));
% Ycross(i)=(dlz(i).*Rx(i))-(dlx(i).*Rz(i));
% Zcross(i)=(dlx(i).*Ry(i))-(dly(i).*Rx(i));
% R(i)=sqrt(Rx(i).^2+Ry(i).^2+Rz(i).^2);
% end
% Bx1=(I*u0./(4*pi*(R.^3))).*Xcross;
% By1=(I*u0./(4*pi*(R.^3))).*Ycross;
% Bz1=(I*u0./(4*pi*(R.^3))).*Zcross;
% BX(a,b,c)=0;       
% BY(a,b,c)=0;
% BZ(a,b,c)=0;
% for i=1:N     
%     BX(a,b,c)=BX(a,b,c)+Bx1(i);
%     BY(a,b,c)=BY(a,b,c)+By1(i);
%     BZ(a,b,c)=BZ(a,b,c)+Bz1(i);
% end
% end
% end
% end
% 

% 
% 
% for a=1:np  
% for b=1:np
% for c=1:np
% for i=1:length(W)-1
%     
% Rx(i)=(X(a,b,c)-(0.5*(Xw(i)+Xw(i+1))));
% Ry(i)=(Y(a,b,c)-(0.5*(Yw(i)+Yw(i+1))));
% Rz(i)=(Z(a,b,c)-(0.5*(Zw(i)+Zw(i+1))));
% dlx(i)=Xw(i+1)-Xw(i);
% dly(i)=Yw(i+1)-Yw(i);
% dlz(i)=Zw(i+1)-Zw(i);
% end
% Rx(N)=(X(a,b,c)-0.5*(Xw(N)+1));
% Ry(N)=(Y(a,b,c)-(0.5*(Yw(N)+1)));
% Rz(N)=(Z(a,b,c)-(0.5*(Zw(N)+1)));
% dlx(N)=-Xw(N)+1;
% dly(N)=-Yw(N)+1;
% dlz(N)=-Zw(N)+1;
% for i=1:N
% Xcross(i)=(dly(i).*Rz(i))-(dlz(i).*Ry(i));
% Ycross(i)=(dlz(i).*Rx(i))-(dlx(i).*Rz(i));
% Zcross(i)=(dlx(i).*Ry(i))-(dly(i).*Rx(i));
% R(i)=sqrt(Rx(i).^2+Ry(i).^2+Rz(i).^2);
% end
% Bx1=(const./(R.^3)).*Xcross;
% By1=(const./(R.^3)).*Ycross;
% Bz1=(const./(R.^3)).*Zcross;
% BX(a,b,c)=0;       
% BY(a,b,c)=0;
% BZ(a,b,c)=0;
% for i=1:N     
%     BX(a,b,c)=BX(a,b,c)+Bx1(i);
%     BY(a,b,c)=BY(a,b,c)+By1(i);
%     BZ(a,b,c)=BZ(a,b,c)+Bz1(i);
% end
% 
% end
% end
% end
% 
% 
% figure(1)
% quiver3(Z,Y,X,BZ,BY,BX,3);
% hold on
% xlabel('X-axis')
% ylabel('Y-axis')
% zlabel('Z-axis')
% title('B-field of a current wire along Z-axis')
% h=gca; 
% fh = figure(1); 
% set(fh, 'color', 'white'); 