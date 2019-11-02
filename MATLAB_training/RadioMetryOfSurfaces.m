close all
clear all
clc
D=500;
r=180;
n=1.5;
v=[0;0;1];
thetai= pi/4;
phii= pi/3;
theta = zeros(D);
phi = zeros(D);
L=[sin(thetai)*cos(phii);sin(thetai)*sin(phii);cos(thetai)];
Z=L+v;
Normm=norm(Z);
H=Z./Normm;
N = zeros(1,3,D,D);
for i=1:D
    for j=1:D
        dd=((i-D/2)^2 +(j-D/2)^2);
        d=sqrt(dd);
        if (d < r)
            theta=asin(d/r);
            phi=atan2((j-D/2) , (i-D/2));
            
            N=[sin(theta).*sin(phi),sin(theta).*cos(phi),cos(theta)];
            
            theta_h = acos(dot(N,H));
            theta_i =acos(dot(L,N));
            theta_iprim =acos(dot(L,H));
            theta_r = acos(dot(v,N));
            theta_g = acos(dot(v,L));
            Nh=cos(theta_h);
            Li=cos(theta_i);
            vr=cos(theta_r);
            vg=cos(theta_g);
            Th=tan(theta_h);
            Ci=cot(theta_i);
            Cr=cot(theta_r);
            Cor=cos(theta_r);
            Coi=cos(theta_i);
            CoP=cos(theta_iprim);
            CP=cot(theta_iprim);
            Dh=exp(-(Th^2)/2)/2*pi*Nh^3;
            Atheta=[sqrt(2)*exp(-Ci^2/2)/sqrt(pi)*Ci-erfc(Ci/sqrt(2))]/2;
            AthetaR=[sqrt(2)*exp(-Cr^2/2)/sqrt(pi)*Cr-erfc(Cr/sqrt(2))]/2;
            if (CoP > 0)
                gTheta(i,j)=1/ Atheta+1;
                gThetaR(i,j)=1/ AthetaR+1;
                G(i,j)=(1/ Atheta+1).*(1/ AthetaR+1);
            else
                gTheta(i,j)=0;
                gThetaR(i,j)=0;
            end
            %C=2*Li*vr-vg;
            %L(theta_i,theta_r,theta_g) = C^n + Li
            I(i,j)=mean(dot(L,N));
        else
            I(i,j)=0;
        end
    end
end
%plot(I(:,:),I(:,:))
%imagesc(I)
%sphere(I,'b')
imshow(I)
title('Lambertian Ball');