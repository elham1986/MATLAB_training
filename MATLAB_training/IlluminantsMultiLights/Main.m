clear all
close all
clc

D=500;
r=180;
ng=1.5;%refrective index of the glass
v=[0;0;1];
thetai= pi/4;
phii= pi/3;
m=0.2;
n = 0.15016;%refrective index of Silver
k = 3.4727;%imaginary part of refrective index of Silver

[I,Imetal,Iglass]=RadioMetryOfSurfaces(D,r,ng,v,thetai,phii,m,n,k);

%plot(I(:,:),I(:,:))
%imagesc(I)
%sphere(I,'b')
figure(1)
imshow(I)
figure(2)
imagesc(I)
title('Lambertian Ball');

figure(3)
imshow(Imetal)

figure(4)
imagesc(Imetal)
title('Metal Ball');

figure(5)
imagesc(Iglass)

figure(6)
imshow(Iglass)
%scatter3(Iglass)
%colormap(Iglass)
title('Glass Ball');