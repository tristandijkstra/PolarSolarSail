Planet_1=3;
ibody_E_vect=3;
mu_Sun=1.3271e11;
ksun=mu_Sun;
mjd2000_dep=2460313.50000-2451544.5 
dth=0.001;
[kep_E_dep,ksun] = uplanet(mjd2000_dep, ibody_E_vect);
[rr_0_E,vv_0_E]=par2car(kep_E_dep(1),kep_E_dep(2),kep_E_dep(3),kep_E_dep(4),kep_E_dep(5),kep_E_dep(6),ksun);
plotOrbit(50*kep_E_dep(1),kep_E_dep(2),kep_E_dep(3),kep_E_dep(4),kep_E_dep(5),kep_E_dep(6),kep_E_dep(6)+2*pi,dth,mu_Sun);
hold on
%plot Earth
plot3(50*rr_0_E(1),50*rr_0_E(2),50*rr_0_E(3), 'bo' ,'LineWidth',3)
hold on
opts_example2.RefPlane = 'ecliptic';
planet3D('Sun',opts_example2);
%plot3(0 , 0 , 0, 'yo', 'Linewidth', 5)
rapo=kep_E_dep(1)*(1+kep_E_dep(3))
%study for the final orbit
%Theta will be defined, aswell as 
eccentricity= linspace(0,0.35, 50);
inclination=linspace(kep_E_dep(3), kep_E_dep(3)+pi/6,50);
r=rr_0_E;
rnorm=rr_0_E/norm(rr_0_E);
for i=1:length(eccentricity)
    a(i)=rapo/(1+eccentricity(i));
    P(i)=a(i)*(1-eccentricity(i)^2);
    VapoMagn(i)=sqrt(mu_Sun/P(i))*(1-eccentricity(i));
    for j=1:length(inclination)
        incl(j)=kep_E_dep(3)+inclination(j);
        if cross(rr_0_E,vv_0_E)~=0
           VPlane= vv_0_E-dot(rnorm,vv_0_E)*rnorm;
           Vrot=VPlane*cos(incl(j))+ cross(rnorm, Vplane)*sin(incl(j))+rnorm*(dot(rnorm,VPlane))*(1-cos(incl(j)));
           Direction=Vrot/norm(Vrot);
           VelVec=VapoMagn(i)*Direction;
        else
            Vrot=vv_0_E*cos(incl(j))+ cross(rnorm, vv_0_E)*sin(incl(j))+rnorm*(dot(rnorm,vv_0_E))*(1-cos(incl(j)));
            Direction=Vrot/norm(Vrot);
            VelVec=VapoMagn(i)*Direction;
        end
        DeltaV=VelVec-vv_0_E;
        DV(i,j)=norm(DeltaV);
        C3(i,j)=DV(i,j)^2;
        if C3(i,j)>100
            C3(i,j)=NaN;
        else
            C3(i,j)=C3(i,j);
        end
    end
end

 figure
 contour(inclination, eccentricity, DV)
 ylabel 'Different Inclination Values'
 xlabel 'Different Eccentricity'
 zlabel 'Characteristic energy [km^2/s^2]'
 colormap
 colorbar
 grid on
 title 'C_3'
 
 figure
 surf(inclination,eccentricity,C3)
 colorbar
 ylabel 'Different eccentricity values'
 xlabel 'Different inclination values'
 zlabel 'C_3 [km^2/s^2]'
 title 'C_3'
% 
% 
% figure
% [X, Y]=meshgrid(eccentricity, inclination);
% contour(Y, X, DV)
% xlabel 'Different Inclination'
% ylabel 'Different Eccentricity'
% zlabel 'Characteristic energy [\km^2/s^2]'
% colormap
% grid on
% title '\DeltaVdep'
% 
% DV(0,10);
% inclination(10);




%That is just for a specific day.. Once we obtain other valus we can
%optimize for specific days around the Earth's year

%%
%Computation for good inclination changes using different days of the year
%setting the eccentricity of the final orbit to be zero
clc
clear all
close all
mjd2000_dep=linspace(2460313.50000-2451544.5,  2460313.50000-2451544.5+365, 365);
mjd2000_dep=linspace(64328-51544.5, 64328+365-51544.5, 364)
Planet_1=3;
ibody_E_vect=3;
mu_Sun=1.3271e11;
dth=0.001;
inclination=linspace(0, pi/6,365);
for i=1:length(mjd2000_dep)
    [kep_E_dep(i,:),ksun] = uplanet(mjd2000_dep(i), ibody_E_vect);
    [rr_0_E(i,:),vv_0_E(i,:)]=par2car(kep_E_dep(i,1),kep_E_dep(i,2),kep_E_dep(i,3),kep_E_dep(i,4),kep_E_dep(i,5),kep_E_dep(i,6),ksun);
    eccentricity=kep_E_dep(i,2);
    r(i,:)=rr_0_E(i,:);
    rnorm(i,:)=r(i,:)/norm(r(i,:));
    rapo(i)=kep_E_dep(i,1)*(1+kep_E_dep(i,3))
    a(i)=rapo(i)/(1+eccentricity);
    P(i)=a(i)*(1-eccentricity^2);
    VapoMagn(i)=sqrt(mu_Sun/P(i))*(1-eccentricity);
    for j=1:length(inclination)
        incl(j)=kep_E_dep(i,3)+inclination(j);
        if cross(rr_0_E(i,:),vv_0_E(i,:))~=0
           VPlane(i,:)= vv_0_E(i,:)-dot(rnorm(i,:),vv_0_E(i,:))*rnorm(i,:);
           Vrot(i,:)=VPlane(i,:)*cos(incl(j))+ cross(rnorm(i,:), Vplane(i,:))*sin(incl(j))+rnorm(i,:)*(dot(rnorm(i,:),VPlane(i,:)))*(1-cos(incl(j)));
           Direction(i,:)=Vrot(i,:)/norm(Vrot(i,:));
           VelVec(i,:)=VapoMagn(i)*Direction(i,:);
        else
            Vrot(i,:)=vv_0_E(i,:)*cos(incl(j))+ cross(rnorm(i,:), vv_0_E(i,:))*sin(incl(j))+rnorm(i,:)*(dot(rnorm(i,:),vv_0_E(i,:)))*(1-cos(incl(j)));
            Direction(i,:)=Vrot(i,:)/norm(Vrot(i,:));
            VelVec(i,:)=VapoMagn(i)*Direction(i,:);
        end
        DeltaV=VelVec(i,:)-vv_0_E(i,:);
        DV(i,j)=norm(DeltaV);
        C3(i,j)=DV(i,j)^2;
        if C3(i,j)>25
            C3(i,j)=NaN;
        else
            C3(i,j)=C3(i,j);
        end
    end
end
% 
figure(1)
[X, Y]=meshgrid(mjd2000_dep(1,:), inclination(1,:));
contour(inclination(1,:),datenum(mjd2000_dep(1,:)), C3)
datetick('y','dd-mmm-yy')
colorbar 
xlabel 'Departure date from Earth'
ylabel 'Inclination Changes'
zlabel 'C3 [km^2/s^2]'
grid on
title 'C3'


figure(2)
surf(inclination(1,:),datenum(mjd2000_dep(1,:)),C3)
datetick('y','dd-mmm-yy')
colorbar
xlabel 'Different Inclination changes'
ylabel 'Different Date of Departure'
zlabel 'C_3 [km^2/s^2]'
title 'C_3'

%THE 183th element of the inclination change is the change of 15 degrees in
%inclination

[m, index]=min(C3(:,118));
[M, index2]=max(C3(:,118));
val=sqrt(m);
VAL=sqrt(M);
VAL-val

day=mjd2000_dep(index)+51544.5

%Plot the orbit launch for the case of the lowest delta V for the highest
%inclination change possible 

figure(3)
[kep_E_dep,ksun] = uplanet(mjd2000_dep(index), ibody_E_vect);
[rr_0_E,vv_0_E]=par2car(kep_E_dep(1),kep_E_dep(2),kep_E_dep(3),kep_E_dep(4),kep_E_dep(5),kep_E_dep(6),ksun);
plotOrbit(50*kep_E_dep(1),kep_E_dep(2),kep_E_dep(3),kep_E_dep(4),kep_E_dep(5),kep_E_dep(6),kep_E_dep(6)+2*pi,dth,mu_Sun);
hold on
opts_example2.RefPlane = 'ecliptic';
planet3D('Sun',opts_example2);
plotOrbit(50*kep_E_dep(1),kep_E_dep(2),kep_E_dep(3)+inclination(116),kep_E_dep(4),kep_E_dep(5),kep_E_dep(6),kep_E_dep(6)+2*pi,dth,mu_Sun);
xlabel 'X axis'
ylabel 'Y axis'
zlabel 'Z axis'
legend('Earth Orbit', 'Sun', 'Injected orbit after the launch')

%%
%Same Code, setting the eccentricity to zero
clc
clear all
close all
mjd2000_dep=linspace(2464328.5-2451544.5,  2464328.5-2451544.5+365, 366);
%mjd2000_dep=linspace(64328-51544.5, 64328+365-51544.5, 364)
Planet_1=3;
ibody_E_vect=3;
mu_Sun=1.3271e11;
dth=0.001;
inclination=linspace(0, pi/6,365);
for i=1:length(mjd2000_dep)
    [kep_E_dep(i,:),ksun] = uplanet(mjd2000_dep(i), ibody_E_vect);
    [rr_0_E(i,:),vv_0_E(i,:)]=par2car(kep_E_dep(i,1),kep_E_dep(i,2),kep_E_dep(i,3),kep_E_dep(i,4),kep_E_dep(i,5),kep_E_dep(i,6),ksun);
    eccentricity=0;
    r(i,:)=rr_0_E(i,:);
    rnorm(i,:)=r(i,:)/norm(r(i,:));
    rapo(i)=kep_E_dep(i,1)*(1+kep_E_dep(i,3))
    a(i)=rapo(i)/(1+eccentricity);
    P(i)=a(i)*(1-eccentricity^2);
    VapoMagn(i)=sqrt(mu_Sun/P(i))*(1-eccentricity);
    for j=1:length(inclination)
        incl(j)=kep_E_dep(i,3)+inclination(j);
        if cross(rr_0_E(i,:),vv_0_E(i,:))~=0
           VPlane(i,:)= vv_0_E(i,:)-dot(rnorm(i,:),vv_0_E(i,:))*rnorm(i,:);
           Vrot(i,:)=VPlane(i,:)*cos(incl(j))+ cross(rnorm(i,:), Vplane(i,:))*sin(incl(j))+rnorm(i,:)*(dot(rnorm(i,:),VPlane(i,:)))*(1-cos(incl(j)));
           Direction(i,:)=Vrot(i,:)/norm(Vrot(i,:));
           VelVec(i,:)=VapoMagn(i)*Direction(i,:);
        else
            Vrot(i,:)=vv_0_E(i,:)*cos(incl(j))+ cross(rnorm(i,:), vv_0_E(i,:))*sin(incl(j))+rnorm(i,:)*(dot(rnorm(i,:),vv_0_E(i,:)))*(1-cos(incl(j)));
            Direction(i,:)=Vrot(i,:)/norm(Vrot(i,:));
            VelVec(i,:)=VapoMagn(i)*Direction(i,:);
        end
        DeltaV=VelVec(i,:)-vv_0_E(i,:);
        DV(i,j)=norm(DeltaV);
        C3(i,j)=DV(i,j)^2;
        if C3(i,j)>25
            C3(i,j)=NaN;
        else
            C3(i,j)=C3(i,j);
        end
    end
end
% 
figure(1)
[X, Y]=meshgrid(mjd2000_dep(1,:), inclination(1,:));
contour(rad2deg(inclination(1,:)),datenum(mjd2000_dep(1,:)), C3)
datetick('y','dd-mmm-yy')
colorbar 
xlabel 'Departure date from Earth'
ylabel '\Delta i [deg]'
zlabel 'C_3 [km^2/s^2]'
grid on
title 'C_3'


figure(2)
surf(rad2deg(inclination(1,:)),datenum(mjd2000_dep(1,:)),C3)
datetick('y','dd-mmm-yy')
colorbar
xlabel '\Delta i [deg]'
ylabel 'Different Date of Departure'
zlabel 'C_3 [km^2/s^2]'
title 'C_3'

%THE 183th element of the inclination change is the change of 15 degrees in
%inclination

[m, index]=min(C3(:,118));
[M, index2]=max(C3(:,118));
val=sqrt(m);
VAL=sqrt(M);
VAL-val

day=mjd2000_dep(index)+51544.5

MJD2000_NEW=mjd2000_dep+51544.5

%Plot the orbit launch for the case of the lowest delta V for the highest
%inclination change possible 

figure(3)
[kep_E_dep,ksun] = uplanet(mjd2000_dep(index), ibody_E_vect);
[rr_0_E,vv_0_E]=par2car(kep_E_dep(1),kep_E_dep(2),kep_E_dep(3),kep_E_dep(4),kep_E_dep(5),kep_E_dep(6),ksun);
plotOrbit(50*kep_E_dep(1),kep_E_dep(2),kep_E_dep(3),kep_E_dep(4),kep_E_dep(5),kep_E_dep(6),kep_E_dep(6)+2*pi,dth,mu_Sun);
hold on
opts_example2.RefPlane = 'ecliptic';
planet3D('Sun',opts_example2);
plotOrbit(50*kep_E_dep(1),kep_E_dep(2),kep_E_dep(3)+inclination(116),kep_E_dep(4),kep_E_dep(5),kep_E_dep(6),kep_E_dep(6)+2*pi,dth,mu_Sun);
xlabel 'X axis'
ylabel 'Y axis'
zlabel 'Z axis'
legend('Earth Orbit', 'Sun', 'Injected orbit after the launch')

%154
inclination(117)

%%
% We can compute the array for the DeltaV
[kep_E_dep,ksun] = uplanet(mjd2000_dep(154), ibody_E_vect);
[rr_0_E,vv_0_E]=par2car(kep_E_dep(1),kep_E_dep(2),kep_E_dep(3),kep_E_dep(4),kep_E_dep(5),kep_E_dep(6),ksun);

% Compute the Array
eccentricity=0;
rapo2=kep_E_dep(1)*(1+kep_E_dep(3))
a2=rapo2/(1+eccentricity);
P2=a2*(1-eccentricity^2);
VapoMagn2=sqrt(mu_Sun/P2)*(1-eccentricity);
rnorm=rr_0_E/norm(rr_0_E);
VPlane= vv_0_E-dot(rnorm,vv_0_E)*rnorm;
Vrot=VPlane*cos(inclination(117))+ cross(rnorm, VPlane)*sin(inclination(117))+rnorm*(dot(rnorm,VPlane))*(1-cos(inclination(117)));
Direction=Vrot/norm(Vrot);
VelVec2=VapoMagn2*Direction;
DeV=VelVec2-vv_0_E
DeV=norm(DeV)
C3=DeV^2

%%
%Approximate Hohmann Transfer
%We exlude the launch part, we start from the Earth SOI
mu_Sun=1.3271e11;
dE=150000000; %km
V_Earth=sqrt(mu_Sun/dE);
finalD=0.48*dE;
eT=(dE-finalD)/(dE+finalD);
aT=(dE+finalD)/2;
PT=aT*(1-eT^2);
V_apoT=sqrt(mu_Sun/PT)*(1-eT);
DeltaV1=abs(V_apoT-V_Earth);
V_periT=sqrt(mu_Sun/PT)*(1+eT);
Vcirc=sqrt(mu_Sun/finalD);
DeltaV2=abs(Vcirc-V_periT);
DeltaV=DeltaV1+DeltaV2

%Voice, My part.. I divided.. Preliminary Design?? 2024-2025?



