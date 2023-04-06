
clc
clear
close all


Planet_1=3;
Planet_2=2;

orbitType=0;
% Physical constrain on the max delta_v that the spacecraft is able to grant:
delta_v_1_MAX=25;
delta_v_3_MAX=25;


%% STEP 0 : COUNTOUR PLOTS
% In order to visualize and study any possible pattern on the delta_v needed for the transfers Earth_Saturn and Saturn_NEO, has been performed
% a computation of delta_v_1 and delta_v_3 where delta_v_1 is the delta_v needed for the transfer from Earth's orbit to the first Lambert leg,
% delta_v_3 is the delta_v needed for the transfer from second Lambert leg to NEO's orbit. ATTENTION: this computation does not consider the flyby.
% In order to observe properly the results, delta_v_1 will be computed for several Earth-Saturn synodic periods (n_1) and delta_v_3 will be computed
% for several Saturn-NEO synodic periods (n_2). The time discretization used is unit_1 (days)

n_1=5;

unit_1=15;


% mjd_departure_fin_1=mjd_departure_in+n_1*T_syn_ES;
% mjd_arrival_in_1=mjd_departure_fin_1;
% mjd_arrival_fin_1=mjd_arrival_in_1 + (n_1+0.1)*T_syn_ES;

mjd_DEP_in=60000;
mjd_DEP_fin=62000;
mjd_ARR_in=60000;
mjd_ARR_fin=62000;
[delta_v_venus,delta_v_dep,v_final_Lambert_1_x_1,v_final_Lambert_1_y_1,v_final_Lambert_1_z_1] = Lambert_optimal_venus(mjd_DEP_in,mjd_DEP_fin,mjd_ARR_in,mjd_ARR_fin,unit_1,Planet_1,Planet_2,orbitType,delta_v_1_MAX);

mjd2000_departure_1=[mjd_DEP_in:unit_1:mjd_DEP_fin] - 51544.5;
mjd2000_arrival_1=[mjd_ARR_in:unit_1:mjd_ARR_fin] - 51544.5;

% figure
% contour(datenum(mjd2000_departure_1),datenum(mjd2000_arrival_1),delta_v_venus',ShowText="off")
% datetick('x','dd-mmm-yy'),datetick('y','dd-mmm-yy');
% colorbar 
%  xlabel 'Departure date from Earth'
% ylabel 'Arrival date on venus'
% zlabel 'delta_v_tot [km/s]'
% grid on
% title '\DeltaVvenus'

figure
contour(datenum(mjd2000_departure_1),datenum(mjd2000_arrival_1),delta_v_dep',ShowText="off")
datetick('x','dd-mmm-yy'),datetick('y','dd-mmm-yy');
colorbar 
 xlabel 'Departure date from Earth'
ylabel 'Arrival date on venus'
zlabel 'delta_v_dep [km/s]'
grid on
title '\DeltaVdep'

% figure
% contour(datenum(mjd2000_departure_1),datenum(mjd2000_arrival_1),delta_v_arr',ShowText="off")
% datetick('x','dd-mmm-yy'),datetick('y','dd-mmm-yy');
% colorbar 
%  xlabel 'Departure date from Earth'
% ylabel 'Arrival date on venus'
% zlabel 'delta_v_arr [km/s]'
% grid on
% title '\DeltaVarr'

% figure
% surf(datenum(mjd2000_departure_1),datenum(mjd2000_arrival_1),delta_v_venus')
% datetick('x','dd-mmm-yy'),datetick('y','dd-mmm-yy');
% colorbar
% xlabel 'Departure date from Earth'
% ylabel 'Arrival date on venus'
% zlabel 'delta_v_tot [km/s]'
% title '\DeltaVvenus'

figure
surf(datenum(mjd2000_departure_1),datenum(mjd2000_arrival_1),delta_v_dep')
datetick('x','dd-mmm-yy'),datetick('y','dd-mmm-yy');
colorbar
xlabel 'Departure date from Earth'
ylabel 'Arrival date on venus'
zlabel 'delta_v_dep [km/s]'
title '\DeltaVdep'

% figure
% surf(datenum(mjd2000_departure_1),datenum(mjd2000_arrival_1),delta_v_arr')
% datetick('x','dd-mmm-yy'),datetick('y','dd-mmm-yy');
% colorbar
% xlabel 'Departure date from Earth'
% ylabel 'Arrival date on venus'
% zlabel 'delta_v_arr [km/s'
% title '\DeltaVarr'

%% PLOT THE ORBITS

%% ATTENZIONE: GIORNI FISSATI A CASO DA UNA REGIONE DI MINIMO
mjd2000_dep=8530.5;
mjd2000_arr_venus=8695.5;

a_EARTH=14964958;
a_VENUS=10821000;
mu_Sun=1.3271e11;
T_E=2*pi*sqrt(a_EARTH.^3/mu_Sun);
T_VEN=2*pi*sqrt(a_VENUS(1).^3/mu_Sun);
T_syn_EV=((T_E*T_VEN)/abs(T_E-T_VEN))/(60*60*24);

% mjd2000_dep=7400.5;
% mjd2000_arr_venus=7565.5;

% mjd2000_dep=7400.5 + 2*T_syn_EV;
% mjd2000_arr_venus=7565.5 + 2*T_syn_EV;



ksun=mu_Sun;
dth=0.1;
% EARTH ORBIT AT REAL DEPARTURE DAY


ibody_E_vect=3;
[kep_E_dep,ksun] = uplanet(mjd2000_dep, ibody_E_vect);
[rr_0_E,vv_0_E]=par2car(kep_E_dep(1),kep_E_dep(2),kep_E_dep(3),kep_E_dep(4),kep_E_dep(5),kep_E_dep(6),ksun);

%Venus ORBIT AT REAL arrival DAY

ibody_Sat_vect=2;

[kep_venus,ksun] = uplanet(mjd2000_arr_venus, ibody_Sat_vect);
[rr_0_venus,vv_0_venus]=par2car(kep_venus(1),kep_venus(2),kep_venus(3),kep_venus(4),kep_venus(5),kep_venus(6),mu_Sun);





% plot Earth orbit
figure
plotOrbit(kep_E_dep(1),kep_E_dep(2),kep_E_dep(3),kep_E_dep(4),kep_E_dep(5),kep_E_dep(6),kep_E_dep(6)+2*pi,dth,mu_Sun);
hold on
%plot Earth
plot3(rr_0_E(1),rr_0_E(2),rr_0_E(3), 'bo' ,'LineWidth',3)
hold on

%plot venus orbit
plotOrbit(kep_venus(1),kep_venus(2),kep_venus(3),kep_venus(4),kep_venus(5),kep_venus(6),kep_venus(6)+2*pi,dth,mu_Sun);
%plot venus
hold on
plot3(rr_0_venus(1),rr_0_venus(2),rr_0_venus(3),'ro','LineWidth',3)
%plot the Sun
hold on
% opts_example2.RefPlane = 'ecliptic';
% planet3D('Sun',opts_example2);
plot3(0 , 0 , 0, 'yo', 'Linewidth', 5)


%% PLOT THE TRANSFERS


Nrev=0; 
optionsLMR=1;
orbitType=0;


% Transfer between Earth and Venus
TOF_T1=(mjd2000_arr_venus-mjd2000_dep)*24*60*60;
[a_T1,p_T1,e_T1,ERROR,v_i_T1,v_f_T1,TPAR_T1,THETA_T1] = lambertMR(rr_0_E,rr_0_venus,TOF_T1,mu_Sun,orbitType,Nrev,optionsLMR);

T_T1=(2*pi)*sqrt(a_T1^3/ksun);
tspan_T1=linspace(0,TOF_T1,1000);

%options
options = odeset('RelTol',1e-13,'AbsTol',1e-14);

[time_T1, state_T1 ] = ode45(@(t,s) twobody_problem_ode(t,s,ksun), tspan_T1, [rr_0_E v_i_T1'],options);

%plot the transfer orbit

hold on
plot3(state_T1(:,1), state_T1(:,2), state_T1(:,3) , 'r', 'LineWidth', 2);


legend('Earth orbit','Earth','Venus orbit','Venus','Sun','1st Leg')
title 'Earth-Venus transfer'



%% FIXING ARRIVAL ORBIT

%parameters at the arrival day on Venus
[kep_2,ksun] = uplanet(mjd2000_arr_venus,Planet_2);
a_2=kep_2(1);
e_2=kep_2(2);
i_2=kep_2(3);
OM_2=kep_2(4);
om_2=kep_2(5);
theta_0_2=kep_2(6);

% initial cartesian parameters
[rr_0_2(:),vv_0_2(:)] = par2car(a_2,e_2,i_2,OM_2,om_2,theta_0_2,ksun);



%% 
[kep_3,ksun] = uplanet(mjd2000_arr_venus+100,Planet_2);
a_3=kep_3(1);
e_3=kep_3(2);
i_3=kep_3(3);
OM_3=kep_3(4);
om_3=kep_3(5);
theta_0_3=kep_3(6);

% initial cartesian parameters
[rr_0_3(:),vv_0_3(:)] = par2car(a_3,e_3,i_3,OM_3,om_3,theta_0_3,ksun);
%%

AS=1.49597870e8;
%arrival orbit
r_apo_arr=norm(rr_0_2);
e_arr=0.055;
a_arr=r_apo_arr/(1+e_arr);

r_p_arr=a_arr*(1-e_arr);
APO=r_apo_arr/AS
PERI=r_p_arr/AS

vv_apo_norm=sqrt(ksun/(a_arr*(1-e_arr^2)))*(1-e_arr);

% define velocity versor 
vec=cross(cross(rr_0_2,rr_0_3),rr_0_2);
ver=vec/norm(vec);
vv_apo_arr=vv_apo_norm*ver;

%[a_arr,e_arr,i_arr,OM_arr,om_arr, theta_0_2 ] = car2par(rr_0_2, vv_apo_arr,mu_Sun);
[a_arr, e_arr, i_arr,igradarr, OM_arr,OMEGAgrad, om_arr,omegagrad, theta_0_2, thetagrad]=car2par(rr_0_2,vv_apo_arr,mu_Sun)

V_minus_SC=v_f_T1; %heliocentric vel BEFORE the flyby
vv_0_VEN_flyby=vv_0_2; %vel of Venus the day of the flyby
V_plus_SC=vv_apo_arr; %heliocentric velocity AFTER the flyby


v_minus_infinite= V_minus_SC - vv_0_VEN_flyby ;
v_plus_infinite= V_plus_SC - vv_0_VEN_flyby;
delta_v=v_plus_infinite - v_minus_infinite;

v_minus_infinite_norm=norm(v_minus_infinite);
v_plus_infinite_norm=norm(v_plus_infinite);

turn_angle=acos(dot(v_minus_infinite,v_plus_infinite)/(v_plus_infinite_norm*v_minus_infinite_norm));

mu_VEN=324859;
e_minus_fun= @(r_P) 1 + (r_P * v_minus_infinite_norm^2)/mu_VEN;
delta_minus_fun= @(r_P) 2*asin(1./e_minus_fun(r_P));

e_plus_fun= @(r_P) 1 + (r_P * v_plus_infinite_norm^2)/mu_VEN;
delta_plus_fun= @(r_P) 2*asin(1./e_plus_fun(r_P));

delta= @(r_P) delta_minus_fun(r_P)./2 + delta_plus_fun(r_P)./2 ; %total turn angle
fun=@(r_P) delta(r_P) - turn_angle;

R_ven=6051.8;
h_atm=250;
R_VEN=R_ven+h_atm;

mass_VEN=4.867 * 1e24;
mass_SUN=1.989 * 1e30;
SOI_VEN=kep_2(1)*(mass_VEN/mass_SUN)^(2/5);
r_P_guess=10000;

r_P_find1=fzero(fun,r_P_guess);
control=fun(r_P_find1);

  if r_P_find1 < R_VEN   %&&  r_P_find > SOI_VEN
     r_P_find1=NaN;
 else
     r_P_find1=r_P_find1;
 end

v_P_minus=sqrt(v_minus_infinite_norm^2 + (2*mu_VEN/r_P_find1)); %vel at pericenter at entry arc
v_P_plus=sqrt(v_plus_infinite_norm^2 + (2*mu_VEN/r_P_find1)); %vel at percinter at exit arc

% delta_v to give at pericentre of flyby hyperbolas
delta_v_P1=abs(v_P_plus - v_P_minus);



%plot final orbit
hold on
plotOrbit(a_arr,e_arr,i_arr,OM_arr,om_arr,theta_0_2,theta_0_2+2*pi,dth,mu_Sun);




 %% PLOT THE HYPERBOLA ARCS OF THE FLYBY

u_flyby=cross(v_minus_infinite,v_plus_infinite)/norm(cross(v_minus_infinite,v_plus_infinite));
i_flyby=acos(u_flyby(3));
N_flyby=cross([0 ;0; 1],u_flyby);
N_flyby_norm=N_flyby/norm(N_flyby);
OM_flyby=acos(N_flyby_norm(1));
if N_flyby_norm(2) < 0
    OM_flyby=2*pi - acos(N_flyby_norm(1));
end
 
%a_T1,p_T1,e_T1,ERROR,v_i_T1,v_f_T1,TPAR_T1,THETA_T1
e_minus_fun_flyby= @(r_P) 1 + (r_P *v_minus_infinite_norm^2)/mu_VEN;
delta_minus_fun_flyby= @(r_P) 2*asin(1./e_minus_fun_flyby(r_P));
e_plus_fun_flyby= @(r_P) 1 + (r_P *v_plus_infinite_norm^2)/mu_VEN;
delta_plus_fun_flyby= @(r_P) 2*asin(1./e_plus_fun_flyby(r_P));

e_minus=e_minus_fun_flyby(r_P_find1);
e_plus=e_plus_fun_flyby(r_P_find1);
a_minus=-mu_VEN/v_minus_infinite_norm^2;
a_plus=-mu_VEN/v_plus_infinite_norm^2;
delta_minus=delta_minus_fun_flyby(r_P_find1);
delta_plus=delta_plus_fun_flyby(r_P_find1);


% PLOT HYPERBOLAS IN PLANETOCENTRIC FRAME
% INCOMING HYPERBOLA
beta_minus=(pi - delta_minus)/2;
i_minus=i_flyby;
om_minus=beta_minus;
OM_minus=OM_flyby;
theta_0_in=-0.5*pi;
theta_f_in=0;
dth=0.5*pi/180;


%PLOT VENUS
opts_example2.RefPlane = 'ecliptic';
figure;
planet3D('Venus',opts_example2);
grid on;

%PLOT INCOMING HYPERBOLA
if(theta_f_in<theta_0_in)
   theta_f_in=theta_0_in+2*pi; 
end
[rr_in] = par2car(a_minus,e_minus,i_minus,OM_minus,om_minus, theta_0_in,mu_VEN);
Posizioni_in=1000*[rr_in]; %from [Km] to [m]
for th=theta_0_in:dth:theta_f_in
[rr_in] = par2car(a_minus,e_minus,i_minus,OM_minus,om_minus, th,mu_VEN);
Posizioni_in=[Posizioni_in,1000*rr_in];
end
hold on
plot3(Posizioni_in(1,:),Posizioni_in(2,:),Posizioni_in(3,:),'r','LineWidth',2)
axis equal;

% OUTCOMING HYPERBOLA
beta_plus=(pi - delta_plus)/2;
i_plus=i_flyby;
om_plus=beta_plus;
OM_plus=OM_flyby;
theta_0_out=0;
theta_f_out=0.5*pi;
dth=0.5*pi/180;

%PLOT OUTCOMING HYPERBOLA
if(theta_f_out<theta_0_out)
   theta_f_out=theta_0_out+2*pi; 
end
[rr_out] = par2car(a_plus,e_plus,i_plus,OM_plus,om_plus, theta_0_out,mu_VEN);
Posizioni_out=1000*[rr_out]; %from [Km] to [m]
for th=theta_0_out:dth:theta_f_out
[rr_out] = par2car(a_plus,e_plus,i_plus,OM_plus,om_plus, th,mu_VEN);
Posizioni_out=[Posizioni_out,1000*rr_out];
end
hold on
plot3(Posizioni_out(1,:),Posizioni_out(2,:),Posizioni_out(3,:),'k','LineWidth',2)
hold on
quiver3(0,0,0,vv_0_VEN_flyby(1)*5e5,vv_0_VEN_flyby(2)*5e5,vv_0_VEN_flyby(3)*5e5,'g','LineWidth',1.5)

xlabel 'x [m]'
ylabel 'y [m]'
zlabel 'z [m]'

legend('','incoming hyperbola','outcoming hyperbola','Venus velocity')
title 'Flyby around Venus'

 %%
% %Period of the final orbit
% T=2*pi*sqrt(a_arr^3/mu_Sun)
% Days=T/86400;
% 
% figure()
% %Setting up an animation for the Earth, Venus and the trajectory
% mjd2000_dep=8530.5;
% mjd2000_arr_venus=8695.5;
% 
% %Compute the final eccentrical orbit
% dth=0.01;
% thf=theta_0_2+2*pi;
% [rr] = par2car(a_arr,e_arr,i_arr,OM_arr,om_arr, theta_0_2,mu_Sun);
% Posizioni=[rr];
% for th=theta_0_2:dth:thf
%    [rr] = par2car(a_arr,e_arr,i_arr,OM_arr,om_arr, th,mu_Sun);
%    Posizioni=[Posizioni,rr];
% end
% state=state_T1(:,1:3)
% state=state';
% TotTransfer=horzcat(state, Posizioni(1:3, 2:end));
% TotTransfer=TotTransfer';
% 
% [rrEarth]=par2car(kep_E_dep(1),kep_E_dep(2),kep_E_dep(3),kep_E_dep(4),kep_E_dep(5),kep_E_dep(6), mu_Sun)
% PosE=[rrEarth]
% for th=kep_E_dep(6):0.001:kep_E_dep(6)+2*pi
%     [rrEarth]=par2car(kep_E_dep(1),kep_E_dep(2),kep_E_dep(3),kep_E_dep(4),kep_E_dep(5),th, mu_Sun)
%     PosE=[PosE, rrEarth];
% end
% 
% [rrVen]=par2car(kep_venus(1),kep_venus(2),kep_venus(3),kep_venus(4),kep_venus(5),kep_venus(6), mu_Sun)
% PosV=[rrVen]
% for th=kep_venus(6):0.001:kep_venus(6)+2*pi
%     [rrVen]=par2car(kep_venus(1),kep_venus(2),kep_venus(3),kep_venus(4),kep_venus(5),th, mu_Sun)
%     PosV=[PosV, rrVen];
% end
% 
% comet3(state(1,:), state(2,:), state(3,:))
% hold on
% comet3(PosV(1,:), PosV(2,:), PosV(3,:))
% comet3(PosE(1,:), PosE(2,:), PosE(3,:))





% % plot Earth orbit
% figure
% plotOrbit(kep_E_dep(1),kep_E_dep(2),kep_E_dep(3),kep_E_dep(4),kep_E_dep(5),kep_E_dep(6),kep_E_dep(6)+2*pi,dth,mu_Sun);
% hold on
% %plot Earth
% plot3(rr_0_E(1),rr_0_E(2),rr_0_E(3), 'bo' ,'LineWidth',3)
% hold on
% 
% %plot venus orbit
% plotOrbit(kep_venus(1),kep_venus(2),kep_venus(3),kep_venus(4),kep_venus(5),kep_venus(6),kep_venus(6)+2*pi,dth,mu_Sun);
% %plot venus
% hold on
% plot3(rr_0_venus(1),rr_0_venus(2),rr_0_venus(3),'ro','LineWidth',3)




