function  [delta_v_tot,delta_v_1,v_final_Lambert_1_x,v_final_Lambert_1_y,v_final_Lambert_1_z] = Lambert_optimal_venus(mjd_departure_in,mjd_departure_fin,mjd_arrival_in,mjd_arrival_fin,unit,Planet_1,Planet_2,orbitType,delta_v_1_MAX)

% MOST EFFICIENT ( MIN delta_v) LAMBERT TRANSFER BETWEEN A STARTING PLANET
% (Planet1) AND AN ARRIVAL PLANET (Planet2) GIVEN A DEPARTURE WINDOW AND AN
% ARRIVAL WINDOW COMPUTING ONLY delta_v_1 TO GO FROM PLANET 1 ORBIT TO
% THE LAMBERT LEG.

%NOTE. this function uses only a zero revolution transfer 
%NOTE: this function gives the CONSTRAIN on the Maximum excess velocity from launcher (constrain on delta_v1)


% PLANET'S LEGEND:
%                   1:   Mercury
%                   2:   Venus
%                   3:   Earth
%                   4:   Mars
%                   5:   Jupiter
%                   6:   Saturn
%                   7:   Uranus
%                   8:   Neptune
%                   9:   Pluto
%                   10:  Sun


% INPUT 
 % departure window in mjd: mjd_departure_in,mjd_departure_fin
 % arrival window in mjd: mjd_arrival_in,mjd_arrival_fin
 % time step between a point and the next in departure and arrival
 % windows , expressed in days : unit

 % the number related to starting planet : Planet1
 % the number related to arrival planet : Planet2
 % the options used in the Lambert solver:orbitType    Logical variable defining whether transfer is
 %                       0: direct transfer from R1 to R2 (counterclockwise)
 %                       1: retrograde transfer from R1 to R2 (clockwise)
 % the constrain on delta_v_1: delta_v_1_MAX


% OUTPUT
 % the delta_v_1 needed for the transfer for each departure point and each arrival point: delta_v_1
 % the components of the velocity at the end of the lambert arc: v_final_Lambert_1_x
 %                                                               v_final_Lambert_1_y
 %                                                               v_final_Lambert_1_z

 
% CONTRIBUTORS
 % Francesco Suffoletta


% DEPARTURE WINDOW

mjd_departure=[mjd_departure_in:unit:mjd_departure_fin];
% ARRIVAL WINDOW

mjd_arrival=[mjd_arrival_in:unit:mjd_arrival_fin];

% POSITION AND VELOCITY OF THE SPACECRAFT ON THE INITIAL ORBIT IN THE MOMENT OF DEPARTURE 

mjd2000_departure=mjd_departure - 51544.5; % mjd2000 day for departure
ibody_1_vect=Planet_1.*ones(1,length(mjd2000_departure));
for i=1:length(mjd2000_departure)
[kep_1(i,:),ksun(i)] = uplanet(mjd2000_departure(i), ibody_1_vect(i));  %orbital parameters wrt the starting point for orbit propagation
a_1(i)=kep_1(i,1);
e_1(i)=kep_1(i,2);
i_1(i)=kep_1(i,3);
OM_1(i)=kep_1(i,4);
om_1(i)=kep_1(i,5);
theta_0_1(i)=kep_1(i,6);

% initial cartesian parameters
[rr_0_1(i,:),vv_0_1(i,:)] = par2car(a_1(i),e_1(i),i_1(i),OM_1(i),om_1(i), theta_0_1(i),ksun(i));
end

% POSITION AND VELOCITY OF THE SPACECRAFT ON THE FINAL ORBIT IN THE MOMENT OF ARRIVAL

mjd2000_arrival=mjd_arrival - 51544.5; % mjd2000 day for arrival
ibody_2_vect=Planet_2.*ones(1,length(mjd2000_arrival));
for j=1:length(mjd2000_arrival)
[kep_2(j,:),ksun(j)] = uplanet(mjd2000_arrival(j), ibody_2_vect(j));  %orbital parameters wrt the arriving point for orbit propagation
a_2(j)=kep_2(j,1);
e_2(j)=kep_2(j,2);
i_2(j)=kep_2(j,3);
OM_2(j)=kep_2(j,4);
om_2(j)=kep_2(j,5);
theta_0_2(j)=kep_2(j,6);

% initial cartesian parameters
[rr_0_2(j,:),vv_0_2(j,:)] = par2car(a_2(j),e_2(j),i_2(j),OM_2(j),om_2(j), theta_0_2(j),ksun(j));
end

% TRANSFER ORBIT
% solve Lambert problem for the 2 given orbits, wrt the departure point
% and the arrival point. Compute delta_v needed to go on the transfer orbit
% from initial orbit on the departure point and compute delta_v needed to go
% from fransfer obit to final orbit on the arrival point

Nrev=0; %number of revolutions
optionsLMR=1;%warnings are displayed only when the algorithm does not converge


for i=1:length(mjd2000_departure)
    for j=1:length(mjd2000_arrival)
        delta_t_matrix(i,j)=(mjd2000_arrival(j) - mjd2000_departure(i))*24*60*60 ;% TIME OF FLIGTH [s]
%following the column departure day changes, arrival day is remains the same
%following the row arrival day changes, departure day is remains the same

%starting from each departing point (i) and arriving in each arrival
%point (j)
       [a_T_matr(i,j),p_T_matr(i,j),e_T_matr(i,j),ERROR_matr(i,j),v_i_matr(i,j,:),v_f_matr(i,j,:),TPAR_matr(i,j),THETA_matr(i,j)] = lambertMR(rr_0_1(i,:),rr_0_2(j,:),delta_t_matrix(i,j),ksun(i),orbitType,Nrev,optionsLMR);

       v_i_x=v_i_matr(:,:,1); 
       v_i_y=v_i_matr(:,:,2); 
       v_i_z=v_i_matr(:,:,3); 

  v_f_matr(abs(v_f_matr)>100) = NaN;


       %components of the delta_v1 matrix, delta_v needed to go from first
       %orbit to transfer orbit changhing departure and arrival point
       delta_v_1_x(i,j)=v_i_x(i,j) - vv_0_1(i,1);
       delta_v_1_y(i,j)=v_i_y(i,j) - vv_0_1(i,2);
       delta_v_1_z(i,j)=v_i_z(i,j) - vv_0_1(i,3);


       delta_v_1(i,j)=sqrt(delta_v_1_x(i,j).^2 + delta_v_1_y(i,j).^2 + delta_v_1_z(i,j).^2);

delta_v_1(delta_v_1>delta_v_1_MAX) = NaN;

         v_f_matr(abs(v_f_matr)>100) = NaN;

       v_final_Lambert_1_x=v_f_matr(:,:,1); 
       v_final_Lambert_1_y=v_f_matr(:,:,2); 
       v_final_Lambert_1_z=v_f_matr(:,:,3); 

       
       delta_v_tot(i,j)=delta_v_1(i,j);     %total delta_v needed for the mission, each column has a fixed arrival day and a changing departure day.
       %each row has a fixed departure day and a changing arrival day.


    end
end
    


% figure
% surf(datenum(mjd2000_departure),datenum(mjd2000_arrival),delta_v_tot')
% datetick('x','dd-mmm-yy'),datetick('y','dd-mmm-yy');
% colorbar
% xlabel 'Departure date from Earth'
% ylabel 'Arrival date on Saturn'
% zlabel 'delta_v_1 [km/s]'

% figure
% contour(datenum(mjd2000_departure),datenum(mjd2000_arrival),delta_v_tot',ShowText="on")
% datetick('x','dd-mmm-yy'),datetick('y','dd-mmm-yy');
% colorbar 
%  xlabel 'Departure date'
% ylabel 'Arrival date'
% zlabel 'delta_v [km/s]'
% grid on