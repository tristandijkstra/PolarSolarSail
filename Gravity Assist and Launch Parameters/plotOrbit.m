function plotOrbit(a,e,i,Om,om,th0,thf,dth,mu)

%      3D orbit plot
%
%--------Parametri di ingresso:
% a                Semiasse maggiore                    [Km]
% e                Eccentricità                         [-]
% i                Inclinazione                         [rad]
% Om               Ascensione retta del nodo ascendente [rad]
% om               Anomalia del pericentro              [rad]
% th0              Anomalia vera iniziale               [rad]
% thf              Anomalia vera finale                 [rad]
% dth              Passo anomalia vera                  [rad]
% mu               Costante gravitazionale planetaria   [Km^3/s^2]
%

if(thf<th0)
   thf=th0+2*pi; 
end
[rr] = par2car(a,e,i,Om,om, th0,mu);
Posizioni=[rr];
for th=th0:dth:thf
[rr] = par2car(a,e,i,Om,om, th,mu);
Posizioni=[Posizioni,rr];
end
%Terra3D(6378)
plot3(Posizioni(1,:),Posizioni(2,:),Posizioni(3,:),'linewidth',1.5)
hold on
%plot3(Posizioni(1,1),Posizioni(2,1),Posizioni(3,1),'o','MarkerEdgeColor','red')
axis equal;

