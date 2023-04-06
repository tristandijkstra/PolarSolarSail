function [a, e, i,igrad, OMEGA,OMEGAgrad, omega,omegagrad, theta, thetagrad]=car2par(rr,vv,mu)
%rr=[3x1] vv=[3x1]
%NB angoli di acos sono in radianti, trasformo in gradi

%Function utile al calcolo dei parametri orbitali dati nel sistema
%geocentrico
%rr,vv sono grandezze vettoriali, r e v sono i moduli associati
r=sqrt(rr(1)^2+rr(2)^2+rr(3)^2);
v=sqrt(vv(1)^2+vv(2)^2+vv(3)^2);
%Ricordando la formula dell'energia e sapendo che assume valore -mu/2a al
%pericentro, restando costante otteniamo il primo parametro, ovvero il
%semiasse maggiore
Energia=1/2*v^2- mu/r;
a=-mu/(2*Energia);
%Calcoliamo ora il momento angolare (vettore e modulo)
hh=cross(rr,vv);
h=norm(hh);
%Calcolo del vettore e modulo eccentricità
ee= (cross(vv,hh)/mu) - (rr/r);
e=norm(ee);
%Il vettore inclinazione può essere calcolato come il coseno fra k ed h
k=[0;0;1];
i=acos(hh(3)/h);
igrad= 180*i/pi;
%Calcolo Linea dei nodi
N= (cross(k,hh)/norm(cross(k,hh)));
%Calcolo retta nodo ascendente
I=[1,0,0];
j=[0,1,0];
if N(2)>=0
    OMEGA=acos(N(1));
else OMEGA= 2*pi- acos(N(1));
end
OMEGAgrad=180*OMEGA/pi;
%Calcolo anomalia del pericentro
if ee(3)>=0
    omega=acos(dot(N,ee)/e);
else
    omega=2*pi-acos(dot(N,ee)/e);
end
omegagrad=180*omega/pi;
%Calcolo infine l'anomalia vera
vr=dot(vv,rr)/r;
if vr>0
    theta=acos(dot(rr,ee)/(r*e));
else
    theta=2*pi-acos(dot(rr,ee)/(r*e));
end

thetagrad=theta*180/pi;



