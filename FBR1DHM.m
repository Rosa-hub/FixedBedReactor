function FBR1DHM
clc
clear all;
close all;
%Riešenie FB reaktora na oxidáciu etylénoxidu pomocou 1 rozmerného
%heterogénneho modelu s medzi?astcovým a vnútro?asticovým gradientom
%A-C2H4
%B-C2H4O
%C-O2
%D-CO2
%E-H2O
%F-CH4
%% ---///---Parametre---///--- 
%Reaktor:
dt=0.04;                                                        %[m]
dp=2.5e-3;                                                      %[m]
Por=0.42;                                                       %[-]
Tort=2.74;                                                      %[-]
ros=1750;                                                       %[kg/m3]
lamse=6.0535;                                                   %[W/m/K]
BVF=0.38+0.073*(1+(dt/dp-2)^2/(dt/dp)^2);                       %[-]
av=6*(1-BVF)/dp;                                                %[m2/m3]

%Proces:
U=270;                                                          %[W/m2/K]
Cps=1250;                                                       %[J/kg/K]
hf=550;                                                         %[W/m2/K]
lamfe=2;                                                        %[W/m/K]
u=1.5;                                                          %[m/s]
XF=[0.5 0 0.07 0 0];                                            %[-]
TF=225+273.15;                                                  %[K]
PF=1.2e6;                                                       %[Pa]
Tc=480;                                                         %[K]

N=20;
M=zeros(7+5*N);
for i=1:7
   M(i,i)=1; 
end
opt=odeset('Mass',M,'Events',@Eventfun);

z=[20 50];
xAi=ones(1,N)*0.25;
xBi=ones(1,N)*0.15;
xCi=ones(1,N)*0.1;
xDi=ones(1,N)*0.15;
xEi=ones(1,N)*0.15;
Ts=ones(1,N)*520;
IC=[XF TF xAi xBi xCi xDi xEi Ts];

[z,C]=ode45(@heterModel,z,IC,opt,PF,av,U,Cps,hf,lamfe,lamse,u,Tc,BVF,ros,dt,dp,Por,Tort,N);

function dXdz=heterModel(z,X,P,av,U,Cps,hf,lamfe,lamse,u,Tc,BVF,ros,dt,dp,Por,Tort,N)

xAf=X(1);
xBf=X(2);
xCf=X(3);
xDf=X(4);
xEf=X(5);
xFf=1-xAf-xBf-xCf-xDf-xEf;

T=X(6);

xAp=X(7:N+6);
xBp=X(N+7:2*N+6);
xCp=X(2*N+7:3*N+6);
xDp=X(3*N+7:4*N+6);
xEp=X(4*N+7:5*N+6);
xFp=1-xAp-xBp-xCp-xDp-xEp;

Ts=X(5*N+7:6*N+6);

f=Property;
f.T=T;
f.P=P;
f.X=[xAf xBf xCf xDf xEf xFf];
[A,B,C,D]=f.visc
