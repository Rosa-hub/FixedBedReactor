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
lamse=2;                                                   %[W/m/K]
BVF=0.38+0.073*(1+(dt/dp-2)^2/(dt/dp)^2);                       %[-]
av=6*(1-BVF)/dp;                                                %[m2/m3]
rob=ros*(1-BVF);

%Proces:
U=270;                                                          %[W/m2/K]
Cps=1250;                                                       %[J/kg/K]
hf=550;                                                         %[W/m2/K]
lamfe=2;                                                        %[W/m/K]
u=1.5;                                                          %[m/s]
XF=[0.5 0 0.07 0 0 0.43];                                       %[-]
TF=225+273.15;                                                  %[K]
PF=1.2e6;                                                       %[Pa]
Tc=480;                                                         %[K]

f=Property;
f.X=XF;
f.P=PF;
f.T=TF;

k=length(XF);
VF=u*pi*dt^2/4;
mt=PF*VF/8.314/TF*f.MWm/1000;

nO_0=PF*VF/8.314/TF*XF(3);
N=20;
M=zeros(k*(N+1));
for i=1:k
   M(i,i)=1; 
end
opt=odeset('Mass',M,'Events',@Eventfun);

z=[0 1000];

Ts=linspace(550,500,N);
s=Property;
s.P=PF;
s.T=Ts;
Cts=s.Cm;
cAi=linspace(0.1,XF(1)*0.5,N).*Cts;
cBi=linspace(0,0.15,N).*Cts;
cCi=linspace(0,0.07,N).*Cts;
cDi=linspace(0.01,0.1,N).*Cts;
cEi=linspace(0.01,0.1,N).*Cts;
IC=[XF(1:k-1)*f.Cm TF [cAi cBi cCi cDi cEi] Ts];

[z,C]=ode15s(@heterModel,z,IC,opt,PF,mt,av,U,Cps,hf,lamfe,lamse,Tc,BVF,rob,dt,dp,Por,Tort,N,nO_0);
cAf=C(:,1);
cBf=C(:,2);
cCf=C(:,3);
cDf=C(:,4);
cEf=C(:,5);
Tf=C(:,6);
cAp=C(:,7:N+6);
cBp=C(:,N+7:2*N+6);
cCp=C(:,2*N+7:3*N+6);
cDp=C(:,3*N+7:4*N+6);
cEp=C(:,4*N+7:5*N+6);
Tp=C(5*N+7:6*N+6);

for i=1:length(z)
f=Property;
f.P=PF;
f.T=Tf(i);
Ct=f.Cm;
cFf=Ct-cAf(i)-cBf(i)-cCf(i)-cDf(i)-cEf(i);
f.X=[cAf(i) cBf(i) cCf(i) cDf(i) cEf(i) cFf]./Ct;
nO=mt/f.MWm*f.X(3)*1000;
XC(i)=(nO_0-nO)/nO_0;

end

z_50=interp1(XC,z,0.5);
cApit=interp1(z,cAp,z_50);
cBpit=interp1(z,cBp,z_50);
cCpit=interp1(z,cCp,z_50);
cDpit=interp1(z,cDp,z_50);
cEpit=interp1(z,cEp,z_50);
Tpit=interp1(z,Tp,z_50);

k1=70.4.*exp(-59860/8.314./Tp);
k2=49.4e3.*exp(-89791/8.314./Tp);

ksi1=k1.*cCpit;
ksi2=k2.*cCpit;

ksis=


figure(1)
plot(z,cAf,z,cBf,z,cCf,z,cDf,z,cEf)
legend('A','B','C','D','E')
figure(2)
plot(z,Tf)
r=linspace(0,dp/2,N);




function dCdz=heterModel(z,X,P,mt,av,U,Cps,hf,lamfe,lamse,Tc,BVF,rob,dt,dp,Por,Tort,N,nO)

cAf=X(1);
cBf=X(2);
cCf=X(3);
cDf=X(4);
cEf=X(5);


T=X(6);
%P=X(7);
cAp=X(7:N+6);
cBp=X(N+7:2*N+6);
cCp=X(2*N+7:3*N+6);
cDp=X(3*N+7:4*N+6);
cEp=X(4*N+7:5*N+6);


Ts=X(5*N+7:6*N+6);
%Ps=(P+PF)/2;

f=Property;
f.T=T;
f.P=P;
Cf=f.Cm;
cFf=Cf-cAf-cBf-cCf-cDf-cEf;
f.X=[cAf cBf cCf cDf cEf cFf]/Cf;
V=mt*8.314*T/P/f.MWm*1000;
u=4*V/pi/dt^2;

kgi=f.kg(hf,lamfe);
Cpg=f.Cpm/f.MWm*1000;
% Re=f.ReN(u,dp,BVF);
% fp=(1-BVF)/BVF^3*(1.75+150*(1-BVF)/Re);
% dP=fp*z/dp*f.Densm*u^2*(1-BVF)/BVF^3;
%G=PF-dP-P;

dcAf=-kgi(1)*av*(cAf-cAp(end))/u;
dcBf=-kgi(2)*av*(cBf-cBp(end))/u;
dcCf=-kgi(3)*av*(cCf-cCp(end))/u;
dcDf=-kgi(4)*av*(cDf-cDp(end))/u;
dcEf=-kgi(5)*av*(cEf-cEp(end))/u;

dT=(-4*U/dt*(T-Tc)+hf*av*(Ts(end)-T))/(u*f.Densm*Cpg);

F=[dcAf dcBf dcCf dcDf dcEf dT];

k1=70.4.*exp(-59860/8.314./Ts);
k2=49.4e3.*exp(-89791/8.314./Ts);
dr=dp/2/(N-1);
r=linspace(0,dp/2,N);
%Def=f.Dim.*Por./Tort;
for i=1:N
    s=Property;
    s.T=Ts(i);
    s.P=P;
    cFp=s.Cm-cAp(i)-cBp(i)-cCp(i)-cDp(i)-cEp(i);
    s.X=[cAp(i) cBp(i) cCp(i) cDp(i) cEp(i) cFp]./s.Cm;
    [dhr1,dhr2]=s.dHr;
    dhr1=dhr1/2;
    dhr2=dhr2*3;
    ksi1=k1.*cCp(i);
    ksi2=k2.*cCp(i);
    Def=s.Dim.*Por./Tort;
    
    
    if i==1
        A(i)=Def(1)*((2*cAp(i+1)-2*cAp(i))/dr^2)+rob*(s.v1(1)*ksi1(i)+s.v2(1)*ksi2(i));
        B(i)=Def(2)*((2*cBp(i+1)-2*cBp(i))/dr^2)+rob*(s.v1(2)*ksi1(i)+s.v2(2)*ksi2(i));
        C(i)=Def(3)*((2*cCp(i+1)-2*cCp(i))/dr^2)+rob*(s.v1(3)*ksi1(i)+s.v2(3)*ksi2(i));
        D(i)=Def(4)*((2*cDp(i+1)-2*cDp(i))/dr^2)+rob*(s.v1(4)*ksi1(i)+s.v2(4)*ksi2(i));
        E(i)=Def(5)*((2*cEp(i+1)-2*cEp(i))/dr^2)+rob*(s.v1(5)*ksi1(i)+s.v2(5)*ksi2(i));
       
        Tp(i)=(lamse*((2*Ts(i+1)-2*Ts(i))/dr^2)+rob*((-dhr1)*ksi1(i)+(-dhr2)*ksi2(i)))/1000;
    elseif i==N
        
        A(i)=(Def(1)*(cAp(i)-cAp(i-1))/dr+kgi(1)*(cAp(i)-cAf));
        B(i)=(Def(2)*(cBp(i)-cBp(i-1))/dr+kgi(2)*(cBp(i)-cBf));
        C(i)=(Def(3)*(cCp(i)-cCp(i-1))/dr+kgi(3)*(cCp(i)-cCf));
        D(i)=(Def(4)*(cDp(i)-cDp(i-1))/dr+kgi(4)*(cDp(i)-cDf));
        E(i)=(Def(5)*(cEp(i)-cEp(i-1))/dr+kgi(5)*(cEp(i)-cEf));
    
        Tp(i)=(lamse*(Ts(i)-Ts(i-1))/dr+hf*(Ts(i)-T));
    else
        
        A(i)=Def(1)*((cAp(i-1)-2*cAp(i)+cAp(i+1))/dr^2+2/r(i)*(cAp(i)-cAp(i-1))/dr)+rob*(s.v1(1)*ksi1(i)+s.v2(1)*ksi2(i));
        B(i)=Def(2)*((cBp(i-1)-2*cBp(i)+cBp(i+1))/dr^2+2/r(i)*(cBp(i)-cBp(i-1))/dr)+rob*(s.v1(2)*ksi1(i)+s.v2(2)*ksi2(i));
        C(i)=Def(3)*((cCp(i-1)-2*cCp(i)+cCp(i+1))/dr^2+2/r(i)*(cCp(i)-cCp(i-1))/dr)+rob*(s.v1(3)*ksi1(i)+s.v2(3)*ksi2(i));
        D(i)=Def(4)*((cDp(i-1)-2*cDp(i)+cDp(i+1))/dr^2+2/r(i)*(cDp(i)-cDp(i-1))/dr)+rob*(s.v1(4)*ksi1(i)+s.v2(4)*ksi2(i));
        E(i)=Def(5)*((cEp(i-1)-2*cEp(i)+cEp(i+1))/dr^2+2/r(i)*(cEp(i)-cEp(i-1))/dr)+rob*(s.v1(5)*ksi1(i)+s.v2(5)*ksi2(i));
        
        Tp(i)=(lamse*((Ts(i-1)-2*Ts(i)+Ts(i+1))/dr^2+2/r(i)*(Ts(i)-Ts(i-1))/dr)+rob*((-dhr1)*ksi1(i)+(-dhr2)*ksi2(i)))/1000;

    end
end

dCdz=[F A B C D E Tp]';


function [val,it,dir]=Eventfun(z,X,P,mt,av,U,Cps,hf,lamfe,lamse,Tc,BVF,rob,dt,dp,Por,Tort,N,nO_0)
cAf=X(1);
cBf=X(2);
cCf=X(3);
cDf=X(4);
cEf=X(5);

T=X(6);

f=Property;
f.P=P;
f.T=T;
Ct=f.Cm;
cFf=Ct-cAf-cBf-cCf-cDf-cEf;
f.X=[cAf cBf cCf cDf cEf cFf]./Ct;
nO=mt/f.MWm*f.X(3)*1000;
XO=(nO_0-nO)/nO_0;

val=0.75-XO;
it=1;
dir=0;