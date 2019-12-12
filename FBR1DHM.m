function FBR1DHM
clc
clear all;
close all;
%Riešenie FB reaktora na oxidáciu etylénoxidu pomocou 1 rozmerného
%heterogénneho modelu s medzicastcovým a vnútrc?asticovým gradientom
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
XF=[0.5 0 0.5 0 0 0.0];                                         %[-]
TF=225+273.15;                                                  %[K]
PF=1.2e6;                                                       %[Pa]
Tc=480;                                                         %[K]

%Inicializácia vlastností pre vstup suroviny
f=Property;
f.X=XF;
f.P=PF;
f.T=TF;
CF=XF*f.Cm;
k=length(XF);
VF=u*pi*dt^2/4;
mt=PF*VF/8.314/TF*f.MWm/1000;
nO_0=PF*VF/8.314/TF*XF(3);
N=20;

%Zadefinovanie matice pre výpo?et sústavy DAE
M=zeros(k*(N+1));
for i=1:k
   M(i,i)=1; 
end
opt=odeset('Mass',M,'Events',@Eventfun);

z=[0 1000];

%Definovanie po?iato?ných podmienok a nástrelov
Ts=linspace(600,500,N);
s=Property;
s.P=PF;
s.T=Ts;
Cts=s.Cm;
cAi=linspace(0.1,XF(1)*0.4,N).*Cts;
cBi=linspace(0,0.2,N).*Cts;
cCi=linspace(0,0.05,N).*Cts;
cDi=linspace(0.01,0.1,N).*Cts;
cEi=linspace(0.01,0.1,N).*Cts;
IC=[CF(1:k-1) TF [cAi cBi cCi cDi cEi] Ts];

[z,C]=ode15s(@heterModel,z,IC,opt,PF,mt,av,U,Cps,hf,lamfe,lamse,Tc,BVF,ros,dt,dp,Por,Tort,N,nO_0);
z=real(z);
C=real(C);
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
Tp=C(:,5*N+7:6*N+6);

%Spracovanie výstupných hodnôt
for i=1:length(z)
f=Property;
f.P=PF;
f.T=Tf(i);
Ct=f.Cm;
cFf(i)=Ct-cAf(i)-cBf(i)-cCf(i)-cDf(i)-cEf(i);
f.X=[cAf(i) cBf(i) cCf(i) cDf(i) cEf(i) cFf(i)]./Ct;
nO=mt/f.MWm*f.X(3)*1000;
XC(i)=(nO_0-nO)/nO_0;

end

z_50=interp1(XC,z,0.5);
cCpit=interp1(z,cCp,z_50);
Tpit=interp1(z,Tp,z_50);

k1_50=70.4.*exp(-59860/8.314./Tpit);
k2_50=49.4e3.*exp(-89791/8.314./Tpit);
k1_0=70.4.*exp(-59860/8.314./Tp(1,:));
k2_0=49.4e3.*exp(-89791/8.314./Tp(1,:));

ksi1_50=k1_50.*cCpit;
ksi2_50=k2_50.*cCpit;
ksi1_0=k1_0.*cCp(1,:);
ksi2_0=k2_0.*cCp(1,:);

%Výpo?et faktora ú?innosti
r=linspace(0,dp/2,N);

I1_50=0;
I2_50=0;
I1_0=0;
I2_0=0;
for i=2:N
    I1_50=I1_50+(ksi1_50(i)*r(i)^2+ksi1_50(i-1)*r(i-1)^2)/2*(r(i)-r(i-1));
    I2_50=I2_50+(ksi2_50(i)*r(i)^2+ksi2_50(i-1)*r(i-1)^2)/2*(r(i)-r(i-1));
    I1_0=I1_0+(ksi1_0(i)*r(i)^2+ksi1_0(i-1)*r(i-1)^2)/2*(r(i)-r(i-1));
    I2_0=I2_0+(ksi2_0(i)*r(i)^2+ksi2_0(i-1)*r(i-1)^2)/2*(r(i)-r(i-1));
    
end
Uc1_50=3/r(end)^3*I1_50/ksi1_50(end)*100
Uc2_50=3/r(end)^3*I2_50/ksi2_50(end)*100
Uc1_0=3/r(end)^3*I1_0/ksi1_0(end)*100
Uc2_0=3/r(end)^3*I2_0/ksi2_0(end)*100

Sel=(cBf(end)-CF(2))/f.v1(2)*f.v1(3)/(cCf(end)-CF(3))*100

%Vykreslenie grafov
figure(1)
plot(z,cAf,z,cBf,z,cCf,z,cDf,z,cEf,'LineWidth',2)
legend('A','B','C','D','E')
grid on
title 'Koncentracný profil v reaktore'
xlabel 'L [m]'
ylabel 'c [mol/m^3]'

figure(2)
plot(z,Tf,'LineWidth',2)
title 'Teplotný profil v reaktore'
xlabel 'L [m]'
ylabel 'T [K]'
grid on


disp(z(end))
disp(z(end)/dt)


function dCdz=heterModel(z,X,P,mt,av,U,Cps,hf,lamfe,lamse,Tc,BVF,rob,dt,dp,Por,Tort,N,nO)

cAf=X(1);
cBf=X(2);
cCf=X(3);
cDf=X(4);
cEf=X(5);


T=X(6);
cAp=X(7:N+6);
cBp=X(N+7:2*N+6);
cCp=X(2*N+7:3*N+6);
cDp=X(3*N+7:4*N+6);
cEp=X(4*N+7:5*N+6);


Ts=X(5*N+7:6*N+6);

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
Def=f.Dim.*Por./Tort;
for i=1:N
    s=Property;
    s.T=Ts(i);
    s.P=P;
    cFp=s.Cm-cAp(i)-cBp(i)-cCp(i)-cDp(i)-cEp(i);
    s.X=[cAp(i) cBp(i) cCp(i) cDp(i) cEp(i) cFp]./s.Cm;
    [dhr1,dhr2]=s.dHr;
    dhr1=dhr1;
    dhr2=dhr2;
    ksi1=k1.*cCp(i);
    ksi2=k2.*cCp(i);
    
    
    
    if i==1
        A(i)=Def(1)*((2*cAp(i+1)-2*cAp(i))/dr^2)+rob*(s.v1(1)*ksi1(i)+s.v2(1)*ksi2(i));
        B(i)=Def(2)*((2*cBp(i+1)-2*cBp(i))/dr^2)+rob*(s.v1(2)*ksi1(i)+s.v2(2)*ksi2(i));
        C(i)=Def(3)*((2*cCp(i+1)-2*cCp(i))/dr^2)+rob*(s.v1(3)*ksi1(i)+s.v2(3)*ksi2(i));
        D(i)=Def(4)*((2*cDp(i+1)-2*cDp(i))/dr^2)+rob*(s.v1(4)*ksi1(i)+s.v2(4)*ksi2(i));
        E(i)=Def(5)*((2*cEp(i+1)-2*cEp(i))/dr^2)+rob*(s.v1(5)*ksi1(i)+s.v2(5)*ksi2(i));
       
        Tp(i)=(lamse*((2*Ts(i+1)-2*Ts(i))/dr^2)+rob*((-dhr1)*ksi1(i)+(-dhr2)*ksi2(i)))/10000;
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
        
        Tp(i)=(lamse*((Ts(i-1)-2*Ts(i)+Ts(i+1))/dr^2+2/r(i)*(Ts(i)-Ts(i-1))/dr)+rob*((-dhr1)*ksi1(i)+(-dhr2)*ksi2(i)))/10000;

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