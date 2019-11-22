function FBR1DHM
%% ---///---Parametre---///--- 
%Reaktor:
dr=0.04;                                                        %[m]
dp=2.5e-3;                                                      %[m]
Por=0.42;                                                       %[-]
Tort=2.74;                                                      %[-]
ros=1750;                                                       %[kg/m3]
BVF=0.390+1.740/(dr/dp+1.140)^2;                                %[-]

%Proces:
U=270;                                                          %[W/m2/K]
Cps=1250;                                                       %[J/kg/K]
hf=550;                                                         %[W/m2/K]
lame=2;                                                         %[W/m/K]
u=1.5;                                                          %[m/s]
xF=[0.5 0 0.07 0 0 0.43];                                       %[-]
TF=225+273.15;                                                  %[K]
P=1.2e6;                                                        %[Pa]
Tc=480;                                                         %[K]