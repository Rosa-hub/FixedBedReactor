classdef Property
    
    properties
        P{mustBeNumeric}
        T{mustBeNumeric}
        X{mustBeNumeric}   
    end
    
    properties(Constant)
      M=[28.0532 44.5026 31.9988 44.0095 18.0153 16.0425];
      A=[-6.38788 -23.25802	31.32234 24.99735 30.092 -0.703029];
      B=[184.4019 275.6997 -20.23531 55.18696 6.832514 108.4773];
      C=[-112.9718 -188.9729 57.86644 -33.69137	6.793435 -42.52157];
      D=[28.49593 51.0335 -36.50624	7.948387 -2.53448 5.862788];
      E=[0.31554 0.38693 -0.007374 -0.136638 0.082139 0.678565];
      Pc=[50.6 71.91 50.43 73.4 220.64 46.1];
      Tc=[282.5	469	154.58 304.35 647 190.6];
      Vc=[0.1311 0.1403 0.0735 0.0919 0.0559 0.0986];
      R=8.314;
      v1=[-2 2 -1 0 0 0];
      v2=[-1/3 0 -1 2/3 2/3 0];
      hf=[52.4 -52.63 0 -393.52 -241.83	-74.87];
      dm=[0	1.94 0 0 1.8546	0];
      
    end
    
    methods
        function output=Cpm(obj)
            t=obj.T/1000;
            a=sum(obj.A.*obj.X);
            b=sum(obj.B.*obj.X);
            c=sum(obj.C.*obj.X);
            d=sum(obj.D.*obj.X);
            e=sum(obj.E.*obj.X);
            
            output=a+b*t+c*t^2+d*t^3+e/t^2;
        end
        
        function output=Densm(obj)
            
            output=obj.P*obj.MWm/1000/obj.R/obj.T; 
        end
        
        function output=MWm(obj)
            output=sum(obj.M.*obj.X);
        end
        
        function output=Cm(obj)
            output=obj.P/obj.R/obj.T;
        end
        
        function [output1,output2,output3,output4]=visc(obj)
           mir=52.46.*obj.dm.^2.*obj.Pc./obj.Tc.^2; 
           Zc=obj.Pc.*obj.Vc./obj.Tc./obj.R;
           Tr=obj.T./obj.Tc;
           
           for i=1:length(mir)
           if mir(i)>=0 && mir(i)<0.022
               Fp0(i)=1;
           elseif mir(i)>=0.075
               Fp0(i)=1+30.55.*(0.292-Zc(i)).^1.72.*abs(0.96+0.1.*(Tr(i)-0.7));
               
           else
               Fp0(i)=1+30.55.*(0.292-Zc(i)).^1.72;
           end
           end
           FQ0=1;
           Eps=0.176.*(obj.Tc./obj.M.^3./obj.Pc.^4).^(1/6);
           viski_LP=(0.807.*Tr.^0.618-0.357.*exp(-0.449.*Tr)+0.340.*exp(-4.058.*Tr)+0.018).*Fp0.*FQ0./Eps;
           
           Tcm=sum(obj.X.*obj.Tc);
           Pcm=obj.R*Tcm*sum(obj.X.*Zc)/sum(obj.X.*obj.Vc);
           Mm=obj.MWm;
           Fp0m=sum(obj.X.*Fp0);
           
           [Mh,h]=max(obj.M);
           [Ml,l]=min(obj.M);
           
           if Mh/Ml>9 && obj.X(h)>0.05 && obj.X(h)<0.7
               a=1-0.01*(Mh/Ml)^0.87;
           else
               a=1;
           end
           
           FQ0m=sum(obj.X*FQ0)*a;
           Trm=obj.T/Tcm;
           Epsm=0.176.*(Tcm./Mm.^3./Pcm.^4).^(1/6);
           
           viskm_LP=(0.807.*Trm.^0.618-0.357.*exp(-0.449.*Trm)+0.340.*exp(-4.058.*Trm)+0.018).*Fp0m.*FQ0m./Epsm;
           
           Pr=obj.P*1e-5./obj.Pc;        
           Z1=viski_LP.*Eps;
           
           for i=1:length(Pr)
               if Tr(i)>1 && Tr(i)<40 && Pr(i)>0 && Pr(i)<100
                   a1=1.245e-3;b1=1.6553;c1=0.4489;d1=1.7368;f1=0.9425;
                   a2=5.1726;b2=1.2723;c2=3.0578;d2=2.2310;f2=-0.1853;
                   gam=-0.3286;del=-37.7332;ep=-7.6351;xsi=0.4489;

                   a=a1./Tr(i).*exp(a2.*Tr(i).^gam);
                   b=a.*(b1.*Tr(i)-b2);
                   c=c1./Tr(i).*exp(c2.*Tr(i).^del);
                   d=d1./Tr(i).*exp(d2.*Tr(i).^ep);
                   e=1.3088;
                   f=f1.*exp(f2.*Tr(i).^xsi);

                   Z2(i)=Z1(i)*(1+a*Pr(i)^e/(b*Pr(i)^f+(1+c*Pr(i)^d)^-1));
           
               else
                   a=3.262+14.98*Pr(i)^5.508;
                   b=1.390+5.746*Pr(i);
                   Z2(i)=0.600+0.760*Pr(i)^a+(6.990*Pr(i)^b-0.6)*(1-Tr(i));
               end
           end
           Y=Z2./Z1;
           Fp=(1+(Fp0-1).*Y.^-3)./Fp0;
           FQ=(1+(FQ0-1).*(Y.^-1-0.007.*(log(Y)).^4))./FQ0;
           viski_HP=Z2.*Fp.*FQ./Eps;
           
           Prm=obj.P*1e-5/Pcm;
           Z1m=viskm_LP.*Epsm;
           
           if Trm>1 && Trm<40 && Prm>0 && Prm<100
                   a1=1.245e-3;b1=1.6553;c1=0.4489;d1=1.7368;f1=0.9425;
                   a2=5.1726;b2=1.2723;c2=3.0578;d2=2.2310;f2=-0.1853;
                   gam=-0.3286;del=-37.7332;ep=-7.6351;xsi=0.4489;

                   a=a1./Trm.*exp(a2.*Trm.^gam);
                   b=a.*(b1.*Trm-b2);
                   c=c1./Trm.*exp(c2.*Trm.^del);
                   d=d1./Trm.*exp(d2.*Trm.^ep);
                   e=1.3088;
                   f=f1.*exp(f2.*Trm.^xsi);

                   Z2m=Z1m*(1+a*Prm^e/(b*Prm^f+(1+c*Prm^d)^-1));
           
               else
                   a=3.262+14.98*Prm^5.508;
                   b=1.390+5.746*Prm;
                   Z2m=0.600+0.760*Prm^a+(6.990*Prm^b-0.6)*(1-Trm);
           end
               
           Ym=Z2m./Z1m;
           Fpm=(1+(Fp0m-1).*Ym.^-3)./Fp0m;
           FQm=(1+(FQ0m-1).*(Ym.^-1-0.007.*(log(Ym)).^4))./FQ0m;
           viskm_HP=Z2m.*Fpm.*FQm./Epsm;
           
           output1=viski_LP*1e-3;
           output2=viskm_LP*1e-3;
           output3=viski_HP*1e-3;
           output4=viskm_HP*1e-3;
        end
        
        function [output1,output2]=dHr(obj)
            T0=298.15;
            dHr1=sum(obj.v1.*obj.hf);
            dHr2=sum(obj.v2.*obj.hf);
            
            if obj.T~=T0
                t=obj.T/1000;
                t0=T0/1000;
                
              a1=sum(obj.v1.*obj.A);
              b1=sum(obj.v1.*obj.B);
              c1=sum(obj.v1.*obj.C);
              d1=sum(obj.v1.*obj.D);
              e1=sum(obj.v1.*obj.E);
              
              a2=sum(obj.v2.*obj.A);
              b2=sum(obj.v2.*obj.B);
              c2=sum(obj.v2.*obj.C);
              d2=sum(obj.v2.*obj.D);
              e2=sum(obj.v2.*obj.E);
              
              dHr1=dHr1+(a1*(t-t0)+b1/2*(t-t0)^2+c1/3*(t-t0)^3+d1/4*(t-t0)^4-e1/(t-t0));
              dHr2=dHr2+(a2*(t-t0)+b2/2*(t-t0)^2+c2/3*(t-t0)^3+d2/4*(t-t0)^4-e2/(t-t0));
                
            end
            
            output1=dHr1;
            output2=dHr2;
        end
    end
end