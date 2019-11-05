function [I,Imetal,Iglass]=RadioMetryOfSurfaces(D,r,ng,v,thetai,phii,m,n,k)

theta = zeros(D);
phi = zeros(D);
L=[sin(thetai)*cos(phii);sin(thetai)*sin(phii);cos(thetai)];
Z=L+v;
Normm=norm(Z);
H=Z./Normm;
for i=1:D
    for j=1:D
        dd=((i-D/2)^2 +(j-D/2)^2);
        d=sqrt(dd);
        if (d < r)
            theta=asin(d/r);
            phi=atan2((j-D/2) , (i-D/2));
            
            N=[sin(theta).*sin(phi),sin(theta).*cos(phi),cos(theta)];
            
            theta_h = acos(dot(N,H));
            theta_i =acos(dot(L,N));
            theta_iprim =acos(dot(L,H));
            theta_r = acos(dot(v,N));
            theta_g = acos(dot(v,L));
            %Nh=cos(theta_h);
            Nh=mean(dot(N,H));
            Li=cos(theta_i);
            vr=cos(theta_r);
            vg=cos(theta_g);
            Th=tan(theta_h);
            Ci=cot(theta_i);
            Cr=cot(theta_r);
            Cor=mean(dot(v,N));
            %Cor=cos(theta_r);
            %Coi=cos(theta_i);
            Coi=mean(dot(L,N));
            %CoP=cos(theta_iprim);
            CoP=mean(dot(L,H));
            TAnP=tan(theta_iprim);
            TAni=tan(theta_i);
            CP=cot(theta_iprim);
            SINp=sin(theta_iprim);
            SINi=sin(theta_i);
            Dh=exp(-(Th^2)/(2*m^2))/(2*pi*m^2*Nh^3);
            Atheta=(sqrt(2)*m*exp(-Ci^2/(2*m^2))/(sqrt(pi)*Ci)-erfc(Ci/(sqrt(2)*m)))/2;
            AthetaR=(sqrt(2)*m*exp(-Cr^2/(2*m^2))/(sqrt(pi)*Cr)-erfc(Cr/(sqrt(2)*m)))/2;
            z=n^2-k^2-SINp^2;
            a=sqrt(z^2+4*n*k);
            zG=ng^2-SINi^2;
            aG=zG;
            
            %RsP=((sqrt(a+z)-sqrt(2)*Coi)^2+a-z)/((sqrt(a+z)+sqrt(2)*Coi)^2+a-z);
            %RpP=RsP*((sqrt(a+z)-sqrt(2)*SINi*TAni)^2+a-z)/((sqrt(a+z)+sqrt(2)*SINi*TAni)^2+a-z);
            
            RsP=((sqrt(a+z)-sqrt(2)*CoP)^2+a-z)/((sqrt(a+z)+sqrt(2)*CoP)^2+a-z);
            RpP=RsP*((sqrt(a+z)-sqrt(2)*SINp*TAnP)^2+a-z)/((sqrt(a+z)+sqrt(2)*SINp*TAnP)^2+a-z);
            
            RsPG=((sqrt(2*zG)-sqrt(2)*CoP)^2)/((sqrt(2*zG)+sqrt(2)*CoP)^2);
            RpPG=RsP*((sqrt(2*zG)-sqrt(2)*SINp*TAnP)^2)/((sqrt(2*zG)+sqrt(2)*SINp*TAnP)^2);
            
            %RsPG=((CoP-sqrt(ng^2-SINp^2))/(CoP+sqrt(ng^2-SINp^2)))^2;
            %RpPG=(((ng^2)*CoP-sqrt(ng^2-SINp^2))/((ng^2)*CoP+sqrt(ng^2-SINp^2)))^2;
            
            Rtotal=(RsP+RpP)/2;
            RtotalG=(RsPG+RpPG)/2;
            if (CoP > 0)
                gTheta(i,j)=1/ (Atheta+1);
                gThetaR(i,j)=1/ (AthetaR+1);
                G(i,j)=(1/( Atheta+1)).*(1/( AthetaR+1));
                Frmetal(i,j)=(Dh.*G(i,j).*Rtotal)./(4*Coi*Cor);
                Frglass(i,j)=(Dh.*G(i,j).*RtotalG)./(4*Coi*Cor);
            else
                gTheta(i,j)=0;
                gThetaR(i,j)=0;
                G(i,j)=0;
            end
            I(i,j)=mean(dot(L,N));
            %Frmetal(i,j)=Dh.*G(i,j).*Rtotal./4*Coi*Cor;
            %Frglass(i,j)=Dh.*G(i,j).*RtotalG./4*Coi*Cor;
            %Imetal(i,j)=Dh.*G(i,j).*I(i,j)./4*VN(i,j);
            %Iglass(i,j)=Dh.*G(i,j).*I(i,j)./4*VN(i,j);
            Imetal(i,j)=Frmetal(i,j).*I(i,j);
            Iglass(i,j)=Frglass(i,j).*I(i,j);
            
            
        else
            I(i,j)=0;
            Imetal(i,j)=0.5;
            Iglass(i,j)=0.5;
        end
    end
end
