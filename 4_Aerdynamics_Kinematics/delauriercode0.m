clc
%clear all
%format short
alphastall=13*pi/180;
Cmac=0.025;
aspectratio=4.7;
alpha_o=0.5*pi/180;
gamma=33.5*pi/180; % maximum flapping angle(rad)
eta=0.98; % leading edge suction efficiency
rho=1.23; % density of air(kg/m^3)
theta_a=9*pi/180; % pitching angle of flapping angle with respect to U(rad)
f=6.15; % flapping frequency (Hz)
U=2.8; % forward velocity (m/s)
mew=1.46e-5; % (coefficient of viscocity(kg/ms))
dt=1/(20*f); %time interval
dy=0.01;%section width
ddy=dy*100;
beta=0:0.25:2.5;
avgliftt=1;
avgthrustt=1;
averagepowerinput=1;
efficiency=1;
for z=1:length(beta)
    theta_o=beta(z)*1*pi/180; % pitching motion amplitude(rad/m)
    % beta(z+1)=zz(z)*1;
    for x = 0:ddy:20
        num=(x/ddy)+1;
        c(num) = (- 0.036*x^2 + 0.13*x + 12)*0.01;
        o=(pi*f*c(num)/U); % reduced frequency
        count=1;
        for t = 0:dt:1/f
            dihedral=gamma*cos(2*pi*f*t);
            theta=-theta_o*(x+(9/2))*sin(2*pi*f*t)+theta_a;
            h=-gamma*(x+(9/2))*cos(2*pi*f*t);
            hdot=gamma*(x+(9/2))*sin(2*pi*f*t)*2*pi*f;
            hdoubledot=-((2*pi*f)^2)*h;
            thetadot=-theta_o*(x+(9/2))*cos(2*pi*f*t)*2*pi*f;
            thetadoubledot=theta_o*(x+(9/2))*sin(2*pi*f*t)*2*pi*f*2*pi*f;
            Vpitch=0.75*c(num)*thetadot;
            Vplunge=hdot*cos(theta-theta_a);
            Downwardvelocity=U*(theta-theta_a);
            Vpitchdot=0.75*c(num)*thetadoubledot;
            Vplungedot=-hdot*sin(theta-theta_a)*thetadot+hdoubledot*cos(theta-theta_a);
            Downwardvelocitydot=U*thetadot;
            temp(num,count)=(Vpitch+Vplunge+Downwardvelocity)/U;
            alpha=((temp(num,count)));
            angle(num,count)=((Vpitchdot+Vplungedot+Downwardvelocitydot)/U)*180/pi;
            alphadot=(Vpitchdot+Vplungedot+Downwardvelocitydot)/U;
            AR=aspectratio;
            C1=0.5*AR/(2.32+AR);
            C2=0.181+(0.772/AR);
            F=1-(C1*o*o/(o*o+C2^2));
            G=-(C1*C2*o)/(o*o+C2^2);
            C=(F*F+G*G)^0.5;
            w=U*2*(alpha_o+theta_a)/(2+AR);
            alphap=((AR/(2+AR))*(F*alpha+(c(num)/(2*U))*(G/o)*alphadot))-w/U;
            b=U*alphadot-0.25*c(num)*thetadoubledot;
            Vx=U*cos(theta)-hdot*sin(theta-theta_a);
            % velocityx(num,count)=U*cos(theta)-hdot*sin(theta-theta_a);
            tempp(num,count)=((alphap+theta_a-(0.75*c(num)*thetadot/U)))*180/pi;
            if (alphap+theta_a-((0.75*c(num)*thetadot)/U)) >= alphastall
                Vn=hdot*cos(theta-theta_a)+0.5*c(num)*thetadot+U*sin(theta);
                Vs=sqrt((Vx*Vx+Vn*Vn));
                dNc(num,count)=1.98*Vs*Vn*c(num)*rho*dy/2;
                dNa(num,count)=rho*pi*c(num)*c(num)*(b/8)*dy;
                dDf(num,count)=0;
                dTs(num,count)=0;
                dDcamb(num,count)=0;
                dPin(num,count)=(dNc(num,count)+dNa(num,count))*(hdot*cos(theta-theta_a)+0.5*c(num)*thetadot);
                dL(num,count)=((dNa(num,count)+dNc(num,count))*cos(theta)+(dTs(num,count)-dDf(num,count)-dDcamb(num,count))*sin(theta))*cos(dihedral);
                dT(num,count)= (dTs(num,count)-dDf(num,count)-dDcamb(num,count))*cos(theta)-(dNa(num,count)+dNc(num,count))*sin(theta);
            else
                Vy=U*(alphap+theta_a)-0.5*c(num)*thetadot;
                % velocityy(num,count)=U*(alphap+theta_a)-0.5*c(num)*thetadot;
                V=sqrt((Vx)^2+(Vy)^2);
                dDcamb(num,count)=-2*pi*alpha_o*(alphap+theta_a)*rho*U*V*c(num)*dy*0.5;
                dNc(num,count)=rho*U*V*pi*(alphap+theta_a+alpha_o )*dy*c(num); %circulatory normal force
                dNa(num,count)=rho*pi*c(num)*c(num)*b/4*dy; %normal force due to apparent mass effect
                dTs(num,count)=eta*2*pi*(alphap+theta_a-((0.25*c(num)*thetadot)/U))^2*rho*U*(V/2).*dy*c(num); %leading edge suction force
                re=U*c(num)/(mew);
                Cdp(num,count)=0.89*(log10(re))^-2.58;
                dDf(num,count)=Cdp(num,count)*rho*(Vx)^2*c(num)*dy*0.5; % friction drag
                dFx=(dTs(num,count)-dDcamb(num,count)-dDf(num,count));
                dN=((dNa(num,count)+dNc(num,count)));
                dL(num,count)=(dN*cos(theta)+dFx*sin(theta))*cos(dihedral);
                dT(num,count)= dFx*cos(theta)-dN*sin(theta);
                dM=-((rho*pi*c(num).*c(num).*c(num)*thetadot*U/16)+(rho*pi*c(num).*c(num).*c(num)*c(num)*thetadoubledot)/128)*dy;
                dmac=Cmac*rho*dy*(c(num).^2)*U*V*0.5;
                dPin(num,count)=((dTs(num,count)-dDf(num,count)-dDcamb(num,count))*hdot*sin(theta- theta_a))+((dNa(num,count)+dNc(num,count))*(hdot*cos(theta-theta_a)+0.25*c(num)*thetadot))+(0.25*c(num)*thetadot*dNa(num,count))-(dM+dmac)*thetadot;
            end
            time(count)=t;
            count=count+1;
        end
    end
    avglift=0;
    avgthrust=0;
    avgpower=0;
    for i=1:count-1
        lift(i)=0;
        thrust(i)=0;
        power(i)=0;
    end
    for i=1:count-1
        for k=1:num
            lift(i)=lift(i)+dL(k,i);
            thrust(i)=thrust(i)+dT(k,i);
            power(i)=power(i)+dPin(k,i);
        end
        avglift=avglift+lift(i)*dt;
        avgthrust=avgthrust+thrust(i)*dt;
        avgpower=avgpower+power(i)*dt;
    end
    avgliftt(z)=avglift*f;
    avgthrustt(z)=avgthrust*f;
    averagepowerinput(z)=avgpower*f;
    efficiency(z)=(avgthrustt(z)*U)/averagepowerinput(z);
end
% plot(time,lift,time,thrust)
% plot(beta,avgthrustt)
[LX,A1,A2]= plotyy(beta,2*averagepowerinput,beta,efficiency,'plot')
set(get(LX(1),'Ylabel'),'String','Avg Input Power(watt)')
set(get(LX(2),'Ylabel'),'String','Propulsive Efficiency')
figure (2)
[AX,H1,H2] = plotyy(beta,2*avgliftt,beta,2*avgthrustt,'plot')
set(get(AX(1),'Ylabel'),'String','Avg Lift(N)')
set(get(AX(2),'Ylabel'),'String','Avg Thrust(N)')
