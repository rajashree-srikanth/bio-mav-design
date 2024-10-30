clc
clear all
format short
alphastall=deg2rad(13);
Cmac=0.025;
aspectratio=7.79;
alpha_o=deg2rad(0.5);
gamma=deg2rad(135.98/2); % maximum flapping
eta=0.98; % leading edge suction efficiency
rho= 1.1033; % density of air(kg/m3) at 500m
theta_a=deg2rad(8); % pitching angle of flapping angle with respect to U(rad)
f=7.01; % flapping frequency (Hz)
U=10; % forward velocity (m/s)
mew=1.8458*10^-5; % (coefficient of viscocity(Ns/m2) or in (kg/ms)) at 500m
dt=1/(20*f); %time interval
dy=0.02308;%section width
ddy=dy*10;
span = [0 0.0186915 0.0381958 0.0577 0.0772042 0.0967085 0.1154 0.134904 0.153596 0.1731 0.192604 0.212108 0.2308]
c = [0.0740984 0.0681967 0.0636066 0.0616393 0.0609836 0.0596721 0.0557377 0.052459 0.0498361 0.0478689 0.02950815 0.01508195 0]
beta=0:0.82021:40; % deg/meter (1 to 6.096 in intervals of 0.25) deg per feet
for z=1:length(beta)
    theta_o=deg2rad(beta(z)); % pitching motion amplitude(rad/feet)
    for x = 1:1:12
        o=((pi*f*c(x))/U); % reduced frequency
        count=1;
        for t = 0:dt:1/f
            dihedral=gamma*cos(2*pi*f*t);
            theta=-(theta_o*span(x)*sin(2*pi*f*t))+theta_a;
            h=-gamma*span(x)*cos(2*pi*f*t);
            hdot=gamma*span(x)*sin(2*pi*f*t)*2*pi*f;
            hdoubledot=-((2*pi*f)^2)*h;
            thetadot=-theta_o*span(x)*cos(2*pi*f*t)*2*pi*f;
            thetadoubledot=theta_o*span(x)*sin(2*pi*f*t)*2*pi*f*2*pi*f;
            Vpitch=0.75*c(x)*thetadot;
            Vplunge=hdot*cos(theta-theta_a);
            Downwardvelocity=U*(theta-theta_a);
            Vpitchdot=0.75*c(x)*thetadoubledot;
            Vplungedot=-hdot*sin(theta-theta_a)*thetadot+hdoubledot*cos(theta-theta_a);
            Downwardvelocitydot=U*thetadot;
            temp(x,count)=(Vpitch+Vplunge+Downwardvelocity)/U;
            alpha=(temp(x,count));
            angle(x,count)= (((Vpitchdot+Vplungedot+Downwardvelocitydot)/U))*180/pi;
            alphadot=(Vpitchdot+Vplungedot+Downwardvelocitydot)/U;
            AR=aspectratio;
            C1=(0.5*AR)/(2.32+AR);
            C2=(0.181)+(0.772/AR);
            F=1-((C1*o^2)/(o^2+C2^2));
            G=-((C1*C2*o)/(o^2+C2^2));
            C=(F^2+G^2)^0.5;
            w=(U*2*(alpha_o+theta_a))/(2+AR);
            alphap=((AR/(2+AR))*(F*alpha+(c(x)/(2*U))*(G/o)*alphadot))-w/U;
            b=(U*alphadot)-(0.25*c(x)*thetadoubledot);
            Vx=(U*cos(theta))-(hdot*sin(theta-theta_a));
            % velocityx(x,count)=U*cos(theta)-hdot*sin(theta-theta_a);
            tempp(x,count)=(((alphap+theta_a-(0.75*c(x)*thetadot/U))))*180/pi;
            if (alphap+theta_a-((0.75*c(x)*thetadot)/U)) >= alphastall
                Vn=hdot*cos(theta-theta_a)+0.5*c(x)*thetadot+U*sin(theta);
                Vs=sqrt((Vx^2+Vn^2));
                dNc(x,count)=(1.98*Vs*Vn*c(x)*rho*dy)/2;
                dNa(x,count)=(rho*pi*c(x)^2*dy*b)/8;
                dDf(x,count)=0;
                dTs(x,count)=0;
                dDcamb(x,count)=0;
                dPin(x,count)=(dNc(x,count)+dNa(x,count))*(hdot*cos(theta-theta_a)+0.5*c(x)*thetadot);
                dL(x,count)=((dNa(x,count)+dNc(x,count))*cos(theta))+((dTs(x,count)-dDf(x,count)-dDcamb(x,count))*sin(theta))*cos(dihedral);
                dT(x,count)= ((dTs(x,count)-dDf(x,count)-dDcamb(x,count))*cos(theta))-((dNa(x,count)+dNc(x,count))*sin(theta));
            else
                Vy=(U*(alphap+theta_a))-(0.5*c(x)*thetadot);
                % velocityy(x,count)=U*(alphap+theta_a)-0.5*c(x)*thetadot;
                V=sqrt((Vx)^2+(Vy)^2);
                dDcamb(x,count)=-2*pi*alpha_o*(alphap+theta_a)*rho*U*V*c(x)*dy*0.5;
                dNc(x,count)=rho*U*V*pi*(alphap+theta_a+alpha_o)*dy*c(x);
                %circulatory normal force
                dNa(x,count)=(rho*pi*(c(x)^2)*b*dy)/4; %normal force due to apparent mass effect
                dTs(x,count)=eta*2*pi*((alphap+theta_a-((0.25*c(x)*thetadot)/U))^2)*rho*U*(V/2).*dy*c(x); %leading edge suction force
                re=(rho*U*c(x))/(mew);
                Cdf(x,count)=0.89/((log10(re))^2.58);
                dDf(x,count)=Cdf(x,count)*rho*((Vx)^2)*c(x)*dy*0.5; %friction drag
                dFx=(dTs(x,count)-dDcamb(x,count)-dDf(x,count));
                dN=((dNa(x,count)+dNc(x,count)));
                dL(x,count)=(dN*cos(theta)+dFx*sin(theta))*cos(dihedral);
                dT(x,count)= (dFx*cos(theta))-(dN*sin(theta));
                dM=-(((rho*pi*(c(x)^3)*thetadot*U)/16)+(rho*pi*(c(x)^4)*thetadoubledot)/128)*dy;
                dmac=Cmac*rho*dy*(c(x).^2)*U*V*0.5;
                dPin(x,count)=((dFx)*hdot*sin(theta-theta_a))+((dN)*(hdot*cos(theta-theta_a)+(0.25*c(x)*thetadot)))+(0.25*c(x)*thetadot*dNa(x,count))-((dM+dmac)*thetadot);
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
        for k=1:x
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

figure(1)
yyaxis left
plot(beta,2*averagepowerinput)
ylabel('Average input power in W')
yyaxis right
plot(beta,efficiency)
ylabel('Propulsive efficincy')
xlabel('beta in deg/meter')
% [LX,A1,A2]= plotyy(beta,averagepowerinput,beta,efficiency,'plot')
% set(get(LX(1),'Ylabel'),'String','Avg Input Power(watt)')
% set(get(LX(2),'Ylabel'),'String','Propulsive Efficiency')

figure(2)
% tt = convforce(avgthrustt,'N','lbf');
% ll = convforce(avgliftt,'N','lbf');
yyaxis left
plot(beta,2*avgthrustt)
ylabel('Thrust in N')
yyaxis right
plot(beta,2*avgliftt)
ylabel('Lift in N')
xlabel('beta in deg/meter')
% [AX,H1,H2] = plotyy(beta,2*ll,beta,2*tt,'plot')

% set(get(AX(1),'Ylabel'),'String','Avg Lift(lbf)')
% set(get(AX(2),'Ylabel'),'String','Avg Thrust(lbf)')
% yyaxis left
% plot(beta,2*avgliftt)
% yyaxis right
% plot(beta,2*avgthrustt)