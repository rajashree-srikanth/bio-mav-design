clc
clear all

%%%%%%%%%%%%%%%%%%%%%%INPUTS%%%%%%%%%%%%%%%%%%%%%%
WbyS = 20; %in N/m2 13
S = 0.10859754; %0.02584;%in m2
U = 8.5; % forward velocity (m/s)
st = 0.26; %(f*2*sin(gamma)*(wingspan/2))/U % 0.294689103
theta_a=deg2rad(5.5); % pitching angle of flapping angle with respect to U(rad)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta=-10:0.82021:10; % deg/meter (1 to 6.096 in intervals of 0.25) deg per feet 56.21621622
    
    g=9.830590516507780; %9.80665; Hassanalian paper  lat and height 
    eta=0.98; % leading edge suction efficiency
    rho= 1.1033; % density of air(kg/m3) at 500m
    mew=1.8458*10^-5; % (coefficient of viscocity(Ns/m2) or in (kg/ms)) at 500m
    W = WbyS*S; % in newton
    Wg = W/g % in kg
    
    aspectratio=8.28;
    wingspan = sqrt(aspectratio*S)
    MAC = wingspan/aspectratio;
    dy=wingspan/20;%section width
    ddy=dy*10;
    
%     gamma = deg2rad(60/2); %gamma by 2 not gamma
    f=1.34*(Wg^(3/8))*(g^(1/2))*(wingspan^(-23/24))*(S^(-1/3))*(rho^(-3/8)) %7.11; % flapping frequency (Hz) (0.3*U)/(wingspan*sin(gamma)) 
    gamma = asin((st*U)/(f*wingspan)) %gamma by 2 not gamma
    disp(rad2deg(gamma*2))
    
    alphastall=deg2rad(13);
    Cmac=0.025;
    alpha_o=deg2rad(0.5);

    span = {};
    c = {};
    for x=0:0.0841:0.37653
        y_1 = 3.3571*x^2 - 2.0751*x + 1.8217;
        c = [c, y_1*MAC];
        span = [span, x*(wingspan/2)];
    end
    for x=0.37653:0.0841:0.75131
        y_2 = 2.875*x^3 - 3.8143*x^2 + 0.5437*x + 1.6983;
        c = [c, y_2*MAC];
        span = [span, x*(wingspan/2)];
    end
    for x=0.75131:0.0841:1
        y_3 = -37.92*x^3 + 104.44*x^2 - 99.838*x + 33.316;
        c = [c, y_3*MAC];
        span = [span, x*(wingspan/2)];
    end
    %c = [c, 0*MAC];
    %span = [span, (wingspan/2)];
    
    
    % span = [0 0.0195423 0.0382702 0.0578125 0.0773548 0.0960827 0.115625 0.135167 0.153895 0.173438 0.19298 0.211708 0.23125];
    % c = [0.1018852 0.0937705 0.087459 0.0847541 0.0838525 0.0829508 0.0766394 0.0721311 0.0685245 0.0658196 0.0405738 0.02073766 0];

    dt=1/(20*f); %time interval

    for z=1:length(beta)
        theta_o=deg2rad(beta(z)); %=1 pitching motion amplitude(rad/feet)
        xcount = 0;
        
        for x = 1:1:12
            o=((pi*f*c{x})/U); % reduced frequency
            tcount = 0;
            count=1;
            xcount = xcount+1;
            for t = 0:dt:1/f
                tcount = tcount+1;
                dihedral=gamma*cos(2*pi*f*t);
                theta=-(theta_o*span{x}*sin(2*pi*f*t))+theta_a;
                h=-gamma*span{x}*cos(2*pi*f*t);
                hdot=gamma*span{x}*sin(2*pi*f*t)*2*pi*f;
                hdoubledot=-((2*pi*f)^2)*h;
                thetadot=-theta_o*span{x}*cos(2*pi*f*t)*2*pi*f;
                thetadoubledot=theta_o*span{x}*sin(2*pi*f*t)*2*pi*f*2*pi*f;
                Vpitch=0.75*c{x}*thetadot;
                Vplunge=hdot*cos(theta-theta_a);
                Downwardvelocity=U*(theta-theta_a);
                Vpitchdot=0.75*c{x}*thetadoubledot;
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
                alphap=((AR/(2+AR))*(F*alpha+(c{x}/(2*U))*(G/o)*alphadot))-w/U;%%%%%%%%%
                alpharr(xcount,tcount) = alphap;
                b=(U*alphadot)-(0.25*c{x}*thetadoubledot);
                Vx=(U*cos(theta))-(hdot*sin(theta-theta_a));
                % velocityx(x,count)=U*cos(theta)-hdot*sin(theta-theta_a);
                tempp(x,count)=(((alphap+theta_a-(0.75*c{x}*thetadot/U))))*180/pi;
                if (alphap+theta_a-((0.75*c{x}*thetadot)/U)) >= alphastall
                    Vn=hdot*cos(theta-theta_a)+0.5*c{x}*thetadot+U*sin(theta);
                    Vs=sqrt((Vx^2+Vn^2));
                    dNc(x,count)=(1.98*Vs*Vn*c{x}*rho*dy)/2;
                    dNa(x,count)=(rho*pi*c{x}^2*dy*b)/8;
                    dDf(x,count)=0;
                    dTs(x,count)=0;
                    dDcamb(x,count)=0;
                    dPin(x,count)=(dNc(x,count)+dNa(x,count))*(hdot*cos(theta-theta_a)+0.5*c{x}*thetadot);
                    dL(x,count)=((dNa(x,count)+dNc(x,count))*cos(theta))+((dTs(x,count)-dDf(x,count)-dDcamb(x,count))*sin(theta))*cos(dihedral);
                    dT(x,count)= ((dTs(x,count)-dDf(x,count)-dDcamb(x,count))*cos(theta))-((dNa(x,count)+dNc(x,count))*sin(theta));
                else
                    Vy=(U*(alphap+theta_a))-(0.5*c{x}*thetadot);
                    % velocityy(x,count)=U*(alphap+theta_a)-0.5*c(x)*thetadot;
                    V=sqrt((Vx)^2+(Vy)^2);
                    dDcamb(x,count)=-2*pi*alpha_o*(alphap+theta_a)*rho*U*V*c{x}*dy*0.5;
                    dNc(x,count)=rho*U*V*pi*(alphap+theta_a+alpha_o)*dy*c{x};
                    %circulatory normal force
                    dNa(x,count)=(rho*pi*(c{x}^2)*b*dy)/4; %normal force due to apparent mass effect
                    dTs(x,count)=eta*2*pi*((alphap+theta_a-((0.25*c{x}*thetadot)/U))^2)*rho*U*(V/2).*dy*c{x}; %leading edge suction force
                    re=(rho*U*c{x})/(mew);
                    Cdf(x,count)=0.89/((log10(re))^2.58);
                    dDf(x,count)=Cdf(x,count)*rho*((Vx)^2)*c{x}*dy*0.5; %friction drag
                    dFx=(dTs(x,count)-dDcamb(x,count)-dDf(x,count));
                    dN=((dNa(x,count)+dNc(x,count)));
                    dL(x,count)=(dN*cos(theta)+dFx*sin(theta))*cos(dihedral);
                    dT(x,count)= (dFx*cos(theta))-(dN*sin(theta));
                    dM=-(((rho*pi*(c{x}^3)*thetadot*U)/16)+(rho*pi*(c{x}^4)*thetadoubledot)/128)*dy;
                    dmac=Cmac*rho*dy*(c{x}.^2)*U*V*0.5;
                    dPin(x,count)=((dFx)*hdot*sin(theta-theta_a))+((dN)*(hdot*cos(theta-theta_a)+(0.25*c{x}*thetadot)))+(0.25*c{x}*thetadot*dNa(x,count))-((dM+dmac)*thetadot);
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
        if beta(z)==-10
            timebeta10=time;
            liftbeta10=lift;
            plot(timebeta10,liftbeta10)
        end
        avgliftt(z)=avglift*f;
        avgthrustt(z)=avgthrust*f;
        averagepowerinput(z)=avgpower*f;
        efficiency(z)=(avgthrustt(z)*U)/averagepowerinput(z);
    end

fig1 = figure
left_color = [0.6350 0.0780 0.1840];
right_color = [0 0.4470 0.7410];
set(fig1,'defaultAxesColorOrder',[left_color; right_color]);
yyaxis left
plot(beta,2*averagepowerinput,'linewidth',0.75)%,'color', '#A2142F', 'linewidth',1)
ylabel('Average input power (in W)')
yyaxis right
plot(beta,efficiency)%,'color', '#77AC30', 'linewidth',1)
ylabel('Propulsive efficiency','linewidth',0.75)
xlabel('\beta in deg/metre')

fig2 = figure
left_color = [0.8500 0.3250 0.0980];
right_color = [0.4940 0.1840 0.5560];
set(fig2,'defaultAxesColorOrder',[left_color; right_color]);
yyaxis left
plot(beta,2*avgthrustt,'linewidth',0.75)
ylabel('Average thrust (in N)')
yyaxis right
plot(beta,2*avgliftt,'linewidth',0.75)
ylabel('Average lift (in N)')
xlabel('\beta in (deg/metre)')