%linkage length (in mm) - for tadarida
lh = 39.91; %normalizing to fit dimensions of tadarida?
lr = 66.08;
lw = 131.99;
b = 7.143; %20/2.8
f = 10.714; %30/2.8
%extended state - angle in degrees - converted to radians
theta_se = deg2rad(51);
theta_ee = deg2rad(110);
theta_we = deg2rad(147);
xae = 14.64;
%tucked state - angle in degrees - converted to radians
theta_st = deg2rad(20);
theta_et = deg2rad(41);
theta_wt = deg2rad(35);
xat = 21.147;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%the equations - extended state%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%eq1 = (xae - (lh-x(2))*cos(theta_se))^2 + ((lh-x(2))*sin(theta_se))^2 - (b+x(1))^2 ;
%eq2 = (x(2)*cos(theta_se) + x(5)*cos(theta_ee - theta_se) - b*(xae-(lh-x(2))*cos(theta_se))/(b+x(1)))^2 + (x(2)*sin(theta_se) - x(5)*sin(theta_ee - theta_se) + b*((lh-x(2))*sin(theta_se))/(b+x(1)))^2 - x(3)^2;
%eq3 = ((f/x(3))*(x(2)*cos(theta_se) + x(5)*cos(theta_ee - theta_se) - b*(xae-(lh-x(2))*cos(theta_se))/(b+x(1))) - x(4)*cos(theta_we - theta_ee + theta_se) + (lr + x(5))*cos(theta_ee -theta_se))^2 + ((f/x(3))*(x(2)*sin(theta_se) - x(5)*sin(theta_ee - theta_se) + b*((lh-x(2))*sin(theta_se))/(b+x(1))) - x(4)*sin(theta_we - theta_ee + theta_se) - (lr + x(5))*sin(theta_ee -theta_se))^2 - x(6)^2 ;
% a d e h i j
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%the equations - tucked state%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%eq4 = (xat - (lh-x(2))*cos(theta_st))^2 + ((lh-x(2))*sin(theta_st))^2 - (b+x(1))^2 ;
%eq5 = (x(2)*cos(theta_st) + x(5)*cos(theta_et - theta_st) - b*(xat-(lh-x(2))*cos(theta_st))/(b+x(1)))^2 + (x(2)*sin(theta_st) - x(5)*sin(theta_et - theta_st) + b*((lh-x(2))*sin(theta_st))/(b+x(1)))^2 - x(3)^2;
%eq6 = ((f/x(3))*(x(2)*cos(theta_st) + x(5)*cos(theta_et - theta_st) - b*(xat-(lh-x(2))*cos(theta_st))/(b+x(1))) - x(4)*cos(theta_wt - theta_et + theta_st) + (lr + x(5))*cos(theta_et -theta_st))^2 + ((f/x(3))*(x(2)*sin(theta_st) - x(5)*sin(theta_et - theta_st) + b*((lh-x(2))*sin(theta_st))/(b+x(1))) - x(4)*sin(theta_wt - theta_et + theta_st) - (lr + x(5))*sin(theta_et -theta_st))^2 - x(6)^2 ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F = @(x) [
    (xae - (lh-x(2))*cos(theta_se))^2 + ((lh-x(2))*sin(theta_se))^2 - (b+x(1))^2;
    (x(2)*cos(theta_se) + x(5)*cos(theta_ee - theta_se) - b*(xae-(lh-x(2))*cos(theta_se))/(b+x(1)))^2 + (x(2)*sin(theta_se) - x(5)*sin(theta_ee - theta_se) + b*((lh-x(2))*sin(theta_se))/(b+x(1)))^2 - x(3)^2;
    ((f/x(3))*(x(2)*cos(theta_se) + x(5)*cos(theta_ee - theta_se) - b*(xae-(lh-x(2))*cos(theta_se))/(b+x(1))) - x(4)*cos(theta_we - theta_ee + theta_se) + (lr + x(5))*cos(theta_ee -theta_se))^2 + ((f/x(3))*(x(2)*sin(theta_se) - x(5)*sin(theta_ee - theta_se) + b*((lh-x(2))*sin(theta_se))/(b+x(1))) - x(4)*sin(theta_we - theta_ee + theta_se) - (lr + x(5))*sin(theta_ee -theta_se))^2 - x(6)^2;
    (xat - (lh-x(2))*cos(theta_st))^2 + ((lh-x(2))*sin(theta_st))^2 - (b+x(1))^2;
    (x(2)*cos(theta_st) + x(5)*cos(theta_et - theta_st) - b*(xat-(lh-x(2))*cos(theta_st))/(b+x(1)))^2 + (x(2)*sin(theta_st) - x(5)*sin(theta_et - theta_st) + b*((lh-x(2))*sin(theta_st))/(b+x(1)))^2 - x(3)^2;
    ((f/x(3))*(x(2)*cos(theta_st) + x(5)*cos(theta_et - theta_st) - b*(xat-(lh-x(2))*cos(theta_st))/(b+x(1))) - x(4)*cos(theta_wt - theta_et + theta_st) + (lr + x(5))*cos(theta_et -theta_st))^2 + ((f/x(3))*(x(2)*sin(theta_st) - x(5)*sin(theta_et - theta_st) + b*((lh-x(2))*sin(theta_st))/(b+x(1))) - x(4)*sin(theta_wt - theta_et + theta_st) - (lr + x(5))*sin(theta_et -theta_st))^2 - x(6)^2
    ];
x0 = [10;10;10;10;10;10];
options = optimoptions('fsolve','Display','iter','PlotFcn',@optimplotfirstorderopt);
[x,fval] = fsolve(F,x0,options)
