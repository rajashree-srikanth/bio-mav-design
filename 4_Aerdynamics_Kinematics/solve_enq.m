%%ERROR IN ALGORITHM - solve() cannot solve non-linear system of equations%% 
%syms x
%A = solve(x^2+2, x)
lh = 110;
lr = 180;
lw = 370;
b = 20;
f = 30;
theta_se = deg2rad(51);
theta_ee = deg2rad(110);
theta_we = deg2rad(147);
xae = 45;
theta_st = deg2rad(20);
theta_et = deg2rad(41);
theta_wt = deg2rad(35);
xat = 65;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms a d e h k m %k in place of i and m in place of j
%%%%%%%%%%%%%%%%%%%%%%%%%%%% extended state %%%%%%%%%%%%%%%%%%%%%%%
eq1 = (xae - (lh-d)*cos(theta_se))^2 + ((lh-d)*sin(theta_se))^2 - (b+a)^2 ;
eq2 = (d*cos(theta_se) + k*cos(theta_ee - theta_se) - b*(xae-(lh-d)*cos(theta_se))/(b+a))^2 + (d*sin(theta_se) - k*sin(theta_ee - theta_se) + b*((lh-d)*sin(theta_se))/(b+a))^2 - e^2;
eq3 = ((f/e)*(d*cos(theta_se) + k*cos(theta_ee - theta_se) - b*(xae-(lh-d)*cos(theta_se))/(b+a)) - h*cos(theta_we - theta_ee + theta_se) + (lr + k)*cos(theta_ee -theta_se))^2 + ((f/e)*(d*sin(theta_se) - k*sin(theta_ee - theta_se) + b*((lh-d)*sin(theta_se))/(b+a)) - h*sin(theta_we - theta_ee + theta_se) - (lr + k)*sin(theta_ee -theta_se))^2 - m^2 ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eq4 = (xat - (lh-d)*cos(theta_st))^2 + ((lh-d)*sin(theta_st))^2 - (b+a)^2 ;
eq5 = (d*cos(theta_st) + k*cos(theta_et - theta_st) - b*(xat-(lh-d)*cos(theta_st))/(b+a))^2 + (d*sin(theta_st) - k*sin(theta_et - theta_st) + b*((lh-d)*sin(theta_st))/(b+a))^2 - e^2;
eq6 = ((f/e)*(d*cos(theta_st) + k*cos(theta_et - theta_st) - b*(xat-(lh-d)*cos(theta_st))/(b+a)) - h*cos(theta_wt - theta_et + theta_st) + (lr + k)*cos(theta_et -theta_st))^2 + ((f/e)*(d*sin(theta_st) - k*sin(theta_et - theta_st) + b*((lh-d)*sin(theta_st))/(b+a)) - h*sin(theta_wt - theta_et + theta_st) - (lr + k)*sin(theta_et -theta_st))^2 - m^2 ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%S = solve(eq1,eq2,eq3,eq4,eq5,eq6,a,d,e,h,k,m)
eqns = [eq1==0,eq2==0,eq3==0,eq4==0,eq5==0,eq6==0];
vars = [a d e h k m];
E = solve(eqns,vars)