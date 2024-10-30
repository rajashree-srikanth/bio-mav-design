%for tadarida australis
AR = 8.28;
U = 10; %m/s
C = 0.055;
S = 1;
freq = 9.27;
omega = 2*pi*freq;
k = (C*omega)/(2*U);
C1 = (0.5*AR)/(2.32+AR);
C2 = 0.181 + (0.772/AR);
F1 = 1 - (C1*k^2)/(k^2 + C2^2);
G1 = -(C1*C2*k)/(k^2 + C2^2);
syms theta t y
theta_a = deg2rad(9);
gamma = asin(0.6927);
h = -gamma*y*cos(omega*t)
h_dot = diff(h)
alpha = (h_dot*cos(theta-theta_a) + (3/4)*C*)
alpha_dash = 