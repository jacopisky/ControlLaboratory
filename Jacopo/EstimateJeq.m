% known parameters
Kt = 7.68e-3;
Rs = 0.5;
N  = 14;
Tsample = 1e-3;
Beq = 7.94e-7;
tausf = 9.09e-3;
deg2rad = pi/180;

load("data/accel_decel.mat");
t = ia_best.time;
n_sample = length(t);
i = ia_best.signals.values;
vs = Rs * i;
wm = thl_u_best.signals(2).values * N * deg2rad;
wm_prime = al_best_meas.signals.values * deg2rad*deg2rad;

% plotting results
nexttile
plot(t, vs)
title('Voltage')
nexttile
plot(t, wm_prime)
title('Acceleration')
nexttile
plot(t, wm)
title('Angular velocity')

Jeq = (1/n_sample)*sum((1/wm_prime)*(Kt*Rs*vs - Beq*wm - (tausf/N)*sign(wm)));
sprintf(num2str(Jeq))