% known parameters
Kt = 7.68e-3;
Rs = 0.5;
N  = 14;
Tsample = 1e-3;
Beq = 7.94e-7;
tausf = 9.09e-3;
deg2rad = pi/180;

threshold_min = 10;
threshold_max = 14;

load("data/accel_decel.mat");
t = ia_best.time;
n_sample = length(t);
i = ia_best.signals.values;
vs = Rs * i;
wm = thl_u_best.signals(2).values * N * deg2rad;
wm_prime = al_best_meas.signals.values * N * deg2rad * deg2rad;
wm_prime = abs(wm_prime);
filtered_wm_prime = [];
filtered_wm = [];
filtered_vs = [];
filtered_t = [];
j = 1;
for i=1:n_sample
    if(wm_prime(i)>threshold_min && wm_prime(i) < threshold_max)
        filtered_wm_prime(j) = wm_prime(i);
        filtered_wm(j) = wm(i);
        filtered_vs(j) = vs(i);
        filtered_t(j) = t(i);
        j = j +1;
    end
end

% plotting results
nexttile
plot(filtered_t, filtered_vs)
title('Voltage')
nexttile
plot(filtered_t, filtered_wm_prime)
title('Acceleration (modulus)')
nexttile
plot(filtered_t, filtered_wm)
title('Angular velocity')

Jeq = (1/n_sample)*sum((1/wm_prime)*(vs*Kt/Rs - Beq*wm - (1/N)*tausf*sign(wm)));
sprintf(num2str(Jeq))