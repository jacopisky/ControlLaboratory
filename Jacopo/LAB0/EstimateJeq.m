run("EstimateBeq_tausf.m");

threshold_min = 500;
threshold_max = 10000;

load("data/accel_decel.mat");
t = ia_best.time;
n_sample = length(t);
i = ia_best.signals.values;
vs = Rs * i;
wm = al_best_meas1.signals.values * N * deg2rad;
am = al_best_meas.signals.values * N * deg2rad;
torque = torquem_best.signals.values - (Beq_est * wm + tausf_est * sign(wm) / N);
% plotting results
figure
nexttile
plot(t, vs)
title('Voltage')
nexttile
plot(t, wm)
title('Angular velocity')
nexttile
plot(t, am)
title('Acceleration')

am = abs(am);
torque = abs(torque);
sum = 0;
n = 0;
for i = 1:n_sample
    if am(i) ~= 0
        sum = sum + torque(i)/am(i);
        n = n + 1;
    end
end
Jeq_est = sum / n;
sprintf(num2str(Jeq_est))