% known parameters
Kt = 7.68e-3;
Rs = 0.5;
N  = 14;
Beq = 6.9524e-7;
tausf = 0.0099;
Tsample = 1e-3;
deg2rad = pi/180;

threshold_min = 500;
threshold_max = 10000;

load("data/accel_decel.mat");
t = ia_best.time;
n_sample = length(t);
i = ia_best.signals.values;
vs = Rs * i;
wm = al_best_meas1.signals.values * N * deg2rad;
am = al_best_meas.signals.values * N * deg2rad;
torque = torquem_best.signals.values - (Beq * wm + tausf * sign(wm) / N);
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
Beq = sum / n;
sprintf(num2str(Beq))