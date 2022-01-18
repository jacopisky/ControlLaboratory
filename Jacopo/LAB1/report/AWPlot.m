load("data/aw_sim.mat")
t = aw_sim.time;
ref = aw_sim.signals(1).values;
no_aw = aw_sim.signals(2).values;
aw = aw_sim.signals(3).values;
figure
hold on
plot(t, ref, "--")
plot(t, no_aw)
plot(t, aw)
legend("ref","no aw","aw")
xlabel("t [s]")
ylabel("\vartheta_l [deg]")

load("../data/NoAwNew.mat");

t = thl_meas_noAw.time(1:2500);
no_aw = thl_meas_noAw.signals.values(1:2500);

load("../data/antiwindup360.mat");
aw = thl_meas.signals.values(1:2500);
ref = ones(size(t))*360;

figure
hold on
plot(t, ref, "--")
plot(t, no_aw)
plot(t, aw)
legend("ref","no aw","aw")
xlabel("t [s]")
ylabel("\vartheta_l [deg]")
