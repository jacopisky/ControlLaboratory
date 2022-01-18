load("../data/FFAW.mat");
t = thl_meas.time;
err = error.signals.values;
figure
hold on
plot(t, err)
legend("aw + ff")
xlabel("t [s]")
ylabel("error [deg]")
