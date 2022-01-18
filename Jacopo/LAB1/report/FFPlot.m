load("../data/NoFF.mat");
t = thl_meas_noFF.time;
no_ff = thl_meas_noFF.signals.values(:,2);
ref = thl_meas_noFF.signals.values(:,1);
err_no = error_NoFF.signals.values;
load("../data/FF.mat");
ff = thl_meas.signals.values(:,2);
err = error.signals.values;
figure
hold on
plot(t, err_no)
plot(t, err)
legend("no ff", "with ff")
xlabel("t [s]")
ylabel("error [deg]")
