load("data/hub_estimate.mat");
t = thd_prova.time;
theta_d = thd_prova.signals(1).values();
plot(t, theta_d);