load('integral5000q01disturb10prova2.mat');
figure
plot(t, theta)
hold on
plot(t, gamma)
hold on
plot(t, disturb)
hold on
plot(t, ref)

legend("theta","gamma","disturb", "gamma ref")