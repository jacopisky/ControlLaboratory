%% Importing experiment data
load("tmp/estimate_hub_params.mat");
t = thd_est.time;
theta_d = abs(thd_est.signals.values());
thd = thh_est.signals.values();
%% Time cut
t = t(6500:9000);
theta_d = theta_d(6500:9000);
thd = thd(6500:9000);

%% Signals analyzer
[peaks, idx] = findpeaks(theta_d,'MinPeakProminence',0.7);
visual = zeros(size(t));
visual(idx) = theta_d(idx);
times = t(idx);

dim = size(peaks);
count = dim(1);

%% Linear Regression
Y = log(peaks);
PHI = zeros(count,2);
for i=1:count
    PHI(i, 1) = i;
    PHI(i, 2) = 1;
end
THETA_HAT = (transpose(PHI)*PHI)\transpose(PHI)*Y;
psi_hat = abs(THETA_HAT(1,1));
delta_hat = psi_hat / sqrt(pi^2 + psi_hat^2);
count2 = 1;
Tk = [];
for i=2:count
    Tk(count2) = times(i) - times(i-1);
    count2 = count2 + 1;
end
for i=1:count-1
    Tk(i) = pi/Tk(i);
end
w_hat = 1/count * sum(Tk);
wn_hat = w_hat / sqrt(1-delta_hat^2);
est_par.k  = abs(mld.Jb*wn_hat^2);
est_par.Bb = 2*sqrt(mld.Jb*est_par.k)*delta_hat;

sigma_hat = -delta_hat*wn_hat;
theta_d0  = theta_d(1);
A_hat     = sqrt(theta_d0^2+sigma_hat^2*theta_d0^2/(w_hat^2));
phi_hat   = atan(sigma_hat/w_hat);
tt = linspace(t(1), t(end), 1000);
yy = A_hat*exp(sigma_hat*(tt-times(1)));


yyy = -A_hat*exp(sigma_hat*(tt-times(1))).*cos(w_hat*(tt-times(1))+phi_hat);

%% Plot detected peaks
figure
plot(t, theta_d)
hold on
plot(t, visual)
plot(tt, yy, "--")
legend("|\vartheta_d|", "peak detect","Ae^{\sigma t}")
xlim([6.5 9])
ylim([0 20])
xlabel("t[s]")
ylabel("|\vartheta_d| [deg]")

figure
plot(t, thd)
hold on
plot(tt, yyy, "--")
legend("\vartheta_d", "reconstructed \vartheta_d")
xlim([6.5 9])
ylim([-20 20])
xlabel("t[s]")
ylabel("\vartheta_d [deg]")

%% Workspace cleanup
clear i;
clear t;
clear theta_b;
clear prev;
clear now;
clear dir;
clear visual;
clear Y;
clear PHI;
clear THETA_HAT;
clear count;
clear peaks;
clear count2;
clear Tk;
clear times;
clear delta_hat;
clear w_hat;
clear wn_hat;
clear psi_hat;
clear t_start;
clear t_stop;
clear offset;
clear t_off;
clear new_t;
clear dim;
clear idx;