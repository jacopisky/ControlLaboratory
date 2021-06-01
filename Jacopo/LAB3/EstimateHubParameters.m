%% Importing experiment data & estimation
load("data/hub_estimate.mat");
t = thd_prova.time;
th = abs(thd_prova.signals(1).values());

[peaks, idx] = findpeaks(th,'MinPeakProminence',0.7);
visual = zeros(size(t));
visual(idx) = th(idx);
times = t(idx);

dim = size(peaks);
count = dim(1);
%% Plot detected peaks
figure
plot(t, abs(th))
hold on
plot(t, visual)

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