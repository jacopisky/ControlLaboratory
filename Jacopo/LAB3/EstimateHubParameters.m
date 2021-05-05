%% Preliminary imports
%% Importing experiment data & estimation
load("data/hub_test.mat");
t = hub_test.time;
theta_b = abs(hub_test.signals(1).values());
prev = theta_b(1);
dir = 0;
count = 1;
peaks = [];
times = [];
visual = zeros(size(t));
for i=2:size(t)
    now = theta_b(i);
    if dir == 1 && now < prev
        % record peak
        peaks(count) = theta_b(i-1);
        times(count) = t(i-1);
        count = count + 1;
        dir = -1;
        visual(i) = theta_b(i);
    elseif dir == -1 && now > prev
        dir = 1;
    elseif dir == 0
        if prev < now
            dir = 1;
        else
            dir = -1;
        end
    end
    prev = now;
end

%% Plot detected peaks
% figure
% plot(t, abs(theta_b))
% hold on
% plot(t, visual)

%% Linear Regression
Y = transpose(log(peaks));
PHI = zeros(count-1,2);
for i=1:count-1
    PHI(i, 1) = i;
    PHI(i, 2) = 1;
end
THETA_HAT = (transpose(PHI)*PHI)\transpose(PHI)*Y;
psi_hat = abs(THETA_HAT(1,1));
delta_hat = psi_hat / sqrt(pi^2 + psi_hat^2);
count2 = 1;
Tk = [];
for i=2:count-1
    Tk(count2) = times(i) - times(i-1);
    count2 = count2 + 1;
end
for i=1:count2-1
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
