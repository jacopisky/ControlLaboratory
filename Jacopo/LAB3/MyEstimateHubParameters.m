%% Importing experiment data & estimation
load("data/hub_estimate.mat");
t = thd_prova.time;
theta_b = abs(thd_prova.signals(1).values());

% filtering parameters
t_start = 3.9;
t_stop = 6.9;
offset = 0.13;
t_off = 0.04; % oscillation period (approx)

%% Peak Finder
prev = theta_b(1);
prev_t = 0;
dir = 0;
count = 1;
peaks = [];
times = [];
visual = zeros(size(t));
for i=2:size(t)
    now = theta_b(i);
    new_t = t(i);
    if t(i) > t_start && t(i) < t_stop  && new_t > prev_t + t_off
        if dir == 1 && now < prev-offset
            % record peak
            peaks(count) = theta_b(i-1);
            times(count) = t(i-1);
            count = count + 1;
            dir = -1;
            visual(i) = theta_b(i);
            prev_t = new_t;
        elseif dir == -1 && now > prev+offset
            dir = 1;
            prev_t = new_t;
        elseif dir == 0
            if prev < now
                dir = 1;
            else
                dir = -1;
            end
        end
    end
    if dir == 0
            if prev < now
                dir = 1;
            else
                dir = -1;
            end
    end
    prev = now;
end

%% Plot detected peaks
figure
plot(t, abs(theta_b))
hold on
plot(t, visual)

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
clear t_start;
clear t_stop;
clear offset;
clear t_off;
clear new_t;