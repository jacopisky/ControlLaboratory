run("LoadParameters.m");

% requests (used or not depending on design method)
request.Mp    = 0.1;  % max overshooting percentage
request.ts    = 0.15; % settling time
request.tr    = 0.01; % rising time
request.fc    = 5;    % crossover frequency
request.phm   = 90;   % phase margin

% definitions
request.sim.alpha = 4;
pid.sim.T_l = 0.002;

% simple model according to LAB0
Km = drv.dcgain*mot.Kt/(mot.Kt*mot.Ke);
Tm = (mot.R+sens.curr.Rs)*(mot.J+mld.J/(gbox.N^2))/(mot.Kt*mot.Ke);
P_num_reduced = [Km/gbox.N^2];
P_den_reduced = [Tm 1 0];

P = tf(P_num_reduced, P_den_reduced);

% manual set (best for LAB0)
pid.Kp = 60;
pid.Ki = 100;
pid.Kd = 1;

% override
[pid.Kp, pid.Ki, pid.Kd] = getPIDBodeFreq(P, 5, 60, request.sim.alpha);

simp_model.Req = mot.R+sens.curr.Rs;
simp_model.Km = Km;
simp_model.Tm = Tm;
simp_model.P = P;

% clear useless variables
clear P;
clear Km;
clear Tm;
clear P_den_reduced;
clear P_num_reduced;

% set of utility routines
function [Kp,Ki,Kd] = getPIDBodeRequests(plant, Mp, ts, tr, alpha)
    [wgc, phim] = getWPhi(tr, ts, Mp);
    [Kp, Ki, Kd] = getPIDBode(plant, wgc, phim, alpha);
end

function [wgc, phim] = getWPhi(tr, ts, Mp)
    wgc = 2/tr;
    phim = rad2deg(1.04 - 0.8*Mp);
end

function [Kp,Ki,Kd] = getPIDBodeFreq(plant, freq, phim, alpha)
    wgc = 2*pi*freq;
    [Kp,Ki,Kd] = getPIDBode(plant, wgc, phim, alpha);
end

function [Kp,Ki,Kd] = getPIDBode(plant, wgc, phim, alpha)
    [mag, phase] = bode(plant, wgc);
    DeltaK = 1 / mag;
    DeltaPhi = -180 + phim - phase;
    tmp = tand(DeltaPhi);
    Td = (tmp + sqrt(tmp*tmp + 4/alpha))/(2*wgc);
    Ti = alpha * Td;
    Kp = DeltaK * cosd(DeltaPhi);
    Ki = Kp / Ti;
    Kd = Kp * Td;
    s = tf('s');
    C = Kp + Ki * 1/s + Kd * s;
    [a,b,c,d] = margin(C*plant);
    delta_wgc = c - wgc;
    delta_phim = d - phim;
    assert(delta_wgc < 0.01, "failed to build PID");
    assert(delta_phim < 0.01, "failed to build PID");
end