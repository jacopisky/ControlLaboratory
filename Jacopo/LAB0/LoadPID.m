run("LoadParameters.m");

% requests
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

% manual set
pid.Kp = 60;
pid.Ki = 100;
pid.Kd = 1;

% override
[pid.Kp, pid.Ki, pid.Kd] = getPIDBode(P, 100, 90, request.sim.alpha);

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
end

function [Kp,Ki,Kd] = getPIDProgram(P, Mp, ts, tr, alpha)
    best_Mp = Inf;
    best_ts = Inf;
    best_tr = Inf;
    bestKp  = 0;
    bestKi  = 0;
    bestKd  = 0;
    for kp=0:1:50
        for ki=0:1:50
            kd = alpha*ki;
            s = tf('s');
            C = kp + ki/s + kd*s;
            W = feedback(C*P,1);
            curr_Mp = stepinfo(W).Overshoot;
            curr_ts = stepinfo(W).SettlingTime;
            curr_tr = stepinfo(W).RiseTime;
            if(curr_Mp <= Mp && curr_ts <= ts && curr_tr <= tr)
                if(curr_Mp <= best_Mp && curr_ts <= best_ts && curr_tr <= best_tr)
                    bestKp = kp;
                    bestKi = ki;
                    bestKd = kd;
                end
            end
        end
    end
    Kp = bestKp;
    Ki = bestKi;
    Kd = bestKd;
end