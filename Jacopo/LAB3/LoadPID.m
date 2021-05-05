run("LoadParameters.m");

% requests (used or not depending on design method)
request.Mp    = 0.3;  % max overshooting percentage
request.ts    = 0.85; % settling time
request.tr    = 0.1; % rising time
request.fc    = 5;    % crossover frequency
request.phm   = 60;   % phase margin

% definitions
request.sim.alpha = 4;

% simple model according to LAB0
s = tf('s');
D_tau_prime = mld.Jeq*mld.Jb*s^3+(mld.Jeq*mld.Bb+mld.Jb*mld.Beq)*s^2+(mld.Beq*mld.Bb+mld.k*(mld.Jeq+mld.Jb/gbox.N/gbox.N))*s+mld.k*(mld.Beq+mld.Bb/gbox.N/gbox.N);
P = 1/gbox.N/s*drv.dcgain*mot.Kt*(mld.Jb*s^2+mld.Bb*s+mld.k)/((sens.curr.Rs+mot.R)*D_tau_prime+mot.Kt*mot.Ke*(mld.Jb*s^2+mld.Bb*s+mld.k));

[request.wgc, request.phm] = getWgcPhim(request.Mp, request.ts);
request.fc = request.wgc/2/pi;

[pid.Kp, pid.Ki, pid.Kd] = getPIDBodeFreq(P, request.fc, request.phm, request.sim.alpha);
pid.sim.T_l = 1/(10*(2*pi*request.fc));
pid.Kw      = 5/request.ts;

% clear useless variables
clear P;
clear s;
clear D_tau_prime;

% set of utility routines
function [Kp,Ki,Kd] = getPIDBodeRequests(plant, Mp, ts, tr, alpha)
    [wgc, phim] = getWPhi(tr, ts, Mp);
    [Kp, Ki, Kd] = getPIDBode(plant, wgc, phim, alpha);
end

function [wgc, phim] = getWgcPhim(Mp, ts)
    delta = log(1/Mp)/sqrt(pi^2+log(1/Mp)^2);
    wgc = 3/delta/ts;
    phim = 180/pi*atan(2*delta/sqrt(sqrt(1+4*delta^4)-2*delta^2)); % Phase margin [deg]
end

function [wgc, phim] = getWPhi(tr, Mp)
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
    Td = (tmp + sqrt(tmp^2 + 4/alpha))/(2*wgc);
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

function [Kp,Ki,Kd] = getPIDProgram(P, Mp, ts, tr)
    best_Mp = Inf;
    best_ts = Inf;
    best_tr = Inf;
    bestKp  = 0;
    bestKi  = 0;
    bestKd  = 0;
    for kp=0:1:50
        for ki=0:1:10
            for kd=0:1:10
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
    end
    Kp = bestKp;
    Ki = bestKi;
    Kd = bestKd;
end