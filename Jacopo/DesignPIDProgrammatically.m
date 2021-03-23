% defining reduced plant (no static friction) (dominant poles)
Beq_simpl = 0;
Req = Ra + Rs;
Km  = Kdrv*Kt/(Req*Beq_simpl+Kt*Ke);
Tm  = Req*Jeq/(Req*Beq_simpl+Kt*Ke);
P_num_reduced = [Km/N];
P_den_reduced = [Tm 1 0];

P_reduced = tf(P_num_reduced, P_den_reduced);
P = P_reduced;

% requests
Mp    = 0.1;
ts    = 0.15;
tr    = 0.01;

% definitions
alpha = 4;

[Kp, Ki, Kd] = getPIDProgram(P,Mp,ts,tr, alpha);
s = tf('s');
C = Kp + Ki/s + Kd*s;
W = feedback(C*P,1);
step(W)

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