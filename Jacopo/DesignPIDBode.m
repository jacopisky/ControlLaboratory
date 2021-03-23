% defining reduced plant (no static friction) (dominant poles)
Beq_simpl = 0;
Req = Ra + Rs;
Km  = Kdrv*Kt/(Req*Beq_simpl+Kt*Ke);
Tm  = Req*Jeq/(Req*Beq_simpl+Kt*Ke);
P_num_reduced = [Km/N];
P_den_reduced = [Tm 1 0];

% defining the plant (no static friction)
Kmec  = 1/Beq;
Tmec  = Jeq/Beq;
P_num = [Ke*Kmec*Keq*Kdrv/N];
P_den = [Tdrv*Teq*Tmec Tdrv*Teq+Tdrv*Tmec+Teq*Tmec Tdrv+Teq+Tmec+Tdrv*Ke*Kt*Kmec*Keq 1+Ke*Kt*Kmec*Keq 0];


P_reduced = tf(P_num_reduced, P_den_reduced);
P_full    = tf(P_num, P_den);
P = P_reduced;

% requests
Mp    = 0.1;
ts    = 0.15;
tr    = 0.01;

% definitions
alpha = 4;

[Kp, Ki, Kd] = getPIDBode(P,100,90,alpha);
s = tf('s');
C = Kp + Ki/s + Kd*s;
W = feedback(C*P,1);
step(W)

function [wgc, phim] = getWPhi(tr, ts, Mp)
    wgc = 2/tr;
    phim = rad2deg(1.04 - 0.8*Mp);
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