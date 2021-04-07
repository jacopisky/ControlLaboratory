clear all;

run("LoadParameters.m");
run("Requirements.m");

% defining reduced plant (no static friction) (dominant poles)
Beq_simpl = 0;
Req = Ra + Rs;
Km  = Kdrv*Kt/(Req*Beq_simpl+Kt*Ke);
Tm  = Req*Jeq/(Req*Beq_simpl+Kt*Ke);
P_num_reduced = [Km/N];
P_den_reduced = [Tm 1 0];

P = tf(P_num_reduced, P_den_reduced);

[Kp, Ki, Kd] = getPIDBode(P,100,90,alpha);
s = tf('s');
C = Kp + Ki/s + Kd*s;
W = feedback(C*P,1);

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