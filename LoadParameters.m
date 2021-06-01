%% Preliminary Imports
run("../LAB0/LoadParameters.m");
run("LoadHubParameters.m");
run("../LAB0/EstimateMechParameters.m");
run("../LAB1/LoadEstimatedParamsToPlant.m");

%% Requirements Specification
request.Mp    = 0.3;  % max overshooting percentage
request.ts    = 0.85; % settling time

% definitions
request.sim.alpha = 4;

%% Estimated Parameters Load
mld.Jeq = mot.J + mld.Jh/gbox.N/gbox.N;
mld.Beq = mech.Beq;

%% Estimated Hub Parameters
run("EstimateHubParameters.m");
%% Set Estimated Hub Parameters
mld.Bb = est_par.Bb;
mld.k = est_par.k;

%%  PID Design
% simple model according to LAB3
s = tf('s');
D_tau_prime = mld.Jeq*mld.Jb*s^3+(mld.Jeq*mld.Bb+mld.Jb*mld.Beq)*s^2+(mld.Beq*mld.Bb+mld.k*(mld.Jeq+mld.Jb/gbox.N/gbox.N))*s+mld.k*(mld.Beq+mld.Bb/gbox.N/gbox.N);
num = drv.dcgain*mot.Kt*(mld.Jb*s^2+mld.Bb*s+mld.k);
den = gbox.N*s*(mot.Req*D_tau_prime+mot.Kt*mot.Ke*(mld.Jb*s^2+mld.Bb*s+mld.k));

P = num/den;

[request.wgc, request.phm] = getWgcPhim(request.Mp, request.ts);
request.fc = request.wgc/2/pi;

[pid.Kp, pid.Ki, pid.Kd] = getPIDBodeFreq(P, request.fc, request.phm, request.sim.alpha);
pid.sim.T_l = 1/(10*request.wgc);
pid.Kw      = 5/request.ts;

C = pid.Kp + pid.Ki/s + pid.Kd*s;
%margin(C*P)

%% State Space Definitions
state_space.A = [
    0 0 1 0; 
    0 0 0 1; 
    -mld.k/gbox.N/gbox.N/mld.Jeq mld.k/gbox.N/gbox.N/mld.Jeq -1/mld.Jeq*(mld.Beq+mot.Kt*mot.Ke/mot.Req) 0; 
    mld.k/mld.Jb -mld.k/mld.Jb 0 -mld.Bb/mld.Jb];
state_space.B = [0;0;mot.Kt*drv.dcgain/gbox.N/mld.Jeq/mot.Req;0];
state_space.Bd = [0;0;-1/gbox.N/gbox.N/mld.Jeq;0];
state_space.T = [1 0 0 0; 1 1 0 0; 0 0 1 0; 0 0 1 1];
state_space.A = state_space.T\state_space.A*state_space.T;
state_space.B = state_space.T\state_space.B;
state_space.Bd = state_space.T\state_space.Bd;
state_space.C = [1 0 0 0]*state_space.T;
state_space.D = 0;

state_space.obs.wc = 2*pi*50;
state_space.obs.delta = 1/sqrt(2);

delta = log(1/request.Mp)/sqrt(pi^2 + log(1/request.Mp)^2);
wn    = 3/delta/request.ts;
phi   = atan(sqrt(1-delta^2)/delta);
p1    = wn*exp(1i*(-pi+phi));
p2    = wn*exp(1i*(-pi-phi));
p3    = wn*exp(1i*(-pi+phi/2));
p4    = wn*exp(1i*(-pi-phi/2));
p5    = -wn;
poles = [p1,p2,p3,p4];
state_space.controller.K = acker(state_space.A,state_space.B,poles);

gains = [state_space.A state_space.B;state_space.C state_space.D]\[0;0;0;0;1];
state_space.controller.Nx = [gains(1,1);gains(2,1);gains(3,1);gains(4,1)];
state_space.controller.Nu = gains(5,1);

poles=[p1,p2,p3,p4,p5];
robust.Ae = [0, state_space.C; [0;0;0;0], state_space.A];
robust.Be = [0; state_space.B];
robust.Ce = [1 0 0 0 0];
robust.De = 0;
robust.controller.K = acker(robust.Ae, robust.Be, poles);

robust.controller.Ki = robust.controller.K(1);
robust.controller.K  = robust.controller.K(1,2:5);

sysG = ss(state_space.A,state_space.B,state_space.C,state_space.D);
sysGp = ss(-state_space.A,-state_space.B,state_space.C,state_space.D);
sigma = wn*delta;
for k=0:1:5000
    r = rlocus(sysG*sysGp, k);
    r = transpose(r);
    counter = 0;
    for p=r
        if real(p) > sigma && angle(p)<phi && angle(p)>-phi
            counter = counter + 1;
        elseif real(p) < -sigma && (angle(p)-pi/2<phi || angle(p)+pi/2>-phi)
            counter = counter + 1;
        end
    end
    if counter == 8
        simple_lqr.R = 1/k;
        break;
    end
end
simple_lqr.Q = [0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
simple_lqr.controller.K = lqr(sysG, transpose(state_space.C)*state_space.C, simple_lqr.R);
bryson_lqr.R = 1/daq.dac.vdd/daq.dac.vdd;
bryson_lqr.Q = [
    1/((0.3*50*deg2rad)^2) 0 0 0;
    0 1/((pi/36)^2) 0 0;
    0 0 0 0;
    0 0 0 0
    ];
bryson_lqr.controller.K = lqr(sysG, bryson_lqr.Q, bryson_lqr.R);

sysG = ss(robust.Ae,robust.Be,robust.Ce,robust.De);
sysGp = ss(-robust.Ae,-robust.Be,robust.Ce,robust.De);
sigma = wn*delta;
for k=0:1:3000
    r = rlocus(sysG*sysGp, k);
    r = transpose(r);
    counter = 0;
    for p=r
        if real(p) > sigma && angle(p)<phi && angle(p)>-phi
            counter = counter + 1;
        elseif real(p) < -sigma && (angle(p)-pi/2<phi || angle(p)+pi/2>-phi)
            counter = counter + 1;
        end
    end
    if counter == 6
        extended_lqr.R = 1/k;
        break;
    end
end
extended_lqr.Q = [0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0];
K = lqr(sysG, transpose(robust.Ce)*robust.Ce, extended_lqr.R);
extended_lqr.controller.K = K(2:5);
extended_lqr.controller.Ki = K(1);

extended_bryson_lqr.R = 1/daq.dac.vdd/daq.dac.vdd;
extended_bryson_lqr.Q = [
    1/100 0 0 0 0;
    0 1/((0.3*50*deg2rad)^2) 0 0 0;
    0 0 1/((pi/36)^2) 0 0;
    0 0 0 0 0;
    0 0 0 0 0
    ];
K = lqr(sysG, extended_bryson_lqr.Q, extended_bryson_lqr.R);

extended_bryson_lqr.controller.K = K(2:5);
extended_bryson_lqr.controller.Ki = K(1);
%%
q22 = 1;
p = findResonantPole(state_space.A);
[wn, zeta] = damp(p);
w0 = wn*sqrt(1-zeta^2);
Aq_prime = [0 1; -w0^2 0];
Bq_prime = [0; 1];
Cq_prime = [sqrt(q22)*w0^2 0];
Dq_prime = 0;
state_space.Aq = Aq_prime;
state_space.Bq = [[0;0] Bq_prime [0 0; 0 0]];
state_space.Cq = [[0 0]; Cq_prime; [0 0; 0 0]];
state_space.Dq = diag([1/(5.0*deg2rad), Dq_prime, 0, 0]);
state_space.Aa = [
    state_space.A [0 0; 0 0; 0 0; 0 0];
    state_space.Bq state_space.Aq
];
state_space.Ba = [state_space.B; 0; 0];
state_space.Ca = [state_space.C 0 0];
state_space.Da = 0;
Qa = [
    transpose(state_space.Dq)*state_space.Dq transpose(state_space.Dq)*state_space.Cq;
    transpose(state_space.Cq)*state_space.Dq transpose(state_space.Cq)*state_space.Cq
];
Ra = 1/daq.dac.vdd/daq.dac.vdd;
sysG = ss(state_space.Aa, state_space.Ba, state_space.Ca, state_space.Da);
state_space.controller.Ka = lqr(sysG, Qa, Ra);
gains = [state_space.Aa state_space.Ba; state_space.Ca state_space.Da]\[0;0;0;0;0;0;1];
state_space.controller.Nxa = [gains(1,1);gains(2,1);gains(3,1);gains(4,1);gains(5,1);gains(6,1)];
state_space.controller.Nua = gains(7,1);

robust.Ae = [
    0 state_space.Ca;
    [0;0;0;0;0;0] state_space.Aa
];
robust.Be = [0; state_space.Ba];
robust.Ce = [0 state_space.Ca];
robust.De = 0;
% first trial was 1/100 --> 10
qi = 1/100;
Qe = [
   qi [0 0 0 0 0 0];
   [0;0;0;0;0;0] Qa
];
Re = Ra;

sysG = ss(robust.Ae, robust.Be, robust.Ce, robust.De);
K = lqr(sysG, Qe, Re);
robust.controller.Ke = K(2:7);
robust.controller.Kie = K(1);
%% Workspace Cleanup
clear p;
clear p1;
clear p2;
clear p3;
clear p4;
clear p5;
clear k;
clear r;
clear poles;
clear gains;
clear wn;
clear phi;
clear delta;
clear sysG;
clear sysGp;
clear sigma;
clear K;
clear P;
clear s;
clear D_tau_prime;
clear C;
clear q22;
clear lam;
clear counter;
clear Qa;
clear Qe;
clear Ra;
clear Re;
clear Aq;
clear Cq;
clear Bq;
clear Dq;
clear Aq_prime;
clear Cq_prime;
clear Bq_prime;
clear Dq_prime;
clear num;
clear den;
clear w0;
clear q22;
clear qi;
clear zeta;

%% Definition of Bode Method
% set of utility routines
function [wgc, phim] = getWgcPhim(Mp, ts)
    delta = log(1/Mp)/sqrt(pi^2+log(1/Mp)^2);
    wgc = 3/delta/ts;
    phim = 180/pi*atan(2*delta/sqrt(sqrt(1+4*delta^4)-2*delta^2)); % Phase margin [deg]
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

%% Resonance Detector
% function to get resonance eigenvalues in a matrix
function p=findResonantPole(A)
    lams = eigs(A);
    for i=1:size(lams)
        lam = lams(i);
        for j=1:size(lams)
            if j~=i && lam==conj(lams(j))
                p = lam;
                break;
            end
        end
    end
end