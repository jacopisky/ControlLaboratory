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

%% Estimated Hub Parameters (override)
%run("EstimateHubParameters.m");
%mld.Bb = est_par.Bb;
%mld.k = est_par.k;

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
robust.Ce = [0 1 0 0 0];
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
        elseif real(p) < -sigma && angle(p)-pi/2<phi && angle(p)+pi/2>-phi
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
for k=0:100:5000
    r = rlocus(sysG*sysGp, k);
    r = transpose(r);
    counter = 0;
    for p=r
        if real(p) > sigma && angle(p)<phi && angle(p)>-phi
            counter = counter + 1;
        elseif real(p) < -sigma && angle(p)-pi/2<phi && angle(p)+pi/2>-phi
            counter = counter + 1;
        end
    end
    if counter == 2
        extended_lqr.R = 1/k;
        break;
    end
end
extended_lqr.Q = [0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0];

extended_bryson_lqr.R = 1/daq.dac.vdd/daq.dac.vdd;
extended_bryson_lqr.Q = [
    1/100 0 0 0 0
    0 1/((0.3*50*deg2rad)^2) 0 0 0;
    0 0 1/((pi/36)^2) 0 0;
    0 0 0 0 0;
    0 0 0 0 0
    ];
K = lqr(sysG, extended_bryson_lqr.Q, extended_bryson_lqr.R);
extended_bryson_lqr.controller.K = K(2:5);
extended_bryson_lqr.controller.Ki = K(1);

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