%% Preliminary imports
run("../LAB0/LoadParameters.m");
run("LoadHubParameters.m");
run("../LAB0/EstimateMechParameters.m");
run("../LAB1/LoadEstimatedParamsToPlant.m");

request.Mp    = 0.3;  % max overshooting percentage
request.ts    = 0.85; % settling time

mld.Jeq = mot.J + mld.Jh/gbox.N/gbox.N;
mld.Beq = mech.Beq;

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
poles = [p1,p2,p3,p4];
state_space.controller.K = acker(state_space.A,state_space.B,poles);

gains = [state_space.A state_space.B;state_space.C state_space.D]\[0;0;0;0;1];
state_space.controller.Nx = [gains(1,1);gains(2,1);gains(3,1);gains(4,1)];
state_space.controller.Nu = gains(5,1);

%% Estimated Hub Parameters
%run("EstimateHubParameters.m");
%mld.Bb = est_par.Bb;
%mld.k = est_par.k;

clear p1;
clear p2;
clear p3;
clear p4;
clear poles;
clear gains;
clear wn;
clear phi;
clear delta;
