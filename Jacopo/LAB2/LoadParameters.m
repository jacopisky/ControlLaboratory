run("../LAB0/LoadParameters.m");
run("../LAB0/LoadPID.m");
run("../LAB0/EstimateMechParameters.m");
run("../LAB1/LoadEstimatedParamsToPlant.m");

% choose the sampling time!
sim.T_samp = 0.001;
%sim.T_samp = 0.01;
%sim.T_samp = 0.05;

%% IO Model Emulation
io.awu.Kw = 5/request.ts;

io.fwd.inertiaCompGain  = gbox.N*(mot.R+sens.curr.Rs)*est_par.J_eq/(drv.dcgain*mot.Kt);
io.fwd.frictionCompGain = mot.Req / (drv.dcgain*mot.Kt*gbox.N);
io.fwd.BEMFCompGain     = gbox.N*mot.Ke/drv.dcgain;

%% State-Space Emulation
simp_model.Km = drv.dcgain*mot.Kt/(mot.Kt*mot.Ke + (mot.R+sens.curr.Rs)*mech.Beq);
simp_model.Tm = (mot.R+sens.curr.Rs)*mech.Jeq/(mot.Kt*mot.Ke + (mot.R+sens.curr.Rs)*mech.Beq);

state_space.plant.A = [0 1; 0 -1/simp_model.Tm];
state_space.plant.B = [0;simp_model.Km/gbox.N/simp_model.Tm];
state_space.plant.C = [1 0];
state_space.plant.D = 0;

gains = [state_space.plant.A state_space.plant.B;state_space.plant.C state_space.plant.D]\[0;0;1];

state_space.controller.Nx = [gains(1,1);gains(2,1)];
state_space.controller.Nu = gains(3,1);

delta = log(1/request.Mp)/sqrt(pi^2 + log(1/request.Mp)^2);
wn    = 3/(delta*request.ts);

p1    = -delta*wn + 1i*wn*sqrt(1-delta^2);
p2    = conj(p1);

state_space.controller.K = place(state_space.plant.A,state_space.plant.B,[p1 p2]);

Ae = [0 state_space.plant.C;[0;0] state_space.plant.A];
Be = [0;state_space.plant.B];

sigma = -delta*wn;
wd = wn*sqrt(1-delta^2);

% Pole placement
poles = [2*sigma + 1i*wd, 2*sigma - 1i*wd, 2*sigma];

% Ke
Ke = place(Ae,Be,poles);
robust.controller.Ki = Ke(1);
robust.controller.Ke = Ke(1,2:3);

s = tf('s');
Cs = pid.Kp + pid.Ki/s + pid.Kd*s/(pid.sim.T_l*s + 1);

L = place(-1/simp_model.Tm,1,5*min(eig(Cs)));

exact_discrete.Cz = c2d(Cs, sim.T_samp, 'zoh');
[exact_discrete.pid_num, exact_discrete.pid_den] = tfdata(exact_discrete.Cz);
exact_discrete.pid_num=cell2mat(exact_discrete.pid_num);
exact_discrete.pid_den=cell2mat(exact_discrete.pid_den);

state_space.A0 = -1/simp_model.Tm - L;
state_space.B0 =[simp_model.Km/gbox.N/simp_model.Tm,(-1/simp_model.Tm - L)*L];
state_space.C0 = [0;1];
state_space.D0 = [0 1;0 L];

forward_euler.Phi = 1+state_space.A0*sim.T_samp;
forward_euler.Gam = state_space.B0*sim.T_samp;
forward_euler.H = state_space.C0;
forward_euler.J = state_space.D0;

backward_euler.Phi = inv(1-state_space.A0*sim.T_samp);
backward_euler.Gam = backward_euler.Phi * state_space.B0*sim.T_samp;
backward_euler.H = state_space.C0*backward_euler.Phi;
backward_euler.J = state_space.D0 + backward_euler.H*state_space.B0*sim.T_samp;

sysC = ss(state_space.A0, state_space.B0, state_space.C0, state_space.D0);
sysCd = c2d(sysC, sim.T_samp,'tustin');
[tustin.Phi, tustin.Gam, tustin.H, tustin.J] = ssdata(sysCd);

sysC = ss(state_space.A0, state_space.B0, state_space.C0, state_space.D0);
sysCd = c2d(sysC, sim.T_samp,'zoh');
[exact_discrete.Phi, exact_discrete.Gam, exact_discrete.H, exact_discrete.J] = ssdata(sysCd);

%% Direct Digital Design
sysC = ss(state_space.plant.A, state_space.plant.B, state_space.plant.C, state_space.plant.D);
sysCd = c2d(sysC, sim.T_samp,'zoh');
[direct_design.Phi, direct_design.Gam, direct_design.H, direct_design.J] = ssdata(sysCd);

L = place(direct_design.Phi(2,2), direct_design.Phi(1,2), exp(5*min(eig(Cs))*sim.T_samp));

direct_design.Phi0 = direct_design.Phi(2,2) - L*direct_design.Phi(1,2);
direct_design.Gam0 = [direct_design.Gam(2,1)-L*direct_design.Gam(1,1) (direct_design.Phi(2,2)-L*direct_design.Phi(1,2))*L+direct_design.Phi(2,1)-L*direct_design.Phi(1,1)];
direct_design.H0 = [0;1];
direct_design.J0 = [0 1;0 L];

X = linsolve([direct_design.Phi - eye(2),direct_design.Gam;direct_design.H 0],[0;0;1]);
direct_design.Nx = X(1:2,:);
direct_design.Nu = X(3,:);

delta = log(1/request.Mp)/(sqrt(pi*pi + log(1/request.Mp)*log(1/request.Mp)));
wn = 3/(delta*request.ts);

p1_cont = -delta*wn + 1i*wn*sqrt(1-delta^2);
p1_disc = exp(p1_cont*sim.T_samp);
p2_cont = conj(p1_cont);
p2_disc = exp(p2_cont*sim.T_samp);
poles = [p1_disc p2_disc];

direct_design.K = place(direct_design.Phi,direct_design.Gam, poles);

direct_design.Phie = [1 direct_design.H;[0;0] direct_design.Phi];
direct_design.Game = [0;direct_design.Gam];

sigma = -delta*wn;
wd = wn*sqrt(1-delta^2);

p1_cont = 2*sigma + 1i*wd;
p1_disc = exp(p1_cont*sim.T_samp);
p2_cont = 2*sigma - 1i*wd;
p2_disc = exp(p2_cont*sim.T_samp);
p3_cont = 2*sigma;
p3_disc = exp(p3_cont*sim.T_samp);

poles = [p1_disc, p2_disc, p3_disc];

Ke = place(direct_design.Phie, direct_design.Game, poles);
direct_design.Ki = Ke(1);
direct_design.Ke = Ke(1,2:3);

clear L;
clear Cs;
clear s;
clear sysC;
clear sysCd;
clear X;
clear delta;
clear wn;
clear poles;
clear p1_cont;
clear p2_cont;
clear p3_cont;
clear p1_disc;
clear p2_disc;
clear p3_disc;
clear Ke;
clear wd;
clear sigma;
clear gains;
clear p1;
clear p2;
clear Ae;
clear Be;