run("../LAB1/LoadParameters.m");

% choose the sampling time %
sim.T_samp = 0.001;
%sim.T_samp = 0.01;
%sim.T_samp = 0.05;

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

sysC = ss(state_space.plant.A, state_space.plant.B, state_space.plant.C, state_space.plant.D);
sysCd = c2d(sysC, sim.T_samp,'zoh');
[direct_design.Phi, direct_design.Gam, direct_design.H, direct_design.J] = ssdata(sysCd);


L = place(direct_design.Phi(2,2), direct_design.Phi(1,2), exp(5*min(eig(Cs))*sim.T_samp));

direct_design.Phi0 = direct_design.Phi(2,2) - L*direct_design.Phi(1,2);
direct_design.Gam0 = [direct_design.Gam(2,1)-L*direct_design.Gam(1,1) (direct_design.Phi(2,2)-L*direct_design.Phi(1,2))*L+direct_design.Phi(2,1)-L*direct_design.Phi(1,1)];
direct_design.H0 = [0;1];
direct_design.J0 = [0 1;0 L];

% Nominal Nx Nu [DISCRETE]
X = linsolve([direct_design.Phi - eye(2),direct_design.Gam;direct_design.H 0],[0;0;1]);
direct_design.Nx = X(1:2,:);
direct_design.Nu = X(3,:);

% Dominant Pole Approximation
delta = log(1/request.Mp)/(sqrt(pi*pi + log(1/request.Mp)*log(1/request.Mp)));
wn = 3/(delta*request.ts);

% Pole placement [DISCRETE]
p1_cont = -delta*wn + 1i*wn*sqrt(1-delta^2);
p1_disc = exp(p1_cont*sim.T_samp);
p2_cont = conj(p1_cont);
p2_disc = exp(p2_cont*sim.T_samp);
poles = [p1_disc p2_disc];

direct_design.K = place(direct_design.Phi,direct_design.Gam, poles);

direct_design.Ki = robust.controller3.Ki;

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
clear p1_disc;
clear p2_disc;