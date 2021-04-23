run("../LAB1/LoadParameters.m");

% choose the sampling time %
sim.T_samp = 0.001;
%sim.T_samp = 0.01;
%sim.T_samp = 0.05;

s = tf('s');
Cs = pid.Kp + pid.Ki/s + pid.Kd*s/(pid.sim.T_l*s + 1);
L = place(-1/simp_model.Tm,1,5*min(eig(Cs)));

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
sysCd = c2d(sysC, sim.T_samp,'zoh');
[exact_discrete.Phi, exact_discrete.Gam, exact_discrete.H, exact_discrete.J] = ssdata(sysCd);
sysCd = c2d(sysC, sim.T_samp,'tustin');
[tustin.Phi, tustin.Gam, tustin.H, tustin.J] = ssdata(sysCd);

clear L;
clear Cs;
clear s;
clear sysC;
clear sysCd;