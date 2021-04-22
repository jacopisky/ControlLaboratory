run("../LAB0/LoadParameters.m");
run("../LAB0/LoadPID.m");
run("../LAB0/EstimateMechParameters.m");
run("../LAB1/LoadEstimatedParamsToPlant.m");

io.awu.Kw = 5/request.ts;

io.fwd.inertiaCompGain  = gbox.N*(mot.R+sens.curr.Rs)*est_par.J_eq/(drv.dcgain*mot.Kt);
io.fwd.frictionCompGain = mot.Req / (drv.dcgain*mot.Kt*gbox.N);
io.fwd.BEMFCompGain     = gbox.N*mot.Ke/drv.dcgain;

sim.T_samp = 0.001;
%sim.T_samp = 0.01;
%sim.T_samp = 0.05;

% state space
Km = drv.dcgain*mot.Kt/(mot.Kt*mot.Ke + (mot.R+sens.curr.Rs)*mech.Beq);
Tm = (mot.R+sens.curr.Rs)*mech.Jeq/(mot.Kt*mot.Ke + (mot.R+sens.curr.Rs)*mech.Beq);

ss.plant.A = [0 1; 0 -1/Tm];
ss.plant.B = [0;Km/gbox.N/Tm];
ss.plant.C = [1 0];
ss.plant.D = 0;

s = tf('s');
Cs = pid.Kp + pid.Ki/s + pid.Kd*s/(pid.sim.T_l*s + 1);

L = place(-1/Tm,1,5*min(eig(Cs)));
red_ord.obs.A0 = -1/Tm - L;
red_ord.obs.B0 =[Km/gbox.N/Tm,(-1/Tm - L)*L];
red_ord.obs.C0 = [0;1];
red_ord.obs.D0 = [0 1;0 L];

delta = log(1/request.Mp)/sqrt(pi^2 + log(1/request.Mp)^2);
wn    = 3/(delta*request.ts);

p1    = -delta*wn + 1i*wn*sqrt(1-delta^2);
p2    = conj(p1);

gains = [ss.plant.A ss.plant.B;ss.plant.C ss.plant.D]\[0;0;1];
ss.controller.Nx = [gains(1,1);gains(2,1)];
ss.controller.Nu = gains(3,1);

ss.controller.K = place(ss.plant.A,ss.plant.B,[p1 p2]);

red_ord.discrete.obs.A0 = 1+red_ord.obs.A0*sim.T_samp;
red_ord.discrete.obs.B0 = red_ord.obs.B0*sim.T_samp;
red_ord.discrete.obs.C0 = red_ord.obs.C0;
red_ord.discrete.obs.D0 = red_ord.obs.D0;

clear Tm;
clear Km;
clear gains;
clear p1;
clear p2;
clear delta;
clear wn;
clear L;
clear Cs;
clear s;