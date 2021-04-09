run("../LAB0/LoadParameters.m");
run("../LAB0/LoadPID.m");
run("../LAB0/EstimateMechParameters.m");

io.awu.Kw = 5/request.ts;

io.fwd.inertiaCompGain  = gbox.N*simp_model.Req*est_par.J_eq/(drv.dcgain*mot.Kt);
io.fwd.frictionCompGain = mot.Req / (drv.dcgain*mot.Kt*gbox.N);
io.fwd.BEMFCompGain     = gbox.N*mot.Ke/drv.dcgain;

a22 = -est_par.B_eq/est_par.J_eq + mot.Kt*mot.Ke/(est_par.J_eq*mot.Req);
b2  = mot.Kt*drv.dcgain/(mot.Req*est_par.J_eq);
c1  = 1/gbox.N;

ss.plant.A = [0 1; 0 a22];
ss.plant.B = [0;b2];
ss.plant.C = [c1 0];
ss.plant.D = 0;

gains = inv([ss.plant.A ss.plant.B;ss.plant.C ss.plant.D])*[0;0;1];

ss.controller.Nx = [gains(1,1);gains(2,1)];
ss.controller.Nu = gains(3,1);

delta = log(1/request.Mp)/sqrt(pi^2 + log(1/request.Mp)^2);
wn    = 3/(delta*request.ts);

p1    = -delta*wn + 1i*wn*sqrt(1-delta^2);
p2    = -delta*wn - 1i*wn*sqrt(1-delta^2);

ss.controller.K = place(ss.plant.A,ss.plant.B,[p1 p2]);

ss.obs.wc    = 2*pi*50;
ss.obs.delta = 1/sqrt(2);

clear a22;
clear b2;
clear c1;
clear delta;
clear wn;
clear p1;
clear p2;
clear gains;
