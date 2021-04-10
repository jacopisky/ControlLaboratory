run("../LAB0/LoadParameters.m");
run("../LAB0/LoadPID.m");
run("../LAB0/EstimateMechParameters.m");
run("LoadEstimatedParamsToPlant.m");

io.awu.Kw = 5/request.ts;

io.fwd.inertiaCompGain  = gbox.N*simp_model.Req*est_par.J_eq/(drv.dcgain*mot.Kt);
io.fwd.frictionCompGain = mot.Req / (drv.dcgain*mot.Kt*gbox.N);
io.fwd.BEMFCompGain     = gbox.N*mot.Ke/drv.dcgain;

a22 = -(est_par.B_eq*mot.Req + mot.Kt*mot.Ke)/(est_par.J_eq*mot.Req);
b2  = mot.Kt*drv.dcgain/(mot.Req*est_par.J_eq*gbox.N);
c1  = 1;

ss.plant.A = [0 1; 0 a22];
ss.plant.B = [0;b2];
ss.plant.C = [c1 0];
ss.plant.D = 0;

gains = [ss.plant.A ss.plant.B;ss.plant.C ss.plant.D]\[0;0;1];

ss.controller.Nx = [gains(1,1);gains(2,1)];
ss.controller.Nu = gains(3,1);

delta = log(1/request.Mp)/sqrt(pi^2 + log(1/request.Mp)^2);
wn    = 3/(delta*request.ts);

p1    = -delta*wn + 1i*wn*sqrt(1-delta^2);
p2    = -delta*wn - 1i*wn*sqrt(1-delta^2);

ss.controller.K = place(ss.plant.A,ss.plant.B,[p1 p2]);
ss.obs.wc    = 2*pi*50;
ss.obs.delta = 1/sqrt(2);

sigma = -delta*wn;
wd    = wn*sqrt(1-delta^2);
opt1 = [sigma+1i*wd, sigma-1i*wd, sigma];
opt2 = [sigma, sigma, sigma];
opt3 = [2*sigma+1i*wd, 2*sigma-1i*wd, 2*sigma];
opt4 = [2*sigma+1i*wd, 2*sigma-1i*wd, 3*sigma];
A    = [0, ss.plant.C; [0;0], ss.plant.A];
B    = [0; ss.plant.B];
C    = [0, ss.plant.C];

robust.controller1.K = place(A, B, opt1);
robust.controller2.K = acker(A, B, opt2);
robust.controller3.K = place(A, B, opt3);
robust.controller4.K = place(A, B, opt4);

robust.controller1.Ki = robust.controller1.K(1);
robust.controller1.K  = robust.controller1.K(1,2:3);

robust.controller2.Ki = robust.controller2.K(1);
robust.controller2.K  = robust.controller2.K(1,2:3);

robust.controller3.Ki = robust.controller3.K(1);
robust.controller3.K  = robust.controller3.K(1,2:3);

robust.controller4.Ki = robust.controller4.K(1);
robust.controller4.K  = robust.controller4.K(1,2:3);

clear a22;
clear b2;
clear c1;
clear delta;
clear wn;
clear p1;
clear p2;
clear gains;
clear sigma;
clear wd;
clear opt1;
clear opt2;
clear opt3;
clear opt4;
clear A;
clear B;
clear C;