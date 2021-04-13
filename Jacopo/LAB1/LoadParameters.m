run("../LAB0/LoadParameters.m");
run("../LAB0/LoadPID.m");
run("../LAB0/EstimateMechParameters.m");
run("LoadEstimatedParamsToPlant.m");

io.awu.Kw = 5/request.ts;

io.fwd.inertiaCompGain  = gbox.N*(mot.R+sens.curr.Rs)*est_par.J_eq/(drv.dcgain*mot.Kt);
io.fwd.frictionCompGain = mot.Req / (drv.dcgain*mot.Kt*gbox.N);
io.fwd.BEMFCompGain     = gbox.N*mot.Ke/drv.dcgain;

% state space
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
p2    = conj(p1);

ss.controller.K = place(ss.plant.A,ss.plant.B,[p1 p2]);
ss.obs.wc    = 2*pi*50;
ss.obs.delta = 1/sqrt(2);

% robust state space tracking
sigma = -delta*wn;
wd    = wn*sqrt(1-delta^2);
opt1 = [sigma+1i*wd, sigma-1i*wd, sigma];
opt2 = [sigma, sigma, sigma];
opt3 = [2*sigma+1i*wd, 2*sigma-1i*wd, 2*sigma];
opt4 = [2*sigma+1i*wd, 2*sigma-1i*wd, 3*sigma];

opt5 = [5*sigma+1i*wd, 5*sigma-1i*wd, 3*sigma];

A    = [0, ss.plant.C; [0;0], ss.plant.A];
B    = [0; ss.plant.B];
C    = [0, ss.plant.C];

robust.controller1.K = place(A, B, opt1);
robust.controller2.K = acker(A, B, opt2);
robust.controller3.K = place(A, B, opt3);
robust.controller4.K = place(A, B, opt4);
robust.controller5.K = place(A, B, opt5);

robust.controller1.Ki = robust.controller1.K(1);
robust.controller1.K  = robust.controller1.K(1,2:3);

robust.controller2.Ki = robust.controller2.K(1);
robust.controller2.K  = robust.controller2.K(1,2:3);

robust.controller3.Ki = robust.controller3.K(1);
robust.controller3.K  = robust.controller3.K(1,2:3);

robust.controller4.Ki = robust.controller4.K(1);
robust.controller4.K  = robust.controller4.K(1,2:3);

robust.controller5.Ki = robust.controller5.K(1);
robust.controller5.K  = robust.controller5.K(1,2:3);

% error space
Tr1 = 0.15;
Tr2 = 0.25;
Tr3 = 0.5;
Tr4 = 1;

p1  = wn*exp(1i*(-pi+pi/4));
p2  = conj(p1);
p3  = wn*exp(1i*(-pi+pi/6));
p4  = conj(p3);
p5  = -wn;

es.controller1.w0 = 2*pi/Tr1;
es.controller2.w0 = 2*pi/Tr2;
es.controller3.w0 = 2*pi/Tr3;
es.controller4.w0 = 2*pi/Tr4;

A1  = [0,1,0;0,0,1;0,-es.controller1.w0^2,0];
A2  = [0,1,0;0,0,1;0,-es.controller2.w0^2,0];
A3  = [0,1,0;0,0,1;0,-es.controller3.w0^2,0];
A4  = [0,1,0;0,0,1;0,-es.controller4.w0^2,0];

Az1   = [A1,[[0,0;0,0];ss.plant.C];[[0,0,0;0,0,0],ss.plant.A]];
Az2   = [A2,[[0,0;0,0];ss.plant.C];[[0,0,0;0,0,0],ss.plant.A]];
Az3   = [A3,[[0,0;0,0];ss.plant.C];[[0,0,0;0,0,0],ss.plant.A]];
Az4   = [A4,[[0,0;0,0];ss.plant.C];[[0,0,0;0,0,0],ss.plant.A]];
Bz    = [0;0;0;ss.plant.B];
poles = [p1,p2,p3,p4,p5];

es.controller1.Kz = place(Az1, Bz, poles);
es.controller2.Kz = place(Az2, Bz, poles);
es.controller3.Kz = place(Az3, Bz, poles);
es.controller4.Kz = place(Az4, Bz, poles);
es.controller1.Kpsi = es.controller1.Kz(1, 4:5);
es.controller1.Kz = es.controller1.Kz(1, 1:3);
es.controller2.Kpsi = es.controller2.Kz(1, 4:5);
es.controller2.Kz = es.controller2.Kz(1, 1:3);
es.controller3.Kpsi = es.controller3.Kz(1, 4:5);
es.controller3.Kz = es.controller3.Kz(1, 1:3);
es.controller4.Kpsi = es.controller4.Kz(1, 4:5);
es.controller4.Kz = es.controller4.Kz(1, 1:3);


tmp = es.controller1.Kz(1,3);
es.controller1.Kz(1,3) = es.controller1.Kz(1,1);
es.controller1.Kz(1,1) = tmp;
tmp = es.controller2.Kz(1,3);
es.controller2.Kz(1,3) = es.controller2.Kz(1,1);
es.controller2.Kz(1,1) = tmp;
tmp = es.controller3.Kz(1,3);
es.controller3.Kz(1,3) = es.controller3.Kz(1,1);
es.controller3.Kz(1,1) = tmp;
tmp = es.controller4.Kz(1,3);
es.controller4.Kz(1,3) = es.controller4.Kz(1,1);
es.controller4.Kz(1,1) = tmp;

% extended error space
pc1 = -delta*wn + 1i*wn*sqrt(1-delta^2);
pc2 = conj(pc1);

p1 = 2*wn*exp(1i*(-pi+pi/3));
p2 = conj(p1);
p3 = 2*wn*exp(1i*(-pi+pi/6));
p4 = conj(p3);
p5 = -2*wn;

ees.controller1.w0 = 2*pi/Tr1;
ees.controller2.w0 = 2*pi/Tr2;
ees.controller3.w0 = 2*pi/Tr3;
ees.controller4.w0 = 2*pi/Tr4;

ees.Croh = [1 0 0];

ees.Ae1   = [A1, [0 0; 0 0; 0 0];ss.plant.B*ees.Croh, ss.plant.A];
ees.Ae2   = [A2, [0 0; 0 0; 0 0];ss.plant.B*ees.Croh, ss.plant.A];
ees.Ae3   = [A3, [0 0; 0 0; 0 0];ss.plant.B*ees.Croh, ss.plant.A];
ees.Ae4   = [A4, [0 0; 0 0; 0 0];ss.plant.B*ees.Croh, ss.plant.A];
ees.Be = [0;0;0;ss.plant.B];
ees.Ce = [0, 0, 0, ss.plant.C];

pcs   = [pc1, pc2];
poles = [p1,p2,p3,p4,p5];

ees.estimator1.Le = transpose(place(transpose(ees.Ae1), transpose(ees.Ce), poles));
ees.estimator2.Le = transpose(place(transpose(ees.Ae2), transpose(ees.Ce), poles));
ees.estimator3.Le = transpose(place(transpose(ees.Ae3), transpose(ees.Ce), poles));
ees.estimator4.Le = transpose(place(transpose(ees.Ae4), transpose(ees.Ce), poles));

ees.controller = place(ss.plant.A, ss.plant.B, pcs);

clear a22;
clear b2;
clear c1;
clear delta;
clear wn;
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
clear poles;
clear p1;
clear p2;
clear p3;
clear p4;
clear p5;
clear Az1;
clear Az2;
clear Az3;
clear Az4;
clear Bz;
clear A1;
clear A2;
clear A3;
clear A4;
clear tmp;
clear Tr1;
clear Tr2;
clear Tr3;
clear Tr4;
clear pc1;
clear pc2;
clear pcs;