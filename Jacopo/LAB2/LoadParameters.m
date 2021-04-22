run("../LAB0/LoadParameters.m");
run("../LAB0/LoadPID.m");
run("../LAB0/EstimateMechParameters.m");
run("../LAB1/LoadEstimatedParamsToPlant.m");

io.awu.Kw = 5/request.ts;

io.fwd.inertiaCompGain  = gbox.N*(mot.R+sens.curr.Rs)*est_par.J_eq/(drv.dcgain*mot.Kt);
io.fwd.frictionCompGain = mot.Req / (drv.dcgain*mot.Kt*gbox.N);
io.fwd.BEMFCompGain     = gbox.N*mot.Ke/drv.dcgain;

%sim.T_samp = 0.001;
sim.T_samp = 0.01;
%sim.T_samp = 0.05;

% state space
Km = drv.dcgain*mot.Kt/(mot.Kt*mot.Ke + (mot.R+sens.curr.Rs)*mech.Beq);
Tm = (mot.R+sens.curr.Rs)*mech.Jeq/(mot.Kt*mot.Ke + (mot.R+sens.curr.Rs)*mech.Beq);

ss.plant.A = [0 1; 0 -1/Tm];
ss.plant.B = [0;Km/gbox.N/Tm];
ss.plant.C = [1 0];
ss.plant.D = 0;

clear Tm;
clear Km;