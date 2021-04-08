run("../LAB0/LoadParameters.m");
run("../LAB0/LoadPID.m");
run("../LAB0/EstimateMechParameters.m");

awu.Kw = 5/request.ts;

fwd.inertiaCompGain  = gbox.N*simp_model.Req*est_par.J_eq/(drv.dcgain*mot.Kt);
fwd.frictionCompGain = simp_model.Req / (drv.dcgain*mot.Kt*gbox.N);
fwd.BEMFCompGain     = gbox.N*mot.Ke/drv.dcgain;
