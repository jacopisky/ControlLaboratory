clear all;

run("../LAB0/LoadParameters.m");
run("../LAB0/Requirements.m");
run("../LAB0/EstimationAndController.m");

run("FromEstimationToSet.m");

Tw = ts/5;
Kw = 1/Tw;

inertiaCompGain  = N*Req*Jeq_est/(Kdrv*Kt);
frictionCompGain = Req / (Kdrv*Kt*N);
BEMFCompGain     = N*Ke/Kdrv;
