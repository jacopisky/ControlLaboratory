%% Preliminary imports
run("../LAB0/LoadParameters.m");
run("LoadHubParameters.m");
run("../LAB0/EstimateMechParameters.m");
run("../LAB1/LoadEstimatedParamsToPlant.m");

mld.Jeq = mot.J + mld.Jh/gbox.N/gbox.N;
mld.Beq = mech.Beq;

%% Estimated Hub Parameters
%run("EstimateHubParameters.m");
%mld.Bb = est_par.Bb;
%mld.k = est_par.k;
