run("../LoadParameters.m");

opt5 = [-60, -98, -99];

A    = [0, ss.plant.C; [0;0], ss.plant.A];
B    = [0; ss.plant.B];
C    = [0, ss.plant.C];

robust.controller5.K = place(A, B, opt5);

robust.controller5.Ki = robust.controller5.K(1);
robust.controller5.K  = robust.controller5.K(1,2:3);
robust.controller5.Gaw = 35.2;


clear A;
clear B;
clear C;
clear opt5;