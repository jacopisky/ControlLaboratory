run("../LoadParameters.m");

delta = log(1/request.Mp)/sqrt(pi^2 + log(1/request.Mp)^2);
wn    = 3/(delta*request.ts);

% robust state space tracking
sigma = -delta*wn;
wd    = wn*sqrt(1-delta^2);

opt5 = [5*sigma+1i*wd, 5*sigma-1i*wd, 3*sigma];

A    = [0, ss.plant.C; [0;0], ss.plant.A];
B    = [0; ss.plant.B];
C    = [0, ss.plant.C];

robust.controller5.K = place(A, B, opt5);

robust.controller5.Ki = robust.controller5.K(1);
robust.controller5.K  = robust.controller5.K(1,2:3);
robust.controller5.Gaw = 0.4;


clear A;
clear B;
clear C;
clear opt5;
clear delta;
clear wn;
clear wd;
clear sigma;