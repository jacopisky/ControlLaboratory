run("LoadParameters.m");

s = tf('s');
Cs = pid.Kp + pid.Ki/s + pid.Kd*s/(pid.sim.T_l*s + 1);

Km = drv.dcgain*mot.Kt/(mot.Kt*mot.Ke + (mot.R+sens.curr.Rs)*mech.Beq);
Tm = (mot.R+sens.curr.Rs)*mech.Jeq/(mot.Kt*mot.Ke + (mot.R+sens.curr.Rs)*mech.Beq);
P_num_reduced = [Km];
P_den_reduced = [Tm*gbox.N gbox.N 0];

Ps = tf(P_num_reduced, P_den_reduced);

Cz_tustin = c2d(Cs, sim.T_samp, 'tustin');
Pz_tustin = c2d(Ps, sim.T_samp, 'tustin');

Cz_zoh = c2d(Cs, sim.T_samp, 'zoh');
Pz_zoh = c2d(Ps, sim.T_samp, 'zoh');

figure
step(feedback(Cz_tustin*Pz_tustin,1))
hold on
step(feedback(Cz_zoh*Pz_zoh,1))

clear Km;
clear Tm;
clear P_num_reduced;
clear P_den_reduced;
clear Cz_tustin;
clear Cz_zoh;
clear Pz_zoh;
clear Pz_tustin;
clear Cs;
clear Ps;
clear s;