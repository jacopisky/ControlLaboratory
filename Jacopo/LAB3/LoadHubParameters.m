%% Mechanical load (mld) nominal parameters
% hub params
mld.Jh = 6.84e-4; % moment of inertia
mld.Bh = 2.5e-4;  % viscous friction coeff

% beam & elastic joint params
mld.Jb = 1.4e-3; % moment of inertia
mld.d = 5.0e-2; % flex joint damping coeff (estimated)
mld.wn = 24.5; % flex joint natural freq (estimated)
mld.Bb = mld.Jb * 2*mld.d*mld.wn; % beam viscous friction
mld.k = mld.Jb * mld.wn^2; % flex joint stiffness

% potentiometer 2 (Spectrol 357−0−0−103) − installed on hub
sens.pot2.range.R = 10e3; % ohmic value range
sens.pot2.range.V = 5; % voltage range
sens.pot2.range.th_deg = 340; % angle range [deg]
sens.pot2.range.th = sens.pot2.range.th_deg * deg2rad; % angle range [rad]
sens.pot2.deg2V = sens.pot2.range.V / sens.pot2.range.th_deg; % sensitivity [V/deg]
sens.pot2.rad2V = sens.pot2.range.V / sens.pot2.range.th; % sensitivity [V/rad]
sens.pot2.V2deg = 1/sens.pot2.deg2V; % conversion gain [V] −> [deg]
sens.pot2.V2rad = 1/sens.pot2.rad2V; % conversion gain [V] −> [rad]
sens.pot2.noise.var = 3.5e-7;