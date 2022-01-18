%% General parameters and conversion gains
% conversion gains
rpm2rads = 2*pi/60;  % [rpm]----->[rad/s]
rads2rpm = 60/2/pi;  % [rad/s]--->[rpm]
rpm2degs = 360/60;   % [rpm]----->[deg/s]
degs2rpm = 60/360;   % [deg/s]--->[rpm]
deg2rad  = pi/180;   % [deg]----->[rad]
rad2deg  = 180/pi;   % [rad]----->[deg]
ozin2Nm  = 0.706e-2; % [oz*inch]->[N*m]

%% Sensors data
% shunt resistor
sens.curr.Rs = 0.5;

% Hewlett-Packard HEDS-5540#A06 optical encoder
sens.enc.ppr       = 1024*4;            % pulses per rotation
sens.enc.pulse2deg = 360/sens.enc.ppr;  % [pulses]->[deg]
sens.enc.pulse2rad = 2*pi/sens.enc.ppr; % [pulses]->[rad]
sens.enc.deg2pulse = sens.enc.ppr/360;  % [deg]---->[pulses]
sens.enc.rad2pulse = sens.enc.ppr/2/pi; % [rad]---->[pulses]
sens.enc.T_s       = 1e-3;              %  sampling time T_s
sens.enc.q         = 360/(1024*4);      % encoder quantization step

% potentiometer 1 (Spectrol 138-0-0-103) installed on motor box
sens.pot1.range.R      = 10e3;                                       % ohmic value range
sens.pot1.range.V      = 5;                                          % voltage range
sens.pot1.range.th_deg = 345;                                        % angle range [deg]
sens.pot1.range.th     = sens.pot1.range.th_deg * deg2rad;           % angle range [rad]
sens.pot1.deg2V        = sens.pot1.range.V / sens.pot1.range.th_deg; % sensitivity [V/deg]
sens.pot1.rad2V        = sens.pot1.range.V / sens.pot1.range.th;     % sensitivity [V/rad]
sens.pot1.V2deg        = 1/sens.pot1.deg2V;                          % conversion gain [V]->[deg]
sens.pot1.V2rad        = 1/sens.pot1.rad2V;                          % conversion gain [V]->[rad]

%% DC motor nominal parameters
% brushed DC-motor Faulhaber 2338S006S
mot.R    = 2.6;                 % armature resistance (R_a)
mot.L    = 180e-6;              % armature inductance (L_a)
mot.Kt   = 1.088 * ozin2Nm;     % torque constant
mot.Ke   = 0.804e-3 * rads2rpm; % back EMF constant
mot.J    = 5.523e-5 * ozin2Nm;  % rotor inertia (J_m)
mot.B    = 0.0;                 % viscous friction coeff (n.a.) (B_m)
mot.eta  = 0.69;                % motor efficiency
mot.PN   = 3.23/mot.eta;        % nominal output power
mot.UN   = 6;                   % nominal voltage
mot.IN   = mot.PN/mot.UN;       % nominal current
mot.tauN = mot.Kt*mot.IN;       % nominal torque
mot.taus = 2.42 * ozin2Nm;      % stall torque
mot.w0   = 7200 * rpm2rads;     % no load speed

mot.Req  = mot.R + sens.curr.Rs;

%% Gearbox nominal parameters
% planetary gearbox Micromotor SA 23/1
gbox.N1   = 14;   % 1st reduction ratio (planetary gearbox)
gbox.eta1 = 0.80; % gearbox efficiency

% external transmission gears
gbox.N2   = 1;      % 2nd reduction ratio (external trasmission gears)
gbox.J72  = 1.4e-6; % inertia of a single external 72 tooth gear
gbox.eta2 = 1;      % external trasmission efficiency (n.a.)

% overall gearbox data
gbox.N   = gbox.N1*gbox.N2;     % total reduction ratio
gbox.eta = gbox.eta1*gbox.eta2; % total efficiency
gbox.J   = 3*gbox.J72;          % total inertia (at gearbox output)

%% Mechanical load nominal parameters

% inertia disc params
mld.JD = 3e-5; % load disc inertia
mld.BD = 0.0;  % load viscous coeff (n.a.)

% overall mech load params
mld.J     = mld.JD + gbox.J; % total inertia (J_l)
mld.B     = 2.5e-4;          % total viscous fric coeff (estimated)
mld.tausf = 1.0e-2;          % total static friction (estimated)

%% Overall mechanical parameters
mech.Jeq = mot.J+mld.J/(gbox.N)^2;
mech.Beq = mot.B+mld.B/(gbox.N)^2;
%override 
mech.Beq = 2e-6; % (estimated)

%% Voltage driver nominal parameters
% op-amp circuit params
drv.R1 = 7.5e3;  % op-amp input resistor (dac to non-inverting in)
drv.R2 = 1.6e3;  % op-amp input resistor (non-inverting in to gnd)
drv.R3 = 1.2e3;  % op-amp feedback resistor (output to inverting in)
drv.R4 = 0.5e3;  % op-amp feedback resistor (inverting in to gnd)
drv.C1 = 100e-9; % op-amp input capacitor
drv.outmax = 12; % op-amp max output voltage
drv.vp     = 12;
drv.vm     = -12;

% voltage driver dc gain
drv.dcgain = drv.R2/(drv.R1+drv.R2) * (1 + drv.R3/drv.R4); % K_drv

% voltage driver time constant
drv.Tc = drv.C1 * drv.R1*drv.R2/(drv.R1+drv.R2); % T_drv

%% Data acquisition board (daq) data
% NI PCI-6221 DAC data
daq.dac.bits = 16;                              % resolution (bits)
daq.dac.fs   = 10;                              % full scale
daq.dac.q    = 2*daq.dac.fs/(2^daq.dac.bits-1); % quantization
daq.dac.T_s  = 1e-3;
daq.dac.vdd  = 10;
daq.dac.vss  = -10;

% NI PCI-6221 ADC data
daq.adc.bits = 16;                              % resolution (bits)
daq.adc.fs   = 10;                              % full scale (as set in SLDRT Analog Input block)
daq.adc.q    = 2*daq.adc.fs/(2^daq.adc.bits-1); % quantization
daq.adc.T_s  = 1e-3;
daq.adc.vdd  = 10;
daq.adc.vss  = -10;
