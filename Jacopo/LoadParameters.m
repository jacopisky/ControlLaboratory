% KNOWN PLANT PARAMETERS

%% Electrical
% input filter & motor driver
R1  = 7.5e3;
R2  = 1.6e3;
R3  = 1.2e3;
R4  = 500;
C1  = 100e-9;
vp  = 12;
vm  = -12;

% DAC parameters
vdd    = 10;
vcc    = -10;
n_bits = 16;
q_step = (vdd-vcc)/(2^n_bits-1);

% shunt current resistor
Rs = 0.5;

% electrical part of motor
Ra = 2.6;
La = 180e-6;

%% Electro-mechanical part of the motor
Kt = 7.68e-3;
Ke = 7.68e-3;

%% Mechanical
% encoder
n_windows = 1024;
enc_step  = 360/(n_windows*4);

% motor
Jm = 3.9e-7;
Bm = 0;

% gear box and load
Bl  = 0;
Jd  = 3.0e-5;
J72 = 1.4e-6;
N   = 14;


%% Calculated parameters
Kdrv = ((R3+R4)*R2)/(R4*(R1+R2));
Tdrv = (R1*R2*C1)/(R1+R2);
Kmot = 1/(Rs+Ra);
Tmot = La/(Rs+Ra);
Jeq  = Jm+(Jd+3*J72)/(N*N);
Beq  = Bm+Bl/(N*N);

%override
Beq   = 2e-6;
tausf = 1e-2;
%end override

Keq  = 1/Beq;
Teq  = Jeq/Beq;

%% Convertions
rad2rpm    = 60/2/pi;
rad2deg    = 180/pi;
deg2rad    = pi/180;
pulse2deg  = enc_step;
deg2pulse  = 1/pulse2deg;

%% Simulation
Tsample = 1e-3;
Tl = 0.002;
