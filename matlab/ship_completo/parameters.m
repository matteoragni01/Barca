%% Main particulars of the KVLCC2 tanker

p.rho = 1025;           % Salt water density

p.Lpp = 7.0;            % Ship length

p.d = 0.46;             % Ship draft

p.nabla = 3.2724;       % Displacement volume of ship

p.mxprimo = 0.022;                          % Added mass of x axis direction (non-dimensionalized)
p.mx = p.mxprimo * (0.5*p.rho*p.Lpp^2*p.d); % Added mass of x axis direction 
p.myprimo = 0.223;                          % Added mass of y axis direction (non-dimensionalized)
p.my = p.myprimo * (0.5*p.rho*p.Lpp^2*p.d); % Added mass of y axis direction 

p.xG = 0.244;           % Longitudinal coordinate of center of gravity of ship

p.Jzprimo = 0.011;                          % Added moment of inertia (non-dimensionalized)
p.Jz = p.Jzprimo * (0.5*p.rho*p.Lpp^4*p.d); % Added moment of inertia

p.Hr = 0.144;           % Rudder span length

p.m = p.rho*p.nabla;    % Mass of the ship

p.IzG = p.m*(0.25*p.Lpp)^2;% Moment of inertia of ship around center of gravity


%% Parameters for hull forces and moments (hydrodynamic derivatives on maneuvering)

p.R0 = 0.022;       % Ship resistant coefficient on straight maneuvering

p.Xvv = -0.040;
p.Xvr = 0.002;
p.Xrr = 0.011;
p.Xvvvv = 0.771;

p.Yv = -0.315;
p.Yr = 0.083;
p.Yvvv = -1.607;
p.Yvvr = 0.379;
p.Yvrr = -0.391;
p.Yrrr = 0.008;

p.Nv = -0.137;
p.Nr = -0.049;
p.Nvvv = -0.030;
p.Nvvr = -0.294;
p.Nvrr = 0.055;
p.Nrrr = -0.013;

%% Parameters for propeller forces and moments

p.Dp = 0.216;           % Propeller diameter
p.xPprimo = -0.48;      % Longitudinal coordiniate of propeller position (non-dimensionalized)
p.xP = p.xPprimo*p.Lpp; % Longitudinal coordiniate of propeller position

% Coefficient representing (KT)
p.k0 = 0.2931;
p.k1 = -0.2753;
p.k2 = -0.1385;

p.tP = 0.220;           % Thrust deduction factor

p.wP0 = 0.40;           % Wake coefficient at propeller position in straight moving

%% Parameters for rudder forces and moments

p.Ar = 0.0539;          % Profile area of movable part of mariner rudder

p.tR = 0.387;           % Steering resistance deduction factor

p.C1 = 2.0;             % Experimental constant representing wake charateristics in maneuvering

p.epsilon = 1.09;       % Ratio of wake coefficient at propeller and rudder positions

p.kappa = 0.50;         

p.lRprimo = -0.710;     % Effective longitudinal coordinate of rudder position in formula of p.betaR (non-dimensionalized)

p.xHprimo = -0.50;          % Longitudinal coordinate of acting point of the additional lateral force (non-dimensionalized)
p.xH = p.xHprimo * p.Lpp;   % Longitudinal coordinate of acting point of the additional lateral force

p.xRprimo = -0.50;          % Longitudinal coordinate of rudder position (non-dimensionalized)
p.xR = p.xRprimo * p.Lpp;   % Longitudinal coordinate of rudder position

p.falpha = 2.747;           % Rudder lift gradient coefficient

p.aH = 0.312;               % Rudder force increase factor

p.eta = p.Dp/p.Hr;          % Ratio of propeller diameter to radder span

%% condizione iniziale
p.x0 = [1.179 0 0 0 0 0]';
p.Delta = 35*pi/180;
p.dDelta = 2.3*pi/180;
%% input costante
p.np = 10.23;


