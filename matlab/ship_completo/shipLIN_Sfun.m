function [sys,x0,str,ts,simStateCompliance] = shipLIN_Sfun(t,x,u,flag,p,X)

switch flag,

  %%%%%%%%%%%%%%%%%%
  % Initialization %
  %%%%%%%%%%%%%%%%%%
  case 0,
    [sys,x0,str,ts,simStateCompliance]=mdlInitializeSizes(X);

  %%%%%%%%%%%%%%%
  % Derivatives %
  %%%%%%%%%%%%%%%
  case 1,
    sys=mdlDerivatives(t,x,u,p);

  %%%%%%%%%%
  % Update %
  %%%%%%%%%%
  case 2,
    sys=mdlUpdate(t,x,u);

  %%%%%%%%%%%
  % Outputs %
  %%%%%%%%%%%
  case 3,
    sys=mdlOutputs(t,x,u);

  %%%%%%%%%%%%%%%%%%%%%%%
  % GetTimeOfNextVarHit %
  %%%%%%%%%%%%%%%%%%%%%%%
  case 4,
    sys=mdlGetTimeOfNextVarHit(t,x,u);

  %%%%%%%%%%%%%
  % Terminate %
  %%%%%%%%%%%%%
  case 9,
    sys=mdlTerminate(t,x,u);

  %%%%%%%%%%%%%%%%%%%%
  % Unexpected flags %
  %%%%%%%%%%%%%%%%%%%%
  otherwise
    DAStudio.error('Simulink:blocks:unhandledFlag', num2str(flag));

end

% end sfuntmpl

%
%=============================================================================
% mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function.
%=============================================================================
%
function [sys,x0,str,ts,simStateCompliance]=mdlInitializeSizes(xInit)

%
% call simsizes for a sizes structure, fill it in and convert it to a
% sizes array.
%
% Note that in this example, the values are hard coded.  This is not a
% recommended practice as the characteristics of the block are typically
% defined by the S-function parameters.
%
sizes = simsizes;

sizes.NumContStates  = 5;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 5;
sizes.NumInputs      = 1;
sizes.DirFeedthrough = 0;
sizes.NumSampleTimes = 1;   % at least one sample time is needed

sys = simsizes(sizes);

%
% initialize the initial conditions
%
x0  = xInit;


%
% str is always an empty matrix
%
str = [];

%
% initialize the array of sample times
%
ts  = [0 0];

% Specify the block simStateCompliance. The allowed values are:
%    'UnknownSimState', < The default setting; warn and assume DefaultSimState
%    'DefaultSimState', < Same sim state as a built-in block
%    'HasNoSimState',   < No sim state
%    'DisallowSimState' < Error out when saving or restoring the model sim state
simStateCompliance = 'UnknownSimState';

% end mdlInitializeSizes

%
%=============================================================================
% mdlDerivatives
% Return the derivatives for the continuous states.
%=============================================================================
%
function sys=mdlDerivatives(t,x,u,p)
    
    U = sqrt(x(1)^2 + x(2)^2);
    vmprimo = x(2) / U;
    rprimo = x(3) * p.Lpp / U;
    
    beta = atan(-x(2) / x(1));
    betaP = beta - p.xPprimo * rprimo;
    
    if(betaP > 0)
        C2 = 1.6;
    else
        C2 = 1.1;
    end
    
    wp = -(1 - exp(-p.C1 * abs(betaP))) * (C2 - 1) * (1 - p.wP0);
    Jp = x(1) * (1 - wp) / (p.Dp * p.np);
    betaR = beta - p.lRprimo * rprimo;
    
    if(betaR > 0)
        gammaR = 0.640;
    else
        gammaR = 0.395;
    end
    
    vr = U * gammaR * betaR;
    
    Kt = p.k2 * Jp^2 + p.k1 * Jp + p.k0;
    ur = p.epsilon * x(1) * (1 - wp) * sqrt(p.eta * (1 + p.kappa * (sqrt(1 + (8 * Kt / (pi * Jp^2))) - 1)) + (1 - p.eta));

    %X
    Xp = (1 - p.tP) * p.rho * p.np^2 * p.Dp^4 * Kt;
    Xr = -(1 - p.tR) * 0.5 * p.rho * p.Ar * (ur^2 + vr^2) * p.falpha * sin(u(1) - (vr / ur)) * sin(u(1));
    Xh = 0.5 * p.rho * p.Lpp * p.d * U^2 * (-p.R0 + p.Xvv * vmprimo^2 + p.Xvr * vmprimo * rprimo + p.Xrr * rprimo^2 + p.Xvvvv * vmprimo^4);
    
    X = Xh + Xr + Xp;
    
    %Y
    
    Yr = -(1 - p.aH) * 0.5 * p.rho * p.Ar * (ur^2 + vr^2)* p.falpha * sin(u(1) - (vr / ur)) * cos(u(1));
    Yh = 0.5 * p.rho * p.Lpp * p.d * U^2 * (p.Yv * vmprimo + p.Yr * rprimo + p.Yvvv* vmprimo^3 + p.Yvvr * vmprimo^2 * rprimo + p.Yvrr * vmprimo * rprimo^2 + p.Yrrr * rprimo^3);
    
    Y = Yh + Yr;
    
    %Nm
    
    Nr = -(p.xR - p.aH * p.xH) * 0.5 * p.rho * p.Ar * (ur^2 + vr^2) * p.falpha * sin(u(1) - (vr / ur)) * cos(u(1));
    Nh = 0.5 * p.rho * p.Lpp * p.d * U^2 * (p.Nv* vmprimo + p.Nr * rprimo + p.Nvvv* vmprimo^3 + p.Nvvr * vmprimo^2 * rprimo + p.Nvrr * vmprimo * rprimo^2 + p.Nrrr * rprimo^3);
    
    Nm = Nh + Nr;
    
    k = (p.IzG + p.xG^2 * p.m + p.Jz);
    
    %EQUAZIONI
    sys(1) = (X + (p.m + p.my) * x(2) * x(3) + p.xG * p.m * x(3)^2) / (p.m + p.mx);
    sys(2) = (Y - (p.m + p.mx) * x(1) * x(3) - p.xG * p.m * (Nm - p.xG * p.m * (x(1) * x(3))) / k ) / (p.m + p.my)*(1-(p.xG^2*p.m^2)/(k*(p.m + p.my)));
    sys(3) = (Nm - p.xG * p.m * (sys(2) + x(1) * x(3))) / k ;
    sys(4) = x(3);
    %sys(5) = x(1) * cos(x(4)) - x(2) * sin(x(4));
    sys(5) = x(2) * cos(x(4)) + x(1) * sin(x(4));


% end mdlDerivatives

%
%=============================================================================
% mdlUpdate
% Handle discrete state updates, sample time hits, and major time step
% requirements.
%=============================================================================
%
function sys=mdlUpdate(t,x,u)

sys = [];

% end mdlUpdate

%
%=============================================================================
% mdlOutputs
% Return the block outputs.
%=============================================================================
%
function sys=mdlOutputs(t,x,u)
    sys(1) = x(1);
    sys(2) = x(2);
    sys(3) = x(3);
    sys(4) = x(4);
    sys(5) = x(5);
    

% end mdlOutputs

%
%=============================================================================
% mdlGetTimeOfNextVarHit
% Return the time of the next hit for this block.  Note that the result is
% absolute time.  Note that this function is only used when you specify a
% variable discrete-time sample time [-2 0] in the sample time array in
% mdlInitializeSizes.
%=============================================================================
%
function sys=mdlGetTimeOfNextVarHit(t,x,u)

sampleTime = 1;    %  Example, set the next hit to be one second later.
sys = t + sampleTime;

% end mdlGetTimeOfNextVarHit

%
%=============================================================================
% mdlTerminate
% Perform any end of simulation tasks.
%=============================================================================
%
function sys=mdlTerminate(t,x,u)

sys = [];

% end mdlTerminate
