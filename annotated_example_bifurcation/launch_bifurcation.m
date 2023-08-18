%example file used for Margaret Hopkins' thesis
%June 17, 2023
%see thesis for a description of the system and its auxiliary variables

global U Section Diagram   % Global variables to export point from the diagram. 

% Path of the Manlab SRC file - edit this to where you have the Manlab SRC
addpath('/your-filepath/MANLAB-4-1-7/MANLAB-4.1.7/SRC')

%% parameters of the system
neq = 16;     % Number of main equations
neq_aux = 5; % Number of auxiliary variables used

% Parameters specific to your system
parameters.rho = 1.2041; %density of air at 20C
parameters.c = 378.11; %speed of sound in air at 20C
parameters.Q_l = 3.2;%lip quality factor
parameters.mu_l = 2; %kg/m^2 lip mass per unit surface area
parameters.y0 = 0.1E-3; %m lip position at rest 
parameters.b = 8E-3; %m lip opening width
parameters.d = 17.23E-3; %diameter of 1.5c Bach mouthpiece in m
parameters.omega_l = 385.84*2*pi; %natural lip frequency for Alberto's trumpet, obtained through LSA
parameters.A = (parameters.d/2)^2*pi; %cross-sectional area of mouthpiece opening
parameters.Z_c = parameters.rho*parameters.c/parameters.A; %characteristic impedance based on size of mouthpiece opening
%variables from making the system dimensionless
parameters.P_M = parameters.mu_l*parameters.omega_l^2*parameters.y0;
parameters.xi = parameters.Z_c*parameters.b*parameters.y0*sqrt(2/(parameters.rho*parameters.P_M));
%modal impedance parameters from a fitting scheme, for Yamaha Xeno trumpet
parameters.sn = 1.0e+03*[-0.0251+1.4780i -0.0332+2.1945i -0.0390+2.9181i -0.0442+3.6797i -0.0491+4.3459i -0.0589+5.0555i -0.0675+5.7585i];
parameters.Cn = 1.0e+03*[1.0462 1.5824 2.9424 3.6585 4.5920 4.3775 3.2783];

%% initialization of the system - Should not be changed
sys=SystAQ(neq,neq_aux,@equations,@point_display,@global_display,parameters);

%% starting point
%%% Main variables u and continuation parameter lambda
%initialize according to the state of your system at equilibrium
%the complex pressure components are set to zero, as the system is not
%moving
%the number of real and imaginary componenets is equal to the number of
%modes fit, in this case 7
R1 = 0;
I1 = 0;
R2 = 0;
I2 = 0;
R3 = 0;
I3 = 0;
R4 = 0;
I4 = 0;
R5 = 0;
I5 = 0;
R6 = 0;
I6 = 0;
R7 = 0;
I7 = 0;
x = 1; %lip starts at resting position, and x=y/y0
z = 0;
lambda = 0.001; %continuation parameter - set slightly above zero to avoid numerical issues

%%% Auxiliary variables v :
%initialize according to how they relate to the main variables
p = 2*(R1+R2+R3+R4+R5+R6+R7);
s = sqrt(x^2+1E-3);
v = sqrt(lambda-p);
w = sqrt(v^2+1E-3);
u = parameters.xi*(s+x)*v/2;

%%% Vector containing all the unknowns in order
U0 = [R1 I1 R2 I2 R3 I3 R4 I4 R5 I5 R6 I6 R7 I7 x z lambda p s v w u]'; 
%U0 = [R1 I1 x z lambda p s v w u]';

%% Launch Manlab
Manlab('sys'               ,sys, ...        % description of your system
         'U0value'         ,U0, ...         % starting point
         'NRthreshold'     ,1e-10, ...
         'ANMthreshold'    ,1e-14, ...
         'NRmethod'        ,2, ...
         'PointDisplay'    ,0, ...      % Point display [on]/off
         'GlobalDisplay'   ,0, ...      % Global display [on]/off
         'NRstart'         ,1, ...
         'NRitemax'        ,20, ...
         'StabilityCheck'  ,1, ...          %stability computation on
         'Amax_max'        ,1, ...          %max radius of convergence of the series
         'StabTol'         ,1E-8, ...
         'displayvariables',[neq+1 neq+2]);     % MANLAB run



% With 'displayvariables',[neq+1 1] ; the projected bifurcation diagram
% (figure 2) will represents the first variable with respect to the
% continuation parameter lambda.

%% Optional arguments of Manlab launching function with their default values :
%          'order'           ,20, ...     % order of the series
%          'ANMthreshold'    ,1e-6, ...   % threshold for the domain of validity of the series
%          'Amax_max'        ,1e6, ...    % maximum value of the domain of validity of the series
%          'NRthreshold'     ,2e-5, ...   % threshold for Newton-Raphson (NR) corrections
%          'NRitemax'        ,10, ...     % Maximum number of iteration of NR algorithm
%          'NRstart'         ,1, ...      % NR corrections for the user-defined starting point [on]/off
%          'NRmethod'        ,0, ...      % NR corrections on/[off]
%          'BifDetection'    ,1, ...      % Detection of bifurcation [on]/off
%          'PointDisplay'    ,1, ...      % Point display [on]/off
%          'GlobalDisplay'   ,1, ...      % Global display [on]/off
%          'StabilityCheck'  ,0, ...      % Stability computation on/[off]
%          'StabTol'         ,1e-6, ...   % Stability tolerance
%
%               The stability requires to write the system of equations in
%               the explicit ODE form X' = f(X). 
%               The residue function is then R(X) = f(X) = 0.
%               In all other cases, it gives wrong results.
