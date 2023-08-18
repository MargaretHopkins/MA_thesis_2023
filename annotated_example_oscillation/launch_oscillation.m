% Empty example.
% Write here the description of your system, the auxiliary variables used,
% and its final quadratic recast

global U Section Diagram   % Global variables to export point from the diagram.

% Path of the SRC file.
addpath('/your-filepath/MANLAB-4-1-7/MANLAB-4.1.7/SRC')

%% Parameters of the system
nz= 16;               % number of main equations of the system of the differential-algebraic system (DAE)
nz_aux =5;           % number of auxiliary equations of the system of the DAE

H = 20;              % number of harmonics used to compute the solution-branch - 20 is a good start

%%% Parameters specific to your system
parameters.rho = 1.2041; %density of air at 20C
parameters.c = 343.21; %speed of sound in air at 20C
parameters.Q_l = 3.2;%lip quality factor
parameters.mu_l = 2; %kg/m^2 lip mass per unit surface area
parameters.y0 = 0.1E-3; %m lip position at rest 
parameters.b = 8E-3; %m lip opening width
parameters.d = 17.23E-3; %diameter of 1.5c Bach mouthpiece in m
parameters.omega_l = 385.84*2*pi; %natural lip frequency ***should actually be obtained through LSA
%parameters.d = 16.85E-3; %diameter of 5c Bach mouthpiece in m
parameters.A = (parameters.d/2)^2*pi; %cross-sectional area of mouthpiece opening
parameters.sn = 1.0e+03*[-0.0245+1.4837i -0.0311 + 2.2032i -0.0376+2.9253i -0.0438+3.6883i -0.0507+4.3520i -0.0616+5.0600i -0.0720+5.7560i];
parameters.Cn = 1.0e+03*[1.0462 1.5824 2.9424 3.6585 4.5920 4.3775 3.2783];
parameters.Z_c = parameters.rho*parameters.c/parameters.A; %characteristic impedance based on size of mouthpiece opening
%variables from making the system dimensionless
parameters.P_M = parameters.mu_l*parameters.omega_l^2*parameters.y0;
parameters.xi = parameters.Z_c*parameters.b*parameters.y0*sqrt(2/(parameters.rho*parameters.P_M));


%% initialization of the system
type = 'autonomous';        % type of system (can be 'forced' or 'autonomous')
writing = 'standard';   % way the equations have been written (can be 'standard' or 'vectorial')

sys=SystODE(nz,nz_aux,H,@equations,@point_display,@global_display,parameters,type);


%load section from wherever you saved the exported section in the previous
%step
%vector U0 contains all the unknowns
load('/your-directory/Section.mat')
U0 = sys.init_Hopf(Section); % Automatic initialization from a Hopf point.

%ML_pointdisplay=1;

%% Launch of Manlab with options
Manlab('sys'       ,sys , ...
    'U0value'         ,U0, ...
    'NRthreshold'     ,1e-10, ...
    'ANMthreshold'    ,1e-14, ...
    'NRmethod'        ,0, ...      % important to keep correction off to start with for the oscillating solution
    'PointDisplay'    ,0, ...      % Point display [on]/off
    'GlobalDisplay'   ,0, ...      % Global display [on]/off
    'NRstart'         ,0, ...
    'NRitemax'        ,20, ...
    'StabilityCheck'  ,0, ...
    'Amax_max'        ,1, ...
    'StabTol'         ,1e-8, ...
    'displayvariables',[sys.neq+1 2]) %displaying the first Fourier coefficient in front of first cosine of Fourier series vs. continuation parameters lambda - need a separate file later to plot p
%    'ML_pointdisplay',  1 )
    % MANLAB run




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
