function [Rf,dRa,Forcing] = equations(sys,t,zf,dzf,d2zf)
% Equations of the system of the form 
% Rf(zf) = C + L0(zf) + lambda L1(zf) + D0(dzf) + lambda D1(dzf) + DD(d2zf) + Q(zf,zf) + f(zf).

%%% parameters of the system
rho = sys.parameters.rho; %density of air at 20C
c = sys.parameters.c; %speed of sound in air at 20C
Q_l = sys.parameters.Q_l;%lip quality factor
mu_l = sys.parameters.mu_l; %kg/m^2 lip mass per unit surface area
y0 = sys.parameters.y0; %m lip position at rest 
b = sys.parameters.b; %m lip opening width
d = sys.parameters.d; %diameter of 1.5c Bach mouthpiece in m
omega_l = sys.parameters.omega_l; %natural lip frequency ***should actually be obtained through LSA
A = sys.parameters.A; %cross-sectional area of mouthpiece opening
s1 = sys.parameters.sn(1); %first pole
C1 = sys.parameters.Cn(1); %first residue
s2 = sys.parameters.sn(2);
C2 = sys.parameters.Cn(2);
s3 = sys.parameters.sn(3);
C3 = sys.parameters.Cn(3);
s4 = sys.parameters.sn(4);
C4 = sys.parameters.Cn(4);
s5 = sys.parameters.sn(5);
C5 = sys.parameters.Cn(5);
s6 = sys.parameters.sn(6);
C6 = sys.parameters.Cn(6);
s7 = sys.parameters.sn(7);
C7 = sys.parameters.Cn(7);
Z_c = sys.parameters.Z_c; %characteristic impedance based on size of mouthpiece opening
%variables from making the system dimensionless
P_M = sys.parameters.P_M;
xi = sys.parameters.xi;
omega_lt = omega_l/imag(s1);
C1t = C1/imag(s1);
s1t = s1/imag(s1);
C2t = C2/imag(s1);
s2t = s2/imag(s1);
C3t = C3/imag(s1);
s3t = s3/imag(s1);
C4t = C4/imag(s1);
s4t = s4/imag(s1);
C5t = C5/imag(s1);
s5t = s5/imag(s1);
C6t = C6/imag(s1);
s6t = s6/imag(s1);
C7t = C7/imag(s1);
s7t = s7/imag(s1);

%% Variables : 
%u      = zf(1:sys.nz);         % Main variables
Rn = zf(1:2:14-1); %real part of modal pressure
In = zf(2:2:14);
x = zf(14+1);
z = zf(14+2);
p = zf(2*7+3);
s = zf(2*7+4);
v = zf(2*7+5);
w = zf(2*7+6);
u = zf(2*7+7);
lambda = zf(14+8);

%Ua     = zf(sys.nz+1:end-1);   % Auxiliary variables
%lambda = zf(end);              % Continuation parameter

dRn = dzf(1:2:14-1);
dIn = dzf(2:2:2*7);
dx = dzf(2*7+1);
dz = dzf(2*7+2);
%du     = dzf(1:sys.nz);        % first order derivated of Main variables
%dUa    = dzf(sys.nz+1:end);    % first order derivated of Auxiliary variables

%d2u    = d2zf(1:sys.nz);       % second order derivated of Main variables
%d2Ua   = d2zf(sys.nz+1:end);   % first order derivated of Auxiliary variables

%% Residues
R     = zeros(sys.nz,1);         % Main residue
Ra  = zeros(sys.nz_aux,1);     % Auxiliary residue

dRa = zeros(sys.nz_aux,1);     % Differential form of non-quadratic part of the auxiliary residue

% Equations of the main system
R(1) =  real(C1t)*u +real(s1t)* Rn(1) -imag(s1t)*In(1) - dRn(1);                                
R(2) =  imag(C1t)*u+imag(s1t)*Rn(1)+real(s1t)*In(1) - dIn(1); 
R(3) =  real(C2t)*u +real(s2t)* Rn(2) -imag(s2t)*In(2) - dRn(2);                                
R(4) =  imag(C2t)*u+imag(s2t)*Rn(2)+real(s2t)*In(2) - dIn(2);
R(5) =  real(C3t)*u +real(s3t)* Rn(3) -imag(s3t)*In(3) - dRn(3);                               
R(6) =  imag(C3t)*u+imag(s3t)*Rn(3)+real(s3t)*In(3) - dIn(3);
R(7) =  real(C4t)*u +real(s4t)* Rn(4) -imag(s4t)*In(4) - dRn(4);                               
R(8) =  imag(C4t)*u+imag(s4t)*Rn(4)+real(s4t)*In(4) - dIn(4);
R(9) =  real(C5t)*u +real(s5t)* Rn(5) -imag(s5t)*In(5) - dRn(5);                               
R(10) =  imag(C5t)*u+imag(s5t)*Rn(5)+real(s5t)*In(5) - dIn(5);
R(11) =  real(C6t)*u +real(s6t)* Rn(6) -imag(s6t)*In(6) - dRn(6);                               
R(12) =  imag(C6t)*u+imag(s6t)*Rn(6)+real(s6t)*In(6) - dIn(6);
R(13) =  real(C7t)*u +real(s7t)* Rn(7) -imag(s7t)*In(7) - dRn(7);                               
R(14) =  imag(C7t)*u+imag(s7t)*Rn(7)+real(s7t)*In(7) - dIn(7);
%R(15) =  real(C8t)*u +real(s8t)* Rn(8) -imag(s8t)*In(8) - dRn(8);                               
%R(16) =  imag(C8t)*u+imag(s8t)*Rn(8)+real(s8t)*In(8) - dIn(8);
R(15) = omega_lt*z-dx;                                                                      
R(16) = omega_lt*(1- x - z/Q_l+lambda-p)-dz;  


% Definition of the auxiliary variables | differentiation of the non-quadratic equations
Ra(1) =  2*sum(Rn(1:7))-p;                          
Ra(2) = x^2 +1E-3 - s^2;
Ra(3) = lambda - p - v*w;
Ra(4) = v^2 + 1E-3 -w^2;
Ra(5) = xi*v*(s+x)/2-u;
    
% Concatenation of the two residues
Rf=[R ; Ra];

%% Forcing terms
% Should be written as if the forcing angular frequency value is 1
% i.e. the forcing period is 2*pi
Forcing = zeros(2*sys.H+1,sys.nz_tot); % DO NOT CHANGE this line.

% if the equation number k is forced, write the forcing in Forced(:,k) as :
%Forcing(:,k) = forcing_function(t); 

end
