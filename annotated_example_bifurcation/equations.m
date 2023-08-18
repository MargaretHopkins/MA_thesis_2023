function [Rf,dRf] = equations(sys,Uf,dUf)
% The auxiliary variables are defined first and then the residu of the
% system considered.

%% Variables :  Uf=[ U  ; Ua ]   with U=[ u ; lambda ]    - no need to change

u      = Uf(1:sys.neq);      % Main variables
lambda = Uf(sys.neq+1);      % Continuation parameter
Ua     = Uf(sys.neq+2:end);  % Auxiliary variables

du     = dUf(1:sys.neq);     % differential of Main variables
dlambda= dUf(sys.neq+1);     % differential of the continuation parameter
dUa    = dUf(sys.neq+2:end); % differential of Auxiliary variables
dRf=zeros(sys.neq_aux, 1);

%% Parameters of the system
% paramters either defined in launch file or derived from those defined in
% launch file

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
%s8 = sys.parameters.s8;
%C8 = sys.parameters.C8;
P_M = mu_l*omega_l^2*y0;
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
%C8t = C8/imag(s1);
%s8t = s8/imag(s1);

%% Residues
R     = zeros(sys.neq,1);      % Main residue
Ra    = zeros(sys.neq_aux,1);  % Auxiliary residue

dR    = zeros(sys.neq,1);     % Differential form of non-quadratic part of the main residue
dRa   = zeros(sys.neq_aux,1);    % Differential form of non-quadratic part of the auxiliary residue


% Equations of the main system with all time derivatives set to zero
%Ua contains the auxiliary variables and u contains the main variables, in
%order
R(1) =  real(C1t)*Ua(5) +real(s1t)* u(1) -imag(s1t)*u(2);                              
R(2) =  imag(C1t)*Ua(5)+imag(s1t)*u(1)+real(s1t)*u(2); 
R(3) =  real(C2t)*Ua(5) +real(s2t)* u(3) -imag(s2t)*u(4);                                
R(4) =  imag(C2t)*Ua(5)+imag(s2t)*u(3)+real(s2t)*u(4);
R(5) =  real(C3t)*Ua(5) +real(s3t)* u(5) -imag(s3t)*u(6);                               
R(6) =  imag(C3t)*Ua(5)+imag(s3t)*u(5)+real(s3t)*u(6);
R(7) =  real(C4t)*Ua(5) +real(s4t)* u(7) -imag(s4t)*u(8);                               
R(8) =  imag(C4t)*Ua(5)+imag(s4t)*u(7)+real(s4t)*u(8);
R(9) =  real(C5t)*Ua(5) +real(s5t)* u(9) -imag(s5t)*u(10);                               
R(10) =  imag(C5t)*Ua(5)+imag(s5t)*u(9)+real(s5t)*u(10);
R(11) =  real(C6t)*Ua(5) +real(s6t)* u(11) -imag(s6t)*u(12);                               
R(12) =  imag(C6t)*Ua(5)+imag(s6t)*u(11)+real(s6t)*u(12);
R(13) =  real(C7t)*Ua(5) +real(s7t)* u(13) -imag(s7t)*u(14);                               
R(14) =  imag(C7t)*Ua(5)+imag(s7t)*u(13)+real(s7t)*u(14);
%R(15) =  real(C8t)*Ua(5) +real(s8t)* u(15) -imag(s8t)*u(16);                               
%R(16) =  imag(C8t)*Ua(5)+imag(s8t)*u(15)+real(s8t)*u(16);
R(15) = omega_lt*u(16);                                                                      
R(16) = omega_lt*(1- u(15) - u(16)/Q_l+lambda-Ua(1));  
%R(3) = omega_lt*u(4)-du(3);
%R(4) = omega_lt*(1 - u(3) - u(4)/Q_l+lambda-Ua(1))-du(4);

...
    
%U0 = [R1 I1 x z lambda p s v w u]';
%U0 = [R1 I1 R2 I2 R3 I3 R4 I4 R5 I5 R6 I6 R7 I7 R8 I8 x z lambda p s v w u]'; 

% Definition of the auxiliary variables 
Ra(1) =  2*(sum(u(1:2:14)))-Ua(1);          
%Ra(1) = 2*u(1)-Ua(1);
Ra(2) =  u(15)^2 + 1e-3 - Ua(2)^2;  
%Ra(2) = u(3)^2+1e-3-Ua(2)^2;
Ra(3) = lambda - Ua(1) -Ua(3)*Ua(4);                    
Ra(4) = Ua(3)^2+1E-3-Ua(4)^2;                          
Ra(5) = xi*Ua(3)*(Ua(2)+u(15))/2-Ua(5); 
%Ra(5) = xi*Ua(3)*(Ua(2)+u(3))/2-Ua(5);
...
    

% Concatenation of the two residues
Rf =[R  ; Ra ];
%dRf=[dR ; dRa];


end
