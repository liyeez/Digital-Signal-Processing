
function [V_BMtd, SS]= BM_passive(input,fs)

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:

% input: input/excitation signal to BM
%        1 x T vector

% fs:    sampling frequency of the signal
%        1 scalar


% Output:

% V_BMtd: BM displacement
%         T X N matrix
%         Row corresponding to time step
%         Column corresponding to each element of BM

% SS: structure contains model parameters
%

%% Number of elements in BM in the simulation
N = 1000;

%%
global SS
SS.U = input; %input signal
SS.fs = fs; % sample frequency
SS.tstep = 1/SS.fs;% time interval
SS.t = [0:length(input)-1]/SS.fs; % time for ode solver
SS.NX = N; % Number of positions along the cochlea

%%  Set simulation parameters
SS.L=35e-3;                    % Length of the cochlea [m]
SS.X = 0:SS.L/(SS.NX-1):SS.L;  % Positions along the cochlea [m]
SS.Delta = SS.X(2)-SS.X(1);    % Spatial increament [m]
SS.Q=5;                        % Q factor
SS.rho = 1000;                 % Density of water
SS.W=1e-3;                     % Cochlear partition width [m]
SS.H=1e-3;                     % Cochlear chamber height  [m]
SS.B=0.3e-3;                   % Basilar membrane width [m]
SS.m=0.28;                     % Assumed BM mass per area [kg/m^2]
SS.OmegaC=20000*exp(-SS.X/7e-3)*2*pi; % Characteristics frequency [rad]
SS.k=(SS.OmegaC.^2).*SS.m;     % Assuemd BM stiffness [Nm^-3]
SS.c=sqrt(SS.k.*SS.m)/SS.Q;    % Assuemd BM damping   [Nsm^-3]
SS.h= pi^2*SS.W*SS.H/(8*SS.B); % Effective height

%% Forming A, B, C and F (fluid coupling) matrices
lA = 2*SS.NX; % length of A
Num_sts = 2;      % number of micromechanical states

%%%%%%%%%%%%%%%%%%%%  Selecting 2nd output of displaceent %%%%%%%%%%
% SS.C = zeros(SS.NX, lA);  % Initialize C
% SS.C(1,1) = 1;         % Middle Ear
% SS.C(SS.NX*2+2 : SS.NX*Num_sts+1 : SS.NX*2+(SS.NX*Num_sts+1)*(SS.NX-2)) = 0.4; % cochlea
% SS.C(end,end-1) = 1;   % Helicotrema

SS.C = zeros(SS.NX, lA);  % Initialize C
SS.C(1,2) = 1;         % Middle Ear
SS.C(SS.NX*3+2 : SS.NX*Num_sts+1 : SS.NX*3+(SS.NX*Num_sts+1)*(SS.NX-2)) = 0.4; % cochlea
SS.C(end,end) = 1;   % Helicotrema

%%%%%%%%%%%%%%%%%%%%  Constructing fluid coupling matrix F %%%%%%%%%%
F_top = [-1 1 zeros(1,SS.NX-2)];
F_mid = [1*diag(ones(SS.NX-2,1),0) zeros(SS.NX-2,2)] ...
    + [zeros(SS.NX-2,1) diag(-2*ones(SS.NX-2,1),0) zeros(SS.NX-2,1)]...
    + [zeros(SS.NX-2,2) 1*diag(ones(SS.NX-2,1),0)];
F_bot = [zeros(1,SS.NX-2) SS.Delta/SS.h -(SS.Delta/SS.h + SS.Delta^2/SS.h^2)];
SS.F = [1/(2*SS.rho*SS.Delta)*F_top ; SS.h/(2*SS.rho*SS.Delta^2)*F_mid ;SS.h/(2*SS.rho*SS.Delta^2)*F_bot];
SS.iF = inv(SS.F);

%  Initialize isolated A and B matrices
AE = zeros(2*SS.NX);       % A matrix
BE = zeros(2*SS.NX,SS.NX); % B matrix

As = 3.2*1e-6;                                     % stapes area
% Middle ear properties
m_ME = 4.4e5*As;        % Mass
c_ME = 10e9*As;         % Damping
k_ME = 8.1e13*As;       % Stiffness

% Helicotrema
CH = SS.c(end);        % Damping
MH = SS.m;             % Mass

% Define AE and BE matrices
AE(1:2,1:2) = [-c_ME/m_ME -k_ME/m_ME; 1 0];% boundry conditions applied

A11 = -SS.c./SS.m;  % Vectorise elements of A for speed
A12 = -SS.k./SS.m;
A21 = ones(1,SS.NX);

AE(2*lA+3:lA*Num_sts+Num_sts:2*lA+3+(lA*Num_sts+Num_sts)*(SS.NX-3)) = A11(2:end-1);
AE(3*lA+3:lA*Num_sts+Num_sts:3*lA+3+(lA*Num_sts+Num_sts)*(SS.NX-3)) = A12(2:end-1);
AE(2*lA+4:lA*Num_sts+Num_sts:2*lA+4+(lA*Num_sts+Num_sts)*(SS.NX-3)) = A21(2:end-1);
AE(end-1:end,end-1:end) = [-(CH/MH) , 0; ...
                             1      , 0      ];

BE(1:2,1) = [1/m_ME 0];
BE(lA+3:lA+Num_sts:lA+3+(lA+Num_sts)*(SS.NX-3)) = 1./SS.m;
BE(end-1:end,end) = [1/(MH) ; 0];

%   Coupled A and B matrices
SS.A = (eye(size(AE)) - BE/SS.F*SS.C)\AE;
SS.B= (eye(size(AE)) - BE/SS.F*SS.C)\BE;

%% Defining simulation
% Initial condition
SS.X0 = zeros(1,2*SS.NX);

%%%%%%%%%%%%% Solve differential equations using ode45%%%%%%%%
SS.AbsTol = 11;
SS.AbsTolVect = ones(1,2*SS.NX)*SS.AbsTol;
SS.AbsTolVect(1:2:end) = SS.AbsTol-4;
SS.RelTol = 4;
SS.options = odeset('RelTol',1*10^(-SS.RelTol),'AbsTol',1*10.^(-SS.AbsTolVect));

% Call ode45 function
[SS.T,SS.BM] = ode45(@ss_formulation,SS.t,SS.X0,SS.options);
V_BMtd = SS.BM(:,1:2:end);  % Save only BM displacement

function xdot = ss_formulation(t,x)
global SS
input_index = floor(t./SS.tstep)+1;
xdot = SS.A*x + SS.B*[SS.U(input_index); zeros(SS.NX-1,1)];
