clc
clear all
addpath('../ROUTINES/ROUTINES')
addpath('../ROUTINES/ROUTINES/FEM')
addpath('../ROUTINES/ROUTINES/CONTACTMODELS')
addpath('../ROUTINES/ROUTINES/QUASISTATIC')
addpath('../ROUTINES/ROUTINES/TRANSIENT')
addpath('../ROUTINES/ROUTINES/SOLVERS')

fname = './PREP/BRB_MATS_NR.mat';
load(fname,'M','K','R','L','Fv','TFM');
Prestress = 11580;

%% Create Mesh Structure
Nds = dlmread('./PREP/Nodes.dat');
Quad = dlmread('./PREP/Elements.dat');
MESH = MESH2D(Nds, 3, [], Quad, 2);
MESH = MESH.SETCFUN(@(u, z, ud, P) ELDRYFRICT(u, z, ud, P, 0), sparse(2, MESH.Ne*MESH.Nq^2));  % Contact Function

%% Interface Parameters
mu = 0.20;  %% Friction coefficient
c = 1.0;  %% Standard Deviation (microns) of surface roughness

nu   = 0.29;                % Poisson's ratio
Aint = sum(sum(MESH.Tm));      % Average interface area (lumped)
Pint = Prestress*3/Aint;           % Average lumped interface pressure
sint = 1e-6*c;              % Average interface asperity SD (exponential distribution)
chi  = 2.0;                 % Mindlin constant
ktkn = chi*(1-nu)/(2-nu);   % Tangential-normal stiffness ratio
kt   = 4*(1-nu)*Pint/(sqrt(pi)*(2-nu)*sint);    % Tangential stiffness
kn   = kt/ktkn;             % Normal Stiffness

Pars = [kt; kt; kn; mu];
pA = repmat(eye(4), MESH.Ne*MESH.Nq^2, 1);
%% Prestress Analysis

% Fully Stuck Initial Guess 
Kstuck = zeros(size(L, 1));
Kstuck(1:3:MESH.Nn*3, 1:3:MESH.Nn*3) = MESH.Tm*MESH.Qm*Pars(1);
Kstuck(2:3:MESH.Nn*3, 2:3:MESH.Nn*3) = MESH.Tm*MESH.Qm*Pars(2);
Kstuck(3:3:MESH.Nn*3, 3:3:MESH.Nn*3) = MESH.Tm*MESH.Qm*Pars(3);

Kstuck = K + L'*Kstuck*L;

U0 = Kstuck\(Fv*Prestress);
Z0 = zeros(2, MESH.Ne*MESH.Nq^2);
%% Nonlinear Prestress Simulation
disp('STATIC PRESTRESS ANALYSIS')
opts = struct('reletol', 1e-6, 'rtol', 1e-6, 'utol', 1e-6, 'etol', ...
              1e-6, 'ITMAX', 100, 'Display', true);
[Ustat, ~, eflag, ~, J0] = NSOLVE(@(U) QS_RESFUN([U; 0], Z0, Pars, L, pA, ...
						 MESH, M, K, Fv*Prestress, Fv*0), U0, opts);
[~, Z, ~, ~, ~, Txyn_p] = MESH.CONTACTEVAL(Ustat, Z0, Ustat*0, Pars, pA, L);  % Evaluate Slider

%% Linearized Modal Analysis
[V, D] = eigs(J0, M, 10, 'SM');
[D,si] = sort(diag(D));
Ws = sqrt(D)/2/pi;
V = V(:,si);
V = V./sqrt(diag(V'*M*V))';

% SET LINEAR DAMPING
zts = [0.002; 0.003];
ab = [1./(2*2*pi*Ws([1 3])) 2*pi*Ws([1 3])/2]\zts;
C = ab(1)*M+ab(2)*J0;

% SET DYNAMIC EXCITATION
freq = 500;
famp = 800;
fdyn = @(t) famp*sin(2*pi*freq*t).^2.*(t<0.5/freq);
fex = @(t) R(3, :)'*(fdyn(t))+Fv*Prestress;

% Initial Conditions
U0 = Ustat;
Ud0 = Ustat*0;
Z0 = Z;
disp('INITIAL CONDITIONS SET');

%% HHTA Hysteretic
ABG = [0, 1/4, 1/2];  % Unconditionally Stable Newmark-Alpha
% ABG = [0, 1/6, 1/2];  % Implicit linear acceleration
% ABG = [-0.1, 1/6, 1/2];  % HHT-Alpha

T0 = 0;
T1 = 2.5;
dT = 1e-4;  % 5000 Hz Nyquist

opts = struct('reletol', 1e-12, 'etol', 1e-6, 'rtol', 1e-6, 'utol', 1e-6, ...
    'Display', true, 'ITMAX', 100, 'waitbar', true);
tic
[Th, Xh, zh, Xdh, Xddh] = HHTA_NONLIN_HYST(M, C, K, fex, ...
					   @(t, x, z, xd) MESH.CONTACTEVAL(x, z, xd, Pars, pA, L), ...
					   U0, Z0, Ud0, T0, T1, dT, ABG(1), ABG(2), ABG(3), opts);
ttk = toc
clf()
plot(Th, R(3, :)*Xddh, '.-', 'LineWidth', 2)
Fh = fex(Th);
Finput = fdyn(Th);

save('./DATS/TRANSIENT_REF.mat', 'Th', 'Xh', 'zh', 'Xdh', 'Xddh', 'Fh', 'Finput', 'famp', 'freq', 'ttk')
