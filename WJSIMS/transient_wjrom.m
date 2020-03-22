clc
clear all
addpath('../ROUTINES/ROUTINES')
addpath('../ROUTINES/ROUTINES/FEM')
addpath('../ROUTINES/ROUTINES/QUASISTATIC')
addpath('../ROUTINES/ROUTINES/TRANSIENT')
addpath('../ROUTINES/ROUTINES/SOLVERS')
addpath('../ROUTINES/ROUTINES/CONTACTMODELS')

% ThreshRed
Prestress = 11580; 
sint = 1e-6; 

sel_method = 'U' % 'P'[2,3,4,5,6,7,11] 'PD' [2,3,4,5,6,7,10] 'U'
                 % [3,5,7,9,11,19,21,30,40,45,50,60,70,80,90]
Nlev = 3;
fname = sprintf('../WHOLEJOINT_ROM_PREPARE/REDMATS/P%d_S%.2f_%s_%dLEV_GRED_WJMAT_prop.mat',Prestress, log10(sint), sel_method, Nlev);
load(fname, 'M', 'K', 'L', 'R', 'Fv', 'Th', 'I2', 'cnum', 'MESH', ...
     'PatchAreas', 'Npatches', 'NTG', 'GTG', 'dofred', 'GTG', 'Pels')

mi = 1;
%% Setup GMDOF class
SYS = GMDOF((1:Npatches)', 6, L(1:Npatches*6, :), L(1:Npatches*6, :)');
SYS = SYS.SETCFUN(@(us, z, uds, P) ELDRYFRICT_WJ(us, z, uds, P, I2, 0), sparse(2, Npatches));

%% Interface Parameters
%% Interface Parameters
mu = 0.20;  %% Friction coefficient
c = 1.0;  %% Standard Deviation (microns) of surface roughness

nu   = 0.29;                % Poisson's ratio
Aint = sum(PatchAreas);     % Total interface area (lumped)
Pint = Prestress*3/Aint;    % Average lumped interface pressure
sint = 1e-6*c;              % Average interface asperity SD (exponential distribution)
chi  = 2.0;                 % Mindlin constant
ktkn = chi*(1-nu)/(2-nu);   % Tangential-normal stiffness ratio
kt   = 4*(1-nu)*Pint/(sqrt(pi)*(2-nu)*sint);    % Tangential stiffness
kn   = kt/ktkn;             % Normal Stiffness

Pars = [kt; kt; kn; mu];
pA = zeros(Npatches*4, 4);
pA(1:4:end, 1) = PatchAreas;
pA(2:4:end, 2) = PatchAreas;
pA(3:4:end, 3) = PatchAreas;
pA(4:4:end, 4) = 1.0;

%% Prestress Analysis

				% Fully Stuck Initial Guess
Kst = zeros(Npatches*6);
Kst(1:6:end, 1:6:end) = diag(PatchAreas*kt);
Kst(2:6:end, 2:6:end) = diag(PatchAreas*kt);
Kst(3:6:end, 3:6:end) = diag(PatchAreas*kn);
Kst(4:6:end, 4:6:end) = diag(I2(1,:)*kt);
Kst(4:6:end, 5:6:end) = diag(I2(2,:)*kt);
Kst(5:6:end, 4:6:end) = diag(I2(2,:)*kt);
Kst(5:6:end, 5:6:end) = diag(I2(3,:)*kt);
Kst(6:6:end, 6:6:end) = diag(I2(4,:)*kn);
K0 = K + SYS.Tm*Kst*SYS.Qm;

U0 = K0\(Prestress*Fv);
Z0 = SYS.z;
%% Nonlinear Prestress Simulation
disp('STATIC PRESTRESS ANALYSIS')
opts = struct('reletol', 1e-6, 'rtol', 1e-6, 'utol', 1e-6, 'etol', ...
              1e-6, 'ITMAX', 100, 'Display', true);
[Ustat, ~, eflag, ~, J0] = NSOLVE(@(U) QS_RESFUN([U; 0], Z0, Pars, speye(size(K)), pA, ...
						 SYS, M, K, Fv*Prestress, Fv*0), U0, opts);
[~, Z, ~, ~, ~, Txyn_p] = SYS.CONTACTEVAL(Ustat, Z0, Ustat*0, Pars, pA);  % Evaluate Slider
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
Z0 = full(Z);
disp('INITIAL CONDITIONS SET');

%% HHTA Hysteretic
ABG = [0, 1/4, 1/2];  % Unconditionally Stable Newmark-Alpha
% ABG = [0, 1/6, 1/2];  % Implicit linear acceleration
% ABG = [-0.1, 1/6, 1/2];  % HHT-Alpha

T0 = 0;
T1 = 1e0;
dT = 1e-5;  % 5000 Hz Nyquist

opts = struct('reletol', 1e-12, 'etol', 1e-6, 'rtol', 1e-6, 'utol', 1e-6, ...
	      'Display', false, 'ITMAX', 100, 'waitbar', true);
tic
[Th, Xh, zh, Xdh, Xddh] = HHTA_NONLIN_HYST(M, C, K, fex, ...
					   @(t, x, z, xd) SYS.CONTACTEVAL(x, z, xd, Pars, pA), ...
					   U0, Z0, Ud0, T0, T1, dT, ABG(1), ABG(2), ABG(3), opts);
ttk = toc

figure(1)
clf()
plot(Th, R(3, :)*Xddh, '.-', 'LineWidth', 1);
Fh = fex(Th);
Finput = fdyn(Th);

save(sprintf('./DATS/TRANSIENT_%s_%dLEV.mat',sel_method,Nlev), 'Th', 'Xh', 'zh', 'Xdh', 'Xddh', 'Fh', 'Finput', 'famp', 'freq', 'ttk')
