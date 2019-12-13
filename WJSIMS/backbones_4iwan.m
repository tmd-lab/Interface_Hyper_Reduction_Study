clc
clear all
addpath('../ROUTINES/')
addpath('../ROUTINES/FEM/')
addpath('../ROUTINES/IWAN/')

% 5 Patches
setid = 5;
fname = sprintf('../WHOLEJOINT_ROM_PREPARE/REDMATS/%d_SET_WJMAT.mat', setid);

% ThreshRed
Prestress = 11580; sint = 5e-6; Nlev=8;
fname = sprintf('../WHOLEJOINT_ROM_PREPARE/REDMATS/P%d_S%.2f_%dLEV_GRED_WJMAT.mat',Prestress, log10(sint), Nlev);

load(fname, 'M', 'K', 'L', 'R', 'Fv', 'Th', 'LamT', 'cnum', 'MESH', 'PatchAreas', 'Npatches', 'dofred', 'GTG')
Prestress = 11580;
mi = 1;
%% Access Matrices and functions
% LamT = LamT*GTG;
Qxyn = L(reshape(((1:Npatches)-1)*6+(1:3)', Npatches*3, 1), :);
Txyn = LamT(:, reshape(((1:Npatches)-1)*6+(1:3)', Npatches*3, 1));

copt.lspci = [1 2 4 5];  % [Fs Kt Chi Bt Kn]  (Chi not log-scale)

copt.x.T{1} = [PatchAreas' zeros(Npatches, 4)];  % Fs
copt.x.T{2} = [zeros(Npatches,1) PatchAreas' zeros(Npatches, 3)];  % Kt
copt.x.T{3} = [zeros(Npatches,2) ones(Npatches,1) zeros(Npatches, 2)];  % Chi
copt.x.T{4} = [zeros(Npatches,3) ones(Npatches,1) zeros(Npatches, 1)];  % Bt

copt.y.T{1} = [PatchAreas' zeros(Npatches, 4)];  % Fs
copt.y.T{2} = [zeros(Npatches,1) PatchAreas' zeros(Npatches, 3)];  % Kt
copt.y.T{3} = [zeros(Npatches,2) ones(Npatches,1) zeros(Npatches, 2)];  % Chi
copt.y.T{4} = [zeros(Npatches,3) ones(Npatches,1) zeros(Npatches, 1)];  % Bt

copt.n.T{1} = [zeros(Npatches, 4) PatchAreas'];  % Fs

CFUN = @(uxyn, pars) SPRING_IWAN4(uxyn, pars, copt);
DFUN = @(uxyn, pars) PENALTY_IWAN4_DISS(uxyn, pars, copt);

%% Contact Model
pars = [9; 20; -0.8; -4; 20];  % [Fs Kt Chi Bt Kn]
pars = [5.5; 13.8; -0.95; -10; 12.35]
lpars = pars; lpars(copt.lspci) = 10.^(lpars(copt.lspci)); 

% %% Fully Stuck
Kst = Txyn*diag(reshape([copt.x.T{2}*lpars copt.y.T{2}*lpars copt.n.T{1}*lpars]', Npatches*3, 1))*Qxyn;
K0 = K + Kst;
X0 = K0\(Prestress*Fv);

% %% Prestress Analysis
opt = optimoptions('fsolve', 'Display', 'iter', 'SpecifyObjectiveGradient', true);

[Xstat, ~, eflg] = fsolve(@(X) WJRES_STR([X; 0.0], pars, K, Fv*0, Fv*Prestress, L, ...
    Txyn, Qxyn, CFUN, Npatches), X0, opt);
[~, dRstat, ~, dRdpstat] = WJRES_STR([Xstat; 0.0], pars, K, Fv*0, Fv*Prestress, L, ...
    Txyn, Qxyn, CFUN, Npatches);
dXdpstat = -dRstat\dRdpstat;
%% Modal Analysis
[Vst, Wst] = eigs(dRstat, M, 10, 'SM');
Wst = sqrt(diag(Wst));
Vst = Vst./sqrt(diag(Vst'*M*Vst)');
%% Load Experimental Data
% exx = load('../EXPERIMENTAL_DATA/BRB_EXP_RED.mat');
% Nx = 20; iNs = fix(linspace(1, length(exx.Q1x), Nx));
% expdat.Q = exx.Q1x(iNs);
% expdat.W = exx.W1x(iNs)*2*pi;
% expdat.Z = exx.Z1x(iNs);
% expdat.D = exx.D1x(iNs);

exx = load('../EXPERIMENTAL_DATA/EXP_14Mar2019.mat');
Nx = 20; iNs = fix(linspace(1, length(exx.Q), Nx));
expdat.Q = exx.Q(iNs)/abs(R(3,:)*Vst(:,mi));
expdat.W = exx.W(iNs)*2*pi;
expdat.Z = exx.Z(iNs);
expdat.D = exx.D(iNs);

%% Backbone function
opt.Display = 'off';

pars = [5.5; 13.8; -0.95; -10; 12.35]
Nqp = 4;
tic
[eobj, dedp, BB] = WJMODEL_BBFUN(pars, 1, expdat, K, M, X0, Fv*Prestress, L, Txyn, Qxyn, CFUN, Npatches, Nqp, opt);
% [eobj, dedp, BB] = WJMODEL_BBFUN_MASING(pars, 1, expdat, K, M, X0, Fv*Prestress, L, Txyn, Qxyn, CFUN, DFUN, Npatches, opt);
toc

save(sprintf('./DATS/P%d_S%.2f_%dLEV_BB_IWAN4.mat',Prestress,log10(sint),Nlev), 'BB');

%% Plotting

figure(10)
% clf()
subplot(2,1,1)
semilogx(BB.Q, BB.W/(2*pi), '.-'); hold on;
% semilogx(expdat.Q, expdat.W/(2*pi), 'k-');
subplot(2,1,2)
loglog(BB.Q, BB.Z+expdat.Z(1), '.-'); hold on;
% semilogx(expdat.Q, expdat.Z, 'k-');
% loglog(BB.Q, BB.D, '.-', expdat.Q, expdat.D, 'k-')
disp('Done')