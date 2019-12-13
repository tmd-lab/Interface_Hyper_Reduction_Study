clc
clear all
addpath('../ROUTINES/')
addpath('../ROUTINES/FEM/')
addpath('../ROUTINES/IWAN/')

setid = 5;
fname = sprintf('../WHOLEJOINT_ROM_PREPARE/REDMATS/%d_SET_WJMAT.mat', setid);
load(fname, 'M', 'K', 'L', 'R', 'Fv', 'Th', 'LamT', 'cnum', 'MESH', 'PatchAreas', 'Npatches', 'dofred', 'GTG')
Prestress = 11580;

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
%% RQNM
mi = 1;
q = 1e-10;
prev.uxyn = zeros(Npatches, 3);
prev.Fxyn = zeros(Npatches, 3);
prev.dFxyn = zeros(Npatches, 3);
prev.dFxyndpar = zeros(Npatches, length(pars), 3);
prev.dFxyndXdpar = zeros(Npatches*3, length(pars));

opt = optimoptions('fsolve', 'Display', 'iter', 'SpecifyObjectiveGradient', true);
Xd = fsolve(@(X) WJRES_RQNM([X; q], pars, prev, K, M, Xstat, dXdpstat, ...
    Fv*Prestress, L, Txyn, Qxyn, CFUN, Npatches), [Xstat+Vst(:,mi)*q; Wst(mi)^2], opt);
[~, dRdX, ~, dRdp] = WJRES_RQNM([Xd; q], pars, prev, K, M, Xstat, dXdpstat, Fv*Prestress, ...
    L, Txyn, Qxyn, CFUN, Npatches);
dXdp = -dRdX\dRdp;
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
% pars = PARSGA(18, :)'
% pars = mean(PARSGA)'
% pars = PARMIN_W;
% pars = PARS_BI(12,1:end-1)';
Nqp = 4;
tic
[eobj, dedp, BB] = WJMODEL_BBFUN(pars, 1, expdat, K, M, X0, Fv*Prestress, L, Txyn, Qxyn, CFUN, Npatches, Nqp, opt);
% [eobj, dedp, BB] = WJMODEL_BBFUN_MASING(pars, 1, expdat, K, M, X0, Fv*Prestress, L, Txyn, Qxyn, CFUN, DFUN, Npatches, opt);
toc

figure(10)
% clf()
subplot(2,1,1)
semilogx(BB.Q, BB.W/(2*pi), '.-', expdat.Q, expdat.W/(2*pi), 'k-'); hold on
subplot(2,1,2)
loglog(BB.Q, BB.Z+expdat.Z(1), '.-', expdat.Q, expdat.Z, 'k-'); hold on
% loglog(BB.Q, BB.D, '.-', expdat.Q, expdat.D, 'k-')
disp('Done')

%% Scalarization (weighted)
opts = struct('method', 'weighted', 'gradients', true, 'nobj', 2, 'npar', 1, 'nvar', length(pars));
opts = struct('method', 'sphericalwgts', 'gradients', true, 'nobj', 2, 'npar', 1, 'nvar', length(pars)); 
opts.rpt = [1; 1];
wgt = 0.5;
[O, dOdp, dOdw] = SCALARIZE(@(pars) WJMODEL_BBFUN(pars, 1, expdat, K, M, X0, Fv*Prestress, L, Txyn, Qxyn, CFUN, Npatches, Nqp, opt), ...
    [pars; wgt], opts);

%% Scalarized Optimization for utopian point
opts = struct('method', 'weighted', 'gradients', true, 'nobj', 2, 'npar', 1, 'nvar', length(pars));
fopts = optimoptions('fmincon', 'SpecifyObjectiveGradient', true, 'Display', 'iter');
% gopts = optimoptions('ga', 'Display', 'iter', 'MaxGenerations', 20);
lb = [3, 3, -1, -4, 3];
ub = [30, 20, 0, 4, 20];

% pars0 = [1, 2, -0.5, 3, 4]';
pars0 = [5; 13.8; -0.5; -4; 12.35];

wgt = 1.0;  % Frequency Alone
% [PARMIN_W, Wmin] = fmincon(@(prs) SCALARIZE(@(pars) WJMODEL_BBFUN(pars, 1, expdat, K, M, X0, Fv*Prestress, L, Txyn, Qxyn, CFUN, Npatches, Nqp, opt), [prs; wgt], opts), ...
%     pars0, [], [], [], [], lb, ub, [], fopts);  % constrained opt
[PARMIN_W, Wmin] = fminunc(@(prs) SCALARIZE(@(pars) WJMODEL_BBFUN(pars, 1, expdat, K, M, X0, Fv*Prestress, L, Txyn, Qxyn, CFUN, Npatches, Nqp, opt), [prs; wgt], opts), ...
    pars0, fopts);  % unconstrained opt
% [PARMIN_W, Wmin] = ga(@(prs) SCALARIZE(@(pars) WJMODEL_BBFUN(pars, 1, expdat, K, M, X0, Fv*Prestress, L, Txyn, Qxyn, CFUN, Npatches, Nqp, opt), [prs; wgt], opts), length(pars), ...
%     [], [], [], [], lb', ub');  % genetic algorithm

% pars0 = PARMIN_W;
wgt = 0.0;  % Damping Alone
% [PARMIN_Z, Zmin] = fmincon(@(prs) SCALARIZE(@(pars) WJMODEL_BBFUN(pars, 1, expdat, K, M, X0, Fv*Prestress, L, Txyn, Qxyn, CFUN, Npatches, Nqp, opt), [prs; wgt], opts), ...
%     pars0, [], [], [], [], lb, ub, [], fopts);  % constrained opt
[PARMIN_Z, Zmin] = fminunc(@(prs) SCALARIZE(@(pars) WJMODEL_BBFUN(pars, 1, expdat, K, M, X0, Fv*Prestress, L, Txyn, Qxyn, CFUN, Npatches, Nqp, opt), [prs; wgt], opts), ...
    pars0, fopts);  % unconstrained opt

% Utopian Point
opts.rpt = [Wmin; Zmin]

%% Spherical Weighting
opts = struct('method', 'sphericalwgts', 'gradients', true, 'nobj', 2, 'npar', 1, 'nvar', length(pars)); 
% opts.rpt = [Wmin; Zmin];
opts.rpt = [0; 0];

THETAS_SW = linspace(0, pi/2, 10);
THETAS_SW = [0 0.001 0.005 0.01 0.015 0.02 0.025 0.035 0.05 0.075 0.10 0.15 0.20 0.25 0.75 1.0]*pi/2;
ERRS_SW = zeros(length(THETAS_SW), 2);
PARS_SW = zeros(length(THETAS_SW), length(pars));

PARS_SW(1,:) = PARMIN_W;
ERRS_SW(1,:) = WJMODEL_BBFUN(PARS_SW(1,:), 1, expdat, K, M, X0, Fv*Prestress, L, Txyn, Qxyn, CFUN, Npatches, Nqp, opt);

pars0 = (PARMIN_W + PARMIN_Z)/2;
for k=1:length(THETAS_SW)
    [PARS_SW(k, :), err] = fminunc(@(prs) SCALARIZE(@(pars) WJMODEL_BBFUN(pars, 1, expdat, K, M, X0, Fv*Prestress, L, Txyn, Qxyn, CFUN, Npatches, Nqp, opt), [prs; THETAS_SW(k)], opts), ...
        pars0, fopts);  % unconstrained opt
    ERRS_SW(k,:) = WJMODEL_BBFUN(PARS_SW(k,:), 1, expdat, K, M, X0, Fv*Prestress, L, Txyn, Qxyn, CFUN, Npatches, Nqp, opt);
    
    fprintf('Done %d/%d\n', k, length(THETAS_SW));
end

figure(10)
% clf()
plot(ERRS_SW(:,1), ERRS_SW(:,2), '.'); hold on

%% Optimization (Genetic)
gopts = optimoptions('gamultiobj', 'Display', 'iter', ...
    'PlotFcn', {@gaplotpareto, @gaplotrankhist}, 'MaxGenerations', 100, ...
    'DistanceMeasureFcn', {@distancecrowding, 'genotype'}, ...
    'PopulationSize', 100); %, 'InitialPopulationMatrix', [PARS]);
opt.Display = 'off';
lb = [0, 0, -1, -4, 3];
ub = [30, 20, 0, 4, 20];
[PARSGA, ERRSGA] = gamultiobj(@(pars) WJMODEL_BBFUN(pars, 1, expdat, K, M, X0, Fv*Prestress, L, Txyn, Qxyn, CFUN, Npatches, Nqp, opt), length(pars), [], [], [], [], lb, ub, [], gopts);

figure(10)
plot(ERRSGA(:,1), ERRSGA(:,2), '*')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
grid on

%% Spherical Boundary Intersection
load('GA.mat');  Wmin = min(ERRSGA(:,1));   Zmin = min(ERRSGA(:,2));
opts = struct('method', 'sphericalbi', 'gradients', true, 'nobj', 2, 'npar', 2, 'nvar', length(pars)); 

% opts.rpt = [Wmin; Zmin];
% opts.rpt = [0; 0];
% opts.rpt = log10([Wmin; Zmin]);
% opts.rpt = log10([4e-5; 1.5e-3]);
opts.rpt = [-4.32; -2.74];

loglog(10^opts.rpt(1), 10^opts.rpt(2), 'kx')

% THETAS_BI = [linspace(90, 87.5, 20), linspace(87, 0, 20)];
% THETAS_BI = (1-[0 logspace(-3, 0, 29)])*pi/2;
THETAS_BI = linspace(0, pi/2, 15);

% THETAS_BI = [0 0.001 0.005 0.01 0.015 0.025 0.05 0.075 0.10 0.15 0.20 0.25 0.75 1.0]*pi/2;
ERRS_BI = zeros(length(THETAS_BI), 2);
PARS_BI = zeros(length(THETAS_BI), length(pars)+1);

lb = [0, 0, -1, -5, 3];
ub = [30, 20, 0, 5, 20];
% lb = [];
% ub = [];
fopts = optimoptions('fmincon', 'SpecifyObjectiveGradient', true, ...
    'Display', 'iter', 'SpecifyConstraintGradient', true, 'MaxIterations', 70);

% pars0 = (PARMIN_W + PARMIN_Z)/2;
pars0 = mean(PARSGA)';
t0 = 2;

% Tracking for how long each iteration is taking to estimate cluster time.
iter_times = zeros(1, length(THETAS_BI));

parfor k=1:length(THETAS_BI)
    
    %start timer
    tic;
    
    [PARS_BI(k,:), err] = fmincon(@(prs) deal(prs(end), [zeros(1,length(prs)-1) 1]), [pars0; t0], [], [], [], [], lb, ub, @(prs) SCALARIZE(@(pars) WJMODEL_BBFUN(pars, 1, expdat, K, M, X0, Fv*Prestress, L, Txyn, Qxyn, CFUN, Npatches, Nqp, opt), ...
        [prs; THETAS_BI(k)], opts), fopts);
    
    %temporary added so that parfor rules not broken. 
    tmp = PARS_BI(k,:);
    
    ERRS_BI(k,:) = WJMODEL_BBFUN(tmp(1,1:end-1), 1, expdat, K, M, X0, Fv*Prestress, L, Txyn, Qxyn, CFUN, Npatches, Nqp, opt);
    
	fprintf('Done %d/%d\n', k, length(THETAS_BI));
    
%     plot(ERRS_BI(k,1), ERRS_BI(k,2), 'ro', 'MarkerFaceColor', 'r');
%     loglog(10^ERRS_BI(k,1), 10^ERRS_BI(k,2), 'ro', 'MarkerFaceColor', 'r');
%     drawnow
    
    %read time
    iter_times(k) = toc;
end
toc

% figure(10)
% % clf()
% loglog(ERRS_BI(:,1), ERRS_BI(:,2), '.-'); hold on