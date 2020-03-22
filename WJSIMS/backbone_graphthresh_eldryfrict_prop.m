clc
clear all
addpath('../ROUTINES/')
addpath('../ROUTINES/FEM/')
addpath('../ROUTINES/IWAN/')
addpath('../ROUTINES/WJSIMS')

% 5 Patches
setid = 5;
fname = sprintf('../WHOLEJOINT_ROM_PREPARE/REDMATS/%d_SET_WJMAT.mat', setid);

% ThreshRed
Prestress = 11580; 
sint = 1e-6; 

sel_method = 'PD' % 'P'[2,3,4,5,6,10] 'PD' [2,3,6,7,8,10,11] 'U' [3 , 5 , 7, 9, 11, 19, 21, 25, 30, 40, 45, 50]
for Nlev=[2,3,6,7,10]

fname = sprintf('../WHOLEJOINT_ROM_PREPARE/REDMATS/P%d_S%.2f_%s_%dLEV_GRED_WJMAT_prop.mat',Prestress, log10(sint), sel_method, Nlev);
load(fname, 'M', 'K', 'L', 'R', 'Fv', 'Th', 'I2', 'cnum', 'MESH', 'PatchAreas', 'Npatches', 'NTG', 'GTG', 'dofred', 'GTG', 'Pels')

Prestress = 11580;
mi = 1;
%% Access Matrices and functions
Qxyn = L(reshape(((1:Npatches)-1)*6+(1:3)', Npatches*3, 1), :);
Txyn = L(reshape(((1:Npatches)-1)*6+(1:3)', Npatches*3, 1), :)';
Qtxyn= L(reshape(((1:Npatches)-1)*6+(4:6)', Npatches*3, 1), :);

copt.lspci = [2 3];  % [Kt Kn]  log scale

copt.x.T{1} = [ones(Npatches,1) zeros(Npatches, 2)];  % mu
copt.x.T{2} = [zeros(Npatches,1) PatchAreas' zeros(Npatches, 1)];  % Kt

copt.y.T{1} = [ones(Npatches,1) zeros(Npatches, 2)];  % mu
copt.y.T{2} = [zeros(Npatches,1) PatchAreas' zeros(Npatches, 1)];  % Kt

copt.n.T{1} = [zeros(Npatches, 2) PatchAreas'];  % Kn

% CFUN = @(uxyn, pars) SPRING_JENKINS(uxyn, pars, copt);
% CFUN = @(uxyn, pars) PENALTY_JENKINS(uxyn, pars, copt);
% CFUN = @(uxyn, pars) PENALTY_ELDRYFRICT(uxyn, pars, copt);
CFUN = @(uxyn, pars) PENALTY_ELDRYFRICT_2D(uxyn, pars, copt);  % Inappropriate as long as rotational DoFs are not reconciled
% DFUN = @(uxyn, pars) PENALTY_IWAN4_DISS(uxyn, pars, copt);

%% Load Experimental Data
% exx = load('../EXPERIMENTAL_DATA/BRB_EXP_RED.mat');
% Nx = 20; iNs = fix(linspace(1, length(exx.Q1x), Nx));
% expdat.Q = exx.Q1x(iNs);
% expdat.W = exx.W1x(iNs)*2*pi;
% expdat.Z = exx.Z1x(iNs);
% expdat.D = exx.D1x(iNs);

% exx = load('../EXPERIMENTAL_DATA/EXP_14Mar2019.mat');
% Nx = 20; iNs = fix(linspace(1, length(exx.Q), Nx));
% expdat.Q = exx.Q(iN   s)/abs(R(3,:)*Vst(:,mi));
% expdat.W = exx.W(iNs)*2*pi;
% expdat.Z = exx.Z(iNs);
% expdat.D = exx.D(iNs);

exx = load('../FULLJOINT_REF/DATS/RUN2.mat', 'QS', 'WS', 'ZS', 'DS','ZSM', 'BBS');
% exx = load('OUTPUTS.mat', 'QS', 'WS', 'ZS', 'DS');
expdat.Q = exx.QS;
expdat.W = exx.WS;
expdat.Z = exx.ZS;
expdat.ZM = exx.ZSM;
expdat.D = exx.DS;

%% Contact Model Parameters
mu = 0.20;
c = 1.0;

nu = 0.29;
Aint = sum(PatchAreas);
Pint = Prestress*3/Aint;
sint = 1e-6*c;
chi = 2.0;
ktkn = chi*(1-nu)/(2-nu);
kt = 4*(1-nu)*Pint/(sqrt(pi)*(2-nu)*sint);
kn = kt/ktkn;

% pars = log10([0.20*Pint; kt; kn]);  % [Fs Kt Kn]
pars = [mu; log10([kt; kn])];

%% Linearly Augmenting theta DoFs
% Ktheta_lin = L(4:6:6*Npatches, :)'*diag(I2(1,:)*kn)*L(4:6:6*Npatches, :)*0 ...
%     + L(4:6:6*Npatches, :)'*diag(I2(2,:)*kn)*L(5:6:6*Npatches, :)*0 ...
%     + L(5:6:6*Npatches, :)'*diag(I2(2,:)*kn)*L(4:6:6*Npatches, :)*0 ...
%     + L(5:6:6*Npatches, :)'*diag(I2(3,:)*kn)*L(5:6:6*Npatches, :)*0 ...
%     + L(6:6:6*Npatches, :)'*diag(I2(4,:)*kt)*L(6:6:6*Npatches, :);
% K = K+Ktheta_lin;

%%
lpars = pars; lpars(copt.lspci) = 10.^(lpars(copt.lspci)); 

% %% Fully Stuck
Kst = Txyn*diag(reshape([copt.x.T{2}*lpars copt.y.T{2}*lpars copt.n.T{1}*lpars]', Npatches*3, 1))*Qxyn;
K0 = K + Kst;
X0 = K0\(Prestress*Fv);
% Lto = load(sprintf('../WHOLEJOINT_ROM_PREPARE/MATS/%d_SET_NULLRED.mat', setid), 'L');
% Xtmp = [Lto.L*exx.BBS{1}.U(1:end-1,1); GTG\(NTG'*Lto.L(1:MESH.Nn*3,:)*exx.BBS{1}.U(1:end-1,1)); zeros(Npatches*6,1)];
% X0 = Th\Xtmp;

% %% Prestress Analysis
opt = optimoptions('fsolve', 'Display', 'iter', 'SpecifyObjectiveGradient', true, ...
    'MaxIterations', 400, 'MaxFunctionEvaluations', 500);

[Xstat, ~, eflg] = fsolve(@(X) WJRES_STR([X; 0.0], pars, K, Fv*0, Fv*Prestress, L, ...
    Txyn, Qxyn, CFUN, Npatches), X0, opt);
[~, dRstat, ~, dRdpstat] = WJRES_STR([Xstat; 0.0], pars, K, Fv*0, Fv*Prestress, L, ...
    Txyn, Qxyn, CFUN, Npatches);
dXdpstat = -dRstat\dRdpstat;
Fxyn = CFUN([Qxyn(1:3:end,:)*Xstat Qxyn(2:3:end,:)*Xstat Qxyn(3:3:end,:)*Xstat], pars);
Pxyn = Fxyn./PatchAreas';  % Tractions on each patch
%% Modal Analysis
[Vst, Wst] = eigs(dRstat, M, 10, 'SM');
Wst = sqrt(diag(Wst));
Vst = Vst./sqrt(diag(Vst'*M*Vst)');

%% Backbone function
opt.Display = 'off';
opt.SpecifyObjectiveGradient = true;
opt.MaxIterations = 1000;

Nqp = 20;
tic
[eobj, dedp, BB] = WJMODEL_BBFUN(pars, 1, expdat, K, M, R, X0, Fv*Prestress, L, Txyn, Qxyn, CFUN, Npatches, Nqp, opt);
% [eobj, dedp, BB] = WJMODEL_BBFUN_MASING(pars, 1, expdat, K, M, X0, Fv*Prestress, L, Txyn, Qxyn, CFUN, DFUN, Npatches, opt);
ttk = toc

save(sprintf('./DATS/P%d_S%.2f_%s_%dLEV_BB_ELDRYFRICT_prop.mat',Prestress,log10(sint),sel_method,Nlev), 'BB', 'ttk', 'Pxyn', 'MESH', 'Pels');

%% Plotting

figure(10)
% clf()
subplot(2,1,1)
semilogx(BB.Q, BB.W/(2*pi), '.-'); hold on;
semilogx(expdat.Q, expdat.W/(2*pi), 'ko-');
% xlim([1e-6 1e-2])
subplot(2,1,2)
loglog(BB.Q, BB.Z, '.-'); hold on;
semilogx(expdat.Q, expdat.ZM, 'ko-');
% loglog(BB.Q, BB.D, '.-', expdat.Q, expdat.D, 'k-')
% xlim([1e-6 1e-2])
ylim([1e-5 1e0])
disp('Done')

% figure(11)
% V_ref = cell2mat(cellfun(@(c) c.MS(:,end), exx.BBS, 'UniformOutput', false));
% V_rom = Th(1:size(V_ref,1), :)*BB.Phi(:, :);
% % br=bar3(1-(V_rom'*V_ref).^2./(diag(V_rom'*V_rom).*diag(V_ref'*V_ref)'))
% % for k=1:size(br,2); br(k).CData = br(k).ZData; br(k).FaceColor = 'interp'; end
% imagesc(1-abs((V_rom'*V_ref).^2./(diag(V_rom'*V_rom).*diag(V_ref'*V_ref)')))
% yy=colorbar; ylabel(yy, '1-MAC');
% set(gca, 'YDir', 'normal')
% xlabel('REF')
% ylabel('ROM')
% axis equal
% xlim([1 length(exx.QS)])
% ylim([1 length(exx.QS)])
end