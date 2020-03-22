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

sel_method = 'P' % 'P'[2,3,4,5,6,7,11] 'PD' [2,3,4,5,6,7,10] 'U' [3,5,7,9,11,19,21,30,40,45,50,60,70,80,90]
for Nlev=[7]
fname = sprintf('../WHOLEJOINT_ROM_PREPARE/REDMATS/P%d_S%.2f_%s_%dLEV_GRED_WJMAT_prop.mat',Prestress, log10(sint), sel_method, Nlev);
load(fname, 'M', 'K', 'L', 'R', 'Fv', 'Th', 'I2', 'cnum', 'MESH', 'PatchAreas', 'Npatches', 'NTG', 'GTG', 'dofred', 'GTG', 'Pels')

Prestress = 11580;
mi = 1;
%% Access Matrices and functions
Qxyntxyn = L(reshape(((1:Npatches)-1)*6+(1:6)', Npatches*6, 1), :);
Txyntxyn = L(reshape(((1:Npatches)-1)*6+(1:6)', Npatches*6, 1), :)';

copt.lspci = [2 3];  % [Kt Kn]  log scale

copt.x.T{1} = [ones(Npatches,1) zeros(Npatches, 2)];  % mu
copt.x.T{2} = [zeros(Npatches,1) PatchAreas' zeros(Npatches, 1)];  % Kt

copt.y.T{1} = [ones(Npatches,1) zeros(Npatches, 2)];  % mu
copt.y.T{2} = [zeros(Npatches,1) PatchAreas' zeros(Npatches, 1)];  % Kt

copt.n.T{1} = [zeros(Npatches, 2) PatchAreas'];  % Kn
copt.tx.T{1} = [zeros(Npatches, 1) I2(1,:)' zeros(Npatches, 1)];  % Ktxx
copt.tx.T{2} = [zeros(Npatches, 1) I2(2,:)' zeros(Npatches, 1)];  % Ktxy

copt.ty.T{1} = [zeros(Npatches, 1) I2(2,:)' zeros(Npatches, 1)];  % Ktyx
copt.ty.T{2} = [zeros(Npatches, 1) I2(3,:)' zeros(Npatches, 1)];  % Ktyy

copt.tn.T{1} = [zeros(Npatches, 2) I2(4,:)'];  % Ktzz
% CFUN = @(uxyn, pars) SPRING_JENKINS(uxyn, pars, copt);
% CFUN = @(uxyn, pars) PENALTY_JENKINS(uxyn, pars, copt);
% CFUN = @(uxyn, pars) PENALTY_ELDRYFRICT(uxyn, pars, copt);
% CFUN = @(uxyn, pars) PENALTY_ELDRYFRICT_2D(uxyn, pars, copt);  % Inappropriate as long as rotational DoFs are not reconciled
CFUN = @(uxyntxyn, pars) PENALTY_ELDRYFRICT_2D_WTH(uxyntxyn, pars, copt);  % kn-reconciliation of rotations
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

% %% Linearly Augmenting theta DoFs
% Ktheta_lin = L(4:6:6*Npatches, :)'*diag(I2(1,:)*kn)*L(4:6:6*Npatches, :)*0 ...
%     + L(4:6:6*Npatches, :)'*diag(I2(2,:)*kn)*L(5:6:6*Npatches, :)*0 ...
%     + L(5:6:6*Npatches, :)'*diag(I2(2,:)*kn)*L(4:6:6*Npatches, :)*0 ...
%     + L(5:6:6*Npatches, :)'*diag(I2(3,:)*kn)*L(5:6:6*Npatches, :)*0 ...
%     + L(6:6:6*Npatches, :)'*diag(I2(4,:)*kt)*L(6:6:6*Npatches, :);
% K = K+Ktheta_lin;

%%
lpars = pars; lpars(copt.lspci) = 10.^(lpars(copt.lspci)); 

% %% Fully Stuck
Kst = zeros(Npatches*6);
Kst(1:6:end,1:6:end) = diag(copt.x.T{2}*lpars); 
Kst(2:6:end,2:6:end) = diag(copt.y.T{2}*lpars);
Kst(3:6:end,3:6:end) = diag(copt.n.T{1}*lpars);
Kst(4:6:end,4:6:end) = diag(copt.tx.T{1}*lpars);Kst(4:6:end,5:6:end) = diag(copt.tx.T{2}*lpars);
Kst(5:6:end,4:6:end) = diag(copt.ty.T{1}*lpars);Kst(5:6:end,5:6:end) = diag(copt.ty.T{2}*lpars);
Kst(6:6:end,6:6:end) = diag(copt.tn.T{1}*lpars);
Kst = Txyntxyn*Kst*Qxyntxyn;
K0 = K + Kst;
X0 = K0\(Prestress*Fv);
% Lto = load(sprintf('../WHOLEJOINT_ROM_PREPARE/MATS/%d_SET_NULLRED.mat', setid), 'L');
% Xtmp = [Lto.L*exx.BBS{1}.U(1:end-1,1); GTG\(NTG'*Lto.L(1:MESH.Nn*3,:)*exx.BBS{1}.U(1:end-1,1)); zeros(Npatches*6,1)];
% X0 = Th\Xtmp;

% %% Prestress Analysis
opt = optimoptions('fsolve', 'Display', 'iter', 'SpecifyObjectiveGradient', true, ...
    'MaxIterations', 400, 'MaxFunctionEvaluations', 500);

[Xstat, ~, eflg,~,dRstat] = fsolve(@(X) WJRES_STR_WTH([X; 0.0], pars, K, Fv*0, Fv*Prestress, L, ...
    Txyntxyn, Qxyntxyn, CFUN, Npatches), X0, opt);
if eflg < 0
    disp('unconv');
end
[~, dRstat, ~, dRdpstat] = WJRES_STR_WTH([Xstat; 0.0], pars, K, Fv*0, Fv*Prestress, L, ...
    Txyntxyn, Qxyntxyn, CFUN, Npatches);
dXdpstat = -dRstat\dRdpstat;
Fxyn = CFUN([Qxyntxyn(1:6:end,:)*Xstat Qxyntxyn(2:6:end,:)*Xstat Qxyntxyn(3:6:end,:)*Xstat Qxyntxyn(4:6:end,:)*Xstat Qxyntxyn(5:6:end,:)*Xstat Qxyntxyn(6:6:end,:)*Xstat], pars);
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
[eobj, dedp, BB] = WJMODEL_BBFUN_WTH(pars, 1, expdat, K, M, R, X0, Fv*Prestress, L, Txyntxyn, Qxyntxyn, CFUN, Npatches, Nqp, opt);
% [eobj, dedp, BB] = WJMODEL_BBFUN_MASING(pars, 1, expdat, K, M, X0, Fv*Prestress, L, Txyn, Qxyn, CFUN, DFUN, Npatches, opt);
ttk = toc

save(sprintf('./DATS/P%d_S%.2f_%s_%dLEV_BB_ELDRYFRICT_prop_wth.mat',Prestress,log10(sint),sel_method,Nlev), 'BB', 'ttk', 'Pxyn', 'MESH', 'Pels', 'Xstat');

%% Plotting

% figure(10)
% % clf()
% subplot(2,1,1)
% semilogx(BB.Q, BB.W/(2*pi), '.-'); hold on;
% semilogx(expdat.Q, expdat.W/(2*pi), 'ko-');
% % xlim([1e-6 1e-2])
% subplot(2,1,2)
% loglog(BB.Q, BB.Z, '.-'); hold on;
% semilogx(expdat.Q, expdat.ZM, 'ko-');
% % loglog(BB.Q, BB.D, '.-', expdat.Q, expdat.D, 'k-')
% % xlim([1e-6 1e-2])
% ylim([1e-5 1e0])
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