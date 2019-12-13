clc
clear all
addpath('../ROUTINES/')
addpath('../ROUTINES/REMESHING/')
addpath('../ROUTINES/FEM/')
addpath('../ROUTINES/IQSMA/')
addpath('../ROUTINES/QSHMA/')

% ROM levels
% sel_method = 'P';
% Nels = 52;  % [52 100 140 204 240 304]

sel_method = 'U';
Nels = 568; % [36 74 122 232 448 568 588]
% for Nels=[36 74 122 232 448 568 588]

% sel_method = 'PD';
% Nels = 252;  % [48 108 152 200 252 292]

% FRICTMODEL STUDY
% sel_method = 'PD';
% Nels = 68;

load(sprintf('../FULLJOINT_ROM_PREPARE/ROMS/ROM_%s_%dELS.mat',sel_method,Nels), 'M', 'K', 'R', 'L', 'Fv', 'TFMf', 'TFMfhcb', 'MESH');
[Q1, T1] = ZTE_ND2QP(MESH,1);

ref = load('../FULLJOINT_ROM_PREPARE/DATS/statsol_P11580_0.000001.mat', 'Tn', 'Qm', 'Tm', 'MESH');

Prestress = 11580;
mds = [1 3 5];
mi = 1;
%% Interface Parameters
mu = 0.20;  %% Unknown parameters
c = 1.0; 

nu   = 0.29;                % Poisson's ratio
Aint = sum(sum(T1));      % Average interface area (lumped)
Pint = Prestress*3/Aint;           % Average lumped interface pressure
sint = 1e-6*c;              % Average interface asperity SD (exponential distribution)
chi  = 2.0;                 % Mindlin constant
ktkn = chi*(1-nu)/(2-nu);   % Tangential-normal stiffness ratio
kt   = 4*(1-nu)*Pint/(sqrt(pi)*(2-nu)*sint);    % Tangential stiffness
kn   = kt/ktkn;             % Normal Stiffness

pars.mu	= mu;             % Friction Coefficient - Unsure Parameter Selection
pars.ktx= kt;
pars.kty= kt;
pars.kn = kn;

%% ZTE Quadrature Matrices
No   = 5;                  % Number of GLQ integration points in each direction
[Qm,Tm] = ZTE_ND2QP(MESH,No);
QuadMats.Q = Qm;
QuadMats.T = Tm;%*Qm*inv(Qm'*Qm)*Qm';
SCP = (Qm'*Qm)\Qm';

%% Anonymous Function Declaration for forcing function                                         
NLFORCINGFUNC    = @(U, prev) ZTECONTACTFORCEJAC(@JENKINS_P, U, ...
    MESH, 0., pars, QuadMats, prev, TFMfhcb);
CONTACTFUNC = @(uxyn,uxyntxynp) JENKINS_P(uxyn,0,pars,uxyntxynp);

%% Stuck interface initial guess
Kstuck = zeros(size(L,1));
Kstuck(1:3:MESH.Nn*3, 1:3:MESH.Nn*3) = Tm*kt*Qm;
Kstuck(2:3:MESH.Nn*3, 2:3:MESH.Nn*3) = Tm*kt*Qm;
Kstuck(3:3:MESH.Nn*3, 3:3:MESH.Nn*3) = Tm*kn*Qm;

% Kstuck = TFMfhcb*Kstuck;

K0 = K+L'*Kstuck*L;
X0 = K0\(Fv*Prestress);

%% Prestress Analysis
opt = optimoptions('fsolve', 'Display', 'iter', 'SpecifyObjectiveGradient', true);
Prestress = 11580;
cnum = 1e0;
Xguess = X0;
Cmat = cnum'*speye(length(Xguess));
Xguess = Cmat\Xguess;
prev.uxyntxyn = zeros(size(QuadMats.Q,1),6);

[Xstat,~,~,~,dR0] = fsolve(@(X) NLRES(NLFORCINGFUNC, [X; 0], K, L, Prestress*Fv, Fv*0, prev, Cmat), Xguess, opt);
Uph = L*Cmat*Xstat;
[Tx, Ty, Tn] = CONTACTFUNC(QuadMats.Q*reshape(Uph(1:(3*MESH.Nn)), 3, MESH.Nn)', prev.uxyntxyn);
Pstat = reshape([Tx Ty Tn]', [No*No*MESH.Ne*3 1]);

%% Plot prestress tractions
figure(Nels);
clf()
SHOW2DMESH(MESH.Nds+[0 0.02], MESH.Tri, MESH.Quad, Qm\Tn, -1, -100);
SHOW2DMESH(ref.MESH.Nds-[0 0.02], ref.MESH.Tri, ref.MESH.Quad, ref.Qm\ref.Tn, -1, -100); 
yy = colorbar('southoutside');
set(yy, 'fontsize', 14);
xlabel(yy, 'Normal Traction (Pa)')
colormap(jet(11))
axis equal; 
axis off;
title(sprintf('Net Force Ratio: %f\n', sum(Tm*Tn)/sum(ref.Tm*ref.Tn)))
caxis([0 2.5e7])

% return
%% Low-Amplitude Mode shape Extraction
[Vm, Dm]  = eigs(dR0, M, mds(mi)*3, 'SM');
[Dm, si] = sort(diag(Dm));
Vm   = Vm(:,mds(mi));
Vm	 = Vm/sqrt(Vm'*M*Vm);
Fl	 = M*Vm;

%% SETUP RESIDUE FUNCTION
mspca = 1e4;
Cmat = ones(length(Vm)+2,1);
Copt.lsrch = 1;

NLRESIDFUNC = @(Uw, qm, Xstart, prev) QSHMARES(NLFORCINGFUNC, [Uw; qm], ...
    M, K, L, Prestress*Fv, Xstart, prev, Cmat, mspca);
%% QSHMA Run
QS = logspace(-6, -1, 20);
RS = QS*0;
WS = QS*0;
ZS = QS*0;
ZSM = QS*0;
DSM = QS*0;

Nhp = 10;
HYSTS = cell(size(QS));
BBS = cell(size(QS));
ULS = cell(size(QS));
RLS = cell(size(QS));
MESH.dpn = 3;
QuadMatsGen.Q = kron(Qm, eye(3));
QuadMatsGen.T = kron(Tm, eye(3));

opt = optimoptions('fsolve', 'Display', 'off', 'SpecifyObjectiveGradient', true);
for k=1:length(QS)
%     [HYSTS{k}, BBS{k}, ULS{k}, RLS{k}] = GENSTEPPED_HYSTERETICNMA(NLRESIDFUNC, ...
%         CONTACTFUNC, Xstat, [Vm; Dm(1)], Cmat, QS(k), Nhp, MESH, ...
%         L, opt, QuadMatsGen, Copt, 1);

    Copt.Nqp = 4;
	[BBS{k}] = GENSTEPPED_HYSTERETICNMA_LGL(NLRESIDFUNC, ...
        CONTACTFUNC, Xstat, [Vm; Dm(1)], Cmat, QS(k), Nhp, MESH, ...
        L, opt, QuadMatsGen, Copt, 1);
    % RESPONSE POINT DISPLACEMENT
    RS(k) = abs(R(3,:)*(BBS{k}.U(1:end-1,end)-Xstat));
    
    % FREQUENCY
    WS(k) = BBS{k}.W(end);

%     % FULL HYSTERESIS DISSIPATION
%     DS(k) = ZTEDISSIPATION(L(1:(3*MESH.Nn),:)*(HYSTS{k}.U(1:end-1,:)-Xstat), HYSTS{k}.P-Pstat, QuadMats, MESH);
%     ZS(k) = DS(k)/(2*pi*(WS(k)*QS(k))^2);

    % Integrated Dissipation
    DS(k) = BBS{k}.D;
    ZS(k) = (-4*max(BBS{k}.Q)*max(BBS{k}.A)+8*BBS{k}.D)/(2*pi*WS(k)*QS(k))^2;
    
    % MASING DISSIPATION
    qb = BBS{k}.Q;        ab = BBS{k}.A;
    qb = qb([1:end end]); ab = ab([1:end 1]);
    DSM(k) = -4*max(qb)*max(ab)-8*sum(0.5*qb.*(ab([2:end 1])-ab([end 1:end-1])));
    ZSM(k) = DSM(k)/(2*pi*WS(k)*QS(k))^2;
    
    fprintf('Done %d/%d: (%e; %.2f; %e; %e)\n', k, length(QS), QS(k), WS(k)/(2*pi), ZS(k), ZSM(k));    
end

save(sprintf('./DATS/RUN_M%d_%s_%dELS.mat',mi,sel_method,Nels), 'HYSTS', 'BBS', 'ULS', 'RLS', 'RS', 'QS', 'WS', 'ZS', 'DS', 'DSM', 'ZSM', 'Tx', 'Ty', 'Tn', 'Xstat')
%% Plot Backbone
figure(30)
clf()
yyaxis left
semilogx(QS, WS/(2*pi), '.-');
set(gca, 'XScale', 'log')
% ylim([80 180])
ylabel('Natural Frequency (Hz)')
xlabel('Modal Amplitude')

yyaxis right
loglog(QS, ZS+1e-4*0, '.-', QS, ZSM, 'o-')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
ylim([1e-5 1e0])
ylabel('Effective Damping Factor')
% end

gong = load('gong.mat');
sound(gong.y)