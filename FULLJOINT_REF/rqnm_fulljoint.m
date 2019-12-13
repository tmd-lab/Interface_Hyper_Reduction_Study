clc
clear all
addpath('../ROUTINES/')
addpath('../ROUTINES/REMESHING/')
addpath('../ROUTINES/FEM/')
addpath('../ROUTINES/IQSMA/')
addpath('../ROUTINES/QSHMA/')

setid = 5;
fname = sprintf('./PREP/%d_SET_NULLRED.mat', setid);
load(fname,'M','K','R','Fv','L','T','MESH');
mds = [1 3 5];
Prestress = 11580;
mi = 1;  % Mode of interest
[Q1,T1] = ZTE_ND2QP(MESH,1);

% %% Detect Holts
% [~,MESH.Ndbs] = CreateBRBInterface(MESH.Nds(:,1), MESH.Nds(:,2));
% MESH.Pstiff = 1e12;
% MESH.clrn = 1e-4;
% 
% Kpin = zeros(size(L,1));
% for i=1:length(MESH.Ndbs)
%     Kpin((MESH.Ndbs{1}-1)*3+1, (MESH.Ndbs{1}-1)*3+1) = eye(length(MESH.Ndbs{1}))*MESH.Pstiff;
%     Kpin((MESH.Ndbs{1}-1)*3+2, (MESH.Ndbs{1}-1)*3+2) = eye(length(MESH.Ndbs{1}))*MESH.Pstiff;
% end
% 
% K = K + L'*Kpin*L;
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
    MESH, 0., pars, QuadMats, prev);
CONTACTFUNC = @(uxyn,uxyntxynp) JENKINS_P(uxyn,0,pars,uxyntxynp);

%% Stuck interface initial guess
Kstuck = zeros(size(L,1));
Kstuck(1:3:MESH.Nn*3, 1:3:MESH.Nn*3) = Tm*kt*Qm;
Kstuck(2:3:MESH.Nn*3, 2:3:MESH.Nn*3) = Tm*kt*Qm;
Kstuck(3:3:MESH.Nn*3, 3:3:MESH.Nn*3) = Tm*kn*Qm;

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
    ZS(k) = BBS{k}.Z;
    
    % MASING DISSIPATION
    qb = BBS{k}.Q;        ab = BBS{k}.A;
    qb = qb([1:end end]); ab = ab([1:end 1]);
    DSM(k) = -4*max(qb)*max(ab)-8*sum(0.5*qb.*(ab([2:end 1])-ab([end 1:end-1])));
    ZSM(k) = DSM(k)/(2*pi*WS(k)*QS(k))^2;
    
    fprintf('Done %d/%d: (%e; %.2f; %e; %e)\n', k, length(QS), QS(k), WS(k)/(2*pi), ZS(k), ZSM(k));    
end

save('./DATS/RUNS.mat', 'HYSTS', 'BBS', 'ULS', 'RLS', 'RS', 'QS', 'WS', 'ZS', 'DS', 'DSM', 'ZSM')

%% Plot Backbone
figure(10)
clf()
yyaxis left
semilogx(QS, WS/(2*pi), '.-');
set(gca, 'XScale', 'log')
ylabel('Natural Frequency (Hz)')
xlabel('Modal Amplitude')

yyaxis right
loglog(QS, ZS+1e-4*0, '.-', QS, ZSM+1e-4*0, 'o-')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
ylabel('Effective Damping Factor')