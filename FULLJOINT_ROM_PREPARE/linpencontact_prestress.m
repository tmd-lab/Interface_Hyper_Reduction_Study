clc
clear all
addpath('../ROUTINES/')
addpath('../ROUTINES/FEM/')

setid = 5;
load(sprintf('./MATS/%d_SET_NULLRED.mat', setid), 'M', 'K', 'R', 'L', 'Fv', 'MESH', 'Krel', 'Fvrel');
Prestress = 11580.0;

%% ZTE Quadrature Points
No = 2;
[Qm, Tm] = ZTE_ND2QP(MESH, No);
QuadMats.Q = Qm;
QuadMats.T = Tm;
Aint = full(sum(sum(Tm)));

[Q1, T1] = ZTE_ND2QP(MESH, 1);
%% Contact Model
nu = 0.29;
Pint = Prestress/Aint;
sint = 1e-6;
chi = 2.0;
ktkn = chi*(1-nu)/(2-nu);
kt = 4*(1-nu)*Pint/(sqrt(pi)*(2-nu)*sint);
kn = kt/ktkn;

pars.mu = 0.20;
pars.ktx = kt;
pars.kty = kt;
pars.kn = kn;

%% Forcing Function
NLFORCINGFUNC = @(U, prev) ZTECONTACTFORCEJAC(@JENKINS_P, U, ...
    MESH, 0, pars, QuadMats, prev);
CONTACTFUNC = @(uxyn, uxyntxynp) JENKINS_P(uxyn, 0, pars, uxyntxynp);

prev.uxyntxyn = zeros(size(QuadMats.Q,1),6);

%% Fixed Interface Initial Guess
Kstuck = sparse(zeros(size(L,1)));
Kstuck(1:3:MESH.Nn*3, 1:3:MESH.Nn*3) = pars.ktx*Tm*Qm;
Kstuck(2:3:MESH.Nn*3, 2:3:MESH.Nn*3) = pars.kty*Tm*Qm;
Kstuck(3:3:MESH.Nn*3, 3:3:MESH.Nn*3) = pars.kn*Tm*Qm;

K0 = K + L'*Kstuck*L;
%% Solution Step
opt = optimoptions('fsolve', 'Display', 'iter', 'SpecifyObjectiveGradient', true);
Xguess = K0\(Fv*Prestress);

[Xs, ~, ~, ~, dR0] = fsolve(@(X) NLRES(NLFORCINGFUNC, [X;Prestress], K, L, Fv*0, Fv, prev), Xguess, opt);
Uph = L*Xs;
[Tx, Ty, Tn] = CONTACTFUNC(QuadMats.Q*reshape(Uph(1:3*MESH.Nn), 3, MESH.Nn)', prev.uxyntxyn);

%% Linear Modal Analysis
[V, D] = eigs(dR0, M, 30, 'SM');
Ws = sqrt(diag(D));
V = V./sqrt(diag(V'*M*V)');
V = L*V;

%% Relative-Absolute Coordinates
Nint = MESH.Nn;
Nrest = size(L,1)-Nint*3*2;
Trel = sparse([eye(Nint*3),  eye(Nint*3), zeros(Nint*3, Nrest);
               zeros(Nint*3), eye(Nint*3), zeros(Nint*3, Nrest);
               zeros(Nrest, Nint*3*2),     eye(Nrest)]);
V = Trel*V;

save(sprintf('./DATS/statsol_P%d_%f.mat',Prestress,sint), 'Xs', 'Qm', 'Tm', 'Tx', 'Ty', 'Tn', 'dR0', 'V', 'Ws', 'MESH')

%% Plot
figure(10)
clf()
SHOW2DMESH(MESH.Nds, MESH.Tri, MESH.Quad, Qm\Tn, -1, -100); axis equal; axis off
colorbar('southoutside')
colormap(jet(11))