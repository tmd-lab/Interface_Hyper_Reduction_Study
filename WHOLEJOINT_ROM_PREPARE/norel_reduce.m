clc
clear all
addpath('../ROUTINES/')
addpath('../ROUTINES/FEM/')

setid = 5;
fil = {'1_AROUNDSET', '2_ABOVESET', '3_SINGELEMABOVESET', '4_INTSET', '5_INTSETNPS'};
Ngens = [26, 26, 26, 26, 19];
fname = sprintf('../MATRIX_EXTRACTION/RUNS/%s/BRB_WOPRES_MAT.mat', fil{setid});
load(fname, 'M', 'K', 'R', 'Fv');
K = 0.5*(K+K');
M = 0.5*(M+M');
Fv = Fv';
R = R';
fname = sprintf('./MATS/%d_SET_NULLRED.mat', setid);
load(fname, 'MESH')

%% Patches - Selection
Npatches = 5;
xmin = min(MESH.Nds(:, 1));  xmax = max(MESH.Nds(:, 1));  xrng = range(MESH.Nds(:, 1));
% xthrs = linspace(xmin, xmax, Npatches+1);  xthrs(2:end-1) = xthrs(2:end-1) + 1.27e-3;
xthrs = xmin+xrng/8*[0 1 3 5 7 8]; % + [0 1 1 1 1 0]*1.27e-3;
Pnds = cell(Npatches, 1);
for i=1:Npatches
    Pnds{i} = find(MESH.Nds(:,1)>=(xthrs(i)-1.27e-3) & MESH.Nds(:,1)<=(xthrs(i+1)+1.27e-3));
end
% Pnds = {1:77, 67:256, 246:435, 425:614, 604:680};
%% PATCHING
[Q1, T1] = ZTE_ND2QP(MESH, 1);

Pels = cell(1, Npatches);
Area = cell(1, Npatches);
qps = cell(1, Npatches);
ctrds = zeros(Npatches, 2);
for n=1:5
    Pels{n} = find(sum(ismember(MESH.Quad(:, 2:end), Pnds{n}), 2)==4);
    Area{n} = T1(Pnds{n}, Pels{n})*ones(length(Pels{n}), 1);
    qps{n} = Q1(Pels{n}, :)*MESH.Nds;
    ctrds(n, :) = sum(Area{n}.*MESH.Nds(Pnds{n}, :))/sum(Area{n});
    
    k=n;SHOW2DMESH(MESH.Nds, [], MESH.Quad(Pels{k},:), k, -1, -100)
end
PatchAreas = cellfun(@(c) sum(c), Area);

%% Weak Form Integral Matrices
NTN = sparse(MESH.Nn*3, MESH.Nn*3);
NTG = sparse(MESH.Nn*3, Npatches*6);
GTG = sparse(Npatches*6, Npatches*6);
BNV = sparse(MESH.Nn, Npatches);  % Force Transformation Matrix
for n=1:Npatches
    [P, Nums, NTNmat, NTGmat, GTGmat] = CONSPATCHMAT(MESH.Nds, [], MESH.Quad(Pels{n}, :), ctrds(n, :));
    
    NTN = NTN + NTNmat;
    NTG(:, (n-1)*6+(1:6)) = NTG(:, (n-1)*6+(1:6)) + NTGmat;
    GTG((n-1)*6+(1:6), (n-1)*6+(1:6)) = GTG((n-1)*6+(1:6), (n-1)*6+(1:6)) + GTGmat;
    
    BNV(Pnds{n}, n) = 1.0/sum(Area{n});
end
GN = GTG\(NTG');

%% Reduced Order Model
nit = 1:(MESH.Nn*3);
nib = MESH.Nn*3+(1:(MESH.Nn*3));
mi = setdiff(1:length(K), [nit nib]);
cnum = 1e3;
cnum = 1.0;
Cmat = sparse(diag([ones(length(K)+2*Npatches*6, 1); ones(2*Npatches*6,1)*cnum]));

% Ka = sparse(blkdiag(Krel, zeros(Npatches*6*2)));
Ka = sparse([K([nit nib],:), zeros(length([nit nit]), 2*Npatches*6), -blkdiag(NTG, NTG);
             K(mi,:), zeros(length(mi), 4*Npatches*6);
             zeros(2*Npatches*6, size(K,2)+2*Npatches*6), blkdiag(GTG,GTG);
             -blkdiag(NTG', NTG'), zeros(2*Npatches*6, length(mi)), blkdiag(GTG,GTG), zeros(2*Npatches*6)]);
Ka = Ka*Cmat;
Ka = 0.5*(Ka+Ka');
Ma = sparse(blkdiag(M, zeros(2*Npatches*6*2)));
Ma = Ma*Cmat;
Ma = 0.5*(Ma+Ma');
Fva = [Fv; zeros(2*Npatches*6*2,1)];
Ftmpa = Ka*sparse([zeros(length(K)+2*Npatches*6, 1); cnum^(-1)*NTG\Fv(nit); cnum^(-1)*NTG\Fv(nib)]);

Ra = sparse([R zeros(size(R,1), 2*Npatches*6*2)]);

Ta = sparse([speye(length(K)+2*Npatches*6); cnum^(-1)*blkdiag(NTG,NTG)\[K([nit nib], :) zeros(length([nit nib]), 2*Npatches*6)]; ]);
% FIRST APPROXIMATION ENCOUNTERED HERE DUE TO NTG^+ TERM TO EXPRESS Lambda^V
Knew = Ta'*Ka*Ta;       Knew = 0.5*(Knew+Knew');
Mnew = Ta'*Ma*Ta;       Mnew = 0.5*(Mnew+Mnew');
Fvnew = Ta'*(Fva+Ftmpa);
Rnew = Ra*Ta;

vdofs = length(K) + reshape(1:(2*Npatches*6), 6, 2*Npatches);
% vdofs = reshape(vdofs, 2*Npatches*6, 1);  % All virtual DoFs
vdofs = reshape(vdofs(1:3, :), 2*Npatches*3, 1);   % Only Displacement DoFs (AVOID)

[Mhcb, Khcb, Thcb] = HCBREDUCE(Mnew, Knew, vdofs, Ngens(setid));
Mhcb = 0.5*(Mhcb+Mhcb');  Khcb = 0.5*(Khcb+Khcb');
Fvhcb = Thcb'*Fvnew;
Rhcb = Rnew*Thcb;
disp('DONE!')


%% HCB With the Lagrange Multipliers
vdofs = length(K) + reshape(1:(4*Npatches*6), 6, 4*Npatches);
vdofs = reshape(vdofs, 4*Npatches*6, 1);  % All virtual DoFs
[MhcbO, KhcbO, ThcbO] = HCBREDUCE(Ma, Ka, vdofs, Ngens(setid));
MhcbO = 0.5*(MhcbO+MhcbO');  KhcbO = 0.5*(KhcbO+KhcbO');
disp('DONE!')


% [Mf, Kf, Tf] = HCBREDUCE(MhcbO, KhcbO, 1:(length(KhcbO)-2*Npatches*6), 1);

%% Much Older Form
nit = 1:(MESH.Nn*3);
nib = MESH.Nn*3+(1:(MESH.Nn*3));
mi = setdiff(1:length(K), [nit nib]);
cnum = 1e3;
cnum = 1.0;
Cmat = sparse(diag([ones(length(K)+2*Npatches*6, 1); ones(2*Npatches*6,1)*cnum]));

Ka = [K([nit nib], [nit nib])+blkdiag(NTN,NTN), zeros(length([nit nib]),length(mi)), -blkdiag(NTG,NTG); ...
      K(mi, :), zeros(length(mi),2*Npatches*6); ...
      -blkdiag(NTG',NTG'), zeros(2*Npatches*6, length(mi)), blkdiag(GTG, GTG)];
Ma = blkdiag(M, zeros(2*Npatches*6));
Fva = [Fv; zeros(2*Npatches*6,1)];
LamTa = [blkdiag(NTG,NTG); zeros(length(mi)+2*Npatches*6,2*Npatches*6)];
Ra = [R zeros(size(R,1), 2*Npatches*6)];

vdofs = length(K) + reshape(1:(2*Npatches*6), 6, 2*Npatches);
vdofs = reshape(vdofs, 2*Npatches*6, 1);  % All virtual DoFs
% vdofs = reshape(vdofs(1:3, :), 2*Npatches*3, 1);   % Only Displacement DoFs (AVOID)

[Mhcb, Khcb, Thcb] = HCBREDUCE(Ma, Ka, vdofs, Ngens(setid));
Mhcb = 0.5*(Mhcb+Mhcb');  Khcb = 0.5*(Khcb+Khcb');
Fvhcb = Thcb'*Fva;
Rhcb = Ra*Thcb;
LamThcb = Thcb'*LamTa;
disp('DONE!')
