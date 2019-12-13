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

%% Reduction (after substitution of Lagrange Multipliers)
nit = 1:(MESH.Nn*3);
nib = MESH.Nn*3+(1:(MESH.Nn*3));
mi = setdiff(1:length(K), [nit nib]);

cnum = 1e9;
Ka = [K([nit nib], [nit nib])+cnum*blkdiag(NTN,NTN), zeros(length([nit nib]),length(mi)), -cnum*blkdiag(NTG,NTG); ...
      K(mi, :), zeros(length(mi),2*Npatches*6); ...
      -cnum*blkdiag(NTG',NTG'), zeros(2*Npatches*6, length(mi)), cnum*blkdiag(GTG, GTG)];
Ma = blkdiag(M, zeros(2*Npatches*6));
Fva = [Fv; zeros(2*Npatches*6,1)];
LamTa = [blkdiag(NTG,NTG); zeros(length(mi)+2*Npatches*6,2*Npatches*6)];
Ra = [R zeros(size(R,1), 2*Npatches*6)];

dofred = 6;
vdofs = length(K) + reshape(1:(2*Npatches*6), 6, 2*Npatches);
vdofs = reshape(vdofs(1:dofred, :), 2*Npatches*dofred, 1);

[Mhcb, Khcb, Thcb] = HCBREDUCE(Ma, Ka, vdofs, Ngens(setid));
Mhcb = 0.5*(Mhcb+Mhcb');  Khcb = 0.5*(Khcb+Khcb');
Fvhcb = Thcb'*Fva;
Rhcb = Ra*Thcb;
LamThcb = Thcb'*LamTa;
disp('DONE!')

%% Relative Transformation
ngenr = length(Mhcb)-2*Npatches*dofred;
Trel = [eye(Npatches*dofred) eye(Npatches*dofred) zeros(Npatches*dofred, ngenr);
        zeros(Npatches*dofred) eye(Npatches*dofred) zeros(Npatches*dofred, ngenr);
        zeros(ngenr, 2*Npatches*dofred) eye(ngenr)];
Mrel = Trel'*Mhcb*Trel;
Krel = Trel'*Khcb*Trel;
Fvrel = Trel'*Fvhcb;
Rrel = Rhcb*Trel;
LamTrel = Trel'*LamThcb;

%% NULL-SPACE REDUCTION
[Vrel, Drelf] = eig(full(Krel(Npatches*dofred+1:end, Npatches*dofred+1:end)), full(Mrel(Npatches*dofred+1:end, Npatches*dofred+1:end))); 
% [Vhcb, Dhcb] = eig(full(Khcb), full(Mhcb));
[Drelf, si] = sort(diag(Drelf));  Vrel = Vrel(:, si);
% Vrbms = Vhcb(:, 1:6);
Vrbms = [zeros(Npatches*dofred, 6); Vrel(:, 1:6)];  Vrbms = Vrbms./sqrt(diag(Vrbms'*Mrel*Vrbms)');
L = null(Vrbms'*Mrel);

M = L'*Mrel*L;
M = 0.5*(M+M');
K = L'*Krel*L;
K = 0.5*(K+K');
Fv = L'*Fvrel;
R = Rrel*L;
Th = Trel*L;
LamT = L'*LamTrel;
[sort(eig(K), 'descend') sort(eig(M), 'descend')]
disp(cnum)

%% Save
fname = sprintf('./REDMATS/%d_SET_WJMAT.mat', setid);
save(fname, 'M', 'K', 'L', 'Fv', 'R', 'Th', 'LamT', 'MESH', 'PatchAreas', 'Npatches', 'dofred', 'NTN', 'GTG', 'NTG', 'cnum')
disp('SAVED!')