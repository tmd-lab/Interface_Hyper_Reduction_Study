clc
clear all
addpath('../ROUTINES/')
addpath('../ROUTINES/FEM/')
addpath('../ROUTINES/GRAPH/')

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
Area = cell(1, Npatches);
qps = cell(1, Npatches);
ctrds = zeros(Npatches, 2);
for n=1:Npatches
    Pels{n} = find(sum(ismember(MESH.Quad(:, 2:end), Pnds{n}), 2)==4);
    Area{n} = T1(Pnds{n}, Pels{n})*ones(length(Pels{n}), 1);
    qps{n} = Q1(Pels{n}, :)*MESH.Nds;
    ctrds(n, :) = sum(Area{n}.*MESH.Nds(Pnds{n}, :))/sum(Area{n});
    
    k=n;SHOW2DMESH(MESH.Nds, [], MESH.Quad(Pels{k},:), k, -1, -100)
end
PatchAreas = cellfun(@(c) sum(c), Area);

%% Relative Coordinate Transformation: [XT-XB; XB; eta]
% Ngen = Ngens(setid) + length(K)-Ngens(setid)-MESH.Nn*3*2;
Ngen = length(K) - MESH.Nn*3*2;
Trel = [speye(MESH.Nn*3), speye(MESH.Nn*3), sparse(MESH.Nn*3, Ngen);
        sparse(MESH.Nn*3, MESH.Nn*3), speye(MESH.Nn*3), sparse(MESH.Nn*3, Ngen);
        sparse(Ngen, MESH.Nn*3*2), speye(Ngen)];
Krel = Trel'*K*Trel;
Krel = 0.5*(Krel+Krel');
Mrel = Trel'*M*Trel;
Mrel = 0.5*(Mrel+Mrel');
Fvrel = Trel'*Fv;
Rrel = R*Trel;

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
ni = 1:(MESH.Nn*3);
mi = setdiff(1:length(Krel), ni);
cnum = 1e10;

Ka = sparse([Krel(ni, ni)+cnum*NTN Krel(ni, mi) -cnum*NTG;
             Krel(mi, ni), Krel(mi, mi), zeros(length(mi),Npatches*6);
             -cnum*NTG' zeros(Npatches*6, length(mi)) cnum*GTG]);
Ka = 0.5*(Ka+Ka');
Ma = sparse(blkdiag(Mrel, zeros(Npatches*6)));
Ma = 0.5*(Ma+Ma');
Fva = [Fvrel; zeros(Npatches*6,1)];
LamTa = [NTG*inv(GTG); zeros(length(mi)+Npatches*6,Npatches*6)];
Ftmpa = Ka*sparse([zeros(length(Krel), 1); cnum^(-1)*NTG\Fvrel(ni)]);

Ra = sparse([Rrel zeros(size(Rrel,1), Npatches*6)]);

dofred = 6;
vdofs = length(K) + reshape(1:(Npatches*6), 6, Npatches);
vdofs = reshape(vdofs(1:dofred,:), Npatches*dofred, 1);  % virtual DOFs

[Mhcb, Khcb, Thcb] = HCBREDUCE(Ma, Ka, vdofs, Ngens(setid));
Mhcb = 0.5*(Mhcb+Mhcb');  Khcb = 0.5*(Khcb+Khcb');
Fvhcb = Thcb'*Fva;
Rhcb = Ra*Thcb;
Lamhcb = Thcb'*LamTa;
disp('DONE!')

%% NULL-SPACE REDUCTION
[Vhcb, Dhcb] = eig(full(Khcb(Npatches*3+1:end, Npatches*3+1:end)), full(Mhcb(Npatches*3+1:end, Npatches*3+1:end))); 
% [Vhcb, Dhcb] = eig(full(Khcb), full(Mhcb)); 
[Dhcb, si] = sort(diag(Dhcb));  Vhcb = Vhcb(:, si);
% Vrbms = Vhcb(:, 1:6);
Vrbms = [zeros(Npatches*3, 6); Vhcb(:, 1:6)];  Vrbms = Vrbms./sqrt(diag(Vrbms'*Mhcb*Vrbms)');
L = null(Vrbms'*Mhcb);

M = L'*Mhcb*L;
M = 0.5*(M+M');
K = L'*Khcb*L;
K = 0.5*(K+K');
Fv = L'*Fvhcb;
R = Rhcb*L;
Th = Thcb*L;
LamT = L'*Lamhcb;
[sort(eig(K), 'descend') sort(eig(M), 'descend')]
disp(cnum)
%% Save
fname = sprintf('./REDMATS/%d_SET_WJMAT.mat', setid);
save(fname, 'M', 'K', 'L', 'Fv', 'R', 'Th', 'LamT', 'MESH', 'PatchAreas', 'Npatches', 'dofred', 'NTN', 'GTG', 'NTG', 'cnum')
disp('SAVED!')