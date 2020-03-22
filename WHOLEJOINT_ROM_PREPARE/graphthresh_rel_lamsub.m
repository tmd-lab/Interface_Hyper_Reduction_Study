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

%% Thresholded Patch - Selection
Prestress = 11580; sint = 1e-6;
dat = load(sprintf('./DATS/statsol_P%d_S%.2f.mat',Prestress,log10(sint)), 'Xs', 'Qm', 'Tm', 'Tx', 'Ty', 'Tn', 'dR0');

[Q1, T1] = ZTE_ND2QP(MESH, 1);
Pobj = Q1*(dat.Qm\dat.Tn);

% Pobj = (Pobj-min(Pobj))/range(Pobj);

% 2(3), 4(11), 5(13), 7(21)
Nlevs = 2;  % 2(3), 3(5), 5(9), 8(15)
Plevels = linspace(min(Pobj), max(Pobj), Nlevs+1); Plevels(1)=[];
LevelIs = cell(Nlevs,1);
AdjMxs = cell(Nlevs,1);
G = cell(Nlevs,1);
LevelNDIDS = cell(Nlevs,1);
LevelELIDS = cell(Nlevs,1);
Nsets = 0;
Npatches = 0;
for i=1:Nlevs
    LevelIs{i} = find(Pobj<=Plevels(i));
    for j=1:(i-1)
        LevelIs{i} = setdiff(LevelIs{i},LevelIs{j});
    end
    AdjMxs{i} = NODEELADJ(MESH.Nds, MESH.Quad(LevelIs{i},2:end));
    G{i} = graph(AdjMxs{i});
    [bins,bincounts] = conncomp(G{i});
    binids = find(bincounts~=1);
    Nsets = length(binids);
    LevelNDIDS{i} = cell(Nsets,1);
    LevelELIDS{i} = cell(Nsets,1);
    for j=1:Nsets
        LevelNDIDS{i}{j} = find(bins==binids(j));
        LevelELIDS{i}{j} = find(sum(ismember(MESH.Quad(:, 2:end), LevelNDIDS{i}{j}), 2)==4);
    end
    
    Npatches = Npatches + Nsets;
end
Pnds = cell(Npatches,1);
Pels = cell(Npatches,1);
Plevs = zeros(Npatches,1);
k=1;
for i=1:Nlevs
    for j=1:length(LevelNDIDS{i})
        Pnds{k} = LevelNDIDS{i}{j};
        Pels{k} = LevelELIDS{i}{j};
        Plevs(k) = Plevels(i);
        
        if length(Pels{k})==1
            error('Inappropriate')
        end
        k = k+1;
    end
end
% %% PATCHING
Area = cell(1, Npatches);
qps = cell(1, Npatches);
ctrds = zeros(Npatches, 2);

figure(1)
clf()
for n=1:Npatches
    Area{n} = T1(Pnds{n}, Pels{n})*ones(length(Pels{n}), 1);
    qps{n} = Q1(Pels{n}, :)*MESH.Nds;
    ctrds(n, :) = sum(Area{n}.*MESH.Nds(Pnds{n}, :))/sum(Area{n});
    
    k=n;SHOW2DMESH(MESH.Nds, [], MESH.Quad(Pels{k},:), Plevs(k), -1, -100)
end
PatchAreas = cellfun(@(c) sum(c), Area);
axis equal; axis off; colormap(jet(Npatches));
title(sprintf('%d Patches',Npatches))
set(gca, 'Position', [-0.5 -0.1 2 1])
print(sprintf('./FIGS/P%d_S%.2f_%dLEV_%dPATCH_MESH.eps',Prestress, log10(sint), Nlevs,Npatches), '-depsc')
pause(1)
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
cnum = 1e11;

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
% [Vhcb, Dhcb] = eig(full(Khcb(Npatches*3+1:end, Npatches*3+1:end)), full(Mhcb(Npatches*3+1:end, Npatches*3+1:end))); 
% [Dhcb, si] = sort(diag(Dhcb));  Vhcb = Vhcb(:, si);
% Vrbms = [zeros(Npatches*3, 6); Vhcb(:, 1:6)];  Vrbms = Vrbms./sqrt(diag(Vrbms'*Mhcb*Vrbms)');
% L = null(Vrbms'*Mhcb);

[Vhcb, Dhcb] = eig(full(Khcb(Npatches*6+1:end, Npatches*6+1:end)), full(Mhcb(Npatches*6+1:end, Npatches*6+1:end)));
[Dhcb, si] = sort(diag(Dhcb));  Vhcb = Vhcb(:, si);
Vrbms = [zeros(Npatches*6, 6); Vhcb(:, 1:6)];  Vrbms = Vrbms./sqrt(diag(Vrbms'*Mhcb*Vrbms)');
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

%% Check
Dred = eigs(K, M, 10, 'SM');
Drel = eigs(Krel, Mrel, 16, 'SM');
[Dred Drel(7:end)]
%% Save
fname = sprintf('./REDMATS/P%d_S%.2f_%dLEV_GRED_WJMAT.mat', Prestress, log10(sint), Nlevs);
% save(fname, 'M', 'K', 'L', 'Fv', 'R', 'Th', 'LamT', 'MESH', 'PatchAreas', 'Npatches', 'dofred', 'NTN', 'GTG', 'NTG', 'cnum', 'Pels', 'Pnds')
disp('SAVED!')