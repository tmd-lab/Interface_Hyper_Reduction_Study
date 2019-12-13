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
load(fname, 'MESH', 'L')
[Q1, T1] = ZTE_ND2QP(MESH, 1);
EAreas = sum(T1,1);

wd = 0;  wp = 1;
sn = 10;
%% Data Import
Prestress = 11580; sint = 1e-6;
dat = load(sprintf('./DATS/statsol_P%d_S%.2f.mat',Prestress,log10(sint)), ...
           'Xs', 'Qm', 'Tm', 'Tx', 'Ty', 'Tn', 'dR0', 'V');
%% Objective Selection
[Q1, T1] = ZTE_ND2QP(MESH, 1);
PObj = dat.Qm\dat.Tn;  % Normal traction
DObj = abs(dat.V(3:3:MESH.Nn*3,1)-dat.V(MESH.Nn*3+(3:3:MESH.Nn*3),1));  % Relative Displacement

PObj = (PObj-min(PObj))/range(PObj);
DObj = (DObj-min(DObj))/range(DObj);
Obj = (wd*DObj+wp*PObj)/(wd+wp);

%% Objective Smoothening
AdjMx = NODEELADJ(MESH.Nds, MESH.Quad(:,2:end)); AdjMx = AdjMx./sum(AdjMx,2);
AdjMx = AdjMx + eye(MESH.Nn); AdjMx = AdjMx./sum(AdjMx,2);
Obj = AdjMx^sn*Obj;

Obje = Q1*Obj;  % Element-wise objectives

%% Parameters
Nlevs = 3;  % [3 , 5 , 7, 9, 11, 19, 21]
sel_method = 'GL';

%% Construct and Solve Graph Laplacian Problem
[Q2, T2] = ZTE_ND2QP(MESH, 2);
[Q2x, Q2y, T2x, T2y] = ZTE_ND2QPD(MESH, 2);

M = T2*Q2;  M = 0.5*(M+M');
K = T2x*Q2x + T2y*Q2y;  K = 0.5*(K+K');

% Eigendecomposition
[V, D] = eig(full(K), full(M));
[D, si] = sort(diag(D));
V = V(:, si);

figure(1)
clf()
semilogy(abs(D), '.')
%%
figure(2)
clf()
for i=1:4*3
    subplot(4,3,i)
    SHOW2DMESH(MESH.Nds, MESH.Tri, MESH.Quad, V(:, i+2), -1, -100);
    axis equal; colormap(jet(255)); axis off
    title(sprintf('Mode %d', i))
end

%% Embedding
figure(3)
clf()
k1 = 2;
k2 = 3;
k3 = 4;
plot3(V(:,k1), V(:,k2), V(:,k3), 'ko', 'MarkerFaceColor', 'k')

%% K-Means Clustering
rng(1)
Nk = 140;
% idx = kmeans([Obj V(:,(end-(Nk-1)):end)], Nk);
% idx = kmeans(V(:,(end-(Nk-1)):end), Nk);

idx1 = kmeans(V(:,(end-(Nk-1)):end), Nk);
idx2 = kmeans(V(:,2:(Nk+1)), Nk);

% idx1 = kmeans(V(:,2:(Nk+1)), Nk);
% idx2 = kmeans(V(:,2:10), Nk);

% idx = kmeans([dat.V(1:3:MESH.Nn*3,2:10)], Nk);

colos = DISTINGUISHABLE_COLORS(Nk);

figure(4)
clf()
% plot(MESH.Nds(:,1), MESH.Nds(:,2), 'k.'); hold on
for i=1:Nk
    idi = find(idx1==i);
    plot(MESH.Nds(idi,1), MESH.Nds(idi,2)+0.03, 'o', 'Color', ...
         colos(i,:), 'MarkerFaceColor', colos(i,:)); hold on
    
    idi = find(idx2==i);
    plot(MESH.Nds(idi,1), MESH.Nds(idi,2)-0.03, 'o', 'Color', ...
         colos(i,:), 'MarkerFaceColor', colos(i,:))
    % chi = idi(convhull(MESH.Nds(idi,:)));
    % fill(MESH.Nds(chi,1), MESH.Nds(chi,2), colos(i,:));
end
axis equal; axis off