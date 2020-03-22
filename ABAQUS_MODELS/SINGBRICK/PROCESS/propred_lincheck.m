clc
clear all
addpath('../../../ROUTINES/')
addpath('../../../ROUTINES/FEM/')
load('./MATS.mat', 'K', 'M', 'Kc', 'Mc')

%% MESH STRUCTURE
MESH.Nds = [-1 -1; 1 -1; 1 1; -1 1];
MESH.Quad = [1 1 2 3 4];
MESH.Tri  = [];
MESH.Nn      = size(MESH.Nds, 1);
MESH.Ne_Quad = size(MESH.Quad, 1);
MESH.Ne_Tri  = size(MESH.Tri, 1);
MESH.Ne      = MESH.Ne_Quad + MESH.Ne_Tri;

[Q1, T1] = ZTE_ND2QP(MESH, 1);

%% Patch Reduction
Area = sum(full(T1));
qps = Q1*MESH.Nds;
ctrd = full(sum(T1.*MESH.Nds)/Area);

[P, Nums, NTN, NTG, GTG] = CONSPATCHMAT(MESH.Nds, [], MESH.Quad, ctrd);
NTGp = (NTG'*NTG)\NTG';  % Pseudo inverse
%% ROM Development (With Lagrange Multipliers)
cnum = 1e6;
M1 = sparse(blkdiag(M, zeros(12)));
K1 = sparse([K(1:12,:) zeros(12, 6) -cnum*NTG;
      K(13:end,:) zeros(12, 12);
      zeros(6, 30) cnum*GTG;
      -cnum*NTG' zeros(6, 12) cnum*GTG zeros(6)]);
disp(cond(K1))  
% HCB Here
[Mh1, Kh1, Th] = HCBREDUCE(M1, K1, 25:30, 18);  

%% Modes
[V, D] = eig(full(K), full(M));
[D, si] = sort(diag(D));
V = V(:, si);
V = V./sqrt(diag(V'*M*V)');

[V1, D1] = eig(full(K1), full(M1));
[D1, si] = sort(diag(D1));
V1 = V1(:, si);
V1 = V1./sqrt(diag(V1'*M1*V1)');
m1s = find(isfinite(D1));

[Vh1, Dh1] = eig(full(Kh1), full(Mh1));
[Dh1, si] = sort(diag(Dh1));
Vh1 = Vh1(:, si);
Vh1 = Vh1./sqrt(diag(Vh1'*Mh1*Vh1)');

Vh = Th*Vh1;
mhs = find(isfinite(Dh1));

figure(1)
br=bar3((V1(1:24, m1s)'*M*V).^2./(diag(V1(1:24,m1s)'*M*V1(1:24,m1s)).*diag(V'*M*V)'));
for k=1:size(br,2); br(k).CData = br(k).ZData; br(k).FaceColor = 'interp'; end

figure(2)
br=bar3((Vh(1:24, mhs)'*M*V).^2./(diag(Vh(1:24,mhs)'*M*Vh(1:24,mhs)).*diag(V'*M*V)'));
for k=1:size(br,2); br(k).CData = br(k).ZData; br(k).FaceColor = 'interp'; end
