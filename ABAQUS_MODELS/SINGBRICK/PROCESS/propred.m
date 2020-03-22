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
M1 = sparse(blkdiag(M, zeros(12)));
K1 = sparse([K(1:12,:) zeros(12, 6) -NTG;
      K(13:end,:) zeros(12, 12);
      zeros(6, 30) GTG;
      -NTG' zeros(6, 12) GTG zeros(6)]);
% HCB Here
[Mh1, Kh1, Th] = HCBREDUCE(M1, K1, 25:26, 10);  

%% ROM Development (Lagrange Multipliers Eliminated)
M2 = blkdiag(M, zeros(6));
K2 = [K(1:12,1:12) K(1:12,13:end)*(eye(12)-NTG*NTGp) K(1:12,1:12)*NTGp'*GTG;
      (eye(12)-NTGp'*NTG')*K(13:end,1:12) K(13:end,13:end) K(13:end,1:12)*NTGp'*GTG;
      GTG*NTGp*K(1:12,1:12) GTG*NTGp*K(1:12,13:end) zeros(6)];
  
T2 = sparse([eye(30);
        NTGp*K(1:12,:) zeros(6)]);
K2 = T2'*K1*T2;
K2 = 0.5*(K2+K2');
% HCB Here
[Mh2, Kh2, Th] = HCBREDUCE(M2, K2, 25:30, 10);
[eigs(K, M, 10, 'SM') eigs(Kh2, Mh2, 10, 'SM') eigs(Kh1, Mh1, 10, 'SM')]