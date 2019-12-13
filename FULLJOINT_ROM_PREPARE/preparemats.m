clc
clear all

MEXPATH = '../MATRIX_EXTRACTION/RUNS/';
SETDIRS = {'1_AROUNDSET', '2_ABOVESET', '3_SINGELEMABOVESET', '4_INTSET', '5_INTSETNPS'};

setid = 5;  % To call SETDIRS
%% Mesh Extract
% MESH.Nds     = dlmread('../MATRIX_EXTRACTION/RUNS/MESHPROCESS/Nodes.dat');
% MESH.Quad    = dlmread('../MATRIX_EXTRACTION/RUNS/MESHPROCESS/Elements.dat');
MESH.Nds     = dlmread([MEXPATH SETDIRS{setid} '/Nodes.dat']);
MESH.Quad    = dlmread([MEXPATH SETDIRS{setid} '/Elements.dat']);
MESH.Tri     = [];
MESH.Nn      = size(MESH.Nds, 1);
MESH.Ne_Quad = size(MESH.Quad, 1);
MESH.Ne_Tri  = size(MESH.Tri, 1);
MESH.Ne      = MESH.Ne_Quad + MESH.Ne_Tri;

Nint = MESH.Nn;  % Number of Interface nodes (first Nint)
disp('MESH Extracted')

%% Load Mat Files
load([MEXPATH SETDIRS{setid} '/BRB_WOPRES_MAT.mat'], 'M', 'K', 'R', ...
     'Fv');
Fv = sparse(Fv');
R = sparse(R');
M = sparse(M); M = 0.5*(M+M');
K = sparse(K); K = 0.5*(K+K');
Nrest = size(M, 1)-Nint*3*2;
disp('Matrices Extracted.');

%% Relative Transformation : [Xt-Xb; Xb; Xi..]
Trel = sparse([eye(Nint*3),  eye(Nint*3), zeros(Nint*3, Nrest);
               zeros(Nint*3), eye(Nint*3), zeros(Nint*3, Nrest);
               zeros(Nrest, Nint*3*2),     eye(Nrest)]);
Mrel  = Trel'*M*Trel; Mrel = 0.5*(Mrel+Mrel');
Krel  = Trel'*K*Trel; Krel = 0.5*(Krel+Krel');
Rrel  = R*Trel;
Fvrel = Trel'*Fv;
disp('Relative Transformation Done.')

%% Null-Space Transformation of total motion
[V, D] = eigs(Krel(MESH.Nn*3+1:end, MESH.Nn*3+1:end), Mrel(MESH.Nn*3+1:end, MESH.Nn*3+1:end), 20, 'SM');
[D, si] = sort(sqrt(abs(diag(D)))/(2*pi));  V = V(:, si);
nrbm = 6;
Vrbm = [zeros(MESH.Nn*3,nrbm); V(:, 1:nrbm)]; Vrbm = Vrbm./sqrt(diag(Vrbm'*Mrel*Vrbm)');
L = null(Vrbm'*Mrel);

%%
M = L'*Mrel*L; M = 0.5*(M+M');
K = L'*Krel*L; K = 0.5*(K+K');
R = Rrel*L;
Fv = L'*Fvrel;
disp('Rigid Body Modes Null-Space Transformed out.')
    
save(sprintf('./MATS/%d_SET_NULLRED.mat', setid), 'M', 'K', 'R', 'L', 'Fv', 'MESH', 'Krel', 'Fvrel');
% save('tmp.mat', 'M', 'K', 'R', 'L', 'Fv', 'MESH')