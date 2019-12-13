clc
clear all

MEXPATH = '../MATRIX_EXTRACTION/RUNS/';
SETDIRS = {'1_AROUNDSET', '2_ABOVESET', '3_SINGELEMABOVESET', '4_INTSET', '5_INTSETNPS'};

setid = 5;  % To call SETDIRS

% ROM levels
plevs = 52;  % 52, 124, 300
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

%% Load Reduced Mesh
fname = sprintf('./Tobias_Mats/PREP_20Nm_PMESHZT_%d.mat',plevs);
red = load(fname, 'MESH');
red.MESH.Nds = (red.MESH.Nds-min(red.MESH.Nds)).*range(MESH.Nds)./range(red.MESH.Nds)+min(MESH.Nds);

%% INTERPOLATION MATRICES
Qn = MESH2PTS(red.MESH, MESH.Nds);
P2F = MESH2D_P2FMAT(MESH);
P2F_n = MESH2D_P2FMAT(red.MESH);

Tmat = Qn;
Tmatf = Qn'*P2F*Qn*inv(P2F_n);

clf(); SHOW2DMESH(red.MESH.Nds, red.MESH.Tri, red.MESH.Quad, 1, -1, -100);
plot(MESH.Nds(:,1), MESH.Nds(:,2), 'k.');
plot(Qn*red.MESH.Nds(:,1), Qn*red.MESH.Nds(:,2), 'ko'); axis equal; axis off
%% ROM preparation
TFM = blkdiag(kron(Tmat,eye(3)), eye(Nrest+Nint*3));
TFMf = blkdiag(kron(Tmatf,eye(3)));

Kred = TFM'*Krel*TFM;  Kred = 0.5*(Kred+Kred');
Mred = TFM'*Mrel*TFM;  Mred = 0.5*(Mred+Mred');
Rred = Rrel*TFM;
Fvred = TFM'*Fvrel;
%% HCB Here
[Mhcb, Khcb, Thcb] = HCBREDUCE(Mred,Kred,1:red.MESH.Nn*3,Nrest);
Mhcb = 0.5*(Mhcb+Mhcb');  Khcb = 0.5*(Khcb+Khcb');
Rhcb = Rred*Thcb;
Fvhcb = Thcb'*Fvred;
TFMfhcb = [Thcb(1:red.MESH.Nn*3,:)'*TFMf, zeros(size(Mhcb,1),Nrest)];

%% Null Space Reduction
[V, D] = eigs(Khcb(red.MESH.Nn*3+1:end,red.MESH.Nn*3+1:end), Mhcb(red.MESH.Nn*3+1:end,red.MESH.Nn*3+1:end), 10, 'SM');
[D,si] = sort(sqrt(abs(diag(D)))/(2*pi));
nrbms = 6;
Vrbm = [zeros(red.MESH.Nn*3,nrbms); V(:,1:nrbms)];  Vrbm = Vrbm./sqrt(diag(Vrbm'*Mhcb*Vrbm)');
L = null(Vrbm'*Mhcb);

M = L'*Mhcb*L;  M = 0.5*(M+M');
K = L'*Khcb*L;  K = 0.5*(K+K');
R = Rhcb*L;
Fv = L'*Fvhcb;
TFMf = L'*TFMfhcb;
MESH = red.MESH;

save(sprintf('ROM_20Nm_PMESHZT_%d.mat',plevs), 'M', 'K', 'R', 'L', 'Fv', 'TFMf', 'TFMfhcb', 'MESH');
disp('Done!')