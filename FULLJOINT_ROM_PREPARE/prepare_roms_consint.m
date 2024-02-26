clc
clear all
addpath('../ROUTINES/')
addpath('../ROUTINES/FEM/')
addpath('../ROUTINES/GRAPH/')

MEXPATH = '../MATRIX_EXTRACTION/RUNS/';
SETDIRS = {'1_AROUNDSET', '2_ABOVESET', '3_SINGELEMABOVESET', '4_INTSET', ...
            '5_INTSETNPS', '6_HBRB_Baseline', '7_HBRB_MoreModes'};

setid = 6;  % To call SETDIRS

% ROM levels
% sel_method = 'P';
% Nels = 304;  % [52 100 140 204 240 304]

% sel_method = 'U';
% Nels = 588; % [36 74 122 232 448 568 588]

% sel_method = 'PD';
% Nels = 152;  % [48 68 108 152 200 252 292]

%%%%%%% Other Settings:
% Eliminate singular values lower than this during null space in HCB
% Specifically if S(i) / S(1) < HCB_null_space_tol
% Look at HCBREDUCE.m for details
HCB_null_space_tol = 1e-10; 
Ncomp_final = 20; % Number of fixed interface modes at final step


sel_method = 'U';
for Nels = [232]
% % FRICTMODEL STUDY
% sel_method = 'PD';
% Nels = 68;
%% Mesh Extract
% MESH.Nds     = dlmread('../MATRIX_EXTRACTION/RUNS/MESHPROCESS/Nodes.dat');
% MESH.Quad    = dlmread('../MATRIX_EXTRACTION/RUNS/MESHPROCESS/Elements.dat');

Nodes = dlmread([MEXPATH SETDIRS{setid} '/Nodes.dat']);

MESH.Nds     = Nodes(:, 1:2); % Need to eliminate potential 3rd column of zcoordinates for this code to work. 
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
Fv = sparse(Fv);
R = sparse(R);
M = sparse(M); M = 0.5*(M+M');
K = sparse(K); K = 0.5*(K+K');
Nrest = size(M, 1)-Nint*3*2;
disp('Matrices Extracted.');

tic
%% Relative Transformation : [Xt-Xb; Xb; Xi..]

% Goal: Xabaqus = Trel * [Xt-Xb; Xb; Xi..]  

% BRB Outputs: Xabaqus = [Xt; Xb; Xi..] = Trel * [Xt-Xb; Xb; Xi..]
Trel = sparse([eye(Nint*3),  eye(Nint*3), zeros(Nint*3, Nrest);
               zeros(Nint*3), eye(Nint*3), zeros(Nint*3, Nrest);
               zeros(Nrest, Nint*3*2),     eye(Nrest)]);

% disp('HBRB outputs Surfaces in different order')
% if setid >= 6 
%     % HBRB Outputs: Xabaqus = [Xb; Xt; Xi..] = Trel * [Xt-Xb; Xb; Xi..]
%     Trel = sparse([zeros(Nint*3),  eye(Nint*3), zeros(Nint*3, Nrest);
%                    eye(Nint*3),    eye(Nint*3), zeros(Nint*3, Nrest);
%                    zeros(Nrest, Nint*3*2),      eye(Nrest)]);
% end

Mrel  = Trel'*M*Trel; Mrel = 0.5*(Mrel+Mrel');
Krel  = Trel'*K*Trel; Krel = 0.5*(Krel+Krel');
Rrel  = R*Trel;
Fvrel = Trel'*Fv;
disp('Relative Transformation Done.')

%% Load Reduced Mesh
% fname = sprintf('./Tobias_Mats/PREP_20Nm_PMESHZT_%d.mat',Nels);
fname = sprintf('./MATS/REDMESH_%s_%dELS.mat', sel_method, Nels);
red = load(fname, 'MESH');
red.MESH.Nds = (red.MESH.Nds-min(red.MESH.Nds)).*range(MESH.Nds)./range(red.MESH.Nds)+min(MESH.Nds);

%% MESH Plotting
figure; 
plot(MESH.Nds(:, 1), MESH.Nds(:, 2), 'o'); 
pbaspect([range(MESH.Nds(:, 1)),range(MESH.Nds(:, 2)),1]); 
hold on; 
plot(red.MESH.Nds(:, 1), red.MESH.Nds(:, 2), 'r+');


% Plot Mesh Elements
color_opts = {'b', 'r'};
mesh_opts = {MESH, red.MESH};
% mesh_opts = {MESH};

% for j = 1:length(mesh_opts)
% 
%     figure; 
%     hold on;
%     cMesh = mesh_opts{j};
%     lc = color_opts{j};
% 
%     for i = 1:size(cMesh.Quad, 1)
% 
%         node_inds = [cMesh.Quad(i, 2:end), cMesh.Quad(i, 2)];
% 
%         plot(cMesh.Nds(node_inds, 1), cMesh.Nds(node_inds, 2), lc)
%         fill(cMesh.Nds(node_inds, 1), cMesh.Nds(node_inds, 2), lc)
% 
% %         drawnow; 
%     end
%     pbaspect([range(cMesh.Nds(:, 1)),range(cMesh.Nds(:, 2)),1]); 
% end


%% INTERPOLATION MATRICES
% Qn = MESH2PTS(red.MESH, MESH.Nds);
Qn = MESHTFM(MESH, red.MESH, 4);
% P2F = MESH2D_P2FMAT(MESH);
% P2F_n = MESH2D_P2FMAT(red.MESH);
[Q,T] = ZTE_ND2QP(MESH, 2);  P2F = T*Q;
[Q,T] = ZTE_ND2QP(red.MESH, 2);  P2F_n = T*Q;

Tmat = Qn;
Tmatf = Qn'*P2F*Qn*inv(P2F_n);

% figure(1);clf(); 
% SHOW2DMESH(red.MESH.Nds+[0 0.02], red.MESH.Tri, red.MESH.Quad, 1, -1, -100);
% SHOW2DMESH(Qn*red.MESH.Nds-[0 0.02], MESH.Tri, MESH.Quad, 1, -1, -100); 
% axis equal; axis off; pause(0.01)
%% ROM preparation
TFM = blkdiag(kron(Tmat,eye(3)), eye(Nrest+Nint*3));
TFMf = blkdiag(kron(Tmatf,eye(3)));

Kred = TFM'*Krel*TFM;  Kred = 0.5*(Kred+Kred');
Mred = TFM'*Mrel*TFM;  Mred = 0.5*(Mred+Mred');
Rred = Rrel*TFM;
Fvred = TFM'*Fvrel;

disp('ROM prepared.')

disp('Stuck Interface, Relative coordinates - 6 RBMs')
[V, D1] = eig(full(Krel(Nint*3+1:end,Nint*3+1:end)), full(Mrel(Nint*3+1:end,Nint*3+1:end)));
A = diag(D1); sqrt(A(1:10))/2/pi

% Call that HBM Should be Doing for stuck:
[V, D2] = eigs((Krel(Nint*3+1:end,Nint*3+1:end)), (Mrel(Nint*3+1:end,Nint*3+1:end)), 2*3*Nrest, 'SM');
A = diag(D2); sqrt(A(1:10))/2/pi

%% Create an HCB of the Relative Coordinates - Debugging

Nb_rom = 3*Nrest;
Nrbm = 6;

% Adding an extra argument here adds a large negative eigenvalue in the HCB
% solution, but then the null space may be able to correctly eliminate 6
% RBMs and leave 1 for the separating (null space code here eliminates
% negative so is not perfect, but close).
% [MhcbRel, KhcbRel, ThcbRel] = HCBREDUCE(Mrel,Krel,1:3*Nint,Nb_rom, true); 
% MhcbRel = 0.5*(MhcbRel+MhcbRel');  KhcbRel = 0.5*(KhcbRel+KhcbRel');
% 
% disp('Verify Stuck modes:')
% [V, D] = eigs(KhcbRel(Nint*3+1:end,Nint*3+1:end), MhcbRel(Nint*3+1:end,Nint*3+1:end), 10, 'SM');
% sqrt(diag(D))/2/pi
% 
% [V, D] = eig(KhcbRel(Nint*3+1:end,Nint*3+1:end), MhcbRel(Nint*3+1:end,Nint*3+1:end));
% A = diag(D); sqrt(A(1:10))/2/pi
% 
% 
% disp('Verify full modes:')
% [V, D] = eigs(KhcbRel, MhcbRel, 10, 'SM');
% sqrt(diag(D))/2/pi
% 
% disp('Verify Null of Full Modes')
% 
% disp('VERY HACKY NULL SPACE!')
% Linds = [1:3*Nint, 3*Nint+((Nrbm+1):Nb_rom)];
% 
% MRel_HCB_Null = MhcbRel(Linds, Linds);
% KRel_HCB_Null = KhcbRel(Linds, Linds);
% 
% MRel_HCB_Null = 0.5*(MRel_HCB_Null+MRel_HCB_Null');  
% KRel_HCB_Null = 0.5*(KRel_HCB_Null+KRel_HCB_Null');
% 
% % This shows 4 zero energy modes, when it should show one based on my
% % understanding.
% [V, D] = eigs(KRel_HCB_Null, MRel_HCB_Null, 10, 'SM');
% sqrt(diag(D))/2/pi
% 
% disp('Full Modes of Null Reduced: - THESE HAVE CHANGED AND NOW LOOK WRONG!')
% [V, D] = eig(KRel_HCB_Null, MRel_HCB_Null);
% A = diag(D); sqrt(A(1:10))/2/pi
% 
% disp('Verify fixed modes of Null model') % These are good. - no zero eigenvalues
% [V, D] = eigs(KRel_HCB_Null(Nint*3+1:end,Nint*3+1:end), MRel_HCB_Null(Nint*3+1:end,Nint*3+1:end), 10, 'SM');
% sqrt(diag(D))/2/pi
% 
% disp('Fundamental problem, KRel_HCB is not P.S.D')
% [V, D] = eig(KRel_HCB_Null);
% A = diag(D); sqrt(A(1:10))/2/pi
% 
% 
% disp('Double Check that M is P.D.')
% [V, D] = eig(MRel_HCB_Null);
% A = diag(D); sqrt(A(1:10))/2/pi

%% HCB Here
[Mhcb, Khcb, Thcb] = HCBREDUCE(Mred,Kred,1:red.MESH.Nn*3,...
                                Ncomp_final, HCB_null_space_tol);
Mhcb = 0.5*(Mhcb+Mhcb');  Khcb = 0.5*(Khcb+Khcb');
Rhcb = Rred*Thcb;
Fvhcb = Thcb'*Fvred;
TFMfhcb = [Thcb(1:red.MESH.Nn*3,:)'*TFMf, zeros(size(Mhcb,1),3*Nrest)];
TFMhcb = TFM*Thcb;
%% Null Space Reduction
[V, D] = eigs(Khcb(red.MESH.Nn*3+1:end,red.MESH.Nn*3+1:end), Mhcb(red.MESH.Nn*3+1:end,red.MESH.Nn*3+1:end), 10, 'SM');
[D,si] = sort(sqrt(abs(diag(D)))/(2*pi));
V = V(:, si);
nrbms = 6;
Vrbm = [zeros(red.MESH.Nn*3,nrbms); V(:,1:nrbms)];  Vrbm = Vrbm./sqrt(diag(Vrbm'*Mhcb*Vrbm)');
L = null(Vrbm'*Mhcb);

% disp(['Eliminating fixed interface modes with modal frequencies of :', mat2str(sqrt(D)/2/pi)])

% [U, S, V] = svd(L, 'econ');
% L = U;

% %%%%%%%%%%%%%%%%%%%%%%%%%
% % Alternative L with Gram-Schmidt 
% disp('Using alternative L, this has not been well tested.')
% tmp = eye(size(V, 1));
% ModeToOrtho = [V(:, 1:nrbms), tmp(:, nrbms+1:end)];
% ModeToOrtho = [zeros(3*red.MESH.Nn, size(ModeToOrtho, 2)); ModeToOrtho];
% 
% for col_curr = 2:size(ModeToOrtho, 2)
%     for col_sub = 1:(col_curr-1)
%         tmp = (ModeToOrtho(:, col_sub)'*Mhcb* ModeToOrtho(:, col_curr)) / (ModeToOrtho(:, col_sub)'*Mhcb* ModeToOrtho(:, col_sub));
%         ModeToOrtho(:, col_curr) = ModeToOrtho(:, col_curr) - tmp * ModeToOrtho(:, col_sub);
%     end
% end
% 
% ModeNull = ModeToOrtho(:, nrbms+1:end);
% 
% L = [[eye(red.MESH.Nn*3) ;
%      zeros(size(ModeNull, 1)- 3*red.MESH.Nn, red.MESH.Nn*3)], ModeToOrtho];
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


M = L'*Mhcb*L;  M = 0.5*(M+M'); 
K = L'*Khcb*L;  K = 0.5*(K+K');
R = Rhcb*L;
Fv = L'*Fvhcb;
TFMh = TFMhcb*L;
TFMf = L'*TFMfhcb;
MESH = red.MESH;
ttk = toc

% save(sprintf('ROM_20Nm_PMESHZT_%d.mat',Nels), 'M', 'K', 'R', 'L', 'Fv', 'TFMf', 'TFMfhcb', 'MESH');
save(sprintf('./ROMS/ROM_%s_%dELS.mat',sel_method,Nels), 'M', 'K', 'R', 'L', 'Fv', 'TFMh', 'TFMf', 'TFMfhcb', 'MESH', 'ttk');
disp('Done!')
end