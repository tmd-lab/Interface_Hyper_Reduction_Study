clc
clear all
addpath('../ROUTINES/FEM/') % For HCB

MEXPATH = '../MATRIX_EXTRACTION/RUNS/';
SETDIRS = {'1_AROUNDSET', '2_ABOVESET', '3_SINGELEMABOVESET', '4_INTSET', ...
            '5_INTSETNPS', '6_HBRB_Baseline', '7_HBRB_MoreModes'};

setid = 7;  % To call SETDIRS
Ncomp_final = 100; % Number of fixed interface modes at final step
HCB_null_space_tol = 1e-10; 

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
Fv = sparse(Fv);
R = sparse(R);
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

%% HCB To Eliminate Xb
eliminateXb = true;
if(eliminateXb)    
    
    n_modesCompare = 19;
    
    %%%%%%%% Without fixing contact
    [V_initial, D_initial] = eigs(Krel, Mrel, n_modesCompare, 'SM'); % For comparisons / verification
    D_initial = sqrt(sort(diag(D_initial)))/2/pi;
    
    %%%%%%%% with fixing contact
    [V_initial, D_initial2] = eigs(Krel(3*MESH.Nn+1:end, 3*MESH.Nn+1:end), ...
                                    Mrel(3*MESH.Nn+1:end, 3*MESH.Nn+1:end), n_modesCompare, 'SM'); % For comparisons / verification
    D_initial2 = sqrt(sort(diag(D_initial2)))/2/pi;
    
    %%%%%%%% Do the Reduction
    [Mhcb, Khcb, Thcb] = HCBREDUCE(Mrel,Krel,1:MESH.Nn*3,Ncomp_final, HCB_null_space_tol);
    Mrel = 0.5*(Mhcb+Mhcb');  Krel = 0.5*(Khcb+Khcb');
    Rrel = Rrel*Thcb;
    Fvrel = Thcb'*Fvrel;
end

%% Verify Frequencies
if(eliminateXb)

    % Without fixing contact
    [V_final, D_final] = eigs(Khcb, Mhcb, n_modesCompare, 'SM'); % For comparisons / verification
    D_final = sqrt(sort(diag(D_final)))/2/pi;

    Diffs = D_initial(8:end) - D_final(8:end);
    log10( abs(Diffs ./ D_initial(8:end) ))'
    
    % with fixing contact
    [V_final, D_final2] = eigs(Khcb(3*MESH.Nn+1:end, 3*MESH.Nn+1:end), ...
                                Mhcb(3*MESH.Nn+1:end, 3*MESH.Nn+1:end), n_modesCompare, 'SM'); % For comparisons / verification
    D_final2 = sqrt(sort(diag(D_final2)))/2/pi;

    Diffs2 = D_initial2(7:end) - D_final2(7:end);
    log10( abs(Diffs2 ./ D_initial2(7:end) ))'
end

%% Null-Space Transformation of total motion
[V, D] = eigs(Krel(MESH.Nn*3+1:end, MESH.Nn*3+1:end), Mrel(MESH.Nn*3+1:end, MESH.Nn*3+1:end), 19, 'SM'); % Fully Fixed Interface for Eigen Analysis
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
    

%% Check Frequencies after Null Transformation
if(eliminateXb)

    % Without fixing contact
    [V_null, D_null] = eigs(K, M, n_modesCompare-nrbm, 'SM'); % For comparisons / verification
    D_null = sqrt(sort(diag(D_null)))/2/pi;

    Diffs_null = D_null(2:end) - D_final(8:end);
    log10( abs(Diffs_null ./ D_final(8:end) ) )'
    
    % with fixing contact
    [V_final, D_null2] = eigs(K(3*MESH.Nn+1:end, 3*MESH.Nn+1:end), ...
                                M(3*MESH.Nn+1:end, 3*MESH.Nn+1:end), n_modesCompare-nrbm, 'SM'); % For comparisons / verification
    D_null2 = sqrt(sort(diag(D_null2)))/2/pi;

    Diffs2_null = D_null2 - D_final2(7:end);
    log10( abs(Diffs2_null ./ D_final2(7:end)) )'
end

%% Save Matrices

if(~eliminateXb)
    save(sprintf('./MATS/%d_SET_NULLRED.mat', setid), 'M', 'K', 'R', 'L', 'Fv', 'MESH', 'Krel', 'Fvrel');
else
    save(sprintf('./MATS/%d_SET_NULLRED_NOREL.mat', setid), 'M', 'K', 'R', 'L', 'Fv', 'MESH', 'Krel', 'Fvrel', 'Thcb');
end
% save('tmp.mat', 'M', 'K', 'R', 'L', 'Fv', 'MESH')
