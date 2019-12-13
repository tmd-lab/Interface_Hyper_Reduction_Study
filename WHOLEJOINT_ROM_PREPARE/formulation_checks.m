clc
clear all
addpath('../ROUTINES/')
addpath('../ROUTINES/FEM/')

setid = 4;
fil = '4_INTSET';
fname = sprintf('../MATRIX_EXTRACTION/RUNS/%s/BRB_WOPRES_MAT.mat', fil);
load(fname, 'M', 'K', 'R', 'Fv');
K = 0.5*(K+K');
M = 0.5*(M+M');
Fv = Fv';
R = R';
fname = sprintf('../SUBSTRUCTURED_PRESTRESS/PREPARE/%d_SET_NULLRED.mat', setid);
load(fname, 'MESH')
%% Patches
% 1 - 1:77
% 2 - 67:256
% 3 - 246:435
% 4 - 425:614
% 5 - 604:680
[Q1, T1] = ZTE_ND2QP(MESH, 1);

Pnds = {1:77, 67:256, 246:435, 425:614, 604:680};
Npatches = 5;
Pels = cell(1, Npatches);
Area = cell(1, Npatches);
qps = cell(1, Npatches);
ctrds = zeros(Npatches, 2);
for n=1:5
    Pels{n} = find(sum(ismember(MESH.Quad(:, 2:end), Pnds{n}), 2)==4);
    Area{n} = T1(Pnds{n}, Pels{n})*ones(length(Pels{n}), 1);
    qps{n} = Q1(Pels{n}, :)*MESH.Nds;
    ctrds(n, :) = sum(Area{n}.*MESH.Nds(Pnds{n}, :))/sum(Area{n});
end

%% Relative Coordinate Transformation: [XB-XT; XB; eta]
Ngen = length(M) - MESH.Nn*3*2;
Trel = [-speye(MESH.Nn*3), speye(MESH.Nn*3), sparse(MESH.Nn*3, Ngen);
        sparse(MESH.Nn*3, MESH.Nn*3), speye(MESH.Nn*3), sparse(MESH.Nn*3, Ngen);
        sparse(Ngen, MESH.Nn*3*2), speye(Ngen)];
Krel = Trel'*K*Trel;
Mrel = Trel'*M*Trel;
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
NF = NTN*kron(BNV, speye(3));

%% Reduced Order Model
% Mnew = sparse(blkdiag(zeros(Npatches*6, Npatches*6), Mrel));
% Pmat = sparse(blkdiag([speye(Npatches*6) -GN;
%                         -GN' GN'*GN], ...
%                         sparse(MESH.Nn*3+Ngen, MESH.Nn*3+Ngen)));
% Knew = sparse(blkdiag(sparse(Npatches*6, Npatches*6), Krel));
% 
% Fvnew = [zeros(Npatches*6,1); Fvrel];
% Rnew = [zeros(size(Rrel,1), Npatches*6) Rrel];
% vdofs = reshape(1:(Npatches*6), 6, Npatches);
% vdofs = reshape(vdofs(1:3, :), Npatches*3, 1);
% % vdofs = 1:(Npatches*6);

Mnew = sparse(blkdiag(Mrel, sparse(Npatches*6, Npatches*6)));
Pmat = [GN'*GN sparse(MESH.Nn*3, MESH.Nn*3+Ngen) -GN';
        sparse(MESH.Nn*3+Ngen, length(K)+Npatches*6);
        -GN sparse(Npatches*6, MESH.Nn*3+Ngen) speye(Npatches*6)];
Knew = sparse(blkdiag(Krel, sparse(Npatches*6, Npatches*6))) + Pmat;

Fvnew = [Fvrel; zeros(Npatches*6,1)];
Rnew = [Rrel zeros(size(Rrel,1), Npatches*6)];
vdofs = length(K) + reshape(1:(Npatches*6), 6, Npatches);
vdofs = reshape(vdofs(1:3, :), Npatches*3, 1);
% vdofs = [3:3:(MESH.Nn*3) reshape(vdofs(1:2, :), 1, Npatches*2)];

[Mhcb, Khcb, Thcb] = HCBREDUCE(Mnew, Knew, vdofs, Ngen);
Fvhcb = Thcb'*Fvnew;
Rhcb = Rnew*Thcb;

% %% ROM 2 - No Considerable Improvement
% GN2 = blkdiag(GN, GN);
% Nps = Npatches*6*2;
% Mnew2 = sparse(blkdiag(sparse(Nps, Nps), M));
% Pmat2 = sparse(blkdiag([speye(Nps) -GN2;
%                         -GN2' GN2'*GN2], sparse(Ngen, Ngen)));
% Knew2 = sparse(blkdiag(sparse(Nps, Nps), K));
% Fvnew2 = [zeros(Nps,1); Fv];
% Rnew2 = [zeros(size(Rrel,1), Nps) R];
% vdofs2 = reshape(1:Nps, 6, Npatches*2);
% vdofs2 = reshape(vdofs2(1:3, :), Npatches*3*2, 1);
% 
% [Mhcb2, Khcb2, Thcb2] = HCBREDUCE(Mnew2, Knew2, vdofs2, Ngen);
% Fvhcb2 = Thcb2'*Fvnew2;
% Rhcb2 = Rnew2*Thcb2;

%% ROM3
Pmat = [-NTN sparse(MESH.Nn*3, MESH.Nn*3+Ngen) NTG;
        sparse(MESH.Nn*3+Ngen, length(K)+Npatches*6);
        NTG' sparse(Npatches*6, MESH.Nn*3+Ngen) speye(Npatches*6)];

Mnew = sparse(blkdiag(Mrel, sparse(Npatches*6, Npatches*6)));
Knew = sparse(blkdiag(Krel, sparse(Npatches*6, Npatches*6))) + Pmat;

Fvnew = [Fvrel; zeros(Npatches*6,1)];
Rnew = [Rrel zeros(size(Rrel,1), Npatches*6)];
vdofs = length(K) + reshape(1:(Npatches*6), 6, Npatches);
vdofs = reshape(vdofs(1:3, :), Npatches*3, 1);
% vdofs = [3:3:(MESH.Nn*3) reshape(vdofs(1:2, :), 1, Npatches*2)];

[Mhcb, Khcb, Thcb] = HCBREDUCE(Mnew, Knew, vdofs, Ngen);
Fvhcb = Thcb'*Fvnew;
Rhcb = Rnew*Thcb;
disp('Done')
%% EIGENCHECK
D = eigs(K, M, 40, 'SM');
Drel = eigs(Krel, Mrel, 40, 'SM');
Dnew = eigs(Knew, Mnew, 40, 'SM');
Dhcb = eigs(Khcb, Mhcb, 40, 'SM');

sqrt(abs([D Drel Dnew Dhcb]))/(2*pi)