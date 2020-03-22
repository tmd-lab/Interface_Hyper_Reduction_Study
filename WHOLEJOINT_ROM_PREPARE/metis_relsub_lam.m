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
[Q1, T1] = ZTE_ND2QP(MESH, 1);
[Q3, T3] = ZTE_ND2QP(MESH, 3);
EAreas = sum(T1,1);

%% Parameters
Nlevs = 50;  % [3 , 5 , 7, 9, 11, 19, 21, 25, 30, 40]
sel_method = 'U';

%% Processing for METIS
% %% Using gpmetis
% % %  In terms of Nodes
% % AdjMx = NODEELADJ(MESH.Nds, MESH.Quad(:, 2:end));
% % fprintf(Gpfile, '%d %d\n', MESH.Nn, length(find(AdjMx))/2);
% % ind = 1;  % In terms of nodes
% % In terms of Elements
% AdjMx = ELELADJ(MESH);
% AdjMx = AdjMx.*(AdjMx>=2 | AdjMx==0);  % Choose only those with two
%                                        % nodes connected
% ind = 0;  % In terms of elements
%
% Gpfile = fopen('./DATS/brb_graph.mti', 'w+');
% fprintf(Gpfile, '%d %d\n', MESH.Ne, length(find(AdjMx))/2);
% for i=1:(MESH.Nn*ind + MESH.Ne*(~ind))
%     fprintf(Gpfile, [int2str(find(AdjMx(i,:))) '\n']);
% end
% fclose(Gpfile);
%
% % Call Metis
% system(sprintf('cd DATS; gpmetis -contig brb_graph.mti %d > /dev/null', ...
%                Nlevs));
% %% Collect DataSHOW2DMESH(MESH.Nds, [], [(1:length(Pels{n}))' MESH.Quad(Pels{n},2:end)], n, -1, -100, MESH.Ne)
% metisret = dlmread(sprintf('./DATS/brb_graph.mti.part.%d', Nlevs));
%
% Pnds = cell(Nlevs, 1);
% Pels = cell(Nlevs, 1);
% for i=1:Nlevs
%     if (ind == 1)
%         Pnds{i} = find(metisret==(i-1));
%         Pels{i} = find(sum(ismember(MESH.Quad(:, 2:end), Pnds{i}), ...
%                            2)==4);
%     else
%         Pels{i} = find(metisret==(i-1));
%         Pnds{i} = reshape(MESH.Quad(Pels{i}, 2:end), [], 1);
%     end
% end
% Npatches = Nlevs;

%% Using mpmetis
EAreasi = fix(EAreas*1e8);
Mpfile = fopen('./DATS/brb_mesh.mti', 'w+');
fprintf(Mpfile, '%d 1\n', MESH.Ne);
for e=1:MESH.Ne
    fprintf(Mpfile, [int2str([EAreasi(e) MESH.Quad(e,2:end)]) '\n']);
end
fclose(Mpfile);

% Call Metis
system(sprintf('cd DATS; mpmetis -contig -ncommon=2 -gtype=dual brb_mesh.mti %d > /dev/null', ...
               Nlevs));

%% Collect Data
metisret_E = dlmread(sprintf('./DATS/brb_mesh.mti.epart.%d', ...
                             Nlevs));
metisret_N = dlmread(sprintf('./DATS/brb_mesh.mti.npart.%d', Nlevs));

Pnds = cell(Nlevs, 1);
Pels = cell(Nlevs, 1);
for i=1:Nlevs
    Pnds{i} = find(metisret_N==(i-1));
    Pels{i} = find(metisret_E==(i-1));
end
Npatches = Nlevs;

%% PATCHING
Prestress = 11580; sint = 1e-6;
Area = cell(1, Npatches);
qps = cell(1, Npatches);
ctrds = zeros(Npatches, 2);
I2 = zeros(4, Npatches);  % Ixx, Ixy, Iyy, Izz

figure(1)
clf()
for n=1:Npatches
    Area{n} = T1(Pnds{n}, Pels{n})*ones(length(Pels{n}), 1);
    qps{n} = Q1(Pels{n}, :)*MESH.Nds;
    ctrds(n, :) = sum(Area{n}.*MESH.Nds(Pnds{n}, :))/sum(Area{n});
    
    qids = reshape((Pels{n}-1)'*9+(1:9)', [], 1);
    I2(1,n) = sum(T3(:, qids)*(Q3(qids,:)*MESH.Nds(:,2)-ctrds(n,2)).^2);  % Ixx
    I2(2,n) = -sum(T3(:, qids)*prod(Q3(qids,:)*MESH.Nds-ctrds(n,:), 2));  % Ixy
    I2(3,n) = sum(T3(:, qids)*(Q3(qids,:)*MESH.Nds(:,1)-ctrds(n,1)).^2);  % Iyy
    I2(4,n) = sum(T3(:, qids)*sum((Q3(qids,:)*MESH.Nds-ctrds(n,:)).^2, 2));  % Izz
    
    SHOW2DMESH(MESH.Nds, [], [(1:length(Pels{n}))' MESH.Quad(Pels{n},2:end)], n, -1, -100, MESH.Ne)
end
PatchAreas = cellfun(@(c) sum(c), Area);
axis equal; axis off; colormap(jet(Npatches));
colormap(prism)
title(sprintf('%d Patches',Npatches))
set(gca, 'Position', [-0.5 -0.1 2 1])
% print(sprintf('./FIGS/P%d_S%.2f_%s_%dLEV_%dPATCH_MESH.eps',Prestress, log10(sint), sel_method, Nlevs,Npatches), '-depsc')
pause(1)
% return
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

%% Check
% Dred = eigs(K, M, 10, 'SM');
% Drel = eigs(Krel, Mrel, 16, 'SM');
% [Dred Drel(7:end)]
%% Save
fname = sprintf('./REDMATS/P%d_S%.2f_%s_%dLEV_GRED_WJMAT.mat', Prestress, log10(sint), sel_method, Nlevs);
save(fname, 'M', 'K', 'L', 'Fv', 'R', 'Th', 'I2', 'LamT', 'MESH', 'PatchAreas', 'Npatches', 'dofred', 'NTN', 'GTG', 'NTG', 'cnum', 'Pels', 'Pnds')
disp('SAVED!')