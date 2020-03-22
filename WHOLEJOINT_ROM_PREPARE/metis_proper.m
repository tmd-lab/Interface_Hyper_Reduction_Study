addpath('../ROUTINES/')
addpath('../ROUTINES/FEM/')
addpath('../ROUTINES/GRAPH/')

setid = 5;
fil = {'1_AROUNDSET', '2_ABOVESET', '3_SINGELEMABOVESET', '4_INTSET', '5_INTSETNPS'};
Ngens = [26, 26, 26, 26, 19];
% [3 , 5 , 7, 9, 11, 19, 21, 25, 30, 40, 45, 50, 60, 70, 80, 100]
for Nlevs = 90
% clc
clearvars -except Nlevs Ngens setid fil

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

tic
%% Parameters
% Nlevs = 3;  % [3 , 5 , 7, 9, 11, 19, 21, 25, 30, 40, 45, 50]

sel_method = 'U';

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
    
    SHOW2DMESH(MESH.Nds, zeros(0,4), [(1:length(Pels{n}))' MESH.Quad(Pels{n},2:end)], n, -1, -100, MESH.Ne)
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
NTN = 0.5*(NTN+NTN');

%% Reduced Order Model
ni = 1:(MESH.Nn*3);
mi = setdiff(1:length(Krel), ni);
cnum = 1e0;

M1 = sparse(blkdiag(Mrel, zeros(2*Npatches*6)));
K1 = sparse([Krel(ni,:) zeros(MESH.Nn*3, Npatches*6) -cnum*NTG;
      Krel(mi,:) zeros(length(mi),2*Npatches*6);
      zeros(Npatches*6, Npatches*6+size(K,1)) cnum*GTG;
      -cnum*NTG' zeros(Npatches*6, length(mi)) cnum*GTG zeros(Npatches*6)]);
R1 = sparse([Rrel zeros(size(R,1), 2*Npatches*6)]);
Fv1 = [Fvrel; zeros(2*Npatches*6,1)];
disp(condest(K1))

dofred = 6;
Ngen = Ngens(setid);
% Ngen = 36;
vdofs = length(K) + reshape(((1:Npatches)-1)*6+(1:dofred)', [],1);
[Mhcb, Khcb, Thcb] = HCBREDUCE(M1, K1, vdofs, Ngen);
Fvhcb = Thcb'*Fv1;
Rhcb = R1*Thcb;
disp('HCB DONE!')

%% NULL-SPACE REDUCTION
rdofs = setdiff(1:size(Khcb,1), reshape(((1:Npatches)-1)*6+(1:dofred)',[],1));
[Vhcb, Dhcb] = eig(full(Khcb(rdofs, rdofs)), full(Mhcb(rdofs, rdofs))); 
% [Vhcb, Dhcb] = eig(full(Khcb), full(Mhcb));
[Dhcb, si] = sort(diag(Dhcb));  Vhcb = Vhcb(:, si);
Vrbms = zeros(size(Khcb,1), 6); Vrbms(rdofs,:) = Vhcb(:,1:6);
Vrbms = Vrbms./sqrt(diag(Vrbms'*Mhcb*Vrbms)');
L = null(Vrbms'*Mhcb);

M = L'*Mhcb*L;
M = 0.5*(M+M');
K = L'*Khcb*L;
K = 0.5*(K+K');
Fv = L'*Fvhcb;
R = Rhcb*L;
Th = Thcb*L;
% [sort(eig(K), 'descend') sort(eig(M), 'descend')]
% %% Check MAC values
% [Vrel, Drel] = eigs(Krel, Mrel, size(Khcb,1), 'SM');
% [Drel, si] = sort(diag(Drel));  Vrel = Vrel(:, si);  Vrel = Vrel./sqrt(diag(Vrel'*Mrel*Vrel)');
% 
% [Vhcb, Dhcb] = eig(full(Khcb), full(Mhcb));
% [Dhcb, si] = sort(diag(Dhcb));  Vhcb = Vhcb(:, si);  Vhcb = Vhcb./sqrt(diag(Vhcb'*Mhcb*Vhcb)');
% VH = Thcb(1:size(Krel,1), :)*Vhcb;
% mhs = find(isfinite(Dhcb));

% [Vr, Dr] = eig(full(K), full(M));
% [Dr, si] = sort(diag(Dr)); Vr = Vr(:, si); Vr = Vr./sqrt(diag(Vr'*M*Vr)');
% VH = Th(1:size(Krel,1), :)*Vr;
% mhs = find(isfinite(Dr));
% %%
% figure(2)
% clf()
% % br=bar3((VH(:, mhs)'*Mrel*Vrel).^2./(diag(VH(:,mhs)'*Mrel*VH(:,mhs)).*diag(Vrel'*Mrel*Vrel)'));
% % for k=1:size(br,2); br(k).CData = br(k).ZData; br(k).FaceColor = 'interp'; end
% imagesc((VH(:, mhs)'*Mrel*Vrel).^2./(diag(VH(:,mhs)'*Mrel*VH(:,mhs)).*diag(Vrel'*Mrel*Vrel)'));
% yy=colorbar; caxis([0 1]);
% ylabel(yy, 'MAC')
% axis equal
% xlim([1 size(Khcb,1)])
% ylim([1 size(Khcb,1)])
% set(gca, 'ydir', 'normal')
% hold on;
% rectangle('position', [0 0 7.5 7.5], 'LineWidth', 2)
% % rectangle('position', [7.5 7.5 7.5+3*Npatches+0.5 7.5+3*Npatches+0.5], 'LineWidth', 2)
% 
% xlabel('Original Model')
% ylabel('ROM')

%% Save
ttk = toc
fname = sprintf('./REDMATS/P%d_S%.2f_%s_%dLEV_GRED_WJMAT_prop.mat', Prestress, log10(sint), sel_method, Nlevs);
save(fname, 'M', 'K', 'L', 'Fv', 'R', 'Th', 'I2', 'MESH', 'PatchAreas', 'Npatches', 'dofred', 'NTN', 'GTG', 'NTG', 'cnum', 'Pels', 'Pnds', 'ttk')
disp('SAVED!')
end