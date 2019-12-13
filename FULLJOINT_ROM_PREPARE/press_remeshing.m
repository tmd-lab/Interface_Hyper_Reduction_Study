% clc
clear all
addpath('../ROUTINES/')
addpath('../ROUTINES/FEM/')
addpath('../ROUTINES/REMESHING/')
addpath('../ROUTINES/GRAPH/')
addpath('../ROUTINES/REMESHING/NewNidish/')

%% Parameters for meshes
% 1. sel_method: Pressure


sel_method = 'PD';

switch sel_method
  case 'P' % 'Pressure'
    % mds = [1];  sn = 10;
    % 52 element mesh: (0.0,1.0); (50,10)
    % 100 element mesh: (0.0,1.0); (100,10)
    % 140 element mesh: (0.0,1.0); (150,10)
    % 204 element mesh: (0.0,1.0); (200,10)
    % 240 element mesh: (0.0,1.0); (250,10)
    % 304 element mesh: (0.0,1.0); (300,10)

    % Modes for weighting
    mds = [1];
    % Weights for adding displacements and pressures
    wd = 0; wp = 1;
    % Number of Kernel Smoothening steps
    sn = 10;
    % Settings for sel_method:Pressure
    Nel = 304;  Nel_tol = 10;
  case 'PD' % 'Pressure'
    % mds = [1];  sn = 10;
    % 48 element mesh: (1.0,1.0); (50,10)
    % 108 element mesh: (1.0,1.0); (100,10)
    % 152 element mesh: (1.0,1.0); (150,10)
    % 200 element mesh: (1.0,1.0); (200,10)
    % 252 element mesh: (1.0,1.0); (250,10)
    % 292 element mesh: (1.0,1.0); (300,10)
    
    % FRICTMODEL STUDY
    % 68 element mesh: (1.0, 1.0); (60,10)
    
    % Modes for weighting
    mds = [1];      
    % Weights for adding displacements and pressures
    wd = 1; wp = 1;
    % Number of Kernel Smoothening steps
    sn = 10;
    % Settings for sel_method:Pressure
    Nel = 292;  Nel_tol = 10;    
  case 'U' % 'Uniform'
end
%% Options
% fname = './DATS/statsol_P11580.mat';
fname = './DATS/statsol_P11580_0.000001.mat';
%% Input section
dat = load(fname, 'Tx', 'Ty', 'Tn', 'Xs', 'Qm', 'Tm', 'MESH','V');

%% Processing old mesh
% [Elements_idx]=CreatElementsFromNodes(dat.MESH.Nds(:,1),dat.MESH.Nds(:,2));

Elements_idx_p = mat2cell(dat.MESH.Quad(:,end:-1:2), ones(dat.MESH.Ne_Quad,1));

if sel_method~='U'
    P2F = MESH2D_P2FMAT(dat.MESH);
    [Qm,Tm] = ZTE_ND2QP(dat.MESH,2);

    % Normalized normal displacement
%     Vn_objective = sum(abs(dat.V(dat.MESH.Nn*3+(3:3:dat.MESH.Nn*3),mds))+abs(dat.V(3:3:dat.MESH.Nn*3,mds)),2);
%     Vn_objective = sum(dat.V(dat.MESH.Nn*3+(3:3:dat.MESH.Nn*3),mds)+dat.V(3:3:dat.MESH.Nn*3,mds),2);
    Vn_objective = sum(dat.V(3:3:dat.MESH.Nn*3,mds),2);
    Vn_objective = (Vn_objective-min(Vn_objective))/range(Vn_objective);

    % Normalize normal traction
    Tn_objective = dat.Qm\dat.Tn;
    Tn_objective(Tn_objective<0) = 0;
    Tn_objective = (Tn_objective-min(Tn_objective))/range(Tn_objective);

    % Add traction and displacement objectives
    objective = (wd*Vn_objective+wp*Tn_objective)/(wd+wp);
    %% Smoothening
    % Kernel Smoothing
    AdjMx = NODEELADJ(dat.MESH.Nds, dat.MESH.Quad(:,2:end));
    AdjMx = AdjMx./sum(AdjMx,2);
    % AdjMx(find(AdjMx)) = 1;
    AdjMx = AdjMx + eye(dat.MESH.Nn);
    AdjMx = AdjMx./sum(AdjMx,2);
    obji = AdjMx^sn*objective;
    obji = (obji-min(obji))/range(obji);

    % figure(1);clf(); 
    % SHOW2DMESH(dat.MESH.Nds+[0 0.02], dat.MESH.Tri, dat.MESH.Quad, objective, -1, -100);
    % SHOW2DMESH(dat.MESH.Nds-[0 0.02], dat.MESH.Tri, dat.MESH.Quad, obji, -1, -100); axis equal; axis off
    % yy = colorbar('southoutside', 'FontSize', 16);
    % xlabel(yy, 'Field Objective')
    % colormap(jet(11))
    % disp('Done')
    % pause(0.0001)

    % keyboard
    objective = obji;

    %% Element Generation (New Mesh)
    [NewElements_idx,NewNodelList] = ...
        PressureBasedDistribution(dat.MESH.Nds(:,1),dat.MESH.Nds(:,2),dat.MESH.Nds(:,1)*0,objective, Nel, Nel_tol);
    else 
        % [NewElements_idx, sym_flag] = Uniform_Distribution_v02(Elements_idx, dat.MESH.Nds(:,1), dat.MESH.Nds(:,2), 10,5,2);
        [NewElements_idx, sym_flag] = Uniform_Distribution_v02(Elements_idx_p, dat.MESH.Nds(:,1), dat.MESH.Nds(:,2));
end

e3=find(cellfun(@(c) length(c), NewElements_idx(:,1))==3);
e4=find(cellfun(@(c) length(c), NewElements_idx(:,1))==4);

% Quality Check
QualityCheckMesh([(1:size(NewNodelList,1))' NewNodelList],[e4 cell2mat(NewElements_idx(e4,1))],[e3 cell2mat(NewElements_idx(e3,1))],1);
% QualityCheckMesh([(1:dat.MESH.Nn)' dat.MESH.Nds], dat.MESH.Quad, dat.MESH.Tri, 1);

%% Reduced MESH structure
MESH.Nds = NewNodelList(:,1:2);
MESH.Quad = [e4 cell2mat(NewElements_idx(e4,1))];
MESH.Tri = [e3 cell2mat(NewElements_idx(e3,1))];
MESH.Nn = size(MESH.Nds,1);
MESH.Ne_Quad = size(MESH.Quad,1);
MESH.Ne_Tri= size(MESH.Tri,1);
MESH.Ne = MESH.Ne_Quad+MESH.Ne_Tri;

Qn = MESH2PTS(MESH, dat.MESH.Nds);
save(sprintf('./MATS/REDMESH_%s_%dELS.mat', sel_method, MESH.Ne), 'MESH', 'objective', 'Qn')

% %% Plotting
% figure(1);clf(); 
% SHOW2DMESH(MESH.Nds+[0 0.02], MESH.Tri, MESH.Quad, 0.5, -1, -100);
% SHOW2DMESH(dat.MESH.Nds-[0 0.02], dat.MESH.Tri, dat.MESH.Quad, objective, -1, -100); axis equal; axis off
% yy = colorbar('southoutside', 'FontSize', 16);
% xlabel(yy, 'Field Objective')
% colormap(jet(11))
% pause(0.0001)
%% Plotting field objective

%%
figure(1);clf();
% SHOW2DMESH(MESH.Nds+[0 0.02], MESH.Tri, MESH.Quad, Qn\objective, -1, -100);
% SHOW2DMESH(dat.MESH.Nds-[0 0.02], dat.MESH.Tri, dat.MESH.Quad, objective, -1, -100); 
SHOW2DMESH(MESH.Nds, MESH.Tri, MESH.Quad, Qn\objective, -1, -100);

yy = colorbar('southoutside', 'FontSize', 14);
xlabel(yy, 'Field Objective')
colormap(jet(11))

set(gcf, 'Position', [595 687 560 200])
set(gca, 'Position', [0.13 0.4 0.78 0.55])
set(gca, 'OuterPosition', [0 0.35 1.0 0.75])
axis equal; axis off

% print(sprintf('./FIGS/REDMESH_%s_%dELS.eps', sel_method, MESH.Ne), '-depsc')