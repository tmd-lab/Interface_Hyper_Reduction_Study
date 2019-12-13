clc
clear all
close all
addpath('../ROUTINES/')
addpath('../ROUTINES/export_fig')
addpath('../ROUTINES/FEM/')
addpath('../ROUTINES/GRAPH/')
addpath('../ROUTINES/REMESHING')

set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex'); 
set(0, 'DefaultLegendInterpreter', 'latex'); 
set(0,'defaultAxesFontSize',13)

setid = 5;
fname = sprintf('./MATS/%d_SET_NULLRED.mat', setid);
load(fname, 'MESH')
[Q1, T1] = ZTE_ND2QP(MESH, 1);

%%
% sel_method = 'P' % 'P'(2,3,4,5,6,10) 'PD' (2,3,6,7,8,11) 'U' (3, 5, 7, 9, 11, 19, 21)
% wd = 0;  wp = 1;
% Nlevs = [2 4 10];
% sn = 10;

sel_method = 'PD' % 'P'(2,3,4,5,6,10) 'PD' (2,3,6,7,8,11) 'U' (3, 5, 7, 9, 11, 19, 21)
wd = 1;  wp = 1;
Nlevs = [2 6 11];
sn = 10;

% sel_method = 'U' % 'P'(2,3,4,5,6,10) 'PD' (2,3,6,7,8,11) 'U' (3, 5, 7, 9, 11, 19, 21)
% Nlevs = [3 11 21];

%% Recreate Objectives
if sel_method ~= 'U'
    % Data
    Prestress = 11580; sint = 1e-6;
    dat = load(sprintf('./DATS/statsol_P%d_S%.2f.mat',Prestress,log10(sint)), 'Xs', 'Qm', 'Tm', 'Tx', 'Ty', 'Tn', 'dR0', 'V');
    %Select
    PObj = dat.Qm\dat.Tn;  % Normal traction
    DObj = abs(dat.V(3:3:MESH.Nn*3,1)-dat.V(MESH.Nn*3+(3:3:MESH.Nn*3),1));  % Relative Displacement
    PObj = (PObj-min(PObj))/range(PObj);
    DObj = (DObj-min(DObj))/range(DObj);
    Obj = (wd*DObj+wp*PObj)/(wd+wp);
    % Smoothen
    AdjMx = NODEELADJ(MESH.Nds, MESH.Quad(:,2:end)); AdjMx = AdjMx./sum(AdjMx,2);
    AdjMx = AdjMx + eye(MESH.Nn); AdjMx = AdjMx./sum(AdjMx,2);
    Obj = AdjMx^sn*Obj;
    % Formulate
    Obje = Obj;
end

%% Plots
xr = range(MESH.Nds(:,1));
yr = range(MESH.Nds(:,2));

[~, bi, ~] = CreateBRBInterface(MESH.Nds(:,1), MESH.Nds(:,2));
bi = cellfun(@(c) c(convhull(MESH.Nds(c,:))), bi, 'UniformOutput', false);

if sel_method ~= 'U'
    figure(1)
    clf()
    set(gcf,'color','white')
    
    SHOW2DMESH(MESH.Nds, MESH.Tri, MESH.Quad, Obje, -1, -100);

    yy = colorbar('southoutside', 'FontSize', 14);
    xlabel(yy, sprintf('%s Field Objective',sel_method))
    colormap(jet(11))
    caxis([0 1])
    set(yy, 'Ticks', 0:0.2:1)

    set(gcf, 'Position', [595 687 560 200])
    set(gca, 'Position', [0.13 0.4 0.78 0.55])
    set(gca, 'OuterPosition', [0 0.35 1.0 0.75])
    axis equal; axis off

    % print(sprintf('./FIGS/WJ_Fobjective_%s.eps',sel_method), '-depsc')
    export_fig(sprintf('./FIGS/WJ_Fobjective_%s.eps',sel_method))
end

return

%%
figure(2)
clf()
set(gcf,'color','white')
ys = [yr 0 -yr]*1.5;
Ctrds = Q1*MESH.Nds;
for i=1:length(Nlevs)
    Nlev = Nlevs(i)
    load(sprintf('./REDMATS/P11580_S-6.00_%s_%dLEV_GRED_WJMAT.mat', sel_method, Nlev), ...
        'Pels', 'Pnds')
    
    for p=1:length(Pels)
        SHOW2DMESH(MESH.Nds+[0 ys(i)], MESH.Tri, MESH.Quad(Pels{p},:), ...
                   p, -1, -101, 0, MESH.Ne);
        [~, AdjE] = NODEELADJ(MESH.Nds, MESH.Quad(Pels{p},2:end));
        % tE = diag(AdjE); tE(tE==0) = 1.0;
        % chi = find(sum(AdjE./tE',2)>0 & sum(AdjE./tE',2)<4);
        % chi = chi(convhull(MESH.Nds(chi,:)));
        % plot(MESH.Nds(chi, 1), MESH.Nds(chi, 2)+ys(i), 'k-')
    end
    chi = convhull(MESH.Nds);
    plot(MESH.Nds(chi, 1), MESH.Nds(chi, 2)+ys(i), 'k-')
    for ii=1:3
        plot(MESH.Nds(bi{ii},1), MESH.Nds(bi{ii},2)+ys(i), 'k-')
    end    
    
    text(0, ys(i)+yr*0.6, sprintf('%d Levels: %d Patch ROM', Nlev, length(Pels)), 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'FontSize', 14)
end
axis equal;
axis off;
colormap(prism(length(Pels)))

% print(sprintf('./FIGS/PATCHMESH_%s.eps', sel_method), '-depsc')
export_fig(sprintf('./FIGS/PATCHMESH_%s.eps', sel_method))
