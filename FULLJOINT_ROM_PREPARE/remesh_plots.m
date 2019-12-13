clc
clear all
close all
addpath('../ROUTINES/')
addpath('../ROUTINES/FEM/')
addpath('../ROUTINES/export_fig/')

set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex'); 
set(0, 'DefaultLegendInterpreter', 'latex'); 
set(0,'defaultAxesFontSize',13)

symbs = {'s-', 'd-', 'p-'};
% ROM levels
% sel_method = 'U';
% % Nels = 36; % [36 74 122 232 448 568 588]
% Nelrange = [74 122 232];

sel_method = 'P';
% Nels = 52;  % [52 100 140 204 240 304]
Nelrange = [52 140 304];

% sel_method = 'PD';
% % Nels = 252;  % [68 108 152 200 252 292]
% Nelrange = [68 152 292];

ref = load('../FULLJOINT_ROM_PREPARE/DATS/statsol_P11580_0.000001.mat', 'MESH'); % Reference Mesh

%% Plots
xr = range(ref.MESH.Nds(:,1));
yr = range(ref.MESH.Nds(:,2));

if sel_method~='U'
    load(sprintf('./MATS/REDMESH_%s_%dELS.mat',sel_method,Nelrange(1)), 'objective');
    figure(1)
    clf()
    set(gcf,'color','white')

    SHOW2DMESH(ref.MESH.Nds, ref.MESH.Tri, ref.MESH.Quad, objective, -1, -100);

    yy = colorbar('southoutside', 'FontSize', 14);
    xlabel(yy, sprintf('%s Field Objective',sel_method))
    colormap(jet(11))

    set(gcf, 'Position', [595 687 560 200])
    set(gca, 'Position', [0.13 0.4 0.78 0.55])
    set(gca, 'OuterPosition', [0 0.35 1.0 0.75])
    axis equal; axis off

%     print(sprintf('./FIGS/Fobjective_%s.eps',sel_method), '-depsc')
    export_fig(sprintf('./FIGS/Fobjective_%s.eps',sel_method))
end
return

figure(2)
clf()
set(gcf,'color','white')

ys = [yr 0 -yr]*1.5;
for i=1:length(Nelrange)
    Nels = Nelrange(i);
    if sel_method~='U'
        load(sprintf('./MATS/REDMESH_%s_%dELS.mat',sel_method,Nels), 'MESH', 'Qn', 'objective');
        Pn = Qn\objective;
    else
        load(sprintf('./MATS/REDMESH_%s_%dELS.mat',sel_method,Nels), 'MESH');
        Pn = 1;
    end
    
    SHOW2DMESH(MESH.Nds+[0 ys(i)], MESH.Tri, MESH.Quad, Pn, -1, -100);
    
    text(0, ys(i)+yr*0.6, sprintf('%d Element ROM', Nels), 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'FontSize', 14)
end
axis equal;
axis off; colormap(jet(11));
if sel_method~='U'
%     xx=colorbar('southoutside', 'FontSize', 14);
    caxis([0 1])
%     xlabel(xx, 'Field Objective')
end
% print(sprintf('./FIGS/REDMESH_%s.eps', sel_method), '-depsc')
export_fig(sprintf('./FIGS/REDMESH_%s.eps', sel_method))