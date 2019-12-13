clc
clear all
addpath('../ROUTINES/')

set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex'); 
set(0, 'DefaultLegendInterpreter', 'latex'); 
set(0,'defaultAxesFontSize',13)

symbs = {'s-', 'd-', 'p-'};
% ROM levels
% sel_method = 'U';
% Nels = 36; % [36 74 122 232 448 568 588]
% Nelrange = [74 122 232];

% sel_method = 'P';
% Nels = 52;  % [52 100 140 204 240 304]
% Nelrange = [52 140 304];

sel_method = 'PD';
% Nels = 252;  % [68 108 152 200 252 292]
Nelrange = [68 152 292];

mi = 1;
%% Load Reference Data
exx = load('../FULLJOINT_ROM_PREPARE/DATS/statsol_P11580_0.000001.mat', 'Tn', 'Qm','MESH', 'V');

%% Plots
xr = range(exx.MESH.Nds(:,1));
yr = range(exx.MESH.Nds(:,2));

figure(1)
clf()
set(gcf,'color','white')

ys = [yr 0 -yr]*1.5;
for i=1:length(Nelrange)
    Nels = Nelrange(i);
    load(sprintf('../FULLJOINT_ROM_PREPARE/ROMS/ROM_%s_%dELS.mat',sel_method,Nels), 'MESH');
    load(sprintf('./DATS/RUN_M%d_%s_%dELS.mat',mi,sel_method,Nels), 'Tn');
    Qm = ZTE_ND2QP(MESH, fix(sqrt(size(Tn,1)/MESH.Ne)));
    Pn = Qm\Tn;
    Pn(Pn<0) = 0;
    
    SHOW2DMESH(MESH.Nds+[0 ys(i)], MESH.Tri, MESH.Quad, Pn, -1, -100);
    
    text(0, ys(i)+yr*0.6, sprintf('%d Element ROM', Nels), 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'FontSize', 14)
end
axis equal;
axis off; colormap(jet(11));
% xx=colorbar('southoutside', 'FontSize', 14);
caxis([0 2.5e7])
% xlabel(xx, 'Normal Pressure (Pa)')

% print(sprintf('./FIGURES/ROM_PRES_%s.eps', sel_method), '-depsc')
export_fig(sprintf('./FIGURES/ROM_PRES_%s.eps', sel_method))
return
%% Reference - Smoothing
AdjMx = NODEELADJ(exx.MESH.Nds, exx.MESH.Quad(:,2:end)); AdjMx = AdjMx./sum(AdjMx,2);
AdjMx = AdjMx + eye(exx.MESH.Nn); AdjMx = AdjMx./sum(AdjMx,2);

%% Reference - Plot
figure(2)
set(gcf,'color','white')

clf()
set(gca, 'DataAspectRatioMode', 'manual')

Pn = exx.Qm\exx.Tn;
Pn(Pn<0) = 0;
    
SHOW2DMESH(exx.MESH.Nds, exx.MESH.Tri, exx.MESH.Quad, AdjMx^0*Pn*1e-6, -1, -100);

colormap(jet(11));
xx=colorbar('southoutside', 'FontSize', 14);
caxis([0 2.5e7]*1e-6)
ylim(xr/2*[-1 1])
ylim(yr/2*[-1 1])
xlabel(xx, 'Normal Pressure (GPa)')

set(gcf, 'Position', [595 687 560 200])
set(gca, 'Position', [0.13 0.4 0.78 0.55])
set(gca, 'OuterPosition', [0 0.35 1.0 0.75])
axis equal
axis off;

% print('./FIGURES/PREST_REF.eps', '-depsc')
export_fig('./FIGURES/PREST_REF.eps')

%% MODE 1
figure(3)
clf()
set(gcf,'color','white')

SHOW3D(exx.MESH.Nds, exx.MESH.Tri, exx.MESH.Quad, exx.V(1:3:exx.MESH.Nn*3,1), exx.V(2:3:exx.MESH.Nn*3,1), -exx.V(3:3:exx.MESH.Nn*3,1), 1, exx.V(3:3:exx.MESH.Nn*3,1))
axis off
set(gcf, 'Position', [595 687 560 200])
set(gca, 'View', [-15 40])
xx=colorbar('south', 'FontSize', 14);
title(xx, 'Displacement')
set(xx, 'TickDirection', 'out');
set(xx, 'AxisLocationMode', 'manual');
set(xx, 'AxisLocation', 'out');

% print('./FIGURES/MODE1_IDISP.eps','-depsc')
export_fig('./FIGURES/MODE1_IDISP.eps')