clc
clear all
% close all
addpath('../ROUTINES/')
addpath('../ROUTINES/REMESHING')
addpath('../ROUTINES/export_fig/')

set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex'); 
set(0, 'DefaultLegendInterpreter', 'latex'); 
set(0,'defaultAxesFontSize',13)

%% Load Reference Data
exx = load('../FULLJOINT_ROM_PREPARE/DATS/statsol_P11580_0.000001.mat', 'Tn', 'Qm','MESH');
% 'P'[2,3,4,5,6,11] 'PD' [2,3,4,5,6,7,10] 'U' [3,5,7,9,11,19,21,30,40,45,50,60,70,80,90]
NLc = containers.Map;
NLc('P') = [2,3,4,5,6,7];
NLc('PD') = [2,3,4,5,6,7,10];
NLc('U') = [3,5,7,9,11,19,21,30,40,45,50,60,70,80,90];

% sel_method = 'P'
% Nlevs = [2 4 7];

% sel_method = 'PD'
% Nlevs = [2 6 10];

sel_method = 'U'
Nlevs = [3 11 19];

%% Plots
xr = range(exx.MESH.Nds(:,1));
yr = range(exx.MESH.Nds(:,2));

[~, bi, ~] = CreateBRBInterface(exx.MESH.Nds(:,1), exx.MESH.Nds(:, ...
                                                  2));
bi = cellfun(@(c) c(convhull(exx.MESH.Nds(c,:))), bi, 'UniformOutput', false);

figure(1)
clf()
set(gcf,'color','white')

ys = [yr 0 -yr]*1.5;
for i=1:length(Nlevs)
    load(sprintf('./DATS/P11580_S-6.00_%s_%dLEV_BB_ELDRYFRICT_prop_wth.mat',sel_method,Nlevs(i)), 'Pxyn', 'Pels','MESH')
    
    Pn = zeros(MESH.Ne,1);
    for k=1:length(Pels)
        Pn(Pels{k}) = Pxyn(k,3)*1e-6;
    end
    
    SHOW2DMESH(MESH.Nds+[0 ys(i)], zeros(0,4), MESH.Quad, Pn, -1, - ...
               101);
    chi = convhull(MESH.Nds);
    plot(MESH.Nds(chi,1), MESH.Nds(chi,2)+ys(i), 'k-')
    for ii=1:3
        plot(MESH.Nds(bi{ii},1), MESH.Nds(bi{ii},2)+ys(i), 'k-')
    end
    
    text(0, ys(i)+yr*0.6, sprintf('%d Patch ROM', length(Pels)), 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'FontSize', 14)
end
axis equal;
axis off; colormap(jet(11));
caxis([0 25])

% print(sprintf('./FIGURES/WJROM_PRES_%s.eps', sel_method), '-depsc')
% export_fig(sprintf('./FIGURES/WJROM_PRES_%s.eps', sel_method))