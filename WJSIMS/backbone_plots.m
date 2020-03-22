clc
clear all
% close all
addpath('../ROUTINES/')
addpath('../ROUTINES/export_fig')

set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex'); 
set(0, 'DefaultLegendInterpreter', 'latex'); 
set(0,'defaultAxesFontSize',13)

%% Load Reference Data
exx = load('../FULLJOINT_REF/DATS/RUN2.mat', 'QS', 'WS', 'ZS', 'DS','ZSM','BBS');
ref = load('../FULLJOINT_REF/PREP/5_SET_NULLRED.mat','Mrel','L','Rrel');
Mref = ref.Mrel;
Lref = ref.L;
Rref = ref.Rrel;
MSref = Lref*cell2mat(cellfun(@(c) c.MS(:,end), exx.BBS, 'UniformOutput', false));
MSref = sign(Rref(3,:)*MSref).*MSref;
% 'P'[2,3,4,5,6,7,11] 'PD' [2,3,4,5,6,7,10] 'U' [3,5,7,9,11,19,21,30,40,45,50,60] 

sel_method = 'P'
Nlevs = [2 4 7];

% sel_method = 'PD'
% Nlevs = [2 6 10];

% sel_method = 'U'
% Nlevs = [3 11 19];

%% Plot Backbones
figure(1)
clf()
set(gcf,'color','white')

colos = DISTINGUISHABLE_COLORS(3);
aa = gobjects(4,1);

ax1=subplot('Position', [0.13 0.5 0.8 0.45]);
set(ax1, 'DataAspectRatioMode', 'manual')
aa(1) = semilogx(exx.QS, exx.WS/2/pi, 'k-', 'LineWidth', 2); hold on
legend(aa(1), sprintf('Reference (%d DoFs)',size(exx.BBS{1}.U,1)-1));

ax2=subplot('Position', [0.13 0.12 0.8 0.33]);
loglog(exx.QS, exx.ZSM, 'k-', 'LineWidth', 2); hold on

symbs = {'o-', 'v-', '^-'};
k=1;
for i=Nlevs
    load(sprintf('./DATS/P11580_S-6.00_%s_%dLEV_BB_ELDRYFRICT_prop_wth.mat',sel_method,i), 'BB')
    load(sprintf('../WHOLEJOINT_ROM_PREPARE/REDMATS/P11580_S-6.00_%s_%dLEV_GRED_WJMAT_prop.mat',sel_method, i),'M')
    
    subplot(ax1)
    aa(k+1) = semilogx(BB.Q, BB.W/2/pi, symbs{k}, 'Color', colos(k,:), 'MarkerFaceColor', colos(k,:));
    legend(aa(k+1), sprintf('%d Level ROM (%d DoFs)',i,size(M,1)))
    
    subplot(ax2)
    loglog(BB.Q, BB.Z, symbs{k}, 'Color', colos(k,:), 'MarkerFaceColor', colos(k,:))
    k=k+1;
end
subplot(ax1)
xticklabels([])
ll = legend(aa(1:end), 'Location', 'southwest');
set(ll, 'visible', 'off')
ylabel('Natural Frequency (Hz)')
xlim([1e-6 1e-1])
ylim([80 200])

subplot(ax2)
xlabel('Modal Amplitude')
ylabel('Damping factor')
xlim([1e-6 1e-1])
ylim([9e-6 1e0])

% set(ax1, 'Position', [0.13 0.5 0.8 0.45])
% set(ax2, 'Position', [0.13 0.12 0.8 0.33])

% ax1.Position = [0.13 0.5 0.8 0.45]
% ax2.Position = [0.13 0.12 0.8 0.33]

% print(sprintf('./FIGURES/%s_ELDRYFRICT_BB_WJROM.eps',sel_method), '-depsc')
export_fig(sprintf('./FIGURES/%s_ELDRYFRICT_BB_WJROM.eps',sel_method))

%% Closer in on the backbones
figure(2)
clf()
set(gcf, 'Color', 'white')

aa = gobjects(4,1);
colos = DISTINGUISHABLE_COLORS(3);

aa(1) = semilogx(exx.QS, exx.WS/2/pi, 'k-', 'LineWidth', 2); hold on
legend(aa(1), sprintf('Reference (%d DoFs)',size(exx.BBS{1}.U,1)-1), 'fontsize', 16);
k=1;
for i=Nlevs
    load(sprintf('./DATS/P11580_S-6.00_%s_%dLEV_BB_ELDRYFRICT_prop_wth.mat',sel_method,i), 'BB')
    load(sprintf('../WHOLEJOINT_ROM_PREPARE/REDMATS/P11580_S-6.00_%s_%dLEV_GRED_WJMAT_prop.mat',sel_method, i),'M')
    
    aa(k+1) = semilogx(BB.Q, BB.W/2/pi, symbs{k}, 'Color', colos(k,:), 'MarkerFaceColor', colos(k,:));
    legend(aa(k+1), sprintf('%d Level ROM (%d DoFs)',i,size(M,1)))
    k=k+1;
end
legend(aa, 'Location', 'southwest')
xlabel('Modal Amplitude')
ylabel('Natural Frequency (Hz)')
ylim([132 172])
xlim([1e-6 2e-3])

% print(sprintf('./FIGURES/%s_ELDRYFRICT_BB_WJROM_ZOOM.eps', sel_method), '-depsc')
export_fig(sprintf('./FIGURES/%s_ELDRYFRICT_BB_WJROM_ZOOM.eps', sel_method))

%% Mode Shape Deviations
figure(3)
clf()
set(gcf, 'Color', 'white')

aa = gobjects(4,1);
colos = DISTINGUISHABLE_COLORS(3);

k=1;
for i=Nlevs
    load(sprintf('./DATS/P11580_S-6.00_%s_%dLEV_BB_ELDRYFRICT_prop_wth.mat',sel_method,i), 'BB')
    load(sprintf('../WHOLEJOINT_ROM_PREPARE/REDMATS/P11580_S-6.00_%s_%dLEV_GRED_WJMAT_prop.mat',sel_method, i),'M', 'Th')
    MSrom = Th(1:size(Lref,1),:)*BB.Phi;
	MSrom = sign(Rref(3,:)*MSrom).*MSrom;
    
    semilogx(BB.Q, sqrt(2*(1-diag(MSrom'*Mref*MSref))), symbs{k}, 'Color', colos(k,:), 'MarkerFaceColor', colos(k,:)); hold on
    k=k+1;
end

xlabel('Modal Amplitude')
ylabel('Mode shape deviation')
% ylim([132 172])
xlim([1e-6 2e-3])

% print(sprintf('./FIGURES/%s_ELDRYFRICT_BB_WJROM_MAC.eps', sel_method), '-depsc')
export_fig(sprintf('./FIGURES/%s_ELDRYFRICT_BB_WJROM_MAC.eps', sel_method))
