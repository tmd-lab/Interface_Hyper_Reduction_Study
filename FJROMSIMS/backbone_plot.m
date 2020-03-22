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
% % Nels = 36; % [36 74 122 232 448 568 588]
% Nelrange = [74 122 232];

% sel_method = 'P';
% % Nels = 52;  % [52 100 140 204 240 304]
% Nelrange = [52 140 304];

sel_method = 'PD';
% Nels = 252;  % [68 108 152 200 252 292]
Nelrange = [68 152 292];

mi = 1;
%% Load Reference Data
exx = load('../FULLJOINT_REF/DATS/RUN2.mat', 'QS', 'WS', 'ZS', 'DS','ZSM','BBS');
ref = load('../FULLJOINT_REF/PREP/5_SET_NULLRED.mat','Mrel','L','Rrel');
Mref = ref.Mrel;
Lref = ref.L;
Rref = ref.Rrel;
MSref = Lref*cell2mat(cellfun(@(c) c.MS(:,end), exx.BBS, 'UniformOutput', false));
MSref = sign(Rref(3,:)*MSref).*MSref;
%% Plot Backbones
figure(10);
clf()
set(gcf,'color','white')

colos = DISTINGUISHABLE_COLORS(7);  colos = colos(end-2:end,:);
aa = gobjects(4,1);

ax1=subplot('Position', [0.13 0.5 0.8 0.45]);
set(ax1, 'DataAspectRatioMode', 'manual');
aa(1) = semilogx(exx.QS, exx.WS/2/pi, 'k-', 'LineWidth', 2); hold on
legend(aa(1), sprintf('Reference (%d DoFs)',size(exx.BBS{1}.U,1)-1))

ax2=subplot('Position', [0.13 0.12 0.8 0.33]);
loglog(exx.QS, exx.ZSM, 'k-', 'LineWidth', 2); hold on

k=1;
for ii=1:length(Nelrange)
    Nels = Nelrange(ii);
%     load(sprintf('./DATS/RUN_%d.mat',p), 'QS', 'WS', 'ZSM', 'BBS');
    load(sprintf('./DATS/RUN_M%d_%s_%dELS.mat',mi,sel_method,Nels), 'QS', 'WS', 'ZSM', 'BBS');
    
    subplot(ax1);
    aa(k+1) = semilogx(QS, WS/2/pi, symbs{k}, 'Color', colos(k,:), 'MarkerFaceColor', colos(k,:));
    legend(aa(k+1), sprintf('%d Element ROM (%d DoFs)',Nels,size(BBS{1}.U,1)-1))
    
    subplot(ax2);
    loglog(QS, ZSM, symbs{k}, 'Color', colos(k,:), 'MarkerFaceColor', colos(k,:))
    
    k = k+1;
end
subplot(ax1)
xticklabels([])
ll=legend(aa(1:end), 'Location', 'southwest');
set(ll, 'visible', 'off')
ylabel('Natural Frequency (Hz)')
xlim([1e-6 1e-1])

subplot(ax2)
xlabel('Modal Amplitude')
ylabel('Damping factor')
xlim([1e-6 1e-1])
ylim([9e-6 1e0])

% print(sprintf('./FIGURES/ROMBB_M%d_%s.eps',mi,sel_method), '-depsc')
export_fig(sprintf('./FIGURES/ROMBB_M%d_%s.eps',mi,sel_method))
%% Zoomed Plot
figure(20)
clf()
set(gcf,'color','white')

aa = gobjects(4,1);
colos = DISTINGUISHABLE_COLORS(7);  colos = colos(end-2:end,:);

aa(1) = semilogx(exx.QS, exx.WS/2/pi, 'k-', 'LineWidth', 2); hold on
legend(aa(1), sprintf('Reference (%d DoFs)',size(exx.BBS{1}.U,1)-1), 'fontsize', 16)
k=1;
for ii=1:length(Nelrange)
    Nels = Nelrange(ii);
%     load(sprintf('./DATS/RUN_%d.mat',p), 'QS', 'WS', 'ZSM', 'BBS');
    load(sprintf('./DATS/RUN_M%d_%s_%dELS.mat',mi,sel_method,Nels), 'QS', 'WS', 'ZSM', 'BBS');
    
    aa(k+1) = semilogx(QS, WS/2/pi, symbs{k}, 'Color', colos(k,:), 'MarkerFaceColor', colos(k,:));
    legend(aa(k+1), sprintf('%d Element ROM (%d DoFs)',Nels,size(BBS{1}.U,1)-1))
    k = k+1;
end
legend(aa(1:end), 'Location', 'southwest')
xlabel('Modal Amplitude')
ylabel('Natural Frequency (Hz)')
ylim([132 172])
xlim([1e-6 2e-3])

% print(sprintf('./FIGURES/ROMBB_M%d_%s_ZOOM.eps',mi,sel_method), '-depsc')
export_fig(sprintf('./FIGURES/ROMBB_M%d_%s_ZOOM.eps',mi,sel_method))

%% Mode Shape Deviations
figure(30)
clf()
set(gcf, 'Color', 'white')

aa = gobjects(4,1);
colos = DISTINGUISHABLE_COLORS(7);  colos = colos(end-2:end,:);

k=1;
for ii=1:length(Nelrange)
    Nels = Nelrange(ii);
    load(sprintf('./DATS/RUN_M%d_%s_%dELS.mat',mi,sel_method,Nels), 'QS', 'WS', 'ZSM', 'BBS');
    load(sprintf('../FULLJOINT_ROM_PREPARE/ROMS/ROM_%s_%dELS.mat',sel_method,Nels), 'TFMh');
    MSrom = TFMh*cell2mat(cellfun(@(c) c.MS(:,end), BBS, 'UniformOutput', false));
    MSrom = sign(Rref(3,:)*MSrom).*MSrom;
    
    semilogx(QS, sqrt(2*(1-diag(MSrom'*Mref*MSref))), symbs{k}, 'Color', colos(k,:), 'MarkerFaceColor', colos(k,:)); hold on
    k=k+1;
end

xlabel('Modal Amplitude')
ylabel('Mode shape deviation')
% ylim([132 172])
xlim([1e-6 2e-3])

% print(sprintf('./FIGURES/ROMBB_M%d_%s_MAC.eps',mi,sel_method), '-depsc')
export_fig(sprintf('./FIGURES/ROMBB_M%d_%s_MAC.eps',mi,sel_method))