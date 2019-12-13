clc
clear all
close all
addpath('./ROUTINES/export_fig')

set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex'); 
set(0, 'DefaultLegendInterpreter', 'latex'); 
set(0,'defaultAxesFontSize',13)

mi = 1;
%% Setup 
% FJROM
FJROM.sel_methods = {'U', 'P', 'PD'};
FJROM.Nelranges = {[36 74 122 232 448 568], [52 100 140 204 240 304], [68 108 152 200 252 292]};

% WJROM
WJROM.sel_methods = {'U', 'P', 'PD'};
WJROM.Nlevs = {[3 5 7 9 11 19 21], [2 3 4 5 6 10], [2 3 6 7 8 11]};

%% Load Reference Data
exx = load('./FULLJOINT_REF/DATS/RUN2.mat', 'QS', 'WS', 'ZS', 'DS','ZSM','BBS');

%% FJROM Errors
FJROM.Errs = cell(size(FJROM.sel_methods));

for m=1:length(FJROM.sel_methods)
    FJROM.Errs{m} = zeros(length(FJROM.Nelranges{m}), 4);
    for i=1:length(FJROM.Nelranges{m})
        load(sprintf('./FJROMSIMS/DATS/RUN_M%d_%s_%dELS.mat',mi,FJROM.sel_methods{m},FJROM.Nelranges{m}(i)), 'QS', 'WS', 'ZSM', 'BBS');
        
        FJROM.Errs{m}(i, 1) = size(BBS{1}.U,1);
        FJROM.Errs{m}(i, 2) = rms((WS-exx.WS)./exx.WS);
        FJROM.Errs{m}(i, 3) = rms((ZSM(5:end)-exx.ZS(5:end))./exx.ZS(5:end));
        FJROM.Errs{m}(i, 4) = rms((log10(ZSM(5:end))-log10(exx.ZS(5:end)))./log10(exx.ZS(5:end)));
    end
end

%% WJROM Errors
WJROM.Errs = cell(size(WJROM.sel_methods));

for m=1:length(WJROM.sel_methods)
    WJROM.Errs{m} = zeros(length(WJROM.Nlevs{m}), 4);
    for i=1:length(WJROM.Nlevs{m})
        load(sprintf('./WJSIMS/DATS/P11580_S-6.00_%s_%dLEV_BB_ELDRYFRICT.mat',WJROM.sel_methods{m},WJROM.Nlevs{m}(i)), 'BB')
        load(sprintf('./WHOLEJOINT_ROM_PREPARE/REDMATS/P11580_S-6.00_%s_%dLEV_GRED_WJMAT.mat',WJROM.sel_methods{m},WJROM.Nlevs{m}(i)),'M')
        
        WJROM.Errs{m}(i, 1) = size(M,1);
        WJROM.Errs{m}(i, 2) = rms((BB.W-exx.WS)./exx.WS);
        WJROM.Errs{m}(i, 3) = rms((BB.Z(5:end)-exx.ZS(5:end))./exx.ZS(5:end));
        WJROM.Errs{m}(i, 4) = rms((log10(BB.Z(5:end))-log10(exx.ZS(5:end)))./log10(exx.ZS(5:end)));
    end
end

%% Plots
figure(1)
clf()
set(gcf,'color','white')

symbs = {'s-', 'd-', 'p-', 'o-', 'v-', '^-'};
aa = gobjects(length(FJROM.sel_methods)+length(WJROM.sel_methods),1);
colos = DISTINGUISHABLE_COLORS(length(aa));

for m=1:length(WJROM.sel_methods)
    aa(m) = plot(WJROM.Errs{m}(:,1), WJROM.Errs{m}(:,2), symbs{m}, 'Color', colos(m,:), 'MarkerFaceColor', colos(m,:), 'LineWidth', 2); hold on
    legend(aa(m), sprintf('WJROM (%s)',WJROM.sel_methods{m}))
end
for m=1:length(FJROM.sel_methods)
    aa(length(WJROM.sel_methods)+m) = plot(FJROM.Errs{m}(:,1), FJROM.Errs{m}(:,2), symbs{length(WJROM.sel_methods)+m}, 'Color', colos(length(WJROM.sel_methods)+m,:), 'MarkerFaceColor', colos(length(WJROM.sel_methods)+m,:), 'LineWidth', 2); hold on
    legend(aa(length(WJROM.sel_methods)+m), sprintf('RMROM (%s)',FJROM.sel_methods{m}))
end
ll=legend(aa);
xlabel('Reduced Order Model Degrees-of-Freedom')
ylabel('Relative Frequency Error (RMS)')
set(gca, 'XScale', 'log')
% print('./CONVPLOT_W.eps', '-depsc')
export_fig('./CONVPLOT_W.eps')

figure(2)
clf()
set(gcf,'color','white')

aa = gobjects(length(FJROM.sel_methods)+length(WJROM.sel_methods),1);
k=1;
for m=1:length(WJROM.sel_methods)
    aa(m) = plot(WJROM.Errs{m}(:,1), WJROM.Errs{m}(:,4), symbs{m}, 'Color', colos(m,:), 'MarkerFaceColor', colos(m,:), 'LineWidth', 2); hold on
    legend(aa(m), sprintf('WJROM (%s)',WJROM.sel_methods{m}))
end
for m=1:length(FJROM.sel_methods)
    aa(length(WJROM.sel_methods)+m) = plot(FJROM.Errs{m}(:,1), FJROM.Errs{m}(:,4), symbs{length(WJROM.sel_methods)+m}, 'Color', colos(length(FJROM.sel_methods)+m,:), 'MarkerFaceColor', colos(length(WJROM.sel_methods)+m,:), 'LineWidth', 2); hold on
    legend(aa(length(WJROM.sel_methods)+m), sprintf('RMROM (%s)',FJROM.sel_methods{m}))
end
ll=legend(aa(1:end));
set(ll, 'Visible', 'off')
xlabel('Reduced Order Model Degrees-of-Freedom')
ylabel('Relative Log-Damping Factor Error (RMS)')
set(gca, 'XScale', 'log')
% print('./CONVPLOT_Z.eps', '-depsc')
export_fig('./CONVPLOT_Z.eps')