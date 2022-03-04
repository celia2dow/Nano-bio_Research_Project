% SCRIPT for creating plots of parameter estimates and their lower and
% upper bounds
clf(figure(27))
clf(figure(28))
clf(figure(29))
clf(figure(30))

% Lambda 1 values
fig27 = figure(27);
set(fig27, 'Visible', 'off');
% MLE
plot(binrng(1:end-1), est_lambda1.MLE, 'k--', 'LineWidth', 1.5)
hold on
plot(binrng(1:end-1), est_lambda1UP.MLE, 'k:')
plot(binrng(1:end-1), est_lambda1LO.MLE, 'k:')
% Differences method
plot(binrng(1:end-1), est_lambda1.using_diffs, 'r--', 'LineWidth', 1.5)
plot(binrng(1:end-1), est_lambda1UP.using_diffs, 'r:')
plot(binrng(1:end-1), est_lambda1LO.using_diffs, 'r:')
hold off
title('\lambda_1 calculated via MLE (black) and via differences (red)')

fig27.Position = [100,100,1300,700];
saveas(fig27, [PARAMETERS.folder_path '/lambda1_eror'], 'eps')
saveas(fig27, [PARAMETERS.folder_path '/lambda1_eror'], 'png')

% Lambda 2 values
fig28 = figure(28);
set(fig28, 'Visible', 'off');
% Differences method
plot(binrng(2:end-1),est_lambda2.using_diffs, 'r--', 'LineWidth', 1.5)
hold on
plot(binrng(2:end-1),est_lambda2UP.using_diffs, 'r:')
plot(binrng(2:end-1),est_lambda2LO.using_diffs, 'r:')
% Mix method
plot(binrng(1:end-1),est_lambda2.mix, 'b--', 'LineWidth', 1.5)
plot(binrng(1:end-1),est_lambda2UP.mix, 'b:')
plot(binrng(1:end-1),est_lambda2LO.mix, 'b:')
% Disribution method
plot(binrng,est_lambda2.distrib, 'g--', 'LineWidth', 1.5)
plot(binrng,est_lambda2UP.distrib, 'g:')
plot(binrng,est_lambda2LO.distrib, 'g:')
ylim([min(est_lambda2LO.distrib), max(est_lambda2UP.distrib)])
hold off
title('\lambda_2 calculated via differences (red), distribution (green) and mix (blue)')

fig28.Position = [100,100,1300,700];
saveas(fig28, [PARAMETERS.folder_path '/lambda2_eror'], 'eps')
saveas(fig28, [PARAMETERS.folder_path '/lambda2_eror'], 'png')

% CCvalues
fig29 = figure(29);
set(fig29, 'Visible', 'off');
if any(PARAMETERS.max_prtcls ~= inf) 
    % Differences method
    plot(binrng(1:end-1),est_CC.using_diffs, 'r--', 'LineWidth', 1.5)
    hold on
    plot(binrng(1:end-1),est_CCUP.using_diffs, 'r:')
    plot(binrng(1:end-1),est_CCLO.using_diffs, 'r:')
    % Mix method
    plot(binrng(1:end-1),est_CC.dynamic_diffs, 'm--', 'LineWidth', 1.5)
    plot(binrng(1:end-1),est_CCUP.dynamic_diffs, 'm:')
    plot(binrng(1:end-1),est_CCLO.dynamic_diffs, 'm:')
    % Disribution method
    plot(binrng(1,end-1),est_CC.dynamic_mix, 'c--', 'LineWidth', 1.5)
    plot(binrng(1:end-1),est_CCUP.dynamic_mix, 'c:')
    plot(binrng(1:end-1),est_CCLO.dynamic_mix, 'c:')
    ylim([0 2*PARAMETERS.max_prtcls(end)])
    hold off
    title('CC calculated via differences (red), dynamic differences (magenta) and dynamic mix (cyan)')

    fig29.Position = [100,100,1300,700];
    saveas(fig29, [PARAMETERS.folder_path '/CC_eror'], 'eps')
    saveas(fig29, [PARAMETERS.folder_path '/CC_eror'], 'png')
end

% Comparing data to estimate rates:
% Data and one standard deviation error:
fig30 = figure(30);
set(fig30, 'Visible', 'off');
p1 = patch([binrng fliplr(binrng)], [upper1(1,:) fliplr(lower1(1,:))], [0 0 1], ...
    'FaceAlpha', 0.2, 'EdgeColor', 'w', 'LineStyle', '-', 'DisplayName','1 standard deviation');
hold on
p2 = patch([binrng fliplr(binrng)], [upper1(2,:) fliplr(lower1(2,:))], [1 0 0], ...
    'FaceAlpha', 0.2, 'EdgeColor', 'w', 'LineStyle', '-', 'DisplayName','" " "');
p3 = patch([binrng fliplr(binrng)], [upper1(3,:) fliplr(lower1(3,:))], [.5 0 .5], ...
    'FaceAlpha', 0.2, 'EdgeColor', 'w', 'LineStyle', '-', 'DisplayName','" " "');
p4 = plot(binrng,means(2,:),'r', 'DisplayName','Mean internalised');
p5 = plot(binrng,means(1,:),'b', 'DisplayName','Mean bound');
p6 = plot(binrng,means(3,:),'Color',[.5 0 .5], 'DisplayName','Mean associated (either)');

% Use MLE for lambda1, Hypoexponential distribution for lambda2, and
% dynamic via differences for CC.
lam1 = est_lambda1.MLE_mean; lam1UP = est_lambda1UP.MLE_mean; lam1LO = est_lambda1LO.MLE_mean;
lam2 = est_lambda2.distrib_mean; lam2UP = est_lambda2UP.distrib_mean; lam2LO = est_lambda2LO.distrib_mean;
% Association curves: mean, lower, upper
assoc = PARAMETERS.prtcls_per_cell .* CDF.exp(lam1,binrng);
assocLO = PARAMETERS.prtcls_per_cell .* CDF.exp(lam1LO,binrng);
assocUP = PARAMETERS.prtcls_per_cell .* CDF.exp(lam1UP,binrng);

if any(PARAMETERS.max_prtcls ~= inf) 
    CC = est_CC.dynamic_diffs_mean; CCUP = est_CCUP.dynamic_diffs_mean; CCLO = est_CCLO.dynamic_diffs_mean;
    interns = analytic_distrib_with_CC(binrng,lam2,CC,PARAMETERS, ...
        CDF.hypoexp_MLEl1,CDF.exp,assoc);
    internsLO = analytic_distrib_with_CC(binrng,lam2LO,CC,PARAMETERS, ...
        CDF.hypoexp_MLEl1UP,CDF.exp,assocUP);
    internsUP = analytic_distrib_with_CC(binrng,lam2UP,CC,PARAMETERS, ...
        CDF.hypoexp_MLEl1LO,CDF.exp,assocLO);
    sub_str = ['Heuristic estimates: \lambda_1: ' num2str(est_lambda1.MLE_mean)...
        ', \lambda_2: ' num2str(est_lambda2.distrib_mean)...
        ', CC: ' num2str(est_CC.dynamic_diffs_mean)];
    tit_str = ['Associated, bound, internalised per cell ~ true values:'...
        ' \lambda_1: ' num2str(l1) ...
        ', \lambda_2: ' num2str(l2) ...
        ', CC: ' num2str(PARAMETERS.max_prtcls(end))];
else 
    % Internalisation curves: mean, lower, upper
    interns = PARAMETERS.prtcls_per_cell .* CDF.hypoexp_MLEl1(lam2,binrng);
    internsLO = PARAMETERS.prtcls_per_cell .* CDF.hypoexp_MLEl1UP(lam2LO,binrng);
    internsUP = PARAMETERS.prtcls_per_cell .* CDF.hypoexp_MLEl1LO(lam2UP,binrng);
    sub_str = ['Heuristic estimates: \lambda_1: ' num2str(est_lambda1.MLE_mean)...
        ', \lambda_2: ' num2str(est_lambda2.distrib_mean)];
    tit_str = ['Associated, bound, internalised per cell ~ true values: '...
        ' \lambda_1: ' num2str(l1) ...
        ', \lambda_2: ' num2str(l2)];
end
% Interacting curves: mean, lower, upper
interacts = assoc - interns;
interactsLO = assoc - internsUP;
interactsUP = assoc - internsLO;

% Plot ALL
p7 = plot(binrng,assoc,'k--', 'LineWidth', 1.2,...
    'DisplayName', 'Mean analytic distributions with heuristic estimates');
p8 = plot(binrng,assocLO,'k:', 'LineWidth', 1.2,...
    'DisplayName', 'Lower/upper bound analytic distributions with heuristic estimates');
p9 = plot(binrng,assocUP,'k:', 'LineWidth', 1.2);
hasbehavior(p9,'legend',false);
p10 = plot(binrng,interns,'k--', 'LineWidth', 1.2);
hasbehavior(p10,'legend',false);
p11 = plot(binrng,internsLO,'k:', 'LineWidth', 1.2);
hasbehavior(p11,'legend',false);
p12 = plot(binrng,internsUP,'k:', 'LineWidth', 1.2);
hasbehavior(p12,'legend',false);
p13 = plot(binrng,interacts,'k--', 'LineWidth', 1.2);
hasbehavior(p13,'legend',false);
p14 = plot(binrng,interactsLO,'k:', 'LineWidth', 1.2);
hasbehavior(p14,'legend',false);
p15 = plot(binrng,interactsUP,'k:', 'LineWidth', 1.2);
hasbehavior(p15,'legend',false);

legend('Location','Best')
ylim([0 max(upper1,[],'all')])
title(tit_str);
subtitle(sub_str);
xlabel('time $t$ hours', 'Interpreter', 'latex');

fig30.Position = [100,100,1300,700];
pbaspect([1 1 1])
saveas(fig30, [PARAMETERS.folder_path '/Analytical_distribs'], 'eps')
saveas(fig30, [PARAMETERS.folder_path '/Analytical_distribs'], 'png')
