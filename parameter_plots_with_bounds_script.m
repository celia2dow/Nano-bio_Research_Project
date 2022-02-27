% SCRIPT for creating plots of parameter estimates and their lower and
% upper bounds
clf(figure(27))
clf(figure(28))
clf(figure(29))
clf(figure(30))

% Lambda 1 values
figure(27)
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

% Lambda 2 values
figure(28)
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

% CCvalues
figure(29)
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

% Comparing data to estimate rates:
% Data:
figure(30);
plot(binrng,means(2,:),'r',binrng,means(1,:),'b');
hold on
plot(binrng,means(3,:),'Color',[.5 0 .5])
patch([binrng fliplr(binrng)], [upper1(1,:) fliplr(lower1(1,:))], [0 0 1], 'FaceAlpha', 0.2, 'EdgeColor', 'w', 'LineStyle', '-')
patch([binrng fliplr(binrng)], [upper1(2,:) fliplr(lower1(2,:))], [1 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'w', 'LineStyle', '-')
patch([binrng fliplr(binrng)], [upper1(3,:) fliplr(lower1(3,:))], [.5 0 .5], 'FaceAlpha', 0.2, 'EdgeColor', 'w', 'LineStyle', '-')
plot(binrng,lower2(1,:),'b--',binrng,upper2(1,:),'b--', ...
    binrng,lower2(2,:),'r--',binrng,upper2(2,:),'r--');
plot(binrng,lower2(3,:),'Color',[.5 0 .5],'LineStyle', '--')
plot(binrng,upper2(3,:),'Color',[.5 0 .5],'LineStyle', '--')

% Use MLE for lambda1, Hypoexponential distribution for lambda2, and
% dynamic via differences for CC.
lam1 = est_lambda1.MLE; lam1UP = est_lambda1UP.MLE; lam1LO = est_lambda1LO.MLE;
lam2 = est_lambda2.MLE; lam2UP = est_lambda2UP.MLE; lam2LO = est_lambda2LO.MLE;
CC = est_CC.MLE; CCUP = est_CCUP.MLE; CCLO = est_CCLO.MLE;
intern = hy
lam2_dyn = lam2*(1-(i-1)/CC);

plot(binrng,CDF.exp(lam1,binrng),'y:', 'LineWidth', 1.2);
plot(binrng,CDF.exp(lam1LO.MLE,binrng),'y:', 'LineWidth', 1.2);
plot(binrng,CDF.exp(lam1UP.MLE,binrng),'y:', 'LineWidth', 1.2);

ylim([0 max(upper1,[],'all')])
legend('Internalised','Interacting','Associated (either)',...
    'Std. dev. btwn all cells & runs', '" "', '" "', ...
    'Std. dev. btwn run means', '','" "','', '" "', 'Location','Best')
title('Average number of particles associated per cell: separate trends')
xlabel('time $t$ hours', 'Interpreter', 'latex');


