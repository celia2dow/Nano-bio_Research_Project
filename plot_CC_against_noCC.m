% PLOT_CC_AGAINST_NOCC is a driver function that loads in data from
% previous runs with and without carrying capacity included, and plots the 
% association/internalisation curves against one another.
%
%   This is the work of Celia Dowling 04/02/22

clear
close all

fig22=figure(22);
set(fig22, 'Visible', 'off');
tol_diffCases = 1E-2;
choose = 'distrib'; %'diffs'; %

% Hypoexponential CDF - Internalised 
F_hypoexp = @(t,lambda1,lambda2,num_prtcls) (1 - 1./(lambda2-lambda1) .* ...
    (lambda2 .* exp(-lambda1 .* t) - ...
    lambda1 .* exp(-lambda2 .* t))) .* num_prtcls;

% Hypoexponential PDF - Internalised
f_hypoexp = @(t,lambda1,lambda2,num_prtcls) lambda1 .* lambda2 .* num_prtcls ...
    ./(lambda1 - lambda2) .* (exp(-lambda2 .* t) - exp(-lambda1 .* t));

% Exponential CDF - Associated
F_exp = @(t,lambda1,num_prtcls) (1 - exp(-lambda1 .* t)) .* num_prtcls;

% Difference CDF - Interacting
F_diff = @(t,lambda1,lambda2,num_prtcls) lambda1 .* num_prtcls ./(lambda2-lambda1) .* ...
    (exp(-lambda1 .* t) - exp(-lambda2 .* t));

% First load in all of the data for the case without CC
load('variables_1000pPerC_CCInf_0.05_0.03_20222484436.mat')
close(figure(1),figure(2),figure(3));

switch choose 
    case {'distrib'}
        l1_noCC = mean(est_lambda1.MLE); 
        l2_noCC = est_lambda2.distrib_mean; 
    case {'diffs'}
        l1_noCC = mean(est_lambda1.using_diffs);
        l2_noCC = est_lambda2.diffs_mean;
end

% Save thhe average internalised numbers without CC
internalised_cases = [means(2,:); zeros(1,length(means(2,:)))];

% Plot the average number of interacting, internalised and associated
% particles per cell for the case without CC
figure(22);
plot(binrng,means(2,:),'r',binrng,means(1,:),'b');
hold on
plot(binrng,means(3,:),'Color',[.5 0 .5])
patch([binrng fliplr(binrng)], [upper1(1,:) fliplr(lower1(1,:))], [0 0 1], 'FaceAlpha', 0.1, 'EdgeColor', 'w', 'LineStyle', ':')
patch([binrng fliplr(binrng)], [upper1(2,:) fliplr(lower1(2,:))], [1 0 0], 'FaceAlpha', 0.1, 'EdgeColor', 'w', 'LineStyle', ':')
patch([binrng fliplr(binrng)], [upper1(3,:) fliplr(lower1(3,:))], [.5 0 .5], 'FaceAlpha', 0.1, 'EdgeColor', 'w', 'LineStyle', ':')
hold off

% Second load in all of the data for the case with CC
load('variables_1000pPerC_CC30_0.05_0.03_202223225522.mat')
close(figure(1),figure(2),figure(3));

switch choose 
    case {'distrib'}
        l1_CC = mean(est_lambda1.MLE); 
        l2_CC = est_lambda2.distrib_mean; 
        CC = mean(est_CC.using_lambda2_diffs(12/PARAMETERS.tstep_duration + 1:end));
    case {'diffs'}
        l1_CC = mean(est_lambda1.using_diffs);
        l2_CC = est_lambda2.diffs_mean;
        CC = mean(est_CC.using_diffs(12/PARAMETERS.tstep_duration + 1:end));
end

% Save thhe average internalised numbers with CC
internalised_cases(2,:) = means(2,:);

% Plot the average number of interacting, internalised and associated
% particles per cell for the case with CC
figure(22)
hold on
plot(binrng,means(2,:),'Color',[0.9 0 0.6],'LineStyle', '--')
plot(binrng,means(1,:),'Color',[0.2 0.5 1],'LineStyle', '--')
plot(binrng,means(3,:),'Color',[0.6 0 0.9],'LineStyle', '--')
patch([binrng fliplr(binrng)], [upper1(1,:) fliplr(lower1(1,:))], [0.2 0.5 1], 'FaceAlpha', 0.1, 'EdgeColor', 'w', 'LineStyle', ':')
patch([binrng fliplr(binrng)], [upper1(2,:) fliplr(lower1(2,:))], [0.9 0 0.6], 'FaceAlpha', 0.1, 'EdgeColor', 'w', 'LineStyle', ':')
patch([binrng fliplr(binrng)], [upper1(3,:) fliplr(lower1(3,:))], [0.6 0 0.9], 'FaceAlpha', 0.1, 'EdgeColor', 'w', 'LineStyle', ':')


% Plot the theoretical distributions without CC
assoc_noCC = F_exp(binrng,l1_noCC,PARAMETERS.prtcls_per_cell);
internal_noCC = F_hypoexp(binrng,l1_noCC,l2_noCC,PARAMETERS.prtcls_per_cell);
interact_noCC = F_diff(binrng,l1_noCC,l2_noCC,PARAMETERS.prtcls_per_cell);
plot(binrng,assoc_noCC,'y');
plot(binrng,internal_noCC,'y');
plot(binrng,interact_noCC,'y');

% Plot the theoretical distributions with CC
assoc_CC = F_exp(binrng,l1_CC,PARAMETERS.prtcls_per_cell);
internal_CC = zeros(1,length(binrng));
pdf = zeros(1,length(binrng));
pdf(1) = F_hypoexp(0,l1_CC,l2_CC,PARAMETERS.prtcls_per_cell);
internal_CC(1) = pdf(1);
for i = 2:length(binrng)
    l2_dyn = l2_CC*(1-internal_CC(i-1)/CC);
    pdf(i) = f_hypoexp(binrng(i),l1_CC,l2_dyn,PARAMETERS.prtcls_per_cell);
    internal_CC(i) = internal_CC(i-1)+pdf(i).*PARAMETERS.tstep_duration;
end
%internal_CC = F_hypoexp(binrng,l1_CC,est_lambda2.distrib,PARAMETERS.prtcls_per_cell);
interact_CC = assoc_CC - internal_CC;
plot(binrng,assoc_CC,'y--');
plot(binrng,internal_CC,'y--');
plot(binrng,interact_CC,'y--');

% Details of the plot
ylim([0 max(upper1,[],'all')])
legend('Internalised without CC','Interacting without CC','Associated (either) without CC',...
    'Std. dev. btwn all cells & runs without CC', '" "', '" "', ...
    'Internalised with CC','Interacting with CC','Associated (either) with CC',...
    'Std. dev. btwn all cells & runs with CC', '" "', '" "', ...
    'Theoretical distributions from estimates without CC', '" "', '" "', ...
    'Theoretical distributions from estimates with CC', '" "', '" "', ...
    'Location','Best')
title('Average number of particles associated per cell')
subtitle('With and without carrying capacity (CC)')
xlabel('time $t$ hours', 'Interpreter', 'latex');
ylabel('Average number of particles per cell');

% Save plot
fig22.Position = [100,100,1300,700];
date_time = num2str(fix(clock));
date_time = date_time(find(~isspace(date_time)));
savefig(fig22, [pwd '\With_or_without CC' num2str(PARAMETERS.max_prtcls(end)) choose date_time], 'compact')
saveas(fig22, [pwd '\With_or_without CC' num2str(PARAMETERS.max_prtcls(end)) choose date_time], 'png')

% Find the timestep at which the internalised curve deviates away from the
% case without CC
diff_btwn_cases = internalised_cases(1,:) - internalised_cases(2,:);
tsteps_similar = binrng(diff_btwn_cases <= tol_diffCases);
tstep_deviate = max(tsteps_similar);
fprintf(['The time (hours) at which the internalised curve deviates away'...
    ' from the case without CC by more than ' num2str(tol_diffCases) ' particles: \n'])
disp(tstep_deviate);

% Print the estimates for lambda and the carrying capacity
format short e
fprintf('Actual: \n\tlambda1 \tlambda2 \t\tCC \n')
disp([l1, l2, PARAMETERS.max_prtcls(end)])
fprintf(['Using method: ' choose '\n\n'])
fprintf('Without CC: \n\tlambda1 \tlambda2 \n')
disp([l1_noCC, l2_noCC])
fprintf('Without CC: \n\tlambda1 \tlambda2 \t\tCC \n')
disp([l1_CC, l2_CC, CC])
