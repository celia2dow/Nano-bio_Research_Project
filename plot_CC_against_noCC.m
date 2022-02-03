% PLOT_CC_AGAINST_NOCC is a driver function that loads in data from
% previous runs with and without carrying capacity included, and plots the 
% association/internalisation curves against one another.
%
%   This is the work of Celia Dowling 19/01/22

clear
close all

fig22=figure(22);
set(fig22, 'Visible', 'off');
tol_diffCases = 1E-2;

% First load in all of the data for the case without CC
load('variables_1000pPerC_noCC_0.3_0.5.mat')
close(figure(1),figure(2),figure(3));

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
load('variables_1000pPerC_CC100_0.3_0.5.mat')
close(figure(1),figure(2),figure(3));

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
hold off

% Details of the plot
ylim([0 max(upper1,[],'all')])
legend('Internalised without CC','Interacting without CC','Associated (either) without CC',...
    'Std. dev. btwn all cells & runs without CC', '" "', '" "', ...
    'Internalised with CC','Interacting with CC','Associated (either) with CC',...
    'Std. dev. btwn all cells & runs with CC', '" "', '" "', ...
    'Location','Best')
title('Average number of particles associated per cell')
subtitle('With and without carrying capacity (CC)')
xlabel('time $t$ hours', 'Interpreter', 'latex');
ylabel('Average number of particles per cell');

% Save plot
fig22.Position = [100,100,1300,700];
date_time = num2str(fix(clock));
date_time = date_time(find(~isspace(date_time)));
savefig(fig22, [pwd '\With_or_without CC' date_time], 'compact')
saveas(fig22, [pwd '\With_or_without CC' date_time], 'png')

% Find the timestep at which the internalised curve deviates away from the
% case without CC
diff_btwn_cases = internalised_cases(1,:) - internalised_cases(2,:);
tsteps_similar = binrng(diff_btwn_cases <= tol_diffCases);
tstep_deviate = max(tsteps_similar);
disp(tstep_deviate);