% SCRIPT plotting the mean Pair Correlation Coefficient over time (the mean 
% PCC at each timestep across all of the run PCCs and the variance in the mean)
fig8=figure(8);
set(fig8, 'Visible', 'off');
yyaxis left
plot(binrng,mean(runs.cell_pair_cor_coef,1))
ylabel('Mean PCC of runs')
up = mean(runs.cell_pair_cor_coef,1) + var(runs.cell_pair_cor_coef,1);
low = mean(runs.cell_pair_cor_coef,1) - var(runs.cell_pair_cor_coef,1);
hold on
patch([binrng fliplr(binrng)], [up fliplr(low)], [0 0 1], 'FaceAlpha', 0.2, 'EdgeColor', 'w', 'LineStyle', '-')
% Plot the average change in confluence over time
yyaxis right
plot(binrng, total.cell_population./(PARAMETERS.culture_dim^2 * num_runs).*100)
ylabel('% confluence')
hold off
xlabel('time $t$ hours', 'Interpreter', 'latex')
title('Pair correlation coefficient against dish confluence over time')
if PARAMETERS.EWT_move == inf
    subtitle('Without motility events')
else
    subtitle('With motility events')
end
fig8.Position = [100,100,1300,700];
saveas(fig8, [PARAMETERS.folder_path '/PCC_and_confluence'], 'eps')
saveas(fig8, [PARAMETERS.folder_path '/PCC_and_confluence'], 'png')