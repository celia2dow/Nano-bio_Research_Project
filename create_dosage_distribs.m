function create_dosage_distribs(X, total_tsteps, PARAMETERS, total)
% CREATE_DOSAGE_DISTRIBS crreates dosage distributions after every X hours.

% Find limits for the x and y-axis: the largest values to be plotted, given
% that plots are being plotted every hour
Xhour_indices = 0:floor(X/PARAMETERS.tstep_duration):total_tsteps;
Xhour_indices = Xhour_indices(2:end) + 1;
hourly_total.cell_c_o_p = total.cell_c_o_p(:,Xhour_indices,:);
interact_max = max(hourly_total.cell_c_o_p(:,:,2),[],'all');
internal_max = max(hourly_total.cell_c_o_p(:,:,3),[],'all');
x_max = max(interact_max,internal_max);
local_max = zeros(1,length(Xhour_indices));
for index = 2:length(Xhour_indices)
    [cells_interact,~] = histcounts(hourly_total.cell_c_o_p(:,index,2));
    [cells_internal,~] = histcounts(hourly_total.cell_c_o_p(:,index,3));
    if ~isempty(cells_interact) && ~isempty(cells_internal)
        local_max(index) = max([cells_interact(2:end), cells_internal(2:end)]);
    elseif ~isempty(cells_interact)
        local_max(index) = max(cells_interact(2:end));
    elseif ~isempty(cells_internal)
        local_max(index) = max(cells_internal(2:end));
    end
end
y_max = max([local_max,PARAMETERS.initial_num_cells]);

fig2 = figure(2);
set(fig2, 'Visible', 'off');

for time_plot = 1:length(Xhour_indices)
    tstep = (time_plot-1)/PARAMETERS.tstep_duration; % equivalent timestep index
    N_tstep = total.cell_population(tstep+1); % total number of cells across runs
    
    % FREQUENCY HISTOGRAMS
    set(0,'CurrentFigure',fig2)
    dim1 = floor(sqrt(length(Xhour_indices)));
    dim2 = ceil(length(Xhour_indices)/dim1);
    subplot(dim1,dim2,time_plot);
    % FREQUENCY OF CELLS WITH NUMS OF PARTICLES INTERACTING OVER TIME
    histogram(hourly_total.cell_c_o_p(1:N_tstep,time_plot,2),...
        'FaceColor', [0,0,1], 'FaceAlpha', 0.2);
    hold on;
    % FREQUENCY OF CELLS WITH NUMS OF PARTICLES INTERNALISED OVER TIME
    histogram(hourly_total.cell_c_o_p(1:N_tstep,time_plot,3),...
        'FaceColor', [1,0,0], 'FaceAlpha', 0.2);
    hold off;
    xlim([0,x_max]);
    %ylim([0,y_max]);
    title(['At ' num2str(time_plot*X) ' hours']);
    xlabel('Number of particles');
    ylabel('Cell frequency');
    pbaspect([1 1 1]);
    if time_plot==1
        legend('Interacting','Internalised');
    end
end

% Save figure
sgtitle(fig2, 'Frequency of cells with certain numbers of interacting/internalised particles over time')
fig2.Position = [100,100,1100,700];
saveas(fig2, [PARAMETERS.folder_path '/Dosage_distribs' num2str(X) 'hours'], 'eps')
saveas(fig2, [PARAMETERS.folder_path '/Dosage_distribs' num2str(X) 'hours'], 'png')
figure(fig2)
