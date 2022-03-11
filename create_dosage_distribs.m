function create_dosage_distribs(X, total_tsteps, PARAMETERS, total, FLUORESC)
% CREATE_DOSAGE_DISTRIBS creates dosage distributions, fluorescence plots, 
% and dosages split by number of cell divisions after every X hours.

% Find the hours on which plots will be made
Xhour_indices = 0:floor(X/PARAMETERS.tstep_duration):total_tsteps;
Xhour_indices = Xhour_indices(2:end) + 1;
dim1 = floor(sqrt(length(Xhour_indices)));
dim2 = ceil(length(Xhour_indices)/dim1);

% Find limits for the x and y-axis: the largest values to be plotted, given
% that plots are being plotted every hour
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

% Find the limits of each stain fluorescence 
S1_max = FLUORESC.stain1_part * ((interact_max + internal_max) + ...
    FLUORESC.alpha*(interact_max + internal_max)) + FLUORESC.stain1_background;
S2_max = FLUORESC.stain2_part * (interact_max) + ...
    FLUORESC.stain2_background;
disp([S1_max, S2_max])

fig2 = figure(2);
set(fig2, 'Visible', 'off');
fig3 = figure(3);
set(fig3, 'Visible', 'off');
fig9 = figure(9);
set(fig9, 'Visible', 'off');

for time_plot = 1:length(Xhour_indices)
    tstep = (time_plot-1)/PARAMETERS.tstep_duration; % equivalent timestep index
    N_tstep = total.cell_population(tstep+1); % total number of cells across runs
    
    % FREQUENCY HISTOGRAMS
    set(0,'CurrentFigure',fig2)
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

    % FLUORESCENCE DOT PLOTS
    set(0,'CurrentFigure',fig3)
    subplot(dim1,dim2,time_plot);
    % PLOT EXTERNAL (STAIN 2) FLUORESCENCE AGAINST ASSOCIATED (STAIN 1)
    % FLUORESCENCE
    n_interact = total.cell_c_o_p(1:N_tstep,time_plot,2);
    n_internal = total.cell_c_o_p(1:N_tstep,time_plot,3);
    mean1 = FLUORESC.stain1_part .* (n_interact + FLUORESC.alpha.*n_internal) + FLUORESC.stain1_background; % true signal from stain 1 (x-axis)
    mean2 = FLUORESC.stain2_part .* n_interact + FLUORESC.stain2_background; % true signal from stain 2 (y-axis)
    S1 = sqrt(mean1) .* FLUORESC.stain1_std_dev .* randn(N_tstep,1) + mean1; % received signal from stain 1
    S2 = sqrt(mean2) .* FLUORESC.stain2_std_dev .* randn(N_tstep,1) + mean2; % received signal from stain 2
    S1c = S1; S2c = S2; % Specifically for working out the colour densities
    S1c(S1c<=0)=realmin; % So that 0 is not passed through log
    S2c(S2c<=0)=realmin; % So that 0 is not passed through log
    c = ksdensity([log10(S1c),log10(S2c)], [log10(S1c),log10(S2c)]);
    dot_size = 3*ones(length(S1),1);
    scatter(S1, S2, dot_size, c, 'filled');
    set(gca, 'YScale', 'log','Xscale', 'log','XMinorTick','on','YMinorTick','on');
    xlim([1,S1_max]);
    ylim([1,S2_max]);
    title(['At ' num2str(time_plot*X) ' hour/s']);
    xlabel('Stain 1 (associated)');
    ylabel('Stain 2 (external)');
    cb = colorbar();
    cb.Label.String = 'Density estimate';

%     % FREQUENCY HISTOGRAMS PER NUM OF CELL DIVISIONS EVERY X HOURS
%     set(0,'CurrentFigure',fig9)
%     cell_division_class = zeros(1,N_tstep);
%     if N_tstep>PARAMETERS.initial_num_cells
%         division_history = EVOLUTION_INFO.cell_lineage_history(PARAMETERS.initial_num_cells+1:N_tstep,1);
%         edges = (0.5:1:max(division_history)+0.5);
%         cell_division_class = cell_division_class(1:max(division_history)) + histcounts(division_history',edges);
%     end
%     rows = unique(cell_division_class);
%     
%     for row = rows % ITERATING THROUGH THE DIVISION NUMBERS PRESENT AT THIS TIMESTEP
%         plot_num = row*(floor(PARAMETERS.simulation_duration/6) + 1) + colmn;
%         subplot(oldest_cell_gen,floor(PARAMETERS.simulation_duration/6) + 1,plot_num);
%         % FREQUENCY OF CELLS WITH NUMS OF PARTICLES INTERACTING IN THIS
%         % CLASS OF CELL DIVISIONS
%         has_this_many_divs = (cell_division_class==row);
%         cells_with_divs = EVOLUTION_INFO.cell_lineage_history(has_this_many_divs,2);
%         histogram(EVOLUTION_INFO.cell_c_o_p(cells_with_divs,time_plot,2),...
%             'FaceColor', [0,0,1], 'FaceAlpha', 0.2);
%         hold on;
%         % FREQUENCY OF CELLS WITH NUMS OF PARTICLES INTERNALISED IN
%         % THIS CLASS OF CELL DIVISIONS
%         histogram(EVOLUTION_INFO.cell_c_o_p(cells_with_divs,time_plot,3),...
%             'FaceColor', [1,0,0], 'FaceAlpha', 0.2);
%         hold off;
%         xlim([0,x_max]);
%         title(['At ' num2str(time_plot-1) ' hour/s, cells with ' num2str(row) ' divisions']);
%         xlabel('Num. of particles');
%         ylabel('Cell frequ.');
%         if colmn==1
%             legend('Interacting','Internalised');
%         end
%     end
end

% Save figures
sgtitle(fig2, 'Frequency of cells with certain numbers of interacting/internalised particles over time')
fig2.Position = [100,100,1100,700];
saveas(fig2, [PARAMETERS.folder_path '/Dosage_distribs' num2str(X) 'hours'], 'eps')
saveas(fig2, [PARAMETERS.folder_path '/Dosage_distribs' num2str(X) 'hours'], 'png')
figure(fig2)

sgtitle(fig3, 'External (stain 2) fluorescence of cells against associated (stain 1) fluorescence')
fig3.Position = [100,100,1300,700];
saveas(fig3, [PARAMETERS.folder_path '\Fluor_dot_plots' num2str(X) 'hours'], 'eps')
saveas(fig3, [PARAMETERS.folder_path '\Fluor_dot_plots' num2str(X) 'hours'], 'png')
figure(3)

sgtitle(fig9, 'Frequency of cells over time ... split into number of cell divisions')
fig9.Position = [20,20,1300,700];
saveas(fig9, [PARAMETERS.folder_path '\Hist_split_divs'], 'eps')
saveas(fig9, [PARAMETERS.folder_path '\Hist_split_divs'], 'png')

figure(fig2)
figure(fig3)
figure(fig9)
end
