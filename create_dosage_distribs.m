function create_dosage_distribs(X, total_tsteps, PARAMETERS, total, FLUORESC, max_divs, ith)
% CREATE_DOSAGE_DISTRIBS creates dosage distributions, fluorescence plots, 
% and dosages split by number of cell divisions after every X hours.

% Find the hours on which plots will be made
increments = floor(X/(PARAMETERS.tstep_duration * ith));
Xhour_indices = 0:increments:total_tsteps;
Xhour_indices = Xhour_indices(2:end) + 1;
dim1 = floor(sqrt(length(Xhour_indices)));
dim2 = ceil(length(Xhour_indices)/dim1);

% Find limits for the x and y-axis: the largest values to be plotted, given
% that plots are being plotted every hour
Xhourly_total.cell_c_o_p = total.cell_c_o_p(:,Xhour_indices,:);
interact_max = max(Xhourly_total.cell_c_o_p(:,:,2),[],'all');
internal_max = max(Xhourly_total.cell_c_o_p(:,:,3),[],'all');
associat_max = max(sum(Xhourly_total.cell_c_o_p(:,:,2:3),3),[],'all');
interact_min = min(Xhourly_total.cell_c_o_p(:,:,2),[],'all');
internal_min = min(Xhourly_total.cell_c_o_p(:,:,3),[],'all');
associat_min = min(sum(Xhourly_total.cell_c_o_p(:,:,2:3),3),[],'all');
x_max = max(interact_max,internal_max);
local_max = zeros(1,length(Xhour_indices));

% Find the largest number of cell divisions to be included
max_cell_divs = min(max(total.cell_lineage(:,3:end),[],'all'),max_divs);

for index = Xhour_indices(1:end-1)
    [cells_interact,~] = histcounts(Xhourly_total.cell_c_o_p(:,index,2));
    [cells_internal,~] = histcounts(Xhourly_total.cell_c_o_p(:,index,3));
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
S1_min = FLUORESC.stain1_part * (FLUORESC.alpha*associat_min) + FLUORESC.stain1_background ...
    - sqrt(associat_min) .* FLUORESC.stain1_std_dev;
S1_max = FLUORESC.stain1_part * (associat_max) + FLUORESC.stain1_background ...
    + sqrt(associat_max) .* FLUORESC.stain1_std_dev;
S2_min = FLUORESC.stain2_part * (interact_min) + FLUORESC.stain2_background ...
    - sqrt(interact_min) .* FLUORESC.stain2_std_dev;
S2_max = FLUORESC.stain2_part * (interact_max) + FLUORESC.stain2_background ...
    + sqrt(interact_max) .* FLUORESC.stain2_std_dev;

fig2 = figure(2);
set(fig2, 'Visible', 'off');
fig3 = figure(3);
set(fig3, 'Visible', 'off');
fig9 = figure(9);
set(fig9, 'Visible', 'off');

for time_plot = 1:length(Xhour_indices)
    Xth_hour_index = Xhour_indices(time_plot); % equivalent timestep index
    N_tstep = total.cell_population(Xth_hour_index); % total number of cells across runs
    n_interact = total.cell_c_o_p_corrected(1:N_tstep,Xth_hour_index,2);
    n_internal = total.cell_c_o_p_corrected(1:N_tstep,Xth_hour_index,3);
    
    % FREQUENCY HISTOGRAMS
    set(0,'CurrentFigure',fig2)
    subplot(dim1,dim2,time_plot);
    if any(Xhourly_total.cell_c_o_p(:,:,2), 'all')
        % FREQUENCY OF CELLS WITH NUMS OF PARTICLES INTERACTING OVER TIME
        histogram(n_interact,'FaceColor', [0,0,1], 'FaceAlpha', 0.2, ...
            'DisplayName','Interacting');
        hold on;
    end
    % FREQUENCY OF CELLS WITH NUMS OF PARTICLES INTERNALISED OVER TIME
    histogram(n_internal, 'FaceColor', [1,0,0], 'FaceAlpha', 0.2, ...
            'DisplayName','Interalised');
    hold off;
    xlim([0,x_max]);
    %ylim([0,y_max]);
    title(['At ' num2str(time_plot*X) ' hour/s']);
    xlabel('Number of particles');
    ylabel('Cell frequency');
    pbaspect([1 1 1]);
    if time_plot==1
        legend
    end

    % FLUORESCENCE DOT PLOTS
    set(0,'CurrentFigure',fig3)
    subplot(dim1,dim2,time_plot);
    % PLOT EXTERNAL (STAIN 2) FLUORESCENCE AGAINST ASSOCIATED (STAIN 1)
    % FLUORESCENCE
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
    xlim([S1_min,S1_max]);
    if any(Xhourly_total.cell_c_o_p(:,:,2), 'all')
        ylim([S2_min,S2_max]);
    end
    title(['At ' num2str(time_plot*X) ' hour/s']);
    xlabel('Stain 1 (associated)');
    ylabel('Stain 2 (external)');
    cb = colorbar();
    cb.Label.String = 'Density estimate';

    % FREQUENCY HISTOGRAMS PER NUM OF CELL DIVISIONS EVERY X HOURS
    set(0,'CurrentFigure',fig9)
    max_cell_divs_tstep = max(total.cell_lineage(:,Xth_hour_index+2));
    min_cell_divs_tstep = min(total.cell_lineage(total.cell_lineage(:,...
        Xth_hour_index+2)>0,Xth_hour_index+2));
    % ITERATING THROUGH THE DIVISION NUMBERS PRESENT AT THIS TIMESTEP
    rows = unique(total.cell_lineage(total.cell_lineage(:, ...
        Xth_hour_index+2)>0,Xth_hour_index+2))';
    rows = rows(rows<=max_divs);
    for row = rows
        plot_num = (row-1)*length(Xhour_indices) + time_plot;
        subplot(max_cell_divs,length(Xhour_indices),plot_num);
        cells_with_divs = total.cell_lineage(total.cell_lineage(:,...
            Xth_hour_index+2)==row,2);
        if any(Xhourly_total.cell_c_o_p(:,:,2),'all')
            % INTERACTING PARTICLE DISTRIBUTIONS IN CELLS HAVING DIVIDED THIS
            % MANY TIMES
            histogram(total.cell_c_o_p(cells_with_divs,Xth_hour_index,2),...
                'FaceColor', [0,0,1], 'FaceAlpha', 0.2, ...
                'DisplayName','Interacting');
            hold on;
        end
        % INTERNALISED PARTICLE DISTRIBUTIONS IN CELLS HAVING DIVIDED THIS
        % MANY TIMES
        histogram(total.cell_c_o_p(cells_with_divs,Xth_hour_index,3),...
            'FaceColor', [1,0,0], 'FaceAlpha', 0.2, ...
            'DisplayName','Internalised');
        hold off;
        xlim([0,x_max]);
        subtitle([num2str(row) ' cell division/s'])
        xlabel('Num. of particles', 'Interpreter', 'latex');
        ylabel('Cell frequ.', 'Interpreter', 'latex');
        if row == min_cell_divs_tstep
            title(['At ' num2str(time_plot*X) ' hour/s']);
            if time_plot == 1
                legend
            end
        end
    end
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
