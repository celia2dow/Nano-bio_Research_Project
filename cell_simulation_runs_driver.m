% CELL_SIMULATION_RUNS_DRIVER is a script that drives the cells_simulation 
% simulation with the specified parameters for a given number of iterations 
% and produces summary statistics from the numerous runs.
%
%   This is the work of Celia Dowling 1/12/21
%
%   The input argument for cells_simulation.m is a structure PARAMETERS 
%   which has the following fields that need to be defined by the user:
%
%       simulation_duration (hours) Time length of the simulation.    
%       tstep_duration (hours)      Recommended to be equal to the shortest
%                                   expected waiting time (EWT), usually 
%                                   EWT_move.
%       initial_num_cells           Initial number of cells in the lattice.
%       prtcls_per_cell             Number of particles added to the cell
%                                   monolayer per cell
%       cell_diam (micrometers)     Average diameter of particlular cell-type. 
%       dim (cell diameters)        Lattice dimensions are dim by dim.
%       EWT_move (hours)            Expected waiting time for one cell to 
%                                   move 1 cell diameter. Must be greater than
%                                   or equal to tstep_duration.
%       EWTs_proliferate (hours)    List of K expected waiting times for 1
%                                   cell to transition out of each phase in
%                                   the cell proliferation cycle (e.g. mean 
%                                   time spent in phase 1, mean time spent 
%                                   in phase 2,..., mean time spent in phase
%                                   K before proliferating). Elements must 
%                                   be greater than or equal to
%                                   tstep_duration. They should add up to
%                                   give the expected waiting time for
%                                   proliferation.
%       EWFs_internalise (hours)    List of L expected waiting times for 1
%                                   particle to transition out of each
%                                   stage in the cell-particle interaction 
%                                   model. Can be 1 EWT per stage (e.g. 
%                                   mean time spent free and unbound, mean 
%                                   time spent in stage 1 (bound), mean time
%                                   spent in stage 2, ..., mean time spent 
%                                   in stage L-1 before transitioning into 
%                                   stage L i.e. being internalised) or can 
%                                   be 1 EWT per stage per cell phase; columns 
%                                   indicating interaction stage and rows 
%                                   indicating cell phase. Elements must be 
%                                   greater than or equal to tstep_duration. 
%                                   For a petri dish saturated with cells:
%                                   Given a certain percentage association 
%                                   over 24 hours (a), choose EWFs_internalise(1) 
%                                   such that the CDF of the exponential
%                                   distribution with that rate satisfies
%                                   F(24)=a. Given a certain percentage
%                                   internalisation over 24 hours (i),
%                                   choose all other EWFs_internalise rates
%                                   such that the CDF of the hypo- 
%                                   exponential distribution with these
%                                   rates satisfy F(24)=i.
%       max_prtcls                  The particle capacities of each of the
%                                   L cell-particle interaction model
%                                   stages (i.e., capacity of stage 1, ...,
%                                   capacity of stage L). If a stage has 
%                                   unlimited capacity, set it to "inf". If
%                                   all stages have been set to "inf", 
%                                   carrying capacity is effectively 
%                                   switched off.
%       prob_inherit                The probability of a daughter cell born 
%                                   into the site of the parent cell 
%                                   inheriting 1 particle from the parent
%                                   cell.
%       rate_prtcl_diffusivity (per hour) The rate at which particles interact  
%                                   with a given cell (fraction of free particles
%                                   per hour that hit), reflecting the diffusivity 
%                                   of the particle type.
%       vid_speed (frames per sec)  Speed of movie frame playback.
%       visual                      The form of visualisation desired which 
%                                   can be 0,1 or 2:
%                                       0 = off
%                                       1 = slower, more comprehensive vid
%                                           + gif
%                                       2 = faster, less comprehensive vid
%                                           + gif
%
%   The number of runs desired to be run, num_runs, must also be specified.
%
%   Additionally and optionally, the random number generator characterising
%   the simualtion can also be seeded.
%
%   The following output is generated by the driver function and saved in
%   the specified folder path:
%
%       Frequency_histograms        On hourly intervals, the frequencies  
%                                   of cells with certain numbers of 
%                                   internalised (red) and interacting 
%                                   (blue) particles are plotted. These
%                                   frequencies are the averages over
%                                   multiple runs.
%       Associated_vs_internalised  Plots the average number of
%                                   internalised, interacting and
%                                   associated (both) particles per cell
%                                   over time, as averaged out over the 
%                                   multiple runs with. Overall standard 
%                                   deviation clouds and run-to-run 
%                                   deviation in the mean is visualised.


% Prepare
clear;
close all;

% Seed for efficiency purposes (optional)
%rng(22)

% Choose number of runs
num_runs = 1;


PARAMETERS = struct( ...
    'simulation_duration',24, ... (hours)
    'tstep_duration',1/6, ... (hours)
    'initial_num_cells', 500, ... 
    'prtcls_per_cell', 100, ...
    'cell_diam', 25, ... (micrometers)
    'dim', 25, ... (cell diameters)
    'EWT_move', 1/6, ... (hours) 
    'EWTs_proliferate', [4,4,4], ... [phase 1, ..., phase K](hours) 
    'EWTs_internalise',[34.62471997,12.52770188], ... [free or stage 0 (hit), stage 1 (bound), ..., stage L-1](hours)
    'max_prtcls', [inf,inf], ... [stage 1, ..., stage L]
    'prob_inherit', 0.7, ...     
    'rate_prtcl_diffusivity',0.4, ...  (per hour)
    'vid_speed', 2, ... (frames per sec)     
    'visual', 0 ...
    );

folder_name = ['Number_Runs' num2str(num_runs) ...
    'Simulation_Duration' num2str(PARAMETERS.simulation_duration) ...
    '_tstepDuration' num2str(PARAMETERS.tstep_duration) ...
    '_N0' num2str(PARAMETERS.initial_num_cells) ...
    '_dim' num2str(PARAMETERS.dim) ...
    '_cellDiam' num2str(PARAMETERS.cell_diam) ...
    '_maxPrtcls' num2str(PARAMETERS.max_prtcls) ...
    '_EWTmove' num2str(PARAMETERS.EWT_move) ...
    '_Pinherit' num2str(PARAMETERS.prob_inherit) ...
    '_EWTprolif' num2str(PARAMETERS.EWTs_proliferate) ...
    '_rateDiffus' num2str(PARAMETERS.rate_prtcl_diffusivity) ...
    '_EWTintern' num2str(PARAMETERS.EWTs_internalise) ...
    ];

PARAMETERS.folder_path = [pwd '\' date '\' folder_name];

if ~exist(PARAMETERS.folder_path, 'dir')
    mkdir(PARAMETERS.folder_path)
end
                
% Run multiple simulations and collect data
total_tsteps = floor(PARAMETERS.simulation_duration/PARAMETERS.tstep_duration);
total.cell_c_o_p = zeros(num_runs*PARAMETERS.dim^2,total_tsteps+1,3);
total.cell_population = zeros(1,total_tsteps+1);
runs.cell_pair_cor_coef = zeros(num_runs,total_tsteps+1);
runs.average_cell_c_o_p = zeros(num_runs,total_tsteps+1,3);
for run = 1:num_runs
    EVOLUTION_INFO = cells_simulation(PARAMETERS);
    total.cell_population = total.cell_population + EVOLUTION_INFO.cell_population;
    N_final = EVOLUTION_INFO.cell_population(end); % final number of cells
    cell_start = (run - 1) * N_final + 1;
    cell_end = run * N_final;
    total.cell_c_o_p(cell_start:cell_end,:,:) = EVOLUTION_INFO.cell_c_o_p;
    runs.average_cell_c_o_p(run,:,:) = ...
        [sum(EVOLUTION_INFO.cell_c_o_p(:,:,1),1)./EVOLUTION_INFO.cell_population; ...
        sum(EVOLUTION_INFO.cell_c_o_p(:,:,2),1)./EVOLUTION_INFO.cell_population; ...
        sum(EVOLUTION_INFO.cell_c_o_p(:,:,3),1)./EVOLUTION_INFO.cell_population]';
    runs.cell_pair_cor_coef(run,:) = EVOLUTION_INFO.cell_pair_cor_coef;
end
fprintf('The average final population of cells in a run:')
disp(total.cell_population(end)/num_runs)

% Find the average pair correlation coefficient
total.cell_pair_cor_coef = [mean(runs.cell_pair_cor_coef(:,1), 'all'),...
        var(runs.cell_pair_cor_coef(:,1), 0, 'all')];
fprintf("The mean and variance of the pair correlation coefficient in each run:")
format SHORTE
disp(total.cell_pair_cor_coef)

% Find limits for the x and y-axis: the largest values to be plotted, given
% that plots are being plotted every hour
hour_indices = 0:floor(1/PARAMETERS.tstep_duration):total_tsteps;
hour_indices = hour_indices + 1;
hourly_total.cell_c_o_p = total.cell_c_o_p(:,hour_indices,:);
interact_max = max(hourly_total.cell_c_o_p(:,:,2),[],'all');
internal_max = max(hourly_total.cell_c_o_p(:,:,3),[],'all');
x_max = max(interact_max,internal_max);
local_max = zeros(1,length(hour_indices));
for index = 2:length(hour_indices)
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

for time_plot = 1:length(hour_indices)
    tstep = (time_plot-1)/PARAMETERS.tstep_duration; % equivalent timestep index
    N_tstep = total.cell_population(tstep+1); % total number of cells across runs
    
    % FREQUENCY HISTOGRAMS
    set(0,'CurrentFigure',fig2)
    subplot(5,floor(PARAMETERS.simulation_duration/5)+1,time_plot);
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
    title(['At ' num2str(time_plot-1) ' hour/s']);
    xlabel('Number of particles');
    ylabel('Cell frequency');
    if time_plot==1
        legend('Interacting','Internalised');
    end
end
sgtitle(fig2, 'Frequency of cells with certain numbers of interacting/internalised particles over time')
fig2.Position = [100,100,1300,700];
savefig(fig2, [PARAMETERS.folder_path '\Frequency_histograms'], 'compact')
saveas(fig2, [PARAMETERS.folder_path '\Frequency_histograms'], 'png')

% Plot associated Vs interacting+internalised cells over time
means = zeros(3,total_tsteps+1);
% Average number interacting particles per cell over all runs in a timestep
means(1,:) = sum(total.cell_c_o_p(:,:,2),1)./total.cell_population; 
% Average number internalised particles per cell over all runs in a timestep
means(2,:) = sum(total.cell_c_o_p(:,:,3),1)./total.cell_population;     
% Average associated particles per cell over time
means(3,:) = means(1,:) + means(2,:);   
% Create bin ranges
binrng = 0:PARAMETERS.tstep_duration:PARAMETERS.simulation_duration; 
% Calculate variances for each timestep
variances = zeros(2,3,total_tsteps+1);
total.cell_associated = sum(total.cell_c_o_p(:,:,2:3),3);
for tstep = 0:total_tsteps
    N_tstep = total.cell_population(tstep+1); 
    
    % BETWEEN ALL CELLS IN ALL RUNS
    % Variances in interacting particles
    variances(1,1,tstep+1) = var(total.cell_c_o_p(1:N_tstep,tstep+1,2)); 
    % Variances in internalised particles
    variances(1,2,tstep+1) = var(total.cell_c_o_p(1:N_tstep,tstep+1,3)); 
    % Variances in associated particles
    variances(1,3,tstep+1) = var(total.cell_associated(1:N_tstep,tstep+1));
    
    % BETWEEN RUN MEANS
    variances(2,1,tstep+1) = var(runs.average_cell_c_o_p(:,tstep+1,2)); 
    % Variances in internalised particles
    variances(2,2,tstep+1) = var(runs.average_cell_c_o_p(:,tstep+1,3)); 
    % Variances in associated particles
    variances(2,3,tstep+1) = var(sum(runs.average_cell_c_o_p(:,tstep+1,2:3),3));
end

% Create arrays for plotting standard error clouds and variance between
% runs
st_dev1 = squeeze(sqrt(variances(1,:,:))); % standard deviations between all
upper1 = means + st_dev1;
lower1 = means - st_dev1;
st_dev2 = squeeze(sqrt(variances(2,:,:))); % standard deviations between runs
upper2 = means + st_dev2;
lower2 = means - st_dev2;

fig7=figure(7);
set(fig7, 'Visible', 'off');
subplot(1,2,1)
bar(binrng, means(3,:), 'r')
hold on
bar(binrng, means(1,:), 'b')
hold off
legend('Internalised','Interacting','Location','Best')
title('Stacked histogram')
xlabel('time $t$ hours', 'Interpreter', 'latex');
ylabel('Average number of particles per cell');

subplot(1,2,2)
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
hold off
ylim([0 max(upper1,[],'all')])
legend('Internalised','Interacting','Associated (either)',...
    'Std. dev. btwn all cells & runs', '" "', '" "', ...
    'Std. dev. btwn run means', '','" "','', '" "', 'Location','Best')
title('Separate trends')
xlabel('time $t$ hours', 'Interpreter', 'latex');

subplot(1,2,1)
ylim([0 max(upper1,[],'all')])
sgtitle('Average number of particles associated per cell')
fig7.Position = [100,100,1300,700];
savefig(fig7, [PARAMETERS.folder_path '\Associated_vs_internalised'], 'compact')
saveas(fig7, [PARAMETERS.folder_path '\Associated_vs_internalised'], 'png')

% Print figures
figure(fig2)
figure(fig7)

% Save the workspace
save([PARAMETERS.folder_path '\variables.mat']);
