% CELL_SIMULATION_RUNS_DRIVER is a script that drives the cells_simulation 
% simulation with the specified parameters for a given number of iterations 
% and produces summary statistics from the numerous runs.
%
%   This is the work of Celia Dowling 23/02/22
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
%       culture_dim (cell diameters) Lattice dimensions are culture_dim by culture_dim.
%       culture_media_height (mm)   The height of the cell culture media
%                                   above the cell monolayer
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
%       EWTs_internalise.input_type Contains information about how to 
%                                   calculate the L expected waiting times 
%                                   for 1 particle to transition out of 
%                                   each stage in the cell-particle 
%                                   interaction model. input_type will be 
%                                   either "fraction" or "EWT"
%       EWTs_internalise.values     if input_type == "fraction":
%                                       [fraction associated, fraction
%                                       internalised, number of hours of 
%                                       incubation after which these
%                                       fractions are desired to be
%                                       observed] ... at confluence without
%                                       carrying capacity
%                                   if input_type == "EWT":
%                                       e.g. [free or stage 0 (hit), stage 
%                                       1 (bound), ..., stage L-1] (hours) 
%                                       ... List of L expected waiting 
%                                       times. Can be 1 EWT per  stage 
%                                       (e.g. mean time spent free and 
%                                       unbound, mean time spent in stage 1  
%                                       (bound), mean time spent in stage 
%                                       2, ..., mean time spent in stage 
%                                       L-1 before transitioning into stage
%                                       L i.e. being internalised) OR can 
%                                       be 1 EWT per stage per cell phase; 
%                                       columns indicating interaction 
%                                       stage and rows indicating cell 
%                                       phase. Elements must be greater 
%                                       than or equal to tstep_duration. 
%                                   For a petri dish saturated with cells:
%                                   Given a certain percentage association 
%                                   over 24 hours (a), choose values(1) 
%                                   such that the CDF of the exponential
%                                   distribution with that rate satisfies
%                                   F(24)=a. Given a certain percentage
%                                   internalisation over 24 hours (i),
%                                   choose all other values rates such that
%                                   the CDF of the hypo-exponential 
%                                   distribution with these rates satisfy 
%                                   F(24)=i.
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
%       temp (degrees celsius)      Temperature of the cell culture
%       viscos (kg / (m*s))         (Dynamic) viscosity of the cell culture
%                                   medium
%       prtcl_rad (nanometers)      The average radius of the nanoparticle 
%                                   type.
%       vid_speed (frames per sec)  Speed of movie frame playback.
%       visual                      The form of visualisation desired which 
%                                   can be 0,1 or 2:
%                                       0 = off
%                                       1 = slower, more comprehensive vid
%                                           + gif
%                                       2 = faster, less comprehensive vid
%                                           + gif
%
%   The following must also be specified:
%       num_runs                    The number of desired runs.
%       tol                         The tolerance (about 0) in the 
%                                   derivative of the (smoothed) data of
%                                   average internalised/ interacting
%                                   particles per cell over time required
%                                   for identifying when when zero/ 
%                                   constant derivative is reached
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
%                                   multiple runs. Overall standard 
%                                   deviation clouds and run-to-run 
%                                   deviation in the mean is visualised.
%       PCC_and_confluence          Plots the pair correlation coefficient
%                                   for all cells across domains in all
%                                   runs against the average trajectory in
%                                   cell population growth.
%   Additionally is printed:
%       estimates for lambda values and EWTs for all stages of the
%       cell-particle interaction model
%   the internalisation process are


% Prepare
clear;
close all;

% Seed for efficiency purposes (optional)
%rng(22)

% Choose number of runs
num_runs = 100;

% Choose a tolerance for gradient matching
tol = 1E-1;
tol_l2 = 5E-4;

PARAMETERS = struct( ...
    'simulation_duration',24, ... (hours)
    'tstep_duration',1/6, ... (hours)
    'initial_num_cells', 100, ... 
    'prtcls_per_cell', 5000, ...
    'cell_diam', 25, ... (micrometers)
    'culture_dim', 10, ... (cell diameters)
    'culture_media_height', 5,... (millimeters)
    'EWT_move', 1/6, ... (hours) 
    'EWTs_proliferate', [4,4,4], ... [4,4,4], ... [phase 1, ..., phase K](hours) 
    'EWTs_internalise', struct('input_type', "fraction", ... "fraction" or "EWT"
    'values', [0.01,0.006,24]), ... see notes on EWTs_internalise [26.19256, 5.36034], ...[34.62471997,12.52770188], ... 
    'max_prtcls', [inf,50], ... [stage 1, ..., stage L]
    'prob_inherit', 0.7, ...     
    'temp', 36, ... (degrees celsius)    
    'viscos', 1.0005E-3,... (kiloggrams / (meter*second))      
    'prtcl_rad', 25, ... (nanometers)    
    'vid_speed', 2, ... (frames per sec)     
    'visual', 0 ...
    );

date_time = num2str(fix(clock));
date_time = date_time(find(~isspace(date_time)));

folder_name = ['Number_Runs' num2str(num_runs) ...
    '_N0' num2str(PARAMETERS.initial_num_cells) ...
    '_dim' num2str(PARAMETERS.culture_dim) ...
    '_maxPrtcls' num2str(PARAMETERS.max_prtcls) ...
    '_pPerC' num2str(PARAMETERS.prtcls_per_cell) ...
    '_EWTintern' num2str(PARAMETERS.EWTs_internalise.values) ...
    '_' date_time];
folder_name = folder_name(find(~isspace(folder_name)));

PARAMETERS.folder_path = [pwd '/' date '/' folder_name];

if ~exist(PARAMETERS.folder_path, 'dir')
    mkdir(PARAMETERS.folder_path)
end
                
% Run multiple simulations and collect data
total_tsteps = floor(PARAMETERS.simulation_duration/PARAMETERS.tstep_duration);
total.cell_c_o_p = zeros(num_runs*PARAMETERS.culture_dim^2,total_tsteps+1,3);
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
    runs.average_cell_c_o_p(run,:,:) = ... % Average number of particles per cell in a class...
        [sum(EVOLUTION_INFO.cell_c_o_p(:,:,1),1)./EVOLUTION_INFO.cell_population; ... % free or hit
        sum(EVOLUTION_INFO.cell_c_o_p(:,:,2),1)./EVOLUTION_INFO.cell_population; ... % Interacting
        sum(EVOLUTION_INFO.cell_c_o_p(:,:,3),1)./EVOLUTION_INFO.cell_population]'; % Internalised
    runs.cell_pair_cor_coef(run,:) = EVOLUTION_INFO.cell_pair_cor_coef;
end
fprintf('The average final population of cells in a run:')
disp(total.cell_population(end)/num_runs)

% Find the average pair correlation coefficient
total.cell_pair_cor_coef = [mean(runs.cell_pair_cor_coef(:,:), 'all'),...
        var(runs.cell_pair_cor_coef(:,:), 0, 'all')];
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
savefig(fig2, [PARAMETERS.folder_path '/Frequency_histograms'], 'compact')
saveas(fig2, [PARAMETERS.folder_path '/Frequency_histograms'], 'png')

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
    % Variance in interacting particles
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
exp_cdf = @(lambda,tdata) 1-exp(-lambda.*tdata);

%plot(binrng(2:end),exp_cdf(est_lambda1.using_diffs,binrng(2:end)),'g:')
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
savefig(fig7, [PARAMETERS.folder_path '/Associated_vs_internalised'], 'compact')
saveas(fig7, [PARAMETERS.folder_path '/Associated_vs_internalised'], 'png')


% Plot the mean Pair Correlation Coefficient over time (the mean PCC at
% each timestep across all of the run PCCs and the variance in the mean)
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
savefig(fig8, [PARAMETERS.folder_path '/PCC_and_confluence'], 'compact')
saveas(fig8, [PARAMETERS.folder_path '/PCC_and_confluence'], 'png')

% HEURISTIC ESTIMATES
%   of tmax_noCC, lambda1, lambda2 and CC
% Use the provided data on the mean numbers of particles per cell that are 
% free, interacting and internalised and estimate the parameters of the 
% exponential and hypoexponential distributions that the inter-event times 
% follow (of particle association and particle internalisation respectively).

fprintf("\nHeuristic estimates for lambda from free to stage 1 and from stage 1 to stage 2: \n(units are per hour)\n")

% The mean number of free particles available to a cell at the beginning of
% each timestep is the number that haven't been associated on the previous
% timestep divided amongst the lattice sites. 
freePrtcls_start_of_t = (PARAMETERS.prtcls_per_cell * PARAMETERS.initial_num_cells .* ... % mean
    ones(1,length(binrng)) - (means(3,:) .* total.cell_population) ./ num_runs) ./ ...
    (PARAMETERS.culture_dim^2);
freePrtcls_start_of_tUP = (PARAMETERS.prtcls_per_cell * PARAMETERS.initial_num_cells .* ... % 1 std.dev. above
    ones(1,length(binrng)) - (upper1(3,:) .* total.cell_population) ./ num_runs) ./ ...
    (PARAMETERS.culture_dim^2);
freePrtcls_start_of_tLO = (PARAMETERS.prtcls_per_cell * PARAMETERS.initial_num_cells .* ... % 1 std.dev. below
    ones(1,length(binrng)) - (lower1(3,:) .* total.cell_population) ./ num_runs) ./ ...
    (PARAMETERS.culture_dim^2);
% No interaction occur prior to 0 hours, so the entry for '0 hours' 
% (timestep 0 or typically index 1) is in fact the free particles available 
% at the beginning of timestep 1 (typically index 2) from 0-0.1667 hours 
% and so on. The last entry is simply the number of free particles 
% available per cell at the end of the simulation.
freePrtcls_start_of_t = freePrtcls_start_of_t(1:(end-1)); % for timestep 1, 2, ...
freePrtcls_start_of_tUP = freePrtcls_start_of_tUP(1:(end-1)); % 1 std.dev. above
freePrtcls_start_of_tLO = freePrtcls_start_of_tLO(1:(end-1)); % 1 std.dev. below

% The mean number of particles interacting with a cell at the beginning of
% each timestep is the number that are recorded in the previous timestep.
% I.e. the interacting particles available at the beginning of timestep 1
% (typically index 2) from 0-0.1667 hours is the number recorded as
% interacting at the end of timestep 0 (typically index 1) at 0 hours.
interactPrtcls_start_of_t = means(1,1:(end-1)); % for timestep 1, 2, ...
interactPrtcls_start_of_tUP = means(3,1:(end-1))-upper1(2,1:(end-1)); % 1 std.dev. above
interactPrtcls_start_of_tLO = means(3,1:(end-1))-lower1(2,1:(end-1)); % 1 std.dev. below

% ESTIMATE LAMBDA 1
%   via DIFFERENCES method
% Estimate lambda_1, the parameter for the exponential distribution that
% describes the inter-association times of particles per cell, by
% calculating the mean fraction of association events that occur (the
% number that do occur in a timestep over the number that could occur -
% i.e., the number of free particles) over the duration of time passed.
est_lambda1.using_diffs = (means(3,2:end)- means(3,1:(end-1)))./... new in tstep 1, 2, ...
    (binrng(2).*freePrtcls_start_of_t); % free at the beginning of tstep 1, 2, ...
est_lambda1.using_diffs_mean = mean(est_lambda1.using_diffs);
% 1 std.dev. above
est_lambda1UP.using_diffs = (upper1(3,2:end)- upper1(3,1:(end-1)))./... 
    (binrng(2).*freePrtcls_start_of_tUP); 
est_lambda1UP.using_diffs_mean = mean(est_lambda1UP.using_diffs);
% 1 std.dev. below
est_lambda1LO.using_diffs = (lower1(3,2:end)- lower1(3,1:(end-1)))./... 
    (binrng(2).*freePrtcls_start_of_tLO); 
est_lambda1LO.using_diffs_mean = mean(est_lambda1LO.using_diffs);

% ESTIMATE LAMBDA 1
%   via MLE method
smth_assoc = smooth(means(3,:))'; % number associated at end of tstep 0, 1, 2, ...
est_lambda1.MLE = smth_assoc(2:end)./ ... % num assoc at end of tstep 1, 2, ...
    (binrng(2:end).*PARAMETERS.prtcls_per_cell); % time at end of tstep 1, 2, ...
est_lambda1.MLE_mean = mean(est_lambda1.MLE);
% 1 std.dev. above
smth_assocUP = smooth(upper1(3,:))'; 
est_lambda1UP.MLE = smth_assocUP(2:end)./ ... 
    (binrng(2:end).*PARAMETERS.prtcls_per_cell); 
est_lambda1UP.MLE_mean = mean(est_lambda1UP.MLE);
% 1 std.dev. below
smth_assocLO = smooth(lower1(3,:))'; 
est_lambda1LO.MLE = smth_assocLO(2:end)./ ...
    (binrng(2:end).*PARAMETERS.prtcls_per_cell); 
est_lambda1LO.MLE_mean = mean(est_lambda1LO.MLE);

% CDF of hypoexponential distribution
CDF.hypoexp = @(l1, l2, t) 1 - 1./(l2-l1) .* (l2 .* exp(-l1 .* t) - l1 .* exp(-l2 .* t));
CDF.hypoexp_MLEl1 = @(l2,t) CDF.hypoexp(est_lambda1.MLE_mean,l2,t);
CDF.hypoexp_MLEl1_tstep = @(l2) CDF.hypoexp_MLEl1(l2,PARAMETERS.tstep_duration);
% CDF of exponential distribution
CDF.exp = @(l, t) 1 - exp(- l .*t);
CDF.exp_tstep = @(l) CDF.exp(l,PARAMETERS.tstep_duration);
% Two guesses to feed into fzero in case one gives lambda2=lambda1
guess = [est_lambda1.MLE_mean/3 est_lambda1.MLE_mean*3]; % mean

% ESTIMATE TIME AT WHICH CC KICKS IN (TMAX_NOCC) and LAMBDA 2
%   via DISTRIBUTION method
[tmax_noCC,est_lambda2.distrib_mean,est_lambda2.distrib] = ... 
    tmax_l2_from_hypoexpCDF(binrng,CDF.hypoexp_MLEl1,means,PARAMETERS,...
    guess,est_lambda1.MLE_mean,tol_l2);
% 1 std.dev. above - assuming that lambda1 is the mean MLE
[~,est_lambda2UP.distrib_mean,est_lambda2UP.distrib] = ... 
    tmax_l2_from_hypoexpCDF(binrng,CDF.hypoexp_MLEl1,upper1,PARAMETERS,...
    guess,est_lambda1UP.MLE_mean,tol_l2);
% 1 std.dev. above - assuming that lambda1 is the mean MLE
[~,est_lambda2LO.distrib_mean,est_lambda2LO.distrib] = ... 
    tmax_l2_from_hypoexpCDF(binrng,CDF.hypoexp_MLEl1,lower1,PARAMETERS,...
    guess,est_lambda1LO.MLE_mean,tol_l2);

% ESTIMATE LAMBDA 2
%   via MIX method
[est_lambda2.mix_mean,est_lambda2.mix] = l2_from_mixMethod(...
    binrng,CDF.hypoexp_MLEl1_tstep,CDF.exp_tstep,...
    means,freePrtcls_start_of_t,interactPrtcls_start_of_t,...
    PARAMETERS,guess,est_lambda1.MLE_mean,tmax_noCC);
% 1 std.dev. above - assuming that lambda1 is the mean MLE
[est_lambda2UP.mix_mean,est_lambda2UP.mix] = l2_from_mixMethod(...
    binrng,CDF.hypoexp_MLEl1_tstep,CDF.exp_tstep,...
    upper1,freePrtcls_start_of_t,interactPrtcls_start_of_tUP,...
    PARAMETERS,guess,est_lambda1.MLE_mean,tmax_noCC);
% 1 std.dev. below - assuming that lambda1 is the mean MLE
[est_lambda2LO.mix_mean,est_lambda2LO.mix] = l2_from_mixMethod(...
    binrng,CDF.hypoexp_MLEl1_tstep,CDF.exp_tstep,...
    lower1,freePrtcls_start_of_t,interactPrtcls_start_of_tLO,...
    PARAMETERS,guess,est_lambda1.MLE_mean,tmax_noCC);

tsteps = 0:PARAMETERS.tstep_duration:tmax_noCC; 

% ESTIMATE LAMBDA 2
%   via DIFFERENCES method
[est_lambda2.using_diffs_mean,est_lambda2.using_diffs] = ...
    l2_from_diffsMethod(tsteps,means,tol,freePrtcls_start_of_t, ...
    interactPrtcls_start_of_t,PARAMETERS,est_lambda1.using_diffs_mean);
% 1 std.dev. above - assuming that lambda1 is the mean from differences
[est_lambda2UP.using_diffs_mean,est_lambda2UP.using_diffs] = ...
    l2_from_diffsMethod(tsteps,means,tol,freePrtcls_start_of_t, ...
    interactPrtcls_start_of_tUP,PARAMETERS,est_lambda1.using_diffs_mean);
% 1 std.dev. below - assuming that lambda1 is the mean from differences
[est_lambda2LO.using_diffs_mean,est_lambda2LO.using_diffs] = ...
    l2_from_diffsMethod(tsteps,means,tol,freePrtcls_start_of_t, ...
    interactPrtcls_start_of_tLO,PARAMETERS,est_lambda1.using_diffs_mean);

% CALCULATE ACTUAL RATES
[l1,l2] = input_EWT_from_fraction(...
    PARAMETERS.EWTs_internalise.values(1),...
    PARAMETERS.EWTs_internalise.values(2),...
    PARAMETERS.EWTs_internalise.values(3));

% PRINT RATES
fprintf("\nLAMBDA 1: \nACTUAL %5.4e",l1)
fprintf('\nMETHOD \t\tMEAN \t\tLOWER ESTIMATE \tUPPER ESTIMATE')
fprintf('\nDifferences \t%5.4e \t%5.4e \t%5.4e',...
    est_lambda1.using_diffs_mean, est_lambda1LO.using_diffs_mean, ...
    est_lambda1UP.using_diffs_mean);
fprintf('\nMLE \t\t%5.4e \t%5.4e \t%5.4e \n',...
    est_lambda1.MLE_mean, est_lambda1LO.MLE_mean, ...
    est_lambda1UP.MLE_mean);

fprintf("\nLAMBDA 2: \nACTUAL %5.4e",l2)
fprintf('\nMETHOD \t\tMEAN \t\tLOWER ESTIMATE \tUPPER ESTIMATE')
fprintf('\nDifferences \t%5.4e \t%5.4e \t%5.4e',...
    est_lambda2.using_diffs_mean, est_lambda2LO.using_diffs_mean, ...
    est_lambda2UP.using_diffs_mean);
fprintf('\nDistribution \t%5.4e \t%5.4e \t%5.4e',...
    est_lambda2.distrib_mean, est_lambda2LO.distrib_mean, ...
    est_lambda2UP.distrib_mean);
fprintf('\nMix \t\t%5.4e \t%5.4e \t%5.4e\n',...
    est_lambda2.mix_mean, est_lambda2LO.mix_mean, ...
    est_lambda2UP.mix_mean);


% ESTIMATE CARRYING CAPACITY (CC) **if there is one**
%   via DIFFERENCES method
%   via DYNAMIC RATE calculated from DIFFERENCES METHOD
%   via DYNAMIC RATE calculated from MIX METHOD
if any(PARAMETERS.max_prtcls ~= inf) 
    % The number of internalised particles per cell at the start of a 
    % timestep is the number of internalised particles per cell at the end
    % of the previous timestep.
    internalPrtcls_start_of_t = means(2,1:(end-1)); % for timestep 1, 2, ...
    
    % Prepare input for CC_GIVEN_FRACTION
    using_diffs = gradient(means(2,2:end),PARAMETERS.tstep_duration) ./ ... 
        (interactPrtcls_start_of_t .* est_lambda2.using_diffs_mean);
    dynamic_diffs= est_lambda2.using_diffs./est_lambda2.using_diffs_mean;
    dynamic_mix= est_lambda2.mix./est_lambda2.mix_mean;
    frction = [using_diffs; dynamic_diffs; dynamic_mix];
    methods = {'using_diffs','dynamic_diffs','dynamic_mix'};

    est_CC = CC_given_fraction(methods,internalPrtcls_start_of_t,frction, PARAMETERS);
    
    % PRINT CARRYING CAPACITIES
    fprintf("\n\nHeuristic estimates for carrying capacity (CC): \n(units are particles)\n")
    
    fprintf("\nCARRYING CAPACITY: \nACTUAL %5.3f",PARAMETERS.max_prtcls(end))
    fprintf('\nMETHOD \t\t\t\tMEAN \tLOWER ESTIMATE \tUPPER ESTIMATE')
    fprintf('\nDifferences \t\t\t%5.3f',est_CC.using_diffs_mean);
    fprintf('\nDynamic using differences \t%5.3f',est_CC.dynamic_diffs_mean);
    fprintf('\nDynamic using mix \t\t%5.3f\n',est_CC.dynamic_mix_mean);
end

% Print figures
figure(fig2)
figure(fig7)
figure(fig8)

% Save the workspace
save([PARAMETERS.folder_path '/variables_' num2str(PARAMETERS.prtcls_per_cell) ...
    'pPerC_CC' num2str(PARAMETERS.max_prtcls(end)) ...
    '_' num2str(PARAMETERS.EWTs_internalise.values(1)) ...
    '_' num2str(PARAMETERS.EWTs_internalise.values(2)) ...
    '_' date_time '.mat']);

