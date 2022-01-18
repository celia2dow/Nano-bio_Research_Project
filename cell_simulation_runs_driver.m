% CELL_SIMULATION_RUNS_DRIVER is a script that drives the cells_simulation 
% simulation with the specified parameters for a given number of iterations 
% and produces summary statistics from the numerous runs.
%
%   This is the work of Celia Dowling 19/01/22
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
%   The following must also be specified:
%       num_runs                    The number of desired runs.
%       tol_intern                  The tolerance (about 0) in the second
%                                   derivative of the smoothed data of
%                                   avergae internalised particles over 
%                                   time required for identifying when the
%                                   carrying capacity interferes with the
%                                   data.
%       tol_interact                The tolerance in the first and second
%                                   derivatives of the smoothed data of
%                                   average interacting particles over time
%                                   required for identifying a turning
%                                   point or period of approximately linear
%                                   behaviour.
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
tol_intern = 5E-2;
tol_interact = 5E-2;

PARAMETERS = struct( ...
    'simulation_duration',24, ... (hours)
    'tstep_duration',1/6, ... (hours)
    'initial_num_cells', 100, ... 
    'prtcls_per_cell', 1000, ...
    'cell_diam', 25, ... (micrometers)
    'dim', 10, ... (cell diameters)
    'EWT_move', 1/6, ... (hours) 
    'EWTs_proliferate', [4,4,4], ... [4,4,4], ... [phase 1, ..., phase K](hours) 
    'EWTs_internalise', struct('input_type', "fraction", ... "fraction" or "EWT"
    'values', [0.05,0.03,24]), ... see notes on EWTs_internalise [26.19256, 5.36034], ...[34.62471997,12.52770188], ... 
    'max_prtcls', [inf,50], ... [stage 1, ..., stage L]
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
    '_EWTintern' num2str(PARAMETERS.EWTs_internalise.values) ...
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

%plot(binrng(2:end),exp_cdf(est_lambda1,binrng(2:end)),'g:')
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
plot(binrng, total.cell_population./(PARAMETERS.dim^2 * num_runs).*100)
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
savefig(fig8, [PARAMETERS.folder_path '\PCC_and_confluence'], 'compact')
saveas(fig8, [PARAMETERS.folder_path '\PCC_and_confluence'], 'png')

% Use the provided data on the mean numbers of particles per cell that are 
% free, interacting and internalised and estimate the parameters of the 
% exponential and hypoexponential distributions that the inter-event times 
% follow (of particle association and particle internalisation respectively).

% The mean number of free particles available to a cell at the beginning of
% each timestep is the number that haven't been associated on the previous
% timestep divided amongst the lattice sites. 
freePrtcls_start_of_t = (PARAMETERS.prtcls_per_cell * PARAMETERS.initial_num_cells .* ...
    ones(1,length(binrng)) - (means(3,:) .* total.cell_population) ./ num_runs) ./ ...
    (PARAMETERS.dim^2);
% No interaction occur on timestep 0, so the entry for 0 is in fact the
% free particles available at the beginning of timestep 1 and so on.
freePrtcls_start_of_t = freePrtcls_start_of_t(1:(end-1));

% The mean number of particles interacting with a cell at the beginning of
% each timestep is the number that are recorded in the previous timestep.
interactPrtcls_start_of_t = means(1,1:(end-1));

% The MLE for lambda_1, the parameter for the exponential distribution that
% describes the inter-association times of particles per cell, is
% calculated as the mean fraction of association events that occur (the
% number that do occur in a timestep over the number that could occur -
% i.e., the number of free particles) over the duration of time passed.
est_lambda1 = (means(3,2:end)- means(3,1:(end-1)))./...
    (binrng(2).*freePrtcls_start_of_t);%PARAMETERS.prtcls_per_cell);

% Find when the carrying capacity kicks in, if one is implemented at all
if any(PARAMETERS.max_prtcls ~= inf) 
    fprintf('Carrying capacity is implemented \n')
    sld_smooth_internalised = smooth(means(2,:)).*...
        freePrtcls_start_of_t./PARAMETERS.prtcls_per_cell;
    deriv1_internalised = smooth(gradient(sld_smooth_internalised,PARAMETERS.tstep_duration));
    deriv2_internalised = gradient(deriv1_internalised,PARAMETERS.tstep_duration);
    zero_deriv2_internalise_tsteps = binrng(abs(deriv2_internalised)<=tol_intern);
    tmax_noCC = max(zero_deriv2_internalise_tsteps);
    % Limit the timesteps investigated to the tsteps up until the end of
    % the internalised distributions linear behaviour
    tsteps = 0:PARAMETERS.tstep_duration:tmax_noCC; 
else
    fprintf('No carrying capacity is implemented \n')
    tsteps = binrng;
end

% Find when the steady state in interacting particles is reached, if it
% exists
deriv1_interacting = gradient(means(1,1:length(tsteps)),PARAMETERS.tstep_duration);
zero_deriv1_interacting_tsteps = tsteps(abs(deriv1_interacting(:))<=tol_interact);
zero_deriv1_interacting_tsteps = zero_deriv1_interacting_tsteps...
    (zero_deriv1_interacting_tsteps > 0);
% Find the timesteps when a constant non-zero gradient in interacting
% particles is reached, if it exists
deriv2_interacting = gradient(deriv1_interacting,PARAMETERS.tstep_duration);
zero_deriv2_interacting_tsteps = tsteps(abs(deriv2_interacting(:))<=tol_interact);
zero_deriv2_interacting_tsteps = zero_deriv2_interacting_tsteps...
    (zero_deriv2_interacting_tsteps > 0);

% If there is a turning point and thus a steady state
if any(zero_deriv1_interacting_tsteps) && deriv1_interacting(end)<0
    % lambda_2 can be estimated by equating lambda_1*(# free particles) and
    % lambda_2*(# interacting particles) when a steady state is reached in
    % the number of iteracting particles per cell
    indices = zeros(1,length(zero_deriv1_interacting_tsteps));
    for i = 1:length(zero_deriv1_interacting_tsteps)
        indices(i) = find(tsteps == zero_deriv1_interacting_tsteps(i));
    end
    % All other arrays span the timesteps 1/12 hours to 24 hours so the
    % indices should be shifted left to account for the first index being
    % 1/12 and not 0.
    indices = indices-1;
    
    % Estimate lambda_2 using 
    %       lambda_1*(# free particles at start of t) = 
    %       lambda_2* (# interacting particles at start of t) 
    est_lambda2 = mean(est_lambda1) .* freePrtcls_start_of_t(indices) ./ ...
        interactPrtcls_start_of_t(indices);

    fprintf("Estimates for lambda and EWT from free to stage 1 and from stage 1 to stage 2: \n")
    disp([mean(est_lambda1), 1/mean(est_lambda1)]);
    disp([mean(est_lambda2), 1/mean(est_lambda2)]);
% If there is an interval of constant gradient
elseif zero_deriv2_interacting_tsteps 
    % lambda_2 can be estimated by equating the gradient of the curve of 
    % average interacting particles over time with the difference in
    % gradients of the curve for average associated particles over time *(# free particles) and
    % lambda_2*(# interacting particles) when a steady state is reached in
    % the number of iteracting particles per cell
    indices = zeros(1,length(zero_deriv2_interacting_tsteps));
    for i = 1:length(zero_deriv2_interacting_tsteps)
        indices(i) = find(tsteps == zero_deriv2_interacting_tsteps(i));
    end
    consecutive_start = diff(indices) == 1;
    if any(consecutive_start)
        consecutive_end = find(diff(consecutive_start)==-1)+1;
        consecutive_indices = [consecutive_start 0] .* indices;
        consecutive_indices(consecutive_end) = indices(consecutive_end);
        consecutive_indices=consecutive_indices(consecutive_indices>0);

        % One average gradient for each run of consecutive indices
        av_deriv1_interacting = zeros(1,length(consecutive_end));
        % Find indices of the final term in each run 
        end_of_run = diff(consecutive_indices)>1;
        end_of_run = [find(end_of_run>0) length(consecutive_indices)];
    else
        fprintf('Less accurate lambda2 esitmate incoming \n')
        consecutive_indices = indices;
        end_of_run = 1:length(indices);
    end
    consecutive_deriv1_interacting = deriv1_interacting(consecutive_indices);
    old_index = 1;
    
    % All other arrays span the timesteps 1/12 hours to 24 hours so the
    % indices should be shifted left to account for the first index being
    % 1/12 and not 0.
    consecutive_indices = consecutive_indices - 1;
    est_lambda2 = zeros(1,length(consecutive_indices));
    for i = 1:length(end_of_run)
        av_deriv1_interacting(i) = ...
            mean(consecutive_deriv1_interacting(old_index:end_of_run(i)));
        
        % Find the average gradient across each consecutive run of 
        % timesteps. Estimate lambda_2 using 
        %       (gradient of interacting curve) * tstep_duration = 
        %       lambda_1*(# free particles at start of t) - 
        %       lambda_2* (# interacting particles at start of t)
        est_lambda2(old_index:end_of_run(i)) = ...
            (mean(est_lambda1) .* freePrtcls_start_of_t(consecutive_indices(i)) ...
            - av_deriv1_interacting(i) .* PARAMETERS.tstep_duration)./ ...
            interactPrtcls_start_of_t(consecutive_indices(i));
    
        old_index = end_of_run(i)+1;
    end
    
    fprintf("Estimates for lambda and EWT from free to stage 1 and from stage 1 to stage 2: \n")
    disp([mean(est_lambda1), 1/mean(est_lambda1)]);
    disp([mean(est_lambda2), 1/mean(est_lambda2)]);
else
    fprintf("lambda1 and lambda2 were not able to be estimated due to there being no steady state or constant gradient \n")
end

% Print actual EWTs for internalisation
fprintf("Actual lambda and EWT from free to stage 1 and from stage 1 to stage 2: \n")
[l1,l2] = input_EWT_from_fraction(...
    PARAMETERS.EWTs_internalise.values(1),...
    PARAMETERS.EWTs_internalise.values(2),...
    PARAMETERS.EWTs_internalise.values(3));
disp([l1, 1/l1]); disp([l2, 1/l2]);

% Print figures
figure(fig2)
figure(fig7)
figure(fig8)
% Save the workspace
save([PARAMETERS.folder_path '\variables.mat']);


%%% SUBFUNCTIONS %%%
% hypo_exp = @(lambda,tdata) 100.*(1 - (lambda(2).*exp(-lambda(1).*tdata) - ...
%     lambda(1).*exp(-lambda(2).*tdata))./(lambda(2)-lambda(1)));

function [lambda1,lambda2] = input_EWT_from_fraction(frac_associated,...
    frac_internalised,num_hours)
% INPUT_EWT_FROM_FRACTION calculates the rates that can be input into
% cells_simulation.m given the desired fractional association
% (frac_associated) and the desired fractional internalisation
% (frac_internalised) after so many hours (num_hours). frac_associated
% and frac_internalised are numbers between 0 and 1.

% Find lambda1 from the CDF of the exponential distribution, F(t), given
% that F(num_hours)=frac_associated
lambda1 = log(1-frac_associated)/(-num_hours); 

% Don't want the root finding function to arrive at lambda1=lambda2
if 1 >= frac_internalised/(frac_associated^2) % lambda1 >= lambda2
    guess = lambda1/3;
else % lambda1 < lambda2
    guess = 3*lambda1;
end

% Create a function from the CDF of the hypoexponential distribution, G(t),
% given that G(num_hours)=frac_internalised and lambda1 is as calculated.
fun = @(lambda2) 1 - frac_internalised - 1/(lambda2-lambda1) * ...
    (lambda2 * exp(-lambda1 * num_hours) - ...
    lambda1 * exp(-lambda2 * num_hours));

% Estimate lambda2 by finding the root of this function (away from lambda1)
lambda2 = fzero(fun,guess);
end
