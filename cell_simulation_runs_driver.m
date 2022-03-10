% CELL_SIMULATION_RUNS_DRIVER is a script that drives the cells_simulation 
% simulation with the specified parameters for a given number of iterations 
% and produces summary statistics from the numerous runs.
%
%   This is the work of Celia Dowling 11/03/22
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
%                                   monolayer per lattice site and thus the
%                                   total number of particles available to
%                                   each cell regardless of confluence
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
%                                   either "fraction", "EWT", or
%                                   "prob_and_rates"
%       EWTs_internalise.values     if input_type == "fraction":
%                                       [fraction associated, fraction
%                                       internalised, number of hours of 
%                                       incubation after which these
%                                       fractions are desired to be
%                                       observed] ... at confluence without
%                                       carrying capacity. 
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
%                                   if input_type == "prob_and_rates":
%                                       [probability of transitioning from
%                                       stage 0 (hit) to stage 1 (bound) 
%                                       (i.e. probability of binding once
%                                       hit), rate from stage 1 to stage 2
%                                       (per hour), rate from stage 2 to 
%                                       stage 3 (per hour), ..., rate from
%                                       stage L-1 to stage L (internalised,
%                                       per hour)] ... this will result in
%                                       a rate from free to stage 1 of
%                                       EWTs_internalise.values(1) *
%                                       rate_diffus (per timestep)
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
num_runs = 10;

% Choose a tolerance for gradient matching
tol = 1E-1;
tol_l2 = 5E-4;

% Parameters relating to the fluorescent intensity dot plot
stain1_part = 1000; % Without intensity transformations = 1
stain1_std_dev = 10; % Without intensity transformations = 0
stain1_background = 100; % Without intensity transformations = 0
stain2_part = 2000; % Without intensity transformations = 1
stain2_std_dev = 10; % Without intensity transformations = 0
stain2_background = 100; % Without intensity transformations = 0
alpha = 0.9; % Without intensity transformations = 1

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
    'EWTs_internalise', struct('input_type', "fraction", ... "fraction" or "EWT" or "prob_and_rates"
    'values', [0.5,0.01,24]),...%[0.01,0.006,24]), ... see notes on EWTs_internalise [26.19256, 5.36034], ...[34.62471997,12.52770188], ... 
    'max_prtcls', [inf,inf], ... [stage 1, ..., stage L]
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

% Check if the input parameters for internalisation are inappropriate
[it_is,rate_diffus] = is_l1_greater_than_rate_diffus(PARAMETERS);
if it_is
    % Change the input parameters
    fprintf("The input probability has been replaced with 1\n")
    if PARAMETERS.EWTs_internalise.input_type == "EWT"
        % Cap the EWT of free to bound at the rate of hitting
        PARAMETERS.EWTs_internalise.values(:,1) = 1/rate_diffus; % hours
    elseif PARAMETERS.EWTs_internalise.input_type == "fraction"
        PARAMETERS.EWTs_internalise.input_type = "prob_and_rates";
        l1 = rate_diffus; % per hour
        guess = 3*l1;
        frac_internalised = PARAMETERS.EWTs_internalise.values(2);
        num_hours = PARAMETERS.EWTs_internalise.values(3);
        % Create a function from the CDF of the hypoexponential distribution, G(t),
        % given that G(num_hours)=frac_internalised and lambda1 is as calculated.
        fun = @(l2) 1 - frac_internalised - 1/(l2-l1) * ...
            (l2 * exp(-l1 * num_hours) - l1 * exp(-l2 * num_hours));
        % Estimate lambda2 by finding the root of this function (away from lambda1)
        l2 = fzero(fun,guess); % per hour
        if l2 > 1/PARAMETERS.tstep_duration || isnan(l2)
            % Just in case the desired fraciton of internalisaiton is
            % impossible with the reduced lambda1, then cap lambda2 at the
            % reciprocal of the timestep duration
            l2 = 1/PARAMETERS.tstep_duration;
            fprintf('Other EWTs have been reduced to the timestep duration\n')
        end
        PARAMETERS.EWTs_internalise.values = [1,l2];
    elseif PARAMETERS.EWTs_internalise.input_type == "prob_and_rates"
        % Cap the probability of binding once hit to 1
        PARAMETERS.EWTs_internalise.values(1) = 1;
    end
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
fprintf('\nThe average final population of cells in a run:')
disp(total.cell_population(end)/num_runs)

% Find the average pair correlation coefficient
total.cell_pair_cor_coef = [mean(runs.cell_pair_cor_coef(:,:), 'all'),...
        var(runs.cell_pair_cor_coef(:,:), 0, 'all')];
fprintf("The mean and variance of the pair correlation coefficient in each run:")
format SHORTE
disp(total.cell_pair_cor_coef)

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


close all

% Plot the dosage distribution (the frequency of having so many particles
% bound/internalised) after every X hours
create_dosage_distribs(4, total_tsteps, PARAMETERS, total)

% Plot the mean Pair Correlation Coefficient over time (the mean PCC at
% each timestep across all of the run PCCs and the variance in the mean)
PCC_against_confluence

[~,~,L]=size(total.cell_c_o_p);

% In the 2 parameter case with or without a carrying capacity on the last 
% transition (internalisation).
if L==3 && ~any(PARAMETERS.max_prtcls(1:end-1)~=inf)
    % Heuristic estimation of lambda1, lambda2 and CC if relevant
    heuristic_estimates_2params_CCintern
    % Plot the mean associated/internalised/interacting particles over all
    % of the runs over time with standard deviation and analytical
    % distributions from heuristic estimates
    parameter_plots_with_bounds_script
    % Print figures
end

% Save the workspace
save([PARAMETERS.folder_path '/variables_' num2str(PARAMETERS.prtcls_per_cell) ...
    'pPerC_CC' num2str(PARAMETERS.max_prtcls(end)) ...
    '_' num2str(PARAMETERS.EWTs_internalise.values(1)) ...
    '_' num2str(PARAMETERS.EWTs_internalise.values(2)) ...
    '_' date_time '.mat']);

