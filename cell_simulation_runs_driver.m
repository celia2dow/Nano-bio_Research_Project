% CELL_SIMULATION_RUNS_DRIVER is a script that drives the cells_simulation 
% simulation with the specified parameters for a given number of iterations 
% and produces summary statistics from the numerous runs.
%
%   This is the work of Celia Dowling 20/05/22
%
%   The input argument for cells_simulation.m is a structure PARAMETERS 
%   which has the following fields that need to be defined by the user:
%
%       simulation_duration (hours) Time length of the simulation.    
%       tstep_duration (hours)      Recommended to be equal to the shortest
%                                   expected waiting time (EWT), usually 
%                                   EWT_move.
%       initial_num_cells           Initial number of cells in the lattice.
%       prtcls_per_site             Number of particles added to the cell
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
%       EWTs_recycle (hours)        Expected waiting time for one particle
%                                   that is interacting with a cell to be
%                                   recycled back into the culture medium.
%                                   Typically particles in stages 1 (bound)
%                                   and L-1 (just prior to internalisation)
%                                   will have the opportunity to exit the
%                                   process (unbind or be recycled before
%                                   therapeutics are delivered). A particle 
%                                   that is internalised (in stage L)
%                                   cannot be recycled. If wanting to
%                                   create this scenario, set all other
%                                   waiting times to be infinity.
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
%   One input argument for create_dosage_distribs.m is a structure FLUORESC
%   which has the following fields that need to be defined by the user:
%
%       stain1_part         The artificual Mean Fluorescence Intensity (MFI) 
%                           value per particle, specific to the primary 
%                           stain (and particular setup) used in the 
%                           equivalent flow cytometry experiment. To 
%                           eliminate intensity transformations, let this 
%                           be equal to 1.
%       stain1_std_dev      The artificial normal standard deviation of an
%                           intensity signal received by flow cytometry 
%                           apparatus. This variance will also scale with 
%                           the square root of the mean signal. To
%                           eliminate variance, let this be equal to 0.
%       stain1_background   The artificial background MFI reading. To 
%                           eliminate background noise, let this be equal
%                           to 0.
%       stain2_part         The artificial MFI value per particle, specific 
%                           to the secondary stain. To eliminate intensity 
%                           transformations, let this be equal to 1.
%       stain2_std_dev      The artificial normal standard deviation of an
%                           intensity signal received by flow cytometry 
%                           apparatus. This variance will also scale with 
%                           the square root of the mean signal. To
%                           eliminate variance, let this be equal to 0.
%       stain2_background   The artificial background MFI reading. To 
%                           eliminate background noise, let this be equal
%                           to 0.
%       alpha               The fraction of the primary stain signal 
%                           (stain1_part) that is received by the apparatus
%                           when a particle has been internalised. It is
%                           assumed that this reduction in signal is is due 
%                           to it being blocked by the cell membrane and
%                           not pH deformation. To eliminate signal
%                           reduction of internalised particles, let this
%                           be equal to 1.
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
%       Dosage_distribsXhours       On X hourly intervals, the frequencies  
%           (.png, .eps)            of cells with certain numbers of 
%                                   internalised (red) and interacting 
%                                   (blue) particles are plotted. These
%                                   frequencies are the averages over
%                                   multiple runs.
%       Analytical_distribs         Plots the average number of
%           (.png, .eps)            internalised, interacting and
%                                   associated (both) particles per cell
%                                   over time, as averaged out over the 
%                                   multiple runs. Overall standard 
%                                   deviation clouds and analytical
%                                   distributions are plotted from
%                                   heuristic estimates of rates.
%       PCC_and_confluence          Plots the pair correlation coefficient
%           (.png, .eps)            for all cells across domains in all
%                                   runs against the average trajectory in
%                                   cell population growth.
%       lambda1_error,              Plots the mean, lower and upper
%       lambda2_error, CC_error     estimates for lambda1, lambda2 and the
%          (.png, .eps)             carrying capacity (if applicable)
%                                   respectively.
%       Fluor_dot_plotsXhours       On X hourly intervals, secondary 
%          (.png, .eps)             fluorescence of each cell 
%                                   (corresponding to the number of 
%                                   particles bound to the cell) is plotted 
%                                   against the primary fluorescence of 
%                                   each cell (corresponding to the number  
%                                   of particles associated to the cell).
%       variables_ ...              A PARAMETERS specific name for a file
%           (.mat)                  containing the MATLAB workspace
%                                   
%   Additionally is printed:
%       estimates for lambda values for all stages of the cell-particle
%       interaction model and for the carrying capacity with upper and
%       lower bounds


% Prepare
clear;
close all;
set(0,'DefaultFigureWindowStyle','docked') % If the setting is 'docked', the frame sizes get weird, also 'normal'

% Seed for efficiency purposes (optional)
%rng(22)

% Choose number of runs
num_runs = 200;

% Parameters relating to the fluorescent intensity dot plot
FLUORESC = struct( ...
    'stain1_part', 1000, ... Without intensity transformations = 1
    'stain1_std_dev', 0, ... Without intensity transformations = 0
    'stain1_background',100, ... Without intensity transformations = 0
    'stain2_part', 2000, ... Without intensity transformations = 1
    'stain2_std_dev', 0, ... Without intensity transformations = 0
    'stain2_background', 100, ... Without intensity transformations = 0
    'alpha', 0.9 ... Without intensity transformations = 1
    );


PARAMETERS = struct( ...
    'simulation_duration',24, ... (hours)
    'tstep_duration',1/6, ... (hours)
    'initial_num_cells', 20, ... 
    'prtcls_per_site', 1000, ...
    'cell_diam', 25, ... (micrometers)
    'culture_dim', 10, ... (cell diameters)
    'culture_media_height', 5,...  (millimeters) 0.5
    'EWT_move', 1/6, ... (hours) 
    'EWTs_proliferate', 24, ... [4,4,4], ... [phase 1, ..., phase K](hours) 
    'EWTs_internalise', struct('input_type', "fraction", ... "fraction" or "EWT" or "prob_and_rates"
    'values', [0.01,24]), ...[0.02,0.01,0.005,24]),... [0.02,0.01,24]), ...[0.01,0.006,24]),... [0.2,0.1,24]), ...% see notes on EWTs_internalise [26.19256, 5.36034], ...[34.62471997,12.52770188], ... 
    'EWTs_recycle', [inf],... (hours)
    'max_prtcls', [inf], ... [stage 1, ..., stage L]
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
    '_pPerC' num2str(PARAMETERS.prtcls_per_site) ...
    '_EWTintern' num2str(PARAMETERS.EWTs_internalise.values) ...
    '_' date_time];
folder_name = folder_name(find(~isspace(folder_name)));

PARAMETERS.folder_path = [pwd '/' date '/' folder_name];

if ~exist(PARAMETERS.folder_path, 'dir')
    mkdir(PARAMETERS.folder_path)
end

% Choose number of hours after which to plot
X = 6;
% Choose to use and ananlyse data on every i(th) timestep (or after every Y
% hours)
Y = 1/6; % hours
ith = Y/PARAMETERS.tstep_duration; % timesteps
% Choose a tolerance for gradient matching
tol = 5E-4;
% Choose maximum number of cell divisions to plot
max_divs = 3;

% Check if the input parameters for internalisation are inappropriate
[it_is,rate_diffus] = is_l1_greater_than_rate_diffus(PARAMETERS);
if it_is
    % Change the input parameters
    warning("The input probability of binding once hit has been replaced with 1")
    if PARAMETERS.EWTs_internalise.input_type == "EWT"
        % Cap the EWT of free to bound at the rate of hitting
        PARAMETERS.EWTs_internalise.values(:,1) = 1/rate_diffus; % hours
    elseif PARAMETERS.EWTs_internalise.input_type == "fraction"
        PARAMETERS.EWTs_internalise.values(:,1) = rate_diffus;
    elseif PARAMETERS.EWTs_internalise.input_type == "prob_and_rates"
        % Cap the probability of binding once hit to 1
        PARAMETERS.EWTs_internalise.values(1) = 1;
    end
end

% Check how many stages are in the cell-particle interaciton model
if PARAMETERS.EWTs_internalise.input_type == "EWT"
    L = length(PARAMETERS.EWTs_internalise.values);
elseif PARAMETERS.EWTs_internalise.input_type == "fraction"
    L = length(PARAMETERS.EWTs_internalise.values) - 1;
elseif PARAMETERS.EWTs_internalise.input_type == "prob_and_rates"
    L = length(PARAMETERS.EWTs_internalise.values);
end

% Check if the number of stages in the cell-particle interaction model are
% consistent across all input
assert(length(PARAMETERS.max_prtcls)==L && ...
    length(PARAMETERS.EWTs_recycle)==L, ...
    ['The number of stages of cell-particle interaction are not consistent ' ...
    'across all inputs (max_prtcls, EWTs_recylce, EWTs_internalise)'])

% Check if the number of hours between plots is greater than the number of
% hours between data points that are to be used.
assert(mod(X,Y)==0, ['The number of hours between plots (X) must be a natural multiple' ...
    ' of the number of hours between used data points (Y)'])

% Run multiple simulations and collect data
total_tsteps = floor(PARAMETERS.simulation_duration/PARAMETERS.tstep_duration);
total.cell_c_o_p = zeros(num_runs*PARAMETERS.culture_dim^2,total_tsteps+1,3);
total.cell_lineage = zeros(num_runs*PARAMETERS.culture_dim^2,total_tsteps+3);
total.cell_population = zeros(1,total_tsteps+1);
runs.cell_pair_cor_coef = zeros(num_runs,total_tsteps+1);
runs.average_cell_c_o_p = zeros(num_runs,total_tsteps+1,3);
cell_end = 0;
for run = 1:num_runs
    EVOLUTION_INFO = cells_simulation(PARAMETERS);
    total.cell_population = total.cell_population + EVOLUTION_INFO.cell_population;
    N_final = EVOLUTION_INFO.cell_population(end); % final number of cells
    cell_start = cell_end + 1;
    cell_end = cell_end + N_final;
    total.cell_c_o_p(cell_start:cell_end,:,:) = EVOLUTION_INFO.cell_c_o_p;
    EVOLUTION_INFO.cell_lineage_history(:,1:2) = ...
        EVOLUTION_INFO.cell_lineage_history(:,1:2) + cell_start - 1;
    EVOLUTION_INFO.cell_lineage_history(1:PARAMETERS.initial_num_cells,1) = ...
        EVOLUTION_INFO.cell_lineage_history(1:PARAMETERS.initial_num_cells,1) - ...
        cell_start + 1;
    total.cell_lineage(cell_start:cell_end,:) = EVOLUTION_INFO.cell_lineage_history;
    runs.average_cell_c_o_p(run,:,:) = ... % Average number of particles per cell in a class...
        [sum(EVOLUTION_INFO.cell_c_o_p(:,:,1),1)./EVOLUTION_INFO.cell_population; ... % free or hit
        sum(EVOLUTION_INFO.cell_c_o_p(:,:,2),1)./EVOLUTION_INFO.cell_population; ... % Interacting
        sum(EVOLUTION_INFO.cell_c_o_p(:,:,3),1)./EVOLUTION_INFO.cell_population]'; % Internalised
    runs.cell_pair_cor_coef(run,:) = EVOLUTION_INFO.cell_pair_cor_coef;
end
%%
% Reduce the data used to being on every i(th) timestep
total.cell_population = total.cell_population(1:ith:(total_tsteps+1));
total.cell_c_o_p = total.cell_c_o_p(1:cell_end,1:ith:(total_tsteps+1),:);
total.cell_lineage = total.cell_lineage(1:cell_end,[1,2,3:ith:(total_tsteps+3)]);
total.confluence = total.cell_population/(PARAMETERS.culture_dim^2 * num_runs);
binrng = 0:ith*PARAMETERS.tstep_duration:PARAMETERS.simulation_duration; % Create bin ranges
fprintf('\nThe average final population of cells in a run:')
disp(total.cell_population(end)/num_runs)

% Find the average pair correlation coefficient
total.cell_pair_cor_coef = [mean(runs.cell_pair_cor_coef(:,1:ith:(total_tsteps+1)), 'all'),...
    var(runs.cell_pair_cor_coef(:,1:ith:(total_tsteps+1)), 0, 'all')];
runs.cell_pair_cor_coef = runs.cell_pair_cor_coef(run,1:ith:(total_tsteps+1));
%%
fprintf("The mean and variance of the pair correlation coefficient in each run:")
format SHORTE
disp(total.cell_pair_cor_coef)

% Plot associated Vs interacting+internalised cells over time
means = zeros(3,floor(total_tsteps/ith)+1);
% Average number interacting particles per cell over all runs in a timestep
means(1,:) = sum(total.cell_c_o_p(:,:,2),1)./total.cell_population; 
% Average number internalised particles per cell over all runs in a timestep
means(2,:) = sum(total.cell_c_o_p(:,:,3),1)./total.cell_population;     
% Average associated particles per cell over time
means(3,:) = means(1,:) + means(2,:);   

% Calculate variances for each timestep
variances = zeros(2,3,floor(total_tsteps/ith)+1);
total.cell_c_o_p_corrected = zeros(cell_end,floor(total_tsteps/ith)+1,2);
for tstep = 0:floor(total_tsteps/ith)
    N_tstep = total.cell_population(tstep+1); 
    
    % Getting rid of the entries corresponding to cells that haven't been
    % created yet, but maintaining the pairing of numbers of interacting/
    % internalised particles for a cell.
    n_interact = total.cell_c_o_p(:,tstep+1,2);
    n_internal = total.cell_c_o_p(:,tstep+1,3);
    with_associated_prtcls = any([n_interact n_internal],2);
    n_interact = n_interact(with_associated_prtcls);
    n_internal = n_internal(with_associated_prtcls);
    num_with = length(n_interact);
    n_interact = [zeros(N_tstep-num_with,1); n_interact];
    n_internal = [zeros(N_tstep-num_with,1); n_internal];
    total.cell_c_o_p_corrected(1:N_tstep,tstep+1,2) = n_interact;
    total.cell_c_o_p_corrected(1:N_tstep,tstep+1,3) = n_internal;

    % VARIANCES BETWEEN ALL CELLS IN ALL RUNS
    % Variances in interacting particles
    variances(1,1,tstep+1) = var(n_interact); 
    % Variances in internalised particles
    variances(1,2,tstep+1) = var(n_internal); 
    % Variances in associated particles
    variances(1,3,tstep+1) = var(n_interact + n_internal);
    
    % VARIANCES BETWEEN RUN MEANS
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
% bound/internalised), secondary fluoerescence against primary
% fluorescence, and deosages split by number of cell divisions after every X 
% hours
create_dosage_distribs(X, floor(total_tsteps/ith), PARAMETERS, total, FLUORESC, max_divs, ith);

% Plot the mean Pair Correlation Coefficient over time (the mean PCC at
% each timestep across all of the run PCCs and the variance in the mean)
PCC_against_confluence;

% When there is no carrying capacity on transitions in between any stages
% other than the last transition (internalisation).
if ~any(PARAMETERS.max_prtcls(1:end-1)~=inf)
    % Heuristic estimation of lambda1, lambda2 and CC if relevant
    heuristic_estimates_2params_CCintern;
    % Plot the mean ass  ociated/internalised/interacting particles over all
    % of the runs over time with standard deviation and analytical
    % distributions from heuristic estimates
    parameter_plots_with_bounds_script;
    % Print figures
end

% Save the workspace
save([PARAMETERS.folder_path '/variables_' num2str(PARAMETERS.prtcls_per_site) ...
    'pPerC_CC' num2str(PARAMETERS.max_prtcls(end)) ...
    '_' num2str(PARAMETERS.EWTs_internalise.values(1)) ...
    '_' num2str(PARAMETERS.EWTs_internalise.values(2)) ...
    '_' date_time '.mat']);

