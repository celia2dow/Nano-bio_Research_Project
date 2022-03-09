function EVOLUTION_INFO = cells_simulation(PARAMETERS)
% CELLS_SIMULATION simulates the random proliferation and movement of cells
% in a 2.D cell monolayer within a culture dish as they interact with and
% internalised particles.
%
%   This is the work of Celia Dowling 23/02/22
%
%   The input argument is a structure PARAMETERS containing the fields:
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
%       folder_path                 The driver function generates this
%                                   depending on the date and the other
%                                   parameters set. It will save all output
%                                   into this path.
%
%   The output argument is a structure EVOLUTION_INFO containing the fields:
%
%       cell_population     A record of the cell population over time.
%       cell_lineage_history A record of cell lineage: 
%                           [parent cell, daughter cell, generation at 0 hours,
%                           generation at 6 hours, generation at 12 hours, etc.]
%       cell_phase_history  A record of the cell phases over time.
%       cell_c_o_p          A 3D array recording the class of particles
%                           (free AND hit, interacting, internalised) on
%                           a cellular basis at each timestep.
%       count_catch         The number of times the catch (ensuring no more
%                           particles are internalised than the carrying
%                           capacity allows for) is used.
%       cell_pair_cor_coef  The coefficient of pair correlation for each
%                           timestep which describes the spatial dependence
%                           of the probability of one lattice site being 
%                           occupied (unity implies independence)
%
%   The simulation also generates and saves a gif of the visual petri dish
%   animation, if one is created, in the folder path specified.


%%% CONVERT EXPECTED WAITING TIMES INTO RATES PER TSTEP %%%

% Rate at which a cell moves.
rate_move = PARAMETERS.tstep_duration/PARAMETERS.EWT_move;

% Rate of particle diffusivity (i.e. of particle hitting).
% Use Stokes-Einstein equation to calculate the diffusion coefficient of
% the particle in the given conditions:
k_B = 1.380649E-23; % Boltzmann constant
T_K = PARAMETERS.temp + 273.15; % Temperature in kelvin
r_m = PARAMETERS.prtcl_rad * (1E-9); % Particle radius in meters
h = PARAMETERS.culture_media_height * (1e-3); % Culture media height in meters
ratio_cell2site_area = pi /4; % Ratio of cell area to the lattice site it sits on
eta = PARAMETERS.viscos;
diffus_coeff = k_B * T_K / (6 * pi * eta * r_m) ... % In m^2 per second
    * (60^2)...                                     % In m^2 per hour
    * PARAMETERS.tstep_duration;                    % In m^2 per timestep
% Divide by the culture height squared to get the rate of transition
% through the bottom of the column on the lattice site, then multiply by
% the ratio of this quare area to the cell area (approximated as a flat
% disk) to get the rate of hitting.
rate_diffus = diffus_coeff * ratio_cell2site_area /(h^2); % In per timestep

% A list of K rates of cells transitioning between phases in the cell 
% proliferation cycle (e.g. phase 1 to phase 2, phase 2 to phase 3, ..., 
% phase K-1 to phase K, phase K to phase 1 having proliferated).
rates_proliferate = PARAMETERS.tstep_duration./PARAMETERS.EWTs_proliferate;
% The number of phases in the cell proliferation cycle.
K = length(rates_proliferate); 

% Check to see if the desired association/internalisation fractions were
% given or the EWTs directly:
if PARAMETERS.EWTs_internalise.input_type == "EWT"
    % A list of L rates of particles transitioning between stages of the cell-
    % particle interaction model unaffected by carrying capacity (e.g. free to 
    % stage 1, stage 1 to stage 2, ..., stage L-1 to stage L or internalised).
    % Can be 1 rate per transition or can be 1 rate per cell phase per
    % transition - columns inicate interaction stage and rows indicate cell
    % phase.
    rates_internalise = PARAMETERS.tstep_duration./...
        PARAMETERS.EWTs_internalise.values;
elseif PARAMETERS.EWTs_internalise.input_type == "fraction"
    % Calculate the rates from the desired observed fraction
    % associated/internalised at confluence without carrying capacity.
    [l1,l2] = input_EWT_from_fraction( ...
        PARAMETERS.EWTs_internalise.values(1), ...
        PARAMETERS.EWTs_internalise.values(2), ...
        PARAMETERS.EWTs_internalise.values(3));
    rates_internalise = PARAMETERS.tstep_duration.*[l1,l2];
end

% Convert the first rate from the rate of binding per timestep to the
% probability of binding once hit by dividing by the rate of hitting.
rates_internalise(:,1) = rates_internalise(:,1)./rate_diffus;
% Check how many stages in the cell-particle interactions model there are
% and whether or not the user has defined different rates for different
% cell-proliferation phases (matrix input) or not (vector input).
[check,L] = size(rates_internalise); 
swtch = 0;
if check > 1
    swtch = 1;
end


%%% INITIALISE VARIABLES, CONSTANTS & VISUAL %%%

N_tstep = PARAMETERS.initial_num_cells; % the number of cells at timestep t
prtcls_initial = PARAMETERS.culture_dim^2 * PARAMETERS.prtcls_per_cell; % the number of particles initially
total_sites = PARAMETERS.culture_dim ^2; % total number of possible positions in petri dish
% Number of timesteps of duration tstep_duration that fit into
% simulation_duration number of hours.
total_tsteps = floor(PARAMETERS.simulation_duration/PARAMETERS.tstep_duration);
tstep = 0; % the current timestep

% Randomise the positions of initial cells on a culture_dim by culture_dim lattice. From 
% now on, all arrays with the name prefix 'cell_' correspond in indexing 
% with each other. I.e. cell_...(i) corresponds to cell i.
culture_dish = zeros(PARAMETERS.culture_dim);
cell_sites = zeros(1, total_sites); % preallocate space
cell_sites(1:N_tstep) = randsample(1:total_sites, N_tstep);
culture_dish(cell_sites(1:N_tstep))=1; % place cells on culture dish
                                               
% Randomise the phases of initial cells.
cell_phases = zeros(1, total_sites); % preallocate space
cell_phases(1:N_tstep) = datasample(1:K, N_tstep);

% Initialise arrays recording the number of particles in different stages
% (per column) of interaction with each cell (per row). The first column is
% to simulate the number of particles hovering around one cell initially.
cell_prtcls = zeros(total_sites, L+1);
% Initialise array for recording the particles that transition each
% timestep.
successes_in_stage = zeros(total_sites, L);

% Initialise arrays/constants for EVOLUTION_INFO fields.
cell_population = [N_tstep zeros(1,total_tsteps)];
cell_lineage_history = [zeros(N_tstep,1) (1:N_tstep)' ones(N_tstep,2) ...
    zeros(N_tstep,ceil(PARAMETERS.simulation_duration*PARAMETERS.tstep_duration)-1); ...
    zeros(total_sites,ceil(PARAMETERS.simulation_duration*PARAMETERS.tstep_duration)+3)];
lineage_colmn = 4;
cell_phase_history = [cell_phases' zeros(total_sites, total_tsteps)]; 
num_neighs_occ = count_occupied_neighs(PARAMETERS.culture_dim, N_tstep, cell_sites, culture_dish);
normalising_coeff = 4 * (PARAMETERS.culture_dim^2) * N_tstep * (N_tstep-1)/...
        ((PARAMETERS.culture_dim^2) * (PARAMETERS.culture_dim^2-1));
cell_pair_cor_coef = [num_neighs_occ/normalising_coeff zeros(1,total_tsteps)];
% Record free particles, interacting particles and internalised particles
% per timestep in system.
tally_prtcls = [prtcls_initial zeros(1, total_tsteps); zeros(2, total_tsteps + 1)];
% Record class of particles on a cell basis per hour.
cell_c_o_p = zeros(total_sites,total_tsteps+1,3);
cell_c_o_p(1:N_tstep,1,1) = PARAMETERS.prtcls_per_cell;
% The number of times the binomial distribution overdraws particles to
% internalise.
count_catch = 0;

if PARAMETERS.visual
    % Prepare movie.
    petri_fig = figure;
    axis tight manual % ensures getframe() returns a consistent size
    %ax = gca;
    %ax.NextPlot = 'replaceChildren';
    petri_movie(total_tsteps + 1) = struct('cdata',[],'colormap',[]);
    petri_fig.Visible = 'off';
    
    % Save first movie frame (the initial culture dish).
    draw_frame(PARAMETERS,N_tstep, tstep, L, cell_sites, cell_phases, cell_prtcls)
    frame0 = getframe(petri_fig);
    petri_movie(tstep+1) = frame0;
    
    % Prepare gif and save first frame.
    gifpath = [PARAMETERS.folder_path '\simulation' num2str(PARAMETERS.visual) '.gif'];
    image0 = frame2im(frame0);
    [imind,cm]=rgb2ind(image0,256);
    imwrite(imind,cm,gifpath,'gif','Loopcount',inf)    
end


% Loop stops when timesteps are up or when culture dish is full of cells
% with their maximum capacity of particles.
while tstep < total_tsteps && ~all(tally_prtcls(:,tstep+1) == [0; ...
        prtcls_initial-min(total_sites*PARAMETERS.max_prtcls(L),prtcls_initial); ...
        min(total_sites*PARAMETERS.max_prtcls(L),prtcls_initial)],'all')
    tstep = tstep+1;
    
    
    %%% CELL-PARTICLE INTERACTIONS %%%
    
    % Assuming that particles are constantly evenly dispersing themselves  
    % around the petri dish/matrix, then the number of particles hovering  
    % over a matrix site at the beginning of each timestep cell_prtcls(:,1)
    % will be the total number of particles that are not internalised or 
    % currently bound to/interacting with cells divided amongst the number 
    % of matrix sites.
    tally_prtcls(:,tstep+1) = tally_prtcls(:,tstep);
    num_per_site = tally_prtcls(1,tstep+1)/total_sites; % This may be a fraction
    cell_prtcls(:,1) = floor(num_per_site); % This may be 0
    if mod(num_per_site, floor(num_per_site)) ~= 0 % Check for remainder
        leftovers = tally_prtcls(1,tstep+1) - total_sites * floor(num_per_site);
        % Randomly distribute the remaining particles around the matrix and
        % if a cell lies on such a site, add to its hovering particle
        % count.
        prtcl_sites = randsample(1:total_sites, leftovers);
        [~, indices] = intersect(cell_sites, prtcl_sites);
        cell_prtcls(indices,1) = cell_prtcls(indices,1) + 1;
    end
    
    % We draw the number of free particles that hit a cell and attempt to
    % bind to it (i.e. that are available in stage 0 of the cell-particle
    % interaction model) from a Binomial distribution where the probability
    % of a successful event = rate_prtcl_diffusivity * tstep_duration and 
    % the number of trials = the number of particles hovering on the lattice 
    % site (i.e. currently in the first column of our cell_prtcls array).
    cell_prtcls(1:N_tstep,1) = binornd(cell_prtcls(1:N_tstep,1),rate_diffus);
   
    % Now that every column of the cell_prtcls array represents the
    % population of particles in a stage of interaction with the indexed
    % cell (stage 0 or hit to stage L or internalised), we can simply use 
    % another binomial distribution to find how many transition to the next 
    % stage, but this time using the dynamic rate of transition dependent 
    % on the number of particles in the following stage (relative to the 
    % carrying capacity of that following stage) to determine how many in 
    % each stage are successful in graduating to the next stage. Note that
    % those that successfully transition from stage 0 (hit) to stage 1
    % (bound) remain in stage 1, and those that don't succeed go back to
    % being free particles. I.e., particles that hit have a probability of
    % binding and otherwise move away from the cell again.
    for stage = 1:L
        current_carrying_total = cell_prtcls(1:N_tstep,stage+1);
        % Check if internalisation rates have been defined per cell 
        % proliferation cycle phase or it they're simply uniform over all 
        % phases. Dynamically affect the rate of the transition into a
        % stage, relative to its carrying capacity
        if swtch % Dependent on phase
            dynamic_prtcl_rates =  rates_internalise(cell_phases(1:N_tstep),...
                stage) .* (1-current_carrying_total./PARAMETERS.max_prtcls(stage));      
        else % Not dependent on phase
            dynamic_prtcl_rates =  rates_internalise(stage) .* (1-...
                current_carrying_total./PARAMETERS.max_prtcls(stage));
        end
        % Due to wanting to keep the total number of particles allowed to
        % transition from a stage in one timestep static in spite of 
        % successful transitions on the same timestep, save the number of 
        % successes without yet updating the cell_prtcl array
        successes_in_stage(1:N_tstep,stage) = binornd(cell_prtcls(1:N_tstep,...
            stage), dynamic_prtcl_rates);
        satisfy = (successes_in_stage(1:N_tstep,stage)+cell_prtcls(1:N_tstep,...
            stage+1)>PARAMETERS.max_prtcls(stage));
        if any(satisfy)
            successes_in_stage(satisfy,stage)=PARAMETERS.max_prtcls(stage)-cell_prtcls(satisfy, stage+1);
            count_catch = count_catch + 1;
        end
    end
    
    % Update the total/per cell number of particles in each stage of
    % the cell-particle interaction model
    cell_prtcls(1:N_tstep,1:L) = cell_prtcls(1:N_tstep,1:L) - successes_in_stage(1:N_tstep,:);
    cell_prtcls(1:N_tstep,2:L+1) = cell_prtcls(1:N_tstep,2:L+1) + successes_in_stage(1:N_tstep,:);
    
    % Update the total class of particles
    tally_prtcls(2,tstep+1) = sum(cell_prtcls(1:N_tstep,2:L),'all'); % interacting
    tally_prtcls(3,tstep+1) = sum(cell_prtcls(1:N_tstep,L+1)); % internalised
    tally_prtcls(1,tstep+1) = prtcls_initial - sum(tally_prtcls(2:3,tstep+1)); % free
    
    
    %%% MOVEMENT %%%
    
    % N_tstep cells are selected with replacement, at random, one at a time
    % and are given a chance to move
    for choice = 1:N_tstep
        % Only a portion of those selected will try to move
        if rand <= rate_move
            % Choose a random cell
            current_site = datasample(cell_sites(1:N_tstep),1);
            
            % Call local function to choose a new site to move to
            new_site = choose_adjacent_site(PARAMETERS.culture_dim, current_site);
        
            % If the new site is vacant, the cell moves
            if culture_dish(new_site) == 0
              culture_dish(current_site) = 0;
              culture_dish(new_site) = 1;
              cell_sites(cell_sites == current_site) = new_site;
            end
        end
    end

    
    %%% CELL PROLIFERATION CYCLE %%%
    
    N_tstep_static = N_tstep; % keep N_tstep fixed for the attempted prolif. events
    
    % N_t_static cells are selected with replacement, at random, one at a 
    % time and are given a chance to transition to the next phase in the
    % cell proliferation cycle and perhaps proliferate.    
    for choice = 1:N_tstep_static
        % Choose a random cell
        parent_site = datasample(cell_sites(1:N_tstep),1);
        old_phase = cell_phases(cell_sites == parent_site);
        
        % Only a portion of those selected will try to transition
        if rand <= rates_proliferate(old_phase)
            % For cells in their final proliferation cycle phase
            if old_phase == K
                % Call local function to choose daughter_site in which to 
                % attempt to proliferate into (create a daughter cell)
                daughter_site = choose_adjacent_site(PARAMETERS.culture_dim, parent_site);
        
                % If the new site is vacant, the cell proliferates into it
                % and returns to the first phase of the cell proliferation
                % cycle as well
                if culture_dish(daughter_site) == 0
                    culture_dish(daughter_site) = 1;
                    N_tstep = N_tstep + 1;
                    cell_sites(N_tstep) = daughter_site;
                    cell_phases(N_tstep) = 1;
                    cell_phases(cell_sites == parent_site) = 1;
                    
                    % Add the [parent cell #, daughter cell #,..., generation #]
                    % to track lineage - note that both the "parent" and
                    % "daughter" cells need to have their generations
                    % updated.
                    parent_cell_num = find(cell_sites == parent_site);
                    gen_num = cell_lineage_history(cell_lineage_history(:,2) == parent_cell_num, lineage_colmn) + 1;
                    cell_lineage_history(N_tstep,1:2) = [parent_cell_num, N_tstep];
                    cell_lineage_history([parent_cell_num,N_tstep],lineage_colmn) = gen_num;                    
                    
                    % The inheritance of internalised/ interacting
                    % nanoparticles from a cell with n_int nanoparticles 
                    % per cell-particle interaction stage follows a 
                    % binomial distribution.
                    for stage = 2:L+1
                        n_int = cell_prtcls(parent_cell_num, stage);
                        if n_int >= 1
                            % Possible inheritance numbers
                            x_pos = 0:n_int;
                            % Define pdf for these possible numbers
                            pmf = zeros(1,n_int + 1);
                            for i = x_pos
                                % Probability of daughter cell 1 inheriting
                                % i nanoparticles from parent cell in stage
                                pmf(i+1) = 0.5 * nchoosek(n_int,i) *...
                                    (PARAMETERS.prob_inherit^i * ...
                                    (1-PARAMETERS.prob_inherit)^(n_int-i) + ...
                                    PARAMETERS.prob_inherit^(n_int-i) * ...
                                    (1-PARAMETERS.prob_inherit)^i);
                            end
                            % Generate an inherited number of interacting/ 
                            % internalised particles for daughter cell 1 and 
                            % derive for daughter cell 2
                            daught_1 = rand_pmf(x_pos, pmf, 1);
                            daught_2 = n_int - daught_1;
                        else
                            daught_1 = 0;
                            daught_2 = 0;
                        end
                        % Update record of particles inside/ interacting 
                        % with cells
                        cell_prtcls(parent_cell_num, stage) = daught_1;
                        cell_prtcls(N_tstep, stage) = daught_2;
                    end
                else 
                    % If the new site isn't vacant, the cell doesn't change
                    % phase or proliferate
                    cell_phases(cell_sites == parent_site) = K;
                end
            else
                % If a cell is not in its final proliferation cycle phase, 
                % it transitions to the next phase.
                cell_phases(cell_sites == parent_site) = old_phase + 1;
            end
        end
    end
    
    % If a new block of 6 hours is begun, start recording cell_lineage in a
    % new column
    hours = tstep*PARAMETERS.tstep_duration;
    if mod(hours,6)==0
        cell_lineage_history(:,lineage_colmn+1)=cell_lineage_history(:,lineage_colmn);
        lineage_colmn = lineage_colmn + 1;
    end
    
    %%% RECORD %%%
    
    if PARAMETERS.visual
        clf
        % Save each figure as a frame in the movie illustrating the evolution
        % of the culture dish
        draw_frame(PARAMETERS,N_tstep, tstep, L, cell_sites, cell_phases, cell_prtcls)
        frameNt = getframe(petri_fig);
        petri_movie(tstep+1) = frameNt;
        % Save the gif frame.
        imageN = frame2im(frameNt);
        [imind,cm]=rgb2ind(imageN,256);
        imwrite(imind,cm,gifpath,'gif','WriteMode','append')   
    end
    
    % Record the class of particles on a cell basis every timestep
    cell_c_o_p(1:N_tstep,tstep+1,1) = ...
        cell_prtcls(1:N_tstep,1); % Particles that are free OR hit
    cell_c_o_p(1:N_tstep,tstep+1,2) = ...
        sum(cell_prtcls(1:N_tstep,2:L),2); % Particles that interact
    cell_c_o_p(1:N_tstep,tstep+1,3) = ...
        cell_prtcls(1:N_tstep,L+1); % Particles that are internalised
    
    % Record the cell phases at each timestep
    cell_phase_history(1:N_tstep,tstep+1) = cell_phases(1:N_tstep)';
    
    % Record the cell population at each timestep
    cell_population(tstep+1) = N_tstep;
    
    % Record the pair correlation coefficient of cells at each timestep
    num_neighs_occ = count_occupied_neighs(PARAMETERS.culture_dim, N_tstep, cell_sites, culture_dish);
    normalising_coeff = 4 * (PARAMETERS.culture_dim^2) * N_tstep * (N_tstep-1)/...
        ((PARAMETERS.culture_dim^2) * (PARAMETERS.culture_dim^2-1));
    cell_pair_cor_coef(tstep+1) = num_neighs_occ/normalising_coeff;
end

% Print cell particles
%disp(cell_prtcls((1:N_t),:));

if PARAMETERS.visual
    % Play movie
    petri_fig.Visible = 'on';
    movie(petri_fig, petri_movie, 1, PARAMETERS.vid_speed);
end

% Save evolution information into a structure
EVOLUTION_INFO = struct('cell_population', cell_population(1:tstep+1), ...
    'cell_lineage_history', cell_lineage_history(1:N_tstep,:), ...
    'cell_phase_history', cell_phase_history(1:N_tstep,1:tstep+1),...
    'cell_c_o_p', cell_c_o_p(1:N_tstep,:,:), ...
    'count_catch', count_catch, ...
    'cell_pair_cor_coef', cell_pair_cor_coef);
end

