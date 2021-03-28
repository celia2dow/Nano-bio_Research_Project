function evolution_info = cells_simulation(total_t, N_initial, DIM, ...
    max_prtcl, P_m, cycle_probs, rate_interacts, base_prtcl_probs, speed)
% CELLS_SIMULATION simulates the random proliferation and movement of cells
% in a 2.D cell monolayer within a culture dish.
%
%   The section of the culture dish illustrated is represented by a DIM by
%   by DIM square matrix/ lattice. The inter-lattice spacing is taken to be
%   the average diameter of a cell (25 micrometers).
%
%   The input arguments are:
%
%       total_t     total number of timesteps
%       N_initial   initial number of cells in the lattice
%       DIM         lattice dimensions (DIM by DIM)
%       max_prtcl   maximum number of particles that a cell can internalise
%       P_m         probability that a cell moves in 1 timestep
%       cycle_probs a list of K probabilities of transitioning between
%                   phases in the cell proliferation cycle in 1 timestep 
%                   (e.g. phase 1 to phase 2, phase 2 to phase 3, ...,  
%                   phase N-1 to phase N,phase N to phase 1 having 
%                   proliferated)
%                       note that these probabilities are discrete 
%                       approximates of the exponential waiting time rates
%       rate_interacts  the poisson rate at which particles interact with a
%                       given cell (number of particles per timestep)
%       base_prtcl_probs    a list of L base probabilities of particles 
%                           transitioning between stages of the 
%                           cell-particle interaction model in 1 timestep
%       speed       speed of movie frame playback (frames per sec)
%
%   The output arguments are:
%       
%       evolution_info  a structure containing:
%           cell_population     a record of the cell population over time
%           cell_lineage        a record of cell lineage 
%                               [parent cell, daughter cell, generation]
%           uninternalised_particles    a record of the number of particles
%                                       that have not been internalised 
%                                       over time
%
%   This is the work of Celia Dowling 22/3/21


%%% PREPARE %%%

close all;

% Seed for efficiency purposes
rng(22)

% Initialise variables and constants
t=0; % the timestep
N_t = N_initial; % the number of cells at timestep t
prtcls_initial = N_t * 100; % the number of particles initially
K = length(cycle_probs); % the number of phases in the cell proliferation cycle
L = length(base_prtcl_probs); % "" in the cell-particle interaction model
total_sites = DIM * DIM; % total number of possible positions in petri dish

% Initialise arrays for evolution_info fields
population = [N_t zeros(1,total_t)]; % cell_population
lineage = [zeros(N_t,1) (1:N_t)' ones(N_t,1); zeros(total_t,3)]; % cell_lineage
free_prtcls = [prtcls_initial zeros(1, total_t)]; % uninternalised_particles

% Randomise the positions of initial cells on a DIM by DIM lattice. From 
% now on, all arrays with the name suffix 'cell_' correspond in indexing 
% with each other. I.e. cell_...(i) corresponds to cell i.
culture_dish = zeros(DIM);
cell_sites = zeros(1, total_sites); % preallocate space
cell_sites(1:N_t) = randsample(1:total_sites, N_t);
                                               
% Randomise the phases of initial cells
cell_phases = zeros(1, total_sites); % preallocate space
cell_phases(1:N_t) = datasample(1:K, N_t);

% Initialise arrays recording the number of particles in different stages
% (per column) of interaction with each cell (per row). The first column is
% to simulate the number of particles hovering around one cell initially.
cell_prtcls = [100*ones(total_sites, 1) zeros(total_sites, L)];

% Prepare movie
petri_fig = figure;
axis tight manual % ensures getframe() returns a consistent size
%ax = gca;
%ax.NextPlot = 'replaceChildren';
petri_movie(total_t + 1) = struct('cdata',[],'colormap',[]);
petri_fig.Visible = 'off';

% Save first movie frame (the initial culture dish)
draw_frame(N_t, DIM, t, K, cell_sites, cell_phases)
petri_movie(t+1) = getframe(petri_fig);

% Loop stops when timesteps are up or when culture dish is full
while t < total_t && all(culture_dish, 'all') == 0 
    t = t+1;
    
    %%% CELL-PARTICLE INTERACTIONS %%%
    
    % Allow each cell to attempt to interact with x particles where x is
    % drawn from a Poisson distribution of rate (mean) equal to
    % rate_interacts.
    cell_num_attempts = poissrnd(rate_interacts, [1,N_t]);
    for cell = 1: N_t
        x = cell_num_attempts(cell);
        prtcls_per_stage = cell_prtcls(cell,:);
        
        % Assuming about 100 particles are available to a cell at once, and
        % given so many may already be in different stages of the
        % cell-particle interaction model, create a weighted distribution 
        % of the cell-particle interaction model stages.
        num_internal = prtcls_per_stage(L+1);
        hovering_prtcls = zeros(1, 100-num_internal);
        num_previous_stage = 0;
        for stage = 1:L
            num_in_stage = prtcls_per_stage(stage);
            hovering_prtcls(num_previous_stage + 1: num_in_stage + ...
                num_previous_stage) = stage .* ones(1, num_in_stage);
            num_previous_stage = num_in_stage;
        end
        
        % Choose the current stages of the x particles which attempt to 
        % progress through the cell-particle interaction model from this 
        % weighted distribution.
        attempts = datasample(hovering_prtcls, x);
        
        % Run a Bernoulli trial for each particle one at a time using the 
        % dynamic transition probability to see whether or not it succeeds 
        % in its attempt.
        Bernoulli_trials = rand(1,x);
        successes = zeros(1,x);
        for prtcl = 1:x
            stage = attempts(prtcl);
            dynamic_prtcl_prob =  base_prtcl_probs(stage) * ...
                (1-num_internal/max_prtcl);
            successes(prtcl) = stage * (Bernoulli_trials(prtcl) ...
                <= dynamic_prtcl_prob);
            if stage == L && successes(prtcl) ~=0
                num_internal = num_internal + 1;
            end
        end
        
        % Update the total number of particles in each stage of the
        % cell-particle interaction model
        for stage = 1:L
            success_in_stage = sum(successes == stage);
            prtcls_per_stage(stage) = prtcls_per_stage(stage) - success_in_stage;
            prtcls_per_stage(stage + 1) = prtcls_per_stage(stage + 1) ...
                + success_in_stage;
        end
        cell_prtcls(cell,:) = prtcls_per_stage;
    end
    
    %%% MOVEMENT %%%
    
    % N_t cells are selected with replacement, at random, one at a time
    % and are given a chance to move
    for choice = 1:N_t
        % Only a portion of those selected will try to move
        if rand <= P_m
            % Choose a random cell
            current_site = datasample(cell_sites(1:N_t),1);
            
            % Call local function to choose a new site to move to
            new_site = choose_adjacent_site(DIM, current_site);
        
            % If the new site is vacant, the cell moves
            if culture_dish(new_site) == 0
              culture_dish(current_site) = 0;
              culture_dish(new_site) = 1;
              cell_sites(cell_sites == current_site) = new_site;
            end
        end
    end

    %%% CELL PROLIFERATION CYCLE %%%
    
    N_t_static = N_t; % keep N_t fixed for the attempted prolif. events
    
    % N_t_static cells are selected with replacement, at random, one at a 
    % time and are given a chance to transition to the next phase in the
    % cell proliferation cycle and perhaps proliferate.    
    for choice = 1:N_t_static
        % Choose a random cell
        parent_site = datasample(cell_sites(1:N_t),1);
        old_phase = cell_phases(cell_sites == parent_site);
        
        % Only a portion of those selected will try to transition
        if rand <= cycle_probs(old_phase)
            % For cells in their final proliferation cycle phase
            if old_phase == K
                % Call local function to choose daughter_site in which to 
                % attempt to proliferate into (create a daughter cell)
                daughter_site = choose_adjacent_site(DIM, parent_site);
        
                % If the new site is vacant, the cell proliferates into it
                % and returns to the first phase of the cell proliferation
                % cycle as well
                if culture_dish(daughter_site) == 0
                    culture_dish(daughter_site) = 1;
                    N_t = N_t + 1;
                    cell_sites(N_t) = daughter_site;
                    cell_phases(N_t) = 1;
                    cell_phases(cell_sites == parent_site) = 1;
                    
                    % Add the [parent cell #, daughter cell #, generation #]
                    % to track lineage
                    parent_cell_num = find(cell_sites == parent_site);
                    gen_num = lineage(lineage(:,2) == parent_cell_num, 3) + 1;
                    lineage(N_t,:) = [parent_cell_num, N_t, gen_num];
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
    
    %%% RECORD %%%
    
    % Save each figure as a frame in the movie illustrating the evolution
    % of the culture dish
    draw_frame(N_t, DIM, t, K, cell_sites, cell_phases)
    petri_movie(t+1) = getframe(petri_fig);
    
    % Record the cell population and free particles at each timestep
    population(t+1) = N_t;
    free_prtcls(t+1) = prtcls_initial - sum(cell_prtcls(:,L+1));
end

% Print cell particles
disp(cell_prtcls((1:N_t),:));

% Play movie
petri_fig.Visible = 'on';
movie(petri_fig, petri_movie, 1, speed);

% Save evolution information into a structure
evolution_info = struct('cell_population', population(1:t+1), ...
    'cell_lineage', lineage(1:N_t,:), ...
    'uninternalised_particles', free_prtcls(1:t+1));
end


function new_site = choose_adjacent_site(DIM, current_site)
% CHOOSE_ADJACENT_SITE Given the current_site in a DIM by DIM lattice, 
% choose a new_site that is directly adjacent (up, down, left, right)

% Convert linear indices to [row,column] coordinates 
[current_row, current_col] = ind2sub(DIM, current_site);
        
% Each selected cell chooses a random direction to move in
delta = randsample([randsample([-1,1],1),0],2);
new_row_col = [current_row, current_col] + delta;
            
% For every cell that pops out of the lattice on one side, another cell
% pops into the lattice on the opposite side
new_row_col(new_row_col > DIM) = 1;
new_row_col(new_row_col < 1) = DIM; 
        
% Convert [row, column] coordinates to linear indices
new_site = sub2ind([DIM,DIM], new_row_col(1), new_row_col(2));
end


function draw_frame(N_t, DIM, t, K, cell_sites, cell_phases)
% DRAW_FRAME Save each figure as a frame in the movie illustrating the 
% evolution of the culture dish as cells move, proliferate and internalise
% particles

% Draw a scatter plot of the culture dish with:
%       size and colour indicating phase of cell proliferation cycle
[rows, cols] = ind2sub(DIM, cell_sites);
c = [cell_phases(1:N_t)'/K zeros(N_t,2)]; % colours for phases
scatter(cols(1:N_t), DIM - rows(1:N_t) + 1, 10*cell_phases(1:N_t), c, ...
    'filled'); % Ensuring the graph represents the lattice positioning
xlim([0.5 DIM+0.5]); ylim([0.5 DIM+0.5]);
title(sprintf('timestep = %d',t));
drawnow
end
