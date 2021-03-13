function evolution_info = cells_simulation2(total_t, N_initial, P_m, ...
    cycle_probs, DIM, speed)
% cells_simulation simulates the random proliferation and movement of cells
% in a 2.D cell monolayer within a culture dish.
%
%   The section of the culture dish illustrated is representation by a DIM
%   by DIM square matrix/ lattice. The inter-lattice spacing is taken to be
%   the average diameter of a cell (25 micrometers).
%
%   The input arguments are:
%
%       total_t     total number of time steps
%       N_initial   initial number of agents (cells) in the lattice
%       P_m         ability of an agent (cell) to move during each timestep
%       cycle_probs a list of K probabilities of transitioning between
%                   phases in the cell cycle (e.g. phase 1 to phase 2, 
%                   phase 2 to phase 3, ..., phase N-1 to phase N, phase N
%                   to phase 1 having proliferated)
%                       note that these probabilities are discrete approximates
%                       of the exponential waiting time rate
%       DIM         lattice dimensions (DIM by DIM)
%       speed       speed of movie frame playback (how many frames per sec)
%
%   The output arguments are:
%       
%       evolution_info  a structure containing:
%           cell_population     a record of the cell population over time
%           cell_lineage        a record of cell lineage 
%                               [parent cell, daughter cell, generation]
%
%   This is the work of Celia Dowling 8/3/21

close all;

% Seed for efficiency purposes
rng(22)

% Initialise
t=0;
N_t = N_initial; % the number of cells at timestep t
K = length(cycle_probs); % the number of phases in the cell cycle
cell_pop = [N_t zeros(1,total_t)]; % array for cell_population field
cell_lin = zeros(total_t,3); % array for cell_lineage field

% Randomise position of cells on lattice and their phases
culture_dish = zeros(DIM);
rand_lattice_sites = randperm(numel(culture_dish));
cell_sites = zeros(1, numel(culture_dish)); % preallocate space for array
                                            % tracking location of cells
cell_sites(1:N_t) = rand_lattice_sites(1:N_t); % lattice sites for initial
                                               % cells
initial_cycle_phases = datasample(1:K, N_t); % cycle phases for initial
                                             % cells
culture_dish(cell_sites(1:N_t)) = initial_cycle_phases;

% Record each cells' phase in an array that corresponds to the cell_sites
% array (i.e, cell number i at position cell_sites(i) will have phase
% cell_phases(i))
cell_phases = zeros(1, numel(culture_dish)); % preallocate space
culture_dish_as_list = culture_dish(1:numel(culture_dish));
cell_phases(1:N_t) = culture_dish_as_list(cell_sites(1:N_t));


% Prepare movie
cell_fig = figure;
axis tight manual % ensures getframe() returns a consistent size
%ax = gca;
%ax.NextPlot = 'replaceChildren';
cell_movie(total_t + 1) = struct('cdata',[],'colormap',[]);
cell_fig.Visible = 'off';

% Save first movie frame (the initial culture dish)
[rows, cols] = ind2sub(DIM, cell_sites);
c = [cell_phases(1:N_t)'/K zeros(N_t,2)]; % colours for phases
scatter(cols(1:N_t), DIM - rows(1:N_t) + 1, 10*cell_phases(1:N_t), c, ...
    'filled'); % Ensuring the graph represents the lattice positioning
xlim([0.5 DIM+0.5]); ylim([0.5 DIM+0.5]);
title(sprintf('timestep = %d',t));
drawnow
cell_movie(t+1) = getframe(cell_fig);

% Loop stops when timesteps are up or when culture dish is full of cells in
% the last phase of the cell cycle
while t < total_t && all(cell_phases == K*ones(1, numel(culture_dish)), 'all') == 0 
    t = t+1;

    %%% MOVEMENT %%%
    
    % N_t agents are selected with replacement, at random, one at a time
    % and are given a chance to move
    for choice = 1:N_t
        if rand <= P_m  % only a portion of those selected will try to move
            current_site = datasample(cell_sites(1:N_t),1);
            phase = cell_phases(cell_sites == current_site);
            
            % Convert linear indices to [row,column] coordinates 
            [current_row, current_col] = ind2sub(DIM, current_site);
        
            % Each selected agent chooses a random direction to move in 
            % (up, down, left, right)
            delta = randsample([randsample([-1,1],1),0],2);
            new_row_col = [current_row, current_col] + delta;
            
            % For every cell that pops out of the lattice on one side, 
            % another cell pops into the lattice on the opposite side
            new_row_col(new_row_col > DIM) = 1;
            new_row_col(new_row_col < 1) = DIM; 
        
            % Convert [row,column] coordinates to linear indices
            new_site = sub2ind([DIM,DIM], new_row_col(1), new_row_col(2));
        
            % If the new site is vacant, the cell moves
            if culture_dish(new_site) == 0
              culture_dish(current_site) = 0;
              culture_dish(new_site) = phase;
              cell_sites(cell_sites == current_site) = new_site;
            end
        end
    end

    %%% CELL CYCLE & PROLIFERATION %%%
    
    % All N_t agents wait for an approximately exponentially distributed 
    % time (of some rate lambda) before transitioning to the next phase in 
    % the cell cycle (and proliferating if indeed at the end of the cycle).
    % This is simulated by giving all N_t agents in existence at the 
    % beginning of each timestep the chance to transition to the next phase 
    % in the cell cycle via a Bernoulli trial. The probability of success  
    % is equal to this exponential waiting time rate lambda.
    Bernoulli_nums = rand(N_t,1);
    for j = 1:length(Bernoulli_nums)
        parent_site = cell_sites(j);
        old_phase = cell_phases(j);
        
        if Bernoulli_nums(j) <= cycle_probs(old_phase)
            % Convert linear indices to [row, column] coordinates
            [parent_row, parent_col] = ind2sub(DIM, parent_site);
            
            if old_phase == K   % A cell at the end of it's cycle (ready to
                                % proliferate)
                % Each parent cell in their final cycle phase attempting to
                % proliferate chooses a random direction in which to
                % proliferate (create a daughter cell)
                delta = randsample([randsample([-1,1],1),0],2);
                daughter_row_col = [parent_row, parent_col] + delta;
                
                % For every parent cell that proliferates out of the
                % lattice on one side, another cell is born into the
                % lattice on the opposite side
                daughter_row_col(daughter_row_col > DIM) = 1;
                daughter_row_col(daughter_row_col < 1) = DIM;
                
                % Convert [row,column] coordinates to linear indices
                daughter_site = sub2ind([DIM,DIM], daughter_row_col(1), ...
                    daughter_row_col(2));
                
                % If the new site is vacant, the cell proliferates in it
                % and returns to the first phase of the cell cycle as well
                if culture_dish(daughter_site) == 0
                    culture_dish(daughter_site) = 1;
                    N_t = N_t + 1;
                    cell_sites(N_t) = daughter_site;
                    cell_phases(N_t) = 1;
                    new_phase = 1;
                    
                    % Add the [parent cell #, daughter cell #, generation #]
                    % to track lineage
                    parent_not_gen_1 = (cell_lin(:,2) == j);
                    if any(parent_not_gen_1)
                        gen_num = cell_lin(parent_not_gen_1, 3) + 1;
                    else
                        gen_num = 2;
                    end
                    cell_lin(N_t - N_initial,:) = [j, N_t, gen_num];
                else
                    % If a cell in phase K doesn't successfully 
                    % proliferate, it remains in state K
                    new_phase = K;
                end
            else
                % If a cell is not in its final cycle phase, it transitions
                % to the next phase.
                new_phase = old_phase + 1;
            end
            
            % Update chosen cell's phase
            cell_phases(j) = new_phase;
            culture_dish(parent_site) = new_phase;
        end
    end

    
    %%% RECORD %%%
    
    % Save each figure as a frame in the movie illustrating the evolution
    % of the culture dish
    [rows, cols] = ind2sub(DIM, cell_sites);
    c = [cell_phases(1:N_t)'/K zeros(N_t,2)]; % colours for phases
    scatter(cols(1:N_t), DIM - rows(1:N_t) + 1, 10*cell_phases(1:N_t), ...
        c, 'filled'); % Ensuring the graph represents the lattice positioning
    xlim([0.5 DIM+0.5]); ylim([0.5 DIM+0.5]);
    title(sprintf('timestep = %d',t));
    drawnow
    cell_movie(t+1) = getframe(cell_fig);
    
    % Record the cell population at each timestep
    cell_pop(t+1) = N_t;
end

% Play movie
cell_fig.Visible = 'on';
movie(cell_fig, cell_movie, 1, speed);

% Save evolution information into a structure
evolution_info = struct('cell_population', cell_pop(1:t+1), ...
    'cell_lineage', cell_lin(1:(N_t-N_initial),:));
end
