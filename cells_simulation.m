function evolution_info = cells_simulation(total_t, N_initial, P_m, ...
    P_p, DIM, speed)
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
%       P_p         " " " to proliferate " " "
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
N_t = N_initial;
cell_pop = [N_initial zeros(1,total_t)];
cell_lin = zeros(total_t,3);

% Randomise position of cells on lattice
culture_dish = zeros(DIM);
rand_lattice_sites = randperm(numel(culture_dish));
cell_sites = zeros(1, numel(culture_dish)); % preallocate space
cell_sites(1:N_initial) = rand_lattice_sites(1:N_initial);
culture_dish(cell_sites(1:N_initial)) = 1;

% Prepare movie
cell_fig = figure;
axis tight manual % ensures getframe() returns a consistent size
cell_movie(total_t + 1) = struct('cdata',[],'colormap',[]);
cell_fig.Visible = 'off';

% Save first movie frame (the initial culture dish)
spy(culture_dish);
title(sprintf('timestep = %d',t));
drawnow
cell_movie(t+1) = getframe;

% Loop stops when timesteps are up or when culture dish is full
while t < total_t && all(culture_dish, 'all') == 0 
    t = t+1;

    %%% MOVEMENT %%%
    
    % N_t agents are selected with replacement, at random, one at a time
    % and are given a chance to move
    selected_to_move = datasample(cell_sites(1:N_t), N_t);
    which_try_to_move = selected_to_move(rand(N_t,1) <= P_m);
    attempt_directions = datasample(1:4, length(which_try_to_move));
    
    for i = 1:length(which_try_to_move)
        % Convert linear indices to [row,column] coordinates 
        current_site = which_try_to_move(i);
        [current_row, current_col] = ind2sub(DIM, current_site);
        
        % Each selected agent chooses a random direction to move in (up,
        % down, left, right)
        direction = attempt_directions(i);
        switch direction
            case 1 % try to move up
                new_col = current_col;
                new_row = current_row - 1;
            case 2 % try to move right
                new_row = current_row;
                new_col = current_col + 1;
            case 3 % try to move down
                new_col = current_col;
                new_row = current_row + 1;
            case 4 % try to move left
                new_row = current_row;
                new_col = current_col - 1;
        end    
        new_row_col = [new_row, new_col];
                      
        % For every cell that pops out of the lattice on one side, 
        % another cell pops into the lattice on the opposite side
        new_row_col(new_row_col > DIM) = 1;
        new_row_col(new_row_col < 1) = DIM; 
        
        % Convert [row,column] coordinates to linear indices
        new_site = sub2ind([DIM,DIM], new_row_col(1), new_row_col(2));
        
        % If the new site is vacant, the cell moves
        if culture_dish(new_site) == 0
            culture_dish(current_site) = 0;
            culture_dish(new_site) = 1;
            cell_sites(cell_sites == current_site) = new_site;
            which_try_to_move(which_try_to_move == current_site) = new_site;
        end
    end

    %%% PROLIFERATION %%%
    
    % N_t agents are selected with replacement, at random, one at a time
    % and are given a chance to proliferate
    for choice = 1:N_t
        selected_to_prolif = datasample(cell_sites(1:N_t),1);
        if rand <= P_p  % only a portion of those selected will try to
                        % proliferate
                        
            % Convert linear indices to [row, column] coordinates
            parent_site = selected_to_prolif;   
            [parent_row, parent_col] = ind2sub(DIM, parent_site);

            % Each parent agent attempting proliferation chooses a random 
            % direction in which to proliferate (create a daughter cell)
            delta = randsample([randsample([-1,1],1),0],2);
            daughter_row_col = [parent_row, parent_col] + delta;
        
            % For every parent cell that proliferates out of the lattice on 
            % one side, another cell is born into the lattice on the 
            % opposite side
            daughter_row_col(daughter_row_col > DIM) = 1;
            daughter_row_col(daughter_row_col < 1) = DIM; 
        
            % Convert [row,column] coordinates to linear indices
            daughter_site = sub2ind([DIM,DIM], daughter_row_col(1), ...
                daughter_row_col(2));
        
            % If the new site is vacant, the cell proliferates in it
            if culture_dish(daughter_site) == 0
                culture_dish(daughter_site) = 1;
                N_t = N_t + 1;
                cell_sites(N_t) = daughter_site;
                
                % Add the [parent cell #, daughter cell #, generation #] to
                % track lineage
                parent_cell_num = find(cell_sites == parent_site);
                parent_not_gen_1 = (cell_lin(:,2) == parent_cell_num);
                if any(parent_not_gen_1)
                    gen_num = cell_lin(parent_not_gen_1, 3) + 1;
                else
                    gen_num = 2;
                end
                cell_lin(N_t - N_initial,:) = [parent_cell_num, N_t, ...
                    gen_num];
            end
        end        
    end
    
    %%% RECORD %%%
    
    % Save each figure as a frame in the movie illustrating the evolution
    % of the culture dish
    spy(culture_dish);
    %title(sprintf('timestep = %d',t)); % 
    drawnow
    cell_movie(t+1) = getframe;
    
    % Record the cell population at each timestep
    cell_pop(t+1) = N_t;
end

% Play movie
cell_fig.Visible = 'on';
movie(cell_movie, 1, speed);

% Save evolution information into a structure
evolution_info = struct('cell_population', cell_pop(1:t+1), ...
    'cell_lineage', cell_lin(1:(N_t-N_initial),:));
end
