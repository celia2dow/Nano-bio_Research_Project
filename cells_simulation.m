function evolution_info = cells_simulation(total_t, N_initial, DIM, ...
    max_prtcl, P_move, P_inherit, cycle_probs, rate_interacts, ...
    base_prtcl_probs, speed)
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
%       P_move      probability that a cell moves in 1 timestep
%       P_inherit   the probability of a daughter cell born into the site
%                   of the parent cell inheriting 1 particle from the
%                   parent cell
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
%           class_of_particles  a record of the number of particles that
%                               have not been internalised by or bound to 
%                               cells (i.e free particles), interacting
%                               particles and internalised particles over
%                               time
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
% Records free particles, interacting particles and internalised particles
% per timestep in system
tally_prtcls = [prtcls_initial zeros(1, total_t); zeros(2, total_t + 1)];

% Randomise the positions of initial cells on a DIM by DIM lattice. From 
% now on, all arrays with the name suffix 'cell_' correspond in indexing 
% with each other. I.e. cell_...(i) corresponds to cell i.
culture_dish = zeros(DIM);
cell_sites = zeros(1, total_sites); % preallocate space
cell_sites(1:N_t) = randsample(1:total_sites, N_t);
culture_dish(cell_sites(1:N_t))=1; % place cells on culture dish
                                               
% Randomise the phases of initial cells
cell_phases = zeros(1, total_sites); % preallocate space
cell_phases(1:N_t) = datasample(1:K, N_t);

% Initialise arrays recording the number of particles in different stages
% (per column) of interaction with each cell (per row). The first column is
% to simulate the number of particles hovering around one cell initially.
cell_prtcls = zeros(total_sites, L+1);

% Prepare movie
petri_fig = figure;
axis tight manual % ensures getframe() returns a consistent size
%ax = gca;
%ax.NextPlot = 'replaceChildren';
petri_movie(total_t + 1) = struct('cdata',[],'colormap',[]);
petri_fig.Visible = 'off';

% Save first movie frame (the initial culture dish)
draw_frame(N_t, DIM, t, L, max_prtcl, cell_sites, cell_phases, cell_prtcls)
petri_movie(t+1) = getframe(petri_fig);

% Loop stops when timesteps are up or when culture dish is full
while t < total_t && ~all(culture_dish, 'all') 
    t = t+1;
    clf
    %%% CELL-PARTICLE INTERACTIONS %%%
    
    % Assuming that particles are constantly evenly dispersing themselves  
    % around the petri dish/matrix, then the number of particles hovering  
    % over a matrix site at the beginning of each timestep (prtcls_hover) 
    % will be the total number of particles that are not internalised or 
    % currently bound to/interacting with cells divided amongst the number 
    % of matrix sites.
    tally_prtcls(:,t+1) = tally_prtcls(:,t);
    num_per_site = tally_prtcls(1,t+1)/total_sites; % This may be a fraction
    cell_prtcls(:,1) = floor(num_per_site); % This may be 0
    if mod(num_per_site, floor(num_per_site)) ~= 0 % Check for remainder
        leftovers = tally_prtcls(1,t+1) - total_sites * floor(num_per_site);
        % Randomly distribute the remaining particles around the matrix and
        % if a cell lies on such a site, add to its hovering particle count
        prtcl_sites = randsample(1:total_sites, leftovers);
        [~, indices] = intersect(cell_sites, prtcl_sites);
        cell_prtcls(indices,1) = cell_prtcls(indices,1) + 1;
    end
    
    % Given that one cell can interact with any of the particles hovering
    % over their matrix site OR any of the particles already bound to/
    % interacting with the cell, the total number of particles available to
    % a cell for interaction is the sum of these two types.
    num_available = sum(cell_prtcls(1:N_t,1:L),2);
    
    % Allow each cell to attempt to interact with x particles where x is
    % drawn from a Poisson distribution with a parameter proportionate to
    % the number of particles available for interaction per cell
    %       lambda = rate_interacts * prtcls_avail
    cell_num_attempts = poissrnd(rate_interacts .* num_available);
    non_zero_attempts = find(cell_num_attempts > 0);
    if non_zero_attempts
        for index = 1:length(non_zero_attempts)
            cell = non_zero_attempts(index);
            % Select x particles from the distribution of particles currently 
            % bound to or hovering over the cell and record their stages.
            %x = cell_num_attempts(cell); %This is sometimes > num_avail
            x = min(num_available(cell),cell_num_attempts(cell));
            prtcls_per_stage = cell_prtcls(cell,:);
            attempts = randsample(1:num_available(cell), x);
            prev_index = 0;
            for stage = 1:L
                num_in_stage = prtcls_per_stage(stage);
                attempts(attempts > prev_index & ...
                    attempts <= prev_index + num_in_stage) = stage; 
                prev_index = prev_index + num_in_stage;
            end

            % Run a Bernoulli trial for each particle and one at a time, 
            % compare to the dynamic transition probability to see whether 
            % or not it succeeds in its attempt.
            num_internalised = prtcls_per_stage(L+1);
            Bernoulli_trials = rand(1,x);
            successes = zeros(1,x);
            for prtcl = 1:x
                stage = attempts(prtcl);
                dynamic_prtcl_prob =  base_prtcl_probs(stage) * ...
                    (1-num_internalised/max_prtcl);
                successes(prtcl) = stage * (Bernoulli_trials(prtcl) ...
                    <= dynamic_prtcl_prob);
                if stage == L && successes(prtcl) ~=0
                    num_internalised = num_internalised + 1;
                end
            end

            % Update the total number of particles in each stage of the
            % cell-particle interaction model
            for stage = 1:L
                success_in_stage = sum(successes == stage);
                if stage == 1
                    % Update total number of non-interacting and 
                    % uninternalised (i.e free) particles in system
                    tally_prtcls(1,t+1) = tally_prtcls(1,t+1) - success_in_stage;
                    % Update total number of particles interacting with cells
                    tally_prtcls(2,t+1) = tally_prtcls(2,t+1) + success_in_stage;
                elseif stage == L
                    % Update total number of particles interacting with
                    % cells in system
                    tally_prtcls(2,t+1) = tally_prtcls(2,t+1) - success_in_stage;
                    % Update total number of internalised particles in
                    % system
                    tally_prtcls(3,t+1) = tally_prtcls(3,t+1) + success_in_stage;
                end
                prtcls_per_stage(stage) = prtcls_per_stage(stage) - ...
                    success_in_stage;
                prtcls_per_stage(stage + 1) = prtcls_per_stage(stage ...
                    + 1) + success_in_stage;
            end
            cell_prtcls(cell,:) = prtcls_per_stage;
        end
    end
    
    %%% MOVEMENT %%%
    
    % N_t cells are selected with replacement, at random, one at a time
    % and are given a chance to move
    for choice = 1:N_t
        % Only a portion of those selected will try to move
        if rand <= P_move
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
                                    (P_inherit^i * (1-P_inherit)^(n_int-i) + ...
                                    P_inherit^(n_int-i) * (1-P_inherit)^i);
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
                        cell_prtcls(N_t, stage) = daught_2;
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
    
    %%% RECORD %%%
    
    % Save each figure as a frame in the movie illustrating the evolution
    % of the culture dish
    draw_frame(N_t, DIM, t, L, max_prtcl, cell_sites, cell_phases, cell_prtcls)
    petri_movie(t+1) = getframe(petri_fig);
    
    % Record the cell population at each timestep
    population(t+1) = N_t;
end

% Print cell particles
%disp(cell_prtcls((1:N_t),:));

% Play movie
petri_fig.Visible = 'on';
movie(petri_fig, petri_movie, 1, speed);

% Save evolution information into a structure
evolution_info = struct('cell_population', population(1:t+1), ...
    'cell_lineage', lineage(1:N_t,:), ...
    'class_of_particles', tally_prtcls(:,1:t+1));
end



%%% SUBFUNCTIONS %%%

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


function draw_frame(N_t, DIM, t, L, max_prtcl, cell_sites, cell_phases, cell_prtcls)
% DRAW_FRAME Save each figure as a frame in the movie illustrating the 
% evolution of the culture dish as cells move, proliferate and internalise
% particles

% Draw a scatter plot of the culture dish with:
%       size of cell indicating phase of cell proliferation cycle
%       opaqueness of cell colour indicating number of internalised particles
%       width of cell outline indicating number of interacting particles
[rows, cols] = ind2sub(DIM, cell_sites);
% Sizes of cell depend on cell phase
siz_cell = 10*cell_phases; 
col_cell = [1 0 0];
% Transparency of cell colour depends on number of internalised particles
transp_cell = cell_prtcls(:,end)/max_prtcl;
for cell = 1:N_t
    % Width of edge of cell depends on number of interacting particles
    if L == 1 || sum(cell_prtcls(cell,2:end - 1)) == 0
        wid_edge = 0.1;
        col_edge = [0 1 0];
    else
        wid_edge = 0.1 * sum(cell_prtcls(cell,2:end - 1),2);
        col_edge = [0 0 1];
    end 
    % Plot point for cell ensuring the graph represents the lattice 
    % positioning
    scatter(cols(cell), DIM - rows(cell) + 1, siz_cell(cell), ...
    'MarkerEdgeColor', col_edge, 'MarkerFaceColor', col_cell, ...
    'MarkerFaceAlpha', transp_cell(cell), 'LineWidth', wid_edge);
    hold on
end
xlim([0.5 DIM+0.5]); ylim([0.5 DIM+0.5]);
title(sprintf('timestep = %d',t));
drawnow
end


function X = rand_pmf(x,pmf,num)
% RAND_PMF Random numbers from a user defined discrete distribution
%
% X=rand_gen(x,pmf,N)
% Input:
% x   : set of the all possible values that the desired random signal can
%       assume
% pmf : vector that cointains the probability of each possible 
%       value of x
% N   : number of random values to be chosen
% output:
% X   : random signal whit the desired pmf
%
% Example: 
% pmf=[1/3 1/3 1/3]
% x=[1 2 3];
% num=100;
% X=rand_gen(x,pmf,num);
a=[0;cumsum((pmf(:)))]*(1./rand(1,num));
b=a>ones(length(pmf)+1,num);
[~,index]=max(b);
x=x(:);
if index>1
    X=x(index-1);
else
    X=0;
end
end
