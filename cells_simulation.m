function EVOLUTION_INFO = cells_simulation(PARAMETERS)
% CELLS_SIMULATION simulates the random proliferation and movement of cells
% in a 2.D cell monolayer within a culture dish.
%
%   The section of the culture dish illustrated is represented by a dim by
%   by dim square matrix/ lattice. The inter-lattice spacing is taken to be
%   the average diameter of a cell (25 micrometers).
%
%   The input argument is a structure PARAMETERS containing the fields:
%
%       total_tsteps    total number of timesteps
%       N_initial   initial number of cells in the lattice
%       dim         lattice dimensions (dim by dim)
%       siz_cell    average diameter of particlular cell-type (micrometers)
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
%       rate_interacts  the binomial rate at which particles interact with 
%                       a given cell, reflecting the dispersion of the 
%                       particle type (number of particles per timestep)
%       base_prtcl_probs    a list of L base probabilities of particles 
%                           transitioning between stages of the 
%                           cell-particle interaction model in 1 timestep -
%                           can be 1 probability per transition or can be 1
%                           probability per cell phase per transition -
%                           columns inicate interaction stage and rows 
%                           indicate cell phase
%       speed       speed of movie frame playback (frames per sec)
%       visual      the form of visualisation desired which can be 0,1,2
%                       0 = off
%                       1 = slower, more comprehensive simulation
%                       2 = faster, less comprehensive simulation
%
%   The output argument is a structure EVOLUTION_INFO containing the fields:
%
%       cell_population     a record of the cell population over time
%       cell_lineage        a record of cell lineage 
%                           [parent cell, daughter cell, generation]
%       cell_phase_history  a record of the cell phases over time
%       class_of_particles  a record of the number of particles that have
%                           not been internalised by or bound to cells (i.e
%                           free particles), interacting particles and
%                           internalised particles over time
%       average_c_o_p       a record of the average class of particles over
%                           time
%       cell_c_o_p          a 3D array recording the class of particles on 
%                           a cellular basis after a particular amount of
%                           time (e.g 1 hour, half an hour)
%       count_catch         the number of times the catch (ensuring no more
%                           particles are internalised than the carrying
%                           capacity allows for) is used.
%
%   This is the work of Celia Dowling 4/6/21

%%% PREPARE %%%

% Assign much-used fields of PARAMETERS to shorter variables
total_tsteps = PARAMETERS.total_tsteps;
dim = PARAMETERS.dim;    
max_prtcl = PARAMETERS.max_prtcl; 

% Initialise variables and constants
t = 0; % the timestep
N_t = PARAMETERS.N_initial; % the number of cells at timestep t
prtcls_initial = N_t * 100; % the number of particles initially
total_sites = dim * dim; % total number of possible positions in petri dish
tsteps_per_hour = floor(total_tsteps/24); % number of timesteps per hour 

% Check how many stages in the cell-particle interactions model there are
% and whether or not the user has defined different probabilities for
% different cell-proliferation phases (matrix input) or not (vector input)
[check,L] = size(PARAMETERS.base_prtcl_probs); 
swtch = 0;
if check > 1
    swtch = 1;
end
% The number of phases in the cell proliferation cycle
K = length(PARAMETERS.cycle_probs); 

% Randomise the positions of initial cells on a dim by dim lattice. From 
% now on, all arrays with the name suffix 'cell_' correspond in indexing 
% with each other. I.e. cell_...(i) corresponds to cell i.
culture_dish = zeros(dim);
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
successes_in_stage = zeros(total_sites, L);

% Initialise arrays/constants for evolution_info fields
population = [N_t zeros(1,total_tsteps)]; % cell_population
lineage = [zeros(N_t,1) (1:N_t)' ones(N_t,1); zeros(total_tsteps,3)]; % cell_lineage
cell_phase_history = [cell_phases' zeros(total_sites, total_tsteps)]; % cell_phase_history
% Record free particles, interacting particles and internalised particles
% per timestep in system (class_of_particles)
tally_prtcls = [prtcls_initial zeros(1, total_tsteps); zeros(2, total_tsteps + 1)];
% Record class of particles on a cell basis (cell_c_o_p) per hour
cell_c_o_p_per_hour = zeros(total_sites,24+1,3);
% The number of times the binomial distribution overdraws particles to
% internalise
count_catch = 0;

if PARAMETERS.visual
    % Prepare movie
    petri_fig = figure;
    axis tight manual % ensures getframe() returns a consistent size
    %ax = gca;
    %ax.NextPlot = 'replaceChildren';
    petri_movie(total_tsteps + 1) = struct('cdata',[],'colormap',[]);
    petri_fig.Visible = 'off';
    
    % Save first movie frame (the initial culture dish)
    draw_frame(N_t, dim, t, L, max_prtcl, cell_sites, cell_phases, ...
        cell_prtcls, PARAMETERS.visual, PARAMETERS.siz_cell)
    petri_movie(t+1) = getframe(petri_fig);
end

% Loop stops when timesteps are up or when culture dish is full of cells
% with their maximum capacity of particles
while t < total_tsteps && ~all(tally_prtcls(:,t+1) == [0; ...
        prtcls_initial-min(total_sites*max_prtcl,prtcls_initial); ...
        min(total_sites*max_prtcl,prtcls_initial)],'all')
    t = t+1;
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
    
    % Given that rate_interacts is a measure of the dispersivity of a 
    % particle type, we will begin by seeing if any hovering particles that 
    % are currently "free" (so in the first column of our cell_prtcls 
    % array) and are not interacting with a cell (not actually in stage 0 
    % in which particles are interacting with but unbound to the cell) will 
    % enter into an interaction with the cell (remain in the population of 
    % cells in stage 0 of the cell-particle internalisation model).
    
    % The way we draw this number - the number of particles we are allowing
    % each cell to attempt to interact with from those hovering on the
    % lattice site - is by drawing from a Binomial distribution where the
    % probability of a successful event = rate_interacts and the number of
    % trials = the number of particles hovering on the lattice site.
    cell_prtcls(1:N_t,1) = binornd(cell_prtcls(1:N_t,1),PARAMETERS.rate_interacts);
   
    for stage = 1:L
        % Now that every column of the cell_prtcls array represents the
        % population of particles in a state of interaction with the
        % indexed cell, we can simply use another binomial distribution,
        % but this time using the dynamic probability of transition, to
        % determine how many in each stage are successful in graduating to
        % the next stage.
        num_internalised = cell_prtcls(1:N_t,L+1);
        % If base probabilities have been defined per cell proliferation 
        % cycle phase.
        if swtch 
            dynamic_prtcl_probs =  PARAMETERS.base_prtcl_probs( ...
                cell_phases(1:N_t),stage) .* (1-num_internalised./max_prtcl);      
        else
            dynamic_prtcl_probs =  PARAMETERS.base_prtcl_probs(stage) .* ...
                        (1-num_internalised./max_prtcl .* ones(N_t,1));
        end
        
        % Due to wanting to keep the total number of particles allowed to
        % transition from a stage in one timestep static in spite of 
        % successful transitions on the same timestep, save the number of 
        % successes without yet updating the cell_prtcl array
        successes_in_stage(1:N_t,stage) = binornd(cell_prtcls(1:N_t, stage), ...
            dynamic_prtcl_probs);
        
        satisfy = (successes_in_stage(1:N_t,stage)+cell_prtcls(1:N_t, stage+1)>max_prtcl);
        if stage == L && any(satisfy)
            successes_in_stage(satisfy,stage)=max_prtcl-cell_prtcls(satisfy, stage+1);
            count_catch = count_catch + 1;
        end
    end
    
    % Update the total/per cell number of particles in each stage of
    % the cell-particle interaction model
    cell_prtcls(1:N_t,1:L) = cell_prtcls(1:N_t,1:L) - successes_in_stage(1:N_t,1:L);
    cell_prtcls(1:N_t,2:L+1) = cell_prtcls(1:N_t,2:L+1) + successes_in_stage(1:N_t,1:L);
    
    % Update the total class of particles
    tally_prtcls(2,t+1) = sum(cell_prtcls(1:N_t,2:L),'all'); % interacting
    tally_prtcls(3,t+1) = sum(cell_prtcls(1:N_t,L+1)); % internalised
    tally_prtcls(1,t+1) = prtcls_initial - sum(tally_prtcls(2:3,t+1)); % free
    
    %%% MOVEMENT %%%
    
    % N_t cells are selected with replacement, at random, one at a time
    % and are given a chance to move
    for choice = 1:N_t
        % Only a portion of those selected will try to move
        if rand <= PARAMETERS.P_move
            % Choose a random cell
            current_site = datasample(cell_sites(1:N_t),1);
            
            % Call local function to choose a new site to move to
            new_site = choose_adjacent_site(dim, current_site);
        
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
        if rand <= PARAMETERS.cycle_probs(old_phase)
            % For cells in their final proliferation cycle phase
            if old_phase == K
                % Call local function to choose daughter_site in which to 
                % attempt to proliferate into (create a daughter cell)
                daughter_site = choose_adjacent_site(dim, parent_site);
        
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
                                    (PARAMETERS.P_inherit^i * ...
                                    (1-PARAMETERS.P_inherit)^(n_int-i) + ...
                                    PARAMETERS.P_inherit^(n_int-i) * ...
                                    (1-PARAMETERS.P_inherit)^i);
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
    
    if PARAMETERS.visual
        clf
        % Save each figure as a frame in the movie illustrating the evolution
        % of the culture dish
        draw_frame(N_t, dim, t, L, max_prtcl, cell_sites, cell_phases, ...
            cell_prtcls, PARAMETERS.visual, PARAMETERS.siz_cell)
        petri_movie(t+1) = getframe(petri_fig);
    end
    
    % Record the class of particles on a cell basis every hour
    if mod(t, tsteps_per_hour) == 0
        cell_c_o_p_per_hour(1:N_t,t/tsteps_per_hour+1,1) = ...
            cell_prtcls(1:N_t,1); % Cells that hit
        cell_c_o_p_per_hour(1:N_t,t/tsteps_per_hour+1,2) = ...
            sum(cell_prtcls(1:N_t,2:L),2); % Cells that interact
        cell_c_o_p_per_hour(1:N_t,t/tsteps_per_hour+1,3) = ...
            cell_prtcls(1:N_t,L+1); % Cells that are internalised
    end
    
    % Record the cell phases at each timestep
    cell_phase_history(1:N_t,t+1) = cell_phases(1:N_t)';
    
    % Record the cell population at each timestep
    population(t+1) = N_t;
    
end

% Print cell particles
%disp(cell_prtcls((1:N_t),:));

if PARAMETERS.visual
    % Play movie
    petri_fig.Visible = 'on';
    movie(petri_fig, petri_movie, 1, PARAMETERS.speed);
end

% Calculate average class of particles per timestep (average_c_o_p)
average_prtcls = zeros(3,t+1);
for row = 1:3
    average_prtcls(row,:) = tally_prtcls(row,1:t+1)./population(1:t+1);
end

% Save evolution information into a structure
EVOLUTION_INFO = struct('cell_population', population(1:t+1), ...
    'cell_lineage', lineage(1:N_t,:), ...
    'cell_phase_history', cell_phase_history(1:N_t,1:t+1),...
    'class_of_particles', tally_prtcls(:,1:t+1), ...
    'average_c_o_p', average_prtcls(:,1:t+1), ...
    'cell_c_o_p', cell_c_o_p_per_hour(1:N_t,:,:), ...
    'count_catch', count_catch, ...
    'cell_prtcl_last', cell_prtcls(1:N_t,:));
end



%%% SUBFUNCTIONS %%%

function NEW_SITE = choose_adjacent_site(DIM, CURRENT_SITE)
% CHOOSE_ADJACENT_SITE Given the CURRENT_SITE in a DIM by DIM lattice, 
% choose a NEW_SITE that is directly adjacent (up, down, left, right)

% Convert linear indices to [row,column] coordinates 
[current_row, current_col] = ind2sub(DIM, CURRENT_SITE);
        
% Each selected cell chooses a random direction to move in
delta = randsample([randsample([-1,1],1),0],2);
new_row_col = [current_row, current_col] + delta;
            
% For every cell that pops out of the lattice on one side, another cell
% pops into the lattice on the opposite side
new_row_col(new_row_col > DIM) = 1;
new_row_col(new_row_col < 1) = DIM; 
        
% Convert [row, column] coordinates to linear indices
NEW_SITE = sub2ind([DIM,DIM], new_row_col(1), new_row_col(2));
end


function draw_frame(N_T, DIM, T, L, MAX_PRTCL, CELL_SITES, CELL_PHASES, ...
    CELL_PRTCLS, VISUAL, SIZ_CELL)
% DRAW_FRAME Save each figure as a frame in the movie illustrating the 
% evolution of the culture dish as cells move, proliferate and internalise
% particles.

[rows, cols] = ind2sub(DIM, CELL_SITES);
rows_um = rows * SIZ_CELL;
cols_um = cols * SIZ_CELL;
dim_um = DIM * SIZ_CELL;
siz_phase = 10*CELL_PHASES; % Sizes of cell depend on cell phase

if VISUAL == 1 
    %%% Slower, more comprehensive visual %%%
    
    % Draw a scatter plot of the culture dish with:
    %       size of cell indicating phase of cell proliferation cycle
    %       opacity of cell colour indicating number of internalised particles
    %       width of cell outline indicating number of interacting particles
    col_cell = [1 0 0];
    % Transparency of cell colour depends on number of internalised particles
    transp_cell = CELL_PRTCLS(:,end)/MAX_PRTCL;
    for cell = 1:N_T
        % Width of edge of cell depends on number of interacting particles
        if L == 1 || sum(CELL_PRTCLS(cell,2:end - 1)) == 0
            wid_edge = 0.1;
            col_edge = [0 1 0];
        else
            wid_edge = 0.1 * sum(CELL_PRTCLS(cell,2:end - 1),2);
            col_edge = [0 0 1];
        end 
        % Plot point for cell ensuring the graph represents the lattice 
        % positioning
        scatter(cols_um(cell), dim_um - rows_um(cell) + 1, siz_phase(cell), ...
        'MarkerEdgeColor', col_edge, 'MarkerFaceColor', col_cell, ...
        'MarkerFaceAlpha', transp_cell(cell), 'LineWidth', wid_edge);
        hold on
    end
else
    %%% Faster, less comprehensive visual %%%
    
    % Draw a scatter plot of the culture dish with:
    %       size of cell indicating phase of cell proliferation cycle
    %       pink opacity of cell colour indicating number of internalised particles
    %       acqua opacity of cell colour indicating number of interacting particles
    c = [1-sum(CELL_PRTCLS(1:N_T,2:end - 1),2)./100 ... % interacting particles
        1-(CELL_PRTCLS(1:N_T,end)./MAX_PRTCL) ...% particles internalised
        ones(N_T,1)]; % to ensure that no particles gives white
    scatter(cols_um(1:N_T), dim_um - rows_um(1:N_T) + 1, ...
        siz_phase(1:N_T), c, 'filled', 'MarkerEdgeColor', 'k'); 
end

xlim([0.5*SIZ_CELL dim_um+0.5*SIZ_CELL]); 
ylim([0.5*SIZ_CELL dim_um+0.5*SIZ_CELL]);
xlabel('$x (\mu m)$', 'Interpreter', 'latex');
ylabel('$y (\mu m)$', 'Interpreter', 'latex');
title(sprintf('timestep = %d',T));
drawnow
end


function X = rand_pmf(x,PMF,NUM)
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
% X   : random signal with the desired pmf
%
% Example: 
% pmf=[1/3 1/3 1/3]
% x=[1 2 3];
% num=100;
% X=rand_gen(x,pmf,num);
a=[0;cumsum((PMF(:)))]*(1./rand(1,NUM));
b=a>ones(length(PMF)+1,NUM);
[~,index]=max(b);
x=x(:);
if index>1
    X=x(index-1);
else
    X=0;
end
end
