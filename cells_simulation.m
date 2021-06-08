function EVOLUTION_INFO = cells_simulation(PARAMETERS)
% CELLS_SIMULATION simulates the random proliferation and movement of cells
% in a 2.D cell monolayer within a culture dish as they interact with and
% internalised particles.
%
%   This is the work of Celia Dowling 8/6/21
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
%       dim (cell diameters)        Lattice dimensions are dim by dim.
%       EWT_move (hours)            Expected waiting time for one cell to 
%                                   move 1 cell diameter. Must begreater than
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
%       EWFs_internalise (hours)    List of L expected waiting times for 1
%                                   particle to transition out of each
%                                   stage in the cell-particle interaction 
%                                   model. Can be 1 EWT per stage (e.g. 
%                                   mean time spent free or in stage 0 
%                                   once having hit a cell, mean time spent 
%                                   in stage 1, ..., mean time spent in stage 
%                                   L-1 before transitioning into stage L i.e. 
%                                   being internalised) or can be 1 EWT per
%                                   stage per cell phase; columns indicating
%                                   interaction stage and rows indicating
%                                   cell phase. Elements must be greater than
%                                   or equal to tstep_duration. They should
%                                   add up, along with 1/rate_prtcl_diffusivity,
%                                   to give the expected waiting time for
%                                   internalisation of 1 particle.
%       prob_inherit                The probability of a daughter cell born 
%                                   into the site of the parent cell 
%                                   inheriting 1 particle from the parent
%                                   cell.
%       rate_prtcl_diffusivity (per tstep) The rate at which particles interact  
%                                   with a given cell (fraction of free particles
%                                   per timestep), reflecting the diffusivity 
%                                   of the particle type.
%       max_prtcl                   Maximum number of particles that a cell 
%                                   can have in the specified cell-particle 
%                                   interaction stage (i.e. carrying
%                                   capacity).
%       limiting_stage              The specified cell-particle interaction
%                                   stage that limits the rate of particle
%                                   internalisation. Must be between 1 and
%                                   L for carrying capacity to be applied.
%                                   Setting to 0 will turn off carrying
%                                   capacity.
%       vid_speed (frames per sec)  Speed of movie frame playback.
%       visual                      The form of visualisation desired which 
%                                   can be 0,1 or 2:
%                                       0 = off
%                                       1 = slower, more comprehensive vid
%                                       2 = faster, less comprehensive vid
%
%   The output argument is a structure EVOLUTION_INFO containing the fields:
%
%       cell_population     A record of the cell population over time.
%       cell_lineage        A record of cell lineage: 
%                           [parent cell, daughter cell, generation]
%       cell_phase_history  A record of the cell phases over time.
%       average_c_o_p       A record of the average class of particles over 
%                           time; the number that have not been internalised 
%                           by or bound to cells (i.e free particles),
%                           the number of interacting particles and that of
%                           internalised particles.
%       cell_c_o_p          A 3D array recording the class of particles
%                           (free AND hit, interacting, internalised) on
%                           a cellular basis after a particular amount of
%                           time (e.g 1 hour, half an hour).
%       count_catch         The number of times the catch (ensuring no more
%                           particles are internalised than the carrying
%                           capacity allows for) is used.



%%% CONVERT EXPECTED WAITING TIMES INTO RATES PER TSTEP %%%

% Rate at which a cell moves.
rate_move = PARAMETERS.tstep_duration/PARAMETERS.EWT_move;

% A list of K rates of cells transitioning between phases in the cell 
% proliferation cycle (e.g. phase 1 to phase 2, phase 2 to phase 3, ..., 
% phase K-1 to phase K, phase K to phase 1 having proliferated).
rates_proliferate = PARAMETERS.tstep_duration./PARAMETERS.EWTs_proliferate;
% The number of phases in the cell proliferation cycle.
K = length(rates_proliferate); 

% A list of L rates of particles transitioning between stages of the cell-
% particle interaction model unaffected by carrying capacity (e.g. free to 
% stage 1, stage 1 to stage 2, ..., stage L-1 to stage L or internalised).
% Can be 1 rate per transition or can be 1 rate per cell phase per
% transition - columns inicate interaction stage and rows indicate cell
% phase.
rates_internalise = PARAMETERS.tstep_duration./PARAMETERS.EWTs_internalise;
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
prtcls_initial = N_tstep * PARAMETERS.prtcls_per_cell; % the number of particles initially
total_sites = PARAMETERS.dim ^2; % total number of possible positions in petri dish
% Number of timesteps of duration tstep_duration that fit into
% simulation_duration number of hours.
total_tsteps = floor(PARAMETERS.simulation_duration/PARAMETERS.tstep_duration);
tsteps_per_hour = 1/PARAMETERS.tstep_duration; % number of timesteps per hour 
tstep = 0; % the current timestep

% Randomise the positions of initial cells on a dim by dim lattice. From 
% now on, all arrays with the name suffix 'cell_' correspond in indexing 
% with each other. I.e. cell_...(i) corresponds to cell i.
culture_dish = zeros(PARAMETERS.dim);
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
cell_lineage = [zeros(N_tstep,1) (1:N_tstep)' ones(N_tstep,1); zeros(total_tsteps,3)];
cell_phase_history = [cell_phases' zeros(total_sites, total_tsteps)]; 
% Record free particles, interacting particles and internalised particles
% per timestep in system.
tally_prtcls = [prtcls_initial zeros(1, total_tsteps); zeros(2, total_tsteps + 1)];
% Record class of particles on a cell basis per hour.
cell_c_o_p = zeros(total_sites,24+1,3);
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
    petri_movie(tstep+1) = getframe(petri_fig);
end


% Loop stops when timesteps are up or when culture dish is full of cells
% with their maximum capacity of particles.
while tstep < total_tsteps && ~all(tally_prtcls(:,tstep+1) == [0; ...
        prtcls_initial-min(total_sites*PARAMETERS.max_prtcl,prtcls_initial); ...
        min(total_sites*PARAMETERS.max_prtcl,prtcls_initial)],'all')
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
    % interact with it (i.e. that are available in stage 0 of the cell-
    % particle interaction model) from a Binomial distribution where the
    % probability of a successful event = rate_prtcl_diffusivity and the 
    % number of trials = the number of particles hovering on the lattice 
    % site (i.e. currently in the first column of our cell_prtcls array).
    cell_prtcls(1:N_tstep,1) = binornd(cell_prtcls(1:N_tstep,1),...
        PARAMETERS.rate_prtcl_diffusivity);
   
    for stage = 1:L
        % Now that every column of the cell_prtcls array represents the
        % population of particles in a stage of interaction with the
        % indexed cell (stage 0 to stage L), we can simply use another 
        % binomial distribution to find how many transition to the next stage, 
        % but this time using the dynamic rate of transition dependent on
        % the number of particles in the limiting, stage to determine how many
        % in each stage are successful in graduating to the next stage.
        current_carrying_total = cell_prtcls(1:N_tstep,PARAMETERS.limiting_stage+1);
        % Check if internalisation rates have been defined per cell 
        % proliferation cycle phase or it they're simply uniform over all 
        % phases. Only dynamically affect the rates of transitions occuring 
        % up to the limiting stage.
        if stage < PARAMETERS.limiting_stage+1 % Affected by carrying capacity
            if swtch % Dependent on phase
                dynamic_prtcl_rates =  rates_internalise(cell_phases(1:N_tstep),...
                    stage) .* (1-current_carrying_total./PARAMETERS.max_prtcl);      
            else % Not dependent on phase
                dynamic_prtcl_rates =  rates_internalise(stage) .* (1-...
                    current_carrying_total./PARAMETERS.max_prtcl .* ones(N_tstep,1));
            end
        else % Unaffected by carrying capacity
            if swtch % Dependent on phase
                dynamic_prtcl_rates =  rates_internalise(cell_phases(1:N_tstep),stage);      
            else % Not dependent on phase
                dynamic_prtcl_rates =  rates_internalise(stage);
            end
        end
        % Due to wanting to keep the total number of particles allowed to
        % transition from a stage in one timestep static in spite of 
        % successful transitions on the same timestep, save the number of 
        % successes without yet updating the cell_prtcl array
        successes_in_stage(1:N_tstep,stage) = binornd(cell_prtcls(1:N_tstep,...
            stage), dynamic_prtcl_rates);
        satisfy = (successes_in_stage(1:N_tstep,stage)+cell_prtcls(1:N_tstep,...
            stage+1)>PARAMETERS.max_prtcl);
        if stage == PARAMETERS.limiting_stage-1 && any(satisfy)
            successes_in_stage(satisfy,stage)=PARAMETERS.max_prtcl-cell_prtcls(satisfy, stage+1);
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
            new_site = choose_adjacent_site(PARAMETERS.dim, current_site);
        
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
                daughter_site = choose_adjacent_site(PARAMETERS.dim, parent_site);
        
                % If the new site is vacant, the cell proliferates into it
                % and returns to the first phase of the cell proliferation
                % cycle as well
                if culture_dish(daughter_site) == 0
                    culture_dish(daughter_site) = 1;
                    N_tstep = N_tstep + 1;
                    cell_sites(N_tstep) = daughter_site;
                    cell_phases(N_tstep) = 1;
                    cell_phases(cell_sites == parent_site) = 1;
                    
                    % Add the [parent cell #, daughter cell #, generation #]
                    % to track lineage
                    parent_cell_num = find(cell_sites == parent_site);
                    gen_num = cell_lineage(cell_lineage(:,2) == parent_cell_num, 3) + 1;
                    cell_lineage(N_tstep,:) = [parent_cell_num, N_tstep, gen_num];
                    
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
    
    
    %%% RECORD %%%
    
    if PARAMETERS.visual
        clf
        % Save each figure as a frame in the movie illustrating the evolution
        % of the culture dish
        draw_frame(PARAMETERS,N_tstep, tstep, L, cell_sites, cell_phases, cell_prtcls)
        petri_movie(tstep+1) = getframe(petri_fig);
    end
    
    % Record the class of particles on a cell basis every hour
    if mod(tstep, tsteps_per_hour) == 0
        cell_c_o_p(1:N_tstep,tstep/tsteps_per_hour+1,1) = ...
            cell_prtcls(1:N_tstep,1); % Particles that are free AND hit
        cell_c_o_p(1:N_tstep,tstep/tsteps_per_hour+1,2) = ...
            sum(cell_prtcls(1:N_tstep,2:L),2); % Particles that interact
        cell_c_o_p(1:N_tstep,tstep/tsteps_per_hour+1,3) = ...
            cell_prtcls(1:N_tstep,L+1); % Particles that are internalised
    end
    
    % Record the cell phases at each timestep
    cell_phase_history(1:N_tstep,tstep+1) = cell_phases(1:N_tstep)';
    
    % Record the cell population at each timestep
    cell_population(tstep+1) = N_tstep;
    
end

% Print cell particles
%disp(cell_prtcls((1:N_t),:));

if PARAMETERS.visual
    % Play movie
    petri_fig.Visible = 'on';
    movie(petri_fig, petri_movie, 1, PARAMETERS.vid_speed);
end

% Calculate average class of particles per timestep (average_c_o_p)
average_c_o_p = zeros(3,tstep+1);
for row = 1:3
    average_c_o_p(row,:) = tally_prtcls(row,1:tstep+1)./cell_population(1:tstep+1);
end

% Save evolution information into a structure
EVOLUTION_INFO = struct('cell_population', cell_population(1:tstep+1), ...
    'cell_lineage', cell_lineage(1:N_tstep,:), ...
    'cell_phase_history', cell_phase_history(1:N_tstep,1:tstep+1),...
    'average_c_o_p', average_c_o_p(:,1:tstep+1), ...
    'cell_c_o_p', cell_c_o_p(1:N_tstep,:,:), ...
    'count_catch', count_catch);
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


function draw_frame(PARAMETERS,N_TSTEP, TSTEP, L, CELL_SITES, CELL_PHASES, CELL_PRTCLS)
% DRAW_FRAME Save each figure as a frame in the movie illustrating the 
% evolution of the culture dish as cells move, proliferate and internalise
% particles.

[rows, cols] = ind2sub(PARAMETERS.dim, CELL_SITES);
rows_um = rows * PARAMETERS.cell_diam; % (micrometers)
cols_um = cols * PARAMETERS.cell_diam; % (micrometers)
dim_um = PARAMETERS.dim * PARAMETERS.cell_diam; % (micrometers)
siz_phase = 10*CELL_PHASES; % Sizes of cell depend on cell phase

if PARAMETERS.visual == 1 
    %%% Slower, more comprehensive visual %%%
    
    % Draw a scatter plot of the culture dish with:
    %       size of cell indicating phase of cell proliferation cycle
    %       opacity of cell colour indicating number of internalised particles
    %       width of cell outline indicating number of interacting particles
    col_cell = [1 0 0];
    % Transparency of cell colour depends on number of internalised particles
    transp_cell = CELL_PRTCLS(:,end)/PARAMETERS.max_prtcl;
    for cell = 1:N_TSTEP
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
    c = [1-sum(CELL_PRTCLS(1:N_TSTEP,2:end - 1),2)./PARAMETERS.prtcls_per_cell ... % interacting particles
        1-(CELL_PRTCLS(1:N_TSTEP,end)./PARAMETERS.prtcls_per_cell) ...% particles internalised
        ones(N_TSTEP,1)]; % to ensure that no particles gives white
    scatter(cols_um(1:N_TSTEP), dim_um - rows_um(1:N_TSTEP) + 1, ...
        siz_phase(1:N_TSTEP), c, 'filled', 'MarkerEdgeColor', 'k'); 
end

xlim([0.5*PARAMETERS.cell_diam dim_um+0.5*PARAMETERS.cell_diam]); 
ylim([0.5*PARAMETERS.cell_diam dim_um+0.5*PARAMETERS.cell_diam]);
xlabel('$x (\mu m)$', 'Interpreter', 'latex');
ylabel('$y (\mu m)$', 'Interpreter', 'latex');
title(sprintf('time = %3.1f hours',TSTEP*PARAMETERS.tstep_duration));
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
