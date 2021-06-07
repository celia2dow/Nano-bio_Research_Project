% CELLS_SIMULATION_DRIVER is a script that drives the cells_simulation 
% simulation with the specified parameters.
%
%   This is the work of Celia Dowling 8/6/21
%
%   The input argument for cells_simulation.m is a structure PARAMETERS 
%   which has the following fields:
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
%                                   add up, along with 1/rate_diffusivity,
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
%   The output argument for cells_simulation.m is a structure 
%   EVOLUTION_INFO containing the fields:
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

% Prepare
clear;
close all;

PARAMETERS = struct( ...
    'simulation_duration',24, ... (hours)
    'tstep_duration',1/6, ... (hours)
    'initial_num_cells', 1000, ... 
    'prtcls_per_cell', 100, ...
    'cell_diam', 25, ... (micrometers)
    'dim', 100, ... (cell diameters)x
    'EWT_move', 1/6, ... (hours) 
    'EWTs_proliferate',[4,4,4], ... (hours)
    'EWTs_internalise',[5/3,5/3,5/3],... (hours)
    'prob_inherit', 0.7, ...     
    'rate_prtcl_diffusivity', 0.1, ...  (per timestep)
    'max_prtcl', 30, ... 
    'limiting_stage', 0,... 
    'vid_speed', 2, ... (frames per sec)     
    'visual', 0 ...
    );

folder_name = ['SimulationDuration' num2str(PARAMETERS.simulation_duration) ...
    '_tstepDuration' num2str(PARAMETERS.tstep_duration) ...
    '_N0' num2str(PARAMETERS.initial_num_cells) ...
    '_dim' num2str(PARAMETERS.dim) ...
    '_cellDiam' num2str(PARAMETERS.cell_diam) ...
    '_maxPrtcl' num2str(PARAMETERS.max_prtcl) ...
    '_EWTmove' num2str(PARAMETERS.EWT_move) ...
    '_Pinherit' num2str(PARAMETERS.prob_inherit) ...
    '_EWTprolif' num2str(PARAMETERS.EWTs_proliferate(1)) num2str(PARAMETERS.EWTs_proliferate(end)) ...
    '_rateDiffus' num2str(PARAMETERS.rate_prtcl_diffusivity) ...
    '_EWTintern' num2str(PARAMETERS.EWTs_internalise(1)) num2str(PARAMETERS.EWTs_internalise(end)) ...
    ];

if ~exist(folder_name, 'dir')
    mkdir(folder_name)
end

% Seed for efficiency purposes
rng(22)
                
% Run simulation and collect data
EVOLUTION_INFO = cells_simulation(PARAMETERS);

% Print results
fprintf("\nCell population per time step: \n");
%disp(EVOLUTION_INFO.cell_population)
fprintf("\nCell lineage [parent cell, daughter cell, generation number]: \n");
%disp(EVOLUTION_INFO.cell_lineage)
fprintf("\nPhases per cell (row) per time step (column): \n");
%disp(EVOLUTION_INFO.cell_phase_history)
fprintf("\nAverage mumber of particles that hit, interacting and internalised per cell per timestep: \n");
%disp(EVOLUTION_INFO.average_c_o_p)
fprintf("\nFree particles per cell (row) per hour (column): \n")
%disp(EVOLUTION_INFO.cell_c_o_p(:,:,1))
fprintf("\nInteracting particles per cell (row) per hour (column): \n")
%disp(EVOLUTION_INFO.cell_c_o_p(:,:,2))
fprintf("\nInternalised particles per cell (row) per hour (column): \n")
%disp(EVOLUTION_INFO.cell_c_o_p(:,:,3))
fprintf("\nThe number of times the binomial distribution overdraws particles to internalise: \n")
%disp(EVOLUTION_INFO.count_catch)

% Plot all of the interacting Vs internalised particles over time
figure(2)
% Find a limit for the y-axis to fix: the largest cell frequency to be
% plotted
interact_max = max(EVOLUTION_INFO.cell_c_o_p(:,:,2),[],'all');
internal_max = max(EVOLUTION_INFO.cell_c_o_p(:,:,3),[],'all');
x_max = max(interact_max,internal_max);
local_max = zeros(1,PARAMETERS.simulation_duration+1);
for index = 2:PARAMETERS.simulation_duration+1
    [cells_interact,~] = histcounts(EVOLUTION_INFO.cell_c_o_p(:,index,2));
    [cells_internal,~] = histcounts(EVOLUTION_INFO.cell_c_o_p(:,index,3));
    local_max(index) = max([cells_interact(2:end), cells_internal(2:end)]);
end
y_max = max(local_max);
for time_plot = 1:PARAMETERS.simulation_duration+1
    subplot(5,ceil(PARAMETERS.simulation_duration/5),time_plot)
    tstep = (time_plot-1)/PARAMETERS.tstep_duration + 1; % equivalent timestep index
    % FREQUENCY OF CELLS WITH NUMS OF PARTICLES INTERACTING OVER TIME
    histogram(EVOLUTION_INFO.cell_c_o_p(1:EVOLUTION_INFO.cell_population(tstep),time_plot,2),...
        'FaceColor', [0,0,1], 'FaceAlpha', 0.2)
    hold on
    % FREQUENCY OF CELLS WITH NUMS OF PARTICLES INTERNALISED OVER TIME
    histogram(EVOLUTION_INFO.cell_c_o_p(1:EVOLUTION_INFO.cell_population(tstep),time_plot,3),...
        'FaceColor', [1,0,0], 'FaceAlpha', 0.2)
    hold off
    ylim([0,y_max])
    xlim([0,x_max])
    title(['At ' num2str(time_plot-1) ' hour/s'])
    xlabel('Number of particles')
    ylabel('Cell frequency')
    legend('Interacting','Internalised')
end
sgtitle('Frequency of cells with certain numbers of interacting/internalised particles over time')


% Plot the density of cells in different generations
figure(3)
oldest_cell_gen = max(EVOLUTION_INFO.cell_lineage(:,3));
histogram(EVOLUTION_INFO.cell_lineage(:,3),oldest_cell_gen);
title('Population of cell generations')
xlabel('Generation')
ylabel('Frequency')

% Plot proliferation events over time
figure(4)
scatter(0:(1/6):24,EVOLUTION_INFO.cell_population,4,'k','filled')
title('Cell population over time')
xlabel('time $t$ hours', 'Interpreter', 'latex');
ylabel('Cell population');

% Plot associated Vs interacting+internalised cells over time
d1 = EVOLUTION_INFO.average_c_o_p(2,:);     % Average interacting particles per cell over time
d2 = EVOLUTION_INFO.average_c_o_p(3,:);     % Average internalised particles per cell over time
d3 = d1 + d2;                               % Histogram Sum ‘d1’+‘d2’
binrng = 0:(1/6):24;                        % Create Bin Ranges

figure(5)
subplot(1,2,1)
bar(binrng, d3, 'r')
hold on
bar(binrng, d1, 'b')
hold off
legend('Internalised','Interacting')
title('Stacked histogram')
xlabel('time $t$ hours', 'Interpreter', 'latex');
ylabel('Average number of particles per cell');

subplot(1,2,2)
plot(binrng,d2,'r--',binrng,d1,'b--');
hold on
plot(binrng,d3,'Color',[.5 0 .5])
hold off
legend('Internalised','Interacting','Associated (internalised or interacting)')
title('Separate trends')
xlabel('time $t$ hours', 'Interpreter', 'latex');

sgtitle('Average number of particles associated per cell')

save([pwd '\' folder_name '\variables.mat']);
