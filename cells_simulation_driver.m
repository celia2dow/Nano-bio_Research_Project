% CELLS_SIMULATION_DRIVER is a script that drives the cells_simulation 
% simulation with the specified parameters.
%
%   This is the work of Celia Dowling 3/6/21
%
%   The input argument for cells_simulation.m is a structure PARAMETERS 
%   which has the following fields:
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
%   The output argument for cells_simulation.m is a structure 
%   EVOLUTION_INFO containing the fields:
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

% Prepare
clear;
close all;

PARAMETERS = struct( ...
    'total_tsteps', 144, ... 
    'N_initial', 1000, ... 
    'dim', 100, ...      
    'siz_cell', 25, ...  
    'max_prtcl', 30, ... 
    'P_move', 0.9, ...   
    'P_inherit', 0.7, ... 
    'cycle_probs',[0.04,0.04,0.04], ...
    'rate_interacts', 0.1, ...  
    'base_prtcl_probs', [0.1,0.1,0.1], ...
    'speed', 2, ...      
    'visual', 0 ...
    );

folder_name = ['Simulation_totalTsteps' num2str(PARAMETERS.total_tsteps) ...
    '_N0' num2str(PARAMETERS.N_initial) ...
    '_dim' num2str(PARAMETERS.dim) ...
    '_sizCell' num2str(PARAMETERS.siz_cell) ...
    '_maxPrtcl' num2str(PARAMETERS.max_prtcl) ...
    '_Pmove' num2str(PARAMETERS.P_move) ...
    '_Pinherit' num2str(PARAMETERS.P_inherit) ...
    '_cycleProbs' num2str(PARAMETERS.cycle_probs(1)) num2str(PARAMETERS.cycle_probs(end)) ...
    '_rateInteracts' num2str(PARAMETERS.rate_interacts) ...
    '_basePrtclProbs' num2str(PARAMETERS.base_prtcl_probs(1)) num2str(PARAMETERS.base_prtcl_probs(end)) ...
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
local_max = zeros(1,25);
for index = 2:25
    [cells_interact,~] = histcounts(EVOLUTION_INFO.cell_c_o_p(:,index,2));
    [cells_internal,~] = histcounts(EVOLUTION_INFO.cell_c_o_p(:,index,3));
    local_max(index) = max([cells_interact(2:end), cells_internal(2:end)]);
end
y_max = max(local_max);
for time_plot = 1:25
    subplot(5,5,time_plot)
    tsteps_per_hour = floor(PARAMETERS.total_tsteps/24); % number of timesteps in an hour
    tstep = (time_plot - 1) * tsteps_per_hour + 1; % equivalent timestep index
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
