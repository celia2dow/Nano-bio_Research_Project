function cells_simulation_driver()
% CELLS_SIMULATION_DRIVER drives the cells_simulation simulation with the 
% specified parameters.
%
%   This is the work of Celia Dowling 22/3/21

total_tsteps = 144;     % total number of time steps
N_initial = 1000; % initial number of agents (cells) in the lattice. At this 
                % point, this program only works for when N_initial * 100 >
                % DIM * DIM
DIM = 100;       % lattice dimensions (DIM by DIM)
siz_cell = 25;  % average diameter of particlular cell-type (micrometers)
max_prtcl = 30; % the maximum number of particles a cell can internalise
P_move = 0.9;   % probability that an agent (cell) moves in 1 timestep
P_inherit = 0.7; % the probability of a daughter cell born into the site of
                % the parent cell inheriting 1 particle from the parent cell
cycle_probs = [0.04,0.04,0.04];
                % a list of K probabilities of transitioning between phases
                % phases in the cell proliferation cycle in 1 timestep 
                % (e.g. phase 1 to phase 2, phase 2 to phase 3, ..., phase 
                % N-1 to phase N, phase N to phase 1 having proliferated)
                %    note that these probabilities are discrete 
                %    approximates of the exponential waiting time rates
rate_interacts = 0.002;   % the binomial rate at which particles interact 
                        % with a given cell, reflecting the dispersion of
                        % the particle type (number of particles per 
                        % timestep)
base_prtcl_probs = [0.002,0.002,0.002]; %0.1,0.1,0.1;0.3,0.3,0.3];  
                % a list of L base probabilities of particles transitioning
                % between stages of the cell-particle interaction model in 
                % 1 timestep - can be 1 probability per transition or can
                % be 1 probability per cell phase per transition - columns
                % inicate interaction stage and rows indicate cell phase
speed = 2;      % speed of movie frame playback (how many frames per sec)
visual = 0;     % the form of visualisation desired which can be 0,1,2
                %     0 = off
                %     1 = slower, more comprehensive simulation
                %     2 = faster, less comprehensive simulation

evolution_info = cells_simulation(total_tsteps, N_initial, DIM, ...
    siz_cell, max_prtcl, P_move, P_inherit, cycle_probs, ...
    rate_interacts, base_prtcl_probs, speed, visual);

% Print results
fprintf("\nCell population per time step: \n");
%disp(evolution_info.cell_population)
fprintf("\nCell lineage [parent cell, daughter cell, generation number]: \n");
disp(evolution_info.cell_lineage)
fprintf("\nPhases per cell (row) per time step (column): \n");
%disp(evolution_info.cell_phase_history)
fprintf("\nNumber of particles that are free, interacting and internalised per timestep: \n");
%disp(evolution_info.class_of_particles)
fprintf("\nAverage mumber of particles that are free, interacting and internalised per timestep: \n");
%disp(evolution_info.average_c_o_p)
fprintf("\nFree particles per cell (row) per hour (column): \n")
%disp(evolution_info.cell_c_o_p(:,:,1))
fprintf("\nInteracting particles per cell (row) per hour (column): \n")
%disp(evolution_info.cell_c_o_p(:,:,2))
fprintf("\nInternalised particles per cell (row) per hour (column): \n")
%disp(evolution_info.cell_c_o_p(:,:,3))
fprintf("\nThe number of times the binomial distribution overdraws particles to internalise: \n")
disp(evolution_info.count_catch)

% Plot all of the interacting particles over time
figure(1)
for time_plot = 1:25
subplot(5,5,time_plot)
tsteps_per_hour = floor(total_tsteps/24); % number of timesteps that 
tstep = (time_plot - 1) * tsteps_per_hour + 1; % equivalent timestep
histogram(evolution_info.cell_c_o_p(1:evolution_info.cell_population(tstep),time_plot,2))
title(['At ' num2str(time_plot-1) ' hour/s'])
xlabel('number interacting')
ylabel('frequency')
end
sgtitle('Interacting particles over time')

% Plot all of the internalised particles over time
figure(2)
for time_plot = 1:25
subplot(5,5,time_plot)
tsteps_per_hour = floor(total_tsteps/24); % number of timesteps that 
tstep = (time_plot - 1) * tsteps_per_hour + 1; % equivalent timestep
histogram(evolution_info.cell_c_o_p(1:evolution_info.cell_population(tstep),time_plot,3))
title(['At ' num2str(time_plot-1) ' hour/s'])
xlabel('number internalised')
ylabel('frequency')
end
sgtitle('Internalised particles over time')

% Plot the density of cells in different generations
figure(3)
oldest_cell_gen = max(evolution_info.cell_lineage(:,3));
histogram(evolution_info.cell_lineage(:,3),oldest_cell_gen);
title('Population of cell generations')
xlabel('gen')
ylabel('frequency')

% Plot the density of particles in different stages of the particle-cell
% model
figure(4)
[~,L] = size(base_prtcl_probs);
num_in_stage = zeros(1,L+1);
for s = 1:L+1
    num_in_stage(s) = sum(evolution_info.cell_prtcl_last(:,s),1);
end
histogram('BinEdges',-0.5:1:L+0.5,'BinCounts',num_in_stage)
end

