function cells_simulation_driver()
% CELLS_SIMULATION_DRIVER drives the cells_simulation simulation with the 
% specified parameters.
%
%   This is the work of Celia Dowling 22/3/21

total_t = 100;   % total number of time steps
N_initial = 10; % initial number of agents (cells) in the lattice. At this 
                % point, this program only works for when N_initial * 100 >
                % DIM * DIM
DIM = 10;       % lattice dimensions (DIM by DIM)
max_prtcl = 30; % the maximum number of particles a cell can internalise
P_move = 0.9;   % probability that an agent (cell) moves in 1 timestep
P_inherit = 0.5; % the probability of a daughter cell born into the site of
                % the parent cell inheriting 1 particle from the parent cell
cycle_probs = [0.1, 0.1, 0.1, 0.1, 0.1];   
                % a list of K probabilities of transitioning between phases
                % phases in the cell proliferation cycle in 1 timestep 
                % (e.g. phase 1 to phase 2, phase 2 to phase 3, ..., phase 
                % N-1 to phase N, phase N to phase 1 having proliferated)
                %    note that these probabilities are discrete 
                %    approximates of the exponential waiting time rates
rate_interacts = 0.3;    % the poisson rate at which particles interact with
                        % a given cell
base_prtcl_probs = [0.2,0.2,0.2];  
                % a list of L base probabilities of particles transitioning
                % between stages of the cell-particle interaction model in 
                % 1 timestep
speed = 2;      % speed of movie frame playback (how many frames per sec)

evolution_info = cells_simulation(total_t, N_initial, DIM, max_prtcl, ...
    P_move, P_inherit, cycle_probs, rate_interacts, base_prtcl_probs, speed);

% Print results
fprintf("\nCell population per time step: \n");
disp(evolution_info.cell_population)
fprintf("\nCell lineage [parent cell, daughter cell, generation number]: \n");
disp(evolution_info.cell_lineage)
fprintf("\nNumber of particles that are free, interacting and internalised per time step: \n");
disp(evolution_info.class_of_particles)
end
