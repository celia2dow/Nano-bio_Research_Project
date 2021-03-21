function cells_simulation_driver()
% cells_simulation_driver drives the cells_simulation simulation with the 
% specified parameters.

total_t = 222;   % total number of time steps
N_initial = 30;  % initial number of agents (cells) in the lattice
P_m = 0.9;      % ability of an agent (cell) to move during each timestep
cycle_probs = [0.1, 0.1, 0.1, 0.1, 0.1];   
                % a list of K probabilities of transitioning between phases
                % in the cell cycle (e.g. phase 1 to phase 2, phase 2 to
                % phase 3, ..., phase N-1 to phase N, phase N to phase 1
                % having proliferated)
DIM = 50;       % lattice dimensions (DIM by DIM)
speed = 3;      % speed of movie frame playback (how many frames per sec)

evolution_info = cells_simulation(total_t, N_initial, P_m, cycle_probs, ...
    DIM, speed);
fprintf("\nCell population per time step: \n");
disp(evolution_info.cell_population)
fprintf("\nCell lineage [parent cell, daughter cell, generation number]: \n");
disp(evolution_info.cell_lineage)
end
