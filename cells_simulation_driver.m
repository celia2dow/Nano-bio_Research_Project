function cells_simulation_driver()
% cells_matrix_driver drives the cells_simulation simulation with the specified
% parameters.

total_t = 70;   % total number of time steps
N_initial = 300;  % initial number of agents (cells) in the lattice
P_m = 0.9;      % ability of an agent (cell) to move during each timestep
P_p = 0.01;     % " " " ; to proliferate " " "
DIM = 50;       % lattice dimensions (DIM by DIM)
speed = 3;      % speed of movie frame playback (how many frames per sec)

evolution_info = cells_simulation(total_t, N_initial, P_m, P_p, DIM, speed);
fprintf("\nCell population per time step: \n");
disp(evolution_info.cell_population)
fprintf("\nCell lineage [parent cell, daughter cell]: \n");
disp(evolution_info.cell_lineage)
end
