function cells_matrix_driver()
% cells_matrix_driver drives the cells_matrix simulation with the specified
% parameters.

total_t = 50;   % total number of time steps
N_initial = 30;  % initial number of agents (cells) in the lattice
P_m = 0.9;      % ability of an agent (cell) to move during each timestep
P_p = 0.01;     % " " " ; to proliferate " " "
DIM = 50;       % lattice dimensions (DIM by DIM)
speed = 2;      % speed of movie frame playback

cells_matrix(total_t, N_initial, P_m, P_p, DIM, speed);

end