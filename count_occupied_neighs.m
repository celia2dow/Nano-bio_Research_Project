function NUM_OCCUPIED_NEIGHS = count_occupied_neighs(CULTURE_DIM, N_TSTEP, CELL_SITES, CULTURE_DISH)
% COUNT_OCCUPIED_NEIGHS Given the CELL_SITES which are occupied by cells in
% a CULTURE_DIM by CULTURE_DIM, lattice, count the number of nearest neighbours to each
% occupied site that are also occupied.
NUM_OCCUPIED_NEIGHS = 0;
    for cell = 1:N_TSTEP
        [current_row, current_col] = ind2sub(CULTURE_DIM, CELL_SITES(cell));
        
        % Find the neighbours of the current site
        neighs_sub = [current_row, current_col + 1; ...
            current_row, current_col - 1; ...
            current_row + 1, current_col; ...
            current_row - 1, current_col];
        
        % Assume the torus shape of the domain and allow for edge and
        % corner sites to have neighbours on the opposite edges/corners.
        neighs_sub(neighs_sub>CULTURE_DIM)=1;
        neighs_sub(neighs_sub<1)=CULTURE_DIM;
        
        % Convert [row, column] coordinates to linear indices
        neighs_ind = sub2ind([CULTURE_DIM,CULTURE_DIM], neighs_sub(:,1), neighs_sub(:,2));
        
        % Add to the total the number of neighbours that are occupied
        cell_neighs = CULTURE_DISH(neighs_ind);
        NUM_OCCUPIED_NEIGHS = NUM_OCCUPIED_NEIGHS + nnz(cell_neighs);
    end
end
