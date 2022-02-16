function NEW_SITE = choose_adjacent_site(CULTURE_DIM, CURRENT_SITE)
% CHOOSE_ADJACENT_SITE Given the CURRENT_SITE in a CULTURE_DIM by CULTURE_DIM lattice, 
% choose a NEW_SITE that is directly adjacent (up, down, left, right)

% Convert linear indices to [row,column] coordinates 
[current_row, current_col] = ind2sub(CULTURE_DIM, CURRENT_SITE);
        
% Each selected cell chooses a random direction to move in
delta = randsample([randsample([-1,1],1),0],2);
new_row_col = [current_row, current_col] + delta;
            
% For every cell that pops out of the lattice on one side, another cell
% pops into the lattice on the opposite side
new_row_col(new_row_col > CULTURE_DIM) = 1;
new_row_col(new_row_col < 1) = CULTURE_DIM; 
        
% Convert [row, column] coordinates to linear indices
NEW_SITE = sub2ind([CULTURE_DIM,CULTURE_DIM], new_row_col(1), new_row_col(2));
end
