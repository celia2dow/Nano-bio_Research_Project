function draw_frame(PARAMETERS,N_TSTEP, TSTEP, L, CELL_SITES, CELL_PHASES, CELL_PRTCLS)
% DRAW_FRAME Save each figure as a frame in the movie illustrating the 
% evolution of the culture dish as cells move, proliferate and internalise
% particles.

[rows, cols] = ind2sub(PARAMETERS.CULTURE_DIM, CELL_SITES);
rows_um = rows * PARAMETERS.cell_diam; % (micrometers)
cols_um = cols * PARAMETERS.cell_diam; % (micrometers)
culture_dim_um = PARAMETERS.culture_dim * PARAMETERS.cell_diam; % (micrometers)
siz_phase = 10*CELL_PHASES; % Sizes of cell depend on cell phase

if PARAMETERS.visual == 1
    %%% Slower, more comprehensive visual %%%
    
    % Draw a scatter plot of the culture dish with:
    %       size of cell indicating phase of cell proliferation cycle
    %       opacity of cell colour indicating number of internalised particles
    %       width of cell outline indicating number of interacting particles
    col_cell = [1 0 0];
    % Transparency of cell colour depends on number of internalised particles
    if PARAMETERS.max_prtcls(L) == inf
        transp_cell = CELL_PRTCLS(:,end) .* 0.015;
    else
        transp_cell = CELL_PRTCLS(:,end)/PARAMETERS.max_prtcls(L);
    end
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
        scatter(cols_um(cell), culture_dim_um - rows_um(cell) + 1, siz_phase(cell), ...
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
    if any(PARAMETERS.max_prtcls == inf)
        max_interact = PARAMETERS.prtcls_per_cell .* 0.5;
        c = [1-sum(CELL_PRTCLS(1:N_TSTEP,2:end - 1),2)./max_interact, ... % interacting particles
            (1-CELL_PRTCLS(1:N_TSTEP,end)) .* 0.015, ...% particles internalised
            ones(N_TSTEP,1)]; % to ensure that no particles gives white
    else
        max_interact = sum(PARAMETERS.max_prtcls(1:L-1),2);
        c = [1-sum(CELL_PRTCLS(1:N_TSTEP,2:end - 1),2)./max_interact, ... % interacting particles
            1-(CELL_PRTCLS(1:N_TSTEP,end)./PARAMETERS.max_prtcls(L)), ...% particles internalised
            ones(N_TSTEP,1)]; % to ensure that no particles gives white
    end
    scatter(cols_um(1:N_TSTEP), culture_dim_um - rows_um(1:N_TSTEP) + 1, ...
        siz_phase(1:N_TSTEP), c, 'filled', 'MarkerEdgeColor', 'k'); 
end

xlim([0.5*PARAMETERS.cell_diam culture_dim_um+0.5*PARAMETERS.cell_diam]); 
ylim([0.5*PARAMETERS.cell_diam culture_dim_um+0.5*PARAMETERS.cell_diam]);
xlabel('$x (\mu m)$', 'Interpreter', 'latex');
ylabel('$y (\mu m)$', 'Interpreter', 'latex');
title(sprintf('time = %3.1f hours',TSTEP*PARAMETERS.tstep_duration));
drawnow
end
