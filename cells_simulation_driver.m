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
%                                   add up, along with 1/rate_prtcl_diffusivity,
%                                   to give the expected waiting time for
%                                   internalisation of 1 particle.
%       max_prtcls                  The particle capacities of each of the
%                                   L cell-particle interaction model
%                                   stages (i.e., capacity of stage 1, ...,
%                                   capacity of stage L). If a stage has 
%                                   unlimited capacity, set it to "inf". If
%                                   all stages have been set to "inf", 
%                                   carrying capacity is effectively 
%                                   switched off.
%       prob_inherit                The probability of a daughter cell born 
%                                   into the site of the parent cell 
%                                   inheriting 1 particle from the parent
%                                   cell.
%       rate_prtcl_diffusivity (per tstep) The rate at which particles interact  
%                                   with a given cell (fraction of free particles
%                                   per timestep), reflecting the diffusivity 
%                                   of the particle type.
%       vid_speed (frames per sec)  Speed of movie frame playback.
%       visual                      The form of visualisation desired which 
%                                   can be 0,1 or 2:
%                                       0 = off
%                                       1 = slower, more comprehensive vid
%                                           + gif
%                                       2 = faster, less comprehensive vid
%                                           + gif
%       folder_path                 The driver function generates this
%                                   depending on the date and the other
%                                   parameters set. It will save all output
%                                   into this path.
%
%   The output argument for cells_simulation.m is a structure 
%   EVOLUTION_INFO containing the fields:
%
%       cell_population     A record of the cell population over time.
%       cell_lineage_history A record of cell lineage: 
%                           [parent cell, daughter cell, generation at 0 hours,
%                           generation at 6 hours, generation at 12 hours, etc.]
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

% Parameters relating to the fluorescent intensity dot plot
stain1_part = 4000; % Without intensity transformations = 1
stain1_std_dev = 0; % Without intensity transformations = 0
stain1_background = 100; % Without intensity transformations = 0
stain2_part = 10000; % Without intensity transformations = 1
stain2_std_dev = 0; % Without intensity transformations = 0
stain2_background = 200; % Without intensity transformations = 0
alpha = 0.9; % Without intensity transformations = 1

PARAMETERS = struct( ...
    'simulation_duration',24, ... (hours)
    'tstep_duration',1/6, ... (hours)
    'initial_num_cells', 1000, ... 
    'prtcls_per_cell', 100, ...
    'cell_diam', 25, ... (micrometers)
    'dim', 100, ... (cell diameters)x
    'EWT_move', 1/6, ... (hours) 
    'EWTs_proliferate',[4,4,4], ... [phase 1, ..., phase K](hours) 
    'EWTs_internalise',[5/3,5/3,5/3],... [stage 0, ..., stage L-1](hours)
    'max_prtcls', [inf,inf,30], ... [stage 1, ..., stage L]
    'prob_inherit', 0.7, ...     
    'rate_prtcl_diffusivity', 0.1, ...  (per timestep)
    'vid_speed', 2, ... (frames per sec)     
    'visual', 0 ...
    );

folder_name = ['Simulation_Duration' num2str(PARAMETERS.simulation_duration) ...
    '_tstepDuration' num2str(PARAMETERS.tstep_duration) ...
    '_N0' num2str(PARAMETERS.initial_num_cells) ...
    '_dim' num2str(PARAMETERS.dim) ...
    '_cellDiam' num2str(PARAMETERS.cell_diam) ...
    '_maxPrtcls' num2str(PARAMETERS.max_prtcls) ...
    '_EWTmove' num2str(PARAMETERS.EWT_move) ...
    '_Pinherit' num2str(PARAMETERS.prob_inherit) ...
    '_EWTprolif' num2str(PARAMETERS.EWTs_proliferate) ...
    '_rateDiffus' num2str(PARAMETERS.rate_prtcl_diffusivity) ...
    '_EWTintern' num2str(PARAMETERS.EWTs_internalise) ...
    '_2x_fluoroPrtclStddevBckgrnd' num2str(stain1_part) num2str(stain1_std_dev) num2str(stain1_background) ...
    num2str(stain2_part) num2str(stain2_std_dev) num2str(stain2_background) ...
    '_fluoroAlpha' num2str(alpha) ...
    ];

PARAMETERS.folder_path = [pwd '\' date '\' folder_name];

if ~exist(PARAMETERS.folder_path, 'dir')
    mkdir(PARAMETERS.folder_path)
end

% Seed for efficiency purposes
%rng(22)
                
% Run simulation and collect data
EVOLUTION_INFO = cells_simulation(PARAMETERS);

% Print results
%{
fprintf("\nCell population per time step: \n");
disp(EVOLUTION_INFO.cell_population)
fprintf("\nCell lineage [parent cell, daughter cell, generation number]: \n");
disp(EVOLUTION_INFO.cell_lineage_history)
fprintf("\nPhases per cell (row) per time step (column): \n");
disp(EVOLUTION_INFO.cell_phase_history)
fprintf("\nAverage mumber of particles that hit, interacting and internalised per cell per timestep: \n");
disp(EVOLUTION_INFO.average_c_o_p)
fprintf("\nFree particles per cell (row) per hour (column): \n")
disp(EVOLUTION_INFO.cell_c_o_p(:,:,1))
fprintf("\nInteracting particles per cell (row) per hour (column): \n")
disp(EVOLUTION_INFO.cell_c_o_p(:,:,2))
fprintf("\nInternalised particles per cell (row) per hour (column): \n")
disp(EVOLUTION_INFO.cell_c_o_p(:,:,3))
%}
fprintf("\nThe number of times the binomial distribution overdraws particles to internalise: \n")
disp(EVOLUTION_INFO.count_catch)

% Find limits for the x and y-axis: the largest values to be plotted
interact_max = max(EVOLUTION_INFO.cell_c_o_p(:,:,2),[],'all');
internal_max = max(EVOLUTION_INFO.cell_c_o_p(:,:,3),[],'all');
x_max = max(interact_max,internal_max);
local_max = zeros(1,floor(PARAMETERS.simulation_duration)+1);
for index = 2:floor(PARAMETERS.simulation_duration)+1
    [cells_interact,~] = histcounts(EVOLUTION_INFO.cell_c_o_p(:,index,2));
    [cells_internal,~] = histcounts(EVOLUTION_INFO.cell_c_o_p(:,index,3));
    local_max(index) = max([cells_interact(2:end), cells_internal(2:end)]);
end
y_max = max(local_max);
S1_max = stain1_part * (PARAMETERS.prtcls_per_cell + alpha*PARAMETERS.prtcls_per_cell) + stain1_background;
S2_max = stain2_part * PARAMETERS.prtcls_per_cell + stain2_background;

% Find the oldest cell generation present at the end of the simulation
oldest_cell_gen = max(EVOLUTION_INFO.cell_lineage_history(:,3:end),[],'all');
N_tstep_old = 1;

fig2 = figure(2);
set(fig2, 'Visible', 'off');
fig3 = figure(3);
set(fig3, 'Visible', 'off');
fig4 = figure(4);
set(fig4, 'Visible', 'off');
fig8 = figure(8);
set(fig8, 'Visible', 'off');
for time_plot = 1:floor(PARAMETERS.simulation_duration)+1
    tstep = (time_plot-1)/PARAMETERS.tstep_duration; % equivalent timestep index
    N_tstep = EVOLUTION_INFO.cell_population(tstep+1);
    
    % FREQUENCY HISTOGRAMS
    set(0,'CurrentFigure',fig2)
    subplot(5,floor(PARAMETERS.simulation_duration/5)+1,time_plot);
    % FREQUENCY OF CELLS WITH NUMS OF PARTICLES INTERACTING OVER TIME
    histogram(EVOLUTION_INFO.cell_c_o_p(1:N_tstep,time_plot,2),...
        'FaceColor', [0,0,1], 'FaceAlpha', 0.2);
    hold on;
    % FREQUENCY OF CELLS WITH NUMS OF PARTICLES INTERNALISED OVER TIME
    histogram(EVOLUTION_INFO.cell_c_o_p(1:N_tstep,time_plot,3),...
        'FaceColor', [1,0,0], 'FaceAlpha', 0.2);
    hold off;
    xlim([0,x_max]);
    ylim([0,y_max]);
    title(['At ' num2str(time_plot-1) ' hour/s']);
    xlabel('Number of particles');
    ylabel('Cell frequency');
    legend('Interacting','Internalised');
    
    % FLUORESCENCE DOT PLOTS
    set(0,'CurrentFigure',fig3)
    subplot(5,floor(PARAMETERS.simulation_duration/5)+1,time_plot);
    % PLOT EXTERNAL (STAIN 2) FLUORESCENCE AGAINST ASSOCIATED (STAIN 1)
    % FLUORESCENCE
    n_interact = EVOLUTION_INFO.cell_c_o_p(1:N_tstep,time_plot,2);
    n_internal = EVOLUTION_INFO.cell_c_o_p(1:N_tstep,time_plot,3);
    S1 = stain1_std_dev .* randn(N_tstep,1) + (stain1_part .* (n_interact + ...
        alpha.*n_internal) + stain1_background); % signal from stain 1 (x-axis)
    S2 = stain2_std_dev .* randn(N_tstep,1) + (stain2_part .* n_interact + ...
        stain2_background); % signal from stain 2 (y-axis)
    S1c = S1; S2c = S2; % Specifically for working out the coour densities
    S1c(S1c<=0)=realmin; % So that 0 is not passed through log
    S2c(S2c<=0)=realmin; % So that 0 is not passed through log
    %S1c = S1c/S1_max; % To normalise so that it's between 0 and 1
    %S2c = S2c/S2_max; % To normalise so that it's between 0 and 1
    c = ksdensity([log10(S1c),log10(S2c)], [log10(S1c),log10(S2c)]);
    dot_size = 3*ones(length(S1),1);
    scatter(S1, S2, dot_size, c, 'filled');
    set(gca, 'YScale', 'log','Xscale', 'log','XMinorTick','on','YMinorTick','on');
    xlim([1,S1_max]);
    ylim([1,S2_max]);
    title(['At ' num2str(time_plot-1) ' hour/s']);
    xlabel('Stain 1 (associated)');
    ylabel('Stain 2 (external)');
    cb = colorbar();
    cb.Label.String = 'Density estimate';
    
    if mod(time_plot-1,6)==0 || time_plot-1==0
        % FREQUENCY HISTOGRAMS PER GENERATION EVERY 6 HOURS
        set(0,'CurrentFigure',fig4)
        colmn = (time_plot-1)/6 + 1;
        rows = unique(EVOLUTION_INFO.cell_lineage_history(1:N_tstep,colmn+2));
        for row = rows' % ITERATING THROUGH THE CELL GENS PRESENT AT THIS TIMESTEP
            plot_num = (row - 1)*(floor(PARAMETERS.simulation_duration/6) + 1) + colmn;
            subplot(oldest_cell_gen,floor(PARAMETERS.simulation_duration/6) + 1,plot_num);
            % FREQUENCY OF CELLS WITH NUMS OF PARTICLES INTERACTING IN THIS
            % GEN
            is_in_gen = (EVOLUTION_INFO.cell_lineage_history(1:N_tstep,colmn+2)==row);
            cells_in_gen = EVOLUTION_INFO.cell_lineage_history(is_in_gen,2);
            histogram(EVOLUTION_INFO.cell_c_o_p(cells_in_gen,time_plot,2),...
                'FaceColor', [0,0,1], 'FaceAlpha', 0.2);
            hold on;
            % FREQUENCY OF CELLS WITH NUMS OF PARTICLES INTERNALISED IN
            % THIS GEN
            histogram(EVOLUTION_INFO.cell_c_o_p(cells_in_gen,time_plot,3),...
                'FaceColor', [1,0,0], 'FaceAlpha', 0.2);
            hold off;
            xlim([0,x_max]);
            title(['At ' num2str(time_plot-1) ' hour/s, generation ' num2str(row) ' cells']);
            xlabel('Num. of particles');
            ylabel('Cell frequ.');
            if colmn==1
                legend('Interacting','Internalised');
            end
        end
        
        % 3D FREQUENCY HISTOGRAMS EVERY 6 HOURS
        set(0,'CurrentFigure',fig8)
        subplot(2,floor(PARAMETERS.simulation_duration/12) + 1,colmn);
        X_3D = [EVOLUTION_INFO.cell_c_o_p(1:N_tstep,time_plot,2),EVOLUTION_INFO.cell_c_o_p(1:N_tstep,time_plot,3)];
        hist3(X_3D);
        title(['At ' num2str(time_plot-1) ' hour/s']);
        xlabel('Bound');
        ylabel('Internalised');
    end
    N_tstep_old = N_tstep+1;
end
sgtitle(fig2, 'Frequency of cells with certain numbers of interacting/internalised particles over time')
fig2.Position = [100,100,1300,700];
savefig(fig2, [PARAMETERS.folder_path '\Frequency_histograms'], 'compact')
saveas(fig2, [PARAMETERS.folder_path '\Frequency_histograms'], 'png')

sgtitle(fig3, 'External (stain 2) fluorescence of cells against associated (stain 1) fluorescence')
fig3.Position = [100,100,1300,700];
savefig(fig3, [PARAMETERS.folder_path '\Fluorescence_dot_plots'], 'compact')
saveas(fig3, [PARAMETERS.folder_path '\Fluorescence_dot_plots'], 'png')

sgtitle(fig4, 'Frequency of cells over time ... split into cell generation')
fig4.Position = [20,20,1300,700];
savefig(fig4, [PARAMETERS.folder_path '\Hist_split_gen'], 'compact')
saveas(fig4, [PARAMETERS.folder_path '\Hist_split_gen'], 'png')

sgtitle(fig8, '3D frequency of cells with certain numbers of interacting/internalised particles over time')
fig8.Position = [20,20,1300,700];
savefig(fig8, [PARAMETERS.folder_path '\3D_Hist'], 'compact')
saveas(fig8, [PARAMETERS.folder_path '\3D_Hist'], 'png')

% Plot the density of cells in different generations
fig5=figure(5);
set(fig5, 'Visible', 'off');
histogram(EVOLUTION_INFO.cell_lineage_history(:,end),oldest_cell_gen);
title(['Population of cell generations after hour' num2str(floor(PARAMETERS.simulation_duration/6) *6)]) 
xlabel('Generation')
ylabel('Frequency')
fig5.Position = [100,100,1300,700];
savefig(fig5, [PARAMETERS.folder_path '\Gen_densities'], 'compact')
saveas(fig5, [PARAMETERS.folder_path '\Gen_densities'], 'png')

% Plot proliferation events over time
fig6=figure(6);
set(fig6, 'Visible', 'off');
scatter(0:(1/6):24,EVOLUTION_INFO.cell_population,4,'k','filled');
title('Cell population over time')
xlabel('time $t$ hours', 'Interpreter', 'latex');
ylabel('Cell population');
fig6.Position = [100,100,1300,700];
savefig(fig6, [PARAMETERS.folder_path '\Cell_popoulation'], 'compact');
saveas(fig6, [PARAMETERS.folder_path '\Cell_popoulation'], 'png');

% Plot associated Vs interacting+internalised cells over time
d1 = EVOLUTION_INFO.average_c_o_p(2,:);     % Average interacting particles per cell over time
d2 = EVOLUTION_INFO.average_c_o_p(3,:);     % Average internalised particles per cell over time
d3 = d1 + d2;                               % Histogram Sum ‘d1’+‘d2’
binrng = 0:(1/6):24;                        % Create Bin Ranges

fig7=figure(7);
set(fig7, 'Visible', 'off');
subplot(1,2,1)
bar(binrng, d3, 'r')
hold on
bar(binrng, d1, 'b')
hold off
legend('Internalised','Interacting','Location','Best')
title('Stacked histogram')
xlabel('time $t$ hours', 'Interpreter', 'latex');
ylabel('Average number of particles per cell');

subplot(1,2,2)
plot(binrng,d2,'r--',binrng,d1,'b--');
hold on
plot(binrng,d3,'Color',[.5 0 .5])
hold off
legend('Internalised','Interacting','Associated (internalised or interacting)','Location','Best')
title('Separate trends')
xlabel('time $t$ hours', 'Interpreter', 'latex');

sgtitle('Average number of particles associated per cell')
fig7.Position = [100,100,1300,700];
savefig(fig7, [PARAMETERS.folder_path '\Associated_vs_internalised'], 'compact')
saveas(fig7, [PARAMETERS.folder_path '\Associated_vs_internalised'], 'png')

% Plot associated Vs interacting+internalised cells over time split by gen
%{
fig8=figure(8);
set(fig7, 'Visible', 'off');
for gen = 1:oldest_cell_gen
    subplot(1,oldest_cell_gen,gen)
end
sgtitle('Average number of particles associated per cell ... split into cell generation')
fig8.Position = [100,100,1300,700];
savefig(fig8, [PARAMETERS.folder_path '\Associated_vs_internalised'], 'compact')
saveas(fig8, [PARAMETERS.folder_path '\Associated_vs_internalised'], 'png')
%}

%save([PARAMETERS.folder_path '\variables.mat']);

figure(fig2)
figure(fig3)
figure(fig4)
figure(fig5)
figure(fig6)
figure(fig7)
figure(fig8)
