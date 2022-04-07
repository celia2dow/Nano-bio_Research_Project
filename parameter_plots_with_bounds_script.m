% SCRIPT for creating plots of parameter estimates and their lower and
% upper bounds

% Lambda 1 values
fig27 = figure(27);
set(fig27, 'Visible', 'off');
% MLE
plot(binrng(1:end-1), est_lambda1.MLE, 'k--', 'LineWidth', 1.5)
hold on
plot(binrng(1:end-1), est_lambda1UP.MLE, 'k:')
plot(binrng(1:end-1), est_lambda1LO.MLE, 'k:')
hold off
title('\lambda_1 calculated via MLE')
pbaspect([1 1 1])

if L == 2
    % Lambda 2 values
    fig28 = figure(28);
    set(fig28, 'Visible', 'off');

    % Differences method
    plot(binrng(2:end-1),est_lambda2.using_diffs, 'r--', 'LineWidth', 1.5)
    hold on
    plot(binrng(2:end-1),est_lambda2UP.using_diffs, 'r:')
    plot(binrng(2:end-1),est_lambda2LO.using_diffs, 'r:')
    % Mix method
    plot(binrng(1:end-1),est_lambda2.mix, 'b--', 'LineWidth', 1.5)
    plot(binrng(1:end-1),est_lambda2UP.mix, 'b:')
    plot(binrng(1:end-1),est_lambda2LO.mix, 'b:')
    % Disribution method
    plot(binrng,est_lambda2.distrib, 'g--', 'LineWidth', 1.5)
    plot(binrng,est_lambda2UP.distrib, 'g:')
    plot(binrng,est_lambda2LO.distrib, 'g:')
    ylim([min(est_lambda2LO.distrib), max(est_lambda2UP.distrib)])
    hold off
    title('\lambda_2 calculated via differences (red), distribution (green) and mix (blue)')
    pbaspect([1 1 1])

    fig28.Position = [100,100,1300,700];
    saveas(fig28, [PARAMETERS.folder_path '/lambda2_error'], 'eps')
    saveas(fig28, [PARAMETERS.folder_path '/lambda2_error'], 'png')

    figure(fig28)
end

% CCvalues
if PARAMETERS.max_prtcls(end) ~= inf
    fig29 = figure(29);
    set(fig29, 'Visible', 'off');
    if L==1
        % MLE Poisson method
        plot(binrng(1:end-1),est_CC.using_MLE_Poisson, 'r--', 'LineWidth', 1.5)
        hold on
        plot(binrng(1:end-1),est_CCUP.using_MLE_Poisson, 'r:')
        plot(binrng(1:end-1),est_CCLO.using_MLE_Poisson, 'r:')
        hold off
    elseif L==2
        % Differences method
        plot(binrng(1:end-1),est_CC.using_diffs, 'r--', 'LineWidth', 1.5)
        hold on
        plot(binrng(1:end-1),est_CCUP.using_diffs, 'r:')
        plot(binrng(1:end-1),est_CCLO.using_diffs, 'r:')
        % Mix method
        plot(binrng(1:end-1),est_CC.dynamic_diffs, 'm--', 'LineWidth', 1.5)
        plot(binrng(1:end-1),est_CCUP.dynamic_diffs, 'm:')
        plot(binrng(1:end-1),est_CCLO.dynamic_diffs, 'm:')
        % Disribution method
        plot(binrng(1,end-1),est_CC.dynamic_mix, 'c--', 'LineWidth', 1.5)
        plot(binrng(1:end-1),est_CCUP.dynamic_mix, 'c:')
        plot(binrng(1:end-1),est_CCLO.dynamic_mix, 'c:')
        ylim([0 2*PARAMETERS.max_prtcls(end)])
        hold off
    end
    title('CC calculated via differences (red), dynamic differences (magenta) and dynamic mix (cyan)')
    pbaspect([1 1 1]);
    fig29.Position = [100,100,1300,700];
    saveas(fig29, [PARAMETERS.folder_path '/CC_error'], 'eps')
    saveas(fig29, [PARAMETERS.folder_path '/CC_error'], 'png')

    figure(fig29)
end

% Comparing data to estimate rates:
% Data and one standard deviation error:
fig30 = figure(30);
set(fig30, 'Visible', 'off');
p1 = patch([binrng fliplr(binrng)], [upper1(2,:) fliplr(lower1(2,:))], [1 0 0], ...
    'FaceAlpha', 0.2, 'EdgeColor', 'w', 'LineStyle', '-', 'DisplayName','1 standard deviation');
hold on
if L>1
    p2 = patch([binrng fliplr(binrng)], [upper1(1,:) fliplr(lower1(1,:))], [0 0 1], ...
        'FaceAlpha', 0.2, 'EdgeColor', 'w', 'LineStyle', '-', 'DisplayName','" " "');
    p3 = patch([binrng fliplr(binrng)], [upper1(3,:) fliplr(lower1(3,:))], [.5 0 .5], ...
        'FaceAlpha', 0.2, 'EdgeColor', 'w', 'LineStyle', '-', 'DisplayName','" " "');
    p4 = plot(binrng,means(1,:),'b', 'DisplayName','Mean bound');
    p5 = plot(binrng,means(3,:),'Color',[.5 0 .5], 'DisplayName','Mean associated (either)');
end
p6 = plot(binrng,means(2,:),'r', 'DisplayName','Mean internalised');

% Use (Poisson inspired) MLE for lambda1
lam1 = est_lambda1.mean; 
lam1UP = est_lambda1UP.mean; 
lam1LO = est_lambda1LO.mean;
if L==2
    % Use hypoexponential distribution for lambda2
    lam2 = est_lambda2.distrib_mean; 
    lam2UP = est_lambda2UP.distrib_mean; 
    lam2LO = est_lambda2LO.distrib_mean;
end

% We want to scale lambda1 by the inverse of what we scaled it by to
% calculate it - i.e., divide by the current total cell population and
% multiply by the previous total cell population. In this way, lambda1 will
% decrease over time to compensate for the change in confluence affecting
% particle uptake
if total.cell_population(end) ~= total.cell_population(1)
    EWT_prolif = sum(PARAMETERS.EWTs_proliferate);
    % Using exponential growth formula
    proliferations = total.cell_population(1) .* 2 .^ ...
        ([binrng,binrng(end)+PARAMETERS.tstep_duration] ./ EWT_prolif);
    scaling = proliferations(1:end-1)./proliferations(2:end);
    lam1 = lam1 .* scaling;
end

% Association curves: mean, lower, upper
if total.cell_population(1)==total.cell_population(end)
    assoc = PARAMETERS.prtcls_per_site .* CDF.exp(...
        total.confluence(1).*lam1,binrng) ./ ...
        total.confluence(1);
    assocLO = PARAMETERS.prtcls_per_site .* CDF.exp(...
        total.confluence(1).*lam1LO,binrng)./ ...
        total.confluence(1);
    assocUP = PARAMETERS.prtcls_per_site .* CDF.exp(...
        total.confluence(1).*lam1UP,binrng)./ ...
        total.confluence(1);
else
    fprintf("\nWARNING confluence not accounted for in estimates\n")
    assoc = PARAMETERS.prtcls_per_site .* CDF.exp(lam1,binrng);
    assocLO = PARAMETERS.prtcls_per_site .* CDF.exp(lam1LO,binrng);
    assocUP = PARAMETERS.prtcls_per_site .* CDF.exp(lam1UP,binrng);
end

if PARAMETERS.max_prtcls(end) ~= inf
    if L==1
        CC = est_CC.using_MLE_Poisson; 
        CCUP = est_CCUP.using_MLE_Poisson; 
        CCLO = est_CCLO.using_MLE_Poisson;
        interns = analytic_distrib_with_CC(binrng,lam1,CC,PARAMETERS, ...
            CDF.hypoexp_l1,CDF.exp,assoc,L);
        internsLO = analytic_distrib_with_CC(binrng,lam1LO,CC,PARAMETERS, ...
            CDF.hypoexp_l1UP,CDF.exp,assocUP,L);
        internsUP = analytic_distrib_with_CC(binrng,lam1UP,CC,PARAMETERS, ...
            CDF.hypoexp_l1LO,CDF.exp,assocLO,L);
        sub_str = ['Heuristic estimates: \lambda_1: ' num2str(est_lambda1.mean)...
            ', CC: ' num2str(CC)];
        tit_str = ['Internalised per cell ~ true values:'...
            ' \lambda_1: ' num2str(lambdas(1)) ...
            ', CC: ' num2str(PARAMETERS.max_prtcls(end))];
    elseif L==2
        % Use differences for CC.
        CC = est_CC.using_diffs_mean; 
        CCUP = est_CCUP.using_diffs_mean; 
        CCLO = est_CCLO.using_diffs_mean;
        interns = analytic_distrib_with_CC(binrng,lam2,CC,PARAMETERS, ...
            CDF.hypoexp_l1,CDF.exp,assoc,L);
        internsLO = analytic_distrib_with_CC(binrng,lam2LO,CC,PARAMETERS, ...
            CDF.hypoexp_l1UP,CDF.exp,assocUP,L);
        internsUP = analytic_distrib_with_CC(binrng,lam2UP,CC,PARAMETERS, ...
            CDF.hypoexp_l1LO,CDF.exp,assocLO,L);
        sub_str = ['Heuristic estimates: \lambda_1: ' num2str(est_lambda1.mean)...
            ', \lambda_2: ' num2str(lam2)...
            ', CC: ' num2str(CC)];
        tit_str = ['Associated, bound, internalised per cell ~ true values:'...
            ' \lambda_1: ' num2str(lambdas(1)) ...
            ', \lambda_2: ' num2str(lambdas(2)) ...
            ', CC: ' num2str(PARAMETERS.max_prtcls(end))];
    end
else 
    if L==1
        interns = assoc;
        internsLO = assocLO;
        internsUP = assocUP;
        sub_str = ['Heuristic estimates: \lambda_1: ' num2str(est_lambda1.mean)];
        tit_str = ['Internalised per cell ~ true values: '...
            ' \lambda_1: ' num2str(lambdas(1))];
    elseif L==2
        % Internalisation curves: mean, lower, upper
        interns = PARAMETERS.prtcls_per_site .* CDF.hypoexp_l1(lam2,binrng);
        internsLO = PARAMETERS.prtcls_per_site .* CDF.hypoexp_l1UP(lam2LO,binrng);
        internsUP = PARAMETERS.prtcls_per_site .* CDF.hypoexp_l1LO(lam2UP,binrng);
        sub_str = ['Heuristic estimates: \lambda_1: ' num2str(est_lambda1.mean)...
            ', \lambda_2: ' num2str(lam2)];
        tit_str = ['Associated, bound, internalised per cell ~ true values: '...
            ' \lambda_1: ' num2str(lambdas(1)) ...
            ', \lambda_2: ' num2str(lambdas(2))];
    else
        sub_str = ['Heuristic estimates: \lambda_1: ' num2str(est_lambda1.mean)];
        tit_str = ['Associated per cell ~ true values: '...
            ' \lambda_1: ' num2str(lambdas(1))];
    end
end

% Plot ALL
if L<=2
    % Interacting curves: mean, lower, upper
    interacts = assoc - interns;
    interactsLO = assoc - internsUP;
    interactsUP = assoc - internsLO;

    p7 = plot(binrng,interns,'k--', 'LineWidth', 1.2, ...
        'DisplayName', 'Mean analytic distributions with heuristic estimates');
    %p8 = plot(binrng,internsLO,'k:', 'LineWidth', 1.2, ...
    %    'DisplayName', 'Lower/upper bound analytic distributions with heuristic estimates');
    %p9 = plot(binrng,internsUP,'k:', 'LineWidth', 1.2);
    %hasbehavior(p9,'legend',false);
    if L==2
        p10 = plot(binrng,assoc,'k--', 'LineWidth', 1.2);
        hasbehavior(p10,'legend',false);
        %p11 = plot(binrng,assocLO,'k:', 'LineWidth', 1.2,...
        %hasbehavior(p11,'legend',false);
        %p12 = plot(binrng,assocUP,'k:', 'LineWidth', 1.2);
        %hasbehavior(p12,'legend',false);
        
        p13 = plot(binrng,interacts,'k--', 'LineWidth', 1.2);
        hasbehavior(p13,'legend',false);
        %p14 = plot(binrng,interactsLO,'k:', 'LineWidth', 1.2);
        %hasbehavior(p14,'legend',false);
        %p15 = plot(binrng,interactsUP,'k:', 'LineWidth', 1.2);
        %hasbehavior(p15,'legend',false);
    end
else
    p10 = plot(binrng,assoc,'k--', 'LineWidth', 1.2, ...
        'DisplayName', 'Mean analytic distributions with heuristic estimates');
    %p11 = plot(binrng,assocLO,'k:', 'LineWidth', 1.2,...
    %hasbehavior(p11,'legend',false);
    %p12 = plot(binrng,assocUP,'k:', 'LineWidth', 1.2);
    %hasbehavior(p12,'legend',false);
end

legend('Location','Best')
ylim([0 max(upper1,[],'all')])
title(tit_str);
subtitle(sub_str);
ylabel('Num. particles', 'Interpreter', 'latex')
xlabel('time $t$ hours', 'Interpreter', 'latex');
pbaspect([1 1 1])

% Printing time series averages split by number of cell divisions
fig31 = figure(31);
set(fig31, 'Visible', 'off');
max_cell_divs = max(total.cell_lineage(:,3:end),[],'all');
store_means = zeros(max_cell_divs,total_tsteps+1,3);
store_stddev = zeros(max_cell_divs,total_tsteps+1,3);
% Iterate through cell divisions
for num_divs = 1:max(total.cell_lineage(:,tstep+2)) 
    for tstep = 1:total_tsteps+1 
        cells_concerned = total.cell_lineage(total.cell_lineage(:,tstep+2)==num_divs,2);
        if cells_concerned
            interacting_at_tstep = mean(total.cell_c_o_p(cells_concerned,tstep,2));
            internalised_at_tstep = mean(total.cell_c_o_p(cells_concerned,tstep,3));
            associated_at_tstep = interacting_at_tstep + internalised_at_tstep;
            var_intera = var(total.cell_c_o_p(cells_concerned,tstep,2));
            var_intern = var(total.cell_c_o_p(cells_concerned,tstep,3));
            var_assoc = var(sum(total.cell_c_o_p(cells_concerned,tstep,2:3),1));
        else
            interacting_at_tstep = 0;
            internalised_at_tstep = 0;
            associated_at_tstep = 0;
            var_intera = 0;
            var_intern = 0;
            var_assoc = 0;
        end
        
        store_means(num_divs,tstep,:) = [interacting_at_tstep, ...
            internalised_at_tstep, associated_at_tstep];
        store_stddev(num_divs,tstep,:) = [sqrt(var_intera), ...
            sqrt(var_intern), sqrt(var_assoc)];
    end
    if any(total.cell_c_o_p(:,:,2)>0,"all")
        % Plot associated particles over time
        subplot(1,2,1)
        hold on
        pl1 = plot(binrng,store_means(num_divs,:,3), 'DisplayName',...
            [num2str(num_divs) ' cell divisions']);
        col1 = pl1.Color;
        %patch([binrng fliplr(binrng)], [squeeze(store_means(num_divs,:,3)+store_stddev(num_divs,:,3)) ...
        %    fliplr(squeeze(store_means(num_divs,:,3)-store_stddev(num_divs,:,3)))], col1, ...
        %    'FaceAlpha', 0.2, 'EdgeColor', 'w', 'LineStyle', '-', 'DisplayName','1 standard deviation');
        hold off
        ylabel('Num. particles', 'Interpreter', 'latex')
        xlabel('time $t$ hours', 'Interpreter', 'latex');
        title("Mean number of paticles associated")
        legend("Location","northwest")
        pbaspect([1 1 1])
        % Plot internalised particles over time
        subplot(1,2,2)
        hold on
        pl2 = plot(binrng,store_means(num_divs,:,2), 'DisplayName',...
            [num2str(num_divs) ' cell divisions']);
        col2 = pl2.Color;
        %patch([binrng fliplr(binrng)], [squeeze(store_means(num_divs,:,2)+store_stddev(num_divs,:,2)) ...
        %    fliplr(squeeze(store_means(num_divs,:,2)-store_stddev(num_divs,:,2)))], col2, ...
        %    'FaceAlpha', 0.2, 'EdgeColor', 'w', 'LineStyle', '-', 'DisplayName','1 standard deviation');
        hold off
        xlabel('time $t$ hours', 'Interpreter', 'latex');
        title("Mean number of paticles internalised")
        legend("Location","northwest")
        pbaspect([1 1 1])

        sgtitle('Particle dosages split by number of cell divisions');
    else
        % Plot internalised particles over time
        hold on
        pl1 = plot(binrng,store_means(num_divs,:,2), 'DisplayName',...
            ['Mean internalised for ' num2str(num_divs) ' cell divisions']);
        col1 = pl1.Color;
        %patch([binrng fliplr(binrng)], [squeeze(store_means(num_divs,:,2)+store_stddev(num_divs,:,2)) ...
        %    fliplr(squeeze(store_means(num_divs,:,2)-store_stddev(num_divs,:,2)))], col1, ...
        %    'FaceAlpha', 0.2, 'EdgeColor', 'w', 'LineStyle', '-', 'DisplayName','1 standard deviation');
        hold off
        ylabel('Num. particles', 'Interpreter', 'latex')
        xlabel('time $t$ hours', 'Interpreter', 'latex');
        title('Internalised particle dosages split by number of cell divisions');
        legend("Location","northwest")
        pbaspect([1 1 1])
    end

end

% Save figures
fig27.Position = [100,100,1300,700];
saveas(fig27, [PARAMETERS.folder_path '/lambda1_error'], 'eps')
saveas(fig27, [PARAMETERS.folder_path '/lambda1_error'], 'png')

fig30.Position = [100,100,1300,700];
saveas(fig30, [PARAMETERS.folder_path '/Analytical_distribs'], 'eps')
saveas(fig30, [PARAMETERS.folder_path '/Analytical_distribs'], 'png')

fig31.Position = [100,100,1300,700];
saveas(fig31, [PARAMETERS.folder_path '/Time_series_split_by_divs'], 'eps')
saveas(fig31, [PARAMETERS.folder_path '/Time_series_split_by_divs'], 'png')

figure(fig27)
figure(fig30)
figure(fig31)
