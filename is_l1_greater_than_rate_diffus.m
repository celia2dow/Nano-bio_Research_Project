function check = is_l1_greater_than_rate_diffus(PARAMETERS)
% IS_L1_GREATER_THAN_RATE_DIFFUS checks if the input parameters for 
% internalisation are appropriate. I.e. it checks that they don't require
% that the probability of binding once hit be greater than 1.

% CALCULATE RATE OF DIFFUSIVITY (PER HOUR)
% Rate of particle diffusivity (i.e. of particle hitting).
% Use Stokes-Einstein equation to calculate the diffusion coefficient of
% the particle in the given conditions:
k_B = 1.380649E-23; % Boltzmann constant
T_K = PARAMETERS.temp + 273.15; % Temperature in kelvin
r_m = PARAMETERS.prtcl_rad * (1E-9); % Particle radius in meters
h = PARAMETERS.culture_media_height * (1e-3); % Culture media height in meters
ratio_cell2site_area = pi /4; % Ratio of cell area to the lattice site it sits on
eta = PARAMETERS.viscos;
diffus_coeff = k_B * T_K / (6 * pi * eta * r_m) ... % In m^2 per second
    * (60^2);                                       % In m^2 per hour
% Divide by the culture height squared to get the rate of transition
% through the bottom of the column on the lattice site, then multiply by
% the ratio of this square area to the cell area (approximated as a flat
% disk) to get the rate of hitting.
rate_diffus = diffus_coeff * ratio_cell2site_area /(h^2); % In per hour

% CALCULATE ACTUAL RATE 1 (PER HOUR)
if PARAMETERS.EWTs_internalise.input_type == "EWT"
    % A list of L rates of particles transitioning between stages of the cell-
    % particle interaction model unaffected by carrying capacity (e.g. free to 
    % stage 1, stage 1 to stage 2, ..., stage L-1 to stage L or internalised).
    % Can be 1 rate per transition or can be 1 rate per cell phase per
    % transition - columns inicate interaction stage and rows indicate cell
    % phase.
    l1 = 1./PARAMETERS.EWTs_internalise.values(1);
elseif PARAMETERS.EWTs_internalise.input_type == "fraction"
    % Calculate the rates from the desired observed fraction
    % associated/internalised at confluence without carrying capacity.
    [l1,~] = input_EWT_from_fraction(...
        PARAMETERS.EWTs_internalise.values(1),...
        PARAMETERS.EWTs_internalise.values(2),...
        PARAMETERS.EWTs_internalise.values(3));
elseif PARAMETERS.EWTs_internalise.input_type == "prob_and_rates"
    % If the input are just the rates (per hour) from one stage to the
    % next, with the first value being the probability of binding once hit
    % (and therefore not needing scaling by the timestep duration)
    l1 = rate_diffus * PARAMETERS.EWTs_internalise.values(1) / ...
        PARAMETERS.tstep_duration;
end

% CHECK IF L1 IS LARGER THAN RATE_DIFFUS
if l1 > rate_diffus
    check = 1;
else
    check = 0;
end
end