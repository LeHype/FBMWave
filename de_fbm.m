clearvars
tic

% parpool(4);
n_dec = 100; % for each wind turbine two coordinates are needed
lb = -ones(n_dec, 1)*10;
ub = ones(n_dec, 1)*10;
%% DE parameters
n_mem                 = 50; % number of members for each generation, 20 <= n_mem <= 100
constraint_compliance = "flip"; % Select if and how to handle wind turbines leaving the square
% Possible values: "none", "border" or "flip"

donor_style = 'rand'; % Select if the first point of the DE mutation is a random point ('rand') or the best
% individuum ('best')

don_probabability     = 0.85; % probability of interchanging values with the donor vector
F                     = @(x) randn( 1, n_mem ).^2; % scaling factor of DE
% F                   = @(x) 0.65*ones(1,n_mem);

use_constr_values = true;
constr_weight = 1;

% stopping criterion is a min number of evolutions, a relative enhancement
% over the last 1000 iterations, a max number of evolutions, and a measure
% how far the population is spread
stopping_criterion    = @(nEvo, b, fit, cv) (nEvo > 1e2 && abs((b(end-99)-b(end))/b(end)) <= 0.001) || nEvo > 5e3 || ...
    (abs((min(fit)-max(fit))/min(fit)) < 1e-6 && min(cv) == 0) || max(fit)-min(fit)<=1e-6;

%% Initialize population
decvars = [unifrnd( lb.*ones(1, n_mem), ub.*ones(1, n_mem), n_dec, n_mem )];
fitness = match_random(decvars);
[constr_violations, constr_value] = constr_fun(decvars);

%%
nEvolutions = 0;

lowest_constr_value_idc = find(constr_violations == min(constr_violations));
[~, lowest_fitness_idx] = min(fitness(lowest_constr_value_idc));
best_idx = lowest_constr_value_idc(lowest_fitness_idx);
best_fitness = [fitness(best_idx)];

%do evolutions until stopping criterion is fulfilled
while ~stopping_criterion(nEvolutions, best_fitness, fitness, constr_violations)
    % randomly select vectors for generating the donor vector
    donor_parents = zeros(3, n_mem);
    for iMember = 1:n_mem
        available = setdiff(1:n_mem, iMember);
        donor_parents(1:3, iMember) = available(randperm(n_mem-1, 3))';
    end

    switch donor_style
        case 'rand'
            donor_points = decvars(:, donor_parents(1,:));
        case 'best'
            donor_points = decvars(:, best_idx);
    end

    % build the donor vector
    donor_vectors = donor_points + F().*(decvars(:, donor_parents(2,:)) -...
        decvars(:, donor_parents(3,:)));

    % select the positions where the donor vector donates his values
    donate = rand(n_dec, n_mem) < don_probabability;
    donate(randi(n_dec), all(donate == 0,1)) = 1;

    % build offspring
    offspring = decvars .* (~donate) + donor_vectors .* (donate);

    % prevent constraint violations if a wind turbine leaves the square
    % push it back into it or make it flip along the borders

    switch constraint_compliance
        case "border"
            offspring = min(max(offspring, lb), ub);
        case "flip"
            ub_mat = ub.*ones(1,n_mem);
            lb_mat = lb.*ones(1,n_mem);
            b_mat  = (ub-lb).*ones(1,n_mem);

            under_lb = offspring < lb_mat;
            over_ub = offspring > ub_mat;

            offspring(under_lb) = ub_mat(under_lb) + mod(offspring(under_lb), -b_mat(under_lb));
            offspring(over_ub) = lb_mat(over_ub) + mod(offspring(over_ub), b_mat(over_ub));
        otherwise
    end

    offspring_fitness = match_random(offspring);
    [offspring_violations, offspring_value] = constr_fun(offspring);

    % Did the offspring reduce the number of violations?
    reduced_violations = offspring_violations < constr_violations;
    same_violations    = offspring_violations == constr_violations;

    if use_constr_values
        % Did the offspring violate the constraints less than its parent?
        reduced_constr_value  = constr_weight*offspring_value < constr_weight*constr_value;
        same_constr_values    = constr_weight*offspring_value == constr_weight*constr_value;

        % Replace the parent if its offspring has either less constraint
        % violations, the same number of violations, but doesn't violate
        % them as much, or, if both are the same, if the offspring has at
        % least the same fitness value.
        replace = reduced_violations | (same_violations & reduced_constr_value) | ...
            (same_violations & same_constr_values & offspring_fitness <= fitness);
    else
        % Replace the parent if its offspring has either less constraint
        % violations or the same and an equal or lower fitness value.
        replace = reduced_violations | (same_violations & offspring_fitness <= fitness);
    end

    decvars = [decvars(:, ~replace) offspring(:, replace)];

    % calculate fitness and number of violated constraints
    fitness = [fitness(:, ~replace) offspring_fitness(:, replace)];
    constr_violations = [constr_violations(:, ~replace) offspring_violations(:, replace)];
    constr_value = [constr_value(:, ~replace) offspring_value(:, replace)];

    % select the best individuum
    lowest_constr_value_idc = find(constr_violations == min(constr_violations));
    [~, lowest_fitness_idx] = min(fitness(lowest_constr_value_idc));
    best_idx = lowest_constr_value_idc(lowest_fitness_idx);
    best_fitness = [best_fitness, fitness(best_idx)];

    nEvolutions = nEvolutions + 1;
end
toc


function fitness = match_random(decvars)
% use parfeval for calculating the histograms and matching them
    fitness = sum(decvars.^2, 1);
    %% Here is your interface. Generate Function for fitness and return [1 n_mem].
end

function [violations, value] = constr_fun(decvars)
% constr function that does not have any constraints (I was too lazy to
% delete all the constraints related code)
    violations = zeros(1, size(decvars, 2));
    value = zeros(1, size(decvars, 2));
end
