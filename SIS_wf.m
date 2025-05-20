function [BestChart, BestFitness, farmlayout, farmlayout_NA] = SIS_wf(wf, popuSize, iteration, run_id, func, algorithmDir)
    %% Initialization for a single optimization run
    BestChart = zeros(iteration, 1);            % Fitness value of the best individual each iteration
    BestFitness = zeros(iteration, 1);          % Number of turbines for the best individual
    farmlayout = zeros(iteration, wf.rows * wf.cols);
    farmlayout_NA = zeros(iteration, wf.rows * wf.cols);
    NA_loc = wf.NA_loc;
    NA_l = length(NA_loc);
    
    % Max function evaluations (FES)
    FES = 30000;
    nfes = 0;
    max_nfes = FES;

    %% Population Initialization
    D = 3;  % Dimension: [# turbines, chaos param1, chaos param2]
    LU = [1, 0, 0; (144 - NA_l), 1, 1];  % Lower/Upper bounds
    popu = repmat(LU(1, :), popuSize, 1) + rand(popuSize, D) .* (repmat(LU(2, :) - LU(1, :), popuSize, 1));

    % Evaluate initial population
    fitness = zeros(popuSize, 1);
    for popindex = 1:popuSize
        M = ceil(popu(popindex, 1));           % Number of turbines
        wf.turbine_num = M;
        Chaos = generate_chaos_sequence(popu(popindex, 2:3), M);
        w = ceil(Chaos * (wf.rows * wf.cols));
        w = windfarm_constraint(w, wf.NA_loc, M, 1, wf.rows * wf.cols);
        power_total = cal_P_rate_total(wf);
        [fitness(popindex), ~] = wf_fitness(wf, w);
    end

    %% SIS Parameters
    l = 0.5;
    R = 0.1 * rand;
    n = D;
    N0 = floor((popuSize)^(1 / n));
    constraints1 = ones(1, D);
    constraints2 = ones(1, D) * wf.rows * wf.cols;
    for i = 1:n
        delta_x_m(i) = (constraints2(i) - constraints1(i)) / N0;
    end

    %% Find best initial solution
    bsf_fit_var = 0;
    bsf_solution = zeros(1, D);
    for i = 1:popuSize
        nfes = nfes + 1;
        if fitness(i) > bsf_fit_var
            bsf_fit_var = fitness(i);
            bsf_solution = popu(i, :);
        end
        if nfes > max_nfes; break; end
    end

    %% Main SIS loop
    for iteration = 1:iteration
        c = floor(popuSize * l);   % Number of exploitation individuals
        cc = popuSize - c;         % Number of exploration individuals

        % --- Exploitation Phase ---
        x_sele1 = generate_sis_population(c, delta_x_m, bsf_solution);
        x_sele1 = BoundaryDetection(x_sele1, LU);
        fl = evaluate_population(x_sele1, wf);

        % --- Exploration Phase ---
        x_sele2 = generate_sis_population(cc, delta_x_m, bsf_solution);
        x_sele2 = BoundaryDetection(x_sele2, LU);
        fc = evaluate_population(x_sele2, wf);

        % --- Compare and Update ---
        [optimalFitl, optIl] = max(fl);
        [optimalFitc, optIc] = max(fc);

        if optimalFitc > optimalFitl
            delta_x_m = delta_x_m - R * delta_x_m;
            bsf_solution = x_sele2(optIc, :);
            best_fitness = optimalFitc;
        else
            delta_x_m = delta_x_m + R * delta_x_m;
            bsf_solution = x_sele1(optIl, :);
            best_fitness = optimalFitl;
        end

        nfes = nfes + popuSize;
        if bsf_fit_var < best_fitness
            bsf_fit_var = best_fitness;
        end
        if nfes > max_nfes; break; end

        %% Record current best
        M = ceil(bsf_solution(1));
        Chaos = generate_chaos_sequence(bsf_solution(2:3), M);
        w = ceil(Chaos * (wf.rows * wf.cols));
        w = windfarm_constraint(w, wf.NA_loc, M, 1, wf.rows * wf.cols);

        BestChart(iteration) = bsf_fit_var;
        BestFitness(iteration) = M;
        [best_farmlayout, best_farmlayout_NA] = gene_layout_by_indices_one(wf, w);
        farmlayout(iteration, :) = best_farmlayout;
        farmlayout_NA(iteration, :) = best_farmlayout_NA;

        fprintf('%s problem %s | Dim %d | Run %d | FES %d --> %.16f\n', ...
            algorithmDir, func, D, run_id, nfes, bsf_fit_var);
    end
end

%% ====== Utility Functions ======

% Calculate total expected power generation
function power_total = cal_P_rate_total(wf)
    f_p = 0.0;
    for ind_t = 1:length(wf.theta)
        for ind_v = 1:length(wf.velocity)
            f_p = f_p + wf.f_theta_v(ind_t, ind_v) * P_i_X(wf.velocity(ind_v));
        end
    end
    power_total = wf.turbine_num * f_p;
end

% Power output at given wind speed
function re = P_i_X(v)
    if v < 2.0
        re = 0;
    elseif v < 12.8
        re = 0.3 * v ^ 3;
    elseif v < 18
        re = 629.1;
    else
        re = 0;
    end
end

% Generate chaos sequence for layout encoding
function Chaos = generate_chaos_sequence(params, M)
    C1 = params(1);
    a1 = params(2);
    Chaos = zeros(1, M);
    Chaos(1) = C1;
    for j = 1:M-1
        if Chaos(j) < a1 && Chaos(j) > 0
            Chaos(j+1) = Chaos(j) / a1;
        elseif Chaos(j) < 1 && Chaos(j) >= a1
            Chaos(j+1) = (1 - Chaos(j)) / (1 - a1);
        end
    end
end

% Generate new population via SIS strategy
function X = generate_sis_population(N, delta_x_m, bsf_solution)
    n = length(bsf_solution);
    neighborhood = rand(N, n);
    X = zeros(N, n);
    for j = 1:n
        X(:, j) = (2 * neighborhood(:, j) - 1) * delta_x_m(j) + bsf_solution(j);
    end
end

% Evaluate population (with fitness)
function fitness = evaluate_population(population, wf)
    N = size(population, 1);
    fitness = zeros(N, 1);
    for i = 1:N
        Chaos = generate_chaos_sequence(population(i, 2:3), ceil(population(i, 1)));
        w = ceil(Chaos * (wf.rows * wf.cols));
        w = windfarm_constraint(w, wf.NA_loc, ceil(population(i, 1)), 1, wf.rows * wf.cols);
        wf.turbine_num = ceil(population(i, 1));
        cal_P_rate_total(wf);  % For side effect
        [fitness(i), ~] = wf_fitness(wf, w);
    end
end
