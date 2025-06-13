%% Main script: SIS-based Wind Farm Layout Optimization
clear; clc; tic;

% Add helper functions to path
addpath('./WindFarmOptimization');

% Algorithm settings
algorithmDir = 'SIS';     % Directory/name used for output files
func = 'WindFarm';        % Objective function name (for display)

% Problem parameters
rows          = 12;               % Grid rows
cols          = 12;               % Grid columns
wind_type     = [1,2,3,4];        % Wind profile scenarios
cell_width    = 77.0 * 3;         % Width of each cell (meters)
NA_type_list  = 0:12;             % Number of blocked cells variants
turbine_num   = [15];             % Number of turbines to test
popuSize      = 100;              % Population size for SIS
iteration     = 200;              % Maximum SIS iterations per run
runTime       = 2;                % Independent runs per configuration

% Only one algorithm variant is used here (type = 1)
for type = 1
    for wt = wind_type
        for tn = turbine_num
            % Prepare paths for box-plot and convergence Excel files
            pathOfBoxPlot      = ['./', algorithmDir, num2str(type), ...
                                  '_WT', num2str(wt), '_TN', num2str(tn), ...
                                  '_Mean&Std_and_Box-Plot.xls'];
            pathOfConvergence  = ['./', algorithmDir, num2str(type), ...
                                  '_WT', num2str(wt), '_TN', num2str(tn), ...
                                  '_Convergence.xls'];
            % Additional files (for fitness instead of objective)
            PathOfBoxPlot      = ['./', algorithmDir, num2str(type), ...
                                  '_WT', num2str(wt), '_M', num2str(tn), ...
                                  '_Mean&Std_and_Box-Plot.xls'];
            PathOfConvergence  = ['./', algorithmDir, num2str(type), ...
                                  '_WT', num2str(wt), '_M', num2str(tn), ...
                                  '_Convergence.xls'];

            % Initialize summary tables
            tableOfBoxPlot      = initTableforBoxPlot(length(NA_type_list), runTime);
            tableOfConvergence  = initTableforConvergence(length(NA_type_list));
            TableOfBoxPlot      = initTableforBoxPlot(length(NA_type_list), runTime);
            TableOfConvergence  = initTableforConvergence(length(NA_type_list));

            % Loop over each NA (blocked-cell) configuration
            for NA_type = NA_type_list
                % Generate blocked locations and wind farm object
                NA_loc_array = gene_NA_loc(NA_type);
                wf = gene_windfram(rows, cols, tn, cell_width, NA_loc_array, wt);

                % Create result folder if needed
                folder = sprintf('./results/%s/wind_profile%d/tn%d_NA%d', ...
                                 algorithmDir, wt, tn, NA_type);
                if ~exist(folder, 'dir')
                    mkdir(folder);
                end

                % Preallocate result matrices
                eta     = [];   % Best objective per iteration across runs
                fitness = [];   % Turbine count per iteration across runs

                % Run the SIS algorithm multiple times
                for t = 1:runTime
                    [BestChart, Bestfitness, farmlayout, farmlayout_NA] = ...
                        SIS_wf(wf, popuSize, iteration, t, func, algorithmDir);
                    eta(:, t)     = BestChart;
                    fitness(:, t) = Bestfitness;
                    % Save each run's layout
                    save_results(farmlayout, farmlayout_NA, t, folder);
                end

                % Save raw .mat data
                save(sprintf('%s/eta.mat',     folder), 'eta');
                save(sprintf('%s/fitness.mat', folder), 'fitness');

                % Extract final‐iteration and mean curves
                boxPlotChart     = eta(end, :);
                convergenceChart = mean(eta, 2);
                BoxPlotChart     = fitness(end, :);
                ConvergenceChart = mean(fitness, 2);

                % (Optional) write raw eta into a master Excel sheet
                path = ['./', algorithmDir, num2str(type), '_WT', ...
                        num2str(wt), '_TN', num2str(tn), '.xls'];
                sheetName = ['WT_' num2str(wt)];
                xlswrite(path, eta, [sheetName, num2str(NA_type+1), 'TN_', num2str(tn)]);

                % Update and save the summary tables
                tableOfBoxPlot     = updateTableforBoxPlot(...
                    tableOfBoxPlot, boxPlotChart,     pathOfBoxPlot,     0, NA_type+1);
                tableOfConvergence = updateTableforConvergence(...
                    tableOfConvergence, convergenceChart, pathOfConvergence, 0, NA_type+1);
                TableOfBoxPlot     = updateTableforBoxPlot(...
                    TableOfBoxPlot, BoxPlotChart,     PathOfBoxPlot,     0, NA_type+1);
                TableOfConvergence = updateTableforConvergence(...
                    TableOfConvergence, ConvergenceChart, PathOfConvergence, 0, NA_type+1);
            end
        end
    end
end

toc


%% initTableforBoxPlot
% Create a cell array for box‐plot statistics.
%   problemNum = number of problem variants (rows)
%   runNum     = number of independent runs (data columns)
function tableB = initTableforBoxPlot(problemNum, runNum)
    row = 3 + problemNum;      % header + problem rows + 2 extra
    col = 1 + runNum;          % index col + data columns
    tableB = cell(row, col);
    tableB{1, 2} = 'mean';     % column header
    tableB{1, 3} = 'std';      % column header
    for i = 1:problemNum
        tableB{i+1, 1} = ['F' num2str(i)];  % label each problem
    end
end

%% initTableforConvergence
% Create a cell array for convergence data.
%   problemNum = number of problem variants (columns)
function tableC = initTableforConvergence(problemNum)
    tableC = cell(1, problemNum);
    for i = 1:problemNum
        tableC{1, i} = ['F' num2str(i)];  % header labels
    end
end

%% updateTableforBoxPlot
% Fill in box‐plot table row and write to Excel.
function tableB = updateTableforBoxPlot(tableB, data, path, rgo, problemIndex)
    [~, col] = size(data);
    rgoVec    = repmat(rgo, 1, col);
    dataAdj   = data - rgoVec;
    dataMean  = mean(dataAdj, 2);
    dataStd   = std(dataAdj, 0, 2);
    % Convert to cells for assignment
    meanCell  = mat2cell(dataMean, ones(1,numel(dataMean)), 1);
    stdCell   = mat2cell(dataStd,  ones(1,numel(dataStd)),  1);
    dataCells = mat2cell(dataAdj, 1, ones(1,col));
    % Place into table (row = problemIndex+1)
    tableB(problemIndex+1, 2)    = meanCell;
    tableB(problemIndex+1, 3)    = stdCell;
    tableB(problemIndex+1, 4:3+col) = dataCells;
    % Write out
    xlswrite(path, tableB);
end

%% updateTableforConvergence
% Fill in convergence table column and write to Excel.
function tableC = updateTableforConvergence(tableC, data, path, rgo, problemIndex)
    [row, ~] = size(data);
    rgoMat    = repmat(rgo', 1, size(data,2));
    dataAdj   = data - rgoMat;
    dataCells = mat2cell(dataAdj, ones(1,row), 1);
    % Place into table (starting at row 2)
    tableC(2:row+1, problemIndex) = dataCells;
    % Write out
    xlswrite(path, tableC);
end

%% SIS_wf
% Execute one SIS optimization run on the wind farm.
function [BestChart, BestFitness, farmlayout, farmlayout_NA] = ...
         SIS_wf(wf, popuSize, iteration, run_id, func, algorithmDir)

    % Preallocate outputs
    BestChart       = zeros(iteration, 1);
    BestFitness     = zeros(iteration, 1);
    farmlayout      = zeros(iteration, wf.rows * wf.cols);
    farmlayout_NA   = zeros(iteration, wf.rows * wf.cols);
    
    NA_l = length(wf.NA_loc);
    FES  = 30000;    % evaluation budget
    nfes = 0;        % counter
    
    D = 3;           % [#turbines, C1, a1]
    % Bounds for each variable
    LU = [1, zeros(1,D-1); (144-NA_l), ones(1,D-1)];
    % Initialize population uniformly
    popu = repmat(LU(1,:), popuSize,1) + ...
           rand(popuSize,D).*repmat(LU(2,:)-LU(1,:), popuSize,1);

    % Evaluate initial population
    fitness = zeros(popuSize,1);
    for i = 1:popuSize
        C1 = popu(i,2);
        Chaos = C1;
        a1    = popu(i,3);
        M     = ceil(popu(i,1));
        wf.turbine_num = M;
        for j = 1:M-1
            if Chaos(j)<a1 && Chaos(j)>0
                Chaos(j+1)=Chaos(j)/a1;
            elseif Chaos(j)<1 && Chaos(j)>=a1
                Chaos(j+1)=(1-Chaos(j))/(1-a1);
            end
        end
        w = ceil(Chaos*(wf.rows*wf.cols));
        w = windfarm_constraint(w, wf.NA_loc, M,1,wf.rows*wf.cols);
        [fitness(i),~] = wf_fitness(wf,w);
    end

    % Best-so-far initialization
    [bsf_fit_var, idx] = max(fitness);
    bsf_solution       = popu(idx,:);

    % SIS parameters
    l = 0.5;               % elite fraction
    R = 0.1 * rand;        % step-size factor
    N = popuSize;          % population size
    N0 = floor(N^(1/D));   % grid points per dim
    constraints1 = ones(1,D);
    constraints2 = ones(1,D)*wf.rows*wf.cols;
    delta_x_m = (constraints2-constraints1)/N0;

    for iter = 1:iteration
        % Elite vs non-elite counts
        c  = floor(popuSize*l);
        cc = popuSize - c;

        % --- Elite generation via chaotic map ---
        z1 = rand(D,1);
        forbidden = [0.25,0.5,0.75];
        for s = 1:D
            while any(abs(z1(s)-forbidden)<eps)
                z1(s) = rand;
            end
        end
        a = 4.0;
        ChaosMat = zeros(D,c);
        ChaosMat(:,1)=z1;
        for h = 1:c-1
            ChaosMat(:,h+1)=a*ChaosMat(:,h)-a*ChaosMat(:,h).*ChaosMat(:,h);
        end
        x_sele1 = (2*ChaosMat'-1).*delta_x_m + bsf_solution;
        x_sele1 = BoundaryDetection(x_sele1,LU);

        fl = zeros(c,1);
        for i = 1:c
            C1 = x_sele1(i,2);
            Chaos = C1;
            a1    = x_sele1(i,3);
            M     = ceil(x_sele1(i,1));
            wf.turbine_num = M;
            for j = 1:M-1
                if Chaos(j)<a1 && Chaos(j)>0
                    Chaos(j+1)=Chaos(j)/a1;
                elseif Chaos(j)<1 && Chaos(j)>=a1
                    Chaos(j+1)=(1-Chaos(j))/(1-a1);
                end
            end
            w = ceil(Chaos*(wf.rows*wf.cols));
            w = windfarm_constraint(w, wf.NA_loc, M,1,wf.rows*wf.cols);
            [fl(i),~] = wf_fitness(wf,w);
        end
        [optFitL,optIl] = max(fl);

        % --- Non-elite generation via random sampling ---
        neighborhoodc = rand(cc,D);
        x_sele2 = (2*neighborhoodc-1).*delta_x_m + bsf_solution;
        x_sele2 = BoundaryDetection(x_sele2,LU);

        fc = zeros(cc,1);
        for i = 1:cc
            C1 = x_sele2(i,2);
            Chaos = C1;
            a1    = x_sele2(i,3);
            M     = ceil(x_sele2(i,1));
            wf.turbine_num = M;
            for j = 1:M-1
                if Chaos(j)<a1 && Chaos(j)>0
                    Chaos(j+1)=Chaos(j)/a1;
                elseif Chaos(j)<1 && Chaos(j)>=a1
                    Chaos(j+1)=(1-Chaos(j))/(1-a1);
                end
            end
            w = ceil(Chaos*(wf.rows*wf.cols));
            w = windfarm_constraint(w, wf.NA_loc, M,1,wf.rows*wf.cols);
            [fc(i),~] = wf_fitness(wf,w);
        end
        [optFitC,optIc] = max(fc);

        % --- Update best-so-far ---
        if optFitC > optFitL
            delta_x_m     = delta_x_m*(1-R);
            bsf_solution  = x_sele2(optIc,:);
            best_fit_val  = optFitC;
        else
            delta_x_m     = delta_x_m*(1+R);
            bsf_solution  = x_sele1(optIl,:);
            best_fit_val  = optFitL;
        end

        % Update FE count and best value
        nfes = nfes + popuSize;
        if best_fit_val > bsf_fit_var
            bsf_fit_var = best_fit_val;
        end
        if nfes > FES
            break;
        end

        % Record iteration results
        BestChart(iter)   = bsf_fit_var;
        BestFitness(iter) = ceil(bsf_solution(1));

        % Generate layout for best-so-far
        Chaos = bsf_solution(2);
        a1    = bsf_solution(3);
        M     = ceil(bsf_solution(1));
        for j = 1:M-1
            if Chaos(j)<a1 && Chaos(j)>0
                Chaos(j+1)=Chaos(j)/a1;
            elseif Chaos(j)<1 && Chaos(j)>=a1
                Chaos(j+1)=(1-Chaos(j))/(1-a1);
            end
        end
        w = ceil(Chaos*(wf.rows*wf.cols));
        w = windfarm_constraint(w, wf.NA_loc, M,1,wf.rows*wf.cols);
        [best_farmlayout, best_farmlayout_NA] = gene_layout_by_indices_one(wf, w);
        farmlayout(iter, :)    = best_farmlayout;
        farmlayout_NA(iter, :) = best_farmlayout_NA;

        % Display progress
        fprintf('%s problem %s | Dim %d | run %d | FES %d -----> %f\n', ...
                algorithmDir, func, D, run_id, nfes, bsf_fit_var);
    end
end

%% cal_P_rate_total
% Compute total power rate for the wind farm configuration
function power_total = cal_P_rate_total(wf)
    f_p = 0.0;
    for ind_t = 1:length(wf.theta)
        for ind_v = 1:length(wf.velocity)
            f_p = f_p + wf.f_theta_v(ind_t,ind_v) * P_i_X(wf.velocity(ind_v));
        end
    end
    power_total = wf.turbine_num * f_p;
end

%% P_i_X
% Turbine power curve: returns power output at wind speed v
function re = P_i_X(v)
    if v < 2.0
        re = 0;
    elseif v < 12.8
        re = 0.3 * v^3;
    elseif v < 18
        re = 629.1;
    else
        re = 0;
    end
end
