clear;
clc;
tic;

% Define function and add function path
func = 'WindFarm';
addpath('./WindFarmOptimization');

algorithmDir = 'SIS';  % Algorithm name or result directory

% Experiment parameters
rows = 12;
cols = 12;
wind_type = [1, 2, 3, 4];       % Types of wind conditions
cell_width = 77.0 * 3;         % Width of each grid cell
NA_type_list = 0:12;           % Types of no-build zones
turbine_num = [15];            % Number of turbines
popuSize = 100;                % Population size for algorithm
iteration = 200;               % Number of optimization iterations
runTime = 2;                   % Repetitions for each configuration

% Main experimental loop
for type = 1
    for wt = wind_type
        for tn = turbine_num
            % Prepare Excel output paths
            pathOfBoxPlot = ['./', algorithmDir, num2str(type), '_WT', num2str(wt), '_TN', num2str(tn), '_Mean&Std_and_Box-Plot.xls'];
            pathOfConvergence = ['./', algorithmDir, num2str(type), '_WT', num2str(wt), '_TN', num2str(tn), '_Convergence.xls'];

            % Initialize tables for logging results
            tableOfBoxPlot = initTableforBoxPlot(length(NA_type_list), runTime);
            tableOfConvergence = initTableforConvergence(length(NA_type_list));

            for NA_type = NA_type_list
                % Generate No-Build zone and wind farm structure
                NA_loc_array = gene_NA_loc(NA_type);
                wf = gene_windfram(rows, cols, tn, cell_width, NA_loc_array, wt);

                % Prepare folder to save results
                folder = sprintf('./results/%s/wind_profile%d/tn%d_NA%d', algorithmDir, wt, tn, NA_type);
                if ~exist(folder, 'dir')
                    mkdir(folder);
                end

                eta = [];
                fitness = [];

                % Run algorithm multiple times for statistical robustness
                for t = 1:runTime
                    [BestChart, Bestfitness, farmlayout, farmlayout_NA] = SIS_wf(wf, popuSize, iteration, t, func, algorithmDir);
                    eta(:,t) = BestChart;
                    fitness(:,t) = Bestfitness;

                    save_results(farmlayout, farmlayout_NA, t, folder);
                end

                % Save per-run eta and fitness values
                save(sprintf('%s/eta.mat', folder), "eta");
                save(sprintf('%s/fitness.mat', folder), "fitness");

                % Compute box plot and convergence data
                boxPlotChart = eta(end, :);                    % Final fitness values
                convergenceChart = mean(eta, 2);               % Average convergence curve

                % Save detailed convergence results
                path = ['./', algorithmDir, num2str(type), '_WT', num2str(wt), '_TN', num2str(tn), '.xls'];
                sheetName = ['WT_' num2str(wt)];
                xlswrite(path, eta, [sheetName, num2str(NA_type + 1), 'TN_', num2str(tn)]);

                % Save summary results into Excel
                tableOfBoxPlot = updateTableforBoxPlot(tableOfBoxPlot, boxPlotChart, pathOfBoxPlot, 0, NA_type + 1);
                tableOfConvergence = updateTableforConvergence(tableOfConvergence, convergenceChart, pathOfConvergence, 0, NA_type + 1);
            end
        end
    end
end
toc;

%% ================= Helper Functions ===================

% Initialize table format for boxplot summary (mean, std, each run)
function tableB = initTableforBoxPlot(problemNum, runNum)
    row = 3 + problemNum;
    col = 1 + runNum;
    tableB = cell(row, col);
    tableB{1,2} = 'mean';
    tableB{1,3} = 'std';
    for i = 1:problemNum
        tableB{i+1,1} = ['F' num2str(i)];
    end
end

% Initialize table format for convergence curve data
function tableC = initTableforConvergence(problemNum)
    col = problemNum;
    tableC = cell(1, col);
    for i = 1:col
        tableC{1,i} = ['F' num2str(i)];
    end
end

% Update boxplot table with new data
function tableB = updateTableforBoxPlot(tableB, data, path, rgo, problemIndex)
    [row, col] = size(data);
    rgo = repmat(rgo, 1, col);
    data = data - rgo;
    dataMean = mean(data, 2);
    dataStd = std(data, 0, 2);

    % Convert to cell and update table
    Dcol = ones(1, col);
    Drow = ones(1, row);
    dataMean = mat2cell(dataMean, Drow, [1]);
    dataStd = mat2cell(dataStd, Drow, [1]);
    tableB(problemIndex+1, 2) = dataMean;
    tableB(problemIndex+1, 3) = dataStd;
    dataCol = mat2cell(data, [1], Dcol);
    tableB(problemIndex+1, 4:3+col) = dataCol;

    xlswrite(path, tableB);
end

% Update convergence table with new curve
function tableC = updateTableforConvergence(tableC, data, path, rgo, problemIndex)
    [row, col] = size(data);
    rgo = repmat(rgo', 1, col);
    data = data - rgo;
    Drow = ones(1, row);
    dataCol = mat2cell(data, Drow, [1]);
    tableC(2:row+1, problemIndex) = dataCol;

    xlswrite(path, tableC);
end
