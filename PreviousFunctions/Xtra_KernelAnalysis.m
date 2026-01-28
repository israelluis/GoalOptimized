% Run the complete analysis:
[names,models]=comprehensiveKernelAnalysis(resultsBayesopt, varNames);

% selected_model=models{3};

function [names,models]=comprehensiveKernelAnalysis(resultsBayesopt, varNames)
% Run all analyses together
    
    fprintf('=== COMPREHENSIVE KERNEL ANALYSIS ===\n\n');
    
    % 1. Quick check
    quickKernelCheck(resultsBayesopt);
    fprintf('\n');
    
    % 2. Landscape complexity
    analyzeLandscapeComplexity(resultsBayesopt, varNames);
    fprintf('\n');
    
    % 3. Full kernel evaluation
    [names,models]=kernelAnalysis(resultsBayesopt, varNames);
end
function analyzeLandscapeComplexity(resultsBayesopt, varNames)
% Understand if your landscape is smooth, rough, or has multiple optima
    
    X = table2array(resultsBayesopt.XTrace);
    y = resultsBayesopt.ObjectiveTrace;
    
    fprintf('=== LANDSCAPE COMPLEXITY ANALYSIS ===\n\n');
    
    % 1. Check for multiple optima
    [best_obj, best_idx] = min(y);
    local_optima = find(islocalmin(y));
    fprintf('Local optima detected: %d\n', length(local_optima));
    
    % 2. Calculate roughness (high frequency variations)
    y_smooth = smoothdata(y, 'movmedian', 5);
    roughness = mean(abs(y - y_smooth)) / std(y);
    fprintf('Landscape roughness: %.4f (0=smooth, 1=rough)\n', roughness);
    
    % 3. Parameter sensitivity
    param_sensitivity = zeros(1, size(X, 2));
    for i = 1:size(X, 2)
        correlation = abs(corr(X(:, i), y));
        param_sensitivity(i) = correlation;
        fprintf('  %s sensitivity: %.4f\n', varNames{i}, correlation);
    end
    
    % 4. Recommendation based on complexity
    if roughness > 0.3
        fprintf('\n🎯 LANDSCAPE: Rough → Use Matern32 or Exponential kernels\n');
    elseif length(local_optima) > 3
        fprintf('\n🎯 LANDSCAPE: Multi-modal → Use RationalQuadratic kernel\n');
    else
        fprintf('\n🎯 LANDSCAPE: Smooth → Use SquaredExponential kernel\n');
    end
end
function visualizeKernelFits(results, X_norm, y_norm, varNames)
    figure('Position', [100, 100, 1400, 800]);
    
    % Plot top 3 kernels
    for i = 1:6 %min(3, height(results))
        if ~isempty(results.GPR{i})
            subplot(2, 6, i);
            
            % Predictions vs actual
            y_pred = resubPredict(results.GPR{i});
            scatter(y_norm, y_pred, 50, 'filled');
            hold on;
            plot([min(y_norm), max(y_norm)], [min(y_norm), max(y_norm)], 'r--');
            xlabel('Actual (normalized)');
            ylabel('Predicted (normalized)');
            title(sprintf('%s\nMSE: %.4f', results.Kernel{i}, results.MSE(i)));
            grid on;
            
            % Residuals
            subplot(2, 6, i+6);
            residuals = y_norm - y_pred;
            scatter(y_pred, residuals, 50, 'filled');
            hold on;
            yline(0, 'r--', 'LineWidth', 2);
            xlabel('Predicted');
            ylabel('Residuals');
            title(sprintf('%s Residuals', results.Kernel{i}));
            grid on;
        end
    end
end
function [kernel_names, gpr_models] = kernelAnalysis(resultsBayesopt, varNames)
% Analyze which kernel function best describes your optimization landscape
    
    fprintf('=== KERNEL FUNCTION ANALYSIS ===\n\n');
    
    % Extract training data
    X_train = table2array(resultsBayesopt.XTrace);
    y_train = resultsBayesopt.ObjectiveTrace;
    
    % Normalize data for better kernel performance
    % X_mean = mean(X_train);
    % X_std = std(X_train);
    % X_norm = (X_train - X_mean) ./ X_std;
    % y_norm = (y_train - mean(y_train)) / std(y_train);
    
    % Define kernels to test
    kernels = {'ardmatern32', 'ardsquaredexponential', 'ardexponential', ...
               'matern32', 'squaredexponential', 'rationalquadratic'};
    
    % Preallocate results table with correct variable names
    kernel_names = {};
    mse_values = [];
    nll_values = [];
    gpr_models = {};
    
    % Test each kernel
    for i = 1:length(kernels)
        try
            fprintf('Testing %s... ', kernels{i});
            
            % Fit GPR model
            gpr = fitrgp(X_train, y_train, 'KernelFunction', kernels{i}, ...
                        'Standardize', false, 'Verbose', 0);
            
            % Calculate metrics
            y_pred = resubPredict(gpr);
            mse = mean((y_train - y_pred).^2);
            nll = loss(gpr, X_train, y_train); % Negative log likelihood
            
            % Store results
            kernel_names{end+1} = kernels{i};
            mse_values(end+1) = mse;
            nll_values(end+1) = nll;
            gpr_models{end+1} = gpr;
            
            fprintf('MSE: %.4f, NLL: %.4f\n', mse, nll);
            
        catch ME
            fprintf('FAILED: %s\n', ME.message);
            kernel_names{end+1} = kernels{i};
            mse_values(end+1) = NaN;
            nll_values(end+1) = NaN;
            gpr_models{end+1} = [];
        end
    end
    
    % Create results table
    results = table(kernel_names', mse_values', nll_values', gpr_models', ...
        'VariableNames', {'Kernel', 'MSE', 'NLL', 'GPR'});
    
    % Remove failed kernels and sort by best performance (lowest NLL)
    valid_results = results(~isnan(results.NLL), :);
    valid_results = sortrows(valid_results, 'NLL');
    
    % Display results
    fprintf('\n=== KERNEL RANKING (Best to Worst) ===\n');
    for i = 1:height(valid_results)
        fprintf('%d. %s: NLL=%.4f, MSE=%.4f\n', ...
            i, valid_results.Kernel{i}, valid_results.NLL(i), valid_results.MSE(i));
    end
    
    % Visualize predictions vs actual
    visualizeKernelFits(valid_results, X_train, y_train, varNames);
    
    % Return best kernel
    if ~isempty(valid_results)
        best_kernel = valid_results.Kernel{1};
        fprintf('\n🎯 RECOMMENDED KERNEL: %s\n', best_kernel);
    else
        fprintf('\n❌ No kernels successfully fitted\n');
    end
end
function quickKernelCheck(resultsBayesopt)
% Fast kernel comparison - just the essentials
    
    X = table2array(resultsBayesopt.XTrace);
    y = resultsBayesopt.ObjectiveTrace;
    
    kernels = {'ardmatern32', 'ardsquaredexponential', 'ardexponential', ...
               'matern32', 'squaredexponential', 'rationalquadratic'};
    
    fprintf('Quick Kernel Check:\n');
    for i = 1:length(kernels)
        try
            gpr = fitrgp(X, y, 'KernelFunction', kernels{i}, ...
                        'Standardize', false, 'Verbose', 0);
            nll = loss(gpr, X, y);
            fprintf('  %s: NLL = %.4f\n', kernels{i}, nll);
        catch
            fprintf('  %s: FAILED\n', kernels{i});
        end
    end
end