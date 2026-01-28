% quickContourPlots(resultsBayesopt, varNames)
plotMetabolicContours(resultsBayesopt, varNames)
% plotEnhancedContours(resultsBayesopt, varNames)
% plot3DMetabolicSurfaces(resultsBayesopt, varNames)
function plot3DMetabolicSurfaces(resultsBayesopt, varNames)
% Create 3D surface plots of metabolic cost for different parameter pairs
    
    X = table2array(resultsBayesopt.XTrace);
    y = resultsBayesopt.ObjectiveTrace;
    
    % Get parameter indices
    idx_tr = find(strcmp(varNames, 't_r'));
    idx_tp = find(strcmp(varNames, 't_p')); 
    idx_tf = find(strcmp(varNames, 't_f'));
    idx_Mp = find(strcmp(varNames, 'M_p'));
    
    figure('Position', [100, 100, 1400, 1000]);
    
    % Your requested combinations
    plot3DCombination(X, y, idx_Mp, idx_tp, 1, 'M_p', 't_p', 'Peak Torque vs Peak Time');
    plot3DCombination(X, y, idx_tp, idx_tr, 2, 't_p', 't_r', 'Peak Time vs Rise Time');
    plot3DCombination(X, y, idx_tr, idx_tf, 3, 't_r', 't_f', 'Rise Time vs Fall Time');
    plot3DCombination(X, y, idx_tf, idx_Mp, 4, 't_f', 'M_p', 'Fall Time vs Peak Torque');
    
    sgtitle('3D Metabolic Cost Landscape', 'FontSize', 16, 'FontWeight', 'bold');
end

function plot3DCombination(X, y, idx_x, idx_y, subplot_idx, x_name, y_name, title_str)
    subplot(2, 2, subplot_idx);
    
    x_vals = X(:, idx_x);
    y_vals = X(:, idx_y);
    z_vals = y;
    
    % Create grid for surface
    xi = linspace(min(x_vals), max(x_vals), 30);
    yi = linspace(min(y_vals), max(y_vals), 30);
    [XI, YI] = meshgrid(xi, yi);
    
    % Interpolate metabolic cost values
    ZI = griddata(x_vals, y_vals, z_vals, XI, YI, 'cubic');
    
    % Create 3D surface plot
    surf(XI, YI, ZI, 'EdgeColor', 'none', 'FaceAlpha', 0.8);
    hold on;
    
    % Add scatter points for actual evaluations
    scatter3(x_vals, y_vals, z_vals, 40, z_vals, 'filled', 'MarkerEdgeColor', 'k');
    
    % Add best point marker
    [best_cost, best_idx] = min(z_vals);
    best_x = x_vals(best_idx);
    best_y = y_vals(best_idx);
    plot3(best_x, best_y, best_cost, 'ro', 'MarkerSize', 10, 'LineWidth', 3, 'MarkerFaceColor', 'red');
    
    % Set isometric view
    axis equal;
    view(3); % 3D view
    rotate3d on; % Enable rotation
    
    % Labels and formatting
    xlabel(x_name);
    ylabel(y_name);
    zlabel('Metabolic Cost');
    title(title_str);
    colorbar;
    grid on;
    
    % Set nice colormap
    colormap(jet);
end

function plotEnhancedContours(resultsBayesopt, varNames)
% Enhanced contour plots with sensitivity analysis
    
    X = table2array(resultsBayesopt.XTrace);
    y = resultsBayesopt.ObjectiveTrace;
    
    % Calculate parameter sensitivities
    sensitivities = zeros(1, length(varNames));
    for i = 1:length(varNames)
        sensitivities(i) = abs(corr(X(:, i), y));
    end
    
    figure('Position', [100, 100, 1600, 1200]);
    
    % Parameter indices
    idx_tr = find(strcmp(varNames, 't_r'));
    idx_tp = find(strcmp(varNames, 't_p')); 
    idx_tf = find(strcmp(varNames, 't_f'));
    idx_Mp = find(strcmp(varNames, 'M_p'));
    
    %% Plot 1: Most sensitive parameter (tp) vs others
    subplot(2, 3, 1);
    createEnhancedContour(X(:, idx_tp), X(:, idx_Mp), y, sensitivities([idx_tp, idx_Mp]), ...
        'Peak Time (t_p)', 'Peak Torque (M_p)', 't_p vs M_p (Highest Sensitivity)');
    
    subplot(2, 3, 2);
    createEnhancedContour(X(:, idx_tp), X(:, idx_tr), y, sensitivities([idx_tp, idx_tr]), ...
        'Peak Time (t_p)', 'Rise Time (t_r)', 't_p vs t_r');
    
    subplot(2, 3, 3);
    createEnhancedContour(X(:, idx_tp), X(:, idx_tf), y, sensitivities([idx_tp, idx_tf]), ...
        'Peak Time (t_p)', 'Fall Time (t_f)', 't_p vs t_f');
    
    %% Plot 2: Other combinations
    subplot(2, 3, 4);
    createEnhancedContour(X(:, idx_Mp), X(:, idx_tr), y, sensitivities([idx_Mp, idx_tr]), ...
        'Peak Torque (M_p)', 'Rise Time (t_r)', 'M_p vs t_r');
    
    subplot(2, 3, 5);
    createEnhancedContour(X(:, idx_tr), X(:, idx_tf), y, sensitivities([idx_tr, idx_tf]), ...
        'Rise Time (t_r)', 'Fall Time (t_f)', 't_r vs t_f');
    
    subplot(2, 3, 6);
    createEnhancedContour(X(:, idx_tf), X(:, idx_Mp), y, sensitivities([idx_tf, idx_Mp]), ...
        'Fall Time (t_f)', 'Peak Torque (M_p)', 't_f vs M_p');
    
    sgtitle('Metabolic Cost Landscape with Parameter Sensitivity', 'FontSize', 16);
end

function createEnhancedContour(x_vals, y_vals, z_vals, sensitivities, x_label, y_label, plot_title)
    % Create interpolated surface
    xi = linspace(min(x_vals), max(x_vals), 50);
    yi = linspace(min(y_vals), max(y_vals), 50);
    [XI, YI] = meshgrid(xi, yi);
    ZI = griddata(x_vals, y_vals, z_vals, XI, YI, 'cubic');
    
    % Contour plot
    contourf(XI, YI, ZI, 10, 'LineColor', 'k');%,"FaceAlpha",0.25
    hold on;
    
    % Scatter plot with size indicating evaluation sequence
    scatter(x_vals, y_vals, 50, z_vals, 'filled', ...
            'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
    
    % Best point
    [best_cost, best_idx] = min(z_vals);
    best_x = x_vals(best_idx);
    best_y = y_vals(best_idx);
    plot(best_x, best_y, 'p', 'MarkerSize', 15, 'LineWidth', 3, ...
         'MarkerFaceColor', 'yellow', 'MarkerEdgeColor', 'red');
    
    % Add sensitivity info to title
    sensitivity_str = sprintf('(Sensitivity: X=%.2f, Y=%.2f)', sensitivities(1), sensitivities(2));
    title([plot_title, sensitivity_str], 'FontSize', 11);
    
    xlabel(x_label);
    ylabel(y_label);
    colorbar;
    grid on;
end

function plotMetabolicContours(resultsBayesopt, varNames)
% Create contour plots of metabolic cost for different parameter pairs
    
    fprintf('Creating metabolic cost contour plots...\n');
    
    % Extract data
    X = table2array(resultsBayesopt.XTrace);
    y = resultsBayesopt.ObjectiveTrace;
    
    % Get parameter indices
    idx_tr = find(strcmp(varNames, 't_r'));
    idx_tp = find(strcmp(varNames, 't_p')); 
    idx_tf = find(strcmp(varNames, 't_f'));
    idx_Mp = find(strcmp(varNames, 'M_p'));
    
    % Create figure
    figure(11);
    set(gcf,'Color','white')
    % figure('color','white'); %'Position', [100, 100, 1400, 1000],
    
    %% Subplot 1: Peak Torque (Mp) vs Peak Time (tp)
    subplot(2, 2, 1);
    createContourPlot(X(:, idx_Mp), X(:, idx_tp), y, ...
        'Peak Torque (M_p)', 'Peak Time (t_p)', ...
        'Metabolic Cost vs M_p and t_p');
    
    %% Subplot 2: Peak Time (tp) vs Rise Time (tr)
    subplot(2, 2, 2);
    createContourPlot(X(:, idx_tp), X(:, idx_tr), y, ...
        'Peak Time (t_p)', 'Rise Time (t_r)', ...
        'Metabolic Cost vs t_p and t_r');
    
    %% Subplot 3: Rise Time (tr) vs Fall Time (tf)
    subplot(2, 2, 3);
    createContourPlot(X(:, idx_tr), X(:, idx_tf), y, ...
        'Rise Time (t_r)', 'Fall Time (t_f)', ...
        'Metabolic Cost vs t_r and t_f');
    
    %% Subplot 4: Fall Time (tf) vs Peak Torque (Mp)
    subplot(2, 2, 4);
    createContourPlot(X(:, idx_tf), X(:, idx_Mp), y, ...
        'Fall Time (t_f)', 'Peak Torque (M_p)', ...
        'Metabolic Cost vs t_f and M_p');
    
    % Add overall title
    sgtitle('Metabolic Cost Landscape Across Parameter Pairs', 'FontSize', 16, 'FontWeight', 'bold');
end

function createContourPlot(x_vals, y_vals, z_vals, x_label, y_label, plot_title)
    % Create grid for contour plot
    xi = linspace(min(x_vals), max(x_vals), 50);
    yi = linspace(min(y_vals), max(y_vals), 50);
    [XI, YI] = meshgrid(xi, yi);
    
    % Interpolate metabolic cost values
    ZI = griddata(x_vals, y_vals, z_vals, XI, YI, 'cubic');
    
    % Create contour plot
    contourf(XI, YI, ZI, 10, 'LineColor', 'none');
    % contour(XI, YI, ZI, 10,'LineWidth',5);
    hold on;
    
    % Add scatter points for actual evaluations
    scatter(x_vals, y_vals, 40, z_vals, 'filled', 'MarkerEdgeColor', 'none');
    
    % Add colorbar and labels
    colorbar;
    xlabel(x_label, 'FontSize', 15);
    ylabel(y_label, 'FontSize', 15);
    title(plot_title, 'FontSize', 15);
    set(gca,'FontSize',15);
    
    % Add best point marker
    [best_cost, best_idx] = min(z_vals);
    best_x = x_vals(best_idx);
    best_y = y_vals(best_idx);
    plot(best_x, best_y, 'ro', 'MarkerSize', 10, 'LineWidth', 3, 'MarkerFaceColor', 'red');
    
    grid on;
end

function quickContourPlots(resultsBayesopt, varNames)
% Quick version with just the 4 requested subplots
    
    X = table2array(resultsBayesopt.XTrace);
    y = resultsBayesopt.ObjectiveTrace;
    
    % Get indices
    idx_tr = find(strcmp(varNames, 't_r'));
    idx_tp = find(strcmp(varNames, 't_p')); 
    idx_tf = find(strcmp(varNames, 't_f'));
    idx_Mp = find(strcmp(varNames, 'M_p'));
    
    figure('Position', [100, 100, 1200, 900]);
    
    % Your requested combinations
    plotCombination(X, y, idx_Mp, idx_tp, 1, 'M_p', 't_p', 'Peak Torque vs Peak Time');
    plotCombination(X, y, idx_tp, idx_tr, 2, 't_p', 't_r', 'Peak Time vs Rise Time');
    plotCombination(X, y, idx_tr, idx_tf, 3, 't_r', 't_f', 'Rise Time vs Fall Time');
    plotCombination(X, y, idx_tf, idx_Mp, 4, 't_f', 'M_p', 'Fall Time vs Peak Torque');
    
    sgtitle('Metabolic Cost Contour Analysis', 'FontSize', 14);
end

function plotCombination(X, y, idx_x, idx_y, subplot_idx, x_name, y_name, title_str)
    subplot(2, 2, subplot_idx);
    
    x_vals = X(:, idx_x);
    y_vals = X(:, idx_y);
    
    % Create contour
    xi = linspace(min(x_vals), max(x_vals), 40);
    yi = linspace(min(y_vals), max(y_vals), 40);
    [XI, YI] = meshgrid(xi, yi);
    ZI = griddata(x_vals, y_vals, y, XI, YI, 'cubic');
    
    contourf(XI, YI, ZI, 20, 'LineColor', 'none');
    hold on;
    
    % Plot evaluations
    scatter(x_vals, y_vals, 30, y, 'filled', 'MarkerEdgeColor', 'k');
    
    % Best point
    [~, best_idx] = min(y);
    plot(x_vals(best_idx), y_vals(best_idx), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
    
    xlabel(x_name);
    ylabel(y_name);
    title(title_str);
    colorbar;
    grid on;
end