function [W,H,MActivation_recons_array,synMetrics]=synergyAnalysis(Results,synergy_list,plot_flag)
%% Compute synergies with proper multi-start NNMF
% Non-negative Matrix Factorization (most common)
rng(100); % for reproducibility

MActivation = Results.MActivation;

% NNMF parameters
n_replicates   = 100;    % Number of random starts - more = better but slower
algorithm_type = 'als';  % 'als' (alternating least squares) or 'mult' (multiplicative)
max_iterations = 2000;   % Maximum iterations per replicate
tolerance      = 1e-6;   % Convergence tolerance

fprintf('Running NNMF with %d replicates, algorithm: %s\n', n_replicates, algorithm_type);

nSynList     = length(synergy_list);
recons_error = zeros(1, nSynList);
W = cell(1, nSynList);
H = cell(1, nSynList);

% Set NNMF options
options = statset('MaxIter', max_iterations, 'Display', 'off', 'TolFun', tolerance);

for iSyn = 1:nSynList
    sSyn = synergy_list(iSyn);
    fprintf('  Testing %d synergies... ', sSyn);
    
    % Run NNMF with multiple replicates - MATLAB automatically picks the best
    [W_temp, H_temp, recons_error(1,iSyn)] = nnmf(...
        MActivation, ...
        sSyn, ...
        'replicates', n_replicates, ...  % Multiple random starts
        'algorithm', algorithm_type, ...  % Chosen algorithm
        'options', options);              % Convergence criteria
    
    for s = 1:sSyn
        col_max = max(W_temp(:, s));
        W_temp(:, s) = W_temp(:, s) / col_max;
        H_temp(s, :) = H_temp(s, :) * col_max;
    end
    
    W{1,iSyn} = W_temp;
    H{1,iSyn} = H_temp;

    fprintf('Best reconstruction error: %.4f\n', recons_error(1,iSyn));
end

%% Reconstructed synergies
full_data_length = length(MActivation);
nMuscles = length(Results.MuscleNames);

MActivation_recons_array = zeros(nSynList, nMuscles, full_data_length);
VAF  = zeros(nSynList,1);
RMSE = zeros(nSynList,1);
RMSE_per_muscle=zeros(nSynList,nMuscles);
R_per_muscle   =zeros(nSynList,nMuscles);

% Calculate Variance Accounted For (VAF)
total_variance = sum(var(MActivation, 0, 2));

for iSyn = 1:nSynList
    MActivation_reconstructed = W{1,iSyn} * H{1,iSyn};
    MActivation_recons_array(iSyn, :, :) = MActivation_reconstructed;

    % metrics of interests
    VAF(iSyn,1) = 1 - (sum(var(MActivation - MActivation_reconstructed, 0, 2)) / total_variance);
    RMSE(iSyn,1) = sqrt(mean((MActivation - MActivation_reconstructed).^2, 'all'));
    RMSE_per_muscle(iSyn,:) = sqrt(mean((MActivation - MActivation_reconstructed).^2, 2));
    for m = 1:nMuscles
        R_per_muscle(iSyn,m) = corr(MActivation(m,:)', MActivation_reconstructed(m,:)');
    end    
end
synMetrics.N=synergy_list;
synMetrics.VAF=VAF;
synMetrics.RMSE=RMSE;
synMetrics.RMSE_per_muscle=RMSE_per_muscle;
synMetrics.R_per_muscle=R_per_muscle;

% Display results summary
fprintf('\n=== NNMF Results Summary ===\n');
fprintf('Synergies  Recon.Error  VAF\n');
fprintf('------------------------------\n');
for iSyn = 1:nSynList
    fprintf('    %d        %.4f     %.3f\n', ...
        synergy_list(iSyn), recons_error(1,iSyn), VAF(iSyn,1));
end

%% Plotting section (unchanged except for variable name clarity)
to_plot_recoSyn = plot_flag(1);
to_plot_optiSyn = plot_flag(2);
to_plot_evalSyn = plot_flag(3);

syn_color_list = {'#000000', '#2C5F7E', '#4A8BB5', '#6AB0DE', '#8AC6F0'};

if to_plot_recoSyn == 1
    figure('WindowState', 'maximized', 'Name', 'Muscle Activation Reconstruction');    
    for j = 1:nMuscles
        subplot(5, 8, j);
        pl(1) = plot(MActivation(j, :), 'Color', syn_color_list{1}, ...
            'LineWidth', 2, 'DisplayName', 'Original');
        hold on
        
        for i = 1:length(synergy_list)
            pl(i+1) = plot(squeeze(MActivation_recons_array(i, j, :)), ...
                'Color', syn_color_list{i+1}, 'LineWidth', 1.5, ...
                'DisplayName', sprintf('%d syn', synergy_list(i)));
        end
        
        legend(pl, 'NumColumns', 2, 'Box', 'off', 'FontSize', 8);
        ylim([0 1]);
        title(Results.MuscleNames{j}, 'FontSize', 10);
        grid on;
    end
end

if to_plot_optiSyn == 1
    figure('WindowState', 'maximized', 'Name', 'NNMF Optimization Plot');
    
    % Choose k where VAF > 90% or elbow point
    optimal_Syn = find(VAF > 0.9, 1);
    
    % Create scatter plot with markers
    scatter(recons_error, VAF, 200, 'filled');
    hold on;
    
    % Add labels for each synergy count
    for i = 1:length(synergy_list)
        text(recons_error(i) + 0.001, VAF(i) + 0.005, ...
            sprintf('%d syn', synergy_list(i)), 'FontSize', 12, 'FontWeight', 'bold');
        
        if i == optimal_Syn
            text(recons_error(i) + 0.001, VAF(i) - 0.015, ...
                'Optimal', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'red');
        end
    end
    
    % Add reference lines and labels
    yline(0.9, '--k', 'LineWidth', 1.5, 'Alpha', 0.7);
    text(min(recons_error), 0.91, 'VAF = 0.9', 'FontSize', 11);
    
    % Plot formatting
    xlabel('Reconstruction Error', 'FontSize', 14);
    ylabel('Variance Accounted For (VAF)', 'FontSize', 14);
    title('NNMF Performance vs. Number of Synergies', 'FontSize', 16);
    grid on;
    box on;
end

if to_plot_evalSyn == 1
    for iSyn=1:length(synergy_list)
        sSyn=synergy_list(iSyn);
        W_opt = W{1, iSyn};
        H_opt = H{1, iSyn};
        figure('WindowState', 'maximized', 'Name', sprintf('Synergy Details (%d Synergies)', sSyn));

        for i = 1:sSyn
            % Muscle weights (bar plot)
            subplot(sSyn, 2, 2*i-1);
            bar(W_opt(:, i), 'FaceColor', syn_color_list{min(i+1, length(syn_color_list))});
            ylim([0 1]);
            title(sprintf('Synergy %d Weights', i), 'FontSize', 12);
            ylabel('Normalized Weight');
            grid on;

            % Synergy activation (line plot)
            subplot(sSyn, 2, 2*i);
            plot(H_opt(i, :), 'Color', syn_color_list{min(i+1, length(syn_color_list))}, 'LineWidth', 2);
            xlim([0 full_data_length]);
            title(sprintf('Synergy %d Activation', i), 'FontSize', 12);
            xlabel('Time (samples)');
            ylabel('Activation');
            grid on;
        end

        % Add overall title
        sgtitle(sprintf('Muscle Synergy Analysis - %d Synergies (VAF = %.3f)', sSyn, VAF(iSyn)), 'FontSize', 14);
    end
end