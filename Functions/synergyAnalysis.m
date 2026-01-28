function [W,H,MActivation_reconstructed_array,synMetrics]=synergyAnalysis(Results,synergy_list,plot_flag)
%% Compute synergies with proper multi-start NNMF
% Non-negative Matrix Factorization (most common)
rng(100); % for reproducibility

MActivation = Results.MActivation.genericMRS;

% NNMF parameters - ADJUST THESE BASED ON YOUR NEEDS
n_replicates   = 100;      % Number of random starts - more = better but slower
algorithm_type = 'als'; % 'als' (alternating least squares) or 'mult' (multiplicative)
max_iterations = 2000;  % Maximum iterations per replicate
tolerance      = 1e-6;       % Convergence tolerance

fprintf('Running NNMF with %d replicates, algorithm: %s\n', n_replicates, algorithm_type);

nSynergies = length(synergy_list);
reconstruction_error = zeros(1, nSynergies);
W = cell(1, nSynergies);
H = cell(1, nSynergies);

% Set NNMF options
options = statset('MaxIter', max_iterations, 'Display', 'off', 'TolFun', tolerance);

for iSyn = 1:nSynergies
    fprintf('  Testing %d synergies... ', synergy_list(iSyn));
    
    % Run NNMF with multiple replicates - MATLAB automatically picks the best
    [W{1,iSyn}, H{1,iSyn}, reconstruction_error(1,iSyn)] = nnmf(...
        MActivation, ...
        synergy_list(iSyn), ...
        'replicates', n_replicates, ...  % Multiple random starts
        'algorithm', algorithm_type, ...  % Chosen algorithm
        'options', options);              % Convergence criteria
    
    fprintf('Best reconstruction error: %.4f\n', reconstruction_error(1,iSyn));
end

%% Reconstructed synergies
full_data_length = length(Results.MActivation.genericMRS);
nMuscles = length(Results.MuscleNames);

MActivation_reconstructed_array = zeros(nSynergies, nMuscles, full_data_length);
VAF = zeros(nSynergies,1);
RMSE = zeros(nSynergies,1);
RMSE_per_muscle=zeros(nSynergies,nMuscles);
R_per_muscle   =zeros(nSynergies,nMuscles);

% Calculate Variance Accounted For (VAF)
total_variance = sum(var(MActivation, 0, 2));

for iSyn = 1:nSynergies
    MActivation_reconstructed = W{1,iSyn} * H{1,iSyn};
    MActivation_reconstructed_array(iSyn, :, :) = MActivation_reconstructed;

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
for iSyn = 1:nSynergies
    fprintf('    %d        %.4f     %.3f\n', ...
        synergy_list(iSyn), reconstruction_error(1,iSyn), VAF(iSyn,1));
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
            pl(i+1) = plot(squeeze(MActivation_reconstructed_array(i, j, :)), ...
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
    scatter(reconstruction_error, VAF, 200, 'filled');
    hold on;
    
    % Add labels for each synergy count
    for i = 1:length(synergy_list)
        text(reconstruction_error(i) + 0.001, VAF(i) + 0.005, ...
            sprintf('%d syn', synergy_list(i)), 'FontSize', 12, 'FontWeight', 'bold');
        
        if i == optimal_Syn
            text(reconstruction_error(i) + 0.001, VAF(i) - 0.015, ...
                'Optimal', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'red');
        end
    end
    
    % Add reference lines and labels
    yline(0.9, '--k', 'LineWidth', 1.5, 'Alpha', 0.7);
    text(min(reconstruction_error), 0.91, 'VAF = 0.9', 'FontSize', 11);
    
    % Plot formatting
    xlabel('Reconstruction Error', 'FontSize', 14);
    ylabel('Variance Accounted For (VAF)', 'FontSize', 14);
    title('NNMF Performance vs. Number of Synergies', 'FontSize', 16);
    grid on;
    box on;
end

% Fixed synergy count for detailed visualization (you might want to make this adaptive)
sSyn = 5;  % You could replace this with: sSyn = synergy_list(optimal_Syn);
oSyn = find(synergy_list == sSyn, 1);

if ~isempty(oSyn)
    W_opt = W{1, oSyn};
    H_opt = H{1, oSyn};
    
    if to_plot_evalSyn == 1
        figure('WindowState', 'maximized', 'Name', sprintf('Synergy Details (%d Synergies)', sSyn));    
        
        % Normalize W for better visualization (optional)
        W_normalized = W_opt ./ max(W_opt);  % Normalize to [0, 1]
        
        for i = 1:sSyn
            % Muscle weights (bar plot)
            subplot(sSyn, 2, 2*i-1);
            bar(W_normalized(:, i), 'FaceColor', syn_color_list{min(i+1, length(syn_color_list))});
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
        sgtitle(sprintf('Muscle Synergy Analysis - %d Synergies (VAF = %.3f)', sSyn, VAF(oSyn)), 'FontSize', 14);
    end
else
    warning('Selected synergy count (%d) not in synergy_list. Skipping detailed plots.', sSyn);
end
end