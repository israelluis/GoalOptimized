clc;
colorBase = {'#A0A0A0', '#7CFC00', '#FF6347', '#87CEEB'};
colorOpti = {'k', '#006400', '#8B0000', '#00008B'};

sub_list=[2];
nSubs=length(sub_list);

sVar='MActivation'; % lMtildeopt MExcitation MActivation
fig=figure(1); clf; set(gcf,'color','w'); fig.WindowState = 'maximized';
for iSub=1:nSubs
    sSub=sub_list(iSub);
    Results=MRS_Je0{sSub,1}.Results;
    Misc   =MRS_Je0{sSub,1}.Misc;
    extra_frames=Misc.extra_frames;

    fSel=1+extra_frames:size(Results.(sVar).genericMRS,2)-1-extra_frames;
    gait_cycle=linspace(0,100,length(fSel));

    SSyn_list=[4];
    nSyns=length(SSyn_list);

    indexDOF=1;

    LineStyle_sel='-'; %-.
    for iMus=1:40
        subplot(6,8,iMus); hold on

        % no synergies
        Results_Baseline=MRS_Je{sSub,1};
        MActivation_N=Results_Baseline.Results.(sVar).genericMRS(iMus,fSel);
        plot(gait_cycle,MActivation_N,'Color',colorBase{1},'LineWidth',3,'LineStyle','-'); hold on %'#71797E'

        sSyn=1;
        Results_Bilevel =MRS_JeSD{sSub,sSyn,indexDOF};
        MActivation_OPT=Results_Bilevel.Results.(sVar).genericMRS(iMus,fSel);
        plot(gait_cycle,MActivation_OPT,'Color',colorOpti{1},'LineWidth',3,'LineStyle',LineStyle_sel)

        % synergies
        for iSyn=1:nSyns
            sSyn=SSyn_list(iSyn);
            Results_Baseline=MRS_JeS{sSub,sSyn-1}; %[4 5 6]
            MActivation_N=Results_Baseline.Results.(sVar).genericMRS(iMus,fSel);
            plot(gait_cycle,MActivation_N,'Color',colorBase{iSyn+1},'LineWidth',3,'LineStyle','-'); hold on %'#71797E'

            Results_Bilevel =MRS_JeSD{sSub,sSyn,indexDOF};
            MActivation_OPT=Results_Bilevel.Results.(sVar).genericMRS(iMus,fSel);
            plot(gait_cycle,MActivation_OPT,'Color',colorOpti{iSyn+1},'LineWidth',3,'LineStyle',LineStyle_sel)
        end
         set(gca,'FontSize',13);

         if strcmp(sVar,'lMtildeopt')
             axis([0 100 0.5 1.5])
         else
            % axis([0 100 0 0.9])
         end
         title(Results.MuscleNames{iMus});
    end
    for iDOF=1:5
        subplot(6,8,40+iDOF); hold on
        
        Results_Baseline=MRS_Je{sSub,1};
        DofNames_Input=Results_Baseline.Misc.DofNames_Input;
        RActivation=Results_Baseline.Results.RActivation.genericMRS(iDOF,fSel);

        plot(gait_cycle,RActivation);

        sSyn=1;
        Results_Bilevel =MRS_JeSD{sSub,sSyn,indexDOF};
        RActivation_OPT=Results_Bilevel.Results.RActivation.genericMRS(iDOF,fSel);
        plot(gait_cycle,RActivation_OPT,'Color',colorOpti{1},'LineWidth',3,'LineStyle',LineStyle_sel)

        for iSyn=1:nSyns
            sSyn=SSyn_list(iSyn);
            Results_Baseline=MRS_JeS{sSub,sSyn-1}; %[4 5 6]
            MActivation_N=Results_Baseline.Results.RActivation.genericMRS(iDOF,fSel);
            plot(gait_cycle,MActivation_N,'Color',colorBase{iSyn+1},'LineWidth',3,'LineStyle','-'); hold on %'#71797E'

            Results_Bilevel =MRS_JeSD{sSub,sSyn,indexDOF};
            MActivation_OPT=Results_Bilevel.Results.RActivation.genericMRS(iDOF,fSel);
            plot(gait_cycle,MActivation_OPT,'Color',colorOpti{iSyn+1},'LineWidth',3,'LineStyle',LineStyle_sel)
        end

        % axis([0 100 -10 10])
        title([DofNames_Input{iDOF}])
    end
    % sgtitle(['sub' num2str(sSub) ' target joint: ' DOFNames_list{indexDOF}{1}])
end
%%
clc;
fig = figure(2); clf; set(gcf, 'color', 'w');  fig.WindowState = 'maximized';

% Colors for plots
var1 = {'#A0A0A0', '#7CFC00', '#FF6347', '#87CEEB'}; % Light colors
var2 = {'#2F4F4F', '#006400', '#8B0000', '#00008B'}; % Dark colors

% Synergy counts (4, 5, 6)
syn_list = 4:6;
sSub = 3; % Single subject

% Loop through synergy counts
for iSyn = 1:length(syn_list)
    sSyn = syn_list(iSyn); % Current synergy count (4, 5, or 6)
    
    % Get results for this subject and synergy number
    Results_Baseline = MRS_JeS{sSub, iSyn}; % Assuming indexing starts at 1
    H = Results_Baseline.Results.SynergyControl.H;
    W = Results_Baseline.Results.SynergyControl.W;
    MuscleNames=Results_Baseline.Results.MuscleNames;

    % Get results for this subject and synergy number
    iDev=1;
    Results_Device = MRS_JeSD{sSub, iSyn+1, iDev}; % Assuming indexing starts at 1
    H_dev = Results_Device.Results.SynergyControl.H;
    W_dev = Results_Device.Results.SynergyControl.W;

    nMuscles = size(W, 1);
    nSynergies = size(W, 2);
    
    % Calculate row indices for this synergy count
    weightRow = (iSyn-1)*2 + 1;
    activationRow = weightRow + 1;
    
    % ==================== ROW 1: WEIGHTS (W matrix) ====================
    % Plot weights for all synergies
    subplot(6, nSynergies, sub2ind([nSynergies, 6], 1:nSynergies, ones(1,nSynergies)*weightRow));
    
    % Create bar plot for each muscle weight across synergies
    bar(W', 'grouped');
    title(['Weights for ' num2str(sSyn) ' Synergies']);
    xlabel('Synergy #');
    ylabel('Weight');
    ylim([0 1]);
    grid on;
    
    % Set x-tick labels
    xticks(1:nSynergies);
    
    % Apply colors to bars (one color per muscle)
    hold on;
    for iMuscle = 1:min(nMuscles, length(var1))
        % Find bars for this muscle
        hBars = findobj(gca, 'Type', 'bar', 'SeriesIndex', iMuscle);
        if ~isempty(hBars)
            set(hBars, 'FaceColor', var1{mod(iMuscle-1, length(var1))+1}, ...
                       'EdgeColor', var2{mod(iMuscle-1, length(var2))+1});
        end
    end
    hold off;
    
    % ==================== ROW 2: ACTIVATIONS (H matrix) ====================
    % Plot activations for all synergies
    for synNum = 1:nSynergies
        subplot(6, nSynergies, sub2ind([nSynergies, 6], synNum, activationRow));
        
        plot(H(synNum, :), 'LineWidth', 2, 'Color', var2{mod(synNum-1, length(var2))+1});
        title(['Synergy ' num2str(synNum) ' Activation']);
        xlabel('Time sample');
        ylabel('Activation');
        grid on;
        
        % Adjust y-axis for better visualization
        yLimits = [min(H(synNum, :)), max(H(synNum, :))];
        if diff(yLimits) > 0
            ylim(yLimits);
        end
    end
end

% Add overall title
sgtitle(['Muscle Synergies Analysis - Subject ' num2str(sSub)], 'FontSize', 16, 'FontWeight', 'bold');

% Adjust layout
set(gcf, 'Position', get(0, 'Screensize')); % Maximize figure
%%
DOFLabels= {'ANKLE PLANTARFLEXION' 'KNEE EXTENSION' 'HIP FLEXION' 'HIP ABDUCTION'};

clc;

% Colors for plots
var1 = {'#A0A0A0', '#7CFC00', '#FF6347', '#87CEEB'}; % Light colors
var2 = {'#2F4F4F', '#006400', '#8B0000', '#00008B'}; % Dark colors

% Synergy counts (4, 5, 6)
syn_list = [4 5 6];

sub_list=[5];

for iSub=1:length(sub_list)
    sSub = sub_list(iSub);
    fig = figure(iSub); clf; set(gcf, 'color', 'w');  fig.WindowState = 'maximized';

    % Loop through synergy counts
    for iSyn = 1:length(syn_list)
        sSyn = syn_list(iSyn); % Current synergy count (4, 5, or 6)

        % Get baseline results
        Results_Baseline = MRS_JeS{sSub, iSyn};
        H = Results_Baseline.Results.SynergyControl.H;
        W = Results_Baseline.Results.SynergyControl.W;
        MuscleNames = Results_Baseline.Results.MuscleNames;

        % Get device results
        iDev = 1;
        Results_Device = MRS_JeSD{sSub, iSyn+1, iDev};
        H_dev = Results_Device.Results.SynergyControl.H;
        W_dev = Results_Device.Results.SynergyControl.W;

        nMuscles = size(W, 1);
        nSynergies = size(W, 2);

        % Calculate row indices for this synergy count
        weightRow = (iSyn-1)*2 + 1;
        activationRow = weightRow + 1;

        % ==================== ROW 1: BASELINE WEIGHTS ONLY ====================
        % Plot each synergy's baseline weights in a separate subplot
        for synNum = 1:nSynergies
            % Calculate subplot position for weights
            subplotPos = (weightRow-1)*nSynergies + synNum;
            subplot(6, nSynergies, subplotPos);

            % Plot baseline weights only
            bar(1:nMuscles, W(:, synNum), 'FaceColor', var1{mod(synNum-1, length(var1))+1}, ...
                'EdgeColor', var2{mod(synNum-1, length(var2))+1}, 'BarWidth', 0.8);

            title(['Synergy ' num2str(synNum) ' Weights (Baseline)']);
            xlabel('Muscle');
            ylabel('Weight');
            ylim([0 3]);
            grid on;

            % Set x-ticks
            if nMuscles <= 15
                xticks(1:nMuscles);
                if synNum == 1 && iSyn == 1
                    % Show muscle names for first plot if not too many
                    if nMuscles <= 10
                        xticklabels(MuscleNames);
                        xtickangle(45);
                    else
                        xticklabels(1:nMuscles);
                    end
                else
                    xticklabels(1:nMuscles);
                end
            else
                xticks(1:5:nMuscles);
                xticklabels(1:5:nMuscles);
            end
        end

        % ==================== ROW 2: ACTIVATIONS - BASELINE vs DEVICE OVERLAY ====================
        % Plot overlapped activations for all synergies
        for synNum = 1:nSynergies
            % Calculate subplot position for activations
            subplotPos = (activationRow-1)*nSynergies + synNum;
            subplot(6, nSynergies, subplotPos);

            % Check if H and H_dev have same length
            nSamples = min(size(H, 2), size(H_dev, 2));

            % Plot baseline activation
            plot(1:nSamples, H(synNum, 1:nSamples), 'LineWidth', 2.5, ...
                'Color', var2{mod(synNum-1, length(var2))+1}, 'DisplayName', 'Baseline');
            hold on;

            % Plot device activation
            plot(1:nSamples, H_dev(synNum, 1:nSamples), 'LineWidth', 2, 'LineStyle', '-.', ...
                'Color', var2{mod(synNum-1, length(var2))+1}, 'DisplayName', 'Device');

            title(['Synergy ' num2str(synNum) ' Activation']);
            xlabel('Time sample');
            ylabel('Activation');
            grid on;
            ylim([0 0.4])

            hold off;
        end
    end

    % Add overall title
    sgtitle(['Muscle Synergies Analysis - Subject: N' num2str(sSub) ' - Conditions: Baseline vs Assisted (' DOFLabels{iDev} ')'], ...
        'FontSize', 16, 'FontWeight', 'bold');

    % Adjust layout for better spacing
    % set(gcf, 'Position', get(0, 'Screensize')); % Maximize figure
end
%%

clc;
fig = figure(2); clf; set(gcf, 'color', 'w');  fig.WindowState = 'maximized';

% Colors for plots
var1 = {'#A0A0A0', '#7CFC00', '#FF6347', '#87CEEB'}; % Light colors
var2 = {'#2F4F4F', '#006400', '#8B0000', '#00008B'}; % Dark colors

% Synergy counts (4, 5, 6)
syn_list = 4:4;
sSub = 2; % Single subject

% Loop through synergy counts
for iSyn = 1:length(syn_list)
    sSyn = syn_list(iSyn); % Current synergy count (4, 5, or 6)
    
    % Get baseline results
    Results_Baseline = MRS_JeS{sSub, iSyn}; 
    H = Results_Baseline.Results.SynergyControl.H;
    W = Results_Baseline.Results.SynergyControl.W;
    MuscleNames = Results_Baseline.Results.MuscleNames;
    
    % Get device results
    iDev = 1;
    Results_Device = MRS_JeSD{sSub, iSyn+1, iDev};
    H_dev = Results_Device.Results.SynergyControl.H;
    W_dev = Results_Device.Results.SynergyControl.W;

    nMuscles = size(W, 1);
    nSynergies = size(W, 2);
    
    % Calculate row indices for this synergy count
    weightRow = (iSyn-1) + 1;
    activationRow = weightRow + 1;
    
    % ==================== ROW 1: WEIGHTS - BASELINE vs DEVICE ====================
    % Plot each synergy's weights in a separate subplot
    for synNum = 1:nSynergies
        % Calculate subplot position for weights
        subplotPos = (weightRow-1)*(nSynergies) + (2*synNum-1);
        subplot(6, nSynergies, subplotPos);
        
        % Create x positions for baseline and device bars
        x = 1:nMuscles;
        
        % Plot baseline weights (solid bars)
        bar(x, W(:, synNum), 'FaceColor', var1{mod(synNum-1, length(var1))+1}, ...
            'EdgeColor', var2{mod(synNum-1, length(var2))+1}, ...
            'BarWidth', 0.4, 'DisplayName', 'Baseline');
        hold on;
        
        % Plot device weights (semi-transparent bars offset to the right)
        bar(x + 0.4, W_dev(:, synNum), 'FaceColor', var1{mod(synNum-1, length(var1))+1}, ...
            'EdgeColor', var2{mod(synNum-1, length(var2))+1}, ...
            'FaceAlpha', 0.5, 'BarWidth', 0.4, 'DisplayName', 'Device');
        
        title(['Synergy ' num2str(synNum) ' Weights']);
        xlabel('Muscle');
        ylabel('Weight');
        ylim([0 1]);
        grid on;
        
        % Set x-ticks to show muscle numbers (or names for first few)
        if synNum == 1 && iSyn == 1
            % For first synergy of first set, show muscle names if there are few
            if nMuscles <= 10
                xticks(1:nMuscles);
                xticklabels(MuscleNames(1:min(nMuscles, 10)));
                xtickangle(45);
            else
                xticks(1:5:nMuscles);
                xticklabels(1:5:nMuscles);
            end
        else
            xticks(1:5:nMuscles);
            xticklabels(1:5:nMuscles);
        end
        
        % Add legend for first plot only to avoid clutter
        if synNum == 1 && iSyn == 1
            legend('Location', 'best');
        end
        
        hold off;
    end
    
    % ==================== ROW 2: ACTIVATIONS - BASELINE vs DEVICE OVERLAY ====================
    % Plot activations for all synergies with overlay
    % for synNum = 1:nSynergies
    %     % Calculate subplot position for activations
    %     subplotPos = (activationRow-1)*(nSynergies) + (2*synNum-1);
    %     subplot(6, nSynergies, subplotPos);
    % 
    %     % Check if H and H_dev have same length
    %     nSamples = min(size(H, 2), size(H_dev, 2));
    % 
    %     % Plot baseline activation
    %     plot(1:nSamples, H(synNum, 1:nSamples), 'LineWidth', 2, ...
    %         'Color', var2{mod(synNum-1, length(var2))+1}, 'DisplayName', 'Baseline');
    %     hold on;
    % 
    %     % Plot device activation
    %     plot(1:nSamples, H_dev(synNum, 1:nSamples), 'LineWidth', 2, 'LineStyle', '--', ...
    %         'Color', var2{mod(synNum-1, length(var2))+1}, 'DisplayName', 'Device');
    % 
    %     title(['Synergy ' num2str(synNum) ' Activation']);
    %     xlabel('Time sample');
    %     ylabel('Activation');
    %     grid on;
    % 
    %     % Add legend for first plot only
    %     if synNum == 1 && iSyn == 1
    %         legend('Location', 'best');
    %     end
    % 
    %     % Adjust y-axis for better visualization
    %     allData = [H(synNum, 1:nSamples), H_dev(synNum, 1:nSamples)];
    %     yLimits = [min(allData), max(allData)];
    %     if diff(yLimits) > 0
    %         ylim(yLimits * 1.1); % Add 10% margin
    %     end
    % 
    %     hold off;
    % end
end

% Add overall title
sgtitle(['Muscle Synergies Analysis - Subject ' num2str(sSub) ' (Baseline vs Device)'], ...
    'FontSize', 16, 'FontWeight', 'bold');

%%


clc;
fig = figure(2); clf; set(gcf, 'color', 'w');  fig.WindowState = 'maximized';

% Colors for plots
var1 = {'#A0A0A0', '#7CFC00', '#FF6347', '#87CEEB'}; % Light colors
var2 = {'#2F4F4F', '#006400', '#8B0000', '#00008B'}; % Dark colors

% Synergy counts (4, 5, 6)
syn_list = 4:6;
sSub = 2; % Single subject

% Loop through synergy counts
for iSyn = 1:length(syn_list)
    sSyn = syn_list(iSyn); % Current synergy count (4, 5, or 6)
    
    % Get baseline results
    Results_Baseline = MRS_JeS{sSub, iSyn}; 
    H = Results_Baseline.Results.SynergyControl.H;
    W = Results_Baseline.Results.SynergyControl.W;
    MuscleNames = Results_Baseline.Results.MuscleNames;
    
    % Get device results
    iDev = 1;
    Results_Device = MRS_JeSD{sSub, iSyn+1, iDev};
    H_dev = Results_Device.Results.SynergyControl.H;
    W_dev = Results_Device.Results.SynergyControl.W;

    nMuscles = size(W, 1);
    nSynergies = size(W, 2);
    
    % Calculate row indices for this synergy count
    weightRow = (iSyn-1)*2 + 1;
    activationRow = weightRow + 1;
    
    % ==================== ROW 1: WEIGHTS - BASELINE vs DEVICE ====================
    % Create subplot for weights comparison
    subplot(6, 2*nSynergies, sub2ind([2*nSynergies, 6], 1:2:2*nSynergies, ones(1,nSynergies)*weightRow));
    
    % Plot baseline weights
    xPos = 1:0.25:1+(nSynergies-1)*0.25;
    for synNum = 1:nSynergies
        hBar = bar(xPos(synNum), W(:, synNum)', 'BarWidth', 0.2);
        set(hBar, 'FaceColor', var1{mod(synNum-1, length(var1))+1}, ...
                  'EdgeColor', var2{mod(synNum-1, length(var2))+1});
        hold on;
    end
    
    % Plot device weights (shifted to the right)
    xPosDev = 1.1:0.25:1.1+(nSynergies-1)*0.25;
    for synNum = 1:nSynergies
        hBar = bar(xPosDev(synNum), W_dev(:, synNum)', 'BarWidth', 0.2, 'FaceAlpha', 0.5);
        set(hBar, 'FaceColor', var1{mod(synNum-1, length(var1))+1}, ...
                  'EdgeColor', var2{mod(synNum-1, length(var2))+1}, 'LineStyle', '--');
    end
    
    title(['Weights for ' num2str(sSyn) ' Synergies']);
    xlabel('Synergy #');
    ylabel('Weight');
    ylim([0 1]);
    grid on;
    
    % Set x-tick labels and position them between baseline and device bars
    xTicks = (xPos + xPosDev) / 2;
    xticks(xTicks);
    xticklabels(arrayfun(@(x) ['Syn' num2str(x)], 1:nSynergies, 'UniformOutput', false));
    
    % Add legend
    legend('Baseline', 'Device', 'Location', 'best');
    hold off;
    
    % ==================== ROW 2: ACTIVATIONS - BASELINE vs DEVICE OVERLAY ====================
    % Plot activations for all synergies with overlay
    for synNum = 1:nSynergies
        subplot(6, 2*nSynergies, sub2ind([2*nSynergies, 6], 2*synNum-1:2*synNum, activationRow));
        
        % Check if H and H_dev have same length
        nSamples = min(size(H, 2), size(H_dev, 2));
        
        % Plot baseline activation
        plot(1:nSamples, H(synNum, 1:nSamples), 'LineWidth', 2, ...
            'Color', var2{mod(synNum-1, length(var2))+1}, 'DisplayName', 'Baseline');
        hold on;
        
        % Plot device activation
        plot(1:nSamples, H_dev(synNum, 1:nSamples), 'LineWidth', 2, 'LineStyle', '--', ...
            'Color', var2{mod(synNum-1, length(var2))+1}, 'DisplayName', 'Device');
        
        title(['Synergy ' num2str(synNum) ' Activation']);
        xlabel('Time sample');
        ylabel('Activation');
        grid on;
        
        % Add legend
        legend('Location', 'best');
        
        % Adjust y-axis for better visualization
        allData = [H(synNum, 1:nSamples), H_dev(synNum, 1:nSamples)];
        yLimits = [min(allData), max(allData)];
        if diff(yLimits) > 0
            ylim(yLimits * 1.1); % Add 10% margin
        end
        
        hold off;
    end
end
