folderPath        = 'C:\Users\movea\Dropbox\PostDoc Work\Simulation\Code\ProjectResults\Codesign\Pilot\sub1\v2_t1\iters200';  % Replace with your folder path
currentFolder=pwd;
addpath(genpath(currentFolder));

muscleGroup = {'AP', 'KE' 'HF' 'HB'};
actuatorType = {'Q', 'A'};
arrangement = [1 2];

analysis_list = generate_analysis_list(muscleGroup, actuatorType, arrangement);
% analysis_list={'KE(Q)' 'KE(A)' 'AP(Q)' 'AP(A)' 'AP(Q)&KE(Q)' 'AP(A)&KE(Q)' 'AP(Q)&KE(A)' 'AP(A)&KE(A)'};
NAnalysis=length(analysis_list);
%% Metabolic cost: eDot_MCLU24
set_optimal_criteria= '_VERX_eDot_MCLU24'; 

set_type_result     = 'MRSOptimal';       
fileNames=get_files_from_pattern(folderPath,[set_type_result set_optimal_criteria])';
MRS_eDot=cell(NAnalysis,1);
for i=1:NAnalysis
    index=contains(fileNames,['_' analysis_list{i} '_']);
    MRS_eDot(i)={load(fullfile(folderPath,fileNames{index}))};
end

set_type_result   = 'Bayesopt';       
fileNames_bayes_eDot=get_files_from_pattern(folderPath,[set_type_result set_optimal_criteria]);
Bayesiopt_eDot=cell(NAnalysis,1);
for i=1:NAnalysis
    index=contains(fileNames,['_' analysis_list{i} '_']);
    Bayesiopt_eDot(i)={load(fullfile(folderPath,fileNames_bayes_eDot{index})).resultsBayesopt};
end
%% Knee contact forces: JRXN_knee
set_optimal_criteria= '_VERX_JRXN_knee'; 

set_type_result     = 'MRSOptimal';       
fileNames=get_files_from_pattern(folderPath,[set_type_result set_optimal_criteria]);
MRS_JRXN=cell(NAnalysis,1);
for i=1:NAnalysis
    index=contains(fileNames,['_' analysis_list{i} '_']);
    MRS_JRXN(i)={load(fullfile(folderPath,fileNames{index}))};
end

set_type_result   = 'Bayesopt';       
fileNames_bayes_JRXN=get_files_from_pattern(folderPath,[set_type_result set_optimal_criteria]);
Bayesiopt_JRXN=cell(NAnalysis,1);
for i=1:NAnalysis
    index=contains(fileNames,['_' analysis_list{i} '_']);
    Bayesiopt_JRXN(i)={load(fullfile(folderPath,fileNames_bayes_JRXN{index})).resultsBayesopt};
end
%% Compute relevant metrics: Positive Power & other objectives
PowerPos_nor_sum_eDot=compute_power(MRS_eDot,NAnalysis);
PowerPos_nor_sum_JRXN=compute_power(MRS_JRXN,NAnalysis);

assistance_goal='eDot_MCLU24';
J_baseline_eDot=computeBaseline(assistance_goal);

JeDot_fromJRXN=zeros(NAnalysis,1);
for iTrial=1:NAnalysis 
    MRS_sel=MRS_JRXN{iTrial};
    [J_sel,Jrel_sel,J_extra_sel] = computeOuterLoopFunction(MRS_sel.Misc,MRS_sel.Results,assistance_goal);
    J_value=(J_sel-J_baseline_eDot)/J_baseline_eDot*100;
    JeDot_fromJRXN(iTrial)=J_value;
end
%%
clc;
assistance_goal='JRXN_knee';
J_baseline_eDot=computeBaseline(assistance_goal);
JJRXN_fromeDot=zeros(NAnalysis,1);

for iTrial=1:NAnalysis 
    MRS_sel=MRS_eDot{iTrial};

    % this is because I am using two computers
    MRS_sel.Misc.OutPath='C:\Users\movea\Dropbox\PostDoc Work\Simulation\Code\ProjectResults\Codesign\Pilot\sub1\v2_t1';
    MRS_sel.Misc.SetupPath='C:\Users\movea\Dropbox\PostDoc Work\Simulation\Code\GenericSetups';

    SubjectSelection='sub1';
    MotionSelection='v2_t1';
    extLoadsInfo.fileName =fullfile(currentFolder,'Database',SubjectSelection,'extLoads',['data_'  MotionSelection '.mot']);
    extLoadsInfo.setupName=fullfile(currentFolder,'Database',SubjectSelection,'extLoads',['setup_' MotionSelection '.xml']);
    MRS_sel.Misc.extLoadsInfo=extLoadsInfo;

    MRS_sel.Misc.IKfile= {fullfile(currentFolder,'Database',SubjectSelection,'IK',['IK_' SubjectSelection '_' MotionSelection '.mot'])};

    ModelVariation='';
    MRS_sel.Misc.model_path=fullfile(currentFolder,'Database',SubjectSelection,'model',['model_rajagopal2022_' SubjectSelection ModelVariation '.osim']);
    % !!!

    [J_sel,Jrel_sel,J_extra_sel] = computeOuterLoopFunction(MRS_sel.Misc,MRS_sel.Results,assistance_goal);
    J_value=(J_sel-J_baseline_eDot)/J_baseline_eDot*100;
    JJRXN_fromeDot(iTrial)=J_value;
end

%%
clc;
MinObjective_ordered_eDot=zeros(NAnalysis,1);
MinObjective_ordered_JRXN=zeros(NAnalysis,1);
for i=1:NAnalysis
    MinObjective_ordered_eDot(i,1)=Bayesiopt_eDot{i}.MinObjective;
    MinObjective_ordered_JRXN(i,1)=Bayesiopt_JRXN{i}.MinObjective;
end

muscleGroup = {'AP' 'KE' 'HF' 'HB'};
% muscleGroup = {'AP' 'KE'};
actuatorType = {'Q', 'A'};
arrangement = [1 2];
analysis_list_ss = generate_analysis_list(muscleGroup, actuatorType, arrangement);
NAnalysis_ss     = length(analysis_list_ss);

index_ss=zeros(NAnalysis_ss,1);
for i=1:NAnalysis_ss
    index_ss(i)         = find(strcmp(analysis_list,analysis_list_ss(i)));
end

choose_color_palette=2;
if choose_color_palette==1
    colors_hex = {
        '#1E3A8A', ... % KE(Q) - Deep navy blue
        '#3B82F6', ... % KE(A) - Bright blue
        '#047857', ... % AP(Q) - Deep green
        '#10B981', ... % AP(A) - Emerald green
        '#7C3AED', ... % AP(Q)&KE(Q) - Royal purple
        '#6366F1', ... % AP(A)&KE(Q) - Indigo
        '#DC2626', ... % AP(Q)&KE(A) - Deep red
        '#EA580C'  ... % AP(A)&KE(A) - Orange
        };
    colors_rgb = hex2rgb(colors_hex);

    visual_marker=['s' 's' 's' 's' 'o' 'o' 'o' 'o'];
elseif choose_color_palette==2
    colors_hex = generate_color_map(muscleGroup, actuatorType, arrangement);
    colors_rgb = hex2rgb(colors_hex);

    visual_marker=generate_visual_markers(muscleGroup, actuatorType, arrangement);
end
%% Plot results
JRXN_burden=PowerPos_nor_sum_JRXN;
eDot_burden=PowerPos_nor_sum_eDot;

% JRXN_burden=PowerPos_nor_sum_JRXN/max(PowerPos_nor_sum_JRXN(index_ss))*10;
% eDot_burden=PowerPos_nor_sum_eDot/max(PowerPos_nor_sum_eDot(index_ss))*10;

% max(PowerPos_nor_sum_eDot)
% max(PowerPos_nor_sum_JRXN)

% Performance and hardware burden
figure(1); clf; set(gcf,'color','w');

text_fontSize=12;
label_x_axis='configurations';

subplot(2,2,1)
X = categorical(analysis_list_ss);
X = reordercats(X,analysis_list_ss);
for j=1:NAnalysis_ss
    index =index_ss(j);
    metric=MinObjective_ordered_eDot(index);
    bar(j,metric,'FaceColor',colors_hex{j}); hold on

    text(j,metric-2,[num2str(metric,'%1.0f')],'HorizontalAlignment','center','FontSize',text_fontSize)
end
xticks(1:length(X));    xticklabels(X); ylim([-40 0])
ylabel('Δ metabolic cost [%]'); xlabel(label_x_axis);
set(gca,'FontSize',12);
title('Performance','FontSize',15);

subplot(2,2,3)
X = categorical(analysis_list_ss);
X = reordercats(X,analysis_list_ss);
for j=1:NAnalysis_ss
    index=index_ss(j);
    metric=eDot_burden(index);
    bar(j,metric,'FaceColor',colors_hex{j}); hold on

    text(j,metric+0.5,[num2str(metric,'%1.1f')],'HorizontalAlignment','center','FontSize',text_fontSize)
end
xticks(1:length(X));    xticklabels(X); ylim([0 12])
ylabel('hardware burden  [ ]'); xlabel(label_x_axis);
set(gca,'FontSize',12);
title('Hardware burden','FontSize',15);

subplot(2,2,2)
X = categorical(analysis_list_ss);
X = reordercats(X,analysis_list_ss);
for j=1:NAnalysis_ss
    index=index_ss(j);
    metric=MinObjective_ordered_JRXN(index);
    bar(j,metric,'FaceColor',colors_hex{j}); hold on

    text(j,metric-3,[num2str(metric,'%1.0f')],'HorizontalAlignment','center','FontSize',text_fontSize)
end
xticks(1:length(X));        xticklabels(X); ylim([-70 0])
ylabel('Δ knee contact force [%]');  xlabel(label_x_axis);
set(gca,'FontSize',12);
title('Performance','FontSize',15);

subplot(2,2,4)
X = categorical(analysis_list_ss);
X = reordercats(X,analysis_list_ss);
for j=1:NAnalysis_ss
    index=index_ss(j);
    metric=JRXN_burden(index);
    bar(j,metric,'FaceColor',colors_hex{j}); hold on

    text(j,metric+0.5,[num2str(metric,'%1.1f')],'HorizontalAlignment','center','FontSize',text_fontSize)
end
xticks(1:length(X));    xticklabels(X); ylim([0 12])
ylabel('hardware burden  [ ]'); xlabel(label_x_axis);
set(gca,'FontSize',12);
title('Hardware burden','FontSize',15);

%%
figure(2); clf; set(gcf,'color','w');

% visual_off=[1 -1 -1 1 1 1 1 -1]; 
visual_off=ones(NAnalysis_ss,1); 
% subplot(2,2,1); hold on
subplot(1,1,1); hold on
for i=1:NAnalysis_ss
    index=index_ss(i);
    x=eDot_burden(index);
    % y=MinObjective_ordered_eDot(index);
    y=-MinObjective_ordered_eDot(index);
    plot(x,y,visual_marker(i),'MarkerSize',15,'MarkerFaceColor',colors_hex{i},'MarkerEdgeColor','none')
    % text(x,y+3.0*visual_off(i),analysis_list{i},'FontSize',12,'HorizontalAlignment','center')
    text(x,y+1.0*visual_off(i),analysis_list{i},'FontSize',12,'HorizontalAlignment','center')
end
% yline(0); xline(0,':k'); xlim([-0.5 10]); ylim([-40 15])
yline(0); xline(0,':k'); xlim([-0.5 10]); ylim([-15 40])
ylabel('Δ metabolic cost [%]');
xlabel('hardware burden [ ]'); %+MechPower/maxMechPower
set(gca,'FontSize',12);
% title('Performance (direct) vs Burden','FontSize',15);
title('Performance vs Hardware Burden','FontSize',15);
%%
% visual_off=[1 -1 -1 1 1 1 1 -1]; 
visual_off=ones(NAnalysis_ss,1);
subplot(2,2,2); hold on
for i=1:NAnalysis_ss
    index=index_ss(i);
    x=eDot_burden(index);
    y=JJRXN_fromeDot(index);
    plot(x,y,visual_marker(i),'MarkerSize',15,'MarkerFaceColor',colors_hex{i},'MarkerEdgeColor','none')
    % text(x,y+1.5*visual_off(i),order{i},'FontSize',12,'HorizontalAlignment','center')
end
yline(0); xline(0,':k'); xlim([-0.5 10]); ylim([-60 0])
ylabel('Δ knee contact force [%]');
xlabel('hardware burden [ ]'); %+MechPower/maxMechPower
set(gca,'FontSize',12);
title('Performance (indirect) vs Burden','FontSize',15);

% visual_off=[-1 1 1 1 1 1 1 -1];
visual_off=ones(NAnalysis_ss,1);
subplot(2,2,4); hold on
for i=1:NAnalysis_ss
    index=index_ss(i);
    x=JRXN_burden(index);
    y=MinObjective_ordered_JRXN(index);
    plot(x,y,visual_marker(i),'MarkerSize',15,'MarkerFaceColor',colors_hex{i},'MarkerEdgeColor','none')
    text(x,y+3*visual_off(i),analysis_list{i},'FontSize',12,'HorizontalAlignment','center')
end
yline(0); xline(0,':k'); xlim([-0.5 10]); ylim([-60 0]);
ylabel('Δ knee contact force [%]');
xlabel('hardware burden [ ]');
set(gca,'FontSize',12);
title('Performance (direct) vs Burden','FontSize',15);


% visual_off=[-1 1 1 1 1 1 1 -1];
visual_off=ones(NAnalysis_ss,1);
subplot(2,2,3); hold on
for i=1:NAnalysis_ss
    index=index_ss(i);
    x=JRXN_burden(index);
    y=JeDot_fromJRXN(index);
    plot(x,y,visual_marker(i),'MarkerSize',15,'MarkerFaceColor',colors_hex{i},'MarkerEdgeColor','none')
    % text(x,y+3*visual_off(i),order{i},'FontSize',12,'HorizontalAlignment','center')
end
yline(0); xline(0,':k'); xlim([-0.5 10]); ylim([-40 15])
ylabel('Δ metabolic cost [%]');
xlabel('hardware burden [ ]');
set(gca,'FontSize',12);
title('Performance (indirect) vs Burden','FontSize',15);

function rgb = hex2rgb(hex_cells)
    rgb = zeros(length(hex_cells), 3);
    for i = 1:length(hex_cells)
        hex = hex_cells{i};
        rgb(i,:) = [hex2dec(hex(2:3)), hex2dec(hex(4:5)), hex2dec(hex(6:7))] / 255;
    end
end

function fileNames=get_files_from_pattern(folderPath,label_pattern_list)
label_pattern=label_pattern_list;
pattern = ['*' label_pattern '*.mat'];

files = dir(fullfile(folderPath, pattern));
% Display the results
for i = 1:length(files)
    fprintf('Found: %s\n', files(i).name);
end

% Get just the file names
fileNames = {files.name};
end

function PowerPos_nor_sum=compute_power(MRS_sel,NAnalysis)
PowerPos_trial=zeros(1,NAnalysis);
PowerNet_trial=zeros(1,NAnalysis);
PowerPos_list_all=zeros(NAnalysis,2); % maximum number of devices;
PowerPos_list_name=strings(NAnalysis,2);
for iTrial=1:NAnalysis
    Device=MRS_sel{iTrial}.Misc.Device;
    DOFNames=MRS_sel{iTrial}.DatStore.DOFNames;
    time=MRS_sel{iTrial}.DatStore.Mesh.t'; % sec

    nDevs=length(Device);
    PowerPos_list=zeros(1,nDevs);
    PowerNet_list=zeros(1,nDevs);
    for iDev=1:nDevs
        Torque=Device{iDev}.Assistance.Profile.Torque; % Nm
        devJoint=Device{iDev}.MuscleGroup{1};
        devDir  =Device{iDev}.MuscleGroup{2};
        actType =Device{iDev}.Type{1};

        q_exp=MRS_sel{iTrial}.DatStore.q_exp(:,strcmp(DOFNames,devJoint)); % deg

        dq_exp=gradient(q_exp)./gradient(time);

        dq_exp_rad=dq_exp*pi/180;
        Power=dq_exp_rad.*Torque*devDir;

        PowerPos = trapz(time(Power>0),Power(Power>0));
        PowerPos_list(1,iDev)=PowerPos;

        PowerNet = trapz(time,Power);
        PowerNet_list(1,iDev)=PowerNet;

        if strcmp(actType,'active')
            PowerPos_list_act=PowerPos;
            abb='A';
        elseif strcmp(actType,'quasi-passive')
            PowerPos_list_act=0;
            abb='Q';
        end
        PowerPos_list_all(iTrial,iDev)=PowerPos_list_act;
        PowerPos_list_name(iTrial,iDev)=abb;

        to_plot=0;
        if to_plot==1
            subplot(3,nDevs,(1-1)*nDevs+iDev);
            plot(time,q_exp);
            ylabel('angle [deg]')
            subplot(3,nDevs,(2-1)*nDevs+iDev);
            plot(time,Torque);
            ylabel('torque [Nm]')
            subplot(3,nDevs,(3-1)*nDevs+iDev);
            plot(time,Power);
            ylabel('power [W]')
            title(['positive Power = ' num2str(PowerPos,'%1.1f') 'J'])

            subplot(3,nDevs,(1-1)*nDevs+iDev);
            title(Device{iDev}.Type{1});
        end

    end
    PowerPos_trial(iTrial)=sum(PowerPos_list);
    PowerNet_trial(iTrial)=sum(PowerNet_list);
end


ind_Q=PowerPos_list_name=="Q";

PowerPos_list_all(ind_Q)=1; % 1 J. Assuming. Solenoid 3.6 W and time = 100 ms (10% gait cycle), 0.36 J (=1) https://www.farnell.com/datasheets/3625674.pdf?_gl=1*10n0fj0*_gcl_aw*R0NMLjE3NjI4NTE0NjcuQ2p3S0NBaUEyc3ZJQmhCLUVpd0FSV0RQam9FUHptYUtHQ09acDBsSmJjVGZYRFBYcngzcE55eERwWDJrMTQ1Y1VBRFVqRXpiSHhfZEtSb0NGblVRQXZEX0J3RQ..*_gcl_au*MjA1MTA3NTE4My4xNzYyODUxNDY4
% https://se.farnell.com/adafruit/412/frame-type-d-type-box-type/dp/2816374?gross_price=true&CMP=KNC-GSE-GEN-Shopping-Pmax-Mid-ROAS&mckv=_dc|pcrid||&gad_source=1&gad_campaignid=22205790637&gbraid=0AAAAAD8yeHnAayAiNFDTwidIwx6m8KhNs&gclid=CjwKCAiA2svIBhB-EiwARWDPjoEPzmaKGCOZp0lJbcTfXDPXrx3pNyxDpX2k145cUADUjEzbHx_dKRoCFnUQAvD_BwE
max_val=max(sum(PowerPos_list_all,2),[],'all');


% PowerPos_list_all(ind_Q)=0.1*max_val;

PowerPos_nor=PowerPos_list_all/max_val*10;
% PowerPos_nor(ind_Q)=1;
PowerPos_nor_sum=sum(PowerPos_nor,2);
end

function J_baseline = computeBaseline(assistance_goal)
currentFolder = pwd;

% select subject and motion
info.SubjectSelection ='sub1';        % select subject
info.MotionSelection  ='v2_t1';       % read as-> v(velocity)== 1[slow], 2[normal], 3[fast]) _ t (trial)== 1st, 2nd, or 3rd trial
info.currentFolder    = currentFolder;
% assistance_goal       ='eDot_MCLU24'; % eDot eDot_MCLU24 gasForces KJMusForces RJXN_knee RJXN_knee_par
OutName_optimal       ='_VERX';

% setup unassisted condition
[Misc,Results_normal,DatStore_normal]=setup_and_run_unassisted_condition(info);
%% Check results
Misc.MotionSelection   ='oneSingleMotion';
Misc.ForceLabel        ='SO_label_rajagopal';
Misc.extra_folder_name ='temp';
Misc.extra_file_name   ='unassisted';

[J_baseline,Jrel_TS,J_extra_normal] = computeOuterLoopFunction(Misc,Results_normal,assistance_goal); % this is my starting point, the Edot at unassisted conditions.
end

function analysis_list = generate_analysis_list(muscleGroup, actuatorType, arrangement)
% GENERATE_ANALYSIS_LIST Generate analysis list based on configuration
% Inputs:
%   muscleGroup: cell array of muscle groups (e.g., {'AP', 'KE', 'HF', 'HA'})
%   actuatorType: cell array of actuator types (e.g., {'Q', 'A'})
%   arrangement: array of arrangements [1, 2] where:
%       1 = single muscle groups
%       2 = paired combinations

    analysis_list = {};
    
    % Single muscle group arrangements (arrangement 1)
    if ismember(1, arrangement)
        for i = 1:length(muscleGroup)
            for j = 1:length(actuatorType)
                config_name = [muscleGroup{i} '(' actuatorType{j} ')'];
                analysis_list{end+1} = config_name;
            end
        end
    end
    
    % Paired combinations (arrangement 2)
    if ismember(2, arrangement)
        % Generate all combinations of muscle groups
        muscle_combinations = nchoosek(muscleGroup, 2);
        
        for i = 1:size(muscle_combinations, 1)
            muscle1 = muscle_combinations{i, 1};
            muscle2 = muscle_combinations{i, 2};
            
            % Generate all actuator type combinations for this muscle pair
            for j = 1:length(actuatorType)
                for k = 1:length(actuatorType)
                    config_name = [muscle1 '(' actuatorType{j} ')&' muscle2 '(' actuatorType{k} ')'];
                    analysis_list{end+1} = config_name;
                end
            end
        end
    end
end

function colors_hex = generate_color_map(muscleGroups, actuatorTypes, arrangements)
% Simple color mapping that maintains relationships
    
    % Base colors for muscle groups
    base_colors = {
        [30, 58, 138],    % AP - Deep Blue
        [5, 120, 87],     % KE - Deep Green  
        [120, 53, 15],    % HF - Brown
        [126, 34, 206]    % HA - Purple
    };
    
    % Actuator type modifiers
    actuator_modifiers = [
        0.8, 0.8, 0.8;    % Q - Original
        1.5, 1.5, 1.5     % A - Brighter
    ];
    
    colors_hex = {};
    idx = 1;
    
    % Generate all configurations
    analysis_list = generate_analysis_list(muscleGroups, actuatorTypes, arrangements);
    
    for i = 1:length(analysis_list)
        config = analysis_list{i};
        
        if contains(config, '&')
            % Paired configuration - blend colors
            parts = strsplit(config, '&');
            color1 = get_base_color(parts{1}, muscleGroups, base_colors, actuator_modifiers, actuatorTypes);
            color2 = get_base_color(parts{2}, muscleGroups, base_colors, actuator_modifiers, actuatorTypes);
            
            % Blend the two colors
            final_color = round((color1 + color2) / 2);
        else
            % Single configuration
            final_color = get_base_color(config, muscleGroups, base_colors, actuator_modifiers, actuatorTypes);
        end
        
        colors_hex{idx} = sprintf('#%02X%02X%02X', final_color(1), final_color(2), final_color(3));
        idx = idx + 1;
    end
end

function color = get_base_color(config, muscleGroups, base_colors, modifiers, actuatorTypes)
% Extract base color and apply actuator modifier
    
    % Find muscle group
    for mg_idx = 1:length(muscleGroups)
        if contains(config, muscleGroups{mg_idx})
            base_color = base_colors{mg_idx};
            break;
        end
    end
    
    % Find actuator type and apply modifier
    for at_idx = 1:length(actuatorTypes)
        if contains(config, ['(' actuatorTypes{at_idx} ')'])
            modifier = modifiers(at_idx, :);
            break;
        end
    end
    
    color = min(round(base_color .* modifier), 255);
end

function visual_marker = generate_visual_markers(muscleGroup, actuatorType, arrangement)
% GENERATE_VISUAL_MARKERS Generate marker styles based on arrangement type
% Inputs:
%   muscleGroup: cell array of muscle groups (e.g., {'AP', 'KE', 'HF', 'HA'})
%   actuatorType: cell array of actuator types (e.g., {'Q', 'A'})
%   arrangement: array of arrangements [1, 2] where:
%       1 = single muscle groups (marker: 's' - square)
%       2 = paired combinations (marker: 'o' - circle)

    visual_marker = {};
    
    % Single muscle group arrangements (arrangement 1)
    if ismember(1, arrangement)
        num_singles = length(muscleGroup) * length(actuatorType);
        singles_markers = repmat({'s'}, 1, num_singles);
        visual_marker = [visual_marker, singles_markers];
    end
    
    % Paired combinations (arrangement 2)
    if ismember(2, arrangement)
        num_pairs = nchoosek(length(muscleGroup), 2) * length(actuatorType)^2;
        pairs_markers = repmat({'o'}, 1, num_pairs);
        visual_marker = [visual_marker, pairs_markers];
    end
    
    % Convert to character array if needed (for backward compatibility)
    if all(cellfun(@(x) ischar(x) && length(x) == 1, visual_marker))
        visual_marker = [visual_marker{:}];
    end
end