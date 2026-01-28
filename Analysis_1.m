%% SET PATHS
mainPath      ='C:\Users\movea\Documents\GitHub\GoalOptimized\ProjectResults\DSE';
trialPath     ='v2_t1';

currentFolder=pwd;
addpath(genpath(currentFolder));

sub_list=1:5;
nSubs=length(sub_list);
%% LOAD MRS Je0
folderSpec   = 'Je0';
specificFile = 'Je0';

MRS_Je0=cell(nSubs,1);
for iSub=1:nSubs
    MRS_Je0(iSub,1)={load(fullfile(mainPath,['sub' num2str(iSub)],trialPath,folderSpec,[specificFile 'Results.mat']))};
end

%% LOAD MRS Je
folderSpec   = 'Je';
specificFile = 'Je';

MRS_Je=cell(nSubs,1);
for iSub=1:nSubs
    MRS_Je(iSub,1)={load(fullfile(mainPath,['sub' num2str(iSub)],trialPath,folderSpec,[specificFile 'Results.mat']))};
end
%% LOAD MRS JeS
folderSpec   = 'JeS0';
specificFile = 'JeS0';

Syn_list=[4 5 6];
nSyns=length(Syn_list);
MRS_JeS=cell(nSubs,nSyns);
metric_JeS =cell(nSubs,1);
for iSub=1:nSubs
    for iSyn=1:nSyns
        sSyn=Syn_list(iSyn);
        MRS_JeS(iSub,iSyn)={load(fullfile(mainPath,['sub' num2str(iSub)],trialPath,folderSpec,[specificFile num2str(sSyn) 'Results.mat']))};
    end
    metric_JeS(iSub)={load(fullfile(mainPath,['sub' num2str(iSub)],trialPath,folderSpec,'synergy_metrics.mat'))};
end
%% Summary
[NSYN, VAF, RMSE, RMSE_m, R_m] = deal(zeros(nSubs, length(Syn_list)));

for iSub=1:nSubs
    metric_sel=metric_JeS{iSub};
    NSYN(iSub,:)  =metric_sel.synMetrics.N';
    VAF(iSub,:)   =metric_sel.synMetrics.VAF;
    RMSE(iSub,:)  =metric_sel.synMetrics.RMSE;
    RMSE_m(iSub,:)=mean(metric_sel.synMetrics.RMSE_per_muscle,2);
    R_m(iSub,:)   =mean(metric_sel.synMetrics.R_per_muscle,2);
end

NSYN_mean=mean(NSYN,1)';
VAF_mean = mean(VAF,1)';   VAF_std = std(VAF,0,1)';
RMSE_mean = mean(RMSE,1)'; RMSE_std = std(RMSE,0,1)';
RMSE_m_mean = mean(RMSE_m,1)'; RMSE_m_std = std(RMSE_m,0,1)';
R_m_mean = mean(R_m,1)';   R_m_std = std(R_m,0,1)';

% Create formatted strings for table
VAF_str = arrayfun(@(m,s) sprintf('%.0f ± %.1f', m*100, s*100), VAF_mean, VAF_std, 'UniformOutput', false);
RMSE_str = arrayfun(@(m,s) sprintf('%.3f ± %.3f', m, s), RMSE_mean, RMSE_std, 'UniformOutput', false);
RMSE_m_str = arrayfun(@(m,s) sprintf('%.3f ± %.3f', m, s), RMSE_m_mean, RMSE_m_std, 'UniformOutput', false);
R_m_str = arrayfun(@(m,s) sprintf('%.2f ± %.2f', m, s), R_m_mean, R_m_std, 'UniformOutput', false);

% Display table
T = table(NSYN_mean, VAF_str, RMSE_str, RMSE_m_str, R_m_str, ...
          'VariableNames', {'Synergies', 'VAF (%)', 'RMSE', 'RMSE_m', 'R_m'});
disp(T)
%% Reconstruction error
delta_VAF =zeros(nSubs,nSyns);
delta_RMSE=zeros(nSubs,nSyns);
for iSub=1:nSubs
    for iSyn=1:nSyns
        delta_VAF(iSub,iSyn) =MRS_JeS{iSub,iSyn}.Results.SynergyControl.VAF;
        delta_RMSE(iSub,iSyn)=MRS_JeS{iSub,iSyn}.Results.SynergyControl.RMSE;
    end
end

% Compute mean and std
d_VAF_mean = mean(delta_VAF,1)';   d_VAF_std = std(delta_VAF,0,1)';
d_RMSE_mean = mean(delta_RMSE,1)'; d_RMSE_std = std(delta_RMSE,0,1)';

% Create formatted strings
d_VAF_str = arrayfun(@(m,s) sprintf('%.0f ± %.0f', m*100, s*100), d_VAF_mean, d_VAF_std, 'UniformOutput', false);
d_RMSE_str = arrayfun(@(m,s) sprintf('%.3f ± %.3f', m, s), d_RMSE_mean, d_RMSE_std, 'UniformOutput', false);

% Display table
T = table(Syn_list', d_VAF_str, d_RMSE_str, ...
          'VariableNames', {'Synergies', 'ΔVAF (%)', 'ΔRMSE'});
disp(T)
%%
EMGfile='C:\Users\movea\Documents\GitHub\GoalOptimized\emgFiles.mot';
EMG_data = ReadMotFile(EMGfile);
%%
EMG_names=EMG_data.names;
EMG_val  =EMG_data.data;
data_lengthEMG=length(EMG_val(:,1));
gait_cycle_EMG=linspace(0,100,data_lengthEMG);
%%
clf;
for i=2:39
    subplot(5,8,i)
    plot(EMG_val(:,i))
    title([EMG_names{i} 'N' num2str(i)])
end
%% Muscle activations
Misc=MRS_Je{1,1}.Misc;
[gait_cycle,~]=computeGC(Misc.time,Misc.extra_frames);
gait_cycle_interp=linspace(0,100,101);
fe=Misc.extra_frames;

MuscleNames=MRS_Je{1,1}.Results.MuscleNames;
MuscleNamesSyn_list     ={{'soleus_r'} {'tibant_r'} {'gasmed_r'} {'vaslat_r'}  {'addlong_r'} {'recfem_r'} {'sart_r'} {'glmed2_r'} {'semimem_r'} {'semiten_r'}};
MuscleNamesSyn_list_full={{'soleus'} {'tib. ant.'} {'gas. med.'} {'vas. lat.'}  {'add. long.'} {'rec. fem.'} {'sart'} {'glut. med.'} {'semimem.'} {'semiten.'}};
nMuscleSyns=length(MuscleNamesSyn_list);
musSim_ind=zeros(1,nMuscleSyns);
musExp_ind=zeros(1,nMuscleSyns);
for iMus=1:nMuscleSyns
    musSim_ind(1,iMus)=find(strcmp(MuscleNames,MuscleNamesSyn_list{iMus}));
    musExp_ind(1,iMus)=find(strcmp(EMG_names  ,MuscleNamesSyn_list{iMus}{1}(1:end-2)));
end

clc;
% Preallocate for all subjects
nSubs = size(MRS_Je,1);
nSyns = size(MRS_JeS,2);
all_MA_Je0= zeros(nSubs, nMuscleSyns, 101);
all_MA_Je = zeros(nSubs, nMuscleSyns, 101);
all_MA_JeS = zeros(nSyns, nSubs, nMuscleSyns, 101);
all_EMG = zeros(nMuscleSyns, data_lengthEMG);

% Collect data from all subjects
for iSub = 1:nSubs
    for iMus = 1:nMuscleSyns

        % Je0 activations
        MA_orig = MRS_Je0{iSub,1}.Results.MActivation.genericMRS(musSim_ind(iMus),1+fe:end-fe-1);
        gait_cycle_orig=linspace(0,100,length(MA_orig));
        all_MA_Je0(iSub,iMus,:) = interp1(gait_cycle_orig, MA_orig, gait_cycle_interp, 'pchip');

        % Je activations
        MA_orig = MRS_Je{iSub,1}.Results.MActivation.genericMRS(musSim_ind(iMus),1+fe:end-fe-1);
        gait_cycle_orig=linspace(0,100,length(MA_orig));
        all_MA_Je(iSub,iMus,:) = interp1(gait_cycle_orig, MA_orig, gait_cycle_interp, 'pchip');
         
        % JeS activations
        for iSyn = 1:nSyns
            MA_syn = MRS_JeS{iSub,iSyn}.Results.MActivation.genericMRS(musSim_ind(iMus),1+fe:end-fe-1);
            gait_cycle_orig=linspace(0,100,length(MA_syn));
            all_MA_JeS(iSyn,iSub,iMus,:) = interp1(gait_cycle_orig, MA_syn, gait_cycle_interp, 'pchip');
        end
        % EMG data
        all_EMG(iMus,:) = EMG_val(:,musExp_ind(iMus));
    end
end
%%
label_list ={'EMG' 'Effort'  'Syn #4' 'Syn #5' 'Syn #6'};
color_list={'#7393B3' '#4c4c4c' '#1f7a8c' '#d80032' '#245501' '#FFC0CB'};
line_list={'-' ':' '--' '-.'};

% Create figure
fig = figure(1); set(gcf,'color','w','Visible','on'); clf; 
for iMus = 1:nMuscleSyns
    subplot(2,5,iMus); hold on;
    EMG_TS=EMG_val(:,musExp_ind(1,iMus));

    % 1. EMG shaded region (mean ± std)
    curve1 = zeros(data_lengthEMG,1);
    curve2 = EMG_TS;
    x2 = [gait_cycle_EMG, fliplr(gait_cycle_EMG)];
    inBetween = [curve1', fliplr(curve2')];
    color_emg=hex2rgb(color_list{1}); 

    pl(1)=fill(x2, inBetween,color_emg,'LineStyle','none','facealpha',.6,'displayName','EMG');
    
    % 2. Effort-minimized (mean ± std band)
    MA_Je_mean = squeeze(mean(all_MA_Je(:,iMus,:), 1));
    MA_Je_std = squeeze(std(all_MA_Je(:,iMus,:), 0, 1));
    fill([gait_cycle_interp, fliplr(gait_cycle_interp)], ...
         [(MA_Je_mean + MA_Je_std)', fliplr((MA_Je_mean - MA_Je_std)')], ...
         hex2rgb(color_list{2}), 'FaceAlpha', 0.1, 'EdgeColor', 'none');
    pl(2)=plot(gait_cycle_interp, MA_Je_mean, 'Color', color_list{2}, 'LineWidth', 3, 'DisplayName', 'Effort');
    
    % 3. Synergy-informed (mean lines only)
    for iSyn = 1:nSyns
        MA_JeS_mean = squeeze(mean(all_MA_JeS(iSyn,:,iMus,:), 2));
        MA_JeS_std = squeeze(std(all_MA_JeS(iSyn,:,iMus,:), 0, 2));
        fill([gait_cycle_interp, fliplr(gait_cycle_interp)], ...
             [(MA_JeS_mean + MA_JeS_std)', fliplr((MA_JeS_mean - MA_JeS_std)')], ...
        hex2rgb(color_list{iSyn+2}), 'FaceAlpha', 0.1, 'EdgeColor', 'none');

        pl(iSyn+2)=plot(gait_cycle_interp, MA_JeS_mean, 'Color', color_list{iSyn+2}, 'LineWidth', 3, ...
             'LineStyle', line_list{iSyn}, 'DisplayName', label_list{iSyn+2});
    end
    
    % Je0
    MA_Je0_mean = squeeze(mean(all_MA_Je0(:,iMus,:), 1));
    pl(6)=plot(gait_cycle_interp, MA_Je0_mean, 'Color', color_list{6}, 'LineWidth', 3, 'DisplayName', 'Je0');

    % Formatting
    legend(pl,'Location', 'best', 'NumColumns', 1); legend boxoff;
    xticks(0:50:100); yticks(0:0.3:1.4);
    set(gca, 'FontSize', 13);
    title(MuscleNamesSyn_list_full{iMus}{1});
    axis([0 100 0 1.15]);
end

% Add global labels
han = axes(fig, 'visible', 'off'); 
han.Title.Visible = 'on'; han.XLabel.Visible = 'on'; han.YLabel.Visible = 'on';
label_y = ylabel(han, 'muscle activations [ ]', 'FontSize', 20); 
label_y.Position(1) = -0.05; label_y.Position(2) = 0.5;
label_x = xlabel(han, 'gait cycle [%]', 'FontSize', 20);   
label_x.Position(1) = 0.5; label_x.Position(2) = -0.05;
%% Metabolic cost
assistiveGoal         ='eDot_MCLU24';
J_avg =zeros(nSubs,1);
J_avg0=zeros(nSubs,1);
JSyn_avg=zeros(nSyns,1);
ExpGrossMetCost=zeros(nSubs,1);

for iSub = 1:nSubs
    selectSpeed=5;

    % experimental
    ExpGrossMetCost(iSub)=round(EE_data(iSub).exp_GEE(selectSpeed),2);

    % effort min
    Misc   =MRS_Je{iSub,1}.Misc;
    Results=MRS_Je{iSub,1}.Results;
    [J_avg(iSub,1),  ~,  ~]   = computeOuterLoopFunction(Misc,Results,assistiveGoal);

    % synergy
    for iSyn=1:nSyns
        Results=MRS_JeS{iSub,iSyn}.Results;
        Misc=MRS_JeS{iSub,iSyn}.Misc;
        [JSyn_avg(iSub,iSyn),  ~,  ~]   = computeOuterLoopFunction(Misc,Results,assistiveGoal);
    end

    % effort min
    Misc   =MRS_Je0{iSub,1}.Misc;
    Results=MRS_Je0{iSub,1}.Results;
    [J_avg0(iSub,1),  ~,  ~]   = computeOuterLoopFunction(Misc,Results,assistiveGoal);

end

GrossMetCost = [J_avg, JSyn_avg, J_avg0] + 1.2;
completeEdot = [ExpGrossMetCost, GrossMetCost];

completeEdot_mean = mean(completeEdot, 1);
completeEdot_std = std(completeEdot, 0, 1);

nEEcon=size(completeEdot,2);
fig = figure(2); set(gcf,'color','w','Visible','on'); clf; 
font_size = 20;

hb = bar(completeEdot_mean, 'FaceColor', 'flat');
for i = 1:nEEcon
    hb.CData(i,:) = hex2rgb(color_list{i});
end

% Add error bars
hold on;
errorbar(1:nEEcon, completeEdot_mean, completeEdot_std, 'k.', 'LineWidth', 2, 'CapSize', 15);

% Formatting
set(gca, 'XTickLabel', {'EXP', 'Effort', 'Syn #4', 'Syn #5', 'Syn #6','Effort0'}, ...
    'FontSize', font_size, 'XTick', 1:5);
ylabel('whole-body average metabolic rate [W/kg]', 'FontSize', font_size);
ylim([0, max(completeEdot_mean + completeEdot_std) * 1.3]);


% Add values on top of bars
for i = 1:nEEcon
    % Absolute value
    text(i, completeEdot_mean(i) + completeEdot_std(i) + 0.5, ...
        sprintf('%.1f', completeEdot_mean(i)), ...
        'HorizontalAlignment', 'center', 'FontSize', font_size-2, 'FontWeight', 'bold');


    % Add mean line
    yline(completeEdot_mean(1), ':k', 'LineWidth', 2, 'Alpha', 0.5);

    % Percentage change RELATIVE TO EFFORT-MINIMIZED (bar #2)
    selData=1;
    if i ~= selData  % Skip the effort-minimized bar itself
        perc_change = ((completeEdot_mean(i) - completeEdot_mean(selData)) / completeEdot_mean(selData)) * 100;
        text(i, completeEdot_mean(i) + completeEdot_std(i) + 1.25, ...
            sprintf('%+.1f%%', perc_change), ...
            'HorizontalAlignment', 'center', 'FontSize', font_size-2, 'Color', 'r', 'FontWeight', 'bold');
    end
end

%%
% Add mean experimental line
yline(mean(ExpGrossMetCost), ':k', 'LineWidth', 2, 'Alpha', 0.5);

% Formatting
set(gca, 'XTickLabel', {'EXP', 'Effort', 'Syn #4', 'Syn #5', 'Syn #6'}, ...
    'FontSize', font_size, 'XTick', 1:5);
ylabel('whole-body average metabolic rate [W/kg]', 'FontSize', font_size);
ylim([0, max(completeEdot_mean + completeEdot_std) * 1.3]);

% Add values on top of bars
for i = 1:5
    % Absolute value
    text(i, completeEdot_mean(i) + completeEdot_std(i) + 0.5, ...
        sprintf('%.1f', completeEdot_mean(i)), ...
        'HorizontalAlignment', 'center', 'FontSize', font_size-2, 'FontWeight', 'bold');
    
    % Percentage change (relative to EXP) for i>1
    if i > 1
        perc_change = ((completeEdot_mean(i) - completeEdot_mean(1)) / completeEdot_mean(1)) * 100;
        text(i, completeEdot_mean(i) + completeEdot_std(i) + 1.5, ...
            sprintf('%+.1f%%', perc_change), ...
            'HorizontalAlignment', 'center', 'FontSize', font_size-2, 'Color', 'r', 'FontWeight', 'bold');
    end
end
%%
label_list ={'EXP' 'Effort'  'Syn #4' 'Syn #5' 'Syn #6'};
fig=figure(2); set(gcf,'color','w','Visible','on'); clf; 

font_size=20;
hb = bar(completeEdot,'FaceColor','flat');
for i=1:5, hb.CData(i,:)=hex2rgb(color_list{i}); end
set(gca,'XTickLabel',label_list); 
set(gca,'FontSize',font_size);
ylabel('gross metabolic cost [W/kg]'); 
yline(ExpGrossMetCost,':k','LineWidth',2)

ylim([0 7])
text(1, completeEdot(1) + 0.25, sprintf('%.1f', completeEdot(1)), ...
    'HorizontalAlignment', 'center','FontSize',font_size);

% Add percentage change for other bars
for i = 2:5
    perc_change = ((completeEdot(i) - completeEdot(1)) / completeEdot(1)) * 100;
    text(i, completeEdot(i) + 0.25, sprintf('%+.1f%%', perc_change), ...
        'HorizontalAlignment', 'center','FontSize',font_size);
end
%% COMPUTE HILO
DatStore_normal=MRS_eDot{1,1}.DatStore;
Misc=MRS_eDot{1,1}.Misc;
subject_mass=Misc.subject_data.subject_mass;

genVal=100;

Device{1}.Mode        = 'prescribed'; % opts: optimized and prescribed
Device{1}.MuscleGroup = {['ankle_angle_' Misc.gait_data.side_sel] -1};
Device{1}.Type        = {'active' 'spline#N3'}; % opts: active, quasi-passive, passive, EMG-driven
Device{1}.Params      = [54.13 27.9 9.59 genVal]; %54.8 14.5 8.1 51.8
[assistanceInfo]=generateTorque(Device{1},DatStore_normal,Misc.time,Misc.extra_frames);
% Ankle plantarflexion
GaitCycle=assistanceInfo.Profile.GaitCycle;
TorqueHILO=assistanceInfo.Profile.Torque/genVal*0.68;

TorqueHILO_list=zeros(4,length(TorqueHILO));
TorqueHILO_list(1,:)=TorqueHILO;

% Knee extension
Device{1}.MuscleGroup = {['knee_angle_' Misc.gait_data.side_sel] -1};
Device{1}.Type        = {'quasi-passive' 'clutchSpring'}; % opts: active, quasi-passive, passive, EMG-driven
Device{1}.Params      = [0.03  genVal]; %54.8 14.5 8.1 51.8
[assistanceInfo]=generateTorque(Device{1},DatStore_normal,Misc.time,Misc.extra_frames);
TorquePeak=max(assistanceInfo.Profile.Torque);
TorqueHILO=assistanceInfo.Profile.Torque/TorquePeak*0.15;
TorqueHILO_list(2,:)=TorqueHILO;

% Hip flexion
Device{1}.Type        = {'active' 'spline#N3'}; % opts: active, quasi-passive, passive, EMG-driven
Device{1}.Params      = [65.9 33.3 24.7 genVal]; %54.8 14.5 8.1 51.8
[assistanceInfo]=generateTorque(Device{1},DatStore_normal,Misc.time,Misc.extra_frames);
TorqueHILO=assistanceInfo.Profile.Torque/genVal*0.29;
TorqueHILO_list(3,:)=TorqueHILO;

% Hip abduction
Device{1}.Type        = {'active' 'spline#N3'}; % opts: active, quasi-passive, passive, EMG-driven
Device{1}.Params      = [48 15 15 genVal]; %54.8 14.5 8.1 51.8
[assistanceInfo]=generateTorque(Device{1},DatStore_normal,Misc.time,Misc.extra_frames);
TorqueHILO=assistanceInfo.Profile.Torque/genVal*10.83/subject_mass; %10.83/subject_mass
TorqueHILO_list(4,:)=TorqueHILO;
%% PLOT 0
figure(10); clf; set(gcf,'color','w'); clc;
SSyn_master=[1 2 3 4];

Misc=MRS_eDot{1,1}.Misc;
[gait_cycle,~]=computeGC(Misc.time,Misc.extra_frames);

MuscleNames=MRS_eDot{1,1}.Results.MuscleNames;
MuscleNamesSyn_list     ={{'soleus_r'} {'tibant_r'} {'gasmed_r'} {'vaslat_r'}  {'addlong_r'} {'recfem_r'} {'psoas_r'} {'glmed1_r'} {'semimem_r'} {'semiten_r'}};
MuscleNamesSyn_list_full={{'soleus'} {'tib. ant.'} {'gas. med.'} {'vas. lat.'}  {'add. long.'} {'rec. fem.'} {'psoas'} {'glut. med.'} {'semimem.'} {'semiten.'}};
nMuscleSyns=length(MuscleNamesSyn_list);
musSim_ind=zeros(1,nMuscleSyns);
for iMus=1:nMuscleSyns
    musSim_ind(1,iMus)=find(strcmp(MuscleNames,MuscleNamesSyn_list{iMus}));
end


color_list={'k' '#1f7a8c' '#d80032' '#245501'};
label_list={'Effort' 'Syn #4' 'Syn #5' 'Syn #6' };
line_list={'-' ':' '--' '-.'};
for iMus=1:nMuscleSyns
    subplot(2,5,iMus)

    for iSyn=1:length(SSyn_master)
        sSyn=SSyn_master(iSyn);

        Results_Baseline=MRS_base{1,sSyn};
        MActivation_N=Results_Baseline.Results.MActivation.genericMRS(musSim_ind(1,iMus),:);
        plot(gait_cycle,MActivation_N,'Color',color_list{iSyn},'LineWidth',4,'LineStyle',line_list{iSyn}); hold on
    end
    legend(label_list)
    legend boxoff
    set(gca,'FontSize',13);
    title([MuscleNamesSyn_list_full{iMus}{1}])
    xlabel('gait cycle [%]'); ylabel('activation [ ]');
    axis([0 100 0 1.4])
end

figure(11); clf; set(gcf,'color','w'); clc;
NMuscleTot=length(MuscleNames);
for iMus=1:NMuscleTot
    subplot(5,8,iMus)

    for iSyn=1:length(SSyn_master)
        sSyn=SSyn_master(iSyn);

        Results_Baseline=MRS_base{1,sSyn};
        MActivation_N=Results_Baseline.Results.MActivation.genericMRS(iMus,:);
        plot(gait_cycle,MActivation_N,'Color',color_list{iSyn},'LineWidth',4,'LineStyle',line_list{iSyn}); hold on
    end
    % legend(syn_name)
    % legend boxoff
    set(gca,'FontSize',13);
    title([MuscleNames{iMus}(1:end-2)])
    % xlabel('gait cycle [%]'); ylabel('activation [ ]');
    axis([0 100 0 1])
end

for iSyn=1:length(SSyn_master)
    sSyn=SSyn_master(iSyn);
    Results_Baseline=MRS_base{1,sSyn};

    if iSyn>1
        VAF=Results_Baseline.Results.SynergyControl.VAF;
        RMSE=Results_Baseline.Results.SynergyControl.RMSE;
        disp([label_list{iSyn} ': VAF =' num2str(VAF) ' and RMSE =' num2str(RMSE)])
    end
end
%% PLOT 1

DOFNames_list={{'ankle_angle_r'}    {'knee_angle_r'}    {'hip_flexion_r'}    {'hip_adduction_r'}};
DOFNames=MRS_eDot{1,1}.DatStore.DOFNames;
NDOFs = length(DOFNames_list);

MuscleNames_list={{{'soleus_r'} {'gasmed_r'}}     {{'vaslat_r'}  {'vasmed_r'}}   {{'recfem_r'} {'psoas_r'}}  {{'glmed1_r'} {'addlong_r'}} };
MuscleNames_list_full={{{'soleus'} {'gas. med.'}} {{'vas. lat.'} {'vas. med.'}}  {{'rec. fem.'} {'psoas'}}   {{'glut. med.'} {'add. longus'}}};
MuscleNames=MRS_eDot{1,1}.Results.MuscleNames;

Misc=MRS_eDot{1,1}.Misc;
Results=MRS_eDot{1,1}.Results;

extra_frames=Misc.extra_frames;
fSel=1+extra_frames:size(Results.MActivation.genericMRS,2)-1-extra_frames;

[gait_cycle,~]=computeGC(Misc.time,Misc.extra_frames);
gait_cycle_sel=gait_cycle(fSel);

synLabels= {'noSyn' 'Syn#4' 'Syn#5' 'Syn#6'};
% DOFLabels= {'Ankle plantarflexion' 'Knee extension' 'Hip flexion' 'Hip abduction'};
DOFLabels= {'ANKLE PLANTARFLEXION' 'KNEE EXTENSION' 'HIP FLEXION' 'HIP ABDUCTION'};
color_syn={'#4169E2' '#FF8888' '#C60000' '#DF1717'};
% color_syn={'k' '#FF8888' '#C60000' '#DF1717'};

SSyn_master=[2 3 4]; % 3 4;
for iFig=1:length(SSyn_master)
    figure(iFig); clf; set(gcf,'color','w');
    Syn_sel=SSyn_master(iFig);
    SSyn_list=[1 Syn_sel];
    % plTorque_pos={[1 2] [4 5]  [7 8]   [10 11]};
    plTorque_pos={[1 2 7 8] [4 5 10 11]  [13 14 19 20]   [16 17 22 23]};
    plMuscle_pos={[3 9] [6 12] [15 21] [18 24]};
    sign_list=[-1 -1 1 -1];
    offs_list=[0 -0.5 -0.5 0];
    nSyns=length(SSyn_list);
    for iDOF=1:NDOFs
        indexDOF=strcmp(DOFNames,DOFNames_list{iDOF});
        jointID =MRS_eDot{1,1}.DatStore.IDinterp(fSel,indexDOF);

        musSim_ind=zeros(1,2);
        for iMus=1:2
            musSim_ind(1,iMus)=find(strcmp(MuscleNames,MuscleNames_list{iDOF}{iMus}));
        end

        name_list={'J_E' 'J_S'};
        % torques
        % subplot(2,6,plTorque_pos{iDOF})
        subplot(4,6,plTorque_pos{iDOF})
        plg(1)=plot(gait_cycle_sel,sign_list(iDOF).*jointID/subject_mass,'k','LineWidth',4,'DisplayName','ID'); hold on
        for iSyn=1:nSyns
            sSyn=SSyn_list(iSyn);
            Results_Bilevel=MRS_eDot{indexDOF,sSyn};
            Torque=Results_Bilevel.Results.Device{1}.Assistance.Profile.Torque;

            plg(2+iSyn)=plot(gait_cycle,Torque/subject_mass,'Color',color_syn{sSyn},'LineWidth',5,'DisplayName',name_list{iSyn});
        end
        plg(2)=plot(GaitCycle,TorqueHILO_list(iDOF,:),'color','#50C878','LineWidth',7,'LineStyle','-.','DisplayName','EXP');

        disp(['summary with Syn' synLabels{SSyn_list(end)}(end) ' at ' DOFLabels{iDOF} ...
              ': VAF' num2str(Results_Bilevel.Results.SynergyControl.VAF) ...
              ' RMSE' num2str(Results_Bilevel.Results.SynergyControl.RMSE)])
        % Results_Bilevel.Results.SynergyControl

        xlim([0 100]); ylim([-0.35 2.00]+offs_list(iDOF))
        xlabel('gait cycle [%]'); ylabel('torque [Nm/kg]')
        set(gca,'FontSize',13);
        title(DOFLabels{iDOF})
        legend(plg,'NumColumns',2)
        legend boxoff

        % muscles
        pm=plMuscle_pos{iDOF};


        for iMus=1:2 %2
            subplot(4,6,pm(iMus))

            sSynBase=1;
            Results_Baseline=MRS_base{1,sSynBase};
            MActivation_N=Results_Baseline.Results.MActivation.genericMRS(musSim_ind(1,iMus),:); 
            plot(gait_cycle,MActivation_N,'Color','k','LineWidth',3); hold on %'#71797E'

            for iSyn=1:nSyns
                sSyn=SSyn_list(iSyn);
                Results_Bilevel =MRS_eDot{indexDOF,sSyn};

                MActivation_OPT=Results_Bilevel.Results.MActivation.genericMRS(musSim_ind(1,iMus),:);
                plot(gait_cycle,MActivation_OPT,'Color',color_syn{sSyn},'LineWidth',3)
                set(gca,'FontSize',13);
                if iMus==2; xlabel('gait cycle [%]'); end 
                ylabel('activation [ ]')
                
                legend boxoff
                axis([0 100 0 1.4])
                % if iMus==iSyn
                title([MuscleNames_list_full{iDOF}{iMus}{1}]);
                % end
            end
            legend({'No Exo','J_E', 'J_S'},'NumColumns',2,'FontSize',10)
        end
    end

    % Then after creating all plots, adjust positions:
    h = findobj(gcf, 'Type', 'axes');
    for i = 1:length(h)
        pos = get(h(i), 'Position');
        % Reduce width to prevent overlap
        pos(1) = pos(1) - 0.01;  % Reduce width by 10%
        pos(3) = pos(3) * 0.90;  % Reduce width by 10%
        pos(4) = pos(4) * 0.90;  % Reduce width by 10%
        set(h(i), 'Position', pos);
    end
    annotation('textbox',[0 0.97 1 0.02],'String',['Analysis with ' synLabels{SSyn_list(end)}(end) ' synergies'],'FontSize',20,'HorizontalAlignment','center','EdgeColor','none');
end
%%
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

function rgb = hex2rgb(hex)
% Convert hex color '#RRGGBB' to RGB [0-1]
if hex(1)=='#'; hex=hex(2:end); end
rgb = [hex2dec(hex(1:2)) hex2dec(hex(3:4)) hex2dec(hex(5:6))]/255;
end