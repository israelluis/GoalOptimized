%% SET PATHS
folderSpecific='iters200x_A1'; % iters200_A1'iters200x_A1'
versionSpecific='VS'; %VERXS VS

folderPath        = ['C:\Users\movea\Documents\GitHub\GoalOptimized\ProjectResults\Codesign\Pilot\sub1\v2_t1\' folderSpecific];  % Replace with your folder path
currentFolder=pwd;
addpath(genpath(currentFolder));
%% LOAD MRS 
syn_list=[0 4 5 6];
NSyn=length(syn_list);

DOFOrder={'AP(A)' 'KE(Q)' 'HF(A)' 'HB(A)'};

clc;
NFiles=4;
MRS_eDot=cell(NFiles,NSyn);
for iSyn=1:NSyn
    set_optimal_criteria= ['_' versionSpecific num2str(syn_list(iSyn)) '_eDot_MCLU24'];

    set_type_result     = 'MRSOptimal';
    fileNames=get_files_from_pattern(folderPath,[set_type_result set_optimal_criteria])';

    for iFile=1:NFiles
        index=contains(fileNames,['_' DOFOrder{iFile}]);
        MRS_eDot(iFile,iSyn)={load(fullfile(folderPath,fileNames{index}))};
    end
end
%%
NFiles=1;
MRS_base=cell(NFiles,iSyn);
for iSyn=1:NSyn
    fileNames=get_files_from_pattern(folderPath,'NormalSyn')';

    for iFile=1:NFiles
        index=contains(fileNames,['NormalSyn' num2str(syn_list(iSyn))]);
        MRS_base(iFile,iSyn)={load(fullfile(folderPath,fileNames{index}))};
    end
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
NMuscleSyn=length(MuscleNamesSyn_list);
indexMus=zeros(1,NMuscleSyn);
for iMus=1:NMuscleSyn
    indexMus(1,iMus)=find(strcmp(MuscleNames,MuscleNamesSyn_list{iMus}));
end


color_list={'k' '#1f7a8c' '#d80032' '#245501'};
syn_name={'Effort' 'Syn #4' 'Syn #5' 'Syn #6' };
line_list={'-' ':' '--' '-.'};
for iMus=1:NMuscleSyn
    subplot(2,5,iMus)

    for iSyn=1:length(SSyn_master)
        sSyn=SSyn_master(iSyn);

        Results_Baseline=MRS_base{1,sSyn};
        MActivation_N=Results_Baseline.Results.MActivation.genericMRS(indexMus(1,iMus),:);
        plot(gait_cycle,MActivation_N,'Color',color_list{iSyn},'LineWidth',4,'LineStyle',line_list{iSyn}); hold on
    end
    legend(syn_name)
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
        disp([syn_name{iSyn} ': VAF =' num2str(VAF) ' and RMSE =' num2str(RMSE)])
    end
end
%% PLOT 1

DOFNames_list={{'ankle_angle_r'}    {'knee_angle_r'}    {'hip_flexion_r'}    {'hip_adduction_r'}};
DOFNames=MRS_eDot{1,1}.DatStore.DOFNames;
NDOFs = length(DOFNames_list);

MuscleNames_list={{{'soleus_r'} {'tibant_r'}}     {{'vaslat_r'}  {'vasmed_r'}}   {{'recfem_r'} {'psoas_r'}}  {{'glmed1_r'} {'addlong_r'}} };
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

SSyn_master=[3]; % 3 4;
for iFig=1:length(SSyn_master)
    figure(iFig); clf; set(gcf,'color','w');
    Syn_sel=SSyn_master(iFig);
    SSyn_list=[1 Syn_sel];
    % plTorque_pos={[1 2] [4 5]  [7 8]   [10 11]};
    plTorque_pos={[1 2 7 8] [4 5 10 11]  [13 14 19 20]   [16 17 22 23]};
    plMuscle_pos={[3 9] [6 12] [15 21] [18 24]};
    sign_list=[-1 -1 1 -1];
    offs_list=[0 -0.5 -0.5 0];
    NSyn=length(SSyn_list);
    for iDOF=1:NDOFs
        indexDOF=strcmp(DOFNames,DOFNames_list{iDOF});
        jointID =MRS_eDot{1,1}.DatStore.IDinterp(fSel,indexDOF);

        indexMus=zeros(1,2);
        for iMus=1:2
            indexMus(1,iMus)=find(strcmp(MuscleNames,MuscleNames_list{iDOF}{iMus}));
        end

        name_list={'J_E' 'J_S'};
        % torques
        % subplot(2,6,plTorque_pos{iDOF})
        subplot(4,6,plTorque_pos{iDOF})
        plg(1)=plot(gait_cycle_sel,sign_list(iDOF).*jointID/subject_mass,'k','LineWidth',4,'DisplayName','ID'); hold on
        for iSyn=1:NSyn
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
            MActivation_N=Results_Baseline.Results.MActivation.genericMRS(indexMus(1,iMus),:); 
            plot(gait_cycle,MActivation_N,'Color','k','LineWidth',3); hold on %'#71797E'

            for iSyn=1:NSyn
                sSyn=SSyn_list(iSyn);
                Results_Bilevel =MRS_eDot{indexDOF,sSyn};

                MActivation_OPT=Results_Bilevel.Results.MActivation.genericMRS(indexMus(1,iMus),:);
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