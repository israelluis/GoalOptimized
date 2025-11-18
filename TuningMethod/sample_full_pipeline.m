%% Load this folder and casadi
clc; close all; clear Misc
currentFolder = pwd;
addpath(genpath(currentFolder));
%% Input information
% paths and folders
ExampleFolder= 'DataExample';            % Select folder where... model, IK, and ID
DataFolder   = 'DataDigitalized';        % Select folder where... ultrasound and moment-angle relationship
OutputFolder = 'ResultsPipeline';        % Name your folder for results

% select the leg's side
Misc.side_sel      = 'r';                % Select leg side

% run muscle analysis
Misc.GetAnalysis = 1;                    % 0-> read a prior analysis, 1-> run muscle analysis
%% Setting paths and conventions (no need to change)
Misc.param_label=load_MTU_names();       % here I name the convention name for MTU parameters

Misc.OutPath   = fullfile(currentFolder,OutputFolder);  
Misc.model_path= fullfile(currentFolder,ExampleFolder,'Model.osim');
           
Misc.EMGfile= {fullfile(DataFolder,'emgFiles.mot')};        % named provided in DataFolder
Misc.USfile = {fullfile(DataFolder,'ultrasoundFiles.mot')}; % named provided in DataFolder
%% Workflow setup: Fiber lenghts
% Here we setup various workflows.

workFlow_gait_with_generic    =1; % simulate gait with generic values
workFlow_gait_tuning_FiberLen =1; % tuning parameters based on FIBER LENGTHS as in Paper
workFlow_gait_with_tunedFiber =1; % simulate gait with tuned fiber parameters only

workFlow_Mvsang_with_generic  =1; % simulate passive moments across joint angles with generic values
workFlow_Mvsang_tuning_PassM  =1; % tuning parameters based on PASSIVE MOMENT as in Paper

workFlow_gait_with_allTuned   =1; % simulate gait with all tuned parameters

%% TUNING - PART 1: Muscle fiber lengths based on ultrasound imaging data
Misc.IKfile = {fullfile(currentFolder,ExampleFolder,'IK.mot')};
Misc.IDfile = {fullfile(currentFolder,ExampleFolder,'ID.sto')};
Misc.DofNames_Input={['hip_flexion_' Misc.side_sel] ['hip_adduction_' Misc.side_sel] ['hip_rotation_' Misc.side_sel] ['knee_angle_' Misc.side_sel]  ['ankle_angle_' Misc.side_sel]}; 
Misc.time = [0 10];

% simulation with generic params & kT =15 (equivalent to 150 N/mm with Rajagoal et al. original values)
if workFlow_gait_with_generic==1
Misc.muscleFiberCal= 0; 
Misc.Param_read    = {'model'};
Misc.lMo_use_ward  = 0;
Misc.OutName       = 'Generic';

PF_kT=15; % plantarflexors
VA_kT=15; % vastii muscle
Misc.Set_kT_ByName = {['soleus_' Misc.side_sel],PF_kT; ['gasmed_' Misc.side_sel], PF_kT; ['gaslat_' Misc.side_sel],PF_kT;...
                      ['vasint_' Misc.side_sel],VA_kT; ['vaslat_' Misc.side_sel], VA_kT; ['vasmed_' Misc.side_sel],VA_kT};

[Results_gen,DatStore,Misc] = solveMuscleRedundancy_TV(Misc);
end

% simulation using digitalized muscle fiber lengths. Parameteres: lMo, lTs, and kT are design variables
if workFlow_gait_tuning_FiberLen==1
Misc.muscleFiberCal= 1; 
Misc.Param_read    = {'model'}; 
Misc.lMo_use_ward  = 1;
Misc.OutName       = 'ParmTun';

Misc.Set_kT_ByName ={}; % assign kT=35 for all muscles

[Results_tun,DatStore,Misc] = solveMuscleRedundancy_TV(Misc);
end

% simulation with previous calibrated parameters: lMo, lTs, and kT. This step does not use digitalized muscle fiber lengths
if workFlow_gait_with_tunedFiber==1
Misc.muscleFiberCal= 0;    
Misc.Param_read    = {'selected'}; % then "Misc.myParams" needs to be specified
Misc.lMo_use_ward  = 0;
Misc.OutName       = 'Validation';

Misc.myParams=labelParams(Results_tun.kT.calibratedMRS,Results_tun.params.calibratedMRS,Misc.param_label);

[Results_val,DatStore,Misc] = solveMuscleRedundancy_TV(Misc);
end

%% TUNING - PART 2: Passive moment-angle relationship based on force load and angles
% model:rajagopal: knee bending positive, 2392: knee bending negative angle
model_conf ={'rajagopal'};               Misc.model_conf=model_conf; 

ank_file =fullfile(DataFolder,['TrialAnk_' Misc.side_sel '_kne0hip0_kneBen15hip0_kneBen60hip0_' model_conf{:} 'Model']);
hip_file =fullfile(DataFolder,['TrialHip_' Misc.side_sel '_kneBen15ank0_kneBen60ank0_' model_conf{:} 'Model']);
kne_file =fullfile(DataFolder,['TrialKne_' Misc.side_sel '_ankDor20hip0_ankPla15hip0_ankDor20hipExt15_' model_conf{:} 'Model']);
neu_file =fullfile(DataFolder,['TrialNeu_' Misc.side_sel '_hipExt7hipAbd2kne0ankPla40_allModels']);

Misc.IKfile = {[ank_file '.mot'] [hip_file '.mot'] [kne_file '.mot'] [neu_file '.mot']};
Misc.IDfile = {[ank_file '.sto'] [hip_file '.sto'] [kne_file '.sto'] [neu_file '.sto']};

Misc.DofNames_Input={{['ankle_angle_' Misc.side_sel]};{['hip_flexion_' Misc.side_sel]};{['knee_angle_' Misc.side_sel]};...        % for passive forces
{['hip_flexion_' Misc.side_sel] ['hip_adduction_' Misc.side_sel] ['hip_rotation_' Misc.side_sel] ['knee_angle_' Misc.side_sel] ['ankle_angle_' Misc.side_sel]}};   % for neutral

time = [0 10; 0 10; 0 10; 0 10];

% IMPORTANT: Misc.model_path and Misc.OutPath as before
model_conf ={'rajagopal'};               Misc.model_conf=model_conf; 

Misc.semimem_adjustment=1; % change to "2" to have the lTs recomputed to be non negative
if Misc.semimem_adjustment==1 
    Misc.updParams.lTs.names      ={'semimem'};
    Misc.updParams.lTs.percentages=0.98; % i.e., 98% of its original value
end

if workFlow_Mvsang_with_generic==1
    Misc.Param_read={'selected'};
    Misc.Mode_opti ={'none'}; % no optimization
    Misc.OutName   ='GenericPassM';
    
    Misc.myParams=labelParams(Results_tun.kT.calibratedMRS,Results_tun.params.calibratedMRS,Misc.param_label);
    
    Misc.add_Passconstraints=0; % no need for this
    [Results_PassMGen,DatStore_PassCal,Misc_PassCal] = solveMuscleRedundancy_calibrationPassive([],time,[],Misc);
end

if workFlow_Mvsang_tuning_PassM==1
    Misc.Param_read={'selected'};
    Misc.Mode_opti ={'passive'}; % passive force-length curve optimization
    Misc.OutName   ='TunedPassM';
    
    Misc.myParams=labelParams(Results_tun.kT.calibratedMRS,Results_tun.params.calibratedMRS,Misc.param_label);
    
    Misc.add_Passconstraints=1;
    [Results_PassMTun,DatStore_PassCal,Misc_PassCal] = solveMuscleRedundancy_calibrationPassive([],time,[],Misc);
end

%% Plot moment-angle relationship
% Here I plot the results of moment-angle relationship between generic and
% tuned parameters. This graph helps to visualize the changes from tuning
nTrial_per_file =[3 2 3];
joint_order     =[1 3 2];
jointLeg_label   = {'ankle','knee','hip'};
description_label= {'knee flexion  0°' 'knee flexion 15°' 'knee flexion 60°';
                   {'ankle dorsiflexion 20°'; 'hip neutral 0°'} {'ankle plantarflexion 15°'; 'hip neutral 0°'} {'ankle dorsiflexion 20°'; 'hip extension 15°'};
                   'knee flexion 15°' 'knee flexion 60°' ' '};
jointAngle_label = {'ankle dorsiflexion [deg]' 'knee flexion [deg]' 'hip extension [deg]'};

lim_x_list={[-30 20]; [0 75]; [-15 40]};
lim_y_list={[-45 45]; [-55 35]; [-50 40]}; %[-45 20]; [-60 30]; [-50 40]
color_pass={'#191919' '#4682B4' '#cc0000'};
frame_per_trial=50;   
fig=figure(1); clf;  set(gcf,'color','w','Visible','on','Position',[50 50 1000 700]);
for joint_opt=1:3
    joint_sel=joint_order(joint_opt);
    for trial_sel=1:nTrial_per_file(joint_sel)
        
        sel_off=(trial_sel-1)*frame_per_trial;
        q_exp_array= DatStore_PassCal(joint_sel).q_exp(1+sel_off:frame_per_trial+sel_off);
        T_exp_array= DatStore_PassCal(joint_sel).T_exp(1+sel_off:frame_per_trial+sel_off);
        T_gen_array= Results_PassMGen.MMoment(joint_sel).genericMRS(1+sel_off:frame_per_trial+sel_off);
        T_cal_array= Results_PassMTun.MMoment(joint_sel).calibratedMRS(1+sel_off:frame_per_trial+sel_off);
        
        rcor_gen = corr(T_exp_array, T_gen_array);
        rmse_gen = sqrt(mean((T_exp_array - T_gen_array).^2));
       
        rcor_cal = corr(T_exp_array, T_cal_array);
        rmse_cal = sqrt(mean((T_exp_array - T_cal_array).^2));
       
        subplot(3,3,joint_opt+(trial_sel-1)*3)
        hold on;
        plot(q_exp_array,T_exp_array,'o','MarkerEdgeColor',"#000000",'MarkerFaceColor',color_pass{1},'markerSize',5);
        plot(q_exp_array,T_gen_array,'o','MarkerEdgeColor',"#000000",'MarkerFaceColor',color_pass{2},'markerSize',5);
        plot(q_exp_array,T_cal_array,'o','MarkerEdgeColor',"#000000",'MarkerFaceColor',color_pass{3},'markerSize',5);
       
        text(lim_x_list{joint_opt}(1)+7,lim_y_list{joint_opt}(2)-2,['GEN r=' num2str(rcor_gen,'%4.2f') ' RMSE=' num2str(rmse_gen,'%4.2f')],'fontSize',8,'color',color_pass{2},'FontWeight','bold','HorizontalAlignment','left')
        text(lim_x_list{joint_opt}(1)+7,lim_y_list{joint_opt}(2)-12,['TUN r=' num2str(rcor_cal,'%4.2f') ' RMSE=' num2str(rmse_cal,'%4.2f')],'fontSize',8,'color',color_pass{3},'FontWeight','bold','HorizontalAlignment','left')
       
        text(lim_x_list{joint_opt}(1)+2,lim_y_list{joint_opt}(1)+20,description_label{joint_opt+(trial_sel-1)*3},'fontSize',10)
        xlim(lim_x_list{joint_opt,:});
        ylim(lim_y_list{joint_opt,:});
       
        set(gca,'FontSize',12);
        xlabel(jointAngle_label{joint_opt},'FontSize',12);
        ylabel('moment [Nm]','FontSize',12)
       
       
        if trial_sel==1 && (joint_opt==1 || joint_opt==2 || joint_opt==3)
            title(jointLeg_label{joint_opt},'FontSize',20);
        end
    end
end

subplot(3,3,9)
hold on;
plot(0,0,'o','MarkerEdgeColor',"#000000",'MarkerFaceColor',color_pass{1},'markerSize',10);
plot(0,0,'o','MarkerEdgeColor',"#000000",'MarkerFaceColor',color_pass{2},'markerSize',10); hold on;
plot(0,0,'o','MarkerEdgeColor',"#000000",'MarkerFaceColor',color_pass{3},'markerSize',10); hold on;

axis([10,11,10,11]) %move dummy points out of view
legend('experimental','generic','tuned','location','southwest','fontSize',15); legend boxoff  
axis off %hide axis
%% Simulate walking with prior tuned lMo, lTs, kT, and passive force-length
Misc.IKfile = {fullfile(currentFolder,ExampleFolder,'IK.mot')};
Misc.IDfile = {fullfile(currentFolder,ExampleFolder,'ID.sto')};
Misc.DofNames_Input={['hip_flexion_' Misc.side_sel] ['hip_adduction_' Misc.side_sel] ['hip_rotation_' Misc.side_sel] ['knee_angle_' Misc.side_sel]  ['ankle_angle_' Misc.side_sel]}; 
Misc.time = [0 10];

if workFlow_gait_with_allTuned==1
Misc.muscleFiberCal= 0;    
Misc.Param_read    = {'selected'};
Misc.lMo_use_ward  = 0;
Misc.OutName       = 'Results_allTuned';

Misc.myParams=labelParams(Results_PassMTun.kT.calibratedMRS,Results_PassMTun.params.calibratedMRS,Misc.param_label);
[Results_all,DatStore,Misc] = solveMuscleRedundancy_TV(Misc);
end
%% Plot all results for walking
state_list={'MActivation' 'lMtildeopt'  'TForce'};
axis_list ={[0 1] [0.5 1.5] [0 1500]};
title_list={'muscle activations' 'normalized fiber lenghts' 'muscle-tendon forces'};
yaxis_list={'muscle activations [ ]' 'normalized lengths [ ]' 'forces [N]'};

MuscleNames= DatStore.MuscleNames;
NMuscle    = length(MuscleNames);
getSide    = MuscleNames{1}(end-1:end);
selectedMus_noSide=strrep(MuscleNames,getSide,'');

if ~exist('EMG_data','var') || ~exist('US_data','var')  
    EMG_data = ReadMotFile(Misc.EMGfile{1});
    US_data  = ReadMotFile(Misc.USfile{1});
end

extra_frames=5;

for state_opt=1:length(state_list)
    fig=figure(1+state_opt); clf;  set(gcf,'color','w','Visible','on','Position',[50 50 1000 700]);

    for mus_sel=1:NMuscle
        muscle_name=MuscleNames{mus_sel};
        state_Normal=Results_gen.(state_list{state_opt}).genericMRS(mus_sel,1+extra_frames:end-extra_frames);
        state_Val   =Results_val.(state_list{state_opt}).genericMRS(mus_sel,1+extra_frames:end-extra_frames);
        state_All   =Results_all.(state_list{state_opt}).genericMRS(mus_sel,1+extra_frames:end-extra_frames);

        x_axis      =linspace(0,100,length(state_Normal)); % or time as Results_gen.Time.genericMRS(1+extra_frames:end-extra_frames);

        subplot(5,8,mus_sel)
        hold on;
        plot(x_axis,state_Normal,'color','g','lineWidth',3,'lineStyle','-');
        plot(x_axis,state_Val,'color','b','lineWidth',3,'lineStyle','--');
        plot(x_axis,state_All,'color','r','lineWidth',3,'lineStyle','-.');
        ylim(axis_list{state_opt});
        set(gca,'FontSize',10);
        title(muscle_name)
    end
    sgtitle(title_list(state_opt),'fontSize',20)

    han=axes(fig,'visible','off');
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    label_y = ylabel(han,yaxis_list{state_opt},'FontSize',17); label_y.Position(1) = -0.05; label_y.Position(2) = 0.5;
    label_x = xlabel(han,'gait cycle [%]'     ,'FontSize',17); label_x.Position(1) =   0.5; label_x.Position(2) = -0.05;
end