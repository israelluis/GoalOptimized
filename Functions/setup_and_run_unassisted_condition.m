function [Misc,Results_normal,DatStore_normal] = setup_and_run_unassisted_condition(info)

SubjectSelection=info.SubjectSelection;
MotionSelection =info.MotionSelection;
currentFolder   =info.currentFolder;

ModelVariation   ='';      % none
FolderOutput     =fullfile('ProjectResults','Codesign','Pilot');

% read metadata
[Misc] = loadMetaData(currentFolder,SubjectSelection,MotionSelection);

% experimental data
Misc.IKfile = {fullfile(currentFolder,'Database',SubjectSelection,'IK',['IK_' SubjectSelection '_' MotionSelection '.mot'])};
Misc.IDfile = {fullfile(currentFolder,'Database',SubjectSelection,'ID',['ID_' SubjectSelection '_' MotionSelection '.sto'])};
Misc.model_path= fullfile(currentFolder,'Database',SubjectSelection,'model',['model_rajagopal2022_' SubjectSelection ModelVariation '.osim']);
Misc.extra_frames=5;

% input dofs - select the DOFs you want to include in the optimization
Misc.DofNames_Input={['ankle_angle_' Misc.gait_data.side_sel] ['knee_angle_' Misc.gait_data.side_sel] ['hip_flexion_' Misc.gait_data.side_sel] ['hip_adduction_' Misc.gait_data.side_sel] ['hip_rotation_' Misc.gait_data.side_sel]};

% run muscle analysis
Misc.GetAnalysis = 0;

% output path
Misc.SetupPath      = fullfile(currentFolder,'GenericSetups');
Misc.OutPath        = fullfile(currentFolder,FolderOutput,SubjectSelection,MotionSelection);

% settings for optimization within the same loop
Misc.Advance.AssistiveDevice   = 0; % exo disenabled
Misc.Advance.TuningMethod_fiber= 0; % no parameter tuning

% select tuned parameters. They need to be previously computed
variable_tuned  = load(fullfile('Database',SubjectSelection,'model',['tunedParameters_' SubjectSelection '.mat']));
Misc.param_label= load_MTU_names();       % here I name the convention name for MTU parameters
Misc.Param_read = {'selected'};
Misc.myParams   = variable_tuned.tunedParams;

% Setup names
muscleNames            = getMuscleNames('rajagopal');
Misc.MuscleNames_Input = appendSide(muscleNames, Misc.gait_data.side_sel);

% to name and save results
Misc.to_save_results = 1;
Misc.OutName= 'normal';

% Workflow setup
to_run_workFlow_Normal = 0; % run:1, read:0

if to_run_workFlow_Normal==1
    [Results_normal,DatStore_normal,Misc] = MRS_Complete(Misc);
else
    SavePath        = fullfile(currentFolder,FolderOutput,SubjectSelection,MotionSelection);
    path_dir_nor    = fullfile(SavePath, 'normalResults.mat');
    R_normalResult  = load(path_dir_nor);
    Results_normal  = R_normalResult.Results;
    DatStore_normal = R_normalResult.DatStore;
    Misc            = R_normalResult.Misc;
end
end


function [Misc] = loadMetaData(currentFolder,SubjectSelection,MotionSelection)
% read metadata
subInfo       =load(fullfile(currentFolder,'Database',SubjectSelection,'model',['subject_information_' SubjectSelection '.mat']));
subject_mass  =subInfo.subject_info.mass;   % [kg]
subject_height=subInfo.subject_info.height; % [cm]

gaitData      =load(fullfile(currentFolder,'Database',SubjectSelection,'gaitData',['gaitFeatureData_' MotionSelection '.mat']));
side_sel      =gaitData.gaitData.side; % either right or left leg side
toeOff_time   =gaitData.gaitData.toeOff;
speed         =gaitData.gaitData.speed; 
cadence       =gaitData.gaitData.cadence;

extLoadsInfo.fileName =fullfile(currentFolder,'Database',SubjectSelection,'extLoads',['data_'  MotionSelection '.mot']);
extLoadsInfo.setupName=fullfile(currentFolder,'Database',SubjectSelection,'extLoads',['setup_' MotionSelection '.xml']);


% loading MISC
% select the leg's side and toeOff event
Misc.gait_data.side_sel     = side_sel;    % [r or l]
Misc.gait_data.toeOff_time  = toeOff_time; % in seconds
Misc.gait_data.speed        = speed;       % in m/s
Misc.gait_data.cadence      = cadence*2;   % in #steps per minute

% subject info
Misc.subject_data.subject_mass   =subject_mass; % [kg]
Misc.subject_data.subject_height =subject_height; % [cm]

Misc.extLoadsInfo=extLoadsInfo;
Misc.MotionSelection=MotionSelection;
end

function muscleNames = getMuscleNames(model_name)
    if strcmp(model_name,'rajagopal')
       muscleNames= {'addbrev_' 'addlong_' 'addmagDist_' 'addmagIsch_' 'addmagMid_' 'addmagProx_' 'bflh_' ...
                          'bfsh_'    'edl_'     'ehl_'        'fdl_'        'fhl_'       'gaslat_'     'gasmed_' 'glmax1_' ...
                          'glmax2_'  'glmax3_'  'glmed1_'     'glmed2_'     'glmed3_'    'glmin1_'     'glmin2_' 'glmin3_' ...
                          'grac_'    'iliacus_' 'perbrev_'    'perlong_'    'piri_'      'psoas_'      'recfem_' 'sart_' ...
                          'semimem_' 'semiten_' 'soleus_'     'tfl_'        'tibant_'    'tibpost_'    'vasint_' 'vaslat_' 'vasmed_'};
    else
        error('Provide a valid model');
    end
end

function newNames = appendSide(muscleNames, side)
    % appendSide appends 'r' or 'l' to each muscle name based on input side
    % Input:
    %   muscleNames - cell array of strings ending in '_'
    %   side        - character or string, 'r' or 'l'
    % Output:
    %   newNames    - cell array of strings with 'r' or 'l' appended

    % % Validate input
    % if ~ismember(side, {'r', 'l'})
    %     error('Side must be either ''r'' or ''l''.');
    % end

    % Append side to each name
    newNames = cellfun(@(name) [name side], muscleNames, 'UniformOutput', false);
end