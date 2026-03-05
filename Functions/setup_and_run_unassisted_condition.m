function [Misc,Results_normal,DatStore_normal] = setup_and_run_unassisted_condition(info,DirF,init_data)

SubjectSelection= info.SubjectSelection;
MotionSelection = info.MotionSelection;
BaseFolder      = info.BaseFolder;
ProjectFolder   = info.ProjectFolder;

BaseDir         = DirF.select_folder_N1;
ModelVariation   ='';      % none

% initialization
if init_data==1;     GetAnalysis=1; setWorkflow=1;
elseif init_data==0; GetAnalysis=0; setWorkflow=0;
end
% GetAnalysis=0; setWorkflow=1;

% read metadata
[Misc] = loadMetaData(BaseFolder,SubjectSelection,MotionSelection);

% number of gait cycles
nGCs=numel(MotionSelection);
if nGCs>1;    MotionSelection_label='multiCycle'; else; MotionSelection_label = info.MotionSelection;
end

% model
Misc.model_path= fullfile(BaseFolder,'Database',SubjectSelection,'model',['model_rajagopal2022_' SubjectSelection ModelVariation '.osim']);

% joint kinematics, kinetics, and DOFs
for iGC=1:nGCs
    LegSide=Misc.gait_data(iGC).side_sel;

    Misc.IKfile(1,iGC) = {fullfile(BaseFolder,'Database',SubjectSelection,'IK',['IK_' SubjectSelection '_' MotionSelection{iGC} '.mot'])};
    Misc.IDfile(1,iGC) = {fullfile(BaseFolder,'Database',SubjectSelection,'ID',['ID_' SubjectSelection '_' MotionSelection{iGC} '.sto'])};
    % input dofs - select the DOFs you want to include in the optimization
    Misc.DofNames_Input(1,iGC)={{['ankle_angle_' LegSide] ['knee_angle_' LegSide] ['hip_flexion_' LegSide] ['hip_adduction_' LegSide] ['hip_rotation_' LegSide]}};
end

Misc.extra_frames=5;

% run muscle analysis
Misc.GetAnalysis = GetAnalysis;

% output path
Misc.SetupPath      = fullfile(BaseFolder,'GenericSetups');
Misc.TrialFolder    = fullfile(ProjectFolder,SubjectSelection,MotionSelection_label);
Misc.OutPath        = fullfile(BaseFolder,Misc.TrialFolder,BaseDir);

% settings for optimization within the same loop
Misc.Advance.AssistiveDevice   = 0; % No device
Misc.Advance.TuningMethod_fiber= 0; % No tuning

% select tuned parameters. They need to be previously computed
variable_tuned  = load(fullfile('Database',SubjectSelection,'model',['tunedParameters_' SubjectSelection '.mat']));
Misc.param_label= load_MTU_names();       % here I name the convention name for MTU parameters
Misc.Param_read = {'selected'};
Misc.myParams   = variable_tuned.tunedParams;

% setup names
muscleNames            = getMuscleNames('rajagopal');
for iGC=1:nGCs
    Misc.MuscleNames_Input(1,iGC) = {appendSide(muscleNames, Misc.gait_data(iGC).side_sel)};
end

% to name and save results
OutName_sel=BaseDir;
Misc.to_save_results = 1;
Misc.OutName= OutName_sel;

% Workflow setup
to_run_workFlow_Normal = setWorkflow; % run:1, read:0

if to_run_workFlow_Normal==1
    [Results_normal,DatStore_normal,Misc] = MRS_Complete(Misc);
else
    LoadPath=fullfile(BaseFolder,Misc.TrialFolder,BaseDir);
    MRS_normal      = load(fullfile(LoadPath, [OutName_sel 'Results.mat']));
    Results_normal  = MRS_normal.Results;
    DatStore_normal = MRS_normal.DatStore;
    Misc            = MRS_normal.Misc;

    Misc.TrialFolder    = fullfile(ProjectFolder,SubjectSelection,MotionSelection_label);
end
end


function [Misc] = loadMetaData(BaseFolder,SubjectSelection,MotionSelection)
nGCs=numel(MotionSelection);

% read metadata
subInfo       =load(fullfile(BaseFolder,'Database',SubjectSelection,'model',['subject_information_' SubjectSelection '.mat']));

% subject info
Misc.subject_data.subject_mass   =subInfo.subject_info.mass;   % [kg]
Misc.subject_data.subject_height =subInfo.subject_info.height; % [cm]

for iGC=1:nGCs
    Misc.MotionSelection(1,iGC)=MotionSelection(iGC);

    extLoadsInfo.fileName(1,iGC) ={fullfile(BaseFolder,'Database',SubjectSelection,'extLoads',['data_'  MotionSelection{iGC} '.mot'])};
    extLoadsInfo.setupName(1,iGC)={fullfile(BaseFolder,'Database',SubjectSelection,'extLoads',['setup_' MotionSelection{iGC} '.xml'])};
    Misc.extLoadsInfo=extLoadsInfo;

    gaitData      =load(fullfile(BaseFolder,'Database',SubjectSelection,'gaitData',['gaitFeatureData_' MotionSelection{iGC} '.mat']));

    % loading MISC
    % select the leg's side and toeOff event
    Misc.gait_data(1,iGC).side_sel     = gaitData.gaitData.side;    % [r or l]
    Misc.gait_data(1,iGC).toeOff_time  = gaitData.gaitData.toeOff; % in seconds
    Misc.gait_data(1,iGC).speed        = gaitData.gaitData.speed;       % in m/s
    Misc.gait_data(1,iGC).cadence      = gaitData.gaitData.cadence*2;   % in #steps per minute
end

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