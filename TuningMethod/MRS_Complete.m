function [Results,DatStore,Misc] = MRS_Complete(Misc)
% -----------------------------------------------------------------------%
% INPUTS:
%           model_path: path to the .osim model
%           time: time window
%           OutPath: path to folder where results will be saved
%           Misc: structure with general input data (see manual for more details)
%
% OUTPUTS:
%           Results:    structure with outputs (states, controls, ...)
%           DatStore:   structure with data used for solving the optimal control problem
%           Misc:       structure with general input data (see manual for more details)
% -----------------------------------------------------------------------%
% This is made from "solveMuscleRedundancy_ExoCal"

% update default settings
Misc = DefaultSettings_upd(Misc);

% number of motion trials
Misc.nTrials = length(Misc.IKfile);

%check if we have to adapt the start and end time so that it corresponds to
%time frames in the IK solution
[time] = Check_TimeIndices(Misc,Misc.time);
Misc.time=time;

OutPath   =Misc.OutPath;
model_path=Misc.model_path;
%% Extract muscle information
% ----------------------------------------------------------------------- %
% Perform muscle analysis for the different selected trials
DatStore = struct;
for i = 1:Misc.nTrials
    % select the IK and ID file
    IK_path_trial = Misc.IKfile{i};
    ID_path_trial = Misc.IDfile{i};
    % Run muscle analysis
    % OPTION 0: Dont run and get the muscle analysis
    % OPTION 1: Run and get muscle analysis
    % OPTION 2: Dont run and get the muscle analysis from selected muscle analysis
    if Misc.GetAnalysis==0    
        OutPath_muscleAnalysis= OutPath;
        Misc.RunAnalysis      = 0;
    elseif Misc.GetAnalysis==1
        OutPath_muscleAnalysis= OutPath;
        Misc.RunAnalysis      = 1;
    elseif Misc.GetAnalysis==2
        Misc.RunAnalysis      = 0;
        dir_opt= [Misc.OutFolder '' Misc.OutFile(1:6) '_exoNN_simCAL_funACT_mMOO_\Results'];
        OutPath_muscleAnalysis= dir_opt;%obtain from SIMCAL
        if ~exist(OutPath,'dir'); mkdir(OutPath); end
    end
    MuscleAnalysisPath=fullfile(OutPath_muscleAnalysis,'MuscleAnalysis'); if ~exist(MuscleAnalysisPath,'dir'); mkdir(MuscleAnalysisPath); end
    if Misc.RunAnalysis
        disp('MuscleAnalysis Running .....');
        OpenSim_Muscle_Analysis(IK_path_trial,model_path,MuscleAnalysisPath,[time(i,1) time(i,end)],Misc.DofNames_Input)
        disp('MuscleAnalysis Finished');
    end
    Misc.MuscleAnalysisPath=MuscleAnalysisPath;
    
    % ----------------------------------------------------------------------- %
    % Extract muscle information -------------------------------------------- %
    % Get number of degrees of freedom (dofs), muscle-tendon lengths and moment
    % arms for the selected muscles.
    [~,Misc.trialName,~]=fileparts(IK_path_trial);
    if ~isfield(Misc,'MuscleNames_Input') || isempty(Misc.MuscleNames_Input)
        Misc=getMuscles4DOFS(Misc);
    end
    
    % get muscle geometry information
    [DatStore] = getMuscleInfo(IK_path_trial,ID_path_trial,Misc,DatStore,i);
    
    % display warnings in muscle selection
    Warnings_MuscleNames(DatStore,Misc,i);
    
    % get indexes of the muscles for which optimal fiber length, tendon stiffness are estimated
    [DatStore] = GetIndices_US(DatStore,Misc,i);
end

% set the tendon stiffness
if ~isfield(Misc,'kT') || isempty(Misc.kT)
    Misc.kT =ones(1,length(Misc.MuscleNames_Input)).*35;
end
if isfield(Misc,'Set_kT_ByName') && ~isempty(Misc.Set_kT_ByName)
    Misc= set_kT_ByName(Misc,DatStore);
end

% Shift tendon force-length curve as a function of the tendon stiffness
Misc.shift = getShift(Misc.kT);

% get the EMG information
[DatStore] = GetEMGInfo(Misc,DatStore);
[DatStore] = GetUSInfo(Misc,DatStore);

% get the number of muscles
NMuscles = length(DatStore(1).MuscleNames);

%% Static optimization
% ----------------------------------------------------------------------- %
% Solve the muscle redundancy problem using static optimization
% NOTE: We do not estimate any parameters here, but these results can serve as
% decent initial guess for the later dynamic optimization
% Extract the muscle-tendon properties
[Misc.params,Misc.lMo,Misc.lTs,Misc.FMo,Misc.alphao,Misc.kT,Misc.shift]=ReadMuscleParameters_upd(model_path,DatStore(1).MuscleNames,Misc);

% Static optimization using IPOPT solver (used as an initial guess)
for trial = 1:Misc.nTrials
    DatStore    = SolveStaticOptimization_IPOPT_CasADi(DatStore,Misc,trial,0);
end

%% Input activation and contraction dynamics
% ----------------------------------------------------------------------- %
tau_act = 0.015;    Misc.tauAct = tau_act * ones(NMuscles, 1);       % activation time constant (activation dynamics)
tau_deact = 0.06;   Misc.tauDeact = tau_deact * ones(NMuscles,1);  % deactivation time constant (activation dynamics)
Misc.b = 0.1;       % tanh coefficient for smooth activation dynamics

%% Descretisation
Misc.Mesh_Frequency=100; %200
% mesh descretisation
for trial = 1:Misc.nTrials
    t0 = DatStore(trial).time(1); tf = DatStore(trial).time(end);
    Mesh(trial).N = round((tf-t0)*Misc.Mesh_Frequency);
    Mesh(trial).step = (tf-t0)/Mesh(trial).N;
    Mesh(trial).t = t0:Mesh(trial).step:tf;
end

%% Evaluate splines at Mesh Points
% ----------------------------------------------------------------------- %
% Get IK, ID, muscle analysis and static opt information at mesh points

for trial = 1:Misc.nTrials
    % Discretization
    N = Mesh(trial).N;
    time_opt = Mesh(trial).t;
    % Spline approximation of muscle-tendon length (LMT), moment arms (MA) and inverse dynamic torques (ID)
    for dof = 1:DatStore(trial).nDOF
        for m = 1:NMuscles
            DatStore(trial).JointMASpline(dof).Muscle(m) = spline(DatStore(trial).time,squeeze(DatStore(trial).dM(:,dof,m)));
        end
        DatStore(trial).JointIDSpline(dof) = spline(DatStore(trial).time,DatStore(trial).T_exp(:,dof));
    end
    
    for m = 1:NMuscles
        DatStore(trial).LMTSpline(m) = spline(DatStore(trial).time,DatStore(trial).LMT(:,m));
    end
    
    % Evaluate LMT, VMT, MA and ID at optimization mesh
    DatStore(trial).LMTinterp = zeros(length(time_opt),NMuscles); % Muscle-tendon length
    for m = 1:NMuscles
        [DatStore(trial).LMTinterp(:,m),~,~] = SplineEval_ppuval(DatStore(trial).LMTSpline(m),time_opt,1);
    end
    DatStore(trial).MAinterp = zeros(length(time_opt),DatStore(trial).nDOF*NMuscles); % Moment arm
    DatStore(trial).IDinterp = zeros(length(time_opt),DatStore(trial).nDOF); % Inverse dynamic torque
    for dof = 1:DatStore(trial).nDOF
        for m = 1:NMuscles
            index_sel=(dof-1)*NMuscles+m;
            DatStore(trial).MAinterp(:,index_sel) = ppval(DatStore(trial).JointMASpline(dof).Muscle(m),time_opt);
        end
        DatStore(trial).IDinterp(:,dof) = ppval(DatStore(trial).JointIDSpline(dof),time_opt);
    end
    
    % Interpolate results of static optimization
    DatStore(trial).SoActInterp = interp1(DatStore(trial).time,DatStore(trial).SoAct,time_opt');
    DatStore(trial).SoRActInterp = interp1(DatStore(trial).time,DatStore(trial).SoRAct,time_opt');
    DatStore(trial).SoForceInterp = interp1(DatStore(trial).time,DatStore(trial).SoForce.*DatStore(trial).cos_alpha./Misc.FMo,time_opt);
    [~,DatStore(trial).lMtildeInterp ] = FiberLength_Ftilde(DatStore(trial).SoForceInterp,Misc.params,DatStore(trial).LMTinterp,Misc.kT,Misc.shift);
    DatStore(trial).vMtildeinterp = zeros(size(DatStore(trial).lMtildeInterp));
    for m = 1:NMuscles
        DatStore(trial).lMtildeSpline = spline(time_opt,DatStore(trial).lMtildeInterp(:,m));
        [~,DatStore(trial).vMtildeinterp_norm,~] = SplineEval_ppuval(DatStore(trial).lMtildeSpline,time_opt,1);
        DatStore(trial).vMtildeinterp(:,m) = DatStore(trial).vMtildeinterp_norm;
    end
    
    % interpolate the joint angles
    DatStore(trial).IKinterp = interp1(DatStore(trial).time,DatStore(trial).q_exp,time_opt');
end

%% Fiber length setup
% Obtain index and interpolate for muscle fiber calibration
% if Misc.Advance.TuningMethod_fiber==1
%     [US_data,nUSdata,ind_US,ind_US_AT,ind_USnone,USDigitalizedInterp]=fiberCalibrationSetup(DatStore,Misc,time_opt); % only works for one trial
% end
DatStore.Mesh=Mesh;
[Results,DatStore,Misc] = MRS_Formulate_and_Solve_NeuroCons(Misc,DatStore);
end


