%% Run sample, self-selected parameters / spline controllers
to_run_example=0;

if to_run_example==1
% run muscle analysis and enable assistive moment
Misc.GetAnalysis = 0;
Misc.Advance.AssistiveDevice  = 1;

% NOTE: Available assistance
% anklePlantarflexion: ankle_angle_   -1
% kneeExtension:       knee_angle_    -1
% hipFlexion:          hip_flexion_    1
% hipAbduction:        hip_adduction_ -1

clear Device
select_simulation_type='active spline#N3';
if strcmp(select_simulation_type,'active spline#N3')
    Device{1}.Mode        = 'prescribed'; % opts: optimized and prescribed
    Device{1}.MuscleGroup = {['ankle_angle_' Misc.gait_data.side_sel] -1};
    Device{1}.Type        = {'active' 'spline#N3'}; % opts: active, quasi-passive, passive, EMG-driven
    Device{1}.Params      = [54.1 22 8.6 77.9]; %54.8 14.5 8.1 51.8
    [assistanceInfo]=generateTorque(Device{1},DatStore,Misc.time,Misc.extra_frames);

    Device{1}.Assistance = assistanceInfo;
elseif strcmp(select_simulation_type,'quasi-passive clutchSpring')
    % generate DeviceSetup
    Device{1}.Mode        = 'prescribed';    % opts: optimized and prescribed
    % Device{1}.MuscleGroup = {['hip_adduction_' Misc.gait_data.side_sel] -1}; 
    % Device{1}.MuscleGroup = {['hip_flexion_' Misc.gait_data.side_sel] 1}; 
    Device{1}.MuscleGroup = {['knee_angle_' Misc.gait_data.side_sel] -1};    
    % Device{1}.MuscleGroup = {['ankle_angle_' Misc.gait_data.side_sel] -1};  
    Device{1}.Type        = {'quasi-passive' 'clutchSpring'}; % opts: active, quasi-passive, passive, EMG-driven
    Device{1}.Params      = [2.99 2.872];
    [assistanceInfo]=generateTorque(Device{1},DatStore,Misc.time,Misc.extra_frames);

    Device{1}.Assistance = assistanceInfo;
    Misc.Device=Device;
elseif strcmp(select_simulation_type,'double')
    % DEVICE 1
    Device{1}.Mode        = 'prescribed'; % opts: optimized and prescribed
    Device{1}.MuscleGroup = {['knee_angle_' Misc.gait_data.side_sel] -1};
    Device{1}.Type        = {'active' 'spline#N3'}; % opts: active, quasi-passive, passive, EMG-driven
    Device{1}.Params      = [12.4 6.3 24 67.3];
    [assistanceInfo]=generateTorque(Device{1},DatStore,Misc.time,Misc.extra_frames);

    Device{1}.Assistance = assistanceInfo;

    % DEVICE 2
    Device{2}.Mode        = 'prescribed'; % opts: optimized and prescribed
    Device{2}.MuscleGroup = {['ankle_angle_' Misc.gait_data.side_sel] -1};
    Device{2}.Type        = {'active' 'spline#N3'}; % opts: active, quasi-passive, passive, EMG-driven
    Device{2}.Params      = [39.2 20.2 20.9 45.1];
    [assistanceInfo]=generateTorque(Device{2},DatStore,Misc.time,Misc.extra_frames);

    Device{2}.Assistance = assistanceInfo;    
end
Misc.Device=Device;

% to name and save results
Misc.to_save_results= 0;
Misc.OutName= 'example';

% time performance of a single loop
ExecutionTime_1=datetime('now');
[Results_assisted,~,Misc]= MRS_Formulate_and_Solve_NeuroCons(Misc,DatStore);
ExecutionTime_2=datetime('now'); duration_loop=seconds(duration(ExecutionTime_2-ExecutionTime_1));
disp(['simulation took ' num2str(duration_loop,'%1.2f') ' secs']);

[J_example_avg,J_example_TS,J_example_extra] = computeOuterLoopFunction(Misc,Results_assisted,assistiveGoal);
J_value=(J_example_avg-J_baseline_avg)/J_baseline_avg*100;
disp(['goal change: ' num2str(J_value,'%+1.2f') '%'])

% Check results
to_plot=1;
if to_plot==1
    figure; clf; % muscle activations
    for i=1:40
        subplot(5,8,i)
        hold on
        plot(Results_normal.MActivation(i,:),'k')
        plot(Results_baseline.MActivation(i,:),'--g')
        plot(Results_assisted.MActivation(i,:),':r')
        ylim([0 1])
        title(Results_assisted.MuscleNames{i})

    end
    sgtitle('muscle activations')

    figure; clf; % device torque
    nDevs=length(Device);
    norMass=Misc.subject_data.subject_mass;
    for iDev=1:nDevs
        gaitCycle=Device{iDev}.Assistance.Profile.GaitCycle;
        iDOF=strcmp(DatStore.DOFNames,Device{iDev}.MuscleGroup{1});
        sID =DatStore.IDinterp(:,iDOF)*Device{iDev}.MuscleGroup{2};
        Torque=Device{iDev}.Assistance.Profile.Torque;
        subplot(nDevs,1,iDev)
        plot(gaitCycle,sID/norMass,'k'); hold on
        plot(gaitCycle,Torque/norMass,'r')
        title([Device{iDev}.MuscleGroup{1}],'Interpreter','none')
    end
    sgtitle('assistive torque')

    figure; clf;
    plot(J_baseline_TS,'k','LineWidth',2); hold on
    plot(J_example_TS,'r','LineWidth',2)
    xlabel('gait cycle[%]'); ylabel([J_example_extra.label ' [' J_example_extra.unit ']']);
    title([ 'baseline (avg value) = ' num2str(J_baseline_avg,'%1.1f') ' ' J_example_extra.unit ...
           ' assisted (avg value) = ' num2str(J_example_avg,'%1.1f')  ' ' J_example_extra.unit ...
           ' - change: ' num2str(J_value,'%+1.1f') '%']);
end
end