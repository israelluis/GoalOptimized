clc;
% conditions in the simulation
mySub_list=[1]; % list of subs
myDev_list=[1]; % list of devices
mySyn_list=[4 5 6]; % list of synergies

% total of conditions
nSubs=length(mySub_list);
nDevs=length(myDev_list);
nSyns=length(mySyn_list);

% track time
ExecutionTime_1=datetime('now');

% condition
myCondition='Bilevel'; %onlyBaseline

for iSub=1:nSubs
    sSub=mySub_list(iSub);
    mySubject=['sub' num2str(sSub)];
    for iDev=1:nDevs
        sDev=myDev_list(iDev);
        myDevice=getDeviceSimple(sDev);

        for iSynConf=1:length(mySyn_list)
            myOutName =['_R1' num2str(mySyn_list(iSynConf))];
            Tool_MainGO(mySubject,myCondition,myDevice,myOutName);
        end
    end
end

% compute total time
ExecutionTime_2=datetime('now'); duration_loop=seconds(duration(ExecutionTime_2-ExecutionTime_1));
disp(['Total simulation time: ' num2str(duration_loop,'%1.2f') ' secs']);
%%

allDevice=0;
if allDevice==1

% for iCon=1:length(condition_list)
for iDev=4:4
    Device_delivered=getDevice(iDev);
    Tool_MainGO(Device_delivered);

    lockFile = 'computeOuterLoopFunction.lock';
    delete(lockFile);
end
end

function [Device]=getDeviceSimple(iCon)
Misc.gait_data.side_sel='r';
switch iCon
    case 1
        Device{1}.Mode       = 'prescribed'; % opts: optimized and prescribed
        Device{1}.Type       = {'active' 'spline#N3'};     % opts: active, quasi-passive, passive, EMG-driven
        Device{1}.MuscleGroup= {['ankle_angle_' Misc.gait_data.side_sel] -1};
    case 2
        Device{1}.Mode       = 'prescribed'; % opts: optimized and prescribed
        Device{1}.Type       = {'active' 'spline#N3'};     % opts: active, quasi-passive, passive, EMG-driven
        Device{1}.MuscleGroup= {['knee_angle_' Misc.gait_data.side_sel] -1};
    case 3
        Device{1}.Mode       = 'prescribed'; % opts: optimized and prescribed
        Device{1}.Type       = {'active' 'spline#N3'};     % opts: active, quasi-passive, passive, EMG-driven
        Device{1}.MuscleGroup= {['hip_flexion_' Misc.gait_data.side_sel]  1};
    case 4
        Device{1}.Mode       = 'prescribed'; % opts: optimized and prescribed
        Device{1}.Type       = {'active' 'spline#N3'};     % opts: active, quasi-passive, passive, EMG-driven
        Device{1}.MuscleGroup= {['hip_adduction_' Misc.gait_data.side_sel] -1};
    case 5
        Device{1}.Mode       = 'prescribed'; % opts: optimized and prescribed
        Device{1}.Type       = {'quasi-passive' 'clutchSpring'};     % opts: active, quasi-passive, passive, EMG-driven
        Device{1}.MuscleGroup= {['knee_angle_' Misc.gait_data.side_sel] -1};
end
end


function [Device]=getDevice(iCon)
Misc.gait_data.side_sel='r';
switch iCon
    case 1
        Device{1}.Mode       = 'prescribed'; % opts: optimized and prescribed
        Device{1}.Type       = {'active' 'spline#N3'};     % opts: active, quasi-passive, passive, EMG-driven
        Device{1}.MuscleGroup= {['hip_flexion_' Misc.gait_data.side_sel]  1};

    case 2
        Device{1}.Mode       = 'prescribed'; % opts: optimized and prescribed
        Device{1}.Type       = {'active' 'spline#N3'};     % opts: active, quasi-passive, passive, EMG-driven
        Device{1}.MuscleGroup= {['hip_adduction_' Misc.gait_data.side_sel] -1};

    case 3
        Device{1}.Mode       = 'prescribed'; % opts: optimized and prescribed
        Device{1}.Type       = {'quasi-passive' 'clutchSpring'};     % opts: active, quasi-passive, passive, EMG-driven
        Device{1}.MuscleGroup= {['hip_flexion_' Misc.gait_data.side_sel]  1};

    case 4
        Device{1}.Mode       = 'prescribed'; % opts: optimized and prescribed
        Device{1}.Type       = {'quasi-passive' 'clutchSpring'};     % opts: active, quasi-passive, passive, EMG-driven
        Device{1}.MuscleGroup= {['hip_adduction_' Misc.gait_data.side_sel] -1};
    case 5
        Device{1}.Mode       = 'prescribed'; % opts: optimized and prescribed
        Device{1}.Type       = {'active' 'spline#N3'};     % opts: active, quasi-passive, passive, EMG-driven
        Device{1}.MuscleGroup= {['ankle_angle_' Misc.gait_data.side_sel] -1};
        Device{2}.Mode       = 'prescribed'; % opts: optimized and prescribed
        Device{2}.Type       = {'active' 'spline#N3'};     % opts: active, quasi-passive, passive, EMG-driven
        Device{2}.MuscleGroup= {['hip_adduction_' Misc.gait_data.side_sel] -1};
end
end
% Device{1}.Mode       = 'prescribed'; % opts: optimized and prescribed
% Device{1}.Type       = {'active' 'spline#N3'};     % opts: active, quasi-passive, passive, EMG-driven
% Device{1}.MuscleGroup= {['ankle_angle_' Misc.gait_data.side_sel] -1};
% Device{2}.Mode       = 'prescribed'; % opts: optimized and prescribed        
% Device{2}.Type       = {'active' 'spline#N3'};     % opts: active, quasi-passive, passive, EMG-driven
% Device{2}.MuscleGroup= {['knee_angle_' Misc.gait_data.side_sel] -1}; 

% Device{1}.Mode       = 'prescribed'; % opts: optimized and prescribed        
% Device{1}.Type       = {'quasi-passive' 'clutchSpring'};     % opts: active, quasi-passive, passive, EMG-driven
% Device{1}.MuscleGroup= {['ankle_angle_' Misc.gait_data.side_sel] -1};
% Device{2}.Mode       = 'prescribed'; % opts: optimized and prescribed        
% Device{2}.Type       = {'quasi-passive' 'clutchSpring'};    % opts: active, quasi-passive, passive, EMG-driven
% Device{2}.MuscleGroup= {['knee_angle_' Misc.gait_data.side_sel] -1};

% Device{1}.Mode       = 'prescribed'; % opts: optimized and prescribed        
% Device{1}.Type       = {'active' 'spline#N3'};     % opts: active, quasi-passive, passive, EMG-driven
% Device{1}.MuscleGroup= {['ankle_angle_' Misc.gait_data.side_sel] -1};
% Device{2}.Mode       = 'prescribed'; % opts: optimized and prescribed        
% Device{2}.Type       = {'quasi-passive' 'clutchSpring'};     % opts: active, quasi-passive, passive, EMG-driven
% Device{2}.MuscleGroup= {['knee_angle_' Misc.gait_data.side_sel] -1}; 

% Device{1}.Mode       = 'prescribed'; % opts: optimized and prescribed        
% Device{1}.Type       = {'quasi-passive' 'clutchSpring'};     % opts: active, quasi-passive, passive, EMG-driven
% Device{1}.MuscleGroup= {['ankle_angle_' Misc.gait_data.side_sel] -1};
% Device{2}.Mode       = 'prescribed'; % opts: optimized and prescribed        
% Device{2}.Type       = {'active' 'spline#N3'};     % opts: active, quasi-passive, passive, EMG-driven
% Device{2}.MuscleGroup= {['knee_angle_' Misc.gait_data.side_sel] -1}; 

% Device{1}.Mode       = 'prescribed'; % opts: optimized and prescribed        
% Device{1}.Type       = {'active' 'spline#N3'};     % opts: active, quasi-passive, passive, EMG-driven
% Device{1}.MuscleGroup= {['knee_angle_' Misc.gait_data.side_sel] -1}; 

% Device{1}.Mode       = 'prescribed'; % opts: optimized and prescribed        
% Device{1}.Type       = {'active' 'spline#N3'};     % opts: active, quasi-passive, passive, EMG-driven
% Device{1}.MuscleGroup= {['ankle_angle_' Misc.gait_data.side_sel] -1}; 
% 
% Device{1}.Mode       = 'prescribed'; % opts: optimized and prescribed        
% Device{1}.Type       = {'quasi-passive' 'clutchSpring'};     % opts: active, quasi-passive, passive, EMG-driven
% Device{1}.MuscleGroup= {['ankle_angle_' Misc.gait_data.side_sel] -1}; 

% Device{1}.Mode       = 'prescribed'; % opts: optimized and prescribed        
% Device{1}.Type       = {'quasi-passive' 'clutchSpring'};     % opts: active, quasi-passive, passive, EMG-driven
% Device{1}.MuscleGroup= {['knee_angle_' Misc.gait_data.side_sel] -1}; 