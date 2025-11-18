

condition_list=2:4;

for iCon=1:length(condition_list)
    Device_delivered=getDevice(iCon);
    Tool_MainGO(Device_delivered);

    lockFile = 'computeOuterLoopFunction.lock';
    delete(lockFile);
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