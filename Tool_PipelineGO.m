clc;
computerPath  ='C:\Users\movea\Documents\GitHub\GoalOptimized';

% conditions in the simulation
mySub_list=[5]; % list of subs
myDev_list=[1 4 5]; % list of devices 3
mySyn_list=[4:6]; % list of synergies [4:6]

% total of conditions
nSubs=length(mySub_list);
nDevs=length(myDev_list);
nSyns=length(mySyn_list);

% condition
setCondition='getBilevel'; %runBaseline getBaseline runBilevel getBilevel

% track time
ExecutionTime_1=datetime('now');

% MRS_list=cell(3,1);

% batch computation
for iSub=1:nSubs
    sSub=mySub_list(iSub);
    mySubject=['sub' num2str(sSub)];

    % fig=figure(iSub); 
    % tabgp = uitabgroup(fig, 'Position', [0.05 0.05 0.95 0.95]);
    % fig.Name = ['Subject ' num2str(sSub) ' - Synergy Analysis'];

    for iDev=1:nDevs
        sDev=myDev_list(iDev);
        myDevice=getDeviceSimple(sDev);

        for iSynConf=1:nSyns
            myOutName =['_R1' num2str(mySyn_list(iSynConf))];
            [MRS]=Tool_MainGO(computerPath,setCondition,mySubject,myDevice,myOutName);
           
            % MRS_list(iSynConf)={MRS};
            % tab = uitab(tabgp, 'Title', ['Synergy ' num2str(mySyn_list(iSynConf))]);
            % axes('Parent', tab);

            % M=MRS{2}.Misc;
            % for i=1:M.SynCon.N; subplot(2,M.SynCon.N,i); bar(M.SynCon.W(:, i)); ylim([0 1]); end
            % for i=1:M.SynCon.N; subplot(2,M.SynCon.N,i+M.SynCon.N); plot(M.SynCon.H(i, :)); ylim([0 1]); end

        end
    end
end

% compute total time
ExecutionTime_2=datetime('now');    sim_duration = ExecutionTime_2 - ExecutionTime_1;
disp(['Total simulation time: ' char(sim_duration)]);
%%
for iSub=1:nSubs
    sSub=mySub_list(iSub);
    fig=figure(iSub);
    tabgp = uitabgroup(fig, 'Position', [0.05 0.05 0.95 0.95]);
    fig.Name = ['Subject ' num2str(sSub) ' - Synergy Analysis'];
    for iSynConf=1:nSyns
        tab = uitab(tabgp, 'Title', ['Synergy ' num2str(mySyn_list(iSynConf))]);
        axes('Parent', tab);
        MRS=MRS_list{iSynConf};
        for i=1:40
            subplot(5,8,i); hold on
            plot(MRS{1}.Results.GaitCycle,MRS{1}.Results.MActivation(i,1:end-1),'k','LineWidth',2);
            plot(MRS{2}.Results.GaitCycle,MRS{2}.Results.MActivation(i,1:end-1),'b','LineWidth',2);
            plot(MRS{3}.Results.GaitCycle,MRS{3}.Results.MActivation(i,1:end-1),'r','LineWidth',2);
            ylim([0 1])
        end
    end
end
%%
clc;
folderSpec   = 'JeSD';
specificFile = 'JeSD'; specificFile_dev=specificFile;
trialPath     ='v2_t1';

dev_list ={'HF(A)'};
iter_list={'200'};

sub_list=5:5;     nSubs=length(sub_list);

syn_list=[0];
nSyns=length(syn_list);
nDevs=length(dev_list);
MRS_JeSD=cell(nSubs,nSyns,nDevs);
for iSub=1:nSubs
    sSub=sub_list(iSub);
    for iSyn=1:nSyns
        sSyn=syn_list(iSyn);
        for iDev=1:nDevs
            MRS_JeSD(iSub,iSyn,iDev)={load(fullfile(computerPath,'\ProjectResults\DSE',['sub' num2str(sSub)],trialPath,folderSpec,...
                ['MRSOptimal_R1' num2str(sSyn) '_eDot_MCLU24_' dev_list{iDev} '_iters' iter_list{iDev} 'Results.mat']))};
        end
    end
end
%%
color_list={'#e86975' '#cc2525' '#78080d'};
for iSub=1:nSubs
    sSub=mySub_list(iSub);
    fig=figure(iSub); clf;
    % tabgp = uitabgroup(fig, 'Position', [0.05 0.05 0.95 0.95]);
    % fig.Name = ['Subject ' num2str(sSub) ' - Synergy Analysis'];


    mass=MRS_JeSD{iSub,1,1}.Misc.subject_data.subject_mass;
    subplot(1,1,1); hold on

    plot(MRS_JeSD{iSub,1,1}.Results.Device{1}.Assistance.Profile.GaitCycle,MRS_JeSD{iSub,1,1}.Results.Device{1}.Assistance.Profile.Torque/mass,'color','k','LineWidth',2);
    for iSynConf=1:3
        % tab = uitab(tabgp, 'Title', ['Synergy ' num2str(mySyn_list(iSynConf))]);
        % axes('Parent', tab);
        MRS=MRS_list{iSynConf};

        plot(MRS{3}.Results.Device{1}.Assistance.Profile.GaitCycle,MRS{3}.Results.Device{1}.Assistance.Profile.Torque/mass,'color',color_list{iSynConf},'LineWidth',2);

    end
end
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
% this assign device conf without (leg) side
switch iCon
    case 1
        Device{1}.Mode       = 'prescribed'; % opts: optimized and prescribed
        Device{1}.Type       = {'active' 'spline#N3'};     % opts: active, quasi-passive, passive, EMG-driven
        Device{1}.MuscleGroup= {'ankle_angle_' -1};
    case 2
        Device{1}.Mode       = 'prescribed'; % opts: optimized and prescribed
        Device{1}.Type       = {'active' 'spline#N3'};     % opts: active, quasi-passive, passive, EMG-driven
        Device{1}.MuscleGroup= {'knee_angle_' -1};
    case 3
        Device{1}.Mode       = 'prescribed'; % opts: optimized and prescribed
        Device{1}.Type       = {'active' 'spline#N3'};     % opts: active, quasi-passive, passive, EMG-driven
        Device{1}.MuscleGroup= {'hip_flexion_' 1};
    case 4
        Device{1}.Mode       = 'prescribed'; % opts: optimized and prescribed
        Device{1}.Type       = {'active' 'spline#N3'};     % opts: active, quasi-passive, passive, EMG-driven
        Device{1}.MuscleGroup= {'hip_adduction_' -1};
    case 5
        Device{1}.Mode       = 'prescribed'; % opts: optimized and prescribed
        Device{1}.Type       = {'quasi-passive' 'clutchSpring'};     % opts: active, quasi-passive, passive, EMG-driven
        Device{1}.MuscleGroup= {'knee_angle_' -1};
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