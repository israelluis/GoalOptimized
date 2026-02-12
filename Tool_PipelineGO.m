clc;
computerPath  ='C:\Users\Israel Luis\Documents\GitHub\GoalOptimized';

% Input simulation
mySub_list=[1]; % list of subs 100 ITERS 5 NO
myDev_list=[1]; % list of devices 3 4 5
% mySyn_list=[0]; % list of synergies [4:6]
mySyn_list=[4:6]; % list of synergies [4:6]

% Total of conditions
[nSubs, nDevs, nSyns] = deal(length(mySub_list), length(myDev_list), length(mySyn_list));
MRS_list=cell(nSubs,nDevs,nSyns);

% condition
setCondition='getBilevel'; %runBaseline getBaseline runBilevel getBilevel

% track time
ExecutionTime_1=datetime('now');

% batch computation
for iSub=1:nSubs
    sSub=mySub_list(iSub);
    mySubject=['sub' num2str(sSub)];

    for iDev=1:nDevs
        sDev=myDev_list(iDev);
        myDevice=getDeviceSimple(sDev);

        for iSynConf=1:nSyns
            myOutName =['_RX' num2str(mySyn_list(iSynConf))]; %_R10 _RNO10 _RXNO10 _R1(X) _RX _RXe
            [MRS]=Tool_MainGO(computerPath,setCondition,mySubject,myDevice,myOutName);
           
            MRS_list(iSub,iDev,iSynConf)={MRS};

        end
    end
end

% compute total time
ExecutionTime_2=datetime('now');    sim_duration = ExecutionTime_2 - ExecutionTime_1;
disp(['Total simulation time: ' char(sim_duration)]);
%% Activations
for iSub=1:nSubs
    sSub=mySub_list(iSub);
    fig=figure(iSub);
    % fig=figure(2);
    tabgp = uitabgroup(fig, 'Position', [0.05 0.05 0.95 0.95]);
    fig.Name = ['now Subject ' num2str(sSub) ' - Synergy Analysis'];
    for iSynConf=1:nSyns
        tab = uitab(tabgp, 'Title', ['Synergy ' num2str(mySyn_list(iSynConf))]);
        axes('Parent', tab);
        MRS=MRS_list{iSub,1,iSynConf};
        for i=1:40
            subplot(5,8,i); hold on
            plot(MRS{1}.Results.GaitCycle,MRS{1}.Results.MActivation(i,1:end-1),'k','LineWidth',2);
            plot(MRS{2}.Results.GaitCycle,MRS{2}.Results.MActivation(i,1:end-1),'b','LineWidth',2);
            % plot(MRS{3}.Results.GaitCycle,MRS{3}.Results.MActivation(i,1:end-1),'r','LineWidth',2);
            ylim([0 1])
        end
    end
end
%% Synergies
nSubs=1;
for iSub=1:nSubs
    sSub=mySub_list(iSub);
    fig=figure(iSub);
    tabgp = uitabgroup(fig, 'Position', [0.05 0.05 0.95 0.95]);
    fig.Name = ['now Subject ' num2str(sSub) ' - Synergy Analysis'];
    for iDev=1:nDevs
        tab = uitab(tabgp, 'Title', ['Device ' num2str(iDev) dev_label{iDev}]);
        axes('Parent', tab);

        for iSynConf=1:nSyns
            MRS=MRS_list{iSub,iDev,iSynConf};
            N=MRS{2}.Misc.SynCon.N;

            for i=1:N
                subplot(3,6,i+6*(iSynConf-1)); hold on
                plot(MRS{2}.Results.SynergyControl.SynergyActivation(i,1:end-1),'b','LineWidth',2);
                % plot(MRS{3}.Results.SynergyControl.SynergyActivation(i,1:end-1),'r','LineWidth',2);
                ylim([0 1])
            end
        end
    end
end
%% No synergies
clc;
folderSpec   = 'JeSD';
specificFile = 'JeSD'; specificFile_dev=specificFile;
trialPath     ='v2_t1';

dev_label ={'AP(A)' 'HF(A)' 'HB(A)' 'KE(Q)'};
iter_label={'200' '200' '200' '100'};

sub_list=[1:5];     nSubs=length(sub_list);
dev_list=[1:4];       nDevs=length(dev_list);
syn_list=[0];
nSyns=length(syn_list);
MRS_JeD=cell(nSubs,nSyns,nDevs);
for iSub=1:nSubs
    sSub=sub_list(iSub);
    for iSyn=1:nSyns
        sSyn=syn_list(iSyn);
        for iDev=1:nDevs
            sDev=dev_list(iDev);
            MRS_JeD(iSub,iDev,iSyn)={load(fullfile(computerPath,'\ProjectResults\DSE',['sub' num2str(sSub)],trialPath,folderSpec,...
                ['MRSOptimal_R1' num2str(sSyn) '_eDot_MCLU24_' dev_label{sDev} '_iters' iter_label{sDev} 'Results.mat']))};
        end
    end
end
%% Torques
color_list={'#e86975' '#cc2525' '#78080d'};
nSubs=5;
for iSub=1:nSubs
    sSub=mySub_list(iSub);
    fig=figure(iSub); clf;

    fig.Name = ['Subject ' num2str(sSub)];
    tabgp = uitabgroup(fig, 'Position', [0.05 0.05 0.95 0.95]);

    mass=MRS_JeD{iSub,1,1}.Misc.subject_data.subject_mass;

    for iDev=1:nDevs
        tab = uitab(tabgp, 'Title', ['Device ' num2str(iDev) dev_label{iDev}]);
        axes('Parent', tab);

        subplot(1,1,1); hold on
        MRS_sel=MRS_JeD{iSub,iDev,1};
        plot(MRS_sel.Results.Device{1}.Assistance.Profile.GaitCycle,MRS_sel.Results.Device{1}.Assistance.Profile.Torque/mass,'color','k','LineWidth',2);
        for iSynConf=1:3

            MRS=MRS_list{iSub,iDev,iSynConf};
            plot(MRS{3}.Results.Device{1}.Assistance.Profile.GaitCycle,MRS{3}.Results.Device{1}.Assistance.Profile.Torque/mass,'color',color_list{iSynConf},'LineWidth',2);
        end
        axis([0 100 -0.2 1.8])
    end
end
%% Muscle activations
color_list={'#e86975' '#cc2525' '#78080d'};
nSubs=5;
for iSub=1:nSubs
    sSub=mySub_list(iSub);
    fig=figure(iSub); clf;

    tabgp = uitabgroup(fig, 'Position', [0.05 0.05 0.95 0.95]);
    fig.Name = ['Subject ' num2str(sSub)];

    mass=MRS_JeD{iSub,1,1}.Misc.subject_data.subject_mass;

    for iDev=1:nDevs
        tab = uitab(tabgp, 'Title', ['Device ' num2str(iDev) dev_list{iDev}]);
        axes('Parent', tab);

        for iMus=1:40
            subplot(5,8,iMus); hold on
            MRS_sel=MRS_JeD{iSub,iDev,1};
            plot(MRS_sel.Results.MActivation.genericMRS(iMus,:),'color','k','LineWidth',2);
            for iSynConf=1:3

                MRS=MRS_list{iSub,iDev,iSynConf};

                plot(MRS{2}.Results.MActivation(iMus,:),'color',color_list{iSynConf},'LineWidth',3);
                plot(MRS{3}.Results.MActivation(iMus,:),'color',color_list{iSynConf},'LineWidth',2,'LineStyle','-.');

            end
            title(MRS_sel.Results.MuscleNames{iMus})
            ylim([0 1])
        end
    end
end
%% Metabolic rates
assistiveGoal='eDot_MCLU24';
color_list={'#e86975' '#cc2525' '#78080d'};
condi_list={'-' '-.'};
nSubs=5;
for iSub=1:nSubs
    sSub=mySub_list(iSub);
    fig=figure(iSub); clf;

    tabgp = uitabgroup(fig, 'Position', [0.01 0.01 0.99 0.99]);
    fig.Name = ['Subject ' num2str(sSub)];

    mass=MRS_JeD{iSub,1,1}.Misc.subject_data.subject_mass;

    for iDev=1:nDevs
        tab = uitab(tabgp, 'Title', ['Device ' num2str(iDev) dev_label{iDev}]);
        axes('Parent', tab);

        for iSynConf=1:3
            % Metabolic rate time-series
            subplot(2,3,iSynConf)
            hold on;
            J_avg_all=zeros(2,1);
            for iCon=2:3
                MRS=MRS_list{iSub,iDev,iSynConf}{iCon};
                Misc=MRS.Misc;
                Results=MRS.Results;
                [J_avg,      J_TS,    ~]   = computeOuterLoopFunction(Misc,Results,assistiveGoal);

                GaitCycle=Results.GaitCycle(1+5:end-5);
                plot(GaitCycle,J_TS,'color',color_list{iSynConf},'LineWidth',3,'LineStyle',condi_list{iCon-1});
                J_avg_all(iCon-1)=J_avg;
            end
            ylim([0 7])
            ylabel('one leg metabolic rate [W/kg]')
            xlabel('gait cycle [%]')

            % Metabolic rate average
            x = ["unassisted" "assisted"];

            metReduction=(J_avg_all(2)-J_avg_all(1))/J_avg_all(1)*100;
            subplot(2,3,iSynConf+3)
            hb = bar(x,[J_avg_all(1) J_avg_all(2)], 'FaceColor', 'flat');
            hb.CData(1,:) = hex2rgb(color_list{iSynConf});
            hb.CData(2,:) = hex2rgb(color_list{iSynConf});

            text(2, J_avg_all(2) + 0.5, ...
            sprintf('%+.1f%%', metReduction), ...
            'HorizontalAlignment', 'center', 'FontSize', 15, 'Color', 'k', 'FontWeight', 'bold');

            ylim([0 8])
            ylabel('net metabolic rate [W/kg]')
            % 
        end
    end
end
%%[~,      J_normal_TS,    ~]   = computeOuterLoopFunction(Misc,Results_normal,assistiveGoal);
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