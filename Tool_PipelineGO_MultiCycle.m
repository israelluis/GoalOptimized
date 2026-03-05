clc;
computerPath  ='C:\Users\Israel Luis\Documents\GitHub\GoalOptimized';

% Input simulation
mySub_list  =[1]; % list of subs 100 ITERS 5 NO         2 3 4 5
myTrial_list=[1 2 3];
myDev_list  =[1 3 4 5]; % list of devices 3 4 5                 3 4 5
mySyn_list  =[0 6]; % list of synergies [4:6]
% mySyn_list  =[4:6]; % list of synergies [4:6] 6 0

% 
list_pip={'_RXCM' '_RXCM'}; %'_RXCM' '_RXCM' '_RXCMa'

% Total of conditions
[nSubs, nTrials, nDevs, nSyns] = deal(length(mySub_list), length(myTrial_list), length(myDev_list), length(mySyn_list));
MRS_list=cell(nSubs,nDevs,nSyns);
SMotion =cell(nTrials,1);

% condition
setCondition='getBilevel'; %runBaseline getBaseline runBilevel getBilevel

% assistive goal
setAssistiveGoal ='JRXN_knee';   % select goal. Options: eDot eDot_MCLU24 gasForces KJMusForces JRXN_knee JRXN_knee_par

% track time
ExecutionTime_1=datetime('now');

% batch computation
for iSub=1:nSubs
    sSub=mySub_list(iSub);
    mySubject=['sub' num2str(sSub)];

    SMotionInd=myTrial_list;

    for iDev=1:nDevs
        sDev=myDev_list(iDev);
        myDevice=getDeviceSimple(sDev);

        for iSynConf=1:nSyns
            % LAST myOutName-> '_RXCM' (only 6 0) & DirSpec='TM';
            % NEW  myOutName-> '_RXCMP' (only 6 0) & DirSpec='TMP'; (100) AP
            % NEW  myOutName-> '_RXCMa' (only 6 0) & DirSpec='TMP'; (200) AP
            % ONLY
            lockFile = 'computeOuterLoopFunction.lock'; delete(lockFile);

            myOutName =[list_pip{iSynConf} num2str(mySyn_list(iSynConf))]; %_R10 _RNO10 _RXNO10 _R1(X) _RX _RXe _RXC(26h)  
            [MRS]=Tool_MainGO_MultiCycle(computerPath,setCondition,setAssistiveGoal,mySubject,SMotionInd,myDevice,myOutName);
            
            MRS_list(iSub,iDev,iSynConf)={MRS};
            
        end
    end
    close all;
end

% compute total time
ExecutionTime_2=datetime('now');    sim_duration = ExecutionTime_2 - ExecutionTime_1;
disp(['Total simulation time: ' char(sim_duration)]);
%% Appendex MetCost
assistiveGoal         ='eDot_MCLU24'; 

MRS_list_eDotMuscles=cell(nSubs,nDevs,nSyns,3,3);
for iSub=1:nSubs
    for iDev=1:nDevs
        for iSynConf=1:nSyns
            MRS=MRS_list{iSub,iDev,iSynConf};
            for iSim=1:3
                Misc=MRS{iSim}.Misc;
                Results=MRS{iSim}.Results;
                for iTrial=1:nTrials
                    [J_bilevel_avg,      J_bilevel_TS,    J_bilevel_extra]   = computeOuterLoopFunction(Misc,Results(iTrial),assistiveGoal); % this is my starting point, the Edot at unassisted conditions.

                    MRS_list_eDotMuscles(iSub,iDev,iSim,iSynConf,iTrial)={J_bilevel_extra.eDot.Muscles};

                    % original_x = linspace(0, 1, length(J_bilevel_TS));
                    % new_x = linspace(0, 1, 101);
                    % J_bilevel_TS_int(iTrial, :) = interp1(original_x, J_bilevel_TS, new_x, 'linear');
                end
            end
        end
    end
end
%% MULTI-GAIT: Activations or MetRates
nSubs=5;
sDev= 2;

color_list= {'k' 'b' 'r'};
for iSub=1:nSubs
    sSub=mySub_list(iSub);
    fig=figure(iSub);
    tabgp = uitabgroup(fig, 'Position', [0.05 0.05 0.95 0.95]);
    fig.Name = ['now Subject ' num2str(sSub) ' - Synergy Analysis'];
    for iSynConf=1:nSyns
        for iTrial=1:nTrials
            tab = uitab(tabgp, 'Title', ['Synergy ' num2str(mySyn_list(iSynConf)) '& Trial' num2str(iTrial)]);
            axes('Parent', tab);
            MRS=MRS_list{iSub,sDev,iSynConf};

            for i=1:40
                subplot(5,8,i); hold on

                type=1;
                if type==1
                    for iSim=1:3
                        plot(MRS{iSim}.Results(iTrial).GaitCycle,MRS{iSim}.Results(iTrial).MActivation(i,1:end-1),color_list{iSim},'LineWidth',2);
                    end
                    ylim([0 1])

                else
                    for iSim=1:3
                        eDotMuscles=MRS_list_eDotMuscles{iSub,sDev,iSim,iSynConf,iTrial};
                        plot(eDotMuscles(i,:),color_list{iSim},'LineWidth',2);
                    end
                    % plot(MRS{1}.Results(iTrial).GaitCycle,MRS{1}.Results(iTrial).MActivation(i,1:end-1),'k','LineWidth',2);
                    % plot(MRS{2}.Results(iTrial).GaitCycle,MRS{2}.Results(iTrial).MActivation(i,1:end-1),'b','LineWidth',2);
                    % plot(MRS{3}.Results(iTrial).GaitCycle,MRS{3}.Results(iTrial).MActivation(i,1:end-1),'r','LineWidth',2);
                    ylim([0 2])
                end
            end
        end
    end
end
%% Synergies
nSubs=2;
for iSub=1:nSubs
    sSub=mySub_list(iSub);
    fig=figure(iSub);
    tabgp = uitabgroup(fig, 'Position', [0.05 0.05 0.95 0.95]);
    fig.Name = ['now Subject ' num2str(sSub) ' - Synergy Analysis'];

    for iTrial=1:nTrials
        tab = uitab(tabgp, 'Title', ['Trial' num2str(iTrial)]);
        axes('Parent', tab);

        for iSynConf=1:nSyns
            MRS=MRS_list{iSub,iDev,iSynConf};
            N=MRS{2}.Misc.SynCon.N;

            for i=1:N
                subplot(6,6,i+6*(iSynConf-1)); hold on
                plot(MRS{2}.Results(iTrial).GaitCycle,MRS{2}.Results(iTrial).SynergyControl.SynergyActivation(i,1:end-1),'b','LineWidth',2);
                % plot(MRS{3}.Results(iTrial).GaitCycle,MRS{3}.Results(iTrial).SynergyControl.SynergyActivation(i,1:end-1),'r','LineWidth',2);
                ylim([0 1])

                subplot(6,6,i+6*(iSynConf-1)+18); hold on
                bar(MRS{2}.Results(iTrial).SynergyControl.W(:,i));
                ylim([0 1])
            end
        end
    end
end

%% SUMMARY: Torques
range_min=[0 2];
range_min2=[0 1.5];
Yrange_list  ={range_min-0.05 range_min2-0.25 range_min-0.25 range_min2-0.25};
color_list={'k' 'b' '#78080d'}; %'#e86975' '#cc2525' '#78080d'
line_list ={'-' '-' ':'};
dev_label ={'Ank. Plantar.' 'Hip Flex.' 'Hip Abd.' 'Knee Ext.'};
nSubs=1;

data_length=size(TorqueHILO_list,2); gait_cycle=linspace(0,100,data_length);

fig=figure; clf; set(gcf,'color','w');
for iSub=1:nSubs
    sSub=mySub_list(iSub);

    simType=1;
    mass=MRS_list{iSub,1,1}{simType}.Misc.subject_data.subject_mass;

    simType=3;
    for iDev=1:nDevs

        sTrial=1;

        subplot(5,4,iDev+(iSub-1)*4); hold on

        plot(gait_cycle,TorqueHILO_list(iDev,:),'color','#097969','LineWidth',4);
        for iSynConf=1:2
            MRS=MRS_list{iSub,iDev,iSynConf};
            plot(MRS{simType}.Results(sTrial).Device{1}.Assistance.Profile.GaitCycle,MRS{simType}.Results(sTrial).Device{1}.Assistance.Profile.Torque/mass,...
                'color',color_list{iSynConf},'LineWidth',3,'LineStyle',line_list{iSynConf});
        end


        axis([0 100 Yrange_list{iDev}])
        xlabel('gait cycle [%]'); ylabel('torque [Nm/kg]');

        title(['Sub #' num2str(iSub) ' ' dev_label{iDev}])
    end

end
%% No synergies
clc;
folderSpec   = 'JeSD';
specificFile = 'JeSD'; specificFile_dev=specificFile;
trialPath     ='v2_t1';

dev_label ={'AP(A)' 'HF(A)' 'HB(A)' 'KE(Q)'};
iter_label={'200' '200' '200' '100'};

sub_list=[1:5];         nSubs=length(sub_list);
dev_list=[1 2 3 4];     nDevs=length(dev_list);
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
Dev2DOF_list=[1 3 4 2];
Yrange_list  ={[-0.5 2.0] [-1.25 1.25] [-0.5 2.0] [-1.25 1.25]};
plotT_list =[1 -1 1 1];
color_list={'k' '#78080d'}; %'#e86975' '#cc2525' '#78080d'
nSubs=5;
for iSub=1:nSubs
    sSub=mySub_list(iSub);
    fig=figure; clf;

    fig.Name = ['Subject ' num2str(sSub)];
    tabgp = uitabgroup(fig, 'Position', [0.05 0.05 0.95 0.95]);

    simType=1;
    mass=MRS_list{iSub,1,1}{simType}.Misc.subject_data.subject_mass;

    simType=3;
    for iDev=1:nDevs
        tab = uitab(tabgp, 'Title', ['Device ' num2str(iDev) dev_label{iDev}]);
        axes('Parent', tab);

        sTrial=1;

        % sTrial=1;
        subplot(1,1,sTrial); hold on
        for iSynConf=1:2
            MRS=MRS_list{iSub,iDev,iSynConf};
            plot(MRS{simType}.Results(sTrial).Device{1}.Assistance.Profile.GaitCycle,MRS{simType}.Results(sTrial).Device{1}.Assistance.Profile.Torque/mass,'color',color_list{iSynConf},'LineWidth',2);
        end
        axis([0 100 Yrange_list{iDev}])
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
        for iTrial=1:nTrials
            tab = uitab(tabgp, 'Title', ['Device ' num2str(iDev) dev_label{iDev} '& Trial' num2str(iTrial)]);
            axes('Parent', tab);


            for iSynConf=1:1
                % Metabolic rate time-series
                subplot(2,3,iSynConf)
                hold on;
                J_avg_all=zeros(2,1);
                for iCon=2:3
                    MRS=MRS_list{iSub,iDev,iSynConf}{iCon};
                    Misc=MRS.Misc;
                    Results=MRS.Results(iTrial);
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
end
%%[~,      J_normal_TS,    ~]   = computeOuterLoopFunction(Misc,Results_normal,assistiveGoal);
%%
%% COMPUTE HILO
clc;
DatStore_normal=MRS_list{1,1,1}{1}.DatStore(1);
Misc=MRS_list{1,1,1}{1}.Misc;
subject_mass=Misc.subject_data.subject_mass;

genVal=100;

Device{1}.Mode        = 'prescribed'; % opts: optimized and prescribed
Device{1}.MuscleGroup = {'ankle_angle' -1};
Device{1}.Type        = {'active' 'spline#N3'}; % opts: active, quasi-passive, passive, EMG-driven
Device{1}.Params      = [54.13 27.9 9.59 genVal]; %54.8 14.5 8.1 51.8
[assistanceInfo]=generateTorque(Device{1},DatStore_normal,Misc.time(1,:),Misc.extra_frames);
% Ankle plantarflexion
GaitCycle=assistanceInfo.Profile.GaitCycle;
TorqueHILO=assistanceInfo.Profile.Torque/genVal*0.68;

TorqueHILO_list=zeros(4,length(TorqueHILO));
TorqueHILO_list(1,:)=TorqueHILO;

% Knee extension
Device{1}.MuscleGroup = {'knee_angle' -1};
Device{1}.Type        = {'quasi-passive' 'clutchSpring'}; % opts: active, quasi-passive, passive, EMG-driven
Device{1}.Params      = [0.05  genVal]; %0.03
[assistanceInfo]=generateTorque(Device{1},DatStore_normal,Misc.time(1,:),Misc.extra_frames);
TorquePeak=max(assistanceInfo.Profile.Torque);
TorqueHILO=assistanceInfo.Profile.Torque/TorquePeak*0.15;
TorqueHILO_list(4,:)=TorqueHILO;

% Hip flexion
Device{1}.Type        = {'active' 'spline#N3'}; % opts: active, quasi-passive, passive, EMG-driven
Device{1}.Params      = [65.9 33.3 24.7 genVal]; %54.8 14.5 8.1 51.8
[assistanceInfo]=generateTorque(Device{1},DatStore_normal,Misc.time(1,:),Misc.extra_frames);
TorqueHILO=assistanceInfo.Profile.Torque/genVal*0.29;
TorqueHILO_list(2,:)=TorqueHILO;

% Hip abduction
Device{1}.Type        = {'active' 'spline#N3'}; % opts: active, quasi-passive, passive, EMG-driven
Device{1}.Params      = [48 15 15 genVal]; %54.8 14.5 8.1 51.8
[assistanceInfo]=generateTorque(Device{1},DatStore_normal,Misc.time(1,:),Misc.extra_frames);
TorqueHILO=assistanceInfo.Profile.Torque/genVal*0.15; %10.83/subject_mass
TorqueHILO_list(3,:)=TorqueHILO;
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
        Device{1}.MuscleGroup= {'ankle_angle' -1};
    case 2
        Device{1}.Mode       = 'prescribed'; % opts: optimized and prescribed
        Device{1}.Type       = {'active' 'spline#N3'};     % opts: active, quasi-passive, passive, EMG-driven
        Device{1}.MuscleGroup= {'knee_angle' -1};
    case 3
        Device{1}.Mode       = 'prescribed'; % opts: optimized and prescribed
        Device{1}.Type       = {'active' 'spline#N3'};     % opts: active, quasi-passive, passive, EMG-driven
        Device{1}.MuscleGroup= {'hip_flexion' 1};
    case 4
        Device{1}.Mode       = 'prescribed'; % opts: optimized and prescribed
        Device{1}.Type       = {'active' 'spline#N3'};     % opts: active, quasi-passive, passive, EMG-driven
        Device{1}.MuscleGroup= {'hip_adduction' -1};
    case 5
        Device{1}.Mode       = 'prescribed'; % opts: optimized and prescribed
        Device{1}.Type       = {'quasi-passive' 'clutchSpring'};     % opts: active, quasi-passive, passive, EMG-driven
        Device{1}.MuscleGroup= {'knee_angle' -1};
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