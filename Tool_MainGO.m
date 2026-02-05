function [MRS]=Tool_MainGO(computerPath,Scondition,Ssubject,SDevice,SOutName)
addpath(genpath(computerPath));
%% Input information
% select subject and motion
info.SubjectSelection = Ssubject;       % select subject e.g., 'sub1'
info.MotionSelection  ='v2_t1';         % select velocity & trial. Options: v(velocity)== 1[slow], 2[normal], 3[fast]) _ t (trial)== 1st, 2nd, or 3rd trial
info.BaseFolder       = computerPath;
info.ProjectFolder    = fullfile('ProjectResults','DSE'); %fullfile('ProjectResults','Codesign','Pilot')

% folder for outputs
DirName='Je'; DirSpec='T';
DirF.select_folder_N1 = DirName;                 % specify folder for simulation based on minimal muscle effort
DirF.select_folder_N2 =[DirName DirSpec 'S'];    % specify folder for simulation based on minimal muscle effort with synergies
DirF.select_folder_N3 =[DirName DirSpec 'SD'];  % specify folder for simulation based on "select_folder_N2" with assistive devices

% create folders
for i=1:numel(fieldnames(DirF))
    dirFile= fullfile(info.BaseFolder,info.ProjectFolder,info.SubjectSelection,info.MotionSelection,[DirF.(['select_folder_N' num2str(i)])]);
    if ~exist(dirFile, 'dir'); mkdir(dirFile); end
end

% setup unassisted/baseline condition
init_data=0;
[Misc,Results_normal,DatStore]=setup_and_run_unassisted_condition(info,DirF,init_data);

%% Synergy setup
if strcmp(Scondition, 'runBaseline')
    sSyn =[];                          conditionAnalysis='all';
elseif strcmp(Scondition, 'getBaseline') || strcmp(Scondition, 'runBilevel') || strcmp(Scondition, 'getBilevel')
    sSyn=str2double(SOutName(end));    conditionAnalysis='load'; 
end
infoSyn.sSyn = sSyn;
%% Synergy computation
[Results_baseline,Misc]=formulation_with_informed_synergist(computerPath,Misc,Results_normal,DatStore,DirF,conditionAnalysis,infoSyn);

% condition loading
if strcmp(Scondition, 'runBaseline')
    MRS{1}.Results=[];               MRS{1}.DatStore=[];       MRS{1}.Misc=[];
elseif strcmp(Scondition, 'getBaseline') || strcmp(Scondition, 'runBilevel') || strcmp(Scondition, 'getBilevel')
    MRS{1}.Results=Results_normal;   MRS{1}.DatStore=DatStore; MRS{1}.Misc=Misc;
    MRS{2}.Results=Results_baseline; MRS{2}.DatStore=DatStore; MRS{2}.Misc=Misc;
end

% condition stopping
if strcmp(Scondition, 'getBaseline')  || strcmp(Scondition, 'runBaseline')
    disp(['Condition:' Scondition '. Stop after getting synergy dimension N' num2str(sSyn) '.']);
    return;  % Exit function early
end
%% Compute and plot assistiveGoal in unassisted/baseline condition
assistiveGoal         ='eDot_MCLU24';   % select goal. Options: eDot eDot_MCLU24 gasForces KJMusForces RJXN_knee RJXN_knee_par

[Misc]=setup_for_JRXN(Misc);
[~,      J_normal_TS,    ~]   = computeOuterLoopFunction(Misc,Results_normal,assistiveGoal); % this is my starting point, the Edot at unassisted conditions.
[J_baseline_avg,  J_baseline_TS,  J_baseline_extra]   = computeOuterLoopFunction(Misc,Results_baseline,assistiveGoal); % this is my starting point, the Edot at unassisted conditions.

to_plot=0;
if to_plot==1
    [gait_cycle,~]=computeGC(Misc.time,Misc.extra_frames);
    gait_cycle_sel=gait_cycle(1+Misc.extra_frames:end-Misc.extra_frames-1);

    figure(1); clf; 
    plot(gait_cycle_sel,J_normal_TS,'k','LineWidth',2); hold on;
    plot(gait_cycle_sel,J_baseline_TS,':r','LineWidth',2); hold on;
    xlabel('gait cycle[%]'); ylabel([J_baseline_extra.label ' [' J_baseline_extra.unit ']']); 
    title_arrange=['J_avg = ' num2str(J_baseline_avg,'%1.1f') ' ' J_baseline_extra.unit];

    ylim([0 6])
    xlim([0 100])
    title(title_arrange,'Interpreter','none'); 

    figure(2); clf;
    for i=1:40
        subplot(5,8,i)
        hold on
        plot(Results_normal.MActivation(i,:),'k','LineWidth',2);
        plot(Results_baseline.MActivation(i,:),':r','LineWidth',2);
        axis([0 100 0 1])
    end
end
%% Bilevel setup

% update device
clear Device
side_sel=Misc.gait_data.side_sel;

for i=1:length(SDevice)
    SDevice{i}.MuscleGroup(1)={[SDevice{i}.MuscleGroup{1} side_sel]};
end

Device=SDevice;

% get bounds and names
[bounds,bounds_stackUp,nVars_tot]=getControlParamBounds_perType(Device);
[label_deviceAbb]=abbreviation_for_results(Device,bounds);

% number of iterations
% MaxObjectiveEvaluations=50*nVars; % heuristic
MaxObjectiveEvaluations=100; % heuristic

% path names
BilevelOutPath= fullfile(computerPath,Misc.TrialFolder,DirF.select_folder_N3);
BilevelOutName= ['MRSOptimal' SOutName '_' assistiveGoal '_' label_deviceAbb '_iters' num2str(MaxObjectiveEvaluations)];

if strcmp(Scondition, 'runBilevel')
%% Bilevel formulation
% Here we will use Bayesian optimization. 

% - OUTER LEVEL
% Optimizer=Bayesian Optimization, Objective= assistive goal

% - INNER LEVEL
% Optimizer=Direct collocation,    Objective= min(a^2) 
%% Framework bilevel
%% Setup Misc & Device
% run muscle analysis and enable assistive moment
Misc.GetAnalysis = 0;
Misc.Advance.AssistiveDevice  = 1;

% to name and save results
Misc.to_save_results= 0;
Misc.OutPath        = BilevelOutPath;
Misc.OutName        = 'iterations';

Misc.Device=SDevice;
%% Setup Bounds & Constraints
% get control parameter bound
lParamUnit=bounds_stackUp.varUnits;
varNames_tot  =bounds_stackUp.varNames;

% pile ranges as in first rangeX and then rangeY
lb_list_tot=bounds_stackUp.range(:,1); % lower bound
ub_list_tot=bounds_stackUp.range(:,end); % upper bound
nVars  =sum(nVars_tot);

var_containers(nVars) = optimizableVariable('temp', [0, 1]); % Preallocate with dummy variable
for i = 1:nVars
    var_name = varNames_tot{i};
    var_containers(i) = optimizableVariable(var_name, [lb_list_tot(i), ub_list_tot(i)]);
end

% setup constraints
[finalConstraint,initialX]=setup_Constraints_and_InitialGuess(Misc,bounds);

% plot initial guesses
to_plot_IGs=1;

if to_plot_IGs==1
    [gait_cycle,~] = computeGC(Misc.time,Misc.extra_frames);
    nDevs=length(Misc.Device);

    for i=1:nVars
        subplot(2,nVars,i)
        plot(initialX.(varNames_tot{i}),'.k','MarkerSize',20)
        title(varNames_tot{i},'Interpreter','none')
        ylim([lb_list_tot(i) ub_list_tot(i)])
        set(gca,'FontSize',15);
        xlabel('#iterations');
    end

    initialX_val=initialX{:,:};
    moment_spline_all=zeros(nDevs,length(gait_cycle),length(initialX_val));

    iVarPrior=0;
    for iDev=1:nDevs
        iVarOrder=1+iVarPrior:nVars_tot(iDev)+iVarPrior;
        for i=1:length(initialX_val)
            Device{iDev}.Params=initialX_val(i,iVarOrder);
            [assistanceInfo]=generateTorque(Device{iDev},DatStore,Misc.time,Misc.extra_frames);
            moment_spline_all(iDev,:,i)=assistanceInfo.Profile.Torque;
            subplot(2,nDevs,nDevs+iDev); hold on;
            plot(gait_cycle,moment_spline_all(iDev,:,i),'LineWidth',2)
        end
        iVarPrior=nVars_tot(iDev);
    end

    for iDev=1:nDevs
        for i=1:length(initialX_val)
            subplot(2,nDevs,nDevs+iDev); hold on;
            plot(gait_cycle,moment_spline_all(iDev,:,i),'LineWidth',2)
        end
        title(Device{iDev}.MuscleGroup{1}(1:end-2),'Interpreter','none');
        ylabel('torque [Nm]'); xlabel('gait cycle [%]'); xlim([0 100]);
        set(gca,'FontSize',15);
    end
end

%% Setup inner loop
% inner loop definition
enable_param_ver=1;
if enable_param_ver && strcmp(assistiveGoal(1:4),'JRXN'); assistance_goal_upd=[assistiveGoal '_par']; else; assistance_goal_upd=assistiveGoal; end
funInner = @(x) innerLoop(Misc,DatStore,varNames_tot,{x},J_baseline_avg,assistance_goal_upd); 

%% Bayesian optimization setup
%  NOTES: 
% 'OutputFcn'              ,@stopNoImprovementSimple ,...   % to stop earlier
% 'AcquisitionFunctionName','expected-improvement-plus',... % the 'plus' is
% to avoid local minima, yet it disincentives explotaition and keeps the iterations 
% going to low torques, finding suboptimal results -0.5% the best case is
% Such problem comes from adding constraints in the Bayesian OPT
% In the paper 2024, adding dynamic constrains have better performance 
% i.e. find optimal solutions <100 iterations, yet, it is wasteful. Optimizer evaluated 
% (fake/masked) set of parameters, building the surrogate under that assumption
% Here, we propose adding a rich set of optimal guesses where tp and Mp are
% deliverally exploited to maximize likelihood of torque diversity. Also these parameters
% are the once with higher sensibility to the outcome parameter
% USE ardexponential

resultsBayesopt = bayesopt(funInner,var_containers,'IsObjectiveDeterministic',true,...
                           'MaxObjectiveEvaluations',MaxObjectiveEvaluations,"UseParallel",true,'ParallelMethod','clipped',...
                           'AcquisitionFunctionName','expected-improvement',...
                           'XConstraintFcn', finalConstraint,...
                           'PlotFcn', {@plotObjectiveModel,@plotMinObjective},... %,@plotTorque
                           'InitialX', initialX,'ExplorationRatio',0.5); % ExplorationRatio is 0.5 by default

disp(['Bayesian optimization took ' num2str(resultsBayesopt.TotalElapsedTime,'%1.2f') ' secs']);
to_save_results=1;
if to_save_results==1; save(fullfile(Misc.OutPath,['Bayesopt' SOutName '_' assistiveGoal '_' label_deviceAbb '_iters' num2str(MaxObjectiveEvaluations)]),'resultsBayesopt'); end
%% Summary of results
ObjMin     = resultsBayesopt.MinObjective;
VarAtIters = resultsBayesopt.XTrace;
ObjAtIters = resultsBayesopt.ObjectiveTrace;
VarBest    = resultsBayesopt.XAtMinObjective{1,:}';
ObjMinTrace = resultsBayesopt.ObjectiveMinimumTrace;
BestIter    = find(ObjAtIters==ObjMin,1);
nVars       = numel(resultsBayesopt.VariableDescriptions);

range_list=ub_list_tot-lb_list_tot;
bestParams_upd_per=(VarBest-lb_list_tot)./range_list*100;

disp(['Objective, ΔE ' num2str(ObjMin,'%1.2f') '%'])
disp('Optimal control parameters:');
for iVar=1:nVars
    disp(['Param#' num2str(iVar) ':' char(varNames_tot(iVar))  '= '  num2str(VarBest(iVar),'%1.1f') ' ' lParamUnit{iVar} ', [min=' num2str(lb_list_tot(iVar),'%1.1f') ', max=' num2str(ub_list_tot(iVar),'%1.1f') ', limits,' num2str(bestParams_upd_per(iVar),'%1.1f') '%]']);
end
%% Get results from best iteration
iVarPrior=0;
for iDev=1:nDevs
    sVar=nVars_tot(iDev);
    iVarOrder=1+iVarPrior:sVar+iVarPrior;

    Device{iDev}.Params=VarBest(iVarOrder); % peak time, rise time, fall time, peak magnitude
    [assistanceInfo]=generateTorque(Device{iDev},DatStore,Misc.time,Misc.extra_frames);
    Device{iDev}.Assistance=assistanceInfo;
    iVarPrior=sVar;
end
Misc.Device = Device;

% to name and save results
Misc.to_save_results= 1;
Misc.OutPath        = BilevelOutPath;
Misc.OutName        = BilevelOutName;
[Results_bilevel,~,Misc]= MRS_Formulate_and_Solve_NeuroCons(Misc,DatStore);
[J_bilevel_avg,J_bilevel_TS,J_bilevel_extra] = computeOuterLoopFunction(Misc,Results_bilevel,assistiveGoal);

% output of the function
MRS{3}.Results=Results_bilevel;  MRS{3}.DatStore=DatStore;    MRS{3}.Misc=Misc; %synergy + bilevel
end

if strcmp(Scondition, 'getBilevel')
    disp(['Condition:' Scondition '. Stop after retrieving bilevel optimal results.']);
    MRS_bilevel      = load(fullfile(BilevelOutPath, [BilevelOutName 'Results.mat']));
    Results_bilevel  = MRS_bilevel.Results;
    DatStore         = MRS_bilevel.DatStore;
    Misc             = MRS_bilevel.Misc;

    MRS{3}.Results=Results_bilevel;  MRS{3}.DatStore=DatStore;    MRS{3}.Misc=Misc; %synergy + bilevel
    return;  % Exit function early
end

%% Plot rapid figure
to_plot_fast=0;
if to_plot_fast==1
    figure(10); clf;
    for i=1:40
        subplot(5,8,i)
        plot( Results_baseline.MActivation(i,:),'k','LineWidth',1); hold on
        plot( Results_bilevel.MActivation(i,:),'r','LineWidth',2)
        ylim([0 1])
        title(Results_bilevel.MuscleNames{i})
    end
end
%% Plot best iteration

to_plot_summary=1;
if to_plot_summary==1
    figure(11); clf; set(gcf,'color','w','Visible','on','WindowState', 'maximized');
    
    nDevs = length(Misc.Device);
    IDinterp=DatStore.IDinterp;
    
    [gait_cycle,time_series] = computeGC(Misc.time,Misc.extra_frames);
    to_gc  = gait_cycle(find(time_series>=Misc.gait_data.toeOff_time,1)); % toe off event
    gait_cycle_ranged=gait_cycle(1+Misc.extra_frames:end-Misc.extra_frames-1);

    all_torque_profiles = zeros(nDevs,size(Results_bilevel.MActivation,2),length(resultsBayesopt.UserDataTrace));
    for i = 1:length(resultsBayesopt.UserDataTrace)
        if ~isempty(resultsBayesopt.UserDataTrace{i})
            all_torque_profiles(:,:,i) = resultsBayesopt.UserDataTrace{i}.Torque;
        end
    end
    
    nRowsTot=nDevs+1;

    % PLOT: GENERAL -------------------------------------------------------
    val_min=min(ObjAtIters); val_max=max(ObjAtIters); range=val_max-val_min;
    nIGs=size(initialX{:,:},1);

    % ALL ITERATIONS
    subplot(nRowsTot,3,1)
    plot(ObjAtIters,'.k','MarkerSize',20); hold on;
    xline(nIGs,'--g','LineWidth',2)
    plot(ObjMinTrace,':b','LineWidth',2,'LineStyle',':')
    plot(BestIter,ObjMin,'.r','MarkerSize',30)
    set(gca,'FontSize',13);
    xlabel('#iterations'); ylabel('obj. fun [%]');
    xlim([0 resultsBayesopt.NumObjectiveEvaluations]); ylim([val_min-0.1*range val_max+0.1*range]);
    title(['Convergence, best at ' num2str(BestIter) ' interation'])
    
    % METRIC: TIME SERIES    
    subplot(nRowsTot,3,2);
    plot(gait_cycle_ranged,J_baseline_TS,'k','LineWidth',2); hold on
    plot(gait_cycle_ranged,J_bilevel_TS,'r','LineWidth',2); hold on
    set(gca,'FontSize',13);
    xlabel('gait cycle [%]')
    ylabel([J_bilevel_extra.label ' [' J_bilevel_extra.unit ']'])
    title('Targeted metric: time series')
    % J_real=(mean(J_bilevel_TS)-mean(J_baseline_TS))/mean(J_baseline_TS)*100; % CHECK YOUR RESULTS (2!)
    
    % METRIC: MEAN VALUE
    subplot(nRowsTot,3,3); hold on;
    x_label={'Normal' 'Bilevel'};
    X = categorical(x_label);
    X = reordercats(X,x_label);

    data_val=[J_baseline_avg J_bilevel_avg];

    color_bar={'k' 'r'};

    for j=1:length(x_label)
        data= data_val(j);
        
        bar(X(j),data,'FaceColor',color_bar{j})

        if j>1 % compute relative change compared to the first one
            E_normal=data_val(1);
            change_per=(data-E_normal)/E_normal*100;
            text(j,E_normal,[num2str(change_per,'%+1.1f') ' %'],'HorizontalAlignment','center','FontSize',15)
        end
    end
    ylim([0 4])
    set(gca,'FontSize',13);
    xlabel('conditions'); ylabel([J_bilevel_extra.label ' [' J_bilevel_extra.unit ']']) %(both legs) 
    title('Targeted metric: mean value')

    % TORQUE PROFILES
    for iDev=1:nDevs
        devJoint=Device{iDev}.MuscleGroup{1}(1:end);
        ID_sel=IDinterp(:,strcmp(DatStore.DOFNames,devJoint))*Device{iDev}.MuscleGroup{2};

        devType=Device{iDev}.Type{1};

        subplot(nRowsTot,3,4+(iDev-1)*3);
        % iterations
        for nIter=1:resultsBayesopt.NumObjectiveEvaluations
            plot(Misc.Device{iDev}.Assistance.Profile.GaitCycle,all_torque_profiles(iDev,:,nIter),'linewidth',0.5,'color','#A9A9A9'); hold on
        end

        % inverse dynamics
        plot(Misc.Device{iDev}.Assistance.Profile.GaitCycle,ID_sel,'k','LineWidth',3);
        
        % optimal from MRS
        plot(Misc.Device{iDev}.Assistance.Profile.GaitCycle,Misc.Device{iDev}.Assistance.Profile.Torque,'r','LineWidth',4);
        
        % optimal from bilevel bestIter
        plot(Misc.Device{iDev}.Assistance.Profile.GaitCycle,all_torque_profiles(iDev,:,BestIter),':c','LineWidth',4);

        xline(to_gc,'--m','LineWidth',2)
        % extras
        if strcmp(Misc.Device{iDev}.Type{1},'active')
            plot(Misc.Device{iDev}.Assistance.Additional.Nodes_x,Misc.Device{iDev}.Assistance.Additional.Nodes_y,'.g','MarkerSize',30)
        end

        set(gca,'FontSize',13);
        xlim([0 100])
        xlabel('gait cycle [%]')
        ylabel('torque [Nm]')
        title(['Joint: ' devJoint(1:end-2) '   Type: ' devType],'Interpreter','none')
        ylim([min(ID_sel) max(ID_sel)])
    end

    counter=0;
    % PARAMETERS
    for iDev=1:nDevs
        sVar=nVars_tot(iDev);
        for i=1:sVar
            counter=counter+1;
            ParamName=Misc.Device{iDev}.Assistance.Param.Names{i};
            ParamUnits=Misc.Device{iDev}.Assistance.Param.Units{i};

            subplot(nRowsTot,6,8+(iDev-1)*6+i); hold on
            range=ub_list_tot(counter)-lb_list_tot(counter);
            ylim([lb_list_tot(counter)-0.1*range ub_list_tot(counter)+0.1*range])
            yline(lb_list_tot(counter),':b','LineWidth',2);
            yline(ub_list_tot(counter),':b','LineWidth',2)
            xline(nIGs,'--g','LineWidth',2)
            plot(VarAtIters.(ParamName),'.k','MarkerSize',15)
            xlabel('#iterations'); 
            title(['Optimal ' ParamName '= ' num2str(resultsBayesopt.XAtMinObjective.(ParamName),'%1.1f')],'Interpreter','none');
            ylabel([ParamName ' [' ParamUnits ']' ],'Interpreter','none')

            plot(BestIter,resultsBayesopt.XAtMinObjective.(ParamName),'.r','MarkerSize',20);
            set(gca,'FontSize',13);
        end
    end

    sgtitle(['SUB:' info.SubjectSelection '   MOT:' info.MotionSelection '   DEVICES:' label_deviceAbb],'interpreter','none') 

    to_save=1;
    if to_save==1
        saveas(gcf,fullfile(Misc.OutPath,['Summary' SOutName '_' assistiveGoal '_' label_deviceAbb '_iters' num2str(MaxObjectiveEvaluations)]),'png');
    end
end



end

% Usage:
% analyzeOptimizerResolution(resultsBayesopt, varNames);



function [objective, coupledconstraints, userdata] = innerLoop(Misc,DatStore,varNames,exoParam,J_normalized,devAim)

nDevs=length(Misc.Device);
iVarPrior=0;
for iDev=1:nDevs
    type=Misc.Device{iDev}.Type{1};
    if strcmp(type,'quasi-passive')
        nVars=2; 
    elseif strcmp(type,'active')
        nVars=4; 
    end

    iVarOrder=1+iVarPrior:nVars+iVarPrior;
    varNames_dev=varNames(iVarOrder);
    Params=zeros(1,nVars);
    for i=1:nVars
        Params(i)=exoParam{1}.(varNames_dev{i});
    end

    Misc.Device{iDev}.Params      = Params;
    [assistanceInfo]=generateTorque(Misc.Device{iDev},DatStore,Misc.time,Misc.extra_frames);

    Misc.Device{iDev}.Assistance = assistanceInfo;

    iVarPrior=nVars;
end


% Run inner loop
try
    [Results,~,Misc] = MRS_Formulate_and_Solve_NeuroCons(Misc,DatStore);
    [Javg,~,~]=computeOuterLoopFunction(Misc,Results,devAim);
    objective=(Javg-J_normalized)/J_normalized*100; % as a percentage of the unassisted condition
catch ME
    objective=0; % assigned value if simulation fails
    warning('error! revise')
end

% resultsBayesopt.UserDataTrace{i}.moment_spline;

coupledconstraints=[];
for iDev=1:nDevs
    Torque=Misc.Device{iDev}.Assistance.Profile.Torque;
    userdata.Torque(iDev,:)=Torque;
end
end


% STRATEGIES
% Stop when, after 100 iterations, there as no improvement
function stop = stopNoImprovementSimple(results, state)
    persistent bestValue noImproveCount
    stop = false;
    
    if strcmp(state, 'iteration') && ~isempty(results.ObjectiveTrace)
        currentBest = min(results.ObjectiveTrace);
        
        if isempty(bestValue) || currentBest < bestValue
            % Improvement found
            bestValue = currentBest;
            noImproveCount = 0;
            fprintf('Improvement! New best: %.6f\n', currentBest);
        else
            % No improvement
            noImproveCount = noImproveCount + 1;
            fprintf('No improvement for %d iterations (best: %.6f)\n', noImproveCount, bestValue);
        end
        
        if noImproveCount >= 100
            stop = true;
            fprintf('Stopped: No improvement for 100 iterations\n');
        end
    end
end

% -> WHEREN THE LAST "5" BEST DISTINCT SOLUTIONS ARE LESS THAN 0.01
% PROBLEM: What if the best is found on the first one? Then there wont be a
% better distinct solutions. Hence, it will seems like we punish an early
% good solution
% Define the stopping condition function
function stop = stopOnSmallImprovementsSimple(results, state)
    persistent bestHistory improvementHistory
    stop = false;
    
    switch state
        case 'initial'
            bestHistory = [];
            improvementHistory = [];
            fprintf('Stopping when last 5 distinct best values improve by < 0.01\n');
            
        case 'iteration'
            if ~isempty(results.ObjectiveTrace)
                currentBest = min(results.ObjectiveTrace);
                
                % Only track when we find a new best value
                if isempty(bestHistory) || currentBest < bestHistory(end)
                    
                    % Calculate improvement from previous best
                    if ~isempty(bestHistory)
                        improvement = bestHistory(end) - currentBest;
                        improvementHistory = [improvementHistory, improvement];
                    else
                        % First value - no improvement to calculate
                        improvementHistory = [];
                    end
                    
                    % Store the new best
                    bestHistory = [bestHistory, currentBest];
                    
                    % Check stopping condition (need at least 5 improvements)
                    if length(improvementHistory) >= 5
                        lastFiveImprovements = improvementHistory(end-4:end);
                        
                        if all(lastFiveImprovements < 0.01)
                            stop = true;
                            fprintf('\n*** STOPPING: Last 5 distinct best values improved by < 0.01 ***\n');
                            fprintf('Improvements: ');
                            fprintf('%.4f ', lastFiveImprovements);
                            fprintf('\nFinal best: %.6f\n', currentBest);
                        else
                            fprintf('Iter %d: Best=%.6f, Improvement=%.6f\n', ...
                                length(results.ObjectiveTrace), currentBest, improvement);
                        end
                    else
                        fprintf('Iter %d: Best=%.6f', length(results.ObjectiveTrace), currentBest);
                        if ~isempty(improvementHistory)
                            fprintf(', Improvement=%.6f', improvementHistory(end));
                        end
                        fprintf(' (%d/5 improvements tracked)\n', length(improvementHistory));
                    end
                else
                    fprintf('Iter %d: No new best (current best: %.6f)\n', ...
                        length(results.ObjectiveTrace), bestHistory(end));
                end
            end
    end
end

function stop = plotTorque(results, state)
persistent hs list_torques best_torques counter
stop = false;
    switch state
        case 'initial'
            hs = figure;  
            counter = 0;
            list_torques=zeros(116,results.Options.MaxObjectiveEvaluations);
            best_torques=zeros(116,results.Options.MaxObjectiveEvaluations);

            figure(hs)
            subplot(2,2,1); hold on
            xlabel 'Gait cycle [%]'; ylabel 'Torque [Nm]'; title 'Candidate torques'

            subplot(2,2,2); hold on
            xlabel 'Gait cycle [%]'; ylabel 'Torque [Nm]'; title 'Best torques'
            
            for i=1:4
                subplot(2,4,4+i); hold on
                title(results.VariableDescriptions(i).Name)
                ylim(results.VariableDescriptions(i).Range)
                xlabel 'Iterations'; ylabel 'Variable';
            end

        case 'iteration'
            
            XTrace=results.XTrace;
           
            % current
            assistanceType='spline#N3';
            assistanceParms=[XTrace{end,1} XTrace{end,2} XTrace{end,3} XTrace{end,4}];
            [assistanceInfo]=generateTorque(assistanceType,assistanceParms,[1.2000 2.3500],5);
            moment_spline=assistanceInfo.Profile.Torque;

            list_torques(:,results.NumObjectiveEvaluations)=moment_spline;

            figure(hs)
            subplot(2,2,1)
            plot(moment_spline,'k','LineWidth',1);


            % best
            if results.ObjectiveTrace(end)<=results.MinObjective
                counter=counter+1;
                best_torques(:,counter)=moment_spline;

                subplot(2,2,2)
                plot(moment_spline,'r','LineWidth',1); 
            end

            for i=1:4
            subplot(2,4,4+i)
            plot(XTrace{:,i},'.k','MarkerSize',15)    
            end
    end
end
% function mergedConstraint = mergeConstraints(constraintFcn_tot)
%     % Remove empty cells
%     validConstraints = constraintFcn_tot(~cellfun(@isempty, constraintFcn_tot));
% 
%     if isempty(validConstraints)
%         mergedConstraint = @(x) true; % Always true if no constraints
%         return;
%     end
% 
%     % Start with first constraint
%     mergedConstraint = validConstraints{1};
% 
%     % Combine with remaining constraints using AND
%     for i = 2:length(validConstraints)
%         currentFcn = validConstraints{i};
%         mergedConstraint = @(x) mergedConstraint(x) & currentFcn(x);
%     end
% 
% end
function mergedConstraint = mergeConstraints(constraintFcn_tot)
    % Remove empty cells
    validConstraints = constraintFcn_tot(~cellfun(@isempty, constraintFcn_tot));
    
    if isempty(validConstraints)
        mergedConstraint = [];
        return;
    end
    
    % If there's only one constraint and it's @(x)true, return empty
    if length(validConstraints) == 1
        try
            testX = struct('AP_tp', {1}, 'AP_tf', {1}, 'AP_tr', {1});
            if validConstraints{1}(testX)
                % Quick test - if it returns true for a simple case
                mergedConstraint = [];
                return;
            end
        catch
            % If testing fails, keep the constraint
        end
    end
    
    % Start with first constraint
    mergedConstraint = validConstraints{1};
    
    % Combine with remaining constraints using AND
    for i = 2:length(validConstraints)
        currentFcn = validConstraints{i};
        mergedConstraint = @(x) mergedConstraint(x) & currentFcn(x);
    end
end

function [finalConstraint,initialX]=setup_Constraints_and_InitialGuess(Misc,bounds)

[gait_cycle,time_series] = computeGC(Misc.time,Misc.extra_frames);
to_gc  = gait_cycle(find(time_series>=Misc.gait_data.toeOff_time,1)); % toe off event
min_gc = 15;                                                          % earliest possible applied torque

Device=Misc.Device;
nDevs=length(Device);
IG_method_list    = zeros(1,nDevs); % to list the initial guess seeding method
constraintFcn_tot = cell(1,nDevs);  % to pile all constraints
for iDev=1:nDevs
    DOFName=Device{iDev}.MuscleGroup{1}(1:end-2);
    if strcmp(Device{iDev}.Type{1},'active')
        switch DOFName
            case 'ankle_angle'
                constraintFcn = @(x)(x.AP_tp + x.AP_tf <= to_gc) & (x.AP_tp - x.AP_tr >= min_gc);
            case 'knee_angle'
                constraintFcn = @(x)(x.KE_tp - x.KE_tr >= 0);
            case 'hip_flexion'
                constraintFcn = @(x)(x.HF_tp + x.HF_tf <= 100) & (x.HF_tp - x.HF_tr >= min_gc) ;
            case 'hip_adduction'
                constraintFcn = @(x)(x.HB_tp + x.HB_tf <= 100) & (x.HB_tp - x.HB_tr >= min_gc) ;
            otherwise
                error('not defined')
        end
        SelectedMethod=5;

    elseif strcmp(Device{iDev}.Type{1},'quasi-passive')
        constraintFcn =[];
        SelectedMethod=0;
    end
    constraintFcn_tot{iDev} = constraintFcn;
    IG_method_list(iDev) = SelectedMethod;
end
finalConstraint = mergeConstraints(constraintFcn_tot);

% set # of initial guesses
list_numInitialPoints=zeros(1,nDevs);
for iDev=1:nDevs
    if strcmp(Device{iDev}.Type{1},'active')
        numInitialPoints_def=3*9; % 9 per each of the 3 regions = 27
    end
    if strcmp(Device{iDev}.Type{1},'quasi-passive')
        numInitialPoints_def=3*9; % same as before
    end
    list_numInitialPoints(iDev)=numInitialPoints_def;
end

% compute initial guesses
initialX_list=cell(nDevs,1);
for iDev=1:nDevs
    SelectedMethod  =IG_method_list(iDev);
    numInitialPoints=list_numInitialPoints(iDev);
    constraintFcn   =constraintFcn_tot{iDev};
    infoDevBounds   =bounds{iDev};

    lb_list =infoDevBounds.range(:,1);
    ub_list =infoDevBounds.range(:,end);
    varNames=infoDevBounds.varNames;

    initialX_dev = initialGuessStrategy(SelectedMethod,numInitialPoints,to_gc, min_gc,lb_list,ub_list,varNames,constraintFcn);
    initialX_list(iDev)={initialX_dev};
end
initialX = horzcat(initialX_list{:});
end

% setup abbreviation for labeling results
function [devAbb_con]=abbreviation_for_results(Device,bounds)
nDevs=length(Device);
devAbb=cell(nDevs,1);
for iDev=1:nDevs
    if strcmp(Device{iDev}.Type{1},'active')
        typeAbb='A';
    elseif strcmp(Device{iDev}.Type{1},'quasi-passive')
        typeAbb='Q';
    else
        error('not defined the type of device assistance')
    end
    devAbb(iDev)={[bounds{iDev}.varNames{1}(1:2) '(' typeAbb ')']};
end
devAbb_con=strjoin(devAbb, '&');
end

function [Misc]=setup_for_JRXN(Misc)
Misc.MotionSelection   ='oneSingleMotion';
Misc.ForceLabel        ='SO_label_rajagopal';
Misc.extra_folder_name ='temp';
Misc.extra_file_name   ='unassisted';
end

function [W,H,MActivation_reconstructed_array,synMetrics]=synergyAnalysisSequential(Results,synergy_list,plot_flag)
%% Compute synergies with SEQUENTIAL NNMF for nested subspaces
% This ensures W_n is a subset of W_{n+1} for stability

rng(100); % for reproducibility

MActivation = Results.MActivation;
[nMuscles, nTimePoints] = size(MActivation);

% NNMF parameters
n_replicates   = 100;      % Number of random starts
algorithm_type = 'als';    % 'als' or 'mult'
max_iterations = 2000;     % Maximum iterations
tolerance      = 1e-6;     % Convergence tolerance

fprintf('Running SEQUENTIAL NNMF with %d replicates, algorithm: %s\n', n_replicates, algorithm_type);

% Sort synergy_list to ensure we go from smallest to largest
synergy_list = sort(synergy_list);
nSynergies = length(synergy_list);

% Initialize storage
reconstruction_error = zeros(1, nSynergies);
W = cell(1, nSynergies);
H = cell(1, nSynergies);
W_sequential = cell(1, nSynergies); % Store sequential W matrices

% Set NNMF options
options = statset('MaxIter', max_iterations, 'Display', 'final', 'TolFun', tolerance);

%% SEQUENTIAL EXTRACTION - Key modification!
for iSyn = 1:nSynergies
    current_syn = synergy_list(iSyn);
    fprintf('  Extracting %d synergies (sequential)... ', current_syn);
    
    if iSyn == 1
        % For the first (smallest) synergy count, run normal NNMF
        [W{1,iSyn}, H{1,iSyn}, reconstruction_error(1,iSyn)] = nnmf(...
            MActivation, ...
            current_syn, ...
            'replicates', n_replicates, ...
            'algorithm', algorithm_type, ...
            'options', options);
        
        % Store for sequential initialization
        W_sequential{1,iSyn} = W{1,iSyn};
        
    else
        % For subsequent synergy counts, use previous W as initialization
        
        % Get previous W matrix
        prev_W = W_sequential{1,iSyn-1};
        prev_n_syn = synergy_list(iSyn-1);
        
        % Calculate how many new synergies to add
        n_new_synergies = current_syn - prev_n_syn;
        
        % Initialize with previous W plus random columns for new synergies
        % Add small random perturbation to previous columns to help convergence
        W_init = [prev_W .* (0.95 + 0.1*rand(size(prev_W))), ...
                  rand(nMuscles, n_new_synergies) * 0.5];
        
        % Run NNMF with warm start
        [W{1,iSyn}, H{1,iSyn}, reconstruction_error(1,iSyn)] = nnmf(...
            MActivation, ...
            current_syn, ...
            'w0', W_init, ...          % Initial W matrix
            'h0', rand(current_syn, nTimePoints), ... % Initial H
            'replicates', 1, ...       % Only need 1 replicate with good init
            'algorithm', algorithm_type, ...
            'options', options);
        
        % Store for next iteration
        W_sequential{1,iSyn} = W{1,iSyn};
        
        % Optional: Check if previous columns are preserved
        % (they should be, but let's verify)
        prev_in_current = pinv(W{1,iSyn}) * prev_W;
        reconstruction_check = norm(prev_W - W{1,iSyn} * prev_in_current, 'fro');
        fprintf('Subspace preservation error: %.4f\n', reconstruction_check);
    end
    
    fprintf('Reconstruction error: %.4f\n', reconstruction_error(1,iSyn));
end

%% Alternative: PCA-based approach (guaranteed nesting)
% Uncomment if you want to try PCA instead
% [coeff, ~, ~] = pca(MActivation');
% for iSyn = 1:nSynergies
%     current_syn = synergy_list(iSyn);
%     W{1,iSyn} = coeff(:,1:current_syn);
%     % For PCA, H = W' * MActivation (but not non-negative)
%     % You might need to run NNMF with PCA initialization
% end

%% Reconstructed synergies and metrics (unchanged from your code)
full_data_length = nTimePoints;
MActivation_reconstructed_array = zeros(nSynergies, nMuscles, full_data_length);
VAF = zeros(nSynergies,1);
RMSE = zeros(nSynergies,1);
RMSE_per_muscle = zeros(nSynergies,nMuscles);
R_per_muscle = zeros(nSynergies,nMuscles);

% Calculate Variance Accounted For (VAF)
total_variance = sum(var(MActivation, 0, 2));

for iSyn = 1:nSynergies
    MActivation_reconstructed = W{1,iSyn} * H{1,iSyn};
    MActivation_reconstructed_array(iSyn, :, :) = MActivation_reconstructed;

    % Metrics
    VAF(iSyn,1) = 1 - (sum(var(MActivation - MActivation_reconstructed, 0, 2)) / total_variance);
    RMSE(iSyn,1) = sqrt(mean((MActivation - MActivation_reconstructed).^2, 'all'));
    RMSE_per_muscle(iSyn,:) = sqrt(mean((MActivation - MActivation_reconstructed).^2, 2));
    for m = 1:nMuscles
        R_per_muscle(iSyn,m) = corr(MActivation(m,:)', MActivation_reconstructed(m,:)');
    end    
end

%% Additional check: Verify subspace nesting
fprintf('\n=== Subspace Nesting Verification ===\n');
for iSyn = 2:nSynergies
    prev_syn = synergy_list(iSyn-1);
    curr_syn = synergy_list(iSyn);
    
    W_prev = W{1,iSyn-1};
    W_curr = W{1,iSyn};
    
    % Project previous W onto current subspace
    W_projected = W_curr * (pinv(W_curr) * W_prev);
    
    % Calculate reconstruction error
    nesting_error = norm(W_prev - W_projected, 'fro') / norm(W_prev, 'fro');
    
    fprintf('W_%d in W_%d subspace: Relative error = %.4f%%', ...
        prev_syn, curr_syn, nesting_error*100);
    
    if nesting_error < 0.05  % 5% threshold
        fprintf(' ✓ (Well nested)\n');
    else
        fprintf(' ✗ (Poor nesting)\n');
    end
end

%% Store metrics (unchanged)
synMetrics.N = synergy_list;
synMetrics.VAF = VAF;
synMetrics.RMSE = RMSE;
synMetrics.RMSE_per_muscle = RMSE_per_muscle;
synMetrics.R_per_muscle = R_per_muscle;
synMetrics.W_sequential = W_sequential; % Store sequential W matrices

% Display results summary
fprintf('\n=== SEQUENTIAL NNMF Results Summary ===\n');
fprintf('Synergies  Recon.Error  VAF\n');
fprintf('------------------------------\n');
for iSyn = 1:nSynergies
    fprintf('    %d        %.4f     %.3f\n', ...
        synergy_list(iSyn), reconstruction_error(1,iSyn), VAF(iSyn,1));
end

%% Plotting section (unchanged - can keep your original plotting code)
% ... [Your plotting code here - same as before] ...

end