addpath(genpath('C:\Users\movea\Documents\GitHub\GoalOptimized'));
% Target objective
assistance_goal='eDot_MCLU24';
close all;
%% baseline
Misc.GetAnalysis = 0;
Misc.Advance.AssistiveDevice  = 0;
Misc.Advance.SynCons          = 0;
Misc.Advance.ApplySynergyConstraints=0;

DatStore=DatStore_normal;
ExecutionTime_1=datetime('now');
[Results_noSyn,~,Misc] = MRS_Formulate_and_Solve_NeuroCons(Misc,DatStore);
ExecutionTime_2=datetime('now'); duration_loop=seconds(duration(ExecutionTime_2-ExecutionTime_1));
disp(['simulation took ' num2str(duration_loop,'%1.2f') ' secs']);

[J_baseline_avg,J_baseline_TS,J_baseline_extra] = computeOuterLoopFunction(Misc,Results_noSyn,assistance_goal); % this is my starting point, the Edot at unassisted conditions.
%%
[W,H,~,synergy_list]=synergyAnalysis(Results_noSyn,[0 0 1]);

synSelected=5; % 3 4 5 6
synInd=find(synergy_list==synSelected);
Wsel=W{synInd};
Hsel=H{synInd};
SynergyActivationComputed=Wsel*Hsel;
%%
Misc.Advance.AssistiveDevice  = 0;
Misc.Advance.ApplySynergyConstraints=1;
Misc.SynCons.W=Wsel;
Misc.SynCons.H=Hsel;
Misc.SynCons.N=synSelected;

ExecutionTime_1=datetime('now');
[Results_SynCons,~,Misc]= MRS_Formulate_and_Solve_NeuroCons(Misc,DatStore);
ExecutionTime_2=datetime('now'); duration_loop=seconds(duration(ExecutionTime_2-ExecutionTime_1));
disp(['simulation took ' num2str(duration_loop,'%1.2f') ' secs']);

[J_baselineUpd_avg,J_baselineUpd_TS,J_baselineUpd_extra] = computeOuterLoopFunction(Misc,Results_SynCons,assistance_goal); % this is my starting point, the Edot at unassisted conditions.

%%
figure; clf;
for iMus=1:40
    subplot(5,8,iMus); hold on
    % plot(Results_normalSynCons.SynergyExcitation(iMus,:),'k')
    plot(Results_noSyn.MActivation.genericMRS(iMus,:),'k','LineWidth',3)    
    % plot(SynergyActivationComputed(iMus,:),'-b','LineWidth',2)

    plot(Results_SynCons.MActivation.genericMRS(iMus,:),':g','LineWidth',2)
    % plot(Results_SynCons.SynergyActivation(iMus,:),':r','LineWidth',2)
    
    title(Results_noSyn.MuscleNames{iMus})
    ylim([0 1])
end

%% exoskeleton

Misc.GetAnalysis = 0;
Misc.Advance.AssistiveDevice  = 1;

MP_list=[15:15:100]; % 100
% MP_list=[7.5:7.5:60];
% MP_list=[10:10:80];
% MP_list=[1:1:8];

nMP=length(MP_list);

data_length=length(Results.MActivation.genericMRS(1,1+Misc.extra_frames:end-Misc.extra_frames-1));

Results_MP=cell(nMP,1);
J_mus_list=cell(nMP,1);
J_avg_list=zeros(nMP,1);
J_TS_list =zeros(nMP,data_length);

ExecutionTime_1=datetime('now');
for iMP=1:nMP
    sMP=MP_list(iMP);
    clear Device
    Device{1}.Mode        = 'prescribed'; % opts: optimized and prescribed
    Device{1}.MuscleGroup = {['ankle_angle_' Misc.gait_data.side_sel] -1};
    Device{1}.Type        = {'active' 'spline#N3'}; % opts: active, quasi-passive, passive, EMG-driven
    Device{1}.Params      = [54.1 22 8.6 sMP]; %55 20 10

    % Device{1}.Mode        = 'prescribed'; % opts: optimized and prescribed
    % Device{1}.MuscleGroup = {['hip_flexion_' Misc.gait_data.side_sel] 1};
    % Device{1}.Type        = {'active' 'spline#N3'}; % opts: active, quasi-passive, passive, EMG-driven
    % Device{1}.Params      = [62 30 35 sMP];
    
    % Device{1}.Mode        = 'prescribed'; % opts: optimized and prescribed
    % Device{1}.MuscleGroup = {['hip_adduction_' Misc.gait_data.side_sel] -1};
    % Device{1}.Type        = {'active' 'spline#N3'}; % opts: active, quasi-passive, passive, EMG-driven
    % Device{1}.Params      = [50 15 10 sMP];

    % Device{1}.Mode        = 'prescribed'; % opts: optimized and prescribed
    % Device{1}.MuscleGroup = {['knee_angle_' Misc.gait_data.side_sel] -1};
    % Device{1}.Type        = {'active' 'spline#N3'}; % opts: active, quasi-passive, passive, EMG-driven
    % Device{1}.Params      = [15 12 17 sMP];

    % Device{1}.Mode        = 'prescribed'; % opts: optimized and prescribed
    % Device{1}.MuscleGroup = {['ankle_angle_' Misc.gait_data.side_sel] -1};
    % Device{1}.Type        = {'quasi-passive' 'clutchSpring'}; % opts: active, quasi-passive, passive, EMG-driven
    % Device{1}.Params      = [20 sMP];

    [assistanceInfo]=generateTorque(Device{1},DatStore_normal,Misc.time,Misc.extra_frames);

    Device{1}.Assistance = assistanceInfo;
    Misc.Device=Device;

    % to name and save results
    Misc.to_save_results= 0;
    Misc.OutName= 'example';

    % time performance of a single loop
    ExecutionTime_1=datetime('now');
    % [Results_assisted,~,Misc] = MRS_Formulate_and_Solve_NeuroCons(Misc,DatStore_normal);
    [Results_assisted,~,Misc]= MRS_Formulate_and_Solve_NeuroCons(Misc,DatStore_normal);
    ExecutionTime_2=datetime('now'); duration_loop=seconds(duration(ExecutionTime_2-ExecutionTime_1));
    disp(['simulation took ' num2str(duration_loop,'%1.2f') ' secs']);

    % Check results
    to_plot=0;
    if to_plot==1
        figure(1); clf; % muscle activations
        for i=1:40
            subplot(5,8,i)
            hold on
            plot(Results_normal.MActivation.genericMRS(i,:),'k')
            plot(Results_assisted.MActivation.genericMRS(i,:),'r')
            ylim([0 1])
            title(Results_assisted.MuscleNames{i})

        end
        sgtitle('muscle activations')

        figure(2); clf; % device torque
        nDevs=length(Device);
        norMass=Misc.subject_data.subject_mass;
        for iDev=1:nDevs
            gaitCycle=Device{iDev}.Assistance.Profile.GaitCycle;
            iDOF=strcmp(DatStore_normal.DOFNames,Device{iDev}.MuscleGroup{1});
            sID =DatStore_normal.IDinterp(:,iDOF)*Device{iDev}.MuscleGroup{2};
            Torque=Device{iDev}.Assistance.Profile.Torque;
            subplot(nDevs,1,iDev)
            plot(gaitCycle,sID/norMass,'k'); hold on
            plot(gaitCycle,Torque/norMass,'r')
            title([Device{iDev}.MuscleGroup{1}],'Interpreter','none')
        end
        sgtitle('assistive torque')

        figure(3); clf;
        plot(J_baseline_TS,'k','LineWidth',2); hold on
        plot(J_example_TS,'r','LineWidth',2)
        xlabel('gait cycle[%]'); ylabel([J_example_extra.label ' [' J_example_extra.unit ']']);
        title([ 'baseline (avg value) = ' num2str(J_baseline_avg,'%1.1f') ' ' J_example_extra.unit ...
            ' assisted (avg value) = ' num2str(J_example_avg,'%1.1f')  ' ' J_example_extra.unit ...
            ' - change: ' num2str(J_value,'%+1.1f') '%']);
    end

    Results_MP(iMP)={Results_assisted};

end

ExecutionTime_2=datetime('now'); duration_loop2=seconds(duration(ExecutionTime_2-ExecutionTime_1));
disp(['simulation took ' num2str(duration_loop2,'%1.2f') ' secs']);
%% compute objective
for iMP=1:nMP
    Results_assisted=Results_MP{iMP};
    [J_example_avg,J_example_TS,J_example_extra] = computeOuterLoopFunction(Misc,Results_assisted,assistance_goal);
    J_value=(J_example_avg-J_baselineUpd_avg)/J_baselineUpd_avg*100;
    disp(['goal change: ' num2str(J_value,'%+1.2f') '%'])
    J_avg_list(iMP)=J_value;
    J_TS_list(iMP,:) =J_example_TS;
    J_mus_list(iMP)={J_example_extra};
end    
%% find optimal torque
% color_list = {'k' '#8B0000' '#B22222' '#DC143C' '#FF4500' '#FF6B6B' '#FFB6C1'};
% color_list = {'k' '#5A0000' '#8B0000' '#B22222' '#DC143C' '#FF4500' '#FFB6C1' '#FFE4E1'};
color_list = {'k' '#330000' '#660000' '#990000' '#CC0000' '#FF0000' '#FF4D4D' '#FF8080' '#FFB3B3' '#FFE6E6'};

figure(50); clf; set(gcf,'WindowState', 'maximized','Color','white');
MP_list_nor   =[0 MP_list/Misc.subject_data.subject_mass];
J_avg_list_com=[0 J_avg_list'];

legend_entries = cell(1, length(MP_list_nor));
legend_entries{1} = 'no exo'; % First entry is for zero torque

for iMP=1:length(MP_list_nor)
    plot(MP_list_nor(iMP), J_avg_list_com(iMP),'marker','.', 'MarkerSize', 60,'Color',color_list{iMP},'LineStyle','none');
    hold on;

    if iMP > 1
        legend_entries{iMP} = [num2str(MP_list(iMP-1)) ' Nm'];
    end
end
% Cubic interpolation
p = polyfit(MP_list_nor, J_avg_list_com, 2); % Fit 3rd degree polynomial

% Create smooth curve for plotting
x_fit = linspace(min(MP_list_nor), max(MP_list_nor), 100);
y_fit = polyval(p, x_fit);

% Plot the cubic fit
% plot(x_fit, y_fit, 'b-', 'LineWidth', 1.5);

% Find optimal point - simple method
[optimal_J, idx] = min(y_fit);
optimal_MP = x_fit(idx);

% Plot the optimal point
% plot(optimal_MP, optimal_J, 'go', 'MarkerSize', 30, 'MarkerFaceColor', 'g');

% real optimal
[optimal_J_real, idx_real] = min(J_avg_list_com);
optimal_MP_real = J_avg_list_com(idx_real);

plot(MP_list_nor(idx_real),J_avg_list_com(idx_real),'ok','MarkerSize',40);   
text(MP_list_nor(idx_real),J_avg_list_com(idx_real)+2,'Best','FontSize',40,'HorizontalAlignment','center');   

Exp_Mp_value=0.68; %0.20 0.15 ;0.29; %0.68;
Exp_label   ='HILO'; %'EXP'; %'HILO'
xline(Exp_Mp_value,':b','LineWidth',10)
text(Exp_Mp_value+0.02,2.5,[{Exp_label};{[ num2str(Exp_Mp_value, '%1.2f') ' Nm/kg']}],'FontSize',40,'Color','b','HorizontalAlignment','left');  

set(gca, 'FontSize', 35)

% Add labels and legend
xlabel('M_P values [Nm/kg]')
ylabel('ΔE [%]')
% legend('no exo', [num2str(MP_list(1)) ' Nm'])
legend(legend_entries, 'Location', 'northeast','FontSize',35,'NumColumns',2);
legend boxoff
% ylim([-20 5])
ylim([-10 15])

title('Metabolic Rate Change vs. Peak Magnitudes')
% legend('Data points', 'Cubic fit', sprintf('Optimal (%.3f, %.3f)', optimal_MP, optimal_J), 'Location', 'best');
%% Time series metabolic cost
figure
for iMP=1:nMP
    hold on;
    plot(J_TS_list(iMP,:),'Color',color_list{iMP},'LineStyle','-','LineWidth',2);
end

%% Time series

% result_list={Results_SynCons Results_MP{1} Results_MP{2} Results_MP{3} Results_MP{4} Results_MP{5} Results_MP{6}};
% J_musAll_list={J_baselineUpd_extra J_mus_list{1} J_mus_list{2} J_mus_list{3} J_mus_list{4} J_mus_list{5} J_mus_list{6}};

result_list{1}=Results_SynCons;
for i=1:length(MP_list)
    result_list(1+i)=Results_MP(i);
end

J_musAll_list{1}=J_baselineUpd_extra;
for i=1:length(MP_list)
    J_musAll_list(1+i)=J_mus_list(i);
end

nRes=length(result_list);

version=1;

if version==1
    muscle_list=1:40; sp_lim = [5,8]; FontSize_sel=10;
elseif version==2
    muscle_list=[13 14 34 36]; sp_lim = [2,2]; FontSize_sel=25;
end

nMuscle_list=length(muscle_list);


MP_values = [0, MP_list]; % Add 0 for baseline
legend_entries = cell(1, length(MP_values) + 1);
for i = 1:length(MP_values)
    if i == 1
        legend_entries{i} = 'no exo';
    else
        legend_entries{i} = [num2str(MP_values(i)) ' Nm'];
    end
end
legend_entries{end} = 'optimal';

figure(51); clf; set(gcf,'WindowState', 'maximized','Color','white');
fSel=1:length(result_list{1}.MActivation.genericMRS)-1;
gait_cycle=linspace(0,100,length(fSel));
for iMus=1:nMuscle_list
    sMus=muscle_list(iMus);
    subplot(sp_lim(1),sp_lim(2),iMus);
    for iRes=1:nRes
        Results=result_list{iRes};
        
        outcome=Results.MActivation.genericMRS;
        % outcome=Results.TForce.genericMRS;

        FMo=Results.Misc.FMo;
        MActivation=Results.MActivation.genericMRS(:,fSel);
        
        FMltilde=Results.FMltilde.genericMRS(:,fSel);
        FMvtilde=Results.FMvtilde.genericMRS(:,fSel);
        % outcome=(FMo.*MActivation')'.*FMltilde.*FMvtilde;

        % outcome=J_musAll_list{iRes}.eDot.Muscles;

        plot(gait_cycle,outcome(sMus,fSel),'Color',color_list{iRes},'LineWidth',3);
        hold on

    end

    Results=result_list{idx_real};
    outcome=Results.MActivation.genericMRS;
    plot(gait_cycle,outcome(sMus,fSel),':b','LineWidth',5);

    title(Results_normal.MuscleNames{sMus}(1:end-2))
    ylim([0 1])
    set(gca, 'FontSize', FontSize_sel)

    if version==2
       if iMus == 1
        legend(legend_entries, 'Location', 'northeast', 'FontSize', 20,'NumColumns',2);
        legend boxoff
       end
    end
end
sgtitle('Muscle Activations across Peak Magnitudes','fontSize',30)

han = axes(figure(51), 'visible', 'off'); 
han.Title.Visible = 'on';
han.XLabel.Visible = 'on';
han.YLabel.Visible = 'on';

label_y = ylabel(han,'activations [ ]','FontSize',30); label_y.Position(1) = -0.05; label_y.Position(2) = 0.5;
label_x = xlabel(han,'gait cycle [%]'       ,'FontSize',30);   label_x.Position(1) = 0.5; label_x.Position(2) = -0.05;

%%
figure(52); clf; set(gcf,'WindowState', 'maximized','Color','white');

MP_list_tot = [0 MP_list];

version=1;

% Version settings
if version == 1
    nMuscles = 40;
    sp_lim = [5, 8];    FontSize_sel = 20;  legend_font_size = 10;
elseif version == 2
    nMuscles = length(muscle_list); % Use the muscle_list from your workspace
    sp_lim = [2, 2];    FontSize_sel = 30;  legend_font_size = 20;
end

% Create legend entries (same as before)
legend_entries = cell(1, length(MP_list_tot));
for i = 1:length(MP_list_tot)
    if i == 1
        legend_entries{i} = 'no exo';
    else
        legend_entries{i} = [num2str(MP_list_tot(i),'%1.0f') ' Nm'];
    end
end
% Add optimal entry
legend_entries{end+1} = 'optimal';

for iMus = 1:nMuscles
    if version == 1
        subplot(sp_lim(1), sp_lim(2), iMus);
        muscle_idx = iMus; % For version 1, use all muscles 1:40
    elseif version == 2
        subplot(sp_lim(1), sp_lim(2), iMus);
        muscle_idx = muscle_list(iMus); % For version 2, use specific muscles
    end
    
    outcome_mean = zeros(nRes, 1);
    
    for iRes = 1:nRes
        outcome = J_musAll_list{iRes}.eDot.Muscles;
        outcome_mean(iRes, 1) = mean(outcome(muscle_idx, :));
        plot(MP_list_tot(iRes)/Misc.subject_data.subject_mass, outcome_mean(iRes, 1), 'color', color_list{iRes}, ...
             'Marker', '.', 'MarkerSize', 50, 'LineWidth', 2,'LineStyle','none'); 
        hold on
    end
    
    % Plot optimal case
    optimal_mp = MP_list_tot(idx_real)/Misc.subject_data.subject_mass;
    optimal_outcome = mean(J_musAll_list{idx_real}.eDot.Muscles(muscle_idx, :));
    plot(optimal_mp, optimal_outcome, 'sb', 'MarkerSize', 30, 'LineWidth', 2);
    
    xline(0.68, ':b', 'LineWidth', 2)
    title(Results_normal.MuscleNames{muscle_idx}(1:end-2), 'FontSize', FontSize_sel)
    ylim([0 0.3])
    set(gca, 'FontSize', FontSize_sel)
    
    % Add legend only to first subplot
    if iMus == 1
        legend(legend_entries, 'Location', 'northeast', 'FontSize', legend_font_size, 'NumColumns', 2);
        legend boxoff
    end
    xlim([0 1.7])
end

% Add common labels
sgtitle('Average Muscle Metabolic Rates vs. Peak Magnitudes', 'FontSize', 30)

han = axes(figure(52), 'visible', 'off'); 
han.Title.Visible = 'on';
han.XLabel.Visible = 'on';
han.YLabel.Visible = 'on';

label_y = ylabel(han, 'metabolic rate [W]', 'FontSize', 30); 
label_y.Position(1) = -0.05; 
label_y.Position(2) = 0.5;

label_x = xlabel(han, 'M_p values [Nm]', 'FontSize', 30);   
label_x.Position(1) = 0.5; 
label_x.Position(2) = -0.05;

%%
figure(53); clf;
for i=1:9
    RTorque=result_list{i}.RActivation.genericMRS;

    for j=1:5
        subplot(2,3,j)
        plot(RTorque(j,:),'color',color_list{i}); hold on
    end
end
%% Compute stiffness

data_length=length(Results_noSyn.MActivation.genericMRS(1,:));
dFt_dlTtilde=zeros(data_length,40,nRes);

kT =Misc.kT;
FMo=Misc.FMo;

for iRes=1:nRes
    Result=result_list{iRes};

    lTtilde=Result.lTtilde.genericMRS;
    TForce =Result.TForce.genericMRS;

    for i=1:40
        dFt_dlTtilde(:,i,iRes) = FMo(i).* (kT(i)/5) * exp(kT(i) .* (lTtilde(i,:) - 0.995));
    end

end

%%
Kt=dFt_dlTtilde;

% Km=gradient(TForce)./gradient(lMtildeopt);
% 
% KMTU=(1/Km+1/Kt).^(-1);
%%
clf;
for i=1:40
    subplot(5,8,i); hold on
    % plot(lTtilde(i,:))
    plot(dFt_dlTtilde(:,i,1),'k','LineWidth',2)
    for iRes=1:nRes-1
        plot(dFt_dlTtilde(:,i,1+iRes),'color',color_list{iRes},'LineWidth',2)
    end
    ylim([0 10*10000])
end
%%
% figure; clf;

% result_list{1}=Results_normalSyn;
for iRes=1:nRes
Results=result_list{iRes};
figure('WindowState', 'maximized'); clf;

DatStore=DatStore_normal;

subject_mass=Results.Misc.subject_data.subject_mass;
NMuscles=length(Results.MuscleNames);
fe=Results.Misc.extra_frames;

fSel=1+fe:size(Results.MActivation.genericMRS,2)-1-fe;
data_length=length(fSel);

alphao=Results.Misc.params(4,:)';
lMo   =Results.Misc.params(2,:)';

w  =sin(alphao).*lMo;
lM =Results.lMtildeopt.genericMRS.*lMo;

alpha=zeros(NMuscles,data_length);
for mus_opt=1:NMuscles
    alpha(mus_opt,:)=asin(w(mus_opt)./lM(mus_opt,fSel));
end

Fpas_projected=zeros(NMuscles,data_length);
Fact_projected=zeros(NMuscles,data_length);
for mus_opt=1:NMuscles
    Fpe  =Results.Fpe.genericMRS(mus_opt,fSel);
    Fpas_projected(mus_opt,:)=Fpe.*cos(alpha(mus_opt,:));

    TForce  =Results.TForce.genericMRS(mus_opt,fSel);
    Fact_projected(mus_opt,:)=TForce-Fpas_projected(mus_opt,:);
end

to_plot=0;
if to_plot==1
    for mus_opt=1:NMuscles
    subplot(5,8,mus_opt)
    hold on
    % plot(Fpe,'k','LineWidth',2)
    % plot(Fpas_projected(mus_opt,:),':r','LineWidth',2)

    TForce=Results.TForce.genericMRS(mus_opt,fSel);

    plot(TForce,'k','LineWidth',2)
    plot(Fact_projected(mus_opt,:),'r','LineWidth',2)
    plot(Fpas_projected(mus_opt,:),'b','LineWidth',2)
    plot(Fact_projected(mus_opt,:)+Fpas_projected(mus_opt,:),':g','LineWidth',2)
    % ylim([0 2000])
    % legend('sum','act','pas')
    % legend boxoff
    title(Results_normal.MuscleNames{mus_opt})
    end
end

NMuscleGroups=10;
MG_Muscles_summed    =zeros(NMuscleGroups,data_length);
MG_Muscles_act_summed=zeros(NMuscleGroups,data_length);
MG_Muscles_pas_summed=zeros(NMuscleGroups,data_length);
MG_Muscles_act_summed_CO=zeros(NMuscleGroups/2,data_length);
MG_Muscles_pas_summed_CO=zeros(NMuscleGroups/2,data_length);

DofNames_Input=Results.Misc.DofNames_Input;
MGNames_Input ={'hipFlex' 'hipExt' 'hipAdd' 'hipAbd' 'hipIntRot' 'hipExtRot' 'kneeFlex' 'kneeExt' 'ankDor' 'ankPlantar'}; %+,-

for dof_opt=1:5
    MAinterp=DatStore.MAinterp(fSel,(dof_opt-1)*(NMuscles)+1:(dof_opt*NMuscles));
    TForce  =Results.TForce.genericMRS(:,fSel);

    % for muscles group -> NET
    MMuscles=MAinterp'.*TForce;
    MMuscles_p=MMuscles.*(tanh(1000*(MMuscles))*0.5+0.5);
    MMuscles_n=MMuscles.*(tanh(1000*(-MMuscles))*0.5+0.5);
    
    MG_Muscles_summed(dof_opt*2-1,:)=sum(MMuscles_p); % first positive
    MG_Muscles_summed(dof_opt*2,:)  =sum(MMuscles_n); % second negative

    to_plot=0;
    if to_plot==1
        for i=1:40
            subplot(5,8,i); hold on;
            plot(MMuscles(i,:),'k','LineWidth',2)
            plot(MMuscles_p(i,:),':r','LineWidth',2)
            plot(MMuscles_n(i,:),':b','LineWidth',2)
            plot(MMuscles_p(i,:)+MMuscles_n(i,:),':g','LineWidth',2)
            title(Results_normal.MuscleNames{i})
        end
    end

    % for muscles group -> ACT
    MG_Muscles_act=MAinterp'.*Fact_projected;
    MG_Muscles_act_p=MG_Muscles_act.*(tanh(1000*(MG_Muscles_act))*0.5+0.5);
    MG_Muscles_act_n=MG_Muscles_act.*(tanh(1000*(-MG_Muscles_act))*0.5+0.5);

    MG_Muscles_act_summed(dof_opt*2-1,:)=sum(MG_Muscles_act_p);
    MG_Muscles_act_summed(dof_opt*2,:)  =sum(MG_Muscles_act_n);

    % for muscles group -> PAS
    MG_Muscles_pas=MAinterp'.*Fpas_projected;
    MG_Muscles_pas_p=MG_Muscles_pas.*(tanh(1000*(MG_Muscles_pas))*0.5+0.5);
    MG_Muscles_pas_n=MG_Muscles_pas.*(tanh(1000*(-MG_Muscles_pas))*0.5+0.5);

    MG_Muscles_pas_summed(dof_opt*2-1,:)=sum(MG_Muscles_pas_p);
    MG_Muscles_pas_summed(dof_opt*2,:)  =sum(MG_Muscles_pas_n);
    % 
    MG_Muscles_act_summed_CO(dof_opt,:)=min([MG_Muscles_act_summed(dof_opt*2-1,:);-MG_Muscles_act_summed(dof_opt*2,:)]);
    MG_Muscles_pas_summed_CO(dof_opt,:)=min([MG_Muscles_pas_summed(dof_opt*2-1,:);-MG_Muscles_pas_summed(dof_opt*2,:)]);
end

gait_cycle = linspace(0,100,data_length)';
off_setlist=[0 0 0 0 -30];
% figure(option_result)
for dof_opt=1:5
    subplot(3,5,dof_opt)
    hold on;
    mg_sel_p=dof_opt*2-1;
    plot(gait_cycle,MG_Muscles_summed(mg_sel_p,:),'k','LineWidth',4)
    plot(gait_cycle,MG_Muscles_act_summed(mg_sel_p,:),'r','LineWidth',2)
    plot(gait_cycle,MG_Muscles_pas_summed(mg_sel_p,:),'b','LineWidth',2)
    % plot(MG_Muscles_act_summed(mg_sel,:)+MG_Muscles_pas_summed(mg_sel,:),':g','LineWidth',2)
    yline(0,'--k','LineWidth',2)
    
    mg_sel_n=dof_opt*2;
    plot(gait_cycle,MG_Muscles_summed(mg_sel_n,:),'k','LineWidth',4)
    plot(gait_cycle,MG_Muscles_act_summed(mg_sel_n,:),'r','LineWidth',2)
    plot(gait_cycle,MG_Muscles_pas_summed(mg_sel_n,:),'b','LineWidth',2)
    % plot(MG_Muscles_act_summed(mg_sel,:)+MG_Muscles_pas_summed(mg_sel,:),':g','LineWidth',2)
    yline(0,'--k','LineWidth',2)

    plot(gait_cycle,MG_Muscles_summed(mg_sel_p,:)+MG_Muscles_summed(mg_sel_n,:),':k','LineWidth',2)

    ylim([-100 100]+off_setlist(dof_opt))
    ylabel('torque [Nm]')
    xlabel('gait cycle [%]')
    title(DofNames_Input{dof_opt},'Interpreter','none')
    
    subplot(3,5,dof_opt+5)
    hold on;
    plot(gait_cycle,MG_Muscles_pas_summed_CO(dof_opt,:),'color',[205, 127, 50]/255,'LineWidth',1)
    plot(gait_cycle,-MG_Muscles_pas_summed_CO(dof_opt,:),'color',[205, 127, 50]/255,'LineWidth',1)
    data_input=MG_Muscles_pas_summed_CO(dof_opt,:);
    fill([gait_cycle' fliplr(gait_cycle')].',[data_input fliplr(-data_input)].',[205, 127, 50]/255, 'edgecolor', 'none', 'facealpha', 0.25);
    
    plot(gait_cycle,MG_Muscles_pas_summed(mg_sel_p,:),'b','LineWidth',2)
    plot(gait_cycle,MG_Muscles_pas_summed(mg_sel_n,:),'b','LineWidth',2)
    
    ylim([-20 20])
    ylabel('passive co-contraction torque [Nm]')
    xlabel('gait cycle [%]')
    title(DofNames_Input{dof_opt},'Interpreter','none')
    
    subplot(3,5,dof_opt+10)
    hold on;
    plot(gait_cycle,MG_Muscles_act_summed_CO(dof_opt,:),'color',[205, 127, 50]/255,'LineWidth',1)
    plot(gait_cycle,-MG_Muscles_act_summed_CO(dof_opt,:),'color',[205, 127, 50]/255,'LineWidth',1)

    data_input=MG_Muscles_act_summed_CO(dof_opt,:);
    fill([gait_cycle' fliplr(gait_cycle')].',[data_input fliplr(-data_input)].',[205, 127, 50]/255, 'edgecolor', 'none', 'facealpha', 0.25);
    
    plot(gait_cycle,MG_Muscles_act_summed(mg_sel_p,:),'r','LineWidth',2)
    plot(gait_cycle,MG_Muscles_act_summed(mg_sel_n,:),'r','LineWidth',2)
    
    ylim([-100 100]+off_setlist(dof_opt))
    ylabel('active co-contraction torque [Nm]')
    xlabel('gait cycle [%]')
    title(DofNames_Input{dof_opt},'Interpreter','none')    
end
set(gcf,'Color','white')
end
% figure(2)

%% Synergist
close all
Results=result_list{1};

[W,H,MActivation_Syns]=synergyAnalysis(Results,[1 1 1]);



function [moment_complete, gait_cycle] = interpolateFromStancePhase2GaitCyle(Misc_time,Misc_toeOff_time,extra_frame,moment_spline,to_plot)

% get time as in one frame equals 0.01 s
time        =Misc_time;
time_series =time(1):0.01:time(end);
data_length =length(time_series);

% conversion to gait cycle
data_GC      = data_length-2*extra_frame;
frames_per_GC= 100/(data_GC-1);
gait_cycle   = 0-frames_per_GC*extra_frame:frames_per_GC:100+frames_per_GC*extra_frame; % case without extra frame-> gait_cycle =linspace(0,100,length(time_series));

% toe off as a percentage of the gait cycle
toeOff_time=Misc_toeOff_time;
toeOff_asGC=(toeOff_time-time(1))/(time(end)-time(1))*100;
ind_toe    =find(gait_cycle>toeOff_asGC,1,'first')-1;

% interporlation
x_original = linspace(0, toeOff_asGC, length(moment_spline));
x_new      = linspace(0, toeOff_asGC, ind_toe-extra_frame);
v_inter = interp1(x_original,moment_spline,x_new,'cubic');
% clf; hold on; plot(x_new,v_inter,'r'); plot(x_original,moment_spline,'ok');

moment_complete=zeros(length(time_series),1);
moment_complete(1+extra_frame:ind_toe)=v_inter;

if to_plot
    figure;
    stance_phase=linspace(0,100,length(moment_spline));
    subplot(1,2,1); plot(stance_phase,moment_spline); xlim([0 100]);                         xlabel('stance phase [%]'); ylabel('moment [Nm]'); title('moment vs. stance phase')
    subplot(1,2,2); plot(gait_cycle,moment_complete); xlim([gait_cycle(1) gait_cycle(end)]); xlabel('gait cycle [%]');   ylabel('moment [Nm]'); title(['moment vs. gait cycle using ' num2str(extra_frame) ' extra frames'])
    xline(0,':k');  xline(gait_cycle(ind_toe),'--r'); xline(100,':k');
end
end
