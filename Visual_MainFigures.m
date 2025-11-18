folderPath = 'C:\Users\movea\Dropbox\PostDoc Work\Simulation\Code\ProjectResults\Codesign\Pilot\sub1\v2_t1';  % Replace with your folder path
pattern = '*Summary_VERX_eDot_MCLU24*';
files = dir(fullfile(folderPath, pattern));

% Display the results
for i = 1:length(files)
    fprintf('Found: %s\n', files(i).name);
end

% Get just the file names
fileNames = {files.name};
%%
close all
for i=1:8
    figure(i)
    img = imread(fileNames{i});
    % Display image with proper scaling
    imshow(img, 'Border', 'tight');
    
    % Maximize the figure window
    set(gcf, 'WindowState', 'maximized'); % Maximize to screen size
    % Alternative: set(gcf, 'WindowState', 'maximized'); % For newer MATLAB versions
    
    % Improve image quality
    axis image;      % Maintain aspect ratio
    axis off;        % Remove axes if it's just an image
    % truesize;        % Display at true resolution
    
    % Add title with filename
    [~, name, ext] = fileparts(fileNames{i});
    title(sprintf('Image %d: %s%s', i, name, ext), 'Interpreter', 'none');
end
%%
folderPath        = 'C:\Users\movea\Dropbox\PostDoc Work\Simulation\Code\ProjectResults\Codesign\Pilot\sub1\v2_t1\iters200';  % Replace with your folder path

% TYPE RESULTS:    'Bayesopt'          'MRSOptimal'
% OPTIMAL CRITERA: '_VERX_eDot_MCLU24' '_VERX_RJXN_knee'
% Metabolic cost
set_optimal_criteria= '_VERX_eDot_MCLU24'; 
set_type_result     = 'MRSOptimal';       
fileNames=get_files_from_pattern(folderPath,[set_type_result set_optimal_criteria]);
MRS_sel_metCost=cell(8);
for i=1:8
    MRS_sel_metCost(i)={load(fileNames{i})};
end

set_type_result   = 'Bayesopt';       
fileNames_bayes_metCost=get_files_from_pattern(folderPath,[set_type_result set_optimal_criteria]);
Bayesiopt_metCost=cell(8);
for i=1:8
    Bayesiopt_metCost(i)={load(fileNames_bayes_metCost{i}).resultsBayesopt};
end

fileNames_bayes_metCost_new = erase(fileNames_bayes_metCost, [set_type_result [set_optimal_criteria '_'] '_iters200' '_knee_' ".mat"]);
order={'KE(Q)' 'KE(A)' 'AP(Q)' 'AP(A)' 'AP(Q)&KE(Q)' 'AP(A)&KE(Q)' 'AP(Q)&KE(A)' 'AP(A)&KE(A)'};
[~, reorder_indices_bayes_metCost] = ismember(order, fileNames_bayes_metCost_new);

% knee contact forces
set_optimal_criteria= '_VERX_RJXN_knee'; 
set_type_result     = 'MRSOptimal';       
fileNames=get_files_from_pattern(folderPath,[set_type_result set_optimal_criteria]);
MRS_sel_kneeContact=cell(8);
for i=1:8
    MRS_sel_kneeContact(i)={load(fileNames{i})};
end

set_type_result   = 'Bayesopt';       
fileNames_bayes_kneeContact=get_files_from_pattern(folderPath,[set_type_result set_optimal_criteria]);
Bayesiopt_kneeContact=cell(8);
for i=1:8
    Bayesiopt_kneeContact(i)={load(fileNames_bayes_kneeContact{i}).resultsBayesopt};
end

fileNames_bayes_metCost_new = erase(fileNames_bayes_kneeContact, [set_type_result [set_optimal_criteria '_'] '_iters200' '_knee_' ".mat"]);
order={'KE(Q)' 'KE(A)' 'AP(Q)' 'AP(A)' 'AP(Q)&KE(Q)' 'AP(A)&KE(Q)' 'AP(Q)&KE(A)' 'AP(A)&KE(A)'};
[~, reorder_indices_bayes_kneeContact] = ismember(order, fileNames_bayes_metCost_new);
%% Compute relevant metrics: Positive Power & other objectives
PowerPos_nor_sum_eDot=compute_power(MRS_sel_metCost);
PowerPos_nor_sum_JRXN=compute_power(MRS_sel_kneeContact);

assistance_goal='eDot_MCLU24';
J_baseline_eDot=computeBaseline(assistance_goal);
J_cost_eDot_fromJRXN=zeros(8,1);
for iTrial=1:8 %8
    MRS_sel=MRS_sel_kneeContact{iTrial};
    [J_sel,Jrel_sel,J_extra_sel] = computeOuterLoopFunction(MRS_sel.Misc,MRS_sel.Results,assistance_goal);
    J_value=(J_sel-J_baseline_eDot)/J_baseline_eDot*100;
    J_cost_eDot_fromJRXN(iTrial)=J_value;
end

assistance_goal='RJXN_knee';
J_baseline_eDot=computeBaseline(assistance_goal);
J_cost_JRXN_fromeDot=zeros(8,1);
for iTrial=1:8 %8
    MRS_sel=MRS_sel_metCost{iTrial};
    [J_sel,Jrel_sel,J_extra_sel] = computeOuterLoopFunction(MRS_sel.Misc,MRS_sel.Results,assistance_goal);
    J_value=(J_sel-J_baseline_eDot)/J_baseline_eDot*100;
    J_cost_JRXN_fromeDot(iTrial)=J_value;
end

% J_value=(J_sel-J_baseline)/J_baseline*100;

%%
MinObjective_ordered_eDot=zeros(8,1);
MinObjective_ordered_JRXN=zeros(8,1);
for i=1:8
    MinObjective_ordered_eDot(i,1)=Bayesiopt_metCost{reorder_indices_bayes_metCost(i)}.MinObjective;
    MinObjective_ordered_JRXN(i,1)=Bayesiopt_kneeContact{reorder_indices_bayes_kneeContact(i)}.MinObjective;
end
colors_hex = {
    '#1E3A8A', ... % KE(Q) - Deep navy blue
    '#3B82F6', ... % KE(A) - Bright blue
    '#047857', ... % AP(Q) - Deep green
    '#10B981', ... % AP(A) - Emerald green
    '#7C3AED', ... % AP(Q)&KE(Q) - Royal purple
    '#6366F1', ... % AP(A)&KE(Q) - Indigo
    '#DC2626', ... % AP(Q)&KE(A) - Deep red
    '#EA580C'  ... % AP(A)&KE(A) - Orange
};
colors_rgb = hex2rgb(colors_hex);
%%
clf; figure(1);
set(gcf,'color','w');
PowerPos_nor_sum_ord_eDot=PowerPos_nor_sum_eDot(reorder_indices_bayes_metCost);
PowerPos_nor_sum_ord_JRXN=PowerPos_nor_sum_JRXN(reorder_indices_bayes_kneeContact);

J_cost_eDot_fromJRXN_ord=J_cost_eDot_fromJRXN(reorder_indices_bayes_kneeContact);
J_cost_JRXN_fromeDot_ord=J_cost_JRXN_fromeDot(reorder_indices_bayes_metCost);

subplot(4,3,1)
X = categorical(order);
X = reordercats(X,order);
for j=1:8
    metric=-MinObjective_ordered_eDot(j);
    bar(j,metric,'FaceColor',colors_hex{j}); hold on

    text(j,metric+4,[num2str(metric,'%1.0f') '%'],'HorizontalAlignment','center','FontSize',15)
end
xticks(1:length(X));    xticklabels(X); ylim([0 35])
ylabel('metabolic cost reduction [%]'); xlabel('conditions');
set(gca,'FontSize',12);
title('Performance','FontSize',15);

subplot(4,3,4)
X = categorical(order);
X = reordercats(X,order);
for j=1:8
    metric=PowerPos_nor_sum_ord_eDot(j);
    bar(j,metric,'FaceColor',colors_hex{j}); hold on

    text(j,metric+2,[num2str(metric,'%1.1f')],'HorizontalAlignment','center','FontSize',15)
end
xticks(1:length(X));    xticklabels(X); ylim([0 15])
ylabel('hardware burden  [ ]'); xlabel('conditions');
set(gca,'FontSize',12);
title('Hardware burden','FontSize',15);

visual_off=[1 -1 -1 1 1 1 1 -1]; visual_marker=['s' 's' 's' 's' 'o' 'o' 'o' 'o'];
subplot(2,3,2); hold on
for i=1:8
    x=PowerPos_nor_sum_ord_eDot(i);
    y=-MinObjective_ordered_eDot(i);
    plot(PowerPos_nor_sum_ord_eDot(i),y,visual_marker(i),'MarkerSize',15,'MarkerFaceColor',colors_hex{i},'MarkerEdgeColor','none')
    text(x,y+1.5*visual_off(i),order{i},'FontSize',12,'HorizontalAlignment','center')
end
yline(0); xline(0,':k'); xlim([-0.5 10]); ylim([-15 30])
ylabel('metabolic cost reduction [%]');
xlabel('hardware burden [ ]'); %+MechPower/maxMechPower
set(gca,'FontSize',12);
title('Performance (direct) vs Burden','FontSize',15);

visual_off=[1 -1 -1 1 1 1 1 -1]; visual_marker=['s' 's' 's' 's' 'o' 'o' 'o' 'o'];
subplot(2,3,3); hold on
for i=1:8
    x=PowerPos_nor_sum_ord_eDot(i);
    y=-J_cost_JRXN_fromeDot_ord(i);
    plot(PowerPos_nor_sum_ord_eDot(i),y,visual_marker(i),'MarkerSize',15,'MarkerFaceColor',colors_hex{i},'MarkerEdgeColor','none')
    % text(x,y+1.5*visual_off(i),order{i},'FontSize',12,'HorizontalAlignment','center')
end
yline(0); xline(0,':k'); xlim([-0.5 10]); ylim([0 60])
ylabel('knee contact force reduction [%]');
xlabel('hardware burden [ ]'); %+MechPower/maxMechPower
set(gca,'FontSize',12);
title('Performance (indirect) vs Burden','FontSize',15);


subplot(4,3,7)
X = categorical(order);
X = reordercats(X,order);
for j=1:length(MinObjective_ordered_JRXN)
    metric=-MinObjective_ordered_JRXN(j);
    bar(j,metric,'FaceColor',colors_hex{j}); hold on

    text(j,metric+7,[num2str(metric,'%1.0f') '%'],'HorizontalAlignment','center','FontSize',15)
end
xticks(1:length(X));        xticklabels(X); ylim([0 75])
ylabel('knee contact force reduction [%]');  xlabel('conditions');
set(gca,'FontSize',12);
title('Performance','FontSize',15);


subplot(4,3,10)
X = categorical(order);
X = reordercats(X,order);
for j=1:8
    metric=PowerPos_nor_sum_ord_JRXN(j);
    bar(j,metric,'FaceColor',colors_hex{j}); hold on

    text(j,metric+2,[num2str(metric,'%1.1f')],'HorizontalAlignment','center','FontSize',15)
end
xticks(1:length(X));    xticklabels(X); ylim([0 15])
ylabel('hardware burden  [ ]'); xlabel('conditions');
set(gca,'FontSize',12);
title('Hardware burden','FontSize',15);

visual_off=[-1 1 1 1 1 1 1 -1]; visual_marker=['s' 's' 's' 's' 'o' 'o' 'o' 'o'];
subplot(2,3,6); hold on
for i=1:8
    x=PowerPos_nor_sum_ord_JRXN(i);
    y=-MinObjective_ordered_JRXN(i);
    plot(PowerPos_nor_sum_ord_JRXN(i),y,visual_marker(i),'MarkerSize',15,'MarkerFaceColor',colors_hex{i},'MarkerEdgeColor','none')
    text(x,y+3*visual_off(i),order{i},'FontSize',12,'HorizontalAlignment','center')
end
yline(0); xline(0,':k'); xlim([-0.5 10])
ylabel('knee contact force reduction [%]');
xlabel('hardware burden [ ]');
set(gca,'FontSize',12);
title('Performance (direct) vs Burden','FontSize',15);


visual_off=[-1 1 1 1 1 1 1 -1]; visual_marker=['s' 's' 's' 's' 'o' 'o' 'o' 'o'];
subplot(2,3,5); hold on
for i=1:8
    x=PowerPos_nor_sum_ord_JRXN(i);
    y=-J_cost_eDot_fromJRXN_ord(i);
    plot(PowerPos_nor_sum_ord_JRXN(i),y,visual_marker(i),'MarkerSize',15,'MarkerFaceColor',colors_hex{i},'MarkerEdgeColor','none')
    % text(x,y+3*visual_off(i),order{i},'FontSize',12,'HorizontalAlignment','center')
end
yline(0); xline(0,':k'); xlim([-0.5 10]); ylim([-15 30])
ylabel('metabolic cost reduction [%]');
xlabel('hardware burden [ ]');
set(gca,'FontSize',12);
title('Performance (indirect) vs Burden','FontSize',15);

function rgb = hex2rgb(hex_cells)
    rgb = zeros(length(hex_cells), 3);
    for i = 1:length(hex_cells)
        hex = hex_cells{i};
        rgb(i,:) = [hex2dec(hex(2:3)), hex2dec(hex(4:5)), hex2dec(hex(6:7))] / 255;
    end
end

function fileNames=get_files_from_pattern(folderPath,label_pattern_list)
label_pattern=label_pattern_list;
pattern = ['*' label_pattern '*.mat'];

files = dir(fullfile(folderPath, pattern));
% Display the results
for i = 1:length(files)
    fprintf('Found: %s\n', files(i).name);
end

% Get just the file names
fileNames = {files.name};
end

function PowerPos_nor_sum=compute_power(MRS_sel)
PowerPos_trial=zeros(1,8);
PowerNet_trial=zeros(1,8);
PowerPos_list_all=zeros(8,2); % maximum number of devices;
PowerPos_list_name=strings(8,2);
for iTrial=1:8
    Device=MRS_sel{iTrial}.Misc.Device;
    DOFNames=MRS_sel{iTrial}.DatStore.DOFNames;
    time=MRS_sel{iTrial}.DatStore.Mesh.t'; % sec

    nDevs=length(Device);
    PowerPos_list=zeros(1,nDevs);
    PowerNet_list=zeros(1,nDevs);
    for iDev=1:nDevs
        Torque=Device{iDev}.Assistance.Profile.Torque; % Nm
        devJoint=Device{iDev}.MuscleGroup{1};
        devDir  =Device{iDev}.MuscleGroup{2};
        actType =Device{iDev}.Type{1};

        q_exp=MRS_sel{iTrial}.DatStore.q_exp(:,strcmp(DOFNames,devJoint)); % deg

        dq_exp=gradient(q_exp)./gradient(time);

        dq_exp_rad=dq_exp*pi/180;
        Power=dq_exp_rad.*Torque*devDir;

        PowerPos = trapz(time(Power>0),Power(Power>0));
        PowerPos_list(1,iDev)=PowerPos;

        PowerNet = trapz(time,Power);
        PowerNet_list(1,iDev)=PowerNet;

        if strcmp(actType,'active')
            PowerPos_list_act=PowerPos;
            abb='A';
        elseif strcmp(actType,'quasi-passive')
            PowerPos_list_act=0;
            abb='Q';
        end
        PowerPos_list_all(iTrial,iDev)=PowerPos_list_act;
        PowerPos_list_name(iTrial,iDev)=abb;

        to_plot=0;
        if to_plot==1
            subplot(3,nDevs,(1-1)*nDevs+iDev);
            plot(time,q_exp);
            ylabel('angle [deg]')
            subplot(3,nDevs,(2-1)*nDevs+iDev);
            plot(time,Torque);
            ylabel('torque [Nm]')
            subplot(3,nDevs,(3-1)*nDevs+iDev);
            plot(time,Power);
            ylabel('power [W]')
            title(['positive Power = ' num2str(PowerPos,'%1.1f') 'J'])

            subplot(3,nDevs,(1-1)*nDevs+iDev);
            title(Device{iDev}.Type{1});
        end

    end
    PowerPos_trial(iTrial)=sum(PowerPos_list);
    PowerNet_trial(iTrial)=sum(PowerNet_list);
end


ind_Q=PowerPos_list_name=="Q";

PowerPos_list_all(ind_Q)=1; % 1 J. Assuming. Solenoid 3.6 W and time = 100 ms (10% gait cycle), 0.36 J (=1) https://www.farnell.com/datasheets/3625674.pdf?_gl=1*10n0fj0*_gcl_aw*R0NMLjE3NjI4NTE0NjcuQ2p3S0NBaUEyc3ZJQmhCLUVpd0FSV0RQam9FUHptYUtHQ09acDBsSmJjVGZYRFBYcngzcE55eERwWDJrMTQ1Y1VBRFVqRXpiSHhfZEtSb0NGblVRQXZEX0J3RQ..*_gcl_au*MjA1MTA3NTE4My4xNzYyODUxNDY4
% https://se.farnell.com/adafruit/412/frame-type-d-type-box-type/dp/2816374?gross_price=true&CMP=KNC-GSE-GEN-Shopping-Pmax-Mid-ROAS&mckv=_dc|pcrid||&gad_source=1&gad_campaignid=22205790637&gbraid=0AAAAAD8yeHnAayAiNFDTwidIwx6m8KhNs&gclid=CjwKCAiA2svIBhB-EiwARWDPjoEPzmaKGCOZp0lJbcTfXDPXrx3pNyxDpX2k145cUADUjEzbHx_dKRoCFnUQAvD_BwE
max_val=max(sum(PowerPos_list_all,2),[],'all');


% PowerPos_list_all(ind_Q)=0.1*max_val;

PowerPos_nor=PowerPos_list_all/max_val*10;
% PowerPos_nor(ind_Q)=1;
PowerPos_nor_sum=sum(PowerPos_nor,2);
end

function J_baseline = computeBaseline(assistance_goal)
currentFolder = pwd;

% select subject and motion
info.SubjectSelection ='sub1';        % select subject
info.MotionSelection  ='v2_t1';       % read as-> v(velocity)== 1[slow], 2[normal], 3[fast]) _ t (trial)== 1st, 2nd, or 3rd trial
info.currentFolder    = currentFolder;
% assistance_goal       ='eDot_MCLU24'; % eDot eDot_MCLU24 gasForces KJMusForces RJXN_knee RJXN_knee_par
OutName_optimal       ='_VERX';

% setup unassisted condition
[Misc,Results_normal,DatStore_normal]=setup_and_run_unassisted_condition(info);
%% Check results
Misc.MotionSelection   ='oneSingleMotion';
Misc.ForceLabel        ='SO_label_rajagopal';
Misc.extra_folder_name ='temp';
Misc.extra_file_name   ='unassisted';

[J_baseline,Jrel_TS,J_extra_normal] = computeOuterLoopFunction(Misc,Results_normal,assistance_goal); % this is my starting point, the Edot at unassisted conditions.
end