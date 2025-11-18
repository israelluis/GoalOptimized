addpath(genpath('C:\Users\movea\Dropbox\Collaboration\Guna'));
import org.opensim.modeling.*

% INPUTS
folderPath='C:\Users\movea\Dropbox\Collaboration\Guna\reactionForceData';

modelPath     = fullfile(folderPath,'subject01_simbody.osim');
outputFilePath= folderPath;
forceFilePath = fullfile(folderPath,'subject01_StaticOptimization_force.sto');
motionFilePath= fullfile(folderPath,'subject01_walk1_ik.mot');
externalLoadSetupPath= fullfile(folderPath,'subject01_walk1_grf.xml');
AnalyzeTool_generic_file=fullfile(folderPath,'AnalysisJrxn_generic_setup.xml');

externalLoadPath_nor=fullfile(folderPath,'subject01_walk1_grf.mot');
externalLoadPath_dev=fullfile(folderPath,'subject01_walk1_grf_withDev.mot');

withDevices=1;

externalLoadPath = externalLoadPath_nor;

if withDevices==1 % update externalLoadPath
    externalLoadPath = externalLoadPath_dev;
end
%%
% SETUP
% select GRF file
generate_GRF_file=1;
if generate_GRF_file==0  % none is generated
else
    % if withDevices==0
        loads_setup.externalForceName      ={'right' 'left'};
        loads_setup.applied_to_body        ={'calcn_r' 'calcn_l'};
        loads_setup.force_expressed_in_body={'ground' 'ground'};
        loads_setup.force_identifier       ={'ground_force_v' '1_ground_force_v'};
        loads_setup.point_expressed_in_body={'ground' 'ground'};
        loads_setup.point_identifier       ={'ground_force_p' '1_ground_force_p'};
        loads_setup.torque_identifier      ={'ground_torque_' '1_ground_torque_'};
    % else
    %     loadsSetup.externalForceName      ={'right'   'left'    'exoFoot' 'exoShank'};
    %     loadsSetup.applied_to_body        ={'calcn_r' 'calcn_l' 'calcn_r' 'tibia_r' };
    %     loadsSetup.force_expressed_in_body={'ground'  'ground'  'calcn_r' 'calcn_r' };
    %     loadsSetup.force_identifier       ={'ground_force_v' '1_ground_force_v' '' ''};
    %     loadsSetup.point_expressed_in_body={'ground' 'ground' 'ground' 'ground'};
    %     loadsSetup.point_identifier       ={'ground_force_p' '1_ground_force_p' '' ''};
    %     loadsSetup.torque_identifier      ={'ground_torque_' '1_ground_torque_' 'devTorque1' 'devTorque2'};
    % end
    loads_setup.externalLoadPath=externalLoadPath;
    % externalLoadSetupPath=generateExternalLoads(loads_setup,folderPath);
    generateExternalLoads(loads_setup,folderPath);

    
end
%% Inverse dynamics
import org.opensim.modeling.*

dataSource = Storage(motionFilePath) ;

ID=InverseDynamicsTool();
ID.setName('ID');
ID.setCoordinatesFileName(motionFilePath);
ID.setStartTime(dataSource.getFirstTime);
ID.setEndTime(dataSource.getLastTime);

model= Model(modelPath);
ID.setModel(model);
ID.setModelFileName(model.getDocumentFileName());
ID.setLowpassCutoffFrequency(6);

%Set forces_to_exclude
excludedForces = ArrayStr();
excludedForces.append('Muscles');
ID.setExcludedForces(excludedForces);

externalLoadSetupPath='C:\Users\movea\Dropbox\Collaboration\Guna\reactionForceData\grf.xml';
ID.setExternalLoadsFileName(externalLoadSetupPath);

finalSetupFile=fullfile(outputFilePath, 'ID_setup.xml');
ID.print(finalSetupFile);

ID = InverseDynamicsTool(finalSetupFile);
ID.run;
%% Generate an assistive torque
clc; clf;
GRF_data=importdata(externalLoadPath_nor);

time=GRF_data.data(:,strcmp(GRF_data.colheaders,'time'));
data_length=length(time);

% assistive torque
t_r= 80;
dur= 60;
M_p=-50;

x=linspace(0,1,dur);             y= sin(pi*x)*M_p;

simple_M(1:t_r)=zeros(1,t_r);
simple_M(t_r+1:t_r+dur)=y;
simple_M(t_r+dur+1:250)=zeros(1,250-(t_r+dur));
simple_t=linspace(time(1),time(end),length(simple_M));

M_inter=interpolateQuick(simple_t,simple_M,time)';
plot(simple_t,simple_M,'k','LineWidth',2); hold on
plot(time,M_inter,':r','LineWidth',2);

zero_values=zeros(data_length,1);
M_inter_full_1= [zero_values zero_values M_inter];
M_inter_full_2=-[zero_values zero_values M_inter];

assistive_torque_names ={'devTorque1_x' 'devTorque1_y' 'devTorque1_z' 'devTorque2_x' 'devTorque2_y' 'devTorque2_z'};
assistive_torque_values=[M_inter_full_1 M_inter_full_2];
% plot(x,y)
generateMotFile_updated([GRF_data.data assistive_torque_values],[GRF_data.colheaders assistive_torque_names],fullfile(outputFilePath,'subject01_walk1_grf_withDev.mot'))

% verify
to_verify_data=0;

if to_verify_data==1
    GRF_data_new=importdata(fullfile(outputFilePath,'subject01_walk1_grf_withDev.mot'));
    parts = split(GRF_data_new.textdata(9), char(9));

    list_eval={'devTorque1_x' 'devTorque1_y' 'devTorque1_z' 'devTorque2_x' 'devTorque2_y' 'devTorque2_z'};

    time_new=GRF_data_new.data(:,1);
    for i=1:length(list_eval)
        ind=find(strcmp(parts,list_eval{i}));
        subplot(2,3,i)
        plot(time_new,GRF_data_new.data(:,ind))
        title(list_eval{i})
    end
end
%%
IDdata_nor=importdata(fullfile(outputFilePath,'inverse_dynamics_normal.sto'));

list_joint={'hip_flexion_r_moment' 'hip_adduction_r_moment'	'hip_rotation_r_moment' 'knee_angle_r_moment' 'ankle_angle_r_moment'};
colheaders=IDdata_nor.colheaders;
mass=1; %72.6

clf; figure(1);
for i=1:length(list_joint)
to_plot_result([1,5,i],list_joint{i},list_joint(i),colheaders,IDdata_nor,mass,'-k');
end

IDdata=importdata(fullfile(outputFilePath,'inverse_dynamics.sto'));

list_joint={'hip_flexion_r_moment' 'hip_adduction_r_moment'	'hip_rotation_r_moment' 'knee_angle_r_moment' 'ankle_angle_r_moment'};
colheaders=IDdata.colheaders;
% mass=72.6;

for i=1:length(list_joint)
to_plot_result([1,5,i],list_joint{i},list_joint(i),colheaders,IDdata,mass,':r');
end
subplot(1,5,5); hold on;
plot(IDdata_nor.data(:,1),IDdata_nor.data(:,19)-IDdata.data(:,19),':b','LineWidth',3);

% simple_t=linspace(IDdata.data(1,1),IDdata.data(end,1),length(simple_M));
% M_inter =interpolateQuick(simple_t,simple_M,IDdata.data(:,1))/mass;
subplot(1,5,5); hold on;
plot(time, M_inter,':g','LineWidth',2);
% plot(time,M_inter/mass,':r','lineWidth',2)
%% Joint reaction
dataSource = Storage(motionFilePath) ;

% Joint Reaction Fields
jRxn.inFrame = 'child' ;
jRxn.onBody  = 'child' ;
jRxn.jointNames = 'all'; % to get a given JointName: osimModel.getJointSet.get(0)

inFrame    = ArrayStr;
onBody     = ArrayStr;
jointNames = ArrayStr;
inFrame.set(0,jRxn.inFrame);
onBody.set(0,jRxn.onBody);
jointNames.set(0,jRxn.jointNames);

tool=AnalyzeTool(AnalyzeTool_generic_file,false); %tool=AnalyzeTool(path_generic_file,false);
tool.setName('AnalysisResult');
tool.setModelFilename(modelPath);
tool.setCoordinatesFileName(motionFilePath);
tool.setResultsDir(outputFilePath);
tool.setReplaceForceSet(0);
tool.setExternalLoadsFileName(externalLoadSetupPath);

AnalysisSet=tool.getAnalysisSet().get(0);
jointRxn = JointReaction.safeDownCast(AnalysisSet);
jointRxn.setName('JointRxn');
jointRxn.setInFrame(inFrame);
jointRxn.setOnBody(onBody);
jointRxn.setJointNames(jointNames);
jointRxn.setStartTime(dataSource.getFirstTime);
jointRxn.setEndTime(dataSource.getLastTime);
jointRxn.setInDegrees(1);
jointRxn.setForcesFileName(forceFilePath); %C:\Users\movea\Dropbox\PhDWork\ExperimentalData\subject1\100
finalSetupFile=fullfile(outputFilePath, 'AnalysisJrxn_updated_setup.xml');
tool.print(finalSetupFile);

% run
tool2 = AnalyzeTool(finalSetupFile);
tool2.run();

%%
clc;
import org.opensim.modeling.*
path_result=fullfile(folderPath,'AnalysisResult_JointRxn_ReactionLoads.sto');
dataSource = Storage(fullfile(folderPath,'AnalysisResult_JointRxn_ReactionLoads.sto')) ;

JRF=importdata(path_result);

list_analysis_pelvisF={'ground_pelvis_on_pelvis_in_pelvis_fx' 'ground_pelvis_on_pelvis_in_pelvis_fy' 'ground_pelvis_on_pelvis_in_pelvis_fz'};
list_analysis_pelvisM={'ground_pelvis_on_pelvis_in_pelvis_mx' 'ground_pelvis_on_pelvis_in_pelvis_my' 'ground_pelvis_on_pelvis_in_pelvis_mz'};
list_analysis_hip    ={'hip_r_on_femur_r_in_femur_r_fx'       'hip_r_on_femur_r_in_femur_r_fy'       'hip_r_on_femur_r_in_femur_r_fz'};
list_analysis_knee   ={'knee_r_on_tibia_r_in_tibia_r_fx'      'knee_r_on_tibia_r_in_tibia_r_fy'      'knee_r_on_tibia_r_in_tibia_r_fz'};
list_analysis_ankle  ={'ankle_r_on_talus_r_in_talus_r_fx'     'ankle_r_on_talus_r_in_talus_r_fy'     'ankle_r_on_talus_r_in_talus_r_fz'};

colheaders=JRF.colheaders;

mass=72.6;
BW  =mass*9.8;

clf;
subplot_pos=[2,2,1];
subplot_title='pelvis force';
subplot_data =list_analysis_pelvisF;
to_plot_result(subplot_pos,subplot_title,subplot_data,colheaders,JRF,BW);

subplot_pos=[2,2,2];
subplot_title='pelvis moment';
subplot_data =list_analysis_pelvisM;
to_plot_result(subplot_pos,subplot_title,subplot_data,colheaders,JRF,BW);

subplot_pos=[2,3,4];
subplot_title='hip forces';
subplot_data =list_analysis_hip;
to_plot_result(subplot_pos,subplot_title,subplot_data,colheaders,JRF,BW);

subplot_pos=[2,3,5];
subplot_title='knee forces';
subplot_data =list_analysis_knee;
to_plot_result(subplot_pos,subplot_title,subplot_data,colheaders,JRF,BW);

subplot_pos=[2,3,6];
subplot_title='ankle forces';
subplot_data =list_analysis_ankle;
to_plot_result(subplot_pos,subplot_title,subplot_data,colheaders,JRF,BW);

function [] = to_plot_result(subplot_pos,subplot_title,subplot_data,colheaders,JRF,BW,colorStyle)
subplot(subplot_pos(1),subplot_pos(2),subplot_pos(3));hold on
time      =JRF.data(:,1);

ind=zeros(length(subplot_data),1);
for i=1:length(subplot_data)
    ind(i)=find(strcmp(colheaders,subplot_data{i}));
    plot(time,JRF.data(:,ind(i))/BW,colorStyle,'LineWidth',2);
end
xlim([time(1) time(end)])
legend(subplot_data,'Interpreter','none')
title(subplot_title)
end

function v_inter=interpolateQuick(simple_t,simple_M,time)
data_length= length(time);
x_original = linspace(simple_t(1), simple_t(end), length(simple_t));
x_new      = linspace(time(1), time(end), data_length);
v_inter    = interp1(x_original,simple_M,x_new,'cubic');
end

function []=generateExternalLoads(loads_Setup,folderPath)
import org.opensim.modeling.*

externalForceName=loads_Setup.externalForceName;
applied_to_body=loads_Setup.applied_to_body;
force_expressed_in_body=loads_Setup.force_expressed_in_body;
force_identifier=loads_Setup.force_identifier;
point_expressed_in_body=loads_Setup.point_expressed_in_body;
point_identifier=loads_Setup.point_identifier;
torque_identifier=loads_Setup.torque_identifier;

externalLoadPath =loads_Setup.externalLoadPath;

myLoads = ExternalLoads();
for i = 1:length(applied_to_body) % nExternalForces
    newForce = ExternalForce() ;
    newForce.setName(externalForceName{i}) ;
    newForce.set_applied_to_body(applied_to_body{i}) ;
    newForce.set_force_expressed_in_body(force_expressed_in_body{i}) ;
    newForce.set_force_identifier(force_identifier{i}) ;
    newForce.set_point_expressed_in_body(point_expressed_in_body{i}) ;
    newForce.set_point_identifier(point_identifier{i}) ;
    newForce.set_torque_identifier(torque_identifier{i}) ;
    newForce.setDataSource(Storage(externalLoadPath)) ;

    newForce_cast = ExternalForce.safeDownCast(newForce);

    myLoads.adoptAndAppend(newForce_cast);
end
myLoads.setDataFileName(externalLoadPath);

% pathChosen=fullfile(folderPath,'grf_withDeviceLoads.xml');
% myLoads.print(pathChosen);

pathChosen=[];

externalLoadSetupPath=pathChosen;
end
%%
% modelFile = strcat(modelPath);
% osimModel = Model(modelFile);
% 
% state = osimModel.initSystem() ; % Re-build the model before printing
% jointRxn = JointReaction() ;
% % jointRxn.setName('JointRxn') ;
% % jointRxn.setInFrame(inFrame) ;
% % jointRxn.setOnBody(onBody) ;
% jointRxn.setJointNames(jointNames) ;
% jointRxn.setStartTime(dataSource.getFirstTime);
% jointRxn.setEndTime(dataSource.getLastTime);
% % osimModel.addAnalysis(jointRxn) ;
% jointRxn.setModel(osimModel) ;
% jointRxn.setForcesFileName(forceFilePath); %C:\Users\movea\Dropbox\PhDWork\ExperimentalData\subject1\100
% % jointRxn.set_statesStore('IK_sub1mod6act4GCMOT1.mot');
% jointRxn.begin(state) ;
% % jointRxn.printResults('results_JointReaction.sto') ;
% jointRxn.end(state) ;
% % jointRxn.print(fullfile(outputFilePath, 'JrxnSetup.xml')) ;
% % jointRxn.printResults(fullfile(outputFilePath,'result4')) ;


%%
% import org.opensim.modeling.*
% 
% GRF_data=Storage(externalLoadPath);
% timeColumn=ArrayDouble();
% timeSpam=GRF_data.getTimeColumn(timeColumn);
% newData1 =linspace(GRF_data.getFirstTime,GRF_data.getLastTime,timeSpam);
% newData2 =linspace(GRF_data.getFirstTime,GRF_data.getLastTime,timeSpam);
% 
% nColumns=GRF_data.getColumnLabels.getSize;
% 
% % labels = ArrayStr();
% % for i=1:nColumns
% %     ColumnLabels(i)={char(GRF_data.getColumnLabels.get(i-1))};
% %     labels.append(ColumnLabels{i});
% % end
% 
% % Convert to ArrayDouble for OpenSim
% newData1_vec = ArrayDouble();
% newData2_vec = ArrayDouble();
% for i = 1:timeSpam
%     newData1_vec.append(newData1(i));
%     newData2_vec.append(newData2(i));
% end
% 
% % v = Vector(newData1_vec);
% 
% % labels.append('newData1');
% % labels.append('newData2');
% % Append columns
% % GRF_data.append('newData1', newData1_vec);
% % GRF_data.append(newData1_vec);
% % GRF_data.append('newData2', newData2_vec);
% 
% % v = Vector(timeSpam,1);
% % GRF_data.append(timeSpam, v);
% % GRF_data.setColumnLabels(labels);
% % 
% 
% % GRF_data.print('newGRF.mot');
% % v = Vector(newTime, newTime2);
% % GRF_data.append(1.0, v);
% 
% % GRF_data.append('newTime');
% % GRF_data.setDataColumn('newTime',newTime)
% 
% % GRF_data.getDataColumn('time')
% % GRF_data_down.getData.(0)
% 
% % https://simtk.org/plugins/phpBB/viewtopicPhpbb.php?f=91&t=8143&p=21980&start=0&view=