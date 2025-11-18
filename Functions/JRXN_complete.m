addpath(genpath('C:\Users\movea\Dropbox\Collaboration\Guna'));
% import org.opensim.modeling.*

% inputs
folderPath='C:\Users\movea\Dropbox\Collaboration\Guna\reactionForceData';

exp_path       = fullfile(folderPath,'expFiles');
setup_INpath   = fullfile(folderPath,'INPUTsetupFiles');
setup_OUTpath  = fullfile(folderPath,'OUTPUTsetupFiles');
result_path    = fullfile(folderPath,'results');

model_path     = fullfile(exp_path,'subject01_simbody.osim');
IK_path        = fullfile(exp_path,'subject01_walk1_ik.mot');
extLoads_path  = fullfile(exp_path,'subject01_walk1_grf.mot');      

% generated extLoads file with device torque
extLoadsDev_outpath  = fullfile(exp_path,'subject01_walk1_grf_device.mot');

% setup paths for unassisted conditions
extLoads_unassistedSetup_path = fullfile(setup_INpath,'extLoads_default.xml'); % for normal condition
extLoads_devGenericSetup_path = fullfile(setup_INpath,'extLoads_generic.xml'); % for assisted condition
SO_defaultSetup_path          = fullfile(setup_INpath,'AnalysisSO_default_setup.xml');
JRXN_defaultSetup_path        = fullfile(setup_INpath,'AnalysisJrxn_default_setup.xml');

% setup paths for device/assisted conditions          
extLoads_deviceSetup_path     = fullfile(setup_OUTpath,'extLoads_device.xml');
SO_deviceSetup_path           = fullfile(setup_OUTpath,'AnalysisSO_device_setup.xml');
JRXN_deviceSetup_path         = fullfile(setup_OUTpath,'AnalysisJrxn_device_setup.xml');

% setup generated without prior template
ID_defaultSetup_outpath       = fullfile(setup_OUTpath,'ID_default_setup.xml');
ID_deviceSetup_outpath        = fullfile(setup_OUTpath,'ID_device_setup.xml');

% results
ID_result_path = fullfile(result_path, 'Result_ID');
SO_result_path = fullfile(result_path, 'Result_SO');
JR_result_path = fullfile(result_path, 'Result_JR');


if ~exist(setup_OUTpath, 'dir');  mkdir(setup_OUTpath);  end
if ~exist(result_path, 'dir');    mkdir(result_path);    end
if ~exist(ID_result_path, 'dir'); mkdir(ID_result_path); end
if ~exist(SO_result_path, 'dir'); mkdir(SO_result_path); end
if ~exist(JR_result_path, 'dir'); mkdir(JR_result_path); end
%% Generate an extLoads file with assistive torque
% generate assistive torque
torque_features.t_r= 100;
torque_features.t_d= 40;
torque_features.M_p=-80;

to_plot.assistiveTorque     =1;
to_plot.extLoadsVerification=1;

% generate files
[device]=generateTorque_and_printExtLoadsFile(extLoads_path,extLoadsDev_outpath,torque_features,to_plot);
generateExtLoadsFile(extLoads_devGenericSetup_path,extLoads_deviceSetup_path,extLoadsDev_outpath);
%% Inverse dynamics
[ID_resultNormal_path]=setup_and_run_ID(model_path,IK_path,0,extLoads_unassistedSetup_path, ID_deviceSetup_outpath,ID_result_path);
[ID_resultDevice_path]=setup_and_run_ID(model_path,IK_path,1,    extLoads_deviceSetup_path,ID_defaultSetup_outpath,ID_result_path);

to_plot_results=1;

if to_plot_results==1
    figure(1); clf;
    list_joint={'hip_flexion_r_moment' 'hip_adduction_r_moment'	'hip_rotation_r_moment' 'knee_angle_r_moment' 'ankle_angle_r_moment'};

    IDdata_unassisted=importdata(ID_resultNormal_path);
    IDdata_device    =importdata(ID_resultDevice_path);

    colheaders=IDdata_unassisted.colheaders;
    mass=72.6; %kg

    for i=1:length(list_joint)
        to_plot_result([1,5,i],list_joint{i},list_joint(i),colheaders,IDdata_unassisted,mass,'-k',[-2 1]);
        to_plot_result([1,5,i],list_joint{i},list_joint(i),colheaders,IDdata_device,mass,':r',[-2 1]);
    end

    estimated_torque=(IDdata_unassisted.data(:,19)-IDdata_device.data(:,19))/mass;
    v_inter=interpolateQuick(device.time,device.torque,IDdata_unassisted.data(:,1))'/mass;

    subplot(1,5,5); hold on;
    plot(IDdata_unassisted.data(:,1),estimated_torque,':b','LineWidth',3, 'DisplayName','devT REAL');
    plot(device.time, device.torque/mass,':g','LineWidth',2, 'DisplayName','devT INPUT');
    plot(IDdata_unassisted.data(:,1), estimated_torque-v_inter,':m','LineWidth',2, 'DisplayName','ERROR');
end
%% Static optimization

to_run_computation=1;

if to_run_computation==1
opt.computeReserveActuator=1;
opt.ReserveActuator_path=fullfile(setup_INpath,'gait2354_CMC_Actuators.xml');

[SO_resultNormal_path]=setup_and_run_SO(opt,model_path,IK_path,0,extLoads_unassistedSetup_path,SO_defaultSetup_path, SO_deviceSetup_path,SO_result_path);
[SO_resultDevice_path]=setup_and_run_SO(opt,model_path,IK_path,1,    extLoads_deviceSetup_path,SO_defaultSetup_path,SO_defaultSetup_path,SO_result_path);
end

to_plot_results=1;

if to_plot_results==1
SOdata_unassisted=importdata(SO_resultNormal_path.force);
SOdata_device    =importdata(SO_resultDevice_path.force);

figure; clf;
muscle_list={'soleus_r' 'med_gas_r' 'bifemlh_r' 'bifemsh_r' 'grac_r' 'sar_r'};
muscle_names=SOdata_unassisted.colheaders;

ind=zeros(1,length(muscle_list));
for i=1:length(muscle_list); ind(i)=find(strcmp(muscle_names,muscle_list(i))); end

time=SOdata_unassisted.data(:,1);

for i=1:length(muscle_list)
    yyaxis left
    subplot(2,3,i); hold on
    plot(time,SOdata_unassisted.data(:,ind(i)),'k','LineWidth',2)
    plot(time,SOdata_device.data(:,ind(i)),':r','LineWidth',2)
    title(muscle_list(i),'Interpreter','none'); xlabel('# frames'); ylabel('activation [ ]'); 
    % ylim([0 1]); 
    xlim([time(1) time(end)])

    yyaxis right
    plot(IDdata_unassisted.data(:,1),estimated_torque,':b','LineWidth',2, 'DisplayName','devT REAL');
end
end
%% Joint reaction
clc;
opt.computeReserveActuator=1;
opt.ReserveActuator_path=fullfile(setup_INpath,'gait2354_CMC_Actuators.xml');

opt.force_file_path   =SO_resultNormal_path.force;
[JR_resultNormal_path]=setup_and_run_JR(opt,model_path,IK_path,0,extLoads_unassistedSetup_path,JRXN_defaultSetup_path, JRXN_defaultSetup_path,JR_result_path);

opt.force_file_path   =SO_resultDevice_path.force;
[JR_resultDevice_path]=setup_and_run_JR(opt,model_path,IK_path,1,    extLoads_deviceSetup_path,JRXN_defaultSetup_path,  JRXN_deviceSetup_path,JR_result_path);

%%
clc;

JRdata_unassisted=importdata(JR_resultNormal_path);
JRdata_device    =importdata(JR_resultDevice_path);

list_analysis_pelvisF={'ground_pelvis_on_pelvis_in_pelvis_fx' 'ground_pelvis_on_pelvis_in_pelvis_fy' 'ground_pelvis_on_pelvis_in_pelvis_fz'};
list_analysis_pelvisM={'ground_pelvis_on_pelvis_in_pelvis_mx' 'ground_pelvis_on_pelvis_in_pelvis_my' 'ground_pelvis_on_pelvis_in_pelvis_mz'};
list_analysis_hip    ={'hip_r_on_femur_r_in_femur_r_fx'       'hip_r_on_femur_r_in_femur_r_fy'       'hip_r_on_femur_r_in_femur_r_fz'};
list_analysis_knee   ={'knee_r_on_tibia_r_in_tibia_r_fx'      'knee_r_on_tibia_r_in_tibia_r_fy'      'knee_r_on_tibia_r_in_tibia_r_fz'};
list_analysis_ankle  ={'ankle_r_on_talus_r_in_talus_r_fx'     'ankle_r_on_talus_r_in_talus_r_fy'     'ankle_r_on_talus_r_in_talus_r_fz'};

colheaders=JRdata_unassisted.colheaders;

mass=72.6;
BW  =mass*9.8;

figure; clf;
subplot_pos=[2,2,1];
subplot_title='pelvis force';
subplot_data =list_analysis_pelvisF;
to_plot_result(subplot_pos,subplot_title,subplot_data,colheaders,JRdata_unassisted,BW,'-k',[0 10]);
to_plot_result(subplot_pos,subplot_title,subplot_data,colheaders,JRdata_device,BW,':r',[0 10]);

subplot_pos=[2,2,2];
subplot_title='pelvis moment';
subplot_data =list_analysis_pelvisM;
to_plot_result(subplot_pos,subplot_title,subplot_data,colheaders,JRdata_unassisted,BW,'-k',[0 10]);
to_plot_result(subplot_pos,subplot_title,subplot_data,colheaders,JRdata_device,BW,'-r',[0 10]);

subplot_pos=[2,3,4];
subplot_title='hip forces';
subplot_data =list_analysis_hip;
to_plot_result(subplot_pos,subplot_title,subplot_data,colheaders,JRdata_unassisted,BW,'-k',[-7 7]);
to_plot_result(subplot_pos,subplot_title,subplot_data,colheaders,JRdata_device,BW,':r',[-7 7]);

subplot_pos=[2,3,5];
subplot_title='knee forces';
subplot_data =list_analysis_knee;
to_plot_result(subplot_pos,subplot_title,subplot_data,colheaders,JRdata_unassisted,BW,'-k',[-7 7]);
to_plot_result(subplot_pos,subplot_title,subplot_data,colheaders,JRdata_device,BW,':r',[-7 7]);

subplot_pos=[2,3,6];
subplot_title='ankle forces';
subplot_data =list_analysis_ankle;
to_plot_result(subplot_pos,subplot_title,subplot_data,colheaders,JRdata_unassisted,BW,'-k',[-7 7]);
to_plot_result(subplot_pos,subplot_title,subplot_data,colheaders,JRdata_device,BW,':r',[-7 7]);

function [] = to_plot_result(subplot_pos,subplot_title,subplot_data,colheaders,JRF,BW,colorStyle,ylim_sel)
subplot(subplot_pos(1),subplot_pos(2),subplot_pos(3));hold on
time      =JRF.data(:,1);

ind=zeros(length(subplot_data),1);
for i=1:length(subplot_data)
    ind(i)=find(strcmp(colheaders,subplot_data{i}));
    plot(time,JRF.data(:,ind(i))/BW,colorStyle,'LineWidth',2);
end
xlim([time(1) time(end)])
ylim(ylim_sel)
legend(subplot_data,'Interpreter','none')
title(subplot_title)
end

function v_inter=interpolateQuick(simple_t,simple_M,time)
data_length= length(time);
x_original = linspace(simple_t(1), simple_t(end), length(simple_t));
x_new      = linspace(time(1), time(end), data_length);
v_inter    = interp1(x_original,simple_M,x_new,'cubic');
end

function myLabel=getDevLabel(flag)
    if flag==0; myLabel='';
    elseif flag==1; myLabel='_withDevice';
    end
end

function generateExtLoadsFile(extLoads_devGenericSetup_path,extLoads_deviceSetup_path,extLoadsDev_path)
% generate extLoads setup with device
extLoads_dev=xmlread(extLoads_devGenericSetup_path);             % open generic setup
extLoads_MOT_FILE=extLoads_dev.getElementsByTagName('datafile');
extLoads_MOT_FILE.item(0).setTextContent(extLoadsDev_path);      % set extLoads device file
xmlwrite(extLoads_deviceSetup_path,extLoads_dev);
end

function [result_outpath_full]=setup_and_run_ID(model_path,IK_path,setFlag,extLoads_deviceSetup_path,ID_deviceSetup_path,result_outpath)
import org.opensim.modeling.*

myLabel=getDevLabel(setFlag);
dataSource = Storage(IK_path);

result_outpath_full=fullfile(result_outpath,['ID' myLabel '.sto']);

ID=InverseDynamicsTool();
ID.setName(['ID' myLabel]);
ID.setCoordinatesFileName(IK_path);
ID.setStartTime(dataSource.getFirstTime);
ID.setEndTime(dataSource.getLastTime);
ID.setOutputGenForceFileName(result_outpath_full);
ID.set_results_directory(result_outpath);

model= Model(model_path);
ID.setModel(model);
ID.setModelFileName(model.getDocumentFileName());
ID.setLowpassCutoffFrequency(6);

%Set forces_to_exclude
excludedForces = ArrayStr();
excludedForces.append('Muscles');
ID.setExcludedForces(excludedForces);

ID.setExternalLoadsFileName(extLoads_deviceSetup_path);
ID.print(ID_deviceSetup_path);

% run
ID = InverseDynamicsTool(ID_deviceSetup_path);
ID.run;
end

function [outpath]=setup_and_run_SO(opt,model_path,IK_path,setFlag,    extLoads_deviceSetup_path,SO_defaultSetup_path,SO_deviceSetup_path,SO_result_path)
import org.opensim.modeling.*

myLabel=getDevLabel(setFlag);
dataSource=Storage(IK_path);

outpath_name=['SO' myLabel];
alg_name='staticOptimization';

tool=AnalyzeTool(SO_defaultSetup_path,false);
tool.setName(outpath_name);
tool.setModelFilename(model_path);
tool.setCoordinatesFileName(IK_path);
tool.setResultsDir(SO_result_path);
tool.setReplaceForceSet(0);
tool.setExternalLoadsFileName(extLoads_deviceSetup_path);

if opt.computeReserveActuator==1
    % to set reserve actuators
    forceSet = ArrayStr();
    forceSet.append(opt.ReserveActuator_path);
    tool.setForceSetFiles(forceSet);
end

AnalysisSet= tool.getAnalysisSet().get(0);
SO = StaticOptimization.safeDownCast(AnalysisSet);
SO.setName(alg_name);
SO.setStartTime(dataSource.getFirstTime);
SO.setEndTime(dataSource.getLastTime);
tool.print(SO_deviceSetup_path);

% run
tool2 = AnalyzeTool(SO_deviceSetup_path);
tool2.run();

outpath.activation=fullfile(SO_result_path,[outpath_name '_' alg_name '_activation.sto']);
outpath.controls  =fullfile(SO_result_path,[outpath_name '_' alg_name '_controls.xml']);
outpath.force     =fullfile(SO_result_path,[outpath_name '_' alg_name '_force.sto']);
end

function [result_outpath_full]=setup_and_run_JR(opt,model_path,IK_path,setFlag,    extLoads_deviceSetup_path,JRXN_defaultSetup_path,  JRXN_deviceSetup_path,JR_result_path)
import org.opensim.modeling.*

myLabel=getDevLabel(setFlag);
dataSource = Storage(IK_path) ;

outpath_name=['JRXN' myLabel];
alg_name='jointReaction';

result_outpath_full=fullfile(JR_result_path,[outpath_name '_' alg_name '_ReactionLoads.sto']);

% select input force file
force_file_path=opt.force_file_path;

% select joint reaction fields
jRxn.inFrame = 'child' ;
jRxn.onBody  = 'child' ;
jRxn.jointNames = 'all'; % to get a given JointName: osimModel.getJointSet.get(0)

% set Joint reaction fields
inFrame    = ArrayStr;
onBody     = ArrayStr;
jointNames = ArrayStr;
inFrame.set(0,jRxn.inFrame);
onBody.set(0,jRxn.onBody);
jointNames.set(0,jRxn.jointNames);

% set and run joint reaction
tool=AnalyzeTool(JRXN_defaultSetup_path,false);
tool.setName(outpath_name);
tool.setModelFilename(model_path);
tool.setCoordinatesFileName(IK_path);
tool.setResultsDir(JR_result_path);
tool.setReplaceForceSet(0);
tool.setExternalLoadsFileName(extLoads_deviceSetup_path);

% to set reserve actuators
if opt.computeReserveActuator==1
    forceSet = ArrayStr();
    forceSet.append(opt.ReserveActuator_path);
    tool.setForceSetFiles(forceSet);
end

AnalysisSet=tool.getAnalysisSet().get(0);
jointRxn = JointReaction.safeDownCast(AnalysisSet);
jointRxn.setName(alg_name);
jointRxn.setInFrame(inFrame);
jointRxn.setOnBody(onBody);
jointRxn.setJointNames(jointNames);
jointRxn.setStartTime(dataSource.getFirstTime);
jointRxn.setEndTime(dataSource.getLastTime);
jointRxn.setInDegrees(1);
jointRxn.setForcesFileName(force_file_path); %C:\Users\movea\Dropbox\PhDWork\ExperimentalData\subject1\100
tool.print(JRXN_deviceSetup_path);

% run
tool2 = AnalyzeTool(JRXN_deviceSetup_path);
tool2.run();
end

function [device] = generateTorque_and_printExtLoadsFile(originalExtLoadsPath,newExtLoadsPath,torque_features,to_plot)

% myLabel=getDevLabel(1); % setFlag = 1;
% [filepath,name,ext] = fileparts(originalExtLoadsPath);
% newExtLoadsPath=fullfile(filepath,[name myLabel ext]);

GRF_data=importdata(originalExtLoadsPath);

time=GRF_data.data(:,strcmp(GRF_data.colheaders,'time'));
data_length=length(time);

% assistive torque
t_r=torque_features.t_r;
dur=torque_features.t_d;
M_p=torque_features.M_p;
x=linspace(0,1,dur);             y= sin(pi*x)*M_p;

simple_M(1:t_r)=zeros(1,t_r);
simple_M(t_r+1:t_r+dur)=y;
simple_M(t_r+dur+1:250)=zeros(1,250-(t_r+dur));
simple_t=linspace(time(1),time(end),length(simple_M));

% assistive torque
if to_plot.assistiveTorque==1
    figure;
    M_inter=interpolateQuick(simple_t,simple_M,time)';
    plot(simple_t,simple_M,'k','LineWidth',2); hold on
    plot(time,M_inter,':r','LineWidth',2);
end
device.time  =time;
device.torque=M_inter;

zero_values=zeros(data_length,1);
M_inter_full_1= [zero_values zero_values M_inter];
M_inter_full_2=-[zero_values zero_values M_inter];

assistive_torque_names ={'devTorque1_x' 'devTorque1_y' 'devTorque1_z' 'devTorque2_x' 'devTorque2_y' 'devTorque2_z'};
assistive_torque_values=[M_inter_full_1 M_inter_full_2];

generateMotFile_updated([GRF_data.data assistive_torque_values],[GRF_data.colheaders assistive_torque_names],newExtLoadsPath)

% verify
if to_plot.extLoadsVerification==1
    GRF_data_new=importdata(newExtLoadsPath);
    parts = split(GRF_data_new.textdata(9), char(9));

    list_eval={'devTorque1_x' 'devTorque1_y' 'devTorque1_z' 'devTorque2_x' 'devTorque2_y' 'devTorque2_z'};

    time_new=GRF_data_new.data(:,1);

    figure;
    for i=1:length(list_eval)
        ind=find(strcmp(parts,list_eval{i}));
        subplot(2,3,i)
        plot(time_new,GRF_data_new.data(:,ind))
        title(list_eval{i})
    end
end
end

function []=generateExternalLoads_v2(loads_Setup,folderPath)
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