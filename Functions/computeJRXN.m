function [JRXN]=computeJRXN(Results_normal,Misc,deviceInput,to_plot)
%% to unload from Misc and Results
model_path=Misc.model_path;
IK_path   =Misc.IKfile;
extLoads_file  =Misc.extLoadsInfo.fileName; 
% no need, I am using GRF as
% in unassisted conditions, TForces are obtained from Results and exoT is
% not needed - "deviceInput" does not do anything. I have to verify this
extLoads_setup =Misc.extLoadsInfo.setupName;
MotionSelection=Misc.MotionSelection;

MuscleNames = Results_normal.MuscleNames;
RANames     = appendName(Misc.DofNames_Input, '_reserve');

Time        = Results_normal.Time.genericMRS;
TForce      = Results_normal.TForce.genericMRS'; % value x muscles - check in case % TForce(:,[13 14])=TForce(:,[13 14])*0;
RATorque    = [Results_normal.RActivation.genericMRS'; zeros(1,length(RANames))];

data_length=length(Time);

if isempty(to_plot)
    to_plot.JRXN_summary=0;
end

extra_folder_name=Misc.extra_folder_name;
extra_file_name  =Misc.extra_file_name;

extra_folder_name_full=fullfile(Misc.OutPath,extra_folder_name);
if ~exist(extra_folder_name_full, 'dir');  mkdir(extra_folder_name_full);  end
%% to create a ForceFile for JRXN setup
allActuator_path=fullfile(Misc.OutPath,extra_folder_name,['MRS_ACT_' extra_file_name '.sto']);

label_name=load(fullfile(Misc.SetupPath,Misc.ForceLabel)); % SO_label_2392 OR SO_label_rajagopal
fullForce_label=label_name.SO_label;

MRS_label=[{'time'} MuscleNames RANames];

ind=zeros(1,length(MRS_label));
for i=1:length(MRS_label)
    ind(i)=find(strcmp(fullForce_label,MRS_label(i)));
end

data_full=zeros(data_length,length(fullForce_label));
data_full(:,ind)=[Time TForce RATorque];

% Update if deviceInput is NOT empty - CURRENTLY IT DOES NOT MATTER: 
% because "Special Types of Forces - TorqueActuator," see
% https://opensimconfluence.atlassian.net/wiki/spaces/OpenSim/pages/53089600/Joint+Reactions+Analysis
if ~ isempty(deviceInput)
    assistiveT=deviceInput.T;
    DOF       =[deviceInput.DOF '_reserve'];

    ind_DOF=strcmp(MRS_label,DOF);
    data_full(:,ind_DOF)=assistiveT;
end

generateMotFile_updated(data_full,fullForce_label,allActuator_path);

to_print_forcesOnly=0;
if to_print_forcesOnly==1
    outPath    =fullfile(Misc.OutPath,'TForce.sto');
    generateMotFile_updated([Time TForce],[{'time'} MuscleNames],outPath);

    TForce_max=max(TForce,[],'all');
    outPath_visual=fullfile(Misc.OutPath,'TForce_visual.sto');
    generateMotFile_updated([Time TForce/TForce_max],[{'time'} MuscleNames],outPath_visual);
end
%% to create "extLoads (_dev) SETUP" using "extLoads SETUP" with new "extLoads FILE"
extLoads_deviceSetup_path = fullfile(Misc.OutPath,extra_folder_name,['setup_updated_' extra_file_name '.xml']);
extLoadsDev_file_path     = extLoads_file;
generateExtLoadsFile(extLoads_setup,extLoads_deviceSetup_path,extLoadsDev_file_path); % generic setup , new (created) setup , new (create) extLoads file
% extLoads_deviceSetup_path=extLoads_setup;
%% to setup and run JRXN 
opt.computeReserveActuator=0; % no needed, it was added in ForceFile
opt.force_file_path       =allActuator_path; 
JRXN_defaultSetup_path    ='C:\Users\movea\Dropbox\Collaboration\Guna\GenericSetups\AnalysisJrxn_default_setup.xml';
JRXN_updatedSetup_path    =fullfile(Misc.OutPath,extra_folder_name,['JRXN_' MotionSelection extra_file_name '_setup.xml']);
JR_result_path            =fullfile(Misc.OutPath);

[result_outpath_full]=setup_and_run_JR(opt,model_path,IK_path,0,    extLoads_deviceSetup_path,JRXN_defaultSetup_path,  JRXN_updatedSetup_path,JR_result_path);
%% to compute joint reaction force (NET) and plot joint reaction force (NET & VECTORS) 
JRXN_data  = importdata(result_outpath_full);
JRXN_label = JRXN_data.colheaders;

side=Misc.gait_data.side_sel;

% variables of interest
JModNames = getJointModelNames(Misc.ForceLabel);
list_analysis_hip   = {[JModNames{1} '_' side '_on_femur_' side '_in_femur_' side '_fx'] ...
                       [JModNames{1} '_' side '_on_femur_' side '_in_femur_' side '_fy'] ...
                       [JModNames{1} '_' side '_on_femur_' side '_in_femur_' side '_fz']};
list_analysis_knee  = {[JModNames{2} '_' side '_on_tibia_' side '_in_tibia_' side '_fx'] ...
                       [JModNames{2} '_' side '_on_tibia_' side '_in_tibia_' side '_fy'] ...
                       [JModNames{2} '_' side '_on_tibia_' side '_in_tibia_' side '_fz']};
list_analysis_ank   = {[JModNames{3} '_' side '_on_talus_' side '_in_talus_' side '_fx'] ...
                       [JModNames{3} '_' side '_on_talus_' side '_in_talus_' side '_fy'] ...
                       [JModNames{3} '_' side '_on_talus_' side '_in_talus_' side '_fz']};

BW=Misc.subject_data.subject_mass*9.81; % body_weight=mass*gravity

% get indexes
ind_hip =zeros(1,3); ind_knee=zeros(1,3);   ind_ank =zeros(1,3);
for i=1:3
    ind_hip(i) =find(strcmp(JRXN_label,list_analysis_hip(i)));
    ind_knee(i)=find(strcmp(JRXN_label,list_analysis_knee(i)));
    ind_ank(i) =find(strcmp(JRXN_label,list_analysis_ank(i)));
end

% compute net force
JRXN_net_hip = sqrt(JRXN_data.data(:,ind_hip(1)).^2  + JRXN_data.data(:,ind_hip(2)).^2  + JRXN_data.data(:,ind_hip(3)).^2);
JRXN_net_knee= sqrt(JRXN_data.data(:,ind_knee(1)).^2 + JRXN_data.data(:,ind_knee(2)).^2 + JRXN_data.data(:,ind_knee(3)).^2);
JRXN_net_ank = sqrt(JRXN_data.data(:,ind_ank(1)).^2  + JRXN_data.data(:,ind_ank(2)).^2  + JRXN_data.data(:,ind_ank(3)).^2);
time=JRXN_data.data(:,1);

JRXN.time     = time;
JRXN.net_hip  = JRXN_net_hip;
JRXN.net_knee = JRXN_net_knee;
JRXN.net_ankle= JRXN_net_ank;

if to_plot.JRXN_summary==1

    y_lim_sel=[-3 10];
    clf; set(gcf,'Color','white')
    subplot(2,3,1)
    for i=1:3
        plot(time,-JRXN_data.data(:,ind_hip(i))/BW,'LineWidth',2); hold on;
    end
    legend(list_analysis_hip,'interpreter','none'); legend boxoff;
    title('hip vector forces','interpreter','none')
    xlabel('time [s]'); ylabel('body weight [ ]'); ylim(y_lim_sel)

    subplot(2,3,2)
    for i=1:3
        plot(time,-JRXN_data.data(:,ind_knee(i))/BW,'LineWidth',2); hold on;
    end
    legend(list_analysis_knee,'interpreter','none'); legend boxoff;
    title('knee vector forces','interpreter','none')
    xlabel('time [s]'); ylabel('body weight [ ]'); ylim(y_lim_sel)

    subplot(2,3,3)
    for i=1:3
        plot(time,-JRXN_data.data(:,ind_ank(i))/BW,'LineWidth',2); hold on;
    end
    legend(list_analysis_ank,'interpreter','none'); legend boxoff;
    title('ankle vector forces','interpreter','none')
    xlabel('time [s]'); ylabel('body weight [ ]'); ylim(y_lim_sel)

    subplot(2,3,4)
    plot(time,JRXN_net_hip/BW,'LineWidth',2);
    title('hip net forces','interpreter','none')
    xlabel('time [s]'); ylabel('body weight [ ]'); ylim(y_lim_sel)

    subplot(2,3,5)
    plot(time,JRXN_net_knee/BW,'LineWidth',2);
    title('knee net forces','interpreter','none')
    xlabel('time [s]'); ylabel('body weight [ ]'); ylim(y_lim_sel)

    subplot(2,3,6)
    plot(time,JRXN_net_ank/BW,'LineWidth',2);
    title('ankle net forces','interpreter','none')
    xlabel('time [s]'); ylabel('body weight [ ]'); ylim(y_lim_sel)
end

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
tool.setInitialTime(dataSource.getFirstTime);
tool.setFinalTime(dataSource.getLastTime);

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

function myLabel=getDevLabel(flag)
    if flag==0; myLabel='';
    elseif flag==1; myLabel='_withDevice';
    end
end

function generateExtLoadsFile(extLoads_devGenericSetup_path,extLoads_deviceSetup_path,extLoads_FILE)
% generate extLoads setup with device
extLoads_dev=xmlread(extLoads_devGenericSetup_path);             % open generic setup
extLoads_MOT_FILE=extLoads_dev.getElementsByTagName('datafile');
extLoads_MOT_FILE.item(0).setTextContent(extLoads_FILE);      % set extLoads device file
xmlwrite(extLoads_deviceSetup_path,extLoads_dev);
end

function newNames = appendName(oldNames, extraName)
    newNames = cellfun(@(name) [name extraName], oldNames, 'UniformOutput', false);
end

function jointModelNames = getJointModelNames(ForceLabel)
selModel=ForceLabel(10:end);
if strcmp(selModel,'2392')
    jointModelNames={'hip' 'knee' 'ankle'};
elseif strcmp(selModel,'rajagopal')
    jointModelNames={'hip' 'walker_knee' 'ankle'};
else
    error(['model not found :' selModel])
end
end