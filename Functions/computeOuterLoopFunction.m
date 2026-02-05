function [J,Jrel_TS,J_extra]=computeOuterLoopFunction(Misc,Results,devGoal)

% window of analysis
extra_frames=Misc.extra_frames;
fSel=1+extra_frames:size(Results.MActivation,2)-1-extra_frames;
t_range=[Results.Time(fSel(1)) Results.Time(fSel(end))];

% number of muscles
nMuscles=length(Results.MuscleNames);

% get the mass of the subject using the function GetModelMass
subject_mass  =Misc.subject_data.subject_mass;       % or modelmass = getModelMass(Misc.model_path);
subject_height=Misc.subject_data.subject_height;

% Order results
[musT_param,musE_param,mus_time,mus_dyn,states_field] = Results_states_params(Results,subject_mass,subject_height,extra_frames);

data_length=length(fSel);

goal_category=devGoal(1:4);
goal_specific=devGoal(6:end);

% Goal: Metabolic cost ----------------------------------------------------
if strcmp(goal_category,'eDot')
    % For Bhargava model
    if strcmp(goal_specific,'MCBH04')
    mode_basalOn=0; % Boolean: Select if you want to add a basal heat
    mode_negWork=1; % Boolean: Select if you want negative heat to be dissipated, thus net energy rate is >=0

    muscle_metRate=zeros(nMuscles,data_length);
    for mus_opt=1:nMuscles
        mus_Dynamics= [mus_time(:)'; squeeze(mus_dyn(:,mus_opt,:))];
        [MC_parameter_BH04,~,~,~]= MC_BH04_R(musT_param(:,mus_opt),musE_param(:,mus_opt),mus_Dynamics,mode_negWork,mode_basalOn);
        muscle_metRate(mus_opt,:)= MC_parameter_BH04(1,:)/subject_mass; % normalized by mass
    end

    % For Bhargava model + adjustment to eccentric contraction
    elseif strcmp(goal_specific,'MCLU24')
    mode_basalOn=0; % Boolean: Select if you want to add a basal heat
    mode_negWork=1; % Boolean: Select if you want negative heat to be dissipated, thus net energy rate is >=0

    muscle_metRate=zeros(nMuscles,data_length);
    for mus_opt=1:nMuscles
        mus_Dynamics= [mus_time(:)'; squeeze(mus_dyn(:,mus_opt,:))];
        [MC_parameter,~,~,~]= MC_LU24_R(musT_param(:,mus_opt),musE_param(:,mus_opt),mus_Dynamics,mode_negWork,mode_basalOn);
        muscle_metRate(mus_opt,:)= MC_parameter(1,:)/subject_mass; % normalized by mass
    end    
    end

    % computation of mean values
    Jrel_TS   =sum(muscle_metRate); % normalized by mass
    Edot_mean =mean(Jrel_TS)*2; % 2 legs
    aim_main  =Edot_mean;

    helper_need=0; % if you want to penalize use of reserve actuators
    if helper_need==1
        NDof      =size(Results.RTorque,1);
        helper_obs=sum(abs(Results.RTorque(:,fSel)))/NDof; % RA can be possitive (agonist) or negative (antagonist)
        aim_helper=sum(helper_obs,'all');
    else
        aim_helper=0;
    end
    J=aim_main+0.1*aim_helper;

    J_extra.eDot.Muscles=muscle_metRate;
    J_extra.label ='metabolic cost';
    J_extra.unit  ='W/kg';

% Goal: Joint contact forces ----------------------------------------------
elseif strcmp(goal_category,'JRXN')

    % enable paralization
    if strcmp(goal_specific(end-2:end),'par')
        % Lock file path (unique to your application)
        lockFile = 'computeOuterLoopFunction.lock';
        maxWaitTime = 60;  % Timeout (seconds)

        % Wait for lock release
        startTime = tic;
        while exist(lockFile, 'file')
            pause(0.1);  % Avoid busy-waiting
            if toc(startTime) > maxWaitTime
                error('Timeout waiting for lock.');
            end
        end

        % Acquire lock
        fid = fopen(lockFile, 'w');
        fclose(fid);

        target_joint=goal_specific(1:end-4);
    else
        target_joint=goal_specific;
    end

    % setup the file 
    if ~isfield(Misc,'DeviceFeatures') || isempty(Misc.DeviceFeatures)
    else
        nParams=length(Misc.DeviceFeatures.Params.values);
        text_label=[];
        for i=1:nParams
            text_label=[text_label 'v' num2str(i) '_' num2str(Misc.DeviceFeatures.Params.values(i),'%1.0f') '_'];
        end
        Misc.extra_file_name  =text_label;
    end

    % compute joint reaction forces
    [JRXN]=computeJRXN(Results,Misc,[],[]);
    pause(0.2); % just in case

    % recompute fSel
    fSel_1=find(JRXN.time>=t_range(1),1,'first');
    fSel_2=find(JRXN.time>=t_range(end),1,'first');
    fSel  =fSel_1:fSel_2;

    JRXN_net_joint=JRXN.(['net_' target_joint])(fSel);
    time        =JRXN.time(fSel);
    BW          =subject_mass*9.81;

    toeOff_event=find(time>Misc.gait_data.toeOff_time,1,'first');

    Jrel_TS=JRXN_net_joint/BW;
    JRXN_net_joint_stance=JRXN_net_joint(1:toeOff_event)/BW;
    J      =mean(JRXN_net_joint_stance);

    J_extra.label =[target_joint ' contact force'];
    J_extra.unit  ='N/BW';

    if strcmp(goal_specific(end-2:end),'par')
        % Release lock
        delete(lockFile);
    end
% elseif strcmp(devGoal,'gasForces')
%     musNames={['gasmed_' Misc.gait_data.side_sel] ['gaslat_' Misc.gait_data.side_sel]};
%     nSelMuscles=length(musNames);
%     ind=zeros(1,nSelMuscles);
%     for i=1:nSelMuscles
%         ind(i)=find(strcmp(Results.MuscleNames,musNames(i)));
%     end
% 
%     tendonForces=Results.TForce(ind,fSel);
%     % for i=1:nGasMuscles; subplot(1,2,i); plot(gasForces(i,:)); end
%     Jrel_TS=tendonForces;
%     J      =sum(mean(tendonForces,2));
% 
%     J_extra.label =' ';
%     J_extra.unit  =' ';
% elseif strcmp(devGoal,'KJMusForces')
%     musNames={['bflh_' Misc.gait_data.side_sel]    ['bfsh_' Misc.gait_data.side_sel] ...
%               ['gasmed_' Misc.gait_data.side_sel]  ['gaslat_' Misc.gait_data.side_sel]...
%               ['grac_' Misc.gait_data.side_sel]    ['sart_' Misc.gait_data.side_sel] ...
%               ['semimem_' Misc.gait_data.side_sel] ['semiten_' Misc.gait_data.side_sel] ...
%               ['recfem_' Misc.gait_data.side_sel]  ['vasint_' Misc.gait_data.side_sel] ...
%               ['vasmed_' Misc.gait_data.side_sel]  ['vaslat_' Misc.gait_data.side_sel]};
%     nSelMuscles=length(musNames);
%     ind=zeros(1,nSelMuscles);
%     for i=1:nSelMuscles
%         ind(i)=find(strcmp(Results.MuscleNames,musNames(i)));
%     end
% 
%     tendonForces=Results.TForce(ind,fSel);
%     % for i=1:nGasMuscles; subplot(1,2,i); plot(gasForces(i,:)); end
%     Jrel_TS=tendonForces;
%     J      =sum(mean(tendonForces,2));
% 
%     J_extra.label =' ';
%     J_extra.unit  =' ';
% elseif strcmp(devGoal,'KJActs')
%     musNames={['bflh_' Misc.side_sel] ['bfsh_' Misc.side_sel] ...
%               ['gasmed_' Misc.side_sel] ['gaslat_' Misc.side_sel]...
%               ['grac_' Misc.side_sel] ['sart_' Misc.side_sel] ...
%               ['semimem_' Misc.side_sel] ['semiten_' Misc.side_sel] ...
%               ['recfem_' Misc.side_sel] ['vasint_' Misc.side_sel] ...
%               ['vasmed_' Misc.side_sel] ['vaslat_' Misc.side_sel]};
% 
%     nSelMuscles=length(musNames);
%     ind=zeros(1,nSelMuscles);
%     for i=1:nSelMuscles
%         ind(i)=find(strcmp(Results.MuscleNames,musNames(i)));
%     end
% 
%     activation=Results.MActivation(ind,fSel);
%     % for i=1:nGasMuscles; subplot(1,2,i); plot(gasForces(i,:)); end
%     Jrel_TS=activation;
%     J      =sum(mean(activation,2));
% 
%     J_extra.label =' ';
%     J_extra.unit  =' ';
% elseif strcmp(devGoal,'RJXN_knee_par')
% 
%     % Lock file path (unique to your application)
%     lockFile = 'computeOuterLoopFunction.lock';
%     maxWaitTime = 60;  % Timeout (seconds)
% 
%     % Wait for lock release
%     startTime = tic;
%     while exist(lockFile, 'file')
%         pause(0.1);  % Avoid busy-waiting
%         if toc(startTime) > maxWaitTime
%             error('Timeout waiting for lock.');
%         end
%     end
% 
%     % Acquire lock
%     fid = fopen(lockFile, 'w');
%     fclose(fid);
% 
%     % as in RJXN_knee
%     if ~isfield(Misc,'DeviceFeatures') || isempty(Misc.DeviceFeatures)
%     else
%         nParams=length(Misc.DeviceFeatures.Params.values);
%         text_label=[];
%         for i=1:nParams
%             text_label=[text_label 'v' num2str(i) '_' num2str(Misc.DeviceFeatures.Params.values(i),'%1.0f') '_'];
%         end
%         Misc.extra_file_name  =text_label;
%     end
% 
%     [JRXN]=computeJRXN(Results,Misc,[],[]);
% 
%     % Release lock
%     delete(lockFile);
% 
%     % recompute fSel
%     fSel_1=find(JRXN.time>=t_range(1),1,'first');
%     fSel_2=find(JRXN.time>=t_range(end),1,'first');
%     fSel  =fSel_1:fSel_2;
% 
%     JRXN_net_joint=JRXN.netKnee(fSel);
%     time        =JRXN.time(fSel);
%     BW          =Misc.subject_data.subject_mass*9.81;
% 
%     toeOff_event=find(time>Misc.gait_data.toeOff_time,1,'first');
% 
%     Jrel_TS=JRXN_net_joint/BW;
%     JRXN_net_joint_stance=JRXN_net_joint(1:toeOff_event)/BW;
%     J      =mean(JRXN_net_joint_stance);
% 
%     J_extra.label ='knee contact force';
%     J_extra.unit  ='N/BW';
% elseif strcmp(devGoal,'RJXN_knee')
% 
%     if ~isfield(Misc,'DeviceFeatures') || isempty(Misc.DeviceFeatures)
%     else
%         nParams=length(Misc.DeviceFeatures.Params.values);
%         text_label=[];
%         for i=1:nParams
%             text_label=[text_label 'v' num2str(i) '_' num2str(Misc.DeviceFeatures.Params.values(i),'%1.0f') '_'];
%         end
%         Misc.extra_file_name  =text_label;
%     end
% 
%     [JRXN]=computeJRXN(Results,Misc,[],[]);
%     pause(0.2); % just in case
% 
%     % recompute fSel
%     fSel_1=find(JRXN.time>=t_range(1),1,'first');
%     fSel_2=find(JRXN.time>=t_range(end),1,'first');
%     fSel  =fSel_1:fSel_2;
% 
%     JRXN_net_joint=JRXN.netKnee(fSel);
%     time        =JRXN.time(fSel);
%     BW          =Misc.subject_data.subject_mass*9.81;
% 
%     toeOff_event=find(time>Misc.gait_data.toeOff_time,1,'first');
% 
%     Jrel_TS=JRXN_net_joint/BW;
%     JRXN_net_joint_stance=JRXN_net_joint(1:toeOff_event)/BW;
%     J      =mean(JRXN_net_joint_stance);
% 
%     J_extra.label ='knee contact force';
%     J_extra.unit  ='N/BW';
else
    warning('no aim was defined')
end

end
