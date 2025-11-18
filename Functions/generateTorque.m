function [assistanceInfo]=generateTorque(Device,DatStore,time,addedFrames)

typeProfile=Device.Type{2};
Params     =Device.Params;
jointSel   =Device.MuscleGroup{1};
jointDir   =Device.MuscleGroup{2};
if strcmp(typeProfile,'spline#N3')
    nNodes      =3;
    nodesParam_x=[Params(2) Params(1) Params(3)];
    nodesParam_y=[Params(4)];

    [moment_spline,nodes_x_all,nodes_y_all,infoNodes]=generateSplineBasedAssistance(nNodes,nodesParam_x,nodesParam_y,0);
    [moment_inter,gaitCycle_inter,time_inter]=interpolate2GaitCycle(time,addedFrames,moment_spline,0);
    
    profile.Torque  =moment_inter';
    profile.GaitCycle=gaitCycle_inter';
    profile.Time   =time_inter';

    param.Values=[infoNodes.nodes_x_upd(2) infoNodes.nodes_x_upd(1) infoNodes.nodes_x_upd(3) infoNodes.nodes_y_upd(1)];

    [varNames,varUnits]=getLabelAssistance(jointSel(1:end-2),jointDir,typeProfile);
    param.Names =varNames;
    param.Units =varUnits;

    additional.Name='paramatric time-based torque control';
    additional.Nodes_x=nodes_x_all;
    additional.Nodes_y=nodes_y_all;
elseif strcmp(typeProfile,'clutchSpring')    
    % DOFSel=Device.MuscleGroup{1};
    sDOF=strcmp(DatStore.DOFNames,jointSel);
    qInter=DatStore.IKinterp(:,sDOF);
    [gait_cycle,~] = computeGC(time,addedFrames);
    
    tC=Params(1);
    K =Params(2);
    
    qInterPos=-qInter.*jointDir; % this '-' sign secures all angles strictly are positive within the region of interets. Verify with this for i=1:5; subplot(2,3,i); plot(DatStore.IKinterp(:,i)); end
    
    itC_engaged=find(gait_cycle>=tC,1);
    qInter_tC_engaged=qInterPos(itC_engaged);

    iLowerThan_qtC=find(qInterPos<qInter_tC_engaged);
    iLowerThan_qtC_laterThanitC=iLowerThan_qtC(iLowerThan_qtC>itC_engaged);


    if any(iLowerThan_qtC_laterThanitC)
        itC_disengaged=iLowerThan_qtC_laterThanitC(1)-1;
    else
        itC_disengaged=itC_engaged;
    end

    qInter_tC_disengaged=qInterPos(itC_disengaged);

    iON=itC_engaged:itC_disengaged;

    qSel=qInterPos(iON);
    qSelOffset=qSel-qSel(1);

    qSelExo=zeros(length(gait_cycle),1);
    qSelExo(iON,1)=qSelOffset;

    qSelExo(qSelExo<0)=0; % this secures that all angles are positive, likely the last angle might be (tiny) negative
    torque=qSelExo.*K';

    to_plot=0;
    if to_plot==1
        subplot(2,1,1)
        plot(gait_cycle,qInterPos); hold on;
        yline(qInter_tC_engaged,':b'); xline(tC,':b');
        xline(gait_cycle(itC_disengaged))
        plot(gait_cycle(iON),qInterPos(iON),':r','LineWidth',2);

        subplot(2,1,2)
        plot(gait_cycle,qSelExo); hold on;
        plot(gait_cycle,torque)
    end

    [gait_cycle,time_series] = computeGC(time,addedFrames);
    profile.Torque  =torque;
    profile.GaitCycle=gait_cycle';
    profile.Time   =time_series';

    param.Values=[tC K];

    % (jointSel,jointDir,assistanceType)

    [varNames,varUnits]=getLabelAssistance(jointSel(1:end-2),jointDir,typeProfile);
    param.Names =varNames;
    param.Units =varUnits;

    additional.Name='paramatric quasi-passive spring-based torque control';
    additional.Angle=qInterPos;
    additional.it_eng=itC_engaged;
    additional.it_dis=itC_disengaged;
end

assistanceInfo.Profile=profile;
assistanceInfo.Param  =param;
assistanceInfo.Additional  =additional;
end

function [moment_complete, gait_cycle,time_series] = interpolate2GaitCycle(Misc_time,extra_frame,moment_spline,to_plot)

[gait_cycle,time_series]=computeGC(Misc_time,extra_frame);

% interporlation
x_original = linspace(0, 100, length(moment_spline));
x_new      = linspace(gait_cycle(1), gait_cycle(end), length(gait_cycle));
v_inter = interp1(x_original,moment_spline,x_new,'cubic');
% clf; hold on; plot(x_new,v_inter,'r'); plot(x_original,moment_spline,'ok');

v_inter(1:1+extra_frame-1)=zeros(1,extra_frame);
v_inter(end-extra_frame+1:end)=zeros(1,extra_frame);

moment_complete=v_inter;

if to_plot
    figure;
    gait_cycle_original=linspace(0,100,length(moment_spline));
    subplot(1,1,1); hold on;
    plot(gait_cycle_original,moment_spline,'k','LineWidth',2); 
    plot(gait_cycle,moment_complete,':r','LineWidth',2); xlim([gait_cycle(1) gait_cycle(end)]); 
    xlabel('gait cycle [%]');   ylabel('torque [Nm]'); title(['torque vs. gait cycle using ' num2str(extra_frame) ' extra frames'])
    % subplot(1,2,2); plot(gait_cycle,moment_complete); xlim([gait_cycle(1) gait_cycle(end)]); xlabel('gait cycle [%]');   ylabel('moment [Nm]'); title(['moment vs. gait cycle using ' num2str(extra_frame) ' extra frames'])
    % xline(0,':k');  xline(gait_cycle(ind_toe),'--r'); xline(100,':k');
end
end