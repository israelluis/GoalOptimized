subject_mass=Misc.subject_data.subject_mass;    

Device{1}.Mode        = 'prescribed'; % opts: optimized and prescribed
Device{1}.MuscleGroup = {['ankle_angle_' Misc.gait_data.side_sel] -1};
Device{1}.Type        = {'active' 'spline#N3'}; % opts: active, quasi-passive, passive, EMG-driven
Device{1}.Params      = [54.9 14.2 8.0 55.7]; %54.8 14.5 8.1 51.8
[assistanceInfo]=generateTorque(Device{1},DatStore_normal,Misc.time,Misc.extra_frames);

GaitCycle=assistanceInfo.Profile.GaitCycle;
TorqueSyn=assistanceInfo.Profile.Torque/subject_mass;

Device{1}.Mode        = 'prescribed'; % opts: optimized and prescribed
Device{1}.MuscleGroup = {['ankle_angle_' Misc.gait_data.side_sel] -1};
Device{1}.Type        = {'active' 'spline#N3'}; % opts: active, quasi-passive, passive, EMG-driven
Device{1}.Params      = [54 27 8.6 100]; %54.8 14.5 8.1 51.8
[assistanceInfo]=generateTorque(Device{1},DatStore_normal,Misc.time,Misc.extra_frames);

GaitCycle=assistanceInfo.Profile.GaitCycle;
TorqueHILO=assistanceInfo.Profile.Torque/100*0.68;


Device{1}.Mode        = 'prescribed'; % opts: optimized and prescribed
Device{1}.MuscleGroup = {['ankle_angle_' Misc.gait_data.side_sel] -1};
Device{1}.Type        = {'active' 'spline#N3'}; % opts: active, quasi-passive, passive, EMG-driven
Device{1}.Params      = [55 24 8.6 90]; %54.8 14.5 8.1 51.8
[assistanceInfo]=generateTorque(Device{1},DatStore_normal,Misc.time,Misc.extra_frames);

GaitCycle=assistanceInfo.Profile.GaitCycle;
TorqueNo=assistanceInfo.Profile.Torque/subject_mass;

clf;
hold on
plot(GaitCycle,TorqueNo,':k','LineWidth',5);
plot(GaitCycle,TorqueSyn,'r','LineWidth',5);
plot(GaitCycle,TorqueHILO,'-.b','LineWidth',5);
legend('No synergies','Synergy informed','HILO')

legend boxoff
ylim([0 1.8])
xlim([0 100])

set(gca,'FontSize',15);

xlabel('gait cycle [%]')
ylabel('metabolic cost reduction [%]')