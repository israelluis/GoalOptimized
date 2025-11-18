function [bounds,bounds_stackUp,nVars_dev]=getControlParamBounds_perType(Device)
nDevs=length(Device);
nVars_dev=zeros(nDevs,1);
% totalVars=0;
for iDev=1:nDevs
    jointSel=Device{iDev}.MuscleGroup{1}(1:end-2);
    jointDir=Device{iDev}.MuscleGroup{2};
    typeProfile=Device{iDev}.Type{2};

    if strcmp(typeProfile,'spline#N3')
            % ankle plantarflexion
        if strcmp(jointSel,'ankle_angle')       & (jointDir==-1)
            range=[30 55; 5 30; 5 30; 10 120];     % Range: 25 25 25 110   E= 55 20 5 100

            % knee extension
        elseif strcmp(jointSel,'knee_angle')    & (jointDir==-1)
            range=[5 30; 5 30; 5 30; 10 100];      % Range: 25 25 25 90    E= 15 10 25 70

            % hip flexion
        elseif strcmp(jointSel,'hip_flexion')   & (jointDir== 1)
            range=[45 75; 5 35; 5 35; 10 100];     % Range: 30 30 30 90    E= 65 32 35 40
            % rangeX=[40 80; 5 45; 5 45];   % Width timing changes results
            % rangeY=[10 100];

            % hip abduction
        elseif strcmp(jointSel,'hip_adduction') & (jointDir==-1)
            range=[40 70; 5 35; 5 35; 10 120];    % Range: 30 30 30 110    E= 50 15 10 100
        else
            error('the selected muscle group for assistance is not defined')
        end

    elseif strcmp(typeProfile,'clutchSpring')
            % ankle plantarflexion
        if strcmp(jointSel,'ankle_angle')       & (jointDir==-1)
            range=[10 30; 0.01 8];     % Range: 20 8   E= 21 7

            % knee extension
        elseif strcmp(jointSel,'knee_angle')    & (jointDir==-1)
            range=[0 10; 0.01 8];      % Range: 10 8    E= 6 4            

            % hip flexion
        elseif strcmp(jointSel,'hip_flexion')   & (jointDir== 1)
            range=[25 45; 0.01 5];     % Range: 20 4    E= 36 3

            % hip abduction
        elseif strcmp(jointSel,'hip_adduction') & (jointDir==-1)
            range=[0 20; 0.01 8];      % Range: 20 8    E= 3 3
        end
    end

   [varNames,varUnits]=getLabelAssistance(jointSel,jointDir,typeProfile);

   bounds{iDev}.range=range;
   bounds{iDev}.varNames=varNames;
   bounds{iDev}.varUnits=varUnits;

   nVars_dev(iDev)=length(varNames);
end

totalVars=sum(nVars_dev);
range_SU =zeros(totalVars,2);
vNames_SU=strings(totalVars,1);
vUnits_SU=strings(totalVars,1);
ind_acc=0;
for iDev=1:nDevs
    nVars=length(bounds{iDev}.varNames);
    ind_temp=1:nVars;
    
    range_SU(ind_temp+ind_acc,:) =bounds{iDev}.range;
    vNames_SU(ind_temp+ind_acc)=bounds{iDev}.varNames;
    vUnits_SU(ind_temp+ind_acc)=bounds{iDev}.varUnits;
    ind_acc=ind_acc+nVars;
end
bounds_stackUp.range=range_SU;
bounds_stackUp.varNames=vNames_SU;
bounds_stackUp.varUnits=vUnits_SU;
end