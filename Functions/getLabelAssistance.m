function [varNames,varUnits]=getLabelAssistance(jointSel,jointDir,assistanceType)
if strcmp(assistanceType,'spline#N3')
    varNames={'tp' 'tr' 'tf' 'Mp'};          
    varUnits={'%GC' '%GC' '%GC' 'Nm'};
elseif strcmp(assistanceType,'clutchSpring')
    varNames={'tc' 'K'};          
    varUnits={'%GC' 'Nm/deg'};
end

if strcmp(jointSel,'ankle_angle')       && jointDir==-1
    MG='AP';
elseif strcmp(jointSel,'knee_angle')    && jointDir==-1
    MG='KE';
elseif strcmp(jointSel,'hip_flexion')   && jointDir== 1
    MG='HF';
elseif strcmp(jointSel,'hip_adduction') && jointDir==-1
    MG='HB';
end

varNames = strcat(MG, '_', varNames);
end