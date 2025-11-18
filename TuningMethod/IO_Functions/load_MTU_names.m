function [param_label] = load_MTU_names()
% Consistent with definition in Muscle Redundandy Solver (MRS) and tuning
% method paper i.e., Luis et al. 2024 Sci. Rep. 

% AS IN MRS
% FMo,    maximum isometric force
% lMo,    optimal fiber length
% lTs,    tendon slack length
% alphao, pennation angle
% vMmax,  maximum contraction velocity

% AS IN PAPER
% kpe,    exponential shape factor for the passive force–length relationship
% so,     normalized fiber length at which the passive force starts to increase
% sM,     normalized fiber length, measured from the optimal fiber length at which maximum force is reached

param_label={'FMo' 'lMo' 'lTs' 'alphao' 'vMmax' 'kpe' 'so' 'sM'};
end