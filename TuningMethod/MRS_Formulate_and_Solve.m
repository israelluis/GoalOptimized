function [Results,DatStore,Misc] = MRS_Formulate_and_Solve(Misc,DatStore)
Mesh=DatStore.Mesh;

%% setup options for the solver
% Create an NLP solver
% output.setup.lM_projecteddata = lM_projecteddata;
output.setup.nlp.solver = 'ipopt';
output.setup.nlp.ipoptoptions.linear_solver = 'mumps';
% Set derivativelevel to 'first' for approximating the Hessian
output.setup.derivatives.derivativelevel = 'second';
output.setup.nlp.ipoptoptions.tolerance = 1e-6;
output.setup.nlp.ipoptoptions.maxiterations = 1000;
if strcmp(output.setup.derivatives.derivativelevel, 'first')
    optionssol.ipopt.hessian_approximation = 'limited-memory';
end
% By default, the barrier parameter update strategy is monotone.
% https://www.coin-or.org/Ipopt/documentation/node46.html#SECTION000116020000000000000
% Uncomment the following line to use an adaptive strategy
% optionssol.ipopt.mu_strategy = 'adaptive'; 
optionssol.ipopt.nlp_scaling_method = 'gradient-based';
optionssol.ipopt.linear_solver = output.setup.nlp.ipoptoptions.linear_solver;
optionssol.ipopt.tol = output.setup.nlp.ipoptoptions.tolerance;
optionssol.ipopt.max_iter = output.setup.nlp.ipoptoptions.maxiterations;

%% Dynamic Optimization - Default parameters
% ----------------------------------------------------------------------- %
% Solve muscle redundancy problem with default parameters
% Problem bounds
e_min = 0; e_max = 1;                   % bounds on muscle excitation
a_min = 0; a_max = 1;                   % bounds on muscle activation
vMtilde_min = -10; vMtilde_max = 10;    % bounds on normalized muscle fiber velocity
lMtilde_min = 0.1; lMtilde_max = 1.7;   % bounds on normalized muscle fiber length

% CasADi setup
import casadi.*
opti    = casadi.Opti();    % create opti structure

% get total number of mesh points
nTrials = Misc.nTrials;
N_tot = sum([Mesh().N]);

% get the number of muscles
NMuscles = length(DatStore(1).MuscleNames);

% set intial guess based on static opt data
SoActGuess = zeros(NMuscles,N_tot);
SoExcGuess = zeros(NMuscles,N_tot-nTrials);
lMtildeGuess = zeros(NMuscles,N_tot);
vMtildeGuess = zeros(NMuscles,N_tot-nTrials);
SoRActGuess = zeros(DatStore(1).nDOF,N_tot-nTrials);
ctx = 1;  ctu= 1;

if strcmp(Misc.initial_guess_source,'noExo')
    dir_opt= [Misc.OutFolder '\' Misc.OutFile(1:6) '_exoNN_simCAL_funACT_ORun\Results']; %exoNN_simCAL_funACT
    store_name='simulation_GenCalValResults';
    sim=load(fullfile(dir_opt,store_name));
for trial = 1:nTrials
    ctx_e = ctx+Mesh(trial).N;      % counter for states
    ctu_e = ctu+Mesh(trial).N-1;    % counter for controls
    SoActGuess(:,ctx:ctx_e) = sim.Results.MActivation.genericMRS;
    SoExcGuess(:,ctu:ctu_e) = sim.Results.MActivation.genericMRS(:,1:end-1);
    lMtildeGuess(:,ctx:ctx_e) = sim.Results.lMtildeopt.genericMRS;
    vMtildeGuess(:,ctu:ctu_e) = sim.Results.vMtilde.genericMRS;
    SoRActGuess(:,ctu:ctu_e) = sim.Results.RActivation.genericMRS; %RA torque, not activation
    ctx = ctx_e+1;
    ctu = ctu_e+1;
end    
elseif strcmp(Misc.initial_guess_source,'SO')
for trial = 1:nTrials
    ctx_e = ctx+Mesh(trial).N;      % counter for states
    ctu_e = ctu+Mesh(trial).N-1;    % counter for controls
    SoActGuess(:,ctx:ctx_e) = DatStore(trial).SoActInterp';
    SoExcGuess(:,ctu:ctu_e) = DatStore(trial).SoActInterp(1:end-1,:)';
    lMtildeGuess(:,ctx:ctx_e) = DatStore(trial).lMtildeInterp';
    vMtildeGuess(:,ctu:ctu_e) = DatStore(trial).vMtildeinterp(1:end-1,:)';
    SoRActGuess(:,ctu:ctu_e) = DatStore(trial).SoRActInterp(1:end-1,:)';
    ctx = ctx_e+1;
    ctu = ctu_e+1;
end
end

% States
%   - muscle activation
a = opti.variable(NMuscles,N_tot+nTrials);      % Variable at mesh points
opti.subject_to(a_min < a < a_max);             % Bounds
opti.set_initial(a,SoActGuess);                 % Initial guess (static optimization)
%   - Muscle fiber lengths
lMtilde = opti.variable(NMuscles,N_tot+nTrials);
opti.subject_to(lMtilde_min < lMtilde < lMtilde_max);
opti.set_initial(lMtilde,lMtildeGuess);
%   - Controls
e = opti.variable(NMuscles,N_tot);
opti.subject_to(e_min < e < e_max);
opti.set_initial(e, SoExcGuess);
%   - Reserve actuators
aT = opti.variable(DatStore(trial).nDOF,N_tot);
opti.subject_to(-1 < aT <1);
% opti.set_initial(aT, SoRActGuess/Misc.Topt); %Not added. Longer to converge
%   - Time derivative of muscle-tendon forces (states)
vMtilde = opti.variable(NMuscles,N_tot);
opti.subject_to(vMtilde_min < vMtilde < vMtilde_max);
opti.set_initial(vMtilde,vMtildeGuess);
%   - Projected muscle fiber length - Auxilary variable to avoid muscle buckling & square root expression in muscle dynamics
lM_projected = opti.variable(NMuscles,N_tot + nTrials);
opti.subject_to(1e-4 < lM_projected(:)); % We impose that projected muscle fiber length has strict positive length

%% Declare and setup muscle fiber optimization
J=0;
if Misc.Advance.TuningMethod_fiber==1
    [opti,J,kT,shift,lMo,lTs] = formulateMuscleFiberCalibration(opti,J,Misc,lMtilde,nUSdata,ind_US,ind_USnone,USDigitalizedInterp,N,NMuscles,DatStore.MuscleNames);
else
    [opti,J,kT,shift,lMo,lTs] = formulateMuscleFiberCalibration(opti,J,Misc,lMtilde,[],[],[],[],[],[],[]);
end

%% Initial guess for this variable is retrieved from lMtilde guess
% and geometric relationship between pennation angle, muscle length
% and width
lMo_default = Misc.params(2,:)'; 
alphao = Misc.params(4,:)';
lMGuess = lMtildeGuess.*lMo_default;
w = lMo_default.*sin(alphao);
lM_projectedGuess = sqrt((lMGuess.^2 - w.^2));
opti.set_initial(lM_projected,lM_projectedGuess);

% constraint on projected fiber length
w = lMo.*sin(alphao);
lM = lMtilde.*lMo;
opti.subject_to(lM.^2 - w.^2 == lM_projected.^2);

% output optimization variables
MuscProperties_params =[Misc.params(1,:)' lMo lTs Misc.params(4:8,:)']';
MuscProperties_kT     =kT';
MuscProperties_shift  =shift';
%% Exoskeleton modelling
% Optimized in the same loop as in the muscle redundancy problem
if Misc.Advance.AssistiveDevice
    Device=Misc.Device;
    nDevices=length(Device);

    % Texo(Misc.nTrials,nExos);
    Texo(Misc.nTrials, nDevices) = struct('T', [], 'Ind', []);
    for trial = 1:Misc.nTrials
        
        % initialize variables
        idealT = opti.variable(nDevices,N_tot);

        % loop for each exoskeleton
        for iExo=1:nDevices
            
            Exo_Type  =Device{iExo}.Type{1};
            Exo_Group =Device{iExo}.MuscleGroup{1};

            % % get kinematics and time
            iSel             = find(ismember(Misc.DofNames_Input,Exo_Group));
            MuscleGroup_sign = Device{iExo}.MuscleGroup{2};
            torque           = Misc.Device{iExo}.Assistance.Profile.Torque;
            
            %interpolate to mesh points
            idx            = 1:length(torque);
            idxq           = linspace(min(idx), max(idx), N_tot);
            Vi             = interp1(idx, torque, idxq, 'pchip');
            idealT(iExo,:)=Vi*MuscleGroup_sign;

            Texo(trial,iExo).T   = idealT(iExo,:);
            Texo(trial,iExo).Ind = iSel;      % index of applied torque

        end
    end
end
%% Implemetation of controls, states and states derivatives
N_acc = 0; % Index that keeps track of trials that are accumulated
% Loop over trials --> one simulation for each trial
for trial = 1:Misc.nTrials
    % Time bounds
    t0 = DatStore(trial).time(1); tf = DatStore(trial).time(end);
    % Discretization
    N = Mesh(trial).N;
    h = Mesh(trial).step;
    
    % Loop over mesh points formulating NLP
    for k=1:N
        % Variables within current mesh interval
        ak = a(:,(N_acc+trial-1) + k); lMtildek = lMtilde(:,(N_acc+trial-1) + k);
        vMtildek = vMtilde(:,N_acc + k); aTk = aT(:,N_acc + k); ek = e(:,N_acc + k);
        lM_projectedk = lM_projected(:,(N_acc+trial-1) + k);
        
        % Euler integration  Uk = (X_(k+1) - X_k)/*dt
        Xk = [ak; lMtildek];
        Zk = [a(:,(N_acc+trial-1) + k + 1);lMtilde(:,(N_acc+trial-1) + k + 1)];
        Uk = [ActivationDynamics(ek,ak,Misc.tauAct,Misc.tauDeact,Misc.b); vMtildek];
        opti.subject_to(eulerIntegrator(Xk,Zk,Uk,h) == 0);
        
        % Get muscle-tendon forces and derive Hill-equilibrium
        [Hilldiffk,FTk] = ForceEquilibrium_lMtildeState_optPassive(ak,lMtildek,vMtildek,lM_projectedk,...
            DatStore(trial).LMTinterp(k,:)',MuscProperties_params',MuscProperties_kT',MuscProperties_shift'); 
        
        % Add path constraints
        % Moment constraints
        for dof = 1:DatStore(trial).nDOF
            T_exp = DatStore(trial).IDinterp(k,dof);
            index_sel = (dof-1)*(NMuscles)+1:(dof*NMuscles); % moment is a vector with the different dofs "below" each other
            T_sim = DatStore(trial).MAinterp(k,index_sel)*FTk + Misc.Topt*aTk(dof);
            
            % subtract exoskeleton moment from ID torque
            % (TID = Tmuscles + Texo)
            if Misc.Advance.AssistiveDevice
                for iExo= 1:nDevices
                    if dof == Texo(trial,iExo).Ind
                        T_exp = T_exp - Texo(trial,iExo).T(k);
                    end
                end
            end
            opti.subject_to(T_exp - T_sim == 0);
        end
        % Hill-equilibrium constraint
        opti.subject_to(Hilldiffk == 0);
    end
    N_acc = N_acc + N;
% Cost function
J = J + ... 
    Misc.wAct*0.5*(sumsqr(e)/N/NMuscles + sumsqr(a)/N/NMuscles) + ...
    Misc.wTres*sumsqr(aT)/N/DatStore(trial).nDOF + ...
    Misc.wVm*sumsqr(vMtilde)/N/NMuscles; % this is faster than sumsqr(e) and sumsqr(a) has excitations fast spikes
end

opti.minimize(J); % Define cost function in opti

% Create an NLP solver
opti.solver(output.setup.nlp.solver,optionssol);

% Solve
diary(fullfile(Misc.OutPath,[Misc.OutName 'GenericMRS.txt']));
tic
sol = opti.solve();
dt = toc;
disp(['Computation time solving OCP: ' num2str(dt) ' s'])
diary off

% Extract results
% Variables at mesh points
% Muscle activations and muscle-tendon forces
a_opt = sol.value(a);
lMtilde_opt = sol.value(lMtilde);
% Muscle excitations
e_opt = sol.value(e);
% Reserve actuators
aT_opt = sol.value(aT);
% Time derivatives of muscle-tendon forces
vMtilde_opt = sol.value(vMtilde);
% Optimal lM_projectedilary variable
lM_projected_opt = sol.value(lM_projected);

% get parameters
MuscProperties_params_opt=sol.value(MuscProperties_params);
MuscProperties_kT_opt    =sol.value(MuscProperties_kT);
MuscProperties_shift_opt =sol.value(MuscProperties_shift);

% store parameters
Results.params.calibratedMRS=MuscProperties_params_opt;
Results.kT.calibratedMRS    =MuscProperties_kT_opt;
Results.params.genericMRS   =Misc.params;
Results.kT.genericMRS       =Misc.kT;

% compute unNormalized tendon values
if Misc.Advance.TuningMethod_fiber == 1
kT_sel =MuscProperties_kT_opt(ind_US_AT);
FMo_sel=MuscProperties_params_opt(1,ind_US_AT);
lTs_sel=MuscProperties_params_opt(3,ind_US_AT);
[kT_abs_computed] = unNormalizedTendon(kT_sel,FMo_sel,lTs_sel);
Results.kT_ATunNormalized.calibratedMRS=sol.value(kT_abs_computed);

kT_sel =Misc.kT(ind_US_AT);
FMo_sel=Misc.params(1,ind_US_AT);
lTs_sel=Misc.params(3,ind_US_AT);
[kT_abs_computed] = unNormalizedTendon(kT_sel,FMo_sel,lTs_sel);
Results.kT_ATunNormalized.genericMRS   =sol.value(kT_abs_computed);
% are they differents? quick to check -> compare_param=MuscProperties_params_opt-Misc.params;
end

% Append results to output structures of exoskeleton
if Misc.Advance.AssistiveDevice
    for trial = 1:Misc.nTrials
        % for iExo= 1:nDevices
            Results.Device=Misc.Device;
        % end
    end
end

% Append results to output structures
Ntot = 0;
for trial = 1:nTrials
    t0 = DatStore(trial).time(1); tf = DatStore(trial).time(end);
    N = round((tf-t0)*Misc.Mesh_Frequency);
    % Time grid
    tgrid = linspace(t0,tf,N+1)';
    % Save results
    Results.Time(trial).genericMRS = tgrid;
    Results.MActivation(trial).genericMRS = a_opt(:,(Ntot + trial - 1) + 1:(Ntot + trial - 1) + N + 1);
    Results.lMtildeopt(trial).genericMRS = lMtilde_opt(:,(Ntot + trial - 1) + 1:(Ntot + trial - 1) + N + 1);
    Results.lM(trial).genericMRS = lMtilde_opt(:,(Ntot + trial - 1) + 1:(Ntot + trial - 1) + N + 1).*repmat(Misc.lMo',1,length(tgrid));
    Results.vMtilde(trial).genericMRS = vMtilde_opt(:,Ntot + 1:Ntot + N);
    Results.lM_projected_opt(trial).genericMRS = lM_projected_opt(:,Ntot + 1:Ntot + N);
    Results.MExcitation(trial).genericMRS = e_opt(:,Ntot + 1:Ntot + N);
    Results.RActivation(trial).genericMRS = aT_opt(:,Ntot + 1:Ntot + N)*Misc.Topt;
    Results.MuscleNames = DatStore.MuscleNames;
    Results.OptInfo = output;
    % Tendon force
    Results.lMTinterp(trial).genericMRS = DatStore(trial).LMTinterp';
    [TForcetilde_,TForce_] = TendonForce_lMtilde(Results.lMtildeopt(trial).genericMRS',MuscProperties_params_opt,Results.lMTinterp(trial).genericMRS',MuscProperties_kT_opt,MuscProperties_shift_opt);
    Results.TForcetilde(trial).genericMRS = TForcetilde_';
    Results.TForce(trial).genericMRS = TForce_';
    % get information F/l and F/v properties - updated with passive forces
    [Fpe_,FMltilde_,FMvtilde_] = getForceLengthVelocityProperties_setPassiveParam(Results.lMtildeopt(trial).genericMRS',Results.vMtilde(trial).genericMRS',MuscProperties_params_opt(5,:),...
                                                                                  MuscProperties_params_opt(6,:),MuscProperties_params_opt(7,:),MuscProperties_params_opt(8,:));
    FMo = ones(N+1,1)*Misc.params(1,:);
    Results.Fpe(trial).genericMRS = (Fpe_.*FMo)';
    Results.FMltilde(trial).genericMRS = FMltilde_';
    Results.FMvtilde(trial).genericMRS = FMvtilde_';
    Ntot = Ntot + N;
end

%% Store Results

% store the Misc structure as well in the results
Results.Misc = Misc;

% add selected muscle names to the output structure
Results.MuscleNames = DatStore.MuscleNames;

%% save the results
% plot states and variables from parameter estimation simulation
if Misc.to_save_results==1
    save(fullfile(Misc.OutPath,[Misc.OutName 'Results.mat']),'Results','DatStore','Misc');
end
end
