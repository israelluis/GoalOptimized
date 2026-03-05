function [Results,DatStore,Misc] = MRS_Formulate_and_Solve_NM(Misc,DatStore)
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
e_min = 0.01; e_max = 1;                % bounds on muscle excitation
a_min = 0.01; a_max = 1;                % bounds on muscle activation
vMtilde_min = -10; vMtilde_max = 10;    % bounds on normalized muscle fiber velocity
lMtilde_min = 0.1; lMtilde_max = 1.7;   % bounds on normalized muscle fiber length

% CasADi setup
import casadi.*
opti    = casadi.Opti();    % create opti structure

nTrials = Misc.nTrials;
Mesh    =[DatStore.Mesh]';

% get total number of mesh points
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
%   - Muscle activation
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
%% Synergy
% Muscle synergy constraints
if Misc.Advance.SynergyControl
    nSynergies = Misc.SynCon.N; % Number of synergies

    % Synergy activation: time-variant
    H = opti.variable(nSynergies, N_tot+nTrials);
    opti.subject_to(0 < H < 1); % Bound synergy activations
    % opti.subject_to(0.01 < H <= 1); % Bound synergy activations
    
    % Synergy weight: time-invariant
    % W = opti.variable(NMuscles, nSynergies);
    % opti.subject_to(0 <= W <= 1);
    Wsel=Misc.SynCon.W;

    % Use unassisted synergy activation as initial guess
    if isfield(Misc.SynCon, 'H') && ~isempty(Misc.SynCon.H)
        opti.set_initial(H, Misc.SynCon.H);
    end
end
%%
if Misc.Advance.ActivationLowerBound  % You'll need to add this field to Misc
    alpha = Misc.Advance.alpha;  % e.g., 0.1 or 0.05
    
    % You need the unassisted activations - load them from your baseline simulation
    % This assumes you have stored them in Misc.UnassistedActivations
    if isfield(Misc, 'UnassistedActivations')
        a_unassisted = Misc.UnassistedActivations;  % size: NMuscles × N_tot+nTrials
    else
        warning('Unassisted activations not found. Activation bound not applied.');
    end
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

    % Tdev(Misc.nTrials, nDevices) = struct('T', [], 'Ind', []);
    for trial = 1:Misc.nTrials
        Device=Misc.Device(trial);
        nDevices=length(Device);

        % initialize variables
        % idealT = opti.variable(nDevices,N_tot);
        N_trial=Mesh(trial).N;
        idealT = MX(nDevices,N_trial);

        % loop for each exoskeleton
        for iDev=1:nDevices

            Exo_Group =Device{iDev}.MuscleGroup{1};

            % % get kinematics and time
            iSel             = find(ismember(Misc.DofNames_Input{trial},[Exo_Group '_' Misc.gait_data(trial).side_sel]));
            MuscleGroup_sign = Device{iDev}.MuscleGroup{2};
            torque           = Misc.Device{iDev}.Assistance.Profile.Torque;
            
            %interpolate to mesh points
            idx            = 1:length(torque);
            idxq           = linspace(min(idx), max(idx), N_trial);
            Vi             = interp1(idx, torque, idxq, 'pchip');
            idealT(iDev,:) =Vi*MuscleGroup_sign;

            Tdev(trial,iDev).T   = idealT(iDev,:);
            Tdev(trial,iDev).Ind = iSel;            % index of applied torque

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
        [Hilldiffk,FTk,~, ~, ~, cos_alpha,fce] = ForceEquilibrium_lMtildeState_optPassive(ak,lMtildek,vMtildek,lM_projectedk,...
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
                Device=Misc.Device(trial);
                nDevices=length(Device);
                for iDev= 1:nDevices
                    if dof == Tdev(trial,iDev).Ind
                        T_exp = T_exp - Tdev(trial,iDev).T(k);
                    end
                end
            end
            opti.subject_to(T_exp - T_sim == 0);

            if Misc.Advance.ActivationLowerBound  % You'll need to add this field to Misc
                % Add constraint at each mesh point
                opti.subject_to(alpha * a_unassisted(:,(N_acc+trial-1) + k) < a(:,(N_acc+trial-1) + k) ); %340
            end

            applyMinCoContraction=0;

            if applyMinCoContraction==1
                if dof==1
                % min active co-contraction torque level
                min_ank=5;

                FMo=MuscProperties_params(1,:);
                ind_dorsi  =find(DatStore(trial).MAinterp(k,index_sel)>0);
                ind_plantar=find(DatStore(trial).MAinterp(k,index_sel)<0);

                Fce_proj_dorsi  =(fce(ind_dorsi)'.*FMo(ind_dorsi).*cos_alpha(ind_dorsi)')';
                Fce_proj_plantar=(fce(ind_plantar)'.*FMo(ind_plantar).*cos_alpha(ind_plantar)')';

                T_sim_dorsi   = DatStore(trial).MAinterp(k,ind_dorsi)*Fce_proj_dorsi; %                
                T_sim_plantar = DatStore(trial).MAinterp(k,ind_plantar)*Fce_proj_plantar; %

                T_sim_dorsi_net(k)=T_sim_dorsi;

                opti.subject_to(T_sim_dorsi >= min_ank);
                opti.subject_to(T_sim_plantar <= -min_ank);
            end
            end
        end
        % Hill-equilibrium constraint
        opti.subject_to(Hilldiffk == 0);
    end
    N_acc = N_acc + N;
    % Cost function
    J = J + ...
        Misc.wAct*0.5*(sum(e(:).^2)/N/NMuscles + sum(a(:).^2)/N/NMuscles) + ... % Misc.wAct*0.5*(sumsqr(e)/N/NMuscles + sumsqr(a)/N/NMuscles) +
        Misc.wTres*sumsqr(aT)/N/DatStore(trial).nDOF + ...
        Misc.wVm*sumsqr(vMtilde)/N/NMuscles; % this is faster than sumsqr(e) and sumsqr(a) has excitations fast spikes

    if Misc.Advance.SynergyControl

        % FORMULATION 1
        % This formulation DOES NOT solve the issue
        % Misc.Advance.SynergyWeight=1;
        % synergy_error = a - Wsel * H;
        % J = J + Misc.Advance.SynergyWeight * sumsqr(synergy_error)/N/nSynergies;
        
        % FORMULATION 2
        % This formulation is the best based on objective function
        % penalization
        Misc.wSc=1.00; % 0.10 yet it has regularization issues (30 secs)
        % Misc.wSc=10.00; % 0.10 yet it has regularization issues (30 secs)

        % wSc=0.01 does not obey synergy pattern (quick)
        % wSc=1.00 poor estimation, regularization issues, glmed reach max, met cost red small (4%) (60 secs)
        
        eps_syn=0.01; % good value based on T&E
        % eps_syn=0.001 make the optimization sensitive to small values 
        % eps_syn=0.1 delutes the relative importance formulation
        
        diff = (a - Wsel * H) ./ (a + eps_syn);
        % diff = (a - Wsel * H);
        J = J + Misc.wSc.*sumsqr(diff)/N/NMuscles;
    end

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

%% Store Results
Results(nTrials) = struct(); % This creates an array of nTrials empty structs

% store parameters
for trial = 1:Misc.nTrials
    Results(trial).Params=Misc.params; % MuscProperties_params_opt
    Results(trial).kT    =Misc.kT;     % MuscProperties_kT_opt
end

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
        Results(trial).Device=Misc.Device;
    end
end

if Misc.Advance.SynergyControl
    H_opt=sol.value(H);

    Ntot = 0;
    for trial = 1:Misc.nTrials
        t0 = DatStore(trial).time(1); tf = DatStore(trial).time(end);
        N = round((tf-t0)*Misc.Mesh_Frequency);

        % Load synergy features
        Results(trial).SynergyControl.H = H_opt(:,(Ntot + trial - 1) + 1:(Ntot + trial - 1) + N + 1);
        Results(trial).SynergyControl.W = Wsel;
        Results(trial).SynergyControl.SynergyActivation = Wsel *  Results(trial).SynergyControl.H;

        % Calculate reconstruction error
        synergy_reconstruction_error = a_opt(:,(Ntot + trial - 1) + 1:(Ntot + trial - 1) + N + 1) - Results(trial).SynergyControl.SynergyActivation;
        Results(trial).SynergyControl.RMSE = sqrt(mean(synergy_reconstruction_error.^2, 'all'));
        Results(trial).SynergyControl.VAF = 1 - (sum(synergy_reconstruction_error.^2, 'all') / sum(a_opt.^2, 'all'));

        Ntot = Ntot + N;
    end
end

% Append results to output structures
Ntot = 0;

for trial = 1:nTrials
    t0 = DatStore(trial).time(1); tf = DatStore(trial).time(end);
    N = round((tf-t0)*Misc.Mesh_Frequency);
    % Time grid
    tgrid = linspace(t0,tf,N+1)';

    % Gait Cycle grid
    initial=0;   final=100;

    data_GC      = N-2*Misc.extra_frames;
    frames_per_GC= final/(data_GC-1);
    extra_times  = frames_per_GC*Misc.extra_frames;
    gait_cycle   = initial-extra_times:frames_per_GC:final+extra_times; % case without extra frame-> gait_cycle =linspace(0,100,length(time_series));

    % Save results
    Results(trial).Time     = tgrid;
    Results(trial).GaitCycle = gait_cycle;
    Results(trial).MActivation = a_opt(:,(Ntot + trial - 1) + 1:(Ntot + trial - 1) + N + 1);
    Results(trial).lMtildeopt = lMtilde_opt(:,(Ntot + trial - 1) + 1:(Ntot + trial - 1) + N + 1);
    Results(trial).lM = lMtilde_opt(:,(Ntot + trial - 1) + 1:(Ntot + trial - 1) + N + 1).*repmat(Misc.lMo',1,length(tgrid));
    Results(trial).vMtilde = vMtilde_opt(:,Ntot + 1:Ntot + N);
    Results(trial).lM_projected_opt = lM_projected_opt(:,Ntot + 1:Ntot + N);
    Results(trial).MExcitation = e_opt(:,Ntot + 1:Ntot + N);
    Results(trial).RTorque     = aT_opt(:,Ntot + 1:Ntot + N)*Misc.Topt;
    Results(trial).MuscleNames = DatStore(trial).MuscleNames;
    Results(trial).OptInfo     = output;
    % Tendon force
    Results(trial).lMTinterp   = DatStore(trial).LMTinterp';
    [TForcetilde_,TForce_,lTtilde_] = TendonForce_lMtilde(Results(trial).lMtildeopt',MuscProperties_params_opt,Results(trial).lMTinterp',MuscProperties_kT_opt,MuscProperties_shift_opt);
    Results(trial).TForcetilde = TForcetilde_';
    Results(trial).TForce      = TForce_';
    Results(trial).lTtilde     = lTtilde_';
    % get information F/l and F/v properties - updated with passive forces
    [Fpe_,FMltilde_,FMvtilde_] = getForceLengthVelocityProperties_setPassiveParam(Results(trial).lMtildeopt',Results(trial).vMtilde',MuscProperties_params_opt(5,:),...
                                                                                  MuscProperties_params_opt(6,:),MuscProperties_params_opt(7,:),MuscProperties_params_opt(8,:));
    FMo = ones(N+1,1)*Misc.params(1,:);
    Results(trial).Fpe        = (Fpe_.*FMo)';
    Results(trial).FMltilde   = FMltilde_';
    Results(trial).FMvtilde   = FMvtilde_';
    % Results(trial).Fce = Misc.params(1,:).*Results(trial).MActivation'.*FMltilde_.*FMvtilde_;
    Ntot = Ntot + N;
end

%% save the results
% plot states and variables from parameter estimation simulation
if Misc.to_save_results==1
    save(fullfile(Misc.OutPath,[Misc.OutName 'Results.mat']),'Results','DatStore','Misc');
end
end
