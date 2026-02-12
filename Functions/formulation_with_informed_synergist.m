function [Results_new,Misc_new]=formulation_with_informed_synergist(computerPath,Misc,Results_noSyn,DatStore,DirF,computationCase,infoSyn)
OutLabel     = DirF.select_folder_N2;
SynergyPath  = fullfile(computerPath,Misc.TrialFolder,DirF.select_folder_N2);
synergy_list = 4:6; % Number of synergies to test

if strcmp(computationCase,'all')
    % synergy analysis
    [W,H,~,synMetrics]=synergyAnalysis(Results_noSyn,synergy_list,[0 0 0]);
    save(fullfile(SynergyPath,'synergy_metrics.mat'),'synMetrics');

    % compute for each synergy number
    for iSyn=1:length(synergy_list)
        sSyn=synergy_list(iSyn);
        Wsel=W{iSyn};
        Hsel=H{iSyn};

        % load to Misc
        Misc.SynCon.N=sSyn;
        Misc.SynCon.W=Wsel;
        Misc.SynCon.H=Hsel;
        
        % enable synergy control and save results
        Misc.Advance.SynergyControl=1;
        Misc.to_save_results= 1;

        % set path and run
        Misc.OutPath  = SynergyPath;
        Misc.OutName  = fullfile([OutLabel num2str(sSyn)]);
        % [Results_new,~,Misc_new]= MRS_Formulate_and_Solve_NeuroCons(Misc,DatStore);
        % [Results_new,~,Misc_new]= MRS_Formulate_and_Solve_NeuroConsV2(Misc,DatStore);
        [Results_new,~,Misc_new]= MRS_Formulate_and_Solve_NeuroConsV3(Misc,DatStore);
    end
elseif strcmp(computationCase,'load')
    sSyn=infoSyn.sSyn;
    if sSyn == 0
        Results_new=Results_noSyn;
        Misc_new   =Misc;
    elseif ~isempty((synergy_list==sSyn))
        MRS_new     = load(fullfile(SynergyPath, [OutLabel num2str(sSyn) 'Results.mat']));
        Results_new = MRS_new.Results;
        Misc_new    = MRS_new.Misc;
    end
end
end