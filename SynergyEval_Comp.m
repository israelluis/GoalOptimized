
mySub_list=[1:5];         nSubs=length(mySub_list);

iDev=1;

color_con={'' 'k' 'r'};
clc
for iSub=1:nSubs
    sSub=mySub_list(iSub);
    fig=figure(iSub); set(fig, 'Position', [1+360*(iSub-1), 50, 360, 900], 'Name', ['Subject' num2str(sSub)],'color','w'); clf;

    sTrial=1;
    iSynConf=1;
    MRS=MRS_list{sSub,iDev,iSynConf};

    line_con={'-' '-' '-.'};
    for iCon=2:3
        MRS_sel=MRS{iCon};
        Misc   =MRS_sel.Misc;
        Results=MRS_sel.Results(sTrial);
        N      =Misc.SynCon.N;

        color_list={'r' 'b' 'm'};
        % mus_list={'gasmed_r' 'soleus_r' 'tibant_r'};
        mus_list={'soleus' 'gasmed'};
        mus_list_label={'SO' 'GM'};
        nMusLis=length(mus_list);
        ind=zeros(nMusLis,1);
        for i=1:nMusLis
            ind(i)=find(strcmp(Results.MuscleNames,[mus_list{i} '_' Misc.gait_data(sTrial).side_sel]));
        end

        dataLength=length(Results.GaitCycle);
        % mus_sel=34;
        Hmus=zeros(nMusLis,N,dataLength);
        for iMus=1:nMusLis
            sMus=ind(iMus);
            for i=1:N

                H=Results.SynergyControl.H(i,1:end-1);
                W=Results.SynergyControl.W(sMus,i);
                Hmus(iMus,i,:)=W*H;

                subplot(N+1,3,1+(i-1)*3); hold on
                plot(Results.GaitCycle,H,'k','LineWidth',2,'LineStyle',line_con{iCon})
                if iMus==nMusLis
                    xlabel('gait cycle [%]'); ylabel('syn act [%]'); ylim([0 1]); title(['S#' num2str(i) ' TIME-SERIES'])
                end

                subplot(N+1,3,2+(i-1)*3); hold on
                bar(iMus,W, 'FaceColor', color_list{iMus})
                ylim([0 1])
                if iMus==nMusLis
                    set(gca, 'XTick', 1:nMusLis, 'XTickLabel', mus_list_label); title(['S#' num2str(i) ' WEIGHTS'])
                end

                subplot(N+1,3,3+(i-1)*3); hold on
                plot(squeeze(Hmus(iMus,i,:)),'Color',color_list{iMus},'LineWidth',2,'LineStyle',line_con{iCon})
                ylim([0 1])
                if iMus==nMusLis
                    xlabel('gait cycle [%]'); ylabel('mus act [%]'); ylim([0 1]); title(['S#' num2str(i) ' CONTRIBUTION'])
                end
            end
            subplot(N+1,3,(N+1)*3)
            hold on;
            plot(squeeze(sum(Hmus(iMus,:,:))),'Color',color_list{iMus},'LineWidth',2,'LineStyle',line_con{iCon})
            if iMus==nMusLis
                xlabel('gait cycle [%]'); ylabel('mus act [%]'); ylim([0 1]); title(['S#' num2str(i) ' TOTAL'])
            end

        end
    end
    sgtitle(['SUB' num2str(sSub)])
end
%% Visual inspection muscle activations and load results
mySub_list=[1 2 3 4 5];         nSubs=length(mySub_list);
myTrial_list=[1 2 3];           nTrials=length(myTrial_list);
iDev=1;

clc

fig=figure(10); set(fig, 'Position', [50, 50, 800, 900]); clf;
color_list={'r' 'g' 'c' 'b' 'k'};

clear Results_all
Results_all(nSubs*nTrials) = struct(); % This creates an array of nTrials empty structs
data_length_all =zeros(nSubs,nTrials);

c=0;
for iSub=1:nSubs
    sSub=mySub_list(iSub);

    iSynConf=1;
    conSyn  =1;

    MRS=MRS_list{sSub,iDev,iSynConf};

    MRS_sel=MRS{conSyn};
    Misc   =MRS_sel.Misc;
    Results=MRS_sel.Results;

    % sSyn=Misc.SynCon.N; % number of synergies
    % Wsel=Misc.SynCon.W; % nMuscles x N
    % Hsel=Misc.SynCon.H; % N x times_series (all trials e.g., [trial_1 trial_2 trial_3])

    c_d=0;
    for iTrial=1:nTrials
        sTrial=myTrial_list(iTrial);

        Results_trial=Results(sTrial);

        c=c+1;
        Results_all(c).MActivation=Results_trial.MActivation;
        Results_all(c).MuscleNames=Results_trial.MuscleNames;

        data_length=size(Results_trial.MActivation,2);
        data_length_all(iSub,iTrial)=data_length;

        GaitCycle=Results_trial.GaitCycle;
        for i=1:40
            subplot(5,8,i)
            activation=Results_trial.MActivation(i,1:end-1);

            plot(GaitCycle,activation,'Color',color_list{sSub}); hold on
            ylim([0 1])
            title(Results_trial.MuscleNames{i}(1:end-2))
        end
    end
    % data_length_all(iSub,[1 2])=
end
%% Compute synergies
synergy_list=[4 5 6];
[W,H,Htrial,~,synMetrics]=synergyAnalysis(Results_all,synergy_list,[0 0 0]);
%% Compute ranges
range_time=zeros(5,2);
acc=0;
for i=1:5
    total_sub_length=sum(data_length_all(i,:));
    range_time(i,:)=[1 total_sub_length]+acc;

    acc=acc+total_sub_length;
end
%%
% Sub
for iSyn=1:length(synergy_list)
    sSyn=synergy_list(iSyn);

    % load to SynPop
    for iSub=1:5
        SynPop(iSub).N=sSyn;
        SynPop(iSub).W(iSyn)=W(iSyn);
        SynPop(iSub).H(iSyn)={H{iSyn}(:,range_time(iSub,1):range_time(iSub,2))};
    end
end

save("SynergyPopulation.mat","SynPop")
%%
mySub_list=[1];         nSubs=length(mySub_list);

iDev=1;

color_con={'' 'k' 'r'};
clc
for iSub=1:nSubs
    sSub=mySub_list(iSub);
    fig=figure(iSub); set(fig, 'Position', [1+360*(iSub-1), 50, 360, 900], 'Name', ['Subject' num2str(sSub)],'color','w'); clf;

    sTrial=1;
    iSynConf=1;
    MRS=MRS_list{sSub,iDev,iSynConf};

    line_con={'-' '-' '-.'};
    for iCon=2:3
        MRS_sel=MRS{iCon};
        Misc   =MRS_sel.Misc;
        Results=MRS_sel.Results(sTrial);
        N      =Misc.SynCon.N;

        color_list={'r' 'b' 'm'};
        % mus_list={'gasmed_r' 'soleus_r' 'tibant_r'};
        mus_list={'soleus' 'gasmed'};
        mus_list_label={'SO' 'GM'};
        nMusLis=length(mus_list);
        ind=zeros(nMusLis,1);
        for i=1:nMusLis
            ind(i)=find(strcmp(Results.MuscleNames,[mus_list{i} '_' Misc.gait_data(sTrial).side_sel]));
        end

        dataLength=length(Results.GaitCycle);
        % mus_sel=34;
        % Hmus=zeros(nMusLis,N,dataLength);
        for iMus=1:nMusLis
            sMus=ind(iMus);
            for i=1:N

                % H=Results.SynergyControl.H(i,1:end-1);
                % W=Results.SynergyControl.W(sMus,i);

                H=SynPop(iSub).H(i,1:end-1);
                W=SynPop(iSub).W(sMus,i);

                Hmus(iMus,i,:)=W*H;

                subplot(N+1,3,1+(i-1)*3); hold on
                plot(H,'k','LineWidth',2,'LineStyle',line_con{iCon})
                if iMus==nMusLis
                    xlabel('gait cycle [%]'); ylabel('syn act [%]'); ylim([0 1]); title(['S#' num2str(i) ' TIME-SERIES'])
                end

                subplot(N+1,3,2+(i-1)*3); hold on
                bar(iMus,W, 'FaceColor', color_list{iMus})
                ylim([0 1])
                if iMus==nMusLis
                    set(gca, 'XTick', 1:nMusLis, 'XTickLabel', mus_list_label); title(['S#' num2str(i) ' WEIGHTS'])
                end

                subplot(N+1,3,3+(i-1)*3); hold on
                plot(squeeze(Hmus(iMus,i,:)),'Color',color_list{iMus},'LineWidth',2,'LineStyle',line_con{iCon})
                ylim([0 1])
                if iMus==nMusLis
                    xlabel('gait cycle [%]'); ylabel('mus act [%]'); ylim([0 1]); title(['S#' num2str(i) ' CONTRIBUTION'])
                end
            end
            subplot(N+1,3,(N+1)*3)
            hold on;
            plot(squeeze(sum(Hmus(iMus,:,:))),'Color',color_list{iMus},'LineWidth',2,'LineStyle',line_con{iCon})
            if iMus==nMusLis
                xlabel('gait cycle [%]'); ylabel('mus act [%]'); ylim([0 1]); title(['S#' num2str(i) ' TOTAL'])
            end

        end
    end
    sgtitle(['SUB' num2str(sSub)])
end
%%
for iSub=1:nSubs
    sSub=mySub_list(iSub);
    fig=figure(iSub); clf;

    iSynConf=1;
    MRS=MRS_list{sSub,iDev,iSynConf};
    for iTrial=1:nTrials
        for iCon=2:2
            MRS_sel=MRS{iCon};
            Results=MRS_sel.Results(iTrial);
            N=MRS_sel.Misc.SynCon.N;

            dataLength=length(Results.GaitCycle);
            mus_sel=34;
            Hmus=zeros(N,dataLength);
            HmusAcc=zeros(N,dataLength);
            for i=1:N
                H=Results.SynergyControl.H(i,1:end-1);
                W=Results.SynergyControl.W(:,i);
                Hmus(i,:)=W(mus_sel)*H;
                % HmusAcc(i,:)=sum(Hmus(1:i,:));
            end

            for iMus=mus_sel:mus_sel
                MActivation_JE=Results.MActivation(iMus,1:end-1);
                MActivation_JS=Results.SynergyControl.SynergyActivation(iMus,1:end-1);
                % subplot(5,8,iMus); hold on
                subplot(1,1,1); hold on
                plot(Results.GaitCycle,MActivation_JE,color_con{iCon},'LineWidth',3,'LineStyle','-')
                plot(Results.GaitCycle,MActivation_JS,color_con{iCon},'LineWidth',3,'LineStyle',':')

                for i=2:N %:N
                    plot(Results.GaitCycle,Hmus(i,:),'m','LineWidth',3,'LineStyle','-')
                end

                axis([0 100 0 1])
                title(Results.MuscleNames{iMus})
            end
        end
    end
end