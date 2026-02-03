MRS_normal=load('C:\Users\Israel Luis\Documents\GitHub\GoalOptimized\ProjectResults\DSE\sub4\v2_t1\Je\JeResults.mat');
Results_normal=MRS_normal.Results;
%%
figure(100); clf;
for i =1:40
    subplot(5,8,i); hold on
    plot(Results_normal.MActivation.genericMRS(i,:),'-k','LineWidth',2)
    plot(Results.MActivation.genericMRS(i,:),':r','LineWidth',2)
    title(Results.MuscleNames{i})
    ylim([0 1])
end