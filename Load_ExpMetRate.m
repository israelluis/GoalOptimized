%% PART I - ENERGETICS
computerPath  ='C:\Users\Israel Luis\Documents\GitHub';
folderPath      =fullfile(computerPath,'GoalOptimized\Database\EE');

sub_list=1:5;
sub_mass=[72.6 63.8 59.2 85.2 64];
readData      = 1;
computeAvgVal = 1;    plotAvgVal = 0;
computeEE     = 1;    plotEE     = 0;

Einfo_data=cell.empty;
Einfo_meta=cell.empty;
Einfo_speedKM_H=double.empty;
EE_rest_kg=double.empty;
EE_rest_kcalMin=double.empty;
Einfo_speedM_S =double.empty;
for sub_opt=1:length(sub_list)
    sub_sel=sub_list(sub_opt);
    
    dir_EE=[folderPath '\SUB' num2str(sub_sel)];
    Einfo_data{sub_sel}= readtable(fullfile(dir_EE,['SUB' num2str(sub_sel) '_CompleteTrials_sumary_COT.xlsx']));
    Einfo_meta{sub_sel}= readtable(fullfile(dir_EE,['SUB' num2str(sub_sel) '_meta.xlsx']));
    
    EE_rest=[115.3076,... %1  %averaged SEL ()  & S(115.3076)
             153.6829,... %2  %averaged SEL (3) & S(134.3177)
             162.4584,... %3  %averaged SEL (3) & S() -
             159.5393,... %4  %averaged SEL (3) & S(167.5838)
             117.7489,... %5  %averaged SEL (3) & S(119.9251)
             166.7143,... %6  %averaged SEL (3) & S(134.1043)
             124.2981,... %7  %averaged SEL (3) & S(101.1094)
             203.9296,... %8  %averaged SEL (3) & S(254.5499) -
             132.8599,... %9  %averaged SEL (3) & S(127.0893)
             103.6485,... %10 %averaged SEL (3) & S(114.2198)
             120.4499,... %11 %averaged SEL (3) & S(124.4214)
             112.9731,... %12 %averaged SEL (3) & S(132.3272)
             132.5861,... %13 %averaged SEL (3) & S(117.6208) -
             113.9725,... %14 %averaged SEL (3) & S(101.2821)
              63.8992,... %15 %averaged SEL (3) & S( 72.3247)
              78.6720,... %16 %averaged SEL (3) & S( 89.8115)
             ]; % J/s
    EE_rest_kg(sub_sel)     =EE_rest(sub_sel)/sub_mass(sub_sel); % J/(s*kg)
    EE_rest_kcalMin(sub_sel)=EE_rest(sub_sel)*(60/4184); % kcal/(min)
    Einfo_speedKM_H(:,sub_sel)  =Einfo_meta{sub_sel}.KM_h;     
    Einfo_speedM_S(:,sub_sel)   =Einfo_meta{sub_sel}.KM_h*1000/3600; % Km/h -> m/s
end
numTrials=8;

%% get variables from tables
if readData==1
    VO2_dot=cell.empty; RER    =cell.empty; VCO2_dot=cell.empty;    time    =cell.empty;    last    =cell.empty;
    for sub_opt=1:length(sub_list)
        sub_sel=sub_list(sub_opt);
        for trial_opt=1:numTrials % for each trial (speed)
            
            lastidx = find(isnan(table2array(Einfo_data{sub_sel}(:,-1+3*trial_opt))),1) - 1;
            if isempty(lastidx)
                lastidx = length(table2array(Einfo_data{sub_sel}(:,-1+3*trial_opt))); 
            end
            
            clear duration
            out = datevec(datetime(table2array(Einfo_data{sub_sel}(1:lastidx,-2+3*trial_opt)),'InputFormat','HH:mm:ss.SSS'));
            out = duration(out(:,4:end));
            out = seconds(out);
    
            VO2_dot{sub_sel}(1:lastidx,trial_opt)=table2array(Einfo_data{sub_sel}(1:lastidx,-1+3*trial_opt));
            RER{sub_sel}(1:lastidx,trial_opt)=table2array(Einfo_data{sub_sel}(1:lastidx,0+3*trial_opt));
            VCO2_dot{sub_sel}(1:lastidx,trial_opt)=RER{sub_sel}(1:lastidx,trial_opt).*VO2_dot{sub_sel}(1:lastidx,trial_opt);
        
        
            time{sub_sel}(1:lastidx,trial_opt)=out-out(1);
            last{sub_sel}(trial_opt)=lastidx;
        end
    end
end

%% Compute average value per trial
if computeAvgVal==1
    % minimum_time=[60; 180; 180; 180; 180; 180; 180; 180];    
    minimum_time=ones(1,numTrials)*(60*3); maximum_time=ones(1,numTrials)*(60*6);
%     minimum_time=ones(1,numTrials)*(60*10); maximum_time=ones(1,numTrials)*(60*12);
        VO2_dot_analyzed=double.empty; VCO2_dot_analyzed=double.empty; frame_detected_min=double.empty; 
        VO2_dot_average=double.empty;  VCO2_dot_average=double.empty;  frame_detected_max=double.empty; 
        
    for sub_opt=1:length(sub_list)
        sub_sel=sub_list(sub_opt);
        
        for trial_opt=1:numTrials
            VO2_dot_analyzed{sub_sel}(:,trial_opt)=VO2_dot{sub_sel}(:,trial_opt);
            VCO2_dot_analyzed{sub_sel}(:,trial_opt)=VCO2_dot{sub_sel}(:,trial_opt);
            
            max_tim=find(time{sub_sel}(:,trial_opt)>maximum_time(trial_opt),1,'first');
            if isempty(max_tim)
                frame_detected_max(sub_sel,trial_opt)=last{sub_sel}(trial_opt); 
            else
                frame_detected_max(sub_sel,trial_opt)=max_tim;
            end
           
            min_tim=find(time{sub_sel}(:,trial_opt)>minimum_time(trial_opt),1,'first');
            if isempty(min_tim)
                frame_detected_min(sub_sel,trial_opt)=frame_detected_max(sub_sel,trial_opt)-1*30; % -1*60
            else
                frame_detected_min(sub_sel,trial_opt)=min_tim;
            end
            
            VO2_dot_average(sub_sel,trial_opt) =mean( VO2_dot{sub_sel}(frame_detected_min(sub_sel,trial_opt):frame_detected_max(sub_sel,trial_opt),trial_opt ) );
            VCO2_dot_average(sub_sel,trial_opt)=mean(VCO2_dot{sub_sel}(frame_detected_min(sub_sel,trial_opt):frame_detected_max(sub_sel,trial_opt),trial_opt ) );

             VO2_dot_analyzed{sub_sel}(frame_detected_min(sub_sel,trial_opt):frame_detected_max(sub_sel,trial_opt),trial_opt)= VO2_dot_average(sub_sel,trial_opt);
            VCO2_dot_analyzed{sub_sel}(frame_detected_min(sub_sel,trial_opt):frame_detected_max(sub_sel,trial_opt),trial_opt)=VCO2_dot_average(sub_sel,trial_opt);
            
        end
    end
end


if plotAvgVal==1
    for sub_opt=1:length(sub_list)
        sub_sel=sub_list(sub_opt);
        for trial_opt=1:numTrials
            figure(trial_opt)
            clf
            hold on

            x_val=squeeze(time{sub_sel}(1:last{sub_sel}(trial_opt),trial_opt));
            plot(x_val, VO2_dot{sub_sel}(1:last{sub_sel}(trial_opt),trial_opt),'Color','b');
            plot(x_val,VCO2_dot{sub_sel}(1:last{sub_sel}(trial_opt),trial_opt),'Color','r');
            
            plot(x_val, VO2_dot_analyzed{sub_sel}(1:last{sub_sel}(trial_opt),trial_opt),'Color','b','LineWidth',2);
            plot(x_val,VCO2_dot_analyzed{sub_sel}(1:last{sub_sel}(trial_opt),trial_opt),'Color','r','LineWidth',2);
            
            xline(time{sub_sel}(frame_detected_min(sub_sel,trial_opt),trial_opt));
            xline(time{sub_sel}(frame_detected_max(sub_sel,trial_opt),trial_opt));

            EavgVal=(16.58*VO2_dot_average(sub_sel,trial_opt)*1000/60 + 4.51*VCO2_dot_average(sub_sel,trial_opt)*1000/60)/opt.subjects_mass(sub_sel);
            
            legend({'VO2-dot','VCO2-dot'})
            xlabel('time [s]'); ylabel('Gas volumen velocity [L/min]');
            title(strcat("Trial ", num2str(Einfo_speedKM_H(trial_opt,sub_sel)*1000/3600,'%4.2f'), " m/s - ", num2str(Einfo_speedKM_H(trial_opt,sub_sel)),'km/h EE:',num2str(EavgVal,'%4.2f'),'J/(s*kg)'))
            ylim([0 2])
            xlim([0 maximum_time(trial_opt)])
        end
    end
end

%%
if computeEE ==1
    for sub_opt=1:length(sub_list)
        sub_sel=sub_list(sub_opt);
        % Compute metabolic cost
        % P_met_gross= 16.58*VO2_dot + 4.51*VCO2_dot [W*s/mL]; thus inputs must be in [mL/s]
        % Metamax inputs: VO2_dot & VCO2_dot are in [L/min]

        EE_GROSS=double.empty; EE_NET=double.empty; COT_NET=double.empty; COT_GROSS=double.empty;
        SUB_SPEED=Einfo_speedM_S(:,sub_sel);
        
%         https://fisiologiadelejercicio.com/wp-content/uploads/2018/02/Differences-of-energy-expenditure-while.pdf
        plot_Kcal= 0;
        if plot_Kcal==1
            factor_con=60/4184; % 1 Kcal -> 4184 J
        else; factor_con=1; 
        end
        
        EE_GROSS(1)= EE_rest(sub_sel)*factor_con;
        COT_NET(1)= 0; COT_GROSS(1)= 0; %Metabolic power per distance (gross)

        for trial_opt=2:numTrials
            VO2_dot_average(sub_sel,trial_opt)= VO2_dot_average(sub_sel,trial_opt)*1000/60;
            VCO2_dot_average(sub_sel,trial_opt)=VCO2_dot_average(sub_sel,trial_opt)*1000/60;
            
            EE_GROSS(trial_opt)= (16.58*VO2_dot_average(sub_sel,trial_opt)+4.51*VCO2_dot_average(sub_sel,trial_opt))*factor_con; %Metabolic power
            EE_NET(trial_opt)  = EE_GROSS(trial_opt)-EE_GROSS(1);
            
            COT_GROSS(trial_opt)= EE_GROSS(trial_opt)/(SUB_SPEED(trial_opt)); %Metabolic power per distance (gross)    
            COT_NET(trial_opt)= (EE_GROSS(trial_opt)-EE_GROSS(1))/(SUB_SPEED(trial_opt));
        end

        normalizedByMassEE=1;
        if normalizedByMassEE==1 
            EE_GROSS=EE_GROSS/sub_mass(sub_sel); 
            EE_NET  =EE_NET/sub_mass(sub_sel); 
        end
        
        normalizedByMassCOT=1;
        if normalizedByMassCOT==1
            COT_GROSS=COT_GROSS/sub_mass(sub_sel); 
            COT_NET  =COT_NET/sub_mass(sub_sel); 
        end

        % organize
        initial=2; final=numTrials;
        
        EE_gross_analyzed(1)=EE_GROSS(1); EE_net_analyzed(1)=EE_NET(1); 
        NET_COT_analyzed(1)=COT_NET(1);   GROSS_COT_analyzed(1)=COT_GROSS(1);
        
        EE_gross_analyzed(2:final-initial+2)=EE_GROSS(initial:final); EE_net_analyzed(2:final-initial+2)=EE_NET(initial:final);
        NET_COT_analyzed(2:final-initial+2)=COT_NET(initial:final);   GROSS_COT_analyzed(2:final-initial+2)=COT_GROSS(initial:final);
        
        
        subject_vel_analyzed(1)=SUB_SPEED(1);
        subject_vel_analyzed(2:final-initial+2)=SUB_SPEED(initial:final); % Km/h -> m/s
        
        [sv_sort_value,sv_sort_index]=sort(subject_vel_analyzed);
        
        EE_gross_SORT=EE_gross_analyzed(sv_sort_index);   EE_net_SORT=EE_net_analyzed(sv_sort_index);
        GROSS_COT_SORT=GROSS_COT_analyzed(sv_sort_index); NET_COT_SORT=NET_COT_analyzed(sv_sort_index);
        
        
        SUB_SPEED=sv_sort_value;
        
        EE_GROSS=EE_gross_SORT;   EE_NET=EE_net_SORT;
        COT_GROSS=GROSS_COT_SORT; COT_NET=NET_COT_SORT;


        final=final-initial+2;
        initial=2;

        % compute regression
        grade=2;
        
        % GROSS ENERGY EXPENDITURE
        [coefficients,rsq,~,~] = applyRegression(grade,SUB_SPEED(initial:final),EE_GROSS(initial:final));
        GEE_x = linspace(0,2.5, 200);       GEE_y = polyval(coefficients, GEE_x);
        coefficients_GEE=coefficients;

        GEE_rsq=rsq; GEE_p=coefficients;
        [GEE_MIN,I]=min(GEE_y); GEE_ES= GEE_x(I);
        
        % NET ENERGY EXPENDITURE
        [coefficients,rsq,~,~] = applyRegression(grade,SUB_SPEED(initial:final),EE_NET(initial:final));
        NEE_x = linspace(0,2.5, 200);       NEE_y = polyval(coefficients, NEE_x);
        
        NEE_rsq=rsq; NEE_p=coefficients;
        [NEE_MIN,I]=min(NEE_y); NEE_ES= NEE_x(I);
        
        % GROSS COT
        [coefficients,rsq,~,~] = applyRegression(grade,SUB_SPEED(initial:final),COT_GROSS(initial:final));
        GCOT_x = linspace(0,2.5, 200);       GCOT_y = polyval(coefficients, GCOT_x);
        
        GCOT_rsq=rsq; GCOT_p=coefficients;
        [GCOT_MIN,I]=min(GCOT_y); GCOT_ES= GCOT_x(I);
        
        % NET COT
        [coefficients,rsq,~,~] = applyRegression(grade,SUB_SPEED(initial:final),COT_NET(initial:final));
        NCOT_x = linspace(0,2.5, 200);       NCOT_y = polyval(coefficients, NCOT_x);
        
        NCOT_rsq=rsq; NCOT_p=coefficients;
        [NCOT_MIN,I]=min(NCOT_y); NCOT_ES= NCOT_x(I);


        if plotEE ==1
            figure;    set(gcf,'color','w');
            
            subplot(2,2,1)
            hold on
            plot_EE(SUB_SPEED(initial:final),EE_GROSS(initial:final),GEE_x,GEE_y,900,GEE_p,GEE_rsq,GEE_ES,GEE_MIN,'gross EE [J/s]');
            axis([0.95*min(SUB_SPEED(initial:final)) 1.05*max(SUB_SPEED(initial:final)) 50*factor_con 1000*factor_con])
            
            subplot(2,2,2)
            hold on
            plot_EE(SUB_SPEED(initial:final),EE_NET(initial:final),NEE_x,NEE_y,900,NEE_p,NEE_rsq,NEE_ES,NEE_MIN,'net EE [J/s]');
            axis([0.95*min(SUB_SPEED(initial:final)) 1.05*max(SUB_SPEED(initial:final)) 50*factor_con 1000*factor_con])
            
            subplot(2,2,3)
            hold on
            plot_EE(SUB_SPEED(initial:final),COT_GROSS(initial:final),GCOT_x,GCOT_y,5.5,GCOT_p,GCOT_rsq,GCOT_ES,GCOT_MIN,'gross COT [J/kg/m]');
            plot(GCOT_ES,GCOT_MIN,'sk','MarkerSize',15)
            axis([0.95*min(SUB_SPEED(initial:final)) 1.05*max(SUB_SPEED(initial:final)) 1 6])
            
            subplot(2,2,4)
            hold on
            plot_EE(SUB_SPEED(initial:final),COT_NET(initial:final),NCOT_x,NCOT_y,5.5,NCOT_p,NCOT_rsq,NCOT_ES,NCOT_MIN,'net COT [J/kg/m]');
            plot(NCOT_ES,NCOT_MIN,'sk','MarkerSize',15)
            axis([0.95*min(SUB_SPEED(initial:final)) 1.05*max(SUB_SPEED(initial:final)) 1 6])
       end
            
            %%% SAVING INFORMATION
            EE_data(sub_sel).speed=SUB_SPEED(1:final);
            EE_data(sub_sel).coeff_GEE=coefficients_GEE;
            EE_data(sub_sel).GCOT_ES = GCOT_ES;
            
            EE_data(sub_sel).extCurve_GEE = GEE_y;
            EE_data(sub_sel).extCurve_NEE = NEE_y;
            EE_data(sub_sel).extCurve_GCOT=GCOT_y;
            EE_data(sub_sel).extCurve_NCOT=NCOT_y;
            
            EE_data(sub_sel).exp_GEE = EE_GROSS(1:final);
            EE_data(sub_sel).exp_NEE =   EE_NET(1:final);
            EE_data(sub_sel).exp_GCOT=COT_GROSS(1:final);
            EE_data(sub_sel).exp_NCOT=  COT_NET(1:final);
            
            EE_data(sub_sel).NEE_coe=NEE_p;
            EE_data(sub_sel).NEE_rsq=NEE_rsq;
            EE_data(sub_sel).GEE_coe=GEE_p;
            EE_data(sub_sel).GEE_rsq=GEE_rsq;
            EE_data(sub_sel).NCOT_coe=NCOT_p;
            EE_data(sub_sel).NCOT_rsq=NCOT_rsq;
            EE_data(sub_sel).GCOT_coe=GCOT_p;
            EE_data(sub_sel).GCOT_rsq=GCOT_rsq;
    end
end
%%
advance_computing=0;
if advance_computing==1
    clear p
    figure('Visible',"on",'color','w'); clf
    Einfo_data=double.empty; speed=double.empty;
    opt.sub_list =[1 2 4 6 9 13 14 15];
    color_EE_sub=["#80ff00" "#ffc000" "#ff0000" "#00c0ff" "#00ff40" "#ff00ff" "#0000ff"...
                  "#b23aee" "#ff0080" "#f43f1a" "#cd1076" "#9400d3" "#8b0a50" "#f4ac1a"...
                  "#8acae7" "#ffc100"];


    for sub_opt=1:length(opt.sub_list)
        sub_sel=opt.sub_list(sub_opt); 
        
        Einfo_data =EE_data(sub_sel).exp_GEE(2:end);
        speed=EE_data(sub_sel).speed(2:end);
        
        coe = EE_data(sub_sel).GEE_coe;
        rsq = EE_data(sub_sel).GEE_rsq;
        
        hold on;
        p(sub_opt)=plot(speed,Einfo_data,'Color',color_EE_sub(sub_sel),"LineStyle","none",'Marker',"o",'MarkerSize',10,'MarkerFaceColor',color_EE_sub(sub_sel)...
            ,'DisplayName', [ num2str(coe(1),'%+4.2f') 'X^2 ' num2str(coe(2),'%+4.2f') 'X^1 ' num2str(coe(3),'%+4.2f') ' -' ' r^2 =' num2str(rsq,'%4.2f')]);
    
        GEE_x = linspace(0,2.5, 200);       GEE_y = polyval(coe, GEE_x);
        
        plot(GEE_x,GEE_y,'Color',color_EE_sub(sub_sel),'LineWidth',2);
        
        axis([0.6 2.3 0 9]);
        xlabel('speed [m/s]','fontSize',30); ylabel('COT [J/Kg/m]','fontSize',30);
        set(gca,'FontSize',20)
    end
    legend(p,"Location","northwest",'FontSize',10); 
end

function[]= plot_EE(x_dis,y_dis,x_ext,y_ext,y_offset,coe_p,rsq_adj,ES,MIN,y_label)
plot(x_ext,y_ext,'Color','#77AC30','LineStyle',"-",'LineWidth',4);  plot(x_dis,y_dis,'or','MarkerSize',8);
% text(min(x_dis),y_offset,[strcat("Y=",num2str(coe_p(1),'%4.2f'),"X^2",num2str(coe_p(2),'%+4.2f'),"X^1",num2str(coe_p(3),'%+4.2f')); strcat("R^2 =",num2str(rsq_adj*100,'%4.2f'),"%")],'FontSize',10,"FontWeight","bold")
text(min(x_dis),y_offset,[strcat("Y=",num2str(coe_p(1),'%4.2f'),"X^2",num2str(coe_p(2),'%+4.2f'),"X^1",num2str(coe_p(3),'%+4.2f')); strcat("R^2 =",num2str(rsq_adj*100,'%4.2f'),"%"); 'Velocity:' num2str(ES,'%4.2f') ', Value:' num2str(MIN,'%4.2f')],'FontSize',8,"FontWeight","normal")
% text(ES-0.1,MIN-0.2,['Velocity:' num2str(ES,'%4.2f') ', Value:' num2str(MIN,'%4.2f')],'FontSize',10)
xlabel('speed [m/s]'); ylabel(y_label);
end