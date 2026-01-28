function [coefficients_flip,rsq,r,predict,rsq_adj,LM] = regression(grade,x,y)

%     coefficients = polyfit(x,y,grade);
%     yfit = polyval(coefficients,x); %yfit =  p(1) * x + p(2);
%     yresid = y - yfit;
%     SSresid = sum(yresid.^2);
%     SStotal = (length(y)-1) * var(y);
%     rsq = 1 - SSresid/SStotal; % r^2 https://se.mathworks.com/help/matlab/data_analysis/linear-regression.html
%     rsq_adj = 1 - SSresid/SStotal * (length(y)-1)/(length(y)-length(coefficients)); % NEW IMPLEMENTATION, WATCH OUT FOR THIS!
    
%     predict = SSresid/(length(y) - grade -1); % https://autarkaw.org/2008/07/05/finding-the-optimum-polynomial-order-to-use-for-regression/
    
    if grade==1; type='linear'; elseif grade==2; type='purequadratic'; elseif grade==3; type='poly3'; end
    
        LM = fitlm(x,y,type);
        rsq=LM.Rsquared.Ordinary; % coefficient of determination
        r  =sqrt(rsq);            % correlation coefficient
        rsq_adj=LM.Rsquared.Adjusted;
        coefficients_flip=flip(LM.Coefficients.Estimate); % now first value is x1 and second value is intercept
        predict=[];
end