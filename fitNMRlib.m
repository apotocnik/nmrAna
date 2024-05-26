%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fitNMRlib
% 
% Function library for fitNMR
%
% - A (offset) is automatically added as the first parameter
% - Functions are decomposed first by "T1" and "T2" prefix
%             and after they are decomposed by a "-" sign
%
% - Maintain the order: B, T1, alpha !!!
%
% Anton Potoènik, 27.2.2011-1.4.2013, F5 @ IJS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [funstr,allcoefs,expN,rules] = fitNMRlib(SIMtype)
    funstr = ''; allcoefs = '';
    expN = 0;
    rules = {};
    parts = regexp(SIMtype,'-','split');
    Lib = parts{1};

        for i=1:numel(parts)-1
            [funstrtmp,allcoefstmp,Nexp,rul] = getSimpleForm(Lib,parts{i+1});
            % Add number suffix to the variables 
            for j=1:numel(allcoefstmp)
                funstrtmp=strrep(funstrtmp,allcoefstmp{j},[allcoefstmp{j} num2str(i)]);
                allcoefstmp{j} = [allcoefstmp{j} num2str(i)];
            end
            funstr = [funstr funstrtmp ' + '];
            allcoefs = [allcoefs allcoefstmp];
            expN = expN + Nexp;
            rules = [rules rul];
        end
    % we need also offset and semicolon
    funstr = [funstr 'A;'];
    allcoefs = ['A' allcoefs];

    if isempty(funstr), error('Wrong SIMtype! \nChoose: T1-1mS,...,T1-7mS,T1-2exp_Tfact,T1-2exp_Bfact_Tfact,T2-exp.'); end;
    
    
    
    function [funstr,allcoefs,expN,rules] = getSimpleForm(Lib,SIMtype)
        
    expN = 1;
    rules = {};
    if strcmpi(Lib,'T1')
             
            if strcmpi('1mS',SIMtype)
                funstr = 'B*(1-exp(-(x/T1).^alpha))';
                allcoefs = {'B','T1','alpha'};
            end
            if strcmpi('3mS',SIMtype)
                funstr = 'B*(1-1/10*exp(-(x/T1)^alpha)-9/10*exp(-(6*x/T1)^alpha))';
                allcoefs = {'B','T1','alpha'};
            end
            if strcmpi('5mS',SIMtype)
                funstr = 'B*(1-1/35*exp(-(x/T1).^alpha)-8/45*exp(-(6*x/T1).^alpha)-50/63*exp(-(15*x/T1).^alpha))';
                allcoefs = {'B','T1','alpha'};
            end
            if strcmpi('7mS',SIMtype)
                funstr = 'B*(1-1/84*exp(-(x/T1).^alpha)-3/44*exp(-(6*x/T1).^alpha)-75/364*exp(-(15*x/T1).^alpha)-1225/1716*exp(-(28*x/T1).^alpha))';
                allcoefs = {'B','T1','alpha'};
            end
            
            if strcmpi('2exp_Tfact',SIMtype)
                funstr = 'B1*(1-exp(-(x/T1).^alpha1)) + B2*(1-exp(-(x/T1/Tfact).^alpha2))';
                allcoefs = {'B1','T1','alpha1','B2','Tfact','alpha2'};
                expN=2;
                rules = {'T1(2)=T1(1)*T1(2)'};
            end
            if strcmpi('2exp_bfact_Tfact',SIMtype)
                funstr = 'B*(1-exp(-(x/T1).^alpha1)) + bfact*B*(1-exp(-(x/T1/Tfact).^alpha2))';
                allcoefs = {'B','T1','alpha1','bfact','Tfact','alpha2'};
                expN=2;
                rules = {'T1(2)=T1(1)*T1(2);','B(2)=B(1)*B(2);'};
            end
            if strcmpi('2exp_bfact_Tfact_1alpha',SIMtype)
                funstr = 'B*(1-exp(-(x/T1).^alpha)) + bfact*B*(1-exp(-(x/T1/Tfact).^alpha))';
                allcoefs = {'B','T1','alpha','bfact','Tfact'};
                expN=2;
                rules = {'T1(2)=T1(1)*T1(2);','B(2)=B(1)*B(2);','alpha(2)=alpha(1);'};
            end
            if strcmpi('3exp_122',SIMtype)
                funstr = 'B*(1-exp(-(x/T11).^alpha)) + 2*B*(1-exp(-(x/T12).^alpha)) + 2*B*(1-exp(-(x/T13).^alpha))';
                allcoefs = {'B','T11','T12','T13','alpha'};
                expN=1;

            end
            
%             if strcmpi('2Exp',SIMtype)
%                 funstr = 'B1*(1-exp(-(x/T1a).^alpha1)) + B2*(1-exp(-(x/T1b).^alpha2))';
%                 allcoefs = {'B1','T1a','alpha1','B2','T1b','alpha2'};
%             end
%             if strcmpi('2Exp_13C_C60',SIMtype) % C60 has 3 C in 1:2:2 ratio.
%                 funstr = 'B*(1-exp(-(x/T1a).^alpha)) + B/3*(1-exp(-x/T1b))';
%                 allcoefs = {'B','T1a','alpha','T1b'};
%             end
%             if strcmpi('2Exp',SIMtype)
%                 funstr = 'B1*(1-exp(-(x/T11).^alpha1)) + B2*(1-exp(-(x/T12).^alpha2))';
%                 allcoefs = {'B1','T11','alpha1','B2','T12','alpha2'};
%             end
%             if strcmpi('3Exp_fraction',SIMtype)
%                 funstr = '(fraction/(1-fraction))*B2*(1-exp(-(x/T11).^alpha1)) + B2*(1-exp(-(x/T12).^alpha2)) + B3*(1-exp(-(x/T13).^alpha3))';
%                 allcoefs = {'fraction','T11','alpha1','B2','T12','alpha2','B3','T13','alpha3'};
%             end
%             if strcmpi('3Exp_fraction_teflon',SIMtype)
%                 funstr = '(fraction/(1-fraction))*(B2+((teflon*B2)/(1-teflon-fraction)))*(1-exp(-(x/T11).^alpha1)) + B2*(1-exp(-(x/T12).^alpha2)) + ((teflon*B2)/(1-teflon-fraction))*(1-exp(-(x/T13).^alpha3))';
%                 allcoefs = {'fraction','T11','alpha1','B2','T12','alpha2','teflon','T13','alpha3'};
%             end
%             if strcmpi('3Exp_teflon',SIMtype)
%                 funstr = 'B1*(1-exp(-(x/T11).^alpha1)) + B2*(1-exp(-(x/T12).^alpha2)) + (fraction*(B1+B2)/(1-fraction))*(1-exp(-(x/T13).^alpha3))';
%                 allcoefs = {'B1','T11','alpha1','B2','T12','alpha2','fraction','T13','alpha3'};
%             end
%             if strcmpi('2Exp_fraction',SIMtype)
%                 funstr = '(fraction/(1-fraction))*B2*(1-exp(-(x/T11).^alpha1)) + B2*(1-exp(-(x/T12).^alpha2))';
%                 allcoefs = {'fraction','T11','alpha1','B2','T12','alpha2'};
%             end

    elseif strcmpi(Lib,'T2')
            if strcmpi('exp',SIMtype)
                funstr = 'B*exp(-(2*x/T2).^alpha)';
                allcoefs = {'B','T2','alpha'};
            end
            if strcmpi('gauss',SIMtype)
                funstr = 'C*exp(-2*(x/W).^2)';
                allcoefs = {'C','W'};
            end
%             if strcmpi('1Exp1Gauss',SIMtype)
%                 funstr = 'B*exp(-(2*x/T2).^alpha) + C/W/sqrt(2*pi)*exp(-2*(x/W).^2)';
%                 allcoefs = {'B','T2','alpha','C','W'};
%             end
%             if strcmpi('2exp',SIMtype)
%                 funstr = 'B1*exp(-(2*x/T21).^alpha1) + B2*exp(-(2*x/T22).^alpha2)';
%                 allcoefs = {'B1','T21','alpha1','B2','T22','alpha2'};
%             end
    elseif strcmpi(Lib,'hypT1')
            if strcmpi('1Gauss',SIMtype)
                funstr = 'B/w/sqrt(pi/2)*exp(-2*(x-xc)*(x-xc)/w/w)*(1-d*exp(-t/T1))';
                allcoefs = {'B','w','xc','d','T1'};
            end
            if strcmpi('TTpGauss',SIMtype)
                funstr = 'Ba/wa/sqrt(pi/2)*exp(-2*(x-xca)*(x-xca)/wa/wa)*(1-d*exp(-t/T1))+Bb/wb/sqrt(pi/2)*exp(-2*(x-xcb)*(x-xcb)/wb/wb)*(1-d*exp(-t/T1))';
                allcoefs = {'Ba','wa','xca','d','T1','Bb','wb','xcb'};
%                 rules = {'T1(3)=T1(2);'};
            end
    end