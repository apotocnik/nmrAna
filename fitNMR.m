%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fitNMR
%
%
% Anton Potoènik, 27.2.2011, F5 @ IJS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fitNMR(handles)
    mode = evalin('base', 'mode');
    TAU = evalin('base', 'TAU');
    D1 = evalin('base', 'D1');
    INTall = evalin('base', 'INTall');
    Niter = evalin('base', 'Niter');
    indFIT = evalin('base', 'indFIT');
    TEM = evalin('base', 'TEM');
    SIMtype = evalin('base', 'SIMtype');
    allinitpar = evalin('base', 'allinitpar');
    indFitpar = evalin('base', 'indFitpar');
    PLOT = evalin('base', 'PLOT');
    

    a4 = handles.axeRes;
%     cla(a4);
    % Get X,Y data
    if strcmp(mode,'D1')
        X = D1';
    else
        X = TAU';
    end
    Y = INTall(:,indFIT);

    [funstr allcoefs expN rules] = fitNMRlib(SIMtype);
    if numel(allcoefs)<numel(indFitpar)
        indFitpar = indFitpar(1:numel(allcoefs));
    end
    if numel(allcoefs)<numel(allinitpar)
        allinitpar = allinitpar(1:numel(allcoefs));
    end
    
%     allinitpar = allfitpar; % Always fresh fit
    allfiterr = zeros(size(allinitpar));    
    allfitpar = allinitpar;
    
    initpar = allinitpar(indFitpar);
    indprob = 1:numel(allinitpar); indprob(indFitpar)=[];
    probpar = allinitpar(indprob);

    coefs = allcoefs(indFitpar);   % names
    alcoefs = allcoefs;
    alcoefs(indFitpar)=[];
    problems = alcoefs;        % names
    
    % Set fit model
    model = fittype(funstr,'problem',problems,'coefficients',coefs);

    opts = fitoptions('Method','Nonlinear');
    opts.Algorithm =  'Trust-Region';            % 'Trust-Region'  'Levenberg-Marquardt'
    opts.Display =    'notify';                  % 'notify'  'off'  'iter'
    opts.MaxIter = Niter;
    opts.MaxFunEvals = 1000;
    opts.TolFun = 1e-10;
    opts.TolX = 1e-10;
    opts.Robust = 'Off';
    opts.StartPoint = initpar;        % Update start value

    if strcmp(opts.Algorithm,'Trust-Region')==1
        opts.Upper = inf*ones(numel(indFitpar),1);         % Update upper limits
        opts.Lower = 0*ones(numel(indFitpar),1);        % Update lower limits
        opts.Lower(1) = -inf;
    else
        opts.Upper = [];         % Delete upper limits
        opts.Lower = [];         % Delete lower limits
    end

    if isempty(allinitpar(indprob))
        [f1 gof] = fit(X, Y, model, opts);
    else
        [f1 gof] = fit(X, Y, model, opts, 'problem', num2cell(probpar));
    end

    par = coeffvalues(f1);        % Get fit coefficients values
    fitY = f1(X);
    sse = gof.sse;
    % Extract errors
    ci = confint(f1,0.95);      % Boundaries within 95%
    err = (ci(2,:)-ci(1,:))/2;  % Absolut error

    for i=1:numel(allfitpar)
       allfitpar(i) = 0;
       allfiterr(i) = 0;
       ind1 = find(indFitpar==i);
       ind2 = find(indprob==i);
       if ~isempty(ind1)
           allfitpar(i) = par(ind1);
           allfiterr(i) = err(ind1);
       elseif ~isempty(ind2)
           allfitpar(i) = probpar(ind2);
           allfiterr(i) = 0;
       else
           allfitpar(i) = 0;
           allfiterr(i) = 0;
       end
    end
    
    
    
    if PLOT
        
        hold(a4,'on');
        Tx = 0; Ty = 0.2; dTy = 0.07; Tc = 'r'; Tc1 = '--k';
        if strcmp(mode,'T2')
            Tx = 0.6; Ty = 0.3;
        end
              
        % Draw exponents:
        B = []; T1 = []; alpha = []; M = []; fitT1 = []; Tc1 = '--k';
        M(1) = allfitpar(1);
        plot(a4,X,M(1)*ones(size(X)),Tc1,'Linewidth',2)
        
%         expN = numel(regexp(funstr,'exp','start'));
%         if strcmp(SIMtype,'T1-1ExpT')
%             alpha(1:2) = allfitpar(4);
%             T1(1) = allfitpar(3);
%             T1(2) = allfitpar(3)/6;
%             B(1) = 0.1*allfitpar(2); 
%             B(2) = 0.9*allfitpar(2);   
%         elseif strcmp(SIMtype,'T2-1Exp1Gauss')
%             B = allfitpar(2);
%             T1 = allfitpar(3);
%             alpha = allfitpar(4);
%             C = allfitpar(5);  
%             W = allfitpar(6);
%         elseif strcmp(SIMtype,'T1-2Exp_13C_C60')
%             B(1) = allfitpar(2);
%             B(2) = B(1)/3;
%             T1(1) = allfitpar(3);
%             alpha(1) = allfitpar(4);
%             T1(2) = allfitpar(5);
%             alpha(2) = 1;
%         else

        % Extract parameters for decomposed plot
        for i=1:expN
          B(i) = allfitpar(3*(i-1)+2);
          T1(i) = allfitpar(3*(i-1)+3);
          if 3*(i-1)+4 > numel(allfitpar)
              alpha(i)=0;
          else
            alpha(i) = allfitpar(3*(i-1)+4);
          end
        end
        for i=1:numel(rules)
            eval(rules{i});
        end
%         end

        allPlot = [T1', alpha', B'];
        allPlot = sortrows(allPlot, 1);

%       if strcmp(mode,'T1')
        for i=1:expN
          M(i+1) = M(i) + allPlot(i,3);
          fitT1{i} = M(i) + allPlot(i,3)*(1-exp(-(X/allPlot(i,1)).^allPlot(i,2)));
          plot(a4,X,M(i+1)*ones(size(X)),Tc1,'Linewidth',2)
          plot(a4,X,fitT1{i},Tc1,'Linewidth',2)
        end
        
        % Determine y-axis limits
        min1 = M(1)-0.05*max(fitY);
        max1 = 1.1*max(fitY);
        mmin = min(INTall); mmax = max(INTall); mdif = mmax-mmin;
        mmmin = mmin-0.05*mdif; mmmax = mmax+0.05*mdif;
        ylim(a4,[min([mmmin min1]) max([mmmax max1])]);
%        ylim(a4,[min1 max2]);
        

%       end
        plot(a4,X,fitY,Tc,'Linewidth',2)
        text(0.05+Tx,Ty,['{\Delta}{^2} = ',num2str(round(1000*sse)/1000)],'Units','Norm','Fontweight','bold','FontSize',12,'Color',Tc,'Parent',a4,'BackgroundColor',[1 1 1])

        % Plot parameter values
        for i=1:numel(allfitpar)
            coef = allcoefs{i};
            par = allfitpar(i);
            err = allfiterr(i);
            partxt = val2txt(par,err,3);
            if ~isempty(strfind(coef,'T1')) || ~isempty(strfind(coef,'T2')) || ~isempty(strfind(coef,'D1'))
                partxt = [partxt 's'];
            end
            text(0.05+Tx,Ty+i*dTy,[coef ' = ', partxt],'Units','Norm','Fontweight','bold','FontSize',12,'Color',Tc,'Parent',a4,'BackgroundColor',[1 1 1]);
        end

        hold(a4,'off');
        
% %% Linearize
% %        cla(a4);
%         hold(a4,'off')
%         expN = numel(regexp(funstr,'exp','start'));
%         for i=1:expN
%             B(i) = allfitpar(3*(i-1)+2);
%             T1(i) = allfitpar(3*(i-1)+3);
%             alpha(i) = allfitpar(3*(i-1)+4);
%         end
% 
%         y0 = allfitpar(1);
%         M = sum(B);
%         plot(a4,TAU,(M-INTall+y0)/M,'o')
% %       set(a4,'YScale','log');
% %       set(a4,'XScale','log');
%         hold(a4,'on')
% 
% %        M = allfitpar(2);
%         plot(a4,X,(M-fitY+y0)/M,'-r')
%         hold(a4,'off')
%         grid(a4,'on')
   
        num2clip([TEM sse reshape([allfitpar;allfiterr],1,[])]);
        txt = [num2str(TEM) 'K  ' num2str(sse) '   '];
        for i=1:numel(allfitpar)
            txt = [txt val2txt(allfitpar(i),allfiterr(i),4) '   '];
        end
        disp(txt);
        gdisp(handles,'Parameters copied to clipboard');
    end
   
    
    
    assignin('base','model',model);
    assignin('base','opts',opts);
    assignin('base','sse',sse);
    assignin('base','fitINTall',fitY);
    assignin('base','f1',f1);
    assignin('base','allfitpar',allfitpar);
    assignin('base','allfiterr',allfiterr);
    assignin('base','fitT1',fitT1);
    
    