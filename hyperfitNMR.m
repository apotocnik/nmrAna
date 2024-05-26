%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% hyperfitNMR
%   fits all T1 data at hte same time
%
% Anton Potoènik, 27.2.2011 - 16.04.2013, F5 @ IJS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function hyperfitNMR(handles)
    mode = evalin('base', 'mode');
    TAU = evalin('base', 'TAU');
    drawIND = evalin('base', 'drawIND');

    Niter = evalin('base', 'Niter');
    SPCall = evalin('base', 'SPCall');
    Fall = evalin('base', 'Fall');
    isPPM = evalin('base', 'isPPM');
    DyF = evalin('base', 'DyF');
%     neglowlim = evalin('base', 'neglowlim');
    TEM = evalin('base', 'TEM');
    SIMtype = evalin('base', 'SIMtype');
    allinitpar = evalin('base', 'allinitpar');
    indFitpar = evalin('base', 'indFitpar');
    PLOT = evalin('base', 'PLOT');
    
    if ~strcmpi(mode,'T1')
        msgbox('This function is intendt only for T1 measurements!!!')
        return
    end
    
    
    SAVE = 0;
    a4 = handles.axeRes;
    
    if isPPM 
        F = toPPM(Fall{1});
    else
        F = Fall{1};
    end
    Z = real(cell2mat(SPCall));
    
    Nf = numel(F);
    Nfdesired = 150;
    dNf = round(Nf/Nfdesired);  if dNf<1; dNf=1; end;
    filter = 1:dNf:Nf;
    F = F(filter);
    Nf = numel(F);
    
    
    
    Nt = numel(TAU);
    X = repmat(reshape(F,[],1),Nt,1);
    Y = repmat(reshape(TAU,[],1),1,Nf);
    Y = reshape(Y',[],1);
    
    
    Z = Z(filter,:);
    Z = reshape(Z,[],1);

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
    model = fittype(funstr,'problem',problems,'coefficients',coefs,'independent',{'x','t'});

    opts = fitoptions('Method','Nonlinear');
    opts.Algorithm =  'Levenberg-Marquardt';            % 'Trust-Region'  'Levenberg-Marquardt'
    opts.Display =    'iter';                  % 'notify'  'off'  'iter'
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
%         if exist('neglowlim','var')
%            opts.Lower(neglowlim) = -inf;
%         end
    else
        opts.Upper = [];         % Delete upper limits
        opts.Lower = [];         % Delete lower limits
    end

    disp('Fit start ...')
    if isempty(allinitpar(indprob))
        [f1 gof] = fit([X,Y], Z, model, opts);
    else
        [f1 gof] = fit([X,Y], Z, model, opts, 'problem', num2cell(probpar));
    end
    disp('done')
    
    par = coeffvalues(f1);        % Get fit coefficients values
    fitZ = f1(X,Y);
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
        cla(a4);
        set(a4,'XScale','lin');
        hold(a4,'on');
        Tx = 0; Ty = 0.2; dTy = 0.07; Tc = 'r';
              
        % Draw exponents:
        B = []; T1 = []; M = []; Tc1 = '--k';

        % Extract parameters for decomposed plot
        for i=1:expN
          B(i) = allfitpar(3*(i-1)+2);
          w(i) = allfitpar(3*(i-1)+3);
          xc(i) = allfitpar(3*(i-1)+4);
          d(i) = allfitpar(3*(i-1)+5);
          T1(i) = allfitpar(3*(i-1)+6);
        end
        for i=1:numel(rules)
            eval(rules{i});
        end

        
        % Plot spectra
        minVal = 1e10;
        maxVal = -1e10;
        for i=drawIND
            Y = real(SPCall{Nt+1-i})+DyF*(i-1); 
            if max(Y) > maxVal, maxVal = max(Y); end;
            if min(Y) < minVal, minVal = min(Y); end;
            if isPPM 
                plot(a4,toPPM(Fall{Nt+1-i}),Y,'k','Linewidth',1.5)
            else
                plot(a4,Fall{Nt+1-i},Y,'k','Linewidth',1.5)
            end
            hold(a4,'on');
%             text(Fmin+(Fmax-Fmin)*0.05,0.25*Dy+Dy*(i-1),[num2str(round(TAU(length(IND)+1-i)*1e6)/1000),' ms ' num2str(IND(length(IND)+1-i))],'FontSize',12,'FontWeight','bold','Parent',a3)
        end
        
        if isPPM 
            xlabel(a4,'\nu [ppm]','FontSize',12);
        else
            xlabel(a4,'\nu [kHz]','FontSize',12);
        end
       
%         set(a4,'YTick',[]);
%         axis(a4,[Fmin,Fmax,minVal-0.05*maxVal,1.05*maxVal]);


        % Plot fits
        for i=1:expN
            for j=drawIND
                Z = f1(F,TAU(Nt-j+1));
                Z = Z + DyF*(j-1); 
                if max(Z) > maxVal, maxVal = max(Z); end;
                if min(Z) < minVal, minVal = min(Z); end;
                plot(a4,F,Z,'r','Linewidth',1.5)
%                 disp(TAU(Nt-j+1))
            end
        end
        

         dmVal = 0.1*(maxVal-minVal);
         ylim(a4,[minVal-dmVal maxVal+dmVal]);
         xlim(a4,[min(F) max(F)])
        
        text(0.05+Tx,Ty,['{\Delta}{^2} = ',num2str(round(1000*sse)/1000)],'Units','Norm','Fontweight','bold','FontSize',12,'Color',Tc,'Parent',a4,'BackgroundColor',[1 1 1])

        % Plot parameter values
        gdisp(handles,'======')
        for i=1:numel(allfitpar)
            coef = allcoefs{i};
            par = allfitpar(i);
            err = allfiterr(i);
            gdisp(handles,[coef ' = ', num2str(round(par*1e4)/1e4) ' +- ' num2str(round(err*1e4)/1e4)]);
%             text(0.05+Tx,Ty+i*dTy,[coef ' = ', num2str(round(par*1e3)/1e3) ' \pm ' num2str(round(err*1e3)/1e3)],'Units','Norm','Fontweight','bold','FontSize',12,'Color',Tc,'Parent',a4,'BackgroundColor',[1 1 1]);            
        end

        hold(a4,'off');
   
%         saveDATA = [TAU',INTall];
         num2clip([TEM sse reshape([allfitpar;allfiterr],1,[])]);
         gdisp(handles,'Parameters copied to clipboard');
    end

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save
    
    if SAVE==1
        saveDATA = [TAU',INTall];
        save([PATHname,nameOUT,num2str(TEM),'KR2.txt'],'saveDATA','-ascii')

        if sumLastN > 0
            saveData = real(SPCall{1});
            for i=2:sumLastN
                saveDATA = saveData + real(SPCall{i});
            end
            saveDATAspc = [Fall{1} saveData/max(saveData)]; % Normalize Summed Spectra
            save([PATHname,sumNameOUT,num2str(TEM),'.txt'],'saveDATAspc','-ascii')
            figure(5);
            plot(saveDATAspc(:,1),saveDATAspc(:,2))
            xlim([Fmin Fmax])
            xlabel('\nu [kHz]')

            if autoPH == 0
                gdisp(handles,'Warning: autoPH = 0, not good for summing spectra!')
            end
        end
    end    
    
    
    assignin('base','model',model);
    assignin('base','opts',opts);
    assignin('base','sse',sse);
    assignin('base','f1',f1);
    assignin('base','allfitpar',allfitpar);
    assignin('base','allfiterr',allfiterr);
    
    