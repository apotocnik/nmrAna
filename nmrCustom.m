function [ list ] = nmrCustom(ind, hObject, eventdata, handles)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
list = {'Edit functions...', 'Get T1s from allfitparM', 'momentLS','T1 line scan','Run time dependence','Save movie as','Complex T1 analyse','Improved Spectral Extraction'};
switch ind
    case 1
        openCustom;
    case 2
        getT1s;
    case 3
        momentLS;
    case 4
        T1linescan;
    case 5
        runTimeDep;
    case 6
        savemovie;
    case 7
        T1onPeaks;
    case 8
        ImpSpectre;
    otherwise
end

function openCustom
    edit nmrCustom;
end

function getT1s
    % This works only when Analyse Run has been previousely done with
    % "T1-2Exp_Afact_Tfact" model, and now "T1-2Exp_Tfact" model is used for
    % linescan
    
    TEM = evalin('base','TEM');
    allfitparM = evalin('base','allfitparM');
    TEs = fliplr(evalin('base','TEs'));
    
    allfitpar = allfitparM(TEs==TEM,:);
    
    % Change Afact to B2=B1*Afact
    allfitpar(5) = allfitpar(2)*allfitpar(5);
    
    % Write to edtAnalyse (as in butCpyParam)
    param = get(handles.edtAnalyse,'String');
    for i=1:numel(param)
        if strfind(param{i},'allinitpar = ')
            str = []; 
            for j=1:numel(allfitpar)
                str = [str num2str(allfitpar(j)) ','];
            end
            str = ['[' str(1:end-1) ']'];
            param{i} = ['allinitpar = ' str ';'];
        end
    end
    set(handles.edtAnalyse,'String',param);

end


function momentLS
    TEM = evalin('base','TEM');
    range = evalin('base','range');
    isPPM = evalin('base','isPPM');
    ppmRange = evalin('base','ppmRange');
    X = evalin('base','Fscan');
    INTscan = evalin('base','INTscan');
    N = size(INTscan,2);
    a4 = handles.axeRes;
    out = TEM;
    colors={'r','b','g','k','c','m'};
    
    for i=1:N
        Y = INTscan(:,i);
        cutoff = max(Y)/1.11111;
        [M0 M1 M2] = nrNMRanalysis(X,Y,0,cutoff);
        out = [out M0 M1 M2];
        plot(a4,X,Y,colors{i})
        hold(a4,'on')
        plot(a4,[M1 M1],[0 max(Y)],[colors{i} '-.'],'LineWidth',2)
        plot(a4,[M1-M2 M1-M2],[0 max(Y)],[colors{i} ':'])
        plot(a4,[M1+M2 M1+M2],[0 max(Y)],[colors{i} ':'])
        plot(a4,[min(X) max(X)],[cutoff cutoff],[colors{i} '--'])
    end
    if isPPM
        range = ppmRange;
        xlabel(a4,'\nu - \nu_r_e_f [ppm]')
    else
        xlabel('\nu [kHz]')
    end
    
    xlim(a4,range);
    hold(a4,'off')
%     num2clip(out)
    assignin('base','out',out);
%     gdisp(handles,'Moments saved to clipboard.')
end



function savemovie
colorNum = 256;

if evalin('base','exist(''fmovie'',''var'')')
    if evalin('base','exist(''movieFPS'',''var'')')
        movieFPS = evalin('base','movieFPS');
    else
        movieFPS = 12;
    end
    [FileName,PathName] = uiputfile({'*.avi';'*.gif'},'Save');
    if ~isequal(FileName, 0)
        if strcmp(FileName(end-2:end),'gif')
            fmovie = evalin('base','fmovie');
            [im,map] = rgb2ind(fmovie(1).cdata,colorNum,'nodither');
            im(1,1,1,numel(fmovie)) = 0;
            for k = 1:numel(fmovie) 
              im(:,:,1,k) = rgb2ind(fmovie(k).cdata,map,'nodither');
            end
            imwrite(im,map,[PathName FileName],'DelayTime',0,'LoopCount',1) %LoopCount - inf or 1
        elseif strcmp(FileName(end-2:end),'avi')
            fmovie = evalin('base','fmovie');
            movie2avi(fmovie, [PathName FileName], 'compression','i420','fps',movieFPS);
        end
        msgbox(['Saved to: ' PathName FileName]);
    end    
end
end

function runTimeDep
    allfiles = get(handles.dirFiles,'String');
    allframes = numel(allfiles);
    hours = 0; sumLR = [];
    REF = evalin('base','REF');
    LIM1 = [52.965 52.975];
    LIM2 = [52.975 52.985];
    set(handles.dirFiles,'Value',1);
    nmrAna('dirFiles_Callback',hObject, eventdata,handles);
    nmrAna('btnLoad_Callback',hObject, eventdata,handles);
    text(0.7,0.85,['Time: +' num2str(round(hours)) ' h'],'Units','Norm','Fontweight','bold','FontSize',16,'Color','b','Parent',handles.axeRes,'BackgroundColor',[1 1 1]);
    fmovie(1) = getframe(handles.axeRes);
    clipboardMtmp = evalin('base','clipboardM');
    clipboardM = zeros(numel(clipboardMtmp(:,1)),allframes+1);
    clipboardM(:,1:2) = clipboardMtmp;
    indL = find(clipboardMtmp(:,1)*REF*1e-6+REF<LIM1(2)&clipboardMtmp(:,1)*REF*1e-6+REF>LIM1(1));
    indR = find(clipboardMtmp(:,1)*REF*1e-6+REF<LIM2(2)&clipboardMtmp(:,1)*REF*1e-6+REF>LIM2(1));
    sumLR = [sumLR; hours trapz(clipboardMtmp(indL,1),clipboardMtmp(indL,2)) trapz(clipboardMtmp(indR,1),clipboardMtmp(indR,2))];
    
    hold(handles.axeRes,'on');
    plot(handles.axeRes,[LIM1(1) LIM1(1)],[0 1],'color','b');plot(handles.axeRes,[LIM1(2) LIM1(2)],[0 1],'color','b');
    plot(handles.axeRes,[LIM2(1) LIM2(1)],[0 1],'color','g');plot(handles.axeRes,[LIM2(2) LIM2(2)],[0 1],'color','g');
    for i=2:allframes
       hours = 1.183*(i-1);
       set(handles.dirFiles,'Value',i);
       nmrAna('dirFiles_Callback',hObject, eventdata,handles);
       nmrAna('btnLoad_Callback',hObject, eventdata,handles);
       text(0.7,0.85,['Time: +' num2str(round(hours)) ' h'],'Units','Norm','Fontweight','bold','FontSize',16,'Color','b','Parent',handles.axeRes,'BackgroundColor',[1 1 1]);
       fmovie(i) = getframe(handles.axeRes);
       clipboardMtmp = evalin('base','clipboardM');
       clipboardM(:,i+1) = clipboardMtmp(:,2);
        indL = find(clipboardMtmp(:,1)*REF*1e-6+REF<LIM1(2)&clipboardMtmp(:,1)*REF*1e-6+REF>LIM1(1));
        indR = find(clipboardMtmp(:,1)*REF*1e-6+REF<LIM2(2)&clipboardMtmp(:,1)*REF*1e-6+REF>LIM2(1));
        sumLR = [sumLR; hours trapz(clipboardMtmp(indL,1),clipboardMtmp(indL,2)) trapz(clipboardMtmp(indR,1),clipboardMtmp(indR,2))];
        
        hold(handles.axeRes,'on');
        plot(handles.axeRes,[LIM1(1) LIM1(1)],[0 1],'color','b');plot(handles.axeRes,[LIM1(2) LIM1(2)],[0 1],'color','b');
        plot(handles.axeRes,[LIM2(1) LIM2(1)],[0 1],'color','g');plot(handles.axeRes,[LIM2(2) LIM2(2)],[0 1],'color','g');
    end
    assignin('base','fmovie',fmovie);
    assignin('base','clipboardM',clipboardM);
    assignin('base','sumLR',sumLR);
    figure(1);clf;plot(sumLR(:,1),sumLR(:,2),'color','b'); hold on; plot(sumLR(:,1),sumLR(:,3),'color','g'); hold off;
    %movie(handles.axeRes,fmovie,1,3);
end

function T1linescan
     param = get(handles.edtFFTparam,'String');
     PH_v_table = evalin('base','PH_v_table');
     datatable = [];
     for j=1:numel(PH_v_table)
         for i=1:numel(param)
                if strfind(param{i},'PH_v = ')
                    str = num2str(PH_v_table(j));
                    %str(str==' ')=',';
                    param{i} = ['PH_v = ' str ';'];
                end
        end
            set(handles.edtFFTparam,'String',param);
            evalin('base',sprintf('%s\n',param{:}));
            
            nmrAna('btnLoad_Callback',hObject, eventdata,handles);
            nmrAna('butAnalyse_Callback',hObject, eventdata, handles);
            nmrAna('butParams_Callback',hObject, eventdata, handles);
            drawnow();
            allfitpar = evalin('base','allfitpar');
            allfiterr = evalin('base','allfiterr');
            sse = evalin('base','sse');
            datatable = [datatable; PH_v_table(j) sse reshape([allfitpar;allfiterr],1,[])];
     end
     
     num2clip(datatable);
     gdisp(handles,'Parameters copied to clipboard');
end

% function T1linescan
%      param = get(handles.edtFFTparam,'String');
%      PH_v_table = evalin('base','PH_v_table');
%      datatable = [];
%      for j=1:numel(PH_v_table)
%          for i=1:numel(param)
%                 if strfind(param{i},'PH_v = ')
%                     str = num2str(PH_v_table(j));
%                     %str(str==' ')=',';
%                     param{i} = ['PH_v = ' str ';'];
%                 end
%         end
%             set(handles.edtFFTparam,'String',param);
%             evalin('base',sprintf('%s\n',param{:}));
%             
%             nmrAna('btnLoad_Callback',hObject, eventdata,handles);
%             nmrAna('butAnalyse_Callback',hObject, eventdata, handles);
%             nmrAna('butParams_Callback',hObject, eventdata, handles);
%             drawnow();
%             allfitpar = evalin('base','allfitpar');
%             allfiterr = evalin('base','allfiterr');
%             sse = evalin('base','sse');
%             datatable = [datatable; PH_v_table(j) sse reshape([allfitpar;allfiterr],1,[])];
%      end
%      
%      num2clip(datatable);
%      gdisp(handles,'Parameters copied to clipboard');
% end

function T1onPeaks
    TEs = fliplr(evalin('base','TEs'));
    param_table = {[],[],[],[]};
    shift_table = [];
    for k=1:numel(TEs)
        set(handles.edtTE,'String',num2str(TEs(k)));
        nmrAna('btnLoad_Callback',hObject, eventdata,handles);
        nmrAna('butSPCfit_Callback',hObject, eventdata,handles);
        shift_table = [shift_table; TEs(k) evalin('base','SPCfit_lastpar')];
        nmrAna('butSPCfitrun_Callback',hObject, eventdata,handles);
        for m=1:4
           evalin('base',['indFIT = ' num2str(m) ';']);
           nmrAna('butAnalyse_Callback',hObject, eventdata,handles);
           allfitpar = evalin('base','allfitpar');
           allfiterr = evalin('base','allfiterr');   
           param_table{m} = [param_table{m}; TEs(k) reshape([allfitpar;allfiterr],1,[])];
        end
        drawnow();
    end
    assignin('base','param_table',param_table);
    assignin('base','shift_table',shift_table);
end

% function D1fractionfit
%     a4 = handles.axeRes; Niter = 10000; TEM = evalin('base','TEM');
%     mode = evalin('base','mode');
%     if strcmp(mode,'D1')
%         INTall = evalin('base','INTall');
%         X = evalin('base','1e6*D1''');
%         X = X(1:end-1);
%         INTall = INTall(1:end-1);
%     elseif strcmp(mode,'P1')
%         INTall = evalin('base','INTall');
%         X = evalin('base','1e6*P1''');
%     end;
%     spinA = 1/2; spinB = 3/2;
%     spinRatio = (1/(spinB+1/2))/(1/(spinA+1/2));
%     expDrop = 0; expFactor = 5;
%     %X = [1:0.1:12]';INTall = 1*abs(sin(X*pi/12))+1*abs(sin(2*X*pi/12));
%     Y = INTall(:,1)/max(max(INTall));
%     %plot(a4,X,Y);
%     if expDrop
%         funstr = ['A*abs(exp(-(x-x0)/(' num2str(expFactor) '*2*D1/pi))*sin((x-x0)*pi/(2*D1)))+B*abs(exp(-' num2str(spinRatio) '*(x-x0)/(' num2str(expFactor) '*2*D1/pi))*sin(' num2str(spinRatio) '*(x-x0)*pi/(2*D1)))'];
%     else
%         funstr = ['A*abs(sin((x-x0)*pi/(2*D1)))+B*abs(sin(' num2str(spinRatio) '*(x-x0)*pi/(2*D1)))'];
%     end
%     allcoefs = {'A','B','D1','x0'};
%     
%     indFitpar=[2,3,4];
%     allinitpar = [0,1,1.8,1];
%     
%         if numel(allcoefs)<numel(indFitpar)
%         indFitpar = indFitpar(1:numel(allcoefs));
%     end
%     if numel(allcoefs)<numel(allinitpar)
%         allinitpar = allinitpar(1:numel(allcoefs));
%     end
%     
% %     allinitpar = allfitpar; % Always fresh fit
%     allfiterr = zeros(size(allinitpar));    
%     allfitpar = allinitpar;
%     
%     initpar = allinitpar(indFitpar);
%     indprob = 1:numel(allinitpar); indprob(indFitpar)=[];
%     probpar = allinitpar(indprob);
% 
%     coefs = allcoefs(indFitpar);   % names
%     alcoefs = allcoefs;
%     alcoefs(indFitpar)=[];
%     problems = alcoefs;        % names
%     
%     model = fittype(funstr,'problem',problems,'coefficients',coefs);
% 
%     opts = fitoptions('Method','Nonlinear');
%     opts.Algorithm =  'Trust-Region';            % 'Trust-Region'  'Levenberg-Marquardt'
%     opts.Display =    'notify';                  % 'notify'  'off'  'iter'
%     opts.MaxIter = Niter;
%     opts.MaxFunEvals = 1000;
%     opts.TolFun = 1e-10;
%     opts.TolX = 1e-10;
%     opts.Robust = 'Off';
%     opts.StartPoint = initpar;        % Update start value
%     
%     if strcmp(opts.Algorithm,'Trust-Region')==1
%         opts.Upper = inf*ones(numel(indFitpar),1);         % Update upper limits
%         opts.Lower = 0*ones(numel(indFitpar),1);        % Update lower limits
%         opts.Lower(1) = -inf;
%     else
%         opts.Upper = [];         % Delete upper limits
%         opts.Lower = [];         % Delete lower limits
%     end
%     
%     if isempty(allinitpar(indprob))
%         [f1 gof] = fit(X, Y, model, opts);
%     else
%         [f1 gof] = fit(X, Y, model, opts, 'problem', num2cell(probpar));
%     end
%     
%     par = coeffvalues(f1);       % Get fit coefficients values
%     Xsmth = [0:0.01:max(X)]';
%     fitY = f1(Xsmth);
%     sse = gof.sse;
%     % Extract errors
%     ci = confint(f1,0.95);      % Boundaries within 95%
%     err = (ci(2,:)-ci(1,:))/2;  % Absolut error
%     
%     
%      for i=1:numel(allfitpar)
%        allfitpar(i) = 0;
%        allfiterr(i) = 0;
%        ind1 = find(indFitpar==i);
%        ind2 = find(indprob==i);
%        if ~isempty(ind1)
%            allfitpar(i) = par(ind1);
%            allfiterr(i) = err(ind1);
%        elseif ~isempty(ind2)
%            allfitpar(i) = probpar(ind2);
%            allfiterr(i) = 0;
%        else
%            allfitpar(i) = 0;
%            allfiterr(i) = 0;
%        end
%      end
%     %end drawing
%     Tx = 0; Ty = 0.2; dTy = 0.07; Tc = 'r'; Tc1 = '--k';
% hold(a4,'on');
% plot(a4,Xsmth,fitY,Tc,'Linewidth',2)
% 
% %plot parts:
% COLOR = {'b','m','c','g','k','r'};
% if expDrop
%     plot(a4,Xsmth,allfitpar(1)*abs(exp(-(Xsmth-allfitpar(4))/(expFactor*2*allfitpar(3)/pi)).*sin((Xsmth-allfitpar(4))*pi/(2*allfitpar(3)))),['--' COLOR{1}],'Linewidth',2);
%     plot(a4,Xsmth,allfitpar(2)*abs(exp(-spinRatio*(Xsmth-allfitpar(4))/(expFactor*2*allfitpar(3)/pi)).*sin(spinRatio*(Xsmth-allfitpar(4))*pi/(2*allfitpar(3)))),['--' COLOR{2}],'Linewidth',2);
% else
%     plot(a4,Xsmth,allfitpar(1)*abs(sin((Xsmth-allfitpar(4))*pi/(2*allfitpar(3)))),['--' COLOR{1}],'Linewidth',2);
%     plot(a4,Xsmth,allfitpar(2)*abs(sin(spinRatio*(Xsmth-allfitpar(4))*pi/(2*allfitpar(3)))),['--' COLOR{2}],'Linewidth',2);
% end
% 
% text(0.05+Tx,Ty,['{\Delta}{^2} = ',num2str(round(1000*sse)/1000)],'Units','Norm','Fontweight','bold','FontSize',12,'Color',Tc,'Parent',a4,'BackgroundColor',[1 1 1])
%         
%         for i=1:numel(allfitpar)
%             coef = allcoefs{i};
%             par = allfitpar(i);
%             err = allfiterr(i);
%             if i==1 || i==2
%                 textStr = [coef ' = ', num2str(round(par*1e3)/1e3) ' \pm ' num2str(round(err*1e3)/1e3) ' (' num2str(round(100*par/(allfitpar(1)+allfitpar(2)))) '%)'];
%             else
%                 textStr = [coef ' = ', num2str(round(par*1e3)/1e3) ' \pm ' num2str(round(err*1e3)/1e3)];
%             end
%             text(0.05+Tx,Ty+i*dTy,textStr,'Units','Norm','Fontweight','bold','FontSize',12,'Color',Tc,'Parent',a4,'BackgroundColor',[1 1 1]);
%         end        
%        
%         
%         hold(a4,'off');
%         
%         num2clip([TEM sse reshape([allfitpar;allfiterr],1,[])]);
%         num2clip([Xsmth fitY allfitpar(1)*abs(exp(-(Xsmth-allfitpar(4))/(expFactor*2*allfitpar(3)/pi)).*sin((Xsmth-allfitpar(4))*pi/(2*allfitpar(3)))) allfitpar(2)*abs(exp(-spinRatio*(Xsmth-allfitpar(4))/(expFactor*2*allfitpar(3)/pi)).*sin(spinRatio*(Xsmth-allfitpar(4))*pi/(2*allfitpar(3))))]);
%         gdisp(handles,'Parameters copied to clipboard');
% end


function SPCintegrate
SPCall = evalin('base','SPCall');
Fall = evalin('base','Fall');
LIM = evalin('base','ppmLIM');
Dy = evalin('base','Dy');
mode = evalin('base','mode');
COLOR = {'b','m','c','g','k','r'};

if sum(strcmp(mode,{'T1','T2'}))
    xname = 'TAU [ms]';
    xpar = evalin('base','TAU');
elseif strcmp(mode,'D1')
    xname = 'D1 [\mus]';
    xpar = evalin('base','D1')*1e6;
elseif strcmp(mode,'P1')
    xname = 'P1 [mus]';
    xpar = evalin('base','P1')*1e6;
elseif strcmp(mode,'SW')
    xname = 'T [K]';
    xpar = evalin('base','TEs');
end

a3 = handles.axeSPC;
a4 = handles.axeRes;

hold(a3,'on');
for j=1:length(LIM)
    plot(a3,[LIM(j),LIM(j)],[-1.05,Dy*(numel(SPCall))+1+0.05],'--b','Linewidth',2)
end
hold(a3,'off');

    
cla(a4);
hold(a4,'on');
maxI = 0;
for j=1:numel(LIM)-1
    Iall = [];
    LEG{j} = [num2str(LIM(j)),':',num2str(LIM(j+1)),' ppm'];
        for i=1:numel(SPCall)
            ind = find(Fall{i}>LIM(j)&Fall{i}<LIM(j+1));
            integral = trapz(Fall{i}(ind),real(SPCall{i}(ind)));
            %plot(Fall,SPCall); hold on; plot(Fall(ind),SPCall(ind),'linewidth',5);
            newI = integral*xpar(i);
            Iall = [Iall; newI];
            if maxI<newI
                maxI = newI;
            end
        end
    plot(a4,xpar, Iall,'o','Color',COLOR{j});
end
xlabel(a4,xname)
ylabel(a4,'Integral*T [arb. units]')
ylim(a4,[0 1.1*maxI]);
legend(a4,LEG,2,'Location','NorthWest');
hold(a4,'off');

end


function deadCode
    pressure = {'0kN','10kN','11kN','11.8kN','12kN','14kN','15kN','17kN','19kN','21kN','25kN','40kN'};
kbars = [0 0.5 0.9 2.3 1.7 2.9 3.8 4.2 6.1 7.8 9.4 14.2];
%pressure = {'21kN','25kN','40kN'};
if 0
   cut = [0:0.01:0.5];
   tabela = [];
   for i=1:numel(pressure)
        pressure2 = pressure{i};
        pressure2(pressure2=='.')='-';
        if strcmp(pressure{i},'0kN')
            set(handles.edtFN,'String',['D:\dropbox\Dropbox\Measurements\AK\AM2_4\' pressure{i} '\T1\Cs3_T1_250K-001.DAT']);
        else
            set(handles.edtFN,'String',['D:\dropbox\Dropbox\Measurements\AK\AM2_4\' pressure{i} '\T1\AM24_' pressure2 '_T1_250K-001.DAT']);
        end
        if strcmp(pressure{i},'10kN')
            set(handles.edtTE,'String','[150]');
        else
            set(handles.edtTE,'String','[250]');
        end
        btnLoad_Callback(hObject, eventdata, handles)
%        for cnum=1:numel(cut)
%             nastavi cut
%             param = get(handles.edtAnalyse,'String');
% 
%             for j=1:numel(param)
%                 if strfind(param{j},'CUT = ')
%                     param{j} = ['CUT = ''' cut(cnum) ''';'];
%                 end
%             end
% 
%             set(handles.edtAnalyse,'String',param);
%             param = get(handles.edtAnalyse,'String');
%             for j=1:numel(param)
%                 evalin('base',param{j});
%             end

            butShift_Callback(hObject, eventdata, handles)
            clipboardM = evalin('base','clipboardM');
            TEs = clipboardM(:,1);
            M1table=clipboardM(:,2);
            M2table=clipboardM(:,5);
     
            ind = find(TEs>=80&TEs<=100);
            M1range=M1table(ind);
            M1 = mean(M1range);
            M1err = sqrt(sum((M1-M1range).^2)/numel(M1range));
            if M1err==0
                M1err = sqrt(sum((M1-M1table(ind-1:ind+1)).^2)/numel(M1table(ind-1:ind+1)));
            end
            M2range=M2table(ind);
            M2 = mean(M2range);
            M2err = sqrt(sum((M2-M2range).^2)/numel(M2range));
            if M2err==0
                M2err = sqrt(sum((M2-M2table(ind-1:ind+1)).^2)/numel(M2table(ind-1:ind+1)));
            end
%        end
tabela = [tabela; M1 M1err M2 M2err];
    end
assignin('base','tabela',tabela);
num2clip(tabela);
msgbox('Konec!');
return; 
end

if 0
M1table = zeros(45,2*numel(pressure));
M2table = zeros(45,2*numel(pressure));
Vmaxtable = zeros(45,2*numel(pressure));
for i=1:numel(pressure)
        pressure2 = pressure{i};
        pressure2(pressure2=='.')='-';
        if strcmp(pressure{i},'0kN')
            set(handles.edtFN,'String',['D:\dropbox\Dropbox\Measurements\AK\AM2_4\' pressure{i} '\T1\Cs3_T1_250K-001.DAT']);
        else
            set(handles.edtFN,'String',['D:\dropbox\Dropbox\Measurements\AK\AM2_4\' pressure{i} '\T1\AM24_' pressure2 '_T1_250K-001.DAT']);
        end
        if strcmp(pressure{i},'10kN')
            set(handles.edtTE,'String','[150]');
        else
            set(handles.edtTE,'String','[250]');
        end
        btnLoad_Callback(hObject, eventdata, handles)
        butShift_Callback(hObject, eventdata, handles)
        clipboardM = evalin('base','clipboardM');
        M1table(1:numel(clipboardM(:,1)),2*i-1)=clipboardM(:,1);
        M2table(1:numel(clipboardM(:,1)),2*i-1)=clipboardM(:,1);
        Vmaxtable(1:numel(clipboardM(:,1)),2*i-1)=clipboardM(:,1);
        M1table(1:numel(clipboardM(:,1)),2*i)=toPPM(clipboardM(:,2));
        M2table(1:numel(clipboardM(:,1)),2*i)=toPPM(clipboardM(:,5));
        Vmaxtable(1:numel(clipboardM(:,1)),2*i)=toPPM(clipboardM(:,4));
end
assignin('base','M1table',M1table);
assignin('base','M2table',M2table);
assignin('base','Vmaxtable',Vmaxtable);
return;
end
if 0
%name = {'3exp_1','3exp_2','3exp_2a','3exp_3','3exp_4','3exp_5','3exp_6','3exp_7','3exp_8','2exp_1','2exp_2','2exp_3','2exp_4','2exp_5'};
name = {'3exp_fraction'};
assignin('base','saveIMG',1);
assignin('base','autoRun',1);
for i=1:numel(pressure)
    for j=1:numel(name)
        pressure2 = pressure{i};
        pressure2(pressure2=='.')='-';
        if strcmp(pressure{i},'0kN')
            set(handles.edtFN,'String',['D:\dropbox\Dropbox\Measurements\AK\AM2_4\' pressure{i} '\T1\Cs3_T1_250K-001.DAT']);
        else
            set(handles.edtFN,'String',['D:\dropbox\Dropbox\Measurements\AK\AM2_4\' pressure{i} '\T1\AM24_' pressure2 '_T1_250K-001.DAT']);
        end
        if strcmp(pressure{i},'10kN')
            set(handles.edtTE,'String','[150]');
        else
            set(handles.edtTE,'String','[250]');
        end
        setName(handles,name{j});
        btnLoad_Callback(hObject, eventdata, handles)
        
%         nameOUT = [getFilePath(handles) 'LS-' name{j}];
%         plotLSsettings(handles,nameOUT);
        butLSana_Callback(hObject, eventdata, handles)
    end
end
assignin('base','saveIMG',0);
assignin('base','autoRun',0);
end

%% subtract spectra
clipboardM = evalin('base','clipboardM');
TEs = evalin('base','TEs');
ranges = [-500 300];
TEs2 = fliplr(TEs);
TEs0 = 300; NORM = 0;
clipboardTMP = [clipboardM(:,1) clipboardM(:,2)-TEs0];
indTMP = find(clipboardTMP(:,1)>ranges(1)&clipboardTMP(:,1)<ranges(2));

for i=1:numel(clipboardM(1,:))/2
   ind = find(clipboardM(:,1+2*(i-1))>ranges(1)&clipboardM(:,1+2*(i-1))<ranges(2));
   if NORM
        clipboardM(ind,2+2*(i-1)) = clipboardM(ind,2+2*(i-1))-clipboardTMP(indTMP,2)*(1+TEs0/TEs2(i))/2;
   else
        clipboardM(ind,2+2*(i-1)) = clipboardM(ind,2+2*(i-1))-clipboardTMP(indTMP,2);
   end
end

        a4 = handles.axeRes;
        plot(a4,clipboardM(:,1),clipboardM(:,2))
        [max_val,location] = max(max(clipboardM(:,1),[],2));
        hold(a4,'on');
        text(max_val,clipboardM(location,2),[num2str(TEs(1)),' K'],'FontSize',10,'FontWeight','bold','Color','b','Parent',a4)
        for j=2:length(TEs)
            plot(a4,clipboardM(:,2*j-1),clipboardM(:,2*j))
            if mod(j,4)==0
                [max_val,location] = max(max(clipboardM(:,2*j-1),[],2));
                text(max_val,clipboardM(location,2*j),[num2str(TEs(j)),' K'],'FontSize',10,'FontWeight','bold','Color','b','Parent',a4)
            end
        end
        hold(a4,'off');
        %xlim(a4,[Fmin Fmax])
        xlabel(a4,'\nu - \nu_r_e_f [ppm]')

num2clip(clipboardM);
gdisp(handles,'Data subtracted and copied to clipboard.');   
        
        
end


function ImpSpectre
%     spcSumInd = evalin('base','spcSumInd');
    SPCall = evalin('base','SPCall');
    Fall = evalin('base','Fall');
    TAU = evalin('base','TAU');
    range = evalin('base','range');
    ppmRange = evalin('base','ppmRange');
    autoPH = evalin('base','autoPH');
    PLOT = evalin('base','PLOT');
    TEM = evalin('base','TEM');
    IND = evalin('base','IND');
    NORM = evalin('base','NORM');
    nParams = evalin('base','nParams');
    isPPM = evalin('base', 'isPPM');
    % reloadSettings(handles);
    param = get(handles.edtAnalyse,'String');
    evalin('base',sprintf('%s\n',param{:}));
    pars = evalin('base', 'allinitpar');
    
    
    prompt = {'Recovery fit parameters [A B T1 alpha]:', 'Omit spectra below this absoute value:'};
    dlg_title = 'Improved Spectrum extraction';
    num_lines = 1;
    def = {mat2str(pars),'0.5'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    if isempty(answer)
        return
    end
    
    par = str2num(answer{1});
    omit = str2num(answer{2});

    
    % Initialize
    SPCsum = zeros(numel(real(SPCall{1})),2);
    SPCsum(:,1) = Fall{1};

    
    for i=1:length(TAU)
        Y = par(1)+par(2)*(1-exp(-(TAU(i)/par(3))^par(4)));
        if abs(Y) > omit
            SPCsum(:,2) = SPCsum(:,2) + Y*real(SPCall{i});
        end
    end

    
%     % Normalize
    SPCsum = dataNorm(SPCsum,NORM,nParams);

    if isPPM 
        SPCsum(:,1) = toPPM(SPCsum(:,1));
    end

%         % Save spectrum to file
%         nameOUT = [getFilePath(handles) mode '_' num2str(TEM) 'K.txt'];
%         save(nameOUT,'SPCsum','-ascii');

    if PLOT
        a4 = handles.axeRes;
        plot(a4,SPCsum(:,1),SPCsum(:,2))
        if isPPM
            range = ppmRange;
            xlabel(a4,'\nu - \nu_r_e_f [ppm]')
        else
            xlabel(a4,'\nu [kHz]') 
        end
        if max(SPCsum(:,1)) < range(1) || min(SPCsum(:,1)) > range(2)
           msgbox('Data outside specified range!','Error!')
           xlim(a4,[min(SPCsum(:,1)) max(SPCsum(:,1))])
        else
           xlim(a4,range)
        end
        grid(a4,'on')
        drawnow;
    end

    if autoPH == 0
        gdisp(handles,'Warning: autoPH = 0, not good for summing spectra!')
    end


    num2clip(SPCsum);
    assignin('base','SPCsum',SPCsum);

    nmrFitData.data{1} = SPCsum;
    nmrFitData.TEs{1} = TEM;
    assignin('base','nmrFitData',nmrFitData);
    gdisp(handles,'Current specter copied to clipboard. Sum of all spectra is in SPCsum. nmrFitData ready.');

end
        


end

