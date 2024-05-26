%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fftNMR
%
%
% Anton Potoènik, 26.2.2011 - 16.12.2013, F5 @ IJS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------------------------------------------------------
% Important variables:
%
% phiall = [xPar FID_auto_phase used_phase];
% SHLall = [xPar FID_auto_SHL FR used_SHL];
% S ... integrated SPC weights calculated after first FFT
%
%-------------------------------------------------------------------------

function ret = fftNMR(handles)

function r = getVar(where,var,alt)
    if evalin(where,['exist(''' var ''',''var'')']);
        eval(['r = evalin(''' where ''',''' var ''');'])
    else
        if ~ischar(alt), alt = num2str(alt); end;
        eval(['r = ' alt]);
    end
end

%======================================================================
% FFT initialization
%======================================================================

ret = 0;
mode = evalin('base', 'mode');
TAU = evalin('base', 'TAU');
D1 = evalin('base', 'D1');
LB = evalin('base', 'LB');
DE = evalin('base', 'DE');
SI = evalin('base', 'SI');
FR = evalin('base', 'FR');
IND = evalin('base', 'IND');
Tall = evalin('base', 'Tall');
FIDall = evalin('base', 'FIDall');
REF = evalin('base', 'REF');
LIM = evalin('base', 'LIM');
CUT = evalin('base','CUT');
avrPoints = getVar('base','avrPoints',1);

indN = numel(IND);

autoCorrInd = evalin('base', 'autoCorrInd');
SHL = evalin('base', 'SHL');
PH = evalin('base', 'PH');
autoSHL = evalin('base', 'autoSHL');
autoPH_SPC = evalin('base','autoPH_SPC');
PH_v = evalin('base','PH_v');
autoPH = evalin('base', 'autoPH');
range = evalin('base', 'range');
Fmin = range(1);
Fmax = range(end);
Dy = evalin('base', 'Dy');

if evalin('base','exist(''PLOT'',''var'')')
    PLOT = evalin('base', 'PLOT');
else
    PLOT = 1;
end

SPCall = cell(1,indN); 
Fall = cell(1,indN); 
phiall_spc = zeros(indN,2);

switch mode
    case {'T1','T2'}
        xPar = TAU;
    case 'D1'
        xPar = D1;
    case 'D2'
        D2 = evalin('base', 'D2');
        xPar = D2;
    case 'SW'
        xPar = FR;
    case 'SPC'
        xPar = TAU;
    otherwise
        error(['Unknown mode ' mode ' !!!']);
end

%----------------------------------------------------------------------
% Shuffle parameters
%----------------------------------------------------------------------

[~,SHF] = sort(xPar);

if sum(strcmp(mode,{'T1'})), SHF = fliplr(SHF); end;

xPar = xPar(SHF);
TAU = TAU(SHF);
D1 = D1(SHF);
FR = FR(SHF);
Tall = Tall(SHF);
FIDall = FIDall(SHF);
IND = IND(SHF);

FIDall0 = FIDall;
Tall0 = Tall;
    
phiall = zeros(indN,3); 
SHLall = zeros(indN,3);

%======================================================================
% Automatic Phase correciton and SHL determination
%======================================================================

% Calculate FID Weights
Exclude = zeros(indN,1);
Weights = ones(indN,1);


%----------------------------------------------------------------------
% Automatic SHL
%----------------------------------------------------------------------


for i=1:indN
    [~,ind] = max(abs(FIDall{i}));
    
    % Fit square function to the avrPoints to determine SHL
%     if avrPoints > 2
%         M = round((avrPoints-1)/2);
%         x = (ind-M:ind+M)';
%         sig = abs(FIDall{i}(x));
%         p = polyfit(x,sig,2);
%         if p(1) < 0 && p(2) > 0
%            ind = round(-p(2)/2/p(1));
%         end
%     end
    
    SHLall(i,:) = [xPar(i),ind,FR(i)];
end


% Suggestions
switch mode
    case {'D1','SW'}
        if autoSHL ~= 3 && autoSHL ~= 0
            disp('autoSHL should be 3 or 0 for D1 or SW mode!!!');
%             autoSHL = 3;
        end
    case {'T1','T2'}
        if autoSHL ~= 1 && autoSHL ~= 0
            disp('autoSHL should be 1 or 0 for T1 or T2 mode!!!');
%             autoSHL = 1;
        end
end

switch autoSHL
    case 0
        SHLs = SHL*ones(indN,1);

    case 1  % T1, T2, Mean SHL from last atuoCorrInd spectra
        Include  = intersect(autoCorrInd,1:numel(IND));
        if isempty(Include), Include = 1:numel(IND); end;
        Exclude = ones(indN,1); Exclude(Include) = 0;
        SHLmean = round(mean(SHLall(Include,2)));
        SHLs = SHLmean*ones(indN,1);

    case 2 
        SHLs = SHLall(:,2); % Individual SHL

    case 3 % D1, SW
        if evalin('base','exist(''S'',''var'')')
            Weights = evalin('base','S'); 
        else
            msgbox('Weights are not yet calculated. Please use different autoSHL!','Error');
            ret=-1; return
        end

        Exclude(Weights < CUT) = 1;
        if sum(Exclude) >= indN-2,
            msgbox('There is not enought points left for SHLmeanW. Please reduce CUT!','Error');
            ret=-1; return
        end
        P = fit(SHLall(:,1),SHLall(:,2),'poly1','Weight', Weights,'Exclude',Exclude);
        SHLmeanW = round(mean(P(SHLall(:,1))));
        SHLs = SHLmeanW*ones(indN,1);

    case 4 % D1, fit SHL
        if evalin('base','exist(''S'',''var'')')
            Weights = evalin('base','S'); 
        else
            msgbox('Weights are not yet calculated. Please use different autoSHL!','Error');
            ret=-1; return
        end

        Exclude(Weights < CUT) = 1;
        if sum(Exclude) >= indN-2,
            msgbox('There is not enought points left for fitting SHL. Please reduce CUT!','Error');
            ret=-1; return
        end
        P = fit(SHLall(:,1),SHLall(:,2),'poly1','Weight', Weights,'Exclude',Exclude);
        SHLs = round(P(SHLall(:,1)));
        
    otherwise
        error('Wrong autoSHL option!')
end

SHLall(:,4) = SHLs; % Store current SHLs to SHLall 4th column



%----------------------------------------------------------------------
% Automatic Phase correction
%----------------------------------------------------------------------

% extract individual phases
for i=1:indN
    M = round((avrPoints-1)/2);
    sig = FIDall{i}(SHLs(i)-M:SHLs(i)+M);
    phi0 = angle(mean(sig));
    phiall(i,:) = [xPar(i),phi0,0];
end

phiall(:,2) = phaseSmooth(mod(phiall(:,2),2*pi));

% Suggestions
switch mode
    case 'D1'
        if autoPH ~= 2
            disp('autoPH should be 2 for D1 mode!!!');
%             autoPH = 2;
        end
end


switch autoPH
    case 0
        phi = PH*ones(indN,1);
        
    case 1
        Include  = intersect(autoCorrInd,1:numel(IND));
        if isempty(Include), Include = 1:numel(IND); end;
        Exclude = ones(indN,1); Exclude(Include) = 0;
        PHImean = mean(phiall(Include,2));
        phi = PHImean*ones(indN,1);
        
    case 2
        phi = phiall(:,2)';
        
    case 3
        if evalin('base','exist(''S'',''var'')')
            Weights = evalin('base','S'); 
        else
            msgbox('Weights are not yet calculated. Please use different autoPH!','Error');
            ret=-1; return
        end

        Exclude(Weights < CUT) = 1;
        if sum(Exclude) >= indN-2,
            msgbox('There is not enought points left for fitting PH. Please reduce CUT!','Error');
            ret=-1; return
        end
        P = fit(phiall(:,1),phiall(:,2),'poly1','Weight', Weights,'Exclude',Exclude);
        PHImeanW = round(mean(P(phiall(:,1))));
        phi = PHImeanW*ones(indN,1);
        
    case 4 % Fit phase from previous results     
        if evalin('base','exist(''S'',''var'')')
            Weights = evalin('base','S'); 
        else
            msgbox('Weights are not yet calculated. Please use different autoPH!','Error');
            ret=-1; return
        end

        Exclude(Weights < CUT) = 1;
        if sum(Exclude) >= indN-2,
            msgbox('There is not enought points left for fitting PH. Please reduce CUT!','Error');
            ret=-1; return
        end
        
        try 
            order = evalin('base','PHorder');
        catch err
        	order = 'poly1'; 
        end
        order = inputdlg('Choose polynomial for phase fitting (poly1, poly2, ...):','Phase fit',1,{order});
        assignin('base','PHorder',order{1});
        P = fit(phiall(:,1),phiall(:,2),order{1},'Weight', Weights,'Exclude',Exclude);
        phi = P(phiall(:,1));

                        
%             Niter = 1000;
%             SIMtype = 'poly2';  initpar = [200,-50,-0.01];
%             ind = 1:indN;
% %             [fitpar,model] = FIT_PH_1(FR(ind)',Weights(ind),180/pi*phiall0(ind,3),REF,initpar,SIMtype,Niter);
%             [fitpar,model] = FIT_PH_1(FR(ind)',S(ind),180/pi*phiall0(ind,3),REF,initpar,SIMtype,Niter);
%             [sse,fitPH] = model(fitpar);
%             phi = fitPH*pi/180;         

    otherwise
        error('Wrong autoPH option!')
end




%======================================================================
% Perform FFT
%======================================================================

for i=1:indN   

    %------------------------------------------------------------------
    % Shift Left & Zero Filling 
    %------------------------------------------------------------------

    DW = Tall{i}(2)-Tall{i}(1);
    TD = length(FIDall{i});
    Tall{i} = (0:1:SI-1)'*DW;

    % Double Echo & Shift Left & Zero Filling
    if DE
        if length(FIDall{i}) >= SI
            FIDall{i} = FIDall{i}(1:SI);
        else
            FIDall{i} = [FIDall{i};zeros(SI-length(FIDall{i}),1)];
        end
        FIDall{i} = [FIDall{i}(1+SHLs(i):SI); FIDall{i}(1:SHLs(i))];
    else
        FIDall{i} = FIDall{i}(1+SHLs(i):TD);
        if length(FIDall{i}) >= SI
            FIDall{i} = FIDall{i}(1:SI);
        else
            FIDall{i} = [FIDall{i}; zeros(SI-length(FIDall{i}),1)];
        end                
    end


    %------------------------------------------------------------------
    % Phase Correction 
    %------------------------------------------------------------------

    FIDall{i} = FIDall{i}*exp(-complex(0,1)*phi(i));
    phiall(i,3) = phiall(i,3) + phi(i);


    %------------------------------------------------------------------
    % Line Broadening 
    %------------------------------------------------------------------
    if 1
        if DE
            LBfun = exp(-Tall{i}*LB/1e6)+exp(-(max(Tall{i})+DW-Tall{i})*LB/1e6);
            LBfun = LBfun/max(LBfun);
        else
            LBfun = exp(-Tall{i}*LB/1e6);
        end
        FIDall{i} = FIDall{i}.*LBfun;
    end

    %------------------------------------------------------------------
    % Fourier Transform 
    %------------------------------------------------------------------
    FID = FIDall{i};
    if DE
        % It is already periodic!!!
    else
        %%%%%%%%%%% 
        FID(1) = (FID(1)+FID(length(FID)))/2;
        %%%%%%%%%%%
    end

    SPC = fftshift(fft(FID));
    DF = 1/DW;
    Fall{i} = DF*(-SI/2:SI/2-1)'/SI*1e3;
    
    ind = find(Fall{i}>=2*Fmin & Fall{i}<=2*Fmax);
    Fall{i} = Fall{i}(ind);
    SPCall{i} = SPC(ind);
    
%     SPCall{i} = SPC;
end


%======================================================================
% Check SPC peak phases and corrects if set
%======================================================================

if autoPH_SPC > 0
    
    for i=1:indN
        [ind,cv] = searchclosest(Fall{i},PH_v);
        phi0 = angle(SPCall{i}(ind));
        phiall_spc(i,:) = [xPar(i),phi0];
    end

    switch autoPH_SPC
        case 1
            PHImean = mean(phiall_spc(intersect(autoCorrInd,IND),2));
            phi_spc = PHImean*ones(indN,1);

        case 2
            phi_spc = phiall_spc(:,2);
        
        otherwise
            error('Wrong autoPH_SPC option!')
    end

    for i=1:indN
        SPCall{i} = SPCall{i}*exp(-complex(0,1)*phi_spc(i));
        phiall(i,3) = phiall(i,3) + phi_spc(i);
    end
end


%======================================================================
% Plot all transformed spectra 
%======================================================================

if PLOT % SHL graph
    a1 = handles.axeSHL;
    % Exclude
    in = Exclude>0;
    if sum(in) > 0
        plot(a1,SHLall(in,1),SHLall(in,2),'o','Color',[0.5 0.5 0.5]); 
        hold(a1,'on');
    end
    % Include
    in = Exclude<1;
    plot(a1,SHLall(in,1),SHLall(in,2),'o'); hold(a1,'on');
    
%     plot(a1,[min(SHLall(:,1)) max(SHLall(:,1))], [SHLmean SHLmean],'g')
    plot(a1,SHLall(:,1), SHLall(:,4),'k');
    hold(a1,'off');
    grid(a1,'on');
    ylabel(a1,'SHL');
    switch mode
        case {'T1','T2'}
            set(a1,'XScale','log');
        case {'D1','SW','D2'}
            set(a1,'XScale','lin');
        otherwise
            set(a1,'XScale','lin');
    end
    if numel(SHLall(:,1)) > 1; 
        xlim(a1,[min(SHLall(:,1)) max(SHLall(:,1))]); 
        ylim(a1,[0.9*min(SHLall(in,2)) 1.1*max(SHLall(in,2))]); 
    end
    set(a1,'XTickLabel',[]);
end



if PLOT % Phase graph 
    a2 = handles.axePhi;
    in = Exclude<1;
    plot(a2,phiall(in,1),phiall(in,2),'o'); hold(a2,'on');
    in = Exclude>0;
    if sum(in) > 0
        plot(a2,phiall(in,1),phiall(in,2),'o','Color',[0.5 0.5 0.5]); 
    end
    
%     plot(a2,[min(phiall(:,1)) max(phiall(:,1))], [PHImean PHImean],'g');
    plot(a2,phiall(:,1), phiall(:,3),'k');
    hold(a2,'off');
    
    ylabel(a2,'Phase [rad]');
    grid(a2,'on');
    switch mode
        case {'T1','T2'}
            set(a2,'XScale','log');
            xlabel(a2,'\tau (ms)');
        case {'D1', 'D2'}
            set(a2,'XScale','lin');
            xlabel(a2,'D1 (\mus)');
        case 'SW'
            set(a2,'XScale','lin');
            xlabel(a2,'\nu (MHz)');
        case 'SPC'
            set(a2,'XScale','lin');
            xlabel(a2,'\tau (ms)');
        otherwise
            set(a2,'XScale','lin');
    end
    if numel(phiall(:,1)) > 1; 
        xlim(a2,[min(phiall(:,1)) max(phiall(:,1))]); 
    end
end



if PLOT
    NN = 50;
    if indN > NN
        step = round(indN/NN);
        if step > 1, disp(['Only every ' num2str(step) ' spc is shown!']); end;
    else
        step = 1;
    end
    ii = 1:step:indN; 
    
    %     clf;
    a3 = handles.axeSPC;
    cla(a3);
    NORM = zeros(1,indN);
    if strcmp(mode,'SW')
        for i=ii
            NORM(i) = max(real(SPCall{i}));
        end
        Dy = 1/step; % edini normalni prikaz za moje pojme
        assignin('base','Dy',Dy);
    end
    minVal = 1e10;
    maxVal = -1e10;
    for i=ii
        if strcmp(mode,'SW')
            ind = find(Fall{i}>=Fmin & Fall{i}<=Fmax);
            Y = real(SPCall{i})/max(NORM)+Dy*(i-1);
            if max(Y) > maxVal, maxVal = max(Y); end;
            if min(Y) < minVal, minVal = min(Y); end;   
            plot(a3,Fall{i}(ind),Y(ind),'k','Linewidth',1.5); %+1e3*FR(i)-1e3*REF
            hold(a3,'on');
            text(Fmin+(Fmax-Fmin)*0.0,0.20*Dy+Dy*(i-1),[num2str(FR(i)),' MHz ' num2str(i)],'Parent',a3)
        else
            Y = real(SPCall{indN+1-i})/max(real(SPCall{1}))+Dy*(i-1); 
            if max(Y) > maxVal, maxVal = max(Y); end;
            if min(Y) < minVal, minVal = min(Y); end;
            plot(a3,Fall{length(IND)+1-i},Y,'k','Linewidth',1.5)
            hold(a3,'on');
            if strcmp(mode,'D1')
                text(Fmin+(Fmax-Fmin)*0.05,0.25*Dy+Dy*(i-1),[num2str(D1(length(IND)+1-i)*1e6),' \mus ' num2str(IND(length(IND)+1-i))],'FontSize',12,'FontWeight','bold','Parent',a3)
            else
                text(Fmin+(Fmax-Fmin)*0.05,0.25*Dy+Dy*(i-1),[num2str(round(TAU(length(IND)+1-i)*1e6)/1000),' ms ' num2str(IND(length(IND)+1-i))],'FontSize',12,'FontWeight','bold','Parent',a3)
            end
        end
    end
    plot(a3,[PH_v PH_v],[minVal-0.05*maxVal,1.05*maxVal],'-r')
    
    % Plot LIM
    for j=1:length(LIM)
        plot(a3,[LIM(j),LIM(j)],[-1.05,Dy*(indN)+1+0.05],'-b','Linewidth',2)
    end
    
    hold(a3,'off'); grid(a3,'on');
    set(a3,'YTick',[]);
    xlabel(a3,'\nu [kHz]');
    xlim(a3,[Fmin,Fmax]);
    ylim(a3,[minVal-0.05*maxVal,1.05*maxVal]);

    drawnow
end



%======================================================================
% Save results and finish 
%======================================================================

SPCsum = zeros(length(real(SPCall{1})),2);

assignin('base','SPCall',SPCall);
assignin('base','SPCsum',SPCsum);
assignin('base','Fall',Fall);
assignin('base','SHLall',SHLall);
assignin('base','phiall',phiall);
assignin('base','phiall_spc',phiall_spc);
assignin('base','SHL',SHL);
assignin('base','SHF',SHF);
assignin('base','phi',phi);
assignin('base','IND',IND);
assignin('base','TAU',TAU);
assignin('base','FR',FR);
assignin('base','D1',D1);
assignin('base','Weights',Weights);
assignin('base','PLOT',PLOT);
assignin('base','FIDall',FIDall0);
assignin('base','Tall',Tall0);
assignin('base','FIDall1',FIDall);
assignin('base','Tall1',Tall);

if std(SHLs)==0
    gdisp(handles,['SHL = ' num2str(SHLs(1))]);
else
    gdisp(handles,['SHL = ' num2str(SHLs')]);
end

end