%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% integSW
%
% Routine for joining frequency sweep spectra
%
% Andraž Krajnc, 24.7.2011, FMF
% Anton Potoènik 19.12.2013, IJS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function integSW(handles)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PATHname = evalin('base', 'PATHname');
REF = evalin('base', 'REF');
FR = evalin('base', 'FR');
IND = evalin('base', 'IND');
Fall = evalin('base', 'Fall');
SPCall = evalin('base', 'SPCall');
PLOT = evalin('base', 'PLOT');
CUT = evalin('base', 'CUT');
isPPM = evalin('base', 'isPPM');
ppmRange = evalin('base', 'ppmRange');
NORM = evalin('base','NORM');
nParams = evalin('base','nParams');
Weights = evalin('base','Weights');
TEM = evalin('base','TEM');

Fmin = evalin('base','Fmin');
Fmax = evalin('base','Fmax');
Lmin = evalin('base','Lmin');
Lmax = evalin('base','Lmax');
SUMdelta = evalin('base','SUMdelta');
if evalin('base','exist(''Rng'',''var'')')
    Rng = evalin('base','Rng');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fmin = -100; Fmax = 100;
% SAVE = 0;
WEIGHT = 1;
% SUMdelta = 0.0005;
indN = numel(IND);  


%======================================================================
% Weighting
%======================================================================
     
if WEIGHT==1
    W = ones(indN); 
end
if WEIGHT==2
    if DeltaW==0
        W = ones(indN); 
    else
        W0 = 1-DeltaW:2*DeltaW/(N2D-1):1+DeltaW;
        W = [];
        for i=1:indN
            W = [W,W0];
        end
    end
end
if WEIGHT==3
    W = ones(indN);
    for i=1:size(indW,1)
        W(indW(i,1):indW(i,2))=W(indW(i,1):indW(i,2))*indW(i,3);
    end
end
if WEIGHT==4
    W = ones(indN);
    for i=1:size(indW,1)
        WfunX = 1:1:indW(i,2)-indW(i,1)+1;
        %WfunY = 1+Wa*(WfunX-(length(WfunX)+1)/2).^2;
        %WfunY = 1+Wa*(WfunX-(length(WfunX)+1)/2).^2;
        WfunY = 1+Wa*(WfunX-(length(WfunX)+1)/2);
        W(indW(i,1):indW(i,2))=WfunY;
    end
end

%======================================================================
% Cut spectra
%======================================================================

% LIM = LIM - PH_v; % get range from LIM without PH_v
% Lmin = min(LIM);
% Lmax = max(LIM);
% Fmin = Lmin;
% Fmax = Lmax;
for i=1:indN
    ind = find(Fall{i}<=Lmin | Fall{i}>=Lmax);
    SPCall{i}(ind) = [];
    Fall{i}(ind) = [];
end
% disp(['numel(ind) = ' num2str(numel(ind))]);



%======================================================================
% Integration 
%======================================================================

S = []; Sabs = [];
for i=1:indN
%     ind = find(Fall{i}>=Lmin & Fall{i}<=Lmax);
    S = [S;W(i)*sum(real(SPCall{i}))];
%     Sabs = [Sabs;W(i)*sum(abs(SPCall{i}))]; 
end

S = S/max(S);
% Sabs = Sabs/max(Sabs);
S(S<0)=0;


%======================================================================
% Summation
%======================================================================

SUMallX = [min(FR)+Fmin/1e3:SUMdelta:max(FR)+Fmax/1e3]';
%     SUMallX = [min(FR):SUMdelta:max(FR)]';
SUMallY = zeros(size(SUMallX));

for i=1:indN
    SUMY = interp1(Fall{i}/1e3+FR(i),real(SPCall{i}),SUMallX,'spline',0);
    SUMallY = SUMallY+W(i)*SUMY;
end

% Normalize
SUMallY = dataNorm([SUMallX SUMallY],NORM,nParams);
SUMallY = SUMallY(:,2);

if isPPM 
    SUMallX = (SUMallX - REF)/REF*1e6;
end


%======================================================================
% Plot
%======================================================================

if PLOT
    a4 = handles.axeRes;
    plot(a4,SUMallX,SUMallY,'k');
    hold(a4,'on');
    if isPPM
        xlabel(a4,'\nu-\nu_{ref} [ppm]'); 
        axis(a4,[ppmRange(1),ppmRange(2),-0.1*ceil(max(SUMallY)),1.1*ceil(max(SUMallY))]);
    else
        xlabel(a4,'\nu [MHz]'); 
        axis(a4,[min(SUMallX),max(SUMallX),-0.1*ceil(max(SUMallY)),1.1*ceil(max(SUMallY))]);
    end
    ylabel(a4,'Intensity');
    
    grid(a4,'on');
    hold(a4,'off');
    set(a4,'XScale','lin');
    
    a1 = handles.axeSHL;
    plot(a1,FR,S,'o',[min(FR) max(FR)],[CUT CUT],'g');
    ylabel(a1,'S weights');
    grid(a1,'on');
    set(a1,'XScale','lin');
    set(a1,'XTickLabel',[]);
    if indN > 1; 
        xlim(a1,[min(FR),max(FR)]);
    end
    
%     a2 = handles.axePhi;
%     plot(a2,FR,S,'o');
%     ylabel(a2,'S');
%     grid(a2,'on');
%     set(a2,'XScale','lin');
%     xlabel(a2,'\nu [MHz]'); 
%     if indN > 1; 
%         xlim(a2,[min(FR),max(FR)]);
%     end
    drawnow
end

%======================================================================
% Test Spectra
%======================================================================

if exist('Rng','var')
    ind = find(SUMallX>Rng(1) & SUMallX<Rng(2));
    if PLOT
        lim = get(a4,'YLim');
        hold(a4,'on');
        plot(a4,[Rng(1) Rng(1)],[lim(1) lim(2)],'r');
        plot(a4,[Rng(2) Rng(2)],[lim(1) lim(2)],'r');
        hold(a4,'off');
    end
    gdisp(handles,['mean: ' num2str(mean(SUMallY(ind))) '  std: ' num2str(std(SUMallY(ind)))])
end

    
%======================================================================
% Save & Export
%======================================================================

if isPPM
    nmrFitData.data{1} = [SUMallX SUMallY];
else
    nmrFitData.data{1} = [SUMallX*1e3 SUMallY];
end
nmrFitData.TEs{1} = TEM;
assignin('base','nmrFitData',nmrFitData);

SPCsum = [SUMallX SUMallY];
% num2clip(SPCsum);
assignin('base','SPCsum',SPCsum);
assignin('base','S',S);

gdisp(handles,'SW specter copied to clipboard. nmrFitData ready.');

