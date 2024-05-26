%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% integNMR
%
%
% Anton Potoènik, 26.2.2011, F5 @ IJS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function integNMR(handles)
    mode = evalin('base', 'mode');
    TAU = evalin('base', 'TAU');
    FR = evalin('base', 'FR');
    TEM = evalin('base', 'TEM');
    D1 = evalin('base', 'D1');
    Fall = evalin('base', 'Fall');
    SPCall = evalin('base', 'SPCall');
    Dy = evalin('base', 'Dy');
    LIM = evalin('base', 'LIM');
    PLOT = evalin('base', 'PLOT');
    
    indN = length(Fall);
      
    INTall = [];
    for k=1:indN
        INT0 = [];
        if size(LIM,1)==1
            INT0 = zeros(1,length(LIM)-1);
            for j=1:length(LIM)-1
                ind = find(Fall{k}>=LIM(j)&Fall{k}<=LIM(j+1));
                INT0(j) = sum(real(SPCall{k}(ind)))/numel(ind);    
            end
        end
        if size(LIM,1)==2
            ind = find((Fall{k}>=LIM(1,1)&Fall{k}<=LIM(1,2))|(Fall{k}>=LIM(2,1)&Fall{k}<=LIM(2,2)));
            INT0 = sum(real(SPCall{k}(ind)))/numel(ind);    
        end            
        INTall = [INTall;INT0];
    end
    
    
    
%======================================================================
% S weight
%======================================================================
    Lmin = min(LIM);
    Lmax = max(LIM);
    S = []; Sabs = [];
    for i=1:indN
        ind = find(Fall{i}>=Lmin & Fall{i}<=Lmax);
        S = [S;sum(real(SPCall{i}(ind)))];
    end
    S = S/max(S);
    S(S<0) = 0;
    
    
%======================================================================
% Plot
%======================================================================    
    if PLOT
    a4 = handles.axeRes;

    LEG = {};
    COLOR = {'k','b','r','g','m','k','b','r','g','m'};
    LW = [2,2,2,2,2,1,1,1,1,1];
    for j=1:length(LIM)-1
        if sum(strcmp(mode,{'T1','T2'}))
            plot(a4,TAU,INTall(:,j),'o','Linewidth',LW(j),'Color',COLOR{j})
        elseif strcmp(mode,'D1')
            plot(a4,D1*1e6,INTall(:,j)/max(max(INTall)),'o','Linewidth',LW(j),'Color',COLOR{j})
        elseif strcmp(mode,'D2')
            D2 = evalin('base','D2');
            plot(a4,D2*1e6,INTall(:,j)/max(max(INTall)),'o','Linewidth',LW(j),'Color',COLOR{j})
        elseif strcmp(mode,'SW')
            plot(a4,FR,INTall(:,j)/max(max(INTall)),'o','Linewidth',LW(j),'Color',COLOR{j})
            % 
        end
        hold(a4,'on');
        text(0.45,0.95,['T = ',num2str(TEM),' K'],'Units','Norm','Fontweight','bold','FontSize',14,'Color','black','Parent',a4,'BackgroundColor',[1 1 1])
        LEG{j} = [num2str(LIM(j)),':',num2str(LIM(j+1)),' kHz'];
    end

    if sum(strcmp(mode,{'T1'}))
        set(a4,'XScale','log');
        legend(a4,LEG,2,'Location','NorthWest');
    else
        legend(a4,LEG,2,'Location','NorthEast');
    end
    if sum(strcmp(mode,{'T2'}))
        set(a4,'YScale','log');
    else
        set(a4,'YScale','lin');
    end
    grid(a4,'on');
    
%     title(TITLE,'FontSize',12,'FontWeight','bold');
        hold(a4,'off');
        if sum(strcmp(mode,{'T1','T2'}))
        	xlabel(a4,'\tau [s]');ylabel(a4,'Normalized magnetization')
        elseif strcmp(mode,'D1')
        	xlabel(a4,'D1 [\mus]');ylabel(a4,'Normalized magnetization')
        elseif strcmp(mode,'SW')
        	xlabel(a4,'\nu [MHz]');ylabel(a4,'Normalized magnetization')
            xlim(a4,[min(FR) max(FR)])
        end
    end
    
    drawnow;
    
    assignin('base','INTall',INTall);
    assignin('base','S',S);
    assignin('base','Sabs',Sabs);
