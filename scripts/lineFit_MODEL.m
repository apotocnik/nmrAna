%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIT_Na-LSX_MODEL_1                  %
% Fitting of Na12/Na-LSX spectra      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fitpar,model]=lineFit_MODEL(f,fEXP,expSPC,MODEL,Niter,NS,peakType,allinitpar,indpar,posAmp)
    SHOWFIT = evalin('base','SHOWFIT');
    if SHOWFIT
        figure(100); clf;
        inum = 1;
        sseall = [];
    end
    
    if ismember(MODEL,[1:1])
        model = @SIMPLE;
        initpar = allinitpar(indpar);
        fitpar = fminsearch(model,initpar,optimset('MaxIter',Niter));
    elseif ismember(MODEL,[2:2])
        model = @SIMPLE;
        initpar = allinitpar(indpar);
        numparams = 3*ismember(peakType,[1,2])+4*ismember(peakType,[3]);
        Ucon = evalin('base','MODEL2_Ucon'); Lcon = evalin('base','MODEL2_Lcon');
        ratio = evalin('base','MODEL2_ratio');
        relativeCon = evalin('base','MODEL2_RelativeCon');
        if relativeCon
            peaknum = 1;
            for h=1:numel(peakType)
                if ismember(peakType(h),[1,2])
                    Lcon(peaknum:peaknum+2) = [0,1,0].*allinitpar(peaknum:peaknum+2)+[0,1,0].*Lcon(peaknum:peaknum+2)+allinitpar(peaknum:peaknum+2).*[1,0,1].*Lcon(peaknum:peaknum+2);
                    Ucon(peaknum:peaknum+2) = [0,1,0].*allinitpar(peaknum:peaknum+2)+[0,1,0].*Ucon(peaknum:peaknum+2)+allinitpar(peaknum:peaknum+2).*[1,0,1].*Ucon(peaknum:peaknum+2);
                    peaknum = peaknum+3;
                elseif ismember(peakType(h),[3])
                    peaknum = peaknum+4;
                end
            end
        end
        ratio = ratio/sum(ratio);
        if isempty(ratio), nonlin = []; else nonlin = @mycon; end;
            
        %Acon = [eye(numel(initpar));-eye(numel(initpar))]; bcon = [Ucon,-Lcon];
        fitpar = fmincon(model,initpar,[],[],[],[],Lcon,Ucon,nonlin,optimset('MaxIter',Niter,'Algorithm','interior-point')); %sqp, interior-point, active-set
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % sse - chi^2
    % fitSPC - 
    function [sse,fitSPC,simSPC,simSPCall] = SIMPLE(par)

        fitSPC = []; simSPC = [];
        allpar = allinitpar;
        allpar(indpar) = par;
 
        % simSPCall = allpar(7*NS+1); this used to be offset!
        simSPCall = 0;
        inum = 0;
        for i=1:NS
            % 1 = lorentz, 2 = gauss, 3 = poly
            if peakType(i)==1
                    A = allpar(1+inum);
                    if (posAmp && A<0), A = inf; end;
                    pos = allpar(2+inum);
                    width = allpar(3+inum);
                    inum = inum+3;
                simSPC{i} = A*lorentzian(f,pos,width);
            elseif peakType(i)==2
                    A = allpar(1+inum);
                    if (posAmp && A<0), A = inf; end;
                    pos = allpar(2+inum);
                    width = allpar(3+inum);
                    inum = inum+3;
                simSPC{i} = A*gaussian(f,pos,width);
            elseif peakType(i)==3
                    A = allpar(1+inum);
                    B = allpar(2+inum);
                    C = allpar(3+inum);
                    x0 = allpar(4+inum);
                    inum = inum+4;
                simSPC{i} = A.*(f-x0).^2+B.*(f-x0)+C;
            end
            simSPCall = simSPCall+simSPC{i};
        end

        sse = 0;
        fitSPC = interp1(f,simSPCall,fEXP);
        inds= find(isnan(fitSPC)); %
        fitSPC(inds)= zeros(size(inds));
        sse = sse + sum((fitSPC-expSPC).^2);
        
       if SHOWFIT
           figure(100);
           subplot(2,1,1);
           %hold on;
           sseall = [sseall, sse];
           plot(sseall,'o');
           %inum = inum+1;
            title('Iterations')
           subplot(2,1,2);
           cla;
           hold on;
            plot(fEXP,expSPC,'Color','r');
            plot(fEXP,fitSPC,'Color','k');
           hold off;
           title('Spectrum')
           drawnow();
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [c,ceq] = mycon(par)
%c = ...     % Compute nonlinear inequalities at x.
%ceq = ...   % Compute nonlinear equalities at x.
    c = 0;
    Aall = [];
    for i=1:NS
       zparam = sum(numparams(1:i));
       Aall = [Aall par(zparam-2)*par(zparam)];
    end
    ceq = 0; Aall_sum = sum(Aall);
    for i=1:numel(Aall)
        ceq = ceq + abs(Aall(i)/Aall_sum-ratio(i));
    end
end

end
