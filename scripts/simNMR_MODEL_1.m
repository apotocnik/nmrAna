%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIT_Na-LSX_MODEL_1                  %
% Fitting of Na12/Na-LSX spectra      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fitpar,model]=simNMR_MODEL_1(f,fEXP,expSPC,freqL,fi,Rfi,tmpREF,Qmom,S,LBtype,MODEL,Niter,allinitpar,indpar,GAUSS,NS);

    e0 = 1.60217733e-19;
    h = 6.6260755e-34;
    
    if ismember(MODEL,[1:9])
        model = @S1_3Q12K1;
        aksfi = fi;
        initpar = allinitpar(indpar);
        fitpar = fminsearch(model,initpar,optimset('Display','iter','MaxIter',Niter));
        % fmincon
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [sse,fitSPC,simSPC,simSPCall] = S1_3Q12K1(par)

        fitSPC = []; simSPC = [];
        allpar = allinitpar;
        allpar(indpar) = par;
        
        simSPCall = allpar(7*NS+1);

        for i=1:NS
            Vzz = allpar(1+(i-1)*7);
            eta = allpar(2+(i-1)*7);
            iso = allpar(3+(i-1)*7)-tmpREF;
            aksq = allpar(4+(i-1)*7);
            akseta = allpar(5+(i-1)*7);
            sig = allpar(6+(i-1)*7);
            A = allpar(7+(i-1)*7);
            
            if A<0, A = inf; end;

            freqQ = 3*e0*Qmom{i}*Vzz/(h*(2*S{i}*(2*S{i}-1)));

            CENTefg =(1/8*(S{i}*(S{i}+1)-3/4))*freqQ.^2./freqL{i}/12.*(6*(1-fi{i}(2,:)).*(1-9*fi{i}(2,:))-4*eta.*fi{i}(1,:).*(1-fi{i}(2,:)).*(9*fi{i}(2,:)+1)+eta.^2.*(-16/3+8*fi{i}(2,:)+6*fi{i}(1,:).^2.*(1-fi{i}(2,:)).^2));
            CENTaks =freqL{i}*(iso+aksq./2.*(3*aksfi{i}(2,:)-1+akseta.*(1-aksfi{i}(2,:)).*aksfi{i}(1,:)));
            HIST = hist(CENTefg+CENTaks,f)';

            SATefg =freqQ./2.*(3*fi{i}(2,:)-1+eta.*(1-fi{i}(2,:)).*fi{i}(1,:));
            HISTr = hist(CENTaks+SATefg,f)';
            HISTl = hist(CENTaks-SATefg,f)';

            %SPEK = 4*HIST; % only central transition
            SPEK = 3*HISTl+4*HIST+3*HISTr; % +satellites

            SPEK(1) = SPEK(2);
            SPEK(length(SPEK)) = SPEK(length(SPEK)-1);

            if LBtype==1
                LB = fft(1./(1+(f/sig).^2));
            end
            if LBtype==2
                LB = fft(exp(-(f/sig).^2/2));
            end

            simSPC{i} = real(ifftshift(ifft(fft(SPEK).*LB)));
            simSPC{i} = A*simSPC{i};
            simSPCall = simSPCall+simSPC{i};
        end

        sse = 0;
        fitSPC = interp1(f,simSPCall,fEXP);
        sse = sse + sum((fitSPC-expSPC).^2);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
