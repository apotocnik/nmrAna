%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIT_PH_1                            %
% Fitting phase in sweep experiments  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fitpar,model]=FIT_PH_1(FR,INT,expPH,REF,initpar,SIMtype,Niter)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if strcmp('poly2',SIMtype)
        model = @poly2;
        fitpar = fminsearch(model,initpar,optimset('Display','iter','MaxIter',Niter));
%         fitpar = fminsearch(model,initpar,optimset('MaxIter',Niter));
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [sse,fitPH] = poly2(par)

        A = par(1);
        B = par(2);
        C = par(3);

        fitPH = mod(A+B*(FR-REF)+C*(FR-REF).^2+180,360)-180;
        sse = sum(INT.*(fitPH-expPH).^2)/sum(INT);
        %sse = sum((fitPH-expPH).^2);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
