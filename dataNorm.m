function [dataOUT factor] = dataNorm(dataIN, Ntype, params)
% for input parameters you get dataIn(:,1:2), Ntype, [additional params]
% output is quite obvious
%
%   0 - do not normalize
%   1 - divide by max and multiply with params(3)
%   2 - within a range: divide by max times params(3)
%   3 - within a range: divide by area and multiply with params(3)
%   4 - multiply with temperature and param(3)
%
%   range = [params(1) params(2)]
%
% Andraž Kranjc @IJS, Anton Potoènik @ IJS 21.11.2012
% -------------------------------------------------------------------------

    if nargin < 3
        params = [-inf inf 1];
    end

    Amp = params(3);
    Fmin = params(1);
    Fmax = params(2);
    X = dataIN(:,1);
    Y = dataIN(:,2);
    ind = find(Y>=Fmin & Y<=Fmax);
  
    switch Ntype
        case 0
          factor = Amp;
            
        case 1
          factor = Amp/max(Y);    

        case 2
          factor = Amp/max(Y(ind));

        case 3
          area = trapz(X(ind),Y(ind));
          factor = Amp/area;

        case 4
          TEM = evalin('base','TEM');
          factor = Amp*TEM;

        otherwise
          factor = 1;
    end
    
    dataOUT = [X Y*factor];
end

