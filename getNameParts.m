%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getNameParts
%
% Extracts name, temperature and ind from filename: AM2_4_10kN_300K-023
% Name = 'AM2_4_10kN', TE = 300, IND = 23;
%
% Anton Potoènik, 26.2.2011, F5 @ IJS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Name, TE, IND] = getNameParts(Name)
    % Nothing is behind IND
    % K is behind temperature (dont understand: 100p5K only 100K)
    
    TE = [];
    IND = [];
    
    % test if really nothing is behind IND 
    if isstrprop(Name(end),'digit')
        idigit = isstrprop(Name,'digit');
        res = find(idigit==0);
        idxIND = res(end)+1;
        IND = str2double(Name(idxIND:end));
        Name(idxIND:end) = [];
    end

    idxK = strfind(Name, 'K');
    if isempty(idxK), return; end;
    
    idxK = idxK(end);
    Name(idxK:end) = [];
    idigit = isstrprop(Name,'digit');
    res = find(idigit==0);
    if strcmp(Name(res(end)),'.') == 1
        res(end) = res(end-1);
    end
    idxTEen = idxK-1;
    idxTEst = res(end)+1;
    TE = str2double(Name(idxTEst:idxTEen));
    Name(idxTEst:end)=[];
%     else
%         TE = [];
%         Name((idxIND-1):end)=[];
%     end;
    
    
    