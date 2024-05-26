%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% txt = val2txt(val,err0,dig)
%
% transform engeneering's type values withh error to text
%
% 0.00423 +/- 0.00022 => 4.23(22) m
%
% Anton Potoènik, 1.1.2013, F5 @ IJS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function txt = val2txt(val,err0,dig)
    isErr = false;
    isNeg = false;
    fin = '';
    
    if nargin > 1
        if ~isempty(err0)
            if err0 > 0
                isErr = true;
            end
        end
    end

    if nargin < 3
        dig = 4;
    end
    
    if val < 0
        isNeg = true;
        val = -val;
    end
    
    if val ~= 0
        t = log(val)/log(10);
    else
        t = 0;
    end

    if t < 0
        a = abs(ceil(t)); % poišèi enice
        b = t+a; % poišèi dele
        if abs(b) > 0.302 % <0.5
            c = a+1;
        else
            c = a;
        end
        units = ceil(c/3);
        switch units
            case 0
                fin = '';
            case 1
                fin = 'm';
            case 2
                fin = 'u';
            case 3
                fin = 'n';
            case 4
                fin = 'p';            
        end
        val = val*realpow(10,3*units);
        d = dig - (3*units - a);
        if isErr
            err = round(err0*realpow(10,3*units+d));
            while err == 0 
                d = d + 1;
                err = round(err0*realpow(10,3*units+d));
            end
        end

    else
        a = abs(floor(t)); % poišèi enice
        b = t-a; % poišèi dele
        if abs(b) > 0.6985 % <50
            c = a+1;
        else
            c = a;
        end
        units = floor(c/3);
        switch units
            case 0
                fin = '';
            case 1
                fin = 'k';
            case 2
                fin = 'M';
            case 3
                fin = 'G';
            case 4
                fin = 'T';            
        end
        val = val/realpow(10,3*units);
        d = dig - ((a+1)-3*units);
        if isErr
            err = round(err0/realpow(10,3*units-d));
            while err == 0 
                d = d + 1;
                err = round(err0/realpow(10,3*units-d));
            end
        end
    end
    
    if isNeg 
        val = -val;
    end
    
    if err0 == 0
       isErr = true;
       err = 0;
    end
    
    form = ['%10.' num2str(d) 'f'];
    if isErr
        txt = [num2str(val,form) '(' num2str(err) ')' fin];
    else
        txt = [num2str(val,form) '' fin];
    end
%     disp(txt)
