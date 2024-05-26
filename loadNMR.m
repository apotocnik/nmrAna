%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loadNMR
%
% FileName is used only for a model, only TE temperature will be loaded
% if IND0 is empty, then automatically obtain IND
%
% Anton Potoènik, 26.2.2011, F5 @ IJS
% - initial version
% Andraž Krajnc, 26.7.2011, FMF 
% - changes for speed improvements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [IND TEs] = loadNMR(FileName, TE, IND0, TURN)
    
    [Path Name Ext] = fileparts(FileName);
    Path = [Path '\'];
%     P1mode = 0;
    [preName, tmpTE, tmpIND] = getNameParts(Name); % Get preName
    %if isempty(TE)
    %    TE = tmpTE;
    %end;
    % Check for other files and extract IND of desired TE
    d = dir([Path preName '*001' Ext]);
    
    if numel(d)>0 %%%%%%%%%%%%%%%%%%%%%%%%
    
        % Check for similar files with different temperatures
        TEs = zeros(1,numel(d));
        for i=1:numel(d)
            [tmpPath Name Ext] = fileparts(d(i).name);
            [preName, tmpTE, tmpIND] = getNameParts(Name);
            if isempty(tmpTE), continue; end; % wrong filename: 150K-G ending
            TEs(i) = tmpTE;
        end
        TEs = sort(TEs);
    
        % Check for index files of specific filename and temperature
        if isempty(TE)
            d = dir([Path preName '0*' Ext]);
        else
            d = dir([Path preName num2str(TE) 'K-0*' Ext]);
        end
        
        if isempty(d)
            msgbox('File not found!','error');
            error(['File not found: ' Path preName num2str(TE) 'K-0*' Ext]);
        end
        indN2 = numel(d);
        dstv = 1;
        while indN2>=99 && dstv<10 %???
            d2 = dir([Path preName num2str(TE) 'K-' num2str(dstv) '*' Ext]);
            indN2 = numel(d2); dstv = dstv+1;
            if indN2>0
                d = [d; d2];
            end
        end
        indN = numel(d);
        fnames = cell(1,indN);
        IND = zeros(1,indN);
        for i=1:indN
            [tmpPath Name Ext] = fileparts(d(i).name);
            [preName, tmpTE, tmpIND] = getNameParts(Name);
            IND(i) = tmpIND;
            fnames{i} = d(i).name;
        end
     
        if numel(d)==99 % ???

        end

        fnames_tmp = {}; % ???
        IND_tmp = [];
        for i=1:numel(IND0)
           if ismember(IND0(i),IND)
               [~,COL] = find(IND == IND0(i));
               fnames_tmp = [fnames_tmp, fnames{COL}];
               IND_tmp = [IND_tmp, IND0(i)];
           end
        end
        if ~isempty(IND_tmp)
            IND = IND_tmp;
            fnames = fnames_tmp;
            indN = numel(IND);
        end
    
    else % obvious we dont need index, so this would be classic spectrum
            d = dir([Path preName '*K' Ext]);
            TEs = zeros(1,numel(d));
            indN = 1;
            fnames = {[preName num2str(TE) 'K' Ext]};
            IND = 1;
        for i=1:numel(d)
            [tmpPath Name Ext] = fileparts(d(i).name);
            [preName, tmpTE, tmpIND] = getNameParts(Name);
            if isempty(tmpTE), continue; end; % wrong filename: 150K-G ending
            TEs(i) = tmpTE;
        end
        TEs = sort(TEs);
    end % end of ..
    
    
    FIDall = cell(1,indN); Tall = cell(1,indN); FR = zeros(1,indN); TAU = zeros(1,indN); D1 = zeros(1,indN);
    M1 = zeros(1,indN); M2 = zeros(1,indN); TRIMLEV = zeros(1,indN);
    %FIDall = cell(length(IND));

    for i=1:indN
        id = fopen([Path fnames{i}]);
        if id == -1
            msgbox([Path fnames{i} ' cannot be found !!!']);
            error([Path fnames{i} ' cannot be found !!!'])
            
%         else
%             gdisp(handles,fnames{i});
        end
        str = fgetl(id);
        while 1-strcmp(str,'[DATA]')
            str = fgetl(id);
            ind = strfind(str,'DW=');
            if ind==1
                DW = str2double(str(ind+3:length(str)));
            end
            ind = strfind(str,'TD=');
            if ind==1
                TD = str2double(str(ind+3:length(str)));
            end
            ind = strfind(str,'FR=');
            if ind==1
                FR(i) = str2double(str(ind+3:length(str)));
            end
            ind = strfind(str,'TAU=');
            if ind==1
                TAU(i) = str2double(str(ind+4:length(str)));
            end
            ind = strfind(str,'D1=');
            if ind==1
                D1(i) = str2double(str(ind+3:length(str)));
            end
            ind = strfind(str,'D2=');
            if ind==1
                D2(i) = str2double(str(ind+3:length(str)));
            end
            ind = strfind(str,'P1=');
            if ind==1
                D1(i) = str2double(str(ind+3:length(str)));
%                 P1mode = 1;
            end
            ind = strfind(str,'_M1=');
            if ind==1
                M1(i) = str2double(str(ind+4:length(str)));
            end
            ind = strfind(str,'_M2=');
            if ind==1
                M2(i) = str2double(str(ind+4:length(str)));
            end
            ind = strfind(str,'_TRIMLEV=');
            if ind==1
                TRIMLEV(i) = str2double(str(ind+9:length(str)));
            end
        end
        FID = (textscan(id,'%f %f', [2,inf]));
        %FID = (fscanf(id, '%g %g', [2,inf]))';
        fclose(id);

        Tall{i} = (0:DW:(TD-1)*DW)*1e6';
        
        splitName = regexp(fnames{i},'_','split');
        
        if TURN 
            FIDall{i} = cell2mat(FID(1)) + complex(0,1)*cell2mat(FID(2));
        else
            FIDall{i} = cell2mat(FID(1)) - complex(0,1)*cell2mat(FID(2));
        end
        % za 380 MHz je FIDall potrebno kompleksno konjugirat!!! (more bit
        % +, pri 400 MHz pa -), TODO TODO !!!
        
    end

    if ~ismember('M1',who)
        M1 = 0;
        M2 = 0;
    end
 
    assignin('base','FIDall',FIDall);
    assignin('base','Tall',Tall);
%     assignin('base','FIDall0',FIDall);
%     assignin('base','Tall0',Tall);
    assignin('base','DW',DW);
    assignin('base','M1',M1);
    assignin('base','M2',M2);
    assignin('base','TRIMLEV',TRIMLEV);
    assignin('base','TD',TD);
    assignin('base','FR',FR);
    assignin('base','TAU',TAU);
    assignin('base','D1',D1);
    assignin('base','D2',D2);
%     if P1mode
%         assignin('base','P1',P1);
%     end
    assignin('base','IND',IND);
    assignin('base','TEM',TE);
    assignin('base','TEs',TEs);
    