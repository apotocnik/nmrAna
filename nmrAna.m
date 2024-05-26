function varargout = nmrAna(varargin)
% NMRANA MATLAB code for nmrAna.fig
%      NMRANA, by itself, creates a new NMRANA or raises the existing
%      singleton*.
%
%      H = NMRANA returns the handle to a new NMRANA or the handle to
%      the existing singleton*.
%
%      NMRANA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NMRANA.M with the given input arguments.
%
%      NMRANA('Property','Value',...) creates a new NMRANA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before nmrAna_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to nmrAna_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help nmrAna

% Last Modified by GUIDE v2.5 13-Apr-2014 09:26:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @nmrAna_OpeningFcn, ...
                   'gui_OutputFcn',  @nmrAna_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before nmrAna is made visible.
function nmrAna_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to nmrAna (see VARARGIN)

% Choose default command line output for nmrAna
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% simply ... needed
PATHname = regexprep([mfilename('fullpath') '.m'], 'nmrAna.m', '');
assignin('base','PATHname',PATHname);
assignin('base','FFT',1);
assignin('base','autoRun',0);

addlistener(handles.dir1,'String','PostSet',@(s,e) dir1CB(handles));
addlistener(handles.dir2,'String','PostSet',@(s,e) dir2CB(handles));
addlistener(handles.dir3,'String','PostSet',@(s,e) dir3CB(handles));
addlistener(handles.dir4,'String','PostSet',@(s,e) dir4CB(handles));
addlistener(handles.dir5,'String','PostSet',@(s,e) dir5CB(handles));
addlistener(handles.edtFN,'String','PostSet',@(s,e) checkDir(handles));

list = nmrCustom(0,hObject, eventdata,handles);
set(handles.popCustom,'String',list);
% UIWAIT makes nmrAna wait for user response (see UIRESUME)
% uiwait(handles.figure1);



function checkDir(handles)
%display('debug');
fp = get(handles.edtFN,'String');
fp1 = regexp(fp,'\','split');

%check n-th folder
for i=1:5
    replstr = [fp1{end}];
    for j=1:(5-i)
       replstr = [fp1{end-j} '\' replstr];
    end
    
    fp3 = strrep(fp,replstr,'');
    if ~exist(fp3)
        %daš prvega najdenega v 5. folder
        replstr = [fp1{end-(6-i)} '\' replstr];
        fp3 = strrep(fp,replstr,'');
        a = dir(fp3);
        for j=3:numel(a)
           if a(j).isdir
               fp1{end-(6-i)} = a(j).name;
                   fp = [fp1{1}];
                   for k=2:numel(fp1)
                       fp = [fp '\' fp1{k}];
                   end
               break;
           end
           if j==numel(a)
               msgbox('error! wrong folders.');
               return;
           end
        end
    end
end



%check file
fp3 = strrep(fp,fp1{end},'');
if ~exist(fp3)
    
end

%
    selected = cell(1,5);
    selected{1} = fp1{end-1};
    selected{2} = fp1{end-2};
    selected{3} = fp1{end-3};
    selected{4} = fp1{end-4};
    selected{5} = fp1{end-5};
    folders = cell(1,5);
    
    for i=1:5
        if i==1
            fp3 = strrep(fp,[fp1{end-1} '\' fp1{end}],'');
        elseif i==2
            fp3 = strrep(fp,[fp1{end-2} '\' fp1{end-1} '\' fp1{end}],'');
        elseif i==3
            fp3 = strrep(fp,[fp1{end-3} '\' fp1{end-2} '\' fp1{end-1} '\' fp1{end}],'');
        elseif i==4
            fp3 = strrep(fp,[fp1{end-4} '\' fp1{end-3} '\' fp1{end-2} '\' fp1{end-1} '\' fp1{end}],'');
        elseif i==5
            fp3 = strrep(fp,[fp1{end-5} '\' fp1{end-4} '\' fp1{end-3} '\' fp1{end-2} '\' fp1{end-1} '\' fp1{end}],'');
        end
    A = dir(fp3);
    folders{i} = {};
    for j=3:numel(A)
       if A(j).isdir
           folders{i} = [folders{i} A(j).name];
       end
    end
    end
    
    fp1 = regexp(fp,'\','split');
    if ~exist(fp)
        fp3 = strrep(fp,[fp1{end}],'');
            A = dir(fp3);
            for j=3:numel(A)
               if ~A(j).isdir && numel(regexpi(A(j).name,'.dat$'))>0 && numel(regexpi(A(j).name,'-G.dat$'))==0
                   fp1{end} = A(j).name;
                       fp = [fp1{1}];
                       for i=2:numel(fp1)
                           fp = [fp '\' fp1{i}];
                       end
                   break;
               end
               if j==numel(A)
                   msgbox('Error! No files with extension .dat!');
                   return;
               end
            end
    end
        fp3 = strrep(fp,['\' fp1{end}],'');
        cd(fp3);
        % Extract Temperature
        [~, Name Ext] = fileparts(fp1{end});
        [~, TE IND] = getNameParts(Name);
        set(handles.edtTE,'String',['[' num2str(TE) ']']);
        set(handles.edtFN,'String',fp);
    
assignin('base','folders',folders);
assignin('base','selected',selected);

%set(handles.edtFN,'String',fp);
%drawnow
updateDirs(handles);


function updateFileList(handles)
fp2 = get(handles.edtFN,'String');
fpindFN = regexp(fp2,'\','start');
fp = fp2(1:fpindFN(end));
fname2 = fp2(fpindFN(end)+1:end);
fpind = regexp(fname2,'_','start');
fname = lower(fname2(1:fpind(end)-1));
fnameind = 0;

allfiles = {};
dirline = {'*K-001.dat','*K.dat'};

for j = 1:numel(dirline)
    files = dir([fp dirline{j}]);
    fplist = cell(1,numel(files));
    fplist2 = cell(1,numel(files));
    foffset = numel(allfiles);
    
    for i = 1:numel(fplist)
       namestr = files(i).name;
       fpind = regexp(namestr,'_','start');
       fplist2{i} = namestr;
       fplist{i} = lower(namestr(1:fpind(end)-1));
    end
    [~,ia,~] = unique(fplist);
    C2 = [fplist(ia)];
    findexC = strcmp(C2, fname);
    findex = find(findexC==1);
    if ~isempty(findex)
        fnameind = foffset+findex;
    end
    allfiles = [allfiles, fplist2(ia)];
end
assignin('base','allfiles',allfiles);
if isempty(allfiles)
    allfiles = {''};
end
set(handles.dirFiles,'String',allfiles);
%now we need to select the right file
if fnameind>0
    set(handles.dirFiles,'Value',fnameind);
else
    set(handles.dirFiles,'Value',1);
    if ~strcmp(allfiles{1},'')
        set(handles.edtFN,'String',[fp allfiles{1}]);
    end
end

function updateDirs(handles)

folders = evalin('base','folders');
selected = evalin('base','selected');

for i=1:5
   folders{i} = lower(folders{i});  
   selected(i) = lower(selected(i));
end

str5 = {}; str4 = {}; str3 = {}; str2 = {}; str1 = {};
str5 = lower([str5,get(handles.dir5,'String')]);
old5 = lower(str5{get(handles.dir5,'Value')});
str4 = lower([str4,get(handles.dir4,'String')]);
old4 = lower(str4{get(handles.dir4,'Value')});
str3 = lower([str3,get(handles.dir3,'String')]);
old3 = lower(str3{get(handles.dir3,'Value')});
str2 = lower([str2, get(handles.dir2,'String')]);
old2 = lower(str2{get(handles.dir2,'Value')});
str1 = lower([str1, get(handles.dir1,'String')]);
old1 = lower(str1{get(handles.dir1,'Value')});

num1 = 1; num2 = 1; num3 = 1; num4 = 1; num5 = 1;
if strcmpi(old1,selected{1})
    num1 = find(ismember(folders{1}, old1)==1);
else
    num1 = find(ismember(folders{1}, selected{1})==1);
end
if strcmpi(old2,selected{2})
    num2 = find(ismember(folders{2}, old2)==1);
else
    num2 = find(ismember(folders{2}, selected{2})==1);    
end
if strcmpi(old3,selected{3})
    num3 = find(ismember(folders{3}, old3)==1);
else
    num3 = find(ismember(folders{3}, selected{3})==1);    
end
if strcmpi(old4,selected{4})
    num4 = find(ismember(folders{4}, old4)==1);
else
    num4 = find(ismember(folders{4}, selected{4})==1);    
end
if strcmpi(old5,selected{5})
    num5 = find(ismember(folders{5}, old5)==1);
else
    num5 = find(ismember(folders{5}, selected{5})==1);    
end
folders = evalin('base','folders');

set(handles.dir1,'Value',num1);
set(handles.dir2,'Value',num2);
set(handles.dir3,'Value',num3);
set(handles.dir4,'Value',num4);
set(handles.dir5,'Value',num5);
set(handles.dir1,'String',folders{1});
set(handles.dir2,'String',folders{2});
set(handles.dir3,'String',folders{3});
set(handles.dir4,'String',folders{4});
set(handles.dir5,'String',folders{5});
updateFileList(handles);

% --- Outputs from this function are returned to the command line.
function varargout = nmrAna_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function edtFN_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

    
function edtFFTparam_Callback(hObject, eventdata, handles)
    reloadFFT(handles);

% --- Executes during object creation, after setting all properties.
function edtFFTparam_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edtAnaparam_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function edtAnaparam_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edtTE_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function edtTE_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edtIND_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function edtIND_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in chkINDauto.
function chkINDauto_Callback(hObject, eventdata, handles)


function edit7_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edtAnalyse_Callback(hObject, eventdata, handles)
    reloadSettings(handles);


% --- Executes during object creation, after setting all properties.
function edtAnalyse_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when uipanel3 is resized.
function uipanel3_ResizeFcn(hObject, eventdata, handles)


function edtFit_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function edtFit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in butBrowse.
function butBrowse_Callback(hObject, eventdata, handles)

FilePath = get(handles.edtFN,'String');
Fp=regexp(FilePath,'\','split');
FileName = Fp{end};
Fp = regexp(FilePath,FileName,'split');
PathName = Fp{1};

if isdir(PathName), cd(PathName); end

    [FileName,PathName] = uigetfile({'*.dat';'*.*'});
    if ~isequal(FileName, 0)
        cd(PathName)
        set(handles.edtFN,'String',[PathName FileName]);
        % Extract Temperature
        [~, Name Ext] = fileparts(FileName);
        [~, TE IND] = getNameParts(Name);
        set(handles.edtTE,'String',['[' num2str(TE) ']']);
    end

    

% --- Executes on button press in btnLoad.
function btnLoad_Callback(hObject, eventdata, handles)

    PLOT = get(handles.chkPLOT,'Value');
    assignin('base','PLOT',PLOT);
    
    autoIND = get(handles.chkINDauto,'Value');
    IND0 = [];
    if autoIND == 0
        IND0 = str2num(get(handles.edtIND,'String'));
    end
    FileName = get(handles.edtFN,'String');
    TE = str2num(get(handles.edtTE,'String'));
    
    if evalin('base','exist(''TURN'',''var'')')
        TURN = evalin('base', 'TURN');
    else
        TURN = 0;
    end
    
    [IND TEs] = loadNMR(FileName, TE, IND0, TURN); % takes only 0.2s
    
    strTEs = '';
    for i=1:numel(TEs)
        strTEs = [strTEs num2str(TEs(i)) ' ' ];
    end
    set(handles.txtTemps,'String',strTEs);
    
    if autoIND
%         set(handles.edtIND,'String',['[' num2str(IND(end)) ':-1:' num2str(IND(1)) ']']);
        set(handles.edtIND,'String',[num2str(IND(1)) ':' num2str(IND(end))]);
    end
    
    [~, Name Ext] = fileparts(FileName);
    if strfind(Name,'T1')
        mode = 'T1';
    elseif strfind(Name,'T2')
        mode = 'T2';
    elseif strfind(Name,'D1')
        mode = 'D1';
    elseif strfind(Name,'D2')
        mode = 'D2';
    elseif strfind(Name,'P1')
        mode = 'D1';
    elseif strfind(Name,'SW')
        mode = 'SW';
%         reloadSettings(handles);
    else
        if numel(IND) > 1
            r = inputdlg('What type of measurement is this? [D1, T1, T2, SW]','Mode',1,{'D1'});
            switch upper(r{1})
                case 'SW'
                    mode = 'SW';
                case 'T1'
                    mode = 'T1';
                case 'T2'
                    mode = 'T2';
                case 'D1'
                    mode = 'D1';     
                otherwise 
                    return;
            end
        else
            mode = 'SPC';  
        end
    end
    assignin('base','mode',mode);
    
    evalin('base','clear S'); % Remove weigths since new files can be loaded

%     reloadSettings(handles);
%     reloadFFT(handles); % read FFT parameters
% 
% %     if evalin('base', 'autoPH') == 4 % For phase fitting
% %         PLOT = evalin('base', 'PLOT'); %save plot
% %         assignin('base', 'autoPH',2);
% %         assignin('base', 'PLOT',0);
% %         fftNMR(handles); 
% %         assignin('base', 'autoPH',4);
% %         assignin('base', 'PLOT',PLOT);
% %     end
%     fftNMR(handles); % 0.15s
% %     if strcmp(mode,'SW')
% %         integNMR(handles);
% %         % do nothing here
% % %         showFID(1,1,1,handles);
% % %         FR = evalin('base', 'FR');
% % %         disp(['FR step is ' num2str(1e3*abs(FR(2)-FR(1))) ' kHz']);
% %     else
%     if strcmp(mode,'SPC')
%         butSpecter_Callback(hObject, eventdata, handles)
%     else
%         integNMR(handles);
%     end

    gdisp(handles,'Load complete');
    
    % Continue with FFT
    butFFT_Callback(hObject, eventdata, handles)
    
    
    
% --- Executes on button press in butFFT.
function butFFT_Callback(hObject, eventdata, handles)
    mode = evalin('base','mode');
    reloadFFT(handles); % read FFT parameters
    reloadSettings(handles); % read Settings parameters

    r=fftNMR(handles); % 0.15s
    if r < 0, return; end;

    if strcmp(mode,'SPC')
        butSpecter_Callback(hObject, eventdata, handles)
    else
        integNMR(handles);
    end

    
% --- Executes on button press in butShowIND.
function butShowIND_Callback(hObject, eventdata, handles)
    reloadFFT(handles); % read FFT parameters
    reloadSettings(handles); % read Settings parameters
    
    IND = evalin('base', 'IND');
    answ = (inputdlg({'0=FID or 1=SPC:','Which IND to show?','Show Line Broadening?','Show Phase Correction?'},'Show IND',1,{'0','1','0','0'}));
    if isempty(answ), return; end;
    isSPC = str2double(answ{1});
    ii = str2num(answ{2});
    isLB = str2double(answ{3});
    isPH = str2double(answ{4});
    
    if isSPC==0
        if numel(ii) > 1
            msgbox('Use only one IND!','Error');
            return
        else
            if intersect(ii,IND)
                showFID(ii,isLB,isPH,handles)
            end
        end
    else
        Fall = evalin('base','Fall');
        SPCall = evalin('base','SPCall');
        Dy = evalin('base','Dy');
        range = evalin('base','range');
        
        if numel(ii) == 1
           a4 = handles.axeRes;
           Y = real(SPCall{ii});
           Y = Y/max(Y);
           plot(a4,Fall{ii},Y);
           text(min(range),0.5,['IND: ' num2str(ii)],'Parent',a4);
           xlim(a4,range);
           grid(a4,'on');
        else
            data = [];
           a3 = handles.axeSPC;
           for j=ii
               Y = real(SPCall{j});
               Y = Y/max(Y);
               plot(a3,Fall{j},Y+Dy*(j-1)); hold(a3,'on');
               text(min(range),Dy*(j-1)+0.5,[num2str(j)],'Parent',a3);
               data = [data Fall{j} Y];
           end
           xlim(a3,range);
           hold(a3,'off');
           grid(a3,'on');
           num2clip(data)
        end
        
    end
    
    function showFID(i,isLB,isPH,handles)
        TEM = evalin('base', 'TEM');
        FR = evalin('base', 'FR');
        LB = evalin('base', 'LB');
        Tall = evalin('base', 'Tall');
        FIDall = evalin('base', 'FIDall');
        SHLall = evalin('base', 'SHLall'); 
        phiall = evalin('base', 'phiall'); 

        x = Tall{i}';
        if isPH
        	FIDall{i} = FIDall{i}*exp(-complex(0,1)*phiall(i,3));
        end
        
        if isLB
           FIDall{i} = FIDall{i}.*exp(-x.*LB/1e6); 
        end
 
        im = imag(FIDall{i});
        re = real(FIDall{i});
        ab = abs(FIDall{i});

        a3 = handles.axeRes;
        plot(a3,x,re,'b',x,im,'r',x,ab,'k','Linewidth',2)
        hold(a3,'on');
        Ylimits = get(a3,'YLim');
        plot(a3,[Tall{i}(SHLall(i,4)) Tall{i}(SHLall(i,4))],[Ylimits(1) Ylimits(2)])
        hold(a3,'off');
        xlabel(a3,'t [\mus]'); ylabel(a3,'FID');
        text(0.95,0.9,['\nu = ',num2str(FR(i)),' MHz, ','{\itT} = ',num2str(TEM),' K'],'FontWeight','bold','HorizontalAl','Right','Units','Norm','Parent',a3)
        xlim(a3,[0,max(Tall{i})]);
        grid(a3,'on');
        

    
% function LSparams(handles) % update initial fillting parameters
% indFitpar = evalin('base','indFitpar');
% allinitpar = evalin('base','allinitpar');
% allfitpar = evalin('base','allfitpar');
% 
% for i=1:numel(allinitpar)
%    if ismember(i,indFitpar)
%       allinitpar(i) = allfitpar(i); 
%    end
% end
% assignin('base','allinitpar',allinitpar);


% --- Executes on button press in butAnalyse.
function butAnalyse_Callback(hObject, eventdata, handles)
    reloadSettings(handles);
    TEM = evalin('base','TEM');
    % integNMR(handles);
    fitNMR(handles);
%     saveFIT=evalin('base','saveFIT');
%     LSname=evalin('base','LSname');
%     if saveFIT
%         saveImageWithinGUI(handles,handles.axeRes,[LSname '_' num2str(TEM) 'K']);
%     end
    

% --- Executes on button press in butSIMtype.
function butSIMtype_Callback(hObject, eventdata, handles)
    edit fitNMRlib


% --- Executes on button press in butScan.
function butScan_Callback(hObject, eventdata, handles)
    range = evalin('base','range');
    TEM = evalin('base','TEM');
    isPPM = evalin('base','isPPM');
    NORM = evalin('base','NORM');
    nParams = evalin('base','nParams');
    PLOT0 = evalin('base','PLOT');
%     PATHname = evalin('base','PATHname');
%     mode = evalin('base','mode');
    
    if evalin('base','exist(''LSpar'',''var'')')
        LSpar = evalin('base', 'LSpar');
    else
        LSpar = cell(1,5);
        LSpar{1} = 0;
    end
%     if evalin('base','exist(''LSname'',''var'')')
%         LSname = evalin('base', 'LSname');
%     else
%         LSname = 'noname';
%     end
    
    if LSpar{1} == 0
        dialog = inputdlg({'Plot parameters [ind]:','Frequency interval [kHz]:','Frequency step [kHz]:','Add plot [sorted ind]:','Spc normalize [1/0]:','Inversion Recovery [1/0]:'},'Line Scan',1,{'[2 5]','3','1','1','1','1'});
        if isempty(dialog), return; end;

        indpar = str2num(dialog{1});
        dF = str2double(dialog{2});
        step =  str2double(dialog{3});
        plotind = str2double(dialog{4});
        isNorm = str2double(dialog{5});
        isIR = str2double(dialog{6});
    else
        indpar = LSpar{2};
        dF = LSpar{3};
        step =  LSpar{4};
        plotind = LSpar{5};
        isNorm = LSpar{6};
        isIR = LSpar{7};
    end

    Fall = evalin('base','Fall');
    SPCall = evalin('base','SPCall');
    F0 = Fall{1}(1);
    N = round((Fall{1}(end)-F0)/step);
    Fscan = zeros(N,1);
    INTscan = zeros(N,numel(indpar));
    
    if plotind > 0
        spc = [Fall{plotind},real(SPCall{plotind})];
    end
    if isNorm && plotind > 0
        [spc norm] = dataNorm(spc,NORM,nParams);
    else
        norm = 1;
    end
    
    if isNorm && isIR % Inversion Recovery gives 2x the intensity
        norm = norm/2;
    end
    
    sses = zeros(1,N);
    for i=1:N
        F1 = F0 + dF;
        LIM = [F0 F1];
        F0 = F0 + step;
        assignin('base','LIM',LIM);
        integNMR(handles);
        fitNMR(handles);
        sses(i) = evalin('base','sse');
        Fscan(i) = mean(LIM);
        allfitpar = evalin('base','allfitpar');
        allfiterr = evalin('base','allfiterr');
        INTscan(i,:) = allfitpar(indpar);
        INTscanErr(i,:) = allfiterr(indpar);
        gdisp(handles,[num2str(i) '/' num2str(N) '   delta^2 = ' num2str(round(1000*sses(i))/1000)]);
        drawnow();
    end
    sse = mean(sses);
    gdisp(handles,['Average delta^2 (' num2str(TEM) ') = ' num2str(round(1000*sse)/1000)]);
    assignin('base','sse',sse);
    
    PLOT = 1;
    if PLOT
        a4 = handles.axeRes;
        LEG = {};
        COLOR = {'r','b','g','m','k','b','r','g','m'};
        INTallSpc = zeros(size(Fscan));
        for i=1:numel(indpar)
            
            % remove bad points
            [ind] = find(INTscanErr(:,i) > max(INTscan(:,i)));
            INTscan(ind,:)=0; INTscanErr(ind,:)=0;
            
            INTallSpc = INTallSpc + INTscan(:,i);
            
%             plot(a4,Fscan,INTscan(:,i)*norm,'k','Linewidth',1.5,'Color',COLOR{i})
            errorbar(a4,Fscan,INTscan(:,i),INTscanErr(:,i),'k','Linewidth',1.5,'Color',COLOR{i})
            hold(a4,'on');
            LEG{i} = ['indpar = ' num2str(indpar(i))];
        end

        if plotind > 0
            plot(a4,spc(:,1),spc(:,2),'k','Linewidth',1.5)
            LEG{i+1} = ['plotind = ' num2str(plotind)];
            plot(a4,Fscan,INTallSpc,'k:')
        end

        hold(a4,'off'); grid(a4,'on');
        legend(a4,LEG,2,'Location','NorthEast');
        xlabel(a4,'\nu [kHz]','FontSize',12); 
        xlim(a4,range); 
    end
   
    data={}; dataErr={}; TEs = {};
    for i=1:numel(indpar)
        if isPPM
            data{i} = [toPPM(Fscan) INTscan(:,i)*norm];
        else
            data{i} = [Fscan INTscan(:,i)*norm];
        end
        dataErr{i} = INTscanErr(:,i)*norm;
        TEs{i} = TEM;
    end
    if isPPM
        Fscan = toPPM(Fscan);
    end
    nmrFitData.data = data;
    nmrFitData.dataErr = dataErr;
    nmrFitData.TEs = TEs;
    assignin('base','nmrFitData',nmrFitData);

    assignin('base','Fscan',Fscan);
    assignin('base','INTscan',INTscan);
    PLOT = PLOT0; % Restore previous value
    assignin('base','PLOT',PLOT);
    set(handles.chkPLOT,'Value',PLOT)
    
    d = Fscan;
    for i=1:numel(indpar)
        d = [d INTscan(:,i) INTscanErr(:,i)];
    end
    num2clip(d);
    gdisp(handles,'Scanned spc copied to clipboard.  nmrFitData ready.');
    
%     for i=1:numel(indpar)
%         saveSUM = [toPPM(Fscan(:,1)), INTscan(:,i)];
%         %SPCsum(:,1) = (SPCsum(:,1)+1000*(FR(1)-REF))/REF*1000; % needs to be reimplemented for different frequences
%         nameOUT = [getFilePath(handles) 'LS-' LSname '-' num2str(i) '_' num2str(TEM) 'K.txt'];
%         save(nameOUT,'saveSUM','-ascii');
%     end


function [nameOUT] = getFilePath(handles)
mode = evalin('base','mode');
PATHname = evalin('base','PATHname');
        FileName = get(handles.edtFN,'String');
        [~, Name Ext] = fileparts(FileName);

        splitName = regexp(Name,mode,'split');
        nameOUT = [PATHname 'lineShift\' char(splitName(1))];
        

function [nameOUT] = getImagesPath(handles)
mode = evalin('base','mode');
PATHname = evalin('base','PATHname');
        FileName = get(handles.edtFN,'String');
        [~, Name Ext] = fileparts(FileName);

        splitName = regexp(Name,mode,'split');
        nameOUT = [PATHname 'images\' char(splitName(1))];
        
        
% --- Executes on button press in butOpen.
function butOpen_Callback(hObject, eventdata, handles)


    
    
% --- Executes on button press in butSave.
function butSave_Callback(hObject, eventdata, handles)


    


% --- Executes on button press in butAnaRun.
function butAnaRun_Callback(hObject, eventdata, handles)

TEs=fliplr(evalin('base','TEs'));
dialog = inputdlg({'Temperatures:','Update parameters'},'Analyse',1,{mat2str(TEs),'0'});
if isempty(dialog), return; end;

% saveFIT = 1;
TEs = str2num(dialog{1});
update = str2num(dialog{2});

% gdisp(handles,['T = ' num2str(TEs(1)) 'K']);
% set(handles.edtTE,'String',['[' num2str(TEs(1)) ']']);
% assignin('base','TEM',TEs(1));
% btnLoad_Callback(hObject, eventdata, handles)
% butAnalyse_Callback(hObject, eventdata, handles)

Trun=evalin('base','Trun');

% Initialize variables
SIMtype=evalin('base','SIMtype');
[funct, params] = fitNMRlib(SIMtype);
nparam = length(params);
allfitparM = zeros(length(TEs),nparam);
allfiterrM = zeros(length(TEs),nparam);
clipboardM = zeros(length(TEs),2*nparam+2);

% allfitparM(1,:)=evalin('base','allfitpar');
% allfiterrM(1,:)=evalin('base','allfiterr');
% sses = evalin('base','sse');
% clipboardM(1,:) = [TEs(1) sses reshape([allfitparM(1,:);allfiterrM(1,:)],1,[])];
% drawnow(); % why drawnow(handles.axeRes) not working!?

for j=1:length(TEs)
    if update && j > 1
        butCpyParams_Callback(hObject, eventdata, handles)
    end
    
    gdisp(handles,['T = ' num2str(TEs(j)) 'K']);
    set(handles.edtTE,'String',['[' num2str(TEs(j)) ']']);
    assignin('base','TEM',TEs(j));
    btnLoad_Callback(hObject, eventdata, handles)
    butAnalyse_Callback(hObject, eventdata, handles)

    allfitparM(j,1:nparam)=evalin('base','allfitpar');
    allfiterrM(j,1:nparam)=evalin('base','allfiterr');
    sses = evalin('base','sse');
    clipboardM(j,:) = [TEs(j) sses reshape([allfitparM(j,:);allfiterrM(j,:)],1,[])];
    drawnow();
end

assignin('base','allfitparM',allfitparM);
assignin('base','allfiterrM',allfiterrM);

num2clip(clipboardM);
gdisp(handles,'Data copied to clipboard.');

% narišem graf
if ~isempty(Trun)
	figure;
	COLS = 2;
	ROWS = ceil(length(Trun)/COLS);

    for i=1:length(Trun)
        TrunStr = Trun{i};
        
        if strfind(TrunStr, 'delta')
            h(i) = subplot(ROWS,COLS,i);
            plot(clipboardM(:,1),clipboardM(:,2),'o','Linewidth',2)
            title(TrunStr,'FontWeight','bold')
            continue
        end
        
        xlog = 0;
        ylog = 0;
        if strfind(TrunStr, 'xlog')
            xlog = 1;
            TrunStr = strrep(TrunStr, ' xlog', '');
        end
        if strfind(TrunStr, 'ylog')
            ylog = 1;
            TrunStr = strrep(TrunStr, ' ylog', '');    
        end

        % Replace parameter names with allfitparM(:,j)'
        for j=1:length(params)
            TrunStr = strrep(TrunStr, params{j}, ['(allfitparM(:,' num2str(j) ')'')']);
        end

        TrunStr = strrep(TrunStr, 'T', 'TEs'); % T1 is already replaced
        TrunStr = strrep(TrunStr, '*', '.*');
        TrunStr = strrep(TrunStr, '^', '.^');
        TrunStr = strrep(TrunStr, '/', './');
        TrunStr = strrep(TrunStr, '\', '.\');

        % PLOT
        h(i) = subplot(ROWS,COLS,i);
        plot(TEs,eval(TrunStr),'o','Linewidth',2)
        title(Trun{i},'FontWeight','bold')
        if xlog == 1, set(gca,'XScale','log'); end;
        if ylog == 1, set(gca,'YScale','log'); end;
    end
end

linkaxes(h, 'x');
splitName = regexp(getFilePath(handles),'\','split');
TITLE = splitName{end};
TITLE(TITLE=='_')=' ';
mtit(TITLE,'fontsize',14,'color','b','xoff',-.1,'yoff',.015);

num2clip(clipboardM);
gdisp(handles,'Data copied to clipboard.');
assignin('base','clipboardM',clipboardM);




% --- Executes on button press in butSpecter.
function butSpecter_Callback(hObject, eventdata, handles)

spcSumInd = evalin('base','spcSumInd');
SPCall = evalin('base','SPCall');
Fall = evalin('base','Fall');
range = evalin('base','range');
ppmRange = evalin('base','ppmRange');
autoPH = evalin('base','autoPH');
PLOT = evalin('base','PLOT');
% SPCsum = evalin('base','SPCsum');
% FR = evalin('base','FR');
% REF = evalin('base','REF');
TEM = evalin('base','TEM');
IND = evalin('base','IND');
NORM = evalin('base','NORM');
nParams = evalin('base','nParams');
mode = evalin('base', 'mode');
isPPM = evalin('base', 'isPPM');

if strcmp(mode,'SW')
   msgbox('Please press Sweep button.','Error');
   return
end

if ~isempty(spcSumInd)
    SPCsum = zeros(numel(real(SPCall{1})),2);
    if numel(IND) < numel(spcSumInd), spcSumInd = IND; end

    SPCsum(:,1) = Fall{1};
    for i=1:length(spcSumInd)
        SPCsum(:,2) = SPCsum(:,2) + real(SPCall{spcSumInd(i)});
    end

    % Normalize
    SPCsum = dataNorm(SPCsum,NORM,nParams);

    if isPPM 
        SPCsum(:,1) = toPPM(SPCsum(:,1));
    end

    if PLOT
        a4 = handles.axeRes;
        plot(a4,SPCsum(:,1),SPCsum(:,2))
        if isPPM
            range = ppmRange;
            xlabel(a4,'\nu - \nu_r_e_f [ppm]')
        else
            xlabel(a4,'\nu [kHz]') 
        end
        if max(SPCsum(:,1)) < range(1) || min(SPCsum(:,1)) > range(2)
           msgbox('Data outside specified range!','Error!')
           xlim(a4,[min(SPCsum(:,1)) max(SPCsum(:,1))])
        else
           xlim(a4,range)
        end
        grid(a4,'on')
        drawnow;
    end

    if autoPH == 0
        gdisp(handles,'Warning: autoPH = 0, not good for summing spectra!')
    end
else
    msgbox('scpSumInd is empty!!!');
    return
end

num2clip(SPCsum);
assignin('base','SPCsum',SPCsum);

nmrFitData.data{1} = SPCsum;
nmrFitData.TEs{1} = TEM;
assignin('base','nmrFitData',nmrFitData);
gdisp(handles,'Current specter copied to clipboard. Sum of all spectra is in SPCsum. nmrFitData ready.');


% --- Executes on button press in butCpyParams.
function butCpyParams_Callback(hObject, eventdata, handles)

    param = get(handles.edtAnalyse,'String');
    allfitpar=evalin('base','allfitpar');
%     SIMtype=evalin('base','SIMtype');
% 
%     if evalin('base','exist(''SIMtype_1'',''var'')')
%         SIMtype_1=evalin('base','SIMtype_1');
%     else
%         SIMtype_1 = '';
%     end
%     if evalin('base','exist(''fraction'',''var'')')
%         fraction=evalin('base','fraction');
%     else
%         fraction = 0.5; 
%     end
    
    %gdisp(handles,allfitpar);
    for i=1:numel(param)
        if strfind(param{i},'allinitpar = ')
%             if strcmp('T1-2Exp_fraction',SIMtype) %problem èe je fraction 1!! hmmm
%                 if allfitpar(2)>0.999
%                     allfitpar(2)=0.9999;
%                 end
%             end
%             %
%             if strcmp('T1-3Exp_fraction',SIMtype_1)
%                 if fraction == 0
%                     allfitpar(2)=0;
%                 elseif fraction >0.99
%                     allfitpar(5)=0;
%                 end
%             end

            str = []; 
            for j=1:numel(allfitpar)
                str = [str num2str(allfitpar(j)) ','];
            end
            str = ['[' str(1:end-1) ']'];
            param{i} = ['allinitpar = ' str ';'];
        end
    end
    set(handles.edtAnalyse,'String',param);

    
    
% --- Executes on button press in butCpyParams.
function updatePH_v(handles,PH_v)
    param = get(handles.edtFFTparam,'String');
    for i=1:numel(param)
        if strfind(param{i},'PH_v = ')
            param{i} = ['PH_v = ' num2str(PH_v) ';'];
        end
    end
    set(handles.edtFFTparam,'String',param);
    
    
    
% function onlyFitParams(hObject, eventdata, handles)
% % Copy only those reset parameters into initparam that are not fitted
%     butParams_Callback(hObject, eventdata, handles);
%     reloadSettings(handles);
%     param = get(handles.edtAnalyse,'String');
%     indFitpar=evalin('base','indFitpar');
%     
%     allinitpar = evalin('base','allinitpar');
%     allresetpar = evalin('base','allresetpar');
% 
%         for i=1:numel(allresetpar) %changed on 14.8.2011 from allinitpar
%            if ~ismember(i,indFitpar)
%               allinitpar(i) = allresetpar(i); 
%            end
%         end
%         
%         if numel(allinitpar)~=numel(allresetpar) %added on 14.8.2011
%            allinitpar_temp = zeros(numel(allresetpar),1);
%            for i=1:numel(allresetpar) 
%                   allinitpar_temp(i) = allinitpar(i); 
%            end
%            allinitpar = allinitpar_temp; 
%         end
%         
%         for i=1:numel(param)
%             if strfind(param{i},'allinitpar = ')
%                 str = mat2str(allinitpar);
%                 str(str==' ')=',';
%                 param{i} = ['allinitpar = ' str ';'];
%             end
%         end
%         assignin('base','allinitpar',allinitpar);
%         set(handles.edtAnalyse,'String',param);

        
        
% --- Executes on button press in butSpecterRun.
function butSpecterRun_Callback(hObject, eventdata, handles)

    TEs = fliplr(evalin('base','TEs')); % Must be the same as in Analyse Run

    dialog = inputdlg({'Temperature shift [K^-1]:','Specter scale','Temperatures:'},'Specter plot',1,{'1','1',mat2str(TEs)});
    if isempty(dialog), return; end;
    
    Tshift = str2double(dialog{1});
    A = str2double(dialog{2});
    TEs = eval(dialog{3});
    %pause(0.1);

    isPPM = evalin('base','isPPM');
    if isPPM 
        range = evalin('base','ppmRange');
    else
        range = evalin('base','range');
    end
    Fmin = min(range); Fmax = max(range);

    if evalin('base','exist(''clipboardM'',''var'')')
        clipboardM = evalin('base','clipboardM');
        if size(clipboardM,2)==numel(TEs)*2
            if strcmp(questdlg('Do you want to reload spectra?','Run Spectra','Yes','No','Yes'),'Yes')
                clipboardM = [];
            end
        else
            clipboardM = [];
        end
    else
        clipboardM = [];
    end


    if isempty(clipboardM)
        for j=1:length(TEs)
            gdisp(handles,['T = ' num2str(TEs(j)) 'K']);
            set(handles.edtTE,'String',['[' num2str(TEs(j)) ']']);
            btnLoad_Callback(hObject, eventdata, handles)
            butSpecter_Callback(hObject, eventdata, handles)
            SPCsum = evalin('base','SPCsum');

            %SPCsum(:,2) = SPCsum(:,2)+TEs(j)*Tshift;
            if isempty(clipboardM)
                clipboardM(:,2*j-1) = SPCsum(:,1);
                clipboardM(:,2*j) = SPCsum(:,2);
            else
                nrows = numel(SPCsum(:,1));
                cpbrows = numel(clipboardM(:,1));

                if nrows>cpbrows
                    clipboardM = [clipboardM; NaN(nrows-cpbrows,numel(clipboardM(1,:)))];
                end

                clipboardM(1:nrows,2*j-1) = SPCsum(:,1);
                clipboardM(1:nrows,2*j) = SPCsum(:,2);
                if nrows<cpbrows
                    clipboardM(nrows+1:nrows,2*j-1) = NaN;
                    clipboardM(1:nrows,2*j) = NaN;
                end
            end
            drawnow();
        end
    end

    clipboardMunscaled = clipboardM;
    for j=1:length(TEs)
        clipboardM(:,2*j) = A*clipboardM(:,2*j)+TEs(j)*Tshift;
    end

    % narišem graf
            a4 = handles.axeRes;
            plot(a4,clipboardM(:,1),clipboardM(:,2))
            [max_val,location] = max(max(clipboardM(:,1),[],2));
            hold(a4,'on');
            max_val = min([Fmax max_val]);
            text(max_val,clipboardM(location,2),[num2str(TEs(1)),' K'],'FontSize',10,'FontWeight','bold','Color','b','Parent',a4)
            for j=2:length(TEs)
                plot(a4,clipboardM(:,2*j-1),clipboardM(:,2*j))
                if mod(j,4)==0
                    [max_val,location] = max(max(clipboardM(:,2*j-1),[],2));
                    max_val = min([Fmax max_val]);
                    text(max_val,clipboardM(location,2*j),[num2str(TEs(j)),' K'],'FontSize',10,'FontWeight','bold','Color','b','Parent',a4)
                end
            end
            hold(a4,'off');
            xlim(a4,[Fmin Fmax])
            if isPPM
                xlabel(a4,'\nu - \nu_r_e_f [ppm]')
            else
                xlabel(a4,'\nu - \nu_r_e_f [kHz]')
            end

    num2clip(clipboardM);
    assignin('base','clipboardM',clipboardMunscaled);
    nmrFitData.data{1} = clipboardMunscaled;
    nmrFitData.TEs{1} = TEs;
    assignin('base','nmrFitData',nmrFitData);
    gdisp(handles,'Data copied to clipboard. nmrFitData ready.');


% --- Executes on button press in butLSRun.
function butLSRun_Callback(hObject, eventdata, handles)
    if ~evalin('base','exist(''allfitparM'',''var'')');
        msgbox('First make a Analyse Run!')
        return
    end
    allfitparM = evalin('base','allfitparM');
    TEs0 = fliplr(evalin('base','TEs'));
    
    dialog = inputdlg({'Plot parameters [ind]:','Frequency interval [kHz]:','Frequency step [kHz]:','Add plot [sorted ind]:','Spc normalize [1/0]:','Inversion Recovery [1/0]:','Temperatures:'},'Line Scan',1,{'[2 5]','3','1','1','1','1',mat2str(TEs0)});
    if isempty(dialog), return; end;
    indpar = str2num(dialog{1});
    dF = str2double(dialog{2});
    step =  str2double(dialog{3});
    plotind = str2double(dialog{4});
    isNorm = str2double(dialog{5});
    isIR = str2double(dialog{6});
    TEs = str2num(dialog{7});
    
    if numel(TEs) > size(allfitparM,1)
    	msgbox('First make a Analyse Run with the same TEs!')
        return
    end 
    data=cell(2,1);
    dataErr=cell(2,1);
    TEss=cell(2,1);
    maximums = cell(2,1);
    
    for i = 1:numel(TEs)
        assignin('base','PLOT',0)
        LSpar{1} = 1;
        LSpar{2} = indpar;
        LSpar{3} = dF;
        LSpar{4} = step;
        LSpar{5} = plotind;
        LSpar{6} = isNorm;
        LSpar{7} = isIR;
        assignin('base','LSpar',LSpar)

        ind = find(TEs0==TEs(i));
        if isempty(ind)
           msgbox('TEM not found!!!')
           return;
        end
        
        allfitpar = allfitparM(ind,:);
        
        % Change Afact to B2=B1*Afact
        allfitpar(5) = allfitpar(2)*allfitpar(5);
        
        assignin('base','allfitpar',allfitpar);
        assignin('base','TEM',TEs(i));
        set(handles.edtTE,'String',['[' num2str(TEs(i)) ']']);
        
        butCpyParams_Callback(hObject, eventdata, handles)
        btnLoad_Callback(hObject, eventdata, handles)
        butAnalyse_Callback(hObject, eventdata, handles)
        butScan_Callback(hObject, eventdata, handles)
        
        a4 = handles.axeRes;
        hold(a4,'on')
        color={'b','r'};
        nmrFiDa = evalin('base','nmrFitData');
        for j=1:numel(indpar)
           data{j} = [data{j} nmrFiDa.data{j}];
           dataErr{j} = [dataErr{j} nmrFiDa.dataErr{j}];
           TEss{j} = [TEss{j} nmrFiDa.TEs{j}];
           
           [~, imax] = max(nmrFiDa.data{j}(:,2));
           fmax = toPPM(nmrFiDa.data{j}(imax,1),1);
           maximums{j} = [maximums{j} fmax];
           plot(a4,[fmax fmax],[0 max(data{j}(:,end))],color{j});
        end
        hold(a4,'off')
    end
    
    nmrFitData.data = data;
    nmrFitData.dataErr = dataErr;
    nmrFitData.TEs = TEss;
    assignin('base','nmrFitData',nmrFitData);
    assignin('base','maximums',maximums);
    evalin('base','clear LSpar');
    
    gdisp(handles,'Scan Line Run finnished. nmrFitData ready.');
    

% --- Executes on button press in butYscale.
function butYscale_Callback(hObject, eventdata, handles)

a4 = handles.axeRes;
if strcmp(get(a4,'YScale'),'linear')
    set(a4,'YScale','log');
else
    set(a4,'YScale','lin');
end
drawnow();


% --- Executes on button press in butXscale.
function butXscale_Callback(hObject, eventdata, handles)

a4 = handles.axeRes;
if strcmp(get(a4,'XScale'),'linear')
    set(a4,'XScale','log');
else
    set(a4,'XScale','lin');
end
drawnow();


% --- Executes on button press in butResetParams.
function butResetParams_Callback(hObject, eventdata, handles)
    reloadSettings(handles);

    param = get(handles.edtAnalyse,'String');
    allresetpar = evalin('base','allresetpar');
    
    %gdisp(handles,allresetpar);
    for i=1:numel(param)
        if strfind(param{i},'allinitpar = ')
            str = [];
            for j=1:numel(allresetpar)
                str = [str num2str(allresetpar(j)) ','];
            end
            str = ['[' str(1:end-1) ']'];
            param{i} = ['allinitpar = ' str ';'];
        end
    end
    
    set(handles.edtAnalyse,'String',param);
    reloadSettings(handles);




function plotSpectra(handles, fname)
plotind = 1; % plot saved specter on top of others if exsists
TEs=fliplr(evalin('base','TEs'));
mode=evalin('base','mode');
range = evalin('base','range');
Dy = evalin('base','Dy');
Fmin = min(range); Fmax = max(range);
indpar = [1]; nameIN = {}; TITLE = '';
    if evalin('base','exist(''LSpar'',''var'')')
        if evalin('base','exist(''LSname'',''var'')')
            LSname = evalin('base','LSname');
        else
            LSname = 'noname';
        end
        LSpar = evalin('base', 'LSpar');
        if LSpar{1} == 1
            indpar = LSpar{4};
            for i=1:numel(indpar)
                nameIN{i} = [fname 'LS-' LSname '-' num2str(indpar(i)) '_'];
            end
            TITLE = [fname 'LS-' LSname];
        end
    elseif strcmp(mode,'SW') 
        if evalin('base','exist(''SWname'',''var'')')
            SWname = evalin('base','SWname');
        else
            SWname = 'noname';
        end
        nameIN{1} = [fname 'SW-' SWname '_'];
    else
        nameIN{1} = [fname mode '_'];
    end
    
    
    PATHname = getFilePath(handles);
    ZOOM = 30;
    GRAPHS = ceil(sqrt(numel(TEs)));
    COLS = ceil(sqrt(numel(TEs)));
    ROWS = 1;
    
    normal=zeros(numel(TEs),numel(nameIN)+1);
    if plotind > 0
        for i=1:length(TEs)
            SPC = load([PATHname mode '_' num2str(TEs(i)) 'K.txt'],'-ascii');
            area(i,numel(nameIN)+2) = trapz(SPC(:,1),SPC(:,2));
            normal(i,numel(nameIN)+1) = ZOOM/area(i,numel(nameIN)+2);%area(i,numel(ind)+2)/sum(area(1,2:(end-1)),2);
        end
    end
    for i=1:length(TEs)
        for j=1:numel(nameIN)
            SPC = load([PATHname 'LS-' LSname '-' num2str(j) '_' num2str(TEs(i)) 'K.txt'],'-ascii');
            area(i,j+1) = trapz(SPC(:,1),SPC(:,2));
        end
        %teflonski shit fixing: v bistvu ga relativno na druge
        %renormaliziramo, saj poznamo obmoèje kjer je 0!!
        if numel(nameIN)==3 && plotind==1
            j=1; %prvi
            SPC = load([PATHname 'LS-' LSname '-' num2str(j) '_' num2str(TEs(i)) 'K.txt'],'-ascii');
            inds = find(SPC(:,1)>=7&SPC(:,1)<=30);
            area1 = trapz(SPC(inds,1),SPC(inds,2));
            j=2; %drugi
            SPC = load([PATHname 'LS-' LSname '-' num2str(j) '_' num2str(TEs(i)) 'K.txt'],'-ascii');
            inds = find(SPC(:,1)>=7&SPC(:,1)<=30);
            area2 = trapz(SPC(inds,1),SPC(inds,2));
            
            %skupaj
            SPC = load([PATHname mode '_' num2str(TEs(i)) 'K.txt'],'-ascii');
            inds = find(SPC(:,1)>=7&SPC(:,1)<=30);
            areaTot = trapz(SPC(inds,1),SPC(inds,2));
            
            deltaS = areaTot/(area1+area2);
            
            norm_teflon = (area(i,5)-(area(i,2)+area(i,3))*deltaS)/(deltaS*area(i,4));
            if (area(i,5)-(area(i,2)+area(i,3))*deltaS)<0 norm_teflon = 0; end;
            if TEs(i)<100 norm_teflon=0; end;
            if TEs(i)>=175 norm_teflon=norm_teflon*1.5; end;
            area(i,4) = area(i,4)*norm_teflon;
        end
            
        for j=1:numel(nameIN)
            normal(i,j) = ZOOM/sum(area(i,1:(end-1)),2);
        end
        if numel(nameIN)==3 && plotind==1
            normal(i,3) = normal(i,3)*norm_teflon;
        end
    end
    
    figure;
    COLOR = {'b','r','g','m','k','b','r','g','m'};
    for i=1:length(TEs)
        if mod(i-1,GRAPHS)==0
            %hold off;
            h(ceil(i/GRAPHS))=subplot(ROWS,COLS,ceil(i/GRAPHS));
            hold on;
            grid on;
            xlabel('\nu-\nu_{ref} [ppm]'); ylabel('Intensity');
            axis([toPPM(Fmin),toPPM(Fmax),-Dy*GRAPHS*0.05,Dy*GRAPHS*1.05]);
            %set(gca,'FontSize',8);
            %set(gca,'YTick',[]);
        end
        for ind=1:numel(nameIN)
            SPC = load([nameIN{ind},num2str(TEs(i)),'K.txt'],'-ascii');
            if ~isempty(normal)
                plot(SPC(:,1),SPC(:,2)*normal(i,ind)+Dy*(GRAPHS-1)-Dy*(mod(i-1,GRAPHS)),'-k','Linewidth',1.5,'Color',COLOR{ind})
            else
                plot(SPC(:,1),SPC(:,2)/max(SPC(:,2))+Dy*(GRAPHS-1)-Dy*(mod(i-1,GRAPHS)),'-k','Linewidth',1.5,'Color',COLOR{ind})
            end
        end
        if plotind > 0
                SPC = load([getFilePath(handles) mode '_' num2str(TEs(i)) 'K.txt'],'-ascii');
                plot(SPC(:,1),SPC(:,2)*normal(i,numel(nameIN)+1)+Dy*(GRAPHS-1)-Dy*(mod(i-1,GRAPHS)),'-k','Linewidth',1.5,'Color','black')
        end
        text(toPPM(Fmin+(Fmax-Fmin))*0.80,0.4*Dy+Dy*(GRAPHS-1)-Dy*(mod(i-1,GRAPHS)),[num2str(TEs(i)),' K'],'FontSize',12,'FontWeight','bold','Parent',gca)
    end

for j=1:numel(indpar)
    LEG{j} = ['indpar = ' num2str(indpar(j))];
end
if plotind > 0
    LEG{numel(indpar)+1} = 'T1 specter';
end

splitName = regexp(TITLE,'\','split');
TITLE = splitName{end};
TITLE(TITLE=='_')=' ';
mtit(TITLE,'fontsize',14,'color','b','xoff',-.1,'yoff',.015);
legend(LEG,2,'Location','NorthEast');
linkaxes(h, 'x');
linkaxes(h, 'y');

hold off; %grid on;
saveImage(handles,gca,[LSname '_spc']);


function saveImage(handles,fig,name)
if evalin('base','exist(''saveIMG'',''var'')')
    saveIMG = evalin('base','saveIMG');
else
    return;
end
if evalin('base','exist(''extension'',''var'')')
    ext = evalin('base','extension');
else
    ext='pdf';
end
if saveIMG
    set(gcf, 'Position', [0 0 1200 800]);
    oldSettings = fillPage(gcf, 'papersize', 'A4');
    saveas(fig,[getImagesPath(handles) name '.' ext]); 
    close(gcf);
end




function setName(handles,name)
    param = get(handles.edtAnalyse,'String');
    for i=1:numel(param)
        if strfind(param{i},'LSname = ')
            param{i} = ['LSname = ''' name ''';'];
        end
    end
    set(handles.edtAnalyse,'String',param);
    reloadSettings(handles);


% --- Executes on button press in butCustom.
function butCustom_Callback(hObject, eventdata, handles)
reloadSettings(handles);reloadFFT(handles);
list = nmrCustom(get(handles.popCustom,'Value'),hObject, eventdata,handles);
set(handles.popCustom,'String',list);


function saveImageWithinGUI(handles,axesObject,name,legendObject)
%this function takes in two arguments
%axesObject is the axes object that will be saved (required input)
%legendObject is the legend object that will be saved (optional input)
%legendObject =  findobj(axesObject,'Type','axes','Tag','legend');
set( 0, 'DefaultFigurePaperOrientation', 'landscape');
if evalin('base','exist(''saveIMG'',''var'')')
    saveIMG = evalin('base','saveIMG');
else
    saveIMG=0;
end
if evalin('base','exist(''extension'',''var'')')
    ext = evalin('base','extension');
else
    ext='pdf';
end
if saveIMG
    filepath = [getImagesPath(handles) name '.' ext];
else
    return;
end

%create a new figure
newFig = figure;

splitName = regexp(filepath,'\','split');
splitName = regexp(splitName{end},'\.','split');
TITLE = splitName{1};
TITLE(TITLE=='_')=' ';
%text(0.5,0.95,TITLE,'FontSize',12,'FontWeight','bold','HorizontalAl','Right');
%get the units and position of the axes object
axes_units = get(axesObject,'Units');
axes_pos = get(axesObject,'Position');
%mtit(newFig,TITLE,'fontsize',14,'color','b');
%copies axesObject onto new figure
axesObject2 = copyobj(axesObject,newFig);

%realign the axes object on the new figure
set(axesObject2,'Units',axes_units);
set(axesObject2,'Position',[15 5 1.8*axes_pos(3) 1.8*axes_pos(4)]);

%if a legendObject was passed to this function . . .
if (exist('legendObject'))
 
    %get the units and position of the legend object
    legend_units = get(legendObject,'Units');
    legend_pos = get(legendObject,'Position');
 
    %copies the legend onto the the new figure
    legendObject2 = copyobj(legendObject,newFig);
 
    %realign the legend object on the new figure
    set(legendObject2,'Units',legend_units);
    set(legendObject2,'Position',[15-axes_pos(1)+legend_pos(1) 5-axes_pos(2)+legend_pos(2) legend_pos(3) legend_pos(4)] );
 
end
 
%adjusts the new figure accordingly
set(newFig,'Units',axes_units);
set(newFig,'Position',[15 10 1.8*(axes_pos(3))+30 1.8*(axes_pos(4))+10]);

title(TITLE);
%set(newFig, 'Position', [0 0 1200 800]);
%saves the plot
oldSettings = fillPage(gcf, 'papersize', 'A4');
saveas(gca,filepath) 
 
%closes the figure
close(newFig)


% --- Executes on button press in butReload.
function butReload_Callback(hObject, eventdata, handles)



function dir5CB(handles)
fp = get(handles.edtFN,'String');
fp1 = regexp(fp,'\','split');

str = {};
str = [str,get(handles.dir5,'String')];
old = str{get(handles.dir5,'Value')};

fp1{end-5} = old;
fp = [fp1{1}];
for i=2:numel(fp1)
   fp = [fp '\' fp1{i}];
end
set(handles.edtFN,'String',fp);


function dir4CB(handles)
fp = get(handles.edtFN,'String');
fp1 = regexp(fp,'\','split');

str = {};
str = [str,get(handles.dir4,'String')];
old = str{get(handles.dir4,'Value')};

fp1{end-4} = old;
fp = [fp1{1}];
for i=2:numel(fp1)
   fp = [fp '\' fp1{i}];
end
set(handles.edtFN,'String',fp);


function dir3CB(handles)
fp = get(handles.edtFN,'String');
fp1 = regexp(fp,'\','split');

str = {};
str = [str,get(handles.dir3,'String')];
old = str{get(handles.dir3,'Value')};

fp1{end-3} = old;
fp = [fp1{1}];
for i=2:numel(fp1)
   fp = [fp '\' fp1{i}];
end
set(handles.edtFN,'String',fp);


function dir2CB(handles)
fp = get(handles.edtFN,'String');
fp1 = regexp(fp,'\','split');

str = {};
str = [str,get(handles.dir2,'String')];
old = str{get(handles.dir2,'Value')};

fp1{end-2} = old;
fp = [fp1{1}];
for i=2:numel(fp1)
   fp = [fp '\' fp1{i}];
end
set(handles.edtFN,'String',fp);


function dir1CB(handles)
fp = get(handles.edtFN,'String');
fp1 = regexp(fp,'\','split');

str = {};
str = [str,get(handles.dir1,'String')];
old = str{get(handles.dir1,'Value')};

fp1{end-1} = old;
fp = [fp1{1}];
for i=2:numel(fp1)
   fp = [fp '\' fp1{i}];
end
set(handles.edtFN,'String',fp);

% --- Executes on selection change in dir3.
function dir3_Callback(hObject, eventdata, handles)
dir3CB(handles);

% --- Executes during object creation, after setting all properties.
function dir3_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in dir2.
function dir2_Callback(hObject, eventdata, handles)
dir2CB(handles);

% --- Executes during object creation, after setting all properties.
function dir2_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in dir1.
function dir1_Callback(hObject, eventdata, handles)
dir1CB(handles);

% --- Executes during object creation, after setting all properties.
function dir1_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in chkDE.
function chkDE_Callback(hObject, eventdata, handles)


% --- Executes on button press in chkAutoSHL.
function chkAutoSHL_Callback(hObject, eventdata, handles)


% --- Executes on button press in chkAutoPH.
function chkAutoPH_Callback(hObject, eventdata, handles)


% --- Executes on button press in chkPLOT.
function chkPLOT_Callback(hObject, eventdata, handles)
PLOT = get(hObject,'Value');
assignin('base','PLOT',PLOT);



function edtLog_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function edtLog_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edtFN_Callback(hObject, eventdata, handles)


% --- Executes on button press in butGo.
function butGo_Callback(hObject, eventdata, handles)


% --- Executes on selection change in dir5.
function dir5_Callback(hObject, eventdata, handles)
dir5CB(handles);

% --- Executes during object creation, after setting all properties.
function dir5_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in dir4.
function dir4_Callback(hObject, eventdata, handles)
dir4CB(handles);

% --- Executes during object creation, after setting all properties.
function dir4_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function reloadSettings(handles)
param = get(handles.edtAnalyse,'String');
evalin('base',sprintf('%s\n',param{:}));


function reloadFFT(handles)
param = get(handles.edtFFTparam,'String');
evalin('base',sprintf('%s\n',param{:}));





% --- Executes on selection change in popCustom.
function popCustom_Callback(hObject, eventdata, handles)

list = nmrCustom(0,hObject, eventdata,handles);
set(handles.popCustom,'String',list);

% --- Executes during object creation, after setting all properties.
function popCustom_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in butSpecterSim.
function butSpecterSim_Callback(hObject, eventdata, handles)


% --- Executes on selection change in dirFiles.
function dirFiles_Callback(hObject, eventdata, handles)

fp2 = get(handles.edtFN,'String');
fpindFN = regexp(fp2,'\','start');
fp = fp2(1:fpindFN(end));
fileslist = get(handles.dirFiles,'String');
filesind = get(handles.dirFiles,'Value');
set(handles.edtFN,'String',[fp fileslist{filesind}]);

% --- Executes during object creation, after setting all properties.
function dirFiles_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when figure1 is resized.
function figure1_ResizeFcn(hObject, eventdata, handles)


% --- Executes on button press in butParamLoad.
function butParamLoad_Callback(hObject, eventdata, handles)


% --- Executes on button press in butMoments.
function butMoments_Callback(hObject, eventdata, handles)
CUT2 = evalin('base','CUT2');
isPPM = evalin('base','isPPM');
ppmRange = evalin('base','ppmRange');
PLOT = evalin('base','PLOT');
nmrFitData = evalin('base','nmrFitData');
range = evalin('base','range');

ns = numel(nmrFitData.TEs);
ni = numel(nmrFitData.TEs{1});
if ns > 1 || ni > 1
    dialog = inputdlg({'Select one field to analyse:','Select one temperature to analyse:'},'Multiple spc found',1,{mat2str(1:ns),mat2str(1:ni)});
    if isempty(dialog), return; end;
    ns = str2num(dialog{1});
    ni = str2num(dialog{2});
    ns = ns(1); % if more are left selected
    ni = ni(1); % if more are left selected
end % otherwise its 1

% dialog = inputdlg({'Specify offset:'},'Moments',1,{'0'});
% if isempty(dialog), return; end;
% offset = str2num(dialog{1});

X = nmrFitData.data{ns}(:,2*ni-1);
Y = nmrFitData.data{ns}(:,2*ni);
TEM = nmrFitData.TEs{ns}(ni);

if isPPM
    range = ppmRange;
    ind = X<range(1) | X>range(2); 
    X(ind)=[];
    Y(ind)=[];
end

[M0 M1 M2] = nrNMRanalysis(X,Y,0,CUT2)

if PLOT
    a4 = handles.axeRes;
    plot(a4,X,Y,'b-')
    hold(a4,'on');
    if CUT2 > -inf
        plot(a4,[min(X) max(X)],[CUT2 CUT2],'k-','LineWidth',2)
    end
    plot(a4,[M1 M1],[0 1.05*max(Y)],'r-.','LineWidth',2)
    plot(a4,[M1-M2 M1-M2],[0 1.05*max(Y)],'g--','LineWidth',2)
    plot(a4,[M1+M2 M1+M2],[0 1.05*max(Y)],'g--','LineWidth',2)
    hold(a4,'off');
    grid(a4,'on');
    
    if isPPM
        range = ppmRange;
        xlabel(a4,'\nu - \nu_r_e_f [ppm]')
    else
        xlabel(a4,'\nu [kHz]') 
    end
    if max(X) < range(1) || min(X) > range(2)
       msgbox('Data outside specified range!','Error!')
       xlim(a4,[min(X) max(X)])
    else
       xlim(a4,range)
    end

    Tx = 0.02; Ty = 0.5; dTy = 0.07;
    text(Tx,Ty+2*dTy,[ 'M0 = ', num2str(M0)],'Units','Norm','Fontweight','bold','FontSize',12,'Color','red','Parent',a4,'BackgroundColor',[1 1 1]);            
    text(Tx,Ty+1*dTy,[ 'M1 = ', num2str(M1)],'Units','Norm','Fontweight','bold','FontSize',12,'Color','red','Parent',a4,'BackgroundColor',[1 1 1]);            
    text(Tx,Ty+0*dTy,[ 'sM2 = ', num2str(2*M2)],'Units','Norm','Fontweight','bold','FontSize',12,'Color','red','Parent',a4,'BackgroundColor',[1 1 1]);            
    text(Tx,Ty-1*dTy,[ 'CUT2 = ', num2str(CUT2)],'Units','Norm','Fontweight','bold','FontSize',12,'Color','red','Parent',a4,'BackgroundColor',[1 1 1]);            

end
    % peak to peak for gaussian is 2*sM2!!!
    % 2*M2 contains 68% of gaussian signal (like 2/3)
               
num2clip([TEM M0 M1 2*M2]);
gdisp(handles,['M0 = ' num2str(M0) ' M1 = ' num2str(M1) ' M2 = ' num2str(2*M2) '.  copied!']);



% --- Executes on key press with focus on edtFFTparam and none of its controls.
function edtFFTparam_KeyPressFcn(hObject, eventdata, handles)


% --- Executes on button press in butCpy2nmrFit.
function butCpy2nmrFit_Callback(hObject, eventdata, handles)



function name = extractName(fileName)
    if strfind(fileName,'/')
        [~, fileName, ~] = fileparts(fileName);
    end
    if strfind(fileName,'\')
        [~, fileName, ~] = fileparts(fileName);
    end
    
    ind = strfind(fileName,'_');
    name = fileName(1:(ind(end)-1));


% --- Executes during object creation, after setting all properties.
function butScan_CreateFcn(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function butMoments_CreateFcn(hObject, eventdata, handles)


% --- Executes on button press in butScanIND.
function butScanIND_Callback(hObject, eventdata, handles)




% --- Executes on button press in nutLeft.
function nutLeft_Callback(hObject, eventdata, handles)
TEs = evalin('base','TEs');
TE = str2num(get(handles.edtTE,'String'));

pos = find(TEs==TE);

if pos > 1
    TE = TEs(pos-1);
    set(handles.edtTE,'String',num2str(TE));
    btnLoad_Callback(hObject, eventdata, handles);
%     butAnalyse_Callback(hObject, eventdata, handles);
end


% --- Executes on button press in cpyAnaResults.
function cpyAnaResults_Callback(hObject, eventdata, handles)




% --- Executes during object creation, after setting all properties.
function butAnalyse_CreateFcn(hObject, eventdata, handles)


% --- Executes on button press in ButAnaRight.
function ButAnaRight_Callback(hObject, eventdata, handles)
TEs = evalin('base','TEs');
TE = str2num(get(handles.edtTE,'String'));

pos = find(TEs==TE);

if pos < numel(TEs)
    TE = TEs(pos+1);
    set(handles.edtTE,'String',num2str(TE));
    btnLoad_Callback(hObject, eventdata, handles);
%     butAnalyse_Callback(hObject, eventdata, handles);
end


% --- Executes on button press in butHyperFit.
function butHyperFit_Callback(hObject, eventdata, handles)
    reloadSettings(handles);
    TEM=evalin('base','TEM');
    %integNMR(handles);
    hyperfitNMR(handles);
    


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function mnuExportSPC_Callback(hObject, eventdata, handles)
SPCall = evalin('base','SPCall');
Fall = evalin('base','Fall');
Fall = evalin('base','Fall');
TAU = evalin('base','TAU');

num2clip([0 TAU; Fall{1} real(cell2mat(SPCall))])


% --------------------------------------------------------------------
function mnuExportMag_Callback(hObject, eventdata, handles)
INTall = evalin('base','INTall');
mode = evalin('base','mode');

switch mode
    case {'T1','T2'}
        xPar = evalin('base','TAU');
    case 'D1'
        xPar = evalin('base','D1');
    case 'D2'
        xPar = evalin('base', 'D2');
    case 'SW'
        xPar = evalin('base', 'FR');
    case 'SPC'
        xPar = evalin('base','TAU');
    otherwise
        error(['Unknown mode ' mode ' !!!']);
end

num2clip([xPar' INTall])


% --------------------------------------------------------------------
function mnuExportFit_Callback(hObject, eventdata, handles)
fitT1 = evalin('base','fitT1');
fitINTall = evalin('base','fitINTall');
INTall = evalin('base','INTall');
mode = evalin('base','mode');

switch mode
    case {'T1','T2'}
        xPar = evalin('base','TAU');
    case 'D1'
        xPar = evalin('base','D1');
    case 'D2'
        xPar = evalin('base', 'D2');
    case 'SW'
        xPar = evalin('base', 'FR');
    case 'SPC'
        xPar = evalin('base','TAU');
    otherwise
        error(['Unknown mode ' mode ' !!!']);
end

% num2clip([xPar' INTall fitINTall cell2mat(fitT1)])
num2clip([xPar' INTall fitINTall])
% TAU = evalin('base','TAU');
% num2clip([TAU' INTall])


% --- Executes on button press in butCpyMagnetization.
function butCpyMagnetization_Callback(hObject, eventdata, handles)


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over butShowIND.
function butShowIND_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to butShowIND (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function mnuCpyFitRes_Callback(hObject, eventdata, handles)
allfitpar = evalin('base','allfitpar');
allfiterr = evalin('base','allfiterr');

res = reshape([allfitpar; allfiterr],1,[]);
    
num2clip(res)
gdisp(handles,'Data copied to clipboard.');


% --------------------------------------------------------------------
function mnuCpySpc2nmrFit_Callback(hObject, eventdata, handles)
    isPPM = evalin('base','isPPM');
    if ~evalin('base','exist(''nmrFitData'',''var'')')
        msgbox('No spectra available!','Error');
        return
    end
    nmrFitData = evalin('base','nmrFitData');
    ns = numel(nmrFitData.TEs);
    if ns > 1
        answer = inputdlg({'Choose one field:'},'Cpy to nmrFit',1,{mat2str(1:ns)});
        ns = str2num(answer{1});
        ns = ns(1);
    end
    ind = 1;
    if evalin('base','exist(''nmr'',''var'')')
        answer = questdlg('Append specter to existing nmrFit?', ...
                          'Cpy to nmrFit','Yes','No','Cancel','Yes');
        nmr = evalin('base','nmr');
        switch answer
            case 'Cancel'
                return
            case 'Yes'
                if isfield(nmr,'N')
                    ind = nmr.N + 1;
                end
            case 'No'
                nmr = [];
        end
    end
    
%     TEs=fliplr(nmrFitData.TEs);
%     dialog = inputdlg({'Temperatures:'},'Choose temperatures',1,{mat2str(TEs)});
%     if isempty(dialog), return; end;
%     TEs = num2str(dialog{1});
    
    [Path, fileName, Ext] = fileparts(get(handles.edtFN,'String'));
    fileName = extractName(fileName);
        
    for i=1:numel(nmrFitData.TEs{ns})
        TEM = nmrFitData.TEs{ns}(i);
        exp.Temperature = TEM;
        exp.Freq = 100;
        exp.DW = 1;
        nmr.data{ind}.T = [0 0 0];
        nmr.data{ind}.signal = [0 0 0];
        nmr.data{ind}.exp = exp;
        
        nmr.data{ind}.fname = fullfile(Path,[fileName '_' num2str(TEM) 'K' Ext]);
        nmr.temp(ind) = TEM;
        nmr.freq(ind) = 100;
        nmr.dates{ind} = datestr(now);

        f = nmrFitData.data{ns}(:,2*i-1);
        spc = nmrFitData.data{ns}(:,2*i);
        if isfield(nmrFitData,'dataErr')
            spcErr = nmrFitData.dataErr{ns}(:,i);
            spc = complex(spc,spcErr);
        end
        
        if ~isPPM
            nmr.fft{ind}.f = 1e3*f;
        else
            nmr.fft{ind}.f = f;
        end
        nmr.fft{ind}.spc = spc;
        nmr.fft{ind}.TD = numel(f);
        nmr.fft{ind}.shl = 0;
        nmr.fft{ind}.phs = 0;
        nmr.fft{ind}.de = 0;
        nmr.fft{ind}.lb = 0;
        nmr.fft{ind}.xlim = ['[' num2str(round(min(f))) ' ' num2str(round(max(f))) ']'];
        nmr.fft{ind}.isPPM = isPPM;
        nmr.N = ind;
        nmr.fft{ind}.autoshl = 0;
        nmr.fft{ind}.isNorm = 0;
        nmr.fft{ind}.REF = exp.Freq;
        ind = ind + 1;
    end

assignin('base','nmr',nmr);
nmrFit


% --- Executes on button press in butSweep.
function butSweep_Callback(hObject, eventdata, handles)

reloadSettings(handles);
integSW(handles);


% --- Executes on button press in butShowW.
function butShowW_Callback(hObject, eventdata, handles)
FR = evalin('base', 'FR');
CUT = evalin('base', 'CUT');
TAU = evalin('base', 'TAU');
D1 = evalin('base', 'D1');
mode = evalin('base', 'mode');

if evalin('base','exist(''S'',''var'')')
    Weights = evalin('base','S'); 
else
    msgbox('Weights are not yet calculated. Please rerun FFT!','Error');
    return
end

switch mode
    case {'T1','T2'}
        xPar = TAU;
    case 'D1'
        xPar = D1;
    case 'SW'
        xPar = FR;
    case 'SPC'
        xPar = TAU;
    otherwise
        error(['Unknown mode ' mode ' !!!']);
end

a1 = handles.axeSHL;
plot(a1,xPar,Weights,'o',[min(xPar) max(xPar)],[CUT CUT],'g');
ylabel(a1,'S weights');
grid(a1,'on');
switch mode
    case {'T1','T2'}
        set(a1,'XScale','log');
    case {'D1','SW'}
        set(a1,'XScale','lin');
    otherwise
        set(a1,'XScale','lin');
end
if numel(xPar) > 1; 
    xlim(a1,[min(xPar) max(xPar)]); 
end
set(a1,'XTickLabel',[]);


% --------------------------------------------------------------------
function mnuTools_Callback(hObject, eventdata, handles)



% --------------------------------------------------------------------
function mnuScnIndTemp_Callback(hObject, eventdata, handles)
    
    FileName = get(handles.edtFN,'String');
    TE = str2num(get(handles.edtTE,'String'));
    
    [Path Name Ext] = fileparts(FileName);
    Path = [Path '\'];
    [preName, tmpTE, tmpIND] = getNameParts(Name); % Get preName
    
    TEs = evalin('base','TEs');
    scanINDresults = zeros(numel(TEs),4);
    for i=1:numel(TEs)
        scanINDresults(i,1) = TEs(i);
        d = dir([Path preName num2str(TEs(i)) 'K-0*' Ext]);
        scanINDresults(i,1) = TEs(i);
        scanINDresults(i,2) = numel(d);
        maxTau = 0;
        maxIND = 0;
        for j=1:numel(d)
            id = fopen([Path d(j).name]);
            if id == -1
                disp([Path d(j).name ' cannot be found !!!']);
                continue
            end;
            str = fgetl(id);
            while 1-strcmp(str,'[DATA]')
                str = fgetl(id);
                ind = strfind(str,'TAU=');
                if ind==1
                    tau = str2double(str(ind+4:length(str)));
                    if tau > maxTau
                        maxTau = tau;
                        maxIND = j;
                        break;
                    end
                end
            end
            fclose(id);
        end
        scanINDresults(i,3) = maxTau;
        scanINDresults(i,4) = maxIND;
    end
    assignin('base','scanINDresults',scanINDresults);
    msgbox('Done. Results are in scanINDresults variable: TE, IND, maxTau,maxIND');
    
    


% --------------------------------------------------------------------
function Untitled_2_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function mnuOpen_Callback(hObject, eventdata, handles)
    
FilePath = get(handles.analyseFN,'String');
Fp=regexp(FilePath,'\','split');
FileName = Fp{end};
Fp = regexp(FilePath,FileName,'split');
PathName = Fp{1};
cd(PathName)

[FileName,PathName] = uigetfile({'*.nma';'*.*'},'Open');

if ~isequal(FileName, 0)
    cd(PathName)
    load([PathName FileName],'-mat','nma')
    set(handles.edtFN,'String',nma.filename);
    set(handles.txtTemps,'String',nma.TEs);
    set(handles.edtTE,'String',nma.TEM);
    set(handles.edtIND,'String',nma.IND);
    %set(handles.chkINDauto,'Value',str2double(nma.chkIND));
    set(handles.chkINDauto,'Value',1); % problem!!! because of bug many saved settings files corrupted, temporary solution
    set(handles.chkPLOT,'Value',1);
    assignin('base','PLOT',1);
    set(handles.edtFFTparam,'String',nma.FFTparam);
    set(handles.edtAnalyse,'String',nma.Analyse);
    set(handles.analyseFN,'String',[PathName FileName]);
end
    


% --------------------------------------------------------------------
function mnuSave_Callback(hObject, eventdata, handles)
FilePath = get(handles.analyseFN,'String');
Fp=regexp(FilePath,'\','split');
FileName = Fp{end};
Fp = regexp(FilePath,FileName,'split');
PathName = Fp{1};
cd(PathName)

[FileName,PathName] = uiputfile({'*.nma';'*.*'},'Save');
if ~isequal(FileName, 0)
%         cd(PathName)
    nma.filename = get(handles.edtFN,'String');
    nma.TEs = get(handles.txtTemps,'String');
    nma.TEM = get(handles.edtTE,'String');
    nma.IND = get(handles.edtIND,'String');
    nma.chkIND = num2str(get(handles.chkINDauto,'Value'));
    nma.chkPLOT = num2str(get(handles.chkPLOT,'Value'));
    nma.FFTparam = get(handles.edtFFTparam,'String');
    nma.Analyse = get(handles.edtAnalyse,'String');       
    save([PathName FileName],'nma')
    msgbox(['Saved to: ' PathName FileName]);
    set(handles.analyseFN,'String',[PathName FileName]);
end
    


% --------------------------------------------------------------------
function mnuExit_Callback(hObject, eventdata, handles)
% hObject    handle to mnuExit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in butSWrun.
function butSWrun_Callback(hObject, eventdata, handles)
    TEs = fliplr(evalin('base','TEs')); % Must be the same as in Analyse Run

    dialog = inputdlg({'Temperature shift [K^-1]:','Specter scale','Temperatures:'},'Specter plot',1,{'1','1',mat2str(TEs)});
    if isempty(dialog), return; end;
    
    Tshift = str2double(dialog{1});
    A = str2double(dialog{2});
    TEs = eval(dialog{3});

    isPPM = evalin('base','isPPM');
    PLOT = evalin('base','PLOT');

    if evalin('base','exist(''clipboardM'',''var'')')
        clipboardM = evalin('base','clipboardM');
        if size(clipboardM,2)==numel(TEs)*2
            if strcmp(questdlg('Do you want to reload spectra?','Run Spectra','Yes','No','Yes'),'Yes')
                clipboardM = [];
            end
        else
            clipboardM = [];
        end
    else
        clipboardM = [];
    end


    if isempty(clipboardM)
        for j=1:length(TEs)
            gdisp(handles,['T = ' num2str(TEs(j)) 'K']);
            set(handles.edtTE,'String',['[' num2str(TEs(j)) ']']);
            btnLoad_Callback(hObject, eventdata, handles)
            butSweep_Callback(hObject, eventdata, handles)
            SPCsum = evalin('base','SPCsum');

            %SPCsum(:,2) = SPCsum(:,2)+TEs(j)*Tshift;
            if isempty(clipboardM)
                clipboardM(:,2*j-1) = SPCsum(:,1);
                clipboardM(:,2*j) = SPCsum(:,2);
            else
                nrows = numel(SPCsum(:,1));
                cpbrows = numel(clipboardM(:,1));

                if nrows>cpbrows
                    clipboardM = [clipboardM; NaN(nrows-cpbrows,numel(clipboardM(1,:)))];
                end

                clipboardM(1:nrows,2*j-1) = SPCsum(:,1);
                clipboardM(1:nrows,2*j) = SPCsum(:,2);
                if nrows<cpbrows
                    clipboardM(nrows+1:cpbrows,2*j-1) = NaN;
                    clipboardM(nrows+1:cpbrows,2*j) = NaN;
                end
            end
            drawnow();
        end
    end

    clipboardMunscaled = clipboardM;
    minX = 1e10;
    maxX = -1e10;
    for j=1:length(TEs)
        clipboardM(:,2*j) = A*clipboardM(:,2*j)+TEs(j)*Tshift;
        mv = min(clipboardM(:,2*j-1)); % Get X limits
        if mv < minX; minX = mv; end;
        mv = max(clipboardM(:,2*j-1));
        if mv > maxX; maxX = mv; end;
    end

    if PLOT
        a4 = handles.axeRes;
        plot(a4,clipboardM(:,1),clipboardM(:,2))
        [max_val,location] = max(max(clipboardM(:,1),[],2));
        hold(a4,'on');
        max_val = min([maxX max_val]);
        text(max_val,clipboardM(location,2),[num2str(TEs(1)),' K'],'FontSize',10,'FontWeight','bold','Color','b','Parent',a4)
        for j=2:length(TEs)
            plot(a4,clipboardM(:,2*j-1),clipboardM(:,2*j))
            if mod(j,4)==0
                [max_val,location] = max(max(clipboardM(:,2*j-1),[],2));
                max_val = min([maxX max_val]);
                text(max_val,clipboardM(location,2*j),[num2str(TEs(j)),' K'],'FontSize',10,'FontWeight','bold','Color','b','Parent',a4)
            end
        end
        hold(a4,'off');
        xlim(a4,[minX maxX])
        if isPPM
            xlabel(a4,'\nu - \nu_r_e_f [ppm]')
        else
            xlabel(a4,'\nu - \nu_r_e_f [kHz]')
        end
    end

    num2clip(clipboardM);
    assignin('base','clipboardM',clipboardMunscaled);
    nmrFitData.data{1} = clipboardMunscaled;
    nmrFitData.TEs{1} = TEs;
    assignin('base','nmrFitData',nmrFitData);
    gdisp(handles,'Data copied to clipboard. nmrFitData ready.');


% --------------------------------------------------------------------
function mnuCpyComSpc_Callback(hObject, eventdata, handles)
    SPCsum = evalin('base','SPCsum');
    num2clip(SPCsum)


% --------------------------------------------------------------------
function mnuDataInfo_Callback(hObject, eventdata, handles)
    iTEM = evalin('base','TEM');
    iIND = evalin('base','IND');
    iFolder = get(handles.edtFN,'String');

    assignin('base','iFolder',iFolder);
    assignin('base','iIND',iIND(1));
    assignin('base','iTEM',iTEM);

    infoGUI;


% --------------------------------------------------------------------
function mnuStack_Callback(hObject, eventdata, handles)
figure(5)
Dx = 3;
Dy = -3;
for i=1:2:numel(Fall)
    ind = Fall{i} > -60 & Fall{i} < 60;
    plot(Fall{i}(ind)+i*Dx,real(SPCall{i}(ind))+i*Dy,'k');
    hold on
end
% xlim([-100 100])
hold off
% ylim([-30 15])
axis off


% --------------------------------------------------------------------
function mnuCpyFig_Callback(hObject, eventdata, handles)
    set(handles.text6,'BackgroundColor','white');
    set(handles.text1,'BackgroundColor','white');
    set(handles.text4,'BackgroundColor','white');
    set(handles.text5,'BackgroundColor','white');
    set(handles.chkINDauto,'BackgroundColor','white');
    print -dmeta % -dbitmap

    set(handles.text6,'BackgroundColor',get(handles.figure1,'Color'));
    set(handles.text1,'BackgroundColor',get(handles.figure1,'Color'));
    set(handles.text4,'BackgroundColor',get(handles.figure1,'Color'));
    set(handles.text5,'BackgroundColor',get(handles.figure1,'Color'));
    set(handles.chkINDauto,'BackgroundColor',get(handles.figure1,'Color'));
    msgbox('Done!','Copy Figure')
