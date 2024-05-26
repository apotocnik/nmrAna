function varargout = infoGUI(varargin)
% INFOGUI MATLAB code for infoGUI.fig
%      INFOGUI, by itself, creates a new INFOGUI or raises the existing
%      singleton*.
%
%      H = INFOGUI returns the handle to a new INFOGUI or the handle to
%      the existing singleton*.
%
%      INFOGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in INFOGUI.M with the given input arguments.
%
%      INFOGUI('Property','Value',...) creates a new INFOGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before infoGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to infoGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help infoGUI

% Last Modified by GUIDE v2.5 09-May-2012 09:54:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @infoGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @infoGUI_OutputFcn, ...
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

% --- Executes just before infoGUI is made visible.
function infoGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to infoGUI (see VARARGIN)

% Choose default command line output for infoGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes infoGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = infoGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

fillData(handles);
butLoad_Callback(hObject, eventdata, handles);

function fillData(handles)
IND = evalin('base','iIND');
TE = evalin('base','iTEM');
Folder = evalin('base','iFolder');

set(handles.edtINDbatch,'String',num2str(IND));
set(handles.edtTEbatch,'String',num2str(TE));
set(handles.edtTE,'String',num2str(TE));
set(handles.edtIND,'String',num2str(IND));
set(handles.txtFolder,'String',Folder);

function edtINDbatch_Callback(hObject, eventdata, handles)
% hObject    handle to edtINDbatch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtINDbatch as text
%        str2double(get(hObject,'String')) returns contents of edtINDbatch as a double


% --- Executes during object creation, after setting all properties.
function edtINDbatch_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtINDbatch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edtStep_Callback(hObject, eventdata, handles)
% hObject    handle to edtStep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtStep as text
%        str2double(get(hObject,'String')) returns contents of edtStep as a double


% --- Executes during object creation, after setting all properties.
function edtStep_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtStep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edtEnd_Callback(hObject, eventdata, handles)
% hObject    handle to edtEnd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtEnd as text
%        str2double(get(hObject,'String')) returns contents of edtEnd as a double


% --- Executes during object creation, after setting all properties.
function edtEnd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtEnd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in butLoad.
function butLoad_Callback(hObject, eventdata, handles)
% hObject    handle to butLoad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
FileName = get(handles.txtFolder,'String');
[Path Name Ext] = fileparts(FileName);
Path = [Path '\'];
TE = get(handles.edtTE,'String');
IND = get(handles.edtIND,'String');
gIND = IND;
prefix = get(handles.edtPrefix,'String');

while numel(IND)<3
    IND = ['0' IND];
end

[preName, tmpTE, tmpIND] = getNameParts(Name); % Get preName

if strcmp(IND, '000')
    filepath = [Path preName TE 'K' Ext];
else
    filepath = [Path preName TE 'K-' IND Ext];
end
id = fopen(filepath);
        if id == -1
            error([filepath ' cannot be found !!!'])
        end
        textline = fgetl(id);
        lines = {};
        set(handles.edtFile,'String','');
        while 1-strcmp(textline,'[DATA]') &&  feof(id) == 0
            lines = [lines; {textline}];
            
            if ~isempty(strfind(textline,'='))
                ind = strfind(textline,'=');
                nmchk = textline(ind+1:end);
                if ~isempty(str2num(nmchk)) && isempty(strfind(nmchk,':'))
                    textline = strrep(textline,'=','=[');
                   evalin('base',[prefix textline '];']); 
                else
                   %textline = strrep(textline,'''','''''');
                   textline=strtrim(textline);
                   if ~(textline(1)=='"')
                    textline = strrep(textline,'=','=''');
                    evalin('base',[prefix textline ''';']);
                   end
                end
            end
            textline = fgetl(id);
        end
        fclose(id);
        evalin('base',[prefix 'IND=' num2str(gIND) ';']);
        evalin('base',[prefix 'TE=' num2str(TE) ';']);
        set(handles.edtFile,'String',lines);
%notes:
fileID = fopen([Path 'notes.txt']);
if fileID ~= -1
    tline = fgetl(fileID);
    lines = {};
    while ischar(tline)
        lines = [lines; tline];
        tline = fgetl(fileID);
    end
    set(handles.edtNotes,'String',lines);
    fclose(fileID);
else
    set(handles.edtNotes,'String','');
end

function edtTEbatch_Callback(hObject, eventdata, handles)
% hObject    handle to edtTEbatch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtTEbatch as text
%        str2double(get(hObject,'String')) returns contents of edtTEbatch as a double


% --- Executes during object creation, after setting all properties.
function edtTEbatch_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtTEbatch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in butAllTE.
function butAllTE_Callback(hObject, eventdata, handles)
% hObject    handle to butAllTE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
FileName = get(handles.txtFolder,'String');
[Path Name Ext] = fileparts(FileName);
Path = [Path '\'];

IND = get(handles.edtINDbatch,'String');

while numel(IND)<3
    IND = ['0' IND];
end

[preName, tmpTE, tmpIND] = getNameParts(Name); % Get preName
if strcmp(IND,'000')
    d = dir([Path preName '*K' Ext]);
else
    d = dir([Path preName '*K-' IND Ext]);
end
    TEs = zeros(1,numel(d));
    for i=1:numel(d)
        [tmpPath Name Ext] = fileparts(d(i).name);
        [preName, tmpTE, tmpIND] = getNameParts(Name);
        if isempty(tmpTE), continue; end; % wrong filename: 150K-G ending
        TEs(i) = tmpTE;
    end
    TEs = sort(TEs);
    set(handles.edtTEbatch,'String',num2str(TEs));

% --- Executes on button press in butAllIND.
function butAllIND_Callback(hObject, eventdata, handles)
% hObject    handle to butAllIND (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
FileName = get(handles.txtFolder,'String');
[Path Name Ext] = fileparts(FileName);
Path = [Path '\'];

TE = get(handles.edtTEbatch,'String');


[preName, tmpTE, tmpIND] = getNameParts(Name); % Get preName

    d = dir([Path preName TE 'K-0*' Ext]);
    IND = zeros(1,numel(d));
    for i=1:numel(d)
        [tmpPath Name Ext] = fileparts(d(i).name);
        [preName, tmpTE, tmpIND] = getNameParts(Name);
        IND(i) = tmpIND;
    end
    IND = sort(IND);
    if isempty(IND)
        set(handles.edtINDbatch,'String','0');
    else
        set(handles.edtINDbatch,'String',num2str(IND));
    end


function edtPrefix_Callback(hObject, eventdata, handles)
% hObject    handle to edtPrefix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtPrefix as text
%        str2double(get(hObject,'String')) returns contents of edtPrefix as a double


% --- Executes during object creation, after setting all properties.
function edtPrefix_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtPrefix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double


% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in butTERun.
function butTERun_Callback(hObject, eventdata, handles)
% hObject    handle to butTERun (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in butINDRun.
function butINDRun_Callback(hObject, eventdata, handles)
% hObject    handle to butINDRun (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edtTE_Callback(hObject, eventdata, handles)
% hObject    handle to edtTE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtTE as text
%        str2double(get(hObject,'String')) returns contents of edtTE as a double


% --- Executes during object creation, after setting all properties.
function edtTE_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtTE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edtIND_Callback(hObject, eventdata, handles)
% hObject    handle to edtIND (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtIND as text
%        str2double(get(hObject,'String')) returns contents of edtIND as a double


% --- Executes during object creation, after setting all properties.
function edtIND_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtIND (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edtFile_Callback(hObject, eventdata, handles)
% hObject    handle to edtFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtFile as text
%        str2double(get(hObject,'String')) returns contents of edtFile as a double


% --- Executes during object creation, after setting all properties.
function edtFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edtNotes_Callback(hObject, eventdata, handles)
% hObject    handle to edtNotes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtNotes as text
%        str2double(get(hObject,'String')) returns contents of edtNotes as a double


% --- Executes during object creation, after setting all properties.
function edtNotes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtNotes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in butSave.
function butSave_Callback(hObject, eventdata, handles)
% hObject    handle to butSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
FileName = get(handles.txtFolder,'String');
[Path Name Ext] = fileparts(FileName);
Path = [Path '\'];

fileID = fopen([Path 'notes.txt'],'w');
param = get(handles.edtNotes,'String');
fprintf(fileID,sprintf('%s\r\n',param{:}));
fclose(fileID);

% --- Executes on button press in butRun.
function butRun_Callback(hObject, eventdata, handles)
% hObject    handle to butRun (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
TE = str2num(get(handles.edtTEbatch,'String'));
IND = str2num(get(handles.edtINDbatch,'String'));
TE2 = TE; IND2 = IND;

if hObject == handles.butRun
   INDall = 0; 
else
   INDall = 1;
end

if numel(TE)>1 && numel(IND)>1
   errordlg('There should be exactly one value in IND or TE field!','Batch run');
   return; 
end

param = get(handles.edtBegin,'String');
evalin('base',sprintf('%s\n',param{:}));

for i=1:numel(TE)
    if INDall
       set(handles.edtTEbatch,'String',num2str(TE(i)));
       butAllIND_Callback(hObject, eventdata, handles);
       IND = str2num(get(handles.edtINDbatch,'String'));
    end
    for j=1:numel(IND)
        set(handles.edtTE,'String',num2str(TE(i)));
        set(handles.edtIND,'String',num2str(IND(j)));
        butLoad_Callback(hObject, eventdata, handles);
        param = get(handles.edtStep,'String');
        evalin('base',sprintf('%s\n',param{:}));
    end
end
        param = get(handles.edtEnd,'String');
        evalin('base',sprintf('%s\n',param{:}));
        
set(handles.edtTEbatch,'String',num2str(TE2));
set(handles.edtINDbatch,'String',num2str(IND2));

% --- Executes on button press in butRunAll.
function butRunAll_Callback(hObject, eventdata, handles)
% hObject    handle to butRunAll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
butRun_Callback(hObject, eventdata, handles);



function edtBegin_Callback(hObject, eventdata, handles)
% hObject    handle to edtBegin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtBegin as text
%        str2double(get(hObject,'String')) returns contents of edtBegin as a double


% --- Executes during object creation, after setting all properties.
function edtBegin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtBegin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in butOpenScript.
function butOpenScript_Callback(hObject, eventdata, handles)
% hObject    handle to butOpenScript (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
FilePath = get(handles.txtScript,'String');
Fp=regexp(FilePath,'\','split');
FileName = Fp{end};
Fp = regexp(FilePath,FileName,'split');
PathName = Fp{1};
cd(PathName)

if hObject == handles.butOpenScript
    [FileName,PathName] = uigetfile({'*.nmi';'*.*'},'Open');
end
    if ~isequal(FileName, 0)
        cd(PathName)
        load([PathName FileName],'-mat','nmi')
        set(handles.edtBegin,'String',nmi.begin);
        set(handles.edtStep,'String',nmi.step);
        set(handles.edtEnd,'String',nmi.end);
        set(handles.txtScript,'String',[PathName FileName]);
    end

% --- Executes on button press in butSaveScript.
function butSaveScript_Callback(hObject, eventdata, handles)
% hObject    handle to butSaveScript (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
FilePath = get(handles.txtScript,'String');
Fp=regexp(FilePath,'\','split');
FileName = Fp{end};
Fp = regexp(FilePath,FileName,'split');
PathName = Fp{1};
cd(PathName)

    [FileName,PathName] = uiputfile({'*.nmi';'*.*'},'Save');
    if ~isequal(FileName, 0)
        nmi.begin = get(handles.edtBegin,'String');
        nmi.step = get(handles.edtStep,'String');
        nmi.end = get(handles.edtEnd,'String');  
        save([PathName FileName],'nmi')
        msgbox(['Saved to: ' PathName FileName]);
        set(handles.txtScript,'String',[PathName FileName]);
    end

% --- Executes on button press in butReloadScript.
function butReloadScript_Callback(hObject, eventdata, handles)
% hObject    handle to butReloadScript (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
butOpenScript_Callback(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function txtFolder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtFolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
