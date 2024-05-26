function gdisp(handles,textline)
% Write text into the Log text box in nmrAna
    txt = get(handles.edtLog,'String');
%     while numel(txt) > 100
%         set(handles.edtLog,'String','');
%         txt='';
%     end
    
    while numel(txt) > 100
        txt(end) = [];
    end
    set(handles.edtLog,'String',[{textline}; txt]);
