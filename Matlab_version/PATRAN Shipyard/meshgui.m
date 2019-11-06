function varargout = meshgui(varargin)
%% Done...
%MESHGUI M-file for meshgui.fig
%      MESHGUI, by itself, creates a new MESHGUI or raises the existing
%      singleton*.
%
%      H = MESHGUI returns the handle to a new MESHGUI or the handle to
%      the existing singleton*.
%
%      MESHGUI('Property','Value',...) creates a new MESHGUI using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to meshgui_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      MESHGUI('CALLBACK') and MESHGUI('CALLBACK',hObject,...) call the
%      local function named CALLBACK in MESHGUI.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help meshgui

% Last Modified by GUIDE v2.5 13-Aug-2015 13:54:53
% 
% 
% Revised:
%   Yijun 01052016: GM_l value
%                   kzz



% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @meshgui_OpeningFcn, ...
                   'gui_OutputFcn',  @meshgui_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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


% --- Executes just before meshgui is made visible.
function meshgui_OpeningFcn(hObject, eventdata, handles, varargin)
    t = table(['Lpp  ';'B    ';'D    ';'T    ';'Disp.';'wn_3 '],[' ';' ';' ';' ';' ';' '],...
    ['A_w ';'GM_l';'GM_t';'I_l ';'I_t ';'wn_4'],[' ';' ';' ';' ';' ';' '],...
    ['Mass';'KG  ';'k_xx';'k_yy';'k_zz';'wn_5'],[' ';' ';' ';' ';' ';' ']);
    set(handles.uitable2,'Data',table2cell(t));
    handles.plan = evalin('base','plan');
    set(handles.edit_kxx,'String',handles.plan.B*0.25);
    set(handles.edit_kyy,'String',handles.plan.L*0.3);
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for meshgui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes meshgui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = meshgui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit_kg_Callback(hObject, eventdata, handles)
    iscomplete(handles);
% hObject    handle to edit_kg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_kg as text
%        str2double(get(hObject,'String')) returns contents of edit_kg as a double


% --- Executes during object creation, after setting all properties.
function edit_kg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_kg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_kxx_Callback(hObject, eventdata, handles)
    iscomplete(handles)
% hObject    handle to edit_kxx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_kxx as text
%        str2double(get(hObject,'String')) returns contents of edit_kxx as a double


% --- Executes during object creation, after setting all properties.
function edit_kxx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_kxx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_kyy_Callback(hObject, eventdata, handles)
    iscomplete(handles)
% hObject    handle to edit_kyy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_kyy as text
%        str2double(get(hObject,'String')) returns contents of edit_kyy as a double


% --- Executes during object creation, after setting all properties.
function edit_kyy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_kyy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_panelsize_Callback(hObject, eventdata, handles)
% hObject    handle to edit_panelsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_panelsize as text
%        str2double(get(hObject,'String')) returns contents of edit_panelsize as a double


% --- Executes during object creation, after setting all properties.
function edit_panelsize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_panelsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_compute.
function pushbutton_compute_Callback(hObject, eventdata, handles)
    KG = str2double(get(handles.edit_kg,'String'));
    kxx = str2double(get(handles.edit_kxx,'String'));
    kyy = str2double(get(handles.edit_kyy,'String'));
    plan = handles.plan;
    
    info = linesplan_panelsize(plan,KG,kxx,kyy,str2double(get(handles.edit_waterdepth,'String')));
    if info.panelsize == plan.T/4
        set(handles.text_info,'String','Panelsize set to 0.25 * draught.')
    else
        set(handles.text_info,'String','Panelsize set to 0.25 * smallest wavelength.')
    end
    handles.t = table(['Lpp  ';'B    ';'D    ';'T    ';'Disp.';'wn_3 '],[plan.L;plan.B;plan.D;plan.T;info.displacement;info.wn_3],...
    ['A_w ';'GM_l';'GM_t';'I_l ';'I_t ';'wn_4'],[info.Aw;info.GM_l;info.GM_t;info.I_l;info.I_t;info.wn_4],...
    ['Mass';'KG  ';'k_xx';'k_yy';'k_zz';'wn_5'],[info.displacement*1.025;KG;kxx;kyy;kyy;info.wn_5]);
    set(handles.uitable2,'Data',table2cell(handles.t))
    
	guidata(hObject,handles);
    set(handles.edit_panelsize,'String',info.panelsize);
    set(handles.edit_panelsize,'Enable','on');
    set(handles.pushbutton_mesh,'Enable','on');
    set(handles.pushbutton_file,'Enable','on');
    
% hObject    handle to pushbutton_compute (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton_mesh.
function pushbutton_mesh_Callback(hObject, eventdata, handles)
    if get(handles.checkbox1,'Value')==1
        mesh_linesplan2(handles.plan,str2double(get(handles.edit_panelsize,'String')),get(handles.slider_spacing,'Value')*pi);
    else
        mesh_linesplan(handles.plan,str2double(get(handles.edit_panelsize,'String')),get(handles.slider_spacing,'Value')*pi);
    end
% hObject    handle to pushbutton_mesh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in pushbutton_file.

function pushbutton_file_Callback(hObject, eventdata, handles)
    [pathstr, name] = fileparts(handles.plan.filepath); 
    writetable(handles.t,[name '_info.csv'],'WriteVariableNames',false)
    set(handles.text_info,'String',['File exported as ' handles.plan.type '_info.csv'])
% hObject    handle to pushbutton_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function [ yes ] = iscomplete(handles)
    yes = 1;
    if isnan(str2double(get(handles.edit_kg,'String')))
        yes = 0;
    elseif isnan(str2double(get(handles.edit_kxx,'String')))
        yes = 0;
    elseif isnan(str2double(get(handles.edit_kyy,'String')))
        yes = 0;
    end
    if yes
        set(handles.pushbutton_compute,'Enable','on');
        set(handles.text_info,'String','Press [Compute] to update the output table.')
    else
        set(handles.pushbutton_compute,'Enable','off');
        set(handles.text_info,'String','Please input the missing properties.')
    end


% --- Executes on slider movement.
function slider_spacing_Callback(hObject, eventdata, handles)
% hObject    handle to slider_spacing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider_spacing_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_spacing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit_waterdepth_Callback(hObject, eventdata, handles)
% hObject    handle to edit_waterdepth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_waterdepth as text
%        str2double(get(hObject,'String')) returns contents of edit_waterdepth as a double


% --- Executes during object creation, after setting all properties.
function edit_waterdepth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_waterdepth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1
