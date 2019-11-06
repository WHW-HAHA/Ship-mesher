function varargout = scalegui(varargin)
% SCALEGUI MATLAB code for scalegui.fig
%      SCALEGUI, by itself, creates a new SCALEGUI or raises the existing
%      singleton*.
%
%      H = SCALEGUI returns the handle to a new SCALEGUI or the handle to
%      the existing singleton*.
%
%      SCALEGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SCALEGUI.M with the given input arguments.
%
%      SCALEGUI('Property','Value',...) creates a new SCALEGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before scalegui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to scalegui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help scalegui

% Last Modified by GUIDE v2.5 15-Sep-2015 10:12:18

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @scalegui_OpeningFcn, ...
                   'gui_OutputFcn',  @scalegui_OutputFcn, ...
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


% --- Executes just before scalegui is made visible.
function scalegui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to scalegui (see VARARGIN)

% Choose default command line output for scalegui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes scalegui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = scalegui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_file.
function pushbutton_file_Callback(hObject, eventdata, handles) % select file
    %scale structure of no scaling
    sf.L = 1;
    sf.B = 1;
    sf.D = 1;
    sf.reduc = 0;

    %select and load file
    [file, path] = uigetfile('*.*');
    ship = load_patranfile([path file],sf);
    handles.file = [path file];
    handles.ship = ship;
	guidata(hObject,handles);
        
     %fill in main dimensions
    set(handles.text_selectedfile,'String',file)
    set(handles.text_orlength,'String',ship.L)
    set(handles.text_orbreadth,'String',ship.B)
    set(handles.text_ordraught,'String',ship.T)
    set(handles.text_ordepth,'String',ship.T)
    
    %enable edit boxes
    set(handles.edit_length,'Enable','on')
    set(handles.edit_breadth,'Enable','on')
    set(handles.edit_draught,'Enable','on')
    set(handles.edit_depth,'Enable','on')
    set(handles.edit_lengthscale,'Enable','on')
    set(handles.edit_breadthscale,'Enable','on')
    set(handles.edit_depthscale,'Enable','on')
    set(handles.edit_draughtscale,'Enable','on')
    
    %update info
    set(handles.text_info,'String','Please complete the scaling details.')   
% hObject    handle to pushbutton_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit_length_Callback(hObject, eventdata, handles)
    newL = str2double(get(handles.edit_length, 'string'));
    oldL = str2double(get(handles.text_orlength, 'string'));
    set(handles.edit_lengthscale,'String',newL/oldL);
    
    if iscomplete(handles)
        set(handles.pushbutton_scale,'Enable','on')
        set(handles.text_info,'String','Press scale to compute scaled variables.')  
    else
        set(handles.pushbutton_scale,'Enable','off')
        set(handles.text_info,'String','Please complete the scaling details.')
    end
% hObject    handle to edit_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_length as text
%        str2double(get(hObject,'String')) returns contents of edit_length as a double


% --- Executes during object creation, after setting all properties.
function edit_length_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_breadth_Callback(hObject, eventdata, handles)
    newB = str2double(get(handles.edit_breadth, 'string'));
    oldB = str2double(get(handles.text_orbreadth, 'string'));
    set(handles.edit_breadthscale,'String',newB/oldB);
    
    if iscomplete(handles)
        set(handles.pushbutton_scale,'Enable','on')
        set(handles.text_info,'String','Press scale to compute scaled variables.')  
    else
        set(handles.pushbutton_scale,'Enable','off')
        set(handles.text_info,'String','Please complete the scaling details.')
    end
% hObject    handle to edit_breadth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_breadth as text
%        str2double(get(hObject,'String')) returns contents of edit_breadth as a double


% --- Executes during object creation, after setting all properties.
function edit_breadth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_breadth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_depth_Callback(hObject, eventdata, handles)
    newD = str2double(get(handles.edit_depth, 'string'));
    oldD = str2double(get(handles.text_ordepth, 'string'));
    set(handles.edit_depthscale,'String',newD/oldD);
    
    if iscomplete(handles)
        set(handles.pushbutton_scale,'Enable','on')
        set(handles.text_info,'String','Press scale to compute scaled variables.')  
    else
        set(handles.pushbutton_scale,'Enable','off')
        set(handles.text_info,'String','Please complete the scaling details.')
    end
% hObject    handle to edit_depth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_depth as text
%        str2double(get(hObject,'String')) returns contents of edit_depth as a double


% --- Executes during object creation, after setting all properties.
function edit_depth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_depth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_draught_Callback(hObject, eventdata, handles)
    newT = str2double(get(handles.edit_draught, 'string'));
    oldT = str2double(get(handles.text_ordraught, 'string'));
    set(handles.edit_draughtscale,'String',newT/oldT);
    
    if iscomplete(handles)
        set(handles.pushbutton_scale,'Enable','on')
        set(handles.text_info,'String','Press scale to compute scaled variables.')  
    else
        set(handles.pushbutton_scale,'Enable','off')
        set(handles.text_info,'String','Please complete the scaling details.')
    end
% hObject    handle to edit_draught (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_draught as text
%        str2double(get(hObject,'String')) returns contents of edit_draught as a double


% --- Executes during object creation, after setting all properties.
function edit_draught_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_draught (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_lengthscale_Callback(hObject, eventdata, handles)
    scale = str2double(get(handles.edit_lengthscale,'String'));
    oldL = str2double(get(handles.text_orlength, 'string'));
    set(handles.edit_length,'String',oldL*scale);
    
    if iscomplete(handles)
        set(handles.pushbutton_scale,'Enable','on')
        set(handles.text_info,'String','Press scale to compute scaled variables.')  
    else
        set(handles.pushbutton_scale,'Enable','off')
        set(handles.text_info,'String','Please complete the scaling details.')
    end
% hObject    handle to edit_lengthscale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_lengthscale as text
%        str2double(get(hObject,'String')) returns contents of edit_lengthscale as a double


% --- Executes during object creation, after setting all properties.
function edit_lengthscale_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_lengthscale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_breadthscale_Callback(hObject, eventdata, handles)
    scale = str2double(get(handles.edit_breadthscale,'String'));
    oldB = str2double(get(handles.text_orbreadth, 'string'));
    set(handles.edit_breadth,'String',oldB*scale);
    
    if iscomplete(handles)
        set(handles.pushbutton_scale,'Enable','on')
        set(handles.text_info,'String','Press scale to compute scaled variables.')  
    else
        set(handles.pushbutton_scale,'Enable','off')
        set(handles.text_info,'String','Please complete the scaling details.')
    end
% hObject    handle to edit_breadthscale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_breadthscale as text
%        str2double(get(hObject,'String')) returns contents of edit_breadthscale as a double


% --- Executes during object creation, after setting all properties.
function edit_breadthscale_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_breadthscale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_depthscale_Callback(hObject, eventdata, handles)
    scale = str2double(get(handles.edit_depthscale,'String'));
    oldD = str2double(get(handles.text_ordepth, 'string'));
    set(handles.edit_depth,'String',oldD*scale);

    if iscomplete(handles)
        set(handles.pushbutton_scale,'Enable','on')
        set(handles.text_info,'String','Press scale to compute scaled variables.')  
    else
        set(handles.pushbutton_scale,'Enable','off')
        set(handles.text_info,'String','Please complete the scaling details.')
    end    
    % hObject    handle to edit_depthscale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_depthscale as text
%        str2double(get(hObject,'String')) returns contents of edit_depthscale as a double


% --- Executes during object creation, after setting all properties.
function edit_depthscale_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_depthscale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_draughtscale_Callback(hObject, eventdata, handles)
    scale = str2double(get(handles.edit_draughtscale,'String'));
    oldT = str2double(get(handles.text_ordraught, 'string'));
    set(handles.edit_draught,'String',oldT*scale);
    
    if iscomplete(handles)
        set(handles.pushbutton_scale,'Enable','on')
        set(handles.text_info,'String','Press scale to compute scaled variables.')  
    else
        set(handles.pushbutton_scale,'Enable','off')
        set(handles.text_info,'String','Please complete the scaling details.')
    end
% hObject    handle to edit_draughtscale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_draughtscale as text
%        str2double(get(hObject,'String')) returns contents of edit_draughtscale as a double


% --- Executes during object creation, after setting all properties.
function edit_draughtscale_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_draughtscale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_scale.
function pushbutton_scale_Callback(hObject, eventdata, handles)
    sf.L = str2double(get(handles.edit_lengthscale,'String'));
    sf.B = str2double(get(handles.edit_breadthscale,'String'));
    sf.D = str2double(get(handles.edit_depthscale,'String'));
    sf.reduc = str2double(get(handles.edit_depth,'String')) - str2double(get(handles.edit_draught,'String'));
    
    handles.sf = sf;
	guidata(hObject,handles);
    ship = handles.ship;
    scaled = load_patranfile(handles.file,sf);
        
    set(handles.text_originalLB,'String',round(ship.L/ship.B,2))
    set(handles.text_originalBT,'String',round(ship.B/ship.T,2))
    set(handles.text_originalCB,'String',round(ship.displacement/(ship.L*ship.B*ship.T),2))
    set(handles.text_originalDisp,'String',round(ship.displacement,0))
    
    set(handles.text_scaledLB,'String',round(scaled.L/scaled.B,2))
    set(handles.text_scaledBT,'String',round(scaled.B/scaled.T,2))
    set(handles.text_scaledCB,'String',round(scaled.displacement/(scaled.L*scaled.B*scaled.T),2))
    set(handles.text_scaledDisp,'String',round(scaled.displacement,0))   
    
    set(handles.text_factorLB,'String',round((scaled.L/scaled.B)/(ship.L/ship.B),2))
    set(handles.text_factorBT,'String',round((scaled.B/scaled.T)/(ship.B/ship.T),2))
    set(handles.text_factorCB,'String',round((scaled.displacement/(scaled.L*scaled.B*scaled.T))/(ship.displacement/(ship.L*ship.B*ship.T)),2))
    set(handles.text_factorDisp,'String',round(scaled.displacement/ship.displacement,2))
    
    set(handles.pushbutton_mesh,'Enable','on')
    set(handles.text_info,'String','Review the output and press [Start Meshing] to accept.')

    handles.scaled = scaled;
    guidata(hObject,handles);
% hObject    handle to pushbutton_scale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function [ yes ] = iscomplete(handles)
    yes = 1;
    if isnan(str2double(get(handles.edit_lengthscale,'String')))
        yes = 0;
    elseif isnan(str2double(get(handles.edit_breadthscale,'String')))
        yes = 0;
    elseif isnan(str2double(get(handles.edit_depthscale,'String')))
        yes = 0;
    elseif isnan(str2double(get(handles.edit_draughtscale,'String')))
        yes = 0;
    end


% --- Executes on button press in pushbutton_mesh.
function pushbutton_mesh_Callback(hObject, eventdata, handles)
    assignin('base','file',handles.file)
    assignin('base','sf',handles.sf)
    assignin('base','ship',handles.scaled)
    delete(gcf)
    run meshgui.m

% hObject    handle to pushbutton_mesh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
