function varargout = OPT_PMOML(varargin)
% OPT_PMOML MATLAB code for OPT_PMOML.fig
%      OPT_PMOML, by itself, creates a_opt new OPT_PMOML or raises the existing
%      singleton*.
%
%      H = OPT_PMOML returns the handle to a_opt new OPT_PMOML or the handle to
%      the existing singleton*.
%
%      OPT_PMOML('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in OPT_PMOML.M with the given input arguments.
%
%      OPT_PMOML('Property','Value',...) creates a_opt new OPT_PMOML or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before OPT_PMOML_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to OPT_PMOML_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's_opt Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help OPT_PMOML

% Last Modified by GUIDE v2.5 15-May-2020 14:19:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @OPT_PMOML_OpeningFcn, ...
                   'gui_OutputFcn',  @OPT_PMOML_OutputFcn, ...
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


% --- Executes just before OPT_PMOML is made visible.
function OPT_PMOML_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a_opt future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to OPT_PMOML (see VARARGIN)

% Choose default command line output for OPT_PMOML
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes OPT_PMOML wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = OPT_PMOML_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a_opt future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in run_button.
function run_button_Callback(hObject, eventdata, handles)
% hObject    handle to run_button (see GCBO)
% eventdata  reserved - to be defined in a_opt future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get user input from GUI
a_min = str2double(get(handles.a_min,'String'))*1e-6;
a_max = str2double(get(handles.a_max,'String'))*1e-6;
s_min = str2double(get(handles.s_min,'String'))*1e-6;
s_max = str2double(get(handles.s_max,'String'))*1e-6;
p_min = str2double(get(handles.p_min,'String'))*1e-6;
p_max = str2double(get(handles.p_max,'String'))*1e-6;
%
t1 = str2double(get(handles.t1,'String'))*1e-6;
t2 = str2double(get(handles.t2,'String'))*1e-6;
g = str2double(get(handles.g,'String'))*1e-6;
ep1 = str2double(get(handles.ep1,'String'));
ep2 = str2double(get(handles.ep2,'String'));
%
M = str2double(get(handles.M,'String'));
N = str2double(get(handles.N,'String'));
Ndot = str2double(get(handles.Ndot,'String'));
Coef_fric = str2double(get(handles.Coef_fric,'String'));

% Optimization
d_opt = (t1+t2+g)^2*ep1./((t1+t2)*ep1+g*ep2);
set(handles.d_opt,'string',num2str(d_opt*1e6,4));
p_opt = 1.263*d_opt+0.667*1e-6;
if p_opt>p_max
    p_opt = p_max;
elseif p_opt<p_min
    p_opt = p_min;
end
set(handles.p_opt,'string',num2str(p_opt*1e6,4));
% ratio_sp = s_min/p_opt;
s_opt = s_min;
a_opt = p_opt-s_opt;
Epslon = ep1;

%
set(handles.a_opt,'string',num2str(a_opt*1e6,4));
set(handles.s_opt,'string',num2str(s_opt*1e6,4));
set(handles.Epslon,'string',num2str(Epslon));

% Calculate by PMOML
[MThetax,MPhi,MFx]=GUIFunc_Fx_PMOML(a_opt,d_opt,s_opt,N,M,Ndot,Coef_fric);
MFx = Epslon*MFx;
[minfx1,minind1] = min(MFx);
[min_Fx,minind2] = min(minfx1);
set(handles.min_Fx,'string',num2str(min_Fx*1e3,4));
set(handles.min_thetax,'string',num2str(MThetax(minind1(minind2),minind2)/pi,2));
set(handles.min_phi,'string',num2str(MPhi(minind1(minind2),minind2)/pi,2));
[maxfx1,maxind1] = max(MFx);
[max_Fx,maxind2] = max(maxfx1);
set(handles.max_Fx,'string',num2str(max_Fx*1e3,4));
set(handles.max_thetax,'string',num2str(MThetax(maxind1(maxind2),maxind2)/pi,2));
set(handles.max_phi,'string',num2str(MPhi(maxind1(maxind2),maxind2)/pi,2));

% [x,y] = meshgrid(-2*p_opt:0.1*p_opt:2*p_opt,-1*d_opt:0.1*d_opt:1*d_opt);
% [Fcoulombian,x,y,U]=GUIFunc_MatrixChargeCoulF_PMOML(a_opt,s_opt,d_opt,dx,phi,N,M);
% [Ex,Ey]=gradient(U);
% Ex = Ex./(0.1*p_opt);
% Ey = Ey./(0.1*d_opt);

% Plot
surf(handles.axes_Fx,MThetax,MPhi,MFx*1e3);
xlabel(handles.axes_Fx,'\theta_x')
ylabel(handles.axes_Fx,'\phi')
zlabel(handles.axes_Fx,'Fx (mN/(V^2\cdot m^2))')
% c = colorbar;
% c.Label.String = 'Fx (mN/(V^2\`m^2))'
box on


% 



function a_opt_Callback(hObject, eventdata, handles)
% hObject    handle to a_opt (see GCBO)
% eventdata  reserved - to be defined in a_opt future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of a_opt as text
%        str2double(get(hObject,'String')) returns contents of a_opt as a_opt double
% Validate that the text in the f2 field converts to a_opt real number




function s_opt_Callback(hObject, eventdata, handles)
% hObject    handle to s_opt (see GCBO)
% eventdata  reserved - to be defined in a_opt future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of s_opt as text
%        str2double(get(hObject,'String')) returns contents of s_opt as a_opt double
% Validate that the text in the f2 field converts to a_opt real number




function p_opt_Callback(hObject, eventdata, handles)
% hObject    handle to p_opt (see GCBO)
% eventdata  reserved - to be defined in a_opt future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of p_opt as text
%        str2double(get(hObject,'String')) returns contents of p_opt as a_opt double
% Validate that the text in the f2 field converts to a_opt real number



function M_Callback(hObject, eventdata, handles)
% hObject    handle to M (see GCBO)
% eventdata  reserved - to be defined in a_opt future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of M as text
%        str2double(get(hObject,'String')) returns contents of M as a_opt double
% Validate that the text in the f2 field converts to a_opt real number

function N_Callback(hObject, eventdata, handles)
% hObject    handle to N (see GCBO)
% eventdata  reserved - to be defined in a_opt future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of N as text
%        str2double(get(hObject,'String')) returns contents of N as a_opt double


function Coef_fric_Callback(hObject, eventdata, handles)
% hObject    handle to Coef_fric (see GCBO)
% eventdata  reserved - to be defined in a_opt future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Coef_fric as text
%        str2double(get(hObject,'String')) returns contents of Coef_fric as a_opt double




function d_opt_Callback(hObject, eventdata, handles)
% hObject    handle to d_opt (see GCBO)
% eventdata  reserved - to be defined in a_opt future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of d_opt as text
%        str2double(get(hObject,'String')) returns contents of d_opt as a_opt double




function Ndot_Callback(hObject, eventdata, handles)
% hObject    handle to Ndot (see GCBO)
% eventdata  reserved - to be defined in a_opt future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Ndot as text
%        str2double(get(hObject,'String')) returns contents of Ndot as a_opt double



function a_min_Callback(hObject, eventdata, handles)
% hObject    handle to a_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of a_min as text
%        str2double(get(hObject,'String')) returns contents of a_min as a double



function a_max_Callback(hObject, eventdata, handles)
% hObject    handle to a_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of a_max as text
%        str2double(get(hObject,'String')) returns contents of a_max as a double



function s_min_Callback(hObject, eventdata, handles)
% hObject    handle to s_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of s_min as text
%        str2double(get(hObject,'String')) returns contents of s_min as a double



function s_max_Callback(hObject, eventdata, handles)
% hObject    handle to s_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of s_max as text
%        str2double(get(hObject,'String')) returns contents of s_max as a double



function p_min_Callback(hObject, eventdata, handles)
% hObject    handle to p_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of p_min as text
%        str2double(get(hObject,'String')) returns contents of p_min as a double



function p_max_Callback(hObject, eventdata, handles)
% hObject    handle to p_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of p_max as text
%        str2double(get(hObject,'String')) returns contents of p_max as a double



function t1_Callback(hObject, eventdata, handles)
% hObject    handle to t1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of t1 as text
%        str2double(get(hObject,'String')) returns contents of t1 as a double



function ep1_Callback(hObject, eventdata, handles)
% hObject    handle to ep1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ep1 as text
%        str2double(get(hObject,'String')) returns contents of ep1 as a double



function t2_Callback(hObject, eventdata, handles)
% hObject    handle to t2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of t2 as text
%        str2double(get(hObject,'String')) returns contents of t2 as a double



function g_Callback(hObject, eventdata, handles)
% hObject    handle to g (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of g as text
%        str2double(get(hObject,'String')) returns contents of g as a double



function ep2_Callback(hObject, eventdata, handles)
% hObject    handle to ep2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ep2 as text
%        str2double(get(hObject,'String')) returns contents of ep2 as a double



function Epslon_Callback(hObject, eventdata, handles)
% hObject    handle to Epslon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Epslon as text
%        str2double(get(hObject,'String')) returns contents of Epslon as a double



function min_Fx_Callback(hObject, eventdata, handles)
% hObject    handle to min_Fx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of min_Fx as text
%        str2double(get(hObject,'String')) returns contents of min_Fx as a double


function max_Fx_Callback(hObject, eventdata, handles)
% hObject    handle to max_Fx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of max_Fx as text
%        str2double(get(hObject,'String')) returns contents of max_Fx as a double



function min_thetax_Callback(hObject, eventdata, handles)
% hObject    handle to min_thetax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of min_thetax as text
%        str2double(get(hObject,'String')) returns contents of min_thetax as a double



function min_phi_Callback(hObject, eventdata, handles)
% hObject    handle to min_phi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of min_phi as text
%        str2double(get(hObject,'String')) returns contents of min_phi as a double



function max_thetax_Callback(hObject, eventdata, handles)
% hObject    handle to max_thetax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of max_thetax as text
%        str2double(get(hObject,'String')) returns contents of max_thetax as a double



function max_phi_Callback(hObject, eventdata, handles)
% hObject    handle to max_phi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of max_phi as text
%        str2double(get(hObject,'String')) returns contents of max_phi as a double
