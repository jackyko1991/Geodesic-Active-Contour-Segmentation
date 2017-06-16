function varargout = T2(varargin)
% T2 MATLAB code for T2.fig
%      T2, by itself, creates a new T2 or raises the existing
%      singleton*.
%
%      H = T2 returns the handle to a new T2 or the handle to
%      the existing singleton*.
%
%      T2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in T2.M with the given input arguments.
%
%      T2('Property','Value',...) creates a new T2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before T2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to T2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help T2

% Last Modified by GUIDE v2.5 21-Aug-2012 17:37:58

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @T2_OpeningFcn, ...
                   'gui_OutputFcn',  @T2_OutputFcn, ...
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

% --- Executes just before T2 is made visible.
function T2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to T2 (see VARARGIN)

% Choose default command line output for T2
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
% global T;


% --- Outputs from this function are returned to the command line.
function varargout = T2_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in segment.
function segment_Callback(hObject, eventdata, handles)
% hObject    handle to segment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes1);%用axes命令设定当前操作的坐标轴是axes_src
load matfile\Dimage;
load matfile\Dset1;
imagesc(Image{1});%用imread读入图片，并用imagesc在axes_src上显示
axis image
colormap(gray);hold on;
contour(u,[0 0],'r');
contour(u_0,[0 0],'r');
% UIWAIT makes T2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Executes on button press in sector.
function sector_Callback(hObject, eventdata, handles)
% hObject    handle to sector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load matfile\Dimage; 
T = zeros(size(Image{1}));
f = handles.figure1;
aH = handles.axes1;
% [x0 y0] = ginput(1)
x0 = 55;
y0 = 128;
r = 50;
i = 0;
for theta = 0:1/4*pi:2*pi
    x = x0 + r*cos(theta);
    y = y0 + r*sin(theta);
    i = i+1;
    h(i) = line([x0 x],[y0, y],'color','g');    
end
% x1 = x0+r*cos(4/3*pi);
% y1 = y0+r*sin(4/3*pi);
% h1 = line([x0,x1],[y0,y1],'color','r');
% x2 = x0+r*cos(1/2*pi);
% y2 = y0+r*sin(1/2*pi);
% h2 = line([x0,x2],[y0,y2],'color','g');
% p0 = [x0 y0];
% DragLine(f,aH,h1);
% DragLine(f,aH,h2);
save matfile\h.mat h
% getAssMatrix(p0,h1,h2,Image{1});

% --- Executes on button press in Analysis.
function Analysis_Callback(hObject, eventdata, handles)
% hObject    handle to Analysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load matfile\info
load matfile\Gset1
% load matfile\T.mat
load matfile\h
load matfile\Dimage
T = {};
j = 1;
t = length(h);
for i=2:t-1
    T{j} = getAssMatrix(h(i),h(i+1),Image{1});
    j = j+1;
end
T{j} = getAssMatrix(h(2),h(i+1),Image{1});
u1 = u;
u2 = u_0;
u1 = standard(u1);
u2 = standard(u2);
u_r = xor(u1,u2);
sector = {};
si = [];
Area = [];
for i = 1:j
    sector{i} = and(u_r,T{i});
    [si(i,:),Area(i,:)] = getImageIntensity(sector{i},Image);
end
% %%

[C,h0] = contour(sector{1},[0 0]);


L = length(C(1,:));
x1 = C(1,2:L);
y1 = C(2,2:L);
fill(x1,y1,'g');

rnames = cell(1,9);
for i=1:8
    rnames{i}=['TE' num2str(info{i}.EchoTime)];   
   % rnames{i}=['TE' num2str(i)]
end
rnames{9}='T2*';

for i=1:8
    te(i) = info{i}.EchoTime;
end
a = [];
for i= 1:j
    f = @(p,te)p(1)*exp(p(2)*te);
    x0 = [100,-1];
    p = lsqcurvefit(f,x0,te,si(i,:));
    a(i) = p(1);
    t2star(i) = -1/p(2);
end
si(:,9) = t2star;
rnames = cell(1,9);
for i=1:8
    rnames{i}=['TE' num2str(info{i}.EchoTime)];   
   % rnames{i}=['TE' num2str(i)]
end

rnames{9}='T2*';
t= uitable(handles.T2Table);
% set(t,'Position',[20 200 300 200])
set(t,'RowName',rnames,'Data',si');
save matfile\sector.mat sector;
save matfile\si.mat si;
save matfile\t2star.mat t2star;
save matfile\a.mat a;
save matfile\te.mat te;


% --- Executes on button press in Curve.
function Curve_Callback(hObject, eventdata, handles)
% hObject    handle to Curve (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load matfile\si
load matfile\t2star
load matfile\a
load matfile\te
c = 1;
SI = si(c,1:8);
t2 = t2star(c);
b = a(c);
axes(handles.axes2);
T2curve(t2,te,SI,b)
