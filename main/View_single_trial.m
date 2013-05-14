function varargout = View_single_trial(varargin)
% VIEW_SINGLE_TRIAL MATLAB code for View_single_trial.fig
%      VIEW_SINGLE_TRIAL, by itself, creates a new VIEW_SINGLE_TRIAL or raises the existing
%      singleton*.
%
%      H = VIEW_SINGLE_TRIAL returns the handle to a new VIEW_SINGLE_TRIAL or the handle to
%      the existing singleton*.
%
%      VIEW_SINGLE_TRIAL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VIEW_SINGLE_TRIAL.M with the given input arguments.
%
%      VIEW_SINGLE_TRIAL('Property','Value',...) creates a new VIEW_SINGLE_TRIAL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before View_single_trial_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to View_single_trial_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help View_single_trial

% Last Modified by GUIDE v2.5 23-Feb-2013 10:53:12

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @View_single_trial_OpeningFcn, ...
                   'gui_OutputFcn',  @View_single_trial_OutputFcn, ...
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


% --- Executes just before View_single_trial is made visible.
function View_single_trial_OpeningFcn(hObject, eventdata, handles, varargin)
global h1
% Choose default command line output for View_single_trial
handles.output = hObject;
scrsz = get(0,'ScreenSize');
h1=figure('Position',[scrsz(4)/2 scrsz(4) scrsz(3)/2.5 scrsz(4)/1.5],'name','View Trials','numbertitle','off','color','w');

% Update handles structure
guidata(hObject, handles);



% UIWAIT makes View_single_trial wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = View_single_trial_OutputFcn(hObject, eventdata, handles) 

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in nexttrial.
function nexttrial_Callback(hObject, eventdata, handles)
currenttrial=str2num(get(handles.currenttrial,'String'));
maxtrials=str2num(get(handles.numtrials,'String'));
if(currenttrial+1<=maxtrials)
    currenttrial =currenttrial +1;
    set(handles.currenttrial,'String',num2str(currenttrial));
%     currenttrial_Callback(handles.currenttrial,handles.rois_toplot);
    currenttrial_Callback(hObject,eventdata,handles);
else
    disp('Cannot increment trials anymore');
end

% --- Executes on button press in prevtrial.
function prevtrial_Callback(hObject, eventdata, handles)
currenttrial=str2num(get(handles.currenttrial,'String'));
if(currenttrial-1>0)
    currenttrial =currenttrial -1;
    set(handles.currenttrial,'String',num2str(currenttrial));
%     currenttrial_Callback(handles.currenttrial,handles.rois_toplot);
    currenttrial_Callback(hObject,eventdata,handles);
else
    disp('Cannot decrement trials anymore');
end

function currenttrial_Callback(hObject, eventdata, handles)
currenttrial = str2num(get(handles.currenttrial,'String'));
roiNo = str2num(get(handles.rois_toplot,'String'));
global obj 
global h1
figure(h1);
% allLickTimes = sesObj.behavArray.get_all_lick_times([]);
% fig_export_dir = 'Single_Trials_sessObj_181851_22';

%%
yl_dff = [-50 350];

numtrials = size(obj,2);
numrois= obj(1).nROIs;
trNo = currenttrial;

     
    roiName = num2str(roiNo);
    fig_export_dir = ['Rois_' roiName '_Trialdata'];
    folder=dir(fig_export_dir);
    if(isempty(folder))
        mkdir(fig_export_dir);
    end
    
    
    % Get Ca signal and whisker variable in trials
% % %         %% temporary fix for the 500ms delay in wSig data : advancing CaSig data by 500ms
% % %         delay = ceil(.5/obj(trNo).FrameTime);
% % %         dffarray = obj(trNo).dff(:,delay+1:end);
% % %         barOnOff = [.5, 2];
% % %          %% remove for data after Dec 28,2012 and uncomment the following statements
    [dffarray] = obj(trNo).dff;
       barOnOff = [1.08, 2.58];
    dff = mean(dffarray(roiNo,:),1);
     ts = obj(trNo).ts{1};
     ts = ts - ts(1);
    % Determine trial type
    soloTrNo =  obj(trNo).TrialName

    trType = obj(trNo).trialtype;
    wskNo = '2';
    ts_wsk = obj(trNo).ts_wsk{1};
      ts_wsk =ts_wsk-ts_wsk(1);  

      % Plot Ca signal trial
     contacts=obj(trNo).contacts{1};
      barOnOff =  barOnOff-(contacts/500)+ 0.5;
     licks=obj(trNo).licks;
     lickTime_trial = licks(licks >3);
     lickTime_trial = lickTime_trial - (contacts/500)+ 0.5;

    ha(1) = subaxis(5,1,1, 'sv', 0.05);
%     setpt = sesObj.wsArray.wl_trials{trNo}.Setpoint{wskNo};
    kappa = obj(trNo).kappa{1};   
    plot(ts_wsk, kappa,'r','linewidth',2);
%     ylabel('Setpoint','fontsize',15)
    ylabel('Kappa','fontsize',15);
    set(gca,'XMinorTick','on','XTick',0:.2:6);
    condir= num2str(obj(trNo).contactdir{1});
    if(isempty(condir))
        condir='nocontact';
        
    end
% %     title(sprintf('Trial %d, ROI %d, wsk %d, %s',trNo, roiNo, wskNo, trType), 'fontsize',20)
     title(sprintf('%d:Trial %s, dFFmean(%s) %s   %s ',trNo,soloTrNo,roiName,trType,condir), 'fontsize',20);
%     title([num2str(trNo) ' Trial' num2str(soloTrNo)  ' dFFmean(135) ' trType], 'fontsize',20)

    ha(2) = subaxis(5,1,2, 'MarginTop', 0.1);
    plot(ts, dff, 'k','linewidth',2)
    line([lickTime_trial'; lickTime_trial'], ...
        [zeros(1,length(lickTime_trial)); ones(1,length(lickTime_trial))*20],'color','m','linewidth',2)
    if(barOnOff(1)>0)
        line([barOnOff(1); barOnOff(1)], [0 0; 100 100], 'color','c','linewidth',2);
    end
    if (barOnOff(2)<3.5)
    line([barOnOff(2); barOnOff(2)], [0 0; 100 100], 'color','b','linewidth',2);
    end
    ylabel('dF/F','fontsize',15)
    ylim(yl_dff);    set(gca,'XMinorTick','on','XTick',0:.2:6);
    
    ha(3) = subaxis(5,1,3, 'sv', 0.05);
    dKappa = obj(trNo).deltaKappa{1};
    plot(ts_wsk,dKappa,'r');

    if(~isempty(contacts))
         numcontacts = size(contacts,2);   
        for i = 1:numcontacts
            if(iscell(contacts))
                ithcontact=contacts{i}-contacts{1}+250; %with respect to the first contact being fixed at .5s
            else
               ithcontact=contacts(i)-contacts(1)+250;
            end
            ithcontact=(ithcontact/500);
            line([ithcontact;ithcontact],[zeros(1,length(ithcontact))+.5;ones(1,length(ithcontact))*5],'color',[.6 .5 0],'linewidth',.5)
   
        end
    end
    ylabel('deltaKappa','fontsize',15);    set(gca,'XMinorTick','on','XTick',0:.2:6);

% Plot whisker velocity trial
    ha(4) = subaxis(5,1,4, 'sv', 0.05);

    vel =  obj(trNo).velocity{1};  
    plot(ts_wsk, vel, 'r','linewidth',2);
    ylabel('Velocity','fontsize',15);    set(gca,'XMinorTick','on','XTick',0:.2:6);

    % Plot whisker theta angle trial
    ha(5) = subaxis(5,1,5, 'sv', 0.05);
%     theta = sesObj.wsArray.wl_trials{trNo}.theta{wskNo};
    theta=obj(trNo).theta{1};
    plot(ts_wsk, theta,'r','linewidth',2);
    ylabel('Theta','fontsize',15)
    xlabel('Time (s)','fontsize',20)
    set(ha, 'box','off', 'xlim',[0 max(ts_wsk)]);    set(gca,'XMinorTick','on','XTick',0:.2:6);
    figure(gcf)
    fname=sprintf('ROIs %s_%d_Trial_%s',roiName,trNo,soloTrNo);
    export_fig(fullfile(fig_export_dir,fname ),gcf,'-png')
%     export_fig(fullfile(fig_export_dir, sprintf('dFFmean(1 3 5)_Trial %s',soloTrNo)),gcf,'-png');

% --- Executes during object creation, after setting all properties.
function currenttrial_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in selectObj.
function selectObj_Callback(hObject, eventdata, handles)
global obj

list = get(hObject,'String');
target=list(get(hObject,'Value'));
if(strcmp(target,'contact_CaTrials'))
    [fname,path]=uigetfile(pwd,'contact_CaTrials obj');
    cd (path);
    load([path fname],'-mat');
    obj=contact_CaTrials;
    numrois=obj.nROIs;
else
    [fname,path]=uigetfile(pwd,'sessobj');
    cd (path);
    load([path fname],'-mat');
    obj=sessObj;
    numrois=obj.nROIs;
end
set(handles.totalrois,'String',num2str(numrois));
set(handles.numtrials,'String',num2str(size(obj,2)));
guidata(hObject,handles);
%currenttrial_Callback(handles.currenttrial);
currenttrial_Callback(hObject,eventdata,handles);

% --- Executes during object creation, after setting all properties.
function selectObj_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function rois_toplot_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function rois_toplot_Callback(hObject, eventdata, handles)

handles.roiNo = str2num(get(hObject,'String'));


% --- Executes on key press with focus on rois_toplot and none of its controls.
function rois_toplot_KeyPressFcn(hObject, eventdata, handles)
handles.roiNo = str2num(get(hObject,'String'));
