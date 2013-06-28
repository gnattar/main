function varargout = Universal_GR_11(varargin)
% UNIVERSAL_GR_11 M-file for Universal_GR_11.fig
%      UNIVERSAL_GR_11, by itself, creates a new UNIVERSAL_GR_11 or raises the existing
%      singleton*.
%
%      H = UNIVERSAL_GR_11 returns the handle to a new UNIVERSAL_GR_11 or the handle to
%      the existing singleton*.
%
%      UNIVERSAL_GR_11('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in UNIVERSAL_GR_11.M with the given input arguments.
% 
%      UNIVERSAL_GR_11('Property','Value',...) creates a new UNIVERSAL_GR_11 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Universal_GR_11_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Universal_GR_11_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Universal_GR_11
 
% Last Modified by GUIDE v2.5 21-May-2013 11:37:41

% Begin initialization code - DO NOT EDIT

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Universal_GR_11_OpeningFcn, ...
                   'gui_OutputFcn',  @Universal_GR_11_OutputFcn, ...
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


% --- Executes just before Universal_GR_11 is made visible.
function Universal_GR_11_OpeningFcn(hObject, eventdata, handles, varargin)
global CaSignal % ROIinfo ICA_ROIs
% Choose default command line output for Universal_GR_11
handles.output = hObject;
set (handles.primProc,'Value',1);
set (handles.secProc,'Value',0);

usrpth = '~/Documents/MATLAB'; usrpth = usrpth(1:end-1);
if exist([usrpth filesep 'gr_Universal.info'],'file')
    load([usrpth filesep 'gr_UNiversal.info'], '-mat');
    set(handles.DataPathEdit, 'String',info.DataPath);
    set(handles.AnimalNameEdit, 'String', info.AnimalName);
    set(handles.ExpDate,'String',info.ExpDate);
    set(handles.SessionName, 'String',info.SessionName);
    
    if isfield(info, 'SoloDataPath')
        set(handles.SoloDataPath, 'String', info.SoloDataPath);
        set(handles.SoloDataFileName, 'String', info.SoloDataFileName);
        set(handles.SoloSessionID, 'String', info.SoloSessionName);
        set(handles.SoloStartTrialNo, 'String', info.SoloStartTrialNo);
        set(handles.SoloEndTrialNo, 'String', info.SoloEndTrialNo);
    end
   
    
else
    set(handles.DataPathEdit, 'String', '/Volumes/GR_Data_02/Data/');
    set(handles.SoloDataPath, 'String', '/Volumes/GR_Data_02/Data/');
    set(handles.baseDataPath, 'String', '/Volumes/GR_Data_02/Data/');
end
    set(handles.ephusDataPath, 'String', '/Volumes/GR_Data_02/Data/');

    % Open and Display section
set(handles.dispModeGreen, 'Value', 1);
set(handles.dispModeRed, 'Value', 0);
set(handles.dispModeImageInfoButton, 'Value', 0);
set(handles.dispModeWithROI, 'Value', 1);
% set(handles.LUTminEdit, 'Value', 0);
% set(handles.LUTmaxEdit, 'Value', 500);
% set(handles.LUTminSlider, 'Value', 0);
% set(handles.LUTmaxSlider, 'Value', 0.5);
set(handles.CurrentImageFilenameText, 'String', 'Current Image Filename');
    % ROI section
set(handles.nROIsText, 'String', '0');
set(handles.CurrentROINoEdit, 'String', '1');
set(handles.ROITypeMenu, 'Value', 1);
    % Analysis mode
set(handles.AnalysisModeDeltaFF, 'Value', 1);
set(handles.AnalysisModeBGsub, 'Value', 0);
set(handles.batchStartTrial, 'String', '1');
set(handles.batchEndTrial, 'String', '1');
set(handles.ROI_modify_button, 'Value', 0);
set(handles.CurrentFrameNoEdit,'String',1);
set(handles.nogopos,'String','16');
set(handles.gopos,'String','0 1.5 3 4.5 6 7.5');
set(handles.SoloStartTrialNo,'Value',1);
set(handles.SoloEndTrialNo,'Value',1);
set(handles.CaSignalrois,'String','0');
set(handles.numblocks,'Value',2);


CaSignal.CaTrials = struct([]);
CaSignal.ROIinfo = struct('ROImask',{}, 'ROIpos',{}, 'ROItype',{},'BGpos',[],...
        'BGmask', [], 'Method','');
% CaSignal.ICA_ROIs = struct('ROImask',{}, 'ROIpos',{}, 'ROItype',{},'Method','ICA');
CaSignal.ImageArray = [];
CaSignal.nFrames = 0;
% handles.userdata.CaTrials = [];
CaSignal.h_info_fig = NaN;
CaSignal.FrameNum = 1;
CaSignal.imSize = [];
CaSignal.h_img = NaN;
CaSignal.Scale = [0 500];
% ROIinfo = {};
% ICA_ROIs = struct;
CaSignal.ROIplot = NaN;
CaSignal.avgCorrCoef_trials = [];
CaSignal.CorrMapTrials = [];
CaSignal.CorrMapROINo = [];
CaSignal.AspectRatio_mode = 'Square';
CaSignal.ICA_figs = nan(1,2);
handles.SoloStartTrial=1;
handles.SoloEndTrial=1;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Universal_GR_11 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Universal_GR_11_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function CaTrial = init_CaTrial(filename, TrialNo, header)
% Initialize the struct data for the current trial
CaTrial.DataPath = pwd;
CaTrial.FileName = filename;
CaTrial.FileName_prefix = filename(1:end-7);

CaTrial.TrialNo = TrialNo;
CaTrial.DaqInfo = header;
if isfield(header, 'acq')
    CaTrial.nFrames = header.acq.numberOfFrames;
    CaTrial.FrameTime = header.acq.msPerLine*header.acq.linesPerFrame;
else
    CaTrial.nFrames = header.n_frame;
    if(~isempty(header.SI4.fastZPeriod)&&(header.SI4.fastZPeriod<1))
        CaTrial.FrameTime = header.SI4.fastZPeriod;
    elseif(~isempty(header.SI4.scanFramePeriod)&&(header.SI4.scanFramePeriod<1))
        CaTrial.FrameTime = header.SI4.scanFramePeriod;
    end
end
% % % if   CaTrial.FrameTime< 1 % some earlier version of ScanImage use sec as unit for msPerLine
% % %     CaTrial.FrameTime = CaTrial.FrameTime*1000;
% % % end
CaTrial.nROIs = 0;
CaTrial.BGmask = []; % Logical matrix for background ROI
CaTrial.AnimalName = '';
CaTrial.ExpDate = '';
CaTrial.SessionName = '';
CaTrial.dff = [];
CaTrial.f_raw = [];
% CaTrial.meanImage = [];
CaTrial.RegTargetFrNo = [];
CaTrial.ROIinfo = struct('ROImask',{}, 'ROIpos',{}, 'ROItype',{},'Method','');
CaTrial.SoloDataPath = '';
CaTrial.SoloFileName = '';
CaTrial.SoloSessionName = '';
CaTrial.SoloTrialNo = [];
CaTrial.SoloStartTrialNo = [];
CaTrial.SoloEndTrialNo = [];
CaTrial.behavTrial = [];
% CaTrial.ROIType = '';

% --- Executes on button press in open_image_file_button.
function open_image_file_button_Callback(hObject, eventdata, handles, filename)

global CaSignal % ROIinfo ICA_ROIs

datapath = get(handles.DataPathEdit,'String');
if exist(datapath, 'dir')
    cd(datapath);
else
clear     warning([datapath ' not exist!'])
    if exist('/Volumes/GR_Data_01/Data','dir')
        cd('/Volumes/GR_Data_01/Data');
    end
end;
if ~exist('filename', 'var')
    [filename, pathName] = uigetfile('*.tif', 'Load Image File');
    if isequal(filename, 0) || isequal(pathName,0)
        return
    end
    cd(pathName);
    FileName_prefix = filename(1:end-7);
    CaSignal.data_files = dir([FileName_prefix '*.tif']);
    CaSignal.data_file_names = {};
    for i = 1:length(CaSignal.data_files)
        CaSignal.data_file_names{i} = CaSignal.data_files(i).name;
    end;
end
datapath = pwd;
set(handles.DataPathEdit,'String',datapath);

FileName_prefix = filename(1:end-7);

TrialNo = find(strcmp(filename, CaSignal.data_file_names));
set(handles.CurrentTrialNo,'String', int2str(TrialNo));
disp(['Loading image file ' filename ' ...']);
set(handles.msgBox, 'String', ['Loading image file ' filename ' ...']);
[im, header] = imread_multi_GR(filename, 'g');
% t_elapsed = toc;
set(handles.msgBox, 'String', ['Loaded file ' filename]);
info = imfinfo(filename);
if isfield(info(1), 'ImageDescription')
    CaSignal.ImageDescription = info(1).ImageDescription; % used by Turboreg
else
    CaSignal.ImageDescription = '';
end
CaSignal.ImageArray = im;
CaSignal.imSize = size(im);
if ~isempty(CaSignal.CaTrials)
    if length(CaSignal.CaTrials)<TrialNo || isempty(CaSignal.CaTrials(TrialNo).FileName)
        CaSignal.CaTrials(TrialNo) = init_CaTrial(filename, TrialNo, header);
    end
    if ~strcmp(CaSignal.CaTrials(TrialNo).FileName_prefix, FileName_prefix)
        CaSignal.CaTrials_INIT = 1;
    else
        CaSignal.CaTrials_INIT = 0;
    end
else
    CaSignal.CaTrials_INIT = 1;
end

% Now we are in data file path. Since analysis results are saved in a separate
% folder, we need to find that folder in order to laod or save analysis
% results. If that folder does not exist, a new folder will be created.

CaSignal.results_path = strrep(datapath,[filesep 'data'],[filesep 'results']);
if ~exist(CaSignal.results_path,'dir')
    mkdir(CaSignal.results_path);
    disp('results dir not exists! A new folder created!');
    disp(CaSignal.results_path);
end

CaSignal.results_fn = [CaSignal.results_path filesep 'CaSignal_CaTrials_' FileName_prefix '.mat'];
CaSignal.results_roiinfo = [CaSignal.results_path filesep 'ROIinfo_', FileName_prefix '.mat'];   

if CaSignal.CaTrials_INIT == 1
   CaSignal.CaTrials = []; % ROIinfo = {};
    if exist(CaSignal.results_fn,'file')
        load(CaSignal.results_fn, '-mat');
        CaSignal.CaTrials = CaTrials;
    else
        A = init_CaTrial(filename, TrialNo, header);
        A(TrialNo) = A;
        if TrialNo ~= 1
            names = fieldnames(A);
            for i = 1:length(names)
                A(1).(names{i})=[];
            end
        end
        CaSignal.CaTrials = A;
    end
    
    if exist(CaSignal.results_roiinfo,'file')
        load(CaSignal.results_roiinfo, '-mat');
        if iscell(ROIinfo)
            f1 = fieldnames(ROIinfo{TrialNo}); f2 = fieldnames(CaSignal.ROIinfo);
            for i = 1:length(ROIinfo)
                for j = 1:length(f1),
                    CaSignal.ROIinfo(i).(f2{strcmpi(f2,f1{j})}) = ROIinfo{i}.(f1{j});
                end
            end
        else
            CaSignal.ROIinfo = ROIinfo;
        end
    end
else
    
    notansferROIinfo =get(handles.donottransferROIinfo,'Value');
    if ((notansferROIinfo)&& ~isempty(CaSignal.ROIinfo))
        CaSignal.ROIinfo=CaSignal.ROIinfo;
    else
        import_ROIinfo_from_trial_Callback(handles.import_ROIinfo_from_trial, eventdata, handles);
    end
end

if exist([CaSignal.results_path filesep FileName_prefix(1:end-7) '[dftShift].mat'],'file')
    load([CaSignal.results_path filesep FileName_prefix(1:end-7) '[dftShift].mat']);
    CaSignal.dftreg_shift = shift;
else
    CaSignal.dftreg_shift = [];
end

% Collect info to be displayed in a separate figure

% if get(handles.dispModeImageInfoButton,'Value') == 1
if isfield(header,'acq')
    CaSignal.info_disp = {sprintf('numFramesPerTrial: %d', header.acq.numberOfFrames), ...
    ['Zoom: ' num2str(header.acq.zoomFactor)],...
    ['numOfChannels: ' num2str(header.acq.numberOfChannelsAcquire)],...
    sprintf('ImageDimXY: %d,  %d', header.acq.pixelsPerLine, header.acq.linesPerFrame),...
    sprintf('Frame Rate: %d', header.acq.frameRate), ...
    ['msPerLine: ' num2str(header.acq.msPerLine)],...
    ['fillFraction: ' num2str(header.acq.fillFraction)],...
    ['motor_absX: ' num2str(header.motor.absXPosition)],...
    ['motor_absY: ' num2str(header.motor.absYPosition)],...
    ['motor_absZ: ' num2str(header.motor.absZPosition)],...
    ['num_zSlice: ' num2str(header.acq.numberOfZSlices)],...
    ['zStep: ' num2str(header.acq.zStepSize)] ...
    ['triggerTime: ' header.internal.triggerTimeString]...
    };
else
    CaSignal.info_disp = [];
end
%     dispModeImageInfoButton_Callback(hObject, eventdata, handles)
% end;
%% Initialize UI values
set(handles.TotTrialNum, 'String', int2str(length(CaSignal.data_file_names)));
set(handles.CurrentImageFilenameText, 'String',  filename);
if CaSignal.CaTrials_INIT == 1
    set(handles.DataPathEdit, 'String', pwd);
    set(handles.AnimalNameEdit, 'String', CaSignal.CaTrials(TrialNo).AnimalName);
    set(handles.ExpDate,'String',CaSignal.CaTrials(TrialNo).ExpDate);
    set(handles.SessionName, 'String',CaSignal.CaTrials(TrialNo).SessionName);
    if isfield(CaSignal.CaTrials(TrialNo), 'SoloDataFileName')
        set(handles.SoloDataPath, 'String', CaSignal.CaTrials(TrialNo).SoloDataPath);
        set(handles.SoloDataFileName, 'String', CaSignal.CaTrials(TrialNo).SoloDataFileName);
        set(handles.SoloSessionID, 'String', CaSignal.CaTrials(TrialNo).SoloSessionName);
        set(handles.SoloStartTrialNo, 'String', num2str(CaSignal.CaTrials(TrialNo).SoloStartTrialNo));
        set(handles.SoloEndTrialNo, 'String', num2str(CaSignal.CaTrials(TrialNo).SoloEndTrialNo));
    end
end

nFrames = size(im, 3);
set(handles.FrameSlider, 'SliderStep', [1/nFrames 1/nFrames]);
set(handles.FrameSlider, 'Value', 1/nFrames);
if ~isempty(CaSignal.ROIinfo)
    set(handles.nROIsText, 'String', int2str(length(CaSignal.ROIinfo(TrialNo).ROIpos)));
end
CaSignal.nFrames = nFrames;
set(handles.batchPrefixEdit, 'String', FileName_prefix);
%    handles = get_exp_info(hObject, eventdata, handles);
% CaSignal.CaTrials(TrialNo).meanImage = mean(im,3);

% update target info for TurboReg
% setTargetCurrentFrame_Callback(handles.setTargetCurrentFrame, eventdata, handles);
% setTargetMaxDelta_Callback(handles.setTargetMaxDelta, eventdata,handles);
% setTargetMean_Callback(handles.setTargetMaxDelta, eventdata, handles);

CaSignal.avgCorrCoef_trials = [];

% The trialNo to load ROIinfo from
notansferROIinfo =get(handles.donottransferROIinfo,'Value');
if ((notansferROIinfo)&& ~isempty(CaSignal.ROIinfo))
    TrialNo_load = TrialNo;
else
    TrialNo_load = str2double(get(handles.import_ROIinfo_from_trial,'String'));
end
if TrialNo_load > 0 && length(CaSignal.ROIinfo)>= TrialNo_load
    CaSignal.ROIinfo(TrialNo) = CaSignal.ROIinfo(TrialNo_load);
    nROIs = length(CaSignal.ROIinfo(TrialNo).ROIpos);
    CaSignal.CaTrials(TrialNo).nROIs = nROIs;
    set(handles.nROIsText, 'String', num2str(nROIs));
end


handles = update_image_axes(handles,im);
update_projection_images(handles);
% set(handles.figure1, 'WindowScrollWheelFcn',{@figScroll, handles.figure1, eventdata, handles});
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%% Start of Independent functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = get_exp_info(hObject, eventdata, handles)
global CaSignal % ROIinfo ICA_ROIs

TrialNo = str2double(get(handles.CurrentTrialNo,'String'));
filename = CaSignal.data_file_names{TrialNo};

if ~isempty(CaSignal.CaTrials(TrialNo).ExpDate)
    ExpDate = CaSignal.CaTrials(TrialNo).ExpDate;
    set(handles.ExpDate, 'String', ExpDate);
else
    CaSignal.CaTrials(TrialNo).ExpDate = get(handles.ExpDate, 'String');
end;


if ~isempty(CaSignal.CaTrials(TrialNo).AnimalName)
    AnimalName = CaSignal.CaTrials(TrialNo).AnimalName;
    set(handles.AnimalNameEdit, 'String', AnimalName);
else
    CaSignal.CaTrials(TrialNo).AnimalName = get(handles.AnimalNameEdit, 'String');
end


if ~isempty(CaSignal.CaTrials(TrialNo).SessionName)
    SessionName = CaSignal.CaTrials(TrialNo).SessionName;
    set(handles.SessionName, 'String', SessionName);
else
    CaSignal.CaTrials(TrialNo).SessionName = get(handles.SessionName, 'String');
end



function handles = update_image_axes(handles,varargin)
% update image display, called by most of call back functions
global CaSignal % ROIinfo ICA_ROIs
TrialNo = str2double(get(handles.CurrentTrialNo,'String')); 
LUTmin = str2double(get(handles.LUTminEdit,'String'));
LUTmax = str2double(get(handles.LUTmaxEdit,'String'));
sc = [LUTmin LUTmax];
cmap = 'gray';
fr = str2double(get(handles.CurrentFrameNoEdit,'String'));
if fr > CaSignal.nFrames && CaSignal.nFrames > 0
    fr = CaSignal.nFrames;
end
if ~isempty(varargin)
    CaSignal.ImageArray = varargin{1};
end
CaSignal.Scale = sc;
CaSignal.FrameNum = fr;

if get(handles.dispCorrMapTrials,'Value')==1
    im = CaSignal.CorrMapTrials;
%     cmap = 'jet';
    sc = [0 max(im(:))];
else
    im = CaSignal.ImageArray;
end;
im_size = size(im);
switch CaSignal.AspectRatio_mode
    case 'Square'
        s1 = im_size(2)/max(im_size(1:2));
        s2 = im_size(1)/max(im_size(1:2));
        asp_ratio = [s1 s2 1]; 
    case 'Image'
        asp_ratio = [1 1 1];
end
axes(handles.Image_disp_axes);
% hold on;
% if (isfield(CaSignal, 'h_img')&& ishandle(CaSignal.h_img))
%     delete(CaSignal.h_img);
% end;
% CaSignal.h_img = imagesc(im(:,:,fr), sc);
CaSignal.h_img = imshow(im(:,:,fr), sc);
set(handles.Image_disp_axes, 'DataAspectRatio', asp_ratio);  %'XTickLabel','','YTickLabel','');
time_str = sprintf('%.3f  sec',CaSignal.CaTrials(1).FrameTime*fr/1000);
set(handles.frame_time_disp, 'String', time_str);
% colormap(gray);
if get(handles.dispModeWithROI,'Value') == 1 && length(CaSignal.ROIinfo) >= TrialNo && ~isempty(CaSignal.ROIinfo(TrialNo).ROIpos)
    update_ROI_plot(handles);
end

% set(handles.figure1, 'WindowScrollWheelFcn',{@figScroll, hObject, eventdata, handles});

guidata(handles.figure1, handles);


function update_ROI_plot(handles)
global CaSignal % ROIinfo ICA_ROIs

CurrentROINo = str2double(get(handles.CurrentROINoEdit,'String'));
TrialNo = str2double(get(handles.CurrentTrialNo,'String')); 

if get(handles.dispModeWithROI,'Value') == 1
    axes(handles.Image_disp_axes);
    % delete existing ROI plots
    if any(ishandle(CaSignal.ROIplot))
        try
            delete(CaSignal.ROIplot(ishandle(CaSignal.ROIplot)));
        end
    end
    CaSignal.ROIplot = plot_ROIs(handles);
end
if isfield(CaSignal.ROIinfo, 'BGpos') && ~isempty(CaSignal.ROIinfo(TrialNo).BGpos)
    BGpos = CaSignal.ROIinfo(TrialNo).BGpos;
    CaSignal.BGplot = line(BGpos(:,1),BGpos(:,2), 'Color', 'b', 'LineWidth', 2);
end
% set(handles.figure1, 'WindowScrollWheelFcn',{@figScroll, handles.figure1, eventdata, handles});
guidata(handles.figure1,handles);

function h_roi_plots = plot_ROIs(handles)
%%
global CaSignal % ROIinfo ICA_ROIs
CurrentROINo = str2double(get(handles.CurrentROINoEdit,'String'));
TrialNo = str2double(get(handles.CurrentTrialNo,'String'));
h_roi_plots = [];
roi_pos = {};
%     if get(handles.ICA_ROI_anal, 'Value') == 1 && isfield(ICA_ROIs, 'ROIpos');
%         roi_pos = ICA_ROIs.ROIpos;
%     elseif length(ROIinfo) >= TrialNo && ~isempty(ROIinfo{TrialNo})
%
%         roi_pos = ROIinfo{TrialNo}.ROIpos;
%     end
if length(CaSignal.ROIinfo) >= TrialNo
    roi_pos = CaSignal.ROIinfo(TrialNo).ROIpos;
end
for i = 1:length(roi_pos) % num ROIs
    if i == CurrentROINo
        lw = 1.2;
        col = [.6 .5 .1];
    else
        lw = .5;
        col = [0.8 0 0];
    end
    if ~isempty(roi_pos{i})
        %             if length(CaSignal.ROIplot)>=i & ~isempty(CaSignal.ROIplot(i))...
        %                     & ishandle(CaSignal.ROIplot(i))
        %                 delete(CaSignal.ROIplot(i));
        %             end
        h_roi_plots(i) = line(roi_pos{i}(:,1),roi_pos{i}(:,2), 'Color', col, 'LineWidth', lw);
        text(median(roi_pos{i}(:,1)+5), median(roi_pos{i}(:,2)+5), num2str(i),'Color',[0 .7 0],'FontSize',12);
        set(h_roi_plots(i), 'LineWidth', lw);
    end
end

function handles = update_projection_images(handles)
global CaSignal % ROIinfo ICA_ROIs

if get(handles.dispMeanMode, 'Value')==1
    if ~isfield(CaSignal, 'h_mean_fig') || ~ishandle(CaSignal.h_mean_fig)
        CaSignal.h_mean_fig = figure('Name','Mean Image','Position',[960   500   480   480]);
    else
        figure(CaSignal.h_mean_fig)
    end
    if get(handles.dispCorrMapTrials,'Value')==1
        mean_im = CaSignal.CorrMapMean;
        sc = [0 max(mean_im(:))];
        colormap('jet');
    else
        im = CaSignal.ImageArray;
        sc = CaSignal.Scale;
        mean_im = mean(im,3);
        colormap(gray); 
    end
    imagesc(mean_im, sc);
    set(gca, 'Position',[0.05 0.05 0.9 0.9], 'Visible','off');
    update_projection_image_ROIs(handles);
end
if get(handles.dispMaxDelta,'Value')==1
    if ~isfield(CaSignal, 'h_maxDelta_fig') || ~ishandle(CaSignal.h_maxDelta_fig)
        CaSignal.h_maxDelta_fig = figure('Name','max Delta Image','Position',[960   40   480   480]);
    else
        figure(CaSignal.h_maxDelta_fig);
    end
    im = CaSignal.ImageArray;
    sc = CaSignal.Scale;
    mean_im = uint16(mean(im,3));
    im = im_mov_avg(im,5);
    max_im = max(im,[],3);
    CaSignal.MaxDelta = max_im - mean_im;
    imagesc(CaSignal.MaxDelta, sc);
    colormap(gray); 
    set(gca, 'Position',[0.05 0.05 0.9 0.9], 'Visible','off');
    
    update_projection_image_ROIs(handles);
end
if get(handles.dispMaxMode,'Value')==1
    if ~isfield(CaSignal, 'h_max_fig') || ~ishandle(CaSignal.h_max_fig)
        CaSignal.h_max_fig = figure('Name','Max Projection Image','Position',[960   180   480   480]);
    else
        figure(CaSignal.h_max_fig)
    end
    im = CaSignal.ImageArray;
    sc = CaSignal.Scale;
    im = im_mov_avg(im,5);
    max_im = max(im,[],3);
    imagesc(max_im, sc);
    colormap(gray); 
    set(gca, 'Position',[0.05 0.05 0.9 0.9], 'Visible','off');
    update_projection_image_ROIs(handles);
end
guidata(handles.figure1,handles);

% update ROI plotting in projecting image figure, called only by updata_projection image
function update_projection_image_ROIs(handles)
% global CaSignal ROIinfo ICA_ROIs
if get(handles.dispModeWithROI,'Value') == 1 
    plot_ROIs(handles);
end

function figScroll(src,evnt, hObject, eventdata, handles)
global CaSignal % ROIinfo ICA_ROIs
% callback function for mouse scroll
% 
im = CaSignal.ImageArray;
fr = str2double(get(handles.CurrentFrameNoEdit, 'String'));
sc = CaSignal.Scale;
nFrames = CaSignal.nFrames;
% axes(handles.Image_disp_axes);
if evnt.VerticalScrollCount > 0
    if fr < nFrames
        fr = fr + 1;
    end
    
elseif evnt.VerticalScrollCount < 0
    if fr > 1
        fr = fr - 1;
    end  
end

set(handles.FrameSlider,'Value', fr/nFrames);

CaSignal.FrameNum = fr;
set(handles.CurrentFrameNoEdit, 'String', num2str(fr));

CaSignal.h_img = imagesc(im(:,:,fr), sc);
colormap(gray);

handles = update_image_axes(handles);

% Update handles structure
% guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% End of Independent functions %%%%%%%%%%%%%%%%%%%%%%%%%


function dispModeWithROI_Callback(hObject, eventdata, handles)
value = get(handles.dispModeWithROI,'Value');
handles = update_image_axes(handles);
handles = update_projection_images(handles);

function DataPathEdit_Callback(hObject, eventdata, handles)
handles.datapath = get(hObject, 'String');
guidata(hObject, handles);

function DataPathEdit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ROI_add_Callback(hObject, eventdata, handles)

% global CaSignal % ROIinfo ICA_ROIs
nROIs = str2num(get(handles.nROIsText, 'String'));
nROIs = nROIs + 1;
set(handles.nROIsText, 'String', num2str(nROIs));

CurrentROINo = get(handles.CurrentROINoEdit,'String');
% if strcmp(CurrentROINo, '0')
%     set(handles.CurrentROINoEdit,'String', '1');
% end;
% Use this instead: automatically go to the last ROI added.
set(handles.CurrentROINoEdit,'String', num2str(nROIs));
guidata(hObject, handles);


function ROI_del_Callback(hObject, eventdata, handles)
global CaSignal % ROIinfo ICA_ROIs
TrialNo = str2double(get(handles.CurrentTrialNo,'String'));
CurrentROINo = str2double(get(handles.CurrentROINoEdit,'String'));
if CurrentROINo > 0
    if length(CaSignal.ROIplot) >= CurrentROINo && ishandle(CaSignal.ROIplot(CurrentROINo))
        try
            delete(CaSignal.ROIplot(CurrentROINo))
            CaSignal.ROIplot(CurrentROINo)=[];
        end
    end
    %     if get(handles.ICA_ROI_anal,'Value') ==1 &&  length(ICA_ROIs.ROIpos) >= CurrentROINo
    %         ICA_ROIs.ROIpos(CurrentROINo) = [];
    %         ICA_ROIs.ROIMask(CurrentROINo) = [];
    %         try
    %         ICA_ROIs.ROIType(CurrentROINo) = [];
    %         end
    %     elseif length(ROIinfo{TrialNo}.ROIpos(CurrentROINo)) >= CurrentROINo
    %         ROIinfo{TrialNo}.ROIpos(CurrentROINo) = [];
    %         ROIinfo{TrialNo}.ROIMask(CurrentROINo) = [];
    %         ROIinfo{TrialNo}.ROIType(CurrentROINo) = [];
    %         CaSignal.CaTrials(TrialNo).nROIs = CaSignal.CaTrials(TrialNo).nROIs - 1;
    %         CaSignal.CaTrials(TrialNo).ROIinfo = ROIinfo{TrialNo};
    %     end
    CaSignal.ROIinfo(TrialNo).ROIpos(CurrentROINo) = [];
    CaSignal.ROIinfo(TrialNo).ROImask(CurrentROINo) = [];
    CaSignal.ROIinfo(TrialNo).ROItype(CurrentROINo) = [];
    CaSignal.CaTrials(TrialNo).nROIs = CaSignal.CaTrials(TrialNo).nROIs - 1;
    CaSignal.CaTrials(TrialNo).ROIinfo = CaSignal.ROIinfo(TrialNo);
    set(handles.nROIsText, 'String', num2str(CaSignal.CaTrials(TrialNo).nROIs));
    set(handles.CurrentROINoEdit,'String', int2str(CurrentROINo - 1));
    % TotROI = get(handles.nROIsText, 'String');
    % if strcmp(TotROI, '0');
    %     set(handles.CurrentROINoEdit,'String', '0');
    % end
    update_ROI_plot(handles);
    update_ROI_numbers(handles);
end
guidata(hObject, handles);



function update_ROI_numbers(handles)
global CaSignal % ROIinfo ICA_ROIs
TrialNo = str2double(get(handles.CurrentTrialNo,'string'));
CurrentROINo = str2double(get(handles.CurrentROINoEdit,'String'));
% if get(handles.ICA_ROI_anal,'value') ==1
%     nd = cellfun(@(x) isempty(x), ICA_ROIs.ROIpos);
%     try
%         ICA_ROIs.ROIpos(nd) = [];
%         ICA_ROIs.ROIMask(nd) = [];
%         ICA_ROIs.ROIType(nd) = [];
%     end
%     nROIs = length(ICA_ROIs.ROIpos);
% else
%     for i = 1:length(ROIinfo{TrialNo}.ROIpos)
%         if isempty(ROIinfo{TrialNo}.ROIpos{i})
%             ROIinfo{TrialNo}.ROIpos(i) = [];
%             ROIinfo{TrialNo}.ROIMask(i) = [];
%             ROIinfo{TrialNo}.ROIType(i) = [];
%         end
%     end
%     nROIs = length(ROIinfo{TrialNo}.ROIpos);
% end
for i = 1:length(CaSignal.ROIinfo(TrialNo).ROIpos)
    if isempty(CaSignal.ROIinfo(TrialNo).ROIpos{i})
        CaSignal.ROIinfo(TrialNo).ROIpos(i) = [];
        CaSignal.ROIinfo(TrialNo).ROImask(i) = [];
        CaSignal.ROIinfo(TrialNo).ROItype(i) = [];
    end
end
    nROIs = length(CaSignal.ROIinfo(TrialNo).ROIpos);
set(handles.nROIsText, 'String', num2str(nROIs));
if CurrentROINo > nROIs
    CurrentROINo = nROIs;
elseif CurrentROINo < 1
    CurrentROINo = 1;
end
set(handles.CurrentROINoEdit, 'String', num2str(CurrentROINo));



function ROI_pre_Callback(hObject, eventdata, handles)
global CaSignal % ROIinfo ICA_ROIs
% update_ROI_numbers(handles);
TrialNo = str2double(get(handles.CurrentTrialNo,'String'));
CurrentROINo = str2double(get(handles.CurrentROINoEdit,'String'));
CurrentROINo = CurrentROINo - 1;
if CurrentROINo <= 0
    CurrentROINo = 1;
end;
set(handles.CurrentROINoEdit,'String',int2str(CurrentROINo));

str_menu = get(handles.ROITypeMenu,'String');
ROIType_str = CaSignal.ROIinfo(TrialNo).ROItype{CurrentROINo};
if ~isempty(ROIType_str)
    ROIType_num = find(strcmp(ROIType_str, str_menu));
    set(handles.ROITypeMenu,'Value', ROIType_num);
else
    ROIType_str = str_menu{get(handles.ROITypeMenu,'Value')};
    CaSignal.ROIinfo(TrialNo).ROItype{CurrentROINo} = ROIType_str;
end
% axes(handles.Image_disp_axes);
update_ROI_plot(handles);
handles = update_projection_images(handles);
guidata(hObject, handles);



function ROI_next_Callback(hObject, eventdata, handles)
global CaSignal % ROIinfo ICA_ROIs
% update_ROI_numbers(handles);
TrialNo = str2double(get(handles.CurrentTrialNo,'String'));
CurrentROINo = str2double(get(handles.CurrentROINoEdit,'String'));
CurrentROINo = CurrentROINo + 1;
if CurrentROINo > str2double(get(handles.nROIsText,'String')) 
    CurrentROINo = str2double(get(handles.nROIsText,'String')) ;
end;
set(handles.CurrentROINoEdit,'String',int2str(CurrentROINo));

str_menu = get(handles.ROITypeMenu,'String');
if length(CaSignal.ROIinfo(TrialNo).ROItype)>= CurrentROINo
    % ~isempty(ROIinfo{TrialNo}.ROIType{CurrentROINo})
    
    ROIType_str = CaSignal.ROIinfo(TrialNo).ROItype{CurrentROINo};
    if ~isempty(ROIType_str)
        ROIType_num = find(strcmp(ROIType_str, str_menu));
        set(handles.ROITypeMenu,'Value', ROIType_num);
    else
        ROIType_str = str_menu{get(handles.ROITypeMenu,'Value')};
        CaSignal.ROIinfo(TrialNo).ROItype{CurrentROINo} = ROIType_str;
    end
else
    CaSignal.ROIinfo(TrialNo).ROItype{CurrentROINo} = str_menu{get(handles.ROITypeMenu,'Value')};
end
update_ROI_plot(handles);
handles = update_projection_images(handles);
guidata(hObject, handles);


function ROI_del_all_Callback(hObject, eventdata, handles)


function CurrentROINoEdit_Callback(hObject, eventdata, handles)
update_ROI_plot(handles);
guidata(hObject, handles);


function CurrentROINoEdit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function ROI_set_poly_Callback(hObject, eventdata, handles)
global CaSignal % ROIinfo ICA_ROIs
TrialNo = str2double(get(handles.CurrentTrialNo,'String'));
CurrentROINo = str2double(get(handles.CurrentROINoEdit, 'String'));
str_menu = get(handles.ROITypeMenu,'String');
ROIType = str_menu{get(handles.ROITypeMenu,'Value')};
% ROI_updated_flag = 0; % to determine if update the trial No of ROI updating.
%% Draw an ROI after mouse press
waitforbuttonpress;
% define the way of drawing, freehand or ploygon
if get(handles.ROI_draw_freehand, 'Value') == 1
    draw = @imfreehand;
elseif get(handles.ROI_draw_poly, 'Value') == 1;
    draw = @impoly;
end
h_roi = feval(draw);
finish_drawing = 0;
while finish_drawing == 0
    choice = questdlg('confirm ROI drawing?','confirm ROI', 'Yes', 'Re-draw', 'Cancel','Yes');
    switch choice
        case'Yes',
            pos = h_roi.getPosition;
            line(pos(:,1), pos(:,2),'color','g')
            BW = createMask(h_roi);
            delete(h_roi);
            finish_drawing = 1;
%             ROI_updated_flag = 1;
        case'Re-draw'
            delete(h_roi);
            h_roi = feval(draw); finish_drawing = 0;
        case'Cancel',
            delete(h_roi); finish_drawing = 1;
%             ROI_updated_flag = 0;
            return
    end
end

CaSignal.ROIinfo(TrialNo).ROIpos{CurrentROINo} = pos;
CaSignal.ROIinfo(TrialNo).ROImask{CurrentROINo} = BW;
CaSignal.ROIinfo(TrialNo).ROItype{CurrentROINo} = ROIType;
CaSignal.CaTrials(TrialNo).nROIs = length(CaSignal.ROIinfo(TrialNo).ROIpos);

set(handles.import_ROIinfo_from_trial, 'String', num2str(TrialNo));

if get(handles.ICA_ROI_anal,'Value') == 1
    CaSignal.ROIinfo(TrialNo).Method = 'ICA';
    CaSignal.rois_by_IC{CaSignal.currentIC} = [CaSignal.rois_by_IC{CaSignal.currentIC}  CurrentROINo];
    for jj = 1:length(CaSignal.ICA_figs)
        if ishandle(CaSignal.ICA_figs(jj))
            figure(CaSignal.ICA_figs(jj)),
            plot_ROIs(handles);
        end
    end
else
    
end
set(handles.nROIsText,'String', num2str(length(CaSignal.ROIinfo(TrialNo).ROIpos)));
axes(handles.Image_disp_axes);
if length(CaSignal.ROIplot) >= CurrentROINo
    if ishandle(CaSignal.ROIplot(CurrentROINo))
        delete(CaSignal.ROIplot(CurrentROINo));
    else
        CaSignal.ROIplot(CurrentROINo) = [];
    end
end
guidata(hObject, handles);
%CaSignal.roi_line(CurrentROINo) = line(pos(:,1),pos(:,2), 'Color', 'r', 'LineWidth', 2);
update_ROI_plot(handles);
update_projection_images(handles);
ROITypeMenu_Callback(hObject, eventdata, handles);
guidata(hObject, handles);



function ROI_modify_button_Callback(hObject, eventdata, handles)
global CaSignal % ROIinfo ICA_ROIs
TrialNo = str2double(get(handles.CurrentTrialNo,'String'));
CurrentROINo = str2double(get(handles.CurrentROINoEdit, 'String'));
pos = CaSignal.ROIinfo(TrialNo).ROIpos{CurrentROINo};
h_axes = handles.Image_disp_axes;

if get(hObject, 'Value')==1
    CaSignal.current_poly_obj = impoly(h_axes, pos);
elseif get(hObject, 'Value')== 0 
    if isa(CaSignal.current_poly_obj, 'imroi')
        pos = getPosition(CaSignal.current_poly_obj);
        BW = createMask(CaSignal.current_poly_obj);
        CaSignal.ROIinfo(TrialNo).ROIpos{CurrentROINo} = pos;
        CaSignal.ROIinfo(TrialNo).ROImask{CurrentROINo} = BW;
        axes(h_axes);
        delete(CaSignal.current_poly_obj); % delete polygon object
        if ishandle(CaSignal.ROIplot(CurrentROINo))
            delete(CaSignal.ROIplot(CurrentROINo));
        end
        CaSignal.ROIplot(CurrentROINo) = [];
        % CaSignal.roi_line(CurrentROINo) = line(pos(:,1),pos(:,2), 'Color', 'r', 'LineWidth', 2);
         update_ROI_plot(handles);
         handles = update_projection_images(handles);
    end;
end;
guidata(hObject, handles);

% --- Executes on button press in import_ROIinfo_from_file.
function import_ROIinfo_from_file_Callback(hObject, eventdata, handles)

choice = questdlg('Import ROIinfo from different file/session?', 'Import ROIs', 'Yes','No','No');
switch choice
    case 'Yes'
        [fn, pth] = uigetfile('*.mat');
        r = load([pth filesep fn]);
        ROIinfo = r.ROIinfo(1);
    case 'No'
        return
end
import_ROIinfo(ROIinfo, handles);


function import_ROIinfo_from_trial_Callback(hObject, eventdata, handles)
% get ROIinfo from the specified trial, and call import_ROIinfo function
global CaSignal

notansferROIinfo =get(handles.donottransferROIinfo,'Value');
if ((notansferROIinfo)&& ~isempty(CaSignal.ROIinfo))
    disp('Not transfering ROIinfo because donottranfer flag is on')
    return
end

% The trialNo to load ROIinfo from
TrialNo_load = str2double(get(handles.import_ROIinfo_from_trial,'String'));
if ~isempty(CaSignal.ROIinfo)
    ROIinfo = CaSignal.ROIinfo(TrialNo_load);
    import_ROIinfo(ROIinfo, handles);% getROIinfoButton_Callback(hObject, eventdata, handles)
else
    warning('No ROIs specified!');
end

function import_ROIinfo(ROIinfo, handles)
% update the ROIs of the current trial with the input "ROIinfo".
global CaSignal % ROIinfo ICA_ROIs
TrialNo = str2double(get(handles.CurrentTrialNo,'String'));

FileName_prefix = CaSignal.CaTrials(TrialNo).FileName_prefix;

CaSignal.ROIinfo(TrialNo) = ROIinfo;
% elseif exist(['ROIinfo_' FileName_prefix '.mat'],'file')
%     load([FileName_prefix 'ROIinfo.mat'], '-mat');
%     if length(CaSignal.ROIinfo)>= TrialNo_load
%         CaSignal.ROIinfo(TrialNo) = CaSignal.ROIinfo{TrialNo_load};
%     end
nROIs = length(CaSignal.ROIinfo(TrialNo).ROIpos);
CaSignal.CaTrials(TrialNo).nROIs = nROIs;
set(handles.nROIsText, 'String', num2str(nROIs));
update_ROI_plot(handles);
handles = update_projection_images(handles);



function import_ROIinfo_from_trial_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ROITypeMenu_Callback(hObject, eventdata, handles)
global CaSignal % ROIinfo ICA_ROIs
TrialNo = str2double(get(handles.CurrentTrialNo,'String'));
CurrentROINo = str2double(get(handles.CurrentROINoEdit, 'String'));
Menu = get(handles.ROITypeMenu,'String');
% CaSignal.CaTrials(TrialNo).ROIType{CurrentROINo} = Menu{get(handles.ROITypeMenu,'Value')};
CaSignal.ROIinfo(TrialNo).ROItype{CurrentROINo} = Menu{get(handles.ROITypeMenu,'Value')};
guidata(hObject, handles);

function CalculatePlotButton_Callback(hObject, eventdata, handles, im, plot_flag)
global CaSignal % ROIinfo ICA_ROIs
TrialNo = str2double(get(handles.CurrentTrialNo,'String'));
% ROIMask = CaSignal.CaTrials(TrialNo).ROIMask;
% if get(handles.ICA_ROI_anal, 'Value') == 1
%     nROI_effective = length(ICA_ROIs.ROIMask);
%     ROImask = ICA_ROIs.ROIMask;
% else
%     nROI_effective = length(CaSignal.ROIinfo(TrialNo).ROIpos);
%     ROImask = CaSignal.ROIinfo(TrialNo).ROIMask;
% end
nROI_effective = length(CaSignal.ROIinfo(TrialNo).ROIpos);
ROImask = CaSignal.ROIinfo(TrialNo).ROImask;
if nargin < 4
    im = CaSignal.ImageArray;
end    
if nargin < 5 %~exist('plot_flag','var')
    plot_flag = 1;
end
if ~isempty(CaSignal.ROIinfo(TrialNo).BGmask)
    BGmask = repmat(CaSignal.ROIinfo(TrialNo).BGmask,[1 1 size(im,3)]) ;
    BG_img = BGmask.*double(im);
    BG_img(BG_img==0) = NaN;
    BG = reshape(nanmean(nanmean(BG_img)),1,[]); % 1-by-nFrames array
else
    BG = 0;
end;

F = zeros(nROI_effective, size(im,3));
dff = zeros(size(F));

for i = 1: nROI_effective
    mask = repmat(ROImask{i}, [1 1 size(im,3)]); % reproduce masks for every frame
    % Using indexing and reshape function to increase speed
    nPix = sum(sum(ROImask{i}));
    % Using reshape to partition into different trials.
    roi_img = reshape(im(mask), nPix, []);
    % Raw intensity averaged from pixels of the ROI in each trial.
    if nPix == 0
        F(i,:) = 0;
    else
        F(i,:) = mean(roi_img, 1);
    end
%%%%%%%%%%%%% Obsolete slower method to compute ROI pixel intensity %%%%%%%   
%     roi_img = mask .* double(im);                                       % 
%                                                                         %      
%     roi_img(roi_img<=0) = NaN;                                          %  
%    % F(:,i) = nanmean(nanmean(roi_img));                                %
%     F(i,:) = nanmean(nanmean(roi_img));                                 %  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if get(handles.AnalysisModeBGsub,'Value') == 1
        F(i,:) = F(i,:) - BG;
    end;
    if get(handles.AnalysisModeDeltaFF,'Value') == 1
        [N,X] = hist(F(i,:));
        F_mode = X(find(N==max(N)));
        baseline = mean(F_mode);
        dff(i,:) = (F(i,:)- baseline)./baseline*100;
    else
%         CaTrace(i,:) = F(i,:);
    end
end;
CaSignal.CaTrials(TrialNo).dff = dff;
CaSignal.CaTrials(TrialNo).f_raw = F;
ts = (1:CaSignal.CaTrials(TrialNo).nFrames).*CaSignal.CaTrials(TrialNo).FrameTime;
if plot_flag == 1
    if get(handles.check_plotAllROIs, 'Value') == 1
        roiNos = [];
    else
        roiNos = str2num(get(handles.roiNo_to_plot, 'String'));
    end
    CaSignal.h_CaTrace_fig = plot_CaTraces_ROIs(dff, ts, roiNos);
% % %     set(handles.Image_disp_axes,'Visible','Off');
% % %     set(handles.Image_disp_axes1,'Visible','On');
% % %     set(handles.Image_disp_axes2,'Visible','On');
%     plot_CaTraces_ROIs(dff, ts, roiNos,handles.Image_disp_axes1,handles.Image_disp_axes2);
end
guidata(handles.figure1, handles);

function doBatchButton_Callback(hObject, eventdata, handles)
global CaSignal % ROIinfo ICA_ROIs
batchPrefix = get(handles.batchPrefixEdit, 'String');
Start_trial = str2double(get(handles.batchStartTrial, 'String'));
End_trial = str2double(get(handles.batchEndTrial,'String'));
% CaSignal.CaTrials = [];
h = waitbar(0, 'Start Batch Analysis ...');
for TrialNo = Start_trial:End_trial
    fname = CaSignal.data_file_names{TrialNo};
    if ~exist(fname,'file')
        [fname, pathname] = uigetfile('*.tif', 'Select Image Data file');
        cd(pathname);
    end;
    msg_str1 = sprintf('Batch analyzing %d of total %d trials with %d ROIs...', ...
        TrialNo, End_trial-Start_trial+1, CaSignal.CaTrials(1).nROIs);  
%     disp(['Batch analyzing ' num2str(TrialNo) ' of total ' num2str(End_trial-Start_trial+1) ' trials...']);
    disp(msg_str1);
    waitbar((TrialNo-Start_trial+1)/(End_trial-Start_trial+1), h, msg_str1);
    set(handles.msgBox, 'String', msg_str1);
% %     [im, header] = imread_multi(fname, 'g');
    [im, header] = imread_multi_GR(fname, 'g'); 
    if (length(CaSignal.CaTrials)<TrialNo || isempty(CaSignal.CaTrials(TrialNo).FileName))
        trial_init = init_CaTrial(fname,TrialNo,header);
        CaSignal.CaTrials(TrialNo) = trial_init;
    end
    set(handles.CurrentTrialNo,'String', int2str(TrialNo));
    
    notansferROIinfo =get(handles.donottransferROIinfo,'Value');

    if (notansferROIinfo)
    %###########################################################################
    % Make sure the ROIinfo of the first trial of the batch is up to date
    if TrialNo > Start_trial && ~isempty(CaSignal.ROIinfo(TrialNo-1).ROIpos)
       if( (size(CaSignal.ROIinfo,2)<TrialNo) || ( size(CaSignal.ROIinfo(TrialNo).ROIpos,2)~= size(CaSignal.ROIinfo(TrialNo-1).ROIpos,2)))
            CaSignal.ROIinfo(TrialNo) = CaSignal.ROIinfo(TrialNo-1);
            CaSignal.CaTrials(TrialNo).nROIs = CaSignal.CaTrials(TrialNo-1).nROIs;
            CaSignal.CaTrials.FrameTime(TrialNo)=CaSignal.CaTrials.FrameTime(TrialNo-1);
       else
            CaSignal.ROIinfo(TrialNo)= CaSignal.ROIinfo(TrialNo);
            CaSignal.CaTrials(TrialNo).nROIs=CaSignal.CaTrials(TrialNo).nROIs;
            
       end
    end
    
    %##########################################################################
    else       
        if TrialNo > Start_trial && ~isempty(CaSignal.ROIinfo(TrialNo-1).ROIpos)
                CaSignal.ROIinfo(TrialNo) = CaSignal.ROIinfo(TrialNo-1);
                CaSignal.CaTrials(TrialNo).nROIs = CaSignal.CaTrials(TrialNo-1).nROIs;
        end
        
    end
%     update_image_axes(handles,im);
    CalculatePlotButton_Callback(handles.figure1, eventdata, handles, im, 0);
%     handles = update_projection_images(handles);
    handles = get_exp_info(hObject, eventdata, handles);
%     CaSignal.CaTrials(TrialNo).meanImage = mean(im,3);
%     close(CaSignal.h_CaTrace_fig);
    set(handles.CurrentTrialNo, 'String', int2str(TrialNo));
    set(handles.CurrentImageFilenameText,'String',fname);
%     set(handles.nROIsText,'String',int2str(length(ROIinfo{TrialNo}.ROIpos)));
    guidata(hObject, handles);
end

SaveResultsButton_Callback(hObject, eventdata, handles);
close(h);
disp(['Batch analysis completed for ' CaSignal.CaTrials(1).FileName_prefix]);
set(handles.msgBox, 'String', ['Batch analysis completed for ' CaSignal.CaTrials(1).FileName_prefix]);

function SaveResultsButton_Callback(hObject, eventdata, handles)
%% Save Results
global CaSignal % ROIinfo ICA_ROIs
TrialNo = str2double(get(handles.CurrentTrialNo,'String'));
FileName_prefix = CaSignal.CaTrials(TrialNo).FileName_prefix;
datapath = get(handles.DataPathEdit,'String');
cd (datapath);
% CaSignal.CaTrials = CaSignal.CaTrials;
% ROIinfo = ROIinfo;
ICA_results.ROIinfo = [];
ICA_results.ICA_components = [];
ICA_results.rois_by_IC = [];
if get(handles.ICA_ROI_anal, 'Value') == 1
    ICA_results.ROIinfo = CaSignal.ROIinfo;
    ICA_results.ICA_components = CaSignal.ICA_components;
    ICA_results.rois_by_IC = CaSignal.rois_by_IC;
end
ROIinfo = CaSignal.ROIinfo;
for i = 1:length(CaSignal.CaTrials)
    if length(CaSignal.ROIinfo) >= i
        CaSignal.CaTrials(i).ROIinfo = CaSignal.ROIinfo(i);
    end
end
CaTrials = CaSignal.CaTrials;
save(CaSignal.results_fn, 'CaTrials','ICA_results','-v7.3');
save(CaSignal.results_roiinfo, 'ROIinfo','-v7.3');
save(fullfile(CaSignal.results_path, ['ICA_ROIs_', FileName_prefix '.mat']), 'ICA_ROIs');
msg_str = sprintf('CaTrials Saved, with %d trials, %d ROIs', length(CaSignal.CaTrials), CaSignal.CaTrials(TrialNo).nROIs);
disp(msg_str);
set(handles.msgBox, 'String', msg_str);
current_dir = pwd;
separators = find(current_dir == filesep);
session_dir = current_dir(1:separators(length(separators)-3));
cd (session_dir);
sessObj_found = dir('sessObj.mat');
if isempty(sessObj_found)
    sessObj = {};
    sessObj.CaTrials = CaSignal.CaTrials;
    save('sessObj','sessObj','-v7.3');
else
    load('sessObj.mat');
    sessObj.CaTrials = CaSignal.CaTrials;
    save('sessObj','sessObj','-v7.3');
end
cd (current_dir);
save_gui_info(handles);


function batchStartTrial_Callback(hObject, eventdata, handles)

function batchStartTrial_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function batchEndTrial_Callback(hObject, eventdata, handles)

function batchEndTrial_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function dispModeGreen_Callback(hObject, eventdata, handles)
global CaSignal % ROIinfo ICA_ROIs
TrialNo = str2double(get(handles.CurrentTrialNo,'String'));
if CaSignal.CaTrials(TrialNo).DaqInfo.acq.numberOfChannelsAcquire == 1
    set(hObject,'Value',1);
end;

function dispModeImageInfoButton_Callback(hObject, eventdata, handles)
global CaSignal % ROIinfo ICA_ROIs
if get(hObject, 'Value') == 1
    CaSignal.h_info_fig = figure; set(gca, 'Visible', 'off');
    f_pos = get(CaSignal.h_info_fig, 'Position'); f_pos(3) = f_pos(3)/2;
    set(CaSignal.h_info_fig, 'Position', f_pos);
    info_disp = CaSignal.info_disp;
    for i = 1: length(info_disp),
        x = 0.01;
        y=1-i/length(info_disp);
        text(x,y,info_disp{i},'Interpreter','none');
    end
    guidata(hObject, handles);
else
    close(CaSignal.h_info_fig);
end


function nROIsText_CreateFcn(hObject, eventdata, handles)

function figure1_DeleteFcn(hObject, eventdata, handles)
global CaSignal % ROIinfo ICA_ROIs
save_gui_info(handles);
clear CaSignal % ROIinfo ICA_ROIs
% close all;

function CurrentTrialNo_Callback(hObject, eventdata, handles)
global CaSignal % ROIinfo ICA_ROIs
TrialNo = str2double(get(handles.CurrentTrialNo,'String'));
if TrialNo>0
    filename = CaSignal.data_file_names{TrialNo};
    if exist(filename,'file')
        open_image_file_button_Callback(hObject, eventdata, handles,filename);
    end
end

function PrevTrialButton_Callback(hObject, eventdata, handles)
global CaSignal % ROIinfo ICA_ROIs
TrialNo = str2double(get(handles.CurrentTrialNo,'String'));
if TrialNo>1
    filename = CaSignal.data_file_names{TrialNo-1};
    if exist(filename,'file')
        open_image_file_button_Callback(hObject, eventdata, handles,filename);
    end
end

function NextTrialButton_Callback(hObject, eventdata, handles)
global CaSignal % ROIinfo ICA_ROIs
TrialNo = str2double(get(handles.CurrentTrialNo,'String'));
if  TrialNo+1 <= length(CaSignal.data_file_names) % exist(filename,'file')
    filename = CaSignal.data_file_names{TrialNo+1};
    open_image_file_button_Callback(hObject, eventdata, handles,filename);
end

function AnalysisModeBGsub_Callback(hObject, eventdata, handles)

function batchPrefixEdit_Callback(hObject, eventdata, handles)

function batchPrefixEdit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function AnimalNameEdit_Callback(hObject, eventdata, handles)
% % % global CaSignal % ROIinfo ICA_ROIs
% % % TrialNo = str2double(get(handles.CurrentTrialNo,'String'));
% % % CaSignal.CaTrials(TrialNo).AnimalName = get(hObject, 'String');
% % % guidata(hObject, handles);
% % % save_gui_info(handles);

function AnimalNameEdit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ExpDate_Callback(hObject, eventdata, handles)
% % % global CaSignal % ROIinfo ICA_ROIs
% % % TrialNo = str2double(get(handles.CurrentTrialNo,'String'));
% % % CaSignal.CaTrials(TrialNo).ExpDate = get(hObject, 'String');
% % % guidata(hObject, handles);
% % % save_gui_info(handles);

function ExpDate_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function SessionName_Callback(hObject, eventdata, handles)
% % % global CaSignal % ROIinfo ICA_ROIs
% % % TrialNo = str2double(get(handles.CurrentTrialNo,'String'));
% % % CaSignal.CaTrials(TrialNo).SessionName = get(hObject, 'String');
% % % guidata(hObject, handles);
% % % save_gui_info(handles);

function SessionName_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function FrameSlider_Callback(hObject, eventdata, handles)
global CaSignal % ROIinfo ICA_ROIs
slider_value = get(hObject,'Value');
nFrames = CaSignal.nFrames;
new_frameNum = ceil(nFrames*slider_value);
if new_frameNum == 0, new_frameNum = 1; end;
set(handles.CurrentFrameNoEdit, 'String', num2str(new_frameNum));
handles = update_image_axes(handles);
guidata(hObject, handles);

function FrameSlider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function CurrentFrameNoEdit_Callback(hObject, eventdata, handles)

handles = update_image_axes(handles);

function CurrentFrameNoEdit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function LUTminEdit_Callback(hObject, eventdata, handles)

update_image_axes(handles);
update_projection_images(handles);

function LUTminEdit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function LUTmaxEdit_Callback(hObject, eventdata, handles)

update_image_axes(handles);
update_projection_images(handles);


function LUTmaxEdit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function LUTminSlider_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
value_min = get(hObject,'Value');
value_max = get(handles.LUTmaxSlider,'Value');
if value_min >= value_max
    value_min = value_max - 0.01;
    set(hObject, 'Value', value_min);
end;
set(handles.LUTminEdit, 'String', num2str(value_min*1000));
update_image_axes(handles);
update_projection_images(handles);
% guidata(hObject, handles);


function LUTminSlider_CreateFcn(hObject, eventdata, handles)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function LUTmaxSlider_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
value_max = get(hObject,'Value');
value_min = get(handles.LUTminSlider, 'Value');
if value_max <= value_min
    value_max = value_min + 0.01;
    set(hObject, 'Value', value_max);
end;
set(handles.LUTmaxEdit, 'String', num2str(value_max*1000));
update_image_axes(handles);
update_projection_images(handles);
% guidata(hObject, handles);


function dispMeanMode_Callback(hObject, eventdata, handles)
global CaSignal % ROIinfo ICA_ROIs
if get(hObject, 'Value')==1
    handles = update_projection_images(handles);
else
    try
        if ishandle(CaSignal.h_mean_fig)
            delete(CaSignal.h_mean_fig);
        end;
    catch ME
    end
end
guidata(hObject, handles);


function dispMaxDelta_Callback(hObject, eventdata, handles)
global CaSignal % ROIinfo ICA_ROIs
if get(hObject, 'Value')==1
    handles = update_projection_images(handles);
else
    try
        if ishandle(CaSignal.h_maxDelta_fig)
            delete(CaSignal.h_maxDelta_fig);
        end;
    catch ME
    end
end
guidata(hObject, handles);

function dispMaxMode_Callback(hObject, eventdata, handles)
if get(hObject, 'Value')==1
    handles = update_projection_images(handles);
else
    try
        if ishandle(CaSignal.h_max_fig)
            delete(CaSignal.h_max_fig);
        end;
    catch ME
    end
end
guidata(hObject, handles);

function ROITypeMenu_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function dispModeGreen_CreateFcn(hObject, eventdata, handles)

function LUTmaxSlider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in SaveFrameButton.
function SaveFrameButton_Callback(hObject, eventdata, handles)
global CaSignal % ROIinfo ICA_ROIs
im = CaSignal.ImageArray;
fr = str2double(get(handles.CurrentFrameNoEdit,'String'));
dataFileName = get(handles.CurrentImageFilenameText, 'String');

[fname, pathName] = uiputfile([dataFileName(1:end-4) '_' int2str(fr) '.tif'], 'Save the current frame as');
if ~isequal(fname, 0)&& ~isequal(pathName, 0)
    imwrite(im(:,:,fr), [pathName fname], 'tif','WriteMode','overwrite','Compression','none');
end




% --- Executes on button press in BG_poly_set.
function BG_poly_set_Callback(hObject, eventdata, handles)
global CaSignal % ROIinfo ICA_ROIs
TrialNo = str2double(get(handles.CurrentTrialNo,'String'));
%    if isempty(CaSignal.CaTrials(TrialNo).BGmask)
waitforbuttonpress;
[BW,xi,yi] = roipoly;
CaSignal.ROIinfo(TrialNo).BGmask = BW;
CaSignal.ROIinfo(TrialNo).BGpos = [xi yi];
% axes(CaSignal.image_disp_gui.Image_disp_axes);
% if isfield(CaSignal, 'BGplot')&& ishandle(CaSignal.BGplot)
%     delete(CaSignal.BGplot);
% end
% CaSignal.BGplot = line(xi, yi, 'Color','b', 'LineWidth',2);
update_ROI_plot(handles);       
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function AnalysisModeBGsub_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AnalysisModeBGsub (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on selection change in RegMethodMenu.
function RegMethodMenu_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function RegMethodMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RegMethodMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in MotionEstmOptions.
function MotionEstmOptions_Callback(hObject, eventdata, handles)
% Hints: contents = get(hObject,'String') returns MotionEstmOptions contents as cell array
%        contents{get(hObject,'Value')} returns selected item from MotionEstmOptions
global CaSignal % ROIinfo ICA_ROIs
TrialNo = str2double(get(handles.CurrentTrialNo,'String'));
switch get(hObject,'Value')
    case 2 % plot cross correlation coef for the current trial
        img = CaSignal.ImageArray;
        xcoef = xcoef_img(img);
        figure('Name', ['xCorr. Coefficient for Trial ' num2str(TrialNo)], 'Position', [1200 300 480 300]);
        plot(xcoef); xlabel('Frame #'); ylabel('Corr. Coeff');
        disp(sprintf(['mean xCorr. Coefficient for trial ' num2str(TrialNo) ': %g'],mean(xcoef)));
    case 3 % Compute cross correlation across all trials
        n_trials = length(CaSignal.data_file_names);
        if isempty(CaSignal.avgCorrCoef_trials)
            xcoef_trials = zeros(n_trials,1);
            h_wait = waitbar(0, 'Calculating cross correlation coefficients for trial 0 ...');
            for i = 1:n_trials
                waitbar(i/n_trials, h_wait, ['Calculating cross correlation coefficients for trial ' num2str(i)]);
                img = imread_multi(CaSignal.data_file_names{i},'g');
                xcoef = xcoef_img(img);
                xcoef_trials(i) = mean(xcoef);
            end
            close(h_wait);
            CaSignal.avgCorrCoef_trials = xcoef_trials;
        else
            xcoef_trials = CaSignal.avgCorrCoef_trials;
        end
        figure('Name', 'xCorr. Coef across all trials', 'Position', [1200 300 480 300]);
        plot(xcoef_trials); xlabel('Trial #'); ylabel('mean Corr. Coeff');
    case 4
        
    case 5
        if ~isempty(CaSignal.dftreg_shift)
            for i = 1:str2num(get(handles.TotTrialNum, 'String'))
                avg_shifts(i) = max(mean(abs(CaSignal.dftreg_shift(:,:,i)),2));
            end
            figure;
            plot(avg_shifts,'LineWidth',2); 
            title('Motion estimation of all trials','FontSize',18);
            xlabel('Trial #', 'FontSize', 15); ylabel('Mean shift of all frames', 'FontSize', 15);
        end
end
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function MotionEstmOptions_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function save_gui_info(handles)
% global CaSignal ROIinfo ICA_ROIs
info.DataPath = pwd;
info.AnimalName = get(handles.AnimalNameEdit,'String');
info.ExpDate = get(handles.ExpDate,'String');
info.SessionName = get(handles.SessionName, 'String');
info.SoloDataPath = get(handles.SoloDataPath, 'String');
info.SoloDataFileName = get(handles.SoloDataFileName, 'String');
info.SoloSessionName = get(handles.SoloSessionID, 'String');
info.SoloStartTrialNo = get(handles.SoloStartTrialNo, 'String');
info.SoloEndTrialNo = get(handles.SoloEndTrialNo, 'String');

usrpth = '~/Documents/MATLAB'; %usrpth = usrpth(1:end-1);
save([usrpth filesep 'nx_CaSingal.info'], 'info');



function SoloStartTrialNo_Callback(hObject, eventdata, handles)
handles.SoloStartTrial=str2num(get(hObject,'String'));
guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
function SoloStartTrialNo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to solostarttrialno (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function SoloEndTrialNo_Callback(hObject, eventdata, handles)
%
handles.SoloEndTrial=str2num( get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function SoloEndTrialNo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to soloendtrialno (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function SoloDataPath_Callback(hObject, eventdata, handles)
%

% --- Executes during object creation, after setting all properties.
function SoloDataPath_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function SoloDataFileName_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function SoloDataFileName_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in addBehavTrials.
function addBehavTrials_Callback(hObject, eventdata, handles)
global CaSignal %  ROIinfo ICA_ROIs

Solopath = get(handles.SoloDataPath,'String');
mouseName = get(handles.AnimalNameEdit, 'String');
sessionName = get(handles.SoloDataFileName, 'String');
trialStartEnd(1) = str2num(get(handles.SoloStartTrialNo, 'String'));
trialStartEnd(2) = str2num(get(handles.SoloEndTrialNo, 'String'));
trailsToBeExcluded = str2num(get(handles.behavTrialNoToBeExcluded, 'String'));

[Solo_data, SoloFileName] = Solo.load_data_gr(mouseName, sessionName,trialStartEnd,Solopath);
set(handles.SoloDataFileName, 'String', SoloFileName);
behavTrialNums = trialStartEnd(1):trialStartEnd(2);
behavTrialNums(trailsToBeExcluded) = [];

if length(behavTrialNums) ~= str2num(get(handles.TotTrialNum, 'String'))
    error('Number of behavior trials NOT equal to Number of Ca Image Trials!')
end

for i = 1:length(behavTrialNums)
    behavTrials(i) = Solo.BehavTrial_gr(Solo_data,behavTrialNums(i),1);
    CaSignal.CaTrials(i).behavTrial = behavTrials(i);
end
disp([num2str(i) ' Behavior Trials added to CaSignal.CaTrials']);
set(handles.msgBox, 'String', [num2str(i) ' Behavior Trials added to CaSignal.CaTrials']);
guidata(hObject, handles)


function SoloSessionID_Callback(hObject, eventdata, handles)

% % % % Solopath = get(handles.SoloDataPath,'String');
% % % % mouseName = get(handles.AnimalNameEdit, 'String');
% % % % sessionID = get(handles.SoloSessionID, 'String');
% % % % sessionName = ['data_@pole_detect_gr_0obj_' mouseName '_' sessionID];
% % % % trialStartEnd(1) = str2num(get(handles.SoloStartTrialNo, 'String'));
% % % % trialStartEnd(2) = str2num(get(handles.SoloEndTrialNo, 'String'));
% % % % 
% % % % [Solo_data, SoloFileName] = Solo.load_data_gr(mouseName, sessionName,trialStartEnd,Solopath);
% % % % set(handles.SoloDataFileName, 'String', SoloFileName);


% --- Executes during object creation, after setting all properties.
function SoloSessionID_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function dispModeImageInfoButton_CreateFcn(hObject, eventdata, handles)


% --- Executes on button press in ROI_move_left.
function ROI_move_left_Callback(hObject, eventdata, handles)
global CaSignal % ROIinfo ICA_ROIs
TrialNo = str2double(get(handles.CurrentTrialNo,'String'));
imsize = size(CaSignal.ImageArray);
aspect_ratio = imsize(2)/imsize(1);
fine=get(handles.movefine,'Value');
if(fine)
    move_unit = 1* max(1/aspect_ratio,1);
else
    move_unit = 5* max(1/aspect_ratio,1);
end
if get(handles.ROI_move_all_check, 'Value') == 1
    roi_num_to_move = 1: length(CaSignal.ROIinfo(TrialNo).ROIpos);
else
    roi_num_to_move = str2num(get(handles.CurrentROINoEdit,'String'));
end
for i = roi_num_to_move
    CaSignal.ROIinfo(TrialNo).ROIpos{i}(:,1) = CaSignal.ROIinfo(TrialNo).ROIpos{i}(:,1)-move_unit;
    x = CaSignal.ROIinfo(TrialNo).ROIpos{i}(:,1);
    y = CaSignal.ROIinfo(TrialNo).ROIpos{i}(:,2);
    CaSignal.ROIinfo(TrialNo).ROImask{i} = poly2mask(x,y,imsize(1),imsize(2));
end;
update_ROI_plot(handles);
handles = update_projection_images(handles);
guidata(hObject, handles);


% --- Executes on button press in ROI_move_right.
function ROI_move_right_Callback(hObject, eventdata, handles)
global CaSignal % ROIinfo ICA_ROIs

TrialNo = str2double(get(handles.CurrentTrialNo,'String'));
imsize = size(CaSignal.ImageArray);
aspect_ratio = imsize(2)/imsize(1);
fine=get(handles.movefine,'Value');
if(fine)
    move_unit = 1* max(1/aspect_ratio,1);
else
    move_unit = 5* max(1/aspect_ratio,1);
end
if get(handles.ROI_move_all_check, 'Value') == 1
    roi_num_to_move = 1: length(CaSignal.ROIinfo(TrialNo).ROIpos);
else
    roi_num_to_move = str2num(get(handles.CurrentROINoEdit,'String'));
end
for i = roi_num_to_move
    CaSignal.ROIinfo(TrialNo).ROIpos{i}(:,1) = CaSignal.ROIinfo(TrialNo).ROIpos{i}(:,1)+move_unit;
    x = CaSignal.ROIinfo(TrialNo).ROIpos{i}(:,1);
    y = CaSignal.ROIinfo(TrialNo).ROIpos{i}(:,2);
    CaSignal.ROIinfo(TrialNo).ROImask{i} = poly2mask(x,y,imsize(1),imsize(2));
end;
update_ROI_plot(handles);
handles = update_projection_images(handles);
guidata(hObject, handles)

% --- Executes on button press in ROI_move_up.
function ROI_move_up_Callback(hObject, eventdata, handles)
global CaSignal % ROIinfo ICA_ROIs

TrialNo = str2double(get(handles.CurrentTrialNo,'String'));
imsize = size(CaSignal.ImageArray);
aspect_ratio = imsize(2)/imsize(1);
fine=get(handles.movefine,'Value');
if(fine)
    move_unit = 1* max(1/aspect_ratio,1);
else
    move_unit = 5* max(1/aspect_ratio,1);
end
if get(handles.ROI_move_all_check, 'Value') == 1
    roi_num_to_move = 1: length(CaSignal.ROIinfo(TrialNo).ROIpos);
else
    roi_num_to_move = str2num(get(handles.CurrentROINoEdit,'String'));
end
for i = roi_num_to_move
    CaSignal.ROIinfo(TrialNo).ROIpos{i}(:,2) = CaSignal.ROIinfo(TrialNo).ROIpos{i}(:,2)-move_unit;
    x = CaSignal.ROIinfo(TrialNo).ROIpos{i}(:,1);
    y = CaSignal.ROIinfo(TrialNo).ROIpos{i}(:,2);
    CaSignal.ROIinfo(TrialNo).ROImask{i} = poly2mask(x,y,imsize(1),imsize(2));
end;
update_ROI_plot(handles);
handles = update_projection_images(handles);
guidata(hObject, handles)


% --- Executes on button press in ROI_move_down.
function ROI_move_down_Callback(hObject, eventdata, handles)

global CaSignal % ROIinfo ICA_ROIs
TrialNo = str2double(get(handles.CurrentTrialNo,'String'));
imsize = size(CaSignal.ImageArray);
aspect_ratio = imsize(2)/imsize(1);
fine=get(handles.movefine,'Value');
if(fine)
    move_unit = 1* max(1/aspect_ratio,1);
else
    move_unit = 5* max(1/aspect_ratio,1);
end

if get(handles.ROI_move_all_check, 'Value') == 1
    roi_num_to_move = 1: length(CaSignal.ROIinfo(TrialNo).ROIpos);
else
    roi_num_to_move = str2num(get(handles.CurrentROINoEdit,'String'));
end
for i = roi_num_to_move
    CaSignal.ROIinfo(TrialNo).ROIpos{i}(:,2) = CaSignal.ROIinfo(TrialNo).ROIpos{i}(:,2)+move_unit;
    x = CaSignal.ROIinfo(TrialNo).ROIpos{i}(:,1);
    y = CaSignal.ROIinfo(TrialNo).ROIpos{i}(:,2);
    CaSignal.ROIinfo(TrialNo).ROImask{i} = poly2mask(x,y,imsize(1),imsize(2));
end

update_ROI_plot(handles);
handles = update_projection_images(handles);
guidata(hObject, handles)



function behavTrialNoToBeExcluded_Callback(hObject, eventdata, handles)
% hObject    handle to behavTrialNoToBeExcluded (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of behavTrialNoToBeExcluded as text
%        str2double(get(hObject,'String')) returns contents of behavTrialNoToBeExcluded as a double


% --- Executes during object creation, after setting all properties.
function behavTrialNoToBeExcluded_CreateFcn(hObject, eventdata, handles)
% hObject    handle to behavTrialNoToBeExcluded (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in dispCorrMapTrials.
function dispCorrMapTrials_Callback(hObject, eventdata, handles)
global CaSignal % ROIinfo ICA_ROIs
currROI = str2num(get(handles.CurrentROINoEdit,'String'));
if get(handles.dispCorrMapTrials,'Value')==1
    if isempty(CaSignal.CorrMapTrials) || ~isequal(CaSignal.CorrMapROINo, currROI)
        [fn, pathname] = uigetfile([CaSignal.results_path filesep '*.tif'],'Select image file of correlation map');
        CaSignal.CorrMapTrials = imread_multi([pathname filesep fn]);
        CaSignal.CorrMapROINo = currROI;
        % get the mean correlation map from trials with variance > 70 percentile
        dFF = get_dFF_roi(CaSignal.CaTrials, currROI);
        v = var(dFF,0,2);
        ind = find(v > prctile(v,70));
        CaSignal.CorrMapMean = mean(CaSignal.CorrMapTrials(:,:,ind),3);
    end
%     im = CaSignal.CorrMapTrials;
    CaSignal.nFrames = size(CaSignal.CorrMapTrials,3);
    set(handles.CurrentFrameNoEdit,'String',num2str(1));
    handles = update_image_axes(handles);
else
    CaSignal.nFrames = size(CaSignal.ImageArray,3);
    handles = update_image_axes(handles);
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function dispCorrMapTrials_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dispCorrMapTrials (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

function dFF_array = get_dFF_roi(CaSignal, roiNo)

nTrials = numel(CaSignal.CaTrials);
dFF_array = [];
for i = 1:nTrials
    if ~isempty(CaSignal.CaTrials(i).dff)
        if size(CaSignal.CaTrials(i).dff,1) < roiNo
            return;
        else
            dFF_array = [dFF_array; CaSignal.CaTrials(i).dff(roiNo,:)];
        end
    end
end
    
% --- Executes on button press in CorrMap_button.
function CorrMap_button_Callback(hObject, eventdata, handles)
global CaSignal % ROIinfo ICA_ROIs
confirm = questdlg('Start computing spatial temporal correlation for all ROIs? It would take a long time!');
if ~strcmp(confirm, 'Yes')
    return
else
    savepath = CaSignal.results_path;
    roiNums = 1:CaSignal.CaTrials(1).nROIs; 
    % fVar = var(dffROI,0,2);
    % ind = find(fVar > prctile(fVar,70)); % trNoROI; %  find(peakROI>80);
    nTr = length(CaSignal.CaTrials);
    width = CaSignal.CaTrials(1).DaqInfo.width;
    height = CaSignal.CaTrials(1).DaqInfo.height;
    
    for n = 1:numel(roiNums)
        crr_map = zeros(height,width,nTr);
        h = waitbar(0, sprintf('Analyzing %d of total %d ROIs',n,numel(roiNums)),... ) %['Analyzing ' num2str(n) ' of total ' num2str(numel(roiNums)) ' ROIs ...'],...
            'CreateCancelBtn', 'setappdata(gcbf,''cancelling'',1)');
        setappdata(h,'canceling',0);
        for k = 1:nTr
            if getappdata(h,'canceling')
                break
            end
            %     tr = ind(k);
            %     disp(['Analyzing ' imobj(tr).FileName]);
            fname = CaSignal.data_file_names{k};
            im = imread_multi_GR(fname,'g'); %%GRchange
            f = CaSignal.CaTrials(k).f_raw(roiNums(n),:);
            for i=1:size(im,1)
                for j = 1:size(im,2)
                    p = double(im(i,j,:));
                    c = corrcoef(f,p);
                    crr_map(i,j,k) = c(2,1);
                end
            end
            
            waitbar(k/nTr,h, sprintf('Analyzing %d of total %d ROIs, in Trial %d...', n,numel(roiNums),k)); %  ['Analyzing ROI # ' num2str(roiNums(n)) ', Trial ' num2str(k) ' ...']);
            imwrite(crr_map(:,:,k),sprintf('%s%cCorrMapROI#%d.tif',savepath,filesep,roiNums(n)),'Compression','none','WriteMode','append');
        end
        delete(h);
        save(sprintf('%s%cCorr_Map_ROI_#%d.mat', savepath,filesep, roiNums(n)),'crr_map');
        %     axes(ax1); plot(dffROI(ind(k),:)); axes(ax2); imagesc(crr_im(:,:,k),[-0.1 1]);
    end
end

    


% --- Executes during object creation, after setting all properties.
function Image_disp_axes_CreateFcn(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function msgBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to msgBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on selection change in dispAspectRatio.
function dispAspectRatio_Callback(hObject, eventdata, handles)
global CaSignal % ROIinfo ICA_ROIs
Str = get(hObject, 'String');
CaSignal.AspectRatio_mode = Str{get(hObject,'Value')};
guidata(hObject, handles);
handles = update_image_axes(handles);



% --- Executes during object creation, after setting all properties.
function dispAspectRatio_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dispAspectRatio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function frame_time_disp_CreateFcn(hObject, eventdata, handles)


% --- Executes on button press in export2avi_button.
function export2avi_button_Callback(hObject, eventdata, handles)

global CaSignal
TrialNo = str2double(get(handles.CurrentTrialNo,'String'));
fname = CaSignal.CaTrials(TrialNo).FileName;
[movieFileName, pathname] = uiputfile([fname(1:end-4) '.avi'], 'Export current trial to an avi movie');
movObj = VideoWriter([pathname filesep movieFileName]);
movObj.FrameRate = 15;

open(movObj);

for i = 1:CaSignal.CaTrials(TrialNo).nFrames
    set(handles.CurrentFrameNoEdit,'String',num2str(i));
    handles = update_image_axes(handles);
    F = getframe(handles.Image_disp_axes);
    writeVideo(movObj, F);
end
close(movObj);
% movie2avi(Mov,[pathname filesep movieFileName],'compression','none');


% --- Executes when selected object is changed in ROI_def.
function ROI_def_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in ROI_def 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function ICA_ROI_anal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ICA_ROI_anal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes on button press in ICA_ROI_anal.
function ICA_ROI_anal_Callback(hObject, eventdata, handles)

global CaSignal % ROIinfo ICA_ROIs    
if get(hObject, 'Value') == 1
        
    %%****************
    
        % Display mean ICA map
    if isfield(CaSignal, 'ICA_components') && ~isempty(CaSignal.ICA_components)
        CaSignal.rois_by_IC = cell(1, size(CaSignal.ICA_components,1)) ;
        IC_to_remove = inputdlg('IC to remove', 'Remove bad ICs');
        if ~isempty(IC_to_remove)
            IC_to_remove = str2num(IC_to_remove{1});
            CaSignal.ICA_components(IC_to_remove,:) = NaN;
        end
        CaSignal.ICA_figs(3) = disp_mean_IC_map(CaSignal.ICA_components);
    end
    
    if isfield(CaSignal, 'ica_data') && isfield(CaSignal.ica_data,'Data')
        usr_confirm = questdlg('Display Max Projection of all data is slow and memory intensive. Continue?');
        if strcmpi(usr_confirm, 'Yes')
%             Data = LoadData(pwd,CaSignal.CaTrials(1).FileName_prefix,1:50);
            [CaSignal.ICA_data_norm_max, CaSignal.ICA_figs(4)] = disp_maxDelta_rawData(CaSignal.ica_data.Data);
            [CaSignal.ICA_data_norm_max, CaSignal.ICA_figs(4)] = disp_maxDelta_rawData(CaSignal.ica_data.Data);
            set(gcf,'Name',sprintf('MaxDelta of raw Data (%d~%d)',CaSignal.ica_data.FileNums(1),CaSignal.ica_data.FileNums(end)));
        end
 
    end
    
    %**************
   
    % load_saved_data_SVD
   if isempty(CaSignal.ImageArray)
        error('No Imaging Data loaded. Do this before running ICA!');
   end
  
   if isfield(CaSignal, 'ica_data') && isfield(CaSignal.ica_data,'Data')%%~isempty(CaSignal.ica_data.Data)
          
         CaSignal. ica_data.FileBaseName = get(handles.batchPrefixEdit, 'String');
         CaSignal.ica_data.DataDir = pwd;
         CaSignal.ica_data.FileNums = [1:50];% using filenums [1:50] for now , need to put a text box
         usr_confirm = questdlg('ica_data.Data exists. Reload data from first 50 data .tif files.Continue? ');
          if strcmpi(usr_confirm, 'Yes')
             % Load the data files
             CaSignal.ica_data.Data = ICA_LoadData(CaSignal.ica_data.DataDir,CaSignal.ica_data.FileBaseName, CaSignal.ica_data.FileNums);
          end 
   
   else      

%          h = msgbox('No ica_data.Data, Loading data from first 50 data .tif files.');
       
         CaSignal.ica_data.FileBaseName = get(handles.batchPrefixEdit, 'String');
         CaSignal.ica_data.DataDir = pwd;
         CaSignal.ica_data.FileNums = [1:50];% using filenums [1:50] for now , need to put a text box
         % Load the data files
         CaSignal.ica_data.Data = ICA_LoadData( CaSignal.ica_data.DataDir,  CaSignal.ica_data.FileBaseName,  CaSignal.ica_data.FileNums); 
      
   end

    %if no svd data : compute svd & save it in same dir
    if ~exist(sprintf('ICA_data_%s.mat',  CaSignal.ica_data.FileBaseName), 'file')
        [U,S,V] = svd(CaSignal.ica_data.Data,'econ');
        save(sprintf('ICA_data_%s.mat',  CaSignal.ica_data.FileBaseName), 'U','S','V','-v7.3');
    end
    
% % %     %now load the svd data to run_ICA callback
% % %     [fn, pth] = uigetfile('*.mat', 'Load DATA SVD');
% % %     if fn == 0
% % %         return
% % %     end   
    usv = load(sprintf('ICA_data_%s.mat', CaSignal.ica_data.FileBaseName));

    CaSignal.ica_data.U = usv.U;
    CaSignal.ica_data.S = usv.S;
    CaSignal.ica_data.V = usv.V;

    CaSignal.currentIC = 1;
%     handles.ica_data = ica_data;
%     handles.ICA_datafile = fullfile(pth,fn);
%     ICA_ROIs = struct;
    guidata(hObject, handles);
    runICA_button_Callback(handles.runICA_button, eventdata, handles);
    

    
end


% --- Executes on button press in prevIC_button.
function prevIC_button_Callback(hObject, eventdata, handles)
global CaSignal % ROIinfo ICA_ROIs
if get(handles.ICA_ROI_anal,'Value') == 1 && CaSignal.currentIC > 1
    CaSignal.currentIC = CaSignal.currentIC - 1;
    set(handles.current_ICnum_text, 'String', num2str(CaSignal.currentIC));
%     guidata(hObject, handles);
    disp_ICA(handles);
end


% --- Executes on button press in nextIC_button.
function nextIC_button_Callback(hObject, eventdata, handles)
global CaSignal  % ROIinfo ICA_ROIs
if get(handles.ICA_ROI_anal,'Value') == 1 && CaSignal.currentIC < size(CaSignal.ICA_components,1)
    CaSignal.currentIC = CaSignal.currentIC + 1;
    set(handles.current_ICnum_text, 'String', num2str(CaSignal.currentIC));
    guidata(hObject, handles);
    disp_ICA(handles);
end

% --- Executes on button press in runICA_button.
function runICA_button_Callback(hObject, eventdata, handles)
global CaSignal % ROIinfo ICA_ROIs
% if ~(exist(CaSignal.ica_data.Data)
%    
%     
% end
data = CaSignal.ica_data.Data;
V = CaSignal.ica_data.V;
S = CaSignal.ica_data.S;
ICnum = str2num(get(handles.IC_num_edit,'String'));
% CaSignal.ICnum = str2num(get(handles.IC_num_edit,'String'));
CaSignal.ICA_components = run_ICA(CaSignal.ica_data.Data, {S, V, 30, ICnum});
CaSignal.rois_by_IC = cell(1,ICnum);
% CaSignal.ICnum_prev = ICnum;

guidata(handles.figure1, handles);
disp_ICA(handles);


function IC_num_edit_Callback(hObject, eventdata, handles)
runICA_button_Callback(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function IC_num_edit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function disp_ICA(handles)
global CaSignal % ROIinfo ICA_ROIs
RowNum = CaSignal.imSize(1);
ColNum = CaSignal.imSize(2);
if ~ishandle(CaSignal.ICA_figs)
    CaSignal.ICA_figs(1) = figure('Position', [123   460   512   512]);
    CaSignal.ICA_figs(2) = figure('Position',[115    28   512   512]);
end
disp_ICAcomponent_and_blobs(CaSignal.ICA_components(CaSignal.currentIC,:),RowNum, ColNum, CaSignal.ICA_figs);
for i = 1:length(CaSignal.ICA_figs)
    figure(CaSignal.ICA_figs(i)),
    plot_ROIs(handles);
    title(sprintf('IC #%d',CaSignal.currentIC),'FontSize',15);
end

function fig = disp_mean_IC_map(IC)
for i=1:size(IC,1), 
    IC_norm(i,:) = (IC(i,:)- nanmean(IC(i,:)))./ nanstd(IC(i,:)); 
end
IC_norm_mean = nanmax(abs(IC),[],1); % mean(abs(IC_norm),1);
clim = [0  max(IC_norm_mean)*0.7];
fig = figure('Position', [123   372   512   512]);
imagesc(reshape(IC_norm_mean, 128, 512), clim); 
axis square;

function [data_norm_max,fig] = disp_maxDelta_rawData(data)
% each image data has to be already transformed to 1D
% normalize
data_cell = mat2cell(data,ones(1,size(data,1)));
clear data
data_cell_norm = cellfun(@(x) (x-mean(x))./std(x), data_cell, 'UniformOutput',false);
clear data_cell
data_norm = cell2mat(data_cell_norm);
% for i = 1:size(data,1)
%     data_norm(i,:) = (data(i,:) - mean(data(i,:)))./std(data(i,:));
% end
data_norm_max = max(data_norm,[],1);
clim = [0  max(data_norm_max)*0.7];
fig = figure('Position', [100   100   512   512]);
imagesc(reshape(data_norm_max, 128, 512), clim);
axis square;
% --- Executes during object creation, after setting all properties.
function current_ICnum_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to current_ICnum_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in ROI_load_button.
function ROI_load_button_Callback(hObject, eventdata, handles)
global CaSignal


% --------------------------------------------------------------------
function Load_Ca_results_Callback(hObject, eventdata, handles)
global CaSignal
[fn pathstr] = uigetfile('*.mat', 'Load Previous CaTrials results');
if ischar(fn)
    prev_results = load(fullfile(pathstr, fn));
    CaSignal.CaTrials = prev_results.CaTrials;
end


% --------------------------------------------------------------------
function Load_Callback(hObject, eventdata, handles)
% hObject    handle to Load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Load_ROIinfo_Callback(hObject, eventdata, handles)
global CaSignal
[fn pathstr] = uigetfile('*.mat', 'Load saved ROI info');
if ischar(fn)
    load(fullfile(pathstr, fn));
    if iscell(ROIinfo)
        f1 = fieldnames(ROIinfo{1}); f2 = fieldnames(CaSignal.ROIinfo);
        for i = 1:length(ROIinfo)
            for j = 1:length(f1),
                CaSignal.ROIinfo(i).(f2{strcmpi(f2,f1{j})}) = ROIinfo{i}.(f1{j});
            end
        end
    else
        CaSignal.ROIinfo = ROIinfo;
    end
end


% --------------------------------------------------------------------
function Load_ICA_results_Callback(hObject, eventdata, handles)
global CaSignal
[fn pathstr] = uigetfile('*.mat','Load saved ICA results');
if ischar(fn)
    load(fullfile(pathstr, fn)); % load ICA_results
    fprintf('ICA_results of %s loaded!\n', ICA_results.FileBaseName);
    CaSignal.ICA_components = ICA_results.ICA_components;
    CaSignal.currentIC = 1;
    disp_ICA(handles)
end



function current_ICnum_text_Callback(hObject, eventdata, handles)
global CaSignal  
newIC_No = str2num(get(hObject, 'String'));
if newIC_No <= size(CaSignal.ICA_components,1)
    CaSignal.currentIC = newIC_No;
    guidata(hObject, handles);
    disp_ICA(handles);
end


% --- Executes on button press in maxDelta_only_button.
function maxDelta_only_button_Callback(hObject, eventdata, handles)
global CaSignal
% if get(hObject,'Value') == 1
%     [fn, pth] = uigetfile('*.mat','Load Max Delta Image Array');
%     
%     
% end


% --- Executes during object creation, after setting all properties.
function ROI_def_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ROI_def (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over ROI_modify_button.
function ROI_modify_button_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to ROI_modify_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function roiNo_to_plot_Callback(hObject, eventdata, handles)
% hObject    handle to roiNo_to_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of roiNo_to_plot as text
%        str2double(get(hObject,'String')) returns contents of roiNo_to_plot as a double


% --- Executes during object creation, after setting all properties.
function roiNo_to_plot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to roiNo_to_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in check_plotAllROIs.
function check_plotAllROIs_Callback(hObject, eventdata, handles)
% hObject    handle to check_plotAllROIs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_plotAllROIs


% --- Executes when selected object is changed in uipanel1.
function uipanel1_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel1 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in whisker_batch_preprocess.
function whisker_batch_preprocess_Callback(hObject, eventdata, handles)
% hObject    handle to whisker_batch_preprocess (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

basedatapath= get(handles.baseDataPath,'String');
coordsatfirst= str2num(get(handles.coords1,'String'));
coordsatnext= str2num(get(handles.coords2,'String'));
barposatfirst= str2num(get(handles.barpos1,'String'));
barposatnext= str2num(get(handles.barpos2,'String'));
gopos = [];

[gopos,nogopos]=Batch_whiskers_preprocess(basedatapath,coordsatfirst,coordsatnext,barposatfirst,barposatnext);

    set(handles.nogopos,'String',num2str(nogopos/10000));
    set(handles.gopos,'String',num2str(gopos/10000));
% --- Executes on button press in whisker_batch_process.
function whisker_batch_process_Callback(hObject, eventdata, handles)
% hObject    handle to whisker_batch_process (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
basedatapath= get(handles.baseDataPath,'String');
Batch_process_multi_session_whisker_data_GR4(basedatapath);


% --- Executes on button press in make_solo_obj.
function make_solo_obj_Callback(hObject, eventdata, handles)

Solopath = get(handles.SoloDataPath,'String');
mouseName = get(handles.AnimalNameEdit, 'String');
sessionName = get(handles.SoloDataFileName, 'String');
sessionID = get(handles.SoloSessionID, 'String');
if((isempty(Solopath))||(isempty(mouseName))||(isempty(sessionName))||(isempty(sessionID)))
    msgbox('Fill in details');
    return;
end
% trialStartEnd = [get(handles.SoloStartTrialNo,'value'), get(handles.SoloEndTrialNo,'value')];
trialStartEnd = [handles.SoloStartTrial,handles.SoloEndTrial];
obj = Solo.BehavTrialArray_gr(mouseName, sessionName,trialStartEnd,Solopath);
save(['solodata_' mouseName '_' sessionID],'obj');
%% adding solo data to sessObj
behavTrialNums=[handles.SoloStartTrial:handles.SoloEndTrial];
[Solo_data, SoloFileName] = Solo.load_data_gr(mouseName, sessionName,trialStartEnd,Solopath);
behavTrials = {};
for i = 1:length(behavTrialNums)
    behavTrials{i} = Solo.BehavTrial_gr(Solo_data,behavTrialNums(i),1);
end

current_dir = pwd;
separators = find(current_dir == filesep);
session_dir = current_dir(1:separators(length(separators)));
cd (session_dir);
sessObj_found = dir('sessObj.mat');
if isempty(sessObj_found)
    sessObj = {};
    sessObj.behavTrials = behavTrials;
    save('sessObj','sessObj','-v7.3');
else
    load('sessObj.mat');
    sessObj.behavTrials = behavTrials;
    save('sessObj','sessObj','-v7.3');
end
cd (current_dir);


function baseDataPath_Callback(hObject, eventdata, handles)
% hObject    handle to baseDataPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of baseDataPath as text
%        str2double(get(hObject,'String')) returns contents of baseDataPath as a double



% --- Executes during object creation, after setting all properties.
function baseDataPath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to baseDataPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nogopos_Callback(hObject, eventdata, handles)
% hObject    handle to nogopos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nogopos as text
%        str2double(get(hObject,'String')) returns contents of nogopos as a double


% --- Executes during object creation, after setting all properties.
function nogopos_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nogopos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function gopos_Callback(hObject, eventdata, handles)
% hObject    handle to gopos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gopos as text
%        str2double(get(hObject,'String')) returns contents of gopos as a double


% --- Executes during object creation, after setting all properties.
function gopos_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gopos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function [gopos,nogopos]=Batch_whiskers_preprocess(basedatapath,coordsatfirst,coordsatnext,barposatfirst,barposatnext)
    
    d= basedatapath;
    %%load solo_data
    cd ([d '/solo_data']);
    list = dir('solodata_*.mat');
    load(list(1).name);
    solo_obj=obj;
    cd ..
    %% make sessionInfo obj
    sessionInfo = struct('pxPerMm', 31.25, 'mouseName', solo_obj.mouseName,'sessionName',solo_obj.sessionName(end-6:end),'SoloTrialStartEnd',obj.trialStartEnd,'whiskerPadOrigin_nx',[375,70],'whiskerImageDim',[520,400],...
        'bar_coords',[],'bar_time_window',[1.0827,2.5828],'whisker_trajIDs',[0],'theta_kappa_roi',[0,20,140,170],'gotrials',[],'nogotrials',[],'gopos',[],'nogopos',[],'gopix',[],'nogopix',[],'goPosition_mean',[]);


    cd([d '/whisker_data']);
    list = dir('20*_*');
    cd (list(1).name);
    d =pwd;

    [barposmat,barposmatall]= prep(d,solo_obj,coordsatfirst,coordsatnext,barposatfirst,barposatnext);
    save ('barposmat.mat','barposmat');
    save ('barposmatall.mat','barposmatall');
    % cd and save
    cd ..
    cd ..
    % save ('barposmat.mat');
    % list = dir('SessionInfo_*.mat');
    % load(list(1).name);

    %%
    sessionInfo.bar_coords = [];
    sessionInfo.bar_coords = barposmat(:,[1,2]);

    sessionInfo.nogotrials = find(solo_obj.trialTypes==0);
    sessionInfo.gotrials = find(solo_obj.trialTypes==1);

    
    sessionInfo.nogopos=unique(solo_obj.polePositions(sessionInfo.nogotrials));
    nogopos = sessionInfo.nogopos;
   
    sessionInfo.gopos=unique(solo_obj.polePositions(sessionInfo.gotrials));
    gopos = sessionInfo.gopos;
    
    sessionInfo.gopix(:,1) =unique(barposmatall(sessionInfo.gotrials));
     sessionInfo.gopix(:,2) =barposmatall(1,2);
    sessionInfo.nogopix(:,1) = unique(barposmatall(sessionInfo.nogotrials));
      sessionInfo.nogopix(:,2) =barposmatall(1,2);
      
      sessionInfo.goPosition_mean(:,1) =solo_obj.goPosition_mean;
      
    name = ['SessionInfo_' solo_obj.mouseName '_' datestr(now,'mmddyy') ];
    save(name,'sessionInfo');


function [barposmat,barposmatall]= prep(d,solo_obj,coordsatfirst,coordsatnext,barposatfirst,barposatnext)
    %% to make the barposmat conversion
    %d =
    %load('data_*'); %load solodata

    %% evaluate tracker files
    files = dir('*.whiskers');
    temp=struct2cell(files);
    trialnames= temp(1,:)';
    names = char(trialnames);
    if sum(ismember(names(1,:),'g'))>0
        trialno_ind = [29:32];
    else
        trialno_ind=[30:33];
    end
    whisknumind = str2num(names(:,trialno_ind));   

    files = dir('*.measurements');
    temp=struct2cell(files);
    trialnames= temp(1,:)';
    names = char(trialnames);
    measnumind = str2num(names(:,trialno_ind )); 

    files = dir('*.mp4');
    temp=struct2cell(files);
    trialnames= temp(1,:)';
    names = char(trialnames);
    mp4numind = str2num(names(:,trialno_ind )); 
    incomplete =0;

    if (length (whisknumind) == length(measnumind))
        if(length(whisknumind)== length (mp4numind))
            'OK go ahead with analysis'
        else               

            err(['incomplete data:' num2str(length(whisknumind)) 'whiskers and' num2str(length(mp4numind)) 'mp4']);            
            incomplete=1;
        end
    else
            lookat = find((~ismember(measnumind, whisknumind)))
            err(['incomplete data:' num2str(length(whisknumind)) 'whiskers and' num2str(length(measnumind)) 'measurements']);
            incomplete =1;
    end
    
    
    if (incomplete)
        %% get list of files to delete
           deletewhisk = zeros(length(whisknumind),1);
        for i=1:length(whisknumind)
           check = find(measnumind==whisknumind(i));
           if isempty(check)
               deletewhisk(i) = whisknumind(i);
           end
        end


       deletemeas = zeros(length(measnumind),1);
        for i=1:length(measnumind)
           check = find(whisknumind==measnumind(i));
           if isempty(check)
               deletemeas(i) = measnumind(i);
           end
        end
        pause(1)
      end   

    %% make barposmat
    files = dir('*.whiskers');
    temp=struct2cell(files);
    trialnames= temp(1,:)';
    names = char(trialnames);
    trialnumind = str2num(names(:,trialno_ind));

%     pos = saved_history.MotorsSection_motor_position;
    barpos = solo_obj.polePositions';
%     barpos = cell2mat(pos);
    barpos = round(barpos/1000)/10;

    coordsdiff = abs(coordsatfirst-coordsatnext);
    barposdiff = abs(barposatfirst-barposatnext);
    barposmatall = zeros(length(barpos),2);
    factor = coordsdiff/barposdiff;
    barposmatall(:,1) = repmat(coordsatfirst(1),length(barpos),1)- round((barpos(:,1)-barposatfirst(1))*factor(1)) ;
    if(factor(2)~=1)
        barposmatall(:,2) = repmat(coordsatfirst(2),length(barpos),1)-round((barpos(:,1)-barposatfirst(1))*factor(2));
    else
        barposmatall(:,2) = coordsatfirst(2);
    end

    barposmatall=round(barposmatall);
    barposmat = barposmatall(trialnumind,:);
    ['no trials =' num2str(length(trialnumind))]
    ['no. whisker files =' num2str(length(whisknumind))]
    ['no. meas files =' num2str(length(measnumind))]
    ['length of barposmat =' num2str(length(barposmat))]





% --- Executes on button press in secProc.
function secProc_Callback(hObject, eventdata, handles)
    secProc=get(handles.secProc,'Value');
    if(secProc==1)
        set(handles.secProcPanel,'Visible','On');
        set(handles.primProcPanel,'Visible','Off');
        set(handles.tertProcPanel,'Visible','Off');
    end


% --- Executes on button press in tertProc.
function tertProc_Callback(hObject, eventdata, handles)
    tertProc=get(handles.tertProc,'Value');
    if(tertProc==1)
        set(handles.tertProcPanel,'Visible','On');
        set(handles.secProcPanel,'Visible','Off');
        set(handles.primProcPanel,'Visible','Off');
    end
% --- Executes on button press in primProc.
function primProc_Callback(hObject, eventdata, handles)
    primProc=get(handles.primProc,'Value');
    if(primProc==1)
        set(handles.primProcPanel,'Visible','On');
        set(handles.secProcPanel,'Visible','Off');
        set(handles.tertProcPanel,'Visible','Off');
    end


% --- Executes on button press in uibasedatapath.
function uibasedatapath_Callback(hObject, eventdata, handles)
    basedatapath = get(handles.baseDataPath,'String');
    if exist(basedatapath, 'dir')
        cd(basedatapath);
    end
    [basedatapath] = uigetdir(basedatapath,'Set base path');
    set(handles.baseDataPath,'String',basedatapath);


% --- Executes during object creation, after setting all properties.
function numsolotrials_CreateFcn(hObject, eventdata, handles)



function CaSignal_datapath_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function CaSignal_datapath_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in load_CaSignal.
function load_CaSignal_Callback(hObject, eventdata, handles)
    global CaTrials
    global sorted_CaTrials
    global contact_CaTrials 
    basedatapath = get(handles.CaSignal_datapath,'String');
    cd (basedatapath);
    
    if(get(handles.CaTrials_select,'Value')==1)   
        [filename,pathName]=uigetfile('CaSignal*.mat','Load CaSignal.mat file')
        if isequal(filename, 0) || isequal(pathName,0)
            return
        end
        
        set(handles.CaSignal_datapath,'String',pathName);
        cd(pathName);
        load( [pathName filesep filename], '-mat');
         
        if (strcmp(filename(1:6),'sorted'))
            msgbox('Change settings before loading sorted_CaTrials');
            return
        end
        numrois = CaTrials(1,1).nROIs;
    elseif (get(handles.sorted_CaTrials_select,'Value')==1)              
        [filename,pathName]=uigetfile('sorted_CaSignal*.mat','Load sorted_CaSignal.mat file')
        if isequal(filename, 0) || isequal(pathName,0)
            return
        end
        set(handles.CaSignal_datapath,'String',pathName);
        cd(pathName);
        load( [pathName filesep filename], '-mat');
      
        [filename,pathName]=uigetfile('CaSignal*.mat','Also Load CaSignal.mat file')
        if isequal(filename, 0) || isequal(pathName,0)
            return
        end
        load( [pathName filesep filename], '-mat');
        numrois = CaTrials(1,1).nROIs;
    elseif(get(handles.contact_CaSignal_select,'Value')==1)
        
        [filename,pathName]=uigetfile('contact_CaSignal*.mat','Load contact_CaSignal.mat file')
        if isequal(filename, 0) || isequal(pathName,0)
            return
        end
        set(handles.CaSignal_datapath,'String',pathName);
        cd(pathName);
        load( [pathName filesep filename], '-mat');
%         numrois = size(contact_CaTrials{1,1},1);
        numrois = contact_CaTrials(1).nROIs;
    end
    

    set(handles.CaSignalrois,'String',num2str(numrois));



function fov_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function fov_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function roislist_Callback(hObject, eventdata, handles)



% --- Executes during object creation, after setting all properties.
function roislist_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in plotroiSignals.
function plotroiSignals_Callback(hObject, eventdata, handles)
global CaTrials 
global sorted_CaTrials
fov = get(handles.fov,'String');
roislist = get(handles.roislist,'String');
rois = str2num(roislist);

if(isempty(CaTrials))   
    [filename,pathName]=uigetfile('CaSignal*.mat','Load CaSignal.mat file')
    if isequal(filename, 0) || isequal(pathName,0)
        'No CaTrial loaded!'
         return
    end
    load( [pathName filesep filename], '-mat');
end

if(isempty(sorted_CaTrials))   
    [filename,pathName]=uigetfile('sorted_CaSignal*.mat','Load sorted_CaSignal.mat file')
    if isequal(filename, 0) || isequal(pathName,0)
        'No sorted_CaTrial loaded!'
         return
    end
    load( [pathName filesep filename], '-mat');
end

if (length(rois) > CaTrials(1).nROIs)
    rois = 1:CaTrials(1).nROIs;
    set(handles.roislist,'String',num2str(rois));
end
if(get(handles.CaTrials_select,'Value') ==1)  
    tag_trialtypes =0;
    sfx='Unsort';
    trialtypes = ones(length(CaTrials),1);
    plot_roiSignals(CaTrials,fov,rois,roislist,tag_trialtypes,trialtypes,sfx);
elseif(get(handles.sorted_CaTrials_select,'Value') ==1)
%     global sorted_CaTrials    
    tag_trialtypes =1;
    sfx = 'Bsort';
    count =0;
    trialorder  = [sorted_CaTrials.hits , sorted_CaTrials.misses , sorted_CaTrials.cr, sorted_CaTrials.fa]; 
    disp('order = hits misses cr  fa')
    trialtypes = zeros(length(trialorder),1);
    trialtypes(count+1:count+length(sorted_CaTrials.hits )) = 1;
    count = count +length(sorted_CaTrials.hits);
  
    if(isempty(sorted_CaTrials.misses))
        'No Miss Trials in this list'
    else
       trialtypes(count+1:count+length(sorted_CaTrials.misses )) = 2;
       count = count +length(sorted_CaTrials.misses );
    end
 
    if(isempty(sorted_CaTrials.cr))
        'No Correct Rejection Trials in this list'
   
    else
       trialtypes(count+1:count+length(sorted_CaTrials.cr )) = 3; 
        count = count +length(sorted_CaTrials.cr );
    end

    
    if(isempty(sorted_CaTrials.fa))
        'No False alarm Trials in this list'
    else
       trialtypes(count+1:count+length(sorted_CaTrials.fa )) = 4; 
    end
    
    plot_roiSignals(CaTrials(trialorder),fov ,rois,roislist,tag_trialtypes,trialtypes,sfx);
    
    if(get(handles.sortwtouchInfo,'Value')==1)
        %% sort by trial time in session
        tag_trialtypes =1;
        sfx = 'TSort';
        count =0;
        trialorder  = [sorted_CaTrials.touch , sorted_CaTrials.notouch]; 
        disp('order = touch notouch')
        trialtypes = zeros(length(trialorder),1);
        trialtypes(count+1:count+length(sorted_CaTrials.touch )) = 1;
        count = count +length(sorted_CaTrials.touch);        
        trialtypes(count+1:count+length(sorted_CaTrials.notouch )) = 3;
        count = count +length(sorted_CaTrials.notouch );   
        plot_roiSignals(CaTrials(trialorder),fov ,rois,roislist,tag_trialtypes,trialtypes,sfx);
        
        %% sort by bar_pos_trial
        tag_trialtypes =1;
        sfx = 'TSort_barpos';
        count =0;
        trialtypes = zeros(length(trialorder),1);
        touchinds = zeros(length(sorted_CaTrials.touch_barpos),1);
        notouchinds = zeros(length(sorted_CaTrials.notouch_barpos),1);
        
        barpositions1 = [ unique(sorted_CaTrials.touch_barpos')];       
        for i=1:length(barpositions1)            
          inds= find(sorted_CaTrials.touch_barpos==barpositions1(i)) ;
          touchinds(count+1:count+length(inds))=inds;
          trialtypes(count+1:count+length(inds)) = i;         
          count = count +length(inds);
        end
        offset = count;
        count=0;
        barpositions2 = unique(sorted_CaTrials.notouch_barpos');       
        for j=1:length(barpositions2)            
          inds= find(sorted_CaTrials.notouch_barpos==barpositions2(j)) ;
          notouchinds(count+1:count+length(inds))=inds;
          trialtypes(offset+count+1:offset+count+length(inds)) = find(barpositions1==barpositions2(j));         
          count = count +length(inds);
        end
      
        trialorder  = [sorted_CaTrials.touch(touchinds) , sorted_CaTrials.notouch(notouchinds)]; 
        disp('order = touch_barpos_Ant(Top)-Post(Bottom) notouch_barpos_Ant(Top)-Post(Bottom)')       
        plot_roiSignals(CaTrials(trialorder),fov ,rois,roislist,tag_trialtypes,trialtypes,sfx);

    end
else(get(handles.contact_CaSignal_select,'Value')==1)
   
    global contact_CaTrials 
    tag_trialtypes =0;
    sfx ='Csort';
    trialtypes = ones(length(contact_CaTrials),1);
    plot_roiSignals(contact_CaTrials,fov,rois,roislist,tag_trialtypes,trialtypes,sfx);
    


     %% sort by bar_pos_trial
    tag_trialtypes =1;
    sfx = 'CSort_barpos';
    count =0;
    trialtypes = ones(length(contact_CaTrials),1);

    touchinds = zeros(length(contact_CaTrials),1);        
    barposall=cell2mat(arrayfun(@(x) x.barpos, contact_CaTrials, 'uniformoutput',false));
    barpositions1 = unique(barposall);
    for i=1:length(barpositions1)            
      inds= find(barposall==barpositions1(i)) ;
      touchinds(count+1:count+length(inds))=inds;
      trialtypes(count+1:count+length(inds)) = i;         
      count = count +length(inds);
    end

%         trialorder  = [sorted_CaTrials.touch(touchinds) ];%, sorted_CaTrials.notouch(notouchinds)]; 
     disp('order = touch_barpos_Ant(Top)-Post(Bottom)');% notouch_barpos_Ant(Top)-Post(Bottom)')      

    plot_roiSignals(contact_CaTrials(touchinds),fov,rois,roislist,tag_trialtypes,trialtypes,sfx);
        
        
end



function wSig_datapath_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function wSig_datapath_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in loadwSigsessionInfo.
function loadwSigsessionInfo_Callback(hObject, eventdata, handles)
    global sessionInfo
    global wSigTrials
    basedatapath = get(handles.wSig_datapath,'String');
    cd (basedatapath);
    [filename,pathName]=uigetfile('wSig*.mat','Load wSig.mat file');
    if isequal(filename, 0) || isequal(pathName,0)
        return
    end
    set(handles.wSig_datapath,'String',pathName);
    cd(pathName);
    load( [pathName filesep filename], '-mat');
    [filename,pathName]=uigetfile('SessionInfo*.mat','Load SessionInfo.mat file');
    load( [pathName filesep filename],'-mat');


% --- Executes on button press in plotwSigData.
function plotwSigData_Callback(hObject, eventdata, handles)
global sessionInfo
global wSigTrials
    gopix = sessionInfo.gopix;
    nogopix =  sessionInfo.nogopix;  
    numgotrials = zeros(size(gopix,1),1);
    numnogotrials = [];
    gotrials = sessionInfo.gotrials;
    nogotrials = sessionInfo.nogotrials;
    wAmp = zeros(4,2);
    barposmat = sessionInfo.bar_coords;

 % get trial names 
    if sum(ismember(wSigTrials{1}.trackerFileName,'gr')) >0
        wSigfilenames =char(cellfun(@(x) x.trackerFileName(29:32),wSigTrials,'uniformoutput',false));
    else
        wSigfilenames =char(cellfun(@(x) x.trackerFileName(30:33),wSigTrials,'uniformoutput',false));

    end
    wSigtrials=str2num(wSigfilenames);
     
    p =  pwd;
 %% get trial blocks
    numblocks = str2num(get(handles.numblocks,'String'));
    trialblocks=strsplit(get(handles.trialswindow,'String'));
    tags =strsplit(get(handles.trialswindow_tag,'String'));
    numblocks = size(tags,2);
    if(numblocks==1)
        tags={tags};
        trialblocks={trialblocks};
    end
%     gotrialstoplot = [];
%     nogotrialstoplot=[];
    blocks = struct('tag','','gotrialnums',[],'nogotrialnums',[],'gotrialnames',[],'nogotrialnames',[]);
   
    for i=1:numblocks
       blocks.tag{i}= tags{i};
       temp = str2num(trialblocks{i});
       ind= find(ismember(gotrials,temp(1,:)));
       blocks.gotrialnames{i}=gotrials(ind);
       ind = find(ismember(nogotrials,temp(1,:))); 
       blocks.nogotrialnames{i} =  nogotrials(ind);  
       
       
    nogotrialinds = zeros(length(blocks.nogotrialnames{i} ),1);
    gotrialinds = zeros(length(blocks.gotrialnames{i}),1);
    for k=1:length(blocks.nogotrialnames{i})
        if (ismember(blocks.nogotrialnames{i}(k),wSigtrials))
            nogotrialinds(k)=find(wSigtrials==blocks.nogotrialnames{i}(k));
        end
    end

     for k=1:length(blocks.gotrialnames{i})
        if (ismember(blocks.gotrialnames{i}(k),wSigtrials))
            gotrialinds(k)=find(wSigtrials==blocks.gotrialnames{i}(k));
        end
     end

        nogotrialinds(nogotrialinds==0)=[];
        gotrialinds(gotrialinds==0)=[];
        blocks.nogotrialnums{i}=nogotrialinds;
        blocks.gotrialnums{i}=gotrialinds;
       
    end

     avg_trials = str2num(get(handles.wh_trialstoavg,'string')); %% these many from last
     restrictTime = str2num(get(handles.timewindow_wSiganal,'String'));
     plot_whiskerfits = get(handles.plot_fittedwhisker,'Value');
     timewindowtag = cell2mat(strcat('Trials',strrep(trialblocks,':','_'),'pole',get(handles.timewindowtag,'String')));
  
% %      [w_setpoint_trials,w_setpoint_trials_med,w_setpoint_trials_medbinned,w_setpoint_trials_width,w_setpoint_trials_dur,w_setpoint_trials_prewhisk,w_setpoint_trials_prewhiskbinned,pvalsetpoint,...
% %          w_amp_trials,w_thetaenv_trials,w_thetaenv_dist,w_amp_trials_med,w_amp_trials_medbinned,w_amp_trials_width,pvalamp]= wdatasummary(wSigTrials,blocks.tag,blocks.nogotrialnums,avg_trials,gopix,nogopix,restrictTime,p,plot_whiskerfits,'nogo',timewindowtag);
        [w_thetaenv] =  wdatasummary(wSigTrials,blocks.tag,blocks.nogotrialnums,avg_trials,gopix,nogopix,restrictTime,p,plot_whiskerfits,'nogo',timewindowtag);

    for i=1:numblocks

%          blocks.nogo_thetaenv_med{i} = w_thetaenv.med{i};
         blocks.nogo_thetaenv_binned{i} = w_thetaenv.binned{i};
%          blocks.nogo_thetaenv_width{i} = w_thetaenv.width{i};
%          blocks.nogo_thetaenv_dur{i} =w_thetaenv.dur{i};
%          blocks.nogo_thetaenv_prepole{i} =w_thetaenv.prepole{i};
%          blocks.nogo_thetaenv_prepolebinned{i} =w_thetaenv.prepolebinned{i};
         blocks.nogo_thetaenv_pval= w_thetaenv_pval;
         blocks.nogo.thetaenv_trials{i} = w_thetaenv.trials{i};
          blocks.nogo.thetaenv_binned_dist{i} = w_thetaenv_binned_dist{i};      
     end
         [w_thetaenv] =  wdatasummary(wSigTrials,blocks.tag,blocks.gotrialnums,avg_trials,gopix,nogopix,restrictTime,p,plot_whiskerfits,'go',timewindowtag);
 
    for i=1:numblocks

%          blocks.go_thetaenv_med{i} = w_thetaenv.med{i};
         blocks.go_thetaenv_binned{i} = w_thetaenv.binned{i};
%          blocks.go_thetaenv_width{i} = w_thetaenv.width{i};
%          blocks.go_thetaenv_dur{i} =w_thetaenv.dur{i};
%          blocks.go_thetaenv_prepole{i} =w_thetaenv.prepole{i};
%          blocks.go_thetaenv_prepolebinned{i} =w_thetaenv.prepolebinned{i};
         blocks.go_thetaenv_pval= w_thetaenv.pval;
         blocks.go.thetaenv_trials{i} = w_thetaenv.trials{i};
         blocks.go.thetaenv_binned_dist{i} = w_thetaenv_binned_dist{i};     
     end
     fpath = pwd;
     fpath = [fpath filesep 'plots' filesep timewindowtag];
     list = dir('wSigTrials*.mat');
     fname = list(1).name;   
     fname = strrep(fname,'wSigTrials','wSigBlocks');
     blocks.fname = fname;
     fname = [fpath filesep fname];
     save(fname,'blocks');

% --- Executes on button press in zstack_avg.
function zstack_avg_Callback(hObject, eventdata, handles)

    [filename,path] = uigetfile('*.tif','Pick the zstack file');
    cd(path);
    finfo = imfinfo(filename); 
    if isfield(finfo, 'ImageDescription')
       [header] = extern_scim_opentif(filename); 
       n_channel = length(header.SI4.channelsSave);
    %     header.width = header.acq.pixelsPerLine;
    %     header.height = header.acq.linesPerFrame;
        header.n_frame = length(finfo);
    end
    n_channel=2;
    header.width = finfo(1).Width;
    header.height = finfo(1).Height;
    channel= get(handles.zstack_channel,'String');
    if n_channel > 1
        if strncmpi(channel, 'g', 1)
            firstframe = 1;
            step = n_channel;
        elseif strncmpi(channel, 'r', 1)
            firstframe = 2;
            step = n_channel;
        else
            error('unknown channel name?')
        end
    else
        firstframe = 1;
        step = 1;
    end

    im = zeros(header.height, header.width, header.n_frame/step, 'uint16');

    count = 0;
    for i = firstframe : step : length(finfo)
        count = count+1;
        im (:,:,count) = imread(filename, i);
    end



    nframes = str2num(get(handles.zstack_nframes,'String') );
    temp = uint16(zeros(size(im,1),size(im,2),size(im,3)/nframes));
    % temp = zeros(size(im,1),size(im,2),nframes);
    im_avg=uint16(zeros(size(im,1),size(im,2),size(im,3)/nframes));
    count=0;
    % for i = 1:size(im,3)/nframes
        for i = 1:nframes
         temp = im(:,:,count+1:count+size(im,3)/nframes);
         im_avg =im_avg+ temp;
         count=count+size(im,3)/nframes;
    % % %    temp = im(:,:,count+1:count+nframes);
    % % %    im_avg(:,:,i) = mean(temp,3);
    % % %     count=count+nframes;
    end
    im_avg(:,:,:)=im_avg(:,:,:)./nframes;
    imwrite(uint16(im_avg(:,:,1)),colormap(gray), ['zstack' filename(length(filename)-7 : length(filename)-4) '.tiff'])
    for i = 2:size(im,3)/nframes
     imwrite(uint16(im_avg(:,:,i)),colormap(gray), ['zstack' filename(length(filename)-7 : length(filename)-4) '.tiff'], 'writemode', 'append')
    end


function zstack_channel_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function zstack_channel_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function zstack_nframes_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function zstack_nframes_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function Image_disp_axes1_CreateFcn(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function export2avi_button_CreateFcn(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function SaveFrameButton_CreateFcn(hObject, eventdata, handles)



% --- Executes during object creation, after setting all properties.
function Image_disp_axes2_CreateFcn(hObject, eventdata, handles)




function barpos1_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function barpos1_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function barpos2_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function barpos2_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function coords1_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function coords1_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function coords2_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function coords2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to coords2 (see GCBO)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function solo_datapath_Callback(hObject, eventdata, handles)


% --- Executes during object creation, af0ter setting all properties.
function solo_datapath_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in load_solo_obj.
function load_solo_obj_Callback(hObject, eventdata, handles)
    global solo_data
    basedatapath = get(handles.solo_datapath,'String');
    cd (basedatapath);
    [filename,pathName]=uigetfile('solo_data*.mat','Load solo_data.mat file')
    if isequal(filename, 0) || isequal(pathName,0)
        return
    end
    set(handles.solo_datapath,'String',pathName);
    cd(pathName);
   load( [pathName filesep filename], '-mat');
   solo_data = obj;
    


% --- Executes on button press in uigetsolo.
function uigetsolo_Callback(hObject, eventdata, handles)

[filename,pathName]=uigetfile('data*.mat','Load behavior data .mat file');
    if isequal(filename, 0) || isequal(pathName,0)
        return
    end
set(handles.SoloDataPath,'String',pathName);
set(handles.SoloDataFileName,'String',filename(1 :length(filename)-4));
set(handles.SoloSessionID,'String',filename(length(filename)-10 :length(filename)-4));
if sum(ismember(filename,'gr'))>2
    set(handles.AnimalNameEdit, 'String',filename(length(filename)-19 :length(filename)-12));
else
    set(handles.AnimalNameEdit, 'String',filename(length(filename)-20 :length(filename)-12));
end
load([pathName filesep filename]);
set(handles.SoloEndTrialNo, 'String', saved.AnalysisSection_NumTrials);

[imaging_datapath] = uigetdir(pathName,'Setimaging data path');
    if isequal(imaging_datapath, 0) 
    else
        cd (imaging_datapath);
        list = dir('Image*.tif');
        filenames =cell2mat(arrayfun(@(x) x.name(length(x.name)-6 :length(x.name)-4),list,'uniformoutput',false));
        trials=str2num(filenames);
        alltrials = 1:saved.AnalysisSection_NumTrials;
        excluded=alltrials;
        excluded(find(ismember(alltrials,trials)))=[];
        set(handles.behavTrialNoToBeExcluded, 'String',num2str(excluded));
    end
cd (pathName);




% --- Executes during object creation, after setting all properties.
function CaSignalrois_CreateFcn(hObject, eventdata, handles)


% --- Executes on button press in sort_CaSignal.
function sort_CaSignal_Callback(hObject, eventdata, handles)
global CaTrials

if(get(handles.sortwtouchInfo,'Value')==1)
    global sessObj
    [filename,pathname]=uigetfile('wSigTrials*.mat','Load wSigTrials.mat file');
    load([pathname filesep filename]);  
end

if(isempty(CaTrials))
   msgbox( 'Load Casignal data');
   return;
elseif~(isfield(CaTrials,'behavTrial'))
   msgbox( 'Add behav trials (Prim panel)');
   return;
else
    trialCorrects = arrayfun(@(x) x.behavTrial.trialCorrect, CaTrials,'uniformoutput',false);
    [sorted,sortedtrialinds] = sort(cell2mat(trialCorrects),'descend');
   

    correcttrials = sortedtrialinds(find(sorted==1)); %% contains hits and correct rejections trial inds
    incorrecttrials = sortedtrialinds(find(sorted==0));%% xontains misses and false positives trial inds
    
    correcttrialstypes = arrayfun(@(x) x.behavTrial.trialType, CaTrials(correcttrials),'uniformoutput',false);
    incorrecttrialstypes = arrayfun(@(x) x.behavTrial.trialType, CaTrials(incorrecttrials),'uniformoutput',false);
    
%     sortedtrialinds = zeros(1,length(sortedtrialinds));   
%     [sorted,sortedtrialinds1] = sort(cell2mat(correcttrialstypes),'descend');     
%     sorted_CaTrials.hits = correcttrials(find(sorted==1));
%     sorted_CaTrials.cr = correcttrials(find(sorted==0));
%      sortedtrialinds(1:length(sortedtrialinds1)) = correcttrials(sortedtrialinds1);
%      [sorted,sortedtrialinds2] = sort(cell2mat(incorrecttrialstypes),'descend');
%     sorted_CaTrials.misses = incorrecttrials(find(sorted==1));
%     sorted_CaTrials.fa = incorrecttrials(find(sorted==0));
%      sortedtrialinds(length(sortedtrialinds1)+1:length(sortedtrialinds1)+length(sortedtrialinds2)) = incorrecttrials(sortedtrialinds2);
%     
     sortedtrialinds(1:end) = 0;
      sorted_CaTrials.hits = correcttrials(find(cell2mat(correcttrialstypes)==1));
    sorted_CaTrials.cr = correcttrials(find(cell2mat(correcttrialstypes)==0));
  
    sorted_CaTrials.misses = incorrecttrials(find(cell2mat(incorrecttrialstypes)==1));
    sorted_CaTrials.fa = incorrecttrials(find(cell2mat(incorrecttrialstypes)==0));
     
if(get(handles.sortwtouchInfo,'Value')==1)
    names=arrayfun(@(x) x.FileName(length(x.FileName)-6:length(x.FileName)-4), CaTrials,'uniformoutput',false);%%FileName_prefix removes everything but the trial counter
    CaSig_trialnums = str2num(char(names));%from trial filenames
    names=cellfun(@(x) x.trackerFileName(length(x.trackerFileName)-21:length(x.trackerFileName)-18),wSigTrials,'uniformoutput',false);
    wSig_trialnums =str2num(char(names)); % from trial filenames
    common_trialnums=CaSig_trialnums(find(ismember(CaSig_trialnums,wSig_trialnums)));
    temp=wSig_trialnums(find(ismember(wSig_trialnums,CaSig_trialnums)));
    if(size(temp,1)<size(common_trialnums,1))
        common_trialnums=temp;
    end
  
    tags = find(ismember(CaSig_trialnums,common_trialnums));
    CaSig_tags = tags;
    ind = zeros(length(tags),1);
    ind (find(ismember(tags,sorted_CaTrials.hits)))=1;
    ind (find(ismember(tags,sorted_CaTrials.cr)))=3;
    ind (find(ismember(tags,sorted_CaTrials.misses)))=2;
    ind (find(ismember(tags,sorted_CaTrials.fa)))=4;
    
    tags = find(ismember(wSig_trialnums,common_trialnums));
    wSig_tags=tags;
    trialnums=common_trialnums;% wrt CaTrials indices
    trialinds = 1:length(trialnums);
    whiskerID =1;
    contacttimes=cellfun(@(x) x.contacts{whiskerID}, wSigTrials(wSig_tags(trialinds)),'uniformoutput',false);
    notouchtrials = find(cellfun(@isempty,contacttimes));    
    %sorted_CaTrials.touch
    sorted_CaTrials.notouch =CaSig_tags(notouchtrials)'; 
    sorted_CaTrials.notouch_barpos =cell2mat(cellfun(@(x) x.bar_pos_trial(1,1), wSigTrials(wSig_tags(notouchtrials)),'uniformoutput',false)); % Y position alone A-P

    touchtrials = find(~cellfun(@isempty,contacttimes));
    sorted_CaTrials.touch =CaSig_tags(touchtrials)';
    sorted_CaTrials.touch_barpos=cell2mat(cellfun(@(x) x.bar_pos_trial(1,1), wSigTrials(wSig_tags(touchtrials)),'uniformoutput',false));

% % %     for i=1:length(trialinds)
% % %     
% % %             TrialName = strrep(CaTrials(CaSig_tags(i)).FileName,CaTrials(CaSig_tags(i)).FileName_prefix,'');
% % %             TrialName=strrep(TrialName,'.tif','');
% % %             sessObj(i).dff = CaTrials(CaSig_tags(i)).dff;
% % %             sessObj(i).ts = {[1:size(CaTrials(CaSig_tags(i)).dff,2)]*CaTrials(CaSig_tags(i)).FrameTime};
% % %             sessObj(i).theta = wSigTrials{wSig_tags(i)}.theta;
% % %             sessObj(i).kappa = wSigTrials{wSig_tags(i)}.kappa;
% %             sessObj(i).velocity = wSigTrials{wSig_tags(i)}.Velocity;
% %             sessObj(i).deltaKappa = wSigTrials{wSig_tags(i)}.deltaKappa;
% %             sessObj(i).ts_wsk = wSigTrials{wSig_tags(i)}.time;
% %             sessObj(i).contactdir = wSigTrials{wSig_tags(i)}.contact_direct;
% %             sessObj(i).FrameTime = CaTrials(CaSig_tags(i)).FrameTime;
% %             sessObj(i).nframes = CaTrials(CaSig_tags(i)).nFrames;
% %             sessObj(i).CaSigTrialind=CaTrials(CaSig_tags(i)).TrialNo;
% %             sessObj(i).TrialName= TrialName;
% %             sessObj(i).nROIs = CaTrials(CaSig_tags(i)).nROIs;
% %             sessObj(i).contacts=wSigTrials{wSig_tags(i)}.contacts;
% %             sessObj(i).barpos = cellfun(@(x) x.bar_pos_trial(1,1), wSigTrials(wSig_tags(touchtrials)),'uniformoutput',false);
% %             switch(ind(i))
% %                 case(1)
% %                     sessObj(i).trialtype = 'Hit';               
%                 case(2)
%                     sessObj(i).trialtype = 'Miss';
%                 case(3)
%                     sessObj(i).trialtype = 'CR';                   
%                 case(4)
%                     sessObj(i).trialtype = 'FA';
%   
%             end
%             if( isfield(CaTrials, 'ephusTrial')) 
%                 sessObj(i).licks = CaTrials(CaSig_tags(i)).ephusTrial.licks;
%                 sessObj(i).poleposition = CaTrials(i).ephusTrial.poleposition;
%                 sessObj(i).ephuststep = CaTrials(i).ephusTrial.ephuststep;
%             end
%     end
    dirname = pwd;
    sessdir = strrep(dirname,'/fov_01001/fluo_batch_out/',''); %% save under roi folder 
% %     roiname = strrep(sessname,'/fov_01001/fluo_batch_out/',''); 
%     save('sessObj','sessObj');
end
     save('sorted_CaTrials','sorted_CaTrials');
end


% --- Executes on button press in CaTrials_select.
function CaTrials_select_Callback(hObject, eventdata, handles)
set(handles.CaTrials_select,'Value',1);
set(handles.sorted_CaTrials_select,'Value',0);
set(handles.contact_CaSignal_select,'Value',0);
% --- Executes on button press in sorted_CaTrials_select.
function sorted_CaTrials_select_Callback(hObject, eventdata, handles)
set(handles.sorted_CaTrials_select,'Value',1);
set(handles.CaTrials_select,'Value',0);
set(handles.contact_CaSignal_select,'Value',0);

% --- Executes on button press in compute_contact_CaSignal.
function compute_contact_CaSignal_Callback(hObject, eventdata, handles)
global wSigTrials
global CaTrials
global sorted_CaTrials
global contact_CaTrials 
global sessionObj
framerate = 500; % Hz
whiskerID =1;
 S = {'Protract';'Retract'};
[Selection,ok] = listdlg('PromptString','Select direction of touch','ListString',S,'SelectionMode','single','ListSize',[160,100])

selected_contact_direct = S{Selection};
if(isempty(wSigTrials))
    [filename,pathname]=uigetfile('wSigTrials*.mat','Load wSigTrials.mat file');
    load([pathname filesep filename]);
end
    
if(isempty(CaTrials))
    [filename,pathname]=uigetfile('CaTrials*.mat','Load CaTrials.mat file');
    load([pathname filesep filename]);
end

if(isempty(sorted_CaTrials))
    [filename,pathname]=uigetfile('sorted_CaTrials*.mat','Load sorted_CaTrials.mat file');
    load([pathname filesep filename]);
end

names=arrayfun(@(x) x.FileName(length(x.FileName)-6:length(x.FileName)-4), CaTrials,'uniformoutput',false);%%FileName_prefix removes everything but the trial counter
CaSig_trialnums = str2num(char(names));%from trial filenames
names=cellfun(@(x) x.trackerFileName(length(x.trackerFileName)-21:length(x.trackerFileName)-18),wSigTrials,'uniformoutput',false);
wSig_trialnums =str2num(char(names)); % from trial filenames
common_trialnums=CaSig_trialnums(find(ismember(CaSig_trialnums,wSig_trialnums)));
if(size(wSig_trialnums,1)<size(common_trialnums,1))
    common_trialnums=CaSig_trialnums(find(ismember(wSig_trialnums,CaSig_trialnums)));
end

tags = find(ismember(CaSig_trialnums,common_trialnums));
ind = zeros(length(tags),1);
ind (find(ismember(tags,sorted_CaTrials.hits)))=1;
ind (find(ismember(tags,sorted_CaTrials.cr)))=3;
ind (find(ismember(tags,sorted_CaTrials.misses)))=2;
ind (find(ismember(tags,sorted_CaTrials.fa)))=4;
[Y,indorder] = sort(ind,'ascend');
% nogotrials = [indorder(Y==3); indorder(Y==4)];
CaSig_tags = tags(indorder);
tags = find(ismember(wSig_trialnums,common_trialnums));
wSig_tags=tags(indorder);

trialnums=common_trialnums(indorder);% wrt CaTrials indices
trialinds = 1:length(trialnums);
%% removing contacttimes outside bar available time
 contacttimes=cellfun(@(x) x.contacts{whiskerID}, wSigTrials(wSig_tags(trialinds)),'uniformoutput',false);


for i = 1: length(trialnums)
           barOntime= CaTrials(CaSig_tags(i)).behavTrial.pinDescentOnsetTime;
           barOfftime=CaTrials(CaSig_tags(i)).behavTrial.pinAscentOnsetTime;
           firsttouches = contacttimes{i};
%            extraneouscontacts = contacttimes
end



%% removing all trials with no contact (try plotting only for this later)
    contacttimes=cellfun(@(x) x.contacts{whiskerID}, wSigTrials(wSig_tags(trialinds)),'uniformoutput',false);
     numtrials = size(contacttimes,2);
    nocontacts = find(cellfun(@isempty,contacttimes));
    trialinds(nocontacts) =[];   
    
    contacttimes=cellfun(@(x) x.contacts{whiskerID}, wSigTrials(wSig_tags(trialinds)),'uniformoutput',false);
    thetavals=cellfun(@(x) x.theta{whiskerID}, wSigTrials(wSig_tags(trialinds)),'uniformoutput',false);  
    kappavals=cellfun(@(x) x.kappa{whiskerID}, wSigTrials(wSig_tags(trialinds)),'uniformoutput',false);  
    Velocity=cellfun(@(x) x.Velocity{whiskerID}, wSigTrials(wSig_tags(trialinds)),'uniformoutput',false);
    Setpoint=cellfun(@(x) x.Setpoint{whiskerID}, wSigTrials(wSig_tags(trialinds)),'uniformoutput',false);  
    Amplitude=cellfun(@(x) x.Amplitude{whiskerID}, wSigTrials(wSig_tags(trialinds)),'uniformoutput',false); 
    deltaKappa=cellfun(@(x) x.deltaKappa{whiskerID}, wSigTrials(wSig_tags(trialinds)),'uniformoutput',false);  
    ts_wsk=cellfun(@(x) x.time{whiskerID}, wSigTrials(wSig_tags(trialinds)),'uniformoutput',false);  
    contactdir=cellfun(@(x) x.contact_direct{whiskerID}, wSigTrials(wSig_tags(trialinds)),'uniformoutput',false); 
    contacts=cellfun(@(x) x.contacts{whiskerID}, wSigTrials(wSig_tags(trialinds)),'uniformoutput',false); 
    trialnums(nocontacts) = [];
    wSig_tags(nocontacts)=[];
    CaSig_tags(nocontacts)=[];
    Y(nocontacts)=[];
 
%% get CaSig data from sorted_CaTrials.CaTRials
    CaTrials_data  = arrayfun(@(x) x.dff,CaTrials(CaSig_tags),'uniformoutput',false);
    numtrials = length(CaTrials_data);
    numrois = size(CaTrials_data{1},1);
    numframes =size(CaTrials_data{1},2);
    CaSig = zeros(numrois,numframes,numtrials);
    for i= 1:numtrials
            tempmat = zeros(size(CaTrials_data,1),size(CaTrials_data,2));
            tempmat = CaTrials_data{i};
          if (size(tempmat,2)>size(CaSig,2))
           CaSig(:,:,i) = tempmat(1:numrois,1:numframes);
          else
            CaSig(:,1:size(tempmat,2),i) = tempmat(:,:);
          end      
    end
 
    Caframetime = CaTrials.FrameTime;
    baseline = 0.5;
    dur = 3.0;
    wSigframerate = 500;
    numpts=ceil((dur+baseline)*wSigframerate);% 3.5seconds worth of data
    numframes = ceil((dur+baseline)/Caframetime);% 2.5seconds worth of data
    
    numcontacts =0;
    contact_CaTrials=struct('solo_trial',[],'dff',{},'ts',{},'FrameTime',{},'nframes',{},'trialtype',[],'trialCorrect',[],'FileName_prefix',{},'FileName',{},...
                                     'TrialName',{},   'licks', {},'poleposition',{},'nROIs',{},'theta',{},'kappa',{},...
                                   'deltaKappa',{},'ts_wsk',{},'contactdir',{},'contacts',{},'barpos',[],'Setpoint',{},'Amplitude',{},'Velocity',{});
    count=0;
    handles.aligned_contact = get(handles.align_to_first_touch,'Value'); % 1 for first 
    
    for i = 1:numtrials
            allcontacts = size(contacttimes{i},2);
            if get(handles.align_to_first_touch,'Value')
                 pickedcontact=1; %first contact
                 contact_CaSig_tag = 'Ftouch';
            elseif get(handles.align_to_last_touch,'Value');
                 pickedcontact=allcontacts;%last contact
                 contact_CaSig_tag = 'Ltouch';
            end
            numcontacts = 1;
            contactind = zeros(numcontacts,1);
            for j=1:numcontacts
                 ind= 1; %first ind
                
%                   ind= length(contacttimes{i}{1,pickedcontact(j)}); %last ind
                temp=contacttimes{i}{1,pickedcontact(j)}(ind);%for now just the first contact after bartime
                if(temp/wSigframerate<=baseline)&&(allcontacts>1)                   
                    contactind(j) = contacttimes{i}{1,(j+1)}(1);
                elseif(temp/wSigframerate<=baseline)&&(allcontacts==1)
                    contactind(j) = ceil(baseline*wSigframerate)+1;
                else
                    contactind(j) =temp;% for now just the first contact
                end
            % for now just the first contact
            end
         extractedCaSig = zeros(numrois,numframes);
         extractedTheta=zeros(1,numpts);
         extractedKappa=zeros(1,numpts);
         extractedVelocity=zeros(1,numpts);
         extractedSetpoint=zeros(1,numpts);
         extractedAmplitude=zeros(1,numpts);
         extracteddeltaKappa=zeros(1,numpts);
         extractedts_wsk=zeros(1,numpts);
         contact_sorted_CaSig = zeros(numrois,numframes,numcontacts);
        for j= 1:numcontacts  
 
            timepoint = contactind(j)/wSigframerate;
            strp = contactind(j)-floor(baseline*wSigframerate);
            endp=contactind(j)+ceil(dur*wSigframerate);
            temp=thetavals{i};
            diff=[0,0];
            
            if(strp<1)
                diff(1,1)=strp*-1+1;
                strp=1;
            elseif(endp>length(temp))
                 diff(1,2)=endp-length(temp);
                 endp=length(temp);               
            end
            
%             extractedTheta(1:endp-str)= temp(contactind(j)-ceil(baseline*wSigframerate):contactind(j)+ceil(dur*wSigframerate));
            extractedTheta(diff(1,1)+1:numpts-diff(1,2))= temp(strp+1:endp);
            temp=kappavals{i};
            extractedKappa(diff(1,1)+1:numpts-diff(1,2))= temp(strp+1:endp);

            temp=Velocity{i};
            extractedVelocity(diff(1,1)+1:numpts-diff(1,2))= temp(strp+1:endp);
            temp=Setpoint{i};
            extractedSetpoint(diff(1,1)+1:numpts-diff(1,2))= temp(strp+1:endp);
            temp=Amplitude{i};
            extractedAmplitude(diff(1,1)+1:numpts-diff(1,2))= temp(strp+1:endp);
            temp=deltaKappa{i};
            extracteddeltaKappa(diff(1,1)+1:numpts-diff(1,2))= temp(strp+1:endp);
            temp=ts_wsk{i};
            extractedts_wsk(diff(1,1)+1:numpts-diff(1,2))= temp(strp+1:endp);
            extractedcontactdir=contactdir{i}; 
            extractedcontacts=contacts{i};
            if( isfield(CaTrials, 'ephusTrial') && ~isempty(CaTrials(i).ephusTrial.ephuststep))
                ephussamplerate = 1/CaTrials(i).ephusTrial.ephuststep(1);
                ephusattimepoint = ceil(timepoint*ephussamplerate);
                if( ephusattimepoint< ceil(.5*ephussamplerate))
                    diffephus=ceil(.5*ephussamplerate) - ephusattimepoint;
                    newbaseline=.5-ephusattimepoint/ephussamplerate;
                    ephussamples(1:diffephus)=1;
                    ephussamples(diffephus+1:diffephus+ceil((newbaseline+dur)*ephussamplerate)+1) =(ephusattimepoint-ceil(newbaseline*ephussamplerate)):( ephusattimepoint+ceil(dur*ephussamplerate));
                else
                    ephussamples =(ephusattimepoint-ceil(baseline*ephussamplerate)):( ephusattimepoint+ceil(dur*ephussamplerate));
                end
            end

            if(timepoint<=.5)
                strframe=1;
                diff=ceil((0.5-timepoint)/Caframetime);
            else
                strframe=ceil((timepoint-baseline)/Caframetime);
            end

            endframe = strframe + numframes-1;
            
            if(endframe>size(CaSig,2))
               extractedCaSig(:,1:size(CaSig,2)-strframe+1) = CaSig(:,strframe:size(CaSig,2),i); 
            elseif(timepoint<=baseline)   
                
               extractedCaSig(:,diff+strframe:diff+endframe)= CaSig(:,strframe:endframe,i);  
            else
               extractedCaSig = CaSig(:,strframe:endframe,i);  
            end
            count=count+1;
% % %             contact_CaTrials{count}=struct{'dff',{extractedCaSig},'FrameTime',CaTrials(CaSig_tags(i)).FrameTime,'nFrames',numframes,...
% % %                                        'Trialind', CaTrials(CaSig_tags(i)).TrialNo,'TrialNo',trialnums(i),'nROIs',numrois};
            TrialName = strrep(CaTrials(CaSig_tags(i)).FileName,CaTrials(CaSig_tags(i)).FileName_prefix,'');
            TrialName=strrep(TrialName,'.tif','');
            contact_CaTrials(count).solo_trial = str2num(TrialName);
            contact_CaTrials(count).dff = extractedCaSig;
            contact_CaTrials(count).ts = [strframe:endframe].*Caframetime;
            contact_CaTrials(count).theta = {extractedTheta};
            contact_CaTrials(count).kappa = {extractedKappa};
            contact_CaTrials(count).Velocity = {extractedVelocity};
            contact_CaTrials(count).Setpoint = {extractedSetpoint};
            contact_CaTrials(count).Amplitude = {extractedAmplitude};
            contact_CaTrials(count).deltaKappa = {extracteddeltaKappa};
            contact_CaTrials(count).ts_wsk={extractedts_wsk};
            contact_CaTrials(count).contactdir = {extractedcontactdir};
             contact_CaTrials(count).contacts = {extractedcontacts};
            contact_CaTrials(count).FrameTime = CaTrials(CaSig_tags(i)).FrameTime;
            contact_CaTrials(count).nFrames = numframes;
            contact_CaTrials(count).CaSigTrialind=CaTrials(CaSig_tags(i)).TrialNo;
            contact_CaTrials(count).TrialName= TrialName;
            contact_CaTrials(count).nROIs = numrois;
            contact_CaTrials(count).FileName_prefix=CaTrials(CaSig_tags(i)).FileName_prefix;
%             contact_CaTrials(count).contacts={contactind};
            contact_CaTrials(count).contacts={horzcat(contacttimes{i}{:})};
            contact_CaTrials(count).barpos = wSigTrials{wSig_tags(i)}.bar_pos_trial(1,1);%cellfun(@(x) x.bar_pos_trial(1,1), wSigTrials(wSig_tags),'uniformoutput',false);

            switch(Y(i))
                case(1)
                    contact_CaTrials(count).trialtype = 'Hit';               
                case(2)
                    contact_CaTrials(count).trialtype = 'Miss';
                case(3)
                    contact_CaTrials(count).trialtype = 'CR';                   
                case(4)
                    contact_CaTrials(count).trialtype = 'FA';
            end

              if( isfield(CaTrials, 'ephusTrial')) 
                  if ~isempty(CaTrials(CaSig_tags(i)).ephusTrial.licks)
                     
                    contact_CaTrials(count).licks = CaTrials(CaSig_tags(i)).ephusTrial.licks(ephussamples);
                    contact_CaTrials(count).poleposition = CaTrials(CaSig_tags(i)).ephusTrial.poleposition(ephussamples);
                    contact_CaTrials(count).ephuststep = CaTrials(CaSig_tags(i)).ephusTrial.ephuststep(ephussamples);  
                    
                  else
                    contact_CaTrials(count).licks = zeros(    (baseline +dur)*10000,1);
                    contact_CaTrials(count).poleposition = zeros((baseline +dur)*10000,1);
                    contact_CaTrials(count).ephuststep = zeros((baseline +dur)*10000,1);
                  end          
                  
              end
        end

    end

sorted_CaTrials.contact_trialnums = trialnums;
sorted_CaTrials.contact_trialtypes = Y;


save(['contact_CaTrials' contact_CaSig_tag],'contact_CaTrials');
save('sorted_CaTrials', 'sorted_CaTrials');



% --- Executes on button press in contact_CaSignal_select.
function contact_CaSignal_select_Callback(hObject, eventdata, handles)
set(handles.contact_CaSignal_select,'Value',1);
set(handles.sorted_CaTrials_select,'Value',0);
set(handles.CaTrials_select,'Value',0);

% --- Executes during object creation, after setting all properties.
function contact_CaSignal_select_CreateFcn(hObject, eventdata, handles)


% --- Executes on button press in contactdetect.
function contactdetect_Callback(hObject, eventdata, handles)
    cd(get(handles.baseDataPath,'String'));
	[filename1,path]= uigetfile('wsArray*.mat', 'Load wsArray.mat file');
    cd(path);
    load([path filesep filename1]);
    [filename2,path]= uigetfile('wSigTrials*.mat', 'Load wSigTrials.mat file');
    load([path filesep filename2]);
    contDet_param.threshDistToBarCenter = [.1   .55];
    contDet_param.thresh_deltaKappa = [-.1	.1];
%      contDet_param.bar_time_window = cellfun(@(x) x.bar_time_win, wsArray.ws_trials,'UniformOutput', false);
    barTimeWindow = [1.0 2.5];
    contDet_param.bar_time_window = barTimeWindow;
    contact_inds = cell(wsArray.nTrials,1);
    contact_direct = cell(wsArray.nTrials,1);
   [contact_inds, contact_direct] = Contact_detection_session_auto(wsArray, contDet_param);
    %
    wsArray = NX_WhiskerSignalTrialArray([],wSigTrials);
    
    for i = 1:wsArray.nTrials
        wsArray.ws_trials{i}.contacts = contact_inds{i};
        wsArray.ws_trials{i}.contact_direct = contact_direct{i};
        wsArray.ws_trials{i}.totalTouchKappaTrial = wsArray.totalTouchKappaTrial{1}(i);
        wsArray.ws_trials{i}.maxTouchKappaTrial = wsArray.maxTouchKappaTrial{1}(i);
    end
   
    wSigTrials =wsArray.ws_trials;
    %     save(sprintf('SessionInfo_%s.mat', sessionName{kk}), 'sessionInfo');
    save(filename1, 'wsArray');

    save (filename2,'wSigTrials');
    
    
        %% adding this to sessObj
    sessObj_found = dir('sessObj.mat');
    if isempty(sessObj_found)
        sessObj = {};
        sessObj.wSigTrials = wSigTrials;
        save('sessObj','sessObj','-v7.3');
    else
        load('sessObj.mat');
        sessObj.wSigTrials = wSigTrials;;
        save('sessObj','sessObj','-v7.3');
    end
    
% end


% --- Executes on button press in movefine.
function movefine_Callback(hObject, eventdata, handles)
    


% --- Executes on button press in donottransferROIinfo.
function donottransferROIinfo_Callback(hObject, eventdata, handles)




function trialswindow_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function trialswindow_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function trialswindow_tag_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function trialswindow_tag_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in plot_solodata.
function plot_solodata_Callback(hObject, eventdata, handles)
global solo_data
name = solo_data.sessionName;


[filename,pathName]=uigetfile('wSigTrials*.mat','Want to add contact info? The load wSigTrials.mat file');
if isequal(filename, 0) || isequal(pathName,0)           
    addcontactinfo =0;
else
    addcontactinfo =1;
    load([pathName filesep filename]);
end



figure;

[AX,H1,H2] = plotyy([1:solo_data.trialStartEnd(2)],solo_data.Dprime,[1:solo_data.trialStartEnd(2)],solo_data.PercentCorrect);

hold(AX(1), 'on');
plot(solo_data.polePositions/100000,'b*','MarkerSize',3);

plot(solo_data.hitTrialNums,zeros(1,length(solo_data.hitTrialNums)),'gd','MarkerSize',10,'MarkerFaceColor','g'); 
plot(solo_data.missTrialNums,zeros(1,length(solo_data.missTrialNums)),'rd','MarkerSize',10,'MarkerFaceColor','r'); 
plot(solo_data.falseAlarmTrialNums,ones(1,length(solo_data.falseAlarmTrialNums)),'rd','MarkerSize',10,'MarkerFaceColor','r'); 
plot(solo_data.correctRejectionTrialNums,ones(1,length(solo_data.correctRejectionTrialNums)),'gd','MarkerSize',10,'MarkerFaceColor','g'); 

if addcontactinfo
    wSigfilenames =cellfun(@(x) x.trackerFileName(29:32),wSigTrials,'uniformoutput',false);
    whisker_trials = str2num(char(wSigfilenames));
% %      whisker_trialinds = whisker_trials(find(ismember(whisker_trials,solo_data.trialNums)));
% %     contacttimes=cellfun(@(x) x.contacts{1}, wSigTrials(whisker_trialinds),'uniformoutput',false);
contacttimes=cellfun(@(x) x.contacts{1}, wSigTrials,'uniformoutput',false)   ;
nocontact = cellfun(@isempty,contacttimes);
%     inds = find(contact==1);
%     wSigTrialinds=str2num(wSigfilenames); 
%     plot(wSigTrialinds(inds),contact(inds)*0,'kd','MarkerSize',10);
    plot(whisker_trials,nocontact,'kd','MarkerSize',10);


end

title(['Performance from ' name(length(name)-16:length(name)-8) ' Session' name(length(name)-6:length(name)-1)]);
axis(AX(1),[0 solo_data.trialStartEnd(2) -.5 3.5]);
set(AX(1),'YTick',[-.5 0 0.5 1 1.5 2 2.5 3 3.5]);
hold(AX(1), 'off');

axis(AX(2),[0 solo_data.trialStartEnd(2) -.1 1]);
set(AX(2),'YTick',[-.1:.1:1]);

set(get(AX(1),'Ylabel'),'String','Dprime');set(get(AX(2),'Ylabel'),'String','PercentCorrect');
set(H1,'markersize',5,'Marker','.');set(H2,'markersize',5,'Marker','*') ;

saveas(gcf,'solo_performance_barpos','tif');







function timewindow_wSiganal_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function timewindow_wSiganal_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in plot_fittedwhisker.
function plot_fittedwhisker_Callback(hObject, eventdata, handles)




function numblocks_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function numblocks_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function wh_trialstoavg_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function wh_trialstoavg_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function wSigSum_datapath_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function wSigSum_datapath_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in load_wSigSum.
function load_wSigSum_Callback(hObject, eventdata, handles)

    global wSigSummary
    wSigSummary = {};
    basedatapath = get(handles.wSigSum_datapath,'String');
    if(length(basedatapath)<10)
        basedatapath = '/Volumes/GR_Data_01/Data/';
    end
    cd (basedatapath);
    count=0;

    while(count>=0)
     [filename,pathName]=uigetfile('wSigBlocks*.mat','Load wSigBlocks*.mat file');
        if isequal(filename, 0) || isequal(pathName,0)           
           break
        end
       count=count+1;
      load( [pathName filesep filename], '-mat'); 
       set(handles.wSigSum_datapath,'String',pathName);
       cd (pathName);
      wSigSummary{count} = blocks;
 
    end
    folder = uigetdir;
    cd (folder);
    save('wSigSummary','wSigSummary');


% --- Executes on selection change in wSigSum_toplot.
function wSigSum_toplot_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function wSigSum_toplot_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in plot_wSigSum.
function plot_wSigSum_Callback(hObject, eventdata, handles)
global wSigSummary
global wSigSum_Sessions

bartheta= str2num(get(handles.current_bartheta,'String'));
baseline_bartheta = str2num(get(handles.unbiased_bartheta,'String'));
plotlist = get(handles.wSigSum_toplot,'String');
datatoplot= plotlist{get(handles.wSigSum_toplot,'Value')};

blocks= get(handles.wSigSum_block,'String');
blocklist = blocks(get(handles.wSigSum_block,'Value'));
block= get(handles.wSigSum_block,'Value');
numblocks = length(blocklist);

numsessions=size(wSigSummary,2);
commentstr = get(handles.plot_wSigSum_title,'String');
col = [.75 .75 .75; 1 .6 0];
lwidth = linspace(0 ,1,numsessions);
transparency =  0.5;
legendstr = cell(numsessions,1);
% datacollected = cell(numblocks,1);
datacollected = zeros(numsessions*70,4,numblocks);

spacing = 1/(numsessions+3.5);
baseline_sessions =2;
% plotting setpoint
for j= 1:numblocks
    block =j;
    sc = get(0,'ScreenSize');
    h_fig1 = figure('position', [1000, sc(4)/10-100, sc(3)*3/10, sc(4)*3/4], 'color','w'); %%raw setpoint
    ah1=axes('Parent',h_fig1); title([commentstr 'Raw Data ' ]);%blocklist{j} 'Data ' datatoplot]);
    
    h_fig2 = figure('position', [300, sc(4)/10-100, sc(3)*3/10, sc(4)*3/4], 'color','w'); %% error from bartheat
    ah2=axes('Parent',h_fig2); title([commentstr 'Data Error  ' ]);%blocklist{j} 'Data ' datatoplot]);
    hold on;
%     figure;
    count =0;
    prev=0;
   
    for i = 1:numsessions
% % %        temp = getfield(wSigSummary{i},datatoplot); 
% % %        xdata = temp{block}(:,1)';
% % %        ydata = temp{block}(:,2)';
      
       
%        binneddata = strrep(datatoplot,'med','medbinned');
       temp = getfield(wSigSummary{i},datatoplot); 
       binnedxdata = temp{block}(:,1)';
       binnedydata = temp{block}(:,mod(get(handles.wSigSum_toplot,'Value'),4)+1)';
       
% % %        minmaxdata = strrep(datatoplot,'med','width');
% % %        temp = getfield(wSigSummary{i},minmaxdata); 
       

     

      axes(ah1);  
%       jbfill([count+1:count+length(data)],upper,lower,col(j,:),col(j,:),1,transparency); hold on;
%       plot([count+1:count+length(data)],data,'color',col(j,:),'linewidth',1.5);
      
%         plot(xdata+count,ydata,'color',col(j,:),'linewidth',1.0); hold on;
        if(i<baseline_sessions+1)
            plot(binnedxdata+count,binnedydata,'color','k','Marker','o','MarkerSize',6,'MarkerFaceColor','k');
            hline(baseline_bartheta,'k--');
        else
            plot(binnedxdata+count,binnedydata,'color','k','Marker','o','MarkerSize',6,'MarkerFaceColor','r');
            hline(bartheta,'r--');
        end
       legendstr(i) = {['session' num2str(i) ' ']};
% % %        wSigSum_Sessions.setpoint{i,j}=[data;lower;upper;xval']';
%         wSigSum_Sessions.setpoint{i,j}=[xdata;ydata]';
       wSigSum_Sessions.databinned{i,j} = [binnedxdata;binnedydata]';
       axes(ah2);
% % %        jbfill([count+1:count+length(data)],upper-bartheta,lower-bartheta,col(j,:),col(j,:),1,transparency); hold on;
% % %        plot([count+1:count+length(data)],data-bartheta,'color',col(j,:),'linewidth',1.5);     
    if(i<baseline_sessions+1)
%        plot(xdata+count,(ydata-baseline_bartheta),'color',col(j,:),'linewidth',1.0);  
       plot(binnedxdata+count,(binnedydata-baseline_bartheta),'color','k','Marker','o','MarkerSize',6,'MarkerFaceColor','k');
    else
        
%        plot(xdata+count,(ydata-bartheta),'color',col(j,:),'linewidth',1.0);  
       plot(binnedxdata+count,(binnedydata-bartheta),'color','k','Marker','o','MarkerSize',6,'MarkerFaceColor','r');
    end   
       legendstr(i) = {['session' num2str(i) ' ']};

          count = count+xdata(end)+10;

    end
    axes(ah1); axis([0 count -20 20]);grid on;
   
   saveas(gcf,['allsessions_theta_bl ' blocklist{j}] ,'tif');  
   saveas(gcf,['allsessions_theta_bl' blocklist{j}],'fig');
   
   axes(ah2);axis([0 count -20 20]);grid on;
   saveas(gcf,['allsessions_error_bl ' blocklist{j}] ,'tif');  
   saveas(gcf,['allsessions_error_bl' blocklist{j}],'fig');
end
% title([commentstr ' Block ' blocklist{j} ]);
 hold off;

% % %  
% % %  % plotting amp
% % % for j= 1:numblocks
% % %     block =j;
% % %     sc = get(0,'ScreenSize');
% % %     figure('position', [1000, sc(4)/10-100, sc(3)*3/10, sc(4)*3/4], 'color','w'); %%raw theta
% % % %     title([commentstr 'Amplitude Block ' blocklist{j} 'Data ' 'Amp_med']);
% % %     title([commentstr 'Whisker Amplitude ' ]);%blocklist{j} 'Data ' datatoplot]);
% % %     hold on;
% % %     count =0;
% % %     prev=0;   
% % %     for i = 1:numsessions
% % %         ampdata = strrep(datatoplot,'setpoint','amp');
% % %        temp = getfield(wSigSummary{i},ampdata); 
% % %        ydata = temp{block}(:,2)';
% % %        xdata = temp{block}(:,1)';
% % %        binneddata = strrep(ampdata,'med','medbinned');
% % %         temp = getfield(wSigSummary{i},binneddata); 
% % %        binnedydata = temp{block}(:,2)';
% % %        binnedxdata = temp{block}(:,1)';      
% % % % % %        xdata = strrep(datatoplot,'_setpoint_trials_med','');
% % % % % %        xdata = ([xdata 'trialnums']);
% % % % % %        temp = getfield(wSigSummary{i},xdata); 
% % % % % %        xval= temp{block};
% % %        
% % %        %binning over 25 trials
% % % % % %        windowSize=25;
% % % % % %        data=filter(ones(1,windowSize)/windowSize,1,data);
% % % % % %        plottemp=data;
% % % % % %         plottemp(1:25)=plottemp(26);
% % % % % %          plottemp(2:24)=nan;
% % %  
% % % % % %       plot([count+1:count+length(plottemp)],plottemp,'color',col(j,:),'linewidth',1.5);
% % % 
% % %        plot(xdata+count,ydata,'color',col(j,:),'linewidth',1.0);hold on;
% % %      if(i<baseline_sessions+1)
% % %          plot(binnedxdata+count,binnedydata,'color','k','Marker','o','MarkerSize',6,'MarkerFaceColor','k');
% % %      else
% % %          plot(binnedxdata+count,binnedydata,'color','k','Marker','o','MarkerSize',6,'MarkerFaceColor','r');
% % %      end
% % %        
% % %         legendstr(i) = {['session' num2str(i) ' ']};
% % %        wSigSum_Sessions.amp{i,j}=[xdata;ydata]';
% % %        wSigSum_Sessions.ampbinned{i,j}= [binnedxdata;binnedydata]';
% % %   
% % %           count = count+xdata(end)+10;
% % % 
% % %     end
% % %     axis([0 count -1 30]);grid on;
% % %    saveas(gcf,['allsessions_amplitude_bl ' blocklist{j}] ,'tif');  
% % %    saveas(gcf,['allsessions_amplitude_bl' blocklist{j}],'fig');
% % % end
% % % % title([commentstr ' Block ' blocklist{j} ]);
% % %  hold off;
% % %  
% % %  
% % %   % plotting prewhisk_setpoint
% % % for j= 1:numblocks
% % %     block =j;
% % %     sc = get(0,'ScreenSize');
% % %   figure('position', [1000, sc(4)/10-100, sc(3)*3/10, sc(4)*3/4], 'color','w'); 
% % %     title([commentstr 'Prewhisk_Setpoinbt' ]); %' blocklist{j}]);
% % % 
% % %     hold on;
% % %     count =0;
% % %     prev=0;   
% % %     for i = 1:numsessions
% % %         prewhiskdata = strrep(datatoplot,'med','prewhisk');
% % %        temp = getfield(wSigSummary{i},prewhiskdata); 
% % %        ydata = temp{block}(:,2)';
% % %        xdata = temp{block}(:,1)';
% % %        binneddata = strrep(prewhiskdata,'prewhisk','prewhiskbinned');
% % %         temp = getfield(wSigSummary{i},binneddata); 
% % %        binnedydata = temp{block}(:,2)';
% % %        binnedxdata = temp{block}(:,1)';      
% % %  
% % %        plot(xdata+count,ydata,'color',col(j,:),'linewidth',1);hold on;
% % %       if(i<baseline_sessions+1)
% % %          plot(binnedxdata+count,binnedydata,'color','k','Marker','o','MarkerSize',6,'MarkerFaceColor','k');
% % %           hline(baseline_bartheta,'k--');
% % %       else
% % %           plot(binnedxdata+count,binnedydata,'color','k','Marker','o','MarkerSize',6,'MarkerFaceColor','r');
% % %           hline(bartheta,'r--');
% % %       end
% % %       
% % %     
% % %       
% % %         legendstr(i) = {['session' num2str(i) ' ']};
% % %        wSigSum_Sessions.setpoint_prewhisk{i,j}=[xdata;ydata]';
% % %        wSigSum_Sessions.setpoint_prewhiskbinned{i,j}= [binnedxdata;binnedydata]';
% % %   
% % %           count = count+xdata(end)+10;
% % % 
% % %     end
% % %     axis([0 count -20 20]);grid on;
% % %    saveas(gcf,['allsessions_prewhisk_setpoint_bl ' blocklist{j}] ,'tif');  
% % %    saveas(gcf,['allsessions_prewhisk_setpoint_bl' blocklist{j}],'fig');
% % % end
% % % % title([commentstr ' Block ' blocklist{j} ]);
% % %  hold off;
% % %  
 
   folder=uigetdir;
   cd(folder);
save('wSigSum_Sessions','wSigSum_Sessions');  
cd ..


%  vline(lasttrial(:,1),'k:');
% legend('Initial','Initial','Biased','B+-');
% title([commentstr  'All sessions ' datatoplot]);
 
% --- Executes on selection change in wSigSum_block.
function wSigSum_block_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function wSigSum_block_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function plot_wSigSum_title_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function plot_wSigSum_title_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function timewindowtag_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function timewindowtag_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in sortwtouchInfo.
function sortwtouchInfo_Callback(hObject, eventdata, handles)




function ephusDataPath_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function ephusDataPath_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in make_ephus_obj.
function make_ephus_obj_Callback(hObject, eventdata, handles)

ephuspath= get(handles.ephusDataPath,'String');
[fname,ephuspath] = uigetfile(ephuspath,'Pick ephus dir');
  if isequal(ephuspath, 0) 
      return
  else
      set(handles.ephusDataPath,'String',ephuspath);
      cd(ephuspath);
  end
mouseName = get(handles.AnimalNameEdit, 'String');
sessionID = get(handles.SoloSessionID, 'String');
if((isempty(mouseName))||(isempty(sessionID)))
    msgbox('Fill in details');
    return;
end
% trialStartEnd = [get(handles.SoloStartTrialNo,'value'), get(handles.SoloEndTrialNo,'value')];
% trialStartEnd = [handles.SoloStartTrial,handles.SoloEndTrial];
% mouseName=
obj = ephusTrialArray_gr(mouseName, sessionID,ephuspath);
cd ..
save(['ephusdata_' mouseName '_' sessionID],'obj');



% --- Executes on button press in addephusTrials.
function addephusTrials_Callback(hObject, eventdata, handles)
global CaSignal %  ROIinfo ICA_ROIs
ephuspath = get(handles.ephusDataPath,'String');

[fname,pathName]=uigetfile(ephuspath,'Load the ehus obj');
if(isequal(pathName, 0)) 
else
    load([pathName filesep fname],'-mat');
end
ephustrials = size(obj,2);
trialStartEnd(1) = str2num(obj(1).trialname);
trialStartEnd(2) = str2num(obj(ephustrials).trialname);

[fname,imaging_datapath] = uigetfile(ephuspath,'Setimaging data path');
    if isequal(imaging_datapath, 0) 
    else
        cd (imaging_datapath);
        list = dir('Image*.tif');
        filenames =cell2mat(arrayfun(@(x) x.name(length(x.name)-6 :length(x.name)-4),list,'uniformoutput',false));
        trials=str2num(filenames);
        filenames=arrayfun(@(x) x.trialname,obj,'uniformoutput',false);
        ephustrials=str2num(cell2mat(filenames'));  
        commontrials=ephustrials(ismember(ephustrials,trials));
        commontrialinds = find(ismember(ephustrials,trials));
    end


if length(commontrials) ~= str2num(get(handles.TotTrialNum, 'String'))
    error('Number of ephus trials NOT equal to Number of Ca Image Trials!')
end

for i = 1:length(commontrials)
    
    CaSignal.CaTrials(i).ephusTrial = obj(commontrialinds(i));
end
disp([num2str(i) ' Ephus Trials added to CaSignal.CaTrials']);
set(handles.msgBox, 'String', [num2str(i) ' Ephus Trials added to CaSignal.CaTrials']);

%% adding ephus to sessObj
cd (pathName);
current_dir = pwd;
separators = find(current_dir == filesep);
session_dir = current_dir(1:separators(length(separators)-0)); % one folder up from ephus_data dir
cd (session_dir);
sessObj_found = dir('sessObj.mat');
if isempty(sessObj_found)
    sessObj = {};
    sessObj.ephusTrials = obj;
    save('sessObj','sessObj','-v7.3');
else
    load('sessObj.mat');
    sessObj.ephusTrials = obj;
    save('sessObj','sessObj','-v7.3');
end
cd (current_dir);

guidata(hObject, handles)




function ephusTrialsToBeExcluded_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function ephusTrialsToBeExcluded_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in load_wSigSummary.
function load_wSigSummary_Callback(hObject, eventdata, handles)
global wSigSummary
wSigSummary ={};
[fname,datapath] = uigetfile(pwd,'Select wSigSummary data');
cd (datapath);
load(fname,'-mat');


% --- Executes on button press in pooleddata.
function pooleddata_Callback(hObject, eventdata, handles)
    global wSigSum_anm
    wSigSum_anm = {};
        basedatapath = '/Volumes/GR_Data_03/Data/';

    cd (basedatapath);
    count=0;

    while(count>=0)
     [filename,pathName]=uigetfile('wSigSum_Sessions*.mat','Load wSigSum_Sessions*.mat file');
        if isequal(filename, 0) || isequal(pathName,0)           
           break
        end
       count=count+1;
      load( [pathName filesep filename], '-mat'); 
       cd (pathName);
      wSigSum_anm{count} = wSigSum_Sessions;
 
    end
    folder = uigetdir;
    cd (folder);
    save('wSigSum_anm','wSigSum_anm');




% --- Executes on button press in plot_pooledbatchdata.
function plot_pooledbatchdata_Callback(hObject, eventdata, handles)

[filename,pathName]=uigetfile('wSigSum_anm.mat','Load wSigSum_anm.mat file');
load( [pathName filesep filename], '-mat'); 
cd (pathName);
numanm = size(wSigSum_anm,2);
sess_count=nan(numanm,2);
maxtrialcount=0;
for i = 1:numanm
     tempobj=wSigSum_anm{i}.setpoint;
     numsessions=size(tempobj,1);
     sess_count(i,1) =numsessions;
     
%      new_wSigSum_anm={numsessions,1};
       %combine blocks

         for j= 1:numsessions
            if(size(tempobj,2)>1)
             new_wSigSum_anm {i}{j} =  [tempobj{j,1} ; tempobj{j,2}];     
            else
             new_wSigSum_anm{i}{j}= tempobj{j,1};
           end
         end 

end
sess_trial_count=zeros(max(sess_count(:,1)),numanm);
for j=1:max(sess_count(:,1))

    for i=1:numanm
         if(j<=sess_count(i,1))
         trialcount =   size(new_wSigSum_anm{i}{j},1);
         sess_trial_count(j,i) = trialcount ;
         maxtrialcount = (maxtrialcount<=trialcount)*trialcount + (maxtrialcount>trialcount)*maxtrialcount;
         end
    end
    temp1=nan(maxtrialcount,numanm);
    for i=1:numanm
         if(j<=sess_count(i,1))
            temp2=new_wSigSum_anm{i}{j}(:,1);       
            temp1(1:length(temp2),i)= temp2;        
         end
    end
    sess{j} = temp1;
     t2=temp1(1:30,:);
     
%     baselinepersesion(j,:) = min(t2,[],1);
    baselinepersesion(j,:) = mean(t2,1);
end   

%% calculate and plot dsetpoint from bar theta

%temp override with mean bar theta
%% miceorder of animls 183484,181851,181056,184169,188197
  baseline = [ 12.78, 5.6,9.12,8.68,7.17]; 
%  baseline = [  5.6,9.12,8.68,7.17]; 

temp=0;temp3=0;temp4=0;
for j=1:max(sess_count(:,1))
   temp=sess{j} ;
   temp2 = repmat(baseline(1,:),size(temp,1),1);
   temp= (temp-temp2);%./temp2;
   sess1{j}=temp;
end

figure;temp=0;temp3=0;temp4=0;count =0;
for j=1:max(sess_count(:,1))
    mintrialcount= min(nonzeros(sess_trial_count(j,:)));
    mean_dSetpoint=zeros(mintrialcount,1);
    sem_dSetpoint=zeros(mintrialcount,1);
    for i = 1:numanm
        if(j<=sess_count(i,1))
            temp=sess1{j}(:,i);
            temp=temp(find(~isnan(temp)));
            dSetpoint{j}(1:mintrialcount,i) = temp(1:mintrialcount);%temp(length(temp)-mintrialcount+1:length(temp));
        else
            dSetpoint{j}(1:mintrialcount,i) = nan;
        end
    end
    temp=dSetpoint{j};
    mean_dSetpoint=nanmean(temp,2);
    
    Xsq=(~isnan(temp).*temp.^2);
    exp_Xsq=nanmean(Xsq,2);
    expX_sq=mean_dSetpoint.*mean_dSetpoint;
    expX_sq=mean_dSetpoint.^2;
    sem_dSetpoint=sqrt(exp_Xsq-expX_sq)/sqrt(numanm+2);
    upper = mean_dSetpoint+ sem_dSetpoint;
    lower = mean_dSetpoint - sem_dSetpoint ;
    upper(isnan(upper))=0;
    lower(isnan(lower))=0;
    jbfill([count+1: count+length(sem_dSetpoint)],upper' ,lower' ,[.5 .5 .5 ],[.5 .5 .5 ],1,.5);hold on;  
     plot([count+1: count+length(sem_dSetpoint)],mean_dSetpoint,'color',[0 0 0 ],'linewidth',2);hold off;title('dSetpoint_frombartheta');hold on;
    count = count+length(sem_dSetpoint)+3;
end
 saveas(gcf,'dSetpoint_frombartheta','fig');
 save('dSetpoint_bartheta','dSetpoint');
 
temp=0;temp3=0;temp4=0;
for j=1:max(sess_count(:,1))
   temp=dSetpoint{j};
   for i=1:numanm
        if(j<=sess_count(i,1))
           temp3 = temp(find(~isnan(temp(:,i))),i);
           segment =floor(length(temp3)/3);
           parts = [1,segment;segment+1,2*segment;2*segment+1,length(temp3)];
           for k= 1:length(parts)
              temp4(k,i) = mean(temp3(parts(k,1):parts(k,2)));
           end
        else
            temp4(:,i) = nan;
        end

   end
   binnedsess1{j}=temp4;
end
 
count =0; 
figure;
for j=1:max(sess_count(:,1))
    binnedsess(:,:,j)=binnedsess1{j};
    binned_dSetpoint=binnedsess1{j};
    mean_binned_dSetpoint=nanmean(binned_dSetpoint,2);
    exp_Xsq=nanmean((~isnan(binned_dSetpoint).*binned_dSetpoint.^2),2);

    expX_sq=mean_binned_dSetpoint.^2;
    sem_binned_dSetpoint=sqrt(exp_Xsq-expX_sq)/sqrt(numanm+2);
    L=mean_binned_dSetpoint-sem_binned_dSetpoint;
    U=mean_binned_dSetpoint+sem_binned_dSetpoint;
    xpnts = [count+1 : count+3];
%     errorbar(xpnts,mean_binned_dSetpoint,L,U,'color',[1 .6 0 ],'o');
    errorbar(xpnts,mean_binned_dSetpoint,sem_binned_dSetpoint,'color',[.5 .5 .5 ],'LineWidth',2,'Marker','o','MarkerSize',20,'MarkerFaceColor',[0 0 0 ]);
    hold on;title('binned_dSetpoint_frombartheta');
    count = count + 3;
end
 saveas(gcf,'dSetpoint_frombartheta_binned','fig');
 save('binned_dSetpoint_bartheta','binnedsess1');
 
 
 %% calculate and plot dsetpoint for within each session for each animal
 temp=0;temp3=0;temp4=0;
for j=1:max(sess_count(:,1))
   temp=sess{j} ;
   temp2 = repmat(baselinepersesion(j,:),size(temp,1),1);
   temp= (temp-temp2);%./temp2;
   sess1{j}=temp;
% % % %    for i=1:numanm
% % % %         if(j<=sess_count(i,1))
% % % %            temp3 = temp(find(~isnan(temp(:,i))),i);
% % % %            segment =floor(length(temp3)/3);
% % % %            parts = [1,segment;segment+1,2*segment;2*segment+1,length(temp3)];
% % % %            for k= 1:length(parts)
% % % %               temp4(k,i) = mean(temp3(parts(k,1):parts(k,2)));
% % % %            end
% % % %         else
% % % %             temp4(:,i) = nan;
% % % %         end
% % % % 
% % % %    end
% % % %    binnedsess1{j}=temp4;
end

figure;temp=0;temp3=0;temp4=0;count =0;
for j=1:max(sess_count(:,1))
    mintrialcount= min(nonzeros(sess_trial_count(j,:)));
    mean_dSetpoint=zeros(mintrialcount,1);
    sem_dSetpoint=zeros(mintrialcount,1);
    for i = 1:numanm
        if(j<=sess_count(i,1))
            temp=sess1{j}(:,i);
            temp=temp(find(~isnan(temp)));
            dSetpoint{j}(1:mintrialcount,i) = temp(1:mintrialcount);%(length(temp)-mintrialcount+1:length(temp));
        else
            dSetpoint{j}(1:mintrialcount,i) = nan;
        end
    end
    temp=dSetpoint{j};
    mean_dSetpoint=nanmean(temp,2);
    
    Xsq=(~isnan(temp).*temp.^2);
    exp_Xsq=nanmean(Xsq,2);
%     expX_sq=mean_dSetpoint.*mean_dSetpoint;
    expX_sq=mean_dSetpoint.^2;
    sem_dSetpoint=sqrt(exp_Xsq-expX_sq)/sqrt(numanm+2);
    upper = mean_dSetpoint+ sem_dSetpoint;
    lower = mean_dSetpoint - sem_dSetpoint ;
    upper(isnan(upper))=0;
    lower(isnan(lower))=0;
    jbfill([count+1: count+length(sem_dSetpoint)],upper' ,lower' ,[.5 .5 .5 ],[.5 .5 .5 ],1,.5);hold on;  
     plot([count+1: count+length(sem_dSetpoint)],mean_dSetpoint,'color',[.0 0 0  ],'linewidth',2);hold off;title('dSetpoint_withinsession');hold on;
    count = count+length(sem_dSetpoint)+3;
end
 saveas(gcf,'dSetpoint_withinsession','fig');
 save('dSetpoint_wsession','dSetpoint');
 
 temp=0;temp3=0;temp4=0;
for j=1:max(sess_count(:,1))
    temp=dSetpoint{j};
   for i=1:numanm
        if(j<=sess_count(i,1))
           temp3 = temp(find(~isnan(temp(:,i))),i);
           segment =floor(length(temp3)/3);
           parts = [1,segment;segment+1,2*segment;2*segment+1,length(temp3)];
           for k= 1:length(parts)
              temp4(k,i) = mean(temp3(parts(k,1):parts(k,2)));
           end
        else
            temp4(:,i) = nan;
        end

   end
   binnedsess1{j}=temp4;
end 
 
count =0; 
figure;
for j=1:max(sess_count(:,1))
    binnedsess(:,:,j)=binnedsess1{j};
    binned_dSetpoint=binnedsess1{j};
    mean_binned_dSetpoint=nanmean(binned_dSetpoint,2);
    exp_Xsq=nanmean((~isnan(binned_dSetpoint).*binned_dSetpoint.^2),2);

    expX_sq=mean_binned_dSetpoint.^2;
    sem_binned_dSetpoint=sqrt(exp_Xsq-expX_sq)/sqrt(numanm+2);
    L=mean_binned_dSetpoint-sem_binned_dSetpoint;
    U=mean_binned_dSetpoint+sem_binned_dSetpoint;
    xpnts = [count+1 : count+3];
%     errorbar(xpnts,mean_binned_dSetpoint,L,U,'color',[1 .6 0 ],'o');
    errorbar(xpnts,mean_binned_dSetpoint,sem_binned_dSetpoint,'color',[.5 .5 .5],'LineWidth',2,'Marker','o','MarkerSize',20,'MarkerFaceColor',[0 0 0 ]);
    hold on;title('binned_dSetpoint_withinsession');
    count = count + 3;
end
 saveas(gcf,'dSetpoint_withinsession_binned','fig');
 save('binned_dSetpoint_wsession','binnedsess1');
 
 
 %%%
  %% calculate and plot amprange accross sessions  for each animal
for i = 1:numanm
     tempobj=wSigSum_anm{i}.amp;
     numsessions=size(tempobj,1);
     sess_count(i,1) =numsessions;
         for j= 1:numsessions
            if(size(tempobj,2)>1)
             new_wSigSum_anm {i}{j} =  [tempobj{j,1} ; tempobj{j,2}];     
            else
             new_wSigSum_anm{i}{j}= tempobj{j,1};
           end
         end 
end

for j=1:max(sess_count(:,1))
    temp1=nan(maxtrialcount,numanm);
    for i=1:numanm
         if(j<=sess_count(i,1))
            temp2=new_wSigSum_anm{i}{j}(:,1);       
            temp1(1:length(temp2),i)= temp2;        
         end
    end
    sess{j} = temp1;
     t2=temp1(1:30,:);
    baselinepersesion(j,:) = mean(t2,1);
end   

%% calculate and plot amprange

temp=0;temp3=0;temp4=0;
for j=1:max(sess_count(:,1))
   temp=sess{j} ;
   sess1{j}=temp;
end
figure;temp=0;temp3=0;temp4=0;count =0;
for j=1:max(sess_count(:,1))
    mintrialcount= min(nonzeros(sess_trial_count(j,:)));
    mean_damp=zeros(mintrialcount,1);
    sem_damp=zeros(mintrialcount,1);
    for i = 1:numanm
        if(j<=sess_count(i,1))
            temp=sess1{j}(:,i);
             temp=temp(find(~isnan(temp)));
            damp{j}(1:mintrialcount,i) = temp(length(temp)-mintrialcount+1:length(temp));
        else
            damp{j}(1:mintrialcount,i) = nan;
        end
    end
    temp=damp{j};
    mean_damp=nanmean(temp,2);
    
    Xsq=(~isnan(temp).*temp.^2);
    exp_Xsq=nanmean(Xsq,2);
    expX_sq=mean_damp.*mean_damp;
    expX_sq=mean_damp.^2;
    sem_damp=sqrt(exp_Xsq-expX_sq)/sqrt(numanm+2);
    upper = mean_damp+ sem_damp;
    lower = mean_damp - sem_damp ;
    ind= ~isnan(mean_damp);
    
    mean_damp=mean_damp(ind);
    sem_damp=sem_damp(ind);
    upper=upper(ind);
    lower=lower(ind);
    
%     upper(isnan(upper))=0;
%     lower(isnan(lower))=0;
    mean_damp(end-5:end)=[];
    sem_damp(end-5:end)=[];
    upper(end-5:end)=[];
    lower(end-5:end)=[];
    jbfill([count+1: count+length(mean_damp)],upper' ,lower' ,[.5 .5 .5 ],[.5 .5 .5 ],1,.5);hold on;  
     plot([count+1: count+length(mean_damp)],mean_damp,'color',[0 0 0 ],'linewidth',2);hold off;title('damp');hold on;
    count = count+length(mean_damp)+3;
end
 saveas(gcf,'damp','fig');
 save('damp','damp');
 
temp=0;temp3=0;temp4=0;
for j=1:max(sess_count(:,1))
   temp=damp{j};
   for i=1:numanm
        if(j<=sess_count(i,1))
           temp3 = temp(find(~isnan(temp(:,i))),i);
           segment =floor(length(temp3)/3);
           parts = [1,segment;segment+1,2*segment;2*segment+1,length(temp3)];
           for k= 1:length(parts)
              temp4(k,i) = mean(temp3(parts(k,1):parts(k,2)));
           end
        else
            temp4(:,i) = nan;
        end

   end
   binnedsess1{j}=temp4;
end
 
count =0; 
figure;
for j=1:max(sess_count(:,1))
    binnedsess(:,:,j)=binnedsess1{j};
    binned_damp=binnedsess1{j};
    mean_binned_damp=nanmean(binned_damp,2);
    exp_Xsq=nanmean((~isnan(binned_damp).*binned_damp.^2),2);
    expX_sq=mean_binned_damp.^2;
    sem_binned_damp=sqrt(exp_Xsq-expX_sq)/sqrt(numanm+2);
    xpnts = [count+1 : count+3];
    errorbar(xpnts,mean_binned_damp,sem_binned_damp,'color',[.5 .5 .5 ],'LineWidth',2,'Marker','o','MarkerSize',20,'MarkerFaceColor',[0 0 0 ]);
    hold on;title('binned_damp');
    count = count + 3;
end
 saveas(gcf,'damp_binned','fig');
 save('binned_damp','binnedsess1');

  
 
 
function current_bartheta_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function current_bartheta_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in align_to_first_touch.
function align_to_first_touch_Callback(hObject, eventdata, handles)
handles.aligned_contact = get(handles.align_to_first_touch,'Value'); 
if handles.aligned_contact
    set(handles.align_to_last_touch,'Value',0);
end

% --- Executes on button press in align_to_last_touch.
function align_to_last_touch_Callback(hObject, eventdata, handles)
handles.aligned_contact = get(handles.align_to_last_touch,'Value'); 
if handles.aligned_contact
    set(handles.align_to_first_touch,'Value',0);
end



function unbiased_bartheta_Callback(hObject, eventdata, handles)
% hObject    handle to unbiased_bartheta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of unbiased_bartheta as text
%        str2double(get(hObject,'String')) returns contents of unbiased_bartheta as a double


% --- Executes during object creation, after setting all properties.
function unbiased_bartheta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to unbiased_bartheta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
