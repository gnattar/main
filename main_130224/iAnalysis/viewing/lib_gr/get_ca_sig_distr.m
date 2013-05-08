function out = get_ca_sig_distr(imTrialArray, varargin)
%
%
% SigProp, signal properties, 'peak','fwhm','riseTime'('rt'), '
% varargin{}, 1, ROInums
%                 2, prop_name, 'peak', 'eventPeak', 'area', 'fwhm', 'decay', 'rise', 
%                 3, trialPhase, 'stim', 'reward', 'trial'.     Default, 'stim'
%                 4, trialType, 'hits', 'miss', 'cr', 'fa'.         Defualt, 'hits'
%                 5, ROIType, 'spine', 'branch', 'trunk', 'all'.    Default, 'all'
%

% if ~iscell(imTrialArray)
%     imTrialArray = {imTrialArray};
% end
% for k = 1:length(imTrialArray)
%     events_prop = imTrialArray{k}.ROI_events_param;
% end

if ~isempty(varargin{1})
    ROInums = varargin{1};
end
if length(varargin)<2
    prop_name = 'peak';
    trialPhase = 'stim';
    trialType = 'hit';
    ROIType = 'all';
else
    prop_name = varargin{2};
end
if length(varargin)<3
     trialPhase = 'stim';
     trialType = 'hit';
     ROIType = 'all';
else
    trialPhase = varargin{3};
end
if length(varargin)<4
     trialType = 'hit';
     ROIType = 'all';
else
    trialType = varargin{4};
end
if length(varargin)<5
    ROIType = 'all';
else
    ROIType = varargin{5};
end

k = find(strcmpi(trialPhase, {'pre_stim', 'stim', 'reward', 'trial'}));
event_struct = imTrialArray.ROI_events_param(k);
IndsStruct = imTrialArray.trialInds_sorted_by_behav;
if strcmpi(trialType,'all')
    trialInds = 1:imTrialArray.nTrials;
else
    trialInds = IndsStruct.(trialType);
end
switch prop_name
    case 'peak'
        out = max();
    case 'eventPeak'
        out = event_struct.peaks(trialInds, ROInums);
    case 'area'
        out = event_struct.areas(trialInds, ROInums);
    case 'fwhm'
        out = event_struct.fwhm(trialsInds, ROInums);
    case 'decay'
        out = event_struct.tauDecay(trialInds, ROInums);
%     case 'rise'
end

        
        
        