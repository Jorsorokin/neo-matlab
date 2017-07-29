function [rH] = PlotRasters(SpikeTimes,start,varargin)
% function [rasterhandle] = PlotRasters(SpikeTimes,start,varargin)
%
%   Plots rasters of spikes based on event times denoted as "Start".
%   Start is an nx1 matrix of start times (or scalar if start time is equal for all sweeps).
%   SpikeTimes can be a vector or matrix, but you must specify this. 
%   Output is a handle to the plot.
%
%               >>> INPUTS >>>
% Required:
%   SpikeTimes = matrix or cell of spiketimes in SECONDS
%           If matrix, will assume columns are trials, rows are spiketimes.
%           If cell, will assume each cell is a channel, and in each cell
%           columns are trials, rows are spiketimes.
%   start = n-element vector containing times of events (in seconds)...if
%       spiketimes are relative (i.e. blocks of spikeitmes, each block relative
%       to stim onset), then user should define "start" as the time of the
%       onset of the stim relative to the start of the block. For instance, if
%       stim is 6 seconds into the start of eah block, then set "start" = 6.
%       * if spikes NOT relative, start = nx1 array of starting times.
% Optional:
%   pre_time = time (in SECONDS) to subtract from starting time...
%           makes plots relative to stim onset (default = 1s);
%   post_time = time (in SECONDS) to add to starting time...
%           (default = 1s);
%   name = name to save figure (default = "raster.pdf")
%   saving = 0 or 1 (default = 0). If 1, saves figure to current directory.
%   
%               <<< OUTPUTS <<<
%   rH = handle to figure
%
% Example:
%   SpikeTimes{1} = sort(randn(100,100));
%   SoujeTunes{2} = sort(randn(100,100));
%   PlotRasters(SpikeTimes,.5,.2,0.5,'control',1)
%       % plots a raster of SpikeTimes for each channel
%           specifying 0.5s as the relative event onset 
%           (0 on x-axis), .2s pre and .5s post this start time.
%           Saves the rasters and appends "control" to the filename
%
%   
% By JMS, 11/13/2015
%-------------------------------------------------------

% check optional
if nargin>2 && ~isempty(varargin{1})
    pre_time = varargin{1};
else pre_time = 1; end
if nargin>3 && ~isempty(varargin{2})
    post_time = varargin{2};
else post_time = 1; end
if nargin>4 && ~isempty(varargin{3})
    name = varargin{3};
end
if nargin>5 && ~isempty(varargin{4})
    saving = varargin{4};
else saving = 0; end

% check if SpikeTimes is cell or matrix 
if iscell(SpikeTimes)
    cell_array = 1; 
    ntrials = size(SpikeTimes{1},2); % number of trials
    chans = numel(SpikeTimes); % number of channels
else 
    cell_array = 0;
    chans = 1; 
    ntrials = size(SpikeTimes,2); % columns of spike matrix
end

% error check if size(start) ~= ntrials....if so, repeat "start" so that it
% has same dimension as ntrials
if numel(start) == 1 && ntrials > 1
    start = ones(ntrials,1)*start;
elseif numel(start) ~= ntrials
    assert('Error: starting times and # trials not same dimension');
end

%set plot parameters
ylim = [0 ntrials+1]; % for plotting
xmin = -pre_time*1000; % for plotting
xmax = post_time*1000; % for plotting
ticks = .2; % size of spike tick

% ---- initiate figure ----
if chans > 6
    row = 3;
    col = 3;
elseif chans > 4
    row = 3;
    col = 2;
elseif chans > 2
    row = 2; 
    col = 2;
elseif chans == 2
    row = 2;
    col = 1;
else
    row = 1;
    col = 1;
end

% ---plot spikes in cell array---
disp('plotting rasters...')
if cell_array > 0
    for ch = 1:chans
        for trial = 1:ntrials
            plotspikes = SpikeTimes{ch}(SpikeTimes{ch}(:,trial) > start(trial)-pre_time & SpikeTimes{ch}(:,trial) < start(trial)+post_time,trial); % extract the spikes in this plotting window for this channel/trial
            plotspikes = (plotspikes - start(ch))*1000; % make relative to the event onset, then change to ms
            
            subplot(row,col,ch); hold on;
            if ~isempty(plotspikes)
                plot([plotspikes plotspikes],[trial-ticks trial+ticks],'k'); % plot spikes as ms
            end
            clear plotspikes
        end % trial
        set(gca,'xlim',[xmin xmax],'ylim',ylim,...
        'box','off','tickdir','out','ytick',[],'yticklabel',[]);
        title(['Ch: ',num2str(ch)])
    end % channel 
% if not cell, plot spikes in trials
else      
    for trial = 1:ntrials
        plotspikes = SpikeTimes(SpikeTimes(:,trial) > start(trial)-pre_time & SpikeTimes(:,trial) < start(trial)+post_time,trial); % extract the spikes in this plotting window for this channel/trial
        plotspikes = (plotspikes - start(trial))*1000; % make relative to the event onset, then change to ms
        if ~isempty(plotspikes)
            plot([plotspikes plotspikes],[trial-ticks trial+ticks],'k'); % plot spikes as ms
            hold on;
        end
        clear plotspikes
    end
    set(gca,'xlim',[xmin xmax],'ylim',ylim,...
    'box','off','tickdir','out','ytick',[],'yticklabel',[]);
ylabel('Trials');
end

% print the raster and let the user look through in detail
if saving>0
    disp('saving figure...')
    if exist('name','var')
        print([name,'_raster.pdf'],'-dpdf');
    else
        print('raster.pdf','-dpdf');
    end
end
rH = gcf; 
    
end
