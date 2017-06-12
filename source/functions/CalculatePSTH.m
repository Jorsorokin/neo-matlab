function [psth,psthAvg,variance,FF,ISI,bw] = CalculatePSTH(SpikeTimes,start,varargin)
% function [psth] = CalculatePSTH(SpikeTimes,EventTimes,varargin)
%
%   Calculates peri-stimulus time histograms (PSTHs) from a matrix or cell
%   array of spiketimes (in seconds), given a user-defined bin-width. Plots
%   PSTHs and Fano-Factor and saves in current directory.
%
%   If the user does not supply a bin width, the function will search for
%   the optimal bin width via the convex optimization problem:
%       
%           min C(bw) = (2K - V) / bw^2
%     
%       where K = average histogram
%             V = MSE of histogram of spike rate of all trials
%   
%   * For more detail, refer to Shimazaki & Shinomoto, Neural Comp. 2007
%   
%
%               >>> INPUTS >>>
%
% Required:
%   SpikeTimes = matrix or cell of spiketimes in SECONDS
%           If matrix, will assume columns are trials, rows are spiketimes.
%           If cell, will assume each cell is a trial, and in each cell
%           columns are channels, rows are spiketimes.
%   start = n-element vector containing times of events (in seconds)...if
%           spiketimes are relative (i.e. blocks of spikeitmes, each block relative
%           to stim onset), then user should define "start" as the time of the
%           onset of the stim relative to the start of the block. For instance, if
%           stim is 6 seconds into the start of eah block, then set "start" = 6.
%           * if spikes NOT relative, start = nx1 array of starting times.
% 
% Optional:
%   pre_time = time (in SECONDS) to subtract from starting time...
%           makes plots relative to stim onset (default = 1s);
%   post_time = time (in SECONDS) to add to starting time...
%           (default = 1s);
%   bin_width = bin (ms) for PSTH calculation. Default = [] (search for optimal bw)
%   name = name to save figure (default = "psth.pdf")
%   saving = 0 or 1 (default 1). If 1, saves figures to current directory.
%   
%
%               <<< OUTPUTS <<<
%
%   psth = bin-counts of spiketimes occuring within specified time range.
%           If SpikeTimes is a cell array, psth is an nxtrialsxchan matrix.
%           If SpikeTimes is a matrix, psth is an nxtrials matrix
%   psthTrialAvg = average of psth across trials per channel
%   variance = average variance of psth across trials per channel
%   FF = fano factor as "variance / psthAvg"
%   ISI = inter-spike intervals for each bin, for each trial
%   bw = the optimized histogram bin width
%
%
% Example:
%   SpikeTimes{1} = sort(rand(100,100)); % fake spiketimes for ch1
%   SpikeTimes{2} = sort(rand(100,100)); % fake spiketimes for ch2
%   [psth,trialAvg,varAvg] = CalculatePSTH(SpikeTimes,.5,.2,.5,10,'control',1)
%       % plots histogram and fano factor for each channel, from -0.2s : 0.5s
%       around the starting point, (here 0.5s into the SpikeTimes). Save
%       the figures and appends "control" to the figure name
%   
% By JMS, 11/13/2015
%-------------------------------------------------------------------------------

% check optionals
if nargin>2 && ~isempty(varargin{1})
    pre_time = varargin{1};
else pre_time = 1; end % default 1s before start
if nargin>3 && ~isempty(varargin{2})
    post_time = varargin{2};
else post_time = 1; end % default 1s after start
if nargin>4 && ~isempty(varargin{3})
    bw = varargin{3}; else bw = []; end % default 10s bin width
if nargin>5 && ~isempty(varargin{4})
    name = varargin{4}; end % default no name for saving figs
if nargin>6 && ~isempty(varargin{5})
    saving = varargin{5};
else saving = 1; end


% check if SpikeTimes is cell or matrix
if iscell(SpikeTimes)
     nchans = max(size(SpikeTimes));
     ntrials = size(SpikeTimes{1},2);
else
    ntrials = size(SpikeTimes,2);
    nchans = 1;
end

% error check if size(start) ~= ntrials....if so, repeat "start" so that it
% has same dimension as ntrials
if numel(start) == 1
    start = ones(ntrials,1)*start;
elseif numel(start) ~= ntrials
    assert('Error: starting times and # trials not same dimension');
end

% lastbin of the histogram segmentation
lastBin = ceil((post_time*1000)); % last bin edge in ms


% SEARCH ALGORITHM IF NO BW IS PROVIDED
if isempty(bw)
    bw = optimize_bw(SpikeTimes,pre_time,post_time,nchans);
    if saving == 1
        if exist('name','var')
            print([name,'_C-vs-BW.pdf'],'-dpdf');
       else
            print('C-vs-BW.pdf','-dpdf');
        end
    end
end


% Now compute the histogram edges with the optimal BW or provided BW
edge = EdgeCalculator(bw,-pre_time*1000,lastBin); % extract edges for psth
xmin = edge(1); % for plotting
xmax = edge(end); % for plotting
xlim = [xmin xmax];


% ======= begin loop =======
disp('calculating PSTHs...');
clear psth
try 
    % COMPUTE PSTH
    if iscell(SpikeTimes)
         psth = zeros(numel(edge),ntrials,nchans);
         ISI = zeros(numel(edge),ntrials,nchans);
         for ch = 1:nchans
             sptm = SpikeTimes{ch};
             
             [psthTemp,isiTemp] = bin_spikes(sptm,ntrials,pre_time,post_time,start,edge,bw);
            
             psth(:,:,ch) = psthTemp; 
             ISI(:,:,ch) = isiTemp;
             clear psthTemp isiTemp
         end
    else % if not a cell array
        [psth,ISI] = bin_spikes(SpikeTimes,ntrials,pre_time,post_time,start,edge,bw);
    end
    
    % COMPUTE FANO FACTOR AND VARIANCE
    [variance,FF,psthAvg] = compute_fanofactor(psth);

   % === PSTH bar plots ===
   psthlim = [0 max(max(psthAvg))]; % for plotting
   FFlim = [0 max(max(FF))];
   ISIlim = [0 max(max(nanmean(ISI,2)))];
   
   if nchans > 6
       row = 3;
       col = 3;
   elseif nchans > 4
       row = 3;
       col = 2;
   elseif nchans > 2
       row = 2; 
       col = 2;
   elseif nchans == 2
       row = 2;
       col = 1;
   elseif nchans == 1
       row = 3;
       col = 1;
   end

   % plot variance for channels in separate figure if Spiketimes is cell array
   if nchans > 1
       f1 = figure;
       for ch = 1:nchans
           subplot(row,col,ch);
           plot_psth(psthAvg,ch,edge,xlim,psthlim);
       end
       f2 = figure;
       for ch = 1:nchans
           subplot(row,col,ch);
           plot_FF(FF,ch,edge,xlim,FFlim);
       end
       f3 = figure;
       for ch = 1:nchans
           subplot(row,col,ch);
           plot_isi(ISI,ch,edge,xlim,ISIlim);
       end
   else
       f1 = figure;
       subplot(row,col,1);
       plot_psth(psthAvg,1,edge,xlim,psthlim);
       title('PSTH');
       subplot(row,col,2);
       plot_FF(FF,1,edge,xlim,FFlim);
       title('Fano Factor');
       subplot(row,col,3);
       plot_isi(ISI,1,edge,xlim,ISIlim);
       title('Mean Inter-spike Interval');
       
       psth = squeeze(psth);
       ISI = squeeze(ISI);
   end
   
   % print the psth, append "name" if it exists
   if saving==1
       if exist('f1','var')
           figure(f1);
           if exist('name','var')
               print([name,'_PSTH.pdf'],'-dpdf');
           else
               print('PSTH.pdf','-dpdf');
           end
       end
       if exist('f2','var')
           figure(f2);
           if exist('name','var')
               print([name,'_FanoFactor.pdf'],'-dpdf');
           else
               print('FanoFactor.pdf','-dpdf');
           end
       end
       if exist('f3','var')
           figure(f3);
           if exist('name','var')
               print([name,'_ISI.pdf'],'-dpdf');
           else
               print('ISI.pdf','-dpdf');
           end
       end
   end   
catch            
    disp('Error in file:...problem computing PSTH');
    msg = lasterr;
    disp(msg);
end % catch loop

end
%% Functions

% FANO FACTOR PLOT
function plot_FF(data,ch,edge,xlim,ylim)
    plot(edge,data(:,ch));
    hold on;
    plot(edge,smooth(data(:,ch)),'r','linewidth',2);
    set(gca,'xlim',xlim,'ylim',ylim,...
        'box','off','tickdir','out');
    title(['Ch: ',num2str(ch)])
    ylabel('var / spike rate');
end

% PSTH PLOT
function plot_psth(data,ch,edge,xlim,ylim)
    h = bar(edge,data(:,ch),'histc');
    set(gca,'xlim',xlim,'ylim',ylim,...
        'box','off','tickdir','out');
    set(h,'facecolor','k','edgecolor','k');
    title(['Ch: ',num2str(ch)])
    ylabel('Spike Rate');
end

% ISI PLOT
function plot_isi(data,ch,edge,xlim,ylim)
    D = nanmean(squeeze(data(:,:,ch)),2)*1000;
    col = [0,.6,.7];
    h = bar(edge,D,'histc');
    set(gca,'xlim',xlim,'ylim',ylim*1000,...
        'box','off','tickdir','out');
    set(h,'facecolor',col,'edgecolor',col);
    title(['Ch: ',num2str(ch)]);
    ylabel('ISI');
end 

% SPIKE HISTOGRAM CALCULATION
function [psth,ISI] = bin_spikes(sptime,ntrials,pre_time,post_time,start,edge,bw)
    psth = zeros(numel(edge),ntrials); 
    ISI = zeros(numel(edge),ntrials);
    for trial = 1:ntrials
        if ~isempty(sptime)
            psthSpikes = sptime(sptime(:,trial) > start(trial)-pre_time & sptime(:,trial) < start(trial)+post_time,trial); % extract spikes occuring within pre/post times of start time
            psthSpikes = (psthSpikes - start(trial))*1000;
            for j = 1:numel(edge)-1
                sp = psthSpikes(psthSpikes >= edge(j) & psthSpikes < edge(j+1));
                psth(j,trial) = numel(sp);
                ISI(j,trial) = nanmean(diff(sp));
                clear sp
            end
            clear psthSpikes
        end
    end
    psth = psth / (bw/1000); % make spikes relative to start and extract spike count per bin, divide by bw in seconds to get firing rate
    ISI = ISI / 1000; % convert ISI's to seconds, not ms
end

% FANO FACTOR CALCULATION VIA SPIKE RATE FLUCTUATION
function [variance,FF,psthAvg] = compute_fanofactor(psth)
    psthAvg = mean(psth,2); % take mean across trials, squeeze into nxchan array
    variance = (bsxfun(@minus,psthAvg,psth)).^2;
    variance = squeeze(mean(variance,2));
    psthAvg = squeeze(psthAvg);
    FF = variance ./ psthAvg;
end

% OPTIMAL HISTOGRAM BIN WIDTH
function optBW = optimize_bw(sptime,pretime,posttime,nchans)
    binstart = 40; % 50 bins 
    binstop = 400; % 400 bins
    minbin = (posttime + pretime) / 200; 
    if minbin < .01 % hard lower limit of 20ms
        minbin = .01;
    end
    if minbin > .04
        minbin = .04; % hard upper limit of 40
    end
    N = binstart:binstop;
    D = (posttime + pretime) ./ N;
    N(D<minbin) = [];
    D(D<minbin) = []; % eliminate really small bin widths
   
    C = zeros(numel(N),nchans);
    
    % === loop through bin widths ===
    %for i = 1:length(binrange) %increment by 2 to increase speed
    disp('Searching for optimal bin width')
    disp('');
    for i = 1:numel(N)
        if mod(i,10) == 0
            fprintf('%s',' . ');
        end
        
        % get the bins and the bin spacing
        bins = linspace(-pretime,posttime,N(i)+1);
        bw = D(i) * 1000;
        
        % get psth for all trials/chans
        for ch = 1:nchans
            if nchans > 1
                S = reshape(sptime{ch},size(sptime{ch},1)*size(sptime{ch},2),1); % reshape to long array
            else
                S = reshape(sptime,size(sptime,1)*size(sptime,2),1);
            end
            S(S==0) = []; % remove all 0 indices
            S = S - pretime; % make relative to pretime
            psth(:,ch) = histc(S,bins);
        end
        
        % get the fanofactor, variance, psth mean
        [~,FF,ki] = compute_fanofactor(psth);
        
        % eliminate the last edge since we incremented by N(i)+1, then take
        % mean count for further calculations   
        FF(end,:) = [];
        ki(end,:) = [];  
        
        % mean count for all bins i = 1:N
        K = mean(ki); 
        
        % find the MSE between ki and K for each channel
        V = bsxfun(@minus,ki,K).^2;
        
        % now compute C(bw) for the current V and kAvg, for each channel
        ki = ki.*FF; % weight by the fano factor for non-poisson process
        C(i,:) = mean(bsxfun(@minus,2*ki,V) / D(i)^2);
        
        clear ki K V psth
    end
    
    % find minimum of C(i), which gives index of optimal bin width
    C = mean(C,2);
    C2 = smooth(C,21);
    C2 = reshape(C2,size(C));
    C2 = C2(1:end-1,:);
    D = D(1:end-1);
    C = C(1:end-1,:);
    [cmin,idx] = min(C2);
    optBW = mean(D(idx))*1000;
    
    % plot the optimized C vs. D
    figure; hold on; 
    scatter(D*1000,C,'k');
    plot(D*1000,C2,'b','linewidth',2);
    scatter(optBW,cmin,'ro','markerfacecolor','r'); 
    axis tight
    ylabel('C');
    xlabel('BW');
    disp(' ');
    fprintf('optimal bin width: %f %s',optBW,'ms');
    disp(' ');
end