function [sptm,spamp,spsnip] = RemoveCorrSpikes(d,Fs,sptm,varargin)
% function [sptm,spamp,spsnip] = RemoveCorrSpikes(d,Fs,sptm,varargin)
%----------------------------------------------------------------
%   Removes spikes from multichannel data that are highy correlated with
%   other channels...usually on order of .85 or so. This technique helps eliminate 
%   spikes due to movement artifacts.
%   Refer to Paralikar et. al. 2010 for more information.
%
%               >>> INPUTS >>>
% Required:
%   d = data matrix. note, if vector, function will end without doing
%           anything. Needs to be in COLUMN format
%   Fs = sampling rate
%   sptm = matrix of spike times in samples. also in column format
% Optional: 
%   corr_thresh = correlation coefficient threshold for removing spikes. 
%           Default = 0.8
%   spamp = matrix of spike amplitudes...doesn't contribute to the
%           function, but if included indexes removed from sptm will also be
%           removed from spamp.
%   spsnip = matrix of spike snips...if provided, will make function faster
%           as snips will not need to be calculated for each sptm point.
%
%               <<< OUTPUTS <<<
%   sptm = spiketimes after eliminating bad spikes...bad spikes will have
%           value of 0. 
%   spamp = spike amps after eliminating bad spikes (OPTIONAL)
%   spsnip = matrix of spikesnips after eliminating bad spikes (optional)
%
% By JMS, 11/12/2015
%---------------------------------------------------------------

% check data dimensions
if min(size(d))==1
    error('Data not matrix...ending function');
end

% check inputs
if nargin>3 && ~isempty(varargin{1})
    corr_thresh = varargin{1}; 
else corr_thresh = 0; end
if nargin>4 && ~isempty(varargin{2})
    spamp = varargin{2}; 
else spamp = 0; end
if nargin>5 && ~isempty(varargin{3})
    spsnip = varargin{3};
end


disp('Removing highly correlated spikes ...')
prezeros = sum(sptm==0);

% compute inter-electrode correlation for each spike and discard
% spike if significant correlation between the electrodes and the
% electrode of interest...see Paralikar et al. 2010.
for ch = 1:size(d,2)
    cross_chans = 1:size(d,2);
    cross_chans(ch) = []; % channels for cross correlation compared to ref chan
    if ~isempty(cross_chans)
        for spike = 1:length(sptm)% loop through spikes
            if sptm(spike) ~= 0
                snip = sptm(spike)/Fs*1000; % convert to ms and take +/- 1ms from spike
                if snip >= 2 && snip < size(d,1)/Fs*1000-1
                    snip = (snip-1:snip+1)/1000*Fs; % convert back to samples
                elseif snip < 2
                    snip = (snip:snip+1)/1000*Fs;
                elseif snip >= size(d,2)/Fs*1000-1
                    snip = (snip-2:snip)/1000*Fs;
                end
                
                % compute correlation across electrodes per spike
                snip = floor(snip);
                [r,~]=corrcoef([d(snip(1):snip(end),ch),d(snip(1):snip(end),cross_chans)]); 
                clear snip

                % check correlation value
                if max(r(2:end,1))>corr_thresh % if correlation large enough, discard spike
                    if exist('spamp','var')
                        spamp(spike,ch) = 0; % remove spike amps if they are given
                    end
                    if exist('spsnip','var') && ~isempty('spsnip')
                        spsnip(:,spike,ch) = nan;
                    end
                    sptm(spike,ch) = 0; % remove spiketimes with high correlations
                end
            end % if sptm(spike) ~= 0 
        end %spike loop
    end % if ~isempty(cross_chan)
end % channel loop


postzeros = sum(sptm==0);
fprintf('# spikes removed: %s\n',num2str(abs(postzeros-prezeros))); % print # of spikes removed in each channel

end