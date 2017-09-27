function [sptm,spsnip,spmin,spmax] = detectSpikes( data,Fs,varargin )
% 
% [sptm,spsnip,spMin,spMax] = detectSpikes( data, Fs, (thresh, sumChans) )
%
% detects spikes from single or multi-channel voltage traces. Uses the 
% multi-teager energy operator (MTEO), sped up by converting to a MEX file
% (MTEO_mex) to filter the raw traces for high-energy, high-frequency 
% deflections. Uses an automatic threshold (4 * median( signal/0.6745 ))
% to detect spikes. 
%
% Inputs:
%   data:
%       an n x m matrix, with n = samples, m = channels. 
%
%   Fs:
%       the sampling rate (in Hz)
%
%   (thresh):
%       scalar threshold to mulitply the noise estimate
%       median( abs( signal ) / 0.6745 ); Default value = 4;
%
%   (sumChans):
%       boolean (0 or 1). If 0, each channel will be treated separately and
%       spike times/snips will be stored into a cell array, with elements
%       in the cell array corresponding to channels. If 1, spikes will be 
%       detected by summing the MTEO-filtered signals together. 
%       This is useful for (a) keeping the same size data matrix of spike 
%       snips/times for each channel, and (b) detecting spikes that may
%       have occurred on multiple channels as in tetrode recordings.
%       Default value = 0;
%
%   (artifact):
%       the voltage level considered an artifact. Anything lower than this
%       level will be rejected. Default = 1000 (uV)
%
% Outputs:
%   sptm:
%       the time (in samples) of the detected spikes
%
%   spsnip:
%       the waveforms of the detected spikes. [-.05ms, 1ms] of data will be
%       extracted surrounding each spike.
%
%   spmin:
%       the minimum value (ie. peak amplitude) of the detected spikes
%   
%   spmax:
%       the maximum value following the spike trough (i.e. the AHP peak)
% 
% By JMS, 5/20/2017

% check inputs
if nargin > 2 && ~isempty( varargin{1} )
    thresh = varargin{1};
else thresh = 4; end
if nargin > 3 && ~isempty( varargin{2} )
    sumChans = varargin{2};
else sumChans = 0; end
if nargin > 4 && ~isempty( varargin{3} )
    artifact = varargin{3};
else artifact = 1000; end

% params and constants
[~,m] = size( data );
PREPTS = round( 0.0015 * Fs ); % samples to extract prior to spike
POSTPTS = round( 0.0015 * Fs ); % samples to extract after spike
MINDIST = round( 0.0005 * Fs ); % avoids detecting close spikes

% perform multi-teager energy filtering
mte = MTEO( data,2,5 ); % 4-point averaging, starting at 2-points

% get the noise estimates
noise = mad( mte,1 ) / 0.6745;
threshold = thresh * noise;

% detect spikes for each channel or for the summed channels
switch sumChans
    case 0
        sptm = cell( 1,m );
        spsnip = sptm;
        for ch = 1:m   
            sptm{ch} = get_spike_times( mte(:,ch),threshold(ch) );
            [spsnip{ch},sptm{ch}] = get_snips( sptm{ch},data(:,ch) );
            [spsnip{ch},sptm{ch},spmin{ch},spmax{ch}]...
                = get_aligned_spikes( spsnip{ch},sptm{ch} );
        end
    case 1
        sptm = get_spike_times( sum( mte,2 ),sum( threshold ) );
        [spsnip,sptm] = get_snips( sptm,data );
        [spsnip,sptm,spmin,spmax] = get_aligned_spikes( spsnip,sptm );
end

if m == 1
    sptm = sptm{1};
    spsnip = spsnip{1};
end

%% helper functions
    function times = get_spike_times( signal,threshold )

        % find the peaks
        [~,times] = findpeaks( signal,'minpeakheight', threshold,...
            'minpeakdistance',MINDIST );

        % remove boundary spikes
        times(times<=PREPTS | times>=size(signal,1)-POSTPTS) = []; 
    end

    function [snips,sptm] = get_snips( sptm,data )

        % pull out the spike snips
        [snips,sptm,spamp] = get_spike_snips( sptm,data,PREPTS,POSTPTS );

        % get rid of spikes with positive peak or too large an amplitude
        maxPk = max(spamp,[],2);
        bad = (maxPk > 0) | (mean( abs( spamp ),2 ) >= artifact);
        snips(:,bad,:) = [];
        sptm(bad) = []; 
    end  

    function [snips,sptm,vMin,vMax] = get_aligned_spikes( spikes,sptm )
        
        % check size
        [n,nSp,nChan] = size( spikes );
        
        % resample the spikes for the alignment / spline interpolation
        upspikes = zeros( n*4,nSp,nChan );
        x = 1:n;
        xx = linspace(1,n,n*4); % the interpolant x-values
        for c = 1:nChan
            upspikes(:,:,c) = csaps( x,spikes(:,:,c)',.4,xx )'; % smooths the interpolant
        end
        upFs = Fs * 4;

        % find the minimums/maximums
        [vMin,vMax,tMin] = spikeHeight( upspikes,upFs,.001,.003 );
        clear upspikes

        % align the spikes
        %usnips = alignSpikes( upspikes,upFs,tMin,.0005,.001);
        tSpike = median( tMin,2 );
        snips = interpolate_spikes( spikes,Fs,tSpike,0.0005,0.001,0.6,1 ); % no smoothing and downsample
        bad = (sum( isnan( snips(1,:,:) ),3 ) > 0) | (sum( sum( abs( snips )>=abs(artifact),3 ) ) > 0);
        snips(:,bad,:) = [];
        vMin(bad,:) = [];
        vMax(bad,:) = [];
        sptm(bad) = [];
        
        % check if any remaining
        if isempty( sptm )
            disp( 'no spikes found' );
            return
        end       
    end
 
end

