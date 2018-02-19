function [snips,sptimes,mask] = double_flood_fill( data,fs,varargin )
% [snips,sptimes,mask] = double_flood_fill( data,fs,(chanMap,lowThresh,highThresh,maxPts,maxChans) );
%
% Detects spikes in the n x m matrix "data" using a double flood-fill
% algorithm. Spikes are detected as a continuous group of spatio-temporally
% conencted sample points that all excede a lower threshold, and in which
% at least one point excedes an upper threshold. Optional arguments are 
% specifying using the name-value pair format.
%
%                   >>> INPUTS >>>
% data:
%   an n x m matrix, where n = number of points, m = number of channels
%
% fs:
%   scalar sampling rate
%
% (chanMap):
%   an m-dimensional vector specifying the channel mapping. By default, 
%   the columns of "data" will be assumed to be spatially continuous.
%
% (lowThresh):
%   scalar specifying the lower threshold to use for the flood fill algorithm.
%   Default = 2 (times the SD of each channel)
%
% (highThresh):
%   scalar specifying the upper threshold to use for the flood fill agorithm.
%   Default = 4 (times the SD of each channel)
%
% (maxPts):
%   scalar specifying the max number of points for any connected component to be
%   considered a valid spike. Anything longer than this is considered an artifact.
%   Default = inf
%
% (maxChans):
%   scalar specifying the maximum number of channels allowed per component.
%   Components with too many connected channels are considered movement-related artifacts.
%   Default = inf
%
% (artifact):
%   scalar specifying the value to consider the ith component an artifact
%   Default = inf
%
%
%                   <<< OUTPUTS <<<
% snips:
%   a k x p x m tensor, with k = # of points extracted for each spike, p = # of spikes,
%   and m = # of channels
%
% sptimes:
%   a p-dimensional vector specifying the time of the detected action potential
%
% mask:
%   an m x p matrix, where each element (i,j) is in the set [0,1]. The values are set
%   depending on the threshold crossing of each ith channel for the jth spikes, as:
%               0 : |(i,j)| < alpha
%               1 : |(i,j)| > beta
%       0 < x < 1 : alpha < |(i,j)| < beta
%
%           where   "alpha" = lowThresh  * SD( channel_i )
%           and     "beta"  = highThresh * SD( channel_i )
%
%   Thus, for each spike j, a mask is computed for each channel depending on the peak 
%   amplitude (i.e. negative voltage) of the ith spike on the ith channel. Channels not 
%   associated with the jth connected component of the flood-fill algorithm (i.e. where
%   the peak voltage is more negative than the lowThreshold AND continuous with at least
%   one region of amplitude more negative than the highThreshold) are completely maseked for 
%   that spike. Conversely, channels with a spike amplitude for the jth spike greater (more negative)
%   than the highThresh are completely unmasked. And anything in between is given a mask between
%   [0,1] depending on the value of the spike amplitude for that channel. 
%        
% This function is part of the Masked-EM algorithm, outlined by Kadir et
% al. 2015, but works for any data set with multiple channels
%
% By JMS, 8/7/2017

% check inputs / preallocate
[n,m] = size( data );
p = check_inputs();
n_sub = fs; % 1-second segments
numSegs = floor( n / n_sub ); % total number of 1-second segments
if numSegs * n_sub < n
    numSegs = numSegs + 1;
end
preTime = 0.0005;
postTime = 0.001;
preSamples = floor( preTime * fs);
postSamples = floor( postTime * fs);
totalSamples = postSamples + preSamples;
snips = cell(1,numSegs);
sptimes = snips;
mask = snips;

% get the standard deviation estimate for each channel
noise = zeros( m,1 );
for j = 1:10
    start = min( randi( numSegs,1 ), numSegs-1 ) * n_sub;
    noise = noise + mad( data(start:start+n_sub,:),1 )' / 0.6745;
end
noise = noise(p.chanMap) / j; % ordered according to our channel map

% create the start/stop vectors for each segment
start = n_sub * (0:numSegs-1) + 1;
stop = start + n_sub - 1;
stop(end) = n; % avoids extracting more data than available

% loop over data segments, find connected components via double ff algorithm
for j = 1:numSegs
    seg = data(start(j):stop(j),p.chanMap)'; % using p.chanMap ensures the channels are continuous
    
    % compute the low-threshold connected & adjacency matrix
    lowCheck = bsxfun( @gt,-seg,p.lowThresh*noise );
    highCheck = bsxfun( @gt,-seg,p.highThresh*noise ); 

    % find the connected components among low-threshold crosses
    components = bwconncomp( lowCheck ); 
    labels = labelmatrix( components ); % creates a matrix for the components

    % find only those components that also have at least one high-threshold cross
    finalLabels = unique( labels(highCheck) );
    finalLabels(finalLabels == 0) = [];
    clear lowCheck highCheck labels

    % now loop over these components & extract the spike waveform surrounding the peak of each
    counter = 1;
    nComponents = numel( finalLabels );
    snips{j} = nan( totalSamples,nComponents,m ); 
    sptimes{j} = nan( 1,nComponents );
    mask{j} = zeros( m,nComponents );
    spikes = nan( totalSamples*2,nComponents,m );
    sz = components.ImageSize;
    
    for i = finalLabels'
        thisLabel = components.PixelIdxList{i}; % equivalent to "find( labels == i )"
        [ch,pt] = ind2sub( sz,thisLabel ); 
        pt = unique( pt );
        ch = unique( ch );
        nPt = numel( pt ); 
        nCh = numel( ch );

        % skip if this component is too small or large, or too many channels
        if nPt <= 2 || nPt >= p.maxPts || nCh >= p.maxChans
            continue
        end

        % pull out the points corresponding to this component
        alpha = p.lowThresh * noise(ch)';
        beta = p.highThresh * noise(ch)';
        psi = seg(ch,pt)';
            
        % now find the spike onset as the center of mass of all connected points
        % of this component:
        %       t_spike = SUM{ t * psi^p } / SUM{ psi^p }
        %   where "p" is a power-weighting that determines alignment on spike peak 
        %   or center of mass (inf or 1, respectively)
        t_spike = sum( sum( bsxfun( @times,psi.^10,pt ) ) ) / sum( sum( psi.^10 ) );
        closestPt = round( t_spike );
         
        % create our mask for the ith component as:
        %   max_t{ min{ (-V(t,c) - alpha) / (beta - alpha), 1 } }
        psi_masked = max( min( bsxfun( @rdivide,bsxfun( @minus,-psi,alpha ),(beta - alpha) ),1 ) );

        % skip if this point is too close to the edge of the recording
        if (closestPt - preSamples*2 < 1) || (closestPt + postSamples*2 > n_sub)
            continue
        end
        
        % pull out the spike waveform
        thisSpike = seg(:,closestPt-preSamples*2:closestPt+postSamples*2-1)';
        
        % skip if this is an artifact
        if any( any( abs( thisSpike ) >= p.artifact ) | (any( thisSpike ) > p.artifact/2) )
            continue
        end
        
        % skip this point if another channel in the connected component
        % has a large deflection within the pre/post time of this spike
        %       i.e. overlapping spikes
        if p.removeOverlap
            if any( abs( thisSpike(preSamples:preSamples+postSamples,ch) ) > max( abs( psi ) ) ) 
                continue
            end
        end
        
        % add to our "spikes" matrix and update spike times and mask
        spikes(:,counter,p.chanMap) = thisSpike; % back into original channel order as supplied
        sptimes{j}(counter) = t_spike; 
        mask{j}(p.chanMap(ch),counter) = psi_masked; % we un-do the channel mapping and only change channels associated with this component
        counter = counter + 1;
    end

    % now we need to align the spike waveforms since t_spike may not be an exact sample point.
    % we use cubic spline interpolation surrounding this point, 
    % extract less data before/after this point than the original data
    snips{j} = interpolate_spikes( spikes,fs,(sptimes{j}-round( sptimes{j} ))/fs + preTime*2,...
                                        preTime,postTime,0.6,1 ); % 1x interpolation with slight smoothing
    sptimes{j} = (sptimes{j} + start(j) - 1); % to make relative to the start of "data", not "seg"

    % finally, remove any nan's in our matrices
    badInds = isnan( sptimes{j} );
    if any( badInds )
        snips{j}(:,badInds,:) = [];
        sptimes{j}(badInds) = [];
        mask{j}(:,badInds) = [];
    end
end

% now concatenate the data
snips = [snips{:}];
sptimes = [sptimes{:}];
mask = [mask{:}];

%% HELPER FUNCTIONS
    function p = check_inputs()
        names = {'chanMap','lowThresh','highThresh','maxPts','maxChans','artifact','removeOverlap'};
        defaults = {1:m,2,4,inf,inf,inf,false};
        
        p = inputParser();
        for j = 1:numel( names )
            p.addParameter( names{j},defaults{j} );
        end
        
        p.parse( varargin{:} );
        p = p.Results;
    end
end
