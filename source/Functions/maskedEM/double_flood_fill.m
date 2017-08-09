function [snips,times,mask] = double_flood_fill( data,fs,varargin )
% [snips,times,mask] = double_flood_fill( data,fs,(chanMap,lowThresh,highThresh) );
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
%
%                   <<< OUTPUTS <<<
% snips:
%   a k x p x m tensor, with k = # of points extracted for each spike, p = # of spikes,
%   and m = # of channels
%
% times:
%   a p-dimensional vector specifying the peak time of the detected action potential
%
% mask:
%   an m x p matrix, where each element (i,j) is in the set [0,1]. The values are set
%   depending on the threshold crossing of each ith channel for the jth spikes, as:
%               0 : |(i,j)| < alpha
%               1 : |(i,j)| > beta
%       0 < x < 1 : alpha < |(i,j)| < beta
%
%           where "alpha" = lowThresh * SD_ch_k
%           and "beta" = highThresh * SD_ch_k
%
%   Thus, for each spike j, a mask is computed for each channels depending on the peak 
%   amplitude (i.e. negative voltage) of the ith spike on the ith channel. Channels not 
%   associated with the jth connected component of the flood-fill algorithm (i.e. where
%   the peak voltage is more negative than the lowThreshold AND continuous with at least
%   one region of amplitude more negative than the highThreshold) are completely maseked for 
%   that spike. Conversely, channels with a spike amplitude for the jth spike greater (more negative)
%   than the highThresh are completely unmasked. And anything in between is given a mask between
%   [0,1] depending on the value of the spike amplitude for that channel. 
%        
% This function is part of the Masked-EM algorithm, outlined by Kadir et al. 2015
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
preSamples = floor( 0.0005 * fs);
postSamples = floor( 0.0015 * fs);
totalSamples = postSamples + preSamples;
snips = cell(1,numSegs);
times = snips;
mask = snips;
kernel = [0 1 0; % eliminates isolated low-threshold crossings
          1 0 1;
          0 1 0];

% get the standard deviation estimate for each channel
sd = zeros( 1,m );
for j = 1:10
    start = min( randi( numSegs,1 ), numSegs-1 ) * n_sub;
    sd = sd + median( abs( data(start:start+n_sub,:) ) ) / 0.6745;
end
sd = sd(p.chanMap) / j; % ordered according to our channel map

% create the start/stop vectors for each segment
start = n_sub * (0:numSegs-1) + 1;
stop = start + n_sub;
stop(end) = n; % avoids extracting more data than available

% loop over data segments, find connected components via double ff algorithm
for j = 1:numSegs
    seg = data(start(j):stop(j),p.chanMap); % using p.chanMap ensures the channels are continuous. if not, isolated "spike islands" will occur

    % compute the low-threshold connected & adjacency matrix
    lowCheck = bsxfun( @gt,-seg,p.lowThresh*sd ); % low memory implementation
    highCheck = bsxfun( @gt,-seg,p.highThresh*sd ); 
    connected = conv2( single( lowCheck ),kernel,'same' );
    connected = lowCheck & connected > 0; % gets rid of self-adjacency (spot noise)

    % find the connected components among low-threshold crosses
    components = bwconncomp( connected ); % finds connected components
    labels = labelmatrix( components ); % creates a matrix for the components

    % find only those components that also have at least one high-threshold cross
    finalLabels = unique( labels(highCheck) );
    finalLabels(finalLabels == 0) = [];
    %clear connected lowCheck highCheck components

    % now loop over these components & extract the spike waveform surrounding the peak of each
    counter = 1;
    nComponents = numel( finalLabels );
    snips{j} = nan( totalSamples,nComponents,m ); 
    times{j} = nan( 1,nComponents );
    mask{j} = zeros( m,nComponents,'single' );
    spikes = nan( totalSamples*2,nComponents,m );
    sz = components.ImageSize;
    for i = finalLabels

        % find points connected to this component
        thisLabel = components.PixelIdxList{i}; % equivalent to "find( labels == i )"
        [pt,ch] = ind2sub( sz,thisLabel ); 
        pt = unique( pt );
        ch = unique( ch );
        nPt = numel( ch ); 
        nCh = numel( pt );

        % create our mask for the ith component as:
        %  max_t{ psi(t,c) = min{ (-V(t,c) - alpha) / (beta - alpha), 1 } }
        alpha = lowThresh * sd(ch);
        beta = highThres * sd(ch);
        psi = zeros( nPt,nCh );
        for c = 1:nCh
            for t = 1:nPt
                psi(t,c) = seg(pt(t),ch(c));
            end
        end
        psi = min( bsxfun( @rdivide,bsxfun( @minus,-psi,alpha ),(beta - alpha) ),1 );
        psi_masked = max( psi ); % maximum across points for each channel

        % now find the spike onset as the center of mass of all connected points
        % of this component:
        %       t_spike = SUM{ t * psi^p } / SUM{ psi^p }
        %   where "p" is a power-weighting that determines alignment on spike peak 
        %   or center of mass (inf or 1, respectively). p = 2 is a good balance. 
        t_spike = sum( sum( pt * psi.^2 ) ) / sum( sum( psi.^2 ) );
        closestPt = round( t_spike ); 

        % skip if this point is too close to the edge of the recording
        if closestPt - (preSamples*2) < 1 || closestPt + (postSamples*2) > n_sub
            continue;
        end

        % now that we have t_spike, use the closest actual sample point to pull out the spike across channels
        spikes(:,counter,p.chanMap) = seg( closestPt-preSamples*2:t_spike+postSamples*2-1,: ); % back into original channel order as supplied
        times{j}(counter) = t_spike; 

        % update our mask matrix for this component
        mask{j}(p.chanMap(ch),counter) = psi_masked; % we un-do the channel mapping and only change channels associated with this component
        counter = counter + 1;
    end

    % now we need to align the spike waveforms according to their center of mass.
    % since t_spike may not be an exact sample point, we use cubic spline interpolation surrounding this 
    % point, extract less data before/after this point than the original data, then down sample
    snips{j} = get_aligned_spikes( spikes,times{j},preSamples,postSamples );
    times{j} = (times{j} + start(j) - 1) * fs; % to make relative to the start of "data", not "seg"

    % finally, remove any nan's in our matrices
    badInds = isnan( snips{j}(1,:,1) );
    if any( badInds )
        snips{j}(:,badInds,:) = [];
        times{j}(badInds) = [];
        mask{j}(:,bad) = [];
    end
end


%% HELPER FUNCTIONS
function p = check_inputs()
    names = {'chanMap','lowThresh','highThresh'};
    defaults = {1:m,2,4};

    p = input_parser;
    for k = 1:numel( names )
        p.addParameter( names{k},defaults{k} );
    end

    p = p.parse( varargin{:} );
    p = p.results;
end


end
