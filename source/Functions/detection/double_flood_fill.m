function [snips,sptimes,mask,overlap] = double_flood_fill( data,fs,varargin )
    % [snips,sptimes,mask,overlap] = double_flood_fill( data,fs,varargin );
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
    % (removeOverlap):
    %   boolean flag. If true, spikes that overlap with another spike 
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
    % overlap:
    %   a logical vector indicating whether any detected spike is overlapping
    %   with another spike. If "removeOverlap" is true, then overlapping spikes
    %   will be eliminated
    %        
    % This function is part of the Masked-EM algorithm, outlined by Kadir et
    % al. 2015, but works for any data set with multiple channels
    %
    % By JMS, 8/7/2017

    % check inputs / preallocate
    [n,m] = size( data );
    p = check_inputs();
    lowThresh = p.lowThresh;
    highThresh = p.highThresh;
    artifact = p.artifact;
    maxChans = p.maxChans;
    maxPts = p.maxPts;
    chanMap = p.chanMap;
    removeOverlap = p.removeOverlap;

    n_sub = fs; % length of each segment
    numSegs = floor( n / n_sub ); % total number of 1-second segments
    if numSegs * n_sub < n
        numSegs = numSegs + 1;
    end
    
    preTime = 0.0005;
    postTime = 0.001;
    preSamples = floor( preTime * fs);
    postSamples = floor( postTime * fs);
    totalSamples = postSamples + preSamples;
    snips = cell( 1,numSegs );
    sptimes = cell( 1,numSegs );
    mask = cell( 1,numSegs );
    contaminationTime = 0.0002 * fs;
    rng(1);


    % get the standard deviation estimate for each channel
    noise = zeros( m,1 );
    for j = 1:10
        start = min( randi( numSegs,1 ), numSegs-1 ) * n_sub;
        noise = noise + mad( data(start:start+n_sub,:),1 )' / 0.6745;
    end
    noise = noise(chanMap) / j; % ordered according to our channel map

    % create the start/stop vectors for each segment
    start = n_sub * (0:numSegs-1) + 1;
    stop = start + n_sub - 1;
    stop(end) = n; % avoids extracting more data than available

    % data = parallel.pool.Constant( data );
    for j = 1:numSegs
        seg = data(start(j):stop(j),chanMap)';
        %seg = data.Value(start(j):stop(j),chanMap)'; 

        % compute the low-threshold connected & adjacency matrix
        lowCheck = bsxfun( @gt,-seg,lowThresh*noise );
        highCheck = bsxfun( @gt,-seg,highThresh*noise ); 

        % find the connected components among low-threshold crosses
        components = bwconncomp( lowCheck ); 
        labels = labelmatrix( components ); % creates a matrix for the components

        % find only those components that also have at least one high-threshold cross
        finalLabels = unique( labels(highCheck) );
        finalLabels(finalLabels == 0) = [];

        % now loop over these components & extract the spike waveform surrounding the peak of each
        counter = 1;
        nComponents = numel( finalLabels );
        sz = components.ImageSize;
        tempSpikes = nan( totalSamples*2,nComponents,m );
        tempTimes = nan( 1,nComponents );
        tempMask = zeros( m,nComponents );

        for i = finalLabels'
            thisLabel = components.PixelIdxList{i}; % equivalent to "find( labels == i )"
            [ch,pt] = ind2sub( sz,thisLabel ); 
            pt = unique( pt );
            ch = unique( ch );
            nPt = numel( pt ); 
            nCh = numel( ch );

            % skip if this component is too small or large, or too many channels
            if nPt <= 2 || nPt >= maxPts || nCh >= maxChans
                continue
            end

            % remove channels where peak is too far in time away from largest
            % peak 
            if nCh > 1 
                [peakVals,peakTimes] = min( seg(ch,pt),[],2 );
                [~,bestChan] = min( peakVals );
                badCh = abs( peakTimes - peakTimes(bestChan) ) > contaminationTime;
                ch(badCh) = [];
            end

            % pull out the spike
            psi = seg(ch,pt)';
            alpha = lowThresh * noise(ch)';
            beta = highThresh * noise(ch)';

            % now find the spike onset as the center of mass of all connected points
            % of this component:
            %       t_spike = SUM{ t * psi^p } / SUM{ psi^p }
            %   where "p" is a power-weighting that determines alignment on spike peak 
            %   or center of mass (inf or 1, respectively)
            t_spike = sum( sum( bsxfun( @times,psi.^2,pt ) ) ) / sum( sum( psi.^2 ) );
            closestPt = round( t_spike );

            % skip if this point is too close to the edge of the recording
            if (closestPt - preSamples*2 < 1) || (closestPt + postSamples*2 > n_sub)
                continue
            end

            % skip if this is an artifact
            thisSpike = seg(:,closestPt-preSamples*2:closestPt+postSamples*2-1)';
            if any( any( abs( thisSpike ) >= artifact ) )
                continue
            end

            % create our mask for the ith component as:
            %   max_t{ min{ (-V(t,c) - alpha) / (beta - alpha), 1 } }
            psi_masked = max( min( bsxfun( @rdivide,bsxfun( @minus,-psi,alpha ),(beta - alpha) ),1 ) );

            % add to our "spikes" matrix and update spike times and mask
            %tempOverlap(counter) = any( abs( thisSpike(preSamples:preSamples+postSamples,ch) ) > max( abs( psi ) ) );
            tempSpikes(:,counter,chanMap) = thisSpike;
            tempTimes(counter) = t_spike;
            tempMask(chanMap(ch),counter) = psi_masked;
            counter = counter + 1;
        end

        % now we need to align the spike waveforms since t_spike may not be an exact sample point.
        % we use cubic spline interpolation surrounding this point, 
        % extract less data before/after this point than the original data
        snips{j} = interpolate_spikes( tempSpikes,fs,(tempTimes - round( tempTimes ))/fs + preTime*2,...
                                            preTime,postTime,0.70,4 ); % 4x interpolation with slight smoothing
        sptimes{j} = (tempTimes + start(j) - 1); % to make relative to the start of "data", not "seg"
        mask{j} = tempMask;
    end

    % now concatenate the data
    snips = [snips{:}];
    sptimes = [sptimes{:}];
    mask = [mask{:}];

    % remove bad inds
    badInds = isnan( sptimes );
    snips(:,badInds,:) = [];
    sptimes(badInds) = [];
    mask(:,badInds) = [];

    % re-order by increasing spiketime
    [sptimes,idx] = sort( sptimes );
    snips = snips(:,idx,:);
    mask = mask(:,idx);

    % check for overlaps
    overlap = find_overlapping_spikes( sptimes,mask(chanMap,:),maxPts );

    % remove overlapping spikes if desired
    if removeOverlap
        snips(:,overlap,:) = [];
        sptimes(overlap) = [];
        mask(:,overlap) = [];
    end

    %% HELPER FUNCTIONS
    function p = check_inputs()
        names = {'chanMap','lowThresh','highThresh','maxPts',...
            'maxChans','artifact','removeOverlap','computeNoiseCov'};
        defaults = {1:m,2,4,inf,inf,inf,false,false};

        p = inputParser();
        for k = 1:numel( names )
            p.addParameter( names{k},defaults{k} );
        end

        p.parse( varargin{:} );
        p = p.Results;
    end
end
