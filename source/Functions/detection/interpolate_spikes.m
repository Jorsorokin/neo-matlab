function snips = interpolate_spikes( spikes,fs,sptm,preTime,postTime,varargin )
    % snips = interpolate_spikes( spikes,fs,sptm,preTime,postTime,(smoothing,resampValue) )
    %
    % interpolates the voltage waveforms contained in "spikes" according to the 
    % spike time specified by "sptm" using cubic spline interpolation.
    % "sptm" does not necessarily refer to the peak of the waveform, 
    % but could be the center of mass, for instance.
    %
    % preTime and postTime are scalars specifying the amount of data (in time)
    % to keep before/after the aligned spike. This should be less than or equal to the 
    % total amount of data of the spike waveform.
    %
    % note that sptm should be provided in TIME
    
    % check size & inputs
    smoothing = [];
    resampValue = 1; % default no upsampling
    if nargin > 4 
        smoothing = varargin{1};
    end
    if nargin > 5 && ~isempty( varargin{2} )
        resampValue = varargin{2};
    end
        
    [n,nSp,nChan] = size( spikes );
    prePts = floor( preTime * fs );
    postPts = floor( postTime * fs );
    spPts = sptm * fs;
    totalPts = prePts + postPts;
    if totalPts > n
        error( 'amount of requested time surrounding spike is larger than supplied data' );
    end
    
    % spline interpolation
    warning( 'off' );
    totalPts = floor( totalPts * resampValue );
    snips = nan( totalPts,nSp,nChan );
    xind = 1:n;
    for sp = 1:nSp
        try
            x = linspace( spPts(sp)-prePts,spPts(sp)+postPts-1,totalPts );
            snips(:,sp,:) = csaps( xind,squeeze( spikes(:,sp,:) )',smoothing,x )'; % smooths the interpolant
            
        catch
        end
    end
    warning( 'on' );
end
