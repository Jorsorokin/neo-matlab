function snips = get_aligned_spikes( spikes,sptm,fs,prePts,postPts )
    % snips = get_aligned_spikes( spikes,sptm,fs,preTime,postTime )
    %
    % align the voltage waveforms contained in "spikes" according to the 
    % spike time specified by "sptm". "sptm" does not necessarily refer to 
    % the peak of the waveform, but could be the center of mass, for instance.
    %
    % preTime and postTime are scalars specifying the amount of data (in samples)
    % to keep before/after the aligned spike. This should be less than or equal to the 
    % total amount of data of the spike waveform.
    %
    % note that sptm should be provided in SAMPLES not time. 
    
    % check size & inputs
    [n,nSp,nChan] = size( spikes );
    uFs = fs * 4;
    dt = 1/fs;
    udt = 1/uFs;
    totalPts = prePts + postPts;
    if totalPts > n
        error( 'amount of requested time surrounding spike is larger than supplied data' );
    end
    
    warning( 'off' );
    % resample the spikes for the alignment / spline interpolation
    uPrePts = prePts * 4;
    uPostPts = postPts * 4;
    snips = nan( totalPts,nSp,nChan );
    xind = 1:n;
    for sp = 1:nSp
        try
            x = linspace( sptm(sp)-prePts,sptm(sp)+postPts-1,totalPts );
            xx = linspace( x(1),x(end),totalPts*4 );
            upspikes = csaps( xind,squeeze( spikes(:,sp,:) )',0.3,xx )'; % smooths the interpolant
            
            % now down-sample back to original sampling rate & extract pre/post time around the spike
            snips(:,sp,:) = spline( xx,upspikes',x )';
        catch
        end
    end
    
    warning( 'on' );
end