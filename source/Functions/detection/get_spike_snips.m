function [snips,sptm,spamp] = get_spike_snips( sptm,data,PREPTS,POSTPTS )
    % [snips,sptm,spamp] = get_spike_snips( sptm,data,PREPTS,POSTPTS )
    % 
    % pulls out a window of data surrounding each sample in "sptm". "data"
    % should be a m x c matrix, with m = points, c = channels
    %
    % sptm = n x p x c tensor, with n = points, p = spikes, c = channels
    % sptm = p x 1 vector
    % spamp = p x c matrix
       
    nChan = size( data,2 );

    % get the snips
    if ~isrow( sptm )
        sptm = sptm';
    end
    totalPts = POSTPTS+PREPTS;
    nSpikes = numel( sptm );
    inds = bsxfun( @plus,sptm-PREPTS,repmat([0:totalPts-1]',1,nSpikes) );
    snips = reshape( data(inds,:),totalPts,nSpikes,nChan ); % makes "snips" into an nPts x nSpikes x nChan matrix
    spamp = data(sptm,:);
    clear inds   
end