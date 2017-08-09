function [snips,sptm,spamp] = get_spike_snips( sptm,data,PREPTS,POSTPTS,artifact )
    % [snips,sptm,spamp] = get_spike_snips( sptm,data,PREPTS,POSTPTS,artifact )
    % 
    % pulls out a window of data surrounding each sample in "sptm"
    % and eliminates any peak with a positive voltage       
    nChan = size( data,2 );

    % get the snips
    snips = zeros( PREPTS+POSTPTS,numel( sptm ),nChan );
    spamp = zeros( numel( sptm ),nChan );
    for j = 1:numel( sptm )
        snips(:,j,:) = data(sptm(j)-PREPTS:sptm(j)+POSTPTS-1,:);
        spamp(j,:) = data(sptm(j),:);
    end

    % get rid of spikes with positive peak or too large an amplitude
    maxPk = max(spamp,[],2);
    bad = maxPk>0 | mean( abs( spamp ),2 ) >= artifact;
    snips(:,bad,:) = [];
    sptm(bad) = []; 
    spamp(bad,:) = [];       
end