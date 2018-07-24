function overlap = find_overlapping_spikes( sptimes,mask,maxPts )
    % overlap = find_overlapping_spikes( sptime,mask,maxPts )
    %
    % finds detected spikes that are likely due to overlapping single units
    overlap = (sptimes(2:end) - sptimes(1:end-1) < maxPts);
    if any( overlap )
        overlap = [overlap(1),overlap];
        overlap(2:end) = overlap(2:end) & any( mask(:,2:end).*mask(:,1:end-1) );
        overlap(1) = overlap(2);
    end
end