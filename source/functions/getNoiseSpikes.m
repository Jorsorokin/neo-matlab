function noiseIDX = getNoiseSpikes( spikes,fs )
%
% noiseIDX = getNoiseSpikes( spikes,fs )
%
% finds indices in an n x m matrix of spike waveforms that are putative
% noise, accidentally picked up by a spike-detecting software

% get the amplitudes
[vMin,vMax,tMin] = spikeHeight( spikes,fs );
amp = abs( vMin ) + vMax;

% get the slopes
dvRatio = abs( (max( diff( spikes ) ) ./ min( diff( spikes ) ) ) );

% get the half width
hw = halfWidth( spikes,tMin,fs );

% get the noise indices
noiseIDX = find( dvRatio >= 1 & nanmean( hw,2 )' >= .002 ...
    & (amp <= 200 | amp >= 800));

end