function [vMin,vMax,tMin,tMax] = spikeHeight( spikes,fs,varargin )
% [spMin,spMax,tMin,tMax] = spikeHeight( spikes,fs,start,stop )
%
% calculate the minimum/maximum values of each column in the "spikes" data
% matrix, as well as the times of onsets. "start" & "stop" should be supplied in
% time, not samples, and determines when to begin * end searching for the minimum
% & maximum in the voltage waveform

start = 1;
stop = size( spikes,1 );
if nargin >= 3
    start = round( varargin{1} * fs );
end
if nargin == 4
    stop = round( varargin{2} * fs );
end

[vMin,tMin] = min( spikes(start:stop,:,:),[],1 );
med = round( mean( squeeze( median( tMin ) ) ) );
[vMax,tMax] = max( spikes(start+med:stop,:,:),[],1 );
tMin = (tMin + start - 1) / fs;
tMax = (tMax + start + med - 1) / fs;

% clean up
if ndims( spikes ) == 3
    vMin = squeeze( vMin );
    vMax = squeeze( vMax );
    tMin = squeeze( tMin );
    tMax = squeeze( tMax );
else
    vMin = vMin';
    vMax = vMax';
    tMin = tMin';
    tMax = tMax';
end