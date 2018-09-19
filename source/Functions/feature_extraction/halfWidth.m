function hw = halfWidth( spikes,peaks,fs )
% hw = halfWidth( spikes,peaks,fs )
%
% calculates the half-width of columns in the data matrix "spikes".
% default is to find negative-deflection half widths. Thus, make sure
% "peaks" and your "spikes" variables reflect downward-facing spikes (ie as
% recorded extracellularly). 
%
% Inputs:
%   spikes:
%       n x m data matrix, n = points, m = spikes
%   peaks:
%       m x 1 vector of spike amplitude (defined as the minimum voltage of
%       the peak)
%   fs:
%       scalar sampling rate
%
% Outputs:
%   hw:
%       m x 1 vector of half-width info for each spike

% sizes and pre-allocation
[nPts,nSp] = size( spikes );
hw = nan( nSp,1 );  
x = linspace(1,nPts,nPts);
xx = linspace(1,nPts,nPts*5);
fsUP = fs*5;
spikeUP = interp1( x,spikes,xx,'spline' );
if nSp == 1
    spikeUP = spikeUP';
end

for s = 1:nSp
    try
        % extract the spike
        sp = spikeUP(:,s);

        % get the points just under the half-point
        first = find( sp <= peaks(s)/2,1 ) / fsUP;
        last = find( sp <= peaks(s)/2,1,'last' ) / fsUP;

        % calculate the distance in time
        hw(s) = last - first;
    catch
    end
end
    