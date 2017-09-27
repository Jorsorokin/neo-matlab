function [start,stop] = find_blockIndex_range( n,m,maxSize )
% [start,stop] = find_indexing_range( n,m,(maxSize) )
%
% find the indexing ranges to perform block looping over points "1:m" for 
% nearest-neighbor, etc. This greatly improves performance (not looping over 
% all points) while avoiding massive (> 1 GB) matrices. Default "maxSize",
% which refers to the maximum matrix size allowed, is 100,000,000 entries.
%
% "n" refers to number of elements in first vector, "m" refers to number of
% elements in second. This is more flexible than simply providing one
% number, as neighborhood graphs etc. may not necessarily be square.

if nargin < 3
    maxSize = 100000000;
end

% get max Pts possible
maxPts = floor( maxSize / n );
if maxPts < m
    remainder = mod( m,maxPts );
    tempStart = 1:maxPts:m-remainder;
    tempStop = maxPts:maxPts:m-remainder;
    if remainder > 0
        start = [tempStart,tempStart(end)+remainder];
        stop = [tempStop,m];
    else
        start = tempStart;
        stop = tempStop;
    end
else
    start = 1;
    stop = m;
end

end