function [xcg,lags] = correlogram( train1,train2,fs,varargin )
%
% [xcg,lags] = correlogram( train1,train2,fs,(maxlag) )
%
% compute the cross-correlogram between two spike trains (of equal length)


% check inputs
N = numel(train1);
M = numel(train2);
if N ~= M
    error( 'spike trains durations must be equal' );
end

if nargin > 3 && ~isempty( varargin{1} )
    maxlag = round( varargin{1} * fs );
else
    maxlag = N;
end

% compute the cross-correlation
[xcg,lags] = xcorr( train1,train2,maxlag );

% compute mean firing rates
T = N;
fr1 = sum( train1 ) / (T / fs);
fr2 = sum( train2 ) / (T / fs);

% compute the cross-correlogram
lags = lags / fs;
O = maxlag / fs - abs(lags);
xcg = xcg ./ (O * sqrt(fr1*fr2));

end
    