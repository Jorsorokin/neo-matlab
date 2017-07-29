function mteSignals = MTEO( X,varargin )
% mteSignals = MTEO( data, (kStart,kMax) )
%
% Computes the multi-teager energy of the signals in "X". The Teager
% energy of a vector x is defined as:
%
%           te(n) = x(n)^2 - [(x-k)*(x+k)]
%
%       where k is determined a-priori. 
%
% This results in a filtering of the original signal x to highlight
% high-frequency pulses, however it requires knowledge on the optimal value
% of k for a particular signal. 
%
% The multi-teager energy uses a range of values of "k" and then computes 
% the maxima across the filtered signals for each sample point. This
% results in a combination of the original signals to better detect
% impulses that may have varying widths/frequencies. 
%
% In addition, this function applies a smoothing hamming window of length
% 4k+1 to each signal to eliminate noise.
% 
% Refer to ___ for details
%
%               >>> INPUTS >>>
% X:
%   an N x M matrix of original signals, oriented column-major
% (kStart):
%   scalar specifying the starting value of k (default = 1)
% (kMax):
%   scalar determining the maximum value of k (default = 6)
%
%               <<< OUTPUTS <<<
% mteSignals:
%   an (N-2*kmax) x M matrix of filtered signals in X. 
%
% By JMS, 5/15/2017

% check inputs
kStart = 1;
kMax = 6;
if nargin == 2 
    kStart = varargin{1};
elseif nargin == 3
    kMax = varargin{2};
end

% get size of data matrix X
[N,M] = size( X );

% pre-allocate 
mteSignals = zeros( N,M ); % the first kMax points will = 0, useful for alignment
hWin = zeros(1,4*kMax+1);
h = zeros(kMax,1);
ind = 0;
for k = kStart:kMax
    ind = ind+1;
    L = 4*k + 1;
    
    % create our hamming window
    for i = 1:L
        hWin(i) = 0.54 - 0.46*cos( 2*pi*(i-1)/(L-1) );
    end
    
    % since we only use the middle point, just store this
    h(ind) = hWin(2*k+1) / norm( hWin ); 
end

% loop over the sample points
for n = kMax+1:N-kMax
    
    % create the filtered points from n - kmin : n+kmax
    % This is the vectorized version that avoids inner loops 
    newPoints = bsxfun( @minus,X(n,:).^2,...
        X(n-kMax:n-kStart,:).*X(n+kMax:-1:n+kStart,:) );

    % compute the max for each channel and convolve with the appropriate
    % hamming window (depending on the median "k" between channels)
    mteSignals(n,:) = max( bsxfun( @times,newPoints,h ) ); % "convolve" with h
end


end
    
    
    
    