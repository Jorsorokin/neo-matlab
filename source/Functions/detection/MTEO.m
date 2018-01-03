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
if nargin >= 2 && ~isempty( varargin{1} )
    kStart = varargin{1};
end
if nargin >= 3 && ~isempty( varargin{2} )
    kMax = varargin{2};
end

% pre-allocate 
[N,M] = size( X );
window = -kMax:-kStart;
windowSize = numel( window );
mteSignals = zeros( N,M ); % the first kMax points will = 0, useful for alignment
h = zeros(1,numel( window ));

% create hamming window
ind = 0;
for k = kStart:kMax
    ind = ind+1;
    L = 4*k + 1;
    hWin = 0.54 - 0.46*cos( 2*pi*([1:L]-1)/(L-1) );
    
    % since we only use the middle point, just store this
    h(ind) = hWin(2*k+1) / norm( hWin ); 
    clear hWin
end

% create our indexing matrices
n = (kMax+1:N-kMax)';
nPts = numel( n );
preInd = bsxfun( @plus,n,repmat( window,nPts,1 ) );
postInd = bsxfun( @plus,n,repmat( -window,nPts,1 ) );

% get the pre/post points surrounding each point "n" for each column of X
for m = 1:M
    prePts = reshape( X(preInd,m),nPts,windowSize );
    postPts = reshape( X(postInd,m),nPts,windowSize );
    mteSignals(n,m) = max( bsxfun( @times,...
        bsxfun( @minus,X(n,m).^2,prePts.*postPts ),h ),[],2 );
end
   
end
    
    
    
    
