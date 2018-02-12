function E_spikeXY = get_spike_XY( X,chanLoc,varargin )
% E_spikeXY = get_spike_XY( X,chanLoc,(mask) )
%
% compute the expectation of the (x,y) location of each spike given the
% channel distance vector "chanLoc". X is an n x m x c matrix, with n =
% observations, m = variables, c = channels. Can optionally provide a c x n
% masking matrix "mask", in which E_spikeXY will take only unmasked channels for
% each spike into account

if nargin > 2 && ~isempty( varargin{1} )
    mask = varargin{1};
    if size( mask,2 ) ~= size( X,1 )
        error( 'number of points in the mask matrix and X do not match' );
    end
end

% find the minimum voltage for each spike/channel. Then normalize by
% the minimum across channels to interpolate between [0,1]
amp = squeeze( min( X,[],2 ) );
[chanAmp,bestChan] = min( amp,[],2 );
%amp = gather( bsxfun( @rdivide,amp,chanAmp ) ); % provides a smoother interpolation compared to the original masking 
if exist( 'mask','var' ) 
    % take the E[(X,Y)_spike_i] - the expecation of the (x,y) location of the ith spike
    % by weighting the channel distances by the mask values for each
    % channel and dividing by the number of channels with mask > 0
    
    
    % JS, 10/18/2017
    % Removed the "smooth interpolation" by the amplitudes to reduce
    % infinitely dense regions in the projection space and spread apart 
    % spurious activity off of the main spike 
    %E_spikeXY = bsxfun( @rdivide,chanLoc' * (amp' .* mask>0),sum( mask>0 ) )' - 1; 
    E_spikeXY = (chanLoc' * mask)'; % <-- new code
else
    % if no mask, take the expectation by simply using the channel with
    % the largest voltage deflection
    E_spikeXY = chanLoc(bestChan,:) - 1;
end


end
