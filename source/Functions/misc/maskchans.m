function X = maskchans( X,mask )
    % X = maskchans( X,mask )
    %
    % mask channels in X by a masking matrix "mask"
    %
    % Inputs:
    %   X - n x m x c matrix, with c = channels
    %   
    %   mask - c x n masking matrix
    %
    % Outputs:
    %   X - n x m x c, with values repalced with 0 if mask(c,i) < 0
    
    c = size( X,3 );
    if size( mask,1 ) ~= c
        error( 'mask and X must have equal number of channels' );
    end
    
    for ch = 1:c
        X(:,:,ch) = bsxfun( @times,X(:,:,ch),mask(ch,:)' ); % weights each ith spike in each jth channel by the mask value (i,j)
    end
end