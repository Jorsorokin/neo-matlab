function [W,X] = whiten_data( X )
    % [W,X] = whiten_data( X )
    %
    % Whitens the data in the N x M matrix X by using the 
    % spatial-covariance computed directly from X, or by using a
    % pre-computed whitening matrix (W...optional input), where 
    % the # of rows / columns of W equals M
    
    [V,E] = eig( cov( X ) + eye(size(X,2))*1e-10 );
    W = V * (E^-0.5 * V'); % ZCA whitening
    
    % whiten the data
    if nargout > 1
        X = X * W';
    end
end