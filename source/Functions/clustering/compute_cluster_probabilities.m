function [probs,priors] = compute_cluster_probabilities( X,labels,mu,sigma )
% [probs,priors] = compute_cluster_probabilities( X,label,mu,sigma )
% 
% compute the probability of point i belonging to cluster j for all points
% 1:N and clusters 1:K. Additionally, compute the prior probabilities (i.e.
% number of points belonging to a cluster) for each cluster. X is an n x m
% matrix, with n = observations, m = variables.
%
% sigma can be a K-dim cell array, or an m x m x k tensor. mu is a k x m
% matrix.

[N,M] = size( X );
if ~isrow( labels )
    labels = labels';
end
K = unique( labels );
probs = zeros( N,numel( K ) );
priors = sum( bsxfun( @minus,repmat( K,N,1 ),labels' )==0 ) / N;

% loop over clusters
for k = 1:numel( K )
    if iscell( sigma )
        s = sigma{k};
    else
        s = sigma(:,:,k);
    end
    %U = chol( s );
    %Q = bsxfun( @minus,X,mu(k,:) ) / U;    
    xhat = bsxfun( @minus,X,mu(k,:) );
    probs(:,k) = priors(k) / ((2*pi)^(M/2) * det(s)^0.5) * exp( -0.5*sum( xhat*pinv(s).*xhat,2 ) ); % gaussian pdf
end

% normalize the probabilities
%probs = bsxfun( @rdivide,probs,sum( probs,2 ) );

end