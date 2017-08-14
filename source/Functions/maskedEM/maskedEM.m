function [label, model, likelihood, R] = maskedEM_GMM( X,mask,virtualX,varargin )
% [label, model, likelihood, R] = maskedEM( X,mask,virtualX,(init,searchForK) ) 
%
% Perform EM algorithm on a masked data matrix X for fitting a 
% Gaussian Mixture virtualX of unkown components.
%
%                   >>> INPUTS >>> 
%   X: 
%       n x d feature matrix, with n = observations, d = variables
%   mask: 
%       n x d sparse mask matrix, with each element (i,j) in the set [0,1]
%   virtualX: 
%       structure from "compute_noise_ensemble()"
%   (init):
%       the initial cluster number. Default = 1
%   (searchForK):
%       boolean flag (1 or 0) indicating on whether to search for the optimal
%       cluster number or not. Default = 0
%
%                   <<< OUTPUTS <<<
%   label: 
%       n x 1 cluster label for each ith data point
%   model: 
%       trained model structure with fields:
%           mu - d x k matrix of means of each cluster
%           sigma - d x d x k tensor of covariances of each cluster
%           w - 1 x k vector weights of each cluster
%   likelihood: 
%       1 x maxiter vector of the log likelihood for all iterations
%   R: 
%       n x K matrix of probabilities, for all clusters 1:K
%   penalty:
%       1 x maxiter of the modified BIC penalty function for all iterations
% 
% original mixGaussEm by Mo Chen
%
% adapted for masked EM by JMS, 8/8/2017

% initialize the parameters
tol = 1e-6;
err = inf
maxiter = 500;
[n,d] = size( X );

% input checking
if nargin > 3 && ~isempty( varargin{1} )
    init = varargin{1};
else 
    init = 1;
end
if nargin > 4 && ~isempty( varargin{2} )
    searchForK = varargin{2};
else
    searchForK = false;
end

% error checking
if size( mask ) ~= [n,d]
    error( 'Mask and data matrices not equal sizes' );
end
if init > d
    warning( ['Number of requested clusters is greater than dimension of X',...
        'Reducing number of clusters to equal size( X,2 )-1'] );
    init = d;
end

% check if searching for cluster number
if searchForK
    Kmax = 20;
    penalty = inf( 1,Kmax-init + 1 );
else
    Kmax = init;
end

% Loop over possible # clusters
% ========================================
counter = 0;
for k = init:KMax
    counter = counter + 1;

    % run the EM algorithm for k gaussian components
    [model,R,likelihood] = run();
    [~,label] = max( R,2 );

    % compute the penalty function
    penalty(counter) = penalization();
end
% ========================================

if searchForK
    [~,k] = min( penalty ); % minimum penalty
    fprintf( 'Optimal # of clusters: %i',k );
    [model,R,likelihood] = run();
    [~,label] = max( R,2 );
end


%% FUNCTIONS
function [model,R,likelihood] = run()
    % Main runtime 
    likelihood = -inf( 1,maxiter );
    iter = 1;
    
    % initialize random labels
    [label,R] = initialization( X,k ); % our first "expectation"

    while iter < maxiter && err > tol
        iter = iter + 1;
        nK = sum( R ) / n; % the prior probabilities of each jth cluster

        % compute the model given current cluster assignments 
        model = maximization( X,virtualX,mask,label,k );
        model.w = nK;

        % now compute the posterior probabilities given the cluster model
        [R,likelihood(iter)] = expectation( virtualX,R,model );

        % compute change in likelihood
        err = abs( likelihood(iter) - likelihood(iter-1) ) / abs( likelihood(iter) );
    end
end


function [label,R] = initialization( X,k )
    label = ceil( k*rand( n,1 ) );
    R = full( sparse( 1:n,label,1,n,k,n ) );
end

function [R,likelihood] = expectation( virtualX,R,model )
    % computes new posterior probabilities for each point given the current model
    for i = 1:k
        R(:,i) = loggausspdf( virtualX.y,model.mu(:,i),model.Sigma(:,:,i) );
    end
    R = bsxfun( @plus,R,log( model.w ) ); % log( P(Y=k | x,theta) ) + log( P(Y==k) ) = log{ P(Y==k | x,theta) * P(Y==k) }
    T = logsumexp( R,2 );
    likelihood = sum( T ) / n; % SUM_1:K{ SUM_1:N{ log{ P(Y=k | x,theta) * P(Y == k) } } }
    R = exp( bsxfun( @minus,R,T ) ); % eliminate the log
end

function mdoel = maximization( X,virtualX,mask,label,K )    
    % Computes the current model of the clusters 1:K given the data points for each
    mu = zeros( d,K );
    Sigma = zeros( d,d,K );
    for i = 1:K
        C = (label == i);                   % all points belonging to this cluster
        L = (repmat( C,1,d ) & mask > 0);   % unmasked points belonging to this cluster 
        M = (repmat( C,1,d ) & mask == 0);  % masked points belonging to this cluster

        % compute cluster mean and covariance
        mu(:,i) = sum( bsxfun( @plus,virtualX.y(L),sum( M ) .* virtualX.mu ) ) / nnz( C ); % the cluster mean ( sum of means of unmasked features + noise mean * # masked features )
        yhat = bsxfun( @minus,virtualX.y(L),mu(:,i) ); % de-meaned, unmasked points
        noiseCov = (sum( M ) / nnz( C )) * ( (virtualX.mu' - mu(:,i)) * (virtualX.mu' - mu(:,i))' ); % covariance of the noise for this cluster
        Sigma(:,:,i) = (yhat * yhat' + noiseCov... % sum of covariances of masked and unmasked data
                        + diag( sum( bsxfun( @plus,virtualX.eta(L),sum( M ) * virtualX.sigma ) ) ))... % variances of masked data
                        / nnz( C ); 
    end
    model = struct;
    mdoel.mu = mu;
    mdoel.Sigma = Sigma;
end

function y = loggausspdf(X, mu, Sigma)
    X = bsxfun( @minus,X,mu );
    y = -0.5 * (d*log( 2*pi ) + log( det( Sigma ) )...
        + X' * inv( Sigma ) * X...
        + sum( bsxfun( @times,virtualX.eta,diag( inv( Sigma ) ) ),2 ));   % normalization constant
end

function penalty = penalization()
    % calculates the penalty using a modified model complexity
    freeParams = 0;
    for i = unique( label )'
        C = (label == i);
        L = (repmat( C,1,d) & mask > 0);

        % compute free parameters for this cluster
        r = sum( L,2 ); % # of unmasked features for each data point in this cluster
        F = sum( (r .* r+1) / 2 + r + 1 ) / nnz( C );
        freeParams = freeParams + F;
    end
    % bayesian information criterion
    penalty = freeParams*log( n ) - 2*log( max( likelihood ) );
end

function splitClusters()

end

function mergeClusters()

end

function rmClusters()

end


end