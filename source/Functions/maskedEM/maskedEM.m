function [label, model, llh, R] = maskedEM( X,mask,model ) 
% Perform EM algorithm on a masked data matrix X for fitting a 
% Gaussian Mixture Model of unkown components.
%
%                   >>> INPUTS >>> 
%   X: d x n data matrix
%   mask: d x n mask matrix
%   model: structure from "compute_noise_ensemble()"
%
%                   <<< OUTPUTS <<<
%   label: 1 x n cluster label
%   sortModel: trained model structure
%   llh: loglikelihood
%   R: the probabilities
% 
% original mixGaussEm by Mo Chen
%
% adapted for masked EM by JMS, 8/8/2017

% initialize the parameters
tol = 1e-6;
maxiter = 500;
llh = -inf( 1,maxiter );
[d,n] = size( X );
[label,R] = initialization( X,init );
[~,k] = size( R,2 );

for iter = 2:maxiter
    [~,label(1,:)] = max(R,[],2);
    R = R(:,unique(label));   % remove empty clusters
    model = maximization(X,R);
    [R, llh(iter)] = expectation(X,model);
    if abs(llh(iter)-llh(iter-1)) < tol*abs(llh(iter)); break; end;
end
llh = llh(2:iter);

function [label,R] = initialization( X,init )
    label = ceil( init*rand(1,n) );
    R = full( sparse( 1:n,label,1,n,init,n ) );
end

function [R, llh] = expectation()
    R(:) = 0;
    for i = 1:k
        R(:,i) = loggausspdf( X,mu(:,i),Sigma(:,:,i) );
    end
    R = bsxfun(@plus,R,log(w));
    T = logsumexp(R,2);
    llh = sum(T)/n; % loglikelihood
    R = exp(bsxfun(@minus,R,T));
end

function sortModel = maximization()
    [d,n] = size( X );
    nk = sum( R,1 );
    w = nk / n; % the probability weights
    
    % compute mu/sigma for the clusters
    mu = zeros( d,k );
    Sigma = zeros( d,d,k );
    r = sqrt( R );
    for i = 1:k
        thisLabel = (label == k);
        L = (repmat( thisLabel,d,1) & mask > 0);
        M = (repmat( thisLabel,d,1 ) & mask == 0);
        mu(:,i) = sum( bsxfun( @plus,model.y(L),sum( M,2 )' .* model.mu ) ) / nnz( thisLabel );
        yhat = bsxfun( @minus,model.y(L),mu(:,i) );
        noiseCov = (sum( M,2 )' / nnz( thisLabel )) * ( (model.mu' - mu(:,i)) * (model.mu' - mu(:,i)) );
        sigma(:,:,i) = (yhat * yhat' + noiseCov... % covariances of masked and unmasked data
                        + diag( sum( bsxfun( @plus,model.eta(L),sum( M,2 )*model.sigma ) ) ))... % variances of masked data
                        / nnz( thisLabel ); 
    end
    sortModel.mu = mu;
    sortModel.Sigma = Sigma;
    sortModel.w = w;
end

function y = loggausspdf(X, mu, Sigma)
    X = bsxfun( @minus,X,mu );
    [U,p]= chol( Sigma );
    if p ~= 0
        error('ERROR: Sigma is not PD.');
    end
    Q = U'\X;
    q = dot(Q,Q,1);  % quadratic term (M distance)
    c = d*log(2*pi)+2*sum(log(diag(U)));   % normalization constant
    y = -(c+q)/2;
end

end