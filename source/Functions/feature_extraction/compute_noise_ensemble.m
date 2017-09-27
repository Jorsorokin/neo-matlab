function model = compute_noise_ensemble( features,mask )
% model = compute_noise_ensemble( features,mask )
% 
% computes the mean (mu) and variance (sigma) of all instances
% of masked features, for each feature separately, and replace 
% features(i,j) by a gaussian distribution of all masked features(:,j) 
% with P( 1-mask(i,j) )
%
%                   >>> INPUTS >>>
% features:
%   an n x p matrix, with n = observations, p = variables
%
% mask:
%   an n x p matrix, with each (i,j) in the set [0,1]
%
%                   <<< OUTPUTS <<<
% model:
%   structure with the following fields:
%       mu:
%           1 x p vector of masked-feature means
%       sigma:
%           1 x p vector of masked-feature variances
%       y:
%           n x p mean expectation matrix of x: mask(i,j) * features(i,j) + (1-mask(i,j)) * mu(j)
%       eta:
%           n x p variance matrix of x: z(i,j) - y(i,j)^2; 
% 
%
% - part of the masked-EM algorithm. See Kadir et al. 2015 for details
%
% by JMS, 8/8/2017

% compute the various model variables
[n,p] = size( features );
model = struct;
isMasked = (mask == 0);
nMasked = sum( isMasked );

model.mu = bsxfun( @rdivide,sum( isMasked .* features ),nMasked );
model.sigma = bsxfun( @rdivide,sum( isMasked .* bsxfun( @minus,features,model.mu ).^2 ),nMasked );
model.y = mask .* features + bsxfun( @times,(1-mask),model.mu );
z = mask .* features.^2 + bsxfun( @times,(1-mask),(model.mu.^2 + model.sigma) );
model.eta = z - model.y.^2;

% loop over each each feature dimension and determine if we should replace points in "ensemble_features"
% for i = 1:p
%     replace = rand( n,1 ) > mask(:,i);
%     if any( replace )
%         ensemble_features(replace,i) = model.mu(i) + randn( nnz( replace ),1 )*sqrt( model.sigma(i) );
%     end
% end

end


    