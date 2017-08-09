function model = compute_noise_emsemble( features,mask )
% model = compute_noise_emsemble( features,mask )
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
%       noiseDist:
%           1 x p vector of gaussian distributions: N( mu(i),sigma(i) )
%       x:
%           n x p matrix with (i,j) in the set [0,1], defining whether feature (i,j) should be replaced by
%           a noise ensemble N( mu(j),sigma(j) )
%       y:
%           n x p mean expectation matrix: mask(i,j) * features(i,j) + (1-mask(i,j)) * mu(j)
%       z:
%           n x p variance expectation matrix: mask(i,j) * features(i,j)^2 + (1-mask(i,j)) * (mu(j)^2 + sigma(j)^2)
%       eta:
%           n x p variance matrix: z(i,j) - y(i,j)^2; 
% 
%
% - part of the masked-EM algorithm. See Kadir et al. 2015 for details
%
% by JMS, 8/8/2017

% compute the various model variables
[n,p] = size( features );
isMasked = (mask == 0);
model = struct;
model.mu = nanmean( features(isMasked),[],1 );
model.sigma = nanvar( features(isMasked),[],1 );
model.noiseDist = gaussify( model.mu,model.sigma );
model.y = (mask .* features) + bsxfun( @times,(1-mask),model.mu );
model.z = (mask .* features.^2) + bsxfun( @times,(1-mask),(model.mu.^2 + model.sigma.^2) );
model.eta = model.z - model.y.^2;

% find which elements of the feature matrix should be replaced with the noise ensemble
noiseDistMat = repmat( model.noiseDist,n,1 );
model.x = rand( n,n ) > mask;

end


    