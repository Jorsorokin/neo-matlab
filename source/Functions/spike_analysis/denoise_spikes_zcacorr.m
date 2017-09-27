function [W,filteredSpikes] = denoise_spikes_zcacorr( spikes,mask )
% [W,filteredSpikes] = denoise_spikes_zcacorr( spikes,mask )
%
% de-noises the spikes in the n x m x c spike matrix via noise-covariance 
% estimation and ZCA-corr whitening
%
% Inputs:
%   spikes: m x n x c spike matrix, with
%           m - # of samples
%           n - # of spikes
%           c - # of channels
%   
%   mask: c x n mask matrix with 0 <= (i,j) <= 1. 
% 
% Outputs:
%   W - m x m x c noise whitening matrix
%
%   filteredSpikes - spikes following denoising from W
%
% each ith spike for each jth channel is included in the whitening matrix 
% estimation if mask(j,i) == 0. Whitening matrix W is then created as:
%   
%       W_ch = P^(-1/2) * V^(-1/2)
%
%   where P is the correlation matrix of all masked spikes of channel j
%   and V is a diagonal matrix of the variances of each k -> [1:m] point of
%   masked spikes from channel j
%
% the final projection to obtain filtered spikes is then computed as:
%   
%       filteredSpikes_ch = W_ch * spikes(:,:,ch)
% 
% the resulting transformed spikes are as maximally correlated to the
% originals as possible (<-- perhaps don't want this?)
%
% by JMS, 9/1/17

% check inputs
[m,n,c] = size( spikes );
assert( all( size( mask ) == [c,n] ) );

% preallocate
W = zeros( m,m,c );
P = zeros( m,m );
V = zeros( m,m );
if nargout > 1
    filteredSpikes = spikes;
end

% loop over channels
for j = 1:c
    noise = mask(j,:)==0;
    P(:) = corr( spikes(:,noise,j)' ); % correlation matrix
    V(:) = diag( var( spikes(:,noise,j),[],2 ) ); % diagonal matrix
    W(:,:,j) = P^(-0.5) * V^(-0.5); % whitening matrix for channel j
    
    if nargout > 1
        filteredSpikes(:,:,j) = W(:,:,j)' * spikes(:,:,j); % de-noised 
    end
end

end
    
    