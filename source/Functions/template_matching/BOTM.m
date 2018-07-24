function [labels,D] = BOTM( X,eta,prior )
    % [labels,P] = BOTM( X,eta )
    %
    % performs bayesian optimal template matching by matching the templates
    % in "eta" against the raw data in X
    %
    % Inputs:
    %   X - an n x m matrix of spike waveforms, n = points, m = spikes
    %
    %   eta - an n x k matrix of spike templates
    %
    %   prior - a 1 x k vector of prior probabilities for the templates
    %
    % Outputs:
    %   labels - a 1 x m vector of class labels
    %
    %   P - a k x m matrix of class probabilities
    
    
    [n,m] = size( X );
    k = size( eta,2 );
    %P = zeros( k,m );
    D = zeros( m,k );
    
    for i = 1:k
        thisTemplate = eta(:,i);
        D(:,i) = (X' * thisTemplate) + (thisTemplate'*thisTemplate);% + log( prior(i) );
    end
    
    [~,labels] = max( D,[],2 );
end
        