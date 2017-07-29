function [L,K] = findClustNum( X,k0,kMax )
% [L,K] = findClustNum( X,k0,kMax )
%
% Searches through cluster numbers k0 -> kMax and discovers the optimal
% cluster number for the data X via kmeans and projections of the data onto
% inter-cluster vectors followed by the adtest

K = k0 + 1;
Khat = K;
fprintf( 'Searching for optimal cluster #' );
while K <= kMax
    fprintf( ' . ' );
    
    % run k-means to cluster the data
    [labels,mu] = mixGaussEm( X,K );
    mu = mu.mu; 
    
    % loop over the cluster numbers, compute the distribution of
    % projections of clusters onto cluster-to-cluster vectors. 
    for i = 1:size( mu,2 )
        % get samples X_i (those in cluster i)
        X_i = X(:,labels==i);

        for j = 1:i-1   
            % get samples X_j
            X_j = X(:,labels==j);

            % project onto vector mu_i - mu_j
            muVec = mu(:,i) - mu(:,j);
            Xproj_i = muVec' * X_i;
            Xproj_j = muVec' * X_j;

            % compute the histograms of each X_i, X_j, and the two combined
            % and the anderson-darling test statistic
            if i ~= j
                [~,~,adstat_i] = adtest( Xproj_i );
                [~,~,adstat_j] = adtest( Xproj_j );
                [~,~,adstat_ij] = adtest( [Xproj_i,Xproj_j] );

                % compare the statistics. If any individual X_i or X_j has a
                % smaller statistic than the combined X_ij, then we have too
                % many clusters
                if (adstat_i > adstat_ij) && (adstat_j > adstat_ij)
                    Khat = Khat - 1;
                end    
            end
        end
    end
    % now check if Khat < K. If so, we have surpased the number of
    % clusters that actually exist in our data
    if Khat < K
        break;
    else
        K = K+1;
        Khat = K;
    end
end
fprintf( '\nOptimal # of clusters: %i\n',Khat );
k0 = Khat;
L = kmeans( X',k0 );

end
