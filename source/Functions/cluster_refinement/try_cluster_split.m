function [labels,chisqStat,pval] = try_cluster_split( X )
    % [labels,chisqStat,pval] = try_cluster_split( X )
    %
    % Attempts to split the data in X into two clusters by modeling the
    % data as a chi-squared-distributed random variable with 
    % nDim-1 DOF. 
    %
    % X is first projected onto an nDim subspace using PCA,
    % which transforms the data to fit the assumptions of a 
    % chi-squared distribution (normally distributed,
    % independent variables). 
    %
    % Following the projection, the data is then modeled as a 
    % chi-squared variable with nDim-1 DOF as:
    %
    %       Z(i) = SUM_j{ (X(i,j) - X_hat(j))^2 }
    %               
    %   where X_hat(j) is the mean of dimension j across points.
    %
    % Under the null hypothesis, each Z(i) is a chi-square
    % random variable, and Z is distributed according to an 
    % nDim-1 chi-squared distribution. Then, a chi-square test statistic
    % can be implemented as:
    %
    %       chisq = SUM_i{ (O(i) - E(i))^2 / E(i) }
    %
    %   where O(i) and E(i) are observed and expected counts in the 
    %   empirical and theoretical CDFs of Z and an nDim-1 chi-squared
    %   distribution, respectively, for each bin i. 
    %
    % If the probability of observing the chisq statistic is less
    % than 0.025 (accounting for a two-tailed test), then the algorithm 
    % tries to split points in X via Spectral Clustering 
    % (more robust to odd shapes in the PCA space). The sub-clusters
    % are kept only if the sum of the mean squared error between the
    % points in each cluster to their own centers is less than the 
    % sum of the mean squared error of the points to the combined
    % cluster center. 
    % 
    % Inputs:
    %   X - an n x m matrix, with n = observations, m = dimensions
    %   
    %   nDim - the # of dimensions of the subspace to project onto
    %
    % Outputs:
    %   labels - a 1 x n vector of cluster labels
    %
    %   chisqStat - the chi squared statistic as described above
    %
    %   pval - the p-value of the chisqStat
    %
    % Written by Jordan Sorokin, 10/20/2017
    
    labels = ones( 1,size( X,1 ));
    
    % compute the chi-square statistic / pval that X is normally distributed
    [chisqStat,pval,V,Y] = chi2normal( X );

    % (d) if the p-value of the chisquared statistic is small, perform
    % a split of the subclusters using spectral clustering, which
    % tends to be robust to oddly-shaped clusters that may occur
    % from splits of the PCA projection space
    if pval <= 0.001 && V > 0.3                

        % perform sepctral clusetering on the projected data
        [~,A] = adjacency_graph( Y,'kNN',10 );
        subClusts = SpectralClustering( A,2,2 );

        % compare new average waveforms to old average by
        % computing the mean squared error betwween points and centers
        Z_sub = zeros(1,2);
        for j = 1:2
            Y_sub = Y(subClusts == j,:); 

            % compute the mean squared error of the
            % points-to-center of the new cluster
            Z_sub(1) = Z_sub(1) + mean( sum( (Y_sub - mean( Y_sub )).^2,2 ) );
            Z_sub(2) = Z_sub(2) + mean( sum( (Y_sub - mean( Y )).^2,2 ) ); 
        end
        
        % keep the clusters as separate clusters if the 
        % mean squared error of the sub clusters is smaller
        if Z_sub(1) < Z_sub(2)
            labels = subClusts;
        end        
    end
end