function [merge,dipStat,pval] = try_cluster_merge( X,labels )
    % [merge,dipStat,pval] = try_cluster_merge( X,labels )
    %
    % attempts to merge clusters in the n x m data matrix X by projecting each cluster pair 
    % onto the line separating their centers, and performing Hartigan's dip test to determine if the
    % projection is unimodal or not. If unimodal, then the i,j cluster pair are
    % considered a single cluster.
    %
    % Inputs:
    %   X - n x m data matrix, with n = points, m = dimensions
    %   labels - an n x 1 vector of cluster labels, for indexing
    %
    % Outputs:
    %   merge - a K x K matrix, with K = number of unique clusters, and each
    %           element a boolean indicating whether the i,j cluster pair
    %           should be merged (true) or not (false)
    %
    %   dipStat - a K x K matrix of hartigan's dip statistic
    %   pval - a K x K matrix of the p-value for hartigan's dip test
    %
    % By Jordan Sorokin, 5/20/18

    % get labels
    uID = unique( labels );

    % pre-allocate
    K = numel( uID );
    dipStat = zeros( K,K );
    pval = eye( K,K );
    tested = logical( eye( K ) );

    % loop over cluster pairs
    counter = 0;
    for i = 1:K
        inds_i = (labels==uID(i));
        cluster_i = X(inds_i,:);
        for j = 1:K
            if tested(i,j)
                continue
            end
            
            inds_j = (labels==uID(j));
            cluster_j = X(inds_j,:);

            % project onto the line separating the means
            L = mean( cluster_i ) - mean( cluster_j );
            proj = [cluster_i;cluster_j] * L';
            %P = histcounts( proj,'Normalization','PDF' );
            %P = ksdensity( proj );
            
            % perform hartigan's dip test on the PDF
            [dipStat(i,j),pval(i,j)] = HartigansDipSignifTest( proj,500 );
            dipStat(j,i) = dipStat(i,j);
            pval(j,i) = pval(i,j);
            
            % update our counting variables
            tested(i,j) = true;
            tested(j,i) = true;
            counter = counter + 1;
        end
    end

    % determine which clusters should be merged
    merge = pval > (0.05 / counter); % correct for # of tests we performed
end