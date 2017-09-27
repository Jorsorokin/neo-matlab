function [ID,kept] = st_refine_cluster( ID,probabilities,minProb )
    % refine labels to discard points with less than a specific probability
    % of belonging to that cluster. This increases the specificity of
    % detected spikes at the cost of the total number of kept spikes...this
    % can be troublesome if the clusters are largely overlapping, as this
    % will throw away most of the points that exist within the overlap.

    % parse the probability matrix into each cluster
    kept = true( 1,numel( ID ) );
    nClust = numel( unique( ID ) );
    for i = 1:nClust
        if nClust == 1 % if only one group specified
            bad = medoutlier( probabilities,4 ); % remove outliers
            ID(bad) = 0;
            kept(bad) = 0;
        else
            index = ID==i;
            if sum( index ) <= 10 % less than 10 points in a cluster assumed to be noise
                ID(index) = 0; % make these equal to zero;
            else
                % if probability of being in cluster "i" is less than
                % minProb, drop the point
                bad = probabilities(index,i) <= minProb; 
                ID(index(bad)) = 0; % change ID to 0 for these bad spikes
                kept(index(bad)) = false; % change kept for these spikes to 0
                clear index bad 
            end
        end
    end
end 