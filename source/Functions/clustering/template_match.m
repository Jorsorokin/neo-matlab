function [labels,P,outliers] = template_match( X,templates )
    % [labels,P,outliers] = template_match( X,templates )
    %
    % given a matrix of template vectors, matches each new point to each template
    % and outputs the label and probability of each point belonging to the
    % closest cluster. Also outputs outlier indices (those points with low
    % probability of belonging to a cluster)

    %X = X ./ std( X,[],2 );
    %templates = templates ./ std( templates,[],2 );

    scores = zeros( size( X,2 ),size( templates,2 ) ); 
    for j = 1:size( templates,2 )
        scores(:,j) = mean( bsxfun( @minus,X,templates(:,j) ).^2 )';
    end

    [bestScores,labels] = min( scores,[],2 );

    % get the probability of each point belonging to its cluster
    allLabels = unique( labels )';
    P = zeros( size( labels ) );
    for j = allLabels
        pts = (labels == j);
        minScore = min( bestScores(pts) );
        P(pts) = minScore ./ bestScores(pts);
    end

    outliers = find( P <= 0.05 ); 

    % scores = X' * templates;
    % [scores,labels] = max( scores,[],2 );
    % uID = unique( labels );
    % P = zeros( size( labels ) );
    % for j = uID'
    %     pts = (labels == j);
    %     maxScore = max( scores(pts) );
    %     P(pts) = scores(pts) ./ maxScore;
    % end
    % 
    % outliers = find( P < 0.05 );

end