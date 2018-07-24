function [chisqStat,pval] = chi2test( O,E,df )
    % [chisqStat,pval] = chi2test( O,E,df )
    %
    % compute the chi-square test statistic:
    % 
    %       sum( (O - E)^2 / E)
    %
    % which can be used for checking the goodness of fit of an observed
    % histogram count vector O with the expected count vector E of a
    % particular distribution, or if the data in O are dependent
    %
    % written by Jordan Sorokin, 7/18/18
    
    chisqStat = sum( (O - E).^2 ./ E );   
    
    if nargout > 1
        pvalVec = 0:.0001:1;
        chisqTable = chi2inv( pvalVec,df );
        pval = 1 - pvalVec(find( chisqStat < chisqTable,1 ));
    end
end