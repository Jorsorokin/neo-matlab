function [chisqStat,pval] = check_chisquare_distribution( X,df )
    % [chisqStat,pval] = check_chisquare_distribution( X,df )
    %
    % computes the chi squared statistic that the values in X
    % are chi-squared distributed wtih "df" degrees of freedom
    % Under the null hypothesis, X is a chi-squared distributed
    % variable with "df" degrees of freedom. 
    %
    % The chi-squared statistic is calculated as:
    %
    %   chisqStat = SUM_i{ (O(i) - E(i))^2 / E(i) }
    %
    %   where O(i) and E(i) are the observed and expected values
    %   of the empirical and theoretical CDFs of X and a df-chi
    %   squared distribution, respetively, at each bin i
    %
    % Inputs:
    %   X - vector or matrix
    %
    %   df - scalar specifying the degrees of freedom of X
    %
    % Outputs:
    %   chisqStat - the test statistic as described above
    %
    %   pval - P( x <= chisqStat ) for a df-chi squared distributed 
    %          random variable x, rounded to the nearest 1e-4 probability
    %
    % Written by Jordan Sorokin, 10/20/2017
    
    % get the optimal histogram bins
    nPts = length( X );
    [O,bins] = histcounts( X,'BinMethod','fd' );    % for heavy-tailed dists
    
    % compute the theoritical cdf of a chi squared variable
    % over the supplied bins
    P = gammainc( bins(2:end)/2,df/2 );             % probability of data in [ 0,bin(i) ) for chisq with DOF = df
    E = P * nPts;                                   % total # of counts expected, given nPts
   
    % check if the data is chi-square distributed, using 
    % the empirical CDF vs PDF
    [chisqStat,pval] = chi2test( cumsum( O ),E,df ); 
end