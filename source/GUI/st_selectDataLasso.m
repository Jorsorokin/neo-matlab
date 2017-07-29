function pts = st_selectDataLasso( plotAxes )
    % function pts = st_selectDataLasso( plotAxes )
    %
    % allows the user to draw a lasso in the axes defined by "plotAxes" and
    % returns the indices of points within that lasso
    disp( 'Select some data points' );
    axes( plotAxes );
    pts = selectdata( 'sel','lasso','Fill','off' );
end