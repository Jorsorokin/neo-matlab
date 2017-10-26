function pts = st_selectDataLasso( plotAxes )
    % function pts = st_selectDataLasso( plotAxes )
    %
    % allows the user to draw a lasso in the axes defined by "plotAxes" and
    % returns the indices of points within that lasso
    
    % check if we have something pressed in the toolbar
    if ~isempty( plotAxes.Parent.WindowKeyPressFcn )
        disp( 'Please un-select the toolbar button that is highlighted' );
        pts = [];
        return
    end
    
    disp( 'Select some data points' );
    axes( plotAxes );
    pts = selectdata( 'sel','lasso','Fill','off' );
end