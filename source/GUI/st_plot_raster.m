function st_plot_raster( handles )
    % function st_plot_raster( handles )
    %
    % plots a raster of the selected spikes, color coded by their labels
    % and separated along the y-axis by trials 
    
    % check if any spike times exist
    if isnan( handles.times )
        disp( 'No spike times available' );
        return
    end
    
    % get the selected data points
    cla( handles.rasterplot );
    handles.rasterplot.Parent.Visible = 'on';
    pts = handles.selectedPoints; % only plots selected (highlighted) points
    if ~any( pts )
        return
    end
    
    % loop over trials, plot a scatter of points along rows for each trial,
    % then color code the scatter points by the labels
    selectedTrials = handles.trials(pts);
    uTrials = unique( selectedTrials );
    for trial = uTrials
        thesePts = (pts & handles.trials == trial);
        h = scatter( handles.rasterplot,...
            handles.times(thesePts),...
            ones( 1,sum( thesePts ),class( trial ) )*trial,...
            30,'.' );
        h.CData = handles.plotcolor(thesePts,:);
    end
end