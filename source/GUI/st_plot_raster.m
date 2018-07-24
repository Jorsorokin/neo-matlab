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
    
    if ~any( handles.selectedClusts )
        return
    end
    
    % get the selected data points
    pts = handles.selectedPoints; % only plots selected (highlighted) points
    IDs = handles.clusterIDs(handles.selectedClusts);
    nID = numel( IDs );
    nColors = size( handles.allPlotColor,1 ) - 1;
    
    % create psth matrix
    bw = 0.02; % 40 ms
    gaussKernel = gausswin( 11,2 ); % 200ms smoothing around peak
    gaussKernel = gaussKernel / sum( gaussKernel );
    xlims = handles.rasterplot.Children(1).XLim;
    nBins = ceil( xlims(2) / bw );
    bins = linspace( xlims(1),xlims(2),nBins+1 );
    psth = zeros( nBins,nID,'single' );
    
    % loop over trials, plot a scatter of points along rows for each trial,
    % then color code the scatter points by the labels
    uTrials = unique( handles.trials(pts) );
    for trial = uTrials
        thesePts = (pts & handles.trials == trial);
        h = scatter( handles.rasterplot.Children(2),...
            handles.times(thesePts),... 
            ones( 1,sum( thesePts ),class( trial ) )*trial,... %+ randn( 1,sum( thesePts ),'single' )*0.01,... 
            50,'.' );
        h.CData = handles.plotcolor(thesePts,:);
        
        for id = 1:nID
            psth(:,id) = psth(:,id) + histcounts( handles.times(thesePts & handles.labels == IDs(id)),bins )';
        end
    end
    
    psth = psth * (1/numel( uTrials )) * (1/bw);
    for id = 1:nID
        thisID = IDs(id);
        thisColor = min( thisID,mod( thisID,nColors ) + uint8(thisID ~= 0) ) + 1;
        plot( handles.rasterplot.Children(1),bins(2:end),conv( psth(:,id),gaussKernel,'same' ),'color',handles.allPlotColor(thisColor,:) );
    end
end