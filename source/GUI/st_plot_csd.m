function st_plot_csd( handles )
    % plots the current source density for selected points
    
    if ~handles.availableData
        fprintf( 'No data available\n' );
        return
    end
    
    if ~any( handles.selectedPoints )
        return
    end
    
    nColors = size( handles.allPlotColor,1 ) - 1;
    dt = max( var( handles.data(:,handles.selectedPoints,:) ) );
    for clust = handles.clusterIDs(handles.selectedClusts)
        idx = handles.selectedPoints & (handles.labels == clust);
        W = diff( diff( squeeze( mean( handles.data(:,idx,:),2 ) ),[],2 ),[],2 );
        thisColor = min( clust,mod( clust,nColors ) + uint8(clust ~= 0) ) + 1;
        multisignalplot( W,[],handles.allPlotColor(thisColor,:),dt );
    end
end