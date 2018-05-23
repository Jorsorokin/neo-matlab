function color = st_getClusterColor( handles,clustID )
    % color = st_get_clustColor( handles,clustID )
    %
    % returns the color [r g b] for the cluster with specified ID
    
    nColors = size( handles.allPlotColor,1 );
    color = handles.allPlotColor( min( clustID,mod( clustID,nColors ) ),: );
end