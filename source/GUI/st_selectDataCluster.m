function [pts,label] = st_selectDataCluster( plotAxes,allLabels,varargin )
    % function pts = st_selectDataCluster( plotAxes,allLabels )
    %
    % allows the user to click on a given region within the axes of
    % "plotAxes", then returns all points that have the same label as the
    % point closest to the [x,y] position of the mouse click. 
    
    disp( 'Please select the cluster of interest' );
    axes( plotAxes );
    closestpt = selectdata( 'sel','closest','pointer','crosshair' );
    if iscell( closestpt )
        closestpt = closestpt{end};
    end
    label = allLabels( closestpt );
    pts = find( ismember( allLabels,label ) );
end
    