function st_update_scatter( handles )
    % updates the scatter plot of the projections to reflect new 
    % view angles of the data
    
    scatterView = handles.projection * handles.plotDims;
    h = handles.mapplot.Children(end);
    h.XData = scatterView(:,1);
    h.YData = scatterView(:,2);
end