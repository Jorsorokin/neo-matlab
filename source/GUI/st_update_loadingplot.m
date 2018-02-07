function st_update_loadingplot( handles,axis )
    % updates the loading plot with each new projection viewpoint
    
    % check if loading plot is open
    if isempty( handles.loadingplot )
        return
    end
    
    % get the current loadings and update the Y-data for each subplot\
    if axis == 1
        child = 2;
    else
        child = 1;
    end
    nDims = size( handles.plotDims,1 );
    ax = handles.loadingplot.Children(child);
    h = ax.Children;
    h.YData = handles.plotDims(:,axis);
    h.XData = 1:numel( h.YData );
    set( ax,'xlim',[0,nDims+1],'xtick',1:nDims );
end