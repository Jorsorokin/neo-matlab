function st_plot_clusterQuality( handles )
    % function st_plot_clusterQuality( handles)
    %
    % plots cluster quality, color-coded by the cluster
    
    % turn on the figure
    handles.qualityplot.Visible = 'on';
    
    % end if no cluster quality available
    if any( isnan( handles.R.clusterQuality ) )
        return
    end
    
    % get the first subplot, plot cluster scores
    ax = handles.qualityplot.Children(2);
    ax.YLabel.String = [handles.sortOptions.clusterMetric,' score'];
    
    cla( ax );
    clusts = unique( handles.labels );
    h = scatter( ax,clusts,handles.R.clusterQuality,...
                repmat( 100,1,numel( handles.R.clusterQuality ) ),'.' );
    h.CData = handles.allPlotColor(clusts+1,:);
    ylim = [min( handles.R.clusterQuality )-1, max( handles.R.clusterQuality )+1];
    xlim = [min( handles.labels )-1, max( handles.labels )];
    set( ax,'ylim',ylim,'xlim',xlim );
    
    % highlight the clusters if any points in the projection plot belonging 
    % to those clusters that are highlighted
    pts = handles.selectedPoints;
    if any( pts )
        clusts = unique( handles.labels(pts) );
        ax.Children.SizeData(pts) = 300;
    end 

    % now plot running cluster averages
    ax = handles.qualityplot.Children(1);
    ax.YLabel.String = [handles.sortOptions.clusterMetric, ' mean score'];
    
    if min( handles.labels ) == 0
        meanQuality = mean( handles.R.clusterQuality(2:end) );
    else
        meanQuality = mean( handles.R.clusterQuality );
    end

    if isempty( ax.Children )
        data = [nan,meanQuality];
        xlim = [0 1];
    else
        data = [ax.Children.YData,meanQuality];
        xlim = [0, numel( data )-1];
    end

    cla( ax )
    plot( ax,xlim(1):xlim(end),data,'wo--','MarkerFaceColor','w' );
    set( ax,'xlim',xlim,'ylim',[min( data )-1, max( data )+1] );
end