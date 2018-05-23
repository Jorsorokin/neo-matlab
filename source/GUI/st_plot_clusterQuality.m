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
    clusts = handles.clusterIDs;
    nClusts = numel( clusts );
    h = scatter( ax,clusts,handles.R.clusterQuality,...
                repmat( 100,1,numel( handles.R.clusterQuality ) ),'.' );
    color = zeros( nClusts,3 );
    for i = 1:nClusts
        color(i,:) = handles.plotcolor(find( handles.labels==clusts(i),1 ),:);
    end
    h.CData = color;
    ylim = [min( handles.R.clusterQuality )-1, max( handles.R.clusterQuality )+1];
    xlim = [min( handles.labels )-1, max( handles.labels )];
    set( ax,'ylim',ylim,'xlim',xlim );
    
    % highlight the clusters if any points in the projection plot belonging 
    % to those clusters that are highlighted
    ax.Children.SizeData(handles.labels(handles.selectedPoints)) = 500;

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