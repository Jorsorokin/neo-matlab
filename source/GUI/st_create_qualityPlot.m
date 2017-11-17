function fig = st_create_qualityPlot( handles )
    % ax = st_create_qualityPlot( handles )
    % creates an instance of the "qualityplot" used by the "sortTool" GUI
    %
    % handles is a structure with variables referenced by the main GUI

    fig = figure( 'Visible','off','color','k','position',[500,500,500,400] );
    set( gcf,'MenuBar','none','Name','Cluster quality' ); 
    
    if any( isnan( handles.R.clusterQuality ) )
        ylim = [-1 1];
        xlim = [0 1];
    else
        ylim = [min( handles.R.clusterQuality )-1,max( handles.R.clusterQuality )+1];
        xlim = [min( handles.labels ),max( handles.labels )];
    end
    
    ax = subplot( 2,1,1 ); hold on;
    set( ax,'tickdir','out','box','off','color','k',...
        'fontsize',8,'xcolor','w','ycolor','w',...
        'ylim',ylim,'xlim',xlim );
    xlabel( 'cluster' );
    ylabel( [handles.sortOptions.clusterMetric, ' score'] );

    ax = subplot( 2,1,2 ); hold on;
    set( ax,'tickdir','out','box','off','color','k',...
        'fontsize',8,'xcolor','w','ycolor','w',...
        'ylim',ylim,'xlim',[0 1] );
    xlabel( 'cluster scheme' );
    ylabel( [handles.sortOptions.clusterMetric, ' mean score'] );
end