function fig = st_create_loadingplot( handles )
    % creates the loading plot for viewing the loadings of the different
    % dimensions onto the 2D scatter projection
    fig = figure( 'Position',[500,600,500,350],'Color','k','menu','none','Visible','off' );
    
    totalDims = handles.nDim + handles.nLocationDim;
    ynames = {'x-axis loading','y-axis loading'};
    for j = 1:2
        ax = subplot( 2,1,j );
        stem( 'parent',ax,handles.plotDims(:,j) );
        set( ax,'tickdir','out','box','off',...
            'xlim',[0,totalDims+1],'ylim',[-1 1],...
            'color','k','xcolor','w','ycolor','w',...
            'xtick',1:totalDims );
        ax.XLabel.String = 'projection dim';
        ax.YLabel.String = ynames{j};
    end 
    
    fig.Visible = 'on';
end