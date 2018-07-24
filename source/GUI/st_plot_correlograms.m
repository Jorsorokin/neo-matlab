function st_plot_correlograms( handles )
    % st_plot_correlograms( handles )
    %
    % plots the cross-correlograms between each i,j cluster, including the i-i auto-correlogram
    
    if isnan( handles.times )
        disp( 'No spike times available' );
        return
    end
    
    if ~any( handles.selectedClusts )
        return
    end
    
    uID = handles.clusterIDs( handles.selectedClusts );
    nID = nnz( handles.selectedClusts );
    nColors = size( handles.allPlotColor,1 ) - 1;
    bw = 0.0002; % .5 ms
    maxLag = 0.150; % 150 ms
    
    marg_h = [0.05 0.05]; marg_w = [0.04 0.04]; gap = [0.01 0.01];
    axh = (1-sum(marg_h)-(nID-1)*gap(1))/nID; % x positions
    axw = (1-sum(marg_w)-(nID-1)*gap(2))/nID; % y positions
    py = 1-marg_h(2)-axh; 
    
    idVec = 1:nID;
    plotCombo = true( nID,nID );
    axesColors = {'k','w'};
    grayScale = [0.65 0.65 0.65];
    
    for ih = idVec
        inds = find( plotCombo(ih,:) );
        for ix = inds
            plotCombo(ih,ix) = false;
            plotCombo(ix,ih) = false;
            px = marg_w(1) + axw*(ix-1) + gap(2)*(ix-1);
            sameClust = (ix==ih);
            xAxisColor = axesColors{sameClust+1}; % avoids an if-statement
            
            ax = axes(handles.xcorrplot,...
                'Units','normalized',...
                'Position',[px py axw axh],...
                'NextPlot','add',...
                'Visible','off' );
            
            set( ax,'xlim',[-maxLag maxLag]*1000,'ylim',[-0.01 1],'tickdir','out','yticklabel','',...
                'fontsize',8,'ycolor','k','xcolor',xAxisColor,'color','k' );
            
            % compute cross-correlogram between i,j cluster and plot
            ids = uID([ih,ix]);
            idx = handles.labels == ids(1);
            train1 = spikevec2mat( handles.times(idx),handles.trials(idx) );
            idx = handles.labels == ids(2);
            train2 = spikevec2mat( handles.times(idx),handles.trials(idx) );
            
            if size( train1,2 ) ~= size( train2,2 )
                continue
            end
            
            [xcg,lags] = correlogram( train1,train2,bw,maxLag );
            colorID = min( ids,mod( ids,nColors ) + uint8(ids ~= 0) ) + 1;
            color = handles.allPlotColor(colorID(1),:)*sameClust + grayScale*~sameClust;
            bar( ax,lags*1000,xcg./max(xcg),'FaceColor',color,'EdgeColor','none','BarWidth',1 );
        end
        py = py-axh-gap(1);
    end
    
    linkaxes( handles.xcorrplot.Children,'xy' ); 
    set( handles.xcorrplot.Children,'Visible', 'on' );
end
    