function st_plot_dvdt( handles )
    % plots the phase-portrait (V vs dVdT) of the mean spike waveforms for
    % each cluster
    
    if any( isnan( handles.data ) )
        fprintf( 'No data available\n' );
        return
    end
    
    if ~any( handles.selectedPoints )
        return
    end
    
    nColors = size( handles.allPlotColor,1 ) - 1;

    % loop over selected clusters       
    if ~any( isnan( handles.mask ) )
        for clust = handles.clusterIDs(handles.selectedClusts)
            idx = handles.selectedPoints & (handles.labels == clust);
            W = squeeze( mean( handles.data(:,idx,:),2 ) );
            W = W .* mean( handles.mask(:,idx),2 )';
            W = W / abs( min( W(:) ) ); % normalizes to maximum amp on best chan
            
            thisColor = min( clust,mod( clust,nColors ) + uint8(clust ~= 0) ) + 1;
            line( handles.dvdtplot.Children,W(2:end,:),diff( W ),...
                'color',handles.allPlotColor(thisColor,:) );
        end
    else
        for clust = handles.clusterIDs(handles.selectedClusts)
            W = squeeze( mean( handles.data(:,handles.selectedPoints & (handles.labels == clust),:),2 ) );
            W = W / abs( min( W(:) ) ); % normalizes to maximum amp on best chan
            
            thisColor = min( clust,mod( clust,nColors ) + uint8(clust ~= 0) ) + 1;
            line( handles.dvdtplot.Children,W(2:end,:),diff( W ),...
                'color',handles.allPlotColor(thisColor,:) );
        end
    end
    axis tight
end
        