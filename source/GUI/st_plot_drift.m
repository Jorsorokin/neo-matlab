function st_plot_drift( handles )
    % function st_plot_drift( handles )
    %
    % plots the amplitude / instantaneous ISI of selected points over time
    
    % check if any data exists
    if ~handles.availableData || all( isnan( handles.times ) )
        fprintf( 'No data or spike times available\n' );
        return
    end

    if ~any( handles.selectedClusts )
        return
    end
        
    % get the selected data points
    pts = handles.selectedPoints; % only plots selected (highlighted) points
    IDs = handles.clusterIDs(handles.selectedClusts);
    maxTime = max( handles.times );
    
    for i = IDs
        thesePts = pts & (handles.labels == i); 
        time = handles.times(thesePts);
        trials = (handles.trials(thesePts) - 1) * maxTime;
        if any( isnan( handles.mask ) )
            amp = squeeze( min( handles.data(:,thesePts,:) ) );
            meanAmp = mean( amp );
            noise = median( meanAmp );
            [~,bestchan] = min( meanAmp );
            nChan = sum( amp < (noise - 4*std( meanAmp )),2 );
            amp = amp(:,bestchan);
        else
            m = handles.mask(:,thesePts);
            [~,bestchan] = max( sum( m,2 ) );
            amp = min( handles.data(:,thesePts,bestchan) );
            nChan = sum( m >= 0.5 );
        end
            
        h = scatter( handles.driftplot.Children(2),...
                 time + trials,amp,50,'.' );
        h.CData = handles.plotcolor(thesePts,:);
        
        h = scatter( handles.driftplot.Children(1),...
                    time + trials,nChan,50,'.' );
        h.CData = handles.plotcolor(thesePts,:);
    end
    
    set( handles.driftplot.Children(1),'YLimSpec','tight' );
    set( handles.driftplot.Children(2),'YLimSpec','tight' );
end