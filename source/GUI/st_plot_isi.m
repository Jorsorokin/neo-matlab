function st_plot_isi( handles )
    % function st_plot_isi( handles )
    %
    % plots the inter-spike interval (isi) of the spikes defined by 
    % handles.selectedPoints onto handles.isiplot
    
    % check if any spike times exist
    if isnan( handles.times )
        disp( 'No spike times available' );
        return
    end

    pts = handles.selectedPoints; % only plots selected (highlighted) points
    if ~any( pts )
        return
    end
    pltLims = handles.isiplot.Children.XLim;
    edges = logspace( log10( pltLims(1) ),log10( pltLims(2) ),round( pltLims(2)/10 ) );
    
    % loop over ids/trials
    uID = handles.clusterIDs(handles.selectedClusts);
    uTrials = unique( handles.trials(pts) );
    nColors = size( handles.allPlotColor,1 ) - 1;
    
    for id = uID
        isi = [];
        for trial = uTrials
            thesePts = (handles.labels == id & handles.trials == trial & pts);
            isi = [isi,diff( handles.times(thesePts)*1000 )];
        end
        
        prob = histcounts( isi(isi<=pltLims(2)),edges,'Normalization','probability' );
        thisColor = min( id,mod( id,nColors ) + uint8(id ~= 0) ) + 1;
        stairs( handles.isiplot.Children,edges(2:end),prob,'color',handles.allPlotColor(thisColor,:) );
    end
end