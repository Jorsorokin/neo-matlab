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
    
    % get the selected data points
    cla( handles.isiplot );
    handles.isiplot.Parent.Visible = 'on';
    pts = handles.selectedPoints; % only plots selected (highlighted) points
    if ~any( pts )
        return
    end
    pltLims = handles.isiplot.XLim;
    edges = linspace( pltLims(1),pltLims(2),pltLims(2)*2000 );
    plotedges = edges(1:end-1);
    
    % loop over ids/trials
    selectedIDs = handles.labels(pts);
    selectedTrials = handles.trials(pts);
    uID = handles.clusterIDs(handles.selectedClusts);
    uTrials = unique( selectedTrials );
    nColors = size( handles.allPlotColors,1 );
    
    for id = uID
        isi = [];
        for trial = uTrials
            thesePts = (selectedIDs == id & selectedTrials == trial);
            isi = [isi,diff( handles.times(thesePts) )];
        end
        prob = histcounts( isi,edges,'Normalization','pdf' );
        h = stairs( handles.isiplot,plotedges,prob );
        hold on;
        set( h,'color',handles.allPlotColor(min( id,mod( id,nColors ) )) );
    end
end