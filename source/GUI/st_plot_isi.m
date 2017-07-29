function st_plot_isi( handles )
    % function st_plot_isi( handles )
    %
    % plots the inter-spike interval (isi) of the spikes defined by 
    % handles.selectedPoints onto handles.isiplot
    
    % check if any spike times exist
    if isnan( handles.times )
        disp( 'No spike times available' );
        return;
    end
    
    % get the selected data points
    cla( handles.isiplot );
    handles.isiplot.Parent.Visible = 'on';
    pts = handles.selectedPoints; % only plots selected (highlighted) points
    if ~any( pts )
        return;
    end
    pltLims = handles.isiplot.XLim;
    edges = linspace( pltLims(1),pltLims(2),pltLims(2)*2000 );
    plotedges = edges(1:end-1);
    
    % loop over ids/trials
    selectedIDs = handles.labels(pts);
    selectedTrials = handles.trials(pts);
    uID = unique( selectedIDs );
    uTrials = unique( selectedTrials );
    for id = uID
        prob = 0;
        counter = 0;
        for trial = uTrials
            thesePts = (selectedIDs == id & selectedTrials == trial);
            prob = prob + histcounts( diff( handles.times(thesePts) ),edges,...
                'Normalization','pdf' );
            counter = counter + 1;
        end
        h = stairs( handles.isiplot,plotedges,prob/counter );
        hold on;
        set( h,'color',handles.allPlotColor(id+1,:) );
    end
end