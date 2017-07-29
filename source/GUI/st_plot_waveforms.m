function st_plot_waveforms( handles ) 
    % function st_plot_waveforms( handles )
    %
    % plots the raw waveforms contained in "handles.data" onto the 
    % "handles.waveformplot" axes. This function is part of the suite
    % belonging to the sortTool GUI
    
    % check if any data exists
    if isnan( handles.data )
        disp( 'No data available' );
        return;
    end
    
    % get IDs of the labels, then loop and plot color coded waveforms
    cla( handles.waveformplot );
    handles.waveformplot.Parent.Visible = 'on';
    pts = handles.selectedPoints; % only plots selected (highlighted) points
    numNans = 1;
    if ~any( pts )
        return;
    end
    uID = unique( handles.labels(pts) );
    for i = uID
        idx = handles.labels == i & pts;
        nSp = sum( idx );
        nPt = size( handles.data,1 );
        totalPts = nPt*nSp + nSp*numNans; % extra points for nans

        % add nans between each spike waveform (make one big line)
        % this helps reduce overhead of having MANY "line" objects
        time = repmat( [1:nPt,nan( 1,numNans )],1,nSp );
        bigline = reshape( [handles.data(:,idx);nan( numNans,nSp )],...   
            totalPts,1 ); 
        plot( handles.waveformplot,time,bigline,...
            'color',handles.allPlotColor(i+1,:),...
            'linewidth',0.5 );
    end
end