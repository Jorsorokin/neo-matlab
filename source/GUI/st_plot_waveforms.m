function st_plot_waveforms( handles ) 
    % function st_plot_waveforms( handles )
    %
    % plots the raw waveforms contained in "handles.data" onto the 
    % "handles.waveformplot" axes. This function is part of the suite
    % belonging to the sortTool GUI
    
    % check if any data exists
    if ~handles.availableData
        disp( 'No data available' );
        return
    end
    
    %handles.waveformplot(1).Parent.Visible = 'on';
    pts = handles.selectedPoints; % only plots selected (highlighted) points
    numNans = 1;
    if ~any( pts )
        return;
    end
    
    % clear previous plot
    for c = 1:numel( handles.waveformplot )
        cla( handles.waveformplot(c) );
    end
    
    uID = unique( handles.labels(pts) );
    [nPt,~,nChan] = size( handles.data );
    
    for i = uID
        idx = (handles.labels == i) & pts; 
        nSp = sum( idx );
        totalPts = nPt*nSp + nSp*numNans; % extra points for nans

        % add nans between each spike waveform (make one big line)
        % this helps reduce overhead of having MANY "line" objects
        time = repmat( [1:nPt,nan( 1,numNans )],1,nSp );

        % loop over channels
        for c = 1:nChan
            bigline = reshape( [handles.data(:,idx,c);nan( numNans,nSp )],totalPts,1 ); 
            line( handles.waveformplot(c),time,bigline,...
                'color',handles.allPlotColor(i+1,:),...
                'linewidth',0.5 );
        end
    end
end