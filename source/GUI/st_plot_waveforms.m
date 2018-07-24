function st_plot_waveforms( handles ) 
    % function st_plot_waveforms( handles )
    %
    % plots the raw waveforms contained in "handles.data" onto the 
    % "handles.waveformplot" axes. This function is part of the suite
    % belonging to the sortTool GUI
    
    % check if any data exists
    if ~handles.availableData
        fprintf( 'No data available\n' );
        return
    end

    pts = handles.selectedPoints; % only plots selected (highlighted) points
    if ~any( pts )
        return
    end
    
    uID = handles.clusterIDs(handles.selectedClusts);
    [nPt,~,nChan] = size( handles.data );
    nAxes = numel( handles.waveformplot.Children );
    axDiff = nAxes - nChan;
    
    % loop over channels, plot waveforms for that channel
    if nnz( pts ) >= 100
        time = 1:nPt;
        for i = uID
            idx = (handles.labels==i) & pts;
            color = handles.plotcolor( find( idx,1 ),: );
            mean_waveform = squeeze( mean( handles.data(:,idx,:),2 ) );
            sd_waveform = squeeze( std( handles.data(:,idx,:),[],2 ) );
            for c = nChan:-1:1
                line( handles.waveformplot.Children(c+axDiff),time,mean_waveform(:,c),...
                    'color',color,'linewidth',2 );
                line( handles.waveformplot.Children(c+axDiff),time,mean_waveform(:,c) + sd_waveform(:,c),...
                    'linestyle','--','color',color,'linewidth',0.5 );
                line( handles.waveformplot.Children(c+axDiff),time,mean_waveform(:,c) - sd_waveform(:,c),...
                    'linestyle','--','color',color,'linewidth',0.5 );
            end
        end
    else
        timevec = [1:nPt,nan];
        for i = uID
            idx = (handles.labels == i) & pts; 
            color = handles.plotcolor( find( idx,1 ),: );
            nSp = nnz( idx );
            totalPts = nPt*nSp + nSp; % extra points for nans
            time = repmat( timevec,1,nSp );
            for c = nChan:-1:1
                bigline = reshape( [handles.data(:,idx,c);nan( 1,nSp )],totalPts,1 ); 
                line( handles.waveformplot.Children(c+axDiff),time,bigline,...
                    'color',color,'linewidth',0.5 );
            end
        end
    end
end