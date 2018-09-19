function block = add_epochData_to_block( block,chans,eventTimes,preTime,postTime )
    % block = add_epochData_to_block( block,chans,eventTimes,preTime,postTime )
    %
    % Assumes the input folder "indir" contains files recorded from
    % OpenEphys, and that one has already created the Epoch objects
    % which contain start / ending times to use for data extraction
    
    % load in first channel
    [data,header] = OpenEphysLoader( block.filepath,'chans',chans(1) );
    
    % create epochs and remove boundary cases
    epochs = create_epochs( eventTimes,preTime,postTime,'event' );
    startPt = floor( [epochs.startTime] * header.fs );
    stopPt = floor( [epochs.stopTime] * header.fs );
    
    badEpochs = startPt < 1 | stopPt-1 > length( data );
    epochs(badEpochs) = [];
    startPt(badEpochs) = [];
    stopPt(badEpochs) = [];
    nEpoch = numel( epochs );
    
    if ~isrow( chans )
        chans = chans';
    end
    
    for ch = chans
        
        % load the data into the workspace
        [data,header] = OpenEphysLoader( block.filepath,'chans',ch );
        electrode = Electrode( ch );
        block.addChild( electrode );

        % extract the Epochs surrounding each event
        for ep = 1:nEpoch

            signal = Signal( data(startPt(ep):stopPt(ep)-1),header.fs ); % creates the Signal object

            % add the Signal to the corresponding parents
            epochs(ep).addChild( signal );
            electrode.addChild( signal );
        end
        clear data
    end % electrode loop
    
    % add to the block
    block.addChild( epochs );
end