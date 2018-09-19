function epochs = create_epochs( eventTimes,preTime,postTime,varargin )
    % epochs = create_epochs( eventTimes,preTime,postTime,(name) )
    
    nEpochs = numel( eventTimes );
    
    for ep = 1:nEpochs
        epochs(ep) = Epoch( eventTimes(ep)-preTime,eventTimes(ep)+postTime,ep );
        epochs(ep).eventTime = eventTimes(ep);
    end
    
    if nargin > 3 && ~isempty( varargin{1} ) 
        for ep = 1:nEpochs
            epochs(ep).name = sprintf( '%s%i',varargin{1},ep );
        end
    end
end