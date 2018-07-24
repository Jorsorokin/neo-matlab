function channelindex = resort_neurons( channelindex,neuronIDs,templates )
    % channelindex = resort_neurons( channelindex,neuronIDs,templates )
    % 
    % resorts specific neurons in the channelindex object by comparing
    % spikes to the templates. This is a hard clustering,
    % with each spike assigned to one of the meanWaveforms.
    
    if channelindex.nUnits <= 1
        return
    end
    
    for id = neuronIDs
        neuron = findobj( channelindex.getChild( 'Neuron' ),'ID',id );
        [snips,~,~,mask] = neuron.getSpikes(); 
        snips = concatenateSpikes( maskchans( snips,mask ) );
        
        [~,labels] = max( templates.W(:,:,1)' * snips );
        channelindex.splitNeuron( id,templates.labelMap( labels ) );
    end
end
    
    
    