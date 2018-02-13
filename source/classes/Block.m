classdef Block < Container
    % self = Block( filename,date,condition,filepath )
    %
    % Initiate an instance of the Block class. A Block is 
    % the top-level container for a given recording session,
    % and contains ChannelIndex and Epoch classes as children.
    % 
    % The main purpose of the Block container is to tie the 
    % data (channels, epochs, signals, neurons, & spikes)
    % together into a unified recording session for storage
    % and organization.
    %
    % Children:
    %   Epoch
    %   ChannelIndex
    %   Electrode
    %
    % Parents:
    %   none
    %
    % Properties:
    %   filename - name of the recorded data file
    %   date - date recorded
    %   condition - any condition specific to this recording (i.e. "nose poke", etc.)
    %   filepath - path of the recorded data file
    %   nEpochs - number of Epochs extracted
    %   nChanInds - number of ChannelIndex objects extracted
    %   nUnits - number of Neurons identified
    %   nSignals - number of raw voltage signals (total across channels/epochs)   
    %
    % Methods:
    %   print
    %   update
    %   write
    %   getNeurons
    %   undoSorting
    %   mergeSimilarNeurons
    %   saveSpikeWaveforms
    %
    %       * see also methods in the Container class
    
    properties
        filename
        date
        condition
        filepath
        nElectrodes = 0;
        nEpochs = 0;
        nChanInds = 0;
        nUnits = 0;
        nSignals = 0;
    end
    
    methods
        
        function self = Block( filename,date,condition,filepath )
            % self = Block( filename,date,condition,filepath )
            %
            % initiates the Block object

            self.filename = filename;
            self.date = date;
            self.condition = condition;
            self.filepath = filepath;
        end
        

        function addChild(self,child)
            switch class( child )
                case 'Epoch'
                    addChild@Container( self,child );
                    self.nEpochs = self.nEpochs + numel( child );
                case 'ChannelIndex'
                    addChild@Container( self,child );
                    self.nChanInds = self.nChanInds + numel( child );
                case 'Electrode'
                    addChild@Container( self,child );
                    self.nElectrodes = self.nElectrodes + numel( child );
                otherwise 
                    error( 'Only ChannelIndex, Electrode, & Epoch objects are valid children' );
            end
        end
        

        function addParent(~,~)
            disp( 'Block objects have no parents' );
        end   


        function print( self )
            % print( self )
            % 
            % displays the file metadata and children
            self.update();
            fprintf('----------------------------------------------\n');
            fprintf( '%s (%s)\nRecorded on %s\n\n',...
                self.filename,self.condition,self.date );
            fprintf( '%i electrodes\n%i signals\n%i channel groups\n%i epochs\n%i units\n',...
                self.nElectrodes, self.nSignals, self.nChanInds, self.nEpochs,self.nUnits );
            fprintf('----------------------------------------------\n');
        end
        
        
        function update( self )
            % update( self )
            %
            % updates the properties contained within this block according
            % to its children. It searches through the Epoch and
            % ChannelIndex children (if any) and updates the number of
            % epochs, channels, neurons, electrodes, and signals found.
            epoch = self.getChild( 'Epoch' );
            chanind = self.getChild( 'ChannelIndex' );
            electrode = self.getChild( 'Electrode' );
            neurons = self.getNeurons();
            %signals = self.getSignals();

            % check for validity & update: ChannelIndex
            % =========================================
            if ~isempty( chanind )
                chanind(~isvalid( chanind )) = []; % eliminates bad channel indices
                if ~isempty( chanind )     
                    self.nElectrodes = 0;
                    chanind(~isvalid( chanind )) = [];
                    hasElectrodes = find( [chanind.nElectrodes] );
                    for j = hasElectrodes
                        chanElectrodes = chanind(j).getChild( 'Electrode' );
                        electrodeID = [chanElectrodes.electrodeNum];
                        chanind(j).nElectrodes = numel( electrodeID );
                        chanind(j).chanIDs = electrodeID; % the actual electrode IDs 
                        chanind(j).nUnits = numel( chanind(j).getChild( 'Neuron' ) );
                    end
                end
            end
            self.nChanInds = numel( chanind );
            % =========================================
            
            % check validity & update: Electrode 
            % ==================================
            if isempty( electrode )
                self.nElectrodes = 0;
            else
                electrode(~isvalid( electrode )) = []; % eliminates bad electrodes
                for j = 1:numel( electrode )
                    chanind = electrode(j).getParent( 'ChannelIndex' );
                    if ~isempty( chanind )
                        electrode(j).chanInd = [chanind.chanIndNum]; % all parent ChannelIndex 
                        electrode(j).nChanInd = numel( chanind );
                    else
                        electrode(j).chanInd = [];
                        electrode(j).nChanInd = 0;
                    end                    

                    % update the signals in this electrode
                    for sig = 1:electrode(j).nSignals
                        signal = electrode(j).getChild( 'Signal',sig );

                        % check validity of Signal object
                        if ~isempty( signal )
                            if ~isvalid( signal )
                                electrode(j).removeChild( 'Signal',sig );
                                clear signal; % removes this invalid signal object
                                continue
                            end
                            signal.epoch = signal.getParent( 'Epoch' ).epochNum;
                            signal.chanInd = electrode(j).chanInd; % the parent ChannelIndex
                            signal.electrode = electrode(j).electrodeNum; % this electrode ID
                        end
                    end
                end
                self.nElectrodes = numel( electrode );
            end
            % =================================

            % check validity & update: Epoch
            % ==============================
            nSig = 0;
            nEp = 0;
            if ~isempty( epoch )
                epoch(~isvalid( epoch )) = [];
                nEp = numel( epoch );
                
                % loop over spike/signal objects, remove those that are not valid
                if nEp > 0
                    for ep = 1:nEp
                        spikes = epoch(ep).getChild( 'Spikes' );
                        signals = epoch(ep).getChild( 'Signal' );
                        if ~isempty( signals )
                            signals(~isvalid(signals)) = [];
                        end
                        nSig = nSig + numel( signals );
                        epoch(ep).nSignals = numel( signals );

                        if ~isempty( spikes )
                            spikes(~isvalid( spikes )) = [];
                            for sp = 1:numel( spikes )
                                if isempty( spikes(sp).getParent( 'Neuron' ) )
                                    spikes(sp).deleteSelf(); 
                                else
                                    spikes(sp).nSpikes = numel( spikes(sp).times );
                                end
                            end
                            spikes(~isvalid( spikes )) = [];
                        end
                    end
                end
            end
            self.nSignals = nSig;
            self.nEpochs = nEp;
            % =========================================

            % check validity & update: Neuron
            % ==============================
            self.nUnits = 0;
            if ~isempty( neurons )
                neurons(~isvalid( neurons )) = [];
                for j = 1:numel( neurons )
                    if isempty( neurons(j).getParent( 'ChannelIndex' ) )
                        neurons(j).deleteSelf();
                    end
                end
                neurons(~isvalid( neurons )) = [];
                self.nUnits = numel( neurons );
                if self.nUnits > 0
                    for n = 1:self.nUnits
                        neurons(n).nSpikes = sum( [neurons(n).getChild( 'Spikes' ).nSpikes] );
                    end
                end  
            end
            % =========================================
        end
        
        
        function write( self,outpath )
            % write( self,outpath )
            %
            % saves the Block and its children into the folder specified by
            % "outpath". The file will be saved as "(self.filename)_extractedData.mat".
            % The actual matlab variable itself will be saved as "block"
            if ~isdir( outpath )
                mkdir( outpath );
                addpath( outpath );
            end
            
            disp( 'Writing the Block to disk...' );
            
            self.update(); % update the num channels/epochs/units
            eval( ['block' '=self'] ); % makes variable "block" reference "self"
            save( [outpath,filesep,self.filename,'_extractedData.mat'],'block','-v7.3' );
        end
        
        
        function neurons = getNeurons( self )
            % neurons = getNeurons( self )
            %
            % get all Neuron objects contained within this block.
            chanind = self.getChild( 'ChannelIndex' );
            neurons = [];
            if ~isempty( chanind )
                for j = numel( chanind ):-1:1
                    neurons = [chanind(j).getChild( 'Neuron' ),neurons]; % negative indexing to pre-allocate max array first
                end
            end
        end
        
        
        function undoSorting( self )
            % undoSorting( self )
            %
            % re-merge all spike objects and Neurons,
            % undoing results of a previous sorting
                       
            epoch = self.getChild( 'Epoch' );
            
            for ep = 1:self.nEpochs
                spikes = epoch(ep).getChild('Spikes');
                epoch(ep).removeChild( 'Spikes' );
                epoch(ep).nSpikes = 0;
                
                chanIndNum = zeros( 1,numel( spikes ) );
                for sp = 1:numel( spikes )
                    chanIndNum(sp) = spikes(sp).getParent( 'Neuron' ).chanInd;
                end
                
                % get the voltage waveforms / times / mask
                for ch = 1:numel( unique( chanIndNum ) )
                    thisChanInd = chanIndNum == ch;
                    oldSpikes = spikes(thisChanInd);
                    newSpikes = Spikes( [oldSpikes.times],[oldSpikes.voltage],oldSpikes(1).fs );
                    [newSpikes.times,idx] = sort( newSpikes.times );
                    if ~isempty( newSpikes.voltage )
                        newSpikes.voltage = newSpikes.voltage(:,idx,:);
                    end
                    newSpikes.mask = [oldSpikes.mask];
                    if ~isempty( newSpikes.mask )
                        newSpikes.mask = newSpikes.mask(:,idx);
                    end
                    newSpikes.chanInd = ch;
                    
                    epoch(ep).addChild( newSpikes );
                end
                clear spikes newSpikes oldSpikes
            end
            
            for ch = 1:self.nChanInds
                chan = self.getChild( 'ChannelIndex',ch );
                chan.removeChild( 'Neuron' );
                neuron = Neuron( 0 );
                
                for ep = 1:numel( epoch )
                    allSpikes = epoch(ep).getChild( 'Spikes' );
                    neuron.addChild( allSpikes.findobj( 'chanInd',ch ) );
                end
                
                chan.addChild( neuron );
                clear neuron
            end
    
            self.update;
        end
        
        function mergeSimilarNeurons( self,corrThresh )
            % mergeSimilarNeurons( self,corrThresh )
            %
            % Merge neurons for each ChannelIndex whos mean waveforms have
            % a correlation coefficient >= "corrThresh"

            % check for any neurons
            if self.nUnits == 0
                disp( 'No neurons in this block' );
                return
            end

            % loop over channel indices
            for ch = 1:self.nChanInds
                fprintf( 'Merging neurons for ChannelIndex %i\n',ch );
                chanind = self.getChild( 'ChannelIndex',ch );
                neurons = chanind.getChild( 'Neuron' );
                nunits = numel( neurons );
                if nunits <= 1
                    continue
                end

                % preallocate vars
                nPts = numel( neurons(1).meanWaveform );
                meanWaveform = zeros( nPts,nunits );
                IDs = [neurons.ID];

                % get the mean waveform
                for j = 1:nunits
                    if isempty( neurons(j).meanWaveform )
                        spikes = neurons(j).getSpikes;
                        neurons(j).meanWaveform = squeeze( mean( spikes,2 ) );
                    end
                    meanWaveform(:,j) = reshape( neurons(j).meanWaveform,nPts,1 );
                end
                clear neurons

                % find similar waveforms and merge those neurons
                corrMat = corr( meanWaveform ) - eye( nunits );
                [i,j] = find( corrMat >= corrThresh );
                unique_i = unique( i' );
                nFinalUnits = size( meanWaveform,2 );
                if ~isempty( i )
                    for k = unique_i
                        inds = (i==k);
                        chanind.mergeNeurons( [IDs(k),IDs(j(inds))] );
                    end

                    % re-order the neuron children
                    [~,idx] = sort( [chanind.getChild('Neuron').ID] );
                    neuronChild = find( contains( cellfun( @class,chanind.children,'un',0 ),'Neuron' ) );
                    chanind.children{neuronChild} = chanind.children{neuronChild}(idx);
                    nFinalUnits = numel( idx );
                end
                
                % add channel locations
                neurons = chanind.getChild( 'Neuron' );
                distances = chanind.chanDistances;
                for k = 1:numel( neurons )
                    [~,neurons(k).bestElectrode] = min( min( neurons(k).meanWaveform ) );
                    neurons(k).location = distances(neurons(k).bestElectrode,:);
                end

                fprintf( 'Merged %i neurons\n',nunits - nFinalUnits );
            end

            self.update();  
        end
        
        function saveSpikeWaveforms( self,spikesDir )
            % saveSpikeWaveforms( self,spikesDir )
            % 
            % save the voltage waveforms of all the spikes contained in
            % this block into the "spikesDir" folder. The saved file is a
            % .mat file, with the spike waveforms in the variable "spikes", the
            % masking matrices in "mask", the trials of the spikes in the
            % variable "epochs", the shank from which the spikes came from in the 
            % variable "unitCH', and the neuron ID in the variable
            % "unitID".
            %
            % Note that Each Spike object will have its voltage
            % waveforms and masking matrix removed, so that the parent
            % Block will use less disk space!!
            
            % pre-allocate
            neurons = self.getNeurons();
            nunits = self.nUnits;

            spikes = cell( 1,nunits );
            epochs = cell( 1,nunits );
            mask = cell( 1,nunits );
            unitID = uint8([neurons.ID]);
            unitCH = uint8([neurons.chanInd]);

            % loop over neurons, pull out spikes
            for n = 1:nunits
                [sp,~,ep,m] = neurons(n).getSpikes();
                spikes{n} = single( sp );
                mask{n} = sparse( double( m ) );
                epochs{n} = uint8( ep );
                neurons(n).meanWaveform = squeeze( mean( sp,2 ) );

                spObj = neurons(n).getChild( 'Spikes' );
                for sp = 1:numel( spObj )
                    spObj(sp).voltage = [];
                    spObj(sp).mask = [];
                end
            end

            % save the spikes
            saveFile = [spikesDir,filesep,self.filename,'_spikeWaveforms.mat'];
            save( saveFile,'spikes','epochs','mask','unitID','unitCH','-v7.3' );   
            
            % update the block
            self.update();
        end
        
    end % methods
end

