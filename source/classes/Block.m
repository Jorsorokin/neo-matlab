classdef Block < Container
    
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
            %
            %       * see also methods in the Container class

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

            % check for validity & update: ChannelIndex
            % =========================================
            if ~isempty( chanind )     
                self.nElectrodes = 0;
                chanind(~isvalid( chanind )) = [];
                for j = 1:numel( chanind )
<<<<<<< HEAD
                    chanElectrodes = chanind(j).getChild( 'Electrode' );
                    if ~isempty( chanElectrodes )
                        electrodeID = [chanElectrodes.electrodeNum];
                        chanind(j).nElectrodes = numel( electrodeID );
                        chanind(j).chanIDs = electrodeID; % the actual electrode IDs 
                    end
                    self.nElectrodes = self.nElectrodes + numel( chanElectrodes );
=======
                    electrodes = chanind(j).getChild( 'Electrode' );
                    electrodeID = [electrodes.electrodeNum];
                    chanind(j).nElectrodes = numel( electrodeID );
                    chanind(j).chanIDs = electrodeID; % the actual electrode IDs 
                    chanind(j).channels = 1:chanind(j).nElectrodes; % their location in the Block parent
>>>>>>> origin/master
                end
            end
            self.nChanInds = numel( chanind );
            % =========================================
            
            % check validity & update: Electrode & Signal
            % ==========================================
            if ~isempty( electrode )
                nSig = 0;
                electrode(~isvalid( electrode )) = [];
                for j = 1:numel( electrode )
                    nSig = nSig + electrode(j).nSignals;
                    chanind = electrode(j).getParent( 'ChannelIndex' );
                    if ~isempty( chanind )
                        electrode(j).chanInd = [electrode(j).getParent( 'ChannelIndex' ).chanIndNum]; % all parent ChannelIndex 
                        electrode(j).nChanInd = numel( chanind );
                    else
                        electrode(j).chanInd = [];
                        electrode(j).nChanInd = 0;
                    end
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
                self.nSignals = nSig;
            else
                self.nSignals = 0;
            end
            self.nElectrodes = numel( electrode );
            % ===========================================

            % check validity & update: Epoch
            % ==============================
            if ~isempty( epoch )
                epoch(~isvalid( epoch )) = [];
                
                % loop over spike objects, remove those that are not valid
                for ep = 1:numel( epoch )
                    spikes = epoch(ep).getChild( 'Spikes' );
                    signals = epoch(ep).getChild( 'Signal' );
                    if ~isempty( signals )
                        signals(~isvalid(signals)) = [];
                    end
                    epoch(ep).nSignals = numel( signals );
                    if ~isempty( spikes )
                        spikes(~isvalid( spikes )) = [];
                        for sp = 1:numel( spikes )
                            if isempty( spikes(sp).getParent( 'Neuron' ) )
                                spikes(sp).deleteSelf(); 
                            end
                        end
                        spikes(~isvalid( spikes )) = [];
                    end
                end
            end
            self.nEpochs = numel( epoch );
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
                for j = 1:numel( chanind )
                    neurons = [neurons,chanind(j).getChild( 'Neuron' )];
                end
            end
        end
                
    end % methods

end

