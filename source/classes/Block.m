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
        chanMap = [];
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
            %   chanMap - a sparse block-diagonal matrix defining channel groups among all 
            %                channels (used for creating ChannelIndex objects).
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
            fprintf( '%i electrodes; %i signals; %i channel groups; %i epochs\n',...
                self.nElectrodes, self.nSignals, self.nChanInds, self.nEpochs );
            fprintf( '%i identified units\n', self.nUnits );
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
                chanind(~isvalid( chanind )) = [];
                for j = 1:numel( chanind )
                    electrodeID = [chanind(j).getChild( 'Electrode' )];
                    chanind(j).nElectrodes = numel( electrodeID );
                    chanind(j).chanIDs = electrodeID;
                end
            end
            self.nChanInds = numel( chanind );
            % =========================================
            
            % check validity & update: Electrode
            % ==================================
            if ~isempty( electrode )
                nSig = 0;
                for j = 1:numel( electrode )
                    nSig = nSig + electrode(j).nSignals;
                    for sig = 1:electrode(j).nSignals
                        signal = electrode(j).getChild( 'Signal',sig );
                        signal.chanInd = electrode(j).chanInd;
                    end
                end
                self.nSignals = nSig;
            else
                self.nSignals = 0;
            end
            self.nElectrodes = numel( electrode );
            % =========================================

            % check validity & update: Epoch
            % ==============================
            if ~isempty( epoch )
                epoch(~isvalid( epoch )) = [];
                
                % loop over spike objects, remove those that are not valid
                for ep = 1:numel( epoch )
                    spikes = epoch(ep).getChild( 'Spikes' );
                    if ~isempty( spikes )
                        spikes(~isvalid( spikes )) = [];
                        for sp = 1:numel( spikes )
                            if isempty( spikes(sp).getParent( 'Neuron' ) )
                                spikes(sp).deleteSelf();
                            end
                        end
                    end
                end
            end
            self.nEpochs = numel( epoch );
            % =========================================

            % check validity & update: Neuron
            % ==============================
            if ~isempty( neurons )
                neurons(~isvalid( neurons )) = [];
                for j = 1:numel( neurons )
                    if isempty( neurons(j).getParent( 'ChannelIndex' ) )
                        neurons(j).deleteSelf();
                    end
                end
                self.nUnits = numel( [neurons.ID] );
                if numel( self.nUnits ) > 0
                    for n = 1:self.nUnits
                        spikes = neurons(n).getChild( 'Spikes' );
                        neurons(n).nSpikes = sum( [spikes.nSpikes] );
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
            save( [outpath,filesep,self.filename,'_extractedData.mat'],'block' );
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

