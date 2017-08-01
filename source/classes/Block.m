classdef Block < Container
    
    properties
        filename
        date
        condition
        filepath
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
            %
            % Examples:
            % 
            % % create a Block instance, add Epoch and ChannelIndex
            % % objects as children
            % file = 'myfile_trial1';
            % date = '12-May-2017';
            % condition = 'no stimulus';
            % filepath = '/Users/shared/EphysFiles';
            % block = Block( file,date,condition,filepath );
            % block.addChild( Epoch( 5,10,1 ) );
            % block.addChild( ChannelIndex( [1,2,3],[5,7,3],1 ) );
            % block.print(); % print a summary to the screen

            self.filename = filename;
            self.date = date;
            self.condition = condition;
            self.filepath = filepath;
        end
        

        function addChild(self,child)
            switch class( child )
                case {'Epoch'}
                    addChild@Container( self,child );
                    self.nEpochs = self.nEpochs + numel( child );
                case {'ChannelIndex'}
                    addChild@Container( self,child );
                    self.nChanInds = self.nChanInds + numel( child );
                otherwise 
                    error( 'Only ChannelIndex and Epoch objects are valid children' );
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
            fprintf('-----------------------------------------\n');
            fprintf( '%s (%s)\nRecorded on %s\n',...
                self.filename,self.condition,self.date );
            fprintf( '\n%i identified units over %i channel groups & %i epochs\n',...
                self.nUnits, self.nChanInds, self.nEpochs );
            fprintf('-----------------------------------------\n');
        end
        
        
        function update( self )
            % update( self )
            %
            % updates the properties contained within this block according
            % to its children. It searches through the Epoch and
            % ChannelIndex children (if any) and updates the number of
            % epochs, channels, neurons, and signals found.
            epoch = self.getChild( 'Epoch' );
            chanind = self.getChild( 'ChannelIndex' );
            if ~isempty( chanind )                
                % check for validity
                chanind(~isvalid( chanind )) = [];
                self.nSignals = sum( [chanind.nSignals] );
            end
            if ~isempty( epoch )
                % check for validity
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
            self.nChanInds = numel( chanind );
            
            % update all the Neuron IDs and the IDs associated with their
            % Spikes children
            neurons = self.getNeurons();
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
            eval( ['block' '=self'] ); % makes "block" reference "self"
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
                
    end
    
end

