classdef ChannelIndex < Container
    
    properties 
        channels
        chanIDs
        chanIndNum
        nChans
        nSignals = 0;
        nUnits = 0;
        sortModel
    end    
        
    methods
        
        function self = ChannelIndex( channels,channelIDs,chanIndNum )
            % self = ChannelIndex( channels,channelIDs,chanIndNum )
            %
            % Initiate an instance of the ChannelIndex class.
            % A ChannelIndex object contains a reference to a set of
            % channels and their actual IDs (not necessarily the 
            % the indices of the channels) that are intended as a group.
            % For instance, a tetrode can be organized as a ChannelIndex
            % object with four channels. 
            %
            % The ChannelIndex is intended to facilitate the organization
            % of recorded signals and spikes to easily extract such signals
            % separately for each group of channels across Epochs. It also
            % maintains a running count of the number of Neuron objects contained
            % within the group of channels over the recording session.
            %
            % Children:
            %   Signals
            %   Neurons
            %
            % Parents:
            %   Block
            %
            % Methods:
            %   detectSpikes
            %   sortSpikes
            %   plotSpikes
            %
            %       * see also methods in the Container class
            
            self.channels = channels;
            self.chanIDs = channelIDs;
            self.chanIndNum = chanIndNum;
            self.nChans = numel( channels );
            self.sortModel = struct( 'mu',[],'Sigma',[],'w',[] );
        end
        
        
        function addChild(self,child)
            % overloaded method to add specific attributes to the
            % ChannelIndex class
            if max( strcmp( class( child ),{'Signal','Neuron'} ) ) == 0
                error( 'Only Signal and Neuron objects are valid children' );
            end
            
            for j = 1:numel( child )
                switch class( child )
                    case {'Signal'}
                        if child(j).nSignals ~= self.nChans
                            error( 'Signal and ChanIndex objects must have equal # of channels' );
                        end
                        self.nSignals = self.nSignals + child(j).nSignals;
                    case {'Neuron'}
                        self.nUnits = self.nUnits + 1;
                end
                addChild@Container( self,child(j) );
                child(j).chanInd = self.chanIndNum;
                child(j).nChan = self.nChans;
            end  
        end
        
        
        function addParent( self,parent )
            % overloaded method to add specific attributes to the
            % ChannelIndex class
            switch class( parent )
                case {'Block','Epoch'}
                    addParent@Container( self,parent );
                otherwise
                    error( 'Only Block or Epoch objects are valid parents' );
            end
            if isa( parent,'Epoch' )
                self.epoch = parent.epochNum;
            end
        end
        
        
        function detectSpikes( self,thresh,artifact )
            % detectSpikes( self,thresh,artifact )
            %
            % find spike waveforms from the analog signals in the child 
            % "Signal" object. If no such child exists, end the function.
            %
            % stores the results of each detection into a "Spikes" object
            % and also creates a "Neuron" object that contains a reference
            % to the current ChannelIndex. The "Spikes" object references
            % the appropriate "Epoch" and "Neuron" objects, and the
            % waveforms contained in the "Spikes" objects across epochs can
            % be sorted through the current object via:
            % self.sortSpikes()
            
            % find all Signal object children of current channel index
            signals = self.getChild( 'Signal' );
            if isempty( signals )
                assert( 'No signals detected in current channel index' );
                return
            end
            
            % create a Neuron class if one doesn't exist for the current
            % ChannelIndex
            N = self.getChild( 'Neuron' );
            nSig = numel( signals ); % number of epochs
            if ~isempty( N )
                % remove previous Neurons from self
                self.removeChild( 'Neuron' );
            end
            clear N
            
            % add a new Neuron object to self, with a NaN identity
            self.addChild( Neuron( NaN ) );

            % find spikes using all chans combined in the "ChannelIndex" 
            for sig = 1:nSig
                
                % detect the spikes
                volt = filtfilt2( signals(sig).voltage,300,0,signals(sig).fs );
                [sptm,spsnip] = detectSpikes( volt,signals(sig).fs,...
                    thresh,1,artifact );
                
                % create a "Spikes" object using the found spikes
                Sp(sig) = Spikes( sptm/signals(sig).fs,spsnip,signals(sig).fs );
                
                % get the Epoch if it exists
                E = signals(sig).getParent( 'Epoch' );
                if ~isempty( E )
                    E.removeChild( 'Spikes' ); % remove previous spikes
                    E.addChild( Sp(sig) ); % add the new spikes
                end
            end
            
            % now add the Spikes object to the Neuron object
            self.getChild( 'Neuron' ).addChild( Sp ); 
        end
        

        function sortSpikes( self,varargin )
            % sortSpikes( self,(method,init,PCs,level,reject,search,neuronInds) );
            %
            % sort spikes referenced by the current ChannelIndex object. For each sorted
            % ID, create a new Neuron instance and store as a child under the
            % current ChannelIndex. Also create the appropriate "Spikes" 
            % object for each neuron ID (extracted from the
            % non-sorted "Spikes" objects)
            %
            % supply optional arguments as name-value pairs. 
            % EX:
            %   self.spikeSort( 'method','pca','level',4,'reject',0.1 )
            %
            % refer to "spikesort.m" for more details
                
            % check the optional inputs
            p = check_inputs( varargin );
                           
            % get all Spikes and the ChannelIndex objects associated wtih
            % this ChanIndex Each child will be a "Spikes" associated with 
            % a different "Epoch" object
            N = self.getChild( 'Neuron' );
            if isempty( N )
                error( 'Must detect spikes first!' );
            end
            if isnan( p.neuronInds ) || max( p.neuronInds ) > numel( N )
                p.neuronInds = 1:numel( N );
            end
            
            % pull out spike times/voltages and corresponding epochs
            sptm = [];
            spsnips = [];
            epochnum = [];
            for ind = p.neuronInds
                [v,t,ep] = N(ind).getSpikes();
                
                % concatenate the spiketimes into one long vector
                t = reshape( t,1,size(t,1)*size(t,2) );
                t(isnan(t)) = [];
                
                % add to sptm/spsnip/epochnum
                sptm = [sptm,t];
                spsnips = [spsnips,v];
                epochnum = [epochnum,ep];
                clear t v ep
            end           
            fs = N(1).getChild('Spikes',1).fs;
               
            % get the epochs for later reference
            epoch = self.getPartner( 'Epoch','Signal' ); % get the epochs
            
            % concatenate the spike waveforms accross channels
            allSnips = concatenateSpikes( spsnips );
            
            % now sort the spikes
            [id,Params] = spikesort( allSnips,fs,'init',p.init,'method',p.method,...
                'search',p.search,'level',p.level,'reject',p.reject,'PCs',p.PCs );
            
            % now extract the relevant IDs for each epoch, and change the
            % IDs according to the number of Neurons currently in existence
            % in the ChannelIndex object that were not re-sorted
            previousNeuronInds = ismember( 1:numel( N ),p.neuronInds );
            previousIDs = [N(previousNeuronInds).ID];
            nonSortedIDs = [N(~previousNeuronInds).ID];
            changeIDs = ismember( id,nonSortedIDs );
            
            % update the IDs to avoid double reference if we have sorted
            % in the past
            if any( changeIDs )
                id(changeIDs) = id(changeIDs) + max( currentIDs );
            end
            
            % remove the previous "Spikes" associated with a previously
            % sorted Neuron that is now resorted
            for ep = 1:numel( epoch )
                oldSpikes = epoch(ep).getChild( 'Spikes' );
                if ~isempty( oldSpikes )
                    oldSpikeIDs = [oldSpikes.unitID];
                    oldSpikeInds = (ismember( oldSpikeIDs,previousIDs )...
                        | isnan( oldSpikeIDs ));
                    if any( oldSpikeInds )
                        epoch(ep).removeChild( 'Spikes',oldSpikeInds );
                    end
                end
            end
            
            % now loop over the unique IDs and create a new Neuron object
            % and Spikes objects across Epoch
            for i = unique( id )
                thisID = id==i;
               
                % loop over epochs and create a new Spikes object. Add to
                % it's appropriate epoch
                for ep = 1:numel( epoch )
                    thisEpoch = epochnum == ep;
                                        
                    % now add the new Spikes based on the sorting
                    Sp(ep) = Spikes( sptm(thisEpoch & thisID),...
                        spsnips(:,thisEpoch & thisID,:),fs );
                    epoch(ep).addChild( Sp(ep) );
                end
                
                % create a new Neuron in the ChannelIndex and add the
                % Spikes associated with different Epochs
                N = Neuron( i );
                N.addChild( Sp );
                
                % add the sorting parameters specific to this neuron
                N.features = Params.features(thisID,:);
                N.probabilities = Params.prob(thisID,:);
                N.keptSpikes = Params.keptSpikes(thisID,:);
                N.featureMethod = Params.featureMethod;
                
                % add the neuron to the channelindex
                self.addChild( N );
            end
            
            % finally, update the sorting Params and remove old neurons          
            self.removeChild( 'Neuron',p.neuronInds );
            self.nUnits = numel( self.getChild( 'Neuron' ) );
            self.sortModel = Params.model; % how to update rather than replace?
            
            
            % HELPER FUNCTION
            function p = check_inputs( inputs )
                pnames = {'method','init','search','level','reject','PCs','neuronInds'};
                defaults = {'raw',2,0,5,.1,0,nan};
                options = {{'raw','ica','pca'},{1:50},{0,1},{1:50},{linspace(0,1,101)},{nan},{nan}};

                p = inputParser;             
                % loop over the rest of the optional inputs
                for j = 1:numel(pnames)
                    if ischar(options{j}{1})
                        p.addParameter( pnames{j},defaults{j},@(x) max(strcmp(x,options{j})) == 1 );
                    else
                        if strcmp( pnames{j},'init' ) || strcmp( pnames{j},'PCs' )
                            p.addParameter( pnames{j},defaults{j},@(x) ~isempty(x) );
                        else
                            p.addParameter( pnames{j},defaults{j},@(x) max(x == cell2mat(options{j})) == 1 );
                        end
                    end
                end
                p.parse(inputs{:});
                p = p.Results;
            end            
        end
        

        function plotSpikes( self,varargin )
            % plotSpikes( self,(epoch) )
            %
            % plot the spike waveforms contained in the current
            % ChannelIndex group. If no "Neuron" objects are associated
            % with this ChannelIndex, then end the function. Optionally
            % specify which epochs to plot. If left blank, will plot all
            % epochs.
            neurons = self.getChild( 'Neuron' );
            nNeurons = numel( neurons );
            if ~isempty( neurons )
                
                % check optional input
                if nargin > 1 && ~isempty( varargin{1} )
                    epochs = varargin{1};
                else
                    epochs = 1:min( cellfun( @numel,[neurons.children] ) );
                end
                
                % loop over neurons, get spikes over epochs, plot
                cmap = colormap( winter(nNeurons) );
                for n = 1:numel( neurons )
                    spikes = neurons(n).getChild( 'Spikes',epochs );
                    for ep = epochs
                        spikes(ep).plot( cmap(n,:) );
                    end
                end
                suptitle( ['Epochs ',num2str(epochs)] );
            end
        end
        
    end % methods
    
end
            
            