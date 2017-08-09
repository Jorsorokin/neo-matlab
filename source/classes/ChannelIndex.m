classdef ChannelIndex < Container
    
    properties 
        chanIDs = [];
        channels = NaN;
        chanIndNum
        nElectrodes = 0;
        nSignals = 0;
        nUnits = 0;
        name
    end      
        
    methods
        
        function self = ChannelIndex( chanIndNum )
            % self = ChannelIndex( chanIndNum )
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
            %   Electrode
            %   Neuron
            %
            % Parents:
            %   Block
            %
            % Properties:
            %   channels - the Electrodes from the list of all Electrodes (i.e.
            %              block.getChild( 'Electrode' ) )
            %   chanIDs - the actual Electrode IDs
            %   chanIndNum - the unique number for this ChannelIndex object
            %   nElectrodes - number of Electrode objects contained in this ChannelIndex
            %   nSignals - number of raw Signals provided to this ChannelIndex
            %   nUnits - number of Neurons provided to this ChannelIndex
            %   name - user-defined name specific to this ChannelIndex 
            %          (i.e. tetrode1, stereotrode3, etc.)
            %
            % Methods:
            %   getSignals
            %   waveDenoise
            %   detectSpikes
            %   sortSpikes
            %   sortGUI
            %   plotSpikes
            %   plotFeatures
            %   mergeNeurons
            %
            %       * see also methods in the Container class
            
            self.chanIndNum = chanIndNum;
        end
        
        
        function addChild(self,child)
            % overloaded method to add specific attributes to the
            % ChannelIndex class
            if max( strcmp( class( child ),{'Electrode','Neuron'} ) ) == 0
                error( 'Only Electrode and Neuron objects are valid children' );
            end
            
            for j = 1:numel( child )
                switch class( child )
                    case 'Electrode'
                        self.nElectrodes = self.nElectrodes + 1;
                        self.nSignals = self.nSignals + child(j).nSignals;
                        self.chanIDs(end+1) = child(j).electrodeNum;
                        child(j).chanInd(end+1) = self.chanIndNum;
                    case 'Neuron'
                        self.nUnits = self.nUnits + 1;
                        child(j).chanInd = self.chanIndNum;
                        child(j).nChan = self.nElectrodes;
                end
                addChild@Container( self,child(j) );
            end  
        end
        
        
        function addParent( self,parent )
            % overloaded method to add specific attributes to the
            % ChannelIndex class
            switch class( parent )
                case 'Block'
                    addParent@Container( self,parent );
                otherwise
                    error( 'Only Block or Epoch objects are valid parents' );
            end
        end


        function signals = getSignals( self )
            % returns all Signal objects contained within all Electrodes
            % belonging to this ChannelIndex object
            signals = [];
            if self.nSignals == 0
                disp( 'No signals available' );
                return;
            end

            electrodes = self.getChild( 'Electrode' );
            for j = 1:self.nElectrodes
                signals = [signals,electrodes(j).getChild( 'Signal' )];
            end
        end


        function waveDenoise( self,wLevel,wType )
            % waveDenoise( self,wLevel,wType )
            %
            % denoises the signals contained within all child Electrodes
            % using multi-signal wavelet denoising.
            % 
            % wLevel equals the wavelet decomposition level desired (lower
            % = less smoothing), and wType equals the wavelet to use.
            %
            % Note, this function will not work if the voltages contained within 
            % the various "Signal" children across Electrodes have different 
            % sampling rates and/or number of points. Thus, for each Epoch,
            % one must ensure the voltages have not been filtered differently 
            % for the various electrodes.
            
            % get all Signal objects associated with this channelindex
            signals = self.getSignals();
            epochs = [signals.epoch]; 
            uniqueEpochs = unique( epochs );

            % loop over unique epochs
            for ep = uniqueEpochs

                % pull out the appropriate signals
                theseSignals = signals.findObj( 'Epoch',ep ); 
                voltage = [theseSignals.voltage];

                % get the multi-signal wavelet decomposition
                dec = mdwtdec( 'c',voltage,wLevel,wType );
                
                % denoise the decomposition using multi-resolution
                voltage = mswden( 'den',dec,'rigrsure','mln','s' );

                % now add voltages back to their parent signals
                for j = 1:numel( theseSignals )
                    theseSignals(j).voltage = voltage(:,j);
                end
            end
        end
        
        
        function detectSpikes( self,thresh,artifact )
            % detectSpikes( self,thresh,artifact )
            %
            % find spike waveforms from the analog signals contained in the child 
            % "Electrode" object. If no such child exists, end the function.
            %
            % stores the results of each detection into a "Spikes" object
            % and also creates a "Neuron" object that contains a reference
            % to the current ChannelIndex. The "Spikes" object references
            % the appropriate "Epoch" and "Neuron" objects, and the
            % waveforms contained in the "Spikes" objects across epochs can
            % be sorted through the current object via:
            % self.sortSpikes()
            
            % find all electrodes in this group
            electrodes = self.getChild( 'Electrode' );
            if isempty( electrodes )
                disp( 'No electrode(s) available' );
                return;
            end

            % pull out the actual Signals & the total Epochs/Electrodes
            signals = [];
            for j = 1:self.nElectrodes
                signals = [signals,electrodes(j).getChild( 'Signal' )];
            end

            % get their epochs 
            epochNums = [signals.epoch]; 
            if any( isnan( epochNums ) ) % i.e. no Epoch class instances
                epochNums = 1:numel( signals );
            end
            uniqueEpochs = unique( epochNums );
            
            % create a Neuron class if one doesn't exist for the current
            N = self.getChild( 'Neuron' );
            if ~isempty( N )
                self.removeChild( 'Neuron' );
            end
            self.addChild( Neuron( 0 ) ); % zero-ID neuron indicates pre-sorting / non-sorted

            % Spike Detection
            % ==================================================================
            fprintf( 'Detecting spikes across epochs' )
            for ep = 1:numel( uniqueEpochs )
                fprintf( '.' );
                
                % get the voltage traces associated with the current epoch & filter for spike
                thisEpoch = ismember( epochNums,uniqueEpochs( ep ));
                fs = signals(find( thisEpoch,1 )).fs;
                volt = filtfilt2( [signals(thisEpoch).voltage],300,0,fs );
                [n,m] = size( volt );
                volt = reshape( smooth( medfilt1( volt,3 ) ),n,m ); % removes spot noise

                % detect the spikes
                [sptm,spsnip] = detectSpikes( volt,fs,thresh,1,artifact );
                
                % create a "Spikes" object using the found spikes
                Sp(ep) = Spikes( sptm/fs,spsnip,fs );
                
                % get the Epoch if it exists
                E = signals(find( thisEpoch,1 )).getParent( 'Epoch' );
                if ~isempty( E )
                    
                    % find the previous spikes associated with all Neuron objects in this ChannelIndex
                    allSpikes = E.getChild( 'Spikes' );
                    for neuron = 1:numel( N )
                        prevSpikeInds = find( ismember( allSpikes,N(neuron).getChild( 'Spikes' ) ) );                        
                        E.removeChild( 'Spikes',prevSpikeInds ); % remove previous spikes
                    end
                    E.addChild( Sp(ep) ); % add the new spikes
                end
            end
            fprintf( '\n' );
            % ==================================================================
            
            % now add the spikes to the Neuron parent
            self.getChild( 'Neuron' ).addChild( Sp ); 
        end
        

        function sortSpikes( self,varargin )
            % sortSpikes( self,(method,init,PCs,features,level,reject,search,neuronID) );
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
                           
            % update the block so that all neurons have different IDs
            block = self.getParent( 'Block' );
            block.update();
            
            % get all Spikes and the ChannelIndex objects associated with
            % this ChanIndex. Each child will be a "Spikes" associated with 
            % a different "Epoch" object
            neurons = self.getChild( 'Neuron' );
            if isempty( neurons )
                error( 'Must detect spikes first!' );
            end
            
            % check to make sure the supplied neuronIDs (if provided) exist
            % within this channelindex. If not, end the function.
            if isnan( p.neuronID )
                p.neuronID = [neurons.ID];
            else
                if ~any( ismember( [neurons.ID],p.neuronID ) )
                    error( 'Provided neuron IDs do not exist within this ChannelIndex ' );
                end
            end
            
            % pull out spike times/voltages and corresponding epochs
            sptm = [];
            spsnips = [];
            epochnum = [];
            for ind = p.neuronID
                [v,t,ep] = neurons.findobj('ID',ind).getSpikes();
                
                % concatenate the spiketimes into one long vector
                t = reshape( t,1,numel( t ) );
                t(isnan(t)) = [];
                
                % add to sptm/spsnip/epochnum
                sptm = [sptm,t];
                spsnips = [spsnips,v];
                epochnum = [epochnum,ep];
                clear t v ep
            end           
            fs = neurons(1).getChild('Spikes',1).fs;
            
            % now sort the spikes
            [id,Params] = spikesort( concatenateSpikes( spsnips ),fs,...
                'init',p.init,'method',p.method,...
                'search',p.search,'level',p.level,'reject',p.reject,...
                'PCs',p.PCs,'features',p.features );
            if isfield( Params,'PCs' ) && numel( Params.PCs ) > 1
                self.sortPCs = Params.PCs;
            elseif isfield( Params,'ICs' )
                self.sortICs = Params.ICs;
            end
                   
            % create our new neurons / eliminate old ones
            self.create_new_neurons( spsnips,sptm,epochnum,p.neuronID,id,...
                Params.features,Params.prob,Params.featMethod,...
                Params.model,Params.mapping );
            
            % update the sorting Params 
            %self.sortModel = Params.model; % how to update rather than replace? 
            %self.projMatrix = Params.mapping;
            
            % HELPER FUNCTION
            function p = check_inputs( inputs )
                pnames = {'method','init','search','level','reject','PCs','features','neuronID'};
                defaults = {'raw',2,0,5,.1,0,nan,nan};
                options = {{'raw','ica','pca'},{1:50},{0,1},{1:50},{linspace(0,1,101)},{nan},{nan},{1:100}};

                p = inputParser;             
                % loop over the rest of the optional inputs
                for j = 1:numel(pnames)
                    if ischar(options{j}{1})
                        p.addParameter( pnames{j},defaults{j},@(x) max(strcmp(x,options{j})) == 1 );
                    else
                        if strcmp( pnames{j},'init' ) || strcmp( pnames{j},'PCs' ) || strcmp( pnames{j},'features' )
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
        
        
        function mergeNeurons( self,neuronIDs )
            % mergeNeurons( self,neuronIDs )
            %
            % merge the Neuron objects identified by the supplied neuronIDs
            % into one Neuron object. This is useful for undoing the
            % results of running "self.sortSpikes()", combining
            % noise-spikes, etc. The parent Block will update so that
            % Neuron IDs and their Spike children are correctly
            % re-established. 
            
            % check the neuronIDs, make sure they exist in this
            % channelindex object
            neurons = self.getChild( 'Neuron' );
            if isempty( neurons )
                error( 'No Neuron objects found' );
            end
            IDs = [neurons.ID];
            if ~all( ismember( neuronIDs,IDs) )
                error( 'Some or all provided neuronIDs are not children of this ChannelIndex' );
            end
                        
            % get the IDs of all neurons within this Block
            allIDs = [self.getParent( 'Block' ).getNeurons().ID];
            
            % create a new Neuron with a temporary ID
            newNeuron = Neuron( allIDs(end)+1 );
            
            % loop over neuronIDs, pull out old Spike children
            oldSpikes = [];
            features = [];
            for id = neuronIDs
                oldInd = find( IDs==id,1 );
                oldNeuron = neurons(oldInd);
                oldSpikes = [oldSpikes,oldNeuron.getChild( 'Spikes' )];
                if ~isnan( oldNeuron.features )
                    features = vertcat( features,oldNeuron.features );
                    newNeuron.featureMethod = oldNeuron.featureMethod;
                end
                oldNeuron.deleteSelf(); % removes references from parents and children
                clear oldNeuron 
            end
            
            % loop over the "spikes" objects. Find any that have the same "Epoch"
            % parent and merge the times/voltages together into a new
            % Spikes object
            allEpochs = ( [oldSpikes.epoch] );
            if any( allEpochs )
                counter = 0;
                possibleEpochs = unique( allEpochs );
                for ep = possibleEpochs
                    counter = counter + 1;
                    epochSpikes = oldSpikes( ismember( allEpochs,ep ) );
                    if numel( epochSpikes ) > 1

                        % get the spike times/voltages and sort
                        sptm = [epochSpikes.times];
                        volt = [epochSpikes.voltage];
                        [sptm,idx] = sort( sptm ); % sort smallest -> largest
                        volt = volt(:,idx,:); % sort corresponding waveforms
                        epoch = epochSpikes(1).getParent( 'Epoch' ); % get the actual Epoch
                        fs = epochSpikes(1).fs;
                        
                        % delete the old spikes of this epoch
                        for sp = 1:numel( epochSpikes )
                            epochSpikes(sp).deleteSelf();
                        end

                        % create a new Spikes object with the concatenated
                        % voltages/spike times. Add the corresponding epoc
                        newSpikes(counter) = Spikes( sptm,volt,fs );
                        if ~isempty( epoch )
                            newSpikes(counter).addParent( epoch );
                        end
                        
                        clear epoch
                    else
                        % just add the current spikes object to the
                        % newSpikes variable
                        newSpikes(counter) = epochSpikes;
                    end
                end
            else
                % if no epochs exist, just merge the spikes of all the
                % oldSpikes objects and create a new Spikes object
                sptm = [oldSpikes.times];
                volt = [oldSpikes.voltage];
                [sptm,idx] = sort( sptm );
                volt = volt(:,idx,:);
                newSpikes = Spikes( sptm,volt,oldSpikes(1).fs );
                for sp = 1:numel( oldSpikes )
                    oldSpikes(sp).deleteSelf();
                end
            end
            
            % now add the spikes to the new Neuron, and the Neuron to the
            % ChannelIndex
            newNeuron.addChild( newSpikes );
            if ~isempty( features )
                newNeuron.features = features;
            end
            self.addChild( newNeuron );
            
            % update the block to re-establish consecutive neuron IDs
            self.getParent( 'Block' ).update();
        end

        
        function plotSpikes( self,varargin )
            % plotSpikes( self,(epoch,neuronIDs) )
            %
            % plot the spike waveforms contained in the current
            % ChannelIndex group. If no "Neuron" objects are associated
            % with this ChannelIndex, then end the function. Optionally
            % specify which epochs and neurons to plot (by IDs). 
            % If left blank, will plot all epochs / neurons
            neurons = self.getChild( 'Neuron' );
            
            if ~isempty( neurons )
                
                % check optional input
                if nargin > 1 && ~isempty( varargin{1} )
                    epochs = varargin{1};
                else
                    epochs = 1:min( cellfun( @numel,[neurons.children] ) );
                end
                if nargin > 2 && ~isempty( varargin{2} )
                    neuronIDs = varargin{2};
                    neurons = neurons(ismember( [neurons.ID],neuronIDs ));
                end

                % loop over neurons, get spikes/features over epochs, plot
                nNeurons = numel( neurons );
                cmap = colormap( jet(nNeurons) );
                for n = 1:nNeurons
                    spikes = neurons(n).getChild( 'Spikes',epochs );
                    for ep = epochs
                        spikes(ep).plot( cmap(n,:) );
                    end
                end
                suptitle( ['Epochs ',num2str(epochs)] );
            end
        end
        
        
        function plotFeatures( self,varargin )
            % plotFeatures( self,(neuronIDs) )
            %
            % plot the features contained within the Neurons of this
            % ChannelIndex (if sorting has been done) using "gplotmatrix",
            % color-coded by the neuron ID. Can optionally specify which neurons
            % (by their IDs, not indices) to plot.
            neurons = self.getChild( 'Neuron' );
            
            if ~isempty( neurons )
                
                % check input
                if nargin > 1 && ~isempty( varargin{1} )
                    neuronIDs = varargin{1};
                    neurons = neurons(ismember( [neurons.ID],neuronIDs ));
                end
                
                % determine if any spikes with this ID exist
                nNeurons = numel( neurons );
                if nNeurons == 0
                    error( 'No Neurons with provided IDs belong to this ChannelIndex' );
                end
                
                % loop over neurons, get features
                nSpikes = [neurons.nSpikes];
                features = nan( sum( nSpikes ),size( self.sortModel.mu,1 ) );
                id = zeros( size( features,1 ),1 );
                counter = 0;
                for n = 1:nNeurons
                    inds = counter+1:counter+nSpikes(n);
                    features(inds,:) = neurons(n).features;
                    id(inds) = neurons(n).ID;
                    counter = counter + nSpikes(n);
                end
                
                % now plot the features 
                figure;
                for j = 1:size( features,2 )
                    axlabel{j} = sprintf( 'feature %i',j );
                end
                [~,ax] = gplotmatrix( features,[],id,[],[],3,[],[],axlabel );
                legend(ax(1,end),strsplit( num2str( [neurons.ID] ),' ' ));
                if ~isempty( self.name )
                    suptitle( sprintf( '%s spike features',self.name ) );
                else
                    suptitle( sprintf( 'ChanInd #%i spike features',self.chanIndNum ) );
                end
            end
        end


        function sortGUI( self,varargin )
            % sortGUI( self,(neuronIDs) )
            %
            % creates an instance of the sortTool.fig GUI for visualizing
            % sorting results in better detail, and providing flexibility 
            % regarding manual clustering and projection methods.
            %
            % By default, all spikes from all Neuron children of the current container
            % will be fed into the sortTool GUI, however one can specify a 
            % subset of neurons via the optional arugment "neuronID", which will 
            % only provide the spike waveforms from Neurons with those IDs (for further sorting).
            neurons = self.getChild( 'Neuron' );
            if isempty( neurons )
                disp( 'Must detect spikes first!' );
                return
            end

            % check input
            if nargin > 1 && ~isempty( varargin{1} )
                neuronIDs = varargin{1};
                neurons = neurons(ismember( [neurons.ID],neuronIDs ));
            else
                neuronIDs = [neurons.ID];
            end
            
            % determine if any spikes with this ID exist
            nNeurons = numel( neurons );
            if nNeurons == 0
                disp( 'No Neurons with provided IDs belong to this ChannelIndex' );
                return
            end

            % determine if any electrodes added
            if self.nElectrodes == 0
                disp( 'Mismatch between # of electrodes among data objects' );
                disp( 'Check your data heirarhcy. Try running "block.update()"' );
                return
            end
            
            % get the spikes from these neurons
            nSpikes = [neurons.nSpikes];
            nPts = size( neurons(1).getChild( 'Spikes',1 ).voltage,1 );
            spsnips = nan( nPts,sum( nSpikes ),self.nElectrodes );   % nPt x nSp x nCh matrix (spike waveforms)
            sptimes = nan( sum( nSpikes ),1 );                       % nSp x 1 vector (spike times)
            epoch = sptimes;                                         % nSp x 1 vector (parent epoch)
            id = sptimes;                                            % nSp x 1 vector (class labels) 
            counter = 0;
            for n = 1:nNeurons
                inds = counter+1:counter+nSpikes(n);
                [spsnips(:,inds,:),times,epoch(inds)] = neurons(n).getSpikes();
                times = reshape( times,numel(times),1 );
                times( isnan( times ) ) = [];
                sptimes(inds) = times;
                id(inds) = neurons(n).ID;
                counter = counter + nSpikes(n);
                clear times
            end

            % initiate the sortTool GUI with the provided waveforms and labels
            [newID,features,model] = sortTool( 'data',concatenateSpikes( spsnips ),...
                'times',sptimes','trials',epoch','labels',id' );

            % create new neurons / eliminate old
            if max( newID ) > 0 
               if numel( newID ) ~= numel( sptimes ) || any( ~ismember( newID,id ) )
                    self.create_new_neurons( spsnips(:,model.keptPts,:),sptimes(model.keptPts),...
                        epoch(model.keptPts),neuronIDs,newID,features,...
                        model.probabilities,model.projectMethod,...
                        model.sortModel,model.W );
                   % self.projMatrix = model.W;
                   % self.sortModel = model.sortModel;
               end
            end
        end
        
    end % public methods
    
    
    methods(Access = 'private');
        
        function create_new_neurons( self,spsnips,sptm,trials,prevNeuronID,newID,varargin )
            % merges / creates new Neuron children based on new labels
            
            % check inputs
            features = nan;
            prob = nan;
            featureMethod = nan;
            sortModel = nan;
            projMatrix = nan;
            if nargin > 5 && ~isempty( varargin{1} )
                features = varargin{1};
            end
            if nargin > 6 && ~isempty( varargin{2} )
                prob = varargin{2};
            end
            if nargin > 7 && ~isempty( varargin{3} )
                featureMethod = varargin{3};
            end
            if nargin > 8 && ~isempty( varargin{4} )
                sortModel = varargin{4};
            end
            if nargin > 9 && ~isempty( varargin{5} )
                projMatrix = varargin{5};
            end
            
            % re format if necessary
            if ~isrow( sptm )
                sptm = sptm';
            end
            if ~isrow( trials )
                trials = trials';
            end
            if ~isrow( newID )
                newID = newID';
            end
            
            % get all Neuron children and Epoch partners
            %allNeurons = self.getParent( 'Block' ).getNeurons();
            neurons = self.getChild( 'Neuron' );
            epoch = self.getPartner( 'Epoch','Signal' );
            fs = epoch(1).getChild( 'Spikes',1 ).fs;
            
            % ============================================================
            % create new neurons & associated spikes
            % ============================================================
            % (a) get all of the previous IDs, and any ID's == 0
            allIDs = [neurons.ID];
            zeroID = (newID==0); % find those with IDs == 0
            
            % (b) get the indices of the neurons referring to current ones
            previousNeuronInds = ismember( allIDs,prevNeuronID );
            nonSortedIDs = allIDs(~previousNeuronInds);
            nonSortedIDs(nonSortedIDs==0) = []; % remove the ID = 0 so that
                                                % we can  simply add to it
            
            % (c) get indices of which IDs we wish to change, ignoring
            % those with IDs == 0
            changeIDs = ismember( newID,nonSortedIDs ); 
            oldIDs = unique( newID(changeIDs) );
            totalNum = numel( allIDs ) + numel( unique( newID ) );
            takenID = false( 1,totalNum ); % update this as we use up IDs
            
            % (d) update the IDs to avoid double reference if we have sorted
            % in the past
            if any( ~isnan( nonSortedIDs ) )
                usedIDs = newID(~changeIDs & ~zeroID);
                if ~isempty( usedIDs )
                    takenID( [nonSortedIDs,unique( usedIDs )] ) = true;
                end
            end
            if any( changeIDs )
                for i = 1:numel( oldIDs ) % loop over IDs that we need to change to avoid overlap
                    currentID = (newID==oldIDs(i));
                    new = find( ~takenID,1 ); % find the first ID in all possible ones that is available to use
                    newID(currentID) = new;
                    takenID(new) = true; % no longer available for the next Neuron
                end
            end
            clear takenID currentID changeIDs nonSortedIDs usedIDs zeroID
            
            % (e) remove the "Spikes" associated with a previously
            % sorted Neuron that is now resorted
            for ep = 1:numel( epoch )
                oldSpikes = epoch(ep).getChild( 'Spikes' );
                if ~isempty( oldSpikes )
                    oldSpikeIDs = [oldSpikes.unitID];
                    oldSpikeChans = zeros( 1,numel( oldSpikes ) );
                    for j = 1:numel( oldSpikes )
                        oldSpikeChans(j) = oldSpikes(j).getParent('Neuron').chanInd;
                    end
                    %oldSpikeChans = [oldSpikes.getParent( 'Neuron' ).chanInd];
                    oldSpikeInds = (ismember( oldSpikeIDs,prevNeuronID )...
                    & ismember( oldSpikeChans,self.chanIndNum ))...
                    | isnan( oldSpikeIDs );
                    if any( oldSpikeInds )
                        epoch(ep).removeChild( 'Spikes',find( oldSpikeInds ) );
                    end
                end
            end
            
            % remove the old neurons from this channel index
            self.removeChild( 'Neuron',find(  previousNeuronInds ) );

            % (f) now loop over the unique IDs and create a new Neuron object
            % and Spikes objects across Epochs
            uID = unique( newID );
            counter = 0;
            for i = uID
                counter = counter+1;
                thisID = (newID==i);
                currentIDs = [self.getChild( 'Neuron' ).ID];
                thisNeuron = ismember( currentIDs,i );
                
                % now add the new Spikes based on the sorting if a
                % Neuron with ID == i does not exist
                if ~any( thisNeuron ) 
                   
                    % loop over epochs and create a new Spikes object. Add to
                    % it's appropriate epoch 
                   for ep = 1:numel( epoch )
                        thisEpoch = ismember( trials,ep );
                        Sp(ep) = Spikes( sptm(thisEpoch & thisID),...
                            spsnips(:,thisEpoch & thisID,:),fs );
                        epoch(ep).addChild( Sp(ep) );
                   end

                    % create a new Neuron in the ChannelIndex and add the
                    % Spikes associated with different Epochs
                    newNeuron = Neuron( i );
                    newNeuron.addChild( Sp );

                    % add the sorting parameters specific to this neuron
                    newNeuron.features = features(thisID,:);
                    if any( ~isnan( prob ) ) % <-- fix this after fixing sortTool manual clust
                        newNeuron.probabilities = prob(thisID,:);
                    end
                    newNeuron.featureMethod = featureMethod;
                    newNeuron.sortModel.mu = sortModel.mu(:,counter);
                    newNeuron.sortModel.sigma = sortModel.Sigma(:,:,counter);
                    newNeuron.sortModel.w = sortModel.w(counter);
                    newNeuron.projMatrix = projMatrix;

                    % add the neuron to the channelindex
                    self.addChild( newNeuron );
                else
                    % if a Neuron with this ID already exists (id == 0), 
                    % get the spikes and just add them to the appropriate
                    % spikes objects
                    counter = counter-1;
                    neuron = self.getChild( 'Neuron',find( thisNeuron ) );
                    prevSpikes = neuron.getChild( 'Spikes' );
                    
                    % loop over spikes children
                    for sp = 1:numel( prevSpikes )
                        thisEpoch = ismember( trials,prevSpikes(sp).epoch );
                        IDX = thisEpoch & thisID;
                        prevSpikes(sp).voltage = [prevSpikes(sp).voltage,spsnips(:,IDX,:)];
                        prevSpikes(sp).times = [prevSpikes(sp).times,sptm(IDX)];
                        prevSpikes(sp).nSpikes = prevSpikes(sp) + sum( IDX );
                    end
                    
                    % remove probabilities, features, and kept spikes from
                    % this Neuron, since these variables change between
                    % different sorting runtimes and it makes no sense to,
                    % say, concatenate features taken from PCA on one sort
                    % runtime, and those from NPE or ICA, etc. from another
                    neuron.features = nan;
                    neuron.probabilities = nan;
                    neuron.sortModel = nan;
                    neuron.projMatrix = nan;
                end
            end
            
            % update self
            self.nUnits = numel( self.getChild( 'Neuron' ) );
            self.getParent( 'Block' ).update();
        end
        
    end % private methods

end
            
            