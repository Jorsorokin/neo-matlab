classdef ChannelIndex < Container
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
    %   channels - the indices from the list of all Electrodes (i.e.
    %              block.getChild( 'Electrode' ) )
    %   chanIDs - the actual Electrode IDs from the recording
    %   chanMap - a 1 x self.nElectrodes vector specifying any corrected ordering of electrodes 
    %   chanDistances - an m x self.nElectrodes vector specifying the m-dimensional 
    %                   distances of each electrode (i.e. (x,y), x-only, etc.)
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
    %   plotSpikes
    %   plotFeatures
    %   mergeNeurons
    %   splitNeurons
    %
    %       * see also methods in the Container class
    
    properties 
        chanIDs = [];
        channels = NaN;
        chanMap = [];
        chanDistances = [];
        chanIndNum
        nElectrodes = 0;
        nSignals = 0;
        nUnits = 0;
        name
    end      
        
    methods
        
        function self = ChannelIndex( chanIndNum )
            % initiate the ChannelIndex object
            
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
            for j = self.nElectrodes:-1:1
                signals = [electrodes(j).getChild( 'Signal' ),signals]; % negative indexing for performance
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
                
                % denoise the decomposition using single-resolution to better preserve spikes
                voltage = mswden( 'den',dec,'sqtwolog','sln','s' );

                % now add voltages back to their parent signals
                for j = 1:numel( theseSignals )
                    theseSignals(j).voltage(:) = voltage(:,j);
                end
            end
        end
        
        
        function detectSpikes( self,varargin )
            % detectSpikes( self,(thresh,artifact,masked_detection) )
            %
            % find spike waveforms from the analog signals contained in the child 
            % "Electrode" object. If no such child exists, end the function.
            %
            % Two different methods of spike detection are possible: summed multi-teager 
            % energy operator (MTEO) across electrodes, or masked threshold crossing. 
            % The former is useful for electrode configurations such as tetrodes,
            % where each electrode group is very well isolated from otheres (and thus
            % all electrodes in a ChannelIndex are assumed to record all detected spikes).
            % 
            % However in the case of dense multi-electrode arrays in which
            % electrode grouping is less obvious, it may be invalid to assume each electrode
            % recorded each spike. This later method uses the masked detection method, where
            % for each ith spike, any channel with a voltage deflection less than a predefined 
            % threshold are considered "masked". 
            %
            % The result of "detectSpikes()" is a "Spikes" object for each Epoch, 
            % as well as a "Neuron" object with ID = 0 that contains a reference to the 
            % current ChannelIndex. The "Spikes" object references the appropriate 
            % "Epoch" and "Neuron" objects, and the waveforms contained in the "Spikes" 
            % objects across epochs can be sorted through the current object via:
            % self.sortSpikes(). In addition, if "masked_detection" is true, a sparse matrix
            % indicating which channels contributed to which spike times will be stored into each 
            % Spikes object under the "mask" variable (i.e. Spikes.mask).
            %
            % optional inputs:
            %   thresh - the threshold for spike detection, as a multiple of the SD of the background
            %            (default = 4)
            %   artifact - the uV level to consider a spike as an artifact 
            %              (default = 1000)
            %   masked_detection - boolean indicating MTEO vs. masked spike detection
            %                      (default = false)
            
            % check inputs
            if nargin > 1 && ~isempty( varargin{1} )
                thresh = varargin{1};
            else
                thresh = 4; % * SD(voltage)
            end
            if nargin > 2 && ~isempty( varargin{2} )
                artifact = varargin{2};
            else
               artifact = 1000; % uV
            end
            if nargin > 3 && ~isempty( varargin{3} )
                masked_detection = varargin{3};
            else
                masked_detection = false; % MTEO vs. masked detection
            end

            % find all electrodes in this group
            electrodes = self.getChild( 'Electrode' );
            if isempty( electrodes )
                disp( 'No electrode(s) available' );
                return
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
            
            % create a Neuron class if one doesn't exist for the current channelindex
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
                
                % detect the spikes
                if ~masked_detection
                    [sptm,spsnip] = detectSpikes( volt,fs,thresh,1,artifact );
                    mask = [];
                else
                    maxPts = floor( 0.001 * fs ); % anything > 1ms is an artifact or overlapping spike
                    maxChan = [];

                    % get channel distance matrix to find # of channels within 150um of one another
                    if ~isempty( self.chanDistances )
                        distMat = pdist2( self.chanDistances,self.chanDistances );
                        maxChan = floor( mean( sum( distMat <= 150 ) ) );
                    end
                    [spsnip,sptm,mask] = double_flood_fill( bsxfun( @minus,volt,mean( volt,2 ) ),fs,...
                        self.chanMap,thresh/2,thresh,maxPts,maxChan,artifact );                
                end
                clear volt
                
                % create a "Spikes" object using the found spikes
                Sp = Spikes( single( sptm/fs ),single( spsnip ),fs );
                if nnz( mask ) / numel( mask ) <= 0.1
                    Sp.mask = sparse( mask ); % memory efficient if sparse enough
                else
                    Sp.mask = single( mask );
                end
                self.getChild( 'Neuron' ).addChild( Sp );
                
                % get the Epoch if it exists
                E = signals(find( thisEpoch,1 )).getParent( 'Epoch' );
                if ~isempty( E )
                    
                    % find the previous spikes associated with all Neuron objects in this ChannelIndex
                    allSpikes = E.getChild( 'Spikes' );
                    for neuron = 1:numel( N )
                        prevSpikeInds = find( ismember( allSpikes,N(neuron).getChild( 'Spikes' ) ) );                        
                        E.removeChild( 'Spikes',prevSpikeInds ); % remove previous spikes
                    end
                    E.addChild( Sp ); % add the new spikes
                end
            end
            
            fprintf( '\n' );
            % ==================================================================
            
            % updates
            clear N
            self.getParent( 'Block' ).update;
        end
        

        function sortSpikes( self,varargin )
            % sortSpikes( self,varargin );
            %
            % sort spikes referenced by the current ChannelIndex object. For each sorted
            % ID, create a new Neuron instance and store as a child under the
            % current ChannelIndex. Also create the appropriate "Spikes" 
            % object for each neuron ID (extracted from the
            % non-sorted "Spikes" objects)
            %
            % supply optional arguments as name-value pairs. 
            %
            % optional inputs:
            %   projMethod - the decomposition method (PCA,ICA,...)
            %   init - the initial cluster number 
            %   features - previously computed features
            %   level - the number of dims to keep (default = 5 )
            %   reject - the outlier threshold for HDBSCAN (default = 0.9)
            %   neuronID - vector or scalar indicating which neurons to sort
            %   sortMethod - method for clustering
            %   nNeighbors - # of points for nearest-neighbor calculation
            %   eps - epsilon radius (for DBSCAN / HDBSCAN / Spectral )
            %   minclustsize - minimum # of points in a cluster to keep
            %                  that cluster
            %   useGUI - boolean flag to use the sortTool GUI for sorting.
            %            If true, an instance of the sortTool GUI will be
            %            initiated, and the results of the sorting will be
            %            automatically merged with this ChannelIndex object
            %   
            % EX:
            %   self.spikeSort( 'method','pca','level',4,'reject',0.8 )
            %
            % refer to "sort_clusters.m" & "compute_spike_features.m" for
            % more details, as well as the sortTool GUI documentation
            
                
            % check the optional inputs
            p = check_inputs( varargin );
                           
            % update the block so that all neurons have different IDs
            block = self.getParent( 'Block' );
            if ~isempty( block )
                block.update();
            end
            
            % get all Spikes and the ChannelIndex objects associated with
            % this ChanIndex. Each child will be a "Spikes" associated with 
            % a different "Epoch" object
            neurons = self.getChild( 'Neuron' );
            if isempty( neurons )
                disp( 'Must detect spikes first!' );
                return
            end
            
            % check to make sure the supplied neuronIDs (if provided) exist
            % within this channelindex. If not, end the function.
            if isempty( p.neuronID )
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
            mask = [];
            for ind = p.neuronID
                [v,t,ep,m] = neurons.findobj( 'ID',ind ).getSpikes();
                
                % concatenate the spiketimes into one long vector
                t = reshape( t,1,numel( t ) );
                t(isnan(t)) = [];
                
                % add to sptm/spsnip/epochnum
                sptm = [sptm,t];
                spsnips = [spsnips,v];
                epochnum = [epochnum,ep];
                mask = [mask,m];
                clear t v ep m
            end           
            
            % Spike sorting
            % =============================================================
            
            % pull out spike features
            if isempty( p.features )
                [features,~,mapping] = compute_spike_features( permute( spsnips,[2,1,3] ),...
                                                               p.level,p.projMethod,...
                                                               mask,self.chanDistances,1 );
            else
                features = p.features;
                mapping = [];
            end
                                                    
            % sort the features
            [labels,model] = sort_clusters( features,p.init,p.sortMethod,...
                        'mask',mask,'eps',p.eps,'neighbors',p.nNeighbors,...
                        'minclustsize',p.minclustsize,'kernelWeighting',true,...
                        'outlierthresh',p.reject );
            
            % call an instance of the sortTool GUI if requested
            if p.useGUI
                [labels,features,R] = sortTool( 'data',single( spsnips ),'projection',single( features ),...
                                                'times',single( sptm ),'mask',single( mask ),...
                                                'trials',single( epochnum ),'labels',labels,...
                                                'location',self.chanDistances );
                
                if isstruct( R.mapping )
                    p.projMethod = R.projectMethod;
                    mapping = R.mapping;
                end
                if isfield( R,'sortModel' )
                   model = R.sortModel;
                else
                    model = [];
                end
            end
            
            % get vars from the model
            if ~isfield( model,'probabilities' )
                model.probabilities = [];
            end
            
            % =============================================================
            
            % create our new neurons / eliminate old ones
            self.create_new_neurons( spsnips,sptm,epochnum,p.neuronID,labels,...
                features,model.probabilities,p.projMethod,model,mapping,mask );
            
            % HELPER FUNCTION
            function p = check_inputs( inputs )
                pnames = {'projMethod','init','level','reject',...
                    'features','neuronID','sortMethod','nNeighbors',...
                    'eps','minclustsize','useGUI'};
                defaults = {'PCA',2,5,.9,[],[],'EM-GMM',5,0.1,5,false};

                p = inputParser;             
                for j = 1:numel(pnames)
                    p.addParameter( pnames{j},defaults{j} );
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
                disp( 'No Neuron objects found' );
                return
            end
            
            IDs = [neurons.ID];
            if ~all( ismember( neuronIDs,IDs) )
                disp( 'Some or all provided neuronIDs are not children of this ChannelIndex' );
                return
            end
            
            % create a new Neuron with the lowest ID of those being merged
            newNeuron = Neuron( min( neuronIDs ) );
            
            % loop over neuronIDs, pull out old Spike children
            oldSpikes = [];
            features = [];
            meanWaveform = [];
            for id = neuronIDs
                oldNeuron = neurons.findobj( 'ID',id );
                oldSpikes = [oldSpikes,oldNeuron.getChild( 'Spikes' )];
                meanWaveform = [meanWaveform,reshape( oldNeuron.meanWaveform,...
                                                      numel(oldNeuron.meanWaveform),1 )];
                if ~isnan( oldNeuron.features )
                    features = vertcat( features,oldNeuron.features );
                    newNeuron.featureMethod = oldNeuron.featureMethod;
                end
                oldNeuron.deleteSelf(); % removes references from parents and children
                clear oldNeuron 
            end
            
            % loop over the "spikes" objects. Find any that have the same "Epoch"
            % parent and merge together into a new Spikes object
            for j = 1:numel( oldSpikes )
                oldSpikes(j).epoch = oldSpikes(j).getParent('Epoch').epochNum;
            end
            
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
                        mask = [epochSpikes.mask]; 
                        [sptm,idx] = sort( sptm ); % sort smallest -> largest
                        if ~isempty( volt )
                            volt = volt(:,idx,:); % sort corresponding waveforms
                        end
                        if ~isempty( mask )
                            mask = mask(:,idx); % sort masked matrices
                        end
                        epoch = epochSpikes(1).getParent( 'Epoch' ); % get the actual Epoch
                        fs = epochSpikes(1).fs;
                        
                        % delete the old spikes of this epoch
                        for sp = 1:numel( epochSpikes )
                            epochSpikes(sp).deleteSelf();
                        end

                        % create a new Spikes object with the concatenated
                        % voltages/spike times. Add the corresponding epoch
                        newSpikes(counter) = Spikes( sptm,volt,fs );
                        newSpikes(counter).mask = mask;
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
                mask = [oldSpikes.mask];
                [sptm,idx] = sort( sptm );
                volt = volt(:,idx,:);
                if ~isempty( mask )
                    mask = mask(:,idx); % sort masked matrices
                end
                newSpikes = Spikes( sptm,volt,oldSpikes(1).fs );
                newSpikes.mask = mask;
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
            if ~isempty( meanWaveform )
                meanWaveform = mean( meanWaveform,2 );
                newNeuron.meanWaveform = reshape( meanWaveform,...
                                                  numel( meanWaveform )/self.nElectrodes,...
                                                  self.nElectrodes );
            end
            self.addChild( newNeuron );
            
            % update the block to re-establish consecutive neuron IDs
            self.getParent( 'Block' ).update();
        end
        
        
        function splitNeuron( self,oldID,labels )
            % splitNeuron( self,oldID,labels )
            %
            % splits a child Neuron object, specified by "oldID", into new 
            % neurons by re-assigning spikes according to the provided
            % label vector "labels"
            
            % get neuron children
            if self.nUnits == 0
                disp( 'ChannelIndex has no Neuron children!' );
                return
            end
            
            % check if the specified neuron exists
            oldNeuron = self.getChild( 'Neuron' ).findobj( 'ID',oldID );
            if isempty( oldNeuron )
                disp( 'No Neuron with specified ID found' );
                return
            end
            
            % check if # of labels provided does not match the # of spikes
            % associated with the old neuron
            if oldNeuron.nSpikes ~= numel( labels )
                disp( 'Number of provided labels does not match number of spikes in this Neuron' );
                return
            end
            
            % get the spike snips / times of this neuron
            [snips,sptimes,epochs,mask] = oldNeuron.getSpikes();
            sptimes = reshape( sptimes,numel( sptimes ),1 );
            sptimes(isnan( sptimes )) = [];
            
            % create the new neuron children given the labels / trials
            self.create_new_neurons( snips,sptimes,epochs,oldID,labels,...
                [],[],[],[],[],mask );
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
              
    end % public methods
    
    
    methods(Access = 'private')
        
        function create_new_neurons( self,spsnips,sptm,trials,prevNeuronID,newID,varargin )
            % merges / creates new Neuron children based on new labels
            
            % check inputs
            features = nan;
            prob = nan;
            featureMethod = nan;
            sortModel = nan;
            projMatrix = nan;
            mask = nan;
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
            if nargin > 10 && ~isempty( varargin{6} )
                mask = varargin{6};
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
            neurons = self.getChild( 'Neuron' );
            changedNeurons = neurons.findobj( 'ID',prevNeuronID );
            epoch = [];
            for j = 1:numel( changedNeurons )
                epoch = [epoch,changedNeurons(j).getPartner( 'Epoch','Spikes' )];
            end
            clear changedNeurons
            fs = epoch(1).getChild( 'Spikes',1 ).fs;
            
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
            %mask = cell( 1,numel( epoch ) );
            for ep = 1:numel( epoch )
                oldSpikes = epoch(ep).getChild( 'Spikes' );
                if ~isempty( oldSpikes )
                    oldSpikeIDs = [oldSpikes.unitID];
                    oldSpikeChans = zeros( 1,numel( oldSpikes ) );
                    for j = 1:numel( oldSpikes )
                        oldSpikeChans(j) = oldSpikes(j).getParent('Neuron').chanInd;
                    end
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
            remainingNeurons = self.getChild( 'Neuron' );
            if any( remainingNeurons )
                currentIDs = [remainingNeurons.ID];
            else
                currentIDs = [];
            end
            for i = uID
                counter = counter+1;
                thisID = (newID==i);
                thisNeuron = ismember( currentIDs,i );
                
                % now add the new Spikes based on the sorting if a
                % Neuron with ID == i does not exist
                if ~any( thisNeuron ) 
                   
                    % loop over epochs and create a new Spikes object. Add to
                    % it's appropriate epoch 
                   for ep = 1:numel( epoch )
                        thisEpoch = ismember( trials,ep );
                        IDX = thisEpoch & thisID;
                        Sp(ep) = Spikes( sptm(IDX),spsnips(:,IDX,:),fs );
                        if ~isnan( mask )
                            Sp(ep).mask = mask(:,IDX);
                        end
                        epoch(ep).addChild( Sp(ep) );
                   end

                    % create a new Neuron in the ChannelIndex and add the
                    % Spikes associated with different Epochs
                    newNeuron = Neuron( i );
                    newNeuron.addChild( Sp );

                    % add the sorting parameters specific to this neuron
                    if ~isnan( features )
                        newNeuron.features = features(thisID,:);
                    end
                    if any( ~isnan( prob ) ) 
                        newNeuron.probabilities = prob(thisID,:);
                    end
                    newNeuron.featureMethod = featureMethod;
                    newNeuron.sortModel = sortModel;
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
                        if ~isnan( mask )
                            prevSpikes(sp).mask = [prevSpikes(sp).mask,mask(:,IDX)];
                        end
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
            
            
