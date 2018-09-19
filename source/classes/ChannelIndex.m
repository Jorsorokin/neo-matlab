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
    %   sortSpikes_using_clusterModel
    %   sortSpikes_using_templates
    %   plotSpikes
    %   plotFeatures
    %   mergeNeurons
    %   splitNeuron
    %
    %       * see also methods in the Container class
    
    properties 
        chanIDs = [];
        channels = NaN;
        chanMap = [];
        chanDistances = [];
        chanIndNum
        noiseCov = [];
        noise = [];
        whiteningMatrix = [];
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
            if ~isa( child,'Electrode' ) && ~isa( child,'Neuron' )
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
            assert( isa( parent,'Block' ),'Only Block or Epoch objects are valid parents' );
            addParent@Container( self,parent );
        end


        function signals = getSignals( self )
            % returns all Signal objects contained within all Electrodes
            % belonging to this ChannelIndex object
            if self.nSignals == 0
                disp( 'No signals available' );
                return
            end

            signals = [self.getChild( 'Electrode' ).children];
            signals = [signals{:}];
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
            % detectSpikes( self,(thresh,artifact,masked_detection,whiten,computeNoiseCov) )
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
            %   thresh:
            %       the threshold for spike detection, as a multiple of the SD of the background
            %           Default = 4
            %
            %   artifact:
            %       the uV level to consider a spike as an artifact 
            %           Default = 1000
            %
            %   masked_detection: 
            %       boolean indicating MTEO vs. masked spike detection
            %           Default = false
            %
            %   whiten:
            %       boolean flag... Data will be whitened across channels
            %       using self.whiteningMatrix (can be pre-provided or will
            %       be computed using a low-RMS segment of data to estimate
            %       the the noise covariance)
            %           Default = false
            
            % check inputs
            p = check_inputs( varargin );

            % find all electrodes in this group
            electrodes = self.getChild( 'Electrode' );
            if isempty( electrodes )
                disp( 'No electrode(s) available' );
                return
            end

            % pull out the actual Signals & the total Epochs/Electrodes
            epochs = self.getSibling( 'Epoch','Block' );
            nEpochs = numel( epochs );
            
            % create a Neuron class if one doesn't exist for the current channelindex
            prevNeuron = self.getChild( 'Neuron' );
            if ~isempty( prevNeuron )
                self.removeChild( 'Neuron' );
            end
            newNeuron = Neuron( 0 ); % zero-ID neuron indicates pre-sorting / non-sorted

            % variables for detection
            fs = epochs(1).getChild( 'Signal',1 ).fs; % assumes signal sampling rates are equal
            chanDist = self.chanDistances;
            map = self.chanMap;
            maxPts = floor( 0.002 * fs );
            if ~isempty( chanDist )
                distMat = pdist2( chanDist,chanDist );
                maxChan = floor( mean( sum( distMat <= 250 ) ) );
            else
                maxChan  = self.nElectrodes;
            end
            
            % get whitening matrix 
            if p.whiten && isempty( self.whiteningMatrix )
                volt = filtfilt2( [epochs(1).getChild('Signal').voltage],300,6000,fs ); 
                
                if size( volt,2 ) > 16
                    volt = volt - median( volt,2 );
                end
                
                if ~p.masked_detection
                    sptm = detectSpikes( volt,fs,p.thresh,1,p.artifact );
                else
                    [~,sptm] = double_flood_fill( volt,fs,...
                        'chanMap',map,'lowThresh',p.thresh/2,'highThresh',p.thresh,...
                        'maxPts',maxPts,'maxChans',maxChan,'artifact',p.artifact );   
                end

                % find regions in recording without spikes
                loc = find( diff( sptm/fs*1000 ) >= 10 );
                C = zeros( self.nElectrodes,class( volt ) );
                for k = 1:numel( loc )
                    C = C + cov( volt(round( sptm(loc(k)):sptm(loc(k)+1) ),:) );
                end
                
                self.noiseCov = C / k;
                self.whiteningMatrix = (self.noiseCov + eye(self.nElectrodes)*1e-10)^-0.5;
                clear sptm C volt loc
            end
                
            % Spike Detection
            % ==================================================================
            fprintf( 'Detecting spikes' )
            W = self.whiteningMatrix';
            chNum = self.chanIndNum;
            whiten = p.whiten;
            useMask = p.masked_detection;
            for ep = 1:nEpochs
                
                fprintf( '.' );
                
                % get the voltage traces associated with the current epoch & filter for spike
                volt = filtfilt2( [epochs(ep).getChild( 'Signal' ).voltage],300,6000,fs );   
                
                if size( volt,2 ) > 16
                    volt = volt - median( volt,2 );
                end

                if whiten
                    volt = volt * W; % performs the actual whitening
                end
                
                % detect the spikes
                if ~useMask
                    [sptm,spsnip] = detectSpikes( volt,fs,p.thresh,1,p.artifact );
                    mask = [];
                else
                    [spsnip,sptm,mask] = double_flood_fill( volt,fs,...
                        'chanMap',map,'lowThresh',p.thresh/2,'highThresh',p.thresh,...
                        'maxPts',maxPts,'maxChans',maxChan,'artifact',p.artifact );                
                end
            
                % create a "Spikes" object using the found spikes
                Sp = Spikes( single( sptm/fs ),single( spsnip ),fs );
                Sp.mask = sparse( mask ); 
                
                % find the previous spikes associated with all Neuron objects in this ChannelIndex
                if ~isempty( prevNeuron )
                    prevSpikes = findobj( epochs(ep).getChild( 'Spikes' ),'chanInd',chNum );
                    for neuron = 1:numel( prevNeuron )
                        prevSpikeInds = find( ismember( prevSpikes,prevNeuron(neuron).getChild( 'Spikes' ) ) );                        
                        epochs(ep).removeChild( 'Spikes',prevSpikeInds ); % remove previous spikes
                    end
                end
                
                % add the new spikes
                epochs(ep).addChild( Sp );
                newNeuron.addChild( Sp );
            end
                   
            fprintf( '\n' );
            % ==================================================================
            
            % updates
            clear N
            self.addChild( newNeuron );
            self.getParent( 'Block' ).update;
            
            % helper function
            function p = check_inputs( inputs )
                names = {'thresh','artifact','masked_detection','whiten'};
                defaults = [4,1000,false,false];
                p = inputParser;
                             
                for j = 1:numel(names)
                    p.addParameter( names{j},defaults(j) );
                end
                p.parse(inputs{:});
                p = p.Results;
            end
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
        
        
        function sortSpikes_using_clusterModel( self,W,mapping,model )
            % sortSpikes_using_clusterModel( self,W,mapping,model )
            %
            % sorts the spike waveforms in the block file using the output from
            % "create_cluster_models_peranimal"
            %
            % Inputs:
            %   W - the projection matrix used to project the training set
            %       created by "compute_spike_features"
            %
            %   mapping - the mapping structure used for the projection
            %             created by "compute_spike_features"
            %
            %   model - the HDBSCAN model used to sort the training set
            %
            % Outputs:
            %   creates new neurons internally such that self.nUnits
            %   reflects the number of sorted neurons as children to
            %   this channelindex object
            
            % check size consistency
            assert( all( size( self.getChild( 'Neuron',1 ).meanWaveform ) == [size( W,1 ),size( W,3 )] ),...
                    'Size of spikes not compatible with size of projection matrix W' );
            
            % re-merge neurons in this ChannelIndex object if necessary 
            if self.nUnits > 1
                self.mergeNeurons( [self.getChild('Neuron').ID] );
            end
            
            if mapping.useMap
                chanmap = self.chanMap;
            else
                chanmap = 1:self.nElectrodes;
            end
            
            % loop over epochs to reduce memory requirement
            nEpochs = numel( self.getSibling( 'Epoch','Block' ) );
            neuron = self.getChild( 'Neuron' );
            features = zeros( neuron.nSpikes,size( model.data,2 ) );
            labels = zeros( nSpikes,1 );
            P = zeros( nSpikes,1 );
            counter = 0;
            
            for j = 1:nEpochs
                [spikes,~,~,mask] = neuron.getSpikes(j);
                inds = counter+1:counter+size( spikes,2 );
                features(inds,:) = map_new_spikes( permute( spikes(:,:,chanmap),[2,1,3] ),W,mapping,full( mask(chanmap,:) ) );
                [labels(inds),P(inds)] = model.predict( features(inds,:) );
                counter = counter + size( spikes,2 );
                clear spikes mask
            end
                      
            % split neuron into isolated units
            [spikes,sptimes,epochs,mask] = neuron.getSpikes(); % all spikes
            sptimes = spikemat2vec( sptimes );
            clear neuron
            
            self.create_new_neurons( spikes,sptimes,epochs,ID,labels,...
                features,P,mapping.method,[],[],mask );
            
            clear spikes mask sptimes epochs
        end

        
        function [ID,A] = sortSpikes_using_templates( self,templates,varargin )
            % [ID,A] = sortSpikes_using_templates( self,templates,(useChanMap,useMask,resolveOverlap,thresh) )
            %
            % sorts the spike waveforms contained in this ChannelIndex
            % object by using template matching. Currently, template
            % matching is performed by minimizing the MSE between each
            % spike and template.
            %
            % Inputs:
            %   templates - the template structure created from
            %               "create_spike_templates"
            %
            %   (useMask) - a flag to mask channels (true) or not
            %               (default = false)
            %
            %   (useChanMap) - a flag to re-map the channels 
            %                  (default = false)
            %
            %   (resolveOverlap) - boolean flag for resolving overlapping spikes or not
            %                      (default = false)
            %   
            %   (thresh) - threshold for assigning spikes as: thresh*clust_amplitude_sd
            %              (default = 3)
            %
            % Outputs:
            %   ID - the m x 3 matrix of parent template IDs for each spike
            %
            %   A - the m x 3 matrix of amplitude coefficients for W (see below)
            %
            %   Template matching is done by greedily solving the equation:
            %           X_i = SUM{ A_ij*W_j + eta }
            %
            %   where A are coefficients describing the weighting of the
            %   jth mean waveform template W_j that, when summed across all,
            %   possible template, gives rise to the voltage wavefrom X_i
            %
            % This function also internally splits the neuron so that
            % self.nUnits reflects the newly created neurons. The outputs
            % are intended for verification of results / storage of
            % coefficients for bookeeping
                     
            % re-merge neurons in this ChannelIndex object if necessary 
            if self.nUnits > 1
                self.mergeNeurons([self.getChild('Neuron').ID]);
            end
            
            % check size consistency
            neuron = self.getChild( 'Neuron',1 );
            nEpoch = numel( neuron.getChild( 'Spikes' ) );
            [n,c] = size( neuron.meanWaveform );
            assert( all( n*c == size( templates.W,1 ) ),...
                    'Size of templates must match size of spikes' );
            
            % optional inputs    
            useChanMap = nargin < 2 || isempty( varargin{1} ) || varargin{1};
            useMask = nargin < 3 || isempty( varargin{2} ) || varargin{2};  
            resolveOverlap = nargin < 4 || isempty( varargin{3} ) || varargin{3};
            if nargin > 4; thresh = varargin{4}; else; thresh = 3; end
            
            chanmap = useChanMap*self.chanMap + ~useChanMap*(1:c)';
            [~,reverse_chanmap] = sort( chanmap );
            
            % create input variables from templates
            W = templates.W;
            M = templates.M;
            p = exp( -(templates.A_norm.sd ./ templates.A_norm.mu) );
            minAmp = templates.A_norm.mu - thresh*templates.A_norm.sd.*p;
            maxAmp = templates.A_norm.mu + thresh*templates.A_norm.sd.*p;
           
            % remove ID = 0 template
            muID = templates.labelMap == 0;
            W(:,muID,:) = [];
            M(:,muID) = [];
            minAmp(muID,:) = [];
            maxAmp(muID,:) = [];
            
            g = gpuDevice;      
            if ~isempty( g )
                W = gpuArray( W );
                minAmp = gpuArray( minAmp );
                maxAmp = gpuArray( maxAmp );
            end
            
            % loop over spike object to avoid huge matrix multiplications
            counter = 0;
            ID = zeros( neuron.nSpikes,3,'uint8' );
            A = zeros( neuron.nSpikes,3,'single' );
            
            for ep = 1:nEpoch
                if useMask
                    [spikes,~,~,mask] = neuron.getSpikes( ep );
                    spikes = maskchans( spikes,mask ); 
                    clear mask
                else
                    spikes = neuron.getSpikes( ep );
                end
                
                spikes = concatenateSpikes( spikes(:,:,chanmap) ); 
                
                if ~isempty( g )
                    [id,a] = template_match( gpuArray( spikes ),W,minAmp,maxAmp );
                    id = gather( id );
                    a = gather( a );
                else
                    [id,a] = template_match( spikes,W,bounds );
                end
                
                % store into our pre-allocated vecs
                m = size( id,1 );
                inds = counter+1:counter+m;
                ID(inds,:) = id;
                A(inds,:) = a;
                
                counter = counter + m;
                clear spikes
            end
            
            if ~isempty( g )
                W = gather( W );
            end

            [spikes,spikeTimes,trials,mask] = neuron.getSpikes();
            spikeTimes = spikemat2vec( spikeTimes );
            labels = ID(:,1)';
            
            if resolveOverlap
                [newSpikes,newID,parentSpikes] = resolve_overlapping_spikes( concatenateSpikes( spikes(:,:,chanmap) ),...
                                                                             spikeTimes,trials,ID,A,W(:,:,1) ./ templates.w_norm(~muID)' );
            
                % update spike properties in place with the newly discovered
                newSpikes = expandSpikes( newSpikes,c );
                newSpikes = newSpikes(:,:,reverse_chanmap);
                [parents,inds] = unique( parentSpikes );   
                spikes(:,parents,:) = newSpikes(:,inds,:);
                if ~isempty( mask )
                    mask(:,parents) = M(reverse_chanmap,newID(inds));
                end
                
                remainingInds = ~ismember( 1:numel( parentSpikes ),inds );
                parents = parentSpikes(remainingInds);
                spikes = [spikes,newSpikes(:,remainingInds,:)];
                spikeTimes = [spikeTimes,spikeTimes(parents)];
                trials = [trials,trials(parents)];
                labels = [labels,newID(remainingInds)];
                if ~isempty( mask )
                    mask = [mask,M(reverse_chanmap,newID(remainingInds))];
                end
                
                clear newSpikes newID parentSpikes
            end
            
            % split the unsorted neuron into sorted ones
            clear neuron W V bounds
            map = templates.labelMap(templates.labelMap>0);
            labels(labels>0) = map(labels(labels>0));
            labels = remove_small_clusters( labels,5 );
            
            self.create_new_neurons( spikes,spikeTimes,trials,0,labels,...
                [],[],[],[],[],mask );
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
                    features = [features;oldNeuron.features];
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
                            epoch.addChild( newSpikes(counter) );
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
            
            % don't split if it's the same labels as before
            if all( ismember( labels,oldID ) )
                return
            end
            
            % get the spike snips / times of this neuron
            [snips,sptimes,epochs,mask] = oldNeuron.getSpikes();
            sptimes = spikemat2vec( sptimes );
            
            % create the new neuron children given the labels / trials
            self.create_new_neurons( snips,sptimes,epochs,oldID,labels,...
                oldNeuron.features,oldNeuron.probabilities,[],[],[],mask );
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
        
        
        function [newLabels,newFeatures,newModel] = inspect_sorting_results( self,varargin )
            % [newLabels,newFeatures,newModel] = inspect_sorting_results( self,(neuronIDs) )
            %
            % opens an instance of the sortTool GUI to visualize the 
            % isolated spikes and inspect sorting quality. Can optionally 
            % supply "neuronIDs" argument to specify which neurons to
            % inspect (default = all)
            
            % check inputs
            neurons = self.getChild( 'Neuron' );
            if numel( neurons ) == 0 && neurons.ID == 0
                fprintf( 'Please sort neurons first\n' );
                return
            end
            
            IDs = [neurons.ID];
            if nargin > 1 && ~isempty( varargin{1} )
                requestedIDs = varargin{1};
                if ~all( ismember( requestedIDs,IDs ) )
                    fprintf( 'Some or all requested IDs do not exist\n' );
                    return
                end
                neurons = neurons(ismember( IDs,requestedIDs )); % only keep those we asked for
            end
            
            % extract spikes / features for the neurons
            nNeurons = numel( neurons );
            if nNeurons > 1
                nSpikes = sum( [neurons.nSpikes] );
                nNeurons = numel( neurons );
                [nPts,nChans] = size( neurons(1).meanWaveform );
                spikes = zeros( nPts,nSpikes,nChans,'single' );
                times = zeros( nSpikes,1,'single' );
                epoch = zeros( nSpikes,1,'single' );
                labels = zeros( nSpikes,1,'uint8' );
                mask = zeros( nChans,nSpikes,'single' );
                useFeatures = false;
                if ~isempty( neurons(1).features ) && all( ~any( isnan( neurons(1).features ) ) )
                    useFeatures = true;
                    features = zeros( nSpikes,size( neurons(1).features,2 ),'single' );
                else
                    features = [];
                end

                startInd = 0;
                for i = 1:nNeurons
                    inds = startInd + 1:startInd + neurons(i).nSpikes;
                    [s,t,ep,m] = neurons(i).getSpikes();
                    spikes(:,inds,:) = single( s );
                    epoch(inds) = uint8( ep );
                    mask(:,inds) = single( full( m ) );
                    t = reshape( t,numel( t ),1 );
                    times(inds) = single( t(~isnan(t)) );
                    labels(inds) = uint8( neurons(i).ID );

                    if useFeatures
                        features(inds,:) = single( neurons(i).features );
                    end
                    clear s t ep m
                    startInd = inds(end);
                end
            else
                [spikes,times,epoch,mask] = neurons.getSpikes();
                spikes = single( spikes );
                times = single( reshape( times,1,numel( times ) ) );
                times(isnan(times)) = [];
                epoch = uint8( epoch );
                mask = single( full( mask ) );
                labels = repmat( uint8( neurons.ID ),1,neurons.nSpikes );
                features = single( neurons.features );
            end
            
            % remap
            if ~isempty( self.chanMap )
                spikes = spikes(:,:,self.chanMap);
                mask = mask(self.chanMap,:);
            end
            
            % create an instance of the GUI
            [newLabels,newFeatures,newModel] = sortTool('data',spikes,...
                                                        'times',times,...
                                                        'trials',epoch,...
                                                        'mask',mask,...
                                                        'projection',features,...
                                                        'labels',labels);
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
            nonSortedIDs(nonSortedIDs==0) = []; % remove the ID = 0 so that we can  simply add to it
            
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
            %oldNeurons = self.getChild( 'Neuron',find( previousNeuronInds ) );
            %for i = 1:numel( oldNeurons )
            %    oldNeurons(i).deleteSelf();
            %end
            clear oldNeurons
            self.removeChild( 'Neuron',find( previousNeuronInds ) );

            % (f) now loop over the unique IDs and create a new Neuron object
            % and Spikes objects across Epochs
            uID = unique( newID );
            counter = 0;
            remainingNeurons = self.getChild( 'Neuron' );
            if ~isempty( remainingNeurons )
                currentIDs = [remainingNeurons.ID];
            else
                currentIDs = [];
            end
            for i = uID
                counter = counter+1;
                thisID = (newID==i);
                thisNeuron = ismember( currentIDs,i );
                clear Sp
                
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
                        prevSpikes(sp).nSpikes = prevSpikes(sp).nSpikes + sum( IDX );
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
            
            
