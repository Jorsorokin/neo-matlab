classdef Epoch < Container
    % self = Epoch( startTime,stopTime,epochNum )
    %
    % Initiate an instance of the Epoch class. 
    % An Epoch contains the starting and ending time
    % of a specific segment of data during a recording 
    % session (although this definition can be extended 
    % as a specific pair of times to hold any segment of data,
    % not specifically from just one recording session).
    %
    % The main purpose of the Epoch is to organize data into 
    % time-delimited events, each of which can refer to a specific
    % set of channels. The "epochs" can be anything with a starting 
    % and ending time in your recording. This makes accessing 
    % specific data based on a subset of channels across epochs
    % extremely easy, and helps organize analysis centered on
    % trial-averaging such as peri-stimulus time histograms (PSTHs).
    %
    % You can add a specific "eventTime" after initiating the Epoch 
    % object, which refers to a single point in time. This could be 
    % useful, for instance, for extracting a window of time
    % (the epoch) surrounding some stimulus (the event)
    %
    % Children:
    %   Signals 
    %   Spikes 
    %          
    % Parents:
    %   Block 
    %
    % Properties:
    %   startTime - starting time of the epoch
    %   stopTime - ending time of the epoch
    %   duration - total duration
    %   eventTime - start time of an event within the epoch
    %   nSignals - number of channels contained in this epoch
    %   nSpikes - number of action potentials detected in this epoch across all channels
    %   name - user-defined name to give to this epoch
    %
    % Methods:
    %   addEvent
    %   car (common average reference)
    %   rmMovement
    %   rmRedundantSpikes
    %   plotSignals
    %   stateSpace
    %   
    %       * see also methods in the Container class   

    properties
        startTime
        stopTime
        duration
        eventTime
        epochNum
        nSignals = 0;
        nSpikes = 0;
        name
    end
    
    methods
        
        function self = Epoch( startTime,stopTime,epochNum )
            self.startTime = startTime;
            self.stopTime = stopTime;
            self.duration = stopTime - startTime;
            self.epochNum = epochNum;
        end
        
        
        function addChild( self,child )
            % overloaded the "addChild" method on the Container class to
            % update specific properties of the Epoch class
            if max( strcmp( class( child ),{'Signal','Spikes'} ) ) == 0
                error( 'Only Signal & Spikes objects are valid children' );
            end
            
            switch class( child )
                case 'Signal'
                    self.nSignals = self.nSignals + numel( child );
                case 'Spikes'
                    self.nSpikes = self.nSpikes + sum( [child.nSpikes] );
            end
            addChild@Container( self,child );
            for j = 1:numel( child )
                child(j).epoch = self.epochNum;   
            end
        end
        
        
        function addParent( self,parent )
            % overload the "addParent" method of the Container class for
            % Epoch class specifically
            switch class( parent )
                case {'Block'}
                    addParent@Container( self,parent );
                otherwise
                    error( 'Only Block object is a valid parent' );
            end
        end
        
        
        function addEvent( self,eventTime )
            % add the onset time (in seconds) of an event during the epoch.
            % Can add multiple events into an array
            self.eventTime(end+1) = eventTime;
        end
        
        
        function car( self,varargin )
            % car( self,(electrodes) )
            % Subtracts the common average reference (CAR) from voltages in
            % each "Signal" object, in place.
            %
            % loop over the "Signal" objects specified by the indices "electrodes"
            % and get the mean of the voltages from all others that do not belong
            % to the same ChannelInde group(s). Then, subtract each mean from each 
            % signal
            %
            % This is a virtual reference method to eliminate common noise among
            % channels for the given Epoch. By default, "electrodes" will include
            % all "Signal" objects that are children of the current Epoch.
            
            % check inputs
            if nargin > 1 && ~isempty( varargin{1} )
                electrodes = varargin{1};
            else
                electrodes = 1:self.nSignals;
            end
            
            % pull out the signals
            signals = self.getChild( 'Signal',electrodes );
            if ~isempty( signals )
                noElectrode = cellfun( @isempty,{signals.electrode} ); % removes signals not from Electrodes (i.e. stims,etc)
                signals(noElectrode) = [];
            end
            if numel( signals ) == 0
                disp( 'No signals from specified electrodes available' );
                return
            end

            % create an index reference for extracting the right signals
            nSigs = numel( signals );
            indices = true( 1,nSigs );
            
            % preallocate our reference matrix
            CAR = zeros( signals(1).nPoints,nSigs );

            % get all the ChannelIndex IDs associated with each signal
            allChInd = {signals.chanInd}; % each cell contains channelindex numbers for that signal
            
            % get voltage from all signals
            allVolt = [signals.voltage];
            
            % loop over the signal objects defined by IDX
            count = 0;
            for i = 1:nSigs
                count = count + 1;

                % find Signals that are associated with a common channelindex
                commonChanInd = cell2mat( cellfun( @(x)...
                    any( ismember( x,signals(i).chanInd ) ),allChInd,'un',0 ) );
                indices(commonChanInd) = false;
                
                % get mean voltage waveforms across all other signals
                CAR(:,count) = mean( allVolt(:,indices),2 );
                
                indices(commonChanInd) = true;
            end 
            
            % subtract the CAR from the signal voltages
            count = 0;
            for i = 1:nSigs
                count = count+1;
                signals(i).voltage = signals(i).voltage - CAR(:,count);
            end
        end


        function rmMovement( self,varargin )
            % rmMovement( self, (corrThresh) )
            %
            % eliminates the spiketimes and waveforms contained within
            % the "spikes" children of this Epoch that are too highly 
            % correlated with others over different channels. 
            % 
            % Movement artifacts are usually present across all channels,
            % and this is an attempt to remove those artifacts 
            % from spike detection. You can set the max correlation 
            % coefficient allowed before elimination (default = 0.8)
            %
            % Caveats: with highly-sampled recordings (i.e. dense 
            % microelectrodes or tetrodes), true spikes may appear on multiple
            % channels. 

            % check inputs
            if nargin < 2 || isempty( varargin{1} )
                corrThresh = 0.8;
            else 
                corrThresh = varargin{1};
            end

            % get the Spikes and Signal objects of this epoch
            newSpikeCount = 0;
            spikes = self.getChild( 'Spikes' );
            signals = self.getChild( 'Signal' );

            % check for any spikes
            if isempty( spikes )
                disp( 'Must detect spikes first. Ending function.' );
            end

            % create a big matrix of the raw signals
            rawdata = [signals.voltage]; 
            FS = signals(1).fs;

            % create our pre/post # of samples to extract for each spike
            pretime = floor( 0.0005 * FS );
            posttime = floor( 0.0015 * FS );

            % loop over Spike objects
            for i = 1:numel( spikes )

                % find signals associated with this spike object
                childSigs = spikes(i).getSibling( 'Signal','Epoch' );
                siblingSigs = ismember( signals,childSigs );

                % pull out the spike times
                sptime = round( spikes(i).times * FS );
                badspikes = false( numel( sptime ),1 );

                % loop over spike times, take average spike snip across sibling channels
                % and compare to the other signals 
                for sp = 1:numel( sptime )
                    thisSnip = mean( rawdata( sptime(sp)-pretime:sptime(sp)+posttime,siblingSigs ),2 );
                    otherSnips = rawdata( sptime(sp)-pretime:sptime(sp)+posttime,~siblingSigs );

                    % get the correlation coefficient across channels for this spike snip
                    [coef,~] = corrcoef( [thisSnip,otherSnips] );

                    % check if any coefficient is larger than our max allowed. If so, indicate as a bad 
                    % spike and discard from this spike object
                    if max( coef(2:end,1) ) > corrThresh
                        badspikes(sp) = true;
                    end
                end

                % remove bad spikes
                spikes(i).times(badspikes) = [];
                spikes(i).voltage(:,badspikes,:) = [];
                spikes(i).nSpikes = sum( ~badspikes );
                newSpikeCount = newSpikeCount + spikes(i).nSpikes;

                % update the number of spikes in the Neuron parent
                neuron = spikes(i).getParent( 'Neuron' );
                if ~isempty( neuron )
                    neuron.nSpikes = neuron.nSpikes - sum( badspikes );
                end
            end

            % update the number of spikes in this epoch 
            self.nSpikes = newSpikeCount;
        end
        

        function rmRedundantSpikes( self )
            % rmRedundantSpikes( self )
            %
            % removes the redundant spike times and waveforms that were 
            % detected multiple times from different ChannelIndex objects
            %
            % This is an issue with dense, multi-electrode arrays where 
            % every K-neighbors to any ith electrode are considered a 
            % channel group. Thus, any electrode may be associated wtih multiple
            % ChannelIndex objects, and thus the same spikes may appear in these
            % more than once. 
            %
            % Spikes are deemed redundant if: 
            %   (a) the times of the action potentials are within +/- 0.5ms  
            %       (i.e. within the refractory period)
            %   (b) they originate from ChannelIndex objects with overlapping electrodes
            %   (c) the shapes of the action potentials alligned to their 
            %       peaks are highly overlapping ( >= 0.9 correlation )
            %
            % If any ith spike is deemed redundant, then by default the spike is
            % is removed from Spike objects except for the one 
            % in which it has the largest absolute voltage

            % get the spikes
            if self.nSpikes == 0
                disp( 'No spikes available' );
                return
            end
            
            % get the associated channelindex & electrode numbers
            spikes = self.getChild( 'Spikes' );
            parentChanInd = [spikes.chanInd];
            channelindex = self.getSibling( 'ChannelIndex','Block' );
            channelindex = channelindex( ismember( [channelindex.chanIndNum],parentChanInd ));
            electrodeNum = {channelindex.electrodes};
            
            % create our anonymous functions for checking redundant spikes
            % ======================================================================
            % (a) checks if any spike time in each Spikes object close to the spike time "t"
            checkA = @(x,t) cellfun( @(c)(find( abs( c-t ) <= 5e-4 )),x,'un',0 ); 
            
            % (b) checks if any other Spikes objects associated with a common electrode 
            checkB = @(x,ch) cellfun( @(c)(find( ismember( c,ch ) )),x,'un',0 );

            % general check for emptyness
            emptyCheck = @(x) ~cellfun( @isempty,x );
            % ======================================================================            
            
            % get our presets for the looping
            waveforms = nan( size( spikes(1).voltage,1 ),numel( spikes ) ); % to avoid creating matrices on each loop
            indices = true( 1,numel( spikes ));
            allSpikeTimes = {spikes.times};
            totalRedundantSpikes = 0;

            % now loop over Spikes objects and the spike times in each
            for sp = 1:numel( spikes )
                indices(sp) = false;
                t = 1;
                while t <= self.nSpikes
                    oldSpikeCount = self.nSpikes;
                    thisSpike = spikes(sp).times(t);

                    % check for any common electrode parents AND nearby spike, 
                    % and continue if none found
                    nearbySpikes = checkA( allSpikeTimes,thisSpike ); 
                    sameElectrode = checkB( electrodeNum,electrodeNum{sp} );
                    redundant = emptyCheck( cellfun( @times,sameElectrode,nearbySpikes ) );
                    if ~any( redundant(indices) )
                        t = t+1; % move onto next spike time
                        continue
                    end

                    % get the spike-waveforms for the potential redundant spike
                    % from all of the associated Spikes objects
                    redundant = find( redundant );
                    for j = redundant
                        waveforms(:,j) = squeeze( spikes(j).voltage(:,nearbySpikes{j},sameElectrode{j}) );
                    end
                    
                    % compute the spike-waveform correlation coefficients and 
                    % check for strong correlations
                    p = corrcoef( waveforms );
                    if ~any( p(:,indices) >= 0.9 )
                        t = t+1;
                        continue
                    end

                    % find largest peak-trough (i.e range) and remove all others
                    [~,keep] = max( range( waveforms ) );
                    for j = redundant
                        if j == keep
                            continue
                        end
                        spikes(j).times(nearbySpikes{j}) = [];
                        spikes(j).voltage(:,nearbySpikes{j},:) = [];
                        spikes(j).nSpikes = spikes(j).nSpikes - 1;
                    end

                    % if we didn't affect current spike object, update "t"
                    if spikes(sp).nSpikes == oldSpikeCount
                        t = t+1; 
                    end

                    totalRedundantSpikes = totalRedundantSpikes + 1;
                    waveforms(:) = nan;
                end
            end

            % finally update the block so everything is synched
            fprintf( 'Eliminated %i redundant spikes across %i Spikes objects',...
                totalRedundantSpikes,numel( spikes ) );
            block.update();
        end 


        function plotSignals( self )
            % plotSignals( self )
            %
            % plot signals for one Epoch (i.e. all channels).
            % each channel does not necessarily have to have the same
            % sampling rate, but has the same start/end time according to
            % the epoch

            % get the signals from the epoch
            sig = self.getChild( 'Signal' );
            if isempty( sig )
                disp( 'no signals available' );
            end
            
            % set up plotting parameters
            nChIdx = numel( sig ); % number of Signal objects (grouped channels)
            figure; hold on;
            %cmap = colormap( parula( nChIdx ) );
            
            % get noise estimate average
            noise = 0;
            for chidx = 1:nChIdx
                noise = noise + mean( sig(chidx).estimateNoise() );
            end
            noise = noise / chidx * 10;

            % plot the signals
            multisignalplot( [sig.voltage],sig(1).fs,[],noise );
            
            % clean up graph
            axis tight
            xlabel( 'time (s)' );
            title( sprintf( 'All signals for epoch %i',self.epochNum ) );
            
            % now plot a vertical line indicating an event if one exists
            if ~isempty( self.eventTime )
                yL = get(gca,'ylim');
                plot( [self.eventTime-self.startTime, self.eventTime-self.startTime],yL,'w--','linewidth',2 );
            end  

            % dark theme
            darkPlot( gcf )
        end


        function [projections,rate,kernel] = stateSpace( self,varargin )
            % [projections,rate,kernel] = stateSpace( self,(kernel,start,stop,nDim) )
            %
            % compute the firing rate of the spikes associated with this Epoch
            % and (assuming each spike train is from an isolated neuron)
            % compute the projections of the firing rates onto their PCs.
            %
            % if "kernel" is left blank, an optimal, adaptive kernel will be calculated
            % for each spike train. Additionally, one can specify "start" and "stop" in 
            % seconds relative to the duration of the epoch event time (if it exists)
            % to only calculate the projections for a subset of time of the entire epoch

            % check inputs
            if nargin > 1 && ~isempty( varargin{1} )
                kernel = varargin{1};
                findKernel = false;
            else
                findKernel = true;
            end
            if nargin > 2 && ~isempty( varargin{2} )
                start = varargin{2};
            else
                start = 0;
            end
            if nargin > 3 && ~isempty( varargin{3} )
                stop = varargin{3};
            else
                stop = self.duration;
            end
            if nargin > 4 && ~isempty( varargin{4} )
                nDim = varargin{4};
            else
                nDim = 3;
            end

            % estimate the kernel
            if findKernel
                % TO DO
            end

            % normalize the kernel
            kernel = kernel / norm( kernel );

            % get the spikes associated with this epoch. 
            spikes = self.getChild( 'Spikes' );
            nSp = numel( spikes );
            if isempty( spikes )
                disp( 'Must detect spikes first' );
                projections = nan;
                return
            end

            % Preallocate our firing rate matrix
            startPt = max( 1,ceil( start * spikes(1).fs ) );
            stopPt = ceil( stop * spikes(1).fs );
            rate = zeros( stopPt-startPt,nSp );

            for sp = 1:nSp
                % pull out the neuron
                neuron = spikes(sp).getParent( 'Neuron' );

                % get the firing rate
                fr = neuron.firingRate( kernel );

                % only keep the rate associated with this epoch
                thisEpoch = ismember( [neuron.getChild( 'Spikes' ).epoch],self.epochNum );
                rate(:,sp) = smooth( fr(startPt:stopPt-1,thisEpoch),floor( 0.02 * spikes(1).fs ) );
            end

            % TO DO:
            %   projections using manifold learning vs. PCA,
            %   or use a robust-version of PCA 

            % project the spike rates onto their PCs 
            [u,s] = svd( rate' );
            PCs = u(:,1:nDim) * s(1:nDim,1:nDim);
            projections = rate * PCs;
        end
        
    end % methods
    
end
        