classdef Epoch < Container
    
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
            %   plotSignals
            %   stateSpace
            %   
            %       * see also methods in the Container class
            %
            % Examples:
            %
            %   % create two Epochs, each a specific stimulus and of
            %   % varying lengths.
            %   startTimes = [12, 44];
            %   endTimes = [20, 48];
            %   events = [16, 46];
            %   names = {'stim1','stim2'};
            %   for j = 1:2
            %       epoch(j) = Epoch( startTimes(j),endTimes(j),j );
            %       epoch(j).name = names{j};
            %       epoch(j).eventTime = events(j);
            %   end

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
            
            for j = 1:numel( child )
                addChild@Container( self,child(j) );
                child(j).epoch = self.epochNum;               

                switch class( child )
                    case {'Signal'}
                        self.nSignals = self.nSignals + child(j).nSignals;
                    case {'Spikes'}
                        self.nSpikes = self.nSpikes + child(j).nSpikes;
                end
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
            % car( self,(IDX) )
            % Subtracts the common average reference (CAR) from voltages in
            % each "Signal" object, in place.
            %
            % loop over the "Signal" objects specified by the indices "IDX"
            % and for each, subtract the mean of all others. This is a
            % virtual reference method to eliminate common noise among
            % channels for the given Epoch. By default, "IDX" will include
            % all "Signal" objects that are children of the current Epoch.

            % get the "Signal" children, if any
            signals = self.getChild( 'Signal' );
            if isempty( signals )
                error( 'No signals found' );
            end
            
            % check inputs
            if nargin > 1 && ~isempty( varargin{1} )
                IDX = varargin{1};
            else
                IDX = 1:numel( signals );
            end
            
            % create an index reference for extracting the right signals
            indices = false( 1,numel( signals ) );
            indices(IDX) = true;
            
            % preallocate our reference matrix
            CAR = zeros( signals(IDX(1)).nPoints,numel( IDX ) );
            
            % loop over the signal objects defined by IDX
            count = 0;
            for i = IDX
                count = count+1;
                indices(i) = false;
                
                % get mean voltage waveforms 
                CAR(:,count) = mean( [signals(indices).voltage],2 );
                
                indices(i) = true;
            end 
            
            % subtract the CAR from the signal voltages
            count = 0;
            for i = IDX
                count = count+1;
                signals(i).voltage = bsxfun( @minus,signals(i).voltage,CAR(:,count) );
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
            rawdata = zeros( size( signals(1).voltage,1 ),self.nSignals );
            counter = 1;
            for i = 1:numel( signals )
                rawdata(:,counter:counter+signals(i).nSignals-1) = signals(i).voltage;
                counter = counter + signals(i).nSignals;
            end
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
            counter = 0;
            figure; hold on;
            cmap = colormap( parula( nChIdx ) );
            
            % get noise estimate average
            noise = 0;
            for chidx = 1:nChIdx
                noise = noise + mean( sig(chidx).estimateNoise() );
            end
            noise = noise / chidx * 10;

            % loop over each channel group and individual signal and plot
            for chidx = 1:nChIdx
                data = sig(chidx).voltage;
                multisignalplot( data - (noise*counter),sig(chidx).fs,cmap(chidx,:),noise );
                counter = counter + sig(chidx).nSignals;
            end 
            
            % clean up graph
            set( gca,'tickdir','out','box','off','ytick',[],'yticklabel',[] );
            axis tight
            xlabel( 'time (s)' );
            title( sprintf( 'All signals for epoch %i',self.epochNum ) );
            
            % now plot a vertical line indicating an event if one exists
            if ~isempty( self.eventTime )
                yL = get(gca,'ylim');
                plot( [self.eventTime-self.startTime, self.eventTime-self.startTime],yL,'k--','linewidth',2 );
            end  
        end


        function [projections,rate,kernel] = stateSpace( self,varargin )
            % [projections,rate,kernel] = stateSpace( self, (kernel,start,stop) )
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


            % estimate the kernel
            % TO DO
            %   To be filled in ...

            % normalize the kernel
            kernel = kernel / norm( kernel );

            % get the spikes associated with this epoch. 
            spikes = self.getChild( 'Spikes' );
            nSp = numel( spikes );
            if isempty( spikes )
                disp( 'Must detect spikes first' );
                projections = nan;
                return;
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
            PCs = u(:,1:3) * s(1:3,1:3);
            projections = rate * PCs;
        end
        
    end % methods
    
end
        