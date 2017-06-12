classdef Epoch < Container
    
    properties
        startTime
        stopTime
        duration
        eventTime
        epochNum
        nSignals = 0;
        nSpikes = 0;
    end
    
    methods
        
        function self = Epoch( startTime,stopTime,epochNum )
            % self = Epoch( startTime,stopTime,epochNum )
            %
            % Initiate an instance of the Epoch class. 
            % An Epoch contains the starting and ending time
            % of a specific segment of data recording during a 
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
            % Methods:
            %   car (common average reference)
            %   plotSignals
            %   
            %       * see also methods in the Container class

            
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
            % 
            % loop over the "Signal" objects specified by the indices "idx"
            % and for each, subtract the mean of all others. This is a
            % virtual reference method to eliminate common noise among
            % channels for the given Epoch. By default, "idx" will include
            % all "Signal" objects that are children of the current Epoch.
            %
            % Subtracts the common average reference (CAR) from voltages in
            % each "Signal" object, in place.
            
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
            
            % loop over the signal objects defined by IDX
            for i = IDX
                indices(i) = false;
                
                % get mean voltage waveforms 
                CAR = mean( [signals(indices).voltage],2 );
                
                % subtract the CAR from the ith signal voltage
                signals(i).voltage = bsxfun( @minus,signals(i).voltage,CAR );
                
                indices(i) = true;
            end 
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
            
            % loop over each channel group and individual signal and plot
            for chidx = 1:nChIdx
                data = sig(chidx).voltage;
                time = linspace( self.startTime,self.stopTime,sig(chidx).nPoints );
                vRange = 10 * mean( sig(chidx).estimateNoise() );
                for i = 1:sig(chidx).nSignals
                    plot( time,data(:,i)-(vRange * counter),'color',cmap(chidx,:) );
                    counter = counter + 1;
                end
            end 
            
            % clean up graph
            set( gca,'tickdir','out','box','off' );
            axis tight
            xlabel( 'time (s)' );
            title( sprintf( 'All signals for epoch %i',self.epochNum ) );
            
            % now plot a vertical line indicating an event if one exists
            if ~isempty( self.eventTime )
                yL = get(gca,'ylim');
                plot( [self.eventTime, self.eventTime],yL,'k--','linewidth',2 );
            end  
        end
            
    end
    
end
        