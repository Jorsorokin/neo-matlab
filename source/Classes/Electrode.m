classdef Electrode < Container
    
    properties
        nSignals = 0;
        nChan
        nSpikes
        chanInd = [];
        electrodeNum
        name
    end
    
    methods
        
        function self = Electrode( voltage,fs,electrode )
            % self = Signal( voltage,fs,electrode )
            %
            % Creates an instance of the Electrode class.
            % An Electrode is a container for both Signal and Spike objects,
            % tying the two different types of data that originate from
            % the same recording channel. 
            %
            % An Electrode object can belong to multiple ChannelIndex
            % objects (as in the case of dense, multi-electrode 
            % probes/arrays), and can also have multiple Signal and 
            % Spike children, which refer to different Epochs. 
            %
            % Children:
            %   Spikes
            %   Signal
            %
            % Parents:
            %   ChannelIndex
            %
            % Properties:
            %   nSignals - the number of channels contained in this Signal
            %   nSpikes - the number of spikes 
            %   chanInd - the parent ChannelIndex number
            %   electrodeNum - a unique number refering to this electrode
            %   name - a user-defined tag for the electrode (i.e. "shank 3, chan 2")
            %   
            % Methods:
            %   plot
            %   getSpikes
            %   getSignals
            %   waveDenoise
            %   
            %       * see also methods in the Container class

            self.voltage = voltage;
            self.nPoints = size( voltage,1 );
            self.duration = self.nPoints / fs;
            self.nSignals = size( voltage,2 );
            self.fs = fs;
        end
        
        
        function addChild( self,child )
            switch class( child )
                case {'Spikes','Signal'}
                    addChild@Container( self,child );
                    if isa( child,'Spikes' )
                        self.nSpikes = sum( [child.nSpikes] );
                    else
                        self.nSignals = self.nSignals + numel( child );
                    end

                    for j = 1:numel( child )
                        child(j).electrode(end+1) = self.electrodeNum;
                    end
                otherwise
                    error( 'Only Spikes or Signal objects are valid children' );
            end
        end 


        function addParent( self,parent )
            switch class( parent )
                case 'ChannelIndex'
                    addParent@Container( self,parent );
                    parent.nSignals = parent.nSignals + self.nSignals;
                    parent.nChan = parent.nChan + 1;
                    self.chanInd = [self.chanInd,parent.chanIndNum]; % add the ChannelIndex 
                otherwise
                    error( 'Only ChannelIndex objects are valid parents' );
            end
        end
        
        
        function [snips,times,epochs] = getSpikes( self,varargin )
            % [snips,times,epochs] = getSpikes( self, (epoch) )
            %
            % pull out the Spikes objects (if any) associated with
            % this Electrode. Can optionally specify which epochs

            % check inputs
            if nargin > 1 && ~isempty( varargin{1} )
                epochs = varargin{1};
            else
                epochs = 1:numel( self(1).getChild( 'Spikes' ) );
            end
            
            % get the "Spikes" child of this Neuron
            child = self.getChild( 'Spikes',epochs );
            if isempty( child )
                disp( 'no spikes found' );
                snips = nan;
                times = nan;
                epochs = nan;
                return;
            end
            
            % preallocate vectors
            npoints = size( child(1).voltage,1 );
            nSp = [child.nSpikes];
            nEpoch = numel( child );
            thisChan = ismember( child(1).electrode,self.electrodeNum );
            times = nan( max( nSp ),nEpoch );
            snips = nan( npoints,self.nSpikes );
            epNum = uint16( zeros( 1,self.nSpikes ) );
            counter = 0;
            
            % loop over Spike objects
            for i = 1:numel( child )
                snips(:,counter+1:counter+nSp(i)) = squeeze( child(i).voltage(:,:,thisChan) );
                times(1:nSp(i),i) = child(i).times;
                epNum(counter+1:counter+nSp(i)) = child(i).epoch;
                counter = counter + nSp(i);
            end
        end


        function voltage = getSignals( self,varargin )
            % voltage = getSignals( self,(epochs) )
            % 
            % pull out the voltage waveforms from the Signal
            % children (if any) associated with this Electrode.
            % 
            % voltage is an n x m matrix, with n = max number of points 
            % across epochs, and m = # of signals (i.e. # of epochs).

            % check inputs
            if nargin > 1 && ~isempty( varargin{1} )
                epochs = varargin{1};
            else
                epochs = 1:self.nSignals;
            end

            % pull out the signal objects
            signals = self.getChild( 'Signal',epochs );
            if ~isempty( signals )

                % preallocate the voltage matrix
                n = max( [signals.nPoints] );
                voltage = nan( n,self.nSignals );

                % loop over signals
                for sig = 1:self.nSignals
                    voltage(1:signals(sig).nPoints) = signals.voltage;
                end
            else
                disp( 'No raw data traces available' );
            end
        
        
        function waveDenoise( self,wLevel,wType )
            % waveDenoise( self,wLevel,wType )
            %
            % denoises the signals using multi-signal wavelet denoising.
            % 
            % wLevel equals the wavelet decomposition level desired (lower
            % = less smoothing), and wType equals the wavelet to use.
            
            % get the multi-signal wavelet decomposition
            dec = mdwtdec( 'c',self.voltage,wLevel,wType );
            
            % denoise the decomposition using multi-resolution
            self.voltage = mswden( 'den',dec,'rigrsure','mln','s' );
        end


        function plot( self )
            % plot( self )
            %
            % plot signals in the current Signal object.
            %
            % get the epoch start/end time if this Signal object has an
            % Epoch parent. Plot relative to this start/stop time. Also
            % plot "." for spike times (if any), color coded by the Neuron ID
            ep = self.getPartner( 'Epoch','Signal' );
            if ~isempty( ep )
                time = linspace( ep.startTime,ep.stopTime,self.nPoints );
            else
                time = linspace( 0,self.nPoints/self.fs,self.nPoints );
            end
            
            % get the voltage traces
            voltage = self.getSignals();

            % find the spikes children, if any, and pull out spike times
            spikes = self.getSpikes();
            if ~isempty( spikes )
                cmap = colormap( jet( numel( spikes ) ) );
            end

            % loop over the signals
            vStep = 10 * self.getChild( 'Signal',1 ).estimateNoise();
            for i = 1:self.nSignals
                volt = voltage(:,i) - (vStep * (i-1));
                plot( time,volt,'color',[.85 .85 .85] ); hold on;

                % check if any "spikes" objects exist. If so, plot 
                % as dots, color-coded by the Neuron ID
                if ~isempty( spikes )
                    for sp = 1:numel( spikes )
                        plot( spikes(sp).times + ep.startTime,...
                            volt(round( spikes(sp).times * spikes(sp).fs )),...
                            '.','color',cmap(sp,:),'markersize',10 );
                    end
                end
            end
            
            % clean up graph
            set( gca,'tickdir','out','box','off','color','k' );
            axis tight
            xlabel( 'time (s)' );
            title( sprintf( 'Electrode: %i, ChanIndex: %i, all Epochs',...
                self.electrodeNum, self.chanInd ) );
            
            % now plot a vertical line indicating an event if one exists
            if ~isempty( ep )
                if ~isempty( ep.eventTime )
                    yL = get(gca,'ylim');
                    plot( [ep.eventTime, ep.eventTime],yL,'w--','linewidth',2 );
                end
            end                
            hold off;
        end
        
    end % methods
    
end