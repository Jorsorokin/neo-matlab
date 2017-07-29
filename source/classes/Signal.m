classdef Signal < Container
    
    properties
        units = 'uV';
        nPoints
        duration
        voltage
        nSignals
        nChan
        nSpikes
        fs
        epoch = NaN;
        chanInd = NaN;
    end
    
    methods
        
        function self = Signal( voltage,fs )
            % self = Signal( voltage,fs )
            %
            % Create an instance of the Signal class.
            % A Signal object contains digitally-sampled
            % voltages from a continuous recording. The voltages
            % should be an n x m matrix, where n = number of samples,
            % and m = number of channels. Thus, a ChannelIndex object
            % may reference a Signal object, where each column in 
            % Signal.voltage reflects one channel referenced by the 
            % ChannelIndex parent. Multiple Signal objects with the same
            % ChannelIndex parent may exist, which could reflect multiple
            % Epochs of data from the same group of channels.
            %
            % From the Signal object, one can extract all of parent Epochs
            % or ChannelIndex, and can also perform filtering from within
            % the Signal object itself.
            %
            % Children:
            %   Spikes
            %
            % Parents:
            %   Epoch
            %   ChannelIndex
            %
            % Properties:
            %   units - the measurement units of this signal (default = 'uV')
            %   nPoints - the number of sample points of this Signal
            %   duration - total duration in time
            %   voltage - the actual voltage values of this Signal
            %   nSignals - the number of channels contained in this Signal
            %   nChan - the (equal) number of channels of the parent ChannelIndex
            %   nSpikes - the number of spikes 
            %   fs - the sampling rate for the voltage traces
            %   epoch - the parent Epoch number 
            %   chanInd - the parent ChannelIndex number
            %
            % Methods:
            %   plot
            %   getSpikes
            %   filter
            %   resample
            %   estimateNoise
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
                case 'Spikes'
                    addChild@Container( self,child );
                    self.nSpikes = sum( [child.nSpikes] );
                otherwise
                    error( 'Only Spikes objects are valid children' );
            end
        end 


        function addParent( self,parent )
            switch class( parent )
                case {'Epoch','ChannelIndex'}
                    addParent@Container( self,parent );
                    parent.nSignals = parent.nSignals + self.nSignals;
                    if isa( class( parent ),'Epoch' )
                        self.epoch = parent.epochNum; % add the Epoch                       
                    else
                        self.chanInd = parent.chanIndNum; % add the ChannelIndex 
                    end
                otherwise
                    error( 'Only Epoch or ChannelIndex objects are valid parents' );
            end
        end
        
        
        function spikes = getSpikes( self )
            % spikes = getSpikes( self )
            %
            % pull out the Spikes objects (if any) associated with
            % this Signal object
            N = self.getParent( 'ChannelIndex' ).getChild( 'Neuron' );
            if ~isempty( N )
                for thisNeuron = 1:numel( N )
                    allSpikes = N(thisNeuron).getChild( 'Spikes' );
                    if ~isempty( allSpikes )
                        spikes(thisNeuron) = allSpikes.findobj( 'epoch',self.epoch );
                    end
                end
            else
                spikes = [];
            end
        end
        
        
        function plot( self )
            % plot( self )
            %
            % plot signals in the current Signal object.
            %
            % get the epoch start/end time if this Signal object has an
            % Epoch parent. Plot relative to this start/stop time. Also
            % plot "." for spike times (if any), color coded by the Neuron ID
            ep = self.getParent('Epoch');
            if ~isempty( ep )
                time = linspace( ep.startTime,ep.stopTime,self.nPoints );
            else
                time = linspace( 0,self.nPoints/self.fs,self.nPoints );
            end
            
            % find the spikes children, if any, and pull out spike times
            spikes = self.getSpikes();
            if ~isempty( spikes )
                cmap = colormap( jet( numel( spikes ) ) );
            end

            % loop over the signals
            vStep = 10 * mean( self.estimateNoise() );
            for i = 1:self.nSignals
                volt = self.voltage(:,i) - (vStep * (i-1));
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
            title( sprintf( 'ChanIndex: %i, Epoch: %i',...
                self.chanInd,self.epoch ) );
            
            % now plot a vertical line indicating an event if one exists
            if ~isempty( ep )
                if ~isempty( ep.eventTime )
                    yL = get(gca,'ylim');
                    plot( [ep.eventTime, ep.eventTime],yL,'w--','linewidth',2 );
                end
            end                
            hold off;
        end
        
        
        function filter( self,hp,lp )
            % filter( self,hp,lp )
            %
            % filter the signals between "hp" and "lp" (given as Hz)
            self.voltage = filtfilt2( self.voltage,hp,lp,self.fs );
        end
        
        
        function resample( self,p,q )
            % resample( self,p,q );
            %
            % resample the data as the ratio p/q
            self.voltage = resample( self.voltage,p,q );
            self.fs = self.fs * (p/q);
        end
        
        
        function noiseLevel = estimateNoise( self )
            % noiseLevel = estimateNoise( self )
            %
            % estimates the noise of each signal in self.voltage as:
            %   median( abs( voltage ) / 0.6745 )
            noiseLevel = median( abs( self.voltage ) / 0.6745 );
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
        
    end % methods
    
end