classdef Signal < Container
    
    properties
        units = 'uV';
        nPoints
        duration
        voltage
        nChan
        nSignals
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
            %   none
            %
            % Parents:
            %   Epoch
            %   ChannelIndex
            %
            % Methods:
            %   plot
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
        
        
        function addChild( ~,~ )
            disp ( 'Only Neuron object is a valid child' );
        end 
        
        
        function plot( self )
            % plot( self )
            %
            % plot signals in the current Signal object.
            %
            % get the epoch start/end time if this Signal object has an
            % Epoch parent. Plot relative to this start/stop time if it
            % Epoch is a parent
            ep = self.getParent('Epoch');
            if ~isempty( ep )
                time = linspace( ep.startTime,ep.stopTime,self.nPoints );
            else
                time = linspace( 0,self.nPoints/self.fs,self.nPoints );
            end
            
            % loop over the signals
            vStep = 10 * mean( self.estimateNoise() );
            for i = 1:self.nSignals
                volt = self.voltage(:,i) - (vStep * (i-1));
                plot( time,volt ); hold on;
            end
            
            % clean up graph
            set( gca,'tickdir','out','box','off' );
            axis tight
            xlabel( 'time (s)' );
            title( sprintf( 'ChanIndex: %i, Epoch: %i',...
                self.chanInd,self.epoch ) );
            
            % now plot a vertical line indicating an event if one exists
            if ~isempty( ep )
                if ~isempty( ep.eventTime )
                    yL = get(gca,'ylim');
                    plot( [ep.eventTime, ep.eventTime],yL,'k--','linewidth',2 );
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