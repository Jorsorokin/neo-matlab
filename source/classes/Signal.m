classdef Signal < Container
    
    properties
        units = 'uV'
        voltage
        duration
        nPoints
        fs
        electrode = [];
        epoch = NaN;
    end

    methods

        function self = Signal( voltage,fs )
            % self = Signal( voltage,fs )
            %
            % Create an instance of the Signal class.
            % A Signal object contains digitally-sampled
            % voltages from a continuous recording.
            %
            % Multiple Signal objects with the same
            % ChannelIndex parent may exist, which could reflect multiple
            % Epochs of data from the same group of channels.
            %
            % Children
            %   none
            %
            % Parents
            %   Electrode
            %   Epoch
            %
            % Properties
            %   units - the measurement units of this signal (default = 'uV')
            %   fs - the sampling rate for the voltage traces
            %   nPoints - the number of sample points of this Signal
            %   duration - total duration in time
            %   voltage - the actual voltage values of this Signal
            %   electrode - the parent electrode
            %   epoch - the parent epoch #
            %
            % methods
            %   filter
            %   resample
            %   estimateNoise
            %
            %       * see also methods in the Container class
            
            self.voltage = voltage;
            self.fs = fs;
            self.nPoints = size( voltage,1 );
            self.duration = self.nPoints / self.fs;
        end


        function addChild( self,child )
            disp( 'Signal objects have no children' );
        end 


        function addParent( self,parent )
            switch class( parent )
                case {'Electrode','Epoch'}
                    addParent@Container( self,parent );
                    parent.nSignals = parent.nSignals + 1;
                    if isa( class( parent ),'Epoch' )
                        self.epoch = parent.epochNum; % add the Epoch                       
                    else
                        self.electrode = parent.electrodeNum; % add the ChannelIndex 
                    end
                otherwise
                    error( 'Only Epoch or ChannelIndex objects are valid parents' );
            end
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

    end % methods
end
