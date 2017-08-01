classdef Electrode < Container
    
    properties
        nSignals = 0;
        nChan
        nSpikes
        epoch = NaN;
        chanInd = NaN;
    end
    
    methods
        
        function self = Electrode( voltage,fs )
            % self = Signal( voltage,fs )
            %
            % Creates an instance of the Electrode class.
            % An Electrode is a container for both Signal and Spike objects,
            % tying the two different types of data that originate from
            % the same recording channel. 
            %
            % An Electrode object can belong to multiple ChannelIndex
            %  objects (as in the case of dense, multi-electrode 
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
            %   nChan - the (equal) number of channels of the parent ChannelIndex
            %   nSpikes - the number of spikes 
            %   chanInd - the parent ChannelIndex number
            %   
            % Methods:
            %   plot
            %   getSpikes
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
                otherwise
                    error( 'Only Spikes or Signal objects are valid children' );
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


        function voltage = getVoltage( self )
            % voltage = getVoltage( self )
            % 
            % pull out the voltage waveforms from the RawData 
            % children (if any) associated with this Signal object
            data = self.getChild( 'Signal' );
            if ~isempty( data )
                voltage = [data.voltage];
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
        
    end % methods
    
end