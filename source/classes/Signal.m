classdef Signal < Container
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
        %   plot
        %   filter
        %   resample
        %   estimateNoise
        %
        %       * see also methods in the Container class
        
    properties
        units = 'uV'
        voltage
        fs
        duration
        nPoints
        chanInd = [];
        electrode = [];
        epoch = NaN;
    end

    methods

        function self = Signal( voltage,fs )
            
            self.voltage = voltage;
            self.fs = fs;
            self.nPoints = numel( voltage );
            self.duration = self.nPoints / self.fs;
        end


        function addChild( ~,~ )
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
                        self.electrode = parent.electrodeNum; % add the Electrode
                        self.chanInd(end+1) = parent.chanInd;
                    end
                otherwise
                    error( 'Only Epoch or ChannelIndex objects are valid parents' );
            end
        end

        function plot( self )
            % plot( self )
            %
            % plots the voltage trace contained in this Signal, and any 
            % spikes associated with the signal as '.' color-coded 
            % by the Neuron

            % get the epoch parent, if any
            if ~isnan( self.epoch )
                ep = self.getParent( 'Epoch' );
                start = ep.startTime;
                stop = ep.stopTime;
                event = ep.eventTime;
            else
                start = 0;
                stop = self.nPoints / self.fs;
                event = [];
            end

            % plot the voltage against time & any event
            time = linspace( start,stop,self.nPoints );
            plot( time,self.voltage,'color',[.85 .85 .85] );
            hold on;
            if ~isempty( event )
                plot( [event,event],[min( self.voltage ),max( self.voltage )],'w--' );
            end

            % plot the spikes on top of the voltage, if any
            if ~isnan( self.epoch )
                [~,times,~,neuronID,chanind] = self.getParent( 'Electrode' ).getSpikes( self.epoch );
                uniqueID = unique( neuronID );
                uniqueChanInd = unique( chanind );
                cmap = colormap( jet(numel( uniqueID )) ); % color signifies neuron ID
                markers = {'.','x','o','^','s'}; % marker type signifies channelindex
                volt = self.voltage;
                for id = uniqueID
                    for ch = 1:numel( uniqueChanInd )
                        spikes = times( ismember( neuronID,id) & ismember( chanind,uniqueChanInd(ch) ) );
                        scatter( spikes + start,volt(round( spikes*self.fs )),...
                            80,markers{ch},'color',cmap(id,:) );
                    end
                end
            end

            % clean up
            darkPlot( gcf );
            ylabel( self.units );
            xlabel( 'time (s) ');
            title( sprintf( 'Signal %i from Electrode %i', self.epoch,self.electrode ) );
            axis tight
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
            self.nPoints = size( self.voltage,1 );
        end
        
        
        function noiseLevel = estimateNoise( self )
            % noiseLevel = estimateNoise( self )
            %
            % estimates the noise of each signal in self.voltage as:
            %   median( abs( voltage ) / 0.6745 )
            noiseLevel = median( abs( self.voltage ) ) / 0.6745;
        end

    end % methods
end
