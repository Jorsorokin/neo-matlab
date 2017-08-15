classdef Spikes < Container
    
    properties
        epoch = NaN;
        unitID = NaN;
        chanInd = NaN;
        timeUnits = 's';
        voltUnits = 'uV';
        fs
        times
        voltage
        nSpikes;
        nChan
        mask = [];
    end
    
    methods
        
        function self = Spikes( times,snips,fs )
            % self = Spikes( times,snips,fs )
            %
            % Create an instance of the Spikes class.
            % A Spikes object contains the spike times (in seconds)
            % and voltage waveforms of an extracellularly-recorded
            % neuron. 
            %
            % The Spikes object references a parent Neuron and Epoch
            % instance, each of which can be used to extract the 
            % spike times/voltages. Moreover, each Neuron object can 
            % extract multiple Spikes objects as a way to facilitate
            % trial (Epoch)-averaged analysis. 
            %
            % Children:
            %   none
            %
            % Parents:
            %   Neuron
            %   Epoch
            %
            % Properties:
            %   epoch - the parent Epoch number
            %   unitID - the parent Neuron ID
            %   chanInd - the associated ChannelIndex
            %   timeUnits - the units of the spike times (default = 's' for seconds)
            %   voltUnits - the units of the spike snips (default = 'uV')
            %   fs - the sampling rate of the spike snips
            %   times - the actual spike times
            %   voltage - an n x m x c matrix of spike snips (n = points, m = # of spikes, c = channels)
            %   nSpikes - the total number of spikes
            %   nChan - number of channels associated with this spike (i.e. those in the associated ChannelIndex)
            %   mask - a sparse matrix with 0 <= (i,j) <= 1, indicating the ith channels that detected each jth spikes 
            %          (see function "double_flood_fill.m")
            %
            % Methods:
            %   plot
            %   smooth
            %
            %       * see also methods in the Container object

            self.times = times;
            self.voltage = snips;
            self.fs = fs;
            self.nSpikes = numel( times );
            self.nChan = size( snips,3 );
        end
        
        
        function addChild( ~,~ )
            disp( 'Spike objects have no children' );
        end
        
        
        function addParent( self,parent )
            switch class( parent )
                case {'Neuron','Epoch'}
                    addParent@Container( self,parent );
                    if isa( class( parent ),'Neuron' )
                        self.unitID = parent.ID; % add the associated unit ID
                        parent.nSpikes = parent.nSpikes + self.nSpikes;
                        self.chanInd = parent.chanInd;
                    elseif isa( class( parent ),'Epoch' )
                        self.epoch = parent.epochNum; % add the epoch 
                    end
                otherwise
                    error( 'Only Epoch and Neuron objects are valid parents' );
            end
        end
        
        
        function plot( self,varargin )
            % plot( self,(col,chans) )
            %
            % plot the spike waveforms. Can optionally specify an input
            % color as a second argument, and channels as third argument.
            %
            % If "self.mask" is populated, each channel will only display spikes
            % that are strongly associated with that channel (i.e. mask(ch,:)==1)
            if nargin < 2 || isempty( varargin{1} )
                col = [0.85 0.85 0.85];
            else
                col = varargin{1}; 
            end
            if nargin < 3
                nchan = self.nChan; 
                chans = 1:nchan;
            else
                chans = varargin{2}; 
                nchan = numel( chans ); 
            end
            
            time = (0:size( self.voltage,1 )-1) / self.fs * 1000;
            
            counter = 0;
            for ch = chans
                counter = counter+1;
                subplot( nchan,1,counter ); hold on;
                %fillPlot( self.voltage(:,:,chans(ch))',time,'sd',[],[],col );
                if ~isempty( self.mask )
                    bestSpikes = find( self.mask(ch,:) == 1 );
                else
                    bestSpikes = 1:self.nSpikes;
                end
                if ~isempty( bestSpikes )
                    plot( time,self.voltage(:,bestSpikes,ch),'color',col );
                    set( gca,'xlim',[time(1) time(end)] );
                    xlabel( 'time (ms)' );
                    ylabel( self.voltUnits );
                    title( sprintf( 'CH %i',ch ) );
                end
            end

            % convert to dark theme
            darkPlot( gcf );
        end
        
        
        function smooth( self,varargin )
            % smooth( self, (amount) )
            %
            % smooth the voltage waveforms using a smoothing spline
            % can optionally input a smoothing parameter that controls how
            % much smoothing is applied. Default = 0.5
            if nargin < 2
                amount = 0.5;
            else
                amount = varargin{1};
            end
            
            % get the interpolant values
            n = size( self.voltage,1 );
            nSp = self.nSpikes;
            nchan = self.nChan;
            spikes = zeros( n*4,nSp,nchan ); % up-samples for interpolation
            x = 1:n;
            
            % loop over channels
            for c = 1:nchan
                spikes(:,:,c) = csaps( x,spikes(:,:,c)',amount,x )'; % smoothing cubic spline
            end
            
            % store back into self
            self.voltage = snips;
        end
                   
    end % methods
    
end
        