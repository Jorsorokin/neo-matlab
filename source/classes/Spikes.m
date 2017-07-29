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
            %   Signal
            %   Epoch
            %
            % Properties:
            %   epoch - the parent Epoch number
            %   unitID - the parent Neuron ID
            %   timeUnits - the units of the spike times (default = 's' for seconds)
            %   voltUnits - the units of the spike snips (default = 'uV')
            %   fs - the sampling rate of the spike snips
            %   times - the actual spike times
            %   voltage - an n x m x c matrix of spike snips (n = points, m = # of spikes, c = channels)
            %   nSpikes - the total number of spikes
            %   nChan - the number of channels associated with these spikes
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
                case {'Epoch','Neuron'}
                    addParent@Container( self,parent );
                    if isa( class( parent ),'Neuron' )
                        self.unitID = parent.ID; % add the associated unit ID
                        parent.nSpikes = parent.nSpikes + self.nSpikes;
                        self.chanInd = parent.chanInd;
                    else
                        self.epoch = parent.epochNum; % add the epoch 
                    end
                otherwise
                    error( 'Only Epoch or Neuron objects are valid parents' );
            end
        end
        
        
        function plot( self,varargin )
            % plot( self,(col,chans) )
            %
            % plot the spike waveforms. Can optionally specify an input
            % color as a second argument, and channels as third argument
            if nargin < 2 || isempty( varargin{1} )
                col = 'w';
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
            
            for ch = 1:nchan
                subplot( nchan,1,ch ); hold on;
                %fillPlot( self.voltage(:,:,chans(ch))',time,'sd',[],[],col );
                plot( time,self.voltage(:,:,chans(ch)),'color',col );
                set( gca,'tickdir','out','box','off','xlim',[time(1) time(end)],'color','k' );
                xlabel( 'time (ms)' );
                ylabel( self.voltUnits );
                title( sprintf( 'CH %i',chans(ch) ) );
            end
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
            x = 1:n; % original time-vec
            xx = linspace(1,n,n*4); % the interpolant x-values
            
            % loop over channels
            for c = 1:nchan
                spikes(:,:,c) = csaps( x,spikes(:,:,c)',amount,xx )'; % smooths the interpolant
            end
            
            % now downsample
            % downsample the spikes
            snips = zeros( size(spikes,1)/4,size(spikes,2),nchan );
            n = size( snips,1 );
            xx = 1:n; % the interpolant x-values
            x = linspace(1,n,n*4); % the original time-vec
            for c = 1:nchan
                snips(:,:,c) = spline( x,spikes(:,:,c)',xx )'; % no smoothing
            end
            
            % store back into self
            self.voltage = snips;
        end
                   
    end % methods
    
end
        