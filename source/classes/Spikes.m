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
            %   denoise
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
            
            cols = ceil( sqrt( nchan ) );
            rows = ceil( nchan / cols );
            time = (0:size( self.voltage,1 )-1) / self.fs * 1000;
            ylimits = [min( self.voltage(:) ),max( self.voltage(:) )];
            
            counter = 0;
            rowCounter = 0;
            for ch = chans
                counter = counter+1;
                subplot( rows,cols,counter ); hold on;
                if ~isempty( self.mask )
                    bestSpikes = find( self.mask(ch,:) == 1 );
                else
                    bestSpikes = 1:self.nSpikes;
                end
                if ~isempty( bestSpikes )
                    %fillPlot( self.voltage(:,bestSpikes,ch)',time,'sd',[],[],col );
                    plot( time,self.voltage(:,bestSpikes,ch),'color',col );
                end

                % change plot appearance
                set( gca,'xlim',[time(1) time(end)],'ylim',ylimits,'color','k','ycolor','k','xcolor','k' );
                title( sprintf( 'CH %i',ch ),'color','w' );
                if mod( counter,cols ) == 1
                    ylabel( self.voltUnits );
                    set( gca,'ycolor','w' );
                    rowCounter = rowCounter + 1;
                end
                if rowCounter == rows
                    xlabel( 'time (ms)' );
                    set( gca,'xcolor','w' );
                end
            end

            % convert to dark theme and link plots
            set( gcf,'color','k' );
            linkaxes( get( gcf,'Children'),'xy' );
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
        
        function svdDenoise( self,varargin )
            % svdDenoise( self,(lastEig,varExp) )
            %
            % uses SVD to denoise the spike waveforms for each channel
            % by setting the singular values > lastEig equal to 0. Default
            % value for "lastEig" is that which explains >= "varExp" of the
            % variance. If "mask" is available in the current Spike object,
            % only spikes with mask > 0 for each channel will be used in
            % the SVD, to avoid excessive contamination of noise in the top
            % singular vectors/values of the SVD decomposition.
            %
            % The other optional argument ('varExp') allows the user to
            % specify how much variance to keep, rather than specifying the
            % exact number of eigenvectors. Default value is 75%
            
            % check if "lastEig" provided
            if nargin > 2 && ~isempty( varargin{1} )
                lastEig = varargin{1};
                if lastEig > size( self.voltage,2 )
                    disp( 'requested # of eigenvectors is greater than those available' );
                    return
                end
            end
            
            % check for varExp
            if nargin > 3 && ~isempty( varargin{2} )
                varExp = varargin{2};
                if varExp > 1 % likely provided as a probability
                    disp( 'varExp must lie within [0,1]' );
                    return
                end
            else
                varExp = 0.75;
            end

            % loop over channels
            for c = 1:size( self.voltage,3 )
                
                % check for mask matrix
                if ~isempty( self.mask )
                    unmasked = self.mask(c,:) > 0;
                    if nnz( unmasked ) < 5
                        continue % too few points for accurate denoising
                    end
                else
                    unmasked = true( 1,self.nSpikes );
                end
                
                % perform svd 
                [u,s,v] = svd( self.voltage(:,unmasked,c),'econ' );
                s = diag( s );
                
                % find final singular vector to keep
                if ~exist( 'lastEig','var' ) 
                    last = find( cumsum( s ) / sum( s ) >= varExp,1 );
                    s(last+1:end) = 0;
                else
                    s(lastEig+1:end) = 0;
                end
                
                % reconstruct using a reduced singular vector subspace
                self.voltage(:,unmasked,c) = u*diag( s )*v';
                clear u s v
            end
        end

    end % methods
    
end
        