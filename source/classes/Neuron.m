classdef Neuron < Container
    
    properties
        ID
        chanInd = NaN;
        nSpikes = 0;
        nChan = 0;
        features = NaN;
        probabilities = NaN;
        keptSpikes = NaN;
        featureMethod = NaN;
    end
    
    methods
        
        function self = Neuron( ID )
            % self = Neuron( ID )
            %
            % Create an instance of the Neuron class. 
            % A Neuron object refers to an identified and 
            % isolated Neuron after sorting voltage AP waveforms
            % detected from a Signal object. 
            % 
            % For each ChannelIndex, each Neuron has a unique
            % identifier within that ChannelIndex that separates it
            % and its spike times/waveforms from other Neurons. 
            % It can have multiple Spikes objects as children, which 
            % represents multiple Epochs in which the Neuron was detected.
            %
            % Children:
            %   Spikes
            %
            % Parents:
            %   ChannelIndex
            %
            % Methods:
            %   getSpikes
            %   plotSpikes
            %   raster
            %   psth
            %   firingRate
            %   plotFeatures
            %   
            %       * see also methods in the Container class

            self.ID = ID;
        end
        
        
        function addChild( self,child )
            switch class( child )
                case 'Spikes'
                    addChild@Container( self,child );
                    for j = 1:numel( child )
                        self.nSpikes = self.nSpikes + child(j).nSpikes;
                        child(j).unitID = self.ID;
                    end
                otherwise
                    error( 'Only Spikes objects are valid children' );
            end
        end 
        
        
        function addParent( self,parent )
            switch class( parent )
                case 'ChannelIndex'
                    addParent@Container( self,parent );
                    self.chanInd = parent.chanIndNum; 
                    self.nChan = parent.nChans;
                otherwise
                    error( 'Only ChannelIndex objects are valid parents' );
            end
        end
        
        
        function [snips,times,epNum] = getSpikes( self )
            % [snips,times,epNum] = getSpikes( self )
            %
            % extracts the spike voltage waveforms and spike times across 
            % epochs for this Neuron object. Returns snips, spike times,
            % and a vector "epNum" refering to which epoch (Spikes obj)
            % each voltage/spiketime came from
            
            % get the "Spikes" child of this Neuron
            child = self.getChild( 'Spikes' );
            if isempty( child )
                disp( 'no spikes found' );
                snips = [];
                return;
            end
            
            % preallocate vectors
            npoints = size( child(1).voltage,1 );
            nSp = [child.nSpikes];
            nEpoch = numel( child );
            times = nan( max( nSp ),nEpoch );
            snips = nan( npoints,self.nSpikes,self.nChan );
            epNum = uint16( zeros( 1,self.nSpikes ) );
            counter = 0;
            
            % loop over Spike objects
            for i = 1:numel( child )
                snips(:,counter+1:counter+nSp(i),:) = child(i).voltage;
                times(1:nSp(i),i) = child(i).times;
                epNum(counter+1:counter+nSp(i)) = child(i).epoch;
                counter = counter + nSp(i);
            end
        end
        
        
        function plotSpikes( self,varargin )
            % plotSpikes( self,(chans) )
            %
            % plot the spike waveforms, color-coded by the epoch. Optionally 
            % specify which channels to plot
            spChild = self.getChild( 'Spikes' );
            nChild = numel( spChild );
            if isempty( spChild )
                disp( 'no spikes found' );
                return;
            end
            
            if nargin < 2 || isempty( varargin{1} )
                chans = 1:spChild(1).nChan;
            else
                chans = varargin{1};
            end
            
            cmap = colormap( hsv(nChild) );
            
            % plot spike waveforms for each "Spikes" child
            for j = 1:nChild
                spChild(j).plot( cmap(j,:),chans );
            end
            suptitle( sprintf( 'Unit %i',self.ID ) );
        end
                
        
        function raster( self,start,pre,post )
            % raster( self,start,pre,post );
            %
            % make a raster plot of the spikes across epochs associated
            % with this Neuron. Must provide "start", "pre", and "post"
            % relative to the spiketimes contained within the
            % Neuron. Assumes the trials (epochs) are aligned for
            % meaningful results
            
            % pull out the spike times across epochs
            [~,sptm,~] = self.getSpikes(); 
            
            % plot the raster according to the pre/post time and starting
            % time provided
            PlotRasters( sptm,start,pre,post );
        end
        
        
        function [count,avgCount,sigma] = psth( self,start,pre,post,varargin )
            % [count,avgCount,sigma] = psth( self,start,pre,post,(bw) )
            %
            % compute the peri-stimulus time histogram (PSTH) from the
            % Neuron object across all epochs. Output is the bin-count
            % (count), trial-averaged count (avgCount), and variance across
            % trials (sigma)
            
            % check for bin width input
            if nargin > 4 && ~isempty( varargin )
                bw = varargin{1};
            else 
                bw = [];
            end
            
            % pull out spike times
            [~,sptm,~] = self.getSpikes();
            
            % compute the PSTH
            [count,avgCount,sigma] = CalculatePSTH( sptm,start,pre,post,bw );
        end


        function rate = firingRate( self,kernel )
            % rate = firingRate( self,kernel )
            %
            % compute the firing rate of the Spikes children of the 
            % current Neuron object for each epoch separately by 
            % convolving with "kernel" 
            % 
            % "rate" will be an n x m matrix, where n = total number of 
            % points of the maximum duration of all Epochs, and m = number 
            % of Epochs. 
            
            % check for spikes
            spikes = N.getChild( 'Spikes' );
            nSp = numel( spikes )
            if isempty( spikes )
                assert( 'No spikes found' );
                rate = [];
            end

            % normalize the kernel
            kernel = kernel / norm( kernel );

            % get the Epoch parents of the Spikes and Signal children of the ChannelIndex
            epochs = self.getPartner( 'Epoch','Spikes' );
            signals = self.getSibling( 'Signal','ChannelIndex' );
            duration = [E.duration]; 
            npoints = ceil( duration .* [signals.fs] ); 

            % pre-allocate our firing-rate matrix based on max number of points
            rate = zeros( max( npoints ),nSp ); 

            % loop over the Spikes children, extract the spike times
            for ep = 1:nSp
                sptm = round( spikes(ep).spikeTimes * signals(ep).fs );
                rate(sptm,ep) = 1;
            end

            % convolve with the kernel
            rate = conv( rate,kernel,'same' );
        end

               
        function plotFeatures( self,column1,column2 )
            % plotFeatures( self,column1,column2 )
            %
            % plots a 2D scatter of the two coloumns of "features"
            % specified by the inputs "columns1" and "column2".
            % Must have already sorted spikes for this function to work
            
            if isnan( self.features )
                error( 'Must run spike sorting first! See "ChannelIndex" for more info.' );
            end
            
            % plot the features
            scatter( self.features(:,column1),self.features(:,column2),'.' );
            xlabel( sprintf( 'feature %i',column1 ) );
            ylabel( sprintf( 'feature %i',column2 ) );
            set( gca,'tickdir','out','box','off' );
        end

    end %methods
    
end
        
        
            