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
        sortModel = struct( 'mu',[],'Sigma',[] );
        projMatrix = nan;
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
            % Each Neuron has a unique identifier specific to its parent
            % ChannelIndex that separates it and its spike times/waveforms 
            % from other Neurons within that parent ChannelIndex.
            % It can have multiple Spikes objects as children, each of which 
            % represents an Epoch from which APs belonging to the Neuron 
            % were detected.
            %
            % Children:
            %   Spikes
            %
            % Parents:
            %   ChannelIndex
            %
            % Properties:
            %   ID - the unique ID for this Neuron in this ChannelIndex group
            %   chanInd - the parent ChannelIndex number
            %   nSpikes - number of spikes across all Epochs
            %   nChan - the number of channels from which this Neuron was detected
            %   features - the n x m matrix of features used when sorting. 
            %              If no sorting has occured (or if the ID = 0 and 
            %              multiple re-sorts occurred), this equals "NaN"
            %   probabilities - the probabilities assigned to each spike of 
            %                   this Neuron to other cluster centers when sorting    
            %   sortModel - a structure containing the mean (mu),
            %               covariances (sigma) of theis neuron's cluster of spikes 
            %   projMatrix - the projection matrix used for the
            %                low-dimensional representation of the spikes
            %   featureMethod - the method used for computing "features" 
            %                   (i.e. kPCA, ICA, NPE ...)
            %
            % Methods:
            %   getSpikes
            %   plotSpikes
            %   raster
            %   psth
            %   getISI
            %   plotISI
            %   estimateKernel
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
                        child(j).chanInd = self.chanInd;
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
              
        
        function [snips,times,epNum] = getSpikes( self,varargin )
            % [snips,times,epNum] = getSpikes( self,(epochs) )
            %
            % extracts the spike voltage waveforms and spike times across 
            % epochs for this Neuron object. Returns snips, spike times,
            % and a vector "epNum" refering to which epoch (Spikes obj)
            % each voltage/spiketime came from. Optionally specify the
            % spike objects (i.e. epochs) to extract.
            
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
                snips = [];
                return;
            end
            
            % preallocate vectors
            npoints = size( child(1).voltage,1 );
            nSp = [child.nSpikes];
            nEpoch = numel( child );
            times = nan( max( nSp ),nEpoch );
            snips = nan( npoints,sum( nSp ),self.nChan );
            epNum = uint16( zeros( 1,sum( nSp ) ) );
            
            % loop over Spike objects
            counter = 0;
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
            
            cmap = colormap( jet(nChild) );
            
            % plot spike waveforms for each "Spikes" child
            for j = 1:nChild
                spChild(j).plot( cmap(j,:),chans );
            end
            suptitle( sprintf( 'Unit %i',self.ID ) );
        end
        
        
        function ISI = getISI( self,varargin )
            % ISI = getISI( self,(epoch) )
            %
            % compute the inter-spike interval (ISI) for all spikes
            % contained within this neuron. ISI will be a large matrix,
            % with each column representing one spike train, and rows
            % representing the difference between consecutive spike times
            % in that spike trian. Differences are given in seconds. Can
            % optionally supply a second argument specifying which spike
            % trains (i.e. epochs) to pull out.
            
            % check inputs
            if nargin > 1 && ~isempty( varargin{1} )
                epochs = varargin{1};
            else
                epochs = 1:numel( self(1).getChild( 'Spikes' ) );
            end
            
            % pull out spikes
            [~,sptm] = self.getSpikes( epochs );
            
            % get the difference in the spike times
            ISI = diff( sptm );
        end
           
        
        function plotISI( self,varargin )
            % plotISI( self,(bw) )
            %
            % plot the inter-spike-interval (ISI) histogram for all spikes
            % contained in this Neuron. Optionally specify a bin width (in
            % seconds).
            
            % check input
            if nargin > 1
                bw = varargin{1};
            end
                
            % get the ISI
            ISI = self.getISI();
            
            % reshape the matrix
            [n,m] = size( ISI );
            ISI = reshape( ISI,n,m );
            ISI(isnan( ISI )) = [];
            
            % plot a histogram of the spike times
            if exist( 'bw','var' )
                [N,edges] = histcounts( ISI,'binwidth',bw,...
                    'Normalization','countdensity' );
            else
                [N,edges] = histcounts( ISI,'Normalization','countdensity' );
            end
            stairs( edges(1:end-1),N );
            set( gca,'tickdir','out','box','off' );
            ylabel( 'count / bin width' );
            xlabel( 'ISI (s)' );
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
            suptitle( sprintf( 'Neuron %i, ChannelIndex %i',...
                                self.ID,self.chanInd ) );
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


        function [width,sigma,kernel] = estimateKernel( self )
            % [width,sigma,kernel] = estimateGaussKernel( self )
            %
            % estimate the appropriate width and sigma for the gaussian kernel
            % for convolving spike times to compute a firing rate.
            
            width = nan;
            sigma = nan;
            kernel = nan;
            
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
            spikes = self.getChild( 'Spikes' );
            nSp = numel( spikes );
            if isempty( spikes )
                assert( 'No spikes found' );
            end

            % normalize the kernel
            kernel = kernel / norm( kernel );

            % get the Epoch parents of the Spikes and Signal children of the ChannelIndex
            epochs = self.getPartner( 'Epoch','Spikes' );
            signals = self.getSibling( 'Signal','ChannelIndex' );
            duration = [epochs.duration]; 
            npoints = ceil( duration .* [signals.fs] ); 

            % pre-allocate our firing-rate matrix based on max number of points
            rate = zeros( max( npoints ),nSp ); 

            % loop over the Spikes children, extract the spike times
            for ep = 1:nSp
                sptm = round( spikes(ep).times * signals(ep).fs );
                rate(sptm,ep) = 1;
            end

            % convolve with the kernel
            rate = conv2( rate,kernel,'same' );
        end

               
        function plotFeatures( self )
            % plotFeatures( self,column1,column2 )
            %
            % plots a gplotmatrix of the columns of "features"
            % specified by the inputs "columns1" and "column2".
            % Must have already sorted spikes for this function to work
            if isnan( self.features )
                disp( 'Must run spike sorting first! See "ChannelIndex" for more info.' );
                return
            end
            
            % plot the features
            for j = 1:size( self.features,2 )
                axlabel{j} = sprintf( 'feature %i',j );
            end
            figure;
            gplotmatrix( self.features,[],[],'k',[],3,[],[],axlabel );
            chind = self.getParent( 'ChannelIndex' );
            name = chind.name; 
            if isempty( name )
                name = ['channelindex ',num2str(chind.chanIndNum)];
            end
            suptitle( sprintf( 'neuron %i, from %s',self.ID,name ) );
        end

    end %methods
    
end
        
        
            