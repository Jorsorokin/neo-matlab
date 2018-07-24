classdef Neuron < Container
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
    %   meanWaveform - average spike waveform of all spike children
    %   bestElectrode - the electrode with the largest voltage deflection
    %                   of the meanWaveform
    %   location - the physical distance of the bestElectrode
    %   region - the brain region from which this neuron was recorded
    %   subRegion - the sub brain region (i.e. VPM of the thalamus)
    %
    % Methods:
    %   getSpikes
    %   getSpikes_subset
    %   plotSpikes
    %   raster
    %   psth
    %   getISI
    %   plotISI
    %   getCorrelogram
    %   estimateKernel
    %   firingRate
    %   plotFeatures
    %   
    %       * see also methods in the Container class
    
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
        meanWaveform
        bestElectrode
        location
        region
        subRegion
    end
    
    methods
        
        function self = Neuron( ID )
            self.ID = ID;
        end
        
        
        function addChild( self,child )
            % adds spikes objects to this Neuron
            assert( isa( child,'Spikes' ),'Only Spikes objects are valid children' );
            
            addChild@Container( self,child );
            self.nSpikes = sum( [child.nSpikes] );
            meanVolt = zeros( size( child(1).voltage,1 ),size( child(1).voltage,3 ) );
            
            for j = 1:numel( child )
                child(j).unitID = self.ID;
                child(j).chanInd = self.chanInd;
                if child(j).nSpikes > 0
                    meanVolt = meanVolt + squeeze( mean( child(j).voltage,2 ) );
                end
            end
            self.meanWaveform = meanVolt / j;
        end 
        
        
        function addParent( self,parent )
            asert( isa( parent,'ChannelIndex' ),'Only ChannelIndex objects are valid parents' );
            
            addParent@Container( self,parent );
            self.chanInd = parent.chanIndNum; 
            self.nChan = parent.nElectrodes;
        end
              
        
        function [snips,times,epNum,mask] = getSpikes( self,varargin )
            % [snips,times,epNum,mask] = getSpikes( self,(epochs,bounds) )
            %
            % extracts the spike voltage waveforms and spike times across 
            % epochs for this Neuron object. Returns snips, spike times,
            % a vector "epNum" refering to which epoch (Spikes obj)
            % each voltage/spiketime came from, and a sparse matrix "mask"
            % refering to the spike mask created if running masked spike detection.
            % 
            % One can optionally specify the spike objects (i.e. epochs) to extract by 
            % supplying an integer or vector of integers referring to the desired epochs.
            %
            % One can also supply a 2-element vector specifying the time bounds (in seconds)
            % to extract. I.e. bounds = [2 5] extracts spikes who's spike times are within 2 and 5 
            % seconds of the recording for each epoch. 
            
            if self.nSpikes == 0
                disp( 'no spikes found' );
                snips = [];
                times = [];
                epNum = [];
                mask = [];
                return
            end

            % get spikes objects
            if nargin > 1 && ~isempty( varargin{1} )
                spikes = self.getChild( 'Spikes',varargin{1} );
            else
                spikes = self.getChild( 'Spikes' );
            end
            
            % pull out actual variables
            nSp = [spikes.nSpikes];
            times = [spikes.times];
            mask = [spikes.mask];
            snips = [spikes.voltage];
            epNum = cell2mat( arrayfun( @(x,y)(ones(1,x)*y),nSp,1:numel( spikes ),'un',0 ) );
            
            % remove spikes outside of bounds
            if nargin > 2 && ~isempty( varargin{2} )
                idx = (times < varargin{2}(1) | times > varargin{2}(2));
                times(idx) = [];
                epNum(idx) = [];
                mask(:,idx) = [];
                snips(:,idx,:) = [];
            end
            
            % convert to spike time matrix
            times = spikevec2mat( times,epNum );
        end
        
        
        function [snips,time,epoch,mask,amp,bestElectrode] = getSpikes_subset( self,nTotal,varargin )
            % [snips,time,epoch,mask,amp,bestElectrode] = getSpikes_subset( nTotal,(epochs,bounds,separateByChannels,minAmp) )
            %
            % extracts a subset of spikes contained within this neuron object. If the optional argument
            % separateByChannels, a boolean, is not provided or is false, then a random subset across all
            % spikes is pulled out. Otherwise, nTotal, the # of spikes to extract, is divided by the # of
            % electrodes recording the spikes, and for each electrode, a random subset of spikes is
            % extracted that have their largest voltage fluctuation on that electrode.
            %
            % the optional argument "epochs" specifies which epochs to pull spikes from. If not
            % specified, all epochs are used. 
            %
            % the optional argument "bounds", a two-element vector, specifies (in seconds) when to extract spikes.
            % If not specified, the random subset is pulled out across the entire dataset. Else, the random subset
            % is only pulled out if: bounds(1) <= spike times <= bounds(2)
            %
            % finally, the optional argument "minAmp" specifies an amplitude threshold for the spike subset.
            % Those spikes with voltages > minAmp (meaning, less of a voltage deflection), are not
            % extracted. Default is minAMp = 0
            %
            % For all optional arguments, use the name-value pair format
                                               
            % get optional inputs
            p = check_inputs( varargin );
            [snips,time,epoch,mask] = self.getSpikes( p.epochs,p.bounds );
            nTotal = min( nTotal,numel( epoch ) );
            time = spikemat2vec( time );
            
            % pull out subset of spikes for each electrode where spikes have largest voltage
            % on that electrode
            if ~isempty( mask )
                [amp,bestElectrode] = min( squeeze( min( maskchans( snips,mask ) ) ),[],2 );
            else
                [amp,bestElectrode] = min( squeeze( min( snips ) ),[],2 );
            end
            
            amp = amp';
            bestElectrode = bestElectrode';
            
            if ~p.separateByChannels
                idx = find( amp <= p.minAmp );
                n = numel( idx );
                idx = idx( sort( randperm( n,min( n,nTotal ) ) ) );
            else
                nElectrodeSpikes = floor( nTotal / self.nChan );
                idx = nan( nElectrodeSpikes,self.nChan );
                for electrode = 1:self.nChan
                    tempidx = find( bestElectrode == electrode & amp <= p.minAmp );
                    n = numel( tempidx );
                    k = min( n,nElectrodeSpikes );
                    idx(1:k,electrode) = tempidx( randperm( n,k ) );
                end
                idx = sort( spikemat2vec( idx ) );
            end
            
            % now extract the subset of spikes
            bestElectrode = bestElectrode(idx);
            amp = amp(idx);
            snips = snips(:,idx,:);
            mask = mask(:,idx);
            epoch = epoch(idx);
            time = spikevec2mat( time(idx),epoch );
            
            function p = check_inputs( inputs )
                names = {'epochs','bounds','separateByChannels','minAmp'};
                defaults = {[],[],false,0};
                p = inputParser;
                for j = 1:numel(names)
                    p.addParameter( names{j},defaults{j} );
                end
                p.parse( inputs{:} );
                p = p.Results;
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
                return
            end
            
            if nargin < 2 || isempty( varargin{1} )
                chans = 1:spChild(1).nChan;
            else
                chans = varargin{1};
            end
            
            cmap = colormap( jet(nChild) );
            
            % plot spike waveforms for each "Spikes" child
            set( gcf,'Visible','off' );
            for j = 1:nChild
                spChild(j).plot( cmap(j,:),chans );
            end
            set( gcf,'Visible','on' );
            suptitle( sprintf( 'Unit %i',self.ID ) );
        end
        
        
        function ISI = getISI( self,varargin )
            % ISI = getISI( self,(epoch) )
            %
            % compute the inter-spike interval (ISI) for all spikes
            % contained within this neuron. ISI will be a large matrix,
            % with each column representing one spike train, and rows
            % representing the difference between consecutive spike times
            % in that spike train (in ms). 
            %
            % Can optionally supply a second argument specifying which spike
            % trains (i.e. epochs) to pull out.
            
            if nargin > 1 && ~isempty( varargin{1} )
                epochs = varargin{1};
            else
                epochs = 1:numel( self(1).getChild( 'Spikes' ) );
            end
            
            [~,sptm] = self.getSpikes( epochs );
            ISI = diff( sptm*1000 );
        end
           
        
        function plotISI( self,varargin )
            % plotISI( self,(epochs,bw) )
            %
            % plot the inter-spike-interval (ISI) histogram for all spikes
            % contained in this Neuron. Optionally specify a bin width (in
            % ms).
            
            % check input
            if nargin > 1
                epochs = varargin{1};
            else
                epochs = 1:numel( self.children{1} );
            end
            
            if nargin > 2
                bw = varargin{2};
            end
                
            % reshape the ISI matrix
            ISI = self.getISI( epochs );
            [n,m] = size( ISI );
            ISI = reshape( ISI,n,m );
            ISI(isnan( ISI )) = [];
            
            % plot a histogram of the spike times
            if exist( 'bw','var' )
                edges = linspace(0,max(ISI),round( max(ISI)/bw ));
                xscale = 'linear';
            else
                edges = logspace( log10(0.1),log10(max(ISI)),round( max(ISI)/20 ) );
                xscale = 'log';
            end
            N = histcounts( ISI,edges,'Normalization','Probability' );
            ax = gca();
            s = stairs( ax,edges(2:end),N );
            set( s,'color',[0.85 0.85 0.85] );
            set( ax,'xscale',xscale );
            ylabel( 'probability' );
            xlabel( 'ISI (ms)' );
            darkPlot( gcf );
        end
        
        
        function [xcg,lags] = getCorrelogram( self,binwidth,maxLag,plotFlag,varargin )
            % [xcg,lags] = getCorrelogram( self,binwidth,maxLag,plotFlag,(neuronID) )
            %
            % gets the auto (or cross) correlogram between this neuron and itself
            % (or another neuron, if another neuron object is provided as an optional input)
            % 
            % the variables binwidth and maxLag should be in seconds
            
            [~,train1] = self.getSpikes();            
            
            if nargin > 4 
                otherNeuron = findobj( self.getParent('ChannelIndex').getChild('Neuron'),'ID',varargin{1} );
                [~,train2] = otherNeuron.getSpikes();
                [xcg,lags] = correlogram( train1,train2,binwidth,maxLag );
            else
                [xcg,lags] = correlogram( train1,train1,binwidth,maxLag );
            end
            
            if plotFlag
                bar( lags*1000,xcg,'FaceColor','k','EdgeColor','None','BarWidth',1 );
                ylabel( 'pdf(x)' );
                xlabel( 'ms' );
                set( gca,'tickdir','out','box','off' );
            end
        end
        
        
        function raster( self,start,pre,post,varargin )
            % raster( self,start,pre,post,(epochs) );
            %
            % make a raster plot of the spikes across epochs associated
            % with this Neuron. Must provide "start", "pre", and "post"
            % relative to the spiketimes contained within the
            % Neuron. Assumes the trials (epochs) are aligned for
            % meaningful results
            
            % check for epochs
            if nargin > 4
                epochs = varargin{1};
            else
                epochs = 1:numel( self.children{1} );
            end

            [~,sptm,~] = self.getSpikes( epochs ); 
            PlotRasters( sptm,start,pre,post );
            title( sprintf( 'Neuron %i, ChannelIndex %i',...
                                self.ID,self.chanInd ) );
        end
        
        
        function [count,avgCount,sigma,time] = psth( self,start,pre,post,varargin )
            % [count,avgCount,sigma,time] = psth( self,start,pre,post,(bw,epochs) )
            %
            % compute the peri-stimulus time histogram (PSTH) from the
            % Neuron object across all epochs. Output is the bin-count
            % (count), trial-averaged count (avgCount), and variance across
            % trials (sigma)
            
            % check for bin width input
            if nargin > 4 && ~isempty( varargin{1} )
                bw = varargin{1};
            else 
                bw = [];
            end

            if nargin > 5 && ~isempty( varargin{2} )
                epochs = varargin{2};
            else
                epochs = 1:numel( self.children{1} );
            end
            
            [~,sptm,~] = self.getSpikes( epochs );
            [count,avgCount,sigma] = CalculatePSTH( sptm,start,pre,post,bw );
            time = linspace( start-pre,start+post,numel( avgCount ) );
        end


        function [sigma,kernel] = estimateKernel( self,varargin )
            % [sigma,kernel] = estimateKernel( self,(t,epochs) )
            %
            % estimate the best sigma for a gaussian kernel
            % for convolving spike times to compute a firing rate.
            % Can optionally supply the times from which to sample 
            % the spike-time vectors when estimating the kernel, as well as
            % which epochs to estimate the kernel from
            
            if nargin > 1 && ~isempty( varargin{1} )
                t = varargin{1};
            else
                t = [];
            end

            if nargin > 2 && ~isempty( varargin{2} )
                epochs = varargin{2};
            else
                epochs = 1:numel( self.children{1} );
            end

            % get spike times
            [~,sptimes] = self.getSpikes( epochs );

            % compute the kernel / sigma
            if isempty( t )
                [~,~,sigma] = sskernel( sptimes );
            else
                [~,~,sigma] = sskernel( sptimes,t );
            end

            % TO DO: create the kernel
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
            % of Epochs total in block
            
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
            electrode = self.getParent( 'ChannelIndex' ).getChild( 'Electrode',1 );
            fs = [electrode.getChild( 'Signal' ).fs];
            duration = [epochs.duration]; 
            npoints = ceil( duration .* fs ); 

            % pre-allocate our firing-rate matrix based on max number of points
            rate = zeros( max( npoints ),nSp ); 

            % loop over the Spikes children, extract the spike times
            for ep = 1:nSp
                sptm = round( spikes(ep).times * fs(ep) );
                rate(sptm,ep) = 1;
            end

            % convolve with the kernel
            rate = conv2( rate,kernel,'same' );
        end

               
        function plotFeatures( self )
            % plotFeatures( self )
            %
            % plots a gplotmatrix of the columns of "features"
            % specified by the inputs "columns1" and "column2".
            % Must have already sorted spikes for this function to work
            if isnan( self.features )
                disp( 'Must run spike sorting first! See "ChannelIndex" for more info.' );
                return
            end
            
            % plot the features
            axlabel = cell( 1,size( self.features,2 ) );
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
        
        
            