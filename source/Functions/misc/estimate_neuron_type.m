function [type,peak_trough,hw] = estimate_neuron_type( X,fs,ID )
	% [type,peak_trough,hw] = estimate_neuron_type( X,fs,ID )
	%
	% given an N x M data matrix X of spike waveforms, sampling rate fs, and 
	% 1 x M vector of neuron labels ID, estimates which neurons are inhibitory vs. excitatory 
	% by clustering the peak-to-trough voltage, and half width (inter neurons should have a 
	% smaller peak-to-trough and half width)

	% get half width and peak-to-trough
	[nPts,~] = size( X );
	[spMin,spMax] = spikeHeight( X,fs,0.0001,nPts / fs);
	hw = halfWidth( X,spMin,fs );
	peak_trough = abs( spMin ./ (spMax - spMin) );

	% cluster 
	labels = kmeans( zscore( [hw,peak_trough] ),2 );

    % now loop over IDs and determine probability of being type 1 or 2
    % as the proportion of spikes within cluster 1 vs 2
    uniqueID = unique( ID );
    type = zeros( numel( uniqueID ),2 );
    for i = 1:numel( uniqueID )
        thisNeuron = ismember( ID,uniqueID(i) );
        type(i,1) = nnz( labels(thisNeuron)==1 ) / nnz( thisNeuron );
        type(i,2) = 1 - type(i,1);
    end
end

