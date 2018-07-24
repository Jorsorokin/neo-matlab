function [newSpikes,newID,parentSpike] = resolve_overlapping_spikes( spikes,sptimes,trials,ID,A,W )
    % [newSpikes,newID,parentSpike] = resolve_overlapping_spikes( spikes,sptimes,ID,A,W )
    %
    % resolves overlapping spikes by estimating the contributions of the
    % templates W and V using coefficients A and B to the overlapping
    % waveforms. W and V are created from "create_spike_templates.m" and 
    % ID, A, and B are the output of "template_match". Overlapping spikes
    % are those with more than 1 non-zero element in ID (i.e. each row is
    % one spike, columns represent contributions from the jth template)
    
    % get # of overlapping parent spikes
    nID = sum( ID > 0,2 );
    isOverlap = find( nID > 1 )';
    nOverlap = sum( nID(isOverlap) );
    
    [n,m] = size( spikes );
    newSpikes = zeros( n,nOverlap,class( spikes ) );
    parentSpike = zeros( 1,nOverlap );
    newID = zeros( 1,nOverlap,class( ID ) );
    refractoryViolation = false( 1,nOverlap );
    counter = 1;
    if isrow( trials )
        trials = trials';
    end
    
    for j = isOverlap
        templateInds = find( ID(j,:) > 0 );
        for k = templateInds
            
            % check if this same neuron has fired within the refractory 
            % period of this spike. If so, we're double-counting
            thisTrial = trials == trials(j);
            thisID = ID(:,1) == ID(j,k);
            idx = thisTrial & thisID;
            idx(j) = false;
            refractoryViolation(counter) = any( abs( sptimes(j) - sptimes(idx) ) < 0.001 );
            
            % extract this new spike by subtracting template: X - sum( A_i*W_i )
            parentSpike(counter) = j;
            otherTemplateInds = templateInds(templateInds ~= k);
            template = W(:,ID(j,otherTemplateInds)) * A(j,otherTemplateInds)';
            newSpikes(:,counter) = spikes(:,j) - template;
            newID(counter) = ID(j,k);
            counter = counter + 1;
        end
    end
    
    % remove those spikes that violated the refractory period (duplicate spikes)
    parentSpike(refractoryViolation) = [];
    newSpikes(:,refractoryViolation) = [];
    newID(refractoryViolation) = [];    
end
            
        
        
        
        