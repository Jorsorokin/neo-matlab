function snips = alignSpikes( spikes,Fs,peaktime,pretime,posttime )
%
% snips = alignSpikes( spikes,Fs,peaktime,pretime,posttime )
%
% Align extracellular spike waveforms to their peaks by extracting
% pretime:posttime amount of data around the spike peaktime. 
%
% Inputs:
%   spikes:
%       n x m x c spike-waveform matrix. "c" represents the channels, "n"
%       the number of points, and "m" the number of spikes
%   Fs:
%       scalar sampling rate (in Hz)
%   peaktime:
%       the time (in seconds) of the peak voltage
%   pretime:
%       the time (in seconds) to extract prior to the peak after aligning
%   posttime:
%       the time (in seconds) to extract following the peak after aligning
%
% Outputs:
%   alignedspikes:
%       the k x m x c aligned spike-waveform matrix, with k <= n
%
% By JMS, 5/20/2017


% check inputs and spike matrix size
[n,m,c] = size( spikes );
peakpts = round( peaktime * Fs );
prepts = floor( pretime * Fs );
postpts = floor( posttime * Fs );
totalpts = prepts + postpts;
if totalpts > n
    error( 'total points requested greater than the amount available' ); 
end

% preallocate
snips = nan( totalpts,m,c );

for j = 1:m
    for ch = 1:c
        try
            snips(:,j,ch) = spikes(peakpts(j,ch)-prepts:peakpts(j,ch)+postpts-1,j,ch);
        catch
        end
    end
end

end