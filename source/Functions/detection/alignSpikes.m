function snips = alignSpikes( spikes,Fs,peakTime,pretime,posttime )
%
% snips = alignSpikes( spikes,Fs,peaktime,pretime,posttime )
%
% Align extracellular spike waveforms to their peaks by extracting
% pretime:posttime amount of data around the spike peaktime
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
[n,nSp,nChan] = size( spikes );
spikes = double( spikes );
prepts = floor( pretime * Fs );
postpts = floor( posttime * Fs );
totalpts = prepts + postpts;
if totalpts > n
    error( 'total points requested greater than the amount available' ); 
end

% preallocate
upspikes = zeros( n*4,nSp,nChan );
upFs = Fs * 4;
uPrePts = prepts * 4;
uPostPts = postpts * 4;
x = 1:n;
xx = linspace(1,n,n*4); % the interpolant x-values
warning( 'off' );
for c = 1:nChan
    try
        upspikes(:,:,c) = csaps( x,spikes(:,:,c)',0.3,xx )'; % smooths the interpolant
    catch
    end
end

% find the minimum across channels/points
%[vMin,tMin] = min( upSpikes,[],1 );
[vMin,~,tMin] = spikeHeight( upspikes,upFs,pretime,posttime );
[~,tmin2] = min( vMin,[],2 ); % across channels
for j = 1:numel( tmin2 )
    tmin2(j) = peakHeight(j,tmin2(j));
end

% align the spikes
usnips = interpolate_spikes( upspikes,round( tmin2*upFs ),uPrePts,uPostPts );

% downsample the spikes
snips = zeros( size(usnips,1)/4,size(usnips,2),nChan );
n = size( snips,1 );
xx = 1:n; % the interpolant x-values
x = linspace(1,n,n*4); 
for c = 1:nChan
    try
        snips(:,:,c) = spline( x,usnips(:,:,c)',xx )'; % no smoothing
    catch
    end
end

warning( 'on' );

end