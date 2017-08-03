function [wIC,A,W,IC] = wICA(data,varargin)
%--------------- function [wIC,A,W] = wICA(data,varargin) -----------------
%
% Performs ICA on data matrix (row vector) and subsequent wavelet
% thresholding to remove low-amplitude activity from the computes ICs.
% This is useful for extracting artifact-only ICs in EEG (for example), and
% then subtracting the artifact-reconstruction from the original data. 
%
%               >>> INPUTS >>>
% Required: 
%   data = data matrix in row format
% Optional:
%   type = "fastica" or "radical"...two different ICA algorithms based on
%       entropy. "fastica" (default) is parametric, "radical" is nonparametric.
%   mult = threshold multiplier...multiplies the computed threshold from
%       "ddencmp" by this number. Higher thresh multipliers = less
%       "background" (or low amp. EEG) is kept in the wICs.
%   plotting = 1 or 0. If 1, plots wIC vs. non-wavelet thresholded ICs
%   Fs = sampling rate, (for plotting...default = 1);
%   L = level set for stationary wavelet transform. Higher levels give
%       better frequency resolution, but less temporal resolution. 
%       Default = 5
%   wavename = wavelet family to use. type "wavenames" to see a list of
%       possible wavelets. (default = "sym5");
%   nIC = the number of IC's to extract
%
%               <<< OUTPUTS <<<
%   wIC = wavelet-thresholded ICs
%   A = mixing matrix (inv(W)) (optional)
%   W = demixing matrix (inv(A)) (optional)
%   IC = non-wavelet ICs (optional)
%   
%       * you can reconstruct the signals as: 
%               signals = A*wIC;
%       - upon reconstruction, you can then subtract the signal from your
%       original data set to remove artifacts, for instance
%
% By JMS, 11/10/2015
%---------------------------------------------------------------------------------------

% check inputs
if nargin>1 && ~isempty(varargin{1})
type=varargin{1}; else type='fastica';end
if nargin>2 && ~isempty(varargin{2})
mult=varargin{2};else mult=1;end
if nargin>3 && ~isempty(varargin{3})
plotting=varargin{3}; else plotting=0;end
if nargin>4 && ~isempty(varargin{4})
Fs=varargin{4};else Fs=1;end
if nargin>5 && ~isempty(varargin{5})
L=varargin{5}; else L=5;end
if nargin>6 && ~isempty(varargin{6})
wavename=varargin{6}; else wavename='sym5';end
if nargin>7 && ~isempty(varargin{7}) && varargin{7} <= size(data,2)
nIC = varargin{7}; else nIC = size(data,2);end

% run ICA using "fastica" or "radical"
if strcmp(type,'fastica')
    [IC,A,W] = fastica(data,'approach','defl','g','pow3','displayMode','off','numofIC',nIC); % fastica for parametric...default "pow3" nonlinearity
elseif strcmp(type,'radical')
    [IC,W] = radical(data); % radical ICA for non-parametric
    A = inv(W);
end

% padding data for proper wavelet transform...data must be divisible by
% 2^L, where L = level set for the stationary wavelet transform
modulus = mod(size(data,2),2^L); %32 = 2^5 = 2^level (level for wavelet)
if modulus ~=0
    extra = zeros(1,(2^L)-modulus);
else
    extra = [];
end
      
% loop through ICs and perform wavelet thresholding
disp('Performing wavelet thresholding');
[n,m] = size(IC);
wIC = zeros(n,m+numel(extra));
for s = 1:n
    if ~isempty(extra)
        sig = [IC(s,:),extra]; % pad with zeros
    else
        sig = IC(s,:);
    end
    [thresh,sorh,~] = ddencmp('den','wv',sig); % get automatic threshold value
    thresh = thresh*mult; % multiply threshold by scalar
    swc = swt(sig,L,wavename); % use stationary wavelet transform (SWT) to wavelet transform the ICs
    Y = wthresh(swc,sorh,thresh); % threshold the wavelet to remove small values
    wIC(s,:) = iswt(Y,wavename); % perform inverse wavelet transform to reconstruct a wavelet IC (wIC)
    clear y sig thresh sorh swc 
end

% remove extra padding
if ~isempty(extra)
    wIC = wIC(:,1:end-numel(extra));
end

% plot the ICs vs. wICs
if plotting>0
    disp('Plotting');
    subplot(3,1,1);
        multisignalplot(IC,Fs,'r');
        title('ICs');
    subplot(3,1,2);
        multisignalplot(wIC,Fs,'r');
        title('wICs')
    subplot(3,1,3);
        multisignalplot(IC-wIC,Fs,'r');
        title('residuals');
end

end