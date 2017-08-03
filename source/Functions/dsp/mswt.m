function [sc,fband] = mswt(X,wtype,level,varargin)
% ========= [sc,fband] = mswt(data,wtype,level) =========
%
% Compute the stationary wavelet transform of each column in the matrix X.
% Stationary wavelet transform is an undecimated version of the discrete
% wavelet transform, which offers a translation-invariant wavelet
% decomposition by elmiinating the downsampling at each level. 
%
%               >>> INPUTS >>>
%   X       :   m x n column-oriented matrix
%   wtype   :   the wavelet to use (see "waveinfo")
%   level   :   the decomposition level
%   Fs      :   OPTIONAL...the sampling rate, if provided, will calculate
%               the actual frequency bands represented by the decomposition
%
%               <<< OUTPUTS <<<
%   sc      :   m x level+1 x n tensor of detail coefficients and
%               approximation coefficients (last column)
%   fband   :   OPTIONAL...frequency bands represented. If Fs is not
%               supplied, will be "nan"
%
% By JMS, 09/26/2016
% =============================================================

%% Check inputs
if nargin > 3 && ~isempty(varargin{1})
    Fs = varargin{1};
else Fs = nan; end

m = size(X,1);
n = size(X,2);

% determine if padding is necessary 
pad = abs(m - ceil(m/(2^level))*(2^level));
if pad > 0
    X = vertcat(X,zeros(pad,n));
end

% preallocate 
sc = zeros(m,level+1,n);

% perform swt over each column
for j = 1:n
    c = swt(X(:,j),level,wtype);
    sc(:,:,j) = c(:,1:m)';
end

% if Fs provided, supply the frequency bands represented by the transform
fband = nan(level+1,2); 
if ~isnan(Fs)
    for i = 1:level
        fband(i,1) = Fs/2/(2^i);
        fband(i,2) = Fs/2/(2^(i-1));
    end
    fband(end,1) = 0;
    fband(end,2) = Fs/2/(2^i);
end