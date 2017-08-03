function [Etot,Erel] = WaveEntropy(C,L,varargin)
% --------------- [Etot,Erel] = WaveEntropy(C,L,varargin) -----------------
%
%   Computes the wavelet shannon entropy (see Rosso et al. 2001) on 
%   wavelet coefficients (from discrete wavelet transform)
%
%                   >>> INPUTS >>>
% Required:
%   C: wavelet coefficients (vector or matrix) assuming the
%       following formats:
%           vector: full vector represents coefficients from one time
%                       series 
%           matrix: each column j represents coefficients from dwt
%                       decomposition of time series j
%   L: bookeeping vector from wavedec
%
% Optional:
%   Cref: cell array of data channels, with each cell containing
%                     same structure as C. Assumes "L" is equal between
%                     refchan data and C.
%
%                   <<< OUTPUTS <<<
%   Etot: total wavelet shannon entropy for each time series
%           Etot = -sum(Energy.*log(Energy))
%               where "Energy" is a vector containing relative energy
%               contributions from the various decompositions in C
%               (calculated by "wenergy")
%
%   Erel: relative wavelet entropy between channels (if "Cref" provided)
%           each row = relative entropy between all columns C and one channel of Cref
%           each column = relative entropy between column j of C/Cref
%
% By JMS, 3/30/2016
%-------------------------------------------------------------------------------

% check inputs
if nargin>2 && ~isempty(varargin{1})
    Cref = varargin{1};
end

% calculate relative wavelet energies for each time series
for i = 1:size(C,2)
    [Ea(i),Ed(:,i)] = wenergy(C(:,i),L);
end

% Concatenate Ea,Ed and calculate total entropy for each time series
Energy = vertcat(Ea,Ed) ./ 100;
Etot = -sum(Energy.*log(eps+Energy));

% repeat for multichan data if it exists
if exist('Cref','var') && ~isempty(Cref)
    for j = 1:numel(Cref)
        for i = 1:size(Cref{j},2)
            [Ea_ref{j}(i),Ed_ref{j}(:,i)] = wenergy(Cref{j}(:,i),L); % calculate energy for each channel/time series
        end
        Energy_ref{j} = vertcat(Ea_ref{j},Ed_ref{j}) ./ 100 ; 

        % now compute the relative entropy between Energy and Energy_ref{j}
        Erel(j,:) = -sum(Energy .* log(Energy ./ Energy_ref{j}));
    end
end

end