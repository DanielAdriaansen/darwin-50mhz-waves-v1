function [invec] = runmean(invec,n)
%%%%%%%%%%%%%%%%%%%%%%%%% Header %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: runmean.m
%
% Description: Computes running mean of 1D (vector) timeseries
% 
% Input: 1. Vector of data
%        2. N (number of points in mean)
%
% Output: Low pass filtered vector of data
%
% Author: Dan Adriaansen
%
% Date: 26 Jul 2010
%
% Notes: This is modeled after runmean_old.m from my thesis code
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% old interp method %
offset = ceil((n-1)/2);
%offset = (n-1)/2;
limit = offset+1;

for d=1:length(invec)
    if d<limit
        invec(d)=invec(d);
    elseif d>length(invec)-limit
        invec(d)=invec(d);
    else
        invec(d)=mean(invec(d-offset:d+offset));
    end
end

return