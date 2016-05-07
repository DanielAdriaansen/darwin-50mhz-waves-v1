%
%
% File: st_single_level.m
%
% Author: D. Adriaansen
%
% Date: 07 May 2016
%
% Purpose: Read in pre-processed 50MHz data and perform the S-transform at a single height
%
% Notes:
%_________________________________________________________________________________________

% Clear the workspace
clear;

%%%%%%%%%%%%%%%%%%%%%%% User Config %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Path to netCDF data
ncpath = '/d1/dadriaan/paper/data/maskedmin';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get the list of files we want to process
flist = dir([ncpath,'/*.nc']);

% Number of files
nfiles = length(flist);
fprintf(['\nEXAMINING: ',num2str(nfiles),' FILES.']);

% Read in the first file and get the dimensions we need, then define a new matrix to hold the data
nz = length(ncread([ncpath,'/',flist(1).name],'pagl'));
nt = length(ncread([ncpath,'/',flist(1).name],'unix_time'));

% Based on the number of files, times, and heights create the correctly sized matrix for the data
var = zeros(nz,nt*nfiles);
fprintf(['\nSIZE OF MATRIX:'])
size(var)

% Loop over each file, open the data and store it
for f=1:nfiles

  fprintf(['\n',ncpath,'/',flist(f).name,'\n'])

  w = ncread([ncpath,'/',flist(f).name],'omegpass2');

end

exit;
