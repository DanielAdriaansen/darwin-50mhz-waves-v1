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

% What level do we want the ST output for?
lev = 3000;

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
fprintf(['\nSIZE OF w MATRIX:'])
size(var)

% Vector to store time
ut = zeros(1,nt*nfiles);
fprintf(['\nSIZE OF time VECTOR:'])
size(ut)

% Vector to store the precip flag
pflag = zeros(1,nt*nfiles);
fprintf(['\nSIZE OF pflag VECTOR:'])
size(pflag)

% Store the height array
agl = ncread([ncpath,'/',flist(1).name],'pagl');

% Find the height index
zindex = find(agl==lev);

% Loop over each file, open the data and store it
for f=1:nfiles

  % What file are we reading?
  fprintf(['\n',ncpath,'/',flist(f).name,'\n'])

  % Read in the data
  w = ncread([ncpath,'/',flist(f).name],'omegpass2');
  pf = ncread([ncpath,'/',flist(f).name],'precipflag');
  t = ncread([ncpath,'/',flist(f).name],'unix_time');
  %size(w)

  % Determine the begin and end of the indexes we're storing
  end_ind = 1440*f;
  beg_ind = end_ind-1439;

  % Store the data in the matrix
  var(:,beg_ind:end_ind) = w;
  pflag(beg_ind:end_ind) = pf;
  time(beg_ind:end_ind) = t;

end

%imagesc(var,[-0.5,0.5]);
%set(gca,'YDir','normal')

% Find the good times we want to do the ST for (i.e. NOT precip)
goodtimes = find(~pflag);

% Loop over the valid periods and find the beginning and end of each
pbeg = goodtimes(1);
for p=2:length(goodtimes)
    dt = goodtimes(p)-goodtimes(p-1);
    if dt > 1
        % We've found the end of a good period. Set the end index and calculate some info about the period
        fprintf(['\n']);
        pend = goodtimes(p-1);
        pdt = pend-pbeg;
        nmin = mod(pdt,60);
        nhrs = int8(pdt/60);
        dbeg = datestr(time(pbeg)/86400+datenum(1970,1,1));
        dend = datestr(time(pend)/86400+datenum(1970,1,1));
        fprintf(['\nLENGTH OF PERIOD = ',num2str(nhrs),' HRS ',num2str(nmin),' MIN'])
        fprintf(['\nPER BEG IDX = ',num2str(pbeg)])
        fprintf(['\nPER END IDX = ',num2str(pend)])
        fprintf(['\nBEG TIME = ',dbeg])
        fprintf(['\nEND TIME = ',dend])
        
        % Before we go on to the next period, extract the vector of data we want to examine in the ST and check it for NAN
        stvec = var(zindex,pbeg:pend);
        nmiss = length(find(isnan(stvec)));
        if nmiss > 0
            fprintf(['\n===========> ERROR! ',num2str(nmiss),' MISSING AT THIS LEVEL'])
        else
            [str,stt,stf] = st(stvec);
            size(str)
            min(stf)
            max(stf)
            %break;
        end
        
        % Set the start of the next period and move forward searching for the end
        pbeg = goodtimes(p);
        %break;
    end
end

fprintf(['\n'])