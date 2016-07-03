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
%ncpath = '/d1/dadriaan/paper/data/maskedmin';
ncpath = '/d1/dadriaan/paper/data/maskedminbad';

% What level do we want the ST output for?
lev = 3000;

% Monsoon days
mbeg = 13; % NOTE- actually day 14, but day 0 = day 1 on zpanel plot
mdays = 20;
%mend = 32; % NOTE- actually day 33, but day 0 = day 1 on zpanel plot

% Break days
bbeg = 36; % NOTE- actually day 37, but day 0 = day 1 on zpanel plot
bdays = 23;
%bend = 58; % NOTE- actually day 59 but day 0 = day 1 on zpanel plot

% What hour is the beginning of a day? In Darwin, we will use 02Z - 02Z, or 1130 - 1130 local time
beghr = 2;

% Are we processing the break or monsoon?
bm = 'break';

% What period do we want to break on?
bper = 15;

% What's the frequency bin interval we want?
fint = 0.5/10.0;

% Calculate the frequency binss

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
  mw = ncread([ncpath,'/',flist(f).name],'mask_w');
  %size(w)

  % Determine the begin and end of the indexes we're storing
  end_ind = 1440*f;
  beg_ind = end_ind-1439;

  % Store the data in the matrix
  var(:,beg_ind:end_ind) = w;
  pflag(beg_ind:end_ind) = pf;
  time(beg_ind:end_ind) = t;
  mask_w(:,beg_ind:end_ind) = mw;

end

% Turn all the bad data to NANs
%miss = find(var==-99);
%var(miss) = nan;

% Determine the index of the beginning and the end of the monsoon and break period
begmonsoon = (beghr*60)+(1440*mbeg)+1;
fprintf(['begmonsoon = ',num2str(begmonsoon)])
fprintf(['\n'])
fprintf(['mbegunix = ',num2str(time(begmonsoon))])
fprintf(['\n'])
endmonsoon = begmonsoon+(1440*(mdays-1));
fprintf(['endmonsoon = ',num2str(endmonsoon)])
fprintf(['\n'])
fprintf(['mendunix = ',num2str(time(endmonsoon))])
fprintf(['\n'])
ndayssoon = ((endmonsoon-begmonsoon)/1440);
fprintf(['ndayssoon = ',num2str(ndayssoon)])
fprintf(['\n'])
begbreak = (beghr*60)+(1440*bbeg)+1;
fprintf(['begbreak = ',num2str(begbreak)])
fprintf(['\n'])
fprintf(['bbegunix = ',num2str(time(begbreak))])
fprintf(['\n'])
endbreak = begbreak+(1440*(bdays-1));
fprintf(['endbreak = ',num2str(endbreak)])
fprintf(['\n'])
fprintf(['bendunix = ',num2str(time(endbreak))])
fprintf(['\n'])
ndaysbreak = ((endbreak-begbreak)/1440);
fprintf(['ndaysbreak = ',num2str(ndaysbreak)])
fprintf(['\n'])

% Set the indexes to subset with based on user arguments
if strcmp(bm,'monsoon')
    sub_beg = begmonsoon;
    sub_end = endmonsoon;
else
    sub_beg = begbreak;
    sub_end = endbreak;
end

%imagesc(var,[-0.5,0.5]);
%set(gca,'YDir','normal')

% At the level that was requested, find all of the good periods for the regime requested
dslice = var(zindex,sub_beg:sub_end);
tslice = time(sub_beg:sub_end);
mslice = mask_w(zindex,sub_beg:sub_end);
%dslice = var(zindex,:);
%tslice = time(:);
%mslice = mask_w(zindex,:);

% Find bad data
%badw = find(mslice>2); %% PRECIP ONLY
badw = find(mslice>1); %% PRECIP + BAD

% Counter for the number of periods processed at the level requested
periodcount = 0;

% Loop over the data and find info about the periods
for p=1:length(badw)-1
    dnt = badw(p+1)-badw(p);
    if dnt > 1
        gbeg = badw(p)+1;
        gend = badw(p+1)-1;
        gdiff = gend-gbeg;
        if gdiff == 0
            %fprintf(['\nFOUND =0'])
            %fprintf(['\nbadf = ',num2str(badw(p))])
            %fprintf(['\nbadf+1 = ',num2str(badw(p+1))])
            %fprintf(['\ngbeg = ',num2str(gbeg)])
            %fprintf(['\ngend = ',num2str(gend)])
            nmin = 1;
            nhrs = 0;
        else
            %fprintf(['\nFOUND >0'])
            %fprintf(['\nbadf = ',num2str(badw(p))])
            %fprintf(['\nbadf+1 = ',num2str(badw(p+1))])
            %fprintf(['\ngbeg = ',num2str(gbeg)])
            %fprintf(['\ngend = ',num2str(gend)])
            nmin = mod(gdiff,60);
            nhrs = floor(gdiff/60);
        end
        fprintf(['\n'])
        fprintf(['\nLENGTH OF PERIOD = ',num2str(nhrs),' HRS ',num2str(nmin),' MIN'])
        fprintf(['\nPER BEG IDX = ',num2str(gbeg)])
        fprintf(['\nPER END IDX = ',num2str(gend)])
        fprintf(['\nBEG TIME = ',datestr(tslice(gbeg)/86400+datenum(1970,1,1))])
        fprintf(['\nEND TIME = ',datestr(tslice(gend)/86400+datenum(1970,1,1))])
        fprintf(['\n'])
        
        % Note that there is a chance, since we're subsetting for monsoon and break independently of the 4hr chunk finder code,
        % that by subsetting we could reduce either the first, last, or both the first and last chunk to some length of time less
        % than four hours. But is this possible? Might want to parameterize the minimum chunk length at the top.
        
        % Before we go on to the next period, extract the vector of data we want to examine in the ST and check it for NAN
        stvec = dslice(gbeg:gend);
        nmiss = length(find(isnan(stvec)));
        if nmiss > 0
            fprintf(['\n=================> ERROR! ',num2str(nmiss),' MISSING AT THIS LEVEL']);
            break
        else
            [str,stt,stf] = st(stvec);
        end
                
        % Advance the period counter
        periodcount = periodcount + 1;
    end
end

fprintf(['\nNUM PERIODS PROCESSED = for ',bm,' is ',num2str(periodcount)])
fprintf(['\n'])
