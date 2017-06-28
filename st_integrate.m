%exi
%
% File: st_integrate.m
%
% Author: D. Adriaansen
%
% Date: 03 Jul 2016
%
% Purpose: Read in pre-processed 50MHz data and perform the S-transform at a single height then
%          integrate across frequency bands.
%
% Notes:
%_________________________________________________________________________________________

% Clear the workspace
%clear;

%%%%%%%%%%%%%%%%%%%%%%% User Config %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Path to netCDF data
%ncpath = '/d1/dadriaan/paper/data/c2/maskedmin';
ncpath = '/d1/dadriaan/paper/data/c2/maskedminbad';
%ncpath = '/d1/dadriaan/paper/data/c3/maskedminbad';

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
bm = 'monsoon';

% What period do we want to break on?
bper = -9;

% Make plots or no? 1 = Yes, 0 = No
pmake = 1;

% Frequency bins for integrating
fbins = [0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5];

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
endmonsoon = begmonsoon+(1440*(mdays))-1;
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
endbreak = length(time);
%endbreak = begbreak+(1440*(bdays-1));
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

% Create a vector of MATLAB datenums for whichever (monsoon or break) we're processing
% NOTE: here tslice contains UNIX times for the entire period
dnslice = datenum(1970,1,1,0,0,tslice);

% Create a two-dimensional matrix the size of the period (x-axis) x number of bins (y-axis) to store power
% Also turn it all to NaN's so arithmetic later isn't messed up
seasonpow = nan(length(fbins)-1,length(dnslice));
%seasonpow(:,:) = nan;

% Find bad data
%badw = find(mslice>2); %% PRECIP ONLY
badw = find(mslice>1); %% PRECIP + BAD

% Counter for the number of periods processed at the level requested
periodcount = 1;

% Vector to hold all the stvecs
stboxplot = [];

% Vector to hold the identifiers for period number
pnums = [];

% Counter for number of bad chunks
badchunk = 0;

% Counter for the number of chunks processed
allchunk = 0;

% Total number of chunks, good or bad
nchunks = 0;

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
        fprintf(['\nLENGTH OF CHUNK = ',num2str(nhrs),' HRS ',num2str(nmin),' MIN'])
        fprintf(['\nCNK BEG IDX = ',num2str(gbeg+sub_beg)])
        fprintf(['\nCNK END IDX = ',num2str(gend+sub_beg)])
        fprintf(['\nBEG TIME = ',datestr(tslice(gbeg)/86400+datenum(1970,1,1))])
        fprintf(['\nEND TIME = ',datestr(tslice(gend)/86400+datenum(1970,1,1))])
        fprintf(['\n'])
        nchunks = nchunks + 1;
        
        % Note that there is a chance, since we're subsetting for monsoon and break independently of the 4hr chunk finder code,
        % that by subsetting we could reduce either the first, last, or both the first and last chunk to some length of time less
        % than four hours. But is this possible? Might want to parameterize the minimum chunk length at the top.
        
        % Before we go on to the next period, extract the vector of data we want to examine in the ST and check it for NAN
        stvec = dslice(gbeg:gend);
        tvec = tslice(gbeg:gend);
        nmiss = length(find(isnan(stvec)));
        if nmiss > 0
            fprintf(['\n=================> ERROR! ',num2str(nmiss),' MISSING AT THIS LEVEL']);
            break
        else
            [str,stt,stf] = st(stvec);
        end
       
        % Call the plotter for the panel plot, and also this script will clean the data and examine for subsets within the chunk ignoring bad data near beginning/end
        % Note that a call to this script modifies the stvec so that when it's used below it has been subset. If that effect is not desired, then
        % comment out the call to this script.
        %plot_filter_time_series;
        
        % Clean out original S-transform stuff
        clear str;
        clear stt;
        clear stf;
        
        % Clean up the data before we do anything. If there's any missing data in the first hour or last hour, just discard that portion
        % of the chunk. However, must still ensure chunk is > 4 hours. If the chunk still has missing data in the middle after checking the edges,
        % then just skip the chunk altogether (or also skip if discarding beginning/end reduces length to less than 4 hours).
        % STEP 1: Take a 10-point (10-minute) running mean
        stmean = runmean(abs(stvec),10);
        % STEP 2: Find outliers (magnitude > 0.5 m/s) and turn them to NaN
        out = find(stmean>0.5);
        stmean(out) = nan;
        % STEP 3: Find out whether the first 60 minutes or last 60 minutes have missing data in them
        if ~isempty(find(isnan(stmean(end-60:end))))
            fprintf(['\nDISCARDING END OF CHUNK'])
            stmean(end-60:end) = nan;
        end
        if ~isempty(find(isnan(stmean(1:61))))
            fprintf(['\nDISCARDING BEG OF CHUNK'])
            stmean(1:61) = nan;
        end
        % STEP 4: Check for any NAN's in the middle now. If there are NAN's continue to next chunk
        if ~isempty(find(isnan(stmean(62:end-61))))
            fprintf(['\nDISCARDING CHUNK. NANs FOUND IN MIDDLE'])
            badchunk = badchunk + 1;
            continue
        end
        % STEP 5: Make sure the chunk is still > 4 hours, and then cut off the nan sections
        goodind = find(~isnan(stmean));
        goodchunk = stvec(goodind);
        goodtime = tvec(goodind);
        if length(goodchunk) < 240
            fprintf(['\nDISCARDING CHUNK BECAUSE < 4 HOURS NOW'])
            badchunk = badchunk + 1;
            continue
        end
        allchunk = allchunk + 1;
        % STEP 6: Now do the S-Transform here using only the good part of the chunk
        [str,stt,stf] = st(goodchunk);
        % STEP 7: Reset variables used below
        clear tvec;
        tvec = goodtime;
        clear stvec;
        stvec = goodchunk;
        
        % Create a vector of MATLAB datenums to use for storing the integrated power
        % NOTE: Here tvec contains the UNIX times associated with only the current "chunk" we are processing
        dnvec = datenum(1970,1,1,0,0,tvec);
               
        % Vector to hold the number of voices per band
        nvoice = zeros(1,length(fbins)-1);
        
        % Matrix to hold the power in each band at each time
        totpow = zeros(length(fbins)-1,length(tvec));
                       
        % Now, loop over each time in the period and compute the total power in each frequency band.
        for t=1:length(tvec)
            %fprintf(['\nTIME = ',datestr(tvec(t)/86400+datenum(1970,1,1))]);
            for b=1:length(fbins)-1
                %fprintf(['\nINTEGRATING FROM ',num2str(fbins(b)),' TO ',num2str(fbins(b+1))]);
                % If it's the first band, ignore the zero frequency (0.0)
                if b==1  
                  group = find(stf>fbins(b) & stf<fbins(b+1));
                elseif b==length(fbins)-1
                  group = find(stf>=fbins(b) & stf<=fbins(b+1));                  
                else
                  group = find(stf>=fbins(b) & stf<fbins(b+1));
                end
                       
                % If it's the first time, store the number of voices per band
                if t==1
                  nvoice(b) = length(group);
                end
                
                % Now calculate total power in this band and normalize by the change in frequency which varies per chunk
                %totpow(b,t) = sum(abs(str(group,t)))*(0.5/length(stf));
                
                % Instead of summing/integrating, try using the average
                totpow(b,t) = sum(abs(str(group,t)))/nvoice(b);
            end
        end
        
        % Find the start end end in the seasonal power matrix of the current chunk
        startind = find(dnslice==dnvec(1));
        endind = find(dnslice==dnvec(end));
        
        % Store the chunk power in the seasonal array
        seasonpow(:,startind:endind) = totpow;
        
        if pmake==1
        
          % Before moving on, create a save a plot for this period
          fw = [0,0,900,700];
          figure('visible','off','position',fw);
        
          %clevs = [0.0,1.0,0.05];
          clims = [0.0 1.0];
        
          %imagesc(stt,stf,(abs(str).*abs(str)),clims);
          imagesc((tvec/86400+datenum(1970,1,1)),stf,(abs(str)),clims);
          %imagesc(totpow,clims);
        
          set(gca,'YDir','normal');
          cbh = colorbar;
          ylabel('Frequency 1/s');
          xlabel('TIme (UTC)');
          title({[bm,' chunk ',num2str(periodcount)],['Begin = ',datestr(tslice(gbeg)/86400+datenum(1970,1,1))],['End = ',datestr(tslice(gend)/86400+datenum(1970,1,1))],['Length =',num2str(nhrs),' HRS ',num2str(nmin),' MIN']})
          datetick('x',15);
          axis tight;
        
          %ylabel(cbh,'abs(str)^2 (m^2/s^2)');
          ylabel(cbh,'abs(str) (m/s)');
          %ylabel(cbh,'Average power (m/s)')
        
          %set(gca,'YTick',[0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5])
          %set(gca,'YTickLabel',[0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5])
        
          %saveas(gcf,['st_',bm,'_chunk_',num2str(periodcount),'_square.png']);
          saveas(gcf,['st_',bm,'_chunk_',num2str(periodcount),'_abs.png']);
          %saveas(gcf,['st_',bm,'_chunk_',num2str(periodcount),'_totpow.png']);
          
          % Time series of the input velocity
          figure('visible','off','position',fw);
          plot((tvec/86400+datenum(1970,1,1)),stvec);
          xlabel('Time (UTC)');
          ylabel('Velocity (m/s)');
          title({[bm,' chunk ',num2str(periodcount)],['Begin = ',datestr(tslice(gbeg)/86400+datenum(1970,1,1))],['End = ',datestr(tslice(gend)/86400+datenum(1970,1,1))],['Length =',num2str(nhrs),' HRS ',num2str(nmin),' MIN']})
          datetick('x',15);
          set(gca,'YLim',[-2.0 2.0]);
          %axis tight;
          saveas(gcf,['ts_',bm,'_chunk_',num2str(periodcount),'.png']);
        
        end
               
        % Break if we've reached the period we want
        if periodcount == bper
            break;
        end
        
        % Store the stvec
        if periodcount == 1
            stboxplot = [stvec'];
        else
            stboxplot = [stboxplot; stvec'];
        end
        %stboxplot = [stboxplot; stvec'];
        
        % Store the period number
        if periodcount == 1
            pnums = [(periodcount)*ones(size(stvec'))];
        else
            pnums = [pnums; ((periodcount)*ones(size(stvec')))];
        end
        %pnums = [pnums; periodcount(size(stvec'))];
        
        % Advance the period counter
        periodcount = periodcount + 1;

    end
end

fprintf(['\nNUM CHUNKS EXAMINED for ',bm,' is ',num2str(nchunks)])
fprintf(['\nNUM CHUNKS DISCARDED for ',bm,' is ',num2str(badchunk)])
fprintf(['\nNUM CHUNKS RETAINED for ',bm,' is ',num2str(allchunk)])
fprintf(['\n'])
