
% A single 24-hr 1-minute matrix to hold the collected data
meanbreak = zeros(length(fbins)-1,1440);

% A single 24 hr 1-minute matrix to hold sample size at each time
samplesize = zeros(1440);

% Create a vector of times to collect data for
dnstart = datenum(2006,1,1,2,0,0);
dnend = dnstart+datenum(0,0,0,0,0,86399);

% Take all of the times from the current season (from main script) and make a date vector
alldvec = datevec(dnslice);

% Loop while the current datenumber is <= the end datenumber, incrementing 60 seconds each time (once per minute)
dncurr = dnstart;
count = 1;
while dncurr <= dnend
    % Convert current datenumber to a date vector
    dvec = datevec(dncurr);
    
    % Store current hours and minutes
    HH = dvec(:,4);
    MM = dvec(:,5);
    
    % Find all the times where the time was exactly HH hours and MM minutes
    times = find(alldvec(:,4)==HH & alldvec(:,5)==MM);
    data = seasonpow(:,times);
    
    % Figure out how many of the times actually have good data at this HH:MM
    powgood = find(~isnan(seasonpow(1,times)));
    
    % Store the number of days with good data at this time
    if ~isempty(powgood)
        samplesize(count) = length(powgood);
        % Loop over the number of bins and do stuff
        for nb=1:length(fbins)-1
            % If there is more than 1, do statistics (mean, median) otherwise just store the value
            if length(powgood) > 1
                meanbreak(nb,count) = mean(data(nb,powgood));
            else
                meanbreak(nb,count) = data(nb,powgood);
            end
        end
    else
        samplesize(count) = 0;
    end
    
    % Go to the next time (increment the current datenumber by 60 seconds)
    dncurr = dncurr + datenum(0,0,0,0,0,60);
    
    % Increment the counter
    count = count + 1;
end

% Now plot the collected power
fw = [0,0,1200,700];
figure('visible','on','position',fw);

% Make the plot
imagesc(dnslice(1:1440),fbins,meanbreak);

% Reverse the Y Axis
set(gca,'YDir','normal');

% Add a colorbar
cbh = colorbar;

% Set the Y-axis Label
ylabel(cbh,'Average power (m/s)')

% Set the X-axis Label
xlabel('Time (UTC)');

% Turn the X-axis to readable times
datetick('x',15);

% Tighten the axes
axis tight;

% Adjust the colormap
caxis([0.0 0.2]);