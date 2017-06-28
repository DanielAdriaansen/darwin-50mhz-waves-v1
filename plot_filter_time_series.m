% Plot a 5 panel plot for each period
% a. Raw st time series
% b. Absolute value or square of time series
% c. 10 point (minute) moving average of (b)
% d. As in (c), but with outliers set to NaN (> 0.5 m/s)
% e. As in (a), but with outliers set to NaN (> 0.5 m/s)


% Convert to a MATLAB datenum to use the datetick for plotting
tdnum = tslice(gbeg:gend)/86400+datenum(1970,1,1);

fw = [0,0,1200,900];
figure('visible','on','position',fw);

% Panel a
subplot(5,1,1)
ph = plot(tdnum,stvec);
datetick('x',15);
axis tight;
xl = get(gca,'xlim'); % Get the original x-limits so that we can use them on all panels
set(gca,'ylim',[-1 1]);
set(gca,'FontSize',8);

% Add time string to title
title(['Chunk ',num2str(periodcount),' ',bm,', starting: ',datestr(tdnum(1)),'. Length = ',num2str(nhrs),' hrs ',num2str(nmin),' minutes'])
%title(datestr(tdnum(1)));

% Panel b
subplot(5,1,2)
ph = plot(tdnum,abs(stvec));
datetick('x',15);
%axis tight;
set(gca,'ylim',[0 1]);
set(gca,'xlim',xl);
set(gca,'FontSize',8);


% Take running mean
stmean = runmean(abs(stvec),10);

% Panel c
subplot(5,1,3)
ph = plot(tdnum,stmean);
datetick('x',15);
%axis tight;
set(gca,'ylim',[0 1]);
set(gca,'xlim',xl);
set(gca,'FontSize',8);
hline = refline([0 0.5]);
hline.Color = 'r';

% Find the outliers and apply
out = find(stmean>0.5);
stmean(out) = nan;

% Panel d
subplot(5,1,4)
ph = plot(tdnum,stmean);
datetick('x',15);
%axis tight;
set(gca,'ylim',[0 1]);
set(gca,'xlim',xl);
set(gca,'FontSize',8);

% Apply outliers to stvec
stvec(out) = nan;

% Panel e
subplot(5,1,5)
ph = plot(tdnum,stvec);
datetick('x',15);
%axis tight;
set(gca,'ylim',[-1 1]);
set(gca,'FontSize',8);
set(gca,'xlim',xl);

% Save to output file
%saveas(ph,[bm,'_chunk_',num2str(periodcount),'_panel'],'png');
