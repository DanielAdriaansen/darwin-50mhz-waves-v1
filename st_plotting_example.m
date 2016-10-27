fw = [0,0,900,700];
figure('visible','on','position',fw);
%imagesc(stt,stf,(abs(str).*abs(str)));
%imagesc(stt,stf,(abs(str)));
imagesc(stt,stf,totpow);
set(gca,'YDir','normal');
cbh = colorbar;
ylabel('Frequency 1/s');
xlabel('Time (minutes)');
title({[bm,' period ',num2str(periodcount)],['Begin = ',datestr(tslice(gbeg)/86400+datenum(1970,1,1))],['End = ',datestr(tslice(gend)/86400+datenum(1970,1,1))],[num2str(nhrs),' HRS ',num2str(nmin),' MIN']})
%ylabel(cbh,'abs(str)^2 (m^2/s^2)');
%ylabel(cbh,'abs(str) (m/s)');
ylabel(cbh,'totpow (m/s)');
%saveas(gcf,['st_',bm,'_period_',num2str(periodcount)]);

% Axis labels for totpow
set(gca,'YTick',[0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5])
set(gca,'YTickLabel',[0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5])