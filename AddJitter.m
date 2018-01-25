%jitteriser for streamlined

rng('shuffle')  % set the random seed

vars.rhoMin = 0;
vars.rhoMax = 300;

vars.JitterRepeats = 10;

% jitter_FFD = zeros(size(data,1),vars.JitterRepeats);
% jitter_AFD = zeros(size(data,1),vars.JitterRepeats);
    
for r = 1:vars.JitterRepeats

    data_jitter = data;

    for eventID = 1:size(data_jitter,1)
        x1 = data(eventID,vars.xCol);
        y1 = data(eventID,vars.yCol);

        %make a uniformly random angle and distance from origin
        rnd_theta = 2*pi.*rand; % angle from origin
%        rnd_rho = randi([vars.rhoMin vars.rhoMax]); % distance from origin -- old version
        rnd_rho = sqrt(randi([vars.rhoMin.^2 vars.rhoMax.^2])); % distance from origin; squaring used to offset centre-clustering effect.
        [offsetX, offsetY] = pol2cart(rnd_theta,rnd_rho);

        data_jitter(eventID,vars.xCol) = x1 + offsetX;
        data_jitter(eventID,vars.yCol) = y1 + offsetY;
    end
    
        data_len = length(data_jitter(:,vars.xCol));
        data_rand_ind = randperm(data_len);
        data_randx = data_jitter(data_rand_ind,vars.xCol);
        data_randy = data_jitter(data_rand_ind,vars.yCol);
        eventsplot = figure;
        scatter(data_randx,data_randy,1,'.');
        axis(vars.AxisLimits);
        axis image
        title([vars.ExptTitle,' : Randomised Events (Round ',num2str(r),')']);
        SaveFileName = [vars.ExptTitle,' - randomised events (Round ',num2str(r),').png'];
        print(eventsplot,'-dpng','-r300',SaveFileName);
        close(eventsplot)

    
    Message=['Running Parallel Pool Loop: ',num2str(r)];
    disp(Message);
    
    [ppFFD_tmp,ppAFD_tmp,~] = DTF2ParaFunc(data_jitter,vars);
    
    loop_fname = ['Blackwattle_JittCorr_',num2str(r),'.mat'];
    save(loop_fname,'ppFFD_tmp','ppAFD_tmp');
    
    jitter_FFD(:,r) = ppFFD_tmp;
    jitter_AFD(:,r) = ppAFD_tmp;
     
    clear ppFFD_tmp ppAFD_tmp   
end

%          jitter_FFD = horzcat(jitter_FFD,ppFFD_tmp);
%          jitter_AFD = horzcat(jitter_AFD,ppAFD_tmp);


% jitter_FFD = jitter_FriendDists(:,1:2:end);
% jitter_AFD = jitter_FriendDists(:,2:2:end);

jitter_AvgFDists = mean(jitter_FFD,2);
jitter_AvgADists = mean(jitter_AFD,2);

vars.meanFFDCol = size(data,2)+1;
data(:,vars.meanFFDCol) = jitter_AvgFDists;

vars.meanAFDCol = size(data,2)+1;
data(:,vars.meanAFDCol) = jitter_AvgADists;


vars.corrFFDCol = size(data,2)+1;
data(:,vars.corrFFDCol) = 0;
data(:,vars.corrFFDCol) = data(:,vars.FFDCol) ./ data(:,vars.meanFFDCol);

vars.corrAFDCol = size(data,2)+1;
data(:,vars.corrAFDCol) = data(:,vars.AFDCol) ./ data(:,vars.meanAFDCol);

vars.corrinvAFDCol = size(data,2)+1;
data(:,vars.corrinvAFDCol) = 1 ./ data(:,vars.corrAFDCol);

%     data(:,vars.invAFDCol) = normalsumdistance(1) ./ data(:,vars.AFDCol);
vars.corrloginvAFDCol = size(data,2)+1;
data(:,vars.corrloginvAFDCol) = log(data(:,vars.corrinvAFDCol));
    
    
% % ColourmapMax = 5*normaldistance;
% fig_invlogAFD_jittcorr = figure;
% title('Normal Sum AFD');
% scatter(data(:,vars.xCol),data(:,vars.yCol),10,data(:,vars.corrloginvAFDCol),'filled');
% axis square tight
% colorbar
% caxis([0 ColourmapMax]);


%             peaksindex = sindexplot(1,sindexplot(2,:)==max(sindexplot(2,2:end-1))); % ignoring first and last bins
%             ColourmapMax = 5*peaksindex;
data_randz = data(data_rand_ind,vars.corrloginvAFDCol);
% data_randz = data(data_rand_ind,vars.corrAFDCol);
fig_jittcorrloginvAFDCol = figure;
scatter(data_randx,data_randy,3,data_randz,'filled');
axis(vars.AxisLimits);
axis square
%             caxis([0 ColourmapMax]);
colorbar
% ColourMapAxis = caxis;
% caxis([0 ColourMapAxis(2)]);
title([vars.ExptTitle,' : Log of Sum of distances to n(1..10)-(Jitter Corrected)']);
SaveFileName = [vars.ExptTitle,' - fig_invlogAFD_jittcorr(n=10)-cmaprangeexpanded.png'];
print(fig_jittcorrloginvAFDCol,'-dpng','-r600',SaveFileName);

% set a dark bg
set(gca,'Color',[0 0 0.5]);
set(gcf,'InvertHardcopy','off')
SaveFileName = [vars.ExptTitle,' - fig_invlogAFD_jittcorr(n=10)-cmaprangeexpanded_darkBG.png'];
print(fig_jittcorrloginvAFDCol,'-dpng','-r600',SaveFileName);


% a histogram
corrloginvAFDColplot = [];
xbins = min(data(:,vars.corrloginvAFDCol)):0.1:max(data(:,vars.corrloginvAFDCol)); % max(data(:,3)); % bin width = 10 (nm)
[corrloginvAFDColplot(2,:),corrloginvAFDColplot(1,:)] = hist(data(:,vars.corrloginvAFDCol),xbins);
plotcollection2 = figure;
subplot(2,3,6);
bar(corrloginvAFDColplot(1,:),corrloginvAFDColplot(2,:));
axis tight
title(['invlogsumd(',num2str(vars.FurthestFriendID),')-JittCorr(',num2str(vars.JitterRepeats),'x)']);
% may need to adjust this by 'in cell' or 'out of cell' discriminator
normaldistance = corrloginvAFDColplot(1,corrloginvAFDColplot(2,:)==max(corrloginvAFDColplot(2,2:end-1))); % ignoring first and last bins
axis([floor(corrloginvAFDColplot(1,1)) ceil(corrloginvAFDColplot(1,end)) 0 1.1*corrloginvAFDColplot(2,corrloginvAFDColplot(1,:) == normaldistance)]);


% a histogram with a threshold shown
hist_thold = 1;
corrloginvAFDColplot = [];
xbins = min(data(:,vars.corrloginvAFDCol)):0.1:max(data(:,vars.corrloginvAFDCol)); % max(data(:,3)); % bin width = 10 (nm)
[corrloginvAFDColplot(2,:),corrloginvAFDColplot(1,:)] = hist(data(:,vars.corrloginvAFDCol),xbins);
hist_thold_col_low = size(find(corrloginvAFDColplot < hist_thold));
set(0, 'CurrentFigure', plotcollection2);
subplot(2,3,5);
bar(corrloginvAFDColplot(1,1:hist_thold_col_low),corrloginvAFDColplot(2,1:hist_thold_col_low),'b');
hold on
bar(corrloginvAFDColplot(1,hist_thold_col_low+1:end),corrloginvAFDColplot(2,hist_thold_col_low+1:end),'r');
axis tight

title(['invlogsumd(',num2str(vars.FurthestFriendID),')-JittCorr(',num2str(vars.JitterRepeats),'x)']);
% may need to adjust this by 'in cell' or 'out of cell' discriminator
normaldistance = corrloginvAFDColplot(1,corrloginvAFDColplot(2,:)==max(corrloginvAFDColplot(2,2:end-1))); % ignoring first and last bins
axis([floor(corrloginvAFDColplot(1,1)) ceil(corrloginvAFDColplot(1,end)) 0 1.1*corrloginvAFDColplot(2,corrloginvAFDColplot(1,:) == normaldistance)]);

%% Save the plot collection

% send the focus back to the plots
set(0, 'CurrentFigure', plotcollection2);
SaveFileName = [vars.ExptTitle,' - JitCorr histograms.png'];
print(plotcollection2,'-dpng','-r300',SaveFileName);

% 
% 
% 
% % % Stupid data randomiser to reduce print times
% % data_jitter_len = length(data_jitter(:,vars.xCol));
% % data_jitter_rand_ind = randperm(data_jitter_len);
% % data_jitter_randx = data_jitter(data_jitter_rand_ind,vars.xCol);
% % data_jitter_randy = data_jitter(data_jitter_rand_ind,vars.yCol);
% % 
% % 
% % jittereventsplot = figure;
% % scatter(data_jitter_randx,data_jitter_randy,1,'.');
% % axis(vars.AxisLimits);
% % axis image
% % title([vars.ExptTitle,' : Events (jittered)']);
% 
% 
% 
% %Bin the distances to find the most common distance range.
% distanceplot = [];
% xbins = 0:0.1:5; % max(data(:,3)); % bin width = 10 (nm)
% [distanceplot(2,:),distanceplot(1,:)] = hist(data(:,vars.corrAFDCol),xbins);
% plotcollection2 = figure;
% subplot(2,3,1);
% bar(distanceplot(1,:),distanceplot(2,:));
% axis tight
% title(['d(',num2str(vars.FurthestFriendID),')']);
% % may need to adjust this by 'in cell' or 'out of cell' discriminator
% normaldistance = distanceplot(1,distanceplot(2,:)==max(distanceplot(2,2:end-1))); % ignoring first and last bins
% axis([0 distanceplot(1,end) 0 1.1*distanceplot(2,distanceplot(1,:) == normaldistance)]);
% 
% % ColourmapMax = 5*normaldistance;
% % figure;
% % title('d(10)');
% % scatter(data(:,vars.xCol),data(:,vars.yCol),10,data(:,vars.FFDCol),'filled');
% % axis square tight
% % colorbar
% % caxis([0 ColourmapMax]);
% 
% vars.invFFDCol = size(data,2)+1;
% 
% % Invert distances to n10th in data col 5
% % for d = 1:size(data,1)
%     data(:,vars.invFFDCol) = normaldistance ./ data(:,vars.FFDCol);
% % end
% 
% % Distance to the 10th nearest neighbour
% dindexplot = [];
% xbins = 0:0.02:floor(max(data(:,vars.invFFDCol)));
% [dindexplot(2,:),dindexplot(1,:)] = hist(data(:,vars.invFFDCol),xbins);
% subplot(2,3,4);
% bar(dindexplot(1,:),dindexplot(2,:));
% axis tight
% title('peakd10 / d(10)');
% % may need to adjust this by 'in cell' or 'out of cell' discriminator
% peakdindex = dindexplot(1,dindexplot(2,:)==max(dindexplot(2,2:end-1))); % ignoring first and last bins
% axis([0 dindexplot(1,end) 0 1.1*dindexplot(2,dindexplot(1,:) == peakdindex)]);
% 
% ColourmapMax = 3*peakdindex;
% %randomised z data
% data_randz = data(data_rand_ind,vars.invFFDCol);
% fig_invFFD = figure;
% scatter(data_randx,data_randy,1,data_randz,'filled');
% % scatter(data(:,vars.xCol),data(:,vars.yCol),1,data(:,vars.invFFDCol),'filled');
% axis(vars.AxisLimits);
% axis square
% caxis([0 ColourmapMax]);
% colorbar
% title([vars.ExptTitle,' : Distance to n(10)']);
% SaveFileName = [vars.ExptTitle,' - fig_invFFD.png'];
% print(fig_invFFD,'-dpng','-r600',SaveFileName);
% 
% % send the focus back to the plots
% set(0, 'CurrentFigure', plotcollection);
% 
% % Sum of distances to nearest 10 neighbours
% sumdistanceplot = [];
% xbins = 0:5:5000; % bin width = 10 (nm)
% [sumdistanceplot(2,:),sumdistanceplot(1,:)] = hist(data(:,vars.AFDCol),xbins);
% subplot(2,3,2);
% bar(sumdistanceplot(1,:),sumdistanceplot(2,:));
% axis tight
% title('sumd(10)');
% normalsumdistance = sumdistanceplot(1,sumdistanceplot(2,:)==max(sumdistanceplot(2,2:end-1))); % ignoring first and last bins
% axis([0 sumdistanceplot(1,end) 0 1.1*sumdistanceplot(2,sumdistanceplot(1,:) == normalsumdistance(1))]);
% 
% % ColourmapMax = 5*normalsumdistance(1);
% % figure;
% % title('sumd(10)');
% % scatter(data(:,vars.xCol),data(:,vars.yCol),10,data(:,vars.AFDCol),'filled');
% % axis square tight
% % colorbar
% % caxis([0 ColourmapMax]);
% 
% vars.invAFDCol = size(data,2)+1;
% vars.logAFDCol = size(data,2)+2;
% 
% % Invert sum of distances to col 6
% % for d = 1:size(data,1)
%     data(:,vars.invAFDCol) = normalsumdistance(1) ./ data(:,vars.AFDCol);
%     data(:,vars.logAFDCol) = log(data(:,vars.AFDCol));
% % end
% 
% sindexplot = [];
% xbins = 0:0.01:2;
% [sindexplot(2,:),sindexplot(1,:)] = hist(data(:,vars.invAFDCol),xbins);
% subplot(2,3,5);
% bar(sindexplot(1,:),sindexplot(2,:));
% axis tight
% title('peaksumd10 / sumd(10)');
% axis([0 2 0 ceil(1.1*max(sindexplot(2,:)))]);
% 
% 
% % may need to adjust this by 'in cell' or 'out of cell' discriminator
% %randomised z data
% data_randz = data(data_rand_ind,vars.invAFDCol);
% peaksindex = sindexplot(1,sindexplot(2,:)==max(sindexplot(2,2:end-1))); % ignoring first and last bins
% ColourmapMax = 5*peaksindex;
% fig_invAFD = figure;
% scatter(data_randx,data_randy,1,data_randz,'filled');
% % scatter(data(:,vars.xCol),data(:,vars.yCol),1,data(:,vars.invAFDCol),'filled');
% axis(vars.AxisLimits);
% axis square
% caxis([0 ColourmapMax]);
% colorbar
% title([vars.ExptTitle,' : Sum of distances to n(1..10)']);
% SaveFileName = [vars.ExptTitle,' - fig_invAFD.png'];
% print(fig_invAFD,'-dpng','-r300',SaveFileName);
% 
% 
% % Log
% 
% % send the focus back to the plots
% set(0, 'CurrentFigure', plotcollection);
% 
% % Log of sum of distances to nearest 10 neighbours
% logsumdistanceplot = [];
% xbins = floor(min(data(:,vars.logAFDCol))):0.1:ceil(max(data(:,vars.logAFDCol))); % bin width = 10 (nm)
% [logsumdistanceplot(2,:),logsumdistanceplot(1,:)] = hist(data(:,vars.logAFDCol),xbins);
% subplot(2,3,3);
% bar(logsumdistanceplot(1,:),logsumdistanceplot(2,:));
% axis tight
% title('logsumd(10)');
% normallogsumdistance = logsumdistanceplot(1,logsumdistanceplot(2,:)==max(logsumdistanceplot(2,2:end-1))); % ignoring first and last bins
% axis([floor(min(data(:,vars.logAFDCol))) ceil(max(data(:,vars.logAFDCol))) 0 1.1*logsumdistanceplot(2,logsumdistanceplot(1,:) == normallogsumdistance)]);
% 
% vars.invlogAFDCol = size(data,2) + 1;
% 
% % flip the data so low distances == high values
% maxd7 = max(data(:,vars.logAFDCol));
% data(:,vars.invlogAFDCol) = maxd7 - data(:,vars.logAFDCol);
% 
% invlogsumdistanceplot = [];
% xbins = floor(min(data(:,vars.invlogAFDCol))):0.1:ceil(max(data(:,vars.invlogAFDCol))); % bin width = 10 (nm)
% [invlogsumdistanceplot(2,:),invlogsumdistanceplot(1,:)] = hist(data(:,vars.invlogAFDCol),xbins);
% subplot(2,3,6);
% bar(invlogsumdistanceplot(1,:),invlogsumdistanceplot(2,:));
% axis tight
% title('invlogsumd(10)');
% normallogsumdistance2 = invlogsumdistanceplot(1,invlogsumdistanceplot(2,:)==max(invlogsumdistanceplot(2,2:end-1))); % ignoring first and last bins
% axis([floor(min(data(:,vars.invlogAFDCol))) ceil(max(data(:,vars.invlogAFDCol))) 0 1.1*invlogsumdistanceplot(2,invlogsumdistanceplot(1,:) == normallogsumdistance2)]);




data_randz = data(data_rand_ind,vars.corrloginvAFDCol);
fig_jittcorrloginvAFDCol_thold = figure;
scatter(data_randx,data_randy,1,data_randz,'filled');
axis(vars.AxisLimits);
axis square
colorbar
% ColourMapAxis = caxis;
% caxis([0 ColourMapAxis(2)]);
title([vars.ExptTitle,' : Log of Sum of distances to n(1..10)-(Jitter Corrected)']);

max_colors = 256;
cbar = zeros(max_colors,3);
cbarmin = min(data(:,vars.corrloginvAFDCol));
cbarmax = max(data(:,vars.corrloginvAFDCol));

% set the axis to desired min/max
caxis([cbarmin cbarmax]);

% threshold
cbarthr1 = 1.0;

cbar_unit = (cbarmax-cbarmin)/max_colors;
cbars_under = 1+floor((cbarthr1 - cbarmin) / cbar_unit);
cbars_over = max_colors - cbars_under;

red = [1 0 0];
green = [0 1 0];
blue = [0 0 1];
yellow = [1 1 0];
magenta = [1 0 1];
cyan = [0 1 1];
white = [1 1 1];
darkblue = [0 0 0.5];

cbar(1:cbars_under,:) = repmat(blue,[cbars_under 1]);
cbar(cbars_under+1:end,:) = repmat(red,[cbars_over 1]);
% cbar(cbars_under+1:cbars_thr2,:) = repmat(magenta,[cbars2 1]);
% cbar(cbars_thr2+1:cbars_thr3,:) = repmat(yellow,[cbars3 1]);
% cbar(cbars_thr3+1:end,:) = repmat(cyan,[cbars_over 1]);

colormap(cbar)
set(gca,'Color','k');

set(gcf,'InvertHardcopy','off')
SaveFileName = [vars.ExptTitle,' - fig_invlogAFD_jittcorr(n=10)-threshold.png'];
print(fig_jittcorrloginvAFDCol_thold,'-dpng','-r600',SaveFileName);

% again with BG to match lower colour index
caxis([(mean(data(:,vars.invlogAFDCol))-(2 * std(data(:,vars.invlogAFDCol)))) (mean(data(:,vars.invlogAFDCol))+(2 * std(data(:,vars.invlogAFDCol))))]);
c = colormap;
set(gca,'Color',c(1,:));
colorbar
set(gcf,'InvertHardcopy','off')
title([vars.ExptTitle,' : Log of Sum of distances to n(1..10)']);
SaveFileName = [vars.ExptTitle,' - fig_invlogAFD_dark_bg.png'];
print(fig_invlogAFD_dark_BG,'-dpng','-r300',SaveFileName);

