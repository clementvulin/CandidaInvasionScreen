%% save table graphs

% 1) radius of colony center (area based)
% 2) filament area
% 3) filament max radius
% 4) number of filament sectors
% 5) star area
% 6) star max radius
% 7) number of star sectors
% 8) total colony area
% 9) max colony radius
% 10) colony center circularity
% 11) colony center X
% 12) colony center Y
% 13) prop. of colony on first round that ended up empty in second
% 14) Exp
% 15) Day
% 16) plate
% 17) colony
% 18) T°C
% 19) [Glu]
% 20) [YP]
% 21) pH (0 is no Buffer)
% 22) strain
% 23) inoc from liq (as of exp 15, all exept exp 6)
% 24) GR D3-7
% 25) first day of filamentation

%% 
%addpath('/Users/clementvulin/Documents/p12_software/MATLAB/tight_subplot.m')
cmap=flipud(colormap('parula')); cmap(1,:)=[0 0 0];
pHs=[5 0 7.4 8];TCs=[30 37 40];
YP=1; glus=[0.1 0.25 0.5 2];
figure;set(gcf,'Position',[440     1   647   797]);
[ha, ~] = tight_subplot(10, 4, [.01 .03],[0.05 .02],[.01 .01]);
smallMatAll=nan(numel(TCs),numel(pHs),40);
for stri=1:10
    for glui=1:4
        smallMat=nan(numel(TCs),numel(pHs));
        for pHi=1:4
            for TCi=1:3
                % 14) Exp 15) Day 16) plate 17) colony
                % 18) T°C
                % 19) [Glu]
                % 20) [YP]
                % 21) pH (0 is no Buffer)
                % 22) strain
                % 24) GR D3-7
                % 25) first day of filamentation
                ToKeep=KeepData(AllDataLong, [15, 18:22], [3, TCs(TCi), glus(glui), YP, pHs(pHi),stri]);
                smallMat(TCi,pHi)=nanmean(AllDataLong(ToKeep,25)); %#ok<NANMEAN> 
            end
        end
        ii=glui+(stri-1)*4;
        axes(ha(ii)) %#ok<LAXES> 
        smallMat(smallMat<1)=1;
        smallMatAll(:,:,end+1)=smallMat; %#ok<SAGROW> 
        imagesc(smallMat,[2 15]);colormap(cmap);
        
        axis("square"); set(ha(ii),'XTickLabel','','YTickLabel','')
        if ii<5
            title(['[Glu]=' num2str(glus(glui))])
        end
        if glui==1
            set(gca,'YTickLabel',{'30°C','37°C','40°C'})
            ylabel(['str.' num2str(stri)])
        end

        if stri==10
            set(gca,'XTick',1:1:4,'XTickLabel',{'5','none','7.4','8'})
            xlabel('buffer pH')
        end
    end
end
c=colormap;colorbar
% waitforbuttonpress
%clear pHs TCs stri YP smallMat p TCi glus glui pHi ToKeep

%% plot a mean
figure;
imagesc(nanmean(smallMatAll,3),[0 100]) %#ok<NANMEAN> 
set(gca,'YTickLabel',{'30°C','37°C','40°C'},'YTick',1:1:3)
            set(gca,'XTick',1:1:4,'XTickLabel',{'5','7.4','none','8'})
set(gca, 'LineWidth',2,'FontSize',15, 'Box', 'off','TickDir','out');
colorbar
title('mean filamentation ring [pxl], all strains')

%% 
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%

%% 
% cmap=flipud(colormap('parula')); cmap(1,:)=[0 0 0];
pHs=[5 0 7.4 8];TCs=[30 37 40];
YP=1; glus=[0.1 0.25 0.5 2];
figure;set(gcf,'Position',[440     1   647   797]);
[ha, pos] = tight_subplot(10, 4, [.01 .03],[0.05 .02],[.01 .01]);
smallMatAll=nan(numel(TCs),numel(pHs),40);
for stri=1:10
    for pHi=1:4
        smallMat=nan(numel(TCs),numel(pHs));
        for glui=1:4
            for TCi=1:3
                % 14) Exp 15) Day 16) plate 17) colony
                % 18) T°C
                % 19) [Glu]
                % 20) [YP]
                % 21) pH (0 is no Buffer)
                % 22) strain
                % 24) GR D3-7
                % 25) first day of filamentation
                ToKeep=KeepData(AllDataLong, [15, 18:22], [3, TCs(TCi), glus(glui), YP, pHs(pHi),stri]);
                smallMat(TCi,glui)=nanmean(AllDataLong(ToKeep,25)); %#ok<NANMEAN> 
            end
        end
        ii=pHi+(stri-1)*4;
        axes(ha(ii)) %#ok<LAXES> 
        smallMat(smallMat<1)=1;
        smallMatAll(:,:,end+1)=smallMat; %#ok<SAGROW> 
        imagesc(smallMat,[2 15]);colormap(cmap);
        
        axis("square"); set(ha(ii),'XTickLabel','','YTickLabel','')
        if ii<5
            title(['pH=' num2str(pHs(pHi))])
        end
        if pHi==1
            set(gca,'YTickLabel',{'30°C','37°C','40°C'})
            ylabel(['str.' num2str(stri)])
        end

        if stri==10
            set(gca,'XTick',1:1:4,'XTickLabel',{'0.1','0.25','0.5','2'})
            xlabel('[glu]')
        end
    end
end
c=colormap;colorbar
% waitforbuttonpress
%clear pHs TCs stri YP smallMat p TCi glus glui pHi ToKeep