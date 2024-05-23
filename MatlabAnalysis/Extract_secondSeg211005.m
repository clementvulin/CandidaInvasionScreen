%% Two stage segmentation using ilastik

%% quick analysis on training files

dirOri='/Users/clementvulin/Desktop/ilastik/Ilastik_colvsAgar_training/01_OriginalTraining';
dirM='/Users/clementvulin/Desktop/ilastik/Ilastik_colvsAgar_training/03_MaskedTraining';
dir2S='/Users/clementvulin/Desktop/ilastik/Ilastik_colvsAgar_training/04_SecondSegTraining';

% using masked files instead of calculated masks, because manually corrected
dsMasked = datastore(dirM,'FileExtensions','.tif','Type','image','ReadFcn',@GetMask);
dsMasked.Files=natsortfiles(dsMasked.Files); 

ds2ndSeg = datastore(dir2S,'FileExtensions','.tif','Type','image','ReadFcn',@SimpleRead);
ds2ndSeg.Files=natsortfiles(ds2ndSeg.Files);

dsOri= datastore(dirOri,'FileExtensions','.tif','Type','image');
dsOri.Files=natsortfiles(dsOri.Files);

%% 
figure; set(gcf,'position',[143 113 1239 685])
ExtractedData=nan(numel(ds2ndSeg.Files),16);
for imi=1:1:numel(ds2ndSeg.Files)
    %[imM,fInfo]=readimage(dsMasked,imi); 
    [imO,fInfo]=readimage(dsOri,imi); 
    
    imMasked=readimage(dsMasked,imi);
    
    im2nd=readimage(ds2ndSeg,imi);
    % proportion of void in mask
    imMaskedSmall=imerode(imMasked,strel('disk',10,8)); %reverse operation of extraction (imdilate)
    VoidMask=imMaskedSmall .* im2nd==2;
    propVoid=sum(VoidMask(:))/sum(imMaskedSmall(:));
    % if proportion is too high, replace with filaments
    if propVoid>0.3
        im2nd(VoidMask)=3;
    end

    % imshowpair(imO,im2nd,'montage')
    %im2nd(1,1:4)=1:4;
    subplot(2,1,1); imshowpair(imO, imMasked,'montage')
    subplot(2,1,2); imshowpair(im2nd,imMasked,'montage')
    
    title(['propVoid ', num2str(propVoid*100,2),'%'])
    %subplot(1,3,3); imshowpair(im2nd,im2nd)

    ExtractedData(imi,1:12)=extract_data(im2nd);
    ExtractedData(imi,1:13)=propVoid;
    ExtractedData(imi,14:17)=extract_NameFile(fInfo);
    
    waitforbuttonpress
end

clear imMaskedSmall imVoid fInfo imMasked im2nd imi imO propVoid VoidMask

%% Batch processing
dirOri='/Users/clementvulin/Desktop/ilastik/all_final/1_MasterOriginals210729_bis';
dirM='/Users/clementvulin/Desktop/ilastik/all_final/3_MaskedImages';
dir2S='/Users/clementvulin/Desktop/ilastik/all_final/4_SecondSeg';

% using masked files instead of calculated masks, because manually corrected
dsMasked = datastore(dirM,'FileExtensions','.tif','Type','image','ReadFcn',@GetMask);
dsMasked.Files=natsortfiles(dsMasked.Files); 

ds2ndSeg = datastore(dir2S,'FileExtensions','.tif','Type','image','ReadFcn',@SimpleRead);
ds2ndSeg.Files=natsortfiles(ds2ndSeg.Files);
%%
ExtractedData=batch_process(ds2ndSeg, dsMasked);
%save(['/Users/clementvulin/Desktop/ilastik/all_final/',datestr(now,1),'SecondSegAnalysis.mat'])

%% functions
function imEnd=SimpleRead(imageFileName)
   im=imread(imageFileName);
   % 1 is colony center
   % 2 is void
   % 3 is filament
   % 4 is star
   
   if sum(im,'all')==0
        imEnd=[]; 
        return;
   end

   % fill the center of the colony
   colcenter=im; 
   colcenter(im~=1)=0; %only keeping the colony
   colcenter2=getcenter(colcenter); % function to get center
   if isempty(colcenter2)
       imEnd=[]; 
       return;
   end
   im(im==1)=3; %replacing all colony by filament
   im(colcenter2==1)=1; %and put the colony back

   % in case the whole image is the colony
   if sum(im==1)>0.9*size(im,1)*size(im,2)
       im=im*0+2; %all is void
   end
    
   imEnd=im;
end

function imEnd=getcenter(im)
colcenter2 = imfill(im,'holes');
   colcenter21=imbinarize(colcenter2);
   colcenter22=imopen(colcenter21,strel('square',2));
    
   % only keep the center
   % get centroids
   labeledImage = bwlabel(colcenter22);

   % first, threshold
   T=regionprops(colcenter22,'Centroid','PixelIdxList','Area');
   %T=T([T.Area]>4000);
   ThrsIm=ismember(labeledImage,find([T.Area]>4000));


   % in case it doesn't work, look for the closest to center
      labeledImage = bwlabel(ThrsIm);
         T=regionprops(ThrsIm,'Centroid','PixelIdxList','Area');
   xy = vertcat(T.Centroid);
   if isempty(xy)
       disp('no center');
       imEnd=[];
       return
   end
   x = xy(:, 1);
   y = xy(:, 2);

   % find the closest to center
   distances = sqrt((size(im,1)/2 - x) .^ 2 + (size(im,2)/2 - y) .^ 2);
   [~, I] = min(distances);

   % keep the closest
   imEnd=ismember(labeledImage,I);
%    imDil=imdilate(imCenter,strel('disk',5,8));
%     imClos=imclose(imDil,strel('disk',5,8));
%     imEnd=imfill(imClos,'holes');
end

function imEnd=GetMask(imageFileName)
    im=imread(imageFileName);
    im2=im(:,:,1); 
    im3=im2; 
    im3(im2~=0)=1;
    imEnd=imbinarize(im3);
end

function values=extract_data(im)
% this function get variables from a segmented images:
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
values=nan(1,9);

% Values on image:
% 1 is colony center
imCenter=im; imCenter(im~=1)=0;
ColonyCenterXY=regionprops(imCenter,'Centroid');ColonyCenterXY=ColonyCenterXY(1).Centroid;
% 2 is void
% 3 is filament
imFil=im; imFil(im~=3)=0; imFil=imbinarize(imFil);
% 4 is star
imStar=im; imStar(im~=4)=0; imStar=imbinarize(imStar);

%% 1) radius of colony center (area based)
values(1)=sqrt(sum(imCenter,'all')/pi);

%% 2) filament area
values(2)=sum(imFil,'all');

%% 3) filament max radius
[a,b]= find(imFil);
Distances=sqrt((ColonyCenterXY(2) - a) .^ 2 + (ColonyCenterXY(1) - b) .^ 2);
values(3)=prctile(Distances,99);
%  figure; imshow(imFil); hold on; scatter(ColonyCenterXY(1),ColonyCenterXY(2));
%  viscircles(ColonyCenterXY,values(2));

%% 4) number of filament sectors
% get pixel values at R=30 pctile of Rmax
Rsmall=prctile(Distances,30);
% create a circular mask line
xc = ColonyCenterXY(1);
yc = ColonyCenterXY(2);
[xx,yy] = meshgrid(1:size(im,1),1:1:size(im,2));
LineSize=2; DilSize=8;
mask = hypot(xx - xc, yy - yc) <= (Rsmall+LineSize) & hypot(xx - xc, yy - yc) >= (Rsmall-LineSize) ;
ImSect=imdilate(imFil & mask, strel('disk',DilSize,8));
RegSect=regionprops(ImSect,'Area');
values(4)=sum([RegSect.Area]/(4*Rsmall*4)>0.5); % this is an arbitrary way to count sectors
% figure; imshow(imFil);
% figure; imshow(ImSect)
% figure; hist([test.Area]/(4*Rsmall*4),100)


%% 5) star area
values(5)=sum(imStar,'all');

%% 6) star max radius
[a,b]= find(imStar);
Distances=sqrt((ColonyCenterXY(2) - a) .^ 2 + (ColonyCenterXY(1) - b) .^ 2);
values(6)=prctile(Distances,99);

%% 7) number of star sectors
Rsmall=prctile(Distances,30);
mask = hypot(xx - xc, yy - yc) <= (Rsmall+LineSize) & hypot(xx - xc, yy - yc) >= (Rsmall-LineSize) ;
ImSect=imdilate(imStar & mask, strel('disk',DilSize,8));
RegSect=regionprops(ImSect,'Area');
values(7)=sum([RegSect.Area]/(4*Rsmall*4)>0.5);

%% 8) total colony area
values(8)=sum(im==1 | im==3 | im==4,"all");

%% 9) max colony radius
values(9)=max(values([1,3,6]));
%% 10) colony center circularity
Circ=regionprops(imCenter,'Circularity');
values(10)=Circ(1).Circularity;

%% 11-12) colony center
values(11:12)=ColonyCenterXY;

end

function values=extract_NameFile(fInfo)
    % extracting exp number et al. from file name
    % 1) Exp
    % 2) Day 
    % 3) plate
    % 4) colony
    PathF=fInfo.Filename;
    Seps=strfind(PathF,'/');
    Name=PathF(Seps(end)+1:end-4);
    values(1)=str2double(extractBetween(Name,'Exp','D'));
    values(2)=str2double(extractBetween(Name,'D','p'));
    values(3)=str2double(extractBetween(Name(4:end),'p','c'));
    values(4)=str2double(extractBetween(Name(4:end),'c','_sec'));
end

function ExtractedData=batch_process(ds2ndSeg, dsMasked)
numFiles=numel(ds2ndSeg.Files);
numGapPctFollow=100;
ExtractedData=nan(numFiles,16);
tic
try
    for imi=1:1:numFiles
        imMasked=readimage(dsMasked,imi);
        if sum(imMasked,"all")==0
            continue;
        end

        [im2nd,fInfo]=readimage(ds2ndSeg,imi);
        if isempty(im2nd)
            continue;
        end

        %% proportion of void in mask
        imMaskedSmall=imerode(imMasked,strel('disk',10,8)); %reverse operation of extraction (imdilate)
        VoidMask=imMaskedSmall .* im2nd==2;
        propVoid=sum(VoidMask(:))/sum(imMaskedSmall(:));
        % if proportion is too high, replace with filaments
        if propVoid>0.3
            im2nd(VoidMask)=3;
        end

        %% extract data from image
        ExtractedData(imi,1:12)=extract_data(im2nd);
        ExtractedData(imi,13)=propVoid;
        ExtractedData(imi,14:17)=extract_NameFile(fInfo);

        if imi/numGapPctFollow==round(imi/numGapPctFollow)
            pct=imi/numFiles;
            disp(['N=', num2str(imi), ' => ',num2str(pct*100, 2) ,'%, t='...
                num2str(toc*(1/pct-1)/60,2),' min remaining'])
        end
    end
catch
    disp(['error in file n°' num2str(imi)])
    return
end

end

