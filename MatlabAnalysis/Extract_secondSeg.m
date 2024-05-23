%% extracting filament data from the 2 steps segmentation.

dirMasked='/Volumes/CLEM_MIC6/C_albicans_LastBckUp200730/C_10_newSeg/MaskedImages/'; %all set
%dirMasked='/Users/clementvulin/Desktop/Ilastik_segondsegmentation/'; %training set
dsM = datastore(dirMasked,'FileExtensions','.tif','Type','image');
dsM.Files=natsortfiles(dsF.Files);


%% find all the empty files, modify them
%figure; 
numEmpty=[];thresharea=2000;
for i=1:numel(dsM.Files)
    [im,fInfo]=readimage(dsM,i);
%     imshow(im)
%     waitforbuttonpress;
    if sum((im(:)>0))<thresharea*3
        numEmpty=[numEmpty,i];
    end
end


%% for all found, rerun
% in the end there was p=231/14845=1.5% failed

% where are pics
dirO='/Users/clementvulin/Desktop/ilastik/MasterOriginals210729_bis';
% dirO='/Users/clementvulin/Desktop/ilastik/Ilastik_colvsAgar'; %training set

dirF='/Users/clementvulin/Desktop/210901_ColVsAgar'; %all
% dirF='/Users/clementvulin/Desktop/210901_ColVsAgar_small'; %training set

dsO = datastore(dirO,'FileExtensions','.tif','Type','image');
dsO.Files=natsortfiles(dsO.Files);
dsF = datastore(dirF,'FileExtensions','.tif','Type','image', 'ReadFcn',@CleanReader);
dsF.Files=natsortfiles(dsF.Files);

dirSave='/Volumes/CLEM_MIC6/C_albicans_LastBckUp200730/C_10_newSeg/MaskedImages/';
%dirSave='/Users/clementvulin/Desktop/Ilastik_segondsegmentation/';

tic
tot=numel(dsO.Files);
for i=numEmpty %5021
    [imO,fInfo]=readimage(dsO,i); [imF,fInfoF]=readimage(dsF,i);
    %fname=extractAfter(fInfo.Filename,'bis/'); %for all
    fname=extractAfter(fInfo.Filename,[dirO,'/']); %valid for any
    
    imNew=[];
    if ~isempty(imF) %is empty if there was no colony found
        for colori=1:3
            imNewC=imO(:,:,colori); imNewC(~imF)=nan;
            imNew(:,:,colori)=imNewC;
        end
    else %no colony found
        imNew=imO*0;
    end
    imNew=uint8(imNew);
    imwrite(imNew,[dirSave,fname])
    
    i
end
for i=1:10; beep; pause(0.5); end

%% some functions
function imEnd=CleanReader(imageFileName)
% read
im=imread(imageFileName);
if range(im(:))==0
    imEnd=[]; 
    return;
end

% binarize
bw=imbinarize(im);bw=~bw; bw=imfill(bw,'holes');

% get centroids
labeledImage = bwlabel(bw);
T=regionprops(bw,'Centroid','PixelIdxList');
xy = vertcat(T.Centroid);
x = xy(:, 1);
y = xy(:, 2);

% find the closest to center
distances = sqrt((size(im,1)/2 - x) .^ 2 + (size(im,2)/2 - y) .^ 2);
[~, I] = min(distances);

% keep the closest
imCenter=ismember(labeledImage,I);
imDil=imdilate(imCenter,strel('disk',5,8));
imClos=imclose(imDil,strel('disk',5,8));
imEnd=imfill(imClos,'holes');

end