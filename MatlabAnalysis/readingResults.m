%% check overlay from Ilastik
% numGap=100; % gap between consecutive pictures
% 
% % where are pics
% %dirO='/Users/clementvulin/Desktop/ilastik/MasterOriginals210729_bis';
% dirO='/Users/clementvulin/Desktop/ilastik/Ilastik_colvsAgar'; %training set
% 
% %dirF='/Users/clementvulin/Desktop/210901_ColVsAgar'; %all
% dirF='/Users/clementvulin/Desktop/210901_ColVsAgar_small'; %training set
% 
% dsO = datastore(dirO,'FileExtensions','.tif','Type','image');
% dsO.Files=natsortfiles(dsO.Files);
% dsF = datastore(dirF,'FileExtensions','.tif','Type','image', 'ReadFcn',@CleanReader);
% dsF.Files=natsortfiles(dsF.Files);

%% for YPD 1 (exp13?)

% check overlay from Ilastik
% numGap=100; % gap between consecutive pictures

% where are pics
dirO='/Users/clementvulin/Desktop/MasterOriginals210927'; % original

dirF='/Users/clementvulin/Desktop/210930_ColVsAgar_YPD1'; % col vs agar set

dsO = datastore(dirO,'FileExtensions','.tif','Type','image');
dsO.Files=natsortfiles(dsO.Files);
dsF = datastore(dirF,'FileExtensions','.tif','Type','image', 'ReadFcn',@CleanReader);
dsF.Files=natsortfiles(dsF.Files);

%% watch
if 0
    f=figure; set(f,'position',[ 38         447        1257         351]);
    for i=1:numGap:numel(dsO.Files)
        
        [imO,fInfo]=readimage(dsO,i); imF=readimage(dsF,i);
        
        subplot(2,2,1)
        imshow(imO); title(extractBetween(fInfo.Filename,'bis/','.tif'))
        subplot(2,2,2)
        imagesc(imF); axis('square')
        subplot(2,2,3);
        imshowpair(imO,imF); title('Mach. Learn.')
        subplot(2,2,4);
        imshowpair(imO,adaptthresh(rgb2gray(imO),0.5)); title('adapt. thresh.')
        waitforbuttonpress
    end
end

%% save masked image

dirSave='/Users/clementvulin/Desktop/210930_MaskedImages/';

tic
tot=numel(dsO.Files);
for i=1:1:tot %5021
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
    
    if round(i/100)==i/100
        done=i/tot;
        left=1-done;
        disp(['done ', num2str(i), '/', num2str(numel(dsO.Files)), '=',num2str(done) ,'% ; remain ', num2str((toc/done*left)/60),' min'])
    end
end

%% second run at masked images


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

