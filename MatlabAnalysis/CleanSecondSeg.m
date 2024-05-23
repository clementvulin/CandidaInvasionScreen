function imEnd=CleanSecondSeg(im)
% this is the extracting the colony that will be used from the second
% segmentation

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