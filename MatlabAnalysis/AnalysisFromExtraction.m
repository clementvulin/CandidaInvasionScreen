%% reran extraction since plates are ill-numbered

%% functions
addpath('/Users/clementvulin/Documents/p12_software/MATLAB/ImageAnalysisCV');


%% init
% This thime only extracting the tif to Originals masterfile directly
%dirM_all='/Volumes/CLEM_MIC6/C_albicans/Julian_Sutter';
dirM_all='/Users/clementvulin/Documents/ZinkernagelLab/z09_Candida_new';

dirDest='/Users/clementvulin/Desktop/MasterOriginals210927/';

list=natsortfiles(dir([dirM_all,'/exp*',]));

regYN=1;downsampleF=0.03;% registration is based on downsampling of images:
showReg=0;

for expt=1
    disp([])
    disp(['Expt', num2str(expt)]);
    %% prepare
    % set to directory
    dirM=[dirM_all,'/',list(expt).name,'/'];
    
    % list sub directories
    ListAllDays=natsortfiles(dir([dirM,'*pics_D*'])); %list of directories
    
    % load colony D3 data from ColTApp
    load([dirM,ListAllDays(1).name,'/',ListAllDays(1).name, '_all.mat']);
    
    % list of files in subfolders
    for ti=1:numel(ListAllDays) %   for all times
            imList.(['t',num2str(ti)])=natsortfiles(dir([dirM,ListAllDays(ti).name,'/*.JPG']));
        if isempty(imList.(['t',num2str(ti)]))
            imList.(['t',num2str(ti)])=natsortfiles(dir([dirM,ListAllDays(ti).name,'/*.jpg']));
        end
    end
    
    %% registration
    if regYN==0
        
        RegVal=struct();
        fprintf(['Registration: 0/',num2str(size(p.counts,1)-1)])
        
        for pics=2:size(p.counts,1) % for all pics in D3 exept first one
            
            %   get D3 pic
            imD3=imread([dirM,ListAllDays(1).name,'/',p.l(pics).name]);
            
            for ti=1:numel(ListAllDays) %   for all times (this implies registration on oneself)
                
                % get new pic
                imDi=imread([dirM,ListAllDays(ti).name,'/',imList.(['t',num2str(ti)])(pics).name]);
                
                %  compute and store the registration
                [OutIm,Tform]=SmallReg(imD3,imDi,downsampleF);
                RegVal.(['p',num2str(pics)]).(['t',num2str(ti)])=Tform;
                
                if showReg
                    imshowpair(OutIm,imD3)
                    nmD=extractBetween(ListAllDays(ti).name,'_','_'); nmD=nmD{1};
                    title([nmD,'p',num2str(pics),'t',num2str(ti)])
                    waitforbuttonpress
                end
                
            end %through time
            
            % user output on command line
            if (pics-1)<10 %first, remove the old charaters
                fprintf(repmat(char(8), 1, 4))
            else
                fprintf(repmat(char(8), 1, 5))
            end
            
            fprintf ([num2str(pics-1),'/',num2str(size(p.counts,1)-1)]);
            
        end %for all D3 pics
        
        % save all in D3 folder
        save([dirM,ListAllDays(1).name,'/RegVal210729.mat'],'RegVal');
    else
        % load the registration values
        load([dirM,ListAllDays(1).name,'/RegVal210729.mat'],'RegVal');
    end %if registration
    
    %% extract the images
    
    % 1) decide size of picture based on Scale factor of image
    SizeWell=34000; %size of wells in micrometer
    if ~isnan(p.umConversion(1))
        SizeImPxl=round(SizeWell/p.umConversion(1));
    else
        disp('um conversion missing')
        disp([dirM,ListAllDays(1).name,'/',ListAllDays(1).name, '_all.mat'])
        return
    end
    
    %2)
    % for each pic in D3
    tic;fprintf([' Extraction of colonies on fr: 0/',num2str(size(p.counts,1)-1)])
    for pics=2:size(p.counts,1) % for all pics in D3
        for ti=1:numel(ListAllDays) %   for all times
            
            % load the picture of 6 cols
            imDi=imread([dirM,ListAllDays(ti).name,'/',imList.(['t',num2str(ti)])(pics).name]);
            nmDay=extractBetween(ListAllDays(ti).name,'_','_'); nmDay=nmDay{1}; %create a day stamp for saving
            
            % apply registration
            if ti==1
                szfix=imref2d(size(imDi)); %create the reference frame for reg on D3 (ti==1)
            else
                tform=RegVal.(['p',num2str(pics)]).(['t',num2str(ti)]);
                imDi=imwarp(imDi, tform,'OutputView',szfix); %registration
            end
            
            colCenters=p.counts{pics,1}; %get col centers
            
            if ~isempty(colCenters)
                % sort the centers to get col number
                [c,I]=sortrows(ceil(colCenters/100),[2,1]); %use rounding to force groups
                colcentersS=colCenters(I,:); %this is the sorter colCenter
                
                % over all colonies
                for col=1:numel(p.counts{pics,2})% for each col
                    % cut to the right size
                    y=round(colcentersS(col,1));x=round(colcentersS(col,2));S=round(SizeImPxl/2*0.7);
                    SmallImage=imDi(x-S:x+S,y-S:y+S,:);
                    % save cut picture as OpXcY (Original plate 1:45 col 1-6)
                    if ~exist(dirDest,'dir'); mkdir(dirDest); end %directory
                    
                    namePic=['p',num2str(pics-1),'c',num2str(col),'.tif'];
                    imwrite(SmallImage,[dirDest,'Exp',num2str(expt),nmDay,namePic]);
                    
                end
            end
        end
        % user output on command line
        if (pics-1)<10 %first, remove the old caraters
            fprintf(repmat(char(8), 1, 4))
        else
            fprintf(repmat(char(8), 1, 5))
        end
        fprintf ([num2str(pics-1),'/',num2str(size(p.counts,1)-1)]);
    end
end


