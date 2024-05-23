%% find missing files

dsSS='/Users/clementvulin/Desktop/ilastik/all_final/SecondSeg/';
%dsOri='/Users/clementvulin/Desktop/ilastik/all_final/MasterOriginals210729_bis/';
dsMasked='/Users/clementvulin/Desktop/ilastik/all_final/MaskedImages/';
destMiss='/Users/clementvulin/Desktop/ilastik/all_final/Masked_missing/';

%filesSS=dir([dsSS,'*.tif']);
%filesOri=dir([dsOri,'*.tif']);
filesMasked=dir([dsMasked,'*.tif']);
missList=struct();   

%find missing
tic
for filei=1:numel(filesMasked)
    fName=filesMasked(filei).name(1:end-4);
    if isempty(dir([dsSS,'*',fName,'*']))
        missList(end+1).name=fName; %#ok<SAGROW> 
        source=[filesMasked(filei).folder,'/',filesMasked(filei).name];
        dest=[destMiss,'/',filesMasked(filei).name];
        copyfile (source, dest)
    end
    if filei/200==round(filei/200)
        pct=filei/numel(filesMasked);
        disp(['N=', num2str(filei), '= ',num2str(pct*100, 2) ,'%, t='...
            num2str(toc*(1/pct-1)/60,2),' min remaining'])
    end
end
