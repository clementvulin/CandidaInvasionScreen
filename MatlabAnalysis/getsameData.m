function numlines=getsameData(numline,whC, LongMat)
% specifically intented to use with the AllDataLong matrix for colonies
% using the columns specified as wh

StartLine=LongMat(numline,:);
wh=ones(size(LongMat,1),1);
for columni=whC
    wh=wh & (LongMat(:,columni)==StartLine(columni));
end
    numlines=find(wh);
end