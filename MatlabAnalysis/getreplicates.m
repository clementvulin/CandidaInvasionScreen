function numlines=getreplicates(numline,LongMat)
% specifically intented to use with the AllDataLong matrix for colonies
% using those columns: 
% 18) T°C 
% 19) [Glu]
% 20) [YP]
% 21) pH (0 is no Buffer)
% 22) strain
% 23) inoc from liq (as of exp 15, all exept exp 6)

StartLine=LongMat(numline,:);
wh=ones(size(LongMat,1),1);
for columni=18:23
    wh=wh & (LongMat(:,columni)==StartLine(columni));
end
    numlines=find(wh);
end