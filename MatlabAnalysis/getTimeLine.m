function numlines=getTimeLine(numline,LongMat)
% specifically intented to use with the AllDataLong matrix for colonies
% using thos columns: 
% 14) Exp
% 15) Day 
% 16) plate
% 17) colony

StartLine=LongMat(numline,:);
SameExp=LongMat(:,14)==StartLine(14);
Sameplate=LongMat(:,16)==StartLine(16);
Samecolony=LongMat(:,17)==StartLine(17);

numlines=find(SameExp&Samecolony&Sameplate);

end