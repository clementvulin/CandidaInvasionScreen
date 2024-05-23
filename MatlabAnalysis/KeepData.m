function ToKeep=KeepData(DataLong, col, val)
ToKeep=ones(size(DataLong,1),1);    
for coli=1:numel(col)
    ToKeep= ToKeep & DataLong(:,col(coli))==val(coli);
end
end