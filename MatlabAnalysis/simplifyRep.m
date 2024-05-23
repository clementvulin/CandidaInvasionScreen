function MeanStd=simplifyRep(whC,LongMat)
    % mean, std and number of rep per condition
    whL=~isnan(LongMat(whC));
    Mean=mean(LongMat(whL,whC));
    Std=std(LongMat(whL,whC));
    MeanStd=[Mean,Std,sum(whL)];
end