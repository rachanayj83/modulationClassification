seqBin = de2bi(seqInt,8,'left-msb')';
seqBin = seqBin(:);

dimSeqBin = numel(seqBin);

nBitExtra = (ceil(dimSeqBin/((N-filDelay)*(C-1)))*((N-filDelay)*(C-1))) - dimSeqBin;   % multiple of sub-carriers
newDimSeqBin = dimSeqBin + nBitExtra;

seqBin = [seqBin; round(rand(newDimSeqBin - dimSeqBin,1))];

seqBinMat = reshape(seqBin,N-filDelay,numel(seqBin)/(N-filDelay));
seqBinMatGuard = [seqBinMat; zeros(filDelay,size(seqBinMat,2))];

seqOutput = seqBinMatGuard(:)';
bitOverHead = numel(seqOutput) - dimSeqBin;