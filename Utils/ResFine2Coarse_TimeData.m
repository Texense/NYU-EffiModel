%% Temporal vector, change resolution:
% Res: New Resolution
% WinSize: Old resolution
% SNWData: The original data

function SNWDataRes = ResFine2Coarse_TimeData(Res, WinSize,  SNWData)
 %change bin resolution to 4ms
NewBinSz = floor(Res/WinSize);
Len = length(SNWData.SpE);
SNWDataRes.SpE = mean(reshape(SNWData.SpE,NewBinSz,floor(Len/NewBinSz)));
SNWDataRes.SpI = mean(reshape(SNWData.SpI,NewBinSz,floor(Len/NewBinSz)));
SNWDataRes.VE = mean(reshape(SNWData.VE,NewBinSz,floor(Len/NewBinSz)));
SNWDataRes.VI = mean(reshape(SNWData.VI,NewBinSz,floor(Len/NewBinSz)));