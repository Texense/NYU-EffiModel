%% From neuron vector to spatial map of pixels
% Input: FrE: Quantity vector for neurons, N*1
%        NnE: Spatial mapping of neurons, struct with X and Y coordinate
%        info
%        NPxX/Y: Total Number of pixel on X dir and Y dir
% 
% Output: SpaPixMat: Pixel Mat, averaging neurons in each pixel
%         NnEPixel:  mapping neuindex to pixel index, struct wit X and Y
%         coordinate
% Zhuo-Cheng Xiao 03/31/2021
function [SpaPixMat,NnEPixel] = NeuVec2Pixel(FrE,NnE,NPxX,NPxY)
% first, transfrom vec into the spatial Map
% NeuVecMap = zeros(max(NnE.Y),max(NnE.X));
% for NeuInd = 1:length(FrE)
%     NeuVecMap(NnE.Y(NeuInd),NnE.X(NeuInd)) = FrE(NeuInd);   
% end

% Now, average Spatial Map into pixels
NnEPixel.X = ceil(NnE.X/(max(NnE.X)/NPxX));
NnEPixel.Y = ceil(NnE.Y/(max(NnE.Y)/NPxY));
NnEPixel.Vec = NnEPixel.Y + (NnEPixel.X-1)*NPxY;
SpaPixMat = zeros(NPxY,NPxX);
for PixelY = 1:NPxY
    for PixelX = 1:NPxX
        SpaPixMat(PixelY,PixelX) = nanmean(FrE(NnEPixel.X == PixelX & ...
                                               NnEPixel.Y == PixelY));
    end
end
end
