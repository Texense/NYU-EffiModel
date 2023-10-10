% For a pixel vector, show how the map is

function [] = ShowField(Vec, Range, XPix, YPix, varargin)
if size(Vec,2) == 1
    Map = reshape(Vec(Range),YPix,XPix);
    imagesc(Map)
    colormap jet
    set(gca,'YDir','normal')
    colorbar
elseif size(Vec,2) == 3
    Map = zeros(YPix,XPix,3);
    for LgnInd = 1:3
        Map(:,:,LgnInd) = reshape(Vec(Range,LgnInd),YPix,XPix);
    end
    image(Map)
    set(gca,'YDir','normal')
end
%% Plot HC boundaries
if nargin > 6
    PixNum = varargin{3};
else
    PixNum = 10;
end

hold on
XHc = floor(XPix/PixNum); YHc = floor(YPix/PixNum);

grayColor = [.7 .7 .7];
if XHc>1
for VerInd = 1:XHc -1
    plot(ones(length(0:YPix))*(VerInd*PixNum+0.5),0:YPix,'Color',grayColor,'LineWidth',1.5)
end
end

if YHc>1
for HorInd = 1:YHc -1
    plot(0:XPix,ones(length(0:XPix))*(HorInd*PixNum+0.5),'Color',grayColor,'LineWidth',1.5)
end
end

if nargin > 4
    title(varargin{1})
end
if nargin > 5 && ~isempty(varargin{2})
    caxis(varargin{2})
end
axis square
end