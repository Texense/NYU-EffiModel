%% Plot accuracy indication plots
% input: XData, YData, CData: 1*n vecter
%        bounds: [lower bound, upper bound] ; can be flexible for 10%
%        colorcode: [good;ok;bad] what color to indicate
%        varargin: xlabels; ylabels
function [CIndi] = PlotAccu(XData, YData, CData, bounds, colorcode,TitText,varargin)
if ~isempty(varargin)
    xlabels = varargin{1};
    ylabels = varargin{2};
else 
    xlabels = 'S^{IE}/S^{II}';
    ylabels = 'S^{EI}/S^{EE}';
end
FlexPor = 0.15;

Good = CData>=bounds(1)             & CData<=bounds(2);
OK   = (CData>=bounds(1)*(1-FlexPor) & CData<=bounds(2)*(1+FlexPor)) & (~Good);
Bad = ~(Good | OK);

CIndi = zeros(size(CData));
CIndi(Good) = 0; % good
CIndi(OK)   = 1; % ok
CIndi(Bad)  = 2; % ok

ax1 = gca;
hold on
scatter(XData, YData, [], CIndi,'.')
colormap(ax1,colorcode)
caxis([0,2])
xlabel(xlabels);ylabel(ylabels)
axis square
title(TitText)
axis([min(XData) max(XData) min(YData) max(YData)])
end