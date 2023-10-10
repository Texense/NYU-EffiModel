% A self contained function to present pixel connectivity
% Input: C_mat: Matrix of pixel connectivity
%        PixId: the pixel I care about
%        PixNum, N_HC, NPixX, NPixY: Basic quantities. Scalars 
%        Crange: range of connnectivity, 1*2 vec
%        AXrange: range of the axis to present, 1* 4 vec
%        CellType: a cell containing the cell type of pre and post. E.G:
%        {'S','C'}

function [] = Paper2FigS1_Conn(C_mat, PixId, PixNum, N_HC, NPixX, NPixY, Crange, AXrange, CellType)
tittext = [CellType{1} ' to ' CellType{2}];

ShowField(C_mat(PixId,:)',...
          1:PixNum,N_HC*NPixX,N_HC*NPixY,...
          tittext, Crange)
colorbar off
axis(AXrange)
axis off

end