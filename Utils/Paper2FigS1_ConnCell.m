% A self contained function to present cell connectivity
% ScatAdj between S and other celltype

function [] = Paper2FigS1_ConnCell(NnS, PreSInd, SX, SY, ScatAdj,...
                                   N_HC, n_S_HC, CellType,Color)
hold on
S = scatter(NnS.X(PreSInd),NnS.Y(PreSInd), '.');
S.MarkerEdgeColor = Color;
scatter(SX*ScatAdj,SY*ScatAdj,50,'k+','LineWidth',3)
grayColor = [.7 .7 .7];
for VerInd = 1:N_HC -1 % HC bounds
    plot(ones(length(0:n_S_HC*N_HC))*(VerInd*n_S_HC),0:n_S_HC*N_HC,'Color',grayColor,'LineWidth',1.5)%,1.5*n_S_HC/n_I_HC)
end
for HorInd = 1:N_HC -1
    plot(0:n_S_HC*N_HC,ones(length(0:n_S_HC*N_HC))*(HorInd*n_S_HC),'Color',grayColor,'LineWidth',1.5)%,1.5*n_S_HC/n_I_HC)
end
title([CellType{1} ' to ' CellType{2}])
axis([0.5*n_S_HC 2.5*n_S_HC 0.5*n_S_HC 2.5*n_S_HC])
axis square; axis off

end