function h = CorrCoefFig(FrETemp, FrITemp, mVETemp, mVITemp)

h = figure('Name','Correlation');

subplot 231; 
scatter(FrETemp, mVETemp, 0.1, 'r.','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2); 
Corr1 = corrcoef(FrETemp, mVETemp);
title(sprintf('Corr = %.2f', Corr1(2,1)))
xlabel('FrE');ylabel('mVE')
subplot 232; 
scatter(FrETemp, mVITemp, 0.1, 'r.','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2); 
Corr1 = corrcoef(FrETemp, mVITemp);
title(sprintf('Corr = %.2f', Corr1(2,1)))
xlabel('FrE');ylabel('mVI')
subplot 233; 
scatter(FrETemp, FrITemp, 0.1, 'r.','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2); 
Corr1 = corrcoef(FrETemp, FrITemp);
title(sprintf('Corr = %.2f', Corr1(2,1)))
xlabel('FrE');ylabel('FrI')
subplot 234; 
scatter(FrITemp, mVETemp, 0.1, 'b.','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2); 
Corr1 = corrcoef(FrITemp, mVETemp);
title(sprintf('Corr = %.2f', Corr1(2,1)))
xlabel('FrI');ylabel('mVE')
subplot 235; 
scatter(FrITemp, mVITemp, 0.1, 'b.','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2); 
Corr1 = corrcoef(FrITemp, mVITemp);
title(sprintf('Corr = %.2f', Corr1(2,1)))
xlabel('FrI');ylabel('mVI')
subplot 236; 
scatter(mVITemp, mVETemp, 0.1, 'g.','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2); 
Corr1 = corrcoef(mVITemp, mVETemp);
title(sprintf('Corr = %.2f', Corr1(2,1)))
xlabel('mVI');ylabel('mVE')
end