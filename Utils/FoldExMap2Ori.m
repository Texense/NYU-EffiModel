% Folding Maps back:
% This only works for 3*3 and 5*5 for now
% Can be easily extended to odd numbers.

function EMap_Ext2Ori = FoldExMap2Ori(NnEex,NnE,EMap_Ori2Ext,N_E,n_E_HC)
% from extended back to original: fold
% E
EMap_Ext2Ori = zeros(size(NnEex.X));
EMap_Ext2Ori(EMap_Ori2Ext) = 1:N_E;
CenterEHC = find(NnE.X>1*n_E_HC & NnE.X<=2*n_E_HC ...
               & NnE.Y>1*n_E_HC & NnE.Y<=2*n_E_HC);
CenterEHCRow = find(NnE.Y>1*n_E_HC & NnE.Y<=2*n_E_HC);
CenterEHCCol = find(NnE.X>1*n_E_HC & NnE.X<=2*n_E_HC);
EMap_Ext2Ori(NnEex.X>0*n_E_HC & NnEex.X<=1*n_E_HC ...
           & NnEex.Y>0*n_E_HC & NnEex.Y<=1*n_E_HC) = CenterEHC; % LL corner
EMap_Ext2Ori(NnEex.X>4*n_E_HC & NnEex.X<=5*n_E_HC ...
           & NnEex.Y>0*n_E_HC & NnEex.Y<=1*n_E_HC) = CenterEHC; % RL corner
EMap_Ext2Ori(NnEex.X>0*n_E_HC & NnEex.X<=1*n_E_HC ...
           & NnEex.Y>4*n_E_HC & NnEex.Y<=5*n_E_HC) = CenterEHC; % LH corner
EMap_Ext2Ori(NnEex.X>4*n_E_HC & NnEex.X<=5*n_E_HC ...
           & NnEex.Y>4*n_E_HC & NnEex.Y<=5*n_E_HC) = CenterEHC; % RH corner
EMap_Ext2Ori(NnEex.X>1*n_E_HC & NnEex.X<=4*n_E_HC ...
           & NnEex.Y>0*n_E_HC & NnEex.Y<=1*n_E_HC) = CenterEHCRow; % Low
EMap_Ext2Ori(NnEex.X>1*n_E_HC & NnEex.X<=4*n_E_HC ...
           & NnEex.Y>4*n_E_HC & NnEex.Y<=5*n_E_HC) = CenterEHCRow;    % High
EMap_Ext2Ori(NnEex.X>0*n_E_HC & NnEex.X<=1*n_E_HC ...
           & NnEex.Y>1*n_E_HC & NnEex.Y<=4*n_E_HC) = CenterEHCCol; % Left
EMap_Ext2Ori(NnEex.X>4*n_E_HC & NnEex.X<=5*n_E_HC ...
           & NnEex.Y>1*n_E_HC & NnEex.Y<=4*n_E_HC) = CenterEHCCol;    % Right


end