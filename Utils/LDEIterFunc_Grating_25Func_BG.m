%% Iterations of LDE: Use precomputed function to determine output of a cell type
% Input: L4EmeshX, L4ImeshY  Domain of precomputed functions
%        LDEFrfunc_Subf      Functions of all input type for certain cell
%        population
%        L4EUse,L4IUse       All input L4E L4I
%        LDEUse              Just for output formality
%        varargin:
%        ORTPix,OBLPix,OPTPix,ORTOBLBd_Pix,OPTOBLBd_Pix   Pixel index of different ODs
%        or
%        PixLGNCtgr: n*3

% Zhuo-Cheng Xiao 04/23/2023
% Find a way to incorporate BG local response functions:

function LDEOutS = LDEIterFunc_Grating_25Func_BG(...
    L4EmeshX, L4ImeshY,...
    LDEFrfunc_Subf,...
    L4EmeshXBG, L4ImeshYBG,...
    LDEFrfunc_SubfBG,...
    L4EUse,L4IUse,...
    PixInptCtgrUse)
    % first look up the library
    LDEOutLIBy = cell(size(LDEFrfunc_Subf));
    LibyAll = zeros(length(L4EUse),size(LDEFrfunc_Subf,1),size(LDEFrfunc_Subf,2));
    for LGNInd = 1:size(LDEFrfunc_Subf,1)
        for L6Ind = 1:size(LDEFrfunc_Subf,2)
        LDEOutLIBy{LGNInd,L6Ind} = ...
            interp2(L4EmeshX, L4ImeshY,LDEFrfunc_Subf{LGNInd,L6Ind}, L4EUse,L4IUse);%,'makima'
        LibyAll(:,LGNInd,L6Ind) = LDEOutLIBy{LGNInd,L6Ind};
        end
    end
    %% compute another output using the BG function:
    % There is always going to be NAN values. Just use the good entries!
    LDEOutBG = interp2(L4EmeshXBG, L4ImeshYBG,...
        LDEFrfunc_SubfBG, L4EUse,L4IUse);%,'makima'
%     if ~any(isnan(LDEOutBG),'all')
%         LibyAll(:,end,end) = LDEOutBG;
%         BGflag = true;
%     else
%         BGflag = false;
%     end
    LibyAll(~isnan(LDEOutBG),end,end) = LDEOutBG(~isnan(LDEOutBG));
    % Then compose
    LDEOutS = sum(PixInptCtgrUse.*LibyAll,[2,3]);




end