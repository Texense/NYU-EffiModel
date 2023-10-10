% ConnectionMat_Fix.m is a function generating connetivity matrices
% Input:
%      N_Post,NnPost,Size_Post: number, index projection, and grid size of
%                               presynaptic neurons
%      N_Pre,NnPre,Size_Pre:    Same for postsynaptic neurons
%      Peak_EE:                 Peak connectivity prob
%      SD_E:                    standard daviation of gaussian
%      Dist_UB:                 The upper bound of distances. Connections
%                               between neuron pairs with distances above this are ignored
%      SameCell:                Are pre and post cells are the same or not
%      NPreSyn:                 The numbers of presynaptic cells, which we will fix for all neurons
%      n_Pre_HC, n_Post_HC:     Number of E/I cells per HC
% Output:C_EE, the connectivity matrix. Only includes 0 and 1. Diagonal
% ignored for the E-E and I-I matrix.
%% This version we fix the presyn cell numbers...
function C_EE = ConnectionMat_Fix(N_Post,NnPost,Size_Post,...
                                  N_Pre,NnPre,Size_Pre,...
                                  Peak_EE,SD_E,Dist_UB,SameCell,...
                                  NPreSyn) 
% get location vectors for all post and pre cells
postNn_Loc = [NnPost.X(1:N_Post);NnPost.Y(1:N_Post)]*Size_Post;
preNn_Loc = [NnPre.X(1:N_Pre);NnPre.Y(1:N_Pre)]*Size_Pre;

% record all effective connection index
Ind1 = zeros(1,N_Post*1000); Ind2 = zeros(1,N_Post*1000); 
tic
subInd = 0;
for postNn = 1:N_Post % for each post synaptic cell, do
    % Distances between, and ignore too-far neurons
    Dist = sqrt(sum((repmat(postNn_Loc(:,postNn),1,N_Pre) - preNn_Loc).^2));
    [~,EffInd2] = find(Dist<=Dist_UB);
    EffDist = Dist(EffInd2);
    
    % Probability of having a projection
%     P_proj = sparse(1,EffInd2,Peak_EE * exp(-EffDist.^2/(2*SD_E^2)),1,N_Pre); 
%     Conn = rand(1,length(P_proj))<P_proj;
%      P_ProjUse = Peak_EE * exp(-EffDist.^2/(2*SD_E^2));
%      ConnPInd = rand(1,length(P_ProjUse))<P_ProjUse;
%      ConnInd2 = EffInd2(ConnPInd);
%      if length(ConnInd2)>=NPreSyn
%          ConnInd2  = randsample(ConnInd2,NPreSyn);
%          
%      end
    P_ProjUse = Peak_EE * exp(-EffDist.^2/(2*SD_E^2)); 
    ConnDeter = unique(randsample(EffInd2,NPreSyn*10,true,P_ProjUse/sum(P_ProjUse)),'stable');
                % get a unique list of presyn neurons
                % Presuambly ConnDeter contains more presyn neurons than we
                % need (but with repetance), so we just choose the first few we need
    ConnInd2 = ConnDeter(1:NPreSyn);
    
    % only valid connection indexes are preserved
    %[~,ConnInd2] = find(Conn);
    Ind1(subInd+1:subInd+length(ConnInd2)) = postNn*ones(size(ConnInd2));
    Ind2(subInd+1:subInd+length(ConnInd2)) = ConnInd2;
    
    subInd = subInd + length(ConnInd2);
end

% eliminate ineffective indexes, and assign connections
Ind1(Ind1 ==0) = [];
Ind2(Ind2 ==0) = [];
C_EE = sparse(Ind1,Ind2,ones(size(Ind1)));

if SameCell == 1
    C_EE = C_EE-diag(diag(C_EE));
end

toc
end