% ConnectionMat_Fix.m is a function generating connetivity matrices
% Input:
%      N_Post,NnPost,Size_Post: number, index projection, and grid size of
%                                 presynaptic neurons
%      N_Pre,NnPre,Size_Pre:    Same for postsynaptic neurons
%      NeEx, Map_Ext:           Extended index projection; maping from
%                                 original to extended

%      Peak_EE:                 Peak connectivity prob
%      SD_E:                    standard daviation of gaussian
%      Dist_UB:                 The upper bound of distances. Connections
%                                 between neuron pairs with distances above this are ignored
%      SameCell:                Are pre and post cells are the same or not
%      NPreSyn:                 The numbers of presynaptic cells, which we will fix for all neurons
%      n_Pre_HC, n_Post_HC:     Number of E/I cells per HC
% Output:
%      C_EE, the connectivity matrix. Only includes 0 and 1. Diagonal
%            signored for the E-E and I-I matrix.
%% This version we fix the presyn cell numbers...
%% Ver2: Modify boundaries: Let's cheat a little - Connect to closest 'right' HC
function C_EE = ConnectionMat_Fix_Boundary(N_Post,NnExPost,Size_Post,  PostMap_Ori2Ext,...
                                           N_Pre, NnExPre, Size_Pre,   PreMap_Ext2Ori,...
                                           Peak_EE,SD_E,Dist_LB,SameCell,...
                                           NPreSyn) 
%% Idea: 
% 1. we start from 5*5 (2 HCs larger), but only care about the center
% 3*3 for postsyn neurons.
% 2. Map "out of bound" presyn neurons to the right place (a generic map for N_HC??)
% 3. Map all cells back to 3*3 map
                                       
% get location vectors for all post and pre cells
postNn_Loc = [NnExPost.X(PostMap_Ori2Ext);
              NnExPost.Y(PostMap_Ori2Ext)]*Size_Post;
preNn_Loc = [NnExPre.X;
             NnExPre.Y]*Size_Pre;

% record all effective connection index
% Ind1 for post, Ind2 for pre
Ind1 = zeros(1,N_Post*1000); Ind2 = zeros(1,N_Post*1000); 
tic
subInd = 0;

for postNn = 1:N_Post % for each post synaptic cell, do
    % Distances between extended maps, 
    % and ignore too-far neurons
    Dist = sqrt(sum((repmat(postNn_Loc(:,postNn),1,size(preNn_Loc,2)) - preNn_Loc).^2));
    [~,EffInd2] = find(Dist<=Dist_LB);
    EffDist = Dist(EffInd2);
    
    % Probability of having a projection
    P_ProjUse = Peak_EE * exp(-EffDist.^2/(2*SD_E^2)); 
    ConnDeter = unique(randsample(EffInd2,NPreSyn*10,true,P_ProjUse/sum(P_ProjUse)),'stable');
                % get a unique list of presyn neurons
                % Presuambly ConnDeter contains more presyn neurons than we
                % need (but with repetance), so we just choose the first few we need
    ConnInd2 = ConnDeter(1:NPreSyn);
    % Now, map pre ind back to original 
    
    % only valid connection indexes are preserved
    Ind1(subInd+1:subInd+length(ConnInd2)) = postNn*ones(size(ConnInd2)); % Use original post Ind
    Ind2(subInd+1:subInd+length(ConnInd2)) = PreMap_Ext2Ori(ConnInd2);
    
    subInd = subInd + length(ConnInd2);
end

% eliminate ineffective indexes, and assign connections
Ind1(Ind1 == 0) = [];
Ind2(Ind2 == 0) = [];
C_EE = sparse(Ind1,Ind2,ones(size(Ind1)));

if SameCell == 1
    C_EE = C_EE-diag(diag(C_EE));
end

toc
end