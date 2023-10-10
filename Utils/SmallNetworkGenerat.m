%% Generate a smaller ring network (dense, sort of)

function S_SWconn = SmallNetworkGenerat(N_pre,N_post,N_conn)

S_SWconn = zeros(N_post,N_pre);
for postInd = 1:N_post
     S_SWconn(postInd,randsample(N_pre,round(N_conn))) = 1;
end

end