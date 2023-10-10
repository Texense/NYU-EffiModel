%% mean-field stuff. We use parameters and mean Vs to estimate firing rates
% Input: C_EE,C_EI,C_IE,C_II connectivity matrices
%        N_EE,N_EI,N_IE,N_II number of presynaptic neurons
%        S_EE,S_EI,S_IE,S_II synaptic strength
%        p_EEFail            Synaptic failure prob
%        lambda_E lambda_I   LGN input
%        rE_amb rI_amb       Ambient input
%        S_Elgn S_Ilgn S_amb Synaptic strength of drives
%        gL_E,gL_I           Leaky time constants
%        Ve,Vi               Reversak potentials
%        mVE,mVI             Mean Vs, collected from simulation
% Output:f_EnI               Estimation of firing rates, E and I

function f_EnI = MeanFieldEst_BkGd(N_EE,N_EI,N_IE,N_II,...
                                   S_EE,S_EI,S_IE,S_II,p_EEFail,...
                                   lambda_E,S_Elgn,rE_amb,S_amb,...
                                   lambda_I,S_Ilgn,rI_amb,...
                                   gL_E,gL_I,Ve,Vi,mVE,mVI)
                               
% E_sideInd = floor(1*n_E_HC+1):2*n_E_HC;
% [E_Ind_X,E_Ind_Y] = meshgrid(E_sideInd,E_sideInd);
% E_Ind = (reshape(E_Ind_X,size(E_Ind_X,1)*size(E_Ind_X,2),1)-1)*n_E_HC*N_HC + reshape(E_Ind_Y,size(E_Ind_X,1)*size(E_Ind_X,2),1);
% 
% I_sideInd = floor(1*n_I_HC+1):2*n_I_HC;
% [I_Ind_X,I_Ind_Y] = meshgrid(I_sideInd,I_sideInd);
% I_Ind = (reshape(I_Ind_X,size(I_Ind_X,1)*size(I_Ind_X,2),1)-1)*n_I_HC*N_HC + reshape(I_Ind_Y,size(I_Ind_X,1)*size(I_Ind_X,2),1);
% 
% % Take averaged number of input neurons 
% % NOTE! Only picking up the middle part
% N_EE = mean(sum(C_EE(E_Ind,:),2)); N_EI = mean(sum(C_EI(E_Ind,:),2)); 
% N_IE = mean(sum(C_IE(I_Ind,:),2)); N_II = mean(sum(C_II(I_Ind,:),2));

% Matrix of firing rates
MatEI = [N_EE*S_EE*(Ve-mVE)*(1-p_EEFail), N_EI*S_EI*(Vi-mVE);
         N_IE*S_IE*(Ve-mVI),              N_II*S_II*(Vi-mVI)];

% External Drive and leaky current. Note that simulaton is by ms, so need *1000
VecInput = [(Ve-mVE) * (lambda_E*S_Elgn + rE_amb*S_amb);
            (Ve-mVI) * (lambda_I*S_Ilgn + rI_amb*S_amb)]*1000;
Leak = -[gL_E*mVE;
        gL_I*mVI]*1000; % leak has opposite signs of mean Vs
    
f_EnI = -(MatEI-eye(2))^-1*(Leak+VecInput);

end