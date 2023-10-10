%% mean-field inverse estimation: We use parameters and desired firing rates to estimate input LGN rate
% Input: C_EE,C_EI,C_IE,C_II connectivity matrices
%        S_EE,S_EI,S_IE,S_II synaptic strength
%        p_EEFail            Synaptic failure prob
%        f_E,f_I             Desired firing rates
%        rE_amb rI_amb       Ambient input
%        S_Elgn S_Ilgn S_amb Synaptic strength of drives
%        gL_E,gL_I           Leaky time constants
%        Ve,Vi               Reversak potentials
%        mVE,mVI             Mean Vs, collected from simulation
% Output:lambda_EnI          LGN Input required, E and I

function lambda_EnI = MeanFieldEstInv_BkGd(C_EE,C_EI,C_IE,C_II,...
                                   S_EE,S_EI,S_IE,S_II,p_EEFail,...
                                   S_Elgn,rE_amb,S_amb,...
                                   S_Ilgn,rI_amb,...
                                   gL_E,gL_I,Ve,Vi,mVE,mVI,...
                                   f_E,f_I)
% Take averaged number of input neurons
N_EE = mean(sum(C_EE,2)); N_EI = mean(sum(C_EI,2)); N_IE = mean(sum(C_IE,2)); N_II = mean(sum(C_II,2));

% Matrix of firing rates
MatEI = [N_EE*S_EE*(Ve-mVE)*(1-p_EEFail), N_EI*S_EI*(Vi-mVE);
         N_IE*S_IE*(Ve-mVI),              N_II*S_II*(Vi-mVI)];

% External Drive and leaky current. Note that simulaton is by ms, so need *1000
Leak = -[gL_E*mVE;
        gL_I*mVI]*1000; % leak has opposite signs of mean Vs
    
f_EnI = [f_E;f_I];
VecInput = -(MatEI-eye(2))*f_EnI - Leak;

lambda_EnI = (VecInput./[(Ve-mVE);(Ve-mVI)] - [rE_amb*S_amb;rI_amb*S_amb])./[S_Elgn;S_Ilgn];

end

%VecInput = [(Ve-mVE) * (lambda_E*S_Elgn + rE_amb*S_amb);
%            (Ve-mVI) * (lambda_I*S_Ilgn + rI_amb*S_amb)]*1000;