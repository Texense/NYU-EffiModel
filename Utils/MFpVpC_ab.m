%% mean-field est With ref. We use parameters and mean Vs to estimate firing rates
% Input: C_EE,C_EI,C_IE,C_II connectivity matrices
%        N_EE,N_EI,N_IE,N_II number of presynaptic neurons
%        S_EE,S_EI,S_IE,S_II synaptic strength
%        S_EL6,S_IL6,rE_L6,rI_L6 L6 Input
%        p_EEFail            Synaptic failure prob
%        lambda_E lambda_I   LGN input
%        rE_amb rI_amb       Ambient input
%        S_Elgn S_Ilgn S_amb Synaptic strength of drives
%        gL_E,gL_I           Leaky time constants
%        Ve,Vi               Reversak potentials
%        mVE,mVI             Mean Vs, collected from simulation
%        tau_ref             ref period, in ms
%        f_pre               from previous
% Output:f_EnI               Estimation of firing rates, E and I

function f_EnI = MFpVpC_ab(N_EE,N_EI,N_IE,N_II,...
                                   S_EE,S_EI,S_IE,S_II,p_EEFail,...
                                   lambda_E,S_Elgn,rE_amb,S_amb,...
                                   lambda_I,S_Ilgn,rI_amb,...
                                   S_EL6,S_IL6,rE_L6,rI_L6,...
                                   gL_E,gL_I,Ve,Vi,mVE,mVI,...
                                   tau_ref,f_pre,...
                                   aE, aI, bE, bI) % The amplitudes for corrections
                               
% Downplay current by a ref factor
ref = diag(1-f_pre*tau_ref/1000); ref_c = diag(-1/2*tau_ref/1000*f_pre.*[aE;aI]);
% Matrix of firing rates
ConEnI = [N_EE*S_EE*(1-p_EEFail), N_EI*S_EI;
          N_IE*S_IE,              N_II*S_II];
Mat = ConEnI .* [(Ve-mVE)-1/2*aE*bE, (Vi-mVE)-1/2*aI*bE;
                 (Ve-mVI)-1/2*aE*bI, (Vi-mVI)-1/2*aI*bI];

Mat_c = ConEnI .* [aE*(Ve-mVE)-bE, aI*(Vi-mVE)-bE;
                   aE*(Ve-mVI)-bI, aI*(Vi-mVI)-bI];

% External Drive and leaky current. Note that simulaton is by ms, so need *1000
ExtInp = [(lambda_E*S_Elgn + rE_amb*S_amb + rE_L6*S_EL6);
          (lambda_I*S_Ilgn + rI_amb*S_amb + rI_L6*S_IL6)]; 
ExtL = [ExtInp(1)*(Ve-mVE) - gL_E*mVE;
        ExtInp(2)*(Ve-mVI) - gL_I*mVI]*1000;
ExtL_c = [-ExtInp(1)*bE - gL_E*bE;
          -ExtInp(2)*bI - gL_I*bI]*1000;
% Solve the linear equations    
f_EnI = (ref*Mat + ref_c*Mat_c - eye(2))^-1 * (-ref*ExtL - ref_c*ExtL_c);    
%f_EnI = -(ref_fac*MatEI-eye(2))^-1*ref_fac*(Leak+VecInput);

end