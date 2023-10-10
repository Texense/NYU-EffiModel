%% Doing Gaussian blur to inputs: L6, lgn, etc...
% For new version of network for paper2
% Input: OD_SMap (n_S_HC * N_HC)^2 matrix, representing the spatial ord dom
%        L6S_Drive: Input to each domain
%        FiltSig: Spatial sigma, should be proportional to mat size
%        FigOn: boolean, export a figure or not
% Output: SpatInputFilt: vector of filtered input
% 10/29/2021

function SpatInputFilt = SpatialGaussianFilt(OD_SMap,L6S_Drive,FiltSig,FigOn)
[a,b] = size(OD_SMap);
if a~=b
    disp("Warning! Spatial map not square")
end

FiltMat = imgaussfilt(L6S_Drive(OD_SMap),FiltSig);
SpatInputFilt = reshape(FiltMat,a*b,1);

if FigOn
    figure
    imagesc(imgaussfilt(L6S_Drive(OD_SMap),FiltSig))
    colorbar
end

end