%% Doing Gaussian blur to inputs: L6, lgn, etc...
% For new version of network for paper2
% Input: OD_SMap (n_S_HC * N_HC)^2 matrix, representing the spatial ord dom
%        L6S_Drive: Input to each domain
%        FiltSig: Spatial sigma, should be proportional to mat size
%        Trunc: Truncating. relative to sigma
%        FigOn: boolean, export a figure or not
% Output: SpatInputFilt: vector of filtered input
% 10/29/2021

function SpatInputFilt = SpatialGaussianFilt_my(OD_SMap,N_HC,n_S_HC,FiltSig,Trunc,FigOn)
[a,b] = size(OD_SMap);
if a~=b
    disp("Warning! Spatial map not square")
end
%% get a truncated Gaussian Kernel
Kersize = floor(2*FiltSig*Trunc);
if Kersize<=0
    Kersize= 2; % at least 1
end
[x,y] = meshgrid(1:Kersize,1:Kersize);
c = Kersize/2;
exponent = ((x-c).^2+(y-c).^2)/(2*FiltSig^2);
TruncGausKer = exp(-exponent);
TruncGausKer(TruncGausKer<exp(-(Trunc*1.0)^2/2))=0; % truncate
TruncGausKer = TruncGausKer + ...
               TruncGausKer(end:-1:1,:) + ...
               TruncGausKer(:,end:-1:1) + ...
               TruncGausKer(end:-1:1,end:-1:1);%symmetrize
TruncGausKer = TruncGausKer/sum(TruncGausKer,'all'); % normalize
%% Extend OD_SMap for one HC each side. Assume that the Gaussian Ker is not too large
% first extract a 2-by-2 pattern from northeast
cut1 = floor(Kersize/2); cut2 = Kersize-cut1-1;

Pattern22 = OD_SMap((N_HC-2)*n_S_HC+1:end,(N_HC-2)*n_S_HC+1:end);
N_expand = ceil(N_HC/2)+1;
Pattern_expand = repmat(Pattern22,N_expand,N_expand);

CtgrList = unique(OD_SMap);
SpatInputFilt = zeros(a*b,length(CtgrList));
for CtgrInd = 1:length(CtgrList)
    Pattern_1ctgr = zeros(size(Pattern_expand));
    Pattern_1ctgr(Pattern_expand == CtgrList(CtgrInd)) = 1;
    Pattern_1ctgr_Gauss = conv2(Pattern_1ctgr,TruncGausKer);
    Pattern_1ctgr_Gauss = Pattern_1ctgr_Gauss(cut1+1:end-cut2,cut1+1:end-cut2); % cut exteriors
    % Only get the original part
    Pattern_1ctgr_Gauss = Pattern_1ctgr_Gauss(n_S_HC+1:(N_HC+1)*n_S_HC,n_S_HC+1:(N_HC+1)*n_S_HC);
    SpatInputFilt(:,CtgrInd) = reshape(Pattern_1ctgr_Gauss,a*b,1);
end

if FigOn
    figure
    for CtgrInd = 1:length(CtgrList)
        subplot(1,length(CtgrList),CtgrInd)
        ShowField(SpatInputFilt(:,CtgrInd), 1:a*b,a,b,'',[],n_S_HC)
%         set(gca,'YDir','Normal')
%         colorbar
%         axis square
    end
end

end