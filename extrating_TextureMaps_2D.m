function TextMap = extrating_TextureMaps_2D(vol2D, mask, Ng, Nkern)

% This function extracts texture maps: Grey level co-occurrence matrix 
% (GLCM)[1] was computed in 2D for each pixel using the Nkern neighboring 
% pixels of each direction around it. The dynamic range of signal intensities
% was reduced to Ng gray levels. 13 Haralick feature maps were extracted:
% - Energy
% - Contrast
% - Entropy
% - Homogeneity
% - Dissimilarity
% - Correlation
% - Variance
% - SumAverage
% - AutoCorrelation
% - ClustTendency
% - SumEntropy
% - DiffVariance
% - DiffEntropy
% See IBSI for texture definitions [2] 

% References:
% [1]: Haralick, R. M.et al. (1973). Textural features for image classification.
%      IEEE Transactions on Systems, Man and Cybernetics, smc 3(6), 610–621.
% [2]: Zwanenburg A et al. Initiative for the Image biomarker 
%      standardisation initiative. arXiv:161207003. 2016 

% This code uses functions developped by M. Valliere (Image crop, Intensity  
% rescaling and GLCM calculation), available at 
% <https://github.com/mvallieres/radiomics/> 


%%%%% Inputs: %%%%%%%
% - vol2D: image to compute, in 2D. Can be easily modified for 3D
%           computation
% - mask: Mask of the ROI to compute. Same size as vol2D. Binary image with 
%           1: pixel to compute ; 0: no to consider          
% - Ng: Number of gray levels for intensity rescaling
% - Nkern: Kernel number: number of pixels per direction around 1 pixel 
%           to use to compute GLCM

% Author:   Marion Tardieu <marion.tardieu@umontpellier.fr>
% Ref:  Tardieu, M et al. Assessing histology structures by ex vivo 
%       MR microscopy and exploring the link between MRM-derived radiomic 
%       features and histopathology in ovarian cancer.


%% Image crop
% Use the function computeBoundingBox(mask) (M. Valliere)

[boxBound] = computeBoundingBox(mask);
maskBox = mask(boxBound(1,1)-Nkern:boxBound(1,2)+Nkern,boxBound(2,1)-Nkern:boxBound(2,2)+Nkern);
ROIonly = vol2D(boxBound(1,1)-Nkern:boxBound(1,2)+Nkern,boxBound(2,1)-Nkern:boxBound(2,2)+Nkern);
ROIonly(~maskBox) = NaN; ROIonly(maskBox<0) = NaN;


%% Intensity rescaling
% uniformQuantization(ROIonly,Ng) (M. Valliere)

ROIonly = double(ROIonly);
maxVal = max(ROIonly(:));
minVal = min(ROIonly(:));
VolRescNaN = round((Ng-1)*(ROIonly-minVal)/(maxVal-minVal))+1;

levels = 1:Ng;

%% GLCM computation

% QUANTIZATION EFFECTS CORRECTION (M. Vallieres)
% In case (for example) we initially wanted to have 64 levels, but due to
% quantization, only 60 resulted.
nLevel = length(levels);
if nLevel > 100
    adjust = 10000;
else
    adjust = 1000;
end
levelTemp = max(levels)+1;
VolResc=VolRescNaN;
VolResc(isnan(VolResc)) = levelTemp;
levels = [levels,levelTemp];

dim = size(VolResc);

qs = round(levels*adjust)/adjust;
lqs = length(qs);

VolResc=round(VolResc*adjust)/adjust;
VolRescNaN=round(VolRescNaN*adjust)/adjust;

boxSize=(Nkern*2+1)^2;

indVect = 1:Ng;
[colGrid,rowGrid] = meshgrid(indVect,indVect);


%% Texture maps extraction - M. Tardieu
TextMapT.Energy = zeros(size(VolResc));
TextMapT.Contrast = zeros(size(VolResc));
TextMapT.Entropy = zeros(size(VolResc));
TextMapT.Homogeneity = zeros(size(VolResc));
TextMapT.Dissimilarity = zeros(size(VolResc));
TextMapT.Correlation = zeros(size(VolResc));
TextMapT.Variance = zeros(size(VolResc));
TextMapT.SumAverage = zeros(size(VolResc));
TextMapT.AutoCorrelation = zeros(size(VolResc));
TextMapT.ClustTendency = zeros(size(VolResc));
TextMapT.SumEntropy = zeros(size(VolResc));
TextMapT.DiffVariance = zeros(size(VolResc));
TextMapT.DiffEntropy = zeros(size(VolResc));


for xpix = 1:size(VolResc,1)
    for ypix = 1:size(VolResc,2)
        
        %%%%%%% GLCM computation %%%%%%% 
        xbmin = max(1,xpix-Nkern);
        xbmax = min(xpix+Nkern,dim(1));
        ybmin = max(1,ypix-Nkern);
        ybmax = min(ypix+Nkern,dim(2));
        
        q2box = VolResc(xbmin:xbmax,ybmin:ybmax);
        q2boxNaN = VolRescNaN(xbmin:xbmax,ybmin:ybmax);
        if length(q2box(:))<boxSize
            continue;
        elseif sum(isnan(q2boxNaN(:)))>round(boxSize/2)
            continue;
        else
            GLCM=measure_GLCM_4map2D(q2box,lqs);
                  
            
            %%%%% average and standard deviation of px and py %%%
            ui = indVect*sum(GLCM,2); % Mean of row GLCMi (px) (or mux or rowMean)
            uj = indVect*sum(GLCM)'; % Mean of column GLCMj (py) (or muy or colMean)
            % as GLCM is symetric, normaly, ui=uj
            
            sigi = sqrt(sum((rowGrid(:)-ui).^2.*GLCM(:))); %(or rowStd)
            sigj = sqrt(sum((colGrid(:)-uj).^2.*GLCM(:))); %(or colStd)
            

            
            %%%%%%% p_x+y %%%%%%
            vstart = -(Ng-1);
            vstop = Ng-1;
            
            % Rotate Matrix 90°
            GLCM90 = rot90(GLCM);
            
            % Initilisiere p_x+y
            p_XplusY = zeros((2*Ng)-1,1);
            
            k = 1;
            for index = vstart : vstop
                p_XplusY(k) = sum(diag(GLCM90,index));
                k = k+1;
            end
            
            %%%%%%% p_x-y %%%%%%
            vstart = 1;
            vstop = Ng-1;
            
            % Initialize p_XminusY
            p_XminusY = zeros(Ng,1);
            p_XminusY(1) = sum (diag(GLCM,0) );
            
            k = 2;
            for index = vstart : vstop
                p_XminusY(k) = sum( [diag(GLCM,index);
                    diag(GLCM,-index)] );
                k = k + 1;
            end
            
            %%%%%%% Texture maps extraction %%%%%%%
            
            % Energy / Uniformity / second angular moment / Joint energy
            TextMapT.Energy(xpix,ypix) = sum(sum(GLCM.^2));
            
            % Contrast
            TextMapT.Contrast(xpix,ypix) = sum(sum((abs(rowGrid-colGrid).^2).*GLCM));
            
            % Entropy / Joint entropy / entropy log2
            TextMapT.Entropy(xpix,ypix) = -sum(sum(GLCM.*log2(GLCM + realmin)));
            
            % Homogeneity / Inverse difference
            TextMapT.Homogeneity(xpix,ypix) = sum(sum(GLCM./(1+abs(rowGrid-colGrid))));
            
            % Dissimilarity
            TextMapT.Dissimilarity(xpix,ypix) = sum(sum(abs(rowGrid-colGrid).*GLCM));
            
            % Correlation
            TextMapT.Correlation(xpix,ypix) = sum((rowGrid(:)-ui).*...
                (colGrid(:)-uj).*GLCM(:))/(sigi*sigj);
            
            % Variance / sum of squares / Joint Variance
            TextMapT.Variance(xpix,ypix) = sum(sum((rowGrid-ui).^2.*GLCM));
            
            % Sum Average
            TextMapT.SumAverage(xpix,ypix) = sum((2:(2*Ng))'.*p_XplusY);
            
            % AutoCorrelation
            TextMapT.AutoCorrelation(xpix,ypix) = sum(sum(rowGrid.*colGrid.*GLCM));
            
            % Cluster tendency - Sum variance
            TextMapT.ClustTendency(xpix,ypix) = sum(sum((rowGrid + colGrid - ui - uj)...
                .^2.*GLCM));
            
            % Sum entropy
            TextMapT.SumEntropy(xpix,ypix) = -sum(p_XplusY.*log2(p_XplusY + realmin));
            
            % Difference Variance
            k=(0:Ng-1)';
            mut=sum(k.*p_XminusY);
            TextMapT.DiffVariance(xpix,ypix) = sum((k-mut).^2.*p_XminusY);
            
            % Difference Entropy
            TextMapT.DiffEntropy(xpix,ypix) = -sum(p_XminusY.*log2(p_XminusY + realmin));
        end
        
    end
end

%% Resize Texture maps
TextMap.Energy = zeros(size(vol2D));
TextMap.Contrast = zeros(size(vol2D));
TextMap.Entropy = zeros(size(vol2D));
TextMap.Homogeneity = zeros(size(vol2D));
TextMap.Dissimilarity = zeros(size(vol2D));
TextMap.Correlation = zeros(size(vol2D));
TextMap.Variance = zeros(size(vol2D));
TextMap.SumAverage = zeros(size(vol2D));
TextMap.AutoCorrelation = zeros(size(vol2D));
TextMap.ClustTendency = zeros(size(vol2D));
TextMap.SumEntropy = zeros(size(vol2D));
TextMap.DiffVariance = zeros(size(vol2D));
TextMap.DiffEntropy = zeros(size(vol2D));



TextMap.Energy(boxBound(1,1)-Nkern:boxBound(1,2)+Nkern,boxBound(2,1)-Nkern:boxBound(2,2)+Nkern) =...
    TextMapT.Energy;
TextMap.Contrast(boxBound(1,1)-Nkern:boxBound(1,2)+Nkern,boxBound(2,1)-Nkern:boxBound(2,2)+Nkern) =...
    TextMapT.Contrast;
TextMap.Entropy(boxBound(1,1)-Nkern:boxBound(1,2)+Nkern,boxBound(2,1)-Nkern:boxBound(2,2)+Nkern) =...
    TextMapT.Entropy;
TextMap.Homogeneity(boxBound(1,1)-Nkern:boxBound(1,2)+Nkern,boxBound(2,1)-Nkern:boxBound(2,2)+Nkern) =...
    TextMapT.Homogeneity;
TextMap.Dissimilarity(boxBound(1,1)-Nkern:boxBound(1,2)+Nkern,boxBound(2,1)-Nkern:boxBound(2,2)+Nkern) =...
    TextMapT.Dissimilarity;
TextMap.Correlation(boxBound(1,1)-Nkern:boxBound(1,2)+Nkern,boxBound(2,1)-Nkern:boxBound(2,2)+Nkern) =...
    TextMapT.Correlation;
TextMap.Variance(boxBound(1,1)-Nkern:boxBound(1,2)+Nkern,boxBound(2,1)-Nkern:boxBound(2,2)+Nkern) =...
    TextMapT.Variance;
TextMap.SumAverage(boxBound(1,1)-Nkern:boxBound(1,2)+Nkern,boxBound(2,1)-Nkern:boxBound(2,2)+Nkern) =...
    TextMapT.SumAverage;
TextMap.AutoCorrelation(boxBound(1,1)-Nkern:boxBound(1,2)+Nkern,boxBound(2,1)-Nkern:boxBound(2,2)+Nkern) =...
    TextMapT.AutoCorrelation;
TextMap.ClustTendency(boxBound(1,1)-Nkern:boxBound(1,2)+Nkern,boxBound(2,1)-Nkern:boxBound(2,2)+Nkern) =...
    TextMapT.ClustTendency;
TextMap.SumEntropy(boxBound(1,1)-Nkern:boxBound(1,2)+Nkern,boxBound(2,1)-Nkern:boxBound(2,2)+Nkern) =...
    TextMapT.SumEntropy;
TextMap.DiffVariance(boxBound(1,1)-Nkern:boxBound(1,2)+Nkern,boxBound(2,1)-Nkern:boxBound(2,2)+Nkern) =...
    TextMapT.DiffVariance;
TextMap.DiffEntropy(boxBound(1,1)-Nkern:boxBound(1,2)+Nkern,boxBound(2,1)-Nkern:boxBound(2,2)+Nkern) =...
    TextMapT.DiffEntropy;

end




