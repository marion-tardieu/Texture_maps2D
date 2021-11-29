# Texture_maps2D
This function extract texture maps from a 2D image. Grey level co-occurrence matrix (GLCM)[1] was computed in 2D for each pixel using the Nkern neighboring pixels of each direction around it. The dynamic range of signal intensities was reduced to Ng gray levels. 13 Haralick feature maps were extracted (See IBSI for texture feature definitions [2]) :
- Energy
- Contrast
- Entropy
- Homogeneity
- Dissimilarity
- Correlation
- Variance
- SumAverage
- AutoCorrelation
- ClustTendency
- SumEntropy
- DiffVariance
- DiffEntropy  


[1]: Haralick, R. M.et al. (1973). Textural features for image classification. IEEE Transactions on Systems, Man and Cybernetics, smc 3(6), 610–621.  
[2]: Zwanenburg A et al. Initiative for the Image biomarker standardisation initiative. arXiv:161207003. 2016   

Author: Marion Tardieu <marion.tardieu@umontpellier.fr>   
Reference: Tardieu, M et al. Assessing histology structures by ex vivo MR microscopy and exploring the link between MRM-derived radiomic features and histopathology in ovarian cancer.

Main function:   
- extrating_TextureMaps_2D.m  
Other functions:   
- measure_GLCM_4map2D.m => computes the Gray-Level Co-occurence Matrix (GLCM) of the region of interest (ROI) of a 2D image. Adapted from the original code of M. Vallieres (2015) 'getGLCM.m' - <https://github.com/mvallieres/radiomics/>  
- computeBoundingBox.m => computes the smallest box containing the whole ROI. From M. Vallieres (2015) available at <https://github.com/mvallieres/radiomics/>  



Aknowledgments:  
M. Vallieres - Radiomics toolbox <https://github.com/mvallieres/radiomics/> 
