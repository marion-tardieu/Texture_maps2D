function GLCM = measure_GLCM_4map2D(vol_box,numbGrayL)

% This function computes the Gray-Level Co-occurence Matrix (GLCM) of the 
% region of interest (ROI) of a 2D image. It is adapted from the original
% code of M. Vallieres (2015) 'getGLCM.m' available at 
% <https://github.com/mvallieres/radiomics/>
% Ref: Haralick, R. M.et al. (1973). Textural features for image classification.
%      IEEE Transactions on Systems, Man and Cybernetics, smc 3(6), 610–621.
%     

% Inputs:
% - vol_box: square whose size depends on kernel number: (2*kern+1)^2. 
%            GLCM will be calculated on this square.
% - numbGrayL: Vector containing the quantized gray-levels in the tumor 
%           region (or reconstruction levels of quantization).
%           
% AUTHOR(S): 
% - Issam El Naqa <ielnaqa@med.umich.edu>       
% - Martin Vallieres <mart.vallieres@gmail.com>
% - Marion Tardieu <marion.tardieu@umontpellier.fr>
% HISTORY:
% - Creation: 2009 (I. El Naqa)
% - Revision 1: January 2013 (M. Vallieres)
% - Revision 2: May 2015 (M. Vallieres)
% - Revision 3: March 2021 (M. Tardieu)
% -------------------------------------------------------------------------
% STATEMENT:
% This file is part of <https://github.com/mvallieres/radiomics/>, 
% a package providing MATLAB programming tools for radiomics analysis.
% --> Copyright (C) 2015  Issam Elnaqa, Martin Vallieres
%
%    This package is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This package is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this package.  If not, see <http://www.gnu.org/licenses/>.
% -------------------------------------------------------------------------


GLCM = double(zeros(numbGrayL));
dim = size(vol_box);
for i = 1:dim(1)
    i_min = max(1,i-1);
    i_max = min(i+1,dim(1));
    for j = 1:dim(2)
        j_min = max(1,j-1);
        j_max = min(j+1,dim(2));

            for I2 = i_min:i_max
                for J2 = j_min:j_max
                        if I2 == i && J2 == j 
                            continue;
                        elseif isnan(vol_box(I2,J2))||isnan(vol_box(i,j))
                            continue;
                        else
                       
                            val_neighbor = vol_box(I2,J2);
                            val_q2 = vol_box(i,j);
                            GLCM(val_q2,val_neighbor) = GLCM(val_q2,val_neighbor) + sqrt(abs(I2-i)+abs(J2-j)); % Discretization length correction (M. Vallieres)
                        end

                end
            end
    end
end

GLCM = GLCM(1:end-1,1:end-1);
GLCM = GLCM/(sum(GLCM(:))); % normalization
end