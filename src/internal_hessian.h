/*************************************************************************************************************
 The MIT License

 Copyright (C) 2006-2019 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

     D.Dubbeldam@uva.nl            http://www.uva.nl/profiel/d/u/d.dubbeldam/d.dubbeldam.html
     scaldia@upo.es                http://www.upo.es/raspa/
     t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
     don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
     snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

 Permission is hereby granted, free of charge, to any person
 obtaining a copy of this software and associated documentation
 files (the "Software"), to deal in the Software without
 restriction, including without limitation the rights to use,
 copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the
 Software is furnished to do so, subject to the following
 conditions:

 The above copyright notice and this permission notice shall be
 included in all copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 OTHER DEALINGS IN THE SOFTWARE.
**************************************************************************************************************/

#ifndef INTERNAL_HESSIAN_H
#define INTERNAL_HESSIAN_H

void CalculateAdsorbateBondHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian);
void CalculateCationBondHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian);
void CalculateAdsorbateBendHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian);
void CalculateCationBendHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian);
void CalculateAdsorbateTorsionHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian);
void CalculateCationTorsionHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian);
void CalculateAdsorbateImproperTorsionHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian);
void CalculateCationImproperTorsionHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian);

void CalculateAdsorbateBondBondHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian);

void CalculateAdsorbateIntraVDWHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian);
void CalculateCationIntraVDWHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian);
void CalculateAdsorbateIntraCoulombHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian);
void CalculateCationIntraCoulombHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian);

void CalculateAdsorbateInversionBendHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivativeTensor,int ComputeGradient,int ComputeHessian);
void CalculateCationInversionBendHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivativeTensor,int ComputeGradient,int ComputeHessian);

void CalculateHarmonicBondConstraintHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian);
void CalculateHarmonicBendConstraintHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian);
void CalculateHarmonicDihedralConstraintHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian);

void HessianBendStrainPosition(INT_VECTOR3 index_i,INT_VECTOR3 index_j,INT_VECTOR3 index_k,REAL_MATRIX HessianMatrix,VECTOR vec_u,VECTOR vec_v,REAL u,REAL v,
           REAL rab,REAL rbc,VECTOR Rab,VECTOR Rbc,VECTOR dtA,VECTOR dtC,REAL DF,REAL DDF,REAL_MATRIX3x3 S,REAL CosTheta);
void HessianBendStrainStrain(REAL_MATRIX HessianMatrix,VECTOR vec_u,VECTOR vec_v,REAL u,REAL v,REAL DF,REAL DDF,REAL_MATRIX3x3 S,REAL CosTheta);

void HessianTorsionStrainPosition(REAL_MATRIX HessianMatrix,INT_VECTOR3 index_i,INT_VECTOR3 index_j,INT_VECTOR3 index_k,INT_VECTOR3 index_l,
             VECTOR vec_u,VECTOR vec_v,VECTOR vec_w,REAL u,REAL v,REAL w,VECTOR fa,VECTOR fc,VECTOR fd,
             VECTOR Dab,VECTOR Dcb,VECTOR Ddc,REAL rbc,REAL_MATRIX3x3 D2I,REAL_MATRIX3x3 D2K,REAL_MATRIX3x3 D2L,
             REAL_MATRIX3x3 D2IJ,REAL_MATRIX3x3 D2IK,REAL_MATRIX3x3 D2IL,REAL_MATRIX3x3 D2JK,REAL_MATRIX3x3 D2JL,REAL_MATRIX3x3 D2KL,
             VECTOR dtA,VECTOR dtB,VECTOR dtC,VECTOR dtD,REAL DDF,REAL_MATRIX3x3 S,REAL CosPhi);
void HessianTorsionStrainStrain(REAL_MATRIX HessianMatrix,VECTOR vec_u,VECTOR vec_v,VECTOR vec_w,REAL u,REAL v,REAL w,REAL DF,REAL DDF,REAL_MATRIX3x3 S,REAL CosPhi);

#endif
