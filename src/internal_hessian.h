/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2017 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'internal_hessian.h' is part of RASPA-2.0

    RASPA is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    RASPA is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *************************************************************************************************************/

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
