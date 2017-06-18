/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2017 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'numerical.h' is part of RASPA-2.0

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

#ifndef NUMERICAL_H
#define NUMERICAL_H

void SaveFrameworkPositionsToReferenceValues(void);
void PlaceFrameworkInBoxFromReferenceValues(void);
void SaveAdsorbateAtomPositionsToReferenceValues(void);
void PlaceAdsorbateAtomsInBoxFromReferenceValues(void);
void SaveCationAtomPositionsToReferenceValues(void);
void PlaceCationAtomsInBoxFromReferenceValues(void);

REAL_MATRIX3x3 ComputeStrainDerivativeNumerically(void);
REAL_MATRIX9x9 ComputeStrainSecondDerivativeNumerically(void);
void ComputeHessianMatrixNumerically(REAL_MATRIX HessianMatrix);
void ConstructCrossTermNumerically(void);
void ComputeCrossTermNumerically(REAL_MATRIX CrossTerm);

void ComputeGradientsNumerically(REAL *Gradients);

void CheckStatusNumerically(void);
REAL_MATRIX9x9 ComputeRelaxationTerm(int NumberOfPositionVariables,int NumberOfBoxVariables);
void ConstrucGeneralizedHessianMatrix(REAL_MATRIX GeneralizedHessianMatrix);
void ComputeNormalModeDerivativeNumerically(REAL_MATRIX GeneralizedHessianMatrix,REAL *Eigenvalues,REAL *Positions,
               CUBIC_SPLINE *splines,REAL_MATRIX3x3 StoredBox,REAL_MATRIX3x3 StoredInverseBox);

void ComputeCrossTermNumericallyMinimalSet(REAL_MATRIX CrossTerm);

void ComputeThirdOrderElasticConstantsNumerically(int NumberOfPositionVariables,int NumberOfBoxVariables,REAL_MATRIX6x6x6 *VoigtMatrixThirdOrder);

void ComputeEnergyGradientHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX3x3 *Strain_Derivative_Tensor,
                                  REAL_MATRIX HessianMatrix,REAL_MATRIX CrossTerm);


void AddRemainderOfCrossTermNumerically(REAL_MATRIX HessianMatrix);
void AddRemainderOfBornTermNumerically(REAL_MATRIX HessianMatrix);

 #endif
