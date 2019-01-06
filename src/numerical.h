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
