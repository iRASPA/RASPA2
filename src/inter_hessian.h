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

#ifndef INTER_HESSIAN_H
#define INTER_HESSIAN_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

void PreComputeRotationDerivatives(void);
void HessianAtomicPositionPosition(REAL_MATRIX HessianMatrix,INT_VECTOR3 index_i,INT_VECTOR3 index_j,
                                                 REAL f1,REAL f2,VECTOR dr,REAL ReplicaFactor);
void HessianCenterOfMassOrientation(REAL_MATRIX HessianMatrix,INT_VECTOR3 index_i,INT_VECTOR3 index_i2,INT_VECTOR3 index_j,INT_VECTOR3 index_j2,
              int index1,int index2,REAL f1,REAL f2,VECTOR dr,REAL ReplicaFactor);
void HessianOrientationOrientation(REAL_MATRIX HessianMatrix,INT_VECTOR3 index_i,INT_VECTOR3 index_i2,INT_VECTOR3 index_j,INT_VECTOR3 index_j2,
             int index1,int index2,REAL f1,REAL f2,VECTOR dr,REAL ReplicaFactor);

void HessianCenterOfMassStrain(REAL_MATRIX HessianMatrix,INT_VECTOR3 index_i,INT_VECTOR3 index_j,
      REAL f1,REAL f2,VECTOR dr,VECTOR posA,VECTOR comA,VECTOR posB,VECTOR comB,int RigidA,int RigidB);
void HessianOrientationStrain(REAL_MATRIX HessianMatrix,INT_VECTOR3 index_i2,INT_VECTOR3 index_j2,int index1,int index2,
                  REAL f1,REAL f2,VECTOR posA,VECTOR comA,VECTOR posB,VECTOR comB,VECTOR dr);

void HessianAtomicStrainStrain(REAL_MATRIX HessianMatrix,INT_VECTOR3 index_i,INT_VECTOR3 index_j,
                       REAL f1,REAL f2,VECTOR dr,VECTOR posA,VECTOR comA,VECTOR posB,VECTOR comB);


void GradientStrainI(REAL *Gradient,REAL f1,VECTOR dr,VECTOR posA,VECTOR comA);
void GradientStrainJ(REAL *Gradient,REAL f1,VECTOR dr,VECTOR posB,VECTOR comB);


void ComputeInterVDWMolecularHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian);
void ComputeInterChargeChargeMolecularHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian);

void CalculateBondConstraintExclusionHessian(REAL *Energy,REAL* Gradient,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 *StrainDerivative,int ComputeGradient,int ComputeHessian);

#endif

