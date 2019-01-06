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

#ifndef RIGID_H
#define RIGID_H

#include "constants.h"

enum{UNCONSTRAINED,CONSTRAIN_COM_POSITION,CONSTRAIN_ORIENTATION,CONSTRAIN_COM_AND_ORIENTATION};

enum{DISTANCE_R,DISTANCE_R_SQUARED};
enum{THETA,COS_THETA,COS_THETA_SQUARED};
enum{PHI,COS_PHI,COS_PHI_SQUARED};
enum{CHI,SIN_CHI,SIN_CHI_SQUARED};

extern int DistanceConstraintType;
extern int BendConstraintType;
extern int DihedralConstraintType;
extern int ImproperDihedralConstraintType;
extern int InversionBendConstraintType;
extern int OutOfPlaneConstraintType;

extern int ComputeRattleSteps;
extern REAL *NumberOfRattleCyclesStage1;
extern int *MaximumNumberOfRattleCyclesStage1;
extern REAL *NumberOfRattleCyclesStage2;
extern int *MaximumNumberOfRattleCyclesStage2;

void BuildRotationMatrix(REAL_MATRIX3x3 *M,QUATERNION q);
void BuildRotationMatrixInverse(REAL_MATRIX3x3 *M,QUATERNION q);
void ComputeQuaternions(void);
void ComputeQuaternionAdsorbate(int m);
void ComputeQuaternionCation(int m);

void CreateComVelocities(void);
void ComputeMolecularVelocitiesAndPositions(void);
void ComputeAngularVelocities(void);
void ComputeQuaternionMomenta(void);

void AtomicToMolecularVelocities(void);

QUATERNION AngularVelocityToQuaternionMomentumAdsorbates(int i,int g);
QUATERNION AngularVelocityToQuaternionMomentumCations(int i,int g);

void ConvertAngularVelocityToAtomicVelocity(void);

VECTOR QuaternionMomentumToAngularVelocityAdsorbates(int i,int g);
VECTOR QuaternionMomentumToAngularVelocityCations(int i,int g);
VECTOR AtomicVelocityToAngularVelocityAdsorbates(int i,int g);
VECTOR AtomicVelocityToAngularVelocityCations(int i,int g);

void ComputeEulerAxisFromQuaternionAdsorbate(int m);
void ComputeEulerAxisFromQuaternionCation(int m);

void AdjustSystemAngularRotationToZero(void);
VECTOR TotalAngularMomentum(void);

void ComputeDegreesOfFreedom(void);

void CreateCartesianPositions(void);
void CreateCartesianVelocities(void);
void NoSquishRotate(int k,REAL dt);
void NoSquishFreeRotorOrderTwo(void);
void NoSquishFreeRotorOrderFour(void);
void CreateQuaternionAndComForces(void);
void NoSquishEvolve0Todt2(void);
void NoSquishEvolvedt2Todt(void);

REAL ReturnConstraintDistance(VECTOR posA,VECTOR posB);
REAL ReturnWilsonVectorsDistanceRATTLE(VECTOR posA,VECTOR posB,VECTOR *wa,VECTOR *wb);
REAL ReturnDistanceConstrainDerivative(VECTOR Rab,VECTOR Vab);

REAL ReturnConstraintBendAngle(VECTOR posA,VECTOR posB,VECTOR posC);
REAL ReturnWilsonVectorsBendRATTLE(VECTOR posA,VECTOR posB,VECTOR posC,VECTOR *wa,VECTOR *wb,VECTOR *wc);
REAL ReturnAngleConstrainDerivative(VECTOR Rij,VECTOR Rjk,VECTOR Vij,VECTOR Vjk);

REAL ReturnConstraintDihedralAngle(VECTOR posA,VECTOR posB,VECTOR posC,VECTOR posD);
REAL ReturnWilsonVectorsTorsionRATTLE(VECTOR posA,VECTOR posB,VECTOR posC,VECTOR posD,VECTOR *wa,VECTOR *wb,VECTOR *wc,VECTOR *wd);
REAL ReturnTorsionConstrainDerivative(VECTOR Rij,VECTOR Rjk,VECTOR Rkl,VECTOR Vij,VECTOR Vjk,VECTOR Vkl);

REAL ReturnConstraintInversionBendAngle(VECTOR posA,VECTOR posB,VECTOR posC,VECTOR posD);
REAL ReturnWilsonVectorsInversionBendRATTLE(VECTOR posA,VECTOR posB,VECTOR posC,VECTOR posD,VECTOR *wa,VECTOR *wb,VECTOR *wc,VECTOR *wd);
REAL ReturnInversionBendConstrainDerivative(VECTOR Rij,VECTOR Rjk,VECTOR Rkl,VECTOR Vij,VECTOR Vjk,VECTOR Vkl);

REAL ReturnConstraintOutOfPlaneDistanceAngle(VECTOR posA,VECTOR posB,VECTOR posC,VECTOR posD);
REAL ReturnWilsonVectorsOutOfPlaneDistanceRATTLE(VECTOR posA,VECTOR posB,VECTOR posC,VECTOR posD,VECTOR *wa,VECTOR *wb,VECTOR *wc,VECTOR *wd);
REAL ReturnOutOfPlaneDistanceConstrainDerivative(VECTOR Rij,VECTOR Rjk,VECTOR Rkl,VECTOR Vij,VECTOR Vjk,VECTOR Vkl);

void ShakeInMinimization(void);

void RattleStageZero(void);
void RattleStageOne(void);
void RattleStageTwo(void);

REAL_MATRIX3x3 ComputeRotationMatrix(VECTOR p);
REAL_MATRIX3x3 ComputeRotationMatrixDerivativeX(VECTOR p);
REAL_MATRIX3x3 ComputeRotationMatrixDerivativeY(VECTOR p);
REAL_MATRIX3x3 ComputeRotationMatrixDerivativeZ(VECTOR p);

REAL_MATRIX3x3 ComputeRotationMatrixSecondDerivativeAX(VECTOR p);
REAL_MATRIX3x3 ComputeRotationMatrixSecondDerivativeBY(VECTOR p);
REAL_MATRIX3x3 ComputeRotationMatrixSecondDerivativeCZ(VECTOR p);
REAL_MATRIX3x3 ComputeRotationMatrixSecondDerivativeAY(VECTOR p);
REAL_MATRIX3x3 ComputeRotationMatrixSecondDerivativeAZ(VECTOR p);
REAL_MATRIX3x3 ComputeRotationMatrixSecondDerivativeBZ(VECTOR p);

void CalculateConstraintsExclusionEnergy(void);
void CalculateConstraintsExclusionForce(void);

#endif
