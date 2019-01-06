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

#ifndef SPECTRA_H
#define SPECTRA_H

enum{DIFFRACTION_GAUSSIAN,DIFFRACTION_LORENTZIAN,DIFFRACTION_PSEUDO_VOIGT};
enum{XRAY_DIFFRACTION,NEUTRON_DIFFRACTION,ELECTRON_DIFFRACTION};
enum{CHROMIUM_RADIATION,IRON_RADIATION,COPPER_RADIATION,MOLYBDENUM_RADIATION,SILVER_RADIATION,SYNCHROTRON_RADIATION};
enum{DIFFRACTION_SINGLE,DIFFRACTION_DOUBLET};
enum{HESSIAN_NUMERICAL,HESSIAN_ANALYTICAL};

extern int ComputePowderDiffractionPattern;

extern int ComputeNormalModes;
extern int MinimumMode;
extern int MaximumMode;
extern int ModeResolution;
extern int CorrectNormalModesForConstraints;

typedef struct diffraction
{
  int Type;
  int RadiationType;
  REAL two_theta_min;
  REAL two_theta_max;
  REAL two_theta_step;
  REAL lambda;
  REAL lambda2;
  int lambda_type;
  REAL w,v,u;
  REAL hmax,kmax,lmax;
  REAL asym;
  int PeakShape;
  int n;
  REAL *spectrum;
} DIFFRACTION;

extern DIFFRACTION Diffraction;

void VibrationalAnalysis(void);
void MassWeightHessianMatrix(int n,REAL_MATRIX Hessian,REAL *Weights);
int RemoveTranslationAndRotationFromHessianMatrix(REAL_MATRIX HessianMatrix,REAL *Weights,REAL *Positions);
int RemoveTranslationFromHessianMatrix(REAL_MATRIX HessianMatrix,REAL *Weights,REAL *Positions);
void SolveEigenValuesAndVectorsHessian(REAL_MATRIX HessianMatrix,REAL* Frequencies);
void SymmetrizeHessianMatrix(REAL_MATRIX Hessian);
void WriteVibrationalData(REAL_MATRIX HessianMatrix,REAL* Frequencies,REAL *Charges);
void NormalModeMonteCarlo(REAL_MATRIX HessianMatrix,REAL *Positions,REAL *Eigenvalues,REAL* Weights,CUBIC_SPLINE *Splines);

void WriteIR(void);
void ConstructNumericalHessianMatrix(REAL_MATRIX HessianMatrix);

int ProjectConstraintsFromHessianMatrix(int n,int np,REAL *Gradient,REAL_MATRIX HessianMatrix,int ComputeGradient,int ComputeHessian);
int ProjectConstraintsFromHessianMatrixMassWeighted(int np,int nb,REAL *Gradient,REAL_MATRIX HessianMatrix,REAL *Weights);

void PowderDiffraction(void);
 #endif
