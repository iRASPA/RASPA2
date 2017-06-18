/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2017 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'spectra.h' is part of RASPA-2.0

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
