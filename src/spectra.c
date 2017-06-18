/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2017 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'spectra.c' is part of RASPA-2.0

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

#ifdef HAVE_CONFIG_H
#include <config.h>
//#include <Accelerate/Accelerate.h>
#endif
#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include "constants.h"
#include "framework.h"
#include "framework_energy.h"
#include "framework_force.h"
#include "framework_hessian.h"
#include "inter_force.h"
#include "internal_force.h"
#include "internal_hessian.h"
#include "inter_hessian.h"
#include "molecule.h"
#include "potentials.h"
#include "spectra.h"
#include "integration.h"
#include "minimization.h"
#include "matrix.h"
#include "thermo_baro_stats.h"
#include "ewald.h"
#include "numerical.h"
#include "statistics.h"
#include "output.h"
#include "movies.h"
#include "cbmc.h"
#include "mc_moves.h"
#include "scattering_factors.h"
#include "rigid.h"

#define MAX_POWDER_DIFFRACTION_PEAKS 500000

extern bool STREAM;

void MakeNormalModeMovie(int NumberOfPositionVariables,int NumberOfBoxVariables,REAL_MATRIX HessianMatrix,REAL *Positions,REAL *frequencies,REAL* Weights);

void TestNormalModeApproximation(REAL_MATRIX HessianMatrix,REAL *Positions,REAL *Eigenvalues,int mode,REAL* Weights,CUBIC_SPLINE *Splines);

REAL LeftFrequencyBoundarySpectrum;
REAL RightFrequencyBoundarySpectrum;
int SpectrumWidth;

static int NumberOfRemovedModes;
static int NumberOfPositionVariables;
static int NumberOfBoxVariables;
static int NumberOfVariables;

int ComputeNormalModes;
int MinimumMode;
int MaximumMode;
int ModeResolution;
int CorrectNormalModesForConstraints;

int ComputePowderDiffractionPattern;
DIFFRACTION Diffraction;

typedef struct powder_diffraction_peaks
{
  int status;
  int h;       // Miller indices h,k,l
  int k;
  int l;
  REAL two_theta;  // scattering angle
  REAL two_theta2; // scattering angle
  REAL sintheta;           // sin(theta)
  REAL tantheta;           // tan(theta)
  int Multiplicity;        // the number of equivalent planes
  REAL d;                  // interplana spacing
  REAL s;                  // sin(theta)/lambda
  REAL ScatteringFactor;   // the scattering factor (including anomalous
  REAL Lp;                 // the Lorentz and Polarization factor
  REAL Intensity;          // the total intensity of the peak
} POWDER_DIFF_PEAK;

int NumberOfPeaks;
POWDER_DIFF_PEAK *PowderDiffractionPeaks;


void MassWeightHessianMatrix(int n,REAL_MATRIX Hessian,REAL *Weights)
{
  int i,j;

  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      Hessian.element[i][j]*=(Weights[i]*Weights[j]);
}

void SymmetrizeHessianMatrix(REAL_MATRIX Hessian)
{
  int i,j;

  for(i=0;i<Hessian.m;i++)
    for(j=i+1;j<Hessian.n;j++)
      Hessian.element[j][i]=Hessian.element[i][j];
}

/*********************************************************************************************************
 * Name       | RemoveShellInteractionsFromHessian                                                       *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Remove the presence of the massless shell from the Hessian matrix                        *
 * Parameters | -                                                                                        *
 * Note       | The number of degrees of freedom is the same as for the rigid ion model. The corrected   *
 *            | Hessian matrix where the shell interactions have been removed can be written as:         *
 *            | H'= H_cc - H_cs H_ss^-1 H_sc                                                             *
 *            | cc is core-core, cs is core-shell, sc is shell-core, and ss is the shell-shell part      *
 *********************************************************************************************************/

void RemoveShellInteractionsFromHessian(REAL_MATRIX Hessian,REAL_MATRIX *CorrectedHessian)
{
  int i,j,k,l;
  REAL_MATRIX ss;

  if(ShellSize>0)
  {
    ss=CreateRealMatrix(ShellSize,ShellSize);

    for(i=0;i<ShellSize;i++)
      for(j=0;j<ShellSize;j++)
        ss.element[i][j]=Hessian.element[i+ShellIndex][j+ShellIndex];

    SingularValueDecompositionMatrixInversion(ss);

    for(i=0;i<CoreSize;i++)
    {
      for(j=0;j<CoreSize;j++)
      {
        CorrectedHessian->element[i][j]=Hessian.element[i][j];
        for(k=0;k<ShellSize;k++)
        {
          for (l=0;l<ShellSize;l++)
            CorrectedHessian->element[i][j]-=Hessian.element[i][ShellIndex+l]*ss.element[l][k]*Hessian.element[ShellIndex+k][j];
        }
      }
    }
    DeleteRealMatrix(ss);
  }
  else
  {
    for(i=0;i<Hessian.m;i++)
      for(j=0;j<Hessian.n;j++)
        CorrectedHessian->element[i][j]=Hessian.element[i][j];
  }
}

void VibrationalAnalysis(void)
{
  int i,j;
  REAL Energy,length;
  REAL *Gradients;
  REAL *Weights;
  REAL *Positions;
  REAL *Charges;
  REAL *Frequencies;
  REAL_MATRIX ReducedGeneralizedHessianMatrix;
  REAL_MATRIX GeneralizedHessianMatrix;
  REAL_MATRIX3x3 StrainDerivativeTensor;


  CurrentSystem=0;
  StoredBox=Box[CurrentSystem];
  StoredReplicaBox=ReplicaBox[CurrentSystem];
  StoredInverseBox=InverseBox[CurrentSystem];

  MinimizationVariables=CARTESIAN;
  Ensemble[CurrentSystem]=NVE;

  // index the molecule-atoms into a single list
  NumberOfPositionVariables=OrderNumberOfMinimiationVariables();
  NumberOfBoxVariables=0;   // the box matrix is fixed
  NumberOfVariables=NumberOfPositionVariables+NumberOfBoxVariables;

  GeneralizedHessianMatrix=CreateRealMatrix(CoreSize+ShellSize,CoreSize+ShellSize);
  ReducedGeneralizedHessianMatrix=CreateRealMatrix(CoreSize,CoreSize);

  Gradients=(REAL*)calloc(NumberOfVariables,sizeof(REAL));
  Weights=(REAL*)calloc(NumberOfVariables,sizeof(REAL));
  Positions=(REAL*)calloc(NumberOfVariables,sizeof(REAL));
  Charges=(REAL*)calloc(NumberOfVariables,sizeof(REAL));
  Frequencies=(REAL*)calloc(NumberOfVariables,sizeof(REAL));

  AllocateMinimizationLocalMemory();

  for(i=0;i<NumberOfVariables;i++)
    Gradients[i]=0.0;

  SetWeights(NumberOfCoordinatesMinimizationVariables,Weights,Charges);
  SetStrainToZero(NumberOfCoordinatesMinimizationVariables,Positions);
  CreateGeneralizedCoordinatesFromPositions(NumberOfCoordinatesMinimizationVariables,NumberOfCellMinimizationVariables,Positions);

  // no weighting for the cell matrix
  for(i=0;i<NumberOfBoxVariables;i++)
    Weights[NumberOfPositionVariables+i]=1.0;

  // compute the generalized Hessian
  ComputeDerivativesSpectra(NumberOfPositionVariables,NumberOfBoxVariables,Positions,&Energy,Gradients,GeneralizedHessianMatrix,&StrainDerivativeTensor);

  RemoveShellInteractionsFromHessian(GeneralizedHessianMatrix,&ReducedGeneralizedHessianMatrix);

  MassWeightHessianMatrix(CoreSize,ReducedGeneralizedHessianMatrix,Weights);

  ProjectConstraintsFromHessianMatrixMassWeighted(CoreSize,0,Gradients,ReducedGeneralizedHessianMatrix,Weights);

  SolveEigenValuesAndVectorsHessian(ReducedGeneralizedHessianMatrix,Frequencies);

  // unweight the eigenvectors
  for(i=0;i<CoreSize;i++)
    for(j=0;j<CoreSize;j++)
      ReducedGeneralizedHessianMatrix.element[i][j]*=Weights[i];

  // normalize again
  for(i=0;i<CoreSize;i++)
  {
    length=0;
    for(j=0;j<CoreSize;j++)
      length+=SQR(ReducedGeneralizedHessianMatrix.element[i][j]);
    for(j=0;j<CoreSize;j++)
      ReducedGeneralizedHessianMatrix.element[i][j]/=sqrt(length);
  }

  if(ComputeNormalModes)
    MakeNormalModeMovie(CoreSize,NumberOfBoxVariables,ReducedGeneralizedHessianMatrix,Positions,Frequencies,Weights);

  WriteVibrationalData(ReducedGeneralizedHessianMatrix,Frequencies,Charges);

  free(Frequencies);
  free(Positions);
  free(Charges);
  free(Weights);
  free(Gradients);
  DeleteRealMatrix(GeneralizedHessianMatrix);
}

int RemoveTranslationAndRotationFromHessianMatrix(REAL_MATRIX HessianMatrix,REAL *Weights,REAL *Positions)
{
  int i,j,k,l,n;
  VECTOR com;
  REAL dot_produkt,norm,sum;
  REAL_MATRIX TransRotationMatrix,HV,VHV;

  n=HessianMatrix.m;
  TransRotationMatrix=CreateRealMatrix(6,n);
  HV=CreateRealMatrix(6,n);
  VHV=CreateRealMatrix(6,6);

  CurrentSystem=0;

  com=GetCenterOfMassCurrentSystem();
  com.x=com.y=com.z=0.0;

  for(i=0;i<n/3;i++)
  {
    TransRotationMatrix.element[0][3*i]=1.0/Weights[3*i];
    TransRotationMatrix.element[1][3*i+1]=1.0/Weights[3*i+1];
    TransRotationMatrix.element[2][3*i+2]=1.0/Weights[3*i+2];

    // rotational modes
    TransRotationMatrix.element[3][3*i]=0.0;
    TransRotationMatrix.element[3][3*i+1]=-(1.0/Weights[3*i+1])*(Positions[3*i+2]-com.z);
    TransRotationMatrix.element[3][3*i+2]=(1.0/Weights[3*i+2])*(Positions[3*i+1]-com.y);

    TransRotationMatrix.element[4][3*i]=(1.0/Weights[3*i])*(Positions[3*i+2]-com.z);
    TransRotationMatrix.element[4][3*i+1]=0.0;
    TransRotationMatrix.element[4][3*i+2]=-(1.0/Weights[3*i+2])*(Positions[3*i]-com.x);

    TransRotationMatrix.element[5][3*i]=-(1.0/Weights[3*i])*(Positions[3*i+1]-com.y);
    TransRotationMatrix.element[5][3*i+1]=(1.0/Weights[3*i+1])*(Positions[3*i]-com.x);
    TransRotationMatrix.element[5][3*i+2]=0.0;
  }

  // create orthonormal basis using simple Gram-Schmidt
  NumberOfRemovedModes=0;
  for(i=0;i<6;i++)
  {
    for(j=0;j<NumberOfRemovedModes;j++)
    {
      dot_produkt=0.0;
      for(k=0;k<n;k++)
        dot_produkt+=TransRotationMatrix.element[i][k]*TransRotationMatrix.element[j][k];
      for(k=0;k<n;k++)
        TransRotationMatrix.element[i][k]-=dot_produkt*TransRotationMatrix.element[j][k];
    }
    dot_produkt=0.0;
    for(k=0;k<n;k++)
      dot_produkt+=SQR(TransRotationMatrix.element[i][k]);
    if(dot_produkt>1e-10)
    {
      norm=sqrt(dot_produkt);
      for(k=0;k<n;k++)
        TransRotationMatrix.element[NumberOfRemovedModes][k]=TransRotationMatrix.element[i][k]/norm;
      NumberOfRemovedModes++;
    }
  }

  for(i=0;i<6;i++)
  {
    for(j=0;j<6;j++)
    {
      dot_produkt=0.0;
      for(k=0;k<n;k++)
        dot_produkt+=TransRotationMatrix.element[i][k]*TransRotationMatrix.element[j][k];
    }
  }

  // project out this set of vectors from the Hessian
  // the projection matrix V = I - F*F^T
  // use VHV (H the Hessian) to project these vectors out
  for(j=0;j<NumberOfRemovedModes;j++)
  {
    for(i=0;i<n;i++)
    {
      dot_produkt=0.0;
      for(k=0;k<n;k++)
        dot_produkt+=TransRotationMatrix.element[j][k]*HessianMatrix.element[i][k];
      HV.element[j][i]=dot_produkt;
    }
  }

  for(j=0;j<NumberOfRemovedModes;j++)
  {
    for(i=0;i<NumberOfRemovedModes;i++)
    {
      dot_produkt=0.0;
      for(k=0;k<n;k++)
        dot_produkt+=TransRotationMatrix.element[j][k]*HV.element[i][k];
      VHV.element[i][j]=VHV.element[j][i]=dot_produkt;
    }
  }

  // modify the Hessian
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
    {
      sum=0.0;
      for(k=0;k<NumberOfRemovedModes;k++)
      {
        sum-=HV.element[k][i]*TransRotationMatrix.element[k][j]+
             HV.element[k][j]*TransRotationMatrix.element[k][i];
        for(l=0;l<NumberOfRemovedModes;l++)
          sum+=TransRotationMatrix.element[l][j]*VHV.element[k][l]*TransRotationMatrix.element[k][i];
      }
      HessianMatrix.element[i][j]+=sum;
    }
  DeleteRealMatrix(VHV);
  DeleteRealMatrix(HV);
  DeleteRealMatrix(TransRotationMatrix);

  return NumberOfRemovedModes;
}


int RemoveTranslationFromHessianMatrix(REAL_MATRIX HessianMatrix,REAL *Weights,REAL *Positions)
{
  int i,j,k,l,n;
  REAL dot_produkt,norm,sum;
  REAL_MATRIX TransRotationMatrix,HV,VHV;

  n=HessianMatrix.m;
  TransRotationMatrix=CreateRealMatrix(3,n);
  HV=CreateRealMatrix(3,n);
  VHV=CreateRealMatrix(3,3);

  for(i=0;i<n/3;i++)
  {
    // translational modes
    TransRotationMatrix.element[0][3*i]=1.0/Weights[3*i];
    TransRotationMatrix.element[1][3*i+1]=1.0/Weights[3*i+1];
    TransRotationMatrix.element[2][3*i+2]=1.0/Weights[3*i+2];
  }

  // create orthonormal basis using simple Gram-Schmidt
  NumberOfRemovedModes=0;
  for(i=0;i<3;i++)
  {
    for(j=0;j<NumberOfRemovedModes;j++)
    {
      dot_produkt=0.0;
      for(k=0;k<n;k++)
        dot_produkt+=TransRotationMatrix.element[i][k]*TransRotationMatrix.element[j][k];
      for(k=0;k<n;k++)
        TransRotationMatrix.element[i][k]-=dot_produkt*TransRotationMatrix.element[j][k];
    }
    dot_produkt=0.0;
    for(k=0;k<n;k++)
      dot_produkt+=SQR(TransRotationMatrix.element[i][k]);
    if(dot_produkt>1e-10)
    {
      norm=sqrt(dot_produkt);
      for(k=0;k<n;k++)
        TransRotationMatrix.element[NumberOfRemovedModes][k]=TransRotationMatrix.element[i][k]/norm;
      NumberOfRemovedModes++;
    }
  }

  for(i=0;i<3;i++)
  {
    for(j=0;j<3;j++)
    {
      dot_produkt=0.0;
      for(k=0;k<n;k++)
        dot_produkt+=TransRotationMatrix.element[i][k]*TransRotationMatrix.element[j][k];
    }
  }

  // project out this set of vectors from the Hessian
  // the projection matrix V = I - F*F^T
  // use VHV (H the Hessian) to project these vectors out
  for(j=0;j<NumberOfRemovedModes;j++)
  {
    for(i=0;i<n;i++)
    {
      dot_produkt=0.0;
      for(k=0;k<n;k++)
        dot_produkt+=TransRotationMatrix.element[j][k]*HessianMatrix.element[i][k];
      HV.element[j][i]=dot_produkt;
    }
  }

  for(j=0;j<NumberOfRemovedModes;j++)
  {
    for(i=0;i<NumberOfRemovedModes;i++)
    {
      dot_produkt=0.0;
      for(k=0;k<n;k++)
        dot_produkt+=TransRotationMatrix.element[j][k]*HV.element[i][k];
      VHV.element[i][j]=VHV.element[j][i]=dot_produkt;
    }
  }

  // modify the Hessian
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
    {
      sum=0.0;
      for(k=0;k<NumberOfRemovedModes;k++)
      {
        sum-=HV.element[k][i]*TransRotationMatrix.element[k][j]+
             HV.element[k][j]*TransRotationMatrix.element[k][i];
        for(l=0;l<NumberOfRemovedModes;l++)
          sum+=TransRotationMatrix.element[l][j]*VHV.element[k][l]*TransRotationMatrix.element[k][i];
      }
      HessianMatrix.element[i][j]+=sum;
    }
  DeleteRealMatrix(VHV);
  DeleteRealMatrix(HV);
  DeleteRealMatrix(TransRotationMatrix);

  return NumberOfRemovedModes;
}


void WriteVibrationalData(REAL_MATRIX HessianMatrix,REAL *Eigenvalues,REAL *Charges)
{
  int i,j,k,n,index;
  REAL dot_produkt,BinSize;
  REAL f,freq,temperature;
  REAL dot_produktx,dot_produkty,dot_produktz;
  REAL *spectrum1,*spectrum2,*spectrum3,*spectrum4,*spectrum5,*Frequencies;
  REAL *intensities1,*intensities2,*intensities3,*intensities4,*intensities5;
  REAL *s1,*s2,*s3,*s4,*s5;
  REAL norm_intensities1,norm_intensities2,norm_intensities3,norm_intensities4,norm_intensities5;
  REAL norm_spectrum1,norm_spectrum2,norm_spectrum3,norm_spectrum4,norm_spectrum5;
  FILE *FilePtr;

  if (STREAM)
  {
    fprintf(stderr, "Streaming not yet supported for this function.");
    exit(0);
  }

  n=HessianMatrix.m;
  temperature=therm_baro_stats.ExternalTemperature[0];
  f=PLANCK_CONSTANT*100.0*SPEED_OF_LIGHT/(BOLTZMANN_CONSTANT*temperature);

  LeftFrequencyBoundarySpectrum=10.0;
  RightFrequencyBoundarySpectrum=3700.0;
  SpectrumWidth=800;

  BinSize=(RightFrequencyBoundarySpectrum-LeftFrequencyBoundarySpectrum)/(REAL)SpectrumWidth;

  spectrum1=(REAL*)calloc(SpectrumWidth,sizeof(REAL));
  spectrum2=(REAL*)calloc(SpectrumWidth,sizeof(REAL));
  spectrum3=(REAL*)calloc(SpectrumWidth,sizeof(REAL));
  spectrum4=(REAL*)calloc(SpectrumWidth,sizeof(REAL));
  spectrum5=(REAL*)calloc(SpectrumWidth,sizeof(REAL));

  intensities1=(REAL*)calloc(SpectrumWidth,sizeof(REAL));
  intensities2=(REAL*)calloc(SpectrumWidth,sizeof(REAL));
  intensities3=(REAL*)calloc(SpectrumWidth,sizeof(REAL));
  intensities4=(REAL*)calloc(SpectrumWidth,sizeof(REAL));
  intensities5=(REAL*)calloc(SpectrumWidth,sizeof(REAL));

  Frequencies=(REAL*)calloc(n,sizeof(REAL));
  s1=(REAL*)calloc(n,sizeof(REAL));
  s2=(REAL*)calloc(n,sizeof(REAL));
  s3=(REAL*)calloc(n,sizeof(REAL));
  s4=(REAL*)calloc(n,sizeof(REAL));
  s5=(REAL*)calloc(n,sizeof(REAL));

  for(i=0;i<n;i++)
    Frequencies[i]=SIGN(sqrt(fabs(Eigenvalues[i]))*TO_WAVENUMBERS,Eigenvalues[i]);

  mkdir("VibrationalData",S_IRWXU);
  mkdir("VibrationalData/Spectrum",S_IRWXU);

  // write frequencies (related to eigenvalues)
  FilePtr=fopen("VibrationalData/Spectrum/Frequencies.data","w");
  fprintf(FilePtr,"# %d frequencies (sorted), first %d are zero,\n", 3*n,NumberOfRemovedModes);
  fprintf(FilePtr,"# (projected out and correspond to translation and rotation of the full system)\n");
  fprintf(FilePtr,"# all positive frequencies signals a true local minimum\n");
  fprintf(FilePtr,"# 'n' negative frequencies signals an 'n'-order saddle point\n");
  fprintf(FilePtr,"# number frequency [cm^-1]\n");
  for(i=0;i<n;i++)
    fprintf(FilePtr,"%d %lf\n",i,(double)Frequencies[i]);
  fclose(FilePtr);

  // write normal modes (related to unweighted eigenvectors)
  FilePtr=fopen("VibrationalData/Spectrum/NormalModes.data","w");
   fprintf(FilePtr,"# %d unweighted eigenmodes (sorted), first %d are projected out\n", n,NumberOfRemovedModes);
  for(j=0;j<n;j++)
  {
    for(i=0;i<n;i++)
      fprintf(FilePtr,"%lf ",(double)HessianMatrix.element[i][j]);
    fprintf(FilePtr,"\n");
  }
  fclose(FilePtr);

  for(i=NumberOfRemovedModes;i<n;i++)
  {
    dot_produkt=0.0;
    dot_produktx=0.0;
    dot_produkty=0.0;
    dot_produktz=0.0;
    for(k=0;k<n/3;k++)
      dot_produkt+=HessianMatrix.element[i][3*k]*Charges[k]+
                   HessianMatrix.element[i][3*k+1]*Charges[k]+
                   HessianMatrix.element[i][3*k+2]*Charges[k];
    for(k=0;k<n/3;k++)
    {
      if(fabs(Charges[k])>1e-8)
      {
        dot_produktx+=HessianMatrix.element[i][3*k]*Charges[k];
        dot_produkty+=HessianMatrix.element[i][3*k+1]*Charges[k];
        dot_produktz+=HessianMatrix.element[i][3*k+2]*Charges[k];
      }
      else
      {
        dot_produktx+=HessianMatrix.element[i][3*k];
        dot_produkty+=HessianMatrix.element[i][3*k+1];
        dot_produktz+=HessianMatrix.element[i][3*k+2];
      }
    }
    index=(int)((Frequencies[i]-LeftFrequencyBoundarySpectrum)/BinSize);
    if((index>0)&&(index<SpectrumWidth))
    {
      intensities1[index]+=1.0;
      intensities2[index]+=SQR(dot_produktx);
      intensities3[index]+=SQR(dot_produkty);
      intensities4[index]+=SQR(dot_produktz);
      intensities5[index]+=SQR(dot_produktx)+SQR(dot_produkty)+SQR(dot_produktz);
    }

    if(fabs(Frequencies[i])>1e-6)
      s1[i]=0.5*fabs(Frequencies[i])/(1.0-exp(-f*fabs(Frequencies[i])));
    else
      s1[i]=0.0;

    if((fabs(Frequencies[i])>1e-6)&&(SQR(dot_produktx))>0.0)
      s2[i]=0.5*SQR(dot_produktx)*fabs(Frequencies[i])/(1.0-exp(-f*fabs(Frequencies[i])));
    else
      s2[i]=0.0;

    if((fabs(Frequencies[i])>1e-6)&&(SQR(dot_produkty))>0.0)
      s3[i]=0.5*SQR(dot_produkty)*fabs(Frequencies[i])/(1.0-exp(-f*fabs(Frequencies[i])));
    else
      s3[i]=0.0;

    if((fabs(Frequencies[i])>1e-6)&&(SQR(dot_produktz))>0.0)
      s4[i]=0.5*SQR(dot_produktz)*fabs(Frequencies[i])/(1.0-exp(-f*fabs(Frequencies[i])));
    else
      s4[i]=0.0;

    if((fabs(Frequencies[i])>1e-6)&&(SQR(dot_produktx))>0.0)
      s5[i]=0.5*(SQR(dot_produktx)+SQR(dot_produkty)+SQR(dot_produktz))*fabs(Frequencies[i])/
         (1.0-exp(-f*fabs(Frequencies[i])));
    else
      s5[i]=0.0;
  }

  for(freq=LeftFrequencyBoundarySpectrum;freq<RightFrequencyBoundarySpectrum;freq+=BinSize)
  {
    index=(int)((freq-LeftFrequencyBoundarySpectrum)/BinSize);

    spectrum1[index]=0.0;
    for(j=NumberOfRemovedModes;j<n;j++)
      spectrum1[index]+=s1[j]*exp(-32.0*SQR(freq-fabs(Frequencies[j]))/fabs(Frequencies[j]));

    spectrum2[index]=0.0;
    for(j=NumberOfRemovedModes;j<n;j++)
      spectrum2[index]+=s2[j]*exp(-32.0*SQR(freq-fabs(Frequencies[j]))/fabs(Frequencies[j]));

    spectrum3[index]=0.0;
    for(j=NumberOfRemovedModes;j<n;j++)
      spectrum3[index]+=s3[j]*exp(-32.0*SQR(freq-fabs(Frequencies[j]))/fabs(Frequencies[j]));

    spectrum4[index]=0.0;
    for(j=NumberOfRemovedModes;j<n;j++)
      spectrum4[index]+=s4[j]*exp(-32.0*SQR(freq-fabs(Frequencies[j]))/fabs(Frequencies[j]));

    spectrum5[index]=0.0;
    for(j=NumberOfRemovedModes;j<n;j++)
      spectrum5[index]+=s5[j]*exp(-32.0*SQR(freq-fabs(Frequencies[j]))/fabs(Frequencies[j]));
  }

  norm_spectrum1=norm_spectrum2=norm_spectrum3=norm_spectrum4=norm_spectrum5=1e-10;
  norm_intensities1=norm_intensities2=norm_intensities3=norm_intensities4=norm_intensities5=1e-10;
  for(i=0;i<SpectrumWidth;i++)
  {
    norm_spectrum1+=spectrum1[i];  norm_intensities1+=intensities1[i];
    norm_spectrum2+=spectrum2[i];  norm_intensities2+=intensities2[i];
    norm_spectrum3+=spectrum3[i];  norm_intensities3+=intensities3[i];
    norm_spectrum4+=spectrum4[i];  norm_intensities4+=intensities4[i];
    norm_spectrum5+=spectrum5[i];  norm_intensities5+=intensities5[i];
  }

  FilePtr=fopen("VibrationalData/Spectrum/Spectrum.data","w");
  for(i=0;i<SpectrumWidth;i++)
  {
    freq=(REAL)i*BinSize+(REAL)LeftFrequencyBoundarySpectrum;
    fprintf(FilePtr,"%lf  %lf %lf  %lf %lf  %lf %lf  %lf %lf %lf %lf\n",
           (double)freq,
           (double)(intensities1[i]*BinSize*100.0/norm_intensities1),(double)(spectrum1[i]*BinSize*100.0/norm_spectrum1),
           (double)(intensities2[i]*BinSize*100.0/norm_intensities2),(double)(spectrum2[i]*BinSize*100.0/norm_spectrum2),
           (double)(intensities3[i]*BinSize*100.0/norm_intensities3),(double)(spectrum3[i]*BinSize*100.0/norm_spectrum3),
           (double)(intensities4[i]*BinSize*100.0/norm_intensities4),(double)(spectrum4[i]*BinSize*100.0/norm_spectrum4),
           (double)(intensities5[i]*BinSize*100.0/norm_intensities5),(double)(spectrum5[i]*BinSize*100.0/norm_spectrum5));
  }
  fclose(FilePtr);

  free(spectrum1);
  free(spectrum2);
  free(spectrum3);
  free(spectrum4);
  free(spectrum5);
  free(Frequencies);
  free(s1);
  free(s2);
  free(s3);
  free(s4);
  free(s5);
  free(intensities1);
  free(intensities2);
  free(intensities3);
  free(intensities4);
  free(intensities5);
}

#ifdef HAVE_LAPACK

// zheev

extern void dsyev_(char *jobz, char *uplo, int *n, double *a,
           int *lda, double *w, double *work, int *lwork,
           int *info);

extern void dsyevd_(char *jobz, char *uplo, int *n, double *a,
           int *lda, double *w, double *work, int *lwork,int* iwork,int *liwork,
           int *info);

void SolveEigenValuesAndVectorsHessian(REAL_MATRIX HessianMatrix,REAL* Frequencies)
{
  char jobz, uplo;
  int n,lda, lwork, liwork,info;
  double  *work;
  int *iwork;

  jobz = 'V';
  uplo = 'U';
  n=HessianMatrix.n;
  lda = n; // The leading dimension of the matrix to be solved.

  work=(double *)malloc(1*sizeof(double));
  iwork=(int *)malloc(1*sizeof(int));
  lwork=-1;
  liwork=-1;

  // Workspace size computation
  dsyevd_(&jobz, &uplo, &n, HessianMatrix.element[0], &lda, Frequencies, work, &lwork, iwork,&liwork,&info);
  lwork=(int)work[0];
  liwork=(int)iwork[0];
  free(work);
  free(iwork);
  work=(double *)malloc(lwork*sizeof(double));
  iwork=(int *)malloc(liwork*sizeof(int));

  // DSYEVD for the "real" computation with the optimal block size
  dsyevd_(&jobz, &uplo, &n, HessianMatrix.element[0], &lda, Frequencies, work, &lwork, iwork,&liwork,&info);
  free(work);
}

#else

void SolveEigenValuesAndVectorsHessian(REAL_MATRIX HessianMatrix,REAL* Frequencies)
{
  int i,n;
  REAL *d,*e;

  n=HessianMatrix.m;
  d=(REAL*)calloc(n,sizeof(REAL));
  e=(REAL*)calloc(n,sizeof(REAL));

  tred2(HessianMatrix.element[0],n,d,e);
  tqli(d,e,n,HessianMatrix.element[0]);
  eigsrt(d,HessianMatrix.element[0],n);

  for(i=0;i<n;i++)
    Frequencies[i]=d[i];

  free(d);
  free(e);
}

#endif

void StoredPositionsToRealPositions(int NumberOfPositionVariables,int NumberOfBoxVariables,REAL *Positions,REAL *displacement,REAL_MATRIX HessianMatrix,REAL_MATRIX3x3 StoredBox,REAL_MATRIX3x3 StoredInverseBox)
{
  int i,j,f1,l,m,n,Type,A;
  INT_VECTOR3 index;
  REAL_MATRIX3x3 Transform,RotationMatrix;
  REAL det,TotalMass,Mass,temp,EulerAngle;
  VECTOR r,pos,com,t,p;

  n=NumberOfPositionVariables;

  // Compute total mass in the cell
  TotalMass=0.0;
  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
  {
    for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
      {
        Mass=PseudoAtoms[Framework[CurrentSystem].Atoms[f1][i].Type].Mass;
        TotalMass+=Mass;
      }
  }
  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
    {
      Type=Adsorbates[CurrentSystem][i].Atoms[j].Type;
      Mass=PseudoAtoms[Type].Mass;
      TotalMass+=Mass;
    }
  }
  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
    {
      Type=Cations[CurrentSystem][i].Atoms[j].Type;
      Mass=PseudoAtoms[Type].Mass;
      TotalMass+=Mass;
    }
  }


  switch(Ensemble[CurrentSystem])
  {
    case NPT:
      Transform.ax=1.0+displacement[n]/sqrt(TotalMass);
      Transform.bx=0.0;
      Transform.cx=0.0;

      Transform.ay=0.0;
      Transform.by=1.0+displacement[n+1]/sqrt(TotalMass);
      Transform.cy=0.0;

      Transform.az=0.0;
      Transform.bz=0.0;
      Transform.cz=1.0+displacement[n+2]/sqrt(TotalMass);
      break;
    case NPTPR:
    case NPHPR:
      switch(NPTPRCellType[CurrentSystem])
      {
        case ISOTROPIC:
          temp=(displacement[n]+displacement[n+1]+displacement[n+2])/3.0;
          Transform.ax=1.0+temp; Transform.bx=0.0;      Transform.cx=0.0;
          Transform.ay=0.0;      Transform.by=1.0+temp; Transform.cy=0.0;
          Transform.az=0.0;      Transform.bz=0.0;      Transform.cz=1.0+temp;
        case ANISOTROPIC:
          Transform.ax=1.0+displacement[n];  Transform.bx=0.0;                   Transform.cx=0.0;
          Transform.ay=0.0;                  Transform.by=1.0+displacement[n+1]; Transform.cy=0.0;
          Transform.az=0.0;                  Transform.bz=0.0;                   Transform.cz=1.0+displacement[n+2];
          break;
        case REGULAR:
          Transform.ax=1.0+displacement[n];   Transform.bx=0.5*displacement[n+3]; Transform.cx=0.5*displacement[n+6];
          Transform.ay=0.5*displacement[n+1]; Transform.by=1.0+displacement[n+4]; Transform.cy=0.5*displacement[n+7];
          Transform.az=0.5*displacement[n+2]; Transform.bz=0.5*displacement[n+5]; Transform.cz=1.0+displacement[n+8];
          break;
        case REGULAR_UPPER_TRIANGLE:
          Transform.ax=1.0+displacement[n];   Transform.bx=0.5*displacement[n+1]; Transform.cx=0.5*displacement[n+2];
          Transform.ay=0.5*displacement[n+1]; Transform.by=1.0+displacement[n+3]; Transform.cy=0.5*displacement[n+4];
          Transform.az=0.5*displacement[n+2]; Transform.bz=0.5*displacement[n+4]; Transform.cz=1.0+displacement[n+5];
          break;
        case MONOCLINIC:
          Transform.ax=1.0+displacement[n];   Transform.bx=0;                     Transform.cx=0.5*displacement[n+3];
          Transform.ay=0.0;                   Transform.by=1.0+displacement[n+2]; Transform.cy=0;
          Transform.az=0.5*displacement[n+1]; Transform.bz=0.0;                   Transform.cz=1.0+displacement[n+4];
          break;
        case MONOCLINIC_UPPER_TRIANGLE:
          Transform.ax=1.0+displacement[n];   Transform.bx=0;                     Transform.cx=displacement[n+2];
          Transform.ay=0.0;                   Transform.by=1.0+displacement[n+1]; Transform.cy=0;
          Transform.az=0.0;                   Transform.bz=0.0;                   Transform.cz=1.0+displacement[n+3];
          break;
        default:
          fprintf(stderr, "Unknown NPTPRCellType\n");
          exit(0);
          break;
      }
      break;
    case NVT:
    case NVE:
    default:
      Transform.ax=1.0;  Transform.bx=0.0; Transform.cx=0.0;
      Transform.ay=0.0;  Transform.by=1.0; Transform.cy=0.0;
      Transform.az=0.0;  Transform.bz=0.0; Transform.cz=1.0;
      break;
  }

  Box[CurrentSystem]=MatrixMatrixMultiplication3x3(Transform,StoredBox);
  CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);
  Invert3x3Matrix(&Box[CurrentSystem],&InverseBox[CurrentSystem],&det);
  CellProperties(&InverseBox[CurrentSystem],&InverseBoxProperties[CurrentSystem],&det);

  AlphaAngle[CurrentSystem]=acos(BoxProperties[CurrentSystem].bx);
  BetaAngle[CurrentSystem]=acos(BoxProperties[CurrentSystem].by);
  GammaAngle[CurrentSystem]=acos(BoxProperties[CurrentSystem].bz);

  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
  {
    for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
      {
        index=Framework[CurrentSystem].Atoms[f1][i].HessianIndex;
        Mass=PseudoAtoms[Framework[CurrentSystem].Atoms[f1][i].Type].Mass;
        r.x=Positions[index.x]+displacement[index.x]/sqrt(Mass);
        r.y=Positions[index.y]+displacement[index.y]/sqrt(Mass);
        r.z=Positions[index.z]+displacement[index.z]/sqrt(Mass);
        Framework[CurrentSystem].Atoms[f1][i].Position.x=Transform.ax*r.x+Transform.bx*r.y+Transform.cx*r.z;
        Framework[CurrentSystem].Atoms[f1][i].Position.y=Transform.ay*r.x+Transform.by*r.y+Transform.cy*r.z;
        Framework[CurrentSystem].Atoms[f1][i].Position.z=Transform.az*r.x+Transform.bz*r.y+Transform.cz*r.z;
      }
    }
  }

  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    Type=Adsorbates[CurrentSystem][m].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      index=Adsorbates[CurrentSystem][m].Groups[l].HessianIndex;
      if(Components[Type].Groups[l].Rigid) // rigid unit
      {
        Mass=Components[Type].Groups[l].Mass;
        pos.x=Positions[index.x]+displacement[index.x]/sqrt(Mass);
        pos.y=Positions[index.y]+displacement[index.y]/sqrt(Mass);
        pos.z=Positions[index.z]+displacement[index.z]/sqrt(Mass);
        Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassPosition.x=Transform.ax*pos.x+Transform.bx*pos.y+Transform.cx*pos.z;
        Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassPosition.y=Transform.ay*pos.x+Transform.by*pos.y+Transform.cy*pos.z;
        Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassPosition.z=Transform.az*pos.x+Transform.bz*pos.y+Transform.cz*pos.z;
        com=Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassPosition;

        Adsorbates[CurrentSystem][m].Groups[l].EulerAxis.x=Positions[index.x+3]+displacement[index.x+3]*sqrt(Components[Type].Groups[l].InverseInertiaVector.x);
        Adsorbates[CurrentSystem][m].Groups[l].EulerAxis.y=Positions[index.y+3]+displacement[index.y+3]*sqrt(Components[Type].Groups[l].InverseInertiaVector.y);
        Adsorbates[CurrentSystem][m].Groups[l].EulerAxis.z=Positions[index.z+3]+displacement[index.z+3]*sqrt(Components[Type].Groups[l].InverseInertiaVector.z);
        p=Adsorbates[CurrentSystem][m].Groups[l].EulerAxis;
        EulerAngle=sqrt(SQR(p.x)+SQR(p.y)+SQR(p.z));

        if(EulerAngle<1e-8)
        {
          RotationMatrix.ax=1.0; RotationMatrix.bx=0.0; RotationMatrix.cx=0.0;
          RotationMatrix.ay=0.0; RotationMatrix.by=1.0; RotationMatrix.cy=0.0;
          RotationMatrix.az=0.0; RotationMatrix.bz=0.0; RotationMatrix.cz=1.0;
        }
        else
        {
          RotationMatrix.ax=1.0+(SQR(p.y)+SQR(p.z))*(cos(EulerAngle)-1.0)/SQR(EulerAngle);
          RotationMatrix.ay=(p.x*p.y - p.x*p.y*cos(EulerAngle) + p.z*EulerAngle*sin(EulerAngle))/SQR(EulerAngle);
          RotationMatrix.az=(p.x*p.z - p.x*p.z*cos(EulerAngle) - p.y*EulerAngle*sin(EulerAngle))/SQR(EulerAngle);

          RotationMatrix.bx=(p.x*p.y - p.x*p.y*cos(EulerAngle) - p.z*EulerAngle*sin(EulerAngle))/SQR(EulerAngle);
          RotationMatrix.by=1.0+((SQR(p.x) + SQR(p.z))*(-1 + cos(EulerAngle)))/SQR(EulerAngle);
          RotationMatrix.bz=(p.y*p.z - p.y*p.z*cos(EulerAngle) + p.x*EulerAngle*sin(EulerAngle))/SQR(EulerAngle);

          RotationMatrix.cx=(p.x*p.z - p.x*p.z*cos(EulerAngle) + p.y*EulerAngle*sin(EulerAngle))/SQR(EulerAngle);
          RotationMatrix.cy=(p.y*p.z - p.y*p.z*cos(EulerAngle) - p.x*EulerAngle*sin(EulerAngle))/SQR(EulerAngle);
          RotationMatrix.cz=1.0+((SQR(p.x) + SQR(p.y))*(-1.0 + cos(EulerAngle)))/SQR(EulerAngle);
        }

       for(i=0;i<Components[Type].Groups[l].NumberOfGroupAtoms;i++)
        {
          A=Components[Type].Groups[l].Atoms[i];
          pos=Components[Type].Positions[A];

          t.x=RotationMatrix.ax*pos.x+RotationMatrix.bx*pos.y+RotationMatrix.cx*pos.z;
          t.y=RotationMatrix.ay*pos.x+RotationMatrix.by*pos.y+RotationMatrix.cy*pos.z;
          t.z=RotationMatrix.az*pos.x+RotationMatrix.bz*pos.y+RotationMatrix.cz*pos.z;
          pos.x=com.x+t.x;
          pos.y=com.y+t.y;
          pos.z=com.z+t.z;
          Adsorbates[CurrentSystem][m].Atoms[A].Position=pos;
        }
      }
      else // flexible unit
      {
        for(i=0;i<Components[Type].Groups[l].NumberOfGroupAtoms;i++)
        {
          A=Components[Type].Groups[l].Atoms[i];
          index=Adsorbates[CurrentSystem][m].Atoms[A].HessianIndex;
          Mass=PseudoAtoms[Adsorbates[CurrentSystem][m].Atoms[A].Type].Mass;

          r.x=Positions[index.x]+displacement[index.x]/sqrt(Mass);
          r.y=Positions[index.y]+displacement[index.y]/sqrt(Mass);
          r.z=Positions[index.z]+displacement[index.z]/sqrt(Mass);
          Adsorbates[CurrentSystem][m].Atoms[A].Position.x=Transform.ax*r.x+Transform.bx*r.y+Transform.cx*r.z;
          Adsorbates[CurrentSystem][m].Atoms[A].Position.y=Transform.ay*r.x+Transform.by*r.y+Transform.cy*r.z;
          Adsorbates[CurrentSystem][m].Atoms[A].Position.z=Transform.az*r.x+Transform.bz*r.y+Transform.cz*r.z;
        }
      }
    }
  }

  for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
  {
    Type=Cations[CurrentSystem][m].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      index=Cations[CurrentSystem][m].Groups[l].HessianIndex;
      if(Components[Type].Groups[l].Rigid) // rigid unit
      {
        Mass=Components[Type].Groups[l].Mass;
        pos.x=Positions[index.x]+displacement[index.x]/sqrt(Mass);
        pos.y=Positions[index.y]+displacement[index.y]/sqrt(Mass);
        pos.z=Positions[index.z]+displacement[index.z]/sqrt(Mass);
        Cations[CurrentSystem][m].Groups[l].CenterOfMassPosition.x=Transform.ax*pos.x+Transform.bx*pos.y+Transform.cx*pos.z;
        Cations[CurrentSystem][m].Groups[l].CenterOfMassPosition.y=Transform.ay*pos.x+Transform.by*pos.y+Transform.cy*pos.z;
        Cations[CurrentSystem][m].Groups[l].CenterOfMassPosition.z=Transform.az*pos.x+Transform.bz*pos.y+Transform.cz*pos.z;
        com=Cations[CurrentSystem][m].Groups[l].CenterOfMassPosition;

        Cations[CurrentSystem][m].Groups[l].EulerAxis.x=Positions[index.x+3]+displacement[index.x+3]*sqrt(Components[Type].Groups[l].InverseInertiaVector.x);
        Cations[CurrentSystem][m].Groups[l].EulerAxis.y=Positions[index.y+3]+displacement[index.y+3]*sqrt(Components[Type].Groups[l].InverseInertiaVector.y);
        Cations[CurrentSystem][m].Groups[l].EulerAxis.z=Positions[index.z+3]+displacement[index.z+3]*sqrt(Components[Type].Groups[l].InverseInertiaVector.z);
        p=Cations[CurrentSystem][m].Groups[l].EulerAxis;
        EulerAngle=sqrt(SQR(p.x)+SQR(p.y)+SQR(p.z));

        if(EulerAngle<1e-8)
        {
          RotationMatrix.ax=1.0; RotationMatrix.bx=0.0; RotationMatrix.cx=0.0;
          RotationMatrix.ay=0.0; RotationMatrix.by=1.0; RotationMatrix.cy=0.0;
          RotationMatrix.az=0.0; RotationMatrix.bz=0.0; RotationMatrix.cz=1.0;
        }
        else
        {
          RotationMatrix.ax=1.0+(SQR(p.y)+SQR(p.z))*(cos(EulerAngle)-1.0)/SQR(EulerAngle);
          RotationMatrix.ay=(p.x*p.y - p.x*p.y*cos(EulerAngle) + p.z*EulerAngle*sin(EulerAngle))/SQR(EulerAngle);
          RotationMatrix.az=(p.x*p.z - p.x*p.z*cos(EulerAngle) - p.y*EulerAngle*sin(EulerAngle))/SQR(EulerAngle);

          RotationMatrix.bx=(p.x*p.y - p.x*p.y*cos(EulerAngle) - p.z*EulerAngle*sin(EulerAngle))/SQR(EulerAngle);
          RotationMatrix.by=1.0+((SQR(p.x) + SQR(p.z))*(-1 + cos(EulerAngle)))/SQR(EulerAngle);
          RotationMatrix.bz=(p.y*p.z - p.y*p.z*cos(EulerAngle) + p.x*EulerAngle*sin(EulerAngle))/SQR(EulerAngle);

          RotationMatrix.cx=(p.x*p.z - p.x*p.z*cos(EulerAngle) + p.y*EulerAngle*sin(EulerAngle))/SQR(EulerAngle);
          RotationMatrix.cy=(p.y*p.z - p.y*p.z*cos(EulerAngle) - p.x*EulerAngle*sin(EulerAngle))/SQR(EulerAngle);
          RotationMatrix.cz=1.0+((SQR(p.x) + SQR(p.y))*(-1.0 + cos(EulerAngle)))/SQR(EulerAngle);
        }

       for(i=0;i<Components[Type].Groups[l].NumberOfGroupAtoms;i++)
        {
          A=Components[Type].Groups[l].Atoms[i];
          pos=Components[Type].Positions[A];

          t.x=RotationMatrix.ax*pos.x+RotationMatrix.bx*pos.y+RotationMatrix.cx*pos.z;
          t.y=RotationMatrix.ay*pos.x+RotationMatrix.by*pos.y+RotationMatrix.cy*pos.z;
          t.z=RotationMatrix.az*pos.x+RotationMatrix.bz*pos.y+RotationMatrix.cz*pos.z;
          pos.x=com.x+t.x;
          pos.y=com.y+t.y;
          pos.z=com.z+t.z;
          Cations[CurrentSystem][m].Atoms[A].Position=pos;
        }
      }
      else // flexible unit
      {
        for(i=0;i<Components[Type].Groups[l].NumberOfGroupAtoms;i++)
        {
          A=Components[Type].Groups[l].Atoms[i];
          index=Cations[CurrentSystem][m].Atoms[A].HessianIndex;
          Mass=PseudoAtoms[Cations[CurrentSystem][m].Atoms[A].Type].Mass;

          r.x=Positions[index.x]+displacement[index.x]/sqrt(Mass);
          r.y=Positions[index.y]+displacement[index.y]/sqrt(Mass);
          r.z=Positions[index.z]+displacement[index.z]/sqrt(Mass);
          Cations[CurrentSystem][m].Atoms[A].Position.x=Transform.ax*r.x+Transform.bx*r.y+Transform.cx*r.z;
          Cations[CurrentSystem][m].Atoms[A].Position.y=Transform.ay*r.x+Transform.by*r.y+Transform.cy*r.z;
          Cations[CurrentSystem][m].Atoms[A].Position.z=Transform.az*r.x+Transform.bz*r.y+Transform.cz*r.z;
        }
      }
    }
  }

}

void WriteIR(void)
{
  int i,j,k;
  REAL Charge;
  VECTOR vel,cor;
  FILE *FilePtr;
  char buffer[256];

  if (STREAM)
  {
    fprintf(stderr, "Streaming not yet supported for this function.");
    exit(0);
  }

  mkdir("InfraRedSpectra",S_IRWXU);
  for(k=0;k<NumberOfSystems;k++)
  {
    sprintf(buffer,"InfraRedSpectra/System_%d",k);
    mkdir(buffer,S_IRWXU);

    cor.x=0.0;
    cor.y=0.0;
    cor.z=0.0;
    for(i=0;i<Framework[CurrentSystem].TotalNumberOfAtoms;i++)
    {
      vel=Framework[CurrentSystem].Atoms[0][i].Velocity;
      Charge=Framework[CurrentSystem].Atoms[0][i].Charge;
      cor.x+=Charge*vel.x;
      cor.y+=Charge*vel.y;
      cor.z+=Charge*vel.z;
    }
    for(i=0;i<NumberOfAdsorbateMolecules[k];i++)
    {
      for(j=0;j<Adsorbates[k][i].NumberOfAtoms;j++)
      {
        vel=Adsorbates[k][i].Atoms[j].Velocity;
        Charge=Adsorbates[k][i].Atoms[j].Charge;
        cor.x+=Charge*vel.x;
        cor.y+=Charge*vel.y;
        cor.z+=Charge*vel.z;
      }
    }
    for(i=0;i<NumberOfCationMolecules[k];i++)
    {
      for(j=0;j<Cations[k][i].NumberOfAtoms;j++)
      {
        vel=Cations[k][i].Atoms[j].Velocity;
        Charge=Cations[k][i].Atoms[j].Charge;
        cor.x+=Charge*vel.x;
        cor.y+=Charge*vel.y;
        cor.z+=Charge*vel.z;
      }
    }

    sprintf(buffer,"InfraRedSpectra/System_%d/IR%s.dat",k,FileNameAppend);
    FilePtr=fopen(buffer,"a");
    fprintf(FilePtr,"%g 0.0\n",(double)((cor.x+cor.y+cor.z)/3.0));
    fclose(FilePtr);
  }
}

void MakeNormalModeMovie(int NumberOfPositionVariables,int NumberOfBoxVariables,REAL_MATRIX HessianMatrix,REAL *Positions,REAL *frequencies,REAL* Weights)
{
  int i,j,f1,n,mode,frames,Type;
  REAL frequency,amplitude,factor,temperature;
  FILE *FilePtr;
  char RecordName[7]="ATOM  ";          // ATOM
  char SerialNumber[6]="     ";         // Atom serial number
  char AtomName[5]="   C";              // Atom Name
  char AltLoc[2]=" ";                   // Alternate location indicator
  char ResIdueName[4]="MOL";            // ResIdue Name
  char ChainId[2]=" ";                  // Chain Identifier
  char ResSeq[5]="   2";                // ResIdue sequence number
  char iCode[2]=" ";                    // code for insertion of resIdues
  REAL Occupancy=1.0;                   // Occupancy
  REAL Temp=0.0;                        // Temperature factor
  char SegID[5]="    ";                 // Segment Identifier, left-justified
  char Element[3]="  ";                 // Element symbol, right-justified
  char Charge[3]="  ";                  // Charge on the Atom
  char buffer[256];
  VECTOR pos;
  REAL_MATRIX3x3 StoredBox,StoredInverseBox;
  REAL *displacement;
  REAL StoredVolume,MaximumVolumeChange,VolumeChange;
  const REAL scale=1.0;
  INT_VECTOR3 index;

  if (STREAM)
  {
    fprintf(stderr, "Streaming not yet supported for this function.");
    exit(0);
  }

  n=NumberOfPositionVariables+NumberOfBoxVariables;

  CurrentSystem=0;
  StoredBox=Box[0];
  StoredInverseBox=InverseBox[0];
  CellProperties(&Box[0],&BoxProperties[0],&StoredVolume);

  displacement=(REAL*)calloc(n,sizeof(REAL));

  mkdir("VibrationalData",S_IRWXU);
  mkdir("VibrationalData/Modes",S_IRWXU);

  MaximumMode=MIN2(MaximumMode,n);

  for(mode=MinimumMode;mode<MaximumMode;mode++)
  {
    frequency=frequencies[mode];
    if(fabs(frequency)>1e-2)
    {
      temperature=therm_baro_stats.ExternalTemperature[0];

      // compute amplitude in Angstrom
      amplitude=sqrt(2.0*K_B*temperature/fabs(frequency));

      sprintf(buffer,"VibrationalData/Modes/mode_%d%s.pdb",mode,FileNameAppend);
      FilePtr=fopen(buffer,"w");

      MaximumVolumeChange=0;
      for(frames=0;frames<ModeResolution;frames++)
      {
        factor=amplitude*sin(2.0*M_PI*frames/ModeResolution);

        for(i=0;i<n;i++)
          displacement[i]=factor*HessianMatrix.element[mode][i];

        StoredPositionsToRealPositions(NumberOfPositionVariables,NumberOfBoxVariables,Positions,displacement,HessianMatrix,StoredBox,StoredInverseBox);

        // correct the positions for constraints
        if(CorrectNormalModesForConstraints)
          ShakeInMinimization();


        VolumeChange=fabs(Volume[0]/StoredVolume)*100.0;
        if(VolumeChange>MaximumVolumeChange) MaximumVolumeChange=VolumeChange;
        fprintf(FilePtr,"MODEL%9d\n",frames);
        fprintf(FilePtr,"REMARK   Raspa-1.0 PDB file\n");
        fprintf(FilePtr,"CRYST1%9.3lf%9.3lf%9.3lf%7.2lf%7.2lf%7.2lf\n",
            (double)(scale*BoxProperties[CurrentSystem].ax),(double)(scale*BoxProperties[CurrentSystem].ay),(double)(scale*BoxProperties[CurrentSystem].az),
            (double)AlphaAngle[CurrentSystem]*RAD2DEG,(double)BetaAngle[CurrentSystem]*RAD2DEG,(double)GammaAngle[CurrentSystem]*RAD2DEG);

        for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
        {
          for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
          {
            pos=Framework[CurrentSystem].Atoms[f1][i].Position;
            sprintf(AtomName,"%2s",PseudoAtoms[Framework[CurrentSystem].Atoms[f1][i].Type].PrintToPDBName);
            sprintf(Element,"%2s",PseudoAtoms[Framework[CurrentSystem].Atoms[f1][i].Type].PrintToPDBName);
            fprintf(FilePtr,"%s%-s %4s%s%s %s%s%s   %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf      %s%s%s\n",
                      RecordName,SerialNumber,AtomName,AltLoc,ResIdueName,
                      ChainId,ResSeq,iCode,(double)(scale*pos.x),(double)(scale*pos.y),(double)(scale*pos.z),(double)Occupancy,(double)Temp,SegID,Element,Charge);
          }
        }
        for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
        {
          for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
          {
            Type=Adsorbates[CurrentSystem][i].Atoms[j].Type;
            if(PseudoAtoms[Type].PrintToPDB)
            {
              pos=Adsorbates[CurrentSystem][i].Atoms[j].Position;
              index=Adsorbates[CurrentSystem][i].Atoms[j].HessianIndex;

              sprintf(AtomName,"%2s",PseudoAtoms[Type].PrintToPDBName);
              sprintf(Element,"%2s",PseudoAtoms[Type].PrintToPDBName);
              fprintf(FilePtr,"%s%-s %4s%s%s %s%s%s   %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf      %s%s%s\n",
                      RecordName,SerialNumber,AtomName,AltLoc,ResIdueName,
                      ChainId,ResSeq,iCode,(double)(scale*pos.x),(double)(scale*pos.y),(double)(scale*pos.z),(double)Occupancy,(double)Temp,SegID,Element,Charge);
            }
          }
        }
        fprintf(FilePtr,"ENDMDL\n");
      }
      fclose(FilePtr);
    }
  }
  free(displacement);
}

int HeapSortMillerPlanes(void)
{
  unsigned long i,ir,j,l;
  POWDER_DIFF_PEAK rra;
  int n;

  n=NumberOfPeaks;
  if(n<2) return 0;
  l=(n>>1)+1;
  ir=n;
  for(;;)
  {
    if(l>1)
      rra=PowderDiffractionPeaks[--l-1];
    else
    {
      rra=PowderDiffractionPeaks[ir-1];
      PowderDiffractionPeaks[ir-1]=PowderDiffractionPeaks[0];
      if(--ir==1)
      {
        PowderDiffractionPeaks[0]=rra;
        break;
      }
    }
    i=l;
    j=l+l;
    while(j<=ir)
    {
      if(j<ir&&PowderDiffractionPeaks[j-1].d<PowderDiffractionPeaks[j].d) j++;
      if(rra.d<PowderDiffractionPeaks[j-1].d)
      {
        PowderDiffractionPeaks[i-1]=PowderDiffractionPeaks[j-1];
        i=j;
        j <<= 1;
      }
      else
        break;
    }
    PowderDiffractionPeaks[i-1]=rra;
  }
  return 0;
}

void FindUniqueMillerPlanes(void)
{
  int i;
  int h,k,l,hmax,kmax,lmax;
  VECTOR a,b,c,k1,k2,k3,dvec;
  REAL sinlambda,sintheta,theta,two_theta,delta_theta;
  REAL d,Lp,scattering_factor,Intensity;
  REAL scat,s,B;
  int PeakIsPresent,Type;
  VECTOR pos;
  REAL max,cs,sn;
  COMPLEX anom_scat,part;

  a.x=Box[CurrentSystem].ax; a.y=Box[CurrentSystem].ay; a.z=Box[CurrentSystem].az;
  b.x=Box[CurrentSystem].bx; b.y=Box[CurrentSystem].by; b.z=Box[CurrentSystem].bz;
  c.x=Box[CurrentSystem].cx; c.y=Box[CurrentSystem].cy; c.z=Box[CurrentSystem].cz;

  k1=CrossProduct(b,c);
  k1.x/=Volume[CurrentSystem]; k1.y/=Volume[CurrentSystem]; k1.z/=Volume[CurrentSystem];

  k2=CrossProduct(c,a);
  k2.x/=Volume[CurrentSystem]; k2.y/=Volume[CurrentSystem]; k2.z/=Volume[CurrentSystem];

  k3=CrossProduct(a,b);
  k3.x/=Volume[CurrentSystem]; k3.y/=Volume[CurrentSystem]; k3.z/=Volume[CurrentSystem];

  hmax=Diffraction.hmax;
  kmax=Diffraction.kmax;
  lmax=Diffraction.lmax;

  NumberOfPeaks=0;
  for(h=-hmax;h<=hmax;h++)
    for(k=-kmax;k<=kmax;k++)
      for(l=-lmax;l<=lmax;l++)
      {
        dvec.x=h*k1.x+k*k2.x+l*k3.x;
        dvec.y=h*k1.y+k*k2.y+l*k3.y;
        dvec.z=h*k1.z+k*k2.z+l*k3.z;

        d=sqrt(SQR(dvec.x)+SQR(dvec.y)+SQR(dvec.z));

        sinlambda=0.5*d;
        sintheta=0.5*Diffraction.lambda*d;

        if(sintheta<=1.0)
          theta=asin(sintheta);
        else theta=100.0;

        two_theta=2.0*theta;
        delta_theta=2.0*asin(0.5*Diffraction.lambda2*d);
        s=sinlambda;

        part.re=part.im=0.0;
        for(CurrentFramework=0;CurrentFramework<Framework[CurrentSystem].NumberOfFrameworks;CurrentFramework++)
        {
          for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];i++)
          {
            pos=ConvertFromXYZtoABC(Framework[CurrentSystem].Atoms[CurrentFramework][i].Position);
            Type=Framework[CurrentSystem].Atoms[CurrentFramework][i].Type;
            scat=ScatteringFactor(PseudoAtoms[Type].ScatteringType,SQR(s));
            if(Diffraction.Type==XRAY_DIFFRACTION)
              anom_scat=GetAnomalousScatteringFactor(PseudoAtoms[Type].AnomalousScatteringType);
            else
              anom_scat.re=anom_scat.im=0.0;
            B=exp(-PseudoAtoms[Type].TemperatureFactor*SQR(s));
            cs=cos(2.0*M_PI*(h*pos.x+k*pos.y+l*pos.z));
            sn=sin(2.0*M_PI*(h*pos.x+k*pos.y+l*pos.z));
            part.re+=B*((scat+anom_scat.re)*cs-anom_scat.im*sn);
            part.im+=B*((scat+anom_scat.re)*sn+anom_scat.im*cs);
          }
        }
        scattering_factor=SQR(part.re)+SQR(part.im);
        if(Diffraction.Type==XRAY_DIFFRACTION)
          Lp=(1.0+SQR(cos(two_theta)))/(SQR(sin(theta))*cos(theta));
        else
          Lp=1.0/(sin(two_theta)*sin(theta));

        Intensity=Lp*scattering_factor;

        PeakIsPresent=FALSE;
        for(i=0;i<NumberOfPeaks;i++)
        {
          // identical plane when |F_hkl|^2 and d_hkl are the same
          if((fabs(PowderDiffractionPeaks[i].ScatteringFactor-scattering_factor)<1e-4)&&
             (fabs(PowderDiffractionPeaks[i].d-d)<1e-4))
          {
            // plane already present, so raise multiplcity and store convienent hkl-values
            // that is, all positive and h>=k and k>=l
            PeakIsPresent=TRUE;
            PowderDiffractionPeaks[i].Multiplicity++;
            if((h>=0)&&(l>=0)&&(l>=0)&&(h>=k)&&(k>=l))
            {
              PowderDiffractionPeaks[i].h=h;
              PowderDiffractionPeaks[i].k=k;
              PowderDiffractionPeaks[i].l=l;
            }
            continue;
          }
        }
        // peak was not already stored, so store peak-data now but skip when h=k=l=0
        if((PeakIsPresent==FALSE)&&(!((h==0)&&(k==0)&&(l==0))))
        {
          PowderDiffractionPeaks[NumberOfPeaks].Multiplicity=1;
          PowderDiffractionPeaks[NumberOfPeaks].two_theta=two_theta;
          PowderDiffractionPeaks[NumberOfPeaks].sintheta=0.5*Diffraction.lambda*d;
          PowderDiffractionPeaks[NumberOfPeaks].tantheta=0.5*Diffraction.lambda*d/cos(theta);
          PowderDiffractionPeaks[NumberOfPeaks].d=d;
          PowderDiffractionPeaks[NumberOfPeaks].s=0.5*d;
          PowderDiffractionPeaks[NumberOfPeaks].h=h;
          PowderDiffractionPeaks[NumberOfPeaks].k=k;
          PowderDiffractionPeaks[NumberOfPeaks].l=l;
          PowderDiffractionPeaks[NumberOfPeaks].Lp=Lp;
          PowderDiffractionPeaks[NumberOfPeaks].ScatteringFactor=scattering_factor;
          NumberOfPeaks++;
          if(NumberOfPeaks>=MAX_POWDER_DIFFRACTION_PEAKS)
          {
            fprintf(stderr, "Maximum number of pxrd-peaks to large\n");
            exit(0);
          }
        }
      }

  // calculate full intensity taking multiplicity into account and also find maximum value
  max=0.0;
  for(i=0;i<NumberOfPeaks;i++)
  {
    Intensity=PowderDiffractionPeaks[i].Multiplicity*
              PowderDiffractionPeaks[i].Lp*
              PowderDiffractionPeaks[i].ScatteringFactor;
    PowderDiffractionPeaks[i].Intensity=Intensity;
    if(Intensity>max) max=Intensity;
  }

  // normalize maximum peak height to 100
  for(i=0;i<NumberOfPeaks;i++)
    PowderDiffractionPeaks[i].Intensity/=0.01*max;
}

void CalculatePowderDiffractionPattern(void)
{
  int i,j,n;
  int cur_peak,old_peak;
  REAL two_theta,theta,tantheta,dr,bf,fwhm;
  char buffer[256];
  FILE *FilePtr;

  if (STREAM)
  {
    fprintf(stderr, "Streaming not yet supported for this function.");
    exit(0);
  }

  n=1+(Diffraction.two_theta_max-Diffraction.two_theta_min)/Diffraction.two_theta_step;
  Diffraction.n=n;

  Diffraction.spectrum=(REAL*)calloc(n,sizeof(REAL));

  for(i=0;i<NumberOfPeaks;i++)
  {
    two_theta=PowderDiffractionPeaks[i].two_theta;
    tantheta=PowderDiffractionPeaks[i].tantheta;
    fwhm=sqrt(Diffraction.w+Diffraction.v*tantheta+Diffraction.u*SQR(tantheta));
    theta=0.5*two_theta;
    dr=(two_theta*RAD2DEG-Diffraction.two_theta_min)/Diffraction.two_theta_step;
    cur_peak=NINT(dr);
    old_peak=cur_peak;

    for(j=0;j<n;j++)
    {
      dr=two_theta*RAD2DEG-Diffraction.two_theta_step*j-Diffraction.two_theta_min;

      dr=(j-cur_peak)*Diffraction.two_theta_step;

      bf=1.0/Diffraction.two_theta_step;

      if(fabs(fwhm)>0.00001)
      {
        switch(Diffraction.PeakShape)
        {
          case DIFFRACTION_GAUSSIAN:
            bf=sqrt(4.0*log(2.0))*exp(-4.0*log(2.0)*SQR(dr)/SQR(fwhm))/(fwhm*sqrt(M_PI));
            break;
          case DIFFRACTION_LORENTZIAN:
            bf=fwhm/(2.0*M_PI*(SQR(dr)+0.25*SQR(fwhm)));
            break;
          case DIFFRACTION_PSEUDO_VOIGT:
            bf=Diffraction.asym*2.0/(fwhm*M_PI*(1.0+4.0*(SQR(dr)/(SQR(fwhm)))))+
               (1.0-Diffraction.asym)*sqrt(4.0*log(2.0))*exp(-4.0*log(2.0)*SQR(dr)/SQR(fwhm))/(fwhm*M_PI);
            break;
        }
      }
      else
      {
        if(cur_peak!=j)
          bf=0.0;
      }

      Diffraction.spectrum[j]+=bf*PowderDiffractionPeaks[i].Intensity;
    }
  }

  mkdir("PowderDiffraction",S_IRWXU);
  sprintf(buffer,"PowderDiffraction/System_%d",CurrentSystem);
  mkdir(buffer,S_IRWXU);

  sprintf(buffer,"PowderDiffraction/System_%d/Spectrum_%s_%d.%d.%d_%lf_%lf%s.dat",CurrentSystem,
                Framework[CurrentSystem].Name[0],
                NumberOfUnitCells[CurrentSystem].x,
                NumberOfUnitCells[CurrentSystem].y,
                NumberOfUnitCells[CurrentSystem].z,
                (double)therm_baro_stats.ExternalTemperature[CurrentSystem],
                (double)(therm_baro_stats.ExternalPressure[CurrentSystem][CurrentIsothermPressure]*PRESSURE_CONVERSION_FACTOR),
                FileNameAppend);
  FilePtr=fopen(buffer,"w");

  n=Diffraction.n;
  for(j=0;j<n;j++)
    fprintf(FilePtr,"%g %g\n",Diffraction.two_theta_step*j+Diffraction.two_theta_min,Diffraction.spectrum[j]);
  fclose(FilePtr);

  sprintf(buffer,"PowderDiffraction/System_%d/PeakInformation_%s_%d.%d.%d_%lf_%lf%s.dat",CurrentSystem,
                Framework[CurrentSystem].Name[0],
                NumberOfUnitCells[CurrentSystem].x,
                NumberOfUnitCells[CurrentSystem].y,
                NumberOfUnitCells[CurrentSystem].z,
                (double)therm_baro_stats.ExternalTemperature[CurrentSystem],
                (double)(therm_baro_stats.ExternalPressure[CurrentSystem][CurrentIsothermPressure]*PRESSURE_CONVERSION_FACTOR),
                FileNameAppend);
  FilePtr=fopen(buffer,"w");
  fprintf(FilePtr,"# 2-theta   d         h  k  l  Mult Lp          Scat. Factor    Intensity\n");
  for(j=0;j<NumberOfPeaks;j++)
  {
    if((fabs(PowderDiffractionPeaks[j].Intensity)>1e-10)&&(PowderDiffractionPeaks[j].two_theta*RAD2DEG<Diffraction.two_theta_max))
      fprintf(FilePtr,"%9.5f %9.5f [ %2d, %2d, %2d] %4d %9.5f %18.10f %12.6f\n",
         PowderDiffractionPeaks[j].two_theta*RAD2DEG,
         PowderDiffractionPeaks[j].d,
         PowderDiffractionPeaks[j].h,PowderDiffractionPeaks[j].k,PowderDiffractionPeaks[j].l,
         PowderDiffractionPeaks[j].Multiplicity,
         PowderDiffractionPeaks[j].Lp,
         PowderDiffractionPeaks[j].ScatteringFactor,
         PowderDiffractionPeaks[j].Intensity);
  }
  fclose(FilePtr);

  free(Diffraction.spectrum);
}

void PowderDiffraction(void)
{
  VECTOR a,b,c;
  VECTOR k1,k2,k3;
  REAL thetamax;

  PowderDiffractionPeaks=(POWDER_DIFF_PEAK*)calloc(MAX_POWDER_DIFFRACTION_PEAKS,sizeof(POWDER_DIFF_PEAK));

  for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
  {
    a.x=Box[CurrentSystem].ax; a.y=Box[CurrentSystem].ay; a.z=Box[CurrentSystem].az;
    b.x=Box[CurrentSystem].bx; b.y=Box[CurrentSystem].by; b.z=Box[CurrentSystem].bz;
    c.x=Box[CurrentSystem].cx; c.y=Box[CurrentSystem].cy; c.z=Box[CurrentSystem].cz;

    k1=CrossProduct(b,c);
    k1.x/=Volume[CurrentSystem]; k1.y/=Volume[CurrentSystem]; k1.z/=Volume[CurrentSystem];

    k2=CrossProduct(c,a);
    k2.x/=Volume[CurrentSystem]; k2.y/=Volume[CurrentSystem]; k2.z/=Volume[CurrentSystem];

    k3=CrossProduct(a,b);
    k3.x/=Volume[CurrentSystem]; k3.y/=Volume[CurrentSystem]; k3.z/=Volume[CurrentSystem];

    // determine range of the Miller indices (h,k,l)
    thetamax=0.5*Diffraction.two_theta_max*DEG2RAD;
    Diffraction.hmax=(int)(2.0*sin(thetamax)/(Diffraction.lambda*sqrt(SQR(k1.x)+SQR(k1.y)+SQR(k1.z))));
    Diffraction.kmax=(int)(2.0*sin(thetamax)/(Diffraction.lambda*sqrt(SQR(k2.x)+SQR(k2.y)+SQR(k2.z))));
    Diffraction.lmax=(int)(2.0*sin(thetamax)/(Diffraction.lambda*sqrt(SQR(k3.x)+SQR(k3.y)+SQR(k3.z))));

    FindUniqueMillerPlanes();
    HeapSortMillerPlanes();
    CalculatePowderDiffractionPattern();
  }
}

// NEW1
int ProjectConstraintsFromHessianMatrixMassWeighted(int np,int nb,REAL *Gradient,REAL_MATRIX HessianMatrix,REAL *Weights)
{
  int i,j,k,m,n,f1,l,ia;
  int A,B,C,D;
  int Type,NumberOfConstraints;
  INT_VECTOR3 index_i,index_j,index_k,index_l;
  int NumberOfBonds,NumberOfBends,NumberOfTorsions,NumberOfImproperTorsions,NumberOfInversionBends;
  REAL dot_produkt,norm,length;
  REAL_MATRIX PrimativeUnitvectors,PrimativeUnitvectorsTranspose,ProjectionOperator,ProjectionOperator2;
  VECTOR posA,posB,posC,posD,com;
  REAL *ModifiedGradient;
  VECTOR ha,hb,hc,hd;

  n=np+nb;
  CurrentSystem=0;

  NumberOfConstraints=0;
  if(RemoveTranslationFromHessian) NumberOfConstraints+=Dimension;
  if(RemoveRotationFromHessian) NumberOfConstraints+=Dimension;

  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    Type=Adsorbates[CurrentSystem][m].Type;
    NumberOfBonds=Components[Type].NumberOfBonds;
    for(i=0;i<NumberOfBonds;i++)
      if(Components[Type].BondType[i]==FIXED_BOND) NumberOfConstraints++;
    NumberOfBends=Components[Type].NumberOfBends;
    for(i=0;i<NumberOfBends;i++)
      if(Components[Type].BendType[i]==FIXED_BEND) NumberOfConstraints++;
    NumberOfTorsions=Components[Type].NumberOfTorsions;
    for(i=0;i<NumberOfTorsions;i++)
      if(Components[Type].TorsionType[i]==FIXED_DIHEDRAL) NumberOfConstraints++;
    NumberOfImproperTorsions=Components[Type].NumberOfImproperTorsions;
    for(i=0;i<NumberOfImproperTorsions;i++)
      if(Components[Type].ImproperTorsionType[i]==FIXED_IMPROPER_DIHEDRAL) NumberOfConstraints++;
    NumberOfInversionBends=Components[Type].NumberOfInversionBends;
    for(i=0;i<NumberOfInversionBends;i++)
      if(Components[Type].InversionBendType[i]==FIXED_INVERSION_BEND) NumberOfConstraints++;
  }
  NumberOfConstraints+=NumberOfDistanceConstraints[CurrentSystem];
  NumberOfConstraints+=NumberOfAngleConstraints[CurrentSystem];
  NumberOfConstraints+=NumberOfDihedralConstraints[CurrentSystem];
  NumberOfConstraints+=NumberOfImproperDihedralConstraints[CurrentSystem];
  NumberOfConstraints+=NumberOfInversionBendConstraints[CurrentSystem];

  // allocate the temporary memory
  PrimativeUnitvectors=CreateRealMatrix(NumberOfConstraints,n);
  PrimativeUnitvectorsTranspose=CreateRealMatrix(n,NumberOfConstraints);
  ProjectionOperator=CreateRealMatrix(n,n);
  ProjectionOperator2=CreateRealMatrix(n,n);
  ModifiedGradient=(REAL*)calloc(n,sizeof(REAL));

  // initialize number of constraints
  NumberOfConstraints=0;

  // get the unweighted center of the system
  com=GetCenterOfMassCurrentSystem();

  // remove translation
  if(RemoveTranslationFromHessian)
  {
    for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
      {
        index_i=Framework[CurrentSystem].Atoms[f1][i].HessianIndex;

        switch(Dimension)
        {
          case 2:
            if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.x]=1.0/Weights[index_i.x];
            if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.y]=0.0;

            if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.x]=0.0;
            if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.y]=1.0/Weights[index_i.y];
            break;
          case 3:
            if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.x]=1.0/Weights[index_i.x];
            if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.y]=0.0;
            if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.z]=0.0;

            if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.x]=0.0;
            if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.y]=1.0/Weights[index_i.y];
            if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.z]=0.0;

            if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints+2][index_i.x]=0.0;
            if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints+2][index_i.y]=0.0;
            if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints+2][index_i.z]=1.0/Weights[index_i.z];
            break;
        }

      }
    }

    for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
    {
      Type=Adsorbates[CurrentSystem][m].Type;
      for(l=0;l<Components[Type].NumberOfGroups;l++)
      {
        if(Components[Type].Groups[l].Rigid) // rigid unit
        {
          index_i=Adsorbates[CurrentSystem][m].Groups[l].HessianIndex;

          switch(Dimension)
          {
            case 2:
              if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.x]=1.0/Weights[index_i.x];
              if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.y]=0.0;

              if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.x]=0.0;
              if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.y]=1.0/Weights[index_i.y];
              break;
            case 3:
              if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.x]=1.0/Weights[index_i.x];
              if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.y]=0.0;
              if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.z]=0.0;

              if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.x]=0.0;
              if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.y]=1.0/Weights[index_i.y];
              if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.z]=0.0;

              if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints+2][index_i.x]=0.0;
              if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints+2][index_i.y]=0.0;
              if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints+2][index_i.z]=1.0/Weights[index_i.z];
              break;
          }
        }
        else
        {
          for(ia=0;ia<Components[Type].Groups[l].NumberOfGroupAtoms;ia++)
          {
            A=Components[Type].Groups[l].Atoms[ia];
            index_i=Adsorbates[CurrentSystem][m].Atoms[A].HessianIndex;

            switch(Dimension)
            {
              case 2:
                if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.x]=1.0/Weights[index_i.x];
                if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.y]=0.0;

                if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.x]=0.0;
                if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.y]=1.0/Weights[index_i.y];
                break;
              case 3:
                if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.x]=1.0/Weights[index_i.x];
                if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.y]=0.0;
                if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.z]=0.0;

                if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.x]=0.0;
                if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.y]=1.0/Weights[index_i.y];
                if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.z]=0.0;

                if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints+2][index_i.x]=0.0;
                if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints+2][index_i.y]=0.0;
                if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints+2][index_i.z]=1.0/Weights[index_i.z];
                break;
            }
          }
        }
      }
    }

    for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
    {
      Type=Cations[CurrentSystem][m].Type;
      for(l=0;l<Components[Type].NumberOfGroups;l++)
      {
        if(Components[Type].Groups[l].Rigid) // rigid unit
        {
          index_i=Cations[CurrentSystem][m].Groups[l].HessianIndex;

          switch(Dimension)
          {
            case 2:
              if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.x]=1.0/Weights[index_i.x];
              if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.y]=0.0;

              if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.x]=0.0;
              if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.y]=1.0/Weights[index_i.y];
              break;
            case 3:
              if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.x]=1.0/Weights[index_i.x];
              if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.y]=0.0;
              if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.z]=0.0;

              if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.x]=0.0;
              if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.y]=1.0/Weights[index_i.y];
              if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.z]=0.0;

              if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints+2][index_i.x]=0.0;
              if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints+2][index_i.y]=0.0;
              if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints+2][index_i.z]=1.0/Weights[index_i.z];
              break;
          }
        }
        else
        {
          for(ia=0;ia<Components[Type].Groups[l].NumberOfGroupAtoms;ia++)
          {
            A=Components[Type].Groups[l].Atoms[ia];
            index_i=Cations[CurrentSystem][m].Atoms[A].HessianIndex;

            switch(Dimension)
            {
              case 2:
                if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.x]=1.0/Weights[index_i.x];
                if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.y]=0.0;

                if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.x]=0.0;
                if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.y]=1.0/Weights[index_i.y];
                break;
              case 3:
                if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.x]=1.0/Weights[index_i.x];
                if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.y]=0.0;
                if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.z]=0.0;

                if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.x]=0.0;
                if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.y]=1.0/Weights[index_i.y];
                if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.z]=0.0;

                if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints+2][index_i.x]=0.0;
                if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints+2][index_i.y]=0.0;
                if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints+2][index_i.z]=1.0/Weights[index_i.z];
                break;
            }
          }
        }
      }
    }

    NumberOfConstraints+=3;
  }


  // remove rotation
  if(RemoveRotationFromHessian)
  {
    for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
      {
        index_i=Framework[CurrentSystem].Atoms[f1][i].HessianIndex;
        posA=Framework[CurrentSystem].Atoms[f1][i].Position;

        switch(Dimension)
        {
          case 2:
            if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.x]=0.0;
            if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.y]=-(posA.z-com.z)/Weights[index_i.y];

            if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.x]=(posA.z-com.z)/Weights[index_i.x];
            if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.y]=0.0;
            break;
          case 3:
            if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.x]=0.0;
            if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.y]=-(posA.z-com.z)/Weights[index_i.y];
            if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.z]=(posA.y-com.y)/Weights[index_i.z];

            if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.x]=(posA.z-com.z)/Weights[index_i.x];
            if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.y]=0.0;
            if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.z]=-(posA.x-com.x)/Weights[index_i.z];

            if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints+2][index_i.x]=-(posA.y-com.y)/Weights[index_i.x];
            if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints+2][index_i.y]=(posA.x-com.x)/Weights[index_i.y];
            if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints+2][index_i.z]=0.0;
            break;
        }
      }
    }

    for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
    {
      Type=Adsorbates[CurrentSystem][m].Type;
      for(l=0;l<Components[Type].NumberOfGroups;l++)
      {
        if(Components[Type].Groups[l].Rigid) // rigid unit
        {
          index_i=Adsorbates[CurrentSystem][m].Groups[l].HessianIndex;
          posA=Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassPosition;

          switch(Dimension)
          {
            case 2:
              if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.x]=0.0;
              if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.y]=-(posA.z-com.z)/Weights[index_i.y];

              if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.x]=(posA.z-com.z)/Weights[index_i.x];
              if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.y]=0.0;
              break;
            case 3:
              if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.x]=0.0;
              if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.y]=-(posA.z-com.z)/Weights[index_i.y];
              if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.z]=(posA.y-com.y)/Weights[index_i.z];

              if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.x]=(posA.z-com.z)/Weights[index_i.x];
              if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.y]=0.0;
              if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.z]=-(posA.x-com.x)/Weights[index_i.z];

              if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints+2][index_i.x]=-(posA.y-com.y)/Weights[index_i.x];
              if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints+2][index_i.y]=(posA.x-com.x)/Weights[index_i.y];
              if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints+2][index_i.z]=0.0;
              break;
          }
        }
        else
        {
          for(ia=0;ia<Components[Type].Groups[l].NumberOfGroupAtoms;ia++)
          {
            A=Components[Type].Groups[l].Atoms[ia];
            index_i=Adsorbates[CurrentSystem][m].Atoms[A].HessianIndex;
            posA=Adsorbates[CurrentSystem][m].Atoms[A].Position;

            switch(Dimension)
            {
              case 2:
                if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.x]=0.0;
                if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.y]=-(posA.z-com.z)/Weights[index_i.y];

                if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.x]=(posA.z-com.z)/Weights[index_i.x];
                if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.y]=0.0;
                break;
              case 3:
                if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.x]=0.0;
                if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.y]=-(posA.z-com.z)/Weights[index_i.y];
                if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.z]=(posA.y-com.y)/Weights[index_i.z];

                if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.x]=(posA.z-com.z)/Weights[index_i.x];
                if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.y]=0.0;
                if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.z]=-(posA.x-com.x)/Weights[index_i.z];

                if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints+2][index_i.x]=-(posA.y-com.y)/Weights[index_i.x];
                if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints+2][index_i.y]=(posA.x-com.x)/Weights[index_i.y];
                if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints+2][index_i.z]=0.0;
                break;
            }
          }
        }
      }
    }

    for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
    {
      Type=Cations[CurrentSystem][m].Type;
      for(l=0;l<Components[Type].NumberOfGroups;l++)
      {
        if(Components[Type].Groups[l].Rigid) // rigid unit
        {
          index_i=Cations[CurrentSystem][m].Groups[l].HessianIndex;
          posA=Cations[CurrentSystem][m].Groups[l].CenterOfMassPosition;

          switch(Dimension)
          {
            case 2:
              if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.x]=0.0;
              if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.y]=-(posA.z-com.z)/Weights[index_i.y];

              if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.x]=(posA.z-com.z)/Weights[index_i.x];
              if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.y]=0.0;
              break;
            case 3:
              if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.x]=0.0;
              if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.y]=-(posA.z-com.z)/Weights[index_i.y];
              if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.z]=(posA.y-com.y)/Weights[index_i.z];

              if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.x]=(posA.z-com.z)/Weights[index_i.x];
              if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.y]=0.0;
              if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.z]=-(posA.x-com.x)/Weights[index_i.z];

              if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints+2][index_i.x]=-(posA.y-com.y)/Weights[index_i.x];
              if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints+2][index_i.y]=(posA.x-com.x)/Weights[index_i.y];
              if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints+2][index_i.z]=0.0;
              break;
          }
        }
        else
        {
          for(ia=0;ia<Components[Type].Groups[l].NumberOfGroupAtoms;ia++)
          {
            A=Components[Type].Groups[l].Atoms[ia];
            index_i=Cations[CurrentSystem][m].Atoms[A].HessianIndex;
            posA=Cations[CurrentSystem][m].Atoms[A].Position;

            switch(Dimension)
            {
              case 2:
                if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.x]=0.0;
                if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.y]=-(posA.z-com.z)/Weights[index_i.y];

                if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.x]=(posA.z-com.z)/Weights[index_i.x];
                if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.y]=0.0;
                break;
              case 3:
                if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.x]=0.0;
                if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.y]=-(posA.z-com.z)/Weights[index_i.y];
                if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.z]=(posA.y-com.y)/Weights[index_i.z];

                if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.x]=(posA.z-com.z)/Weights[index_i.x];
                if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.y]=0.0;
                if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.z]=-(posA.x-com.x)/Weights[index_i.z];

                if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints+2][index_i.x]=-(posA.y-com.y)/Weights[index_i.x];
                if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints+2][index_i.y]=(posA.x-com.x)/Weights[index_i.y];
                if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints+2][index_i.z]=0.0;
                break;
            }
          }
        }
      }
    }


    NumberOfConstraints+=3;
  }

  // add bond contraints
  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    Type=Adsorbates[CurrentSystem][m].Type;
    NumberOfBonds=Components[Type].NumberOfBonds;
    for(i=0;i<NumberOfBonds;i++)
    {
      if(Components[Type].BondType[i]==FIXED_BOND)
      {
        A=Components[Type].Bonds[i].A;
        B=Components[Type].Bonds[i].B;

        index_i=Adsorbates[CurrentSystem][m].Atoms[A].HessianIndex;
        index_j=Adsorbates[CurrentSystem][m].Atoms[B].HessianIndex;

        //if((index_i<0)&&(index_j<0)) continue;

        posA=Adsorbates[CurrentSystem][m].Atoms[A].Position;
        posB=Adsorbates[CurrentSystem][m].Atoms[B].Position;

        ReturnWilsonVectorsBond(posA,posB,&ha,&hb);

        switch(Dimension)
        {
          case 3:
            if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.z]=Weights[index_i.z]*ha.z;
          case 2:
            if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.y]=Weights[index_i.y]*ha.y;
          case 1:
            if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.x]=Weights[index_i.x]*ha.x;
            break;
        }

        switch(Dimension)
        {
          case 3:
            if(index_j.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_j.z]=Weights[index_j.z]*hb.z;
          case 2:
            if(index_j.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_j.y]=Weights[index_j.y]*hb.y;
          case 1:
            if(index_j.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_j.x]=Weights[index_j.x]*hb.x;
            break;
        }

        length=0.0;
        for(j=0;j<n;j++)
          length+=SQR(PrimativeUnitvectors.element[NumberOfConstraints][j]);
        length=1.0/sqrt(length);
        for(j=0;j<n;j++)
          PrimativeUnitvectors.element[NumberOfConstraints][j]*=length;

        NumberOfConstraints++;
      }
    }
  }

  // add distance contraints
  for(m=0;m<NumberOfDistanceConstraints[CurrentSystem];m++)
  {
    index_i=DistanceConstraints[CurrentSystem][m][0]->HessianIndex;
    index_j=DistanceConstraints[CurrentSystem][m][1]->HessianIndex;

    posA=DistanceConstraints[CurrentSystem][m][0]->Position;
    posB=DistanceConstraints[CurrentSystem][m][1]->Position;

    ReturnWilsonVectorsBond(posA,posB,&ha,&hb);

    switch(Dimension)
    {
      case 3:
        if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.z]=Weights[index_i.z]*ha.z;
      case 2:
        if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.y]=Weights[index_i.y]*ha.y;
      case 1:
        if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.x]=Weights[index_i.x]*ha.x;
        break;
    }

    switch(Dimension)
    {
      case 3:
        if(index_j.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_j.z]=Weights[index_j.z]*hb.z;
      case 2:
        if(index_j.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_j.y]=Weights[index_j.y]*hb.y;
      case 1:
        if(index_j.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_j.x]=Weights[index_j.x]*hb.x;
        break;
    }

    length=0.0;
    for(j=0;j<n;j++)
      length+=SQR(PrimativeUnitvectors.element[NumberOfConstraints][j]);
    length=1.0/sqrt(length);
    for(j=0;j<n;j++)
      PrimativeUnitvectors.element[NumberOfConstraints][j]*=length;

    NumberOfConstraints++;
  }


  // add bond-angle contraints
  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    Type=Adsorbates[CurrentSystem][m].Type;
    NumberOfBends=Components[Type].NumberOfBends;
    for(i=0;i<NumberOfBends;i++)
    {
      if(Components[Type].BendType[i]==FIXED_BEND)
      {
        A=Components[Type].Bends[i].A;
        B=Components[Type].Bends[i].B;
        C=Components[Type].Bends[i].C;

        index_i=Adsorbates[CurrentSystem][m].Atoms[A].HessianIndex;
        index_j=Adsorbates[CurrentSystem][m].Atoms[B].HessianIndex;
        index_k=Adsorbates[CurrentSystem][m].Atoms[C].HessianIndex;

        posA=Adsorbates[CurrentSystem][m].Atoms[A].Position;
        posB=Adsorbates[CurrentSystem][m].Atoms[B].Position;
        posC=Adsorbates[CurrentSystem][m].Atoms[C].Position;

        ReturnWilsonVectorsBend(posA,posB,posC,&ha,&hb,&hc);

        switch(Dimension)
        {
          case 3:
            if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.z]=Weights[index_i.z]*ha.z;
          case 2:
            if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.y]=Weights[index_i.y]*ha.y;
          case 1:
            if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.x]=Weights[index_i.x]*ha.x;
            break;
        }

        switch(Dimension)
        {
          case 3:
            if(index_j.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_j.z]=Weights[index_j.z]*hb.z;
          case 2:
            if(index_j.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_j.y]=Weights[index_j.y]*hb.y;
          case 1:
            if(index_j.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_j.x]=Weights[index_j.x]*hb.x;
            break;
        }

        switch(Dimension)
        {
          case 3:
            if(index_k.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_k.z]=Weights[index_k.z]*hc.z;
          case 2:
            if(index_k.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_k.y]=Weights[index_k.y]*hc.y;
          case 1:
            if(index_k.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_k.x]=Weights[index_k.x]*hc.x;
            break;
        }

        NumberOfConstraints++;
      }
    }
  }


  // add general angle contraints
  for(m=0;m<NumberOfAngleConstraints[CurrentSystem];m++)
  {
    index_i=AngleConstraints[CurrentSystem][m][0]->HessianIndex;
    index_j=AngleConstraints[CurrentSystem][m][1]->HessianIndex;
    index_k=AngleConstraints[CurrentSystem][m][2]->HessianIndex;

    posA=AngleConstraints[CurrentSystem][m][0]->Position;
    posB=AngleConstraints[CurrentSystem][m][1]->Position;
    posC=AngleConstraints[CurrentSystem][m][2]->Position;

    ReturnWilsonVectorsBend(posA,posB,posC,&ha,&hb,&hc);

    switch(Dimension)
    {
      case 3:
        if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.z]=Weights[index_i.z]*ha.z;
      case 2:
        if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.y]=Weights[index_i.y]*ha.y;
      case 1:
        if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.x]=Weights[index_i.x]*ha.x;
        break;
    }

    switch(Dimension)
    {
      case 3:
        if(index_j.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_j.z]=Weights[index_j.z]*hb.z;
      case 2:
        if(index_j.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_j.y]=Weights[index_j.y]*hb.y;
      case 1:
        if(index_j.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_j.x]=Weights[index_j.x]*hb.x;
        break;
    }

    switch(Dimension)
    {
      case 3:
        if(index_k.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_k.z]=Weights[index_k.z]*hc.z;
      case 2:
        if(index_k.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_k.y]=Weights[index_k.y]*hc.y;
      case 1:
        if(index_k.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_k.x]=Weights[index_k.x]*hc.x;
        break;
    }

    NumberOfConstraints++;
  }




  // add dihedral contraints
  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    Type=Adsorbates[CurrentSystem][m].Type;
    NumberOfTorsions=Components[Type].NumberOfTorsions;
    for(i=0;i<NumberOfTorsions;i++)
    {
      if(Components[Type].TorsionType[i]==FIXED_DIHEDRAL)
      {
        A=Components[Type].Torsions[i].A;
        B=Components[Type].Torsions[i].B;
        C=Components[Type].Torsions[i].C;
        D=Components[Type].Torsions[i].D;

        index_i=Adsorbates[CurrentSystem][m].Atoms[A].HessianIndex;
        index_j=Adsorbates[CurrentSystem][m].Atoms[B].HessianIndex;
        index_k=Adsorbates[CurrentSystem][m].Atoms[C].HessianIndex;
        index_l=Adsorbates[CurrentSystem][m].Atoms[D].HessianIndex;

        posA=Adsorbates[CurrentSystem][m].Atoms[A].Position;
        posB=Adsorbates[CurrentSystem][m].Atoms[B].Position;
        posC=Adsorbates[CurrentSystem][m].Atoms[C].Position;
        posD=Adsorbates[CurrentSystem][m].Atoms[D].Position;

        ReturnWilsonVectorsTorsion(posA,posB,posC,posD,&ha,&hb,&hc,&hd);

        if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.x]=Weights[index_i.x]*ha.x;
        if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.y]=Weights[index_i.y]*ha.y;
        if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.z]=Weights[index_i.z]*ha.z;

        if(index_j.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_j.x]=Weights[index_j.x]*hb.x;
        if(index_j.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_j.y]=Weights[index_j.y]*hb.y;
        if(index_j.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_j.z]=Weights[index_j.z]*hb.z;

        if(index_k.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_k.x]=Weights[index_k.x]*hc.x;
        if(index_k.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_k.y]=Weights[index_k.y]*hc.y;
        if(index_k.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_k.z]=Weights[index_k.z]*hc.z;

        if(index_l.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_l.x]=Weights[index_l.x]*hd.x;
        if(index_l.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_l.y]=Weights[index_l.y]*hd.y;
        if(index_l.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_l.z]=Weights[index_l.z]*hd.z;

        NumberOfConstraints++;
      }
    }
  }

  // add dihedral contraints
  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    Type=Adsorbates[CurrentSystem][m].Type;
    NumberOfImproperTorsions=Components[Type].NumberOfImproperTorsions;
    for(i=0;i<NumberOfImproperTorsions;i++)
    {
      if(Components[Type].ImproperTorsionType[i]==FIXED_IMPROPER_DIHEDRAL)
      {
        A=Components[Type].ImproperTorsions[i].A;
        B=Components[Type].ImproperTorsions[i].B;
        C=Components[Type].ImproperTorsions[i].C;
        D=Components[Type].ImproperTorsions[i].D;

        index_i=Adsorbates[CurrentSystem][m].Atoms[A].HessianIndex;
        index_j=Adsorbates[CurrentSystem][m].Atoms[B].HessianIndex;
        index_k=Adsorbates[CurrentSystem][m].Atoms[C].HessianIndex;
        index_l=Adsorbates[CurrentSystem][m].Atoms[D].HessianIndex;

        posA=Adsorbates[CurrentSystem][m].Atoms[A].Position;
        posB=Adsorbates[CurrentSystem][m].Atoms[B].Position;
        posC=Adsorbates[CurrentSystem][m].Atoms[C].Position;
        posD=Adsorbates[CurrentSystem][m].Atoms[D].Position;

        ReturnWilsonVectorsTorsion(posA,posB,posC,posD,&ha,&hb,&hc,&hd);

        if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.x]=Weights[index_i.x]*ha.x;
        if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.y]=Weights[index_i.y]*ha.y;
        if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.z]=Weights[index_i.z]*ha.z;

        if(index_j.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_j.x]=Weights[index_j.x]*hb.x;
        if(index_j.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_j.y]=Weights[index_j.y]*hb.y;
        if(index_j.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_j.z]=Weights[index_j.z]*hb.z;

        if(index_k.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_k.x]=Weights[index_k.x]*hc.x;
        if(index_k.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_k.y]=Weights[index_k.y]*hc.y;
        if(index_k.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_k.z]=Weights[index_k.z]*hc.z;

        if(index_l.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_l.x]=Weights[index_l.x]*hd.x;
        if(index_l.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_l.y]=Weights[index_l.y]*hd.y;
        if(index_l.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_l.z]=Weights[index_l.z]*hd.z;

        NumberOfConstraints++;
      }
    }
  }


  // add general dihedral angle contraints
  for(m=0;m<NumberOfDihedralConstraints[CurrentSystem];m++)
  {
    index_i=DihedralConstraints[CurrentSystem][m][0]->HessianIndex;
    index_j=DihedralConstraints[CurrentSystem][m][1]->HessianIndex;
    index_k=DihedralConstraints[CurrentSystem][m][2]->HessianIndex;
    index_l=DihedralConstraints[CurrentSystem][m][3]->HessianIndex;

    posA=DihedralConstraints[CurrentSystem][m][0]->Position;
    posB=DihedralConstraints[CurrentSystem][m][1]->Position;
    posC=DihedralConstraints[CurrentSystem][m][2]->Position;
    posD=DihedralConstraints[CurrentSystem][m][3]->Position;

    ReturnWilsonVectorsTorsion(posA,posB,posC,posD,&ha,&hb,&hc,&hd);

    if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.x]=Weights[index_i.x]*ha.x;
    if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.y]=Weights[index_i.y]*ha.y;
    if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.z]=Weights[index_i.z]*ha.z;

    if(index_j.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_j.x]=Weights[index_j.x]*hb.x;
    if(index_j.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_j.y]=Weights[index_j.y]*hb.y;
    if(index_j.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_j.z]=Weights[index_j.z]*hb.z;

    if(index_k.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_k.x]=Weights[index_k.x]*hc.x;
    if(index_k.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_k.y]=Weights[index_k.y]*hc.y;
    if(index_k.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_k.z]=Weights[index_k.z]*hc.z;

    if(index_l.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_l.x]=Weights[index_l.x]*hd.x;
    if(index_l.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_l.y]=Weights[index_l.y]*hd.y;
    if(index_l.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_l.z]=Weights[index_l.z]*hd.z;

    NumberOfConstraints++;
  }

  // add general improper dihedral angle contraints
  for(m=0;m<NumberOfImproperDihedralConstraints[CurrentSystem];m++)
  {
    index_i=ImproperDihedralConstraints[CurrentSystem][m][0]->HessianIndex;
    index_j=ImproperDihedralConstraints[CurrentSystem][m][1]->HessianIndex;
    index_k=ImproperDihedralConstraints[CurrentSystem][m][2]->HessianIndex;
    index_l=ImproperDihedralConstraints[CurrentSystem][m][3]->HessianIndex;

    posA=ImproperDihedralConstraints[CurrentSystem][m][0]->Position;
    posB=ImproperDihedralConstraints[CurrentSystem][m][1]->Position;
    posC=ImproperDihedralConstraints[CurrentSystem][m][2]->Position;
    posD=ImproperDihedralConstraints[CurrentSystem][m][3]->Position;

    ReturnWilsonVectorsTorsion(posA,posB,posC,posD,&ha,&hb,&hc,&hd);

    if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.x]=Weights[index_i.x]*ha.x;
    if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.y]=Weights[index_i.y]*ha.y;
    if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.z]=Weights[index_i.z]*ha.z;

    if(index_j.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_j.x]=Weights[index_j.x]*hb.x;
    if(index_j.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_j.y]=Weights[index_j.y]*hb.y;
    if(index_j.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_j.z]=Weights[index_j.z]*hb.z;

    if(index_k.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_k.x]=Weights[index_k.x]*hc.x;
    if(index_k.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_k.y]=Weights[index_k.y]*hc.y;
    if(index_k.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_k.z]=Weights[index_k.z]*hc.z;

    if(index_l.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_l.x]=Weights[index_l.x]*hd.x;
    if(index_l.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_l.y]=Weights[index_l.y]*hd.y;
    if(index_l.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_l.z]=Weights[index_l.z]*hd.z;

    NumberOfConstraints++;
  }



  // add inversion bend contraints
  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    Type=Adsorbates[CurrentSystem][m].Type;
    NumberOfInversionBends=Components[Type].NumberOfInversionBends;
    for(i=0;i<NumberOfInversionBends;i++)
    {
      if(Components[Type].InversionBendType[i]==FIXED_INVERSION_BEND)
      {
        A=Components[Type].InversionBends[i].A;
        B=Components[Type].InversionBends[i].B;
        C=Components[Type].InversionBends[i].C;
        D=Components[Type].InversionBends[i].D;

        index_i=Adsorbates[CurrentSystem][m].Atoms[A].HessianIndex;
        index_j=Adsorbates[CurrentSystem][m].Atoms[B].HessianIndex;
        index_k=Adsorbates[CurrentSystem][m].Atoms[C].HessianIndex;
        index_l=Adsorbates[CurrentSystem][m].Atoms[D].HessianIndex;

        posA=Adsorbates[CurrentSystem][m].Atoms[A].Position;
        posB=Adsorbates[CurrentSystem][m].Atoms[B].Position;
        posC=Adsorbates[CurrentSystem][m].Atoms[C].Position;
        posD=Adsorbates[CurrentSystem][m].Atoms[D].Position;

        ReturnWilsonVectorsInversionBend(posA,posB,posC,posD,&ha,&hb,&hc,&hd);

        if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.x]=Weights[index_i.x]*ha.x;
        if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.y]=Weights[index_i.y]*ha.y;
        if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.z]=Weights[index_i.z]*ha.z;

        if(index_j.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_j.x]=Weights[index_j.x]*hb.x;
        if(index_j.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_j.y]=Weights[index_j.y]*hb.y;
        if(index_j.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_j.z]=Weights[index_j.z]*hb.z;

        if(index_k.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_k.x]=Weights[index_k.x]*hc.x;
        if(index_k.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_k.y]=Weights[index_k.y]*hc.y;
        if(index_k.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_k.z]=Weights[index_k.z]*hc.z;

        if(index_l.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_l.x]=Weights[index_l.x]*hd.x;
        if(index_l.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_l.y]=Weights[index_l.y]*hd.y;
        if(index_l.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_l.z]=Weights[index_l.z]*hd.z;

        NumberOfConstraints++;
      }
    }
  }

  // add general inversion-bend angle contraints
  for(m=0;m<NumberOfInversionBendConstraints[CurrentSystem];m++)
  {
    index_i=InversionBendConstraints[CurrentSystem][m][0]->HessianIndex;
    index_j=InversionBendConstraints[CurrentSystem][m][1]->HessianIndex;
    index_k=InversionBendConstraints[CurrentSystem][m][2]->HessianIndex;
    index_l=InversionBendConstraints[CurrentSystem][m][3]->HessianIndex;

    //if((index_i<0)&&(index_j<0)&&(index_k<0)&&(index_l<0)) continue;

    posA=InversionBendConstraints[CurrentSystem][m][0]->Position;
    posB=InversionBendConstraints[CurrentSystem][m][1]->Position;
    posC=InversionBendConstraints[CurrentSystem][m][2]->Position;
    posD=InversionBendConstraints[CurrentSystem][m][3]->Position;

    ReturnWilsonVectorsInversionBend(posA,posB,posC,posD,&ha,&hb,&hc,&hd);

    if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.x]=Weights[index_i.x]*ha.x;
    if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.y]=Weights[index_i.y]*ha.y;
    if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.z]=Weights[index_i.z]*ha.z;

    if(index_j.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_j.x]=Weights[index_j.x]*hb.x;
    if(index_j.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_j.y]=Weights[index_j.y]*hb.y;
    if(index_j.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_j.z]=Weights[index_j.z]*hb.z;

    if(index_k.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_k.x]=Weights[index_k.x]*hc.x;
    if(index_k.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_k.y]=Weights[index_k.y]*hc.y;
    if(index_k.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_k.z]=Weights[index_k.z]*hc.z;

    if(index_l.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_l.x]=Weights[index_l.x]*hd.x;
    if(index_l.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_l.y]=Weights[index_l.y]*hd.y;
    if(index_l.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_l.z]=Weights[index_l.z]*hd.z;

    NumberOfConstraints++;
  }


  // create orthonormal basis using simple Gram-Schmidt
  NumberOfRemovedModes=0;
  for(i=0;i<NumberOfConstraints;i++)
  {
    for(j=0;j<NumberOfRemovedModes;j++)
    {
      dot_produkt=0.0;
      for(k=0;k<n;k++)
        dot_produkt+=PrimativeUnitvectors.element[i][k]*PrimativeUnitvectors.element[j][k];
      for(k=0;k<n;k++)
        PrimativeUnitvectors.element[i][k]-=dot_produkt*PrimativeUnitvectors.element[j][k];
    }
    dot_produkt=0.0;
    for(k=0;k<n;k++)
      dot_produkt+=SQR(PrimativeUnitvectors.element[i][k]);
    if(dot_produkt>1e-10)
    {
      norm=sqrt(dot_produkt);
      for(k=0;k<n;k++)
        PrimativeUnitvectors.element[NumberOfRemovedModes][k]=PrimativeUnitvectors.element[i][k]/norm;
      NumberOfRemovedModes++;
    }
  }

  TransposeRealMatrix(PrimativeUnitvectorsTranspose,PrimativeUnitvectors);
  MultiplyRealMatrix(ProjectionOperator,PrimativeUnitvectorsTranspose,PrimativeUnitvectors);

  for(i=0;i<n;i++)
  {
    for(j=0;j<n;j++)
      ProjectionOperator.element[i][j]=-ProjectionOperator.element[i][j];
    ProjectionOperator.element[i][i]+=1.0;
  }

  MultiplyRealMatrixVector(ModifiedGradient,ProjectionOperator,Gradient);
  for(i=0;i<n;i++)
    Gradient[i]=ModifiedGradient[i];

  MultiplyRealMatrix(ProjectionOperator2,HessianMatrix,ProjectionOperator);
  MultiplyRealMatrix(HessianMatrix,ProjectionOperator,ProjectionOperator2);

  DeleteRealMatrix(PrimativeUnitvectors);
  DeleteRealMatrix(PrimativeUnitvectorsTranspose);
  DeleteRealMatrix(ProjectionOperator);
  DeleteRealMatrix(ProjectionOperator2);
  free(ModifiedGradient);

  return NumberOfRemovedModes;
}

// NEW2
int ProjectConstraintsFromHessianMatrix(int np,int nb,REAL *Gradient,REAL_MATRIX HessianMatrix,int ComputeGradient,int ComputeHessian)
{
  int i,j,k,m,n,f1,l,ia;
  int A,B,C,D;
  int Type,NumberOfConstraints;
  INT_VECTOR3 index_i,index_j,index_k,index_l;
  int NumberOfBonds,NumberOfBends,NumberOfTorsions,NumberOfImproperTorsions,NumberOfInversionBends;
  REAL dot_produkt,norm,length;
  REAL_MATRIX PrimativeUnitvectors,PrimativeUnitvectorsTranspose,ProjectionOperator,ProjectionOperator2;
  VECTOR posA,posB,posC,posD,com;
  REAL *ModifiedGradient;
  VECTOR ha,hb,hc,hd;

  n=np+nb;
  CurrentSystem=0;

  NumberOfConstraints=0;
  if(RemoveTranslationFromHessian) NumberOfConstraints+=3;
  if(RemoveRotationFromHessian) NumberOfConstraints+=3;
  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    Type=Adsorbates[CurrentSystem][m].Type;
    NumberOfBonds=Components[Type].NumberOfBonds;
    for(i=0;i<NumberOfBonds;i++)
      if(Components[Type].BondType[i]==FIXED_BOND) NumberOfConstraints++;
    NumberOfBends=Components[Type].NumberOfBends;
    for(i=0;i<NumberOfBends;i++)
      if(Components[Type].BendType[i]==FIXED_BEND) NumberOfConstraints++;
    NumberOfTorsions=Components[Type].NumberOfTorsions;
    for(i=0;i<NumberOfTorsions;i++)
      if(Components[Type].TorsionType[i]==FIXED_DIHEDRAL) NumberOfConstraints++;
    NumberOfImproperTorsions=Components[Type].NumberOfImproperTorsions;
    for(i=0;i<NumberOfImproperTorsions;i++)
      if(Components[Type].ImproperTorsionType[i]==FIXED_IMPROPER_DIHEDRAL) NumberOfConstraints++;
    NumberOfInversionBends=Components[Type].NumberOfInversionBends;
    for(i=0;i<NumberOfInversionBends;i++)
      if(Components[Type].InversionBendType[i]==FIXED_INVERSION_BEND) NumberOfConstraints++;
  }
  NumberOfConstraints+=NumberOfDistanceConstraints[CurrentSystem];
  NumberOfConstraints+=NumberOfAngleConstraints[CurrentSystem];
  NumberOfConstraints+=NumberOfDihedralConstraints[CurrentSystem];
  NumberOfConstraints+=NumberOfImproperDihedralConstraints[CurrentSystem];
  NumberOfConstraints+=NumberOfInversionBendConstraints[CurrentSystem];

  // allocate the temporary memory
  PrimativeUnitvectors=CreateRealMatrix(NumberOfConstraints,n);
  PrimativeUnitvectorsTranspose=CreateRealMatrix(n,NumberOfConstraints);
  ProjectionOperator=CreateRealMatrix(n,n);
  ProjectionOperator2=CreateRealMatrix(n,n);
  ModifiedGradient=(REAL*)calloc(n,sizeof(REAL));

  // initialize number of constraints
  NumberOfConstraints=0;

  // get the unweighted center of the system
  com=GetCenterOfCurrentSystem();

  // remove translation
  if(RemoveTranslationFromHessian)
  {
    for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
      {
        index_i=Framework[CurrentSystem].Atoms[f1][i].HessianIndex;

        switch(Dimension)
        {
          case 2:
            if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.x]=1.0;
            if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.y]=0.0;

            if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.x]=0.0;
            if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.y]=1.0;
            break;
          case 3:
            if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.x]=1.0;
            if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.y]=0.0;
            if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.z]=0.0;

            if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.x]=0.0;
            if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.y]=1.0;
            if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.z]=0.0;

            if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints+2][index_i.x]=0.0;
            if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints+2][index_i.y]=0.0;
            if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints+2][index_i.z]=1.0;
            break;
        }

      }
    }

    for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
    {
      Type=Adsorbates[CurrentSystem][m].Type;
      for(l=0;l<Components[Type].NumberOfGroups;l++)
      {
        if(Components[Type].Groups[l].Rigid) // rigid unit
        {
          index_i=Adsorbates[CurrentSystem][m].Groups[l].HessianIndex;

          switch(Dimension)
          {
            case 2:
              if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.x]=1.0;
              if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.y]=0.0;

              if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.x]=0.0;
              if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.y]=1.0;
              break;
            case 3:
              if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.x]=1.0;
              if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.y]=0.0;
              if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.z]=0.0;

              if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.x]=0.0;
              if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.y]=1.0;
              if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.z]=0.0;

              if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints+2][index_i.x]=0.0;
              if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints+2][index_i.y]=0.0;
              if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints+2][index_i.z]=1.0;
              break;
          }

        }
        else
        {
          for(ia=0;ia<Components[Type].Groups[l].NumberOfGroupAtoms;ia++)
          {
            A=Components[Type].Groups[l].Atoms[ia];
            index_i=Adsorbates[CurrentSystem][m].Atoms[A].HessianIndex;

            switch(Dimension)
            {
              case 2:
                if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.x]=1.0;
                if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.y]=0.0;

                if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.x]=0.0;
                if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.y]=1.0;
                break;
              case 3:
                if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.x]=1.0;
                if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.y]=0.0;
                if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.z]=0.0;

                if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.x]=0.0;
                if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.y]=1.0;
                if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.z]=0.0;

                if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints+2][index_i.x]=0.0;
                if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints+2][index_i.y]=0.0;
                if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints+2][index_i.z]=1.0;
                break;
            }
          }
        }
      }
    }

    for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
    {
      Type=Cations[CurrentSystem][m].Type;
      for(l=0;l<Components[Type].NumberOfGroups;l++)
      {
        if(Components[Type].Groups[l].Rigid) // rigid unit
        {
          index_i=Cations[CurrentSystem][m].Groups[l].HessianIndex;

          switch(Dimension)
          {
            case 2:
              if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.x]=1.0;
              if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.y]=0.0;

              if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.x]=0.0;
              if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.y]=1.0;
              break;
            case 3:
              if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.x]=1.0;
              if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.y]=0.0;
              if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.z]=0.0;

              if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.x]=0.0;
              if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.y]=1.0;
              if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.z]=0.0;

              if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints+2][index_i.x]=0.0;
              if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints+2][index_i.y]=0.0;
              if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints+2][index_i.z]=1.0;
              break;
          }

        }
        else
        {
          for(ia=0;ia<Components[Type].Groups[l].NumberOfGroupAtoms;ia++)
          {
            A=Components[Type].Groups[l].Atoms[ia];
            index_i=Cations[CurrentSystem][m].Atoms[A].HessianIndex;

            switch(Dimension)
            {
              case 2:
                if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.x]=1.0;
                if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.y]=0.0;

                if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.x]=0.0;
                if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.y]=1.0;
                break;
              case 3:
                if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.x]=1.0;
                if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.y]=0.0;
                if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.z]=0.0;

                if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.x]=0.0;
                if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.y]=1.0;
                if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.z]=0.0;

                if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints+2][index_i.x]=0.0;
                if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints+2][index_i.y]=0.0;
                if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints+2][index_i.z]=1.0;
                break;
            }
          }
        }
      }
    }

    NumberOfConstraints+=3;
  }

  // remove rotation
  if(RemoveRotationFromHessian)
  {
    for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
      {
        index_i=Framework[CurrentSystem].Atoms[f1][i].HessianIndex;
        posA=Framework[CurrentSystem].Atoms[f1][i].Position;

        switch(Dimension)
        {
          case 2:
            if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.x]=0.0;
            if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.y]=-(posA.z-com.z);

            if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.x]=posA.z-com.z;
            if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.y]=0.0;
            break;
          case 3:
            if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.x]=0.0;
            if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.y]=-(posA.z-com.z);
            if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.z]=posA.y-com.y;

            if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.x]=posA.z-com.z;
            if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.y]=0.0;
            if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.z]=-(posA.x-com.x);

            if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints+2][index_i.x]=-(posA.y-com.y);
            if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints+2][index_i.y]=posA.x-com.x;
            if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints+2][index_i.z]=0.0;
            break;
        }
      }
    }

    for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
    {
      Type=Adsorbates[CurrentSystem][m].Type;
      for(l=0;l<Components[Type].NumberOfGroups;l++)
      {
        if(Components[Type].Groups[l].Rigid) // rigid unit
        {
          index_i=Adsorbates[CurrentSystem][m].Groups[l].HessianIndex;
          posA=Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassPosition;

          switch(Dimension)
          {
            case 2:
              if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.x]=0.0;
              if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.y]=-(posA.z-com.z);

              if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.x]=posA.z-com.z;
              if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.y]=0.0;
              break;
            case 3:
              if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.x]=0.0;
              if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.y]=-(posA.z-com.z);
              if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.z]=posA.y-com.y;

              if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.x]=posA.z-com.z;
              if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.y]=0.0;
              if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.z]=-(posA.x-com.x);

              if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints+2][index_i.x]=-(posA.y-com.y);
              if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints+2][index_i.y]=posA.x-com.x;
              if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints+2][index_i.z]=0.0;
              break;
          }
        }
        else
        {
          for(ia=0;ia<Components[Type].Groups[l].NumberOfGroupAtoms;ia++)
          {
            A=Components[Type].Groups[l].Atoms[ia];
            index_i=Adsorbates[CurrentSystem][m].Atoms[A].HessianIndex;
            posA=Adsorbates[CurrentSystem][m].Atoms[A].Position;

            switch(Dimension)
            {
              case 2:
                if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.x]=0.0;
                if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.y]=-(posA.z-com.z);

                if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.x]=posA.z-com.z;
                if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.y]=0.0;
                break;
              case 3:
                if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.x]=0.0;
                if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.y]=-(posA.z-com.z);
                if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.z]=posA.y-com.y;

                if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.x]=posA.z-com.z;
                if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.y]=0.0;
                if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.z]=-(posA.x-com.x);

                if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints+2][index_i.x]=-(posA.y-com.y);
                if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints+2][index_i.y]=posA.x-com.x;
                if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints+2][index_i.z]=0.0;
                break;
            }
          }
        }
      }
    }

    for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
    {
      Type=Cations[CurrentSystem][m].Type;
      for(l=0;l<Components[Type].NumberOfGroups;l++)
      {
        if(Components[Type].Groups[l].Rigid) // rigid unit
        {
          index_i=Cations[CurrentSystem][m].Groups[l].HessianIndex;
          posA=Cations[CurrentSystem][m].Groups[l].CenterOfMassPosition;

          switch(Dimension)
          {
            case 2:
              if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.x]=0.0;
              if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.y]=-(posA.z-com.z);

              if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.x]=posA.z-com.z;
              if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.y]=0.0;
              break;
            case 3:
              if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.x]=0.0;
              if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.y]=-(posA.z-com.z);
              if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.z]=posA.y-com.y;

              if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.x]=posA.z-com.z;
              if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.y]=0.0;
              if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.z]=-(posA.x-com.x);

              if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints+2][index_i.x]=-(posA.y-com.y);
              if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints+2][index_i.y]=posA.x-com.x;
              if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints+2][index_i.z]=0.0;
              break;
          }

        }
        else
        {
          for(ia=0;ia<Components[Type].Groups[l].NumberOfGroupAtoms;ia++)
          {
            A=Components[Type].Groups[l].Atoms[ia];
            index_i=Cations[CurrentSystem][m].Atoms[A].HessianIndex;
            posA=Cations[CurrentSystem][m].Atoms[A].Position;

            switch(Dimension)
            {
              case 2:
                if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.x]=0.0;
                if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.y]=-(posA.z-com.z);

                if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.x]=posA.z-com.z;
                if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.y]=0.0;
                break;
              case 3:
                if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.x]=0.0;
                if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.y]=-(posA.z-com.z);
                if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.z]=posA.y-com.y;

                if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.x]=posA.z-com.z;
                if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.y]=0.0;
                if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints+1][index_i.z]=-(posA.x-com.x);

                if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints+2][index_i.x]=-(posA.y-com.y);
                if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints+2][index_i.y]=posA.x-com.x;
                if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints+2][index_i.z]=0.0;
                break;
            }
          }
        }
      }
    }


    NumberOfConstraints+=3;
  }

  // add bond contraints
  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    Type=Adsorbates[CurrentSystem][m].Type;
    NumberOfBonds=Components[Type].NumberOfBonds;
    for(i=0;i<NumberOfBonds;i++)
    {
      if(Components[Type].BondType[i]==FIXED_BOND)
      {
        A=Components[Type].Bonds[i].A;
        B=Components[Type].Bonds[i].B;

        index_i=Adsorbates[CurrentSystem][m].Atoms[A].HessianIndex;
        index_j=Adsorbates[CurrentSystem][m].Atoms[B].HessianIndex;

        posA=Adsorbates[CurrentSystem][m].Atoms[A].Position;
        posB=Adsorbates[CurrentSystem][m].Atoms[B].Position;

        ReturnWilsonVectorsBond(posA,posB,&ha,&hb);


        switch(Dimension)
        {
          case 3:
            if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.z]=ha.z;
          case 2:
            if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.y]=ha.y;
          case 1:
            if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.x]=ha.x;
            break;
        }

        switch(Dimension)
        {
          case 3:
            if(index_j.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_j.z]=hb.z;
          case 2:
            if(index_j.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_j.y]=hb.y;
          case 1:
            if(index_j.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_j.x]=hb.x;
            break;
        }

        length=0.0;
        for(j=0;j<n;j++)
          length+=SQR(PrimativeUnitvectors.element[NumberOfConstraints][j]);
        length=1.0/sqrt(length);
        for(j=0;j<n;j++)
          PrimativeUnitvectors.element[NumberOfConstraints][j]*=length;


        NumberOfConstraints++;
      }
    }
  }

  // add distance contraints
  for(m=0;m<NumberOfDistanceConstraints[CurrentSystem];m++)
  {
    index_i=DistanceConstraints[CurrentSystem][m][0]->HessianIndex;
    index_j=DistanceConstraints[CurrentSystem][m][1]->HessianIndex;

    posA=DistanceConstraints[CurrentSystem][m][0]->Position;
    posB=DistanceConstraints[CurrentSystem][m][1]->Position;

    ReturnWilsonVectorsBond(posA,posB,&ha,&hb);

    switch(Dimension)
    {
      case 3:
        if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.z]=ha.z;
      case 2:
        if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.y]=ha.y;
      case 1:
        if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.x]=ha.x;
        break;
    }

    switch(Dimension)
    {
      case 3:
        if(index_j.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_j.z]=hb.z;
      case 2:
        if(index_j.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_j.y]=hb.y;
      case 1:
        if(index_j.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_j.x]=hb.x;
        break;
    }

    length=0.0;
    for(j=0;j<n;j++)
      length+=SQR(PrimativeUnitvectors.element[NumberOfConstraints][j]);
    length=1.0/sqrt(length);
    for(j=0;j<n;j++)
      PrimativeUnitvectors.element[NumberOfConstraints][j]*=length;

    NumberOfConstraints++;
  }


  // add bond-angle contraints
  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    Type=Adsorbates[CurrentSystem][m].Type;
    NumberOfBends=Components[Type].NumberOfBends;
    for(i=0;i<NumberOfBends;i++)
    {
      if(Components[Type].BendType[i]==FIXED_BEND)
      {
        A=Components[Type].Bends[i].A;
        B=Components[Type].Bends[i].B;
        C=Components[Type].Bends[i].C;

        index_i=Adsorbates[CurrentSystem][m].Atoms[A].HessianIndex;
        index_j=Adsorbates[CurrentSystem][m].Atoms[B].HessianIndex;
        index_k=Adsorbates[CurrentSystem][m].Atoms[C].HessianIndex;

        posA=Adsorbates[CurrentSystem][m].Atoms[A].Position;
        posB=Adsorbates[CurrentSystem][m].Atoms[B].Position;
        posC=Adsorbates[CurrentSystem][m].Atoms[C].Position;

        ReturnWilsonVectorsBend(posA,posB,posC,&ha,&hb,&hc);

        switch(Dimension)
        {
          case 3:
            if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.z]=ha.z;
          case 2:
            if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.y]=ha.y;
          case 1:
            if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.x]=ha.x;
            break;
        }

        switch(Dimension)
        {
          case 3:
            if(index_j.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_j.z]=hb.z;
          case 2:
            if(index_j.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_j.y]=hb.y;
          case 1:
            if(index_j.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_j.x]=hb.x;
            break;
        }

        switch(Dimension)
        {
          case 3:
            if(index_k.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_k.z]=hc.z;
          case 2:
            if(index_k.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_k.y]=hc.y;
          case 1:
            if(index_k.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_k.x]=hc.x;
            break;
        }

        NumberOfConstraints++;
      }
    }
  }


  // add general angle contraints
  for(m=0;m<NumberOfAngleConstraints[CurrentSystem];m++)
  {
    index_i=AngleConstraints[CurrentSystem][m][0]->HessianIndex;
    index_j=AngleConstraints[CurrentSystem][m][1]->HessianIndex;
    index_k=AngleConstraints[CurrentSystem][m][2]->HessianIndex;

    posA=AngleConstraints[CurrentSystem][m][0]->Position;
    posB=AngleConstraints[CurrentSystem][m][1]->Position;
    posC=AngleConstraints[CurrentSystem][m][2]->Position;

    ReturnWilsonVectorsBend(posA,posB,posC,&ha,&hb,&hc);

    switch(Dimension)
    {
      case 3:
        if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.z]=ha.z;
      case 2:
        if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.y]=ha.y;
      case 1:
        if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.x]=ha.x;
        break;
    }

    switch(Dimension)
    {
      case 3:
        if(index_j.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_j.z]=hb.z;
      case 2:
        if(index_j.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_j.y]=hb.y;
      case 1:
        if(index_j.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_j.x]=hb.x;
        break;
    }

    switch(Dimension)
    {
      case 3:
        if(index_k.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_k.z]=hc.z;
      case 2:
        if(index_k.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_k.y]=hc.y;
      case 1:
        if(index_k.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_k.x]=hc.x;
        break;
    }

    NumberOfConstraints++;
  }




  // add dihedral contraints
  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    Type=Adsorbates[CurrentSystem][m].Type;
    NumberOfTorsions=Components[Type].NumberOfTorsions;
    for(i=0;i<NumberOfTorsions;i++)
    {
      if(Components[Type].TorsionType[i]==FIXED_DIHEDRAL)
      {
        A=Components[Type].Torsions[i].A;
        B=Components[Type].Torsions[i].B;
        C=Components[Type].Torsions[i].C;
        D=Components[Type].Torsions[i].D;

        index_i=Adsorbates[CurrentSystem][m].Atoms[A].HessianIndex;
        index_j=Adsorbates[CurrentSystem][m].Atoms[B].HessianIndex;
        index_k=Adsorbates[CurrentSystem][m].Atoms[C].HessianIndex;
        index_l=Adsorbates[CurrentSystem][m].Atoms[D].HessianIndex;

        posA=Adsorbates[CurrentSystem][m].Atoms[A].Position;
        posB=Adsorbates[CurrentSystem][m].Atoms[B].Position;
        posC=Adsorbates[CurrentSystem][m].Atoms[C].Position;
        posD=Adsorbates[CurrentSystem][m].Atoms[D].Position;

        ReturnWilsonVectorsTorsion(posA,posB,posC,posD,&ha,&hb,&hc,&hd);

        if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.x]=ha.x;
        if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.y]=ha.y;
        if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.z]=ha.z;

        if(index_j.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_j.x]=hb.x;
        if(index_j.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_j.y]=hb.y;
        if(index_j.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_j.z]=hb.z;

        if(index_k.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_k.x]=hc.x;
        if(index_k.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_k.y]=hc.y;
        if(index_k.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_k.z]=hc.z;

        if(index_l.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_l.x]=hd.x;
        if(index_l.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_l.y]=hd.y;
        if(index_l.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_l.z]=hd.z;

        NumberOfConstraints++;
      }
    }
  }

 // add dihedral contraints
  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    Type=Adsorbates[CurrentSystem][m].Type;
    NumberOfImproperTorsions=Components[Type].NumberOfImproperTorsions;
    for(i=0;i<NumberOfImproperTorsions;i++)
    {
      if(Components[Type].ImproperTorsionType[i]==FIXED_IMPROPER_DIHEDRAL)
      {
        A=Components[Type].ImproperTorsions[i].A;
        B=Components[Type].ImproperTorsions[i].B;
        C=Components[Type].ImproperTorsions[i].C;
        D=Components[Type].ImproperTorsions[i].D;

        index_i=Adsorbates[CurrentSystem][m].Atoms[A].HessianIndex;
        index_j=Adsorbates[CurrentSystem][m].Atoms[B].HessianIndex;
        index_k=Adsorbates[CurrentSystem][m].Atoms[C].HessianIndex;
        index_l=Adsorbates[CurrentSystem][m].Atoms[D].HessianIndex;

        posA=Adsorbates[CurrentSystem][m].Atoms[A].Position;
        posB=Adsorbates[CurrentSystem][m].Atoms[B].Position;
        posC=Adsorbates[CurrentSystem][m].Atoms[C].Position;
        posD=Adsorbates[CurrentSystem][m].Atoms[D].Position;

        ReturnWilsonVectorsTorsion(posA,posB,posC,posD,&ha,&hb,&hc,&hd);

        if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.x]=ha.x;
        if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.y]=ha.y;
        if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.z]=ha.z;

        if(index_j.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_j.x]=hb.x;
        if(index_j.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_j.y]=hb.y;
        if(index_j.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_j.z]=hb.z;

        if(index_k.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_k.x]=hc.x;
        if(index_k.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_k.y]=hc.y;
        if(index_k.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_k.z]=hc.z;

        if(index_l.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_l.x]=hd.x;
        if(index_l.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_l.y]=hd.y;
        if(index_l.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_l.z]=hd.z;

        NumberOfConstraints++;
      }
    }
  }


  // add general dihedral angle contraints
  for(m=0;m<NumberOfDihedralConstraints[CurrentSystem];m++)
  {
    index_i=DihedralConstraints[CurrentSystem][m][0]->HessianIndex;
    index_j=DihedralConstraints[CurrentSystem][m][1]->HessianIndex;
    index_k=DihedralConstraints[CurrentSystem][m][2]->HessianIndex;
    index_l=DihedralConstraints[CurrentSystem][m][3]->HessianIndex;

    posA=DihedralConstraints[CurrentSystem][m][0]->Position;
    posB=DihedralConstraints[CurrentSystem][m][1]->Position;
    posC=DihedralConstraints[CurrentSystem][m][2]->Position;
    posD=DihedralConstraints[CurrentSystem][m][3]->Position;

    ReturnWilsonVectorsTorsion(posA,posB,posC,posD,&ha,&hb,&hc,&hd);

    if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.x]=ha.x;
    if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.y]=ha.y;
    if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.z]=ha.z;

    if(index_j.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_j.x]=hb.x;
    if(index_j.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_j.y]=hb.y;
    if(index_j.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_j.z]=hb.z;

    if(index_k.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_k.x]=hc.x;
    if(index_k.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_k.y]=hc.y;
    if(index_k.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_k.z]=hc.z;

    if(index_k.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_l.x]=hd.x;
    if(index_k.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_l.y]=hd.y;
    if(index_k.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_l.z]=hd.z;

    NumberOfConstraints++;
  }

  // add general improper dihedral angle contraints
  for(m=0;m<NumberOfImproperDihedralConstraints[CurrentSystem];m++)
  {
    index_i=ImproperDihedralConstraints[CurrentSystem][m][0]->HessianIndex;
    index_j=ImproperDihedralConstraints[CurrentSystem][m][1]->HessianIndex;
    index_k=ImproperDihedralConstraints[CurrentSystem][m][2]->HessianIndex;
    index_l=ImproperDihedralConstraints[CurrentSystem][m][3]->HessianIndex;

    posA=ImproperDihedralConstraints[CurrentSystem][m][0]->Position;
    posB=ImproperDihedralConstraints[CurrentSystem][m][1]->Position;
    posC=ImproperDihedralConstraints[CurrentSystem][m][2]->Position;
    posD=ImproperDihedralConstraints[CurrentSystem][m][3]->Position;

    ReturnWilsonVectorsTorsion(posA,posB,posC,posD,&ha,&hb,&hc,&hd);

    if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.x]=ha.x;
    if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.y]=ha.y;
    if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.z]=ha.z;

    if(index_j.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_j.x]=hb.x;
    if(index_j.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_j.y]=hb.y;
    if(index_j.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_j.z]=hb.z;

    if(index_k.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_k.x]=hc.x;
    if(index_k.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_k.y]=hc.y;
    if(index_k.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_k.z]=hc.z;

    if(index_l.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_l.x]=hd.x;
    if(index_l.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_l.y]=hd.y;
    if(index_l.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_l.z]=hd.z;

    NumberOfConstraints++;
  }



  // add inversion bend contraints
  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    Type=Adsorbates[CurrentSystem][m].Type;
    NumberOfInversionBends=Components[Type].NumberOfInversionBends;
    for(i=0;i<NumberOfInversionBends;i++)
    {
      if(Components[Type].InversionBendType[i]==FIXED_INVERSION_BEND)
      {
        A=Components[Type].InversionBends[i].A;
        B=Components[Type].InversionBends[i].B;
        C=Components[Type].InversionBends[i].C;
        D=Components[Type].InversionBends[i].D;

        index_i=Adsorbates[CurrentSystem][m].Atoms[A].HessianIndex;
        index_j=Adsorbates[CurrentSystem][m].Atoms[B].HessianIndex;
        index_k=Adsorbates[CurrentSystem][m].Atoms[C].HessianIndex;
        index_l=Adsorbates[CurrentSystem][m].Atoms[D].HessianIndex;

        posA=Adsorbates[CurrentSystem][m].Atoms[A].Position;
        posB=Adsorbates[CurrentSystem][m].Atoms[B].Position;
        posC=Adsorbates[CurrentSystem][m].Atoms[C].Position;
        posD=Adsorbates[CurrentSystem][m].Atoms[D].Position;

        ReturnWilsonVectorsInversionBend(posA,posB,posC,posD,&ha,&hb,&hc,&hd);

        if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.x]=ha.x;
        if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.y]=ha.y;
        if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.z]=ha.z;

        if(index_j.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_j.x]=hb.x;
        if(index_j.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_j.y]=hb.y;
        if(index_j.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_j.z]=hb.z;

        if(index_k.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_k.x]=hc.x;
        if(index_k.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_k.y]=hc.y;
        if(index_k.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_k.z]=hc.z;

        if(index_l.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_l.x]=hd.x;
        if(index_l.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_l.y]=hd.y;
        if(index_l.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_l.z]=hd.z;

        NumberOfConstraints++;
      }
    }
  }

  // add general inversion-bend angle contraints
  for(m=0;m<NumberOfInversionBendConstraints[CurrentSystem];m++)
  {
    index_i=InversionBendConstraints[CurrentSystem][m][0]->HessianIndex;
    index_j=InversionBendConstraints[CurrentSystem][m][1]->HessianIndex;
    index_k=InversionBendConstraints[CurrentSystem][m][2]->HessianIndex;
    index_l=InversionBendConstraints[CurrentSystem][m][3]->HessianIndex;

    posA=InversionBendConstraints[CurrentSystem][m][0]->Position;
    posB=InversionBendConstraints[CurrentSystem][m][1]->Position;
    posC=InversionBendConstraints[CurrentSystem][m][2]->Position;
    posD=InversionBendConstraints[CurrentSystem][m][3]->Position;

    ReturnWilsonVectorsInversionBend(posA,posB,posC,posD,&ha,&hb,&hc,&hd);

    if(index_i.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.x]=ha.x;
    if(index_i.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.y]=ha.y;
    if(index_i.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_i.z]=ha.z;

    if(index_j.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_j.x]=hb.x;
    if(index_j.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_j.y]=hb.y;
    if(index_j.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_j.z]=hb.z;

    if(index_k.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_k.x]=hc.x;
    if(index_k.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_k.y]=hc.y;
    if(index_k.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_k.z]=hc.z;

    if(index_l.x>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_l.x]=hd.x;
    if(index_l.y>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_l.y]=hd.y;
    if(index_l.z>=0) PrimativeUnitvectors.element[NumberOfConstraints][index_l.z]=hd.z;

    NumberOfConstraints++;
  }


  // create orthonormal basis using simple Gram-Schmidt
  NumberOfRemovedModes=0;
  for(i=0;i<NumberOfConstraints;i++)
  {
    for(j=0;j<NumberOfRemovedModes;j++)
    {
      dot_produkt=0.0;
      for(k=0;k<n;k++)
        dot_produkt+=PrimativeUnitvectors.element[i][k]*PrimativeUnitvectors.element[j][k];
      for(k=0;k<n;k++)
        PrimativeUnitvectors.element[i][k]-=dot_produkt*PrimativeUnitvectors.element[j][k];
    }
    dot_produkt=0.0;
    for(k=0;k<n;k++)
      dot_produkt+=SQR(PrimativeUnitvectors.element[i][k]);
    if(dot_produkt>1e-10)
    {
      norm=sqrt(dot_produkt);
      for(k=0;k<n;k++)
        PrimativeUnitvectors.element[NumberOfRemovedModes][k]=PrimativeUnitvectors.element[i][k]/norm;
      NumberOfRemovedModes++;
    }
  }

  if(NumberOfRemovedModes==0) return 0;

  TransposeRealMatrix(PrimativeUnitvectorsTranspose,PrimativeUnitvectors);
  MultiplyRealMatrix(ProjectionOperator,PrimativeUnitvectorsTranspose,PrimativeUnitvectors);

  for(i=0;i<n;i++)
  {
    for(j=0;j<n;j++)
      ProjectionOperator.element[i][j]=-ProjectionOperator.element[i][j];
    ProjectionOperator.element[i][i]+=1.0;
  }

  if(ComputeGradient)
  {
    MultiplyRealMatrixVector(ModifiedGradient,ProjectionOperator,Gradient);
    for(i=0;i<n;i++)
      Gradient[i]=ModifiedGradient[i];
  }

  if(ComputeHessian)
  {
    MultiplyRealMatrix(ProjectionOperator2,HessianMatrix,ProjectionOperator);
    MultiplyRealMatrix(HessianMatrix,ProjectionOperator,ProjectionOperator2);
  }

  DeleteRealMatrix(PrimativeUnitvectors);
  DeleteRealMatrix(PrimativeUnitvectorsTranspose);
  DeleteRealMatrix(ProjectionOperator);
  DeleteRealMatrix(ProjectionOperator2);
  free(ModifiedGradient);

  return NumberOfRemovedModes;
}

