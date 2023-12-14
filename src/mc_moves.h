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

#ifndef MC_MOVES_H
#define MC_MOVES_H

#include "simulation.h"
#include "molecule.h"
#include "framework.h"
#include "input.h"

// delta energies
// Modified by Ambroise de Izarra
//-------------------------------------------------------------------
// switched from static to extern for use in alchemical transformation.
int ComputeNewPolarizationEnergy(int New,int excl_ads,int excl_cation);
extern REAL UDeltaPolarization;
extern REAL *UHostPolarizationNew;
extern REAL *UAdsorbatePolarizationNew;
extern REAL *UCationPolarizationNew;
extern REAL *UPolarizationNew;

extern REAL *UHostBackPolarizationNew;
extern REAL *UAdsorbateBackPolarizationNew;
extern REAL *UCationBackPolarizationNew;
extern REAL *UBackPolarizationNew;
//-------------------------------------------------------------------

extern REAL *UHostVDWDelta;
extern REAL *UHostChargeChargeRealDelta;
extern REAL *UHostChargeBondDipoleRealDelta;
extern REAL *UHostBondDipoleBondDipoleRealDelta;
extern REAL *UHostFourierDelta;

extern REAL *UCationVDWDelta;
extern REAL *UAdsorbateVDWDelta;

extern REAL *UCationChargeChargeRealDelta;
extern REAL *UAdsorbateChargeChargeRealDelta;

extern REAL *UCationChargeBondDipoleRealDelta;
extern REAL *UAdsorbateChargeBondDipoleRealDelta;

extern REAL *UCationBondDipoleBondDipoleRealDelta;
extern REAL *UAdsorbateBondDipoleBondDipoleRealDelta;

extern REAL *UCationFourierDelta;
extern REAL *UAdsorbateFourierDelta;


extern VECTOR **MaximumTranslation;
extern VECTOR **MaximumTranslationInPlane;
extern VECTOR **MaximumRotation;

extern VECTOR **TrialPosition;
extern VECTOR **TrialAnisotropicPosition;
extern VECTOR ***RXMCTrialAnisotropicPositions;
extern VECTOR ***RXMCTrialAnisotropicPositions2;
extern VECTOR **ElectricFieldAtTrialPosition;
extern VECTOR **InducedDipoleAtTrialPosition;
extern REAL TargetAccRatioTranslation;
extern REAL TargetAccRatioRotation;

extern REAL *MaximumVolumeChange;
extern REAL TargetAccRatioVolumeChange;

extern REAL_MATRIX3x3 *MaximumBoxShapeChange;
extern REAL TargetAccRatioBoxShapeChange;

extern REAL *MaximumGibbsVolumeChange;
extern REAL TargetAccRatioGibbsVolumeChange;

extern REAL **FrameworkMaximumTranslation;
extern VECTOR **FrameworkMaximumShiftTranslation;

extern int ParallelMolFractionComponentA;
extern int ParallelMolFractionComponentB;

extern int NumberOfHybridNVESteps;
extern int NumberOfHybridNPHSteps;
extern int NumberOfHybridNPHPRSteps;


//----------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------
// CFC-RXMC Parameters
//----------------------------------------------------------------------------------------

extern REAL **MaximumReactionLambdaChange;
extern REAL TargetAccRatioReactionLambdaChange;

extern int CFWangLandauEvery;
extern REAL **MaximumCFLambdaChange;
extern REAL **MaximumCFWidomLambdaChange;
extern REAL **MaximumCBCFLambdaChange;
extern REAL TargetAccRatioLambdaChange;

//----------------------------------------------------------------------------------------

enum {CF_INSERT_MOVE,CF_DELETE_MOVE,CF_MOVE};

void InitializeMCMovesStatisticsAllSystems(void);


REAL ComputeEnergyOfFractionalMoleculesAdsorbates();
REAL ComputeEnergyOfFractionalMoleculesCations();

void TranslationMove(void);
void RandomTranslationMove(void);
void RotationMove(void);
void RandomRotationMove(void);
void PartialReinsertionMove(void);
void ReinsertionMove(void);
void ReinsertionInPlaceMove(void);
void ReinsertionInPlaneMove(void);
void SwapAddMove(void);
void SwapRemoveMove(void);
void IdentityChangeMove(void);

void WidomMove(void);
void CFWidomLambaMove(void);
int CFWidomLambaAdsorbateMove(void);
int CFWidomLambaCationMove(void);
void OptimizeCFWidomAcceptence(void);
void PrintCFWidomLambdaStatistics(FILE *FilePtr);

void GibbsWidomMove(void);
void GibbsParticleTransferMove(void);
void GibbsIdentityChangeMove(void);

int TranslationMoveAdsorbate(void);
int RandomTranslationMoveAdsorbate(void);
int RotationMoveAdsorbate(void);
int SwapAddAdsorbateMove(void);
int SwapRemoveAdsorbateMove(void);
int ReinsertionAdsorbateMove(void);
int ReinsertionInPlaceAdsorbateMove(void);
int ReinsertionInPlaneAdsorbateMove(void);
int IdentityChangeAdsorbateMove(void);
// Added by Ambroise de Izarra
//-------------------------------------------------------------------
int AlchemicalChangeAdsorbateMove(void);
int WidomOsmostatCalculation(void);
//-------------------------------------------------------------------
int PartialReinsertionAdsorbateMove(void);

void OptimizeTranslationAcceptence(void);
void OptimizeRotationAcceptence(void);
void OptimizeTranslationInPlaneAcceptence(void);
void OptimizeGibbsVolumeChangeAcceptence(void);
void OptimizeVolumeChangeAcceptence(void);
void OptimizeFrameworkChangeAcceptence(void);
void OptimizeFrameworkShiftAcceptence(void);

int TranslationMoveCation(void);
int RandomTranslationMoveCation(void);
int RotationMoveCation(void);
int SwapAddCationMove(void);
int SwapRemoveCationMove(void);
int ReinsertionCationMove(void);
int ReinsertionInPlaceCationMove(void);
int ReinsertionInPlaneCationMove(void);
int PartialReinsertionCationMove(void);
int IdentityChangeCationMove(void);

int ParallelTemperingMove(void);
int HyperParallelTemperingMove(void);
int ParallelMolFractionMove(void);
int ChiralInversionMove(void);
int VolumeMove(void);
int BoxShapeChangeMove(void);
int GibbsParticleTransferAdsorbateMove(void);
int GibbsParticleTransferCationMove(void);
int GibbsVolumeMove(void);
int GibbsIdentityChangeAdsorbateMove(void);
int FrameworkChangeMove(void);
int FrameworkShiftMove(void);
int ExchangeFractionalParticleMove(void);

void PrintGibbsSwapStatistics(FILE *FilePtr);
void PrintTranslationStatistics(FILE *FilePtr);
void PrintRandomTranslationStatistics(FILE *FilePtr);
void PrintRotationStatistics(FILE *FilePtr);
void PrintRandomRotationStatistics(FILE *FilePtr);
void PrintSwapAddStatistics(FILE *FilePtr);
void PrintSwapRemoveStatistics(FILE *FilePtr);
void PrintReinsertionStatistics(FILE *FilePtr);
void PrintReinsertionInPlaneStatistics(FILE *FilePtr);
void PrintReinsertionInPlaceStatistics(FILE *FilePtr);
void PrintPartialReinsertionStatistics(FILE *FilePtr);
void PrintIdentityChangeStatistics(FILE *FilePtr);
void PrintParallelTemperingStatistics(FILE *FilePtr);
void PrintHyperParallelTemperingStatistics(FILE *FilePtr);

// Added by Ambroise de Izarra
//-------------------------------------------------------------------
void PrintAlchemicalChangeStatistics(FILE *FilePtr);
void PrintwidomOsmostatStatistics(FILE *FilePtr);

extern REAL **AlchemicalChangeAttempts;
extern REAL (**AlchemicalChangeAccepted)[2];

extern REAL **WidomOsmostat;
//-------------------------------------------------------------------
void PrintParallelMolFractionStatistics(FILE *FilePtr);
void PrintChiralInversionStatistics(FILE *FilePtr);
void PrintVolumeChangeStatistics(FILE *FilePtr);
void PrintBoxShapeChangeStatistics(FILE *FilePtr);
void PrintGibbsVolumeChangeStatistics(FILE *FilePtr);
void PrintFrameworkStatistics(FILE *FilePtr);
void PrintFrameworkShiftStatistics(FILE *FilePtr);
void PrintGibbsIdentityChangeStatistics(FILE *FilePtr);
void PrintCFSwapLambdaStatistics(FILE *FilePtr);
void PrintCBCFSwapLambdaStatistics(FILE *FilePtr);
void PrintCFGibbsLambdaStatistics(FILE *FilePtr);
void PrintCBCFGibbsLambdaStatistics(FILE *FilePtr);
void PrintRXMCStatistics(FILE *FilePtr);
void PrintExchangeFractionalParticleStatistics(FILE *FilePtr);
void PrintCFGibbsLambdaChangeStatistics(FILE *FilePtr);
void PrintCFGibbsSwapFractionalMoleculeToOtherBoxStatistics(FILE *FilePtr);
void PrintCFGibbsFractionalToIntegerStatistics(FILE *FilePtr);
void PrintCFGibbsWidomStatistics(FILE *FilePtr);

REAL WidomAdsorbateMove(void);
REAL WidomCationMove(void);

void HybridNVEMove(void);
void PrintHybridNVEStatistics(FILE *FilePtr);

void HybridNPHMove(void);
void PrintHybridNPHStatistics(FILE *FilePtr);

void HybridNPHPRMove(void);
void PrintHybridNPHPRStatistics(FILE *FilePtr);

void SurfaceAreaMove(void);

void CFWangLandauIteration(int Switch);
void CFSwapLambaMove(void);
int CFSwapLambaAdsorbateMove(void);
int CFSwapLambaCationMove(void);
void OptimizeCFLambdaAcceptence(void);
void CBCFSwapLambaMove(void);
int CBCFSwapLambaAdsorbateMove(void);
int CBCFSwapLambaCationMove(void);
void OptimizeCBCFLambdaChangeAcceptence(void);

void CFGibbsParticleTransferMove(void);
int CFGibbsParticleTransferAdsorbateMove(void);
int CFGibbsParticleTransferCationMove(void);
void OptimizeCFGibbsLambdaAcceptence(void);
void CBCFGibbsParticleTransferMove(void);
int CBCFGibbsParticleTransferAdsorbateMove(void);
int CBCFGibbsParticleTransferCationMove(void);
void OptimizeCBCFGibbsLambdaChangeAcceptence(void);

void CFCRXMCLambdaChangeMove(void);
void OptimizeRXMCLambdaChangeAcceptence(void);
void CFRXMCWangLandauIteration(int Switch);

int CFGibbsSwapFractionalMoleculeToOtherBoxMove(void);
int CFGibbsLambdaChangeMove(void);
void SampleLambdaHistogram();
void ClearLambdaHistogram();
int CFGibbsFractionalToIntegerMove(void);
void OptimizeCFGibbsLambdaChangeAcceptence(void);

void WriteRestartMcMoves(FILE *FilePtr);
void AllocateMCMovesMemory(void);
void ReadRestartMcMoves(FILE *FilePtr);


#endif
