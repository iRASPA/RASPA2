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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "constants.h"
#include "cbmc.h"
#include "simulation.h"
#include "utils.h"
#include "molecule.h"
#include "potentials.h"
#include "internal_energy.h"
#include "inter_energy.h"
#include "framework.h"
#include "framework_energy.h"
#include "alchemical_transformation.h"
#include "ewald.h"
#include "sample.h"
#include "utils.h"
#include "grids.h"
#include "rigid.h"
#include "output.h"
#include "statistics.h"
#include "spacegroup.h"
#include "thermo_baro_stats.h"
#include "molecular_dynamics.h"
#include "warnings.h"
#include "integration.h"
#include "equations_of_state.h"
#include "mc_moves.h"

long long *BlockCycle;
REAL **BlockWeightedCount;
REAL **BlockCount;
int Block;

REAL MeasureLambdaBelow;

// average energies
static REAL **UHostHostAccumulated;
static REAL **UAdsorbateAdsorbateAccumulated;
static REAL **UCationCationAccumulated;
static REAL **UHostAdsorbateAccumulated;
static REAL **UHostCationAccumulated;
static REAL **UAdsorbateCationAccumulated;

static REAL **UHostHostVDWAccumulated;
static REAL **UAdsorbateAdsorbateVDWAccumulated;
static REAL **UCationCationVDWAccumulated;
static REAL **UHostAdsorbateVDWAccumulated;
static REAL **UHostCationVDWAccumulated;
static REAL **UAdsorbateCationVDWAccumulated;

static REAL **UHostHostCoulombAccumulated;
static REAL **UAdsorbateAdsorbateCoulombAccumulated;
static REAL **UCationCationCoulombAccumulated;
static REAL **UHostAdsorbateCoulombAccumulated;
static REAL **UHostCationCoulombAccumulated;
static REAL **UAdsorbateCationCoulombAccumulated;

static REAL **UHostBondAccumulated;
static REAL **UHostUreyBradleyAccumulated;
static REAL **UHostBendAccumulated;
static REAL **UHostInversionBendAccumulated;
static REAL **UHostTorsionAccumulated;
static REAL **UHostImproperTorsionAccumulated;
static REAL **UHostBondBondAccumulated;
static REAL **UHostBendBendAccumulated;
static REAL **UHostBondBendAccumulated;
static REAL **UHostBondTorsionAccumulated;
static REAL **UHostBendTorsionAccumulated;

static REAL **UAdsorbateBondAccumulated;
static REAL **UAdsorbateUreyBradleyAccumulated;
static REAL **UAdsorbateBondBondAccumulated;
static REAL **UAdsorbateBendAccumulated;
static REAL **UAdsorbateInversionBendAccumulated;
static REAL **UAdsorbateTorsionAccumulated;
static REAL **UAdsorbateImproperTorsionAccumulated;
static REAL **UAdsorbateBondBondAccumulated;
static REAL **UAdsorbateBendBendAccumulated;
static REAL **UAdsorbateBondBendAccumulated;
static REAL **UAdsorbateBondTorsionAccumulated;
static REAL **UAdsorbateBendTorsionAccumulated;
static REAL **UAdsorbateIntraVDWAccumulated;
static REAL **UAdsorbateIntraChargeChargeAccumulated;
static REAL **UAdsorbateIntraChargeBondDipoleAccumulated;
static REAL **UAdsorbateIntraBondDipoleBondDipoleAccumulated;

static REAL **UCationBondAccumulated;
static REAL **UCationUreyBradleyAccumulated;
static REAL **UCationBendAccumulated;
static REAL **UCationInversionBendAccumulated;
static REAL **UCationTorsionAccumulated;
static REAL **UCationImproperTorsionAccumulated;
static REAL **UCationBondBondAccumulated;
static REAL **UCationBendBendAccumulated;
static REAL **UCationBondBendAccumulated;
static REAL **UCationBondTorsionAccumulated;
static REAL **UCationBendTorsionAccumulated;
static REAL **UCationIntraVDWAccumulated;
static REAL **UCationIntraChargeChargeAccumulated;
static REAL **UCationIntraChargeBondDipoleAccumulated;
static REAL **UCationIntraBondDipoleBondDipoleAccumulated;

static REAL **UHostPolarizationAccumulated;
static REAL **UAdsorbatePolarizationAccumulated;
static REAL **UCationPolarizationAccumulated;
static REAL **UHostBackPolarizationAccumulated;
static REAL **UAdsorbateBackPolarizationAccumulated;
static REAL **UCationBackPolarizationAccumulated;

static REAL **UTailCorrectionAccumulated;

static REAL **UDistanceConstraintsAccumulated;
static REAL **UAngleConstraintsAccumulated;
static REAL **UDihedralConstraintsAccumulated;
static REAL **UInversionBendConstraintsAccumulated;
static REAL **UOutOfPlaneDistanceConstraintsAccumulated;
static REAL **UExclusionConstraintsAccumulated;

static REAL **UTotalAccumulated;

static REAL ***NumberOfIntegerMoleculesPerComponentAccumulated;
static REAL ***NumberOfFractionalMoleculesPerComponentAccumulated;
static REAL ***NumberOfExcessMoleculesPerComponentAccumulated;
static REAL ***DensityPerComponentAccumulated;
static REAL **TotalEnergyTimesNumberOfMoleculesAccumulated;
static REAL ***TotalEnergyTimesNumberOfMoleculesPerComponentAccumulated;
static REAL ***HostAdsorbateEnergyTimesNumberOfMoleculesAccumulated;
static REAL ***AdsorbateAdsorbateEnergyTimesNumberOfMoleculesAccumulated;


static VECTOR **TotalSystemDipoleAccumulated;
static VECTOR **TotalSystemDipoleSquaredAccumulated;
static REAL **TotalSystemNormDipoleAccumulated;
static REAL **TotalSystemNormDipoleSquaredAccumulated;

static REAL **NumberOfIntegerMoleculesAccumulated;
static REAL **NumberOfMoleculesSquaredAccumulated;
static REAL ****NumberOfMoleculesPerComponentSquaredAccumulated;
static REAL **DensityAccumulated;
static REAL **InverseDensityAccumulated;
static REAL **GibbsInverseDensityAccumulated;
static REAL ***GibbsInverseDensityPerComponentAccumulated;

static VECTOR **BoxAccumulated;
static REAL **BoxAXAccumulated;
static REAL **BoxAYAccumulated;
static REAL **BoxAZAccumulated;
static REAL **BoxBXAccumulated;
static REAL **BoxBYAccumulated;
static REAL **BoxBZAccumulated;
static REAL **BoxCXAccumulated;
static REAL **BoxCYAccumulated;
static REAL **BoxCZAccumulated;
static VECTOR **BoxLengthAccumulated;
static REAL **AlphaAngleAccumulated;
static REAL **BetaAngleAccumulated;
static REAL **GammaAngleAccumulated;
static REAL **VolumeAccumulated;
static REAL **VolumeSquaredAccumulated;

static REAL **TotalEnergyAccumulated;
static REAL **TotalEnergySquaredAccumulated;
static REAL **EnthalpyAccumulated;
static REAL **EnthalpySquaredAccumulated;
static REAL **EnthalpyTimesVolumeAccumulated;
static REAL **EnthalpyTimesEnergyAccumulated;

static REAL **TemperatureAccumulated;
static REAL **TemperatureCellAccumulated;
static REAL **TemperatureTranslationAccumulated;
static REAL **TemperatureRotationAccumulated;
static REAL **TemperatureRotationAdsorbateAccumulated;
static REAL **TemperatureTranslationAdsorbateAccumulated;

static REAL **TemperatureAdsorbatesAccumulated;
static REAL **TemperatureCationsAccumulated;
static REAL **TemperatureFrameworkAccumulated;

static REAL **MolecularPressureAccumulated;
static REAL_MATRIX3x3 **MolecularStressTensorAccumulated;
static REAL **PressureIdealGasPartAccumulated;
static REAL **PressureExcessPartAccumulated;
static REAL **PressureTailCorrectionAccumulated;
static REAL **PressureAccumulated;
static REAL **UNoseHooverAccumulated;

static REAL **HeatOfVaporization;
static REAL **EnergyPerMolecule;
static REAL **VolumePerMolecule;
static REAL **CompressibilityAccumulated;

static REAL_MATRIX9x9 **BornTermAccumulated;
static REAL_MATRIX9x9 **ConfigurationalStressFluctuationTermAccumulated;
static REAL_MATRIX3x3 **ConfigurationalStressTensorAccumulated;
static REAL_MATRIX3x3 **StressTensorAccumulated;

REAL ***WidomRosenbluthFactorAccumulated;
REAL ***WidomIdealGasAccumulated;
REAL ***WidomRosenbluthFactorCount;
REAL ***WidomEnergyDifferenceAccumulated;
REAL ***WidomEnergyFrameworkAccumulated;
REAL ***WidomEnergyFrameworkCount;

REAL ***GibbsWidomRosenbluthFactorAccumulated;
REAL ***GibbsWidomIdealGasAccumulated;
REAL ***GibbsWidomRosenbluthFactorCount;
REAL ***GibbsWidomEnergyDifferenceAccumulated;
REAL ***GibbsWidomEnergyFrameworkAccumulated;
REAL ***GibbsWidomEnergyFrameworkCount;

REAL **SurfaceAreaFrameworkAccumulated;
REAL ***SurfaceAreaFrameworksAccumulated;
REAL **SurfaceAreaCationsAccumulated;
REAL **SurfaceAreaCount;

VECTOR ***PrincipleMomentsOfInertiaAccumulated;
REAL ***PrincipleMomentsOfInertiaCount;

void AddBornTermToAverages(void)
{

  // BornTerm
  BornTermAccumulated[CurrentSystem][Block].xxxx+=BornTerm[CurrentSystem].xxxx/Volume[CurrentSystem];
  BornTermAccumulated[CurrentSystem][Block].yxxx+=BornTerm[CurrentSystem].yxxx/Volume[CurrentSystem];
  BornTermAccumulated[CurrentSystem][Block].zxxx+=BornTerm[CurrentSystem].zxxx/Volume[CurrentSystem];
  BornTermAccumulated[CurrentSystem][Block].xxyx+=BornTerm[CurrentSystem].xxyx/Volume[CurrentSystem];
  BornTermAccumulated[CurrentSystem][Block].yxyx+=BornTerm[CurrentSystem].yxyx/Volume[CurrentSystem];
  BornTermAccumulated[CurrentSystem][Block].zxyx+=BornTerm[CurrentSystem].zxyx/Volume[CurrentSystem];
  BornTermAccumulated[CurrentSystem][Block].xxzx+=BornTerm[CurrentSystem].xxzx/Volume[CurrentSystem];
  BornTermAccumulated[CurrentSystem][Block].yxzx+=BornTerm[CurrentSystem].yxzx/Volume[CurrentSystem];
  BornTermAccumulated[CurrentSystem][Block].zxzx+=BornTerm[CurrentSystem].zxzx/Volume[CurrentSystem];

  BornTermAccumulated[CurrentSystem][Block].xyxx+=BornTerm[CurrentSystem].xyxx/Volume[CurrentSystem];
  BornTermAccumulated[CurrentSystem][Block].yyxx+=BornTerm[CurrentSystem].yyxx/Volume[CurrentSystem];
  BornTermAccumulated[CurrentSystem][Block].zyxx+=BornTerm[CurrentSystem].zyxx/Volume[CurrentSystem];
  BornTermAccumulated[CurrentSystem][Block].xyyx+=BornTerm[CurrentSystem].xyyx/Volume[CurrentSystem];
  BornTermAccumulated[CurrentSystem][Block].yyyx+=BornTerm[CurrentSystem].yyyx/Volume[CurrentSystem];
  BornTermAccumulated[CurrentSystem][Block].zyyx+=BornTerm[CurrentSystem].zyyx/Volume[CurrentSystem];
  BornTermAccumulated[CurrentSystem][Block].xyzx+=BornTerm[CurrentSystem].xyzx/Volume[CurrentSystem];
  BornTermAccumulated[CurrentSystem][Block].yyzx+=BornTerm[CurrentSystem].yyzx/Volume[CurrentSystem];
  BornTermAccumulated[CurrentSystem][Block].zyzx+=BornTerm[CurrentSystem].zyzx/Volume[CurrentSystem];

  BornTermAccumulated[CurrentSystem][Block].xzxx+=BornTerm[CurrentSystem].xzxx/Volume[CurrentSystem];
  BornTermAccumulated[CurrentSystem][Block].yzxx+=BornTerm[CurrentSystem].yzxx/Volume[CurrentSystem];
  BornTermAccumulated[CurrentSystem][Block].zzxx+=BornTerm[CurrentSystem].zzxx/Volume[CurrentSystem];
  BornTermAccumulated[CurrentSystem][Block].xzyx+=BornTerm[CurrentSystem].xzyx/Volume[CurrentSystem];
  BornTermAccumulated[CurrentSystem][Block].yzyx+=BornTerm[CurrentSystem].yzyx/Volume[CurrentSystem];
  BornTermAccumulated[CurrentSystem][Block].zzyx+=BornTerm[CurrentSystem].zzyx/Volume[CurrentSystem];
  BornTermAccumulated[CurrentSystem][Block].xzzx+=BornTerm[CurrentSystem].xzzx/Volume[CurrentSystem];
  BornTermAccumulated[CurrentSystem][Block].yzzx+=BornTerm[CurrentSystem].yzzx/Volume[CurrentSystem];
  BornTermAccumulated[CurrentSystem][Block].zzzx+=BornTerm[CurrentSystem].zzzx/Volume[CurrentSystem];

  BornTermAccumulated[CurrentSystem][Block].xxxy+=BornTerm[CurrentSystem].xxxy/Volume[CurrentSystem];
  BornTermAccumulated[CurrentSystem][Block].yxxy+=BornTerm[CurrentSystem].yxxy/Volume[CurrentSystem];
  BornTermAccumulated[CurrentSystem][Block].zxxy+=BornTerm[CurrentSystem].zxxy/Volume[CurrentSystem];
  BornTermAccumulated[CurrentSystem][Block].xxyy+=BornTerm[CurrentSystem].xxyy/Volume[CurrentSystem];
  BornTermAccumulated[CurrentSystem][Block].yxyy+=BornTerm[CurrentSystem].yxyy/Volume[CurrentSystem];
  BornTermAccumulated[CurrentSystem][Block].zxyy+=BornTerm[CurrentSystem].zxyy/Volume[CurrentSystem];
  BornTermAccumulated[CurrentSystem][Block].xxzy+=BornTerm[CurrentSystem].xxzy/Volume[CurrentSystem];
  BornTermAccumulated[CurrentSystem][Block].yxzy+=BornTerm[CurrentSystem].yxzy/Volume[CurrentSystem];
  BornTermAccumulated[CurrentSystem][Block].zxzy+=BornTerm[CurrentSystem].zxzy/Volume[CurrentSystem];

  BornTermAccumulated[CurrentSystem][Block].xyxy+=BornTerm[CurrentSystem].xyxy/Volume[CurrentSystem];
  BornTermAccumulated[CurrentSystem][Block].yyxy+=BornTerm[CurrentSystem].yyxy/Volume[CurrentSystem];
  BornTermAccumulated[CurrentSystem][Block].zyxy+=BornTerm[CurrentSystem].zyxy/Volume[CurrentSystem];
  BornTermAccumulated[CurrentSystem][Block].xyyy+=BornTerm[CurrentSystem].xyyy/Volume[CurrentSystem];
  BornTermAccumulated[CurrentSystem][Block].yyyy+=BornTerm[CurrentSystem].yyyy/Volume[CurrentSystem];
  BornTermAccumulated[CurrentSystem][Block].zyyy+=BornTerm[CurrentSystem].zyyy/Volume[CurrentSystem];
  BornTermAccumulated[CurrentSystem][Block].xyzy+=BornTerm[CurrentSystem].xyzy/Volume[CurrentSystem];
  BornTermAccumulated[CurrentSystem][Block].yyzy+=BornTerm[CurrentSystem].yyzy/Volume[CurrentSystem];
  BornTermAccumulated[CurrentSystem][Block].zyzy+=BornTerm[CurrentSystem].zyzy/Volume[CurrentSystem];

  BornTermAccumulated[CurrentSystem][Block].xzxy+=BornTerm[CurrentSystem].xzxy/Volume[CurrentSystem];
  BornTermAccumulated[CurrentSystem][Block].yzxy+=BornTerm[CurrentSystem].yzxy/Volume[CurrentSystem];
  BornTermAccumulated[CurrentSystem][Block].zzxy+=BornTerm[CurrentSystem].zzxy/Volume[CurrentSystem];
  BornTermAccumulated[CurrentSystem][Block].xzyy+=BornTerm[CurrentSystem].xzyy/Volume[CurrentSystem];
  BornTermAccumulated[CurrentSystem][Block].yzyy+=BornTerm[CurrentSystem].yzyy/Volume[CurrentSystem];
  BornTermAccumulated[CurrentSystem][Block].zzyy+=BornTerm[CurrentSystem].zzyy/Volume[CurrentSystem];
  BornTermAccumulated[CurrentSystem][Block].xzzy+=BornTerm[CurrentSystem].xzzy/Volume[CurrentSystem];
  BornTermAccumulated[CurrentSystem][Block].yzzy+=BornTerm[CurrentSystem].yzzy/Volume[CurrentSystem];
  BornTermAccumulated[CurrentSystem][Block].zzzy+=BornTerm[CurrentSystem].zzzy/Volume[CurrentSystem];

  BornTermAccumulated[CurrentSystem][Block].xxxz+=BornTerm[CurrentSystem].xxxz/Volume[CurrentSystem];
  BornTermAccumulated[CurrentSystem][Block].yxxz+=BornTerm[CurrentSystem].yxxz/Volume[CurrentSystem];
  BornTermAccumulated[CurrentSystem][Block].zxxz+=BornTerm[CurrentSystem].zxxz/Volume[CurrentSystem];
  BornTermAccumulated[CurrentSystem][Block].xxyz+=BornTerm[CurrentSystem].xxyz/Volume[CurrentSystem];
  BornTermAccumulated[CurrentSystem][Block].yxyz+=BornTerm[CurrentSystem].yxyz/Volume[CurrentSystem];
  BornTermAccumulated[CurrentSystem][Block].zxyz+=BornTerm[CurrentSystem].zxyz/Volume[CurrentSystem];
  BornTermAccumulated[CurrentSystem][Block].xxzz+=BornTerm[CurrentSystem].xxzz/Volume[CurrentSystem];
  BornTermAccumulated[CurrentSystem][Block].yxzz+=BornTerm[CurrentSystem].yxzz/Volume[CurrentSystem];
  BornTermAccumulated[CurrentSystem][Block].zxzz+=BornTerm[CurrentSystem].zxzz/Volume[CurrentSystem];

  BornTermAccumulated[CurrentSystem][Block].xyxz+=BornTerm[CurrentSystem].xyxz/Volume[CurrentSystem];
  BornTermAccumulated[CurrentSystem][Block].yyxz+=BornTerm[CurrentSystem].yyxz/Volume[CurrentSystem];
  BornTermAccumulated[CurrentSystem][Block].zyxz+=BornTerm[CurrentSystem].zyxz/Volume[CurrentSystem];
  BornTermAccumulated[CurrentSystem][Block].xyyz+=BornTerm[CurrentSystem].xyyz/Volume[CurrentSystem];
  BornTermAccumulated[CurrentSystem][Block].yyyz+=BornTerm[CurrentSystem].yyyz/Volume[CurrentSystem];
  BornTermAccumulated[CurrentSystem][Block].zyyz+=BornTerm[CurrentSystem].zyyz/Volume[CurrentSystem];
  BornTermAccumulated[CurrentSystem][Block].xyzz+=BornTerm[CurrentSystem].xyzz/Volume[CurrentSystem];
  BornTermAccumulated[CurrentSystem][Block].yyzz+=BornTerm[CurrentSystem].yyzz/Volume[CurrentSystem];
  BornTermAccumulated[CurrentSystem][Block].zyzz+=BornTerm[CurrentSystem].zyzz/Volume[CurrentSystem];

  BornTermAccumulated[CurrentSystem][Block].xzxz+=BornTerm[CurrentSystem].xzxz/Volume[CurrentSystem];
  BornTermAccumulated[CurrentSystem][Block].yzxz+=BornTerm[CurrentSystem].yzxz/Volume[CurrentSystem];
  BornTermAccumulated[CurrentSystem][Block].zzxz+=BornTerm[CurrentSystem].zzxz/Volume[CurrentSystem];
  BornTermAccumulated[CurrentSystem][Block].xzyz+=BornTerm[CurrentSystem].xzyz/Volume[CurrentSystem];
  BornTermAccumulated[CurrentSystem][Block].yzyz+=BornTerm[CurrentSystem].yzyz/Volume[CurrentSystem];
  BornTermAccumulated[CurrentSystem][Block].zzyz+=BornTerm[CurrentSystem].zzyz/Volume[CurrentSystem];
  BornTermAccumulated[CurrentSystem][Block].xzzz+=BornTerm[CurrentSystem].xzzz/Volume[CurrentSystem];
  BornTermAccumulated[CurrentSystem][Block].yzzz+=BornTerm[CurrentSystem].yzzz/Volume[CurrentSystem];
  BornTermAccumulated[CurrentSystem][Block].zzzz+=BornTerm[CurrentSystem].zzzz/Volume[CurrentSystem];

  // stress fluctuations
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].xxxx+=ConfigurationalStressTensor[CurrentSystem].ax*ConfigurationalStressTensor[CurrentSystem].ax;
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].yxxx+=ConfigurationalStressTensor[CurrentSystem].bx*ConfigurationalStressTensor[CurrentSystem].ax;
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].zxxx+=ConfigurationalStressTensor[CurrentSystem].cx*ConfigurationalStressTensor[CurrentSystem].ax;
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].xxyx+=ConfigurationalStressTensor[CurrentSystem].ax*ConfigurationalStressTensor[CurrentSystem].bx;
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].yxyx+=ConfigurationalStressTensor[CurrentSystem].bx*ConfigurationalStressTensor[CurrentSystem].bx;
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].zxyx+=ConfigurationalStressTensor[CurrentSystem].cx*ConfigurationalStressTensor[CurrentSystem].bx;
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].xxzx+=ConfigurationalStressTensor[CurrentSystem].ax*ConfigurationalStressTensor[CurrentSystem].cx;
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].yxzx+=ConfigurationalStressTensor[CurrentSystem].bx*ConfigurationalStressTensor[CurrentSystem].cx;
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].zxzx+=ConfigurationalStressTensor[CurrentSystem].cx*ConfigurationalStressTensor[CurrentSystem].cx;

  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].xyxx+=ConfigurationalStressTensor[CurrentSystem].ay*ConfigurationalStressTensor[CurrentSystem].ax;
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].yyxx+=ConfigurationalStressTensor[CurrentSystem].by*ConfigurationalStressTensor[CurrentSystem].ax;
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].zyxx+=ConfigurationalStressTensor[CurrentSystem].cy*ConfigurationalStressTensor[CurrentSystem].ax;
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].xyyx+=ConfigurationalStressTensor[CurrentSystem].ay*ConfigurationalStressTensor[CurrentSystem].bx;
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].yyyx+=ConfigurationalStressTensor[CurrentSystem].by*ConfigurationalStressTensor[CurrentSystem].bx;
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].zyyx+=ConfigurationalStressTensor[CurrentSystem].cy*ConfigurationalStressTensor[CurrentSystem].bx;
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].xyzx+=ConfigurationalStressTensor[CurrentSystem].ay*ConfigurationalStressTensor[CurrentSystem].cx;
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].yyzx+=ConfigurationalStressTensor[CurrentSystem].by*ConfigurationalStressTensor[CurrentSystem].cx;
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].zyzx+=ConfigurationalStressTensor[CurrentSystem].cy*ConfigurationalStressTensor[CurrentSystem].cx;

  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].xzxx+=ConfigurationalStressTensor[CurrentSystem].az*ConfigurationalStressTensor[CurrentSystem].ax;
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].yzxx+=ConfigurationalStressTensor[CurrentSystem].bz*ConfigurationalStressTensor[CurrentSystem].ax;
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].zzxx+=ConfigurationalStressTensor[CurrentSystem].cz*ConfigurationalStressTensor[CurrentSystem].ax;
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].xzyx+=ConfigurationalStressTensor[CurrentSystem].az*ConfigurationalStressTensor[CurrentSystem].bx;
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].yzyx+=ConfigurationalStressTensor[CurrentSystem].bz*ConfigurationalStressTensor[CurrentSystem].bx;
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].zzyx+=ConfigurationalStressTensor[CurrentSystem].cz*ConfigurationalStressTensor[CurrentSystem].bx;
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].xzzx+=ConfigurationalStressTensor[CurrentSystem].az*ConfigurationalStressTensor[CurrentSystem].cx;
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].yzzx+=ConfigurationalStressTensor[CurrentSystem].bz*ConfigurationalStressTensor[CurrentSystem].cx;
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].zzzx+=ConfigurationalStressTensor[CurrentSystem].cz*ConfigurationalStressTensor[CurrentSystem].cx;


  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].xxxy+=ConfigurationalStressTensor[CurrentSystem].ax*ConfigurationalStressTensor[CurrentSystem].ay;
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].yxxy+=ConfigurationalStressTensor[CurrentSystem].bx*ConfigurationalStressTensor[CurrentSystem].ay;
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].zxxy+=ConfigurationalStressTensor[CurrentSystem].cx*ConfigurationalStressTensor[CurrentSystem].ay;
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].xxyy+=ConfigurationalStressTensor[CurrentSystem].ax*ConfigurationalStressTensor[CurrentSystem].by;
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].yxyy+=ConfigurationalStressTensor[CurrentSystem].bx*ConfigurationalStressTensor[CurrentSystem].by;
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].zxyy+=ConfigurationalStressTensor[CurrentSystem].cx*ConfigurationalStressTensor[CurrentSystem].by;
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].xxzy+=ConfigurationalStressTensor[CurrentSystem].ax*ConfigurationalStressTensor[CurrentSystem].cy;
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].yxzy+=ConfigurationalStressTensor[CurrentSystem].bx*ConfigurationalStressTensor[CurrentSystem].cy;
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].zxzy+=ConfigurationalStressTensor[CurrentSystem].cx*ConfigurationalStressTensor[CurrentSystem].cy;

  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].xyxy+=ConfigurationalStressTensor[CurrentSystem].ay*ConfigurationalStressTensor[CurrentSystem].ay;
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].yyxy+=ConfigurationalStressTensor[CurrentSystem].by*ConfigurationalStressTensor[CurrentSystem].ay;
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].zyxy+=ConfigurationalStressTensor[CurrentSystem].cy*ConfigurationalStressTensor[CurrentSystem].ay;
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].xyyy+=ConfigurationalStressTensor[CurrentSystem].ay*ConfigurationalStressTensor[CurrentSystem].by;
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].yyyy+=ConfigurationalStressTensor[CurrentSystem].by*ConfigurationalStressTensor[CurrentSystem].by;
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].zyyy+=ConfigurationalStressTensor[CurrentSystem].cy*ConfigurationalStressTensor[CurrentSystem].by;
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].xyzy+=ConfigurationalStressTensor[CurrentSystem].ay*ConfigurationalStressTensor[CurrentSystem].cy;
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].yyzy+=ConfigurationalStressTensor[CurrentSystem].by*ConfigurationalStressTensor[CurrentSystem].cy;
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].zyzy+=ConfigurationalStressTensor[CurrentSystem].cy*ConfigurationalStressTensor[CurrentSystem].cy;

  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].xzxy+=ConfigurationalStressTensor[CurrentSystem].az*ConfigurationalStressTensor[CurrentSystem].ay;
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].yzxy+=ConfigurationalStressTensor[CurrentSystem].bz*ConfigurationalStressTensor[CurrentSystem].ay;
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].zzxy+=ConfigurationalStressTensor[CurrentSystem].cz*ConfigurationalStressTensor[CurrentSystem].ay;
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].xzyy+=ConfigurationalStressTensor[CurrentSystem].az*ConfigurationalStressTensor[CurrentSystem].by;
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].yzyy+=ConfigurationalStressTensor[CurrentSystem].bz*ConfigurationalStressTensor[CurrentSystem].by;
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].zzyy+=ConfigurationalStressTensor[CurrentSystem].cz*ConfigurationalStressTensor[CurrentSystem].by;
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].xzzy+=ConfigurationalStressTensor[CurrentSystem].az*ConfigurationalStressTensor[CurrentSystem].cy;
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].yzzy+=ConfigurationalStressTensor[CurrentSystem].bz*ConfigurationalStressTensor[CurrentSystem].cy;
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].zzzy+=ConfigurationalStressTensor[CurrentSystem].cz*ConfigurationalStressTensor[CurrentSystem].cy;


  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].xxxz+=ConfigurationalStressTensor[CurrentSystem].ax*ConfigurationalStressTensor[CurrentSystem].az;
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].yxxz+=ConfigurationalStressTensor[CurrentSystem].bx*ConfigurationalStressTensor[CurrentSystem].az;
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].zxxz+=ConfigurationalStressTensor[CurrentSystem].cx*ConfigurationalStressTensor[CurrentSystem].az;
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].xxyz+=ConfigurationalStressTensor[CurrentSystem].ax*ConfigurationalStressTensor[CurrentSystem].bz;
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].yxyz+=ConfigurationalStressTensor[CurrentSystem].bx*ConfigurationalStressTensor[CurrentSystem].bz;
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].zxyz+=ConfigurationalStressTensor[CurrentSystem].cx*ConfigurationalStressTensor[CurrentSystem].bz;
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].xxzz+=ConfigurationalStressTensor[CurrentSystem].ax*ConfigurationalStressTensor[CurrentSystem].cz;
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].yxzz+=ConfigurationalStressTensor[CurrentSystem].bx*ConfigurationalStressTensor[CurrentSystem].cz;
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].zxzz+=ConfigurationalStressTensor[CurrentSystem].cx*ConfigurationalStressTensor[CurrentSystem].cz;

  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].xyxz+=ConfigurationalStressTensor[CurrentSystem].ay*ConfigurationalStressTensor[CurrentSystem].az;
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].yyxz+=ConfigurationalStressTensor[CurrentSystem].by*ConfigurationalStressTensor[CurrentSystem].az;
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].zyxz+=ConfigurationalStressTensor[CurrentSystem].cy*ConfigurationalStressTensor[CurrentSystem].az;
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].xyyz+=ConfigurationalStressTensor[CurrentSystem].ay*ConfigurationalStressTensor[CurrentSystem].bz;
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].yyyz+=ConfigurationalStressTensor[CurrentSystem].by*ConfigurationalStressTensor[CurrentSystem].bz;
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].zyyz+=ConfigurationalStressTensor[CurrentSystem].cy*ConfigurationalStressTensor[CurrentSystem].bz;
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].xyzz+=ConfigurationalStressTensor[CurrentSystem].ay*ConfigurationalStressTensor[CurrentSystem].cz;
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].yyzz+=ConfigurationalStressTensor[CurrentSystem].by*ConfigurationalStressTensor[CurrentSystem].cz;
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].zyzz+=ConfigurationalStressTensor[CurrentSystem].cy*ConfigurationalStressTensor[CurrentSystem].cz;

  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].xzxz+=ConfigurationalStressTensor[CurrentSystem].az*ConfigurationalStressTensor[CurrentSystem].az;
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].yzxz+=ConfigurationalStressTensor[CurrentSystem].bz*ConfigurationalStressTensor[CurrentSystem].az;
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].zzxz+=ConfigurationalStressTensor[CurrentSystem].cz*ConfigurationalStressTensor[CurrentSystem].az;
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].xzyz+=ConfigurationalStressTensor[CurrentSystem].az*ConfigurationalStressTensor[CurrentSystem].bz;
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].yzyz+=ConfigurationalStressTensor[CurrentSystem].bz*ConfigurationalStressTensor[CurrentSystem].bz;
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].zzyz+=ConfigurationalStressTensor[CurrentSystem].cz*ConfigurationalStressTensor[CurrentSystem].bz;
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].xzzz+=ConfigurationalStressTensor[CurrentSystem].az*ConfigurationalStressTensor[CurrentSystem].cz;
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].yzzz+=ConfigurationalStressTensor[CurrentSystem].bz*ConfigurationalStressTensor[CurrentSystem].cz;
  ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][Block].zzzz+=ConfigurationalStressTensor[CurrentSystem].cz*ConfigurationalStressTensor[CurrentSystem].cz;
}

REAL_MATRIX3x3 ComputeCurrentAverageStressTensor(void)
{
  int i;
  REAL count;
  REAL_MATRIX3x3 Average;

  count=0.0;
  Average.ax=0.0; Average.bx=0.0; Average.cx=0.0;
  Average.ay=0.0; Average.by=0.0; Average.cy=0.0;
  Average.az=0.0; Average.bz=0.0; Average.cz=0.0;
  for(i=0;i<NR_BLOCKS;i++)
  {
    if(BlockWeightedCount[CurrentSystem][i]>0.0)
    {
      Average.ax+=StressTensorAccumulated[CurrentSystem][i].ax;
      Average.ay+=StressTensorAccumulated[CurrentSystem][i].ay;
      Average.az+=StressTensorAccumulated[CurrentSystem][i].az;

      Average.bx+=StressTensorAccumulated[CurrentSystem][i].bx;
      Average.by+=StressTensorAccumulated[CurrentSystem][i].by;
      Average.bz+=StressTensorAccumulated[CurrentSystem][i].bz;

      Average.cx+=StressTensorAccumulated[CurrentSystem][i].cx;
      Average.cy+=StressTensorAccumulated[CurrentSystem][i].cy;
      Average.cz+=StressTensorAccumulated[CurrentSystem][i].cz;

      count+=BlockWeightedCount[CurrentSystem][i];
    }
  }
  if(count>0)
  {
    Average.ax/=count; Average.bx/=count; Average.cx/=count;
    Average.ay/=count; Average.by/=count; Average.cy/=count;
    Average.az/=count; Average.bz/=count; Average.cz/=count;
  }
  return Average;
}

REAL_MATRIX3x3 ComputeCurrentAverageConfigurationalStressTensor(void)
{
  int i;
  REAL count;
  REAL_MATRIX3x3 Average;

  count=0.0;
  Average.ax=0.0; Average.bx=0.0; Average.cx=0.0;
  Average.ay=0.0; Average.by=0.0; Average.cy=0.0;
  Average.az=0.0; Average.bz=0.0; Average.cz=0.0;
  for(i=0;i<NR_BLOCKS;i++)
  {
    if(BlockWeightedCount[CurrentSystem][i]>0.0)
    {
      Average.ax+=ConfigurationalStressTensorAccumulated[CurrentSystem][i].ax;
      Average.ay+=ConfigurationalStressTensorAccumulated[CurrentSystem][i].ay;
      Average.az+=ConfigurationalStressTensorAccumulated[CurrentSystem][i].az;

      Average.bx+=ConfigurationalStressTensorAccumulated[CurrentSystem][i].bx;
      Average.by+=ConfigurationalStressTensorAccumulated[CurrentSystem][i].by;
      Average.bz+=ConfigurationalStressTensorAccumulated[CurrentSystem][i].bz;

      Average.cx+=ConfigurationalStressTensorAccumulated[CurrentSystem][i].cx;
      Average.cy+=ConfigurationalStressTensorAccumulated[CurrentSystem][i].cy;
      Average.cz+=ConfigurationalStressTensorAccumulated[CurrentSystem][i].cz;

      count+=BlockWeightedCount[CurrentSystem][i];
    }
  }
  if(count>0)
  {
    Average.ax/=count; Average.bx/=count; Average.cx/=count;
    Average.ay/=count; Average.by/=count; Average.cy/=count;
    Average.az/=count; Average.bz/=count; Average.cz/=count;
  }
  return Average;
}


REAL GetAverageTemperature(void)
{
  int i;
  REAL sum1,sum2;

  sum1=sum2=0.0;
  for(i=0;i<NR_BLOCKS;i++)
    if(BlockWeightedCount[CurrentSystem][i]>0.0)
    {
      sum1+=TemperatureAccumulated[CurrentSystem][i];
      sum2+=BlockWeightedCount[CurrentSystem][i];
    }
  if((sum1>0.0)&&(sum2>0.0))
    return sum1/sum2;
  else
    return 0.0;
}

REAL GetAverageCellTemperature(void)
{
  int i;
  REAL sum1,sum2;

  sum1=sum2=0.0;
  for(i=0;i<NR_BLOCKS;i++)
    if(BlockWeightedCount[CurrentSystem][i]>0.0)
    {
      sum1+=TemperatureCellAccumulated[CurrentSystem][i];
      sum2+=BlockWeightedCount[CurrentSystem][i];
    }
  if((sum1>0.0)&&(sum2>0.0))
    return sum1/sum2;
  else
    return 0.0;
}


REAL GetAverageVolume(void)
{
  int i;
  REAL count,sum;

  count=sum=0.0;
  for(i=0;i<NR_BLOCKS;i++)
  {
    if(BlockWeightedCount[CurrentSystem][i]>0.0)
    {
      sum+=VolumeAccumulated[CurrentSystem][i];
      count+=BlockWeightedCount[CurrentSystem][i];
    }
  }
  return sum/count;
}


REAL_MATRIX6x6 ComputeCurrentAverageElasticConstants(REAL_MATRIX6x6 *ElasticBornTerm, REAL_MATRIX6x6 *ElasticStressFluctuationTerm,REAL_MATRIX6x6 *ElasticIdealGasTerm)
{
  int i;
  REAL count;
  REAL_MATRIX3x3 Stress;
  REAL_MATRIX9x9 BornTerm;
  REAL_MATRIX9x9 ElasticConstants;
  REAL_MATRIX9x9 StressFluctuationTerm;
  REAL_MATRIX9x9 IdealGasTerm;
  REAL_MATRIX6x6 ElasticConstantsVoight;
  REAL V,T,rho;

  count=0.0;
  InitializeMatrix6x6(&ElasticConstantsVoight);
  InitializeMatrix9x9(&BornTerm);
  InitializeMatrix9x9(&StressFluctuationTerm);
  InitializeMatrix9x9(&IdealGasTerm);
  InitializeMatrix9x9(&ElasticConstants);
  Stress=ComputeCurrentAverageConfigurationalStressTensor();
  V=GetAverageVolume();
  T=GetAverageTemperature();
  rho=Framework[CurrentSystem].TotalNumberOfAtoms/V;

  Stress.ax=Stress.bx=Stress.cx=0.0;
  Stress.ay=Stress.by=Stress.cy=0.0;
  Stress.az=Stress.bz=Stress.cz=0.0;

  for(i=0;i<NR_BLOCKS;i++)
    if(BlockWeightedCount[CurrentSystem][i]>0.0)
    {
      AddRealMatrix9x9(&BornTerm,BornTerm,BornTermAccumulated[CurrentSystem][i]);
      AddRealMatrix9x9(&StressFluctuationTerm,StressFluctuationTerm,ConfigurationalStressFluctuationTermAccumulated[CurrentSystem][i]);

      Stress.ax+=ConfigurationalStressTensorAccumulated[CurrentSystem][i].ax;
      Stress.ay+=ConfigurationalStressTensorAccumulated[CurrentSystem][i].ay;
      Stress.az+=ConfigurationalStressTensorAccumulated[CurrentSystem][i].az;

      Stress.bx+=ConfigurationalStressTensorAccumulated[CurrentSystem][i].bx;
      Stress.by+=ConfigurationalStressTensorAccumulated[CurrentSystem][i].by;
      Stress.bz+=ConfigurationalStressTensorAccumulated[CurrentSystem][i].bz;

      Stress.cx+=ConfigurationalStressTensorAccumulated[CurrentSystem][i].cx;
      Stress.cy+=ConfigurationalStressTensorAccumulated[CurrentSystem][i].cy;
      Stress.cz+=ConfigurationalStressTensorAccumulated[CurrentSystem][i].cz;

      count+=BlockWeightedCount[CurrentSystem][i];
    }

  DivideRealMatrix9x9ByReal(&BornTerm,BornTerm,count);
  DivideRealMatrix9x9ByReal(&StressFluctuationTerm,StressFluctuationTerm,count);

  Stress.ax/=count; Stress.bx/=count; Stress.cx/=count;
  Stress.ay/=count; Stress.by/=count; Stress.cy/=count;
  Stress.az/=count; Stress.bz/=count; Stress.cz/=count;

  StressFluctuationTerm.xxxx-=Stress.ax*Stress.ax;
  StressFluctuationTerm.yxxx-=Stress.bx*Stress.ax;
  StressFluctuationTerm.zxxx-=Stress.cx*Stress.ax;
  StressFluctuationTerm.xxyx-=Stress.ax*Stress.bx;
  StressFluctuationTerm.yxyx-=Stress.bx*Stress.bx;
  StressFluctuationTerm.zxyx-=Stress.cx*Stress.bx;
  StressFluctuationTerm.xxzx-=Stress.ax*Stress.cx;
  StressFluctuationTerm.yxzx-=Stress.bx*Stress.cx;
  StressFluctuationTerm.zxzx-=Stress.cx*Stress.cx;

  StressFluctuationTerm.xyxx-=Stress.ay*Stress.ax;
  StressFluctuationTerm.yyxx-=Stress.by*Stress.ax;
  StressFluctuationTerm.zyxx-=Stress.cy*Stress.ax;
  StressFluctuationTerm.xyyx-=Stress.ay*Stress.bx;
  StressFluctuationTerm.yyyx-=Stress.by*Stress.bx;
  StressFluctuationTerm.zyyx-=Stress.cy*Stress.bx;
  StressFluctuationTerm.xyzx-=Stress.ay*Stress.cx;
  StressFluctuationTerm.yyzx-=Stress.by*Stress.cx;
  StressFluctuationTerm.zyzx-=Stress.cy*Stress.cx;

  StressFluctuationTerm.xzxx-=Stress.az*Stress.ax;
  StressFluctuationTerm.yzxx-=Stress.bz*Stress.ax;
  StressFluctuationTerm.zzxx-=Stress.cz*Stress.ax;
  StressFluctuationTerm.xzyx-=Stress.az*Stress.bx;
  StressFluctuationTerm.yzyx-=Stress.bz*Stress.bx;
  StressFluctuationTerm.zzyx-=Stress.cz*Stress.bx;
  StressFluctuationTerm.xzzx-=Stress.az*Stress.cx;
  StressFluctuationTerm.yzzx-=Stress.bz*Stress.cx;
  StressFluctuationTerm.zzzx-=Stress.cz*Stress.cx;


  StressFluctuationTerm.xxxy-=Stress.ax*Stress.ay;
  StressFluctuationTerm.yxxy-=Stress.bx*Stress.ay;
  StressFluctuationTerm.zxxy-=Stress.cx*Stress.ay;
  StressFluctuationTerm.xxyy-=Stress.ax*Stress.by;
  StressFluctuationTerm.yxyy-=Stress.bx*Stress.by;
  StressFluctuationTerm.zxyy-=Stress.cx*Stress.by;
  StressFluctuationTerm.xxzy-=Stress.ax*Stress.cy;
  StressFluctuationTerm.yxzy-=Stress.bx*Stress.cy;
  StressFluctuationTerm.zxzy-=Stress.cx*Stress.cy;

  StressFluctuationTerm.xyxy-=Stress.ay*Stress.ay;
  StressFluctuationTerm.yyxy-=Stress.by*Stress.ay;
  StressFluctuationTerm.zyxy-=Stress.cy*Stress.ay;
  StressFluctuationTerm.xyyy-=Stress.ay*Stress.by;
  StressFluctuationTerm.yyyy-=Stress.by*Stress.by;
  StressFluctuationTerm.zyyy-=Stress.cy*Stress.by;
  StressFluctuationTerm.xyzy-=Stress.ay*Stress.cy;
  StressFluctuationTerm.yyzy-=Stress.by*Stress.cy;
  StressFluctuationTerm.zyzy-=Stress.cy*Stress.cy;

  StressFluctuationTerm.xzxy-=Stress.az*Stress.ay;
  StressFluctuationTerm.yzxy-=Stress.bz*Stress.ay;
  StressFluctuationTerm.zzxy-=Stress.cz*Stress.ay;
  StressFluctuationTerm.xzyy-=Stress.az*Stress.by;
  StressFluctuationTerm.yzyy-=Stress.bz*Stress.by;
  StressFluctuationTerm.zzyy-=Stress.cz*Stress.by;
  StressFluctuationTerm.xzzy-=Stress.az*Stress.cy;
  StressFluctuationTerm.yzzy-=Stress.bz*Stress.cy;
  StressFluctuationTerm.zzzy-=Stress.cz*Stress.cy;


  StressFluctuationTerm.xxxz-=Stress.ax*Stress.az;
  StressFluctuationTerm.yxxz-=Stress.bx*Stress.az;
  StressFluctuationTerm.zxxz-=Stress.cx*Stress.az;
  StressFluctuationTerm.xxyz-=Stress.ax*Stress.bz;
  StressFluctuationTerm.yxyz-=Stress.bx*Stress.bz;
  StressFluctuationTerm.zxyz-=Stress.cx*Stress.bz;
  StressFluctuationTerm.xxzz-=Stress.ax*Stress.cz;
  StressFluctuationTerm.yxzz-=Stress.bx*Stress.cz;
  StressFluctuationTerm.zxzz-=Stress.cx*Stress.cz;

  StressFluctuationTerm.xyxz-=Stress.ay*Stress.az;
  StressFluctuationTerm.yyxz-=Stress.by*Stress.az;
  StressFluctuationTerm.zyxz-=Stress.cy*Stress.az;
  StressFluctuationTerm.xyyz-=Stress.ay*Stress.bz;
  StressFluctuationTerm.yyyz-=Stress.by*Stress.bz;
  StressFluctuationTerm.zyyz-=Stress.cy*Stress.bz;
  StressFluctuationTerm.xyzz-=Stress.ay*Stress.cz;
  StressFluctuationTerm.yyzz-=Stress.by*Stress.cz;
  StressFluctuationTerm.zyzz-=Stress.cy*Stress.cz;

  StressFluctuationTerm.xzxz-=Stress.az*Stress.az;
  StressFluctuationTerm.yzxz-=Stress.bz*Stress.az;
  StressFluctuationTerm.zzxz-=Stress.cz*Stress.az;
  StressFluctuationTerm.xzyz-=Stress.az*Stress.bz;
  StressFluctuationTerm.yzyz-=Stress.bz*Stress.bz;
  StressFluctuationTerm.zzyz-=Stress.cz*Stress.bz;
  StressFluctuationTerm.xzzz-=Stress.az*Stress.cz;
  StressFluctuationTerm.yzzz-=Stress.bz*Stress.cz;
  StressFluctuationTerm.zzzz-=Stress.cz*Stress.cz;

  DivideRealMatrix9x9ByReal(&StressFluctuationTerm,StressFluctuationTerm,(K_B*T)/V);

  IdealGasTerm.xxxx=2.0*rho*K_B*T;
  IdealGasTerm.yxyx=rho*K_B*T;
  IdealGasTerm.zxzx=rho*K_B*T;
  IdealGasTerm.xyyx=rho*K_B*T;
  IdealGasTerm.xzzx=rho*K_B*T;
  IdealGasTerm.yxxy=rho*K_B*T;
  IdealGasTerm.xyxy=rho*K_B*T;
  IdealGasTerm.yyyy=2.0*rho*K_B*T;
  IdealGasTerm.zyzy=rho*K_B*T;
  IdealGasTerm.yzzy=rho*K_B*T;
  IdealGasTerm.zxxz=rho*K_B*T;
  IdealGasTerm.zyyz=rho*K_B*T;
  IdealGasTerm.xzxz=rho*K_B*T;
  IdealGasTerm.yzyz=rho*K_B*T;
  IdealGasTerm.zzzz=2.0*rho*K_B*T;

  // convert to Voight notation
  // 1=xx 2=yy 3=zz 4=yz 5=zx 6=xy

  *ElasticBornTerm=ConvertToVoigt3D(BornTerm);
  *ElasticStressFluctuationTerm=ConvertToVoigt3D(StressFluctuationTerm);
  *ElasticIdealGasTerm=ConvertToVoigt3D(IdealGasTerm);

  AddRealMatrix6x6(&ElasticConstantsVoight,ElasticConstantsVoight,*ElasticBornTerm);
  SubtractRealMatrix6x6(&ElasticConstantsVoight,ElasticConstantsVoight,*ElasticStressFluctuationTerm);
  AddRealMatrix6x6(&ElasticConstantsVoight,ElasticConstantsVoight,*ElasticIdealGasTerm);

  return ElasticConstantsVoight;
}

void InitializesEnergiesAllSystems(void)
{
  for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
    InitializesEnergiesCurrentSystem();
  CurrentSystem=0;
}

void InitializesEnergiesCurrentSystem(void)
{
  // initialize energies
  UHostBond[CurrentSystem]=0.0;
  UHostUreyBradley[CurrentSystem]=0.0;
  UHostBend[CurrentSystem]=0.0;
  UHostInversionBend[CurrentSystem]=0.0;
  UHostTorsion[CurrentSystem]=0.0;
  UHostImproperTorsion[CurrentSystem]=0.0;
  UHostBondBond[CurrentSystem]=0.0;
  UHostBendBend[CurrentSystem]=0.0;
  UHostBondBend[CurrentSystem]=0.0;
  UHostBondTorsion[CurrentSystem]=0.0;
  UHostBendTorsion[CurrentSystem]=0.0;

  UCationBond[CurrentSystem]=0.0;
  UCationUreyBradley[CurrentSystem]=0.0;
  UCationBend[CurrentSystem]=0.0;
  UCationInversionBend[CurrentSystem]=0.0;
  UCationTorsion[CurrentSystem]=0.0;
  UCationImproperTorsion[CurrentSystem]=0.0;
  UCationBondBond[CurrentSystem]=0.0;
  UCationBendBend[CurrentSystem]=0.0;
  UCationBondBend[CurrentSystem]=0.0;
  UCationBondTorsion[CurrentSystem]=0.0;
  UCationBendTorsion[CurrentSystem]=0.0;
  UCationIntraVDW[CurrentSystem]=0.0;
  UCationIntraChargeCharge[CurrentSystem]=0.0;
  UCationIntraChargeBondDipole[CurrentSystem]=0.0;
  UCationIntraBondDipoleBondDipole[CurrentSystem]=0.0;

  UAdsorbateBond[CurrentSystem]=0.0;
  UAdsorbateUreyBradley[CurrentSystem]=0.0;
  UAdsorbateBend[CurrentSystem]=0.0;
  UAdsorbateInversionBend[CurrentSystem]=0.0;
  UAdsorbateTorsion[CurrentSystem]=0.0;
  UAdsorbateImproperTorsion[CurrentSystem]=0.0;
  UAdsorbateBondBond[CurrentSystem]=0.0;
  UAdsorbateBendBend[CurrentSystem]=0.0;
  UAdsorbateBondBend[CurrentSystem]=0.0;
  UAdsorbateBondTorsion[CurrentSystem]=0.0;
  UAdsorbateBendTorsion[CurrentSystem]=0.0;
  UAdsorbateIntraVDW[CurrentSystem]=0.0;
  UAdsorbateIntraChargeCharge[CurrentSystem]=0.0;
  UAdsorbateIntraChargeBondDipole[CurrentSystem]=0.0;
  UAdsorbateIntraBondDipoleBondDipole[CurrentSystem]=0.0;

  UHostHost[CurrentSystem]=0.0;
  UAdsorbateAdsorbate[CurrentSystem]=0.0;
  UCationCation[CurrentSystem]=0.0;
  UHostAdsorbate[CurrentSystem]=0.0;
  UHostCation[CurrentSystem]=0.0;
  UAdsorbateCation[CurrentSystem]=0.0;

  UHostHostVDW[CurrentSystem]=0.0;
  UAdsorbateAdsorbateVDW[CurrentSystem]=0.0;
  UCationCationVDW[CurrentSystem]=0.0;
  UHostAdsorbateVDW[CurrentSystem]=0.0;
  UHostCationVDW[CurrentSystem]=0.0;
  UAdsorbateCationVDW[CurrentSystem]=0.0;

  UHostHostCoulomb[CurrentSystem]=0.0;
  UAdsorbateAdsorbateCoulomb[CurrentSystem]=0.0;
  UCationCationCoulomb[CurrentSystem]=0.0;
  UHostAdsorbateCoulomb[CurrentSystem]=0.0;
  UHostCationCoulomb[CurrentSystem]=0.0;
  UAdsorbateCationCoulomb[CurrentSystem]=0.0;

  UHostHostChargeChargeReal[CurrentSystem]=0.0;
  UAdsorbateAdsorbateChargeChargeReal[CurrentSystem]=0.0;
  UCationCationChargeChargeReal[CurrentSystem]=0.0;
  UHostAdsorbateChargeChargeReal[CurrentSystem]=0.0;
  UHostCationChargeChargeReal[CurrentSystem]=0.0;
  UAdsorbateCationChargeChargeReal[CurrentSystem]=0.0;

  UHostHostChargeBondDipoleReal[CurrentSystem]=0.0;
  UAdsorbateAdsorbateChargeBondDipoleReal[CurrentSystem]=0.0;
  UCationCationChargeBondDipoleReal[CurrentSystem]=0.0;
  UHostAdsorbateChargeBondDipoleReal[CurrentSystem]=0.0;
  UHostCationChargeBondDipoleReal[CurrentSystem]=0.0;
  UAdsorbateCationChargeBondDipoleReal[CurrentSystem]=0.0;

  UHostHostBondDipoleBondDipoleReal[CurrentSystem]=0.0;
  UAdsorbateAdsorbateBondDipoleBondDipoleReal[CurrentSystem]=0.0;
  UCationCationBondDipoleBondDipoleReal[CurrentSystem]=0.0;
  UHostAdsorbateBondDipoleBondDipoleReal[CurrentSystem]=0.0;
  UHostCationBondDipoleBondDipoleReal[CurrentSystem]=0.0;
  UAdsorbateCationBondDipoleBondDipoleReal[CurrentSystem]=0.0;

  UHostHostChargeChargeFourier[CurrentSystem]=0.0;
  UAdsorbateAdsorbateChargeChargeFourier[CurrentSystem]=0.0;
  UCationCationChargeChargeFourier[CurrentSystem]=0.0;
  UHostAdsorbateChargeChargeFourier[CurrentSystem]=0.0;
  UHostCationChargeChargeFourier[CurrentSystem]=0.0;
  UAdsorbateCationChargeChargeFourier[CurrentSystem]=0.0;

  UHostHostChargeBondDipoleFourier[CurrentSystem]=0.0;
  UAdsorbateAdsorbateChargeBondDipoleFourier[CurrentSystem]=0.0;
  UCationCationChargeBondDipoleFourier[CurrentSystem]=0.0;
  UHostAdsorbateChargeBondDipoleFourier[CurrentSystem]=0.0;
  UHostCationChargeBondDipoleFourier[CurrentSystem]=0.0;
  UAdsorbateCationChargeBondDipoleFourier[CurrentSystem]=0.0;

  UHostHostBondDipoleBondDipoleFourier[CurrentSystem]=0.0;
  UAdsorbateAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]=0.0;
  UCationCationBondDipoleBondDipoleFourier[CurrentSystem]=0.0;
  UHostAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]=0.0;
  UHostCationBondDipoleBondDipoleFourier[CurrentSystem]=0.0;
  UAdsorbateCationBondDipoleBondDipoleFourier[CurrentSystem]=0.0;

  UHostPolarization[CurrentSystem]=0.0;
  UAdsorbatePolarization[CurrentSystem]=0.0;
  UCationPolarization[CurrentSystem]=0.0;
  UHostBackPolarization[CurrentSystem]=0.0;
  UAdsorbateBackPolarization[CurrentSystem]=0.0;
  UCationBackPolarization[CurrentSystem]=0.0;

  UTailCorrection[CurrentSystem]=0.0;
  UDistanceConstraints[CurrentSystem]=0.0;
  UAngleConstraints[CurrentSystem]=0.0;
  UDihedralConstraints[CurrentSystem]=0.0;
  UInversionBendConstraints[CurrentSystem]=0.0;
  UOutOfPlaneDistanceConstraints[CurrentSystem]=0.0;
  UExclusionConstraints[CurrentSystem]=0.0;
  UNoseHoover[CurrentSystem]=0.0;
  UTotal[CurrentSystem]=0.0;

  NumberOfRattleCyclesStage1[CurrentSystem]=0.0;
  MaximumNumberOfRattleCyclesStage1[CurrentSystem]=0;
  NumberOfRattleCyclesStage2[CurrentSystem]=0.0;
  MaximumNumberOfRattleCyclesStage2[CurrentSystem]=0;
}

void InitializesEnergyAveragesAllSystems(void)
{
  int i,j,k,l;

  Block=0;

  for(i=0;i<NR_BLOCKS;i++)
  {
    BlockCycle[i]=(long long)((i+1)*NumberOfCycles/(double)NR_BLOCKS);
    for(j=0;j<NumberOfSystems;j++)
    {
      BlockWeightedCount[j][i]=0.0;
      BlockCount[j][i]=0.0;
    }
  }

  for(k=0;k<NumberOfSystems;k++)
  {
    for(i=0;i<NR_BLOCKS;i++)
    {
      ConfigurationalStressTensorAccumulated[k][i].ax=0.0; ConfigurationalStressTensorAccumulated[k][i].bx=0.0; ConfigurationalStressTensorAccumulated[k][i].cx=0.0;
      ConfigurationalStressTensorAccumulated[k][i].ay=0.0; ConfigurationalStressTensorAccumulated[k][i].by=0.0; ConfigurationalStressTensorAccumulated[k][i].cy=0.0;
      ConfigurationalStressTensorAccumulated[k][i].az=0.0; ConfigurationalStressTensorAccumulated[k][i].bz=0.0; ConfigurationalStressTensorAccumulated[k][i].cz=0.0;

      // initialize energies
      UHostBondAccumulated[k][i]=0.0;
      UHostUreyBradleyAccumulated[k][i]=0.0;
      UHostBendAccumulated[k][i]=0.0;
      UHostInversionBendAccumulated[k][i]=0.0;
      UHostTorsionAccumulated[k][i]=0.0;
      UHostImproperTorsionAccumulated[k][i]=0.0;
      UHostBondBondAccumulated[k][i]=0.0;
      UHostBendBendAccumulated[k][i]=0.0;
      UHostBondBendAccumulated[k][i]=0.0;
      UHostBondTorsionAccumulated[k][i]=0.0;
      UHostBendTorsionAccumulated[k][i]=0.0;

      UCationBondAccumulated[k][i]=0.0;
      UCationUreyBradleyAccumulated[k][i]=0.0;
      UCationBendAccumulated[k][i]=0.0;
      UCationInversionBendAccumulated[k][i]=0.0;
      UCationTorsionAccumulated[k][i]=0.0;
      UCationImproperTorsionAccumulated[k][i]=0.0;
      UCationBondBondAccumulated[k][i]=0.0;
      UCationBendBendAccumulated[k][i]=0.0;
      UCationBondBendAccumulated[k][i]=0.0;
      UCationBondTorsionAccumulated[k][i]=0.0;
      UCationBendTorsionAccumulated[k][i]=0.0;
      UCationIntraVDWAccumulated[k][i]=0.0;
      UCationIntraChargeChargeAccumulated[k][i]=0.0;
      UCationIntraChargeBondDipoleAccumulated[k][i]=0.0;
      UCationIntraBondDipoleBondDipoleAccumulated[k][i]=0.0;

      UAdsorbateBondAccumulated[k][i]=0.0;
      UAdsorbateUreyBradleyAccumulated[k][i]=0.0;
      UAdsorbateBondBondAccumulated[k][i]=0.0;
      UAdsorbateBendAccumulated[k][i]=0.0;
      UAdsorbateInversionBendAccumulated[k][i]=0.0;
      UAdsorbateTorsionAccumulated[k][i]=0.0;
      UAdsorbateImproperTorsionAccumulated[k][i]=0.0;
      UAdsorbateBondBondAccumulated[k][i]=0.0;
      UAdsorbateBendBendAccumulated[k][i]=0.0;
      UAdsorbateBondBendAccumulated[k][i]=0.0;
      UAdsorbateBondTorsionAccumulated[k][i]=0.0;
      UAdsorbateBendTorsionAccumulated[k][i]=0.0;
      UAdsorbateIntraVDWAccumulated[k][i]=0.0;
      UAdsorbateIntraChargeChargeAccumulated[k][i]=0.0;
      UAdsorbateIntraChargeBondDipoleAccumulated[k][i]=0.0;
      UAdsorbateIntraBondDipoleBondDipoleAccumulated[k][i]=0.0;

      UHostHostAccumulated[k][i]=0.0;
      UAdsorbateAdsorbateAccumulated[k][i]=0.0;
      UCationCationAccumulated[k][i]=0.0;
      UHostAdsorbateAccumulated[k][i]=0.0;
      UHostCationAccumulated[k][i]=0.0;
      UAdsorbateCationAccumulated[k][i]=0.0;

      UHostHostVDWAccumulated[k][i]=0.0;
      UAdsorbateAdsorbateVDWAccumulated[k][i]=0.0;
      UCationCationVDWAccumulated[k][i]=0.0;
      UHostAdsorbateVDWAccumulated[k][i]=0.0;
      UHostCationVDWAccumulated[k][i]=0.0;
      UAdsorbateCationVDWAccumulated[k][i]=0.0;

      UHostHostCoulombAccumulated[k][i]=0.0;
      UAdsorbateAdsorbateCoulombAccumulated[k][i]=0.0;
      UCationCationCoulombAccumulated[k][i]=0.0;
      UHostAdsorbateCoulombAccumulated[k][i]=0.0;
      UHostCationCoulombAccumulated[k][i]=0.0;
      UAdsorbateCationCoulombAccumulated[k][i]=0.0;

      TotalSystemDipoleAccumulated[k][i].x=0.0;
      TotalSystemDipoleAccumulated[k][i].y=0.0;
      TotalSystemDipoleAccumulated[k][i].z=0.0;
      TotalSystemDipoleSquaredAccumulated[k][i].x=0.0;
      TotalSystemDipoleSquaredAccumulated[k][i].y=0.0;
      TotalSystemDipoleSquaredAccumulated[k][i].z=0.0;
      TotalSystemNormDipoleAccumulated[k][i]=0.0;
      TotalSystemNormDipoleSquaredAccumulated[k][i]=0.0;

      UTailCorrectionAccumulated[k][i]=0.0;
      UDistanceConstraintsAccumulated[k][i]=0.0;
      UAngleConstraintsAccumulated[k][i]=0.0;
      UDihedralConstraintsAccumulated[k][i]=0.0;
      UInversionBendConstraintsAccumulated[k][i]=0.0;
      UOutOfPlaneDistanceConstraintsAccumulated[k][i]=0.0;
      UExclusionConstraintsAccumulated[k][i]=0.0;

      UHostPolarizationAccumulated[k][i]=0.0;
      UAdsorbatePolarizationAccumulated[k][i]=0.0;
      UCationPolarizationAccumulated[k][i]=0.0;
      UHostBackPolarizationAccumulated[k][i]=0.0;
      UAdsorbateBackPolarizationAccumulated[k][i]=0.0;
      UCationBackPolarizationAccumulated[k][i]=0.0;

      UTotalAccumulated[k][i]=0.0;
      NumberOfIntegerMoleculesAccumulated[k][i]=0.0;
      DensityAccumulated[k][i]=0.0;
      InverseDensityAccumulated[k][i]=0.0;
      GibbsInverseDensityAccumulated[k][i]=0.0;
      for(j=0;j<NumberOfComponents;j++)
      {
        GibbsInverseDensityPerComponentAccumulated[k][j][i]=0.0;
        NumberOfIntegerMoleculesPerComponentAccumulated[k][j][i]=0.0;
        NumberOfFractionalMoleculesPerComponentAccumulated[k][j][i]=0.0;
        NumberOfExcessMoleculesPerComponentAccumulated[k][j][i]=0.0;
        DensityPerComponentAccumulated[k][j][i]=0.0;

        WidomRosenbluthFactorAccumulated[k][j][i]=0.0;
        WidomIdealGasAccumulated[k][j][i]=0.0;
        WidomRosenbluthFactorCount[k][j][i]=0.0;
        WidomEnergyDifferenceAccumulated[k][j][i]=0.0;
        WidomEnergyFrameworkAccumulated[k][j][i]=0.0;
        WidomEnergyFrameworkCount[k][j][i]=0.0;

        GibbsWidomRosenbluthFactorAccumulated[k][j][i]=0.0;
        GibbsWidomIdealGasAccumulated[k][j][i]=0.0;
        GibbsWidomRosenbluthFactorCount[k][j][i]=0.0;
        GibbsWidomEnergyDifferenceAccumulated[k][j][i]=0.0;
        GibbsWidomEnergyFrameworkAccumulated[k][j][i]=0.0;
        GibbsWidomEnergyFrameworkCount[k][j][i]=0.0;

        TotalEnergyTimesNumberOfMoleculesPerComponentAccumulated[k][j][i]=0.0;
        HostAdsorbateEnergyTimesNumberOfMoleculesAccumulated[k][j][i]=0.0;
        AdsorbateAdsorbateEnergyTimesNumberOfMoleculesAccumulated[k][j][i]=0.0;

        for(l=0;l<NumberOfComponents;l++)
          NumberOfMoleculesPerComponentSquaredAccumulated[k][j][l][i]=0.0;
      }
      TotalEnergyTimesNumberOfMoleculesAccumulated[k][i]=0.0;
      NumberOfMoleculesSquaredAccumulated[k][i]=0.0;

      TemperatureAccumulated[k][i]=0.0;
      TemperatureCellAccumulated[k][i]=0.0;
      MolecularPressureAccumulated[k][i]=0.0;
      CompressibilityAccumulated[k][i]=0.0;

      MolecularStressTensorAccumulated[k][i].ax=0.0;
      MolecularStressTensorAccumulated[k][i].ay=0.0;
      MolecularStressTensorAccumulated[k][i].az=0.0;
      MolecularStressTensorAccumulated[k][i].bx=0.0;
      MolecularStressTensorAccumulated[k][i].by=0.0;
      MolecularStressTensorAccumulated[k][i].bz=0.0;
      MolecularStressTensorAccumulated[k][i].cx=0.0;
      MolecularStressTensorAccumulated[k][i].cy=0.0;
      MolecularStressTensorAccumulated[k][i].cz=0.0;
      PressureIdealGasPartAccumulated[k][i]=0.0;
      PressureExcessPartAccumulated[k][i]=0.0;
      PressureTailCorrectionAccumulated[k][i]=0.0;
      PressureAccumulated[k][i]=0.0;

      BoxAccumulated[k][i].x=0.0;
      BoxAccumulated[k][i].y=0.0;
      BoxAccumulated[k][i].z=0.0;

      BoxAXAccumulated[k][i]=0.0;
      BoxAYAccumulated[k][i]=0.0;
      BoxAZAccumulated[k][i]=0.0;
      BoxBXAccumulated[k][i]=0.0;
      BoxBYAccumulated[k][i]=0.0;
      BoxBZAccumulated[k][i]=0.0;
      BoxCXAccumulated[k][i]=0.0;
      BoxCYAccumulated[k][i]=0.0;
      BoxCZAccumulated[k][i]=0.0;

      BoxLengthAccumulated[k][i].x=0.0;
      BoxLengthAccumulated[k][i].y=0.0;
      BoxLengthAccumulated[k][i].z=0.0;
      AlphaAngleAccumulated[k][i]=0.0;
      BetaAngleAccumulated[k][i]=0.0;
      GammaAngleAccumulated[k][i]=0.0;

      VolumeAccumulated[k][i]=0.0;
      VolumeSquaredAccumulated[k][i]=0.0;

      TotalEnergyAccumulated[k][i]=0.0;
      TotalEnergySquaredAccumulated[k][i]=0.0;

      EnthalpyAccumulated[k][i]=0.0;
      EnthalpySquaredAccumulated[k][i]=0.0;

      EnthalpyTimesVolumeAccumulated[k][i]=0.0;
      EnthalpyTimesEnergyAccumulated[k][i]=0.0;

      UNoseHooverAccumulated[k][i]=0.0;

      HeatOfVaporization[k][i]=0.0;
      EnergyPerMolecule[k][i]=0.0;
      VolumePerMolecule[k][i]=0.0;

      InitializeMatrix9x9(&BornTermAccumulated[k][i]);
    }
  }

}

void UpdateEnergyAveragesCurrentSystem(void)
{
  int i,j,nr;
  VECTOR dipole_adsorbates,dipole_cations;
  REAL Enthalpy,densityForComponent;
  REAL NumberOfMolecules;
  REAL PressureIdealGas;
  REAL PressureTail;
  REAL UCFMCAdsorbate,UCFMCCation;
  REAL weight,lambda;
  int numberOfAbsoluteIntegerMoleculesForComponent;
  double numberOfExcessIntegerMoleculesForComponent;

  // check for new block
  if(CurrentCycle==BlockCycle[Block])
    Block++;

  // get the Lambda-weight for CFCMC; the weight=1.0 for CBMC
  weight=CFBiasingWeight();
  lambda=CFBiasingLambda(CurrentSystem,0);

  if (lambda>=0.0 && lambda <= MeasureLambdaBelow)
  {
    BlockWeightedCount[CurrentSystem][Block]+=weight;
    BlockCount[CurrentSystem][Block]+=1.0;

    if(ComputeBornTerm)
      AddBornTermToAverages();

    ConfigurationalStressTensorAccumulated[CurrentSystem][Block].ax+=weight*ConfigurationalStressTensor[CurrentSystem].ax;
    ConfigurationalStressTensorAccumulated[CurrentSystem][Block].bx+=weight*ConfigurationalStressTensor[CurrentSystem].bx;
    ConfigurationalStressTensorAccumulated[CurrentSystem][Block].cx+=weight*ConfigurationalStressTensor[CurrentSystem].cx;

    ConfigurationalStressTensorAccumulated[CurrentSystem][Block].ay+=weight*ConfigurationalStressTensor[CurrentSystem].ay;
    ConfigurationalStressTensorAccumulated[CurrentSystem][Block].by+=weight*ConfigurationalStressTensor[CurrentSystem].by;
    ConfigurationalStressTensorAccumulated[CurrentSystem][Block].cy+=weight*ConfigurationalStressTensor[CurrentSystem].cy;

    ConfigurationalStressTensorAccumulated[CurrentSystem][Block].az+=weight*ConfigurationalStressTensor[CurrentSystem].az;
    ConfigurationalStressTensorAccumulated[CurrentSystem][Block].bz+=weight*ConfigurationalStressTensor[CurrentSystem].bz;
    ConfigurationalStressTensorAccumulated[CurrentSystem][Block].cz+=weight*ConfigurationalStressTensor[CurrentSystem].cz;

    StressTensorAccumulated[CurrentSystem][Block].ax+=weight*StressTensor[CurrentSystem].ax;
    StressTensorAccumulated[CurrentSystem][Block].bx+=weight*StressTensor[CurrentSystem].bx;
    StressTensorAccumulated[CurrentSystem][Block].cx+=weight*StressTensor[CurrentSystem].cx;

    StressTensorAccumulated[CurrentSystem][Block].ay+=weight*StressTensor[CurrentSystem].ay;
    StressTensorAccumulated[CurrentSystem][Block].by+=weight*StressTensor[CurrentSystem].by;
    StressTensorAccumulated[CurrentSystem][Block].cy+=weight*StressTensor[CurrentSystem].cy;

    StressTensorAccumulated[CurrentSystem][Block].az+=weight*StressTensor[CurrentSystem].az;
    StressTensorAccumulated[CurrentSystem][Block].bz+=weight*StressTensor[CurrentSystem].bz;
    StressTensorAccumulated[CurrentSystem][Block].cz+=weight*StressTensor[CurrentSystem].cz;



    UHostBondAccumulated[CurrentSystem][Block]+=weight*UHostBond[CurrentSystem];
    UHostUreyBradleyAccumulated[CurrentSystem][Block]+=weight*UHostUreyBradley[CurrentSystem];
    UHostBendAccumulated[CurrentSystem][Block]+=weight*UHostBend[CurrentSystem];
    UHostInversionBendAccumulated[CurrentSystem][Block]+=weight*UHostInversionBend[CurrentSystem];
    UHostTorsionAccumulated[CurrentSystem][Block]+=weight*UHostTorsion[CurrentSystem];
    UHostImproperTorsionAccumulated[CurrentSystem][Block]+=weight*UHostImproperTorsion[CurrentSystem];
    UHostBondBondAccumulated[CurrentSystem][Block]+=weight*UHostBondBond[CurrentSystem];
    UHostBendBendAccumulated[CurrentSystem][Block]+=weight*UHostBendBend[CurrentSystem];
    UHostBondBendAccumulated[CurrentSystem][Block]+=weight*UHostBondBend[CurrentSystem];
    UHostBondTorsionAccumulated[CurrentSystem][Block]+=weight*UHostBondTorsion[CurrentSystem];
    UHostBendTorsionAccumulated[CurrentSystem][Block]+=weight*UHostBendTorsion[CurrentSystem];

    UCationBondAccumulated[CurrentSystem][Block]+=weight*UCationBond[CurrentSystem];
    UCationUreyBradleyAccumulated[CurrentSystem][Block]+=weight*UCationUreyBradley[CurrentSystem];
    UCationBendAccumulated[CurrentSystem][Block]+=weight*UCationBend[CurrentSystem];
    UCationInversionBendAccumulated[CurrentSystem][Block]+=weight*UCationInversionBend[CurrentSystem];
    UCationTorsionAccumulated[CurrentSystem][Block]+=weight*UCationTorsion[CurrentSystem];
    UCationImproperTorsionAccumulated[CurrentSystem][Block]+=weight*UCationImproperTorsion[CurrentSystem];
    UCationBondBondAccumulated[CurrentSystem][Block]+=weight*UCationBondBond[CurrentSystem];
    UCationBendBendAccumulated[CurrentSystem][Block]+=weight*UCationBendBend[CurrentSystem];
    UCationBondBendAccumulated[CurrentSystem][Block]+=weight*UCationBondBend[CurrentSystem];
    UCationBondTorsionAccumulated[CurrentSystem][Block]+=weight*UCationBondTorsion[CurrentSystem];
    UCationBendTorsionAccumulated[CurrentSystem][Block]+=weight*UCationBendTorsion[CurrentSystem];
    UCationIntraVDWAccumulated[CurrentSystem][Block]+=weight*UCationIntraVDW[CurrentSystem];
    UCationIntraChargeChargeAccumulated[CurrentSystem][Block]+=weight*UCationIntraChargeCharge[CurrentSystem];
    UCationIntraChargeBondDipoleAccumulated[CurrentSystem][Block]+=weight*UCationIntraChargeBondDipole[CurrentSystem];
    UCationIntraBondDipoleBondDipoleAccumulated[CurrentSystem][Block]+=weight*UCationIntraBondDipoleBondDipole[CurrentSystem];

    UAdsorbateBondAccumulated[CurrentSystem][Block]+=weight*UAdsorbateBond[CurrentSystem];
    UAdsorbateUreyBradleyAccumulated[CurrentSystem][Block]+=weight*UAdsorbateUreyBradley[CurrentSystem];
    UAdsorbateBendAccumulated[CurrentSystem][Block]+=weight*UAdsorbateBend[CurrentSystem];
    UAdsorbateInversionBendAccumulated[CurrentSystem][Block]+=weight*UAdsorbateInversionBend[CurrentSystem];
    UAdsorbateTorsionAccumulated[CurrentSystem][Block]+=weight*UAdsorbateTorsion[CurrentSystem];
    UAdsorbateImproperTorsionAccumulated[CurrentSystem][Block]+=weight*UAdsorbateImproperTorsion[CurrentSystem];
    UAdsorbateBondBondAccumulated[CurrentSystem][Block]+=weight*UAdsorbateBondBond[CurrentSystem];
    UAdsorbateBendBendAccumulated[CurrentSystem][Block]+=weight*UAdsorbateBendBend[CurrentSystem];
    UAdsorbateBondBendAccumulated[CurrentSystem][Block]+=weight*UAdsorbateBondBend[CurrentSystem];
    UAdsorbateBondTorsionAccumulated[CurrentSystem][Block]+=weight*UAdsorbateBondTorsion[CurrentSystem];
    UAdsorbateBendTorsionAccumulated[CurrentSystem][Block]+=weight*UAdsorbateBendTorsion[CurrentSystem];
    UAdsorbateIntraVDWAccumulated[CurrentSystem][Block]+=weight*UAdsorbateIntraVDW[CurrentSystem];
    UAdsorbateIntraChargeChargeAccumulated[CurrentSystem][Block]+=weight*UAdsorbateIntraChargeCharge[CurrentSystem];
    UAdsorbateIntraChargeBondDipoleAccumulated[CurrentSystem][Block]+=weight*UAdsorbateIntraChargeBondDipole[CurrentSystem];
    UAdsorbateIntraBondDipoleBondDipoleAccumulated[CurrentSystem][Block]+=weight*UAdsorbateIntraBondDipoleBondDipole[CurrentSystem];

    UHostHostAccumulated[CurrentSystem][Block]+=weight*UHostHost[CurrentSystem];
    UAdsorbateAdsorbateAccumulated[CurrentSystem][Block]+=weight*UAdsorbateAdsorbate[CurrentSystem];
    UCationCationAccumulated[CurrentSystem][Block]+=weight*UCationCation[CurrentSystem];
    UHostAdsorbateAccumulated[CurrentSystem][Block]+=weight*UHostAdsorbate[CurrentSystem];
    UHostCationAccumulated[CurrentSystem][Block]+=weight*UHostCation[CurrentSystem];
    UAdsorbateCationAccumulated[CurrentSystem][Block]+=weight*UAdsorbateCation[CurrentSystem];

    UHostHostVDWAccumulated[CurrentSystem][Block]+=weight*UHostHostVDW[CurrentSystem];
    UAdsorbateAdsorbateVDWAccumulated[CurrentSystem][Block]+=weight*UAdsorbateAdsorbateVDW[CurrentSystem];
    UCationCationVDWAccumulated[CurrentSystem][Block]+=weight*UCationCationVDW[CurrentSystem];
    UHostAdsorbateVDWAccumulated[CurrentSystem][Block]+=weight*UHostAdsorbateVDW[CurrentSystem];
    UHostCationVDWAccumulated[CurrentSystem][Block]+=weight*UHostCationVDW[CurrentSystem];
    UAdsorbateCationVDWAccumulated[CurrentSystem][Block]+=weight*UAdsorbateCationVDW[CurrentSystem];

    UHostHostCoulombAccumulated[CurrentSystem][Block]+=weight*UHostHostCoulomb[CurrentSystem];
    UAdsorbateAdsorbateCoulombAccumulated[CurrentSystem][Block]+=weight*UAdsorbateAdsorbateCoulomb[CurrentSystem];
    UCationCationCoulombAccumulated[CurrentSystem][Block]+=weight*UCationCationCoulomb[CurrentSystem];
    UHostAdsorbateCoulombAccumulated[CurrentSystem][Block]+=weight*UHostAdsorbateCoulomb[CurrentSystem];
    UHostCationCoulombAccumulated[CurrentSystem][Block]+=weight*UHostCationCoulomb[CurrentSystem];
    UAdsorbateCationCoulombAccumulated[CurrentSystem][Block]+=weight*UAdsorbateCationCoulomb[CurrentSystem];

    dipole_adsorbates=ComputeTotalDipoleMomentSystemAdsorbates();
    dipole_cations=ComputeTotalDipoleMomentSystemCations();
    dipole_cations.x=dipole_cations.y=dipole_cations.z=0.0;

    TotalSystemDipoleAccumulated[CurrentSystem][Block].x+=weight*(dipole_adsorbates.x+dipole_cations.x);
    TotalSystemDipoleAccumulated[CurrentSystem][Block].y+=weight*(dipole_adsorbates.y+dipole_cations.y);
    TotalSystemDipoleAccumulated[CurrentSystem][Block].z+=weight*(dipole_adsorbates.z+dipole_cations.z);
    TotalSystemNormDipoleAccumulated[CurrentSystem][Block]+=weight*(sqrt(SQR(dipole_adsorbates.x+dipole_cations.x)+
          SQR(dipole_adsorbates.y+dipole_cations.y)+SQR(dipole_adsorbates.z+dipole_cations.z)));

    TotalSystemDipoleSquaredAccumulated[CurrentSystem][Block].x+=weight*SQR(dipole_adsorbates.x+dipole_cations.x);
    TotalSystemDipoleSquaredAccumulated[CurrentSystem][Block].y+=weight*SQR(dipole_adsorbates.y+dipole_cations.y);
    TotalSystemDipoleSquaredAccumulated[CurrentSystem][Block].z+=weight*SQR(dipole_adsorbates.z+dipole_cations.z);
    TotalSystemNormDipoleSquaredAccumulated[CurrentSystem][Block]+=weight*(SQR(dipole_adsorbates.x+dipole_cations.x)+
          SQR(dipole_adsorbates.y+dipole_cations.y)+SQR(dipole_adsorbates.z+dipole_cations.z));

    UHostPolarizationAccumulated[CurrentSystem][Block]+=weight*UHostPolarization[CurrentSystem];
    UAdsorbatePolarizationAccumulated[CurrentSystem][Block]+=weight*UAdsorbatePolarization[CurrentSystem];
    UCationPolarizationAccumulated[CurrentSystem][Block]+=weight*UCationPolarization[CurrentSystem];
    UHostBackPolarizationAccumulated[CurrentSystem][Block]+=weight*UHostBackPolarization[CurrentSystem];
    UAdsorbateBackPolarizationAccumulated[CurrentSystem][Block]+=weight*UAdsorbateBackPolarization[CurrentSystem];
    UCationBackPolarizationAccumulated[CurrentSystem][Block]+=weight*UCationBackPolarization[CurrentSystem];

    UTailCorrectionAccumulated[CurrentSystem][Block]+=weight*UTailCorrection[CurrentSystem];
    UDistanceConstraintsAccumulated[CurrentSystem][Block]+=weight*UDistanceConstraints[CurrentSystem];
    UAngleConstraintsAccumulated[CurrentSystem][Block]+=weight*UAngleConstraints[CurrentSystem];
    UDihedralConstraintsAccumulated[CurrentSystem][Block]+=weight*UDihedralConstraints[CurrentSystem];
    UInversionBendConstraintsAccumulated[CurrentSystem][Block]+=weight*UInversionBendConstraints[CurrentSystem];
    UOutOfPlaneDistanceConstraintsAccumulated[CurrentSystem][Block]+=weight*UOutOfPlaneDistanceConstraints[CurrentSystem];
    UExclusionConstraintsAccumulated[CurrentSystem][Block]+=weight*UExclusionConstraints[CurrentSystem];

    UCFMCAdsorbate = ComputeEnergyOfFractionalMoleculesAdsorbates();
    UCFMCCation = ComputeEnergyOfFractionalMoleculesCations();
    UTotalAccumulated[CurrentSystem][Block]+=weight*(UTotal[CurrentSystem]-UCFMCAdsorbate-UCFMCCation);

    NumberOfIntegerMoleculesAccumulated[CurrentSystem][Block]+=weight*TotalNumberOfIntegerAdsorbates();

    NumberOfMoleculesSquaredAccumulated[CurrentSystem][Block]+=weight*(SQR(TotalNumberOfIntegerAdsorbates()));
    TotalEnergyTimesNumberOfMoleculesAccumulated[CurrentSystem][Block]+=weight*((UTotal[CurrentSystem]-UCFMCAdsorbate-UCFMCCation)*
                   TotalNumberOfIntegerAdsorbates());

    nr=NumberOfUnitCells[0].x*NumberOfUnitCells[0].y*NumberOfUnitCells[0].z;
    for(i=0;i<NumberOfComponents;i++)
    {
      REAL loading_i = (Components[i].NumberOfMolecules[CurrentSystem]
                      -(Components[i].FractionalMolecule[CurrentSystem]>=0?1:0)
                      -Components[i].NumberOfRXMCMoleculesPresent[CurrentSystem]);
      TotalEnergyTimesNumberOfMoleculesPerComponentAccumulated[CurrentSystem][i][Block]+=weight*((UTotal[CurrentSystem]-UCFMCAdsorbate-UCFMCCation)*loading_i);
      HostAdsorbateEnergyTimesNumberOfMoleculesAccumulated[CurrentSystem][i][Block]+=weight*(UHostAdsorbate[CurrentSystem]*loading_i);
      AdsorbateAdsorbateEnergyTimesNumberOfMoleculesAccumulated[CurrentSystem][i][Block]+=weight*(UAdsorbateAdsorbate[CurrentSystem]*loading_i);
      for(j=0;j<NumberOfComponents;j++)
      {
        REAL loading_j = (Components[j].NumberOfMolecules[CurrentSystem]
                        -(Components[i].FractionalMolecule[CurrentSystem]>=0?1:0)
                        -Components[j].NumberOfRXMCMoleculesPresent[CurrentSystem]);
        NumberOfMoleculesPerComponentSquaredAccumulated[CurrentSystem][i][j][Block]+=weight*loading_i*loading_j;
      }

      numberOfAbsoluteIntegerMoleculesForComponent=(Components[i].NumberOfMolecules[CurrentSystem]
                                                   -(Components[i].FractionalMolecule[CurrentSystem]>=0?1:0)
                                                   -Components[i].NumberOfRXMCMoleculesPresent[CurrentSystem]);
      numberOfExcessIntegerMoleculesForComponent=numberOfAbsoluteIntegerMoleculesForComponent-Components[i].AmountOfExcessMolecules[CurrentSystem];

      NumberOfIntegerMoleculesPerComponentAccumulated[CurrentSystem][i][Block]+=weight*numberOfAbsoluteIntegerMoleculesForComponent;
      NumberOfFractionalMoleculesPerComponentAccumulated[CurrentSystem][i][Block]+=(Components[i].FractionalMolecule[CurrentSystem]>=0?1.0:0.0);
      NumberOfExcessMoleculesPerComponentAccumulated[CurrentSystem][i][Block]+=
               weight*numberOfExcessIntegerMoleculesForComponent;

      densityForComponent=Components[i].Mass*(REAL)numberOfAbsoluteIntegerMoleculesForComponent/Volume[CurrentSystem];
      DensityPerComponentAccumulated[CurrentSystem][i][Block]+=weight*densityForComponent;
      DensityAccumulated[CurrentSystem][Block]+=weight*densityForComponent;
      InverseDensityAccumulated[CurrentSystem][Block]+=weight*(Volume[CurrentSystem]/((REAL)TotalNumberOfIntegerMolecules()+1.0));
      GibbsInverseDensityAccumulated[CurrentSystem][Block]+=weight*(Volume[CurrentSystem]/(TotalNumberOfIntegerMolecules()+1.0));
      GibbsInverseDensityPerComponentAccumulated[CurrentSystem][i][Block]+=weight*(Volume[CurrentSystem]/(numberOfAbsoluteIntegerMoleculesForComponent+1.0));
    }

    TemperatureAccumulated[CurrentSystem][Block]+=weight*(2.0*UKinetic[CurrentSystem]/(K_B*DegreesOfFreedom[CurrentSystem]));
    TemperatureCellAccumulated[CurrentSystem][Block]+=weight*GetCellTemperature();
    TemperatureTranslationAccumulated[CurrentSystem][Block]+=weight*(2.0*(UAdsorbateTranslationalKinetic[CurrentSystem]+
           UCationTranslationalKinetic[CurrentSystem]+UHostKinetic[CurrentSystem])/(K_B*DegreesOfFreedomTranslation[CurrentSystem]));
    TemperatureRotationAccumulated[CurrentSystem][Block]+=weight*(2.0*(UAdsorbateRotationalKinetic[CurrentSystem]+
           UCationRotationalKinetic[CurrentSystem])/(K_B*DegreesOfFreedomRotation[CurrentSystem]));
    TemperatureRotationAdsorbateAccumulated[CurrentSystem][Block]+=weight*(2.0*(UAdsorbateRotationalKinetic[CurrentSystem])/
           (K_B*DegreesOfFreedomRotationalAdsorbates[CurrentSystem]));
    TemperatureTranslationAdsorbateAccumulated[CurrentSystem][Block]+=weight*(2.0*(UAdsorbateTranslationalKinetic[CurrentSystem])/
           (K_B*DegreesOfFreedomTranslationalAdsorbates[CurrentSystem]));

    TemperatureAdsorbatesAccumulated[CurrentSystem][Block]+=weight*(2.0*(UAdsorbateTranslationalKinetic[CurrentSystem]+
                UAdsorbateRotationalKinetic[CurrentSystem])/(K_B*DegreesOfFreedomAdsorbates[CurrentSystem]));
    TemperatureCationsAccumulated[CurrentSystem][Block]+=weight*(2.0*UCationKinetic[CurrentSystem]/
                                               (K_B*DegreesOfFreedomCations[CurrentSystem]));
    TemperatureFrameworkAccumulated[CurrentSystem][Block]+=weight*(2.0*UHostKinetic[CurrentSystem]/
                                               (K_B*DegreesOfFreedomFramework[CurrentSystem]));

    NumberOfMolecules=TotalNumberOfIntegerMolecules();

    if(ComputeMolecularPressure[CurrentSystem])
    {
      ComputeMolecularPressureTensor(&MolecularStressTensor[CurrentSystem],&PressureIdealGas,&PressureTail);

      PressureIdealGasPartAccumulated[CurrentSystem][Block]+=weight*PressureIdealGas;
      PressureExcessPartAccumulated[CurrentSystem][Block]-=weight*((MolecularStressTensor[CurrentSystem].ax+MolecularStressTensor[CurrentSystem].by+
                  MolecularStressTensor[CurrentSystem].cz)/(3.0*Volume[CurrentSystem]));
      PressureTailCorrectionAccumulated[CurrentSystem][Block]+=weight*PressureTail;

      MolecularStressTensor[CurrentSystem].ax=(PressureIdealGas+PressureTail-MolecularStressTensor[CurrentSystem].ax/Volume[CurrentSystem]);
      MolecularStressTensor[CurrentSystem].ay=(-MolecularStressTensor[CurrentSystem].ay/Volume[CurrentSystem]);
      MolecularStressTensor[CurrentSystem].az=(-MolecularStressTensor[CurrentSystem].az/Volume[CurrentSystem]);

      MolecularStressTensor[CurrentSystem].bx=(-MolecularStressTensor[CurrentSystem].bx/Volume[CurrentSystem]);
      MolecularStressTensor[CurrentSystem].by=(PressureIdealGas+PressureTail-MolecularStressTensor[CurrentSystem].by/Volume[CurrentSystem]);
      MolecularStressTensor[CurrentSystem].bz=(-MolecularStressTensor[CurrentSystem].bz/Volume[CurrentSystem]);

      MolecularStressTensor[CurrentSystem].cx=(-MolecularStressTensor[CurrentSystem].cx/Volume[CurrentSystem]);
      MolecularStressTensor[CurrentSystem].cy=(-MolecularStressTensor[CurrentSystem].cy/Volume[CurrentSystem]);
      MolecularStressTensor[CurrentSystem].cz=(PressureIdealGas+PressureTail-MolecularStressTensor[CurrentSystem].cz/Volume[CurrentSystem]);

      MolecularStressTensorAccumulated[CurrentSystem][Block].ax+=weight*MolecularStressTensor[CurrentSystem].ax;
      MolecularStressTensorAccumulated[CurrentSystem][Block].ay+=weight*MolecularStressTensor[CurrentSystem].ay;
      MolecularStressTensorAccumulated[CurrentSystem][Block].az+=weight*MolecularStressTensor[CurrentSystem].az;

      MolecularStressTensorAccumulated[CurrentSystem][Block].bx+=weight*MolecularStressTensor[CurrentSystem].bx;
      MolecularStressTensorAccumulated[CurrentSystem][Block].by+=weight*MolecularStressTensor[CurrentSystem].by;
      MolecularStressTensorAccumulated[CurrentSystem][Block].bz+=weight*MolecularStressTensor[CurrentSystem].bz;

      MolecularStressTensorAccumulated[CurrentSystem][Block].cx+=weight*MolecularStressTensor[CurrentSystem].cx;
      MolecularStressTensorAccumulated[CurrentSystem][Block].cy+=weight*MolecularStressTensor[CurrentSystem].cy;
      MolecularStressTensorAccumulated[CurrentSystem][Block].cz+=weight*MolecularStressTensor[CurrentSystem].cz;
    }

    UNoseHooverAccumulated[CurrentSystem][Block]+=weight*UNoseHoover[CurrentSystem];

    BoxAccumulated[CurrentSystem][Block].x+=weight*Box[CurrentSystem].ax;
    BoxAccumulated[CurrentSystem][Block].y+=weight*Box[CurrentSystem].by;
    BoxAccumulated[CurrentSystem][Block].z+=weight*Box[CurrentSystem].cz;

    BoxAXAccumulated[CurrentSystem][Block]+=weight*Box[CurrentSystem].ax;
    BoxAYAccumulated[CurrentSystem][Block]+=weight*Box[CurrentSystem].ay;
    BoxAZAccumulated[CurrentSystem][Block]+=weight*Box[CurrentSystem].az;
    BoxBXAccumulated[CurrentSystem][Block]+=weight*Box[CurrentSystem].bx;
    BoxBYAccumulated[CurrentSystem][Block]+=weight*Box[CurrentSystem].by;
    BoxBZAccumulated[CurrentSystem][Block]+=weight*Box[CurrentSystem].bz;
    BoxCXAccumulated[CurrentSystem][Block]+=weight*Box[CurrentSystem].cx;
    BoxCYAccumulated[CurrentSystem][Block]+=weight*Box[CurrentSystem].cy;
    BoxCZAccumulated[CurrentSystem][Block]+=weight*Box[CurrentSystem].cz;

    BoxLengthAccumulated[CurrentSystem][Block].x+=weight*BoxProperties[CurrentSystem].ax;
    BoxLengthAccumulated[CurrentSystem][Block].y+=weight*BoxProperties[CurrentSystem].ay;
    BoxLengthAccumulated[CurrentSystem][Block].z+=weight*BoxProperties[CurrentSystem].az;

    AlphaAngleAccumulated[CurrentSystem][Block]+=weight*AlphaAngle[CurrentSystem]*RAD2DEG;
    BetaAngleAccumulated[CurrentSystem][Block]+=weight*BetaAngle[CurrentSystem]*RAD2DEG;
    GammaAngleAccumulated[CurrentSystem][Block]+=weight*GammaAngle[CurrentSystem]*RAD2DEG;

    VolumeAccumulated[CurrentSystem][Block]+=weight*Volume[CurrentSystem];

    VolumeSquaredAccumulated[CurrentSystem][Block]+=weight*SQR(Volume[CurrentSystem]);

    TotalEnergyAccumulated[CurrentSystem][Block]+=weight*UTotal[CurrentSystem];
    TotalEnergySquaredAccumulated[CurrentSystem][Block]+=weight*SQR(UTotal[CurrentSystem]);

    Enthalpy=(UTotal[CurrentSystem]-UCFMCAdsorbate-UCFMCCation)+Volume[CurrentSystem]*therm_baro_stats.ExternalPressure[CurrentSystem][0];
    EnthalpyAccumulated[CurrentSystem][Block]+=weight*Enthalpy;
    EnthalpySquaredAccumulated[CurrentSystem][Block]+=weight*SQR(Enthalpy);

    EnthalpyTimesVolumeAccumulated[CurrentSystem][Block]+=weight*Enthalpy*Volume[CurrentSystem];
    EnthalpyTimesEnergyAccumulated[CurrentSystem][Block]+=weight*Enthalpy*(UTotal[CurrentSystem]-UCFMCAdsorbate-UCFMCCation);

    HeatOfVaporization[CurrentSystem][Block]+=weight*(therm_baro_stats.ExternalTemperature[CurrentSystem]-
                                (UAdsorbateAdsorbate[CurrentSystem]+UCationCation[CurrentSystem])/
                                (NumberOfAdsorbateMolecules[CurrentSystem]+NumberOfCationMolecules[CurrentSystem]));
    if(NumberOfMolecules>0)
    {
      EnergyPerMolecule[CurrentSystem][Block]+=weight*UTotal[CurrentSystem]/NumberOfMolecules;
      VolumePerMolecule[CurrentSystem][Block]+=weight*Volume[CurrentSystem]/NumberOfMolecules;
    }

    if(ComputePrincipleMomentsOfInertia)
      MeasurePrincipleMomentsOfInertia();

    CompressibilityAccumulated[CurrentSystem][Block]+=weight*(((MolecularStressTensor[CurrentSystem].ax+MolecularStressTensor[CurrentSystem].by+MolecularStressTensor[CurrentSystem].cz)/3.0)*
             Volume[CurrentSystem]*Beta[CurrentSystem]/NumberOfMolecules);

  }
}

REAL GetAverageProperty(REAL **Property)
{
  int i;
  REAL sum1,sum2;

  sum1=sum2=0.0;
  for(i=0;i<NR_BLOCKS;i++)
    if(BlockCount[CurrentSystem][i]>0.0)
    {
      sum1+=Property[CurrentSystem][i];
      sum2+=BlockCount[CurrentSystem][i];
    }
  if(sum2>0.0)
    return sum1/sum2;
  else
    return 0.0;
}

REAL GetAverageWeightedProperty(REAL **Property)
{
  int i;
  REAL sum1,sum2;

  sum1=sum2=0.0;
  for(i=0;i<NR_BLOCKS;i++)
    if(BlockWeightedCount[CurrentSystem][i]>0.0)
    {
      sum1+=Property[CurrentSystem][i];
      sum2+=BlockWeightedCount[CurrentSystem][i];
    }
  if(sum2>0.0)
    return sum1/sum2;
  else
    return 0.0;
}

VECTOR GetAverageVectorProperty(VECTOR **Property)
{
  int i;
  VECTOR sum1;
  REAL sum2;

  sum1.x=sum1.y=sum1.z=sum2=0.0;
  for(i=0;i<NR_BLOCKS;i++)
    if(BlockWeightedCount[CurrentSystem][i]>0.0)
    {
      sum1.x+=Property[CurrentSystem][i].x;
      sum1.y+=Property[CurrentSystem][i].y;
      sum1.z+=Property[CurrentSystem][i].z;
      sum2+=BlockWeightedCount[CurrentSystem][i];
    }
  if(sum2>0.0)
  {
    sum1.x/=sum2;
    sum1.y/=sum2;
    sum1.z/=sum2;
    return sum1;
  }
  else
  {
    sum1.x=sum1.y=sum1.z=0.0;
    return sum1;
  }
}

REAL_MATRIX9x9 GetAverageRealMatrix9x9Property(REAL_MATRIX9x9 **Property)
{
  int i;
  REAL_MATRIX9x9 sum1;
  REAL sum2;

  sum2=0.0;
  InitializeMatrix9x9(&sum1);
  for(i=0;i<NR_BLOCKS;i++)
    if(BlockWeightedCount[CurrentSystem][i]>0.0)
    {
      AddRealMatrix9x9(&sum1,sum1,Property[CurrentSystem][i]);
      sum2+=BlockWeightedCount[CurrentSystem][i];
    }
  if(sum2>0.0)
    DivideRealMatrix9x9ByReal(&sum1,sum1,sum2);
  return sum1;
}

REAL GetAverageComponentProperty(REAL ***Property,int comp)
{
  int i;
  REAL sum1,sum2;

  sum1=sum2=0.0;
  for(i=0;i<NR_BLOCKS;i++)
    if(BlockCount[CurrentSystem][i]>0.0)
    {
      sum1+=Property[CurrentSystem][comp][i];
      sum2+=BlockCount[CurrentSystem][i];
    }
  if(sum2>0.0)
    return sum1/sum2;
  else
    return 0.0;
}


REAL GetAverageWeightedComponentProperty(REAL ***Property,int comp)
{
  int i;
  REAL sum1,sum2;

  sum1=sum2=0.0;
  for(i=0;i<NR_BLOCKS;i++)
    if(BlockWeightedCount[CurrentSystem][i]>0.0)
    {
      sum1+=Property[CurrentSystem][comp][i];
      sum2+=BlockWeightedCount[CurrentSystem][i];
    }
  if(sum2>0.0)
    return sum1/sum2;
  else
    return 0.0;
}

VECTOR GetAverageComponentVectorProperty(VECTOR ***Property,int comp)
{
  int i;
  VECTOR sum1;
  REAL sum2;

  sum1.x=sum1.y=sum1.z=0.0;
  sum2=0.0;
  for(i=0;i<NR_BLOCKS;i++)
    if(BlockWeightedCount[CurrentSystem][i]>0.0)
    {
      sum1.x+=Property[CurrentSystem][comp][i].x;
      sum1.y+=Property[CurrentSystem][comp][i].y;
      sum1.z+=Property[CurrentSystem][comp][i].z;
      sum2+=BlockWeightedCount[CurrentSystem][i];
    }
  if(sum2>0.0)
  {
    sum1.x/=sum2;
    sum1.y/=sum2;
    sum1.z/=sum2;
    return sum1;
  }
  else
  {
    sum1.x=sum1.y=sum1.z=0.0;
    return sum1;
  }
}



void PrintIntervalStatusInit(long long CurrentCycle,long long NumberOfCycles,FILE *FilePtr)
{
  int i,j,k;
  REAL number_of_unit_cells;
  int FractionalMolecule;
  REAL loading,average_loading,Lambda,shift;
  VECTOR com;

  fprintf(FilePtr,"[Init] Current cycle: %lld out of %lld\n",CurrentCycle,NumberOfCycles);
  fprintf(FilePtr,"========================================================================================================\n\n");
  fprintf(FilePtr,"Net charge: %g (F: %g, A: %g, C: %g)\n",NetChargeSystem[CurrentSystem],NetChargeFramework[CurrentSystem],
                                                           NetChargeAdsorbates[CurrentSystem],NetChargeCations[CurrentSystem]);

  // write out the boxlengths
  fprintf(FilePtr,"Current Box: %9.5lf %9.5lf %9.5lf [A]\n",
          (double)Box[CurrentSystem].ax,(double)Box[CurrentSystem].bx,(double)Box[CurrentSystem].cx);
  fprintf(FilePtr,"             %9.5lf %9.5lf %9.5lf [A]\n",
          (double)Box[CurrentSystem].ay,(double)Box[CurrentSystem].by,(double)Box[CurrentSystem].cy);
  fprintf(FilePtr,"             %9.5lf %9.5lf %9.5lf [A]\n",
          (double)Box[CurrentSystem].az,(double)Box[CurrentSystem].bz,(double)Box[CurrentSystem].cz);
  fprintf(FilePtr,"Box-lengths: %9.5lf %9.5lf %9.5lf Box-angles:  %9.5lf %9.5lf %9.5lf [degrees]\n",
          (double)(BoxProperties[CurrentSystem].ax),
          (double)(BoxProperties[CurrentSystem].ay),
          (double)(BoxProperties[CurrentSystem].az),
          (double)(AlphaAngle[CurrentSystem]*RAD2DEG),
          (double)(BetaAngle[CurrentSystem]*RAD2DEG),
          (double)(GammaAngle[CurrentSystem]*RAD2DEG));
  fprintf(FilePtr,"Volume: %9.5lf [A^3]\n\n",(double)Volume[CurrentSystem]);
  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
  {
    com=GetFrameworkCenterOfMass();
    fprintf(FilePtr,"Framework center-of-mass drift: %18.10lf [A] %18.10lf [A] %18.10lf [A]\n",
      (double)(com.x-Framework[CurrentSystem].IntialCenterOfMassPosition.x),
      (double)(com.y-Framework[CurrentSystem].IntialCenterOfMassPosition.y),
      (double)(com.z-Framework[CurrentSystem].IntialCenterOfMassPosition.z));
  }

  if(NumberOfReactions>0)
  {
    fprintf(FilePtr,"Reactions:\n");
    fprintf(FilePtr,"----------------------------------------------------------------------------------------------------------------------------------------------------\n");
    for(i=0;i<NumberOfReactions;i++)
    {
      fprintf(FilePtr,"Reaction %d, current Lambda: %18.10f, maximum Lambda-change: %18.10f\n",i,CFCRXMCLambda[CurrentSystem][i],MaximumReactionLambdaChange[CurrentSystem][i]);
      fprintf(FilePtr,"\tFractional molecules: ");
      for(j=0;j<NumberOfComponents;j++)
      {
        if(ProductsStoichiometry[i][j]>0)
        {
          for(k=0;k<ProductsStoichiometry[i][j];k++)
            fprintf(FilePtr," %d (%s)",Components[j].ProductFractionalMolecules[CurrentSystem][i][k],Components[j].Name);
        }
      }
      fprintf(FilePtr," <-->");
      for(j=0;j<NumberOfComponents;j++)
      {
        if(ReactantsStoichiometry[i][j]>0)
        {
          for(k=0;k<ReactantsStoichiometry[i][j];k++)
            fprintf(FilePtr," %d (%s)",Components[j].ReactantFractionalMolecules[CurrentSystem][i][k],Components[j].Name);
        }
      }
      fprintf(FilePtr,"\n");

      shift=RXMCBiasingFactors[CurrentSystem][i][0];
      fprintf(FilePtr,"\tBiasing Factors: ");
      for(k=0;k<RXMCLambdaHistogramSize;k++)
      {
        fprintf(FilePtr,"%4f ",RXMCBiasingFactors[CurrentSystem][i][k]-shift);
        if((k+1)%10==0&&(k+1)!=RXMCLambdaHistogramSize)  fprintf(FilePtr,"\n\t                 ");
      }
      fprintf(FilePtr,"\n");
    }
    fprintf(FilePtr,"\n");
  }

  number_of_unit_cells=NumberOfUnitCells[CurrentSystem].x*NumberOfUnitCells[CurrentSystem].y*NumberOfUnitCells[CurrentSystem].z;
  if(Framework[CurrentSystem].FrameworkModel==NONE)
  {
    fprintf(FilePtr,"Amount of molecules per component:\n");
    fprintf(FilePtr,"----------------------------------------------------------------------------------------------------------------------------------------------------\n");
    for(i=0;i<NumberOfComponents;i++)
    {
      fprintf(FilePtr,"Component %d (%s), current number of integer/fractional/reaction molecules: %d/%d/%d, density: %9.5lf [kg/m^3]\n",
        i,
        Components[i].Name,
        Components[i].NumberOfMolecules[CurrentSystem]-(Components[i].FractionalMolecule[CurrentSystem]>=0?1:0)-Components[i].NumberOfRXMCMoleculesPresent[CurrentSystem],
        (Components[i].FractionalMolecule[CurrentSystem]>=0?1:0),
        Components[i].NumberOfRXMCMoleculesPresent[CurrentSystem],
        (double)((Components[i].Mass*(REAL)(Components[i].NumberOfMolecules[CurrentSystem]
                 -(Components[i].FractionalMolecule[CurrentSystem]>=0?1:0)
                 -Components[i].NumberOfRXMCMoleculesPresent[CurrentSystem])/
                  Volume[CurrentSystem])*DENSITY_CONVERSION_FACTOR));

      if(Components[i].CFMoleculePresent[CurrentSystem])
      {
        char FractionalMoleculeString[6];
        FractionalMolecule=Components[i].FractionalMolecule[CurrentSystem];
        if (FractionalMolecule<0)
          strcpy(FractionalMoleculeString,"-");
        else
          sprintf(FractionalMoleculeString,"%d", FractionalMolecule);

        fprintf(FilePtr,"\tFractional molecule-id: %s, max. Lambda-change: %5lf (CFMC) %5lf (CB/CFMC)\n",FractionalMoleculeString,MaximumCFLambdaChange[CurrentSystem][i],MaximumCBCFLambdaChange[CurrentSystem][i]);
        fprintf(FilePtr,"\tLambda factors: ");
        for(k=0;k<Components[i].NumberOfAtoms;k++)
        {
          Lambda=0.0;
          if(FractionalMolecule>=0)
          {
            if(Components[i].ExtraFrameworkMolecule)
              Lambda=Cations[CurrentSystem][FractionalMolecule].Atoms[k].CFVDWScalingParameter;
            else
              Lambda=Adsorbates[CurrentSystem][FractionalMolecule].Atoms[k].CFVDWScalingParameter;
          }
          fprintf(FilePtr,"%4f ",Lambda);
          if((k+1)%10==0&&(k+1)!=Components[i].NumberOfAtoms)  fprintf(FilePtr,"\n\t                ");
        }
        fprintf(FilePtr,"\n");
        shift=Components[i].CFBiasingFactors[CurrentSystem][0];
        shift=0.0;
        fprintf(FilePtr,"\tBiasing Factors: ");
        for(k=0;k<Components[i].CFLambdaHistogramSize;k++)
        {
          fprintf(FilePtr,"%4f ",Components[i].CFBiasingFactors[CurrentSystem][k]-shift);
          if((k+1)%10==0&&(k+1)!=Components[i].CFLambdaHistogramSize)  fprintf(FilePtr,"\n\t                 ");
        }
        fprintf(FilePtr,"\n");
      }
    }
    fprintf(FilePtr,"----------------------------------------------------------------------------------------------------------------------------------------------------\n");
  }
  else
  {
    fprintf(FilePtr,"Loadings per component:\n");
    fprintf(FilePtr,"----------------------------------------------------------------------------------------------------------------------------------------------------\n");
    for(i=0;i<NumberOfComponents;i++)
    {
      fprintf(FilePtr,"Component %d (%s), current number of integer/fractional/reaction molecules: %d/%d/%d, density: %9.5lf [kg/m^3]\n",
        i,
        Components[i].Name,
        Components[i].NumberOfMolecules[CurrentSystem]-(Components[i].FractionalMolecule[CurrentSystem]>=0?1:0)-Components[i].NumberOfRXMCMoleculesPresent[CurrentSystem],
        (Components[i].FractionalMolecule[CurrentSystem]>=0?1:0),
        Components[i].NumberOfRXMCMoleculesPresent[CurrentSystem],
        (double)(Components[i].Mass*(REAL)(Components[i].NumberOfMolecules[CurrentSystem]
                                           -(Components[i].FractionalMolecule[CurrentSystem]>=0?1:0)
                                           -Components[i].NumberOfRXMCMoleculesPresent[CurrentSystem])
                                           /Volume[CurrentSystem])*DENSITY_CONVERSION_FACTOR);
      if(Components[i].CFMoleculePresent[CurrentSystem])
      {
        FractionalMolecule=Components[i].FractionalMolecule[CurrentSystem];
        fprintf(FilePtr,"\tFractional molecule-id: %d, max. Lambda-change: %5lf (CFMC) %5lf (CB/CFMC)\n",FractionalMolecule,MaximumCFLambdaChange[CurrentSystem][i],MaximumCBCFLambdaChange[CurrentSystem][i]);
        fprintf(FilePtr,"\tLambda factors: ");
        for(k=0;k<Components[i].NumberOfAtoms;k++)
        {
          if(Components[i].ExtraFrameworkMolecule)
            Lambda=Cations[CurrentSystem][FractionalMolecule].Atoms[k].CFVDWScalingParameter;
          else
            Lambda=Adsorbates[CurrentSystem][FractionalMolecule].Atoms[k].CFVDWScalingParameter;
          fprintf(FilePtr,"%4f ",Lambda);
          if((k+1)%10==0&&(k+1)!=Components[i].NumberOfAtoms)  fprintf(FilePtr,"\n\t                ");
        }
        fprintf(FilePtr,"\n");
        shift=Components[i].CFBiasingFactors[CurrentSystem][0];
        shift=0.0;     
        fprintf(FilePtr,"\tBiasing Factors: ");
        for(k=0;k<Components[i].CFLambdaHistogramSize;k++)
        {
          fprintf(FilePtr,"%4f ",Components[i].CFBiasingFactors[CurrentSystem][k]-shift);
          if((k+1)%10==0&&(k+1)!=Components[i].CFLambdaHistogramSize)  fprintf(FilePtr,"\n\t                 ");
        }
        fprintf(FilePtr,"\n");
      }

      loading=(REAL)(Components[i].NumberOfMolecules[CurrentSystem]
                    -(Components[i].FractionalMolecule[CurrentSystem]>=0?1:0)
                    -Components[i].NumberOfRXMCMoleculesPresent[CurrentSystem])/(REAL)number_of_unit_cells;
      average_loading=GetAverageWeightedComponentProperty(NumberOfIntegerMoleculesPerComponentAccumulated,i)/(REAL)number_of_unit_cells;
      fprintf(FilePtr,"\tabsolute adsorption: %9.5lf [mol/uc], %14.4lf [mol/kg],      %14.4lf [mg/g]\n",
         (double)loading,
         (double)(Components[i].MOLEC_PER_UC_TO_MOL_PER_KG[CurrentSystem]*loading),
         (double)(Components[i].MOLEC_PER_UC_TO_MILLIGRAM_PER_GRAM_OF_FRAMEWORK[CurrentSystem]*loading));
      fprintf(FilePtr,"\t                                         %14.4lf [cm^3 STP/g],  %14.4lf [cm^3 STP/cm^3]\n",
         (double)(Components[i].MOLEC_PER_UC_TO_CC_STP_G[CurrentSystem]*loading),
         (double)(Components[i].MOLEC_PER_UC_TO_CC_STP_CC[CurrentSystem]*loading));

      loading=(REAL)(Components[i].NumberOfMolecules[CurrentSystem]
                    -Components[i].AmountOfExcessMolecules[CurrentSystem]
                    -(Components[i].FractionalMolecule[CurrentSystem]>=0?1:0)
                    -Components[i].NumberOfRXMCMoleculesPresent[CurrentSystem])/(REAL)number_of_unit_cells;
      average_loading=GetAverageWeightedComponentProperty(NumberOfExcessMoleculesPerComponentAccumulated,i)/(REAL)number_of_unit_cells;
      fprintf(FilePtr,"\texcess adsorption:   %9.5lf [mol/uc], %14.4lf [mol/kg],      %14.4lf [mg/g]\n",
         (double)loading,
         (double)(Components[i].MOLEC_PER_UC_TO_MOL_PER_KG[CurrentSystem]*loading),
         (double)(Components[i].MOLEC_PER_UC_TO_MILLIGRAM_PER_GRAM_OF_FRAMEWORK[CurrentSystem]*loading));
      fprintf(FilePtr,"\t                                         %14.4lf [cm^3 STP/g],  %14.4lf [cm^3 STP/cm^3]\n",
         (double)(Components[i].MOLEC_PER_UC_TO_CC_STP_G[CurrentSystem]*loading),
         (double)(Components[i].MOLEC_PER_UC_TO_CC_STP_CC[CurrentSystem]*loading));

      fprintf(FilePtr,"----------------------------------------------------------------------------------------------------------------------------------------------------\n");
    }
  }

  fprintf(FilePtr,"Degrees of freedom: %d %d %d %d\n",DegreesOfFreedom[CurrentSystem],DegreesOfFreedomFramework[CurrentSystem],
          DegreesOfFreedomAdsorbates[CurrentSystem],DegreesOfFreedomCations[CurrentSystem]);
  fprintf(FilePtr,"Number of Framework-atoms: %6d\n",Framework[CurrentSystem].TotalNumberOfAtoms);
  fprintf(FilePtr,"Number of Adsorbates:      %6d (%d integer, %d fractional, %d reaction)\n",
          NumberOfAdsorbateMolecules[CurrentSystem],
          TotalNumberOfIntegerAdsorbates(),
          TotalNumberOfFractionalAdsorbates(),
          NumberOfReactionAdsorbateMolecules[CurrentSystem]);
  fprintf(FilePtr,"Number of Cations:         %6d (%d integer, %d fractional, %d reaction\n",
      NumberOfCationMolecules[CurrentSystem],
      TotalNumberOfIntegerCations(),
      TotalNumberOfFractionalCations(),
      NumberOfReactionCationMolecules[CurrentSystem]);
  fprintf(FilePtr,"\n");

  fprintf(FilePtr,"Current total potential energy:       % 22.10lf [K]\n",
      (double)(UTotal[CurrentSystem]*ENERGY_TO_KELVIN));

  fprintf(FilePtr,"\tCurrent Host-Host energy:           % 22.10lf [K]\n",(double)(UHostHost[CurrentSystem]*ENERGY_TO_KELVIN));
  fprintf(FilePtr,"\tCurrent Host-Adsorbate energy:      % 22.10lf [K]\n",(double)(UHostAdsorbate[CurrentSystem]*ENERGY_TO_KELVIN));
  fprintf(FilePtr,"\tCurrent Host-Cation energy:         % 22.10lf [K]\n",(double)(UHostCation[CurrentSystem]*ENERGY_TO_KELVIN));
  fprintf(FilePtr,"\tCurrent Adsorbate-Adsorbate energy: % 22.10lf [K]\n",(double)(UAdsorbateAdsorbate[CurrentSystem]*ENERGY_TO_KELVIN));
  fprintf(FilePtr,"\tCurrent Cation-Cation energy:       % 22.10lf [K]\n",(double)(UCationCation[CurrentSystem]*ENERGY_TO_KELVIN));
  fprintf(FilePtr,"\tCurrent Adsorbate-Cation energy:    % 22.10lf [K]\n",(double)(UAdsorbateCation[CurrentSystem]*ENERGY_TO_KELVIN));

  if(ComputePolarization)
  {
    fprintf(FilePtr,"\tCurrent polarization energy:        % 22.10lf [K]\n",(double)(UHostPolarization[CurrentSystem]+
         UAdsorbatePolarization[CurrentSystem]+UCationPolarization[CurrentSystem])*ENERGY_TO_KELVIN*ENERGY_TO_KELVIN);
    fprintf(FilePtr,"\tCurrent back-polarization energy:   % 22.10lf [K]\n",(double)(UHostBackPolarization[CurrentSystem]+
         UAdsorbateBackPolarization[CurrentSystem]+UCationBackPolarization[CurrentSystem])*ENERGY_TO_KELVIN*ENERGY_TO_KELVIN);
  }

  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
  {
    if(Framework[CurrentSystem].NumberOfBondsDefinitions>0)
      fprintf(FilePtr,"\t\tCurrent Host-Bond energy:             % 22.10lf [K]\n",
          (double)UHostBond[CurrentSystem]*ENERGY_TO_KELVIN);
    if(Framework[CurrentSystem].NumberOfUreyBradleyDefinitions>0)
      fprintf(FilePtr,"\t\tCurrent Host-UreyBradly energy:       % 22.10lf [K]\n",
          (double)UHostUreyBradley[CurrentSystem]*ENERGY_TO_KELVIN);
    if(Framework[CurrentSystem].NumberOfBendDefinitions>0)
      fprintf(FilePtr,"\t\tCurrent Host-Bend energy:             % 22.10lf [K]\n",
          (double)UHostBend[CurrentSystem]*ENERGY_TO_KELVIN);
    if(Framework[CurrentSystem].NumberOfInversionBendDefinitions>0)
      fprintf(FilePtr,"\t\tCurrent Host-Inversion Bend energy:   % 22.10lf [K]\n",
          (double)UHostInversionBend[CurrentSystem]*ENERGY_TO_KELVIN);
    if(Framework[CurrentSystem].NumberOfTorsionDefinitions>0)
      fprintf(FilePtr,"\t\tCurrent Host-Torsion energy:          % 22.10lf [K]\n",
          (double)UHostTorsion[CurrentSystem]*ENERGY_TO_KELVIN);
    if(Framework[CurrentSystem].NumberOfImproperTorsionDefinitions>0)
      fprintf(FilePtr,"\t\tCurrent Host-Improper torsion energy: % 22.10lf [K]\n",
          (double)UHostImproperTorsion[CurrentSystem]*ENERGY_TO_KELVIN);
    if(Framework[CurrentSystem].NumberOfBondBondDefinitions>0)
      fprintf(FilePtr,"\t\tCurrent Host-Bond/Bond energy:        % 22.10lf [K]\n",
          (double)UHostBondBond[CurrentSystem]*ENERGY_TO_KELVIN);
    if(Framework[CurrentSystem].NumberOfBendBendDefinitions>0)
      fprintf(FilePtr,"\t\tCurrent Host-Bend/Bend energy:        % 22.10lf [K]\n",
          (double)UHostBendBend[CurrentSystem]*ENERGY_TO_KELVIN);
    if(Framework[CurrentSystem].NumberOfBondBendDefinitions>0)
      fprintf(FilePtr,"\t\tCurrent Host-Bond/Bend energy:        % 22.10lf [K]\n",
          (double)UHostBondBend[CurrentSystem]*ENERGY_TO_KELVIN);
    if(Framework[CurrentSystem].NumberOfBondTorsionDefinitions>0)
      fprintf(FilePtr,"\t\tCurrent Host-Bond/Torsion energy:     % 22.10lf [K]\n",
          (double)UHostBondTorsion[CurrentSystem]*ENERGY_TO_KELVIN);
    if(Framework[CurrentSystem].NumberOfBendTorsionDefinitions>0)
      fprintf(FilePtr,"\t\tCurrent Host-Bend/Torsion energy:     % 22.10lf [K]\n",
          (double)UHostBendTorsion[CurrentSystem]*ENERGY_TO_KELVIN);
  }
  fprintf(FilePtr,"\n");
  PrintWarningStatus();
  fprintf(FilePtr,"\n\n");
  fflush(FilePtr);
}

void PrintIntervalStatusEquilibration(long long CurrentCycle,long long NumberOfCycles,FILE *FilePtr)
{
  int i,j,k;
  REAL number_of_unit_cells;
  int FractionalMolecule;
  REAL loading,average_loading,shift,Lambda;
  VECTOR com;

  fprintf(FilePtr,"[Equilibration] Current cycle: %lld out of %lld\n",CurrentCycle,NumberOfCycles);
  fprintf(FilePtr,"========================================================================================================\n\n");
  fprintf(FilePtr,"Net charge: %g (F: %g, A: %g, C: %g)\n",NetChargeSystem[CurrentSystem],NetChargeFramework[CurrentSystem],
                                                           NetChargeAdsorbates[CurrentSystem],NetChargeCations[CurrentSystem]);

  // write out the boxlengths
  fprintf(FilePtr,"Current Box: %9.5lf %9.5lf %9.5lf [A]\n",
          (double)Box[CurrentSystem].ax,(double)Box[CurrentSystem].bx,(double)Box[CurrentSystem].cx);
  fprintf(FilePtr,"             %9.5lf %9.5lf %9.5lf [A]\n",
          (double)Box[CurrentSystem].ay,(double)Box[CurrentSystem].by,(double)Box[CurrentSystem].cy);
  fprintf(FilePtr,"             %9.5lf %9.5lf %9.5lf [A]\n",
          (double)Box[CurrentSystem].az,(double)Box[CurrentSystem].bz,(double)Box[CurrentSystem].cz);
  fprintf(FilePtr,"Box-lengths: %9.5lf %9.5lf %9.5lf Box-angles:  %9.5lf %9.5lf %9.5lf [degrees]\n",
          (double)(BoxProperties[CurrentSystem].ax),
          (double)(BoxProperties[CurrentSystem].ay),
          (double)(BoxProperties[CurrentSystem].az),
          (double)(AlphaAngle[CurrentSystem]*RAD2DEG),
          (double)(BetaAngle[CurrentSystem]*RAD2DEG),
          (double)(GammaAngle[CurrentSystem]*RAD2DEG));
  fprintf(FilePtr,"Volume: %9.5lf [A^3]\n\n",(double)Volume[CurrentSystem]);
  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
  {
    com=GetFrameworkCenterOfMass();
    fprintf(FilePtr,"Framework center-of-mass drift: %18.10lf [A] %18.10lf [A] %18.10lf [A]\n",
      (double)(com.x-Framework[CurrentSystem].IntialCenterOfMassPosition.x),
      (double)(com.y-Framework[CurrentSystem].IntialCenterOfMassPosition.y),
      (double)(com.z-Framework[CurrentSystem].IntialCenterOfMassPosition.z));
  }
  if(SimulationType==MOLECULAR_DYNAMICS)
  {
    com=MeasureVelocityDrift();
    fprintf(FilePtr,"Total linear momentum: %18.10lf [A] %18.10lf [A] %18.10lf [A]\n",
           (double)com.x,(double)com.y,(double)com.z);
    if(Framework[CurrentSystem].NumberOfCoreShellDefinitions>0)
      fprintf(FilePtr,"Core/Shell temperature: %g\n",GetCoreShellTemperature());
    fprintf(FilePtr,"\n");

  }

  if(NumberOfReactions>0)
  {
    fprintf(FilePtr,"Reactions:\n");
    fprintf(FilePtr,"----------------------------------------------------------------------------------------------------------------------------------------------------\n");
    for(i=0;i<NumberOfReactions;i++)
    {
      fprintf(FilePtr,"Reaction %d, current Lambda: %18.10f, maximum Lambda-change: %18.10f\n",i,CFCRXMCLambda[CurrentSystem][i],MaximumReactionLambdaChange[CurrentSystem][i]);
      fprintf(FilePtr,"\tFractional molecules: ");
      for(j=0;j<NumberOfComponents;j++)
      {
        if(ProductsStoichiometry[i][j]>0)
        {
          for(k=0;k<ProductsStoichiometry[i][j];k++)
            fprintf(FilePtr," %d (%s)",Components[j].ProductFractionalMolecules[CurrentSystem][i][k],Components[j].Name);
        }
      }
      fprintf(FilePtr," <-->");
      for(j=0;j<NumberOfComponents;j++)
      {
        if(ReactantsStoichiometry[i][j]>0)
        {
          for(k=0;k<ReactantsStoichiometry[i][j];k++)
            fprintf(FilePtr," %d (%s)",Components[j].ReactantFractionalMolecules[CurrentSystem][i][k],Components[j].Name);
        }
      }
      fprintf(FilePtr,"\n");

      shift=RXMCBiasingFactors[CurrentSystem][i][0];
      fprintf(FilePtr,"\tBiasing Factors: ");
      for(k=0;k<RXMCLambdaHistogramSize;k++)
      {
        fprintf(FilePtr,"%4f ",RXMCBiasingFactors[CurrentSystem][i][k]-shift);
        if((k+1)%10==0&&(k+1)!=RXMCLambdaHistogramSize)  fprintf(FilePtr,"\n\t                 ");
      }
      fprintf(FilePtr,"\n");
    }
    fprintf(FilePtr,"\n");
  }


  number_of_unit_cells=NumberOfUnitCells[CurrentSystem].x*NumberOfUnitCells[CurrentSystem].y*NumberOfUnitCells[CurrentSystem].z;
  if(Framework[CurrentSystem].FrameworkModel==NONE)
  {
    fprintf(FilePtr,"Amount of molecules per component:\n");
    fprintf(FilePtr,"----------------------------------------------------------------------------------------------------------------------------------------------------\n");
    for(i=0;i<NumberOfComponents;i++)
    {
      fprintf(FilePtr,"Component %d (%s), current number of integer/fractional/reaction molecules: %d/%d/%d, density: %9.5lf [kg/m^3]\n",
        i,
        Components[i].Name,
        Components[i].NumberOfMolecules[CurrentSystem]-(Components[i].FractionalMolecule[CurrentSystem]>=0?1:0)-Components[i].NumberOfRXMCMoleculesPresent[CurrentSystem],
        (Components[i].FractionalMolecule[CurrentSystem]>=0?1:0),
        Components[i].NumberOfRXMCMoleculesPresent[CurrentSystem],
        (double)((Components[i].Mass*(REAL)(Components[i].NumberOfMolecules[CurrentSystem]
                       -(Components[i].FractionalMolecule[CurrentSystem]>=0?1:0)
                       -Components[i].NumberOfRXMCMoleculesPresent[CurrentSystem])/
                        Volume[CurrentSystem])*DENSITY_CONVERSION_FACTOR));
      if(Components[i].CFMoleculePresent[CurrentSystem])
      {
        char FractionalMoleculeString[6];
        FractionalMolecule=Components[i].FractionalMolecule[CurrentSystem];
        if (FractionalMolecule<0)
          strcpy(FractionalMoleculeString,"-");
        else
          sprintf(FractionalMoleculeString,"%d", FractionalMolecule);
        fprintf(FilePtr,"\tFractional molecule-id: %s, max. Lambda-change: %5lf (CFMC) %5lf (CB/CFMC)\n",FractionalMoleculeString,MaximumCFLambdaChange[CurrentSystem][i],MaximumCBCFLambdaChange[CurrentSystem][i]);
        fprintf(FilePtr,"\tLambda factors: ");
        for(k=0;k<Components[i].NumberOfAtoms;k++)
        {
          Lambda=0.0;
          if(FractionalMolecule>=0)
          {
            if(Components[i].ExtraFrameworkMolecule)
              Lambda=Cations[CurrentSystem][FractionalMolecule].Atoms[k].CFVDWScalingParameter;
            else
              Lambda=Adsorbates[CurrentSystem][FractionalMolecule].Atoms[k].CFVDWScalingParameter;
          }
          fprintf(FilePtr,"%4f ",Lambda);
          if((k+1)%10==0&&(k+1)!=Components[i].NumberOfAtoms)  fprintf(FilePtr,"\n\t                ");
        }
        fprintf(FilePtr,"\n");
        shift=Components[i].CFBiasingFactors[CurrentSystem][0];
        shift=0.0;
        fprintf(FilePtr,"\tBiasing Factors: ");
        for(k=0;k<Components[i].CFLambdaHistogramSize;k++)
        {
          fprintf(FilePtr,"%4f ",Components[i].CFBiasingFactors[CurrentSystem][k]-shift);
          if((k+1)%10==0&&(k+1)!=Components[i].CFLambdaHistogramSize)  fprintf(FilePtr,"\n\t                 ");
        }
        fprintf(FilePtr,"\n");
      }
    }
    fprintf(FilePtr,"----------------------------------------------------------------------------------------------------------------------------------------------------\n");
  }
  else
  {
    fprintf(FilePtr,"Loadings per component:\n");
    fprintf(FilePtr,"----------------------------------------------------------------------------------------------------------------------------------------------------\n");
    for(i=0;i<NumberOfComponents;i++)
    {
      fprintf(FilePtr,"Component %d (%s), current number of integer/fractional/reaction molecules: %d/%d/%d, density: %9.5lf [kg/m^3]\n",
        i,
        Components[i].Name,
        Components[i].NumberOfMolecules[CurrentSystem]-(Components[i].FractionalMolecule[CurrentSystem]>=0?1:0)-Components[i].NumberOfRXMCMoleculesPresent[CurrentSystem],
        (Components[i].FractionalMolecule[CurrentSystem]>=0?1:0),
        Components[i].NumberOfRXMCMoleculesPresent[CurrentSystem],
        (double)(Components[i].Mass*(REAL)(Components[i].NumberOfMolecules[CurrentSystem]
                                         -(Components[i].FractionalMolecule[CurrentSystem]>=0?1:0)
                                         -Components[i].NumberOfRXMCMoleculesPresent[CurrentSystem])/
                                          Volume[CurrentSystem])*DENSITY_CONVERSION_FACTOR);
      if(Components[i].CFMoleculePresent[CurrentSystem])
      {
        char FractionalMoleculeString[6];
        FractionalMolecule=Components[i].FractionalMolecule[CurrentSystem];
        if (FractionalMolecule<0)
          strcpy(FractionalMoleculeString,"-");
        else
          sprintf(FractionalMoleculeString,"%d", FractionalMolecule);
        fprintf(FilePtr,"\tFractional molecule-id: %s, max. Lambda-change: %5lf (CFMC) %5lf (CB/CFMC)\n",FractionalMoleculeString,MaximumCFLambdaChange[CurrentSystem][i],MaximumCBCFLambdaChange[CurrentSystem][i]);
        fprintf(FilePtr,"\tLambda factors: ");
        for(k=0;k<Components[i].NumberOfAtoms;k++)
        {
          Lambda=0.0;
          if(FractionalMolecule>=0)
          {
            if(Components[i].ExtraFrameworkMolecule)
              Lambda=Cations[CurrentSystem][FractionalMolecule].Atoms[k].CFVDWScalingParameter;
            else
              Lambda=Adsorbates[CurrentSystem][FractionalMolecule].Atoms[k].CFVDWScalingParameter;
          }
          fprintf(FilePtr,"%4f ",Lambda);
          if((k+1)%10==0&&(k+1)!=Components[i].NumberOfAtoms)  fprintf(FilePtr,"\n\t                ");
        }
        fprintf(FilePtr,"\n");
        shift=Components[i].CFBiasingFactors[CurrentSystem][0];
        shift=0.0;
        fprintf(FilePtr,"\tBiasing Factors: ");
        for(k=0;k<Components[i].CFLambdaHistogramSize;k++)
        {
          fprintf(FilePtr,"%4f ",Components[i].CFBiasingFactors[CurrentSystem][k]-shift);
          if((k+1)%10==0&&(k+1)!=Components[i].CFLambdaHistogramSize)  fprintf(FilePtr,"\n\t                 ");
        }
        fprintf(FilePtr,"\n");
      }

      loading=(REAL)(Components[i].NumberOfMolecules[CurrentSystem]
                   -(Components[i].FractionalMolecule[CurrentSystem]>=0?1:0)
                   -Components[i].NumberOfRXMCMoleculesPresent[CurrentSystem])/(REAL)number_of_unit_cells;
      average_loading=GetAverageWeightedComponentProperty(NumberOfIntegerMoleculesPerComponentAccumulated,i)/(REAL)number_of_unit_cells;
      fprintf(FilePtr,"\tabsolute adsorption: %9.5lf [mol/uc], %14.4lf [mol/kg],      %14.4lf [mg/g]\n",
         (double)loading,
         (double)(Components[i].MOLEC_PER_UC_TO_MOL_PER_KG[CurrentSystem]*loading),
         (double)(Components[i].MOLEC_PER_UC_TO_MILLIGRAM_PER_GRAM_OF_FRAMEWORK[CurrentSystem]*loading));
      fprintf(FilePtr,"\t                                          %14.4lf [cm^3 STP/g],  %14.4lf [cm^3 STP/cm^3]\n",
         (double)(Components[i].MOLEC_PER_UC_TO_CC_STP_G[CurrentSystem]*loading),
         (double)(Components[i].MOLEC_PER_UC_TO_CC_STP_CC[CurrentSystem]*loading));

      loading=((REAL)Components[i].NumberOfMolecules[CurrentSystem]
                    -Components[i].AmountOfExcessMolecules[CurrentSystem]
                    -(Components[i].FractionalMolecule[CurrentSystem]>=0?1:0)
                    -Components[i].NumberOfRXMCMoleculesPresent[CurrentSystem])/(REAL)number_of_unit_cells;
      average_loading=GetAverageWeightedComponentProperty(NumberOfExcessMoleculesPerComponentAccumulated,i)/(REAL)number_of_unit_cells;
      fprintf(FilePtr,"\texcess adsorption:   %9.5lf [mol/uc], %14.4lf [mol/kg],      %14.4lf [mg/g]\n",
         (double)loading,
         (double)(Components[i].MOLEC_PER_UC_TO_MOL_PER_KG[CurrentSystem]*loading),
         (double)(Components[i].MOLEC_PER_UC_TO_MILLIGRAM_PER_GRAM_OF_FRAMEWORK[CurrentSystem]*loading));
      fprintf(FilePtr,"\t                                          %14.4lf [cm^3 STP/g],  %14.4lf [cm^3 STP/cm^3]\n",
         (double)(Components[i].MOLEC_PER_UC_TO_CC_STP_G[CurrentSystem]*loading),
         (double)(Components[i].MOLEC_PER_UC_TO_CC_STP_CC[CurrentSystem]*loading));

      fprintf(FilePtr,"----------------------------------------------------------------------------------------------------------------------------------------------------\n");
    }
  }
  fprintf(FilePtr,"Degrees of freedom: %d %d %d %d\n",DegreesOfFreedom[CurrentSystem],DegreesOfFreedomFramework[CurrentSystem],
          DegreesOfFreedomAdsorbates[CurrentSystem],DegreesOfFreedomCations[CurrentSystem]);
  fprintf(FilePtr,"Number of Framework-atoms: %6d\n",Framework[CurrentSystem].TotalNumberOfAtoms);
  fprintf(FilePtr,"Number of Adsorbates:      %6d (%d integer, %d fractional, %d reaction)\n",
          NumberOfAdsorbateMolecules[CurrentSystem],
          TotalNumberOfIntegerAdsorbates(),
          TotalNumberOfFractionalAdsorbates(),
          NumberOfReactionAdsorbateMolecules[CurrentSystem]);
  fprintf(FilePtr,"Number of Cations:         %6d (%d integer, %d fractional, %d reaction)\n",
          NumberOfCationMolecules[CurrentSystem],
          TotalNumberOfIntegerCations(),
          TotalNumberOfFractionalCations(),
          NumberOfReactionCationMolecules[CurrentSystem]);
  fprintf(FilePtr,"\n");

  if(SimulationType==MOLECULAR_DYNAMICS)
  {
    fprintf(FilePtr,"Conserved energy: % 22.10lf Energy drifts: % 11.10lf % 22.10lf\n",
        (double)(ConservedEnergy[CurrentSystem]*ENERGY_TO_KELVIN),
        (double)(fabs((ConservedEnergy[CurrentSystem]-ReferenceEnergy[CurrentSystem])/ReferenceEnergy[CurrentSystem])*ENERGY_TO_KELVIN),
        (double)(Drift[CurrentSystem]/(CurrentCycle+1.0)));

    if(DegreesOfFreedom[CurrentSystem]>0)
      fprintf(FilePtr,"Temperature:            % 22.10lf\n",
        (double)(2.0*UKinetic[CurrentSystem]/(K_B*DegreesOfFreedom[CurrentSystem])));
    if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
      fprintf(FilePtr,"Temperature Framework:  % 22.10lf\n",
        (double)(2.0*UHostKinetic[CurrentSystem]/(K_B*DegreesOfFreedomFramework[CurrentSystem])));
    if(DegreesOfFreedomAdsorbates[CurrentSystem]>0)
      fprintf(FilePtr,"Temperature Adsorbates: % 22.10lf\n",
        (double)(2.0*UAdsorbateKinetic[CurrentSystem]/(K_B*DegreesOfFreedomAdsorbates[CurrentSystem])));
    if(DegreesOfFreedomCations[CurrentSystem]>0)
      fprintf(FilePtr,"Temperature Cations:    % 22.10lf\n",
        (double)(2.0*UCationKinetic[CurrentSystem]/(K_B*DegreesOfFreedomCations[CurrentSystem])));
    fprintf(FilePtr,"Cell temperature: % 8.3lf\n",GetCellTemperature());
    fprintf(FilePtr,"Current total kinetic energy:       % 22.10lf [K]\n",(double)UKinetic[CurrentSystem]);
    fprintf(FilePtr,"Current total Nose-Hoover energy:   % 22.10lf [K]\n",(double)UNoseHoover[CurrentSystem]);
  }
  fprintf(FilePtr,"Current total potential energy:       % 22.10lf [K]\n",
      (double)(UTotal[CurrentSystem]*ENERGY_TO_KELVIN));

  fprintf(FilePtr,"\tCurrent Host-Host energy:           % 22.10lf [K]\n",(double)(UHostHost[CurrentSystem]*ENERGY_TO_KELVIN));
  fprintf(FilePtr,"\tCurrent Host-Adsorbate energy:      % 22.10lf [K]\n",(double)(UHostAdsorbate[CurrentSystem]*ENERGY_TO_KELVIN));
  fprintf(FilePtr,"\tCurrent Host-Cation energy:         % 22.10lf [K]\n",(double)(UHostCation[CurrentSystem]*ENERGY_TO_KELVIN));
  fprintf(FilePtr,"\tCurrent Adsorbate-Adsorbate energy: % 22.10lf [K]\n",(double)(UAdsorbateAdsorbate[CurrentSystem]*ENERGY_TO_KELVIN));
  fprintf(FilePtr,"\tCurrent Cation-Cation energy:       % 22.10lf [K]\n",(double)(UCationCation[CurrentSystem]*ENERGY_TO_KELVIN));
  fprintf(FilePtr,"\tCurrent Adsorbate-Cation energy:    % 22.10lf [K]\n",(double)(UAdsorbateCation[CurrentSystem]*ENERGY_TO_KELVIN));

  if(ComputePolarization)
  {
    fprintf(FilePtr,"\tCurrent polarization energy:        % 22.10lf [K]\n",(double)(UHostPolarization[CurrentSystem]+
         UAdsorbatePolarization[CurrentSystem]+UCationPolarization[CurrentSystem])*ENERGY_TO_KELVIN*ENERGY_TO_KELVIN);
    fprintf(FilePtr,"\tCurrent back-polarization energy:   % 22.10lf [K]\n",(double)(UHostBackPolarization[CurrentSystem]+
         UAdsorbateBackPolarization[CurrentSystem]+UCationBackPolarization[CurrentSystem])*ENERGY_TO_KELVIN*ENERGY_TO_KELVIN);
  }

  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
  {
    if(Framework[CurrentSystem].NumberOfBondsDefinitions>0)
      fprintf(FilePtr,"\t\tCurrent Host-Bond energy:             % 22.10lf [K]\n",
          (double)UHostBond[CurrentSystem]*ENERGY_TO_KELVIN);
    if(Framework[CurrentSystem].NumberOfUreyBradleyDefinitions>0)
      fprintf(FilePtr,"\t\tCurrent Host-UreyBradly energy:       % 22.10lf [K]\n",
          (double)UHostUreyBradley[CurrentSystem]*ENERGY_TO_KELVIN);
    if(Framework[CurrentSystem].NumberOfBendDefinitions>0)
      fprintf(FilePtr,"\t\tCurrent Host-Bend energy:             % 22.10lf [K]\n",
          (double)UHostBend[CurrentSystem]*ENERGY_TO_KELVIN);
    if(Framework[CurrentSystem].NumberOfInversionBendDefinitions>0)
      fprintf(FilePtr,"\t\tCurrent Host-Inversion Bend energy:   % 22.10lf [K]\n",
          (double)UHostInversionBend[CurrentSystem]*ENERGY_TO_KELVIN);
    if(Framework[CurrentSystem].NumberOfTorsionDefinitions>0)
      fprintf(FilePtr,"\t\tCurrent Host-Torsion energy:          % 22.10lf [K]\n",
          (double)UHostTorsion[CurrentSystem]*ENERGY_TO_KELVIN);
    if(Framework[CurrentSystem].NumberOfImproperTorsionDefinitions>0)
      fprintf(FilePtr,"\t\tCurrent Host-Improper torsion energy: % 22.10lf [K]\n",
          (double)UHostImproperTorsion[CurrentSystem]*ENERGY_TO_KELVIN);
    if(Framework[CurrentSystem].NumberOfBondBondDefinitions>0)
      fprintf(FilePtr,"\t\tCurrent Host-Bond/Bond energy:        % 22.10lf [K]\n",
          (double)UHostBondBond[CurrentSystem]*ENERGY_TO_KELVIN);
    if(Framework[CurrentSystem].NumberOfBendBendDefinitions>0)
      fprintf(FilePtr,"\t\tCurrent Host-Bend/Bend energy:        % 22.10lf [K]\n",
          (double)UHostBendBend[CurrentSystem]*ENERGY_TO_KELVIN);
    if(Framework[CurrentSystem].NumberOfBondBendDefinitions>0)
      fprintf(FilePtr,"\t\tCurrent Host-Bond/Bend energy:        % 22.10lf [K]\n",
          (double)UHostBondBend[CurrentSystem]*ENERGY_TO_KELVIN);
    if(Framework[CurrentSystem].NumberOfBondTorsionDefinitions>0)
      fprintf(FilePtr,"\t\tCurrent Host-Bond/Torsion energy:     % 22.10lf [K]\n",
          (double)UHostBondTorsion[CurrentSystem]*ENERGY_TO_KELVIN);
    if(Framework[CurrentSystem].NumberOfBendTorsionDefinitions>0)
      fprintf(FilePtr,"\t\tCurrent Host-Bend/Torsion energy:     % 22.10lf [K]\n",
          (double)UHostBendTorsion[CurrentSystem]*ENERGY_TO_KELVIN);
  }

  fprintf(FilePtr,"\n");
  PrintWarningStatus();
  fprintf(FilePtr,"\n\n");
  fflush(FilePtr);
}

REAL GetAverageAdsorbatesTemperature(void)
{
  int i;
  REAL sum1,sum2;

  sum1=sum2=0.0;
  for(i=0;i<NR_BLOCKS;i++)
    if(BlockWeightedCount[CurrentSystem][i]>0.0)
    {
      sum1+=TemperatureAdsorbatesAccumulated[CurrentSystem][i];
      sum2+=BlockWeightedCount[CurrentSystem][i];
    }
  if((sum1>0.0)&&(sum2>0.0))
    return sum1/sum2;
  else
    return 0.0;
}

REAL GetAverageCationsTemperature(void)
{
  int i;
  REAL sum1,sum2;

  sum1=sum2=0.0;
  for(i=0;i<NR_BLOCKS;i++)
    if(BlockWeightedCount[CurrentSystem][i]>0.0)
    {
      sum1+=TemperatureCationsAccumulated[CurrentSystem][i];
      sum2+=BlockWeightedCount[CurrentSystem][i];
    }
  if((sum1>0.0)&&(sum2>0.0))
    return sum1/sum2;
  else
    return 0.0;
}

REAL GetAverageFrameworkTemperature(void)
{
  int i;
  REAL sum1,sum2;

  sum1=sum2=0.0;
  for(i=0;i<NR_BLOCKS;i++)
    if(BlockWeightedCount[CurrentSystem][i]>0.0)
    {
      sum1+=TemperatureFrameworkAccumulated[CurrentSystem][i];
      sum2+=BlockWeightedCount[CurrentSystem][i];
    }
  if((sum1>0.0)&&(sum2>0.0))
    return sum1/sum2;
  else
    return 0.0;
}

REAL GetAverageMolecularPressure(void)
{
  int i;
  REAL sum1,sum2;

  sum1=sum2=0.0;
  for(i=0;i<NR_BLOCKS;i++)
    if(BlockWeightedCount[CurrentSystem][i]>0.0)
    {
      sum1+=MolecularPressureAccumulated[CurrentSystem][i];
      sum2+=BlockWeightedCount[CurrentSystem][i];
    }
  if(sum2>0.0)
    return sum1/sum2;
  else
    return 0.0;
}

REAL GetAverageCompressibility(void)
{
  int i;
  REAL sum1,sum2;

  sum1=sum2=0.0;
  for(i=0;i<NR_BLOCKS;i++)
    if(BlockWeightedCount[CurrentSystem][i]>0.0)
    {
      sum1+=CompressibilityAccumulated[CurrentSystem][i];
      sum2+=BlockWeightedCount[CurrentSystem][i];
    }
  if(sum2>0.0)
    return sum1/sum2;
  else
    return 0.0;
}

REAL GetAveragePressure(void)
{
  int i;
  REAL sum1,sum2;

  sum1=sum2=0.0;
  for(i=0;i<NR_BLOCKS;i++)
    if(BlockWeightedCount[CurrentSystem][i]>0.0)
    {
      sum1+=MolecularPressureAccumulated[CurrentSystem][i];
      sum2+=BlockWeightedCount[CurrentSystem][i];
    }
  if(sum2>0.0)
    return sum1/sum2;
  else
    return 0.0;
}

REAL_MATRIX3x3 GetAverageMolecularStressTensor(void)
{
  int i;
  REAL sum;
  REAL_MATRIX3x3 Stress;

  sum=0.0;
  Stress.ax=Stress.bx=Stress.cx=0.0;
  Stress.ay=Stress.by=Stress.cy=0.0;
  Stress.az=Stress.bz=Stress.cz=0.0;
  for(i=0;i<NR_BLOCKS;i++)
    if(BlockWeightedCount[CurrentSystem][i]>0.0)
    {
      Stress.ax+=MolecularStressTensorAccumulated[CurrentSystem][i].ax;
      Stress.ay+=MolecularStressTensorAccumulated[CurrentSystem][i].ay;
      Stress.az+=MolecularStressTensorAccumulated[CurrentSystem][i].az;
      Stress.bx+=MolecularStressTensorAccumulated[CurrentSystem][i].bx;
      Stress.by+=MolecularStressTensorAccumulated[CurrentSystem][i].by;
      Stress.bz+=MolecularStressTensorAccumulated[CurrentSystem][i].bz;
      Stress.cx+=MolecularStressTensorAccumulated[CurrentSystem][i].cx;
      Stress.cy+=MolecularStressTensorAccumulated[CurrentSystem][i].cy;
      Stress.cz+=MolecularStressTensorAccumulated[CurrentSystem][i].cz;
      sum+=BlockWeightedCount[CurrentSystem][i];
    }
  if(sum>0.0)
  {
    Stress.ax/=sum; Stress.bx/=sum; Stress.cx/=sum;
    Stress.ay/=sum; Stress.by/=sum; Stress.cy/=sum;
    Stress.az/=sum; Stress.bz/=sum; Stress.cz/=sum;
  }
  return Stress;
}

REAL GetAverageIdealGasPartPressure(void)
{
  int i;
  REAL count,sum;

  count=sum=0.0;
  for(i=0;i<NR_BLOCKS;i++)
  {
    if(BlockWeightedCount[CurrentSystem][i]>0.0)
    {
      sum+=PressureIdealGasPartAccumulated[CurrentSystem][i];
      count+=BlockWeightedCount[CurrentSystem][i];
    }
  }
  return sum/count;
}

REAL GetAverageExcessPartPressure(void)
{
  int i;
  REAL count,sum;

  count=sum=0.0;
  for(i=0;i<NR_BLOCKS;i++)
  {
    if(BlockWeightedCount[CurrentSystem][i]>0.0)
    {
      sum+=PressureExcessPartAccumulated[CurrentSystem][i];
      count+=BlockWeightedCount[CurrentSystem][i];
    }
  }
  return sum/count;
}


REAL GetAverageTailCorrectionPressure(void)
{
  int i;
  REAL count,sum;

  count=sum=0.0;
  for(i=0;i<NR_BLOCKS;i++)
  {
    if(BlockWeightedCount[CurrentSystem][i]>0.0)
    {
      sum+=PressureTailCorrectionAccumulated[CurrentSystem][i];
      count+=BlockWeightedCount[CurrentSystem][i];
    }
  }
  return sum/count;
}


REAL GetAverageVolumeSquared(void)
{
  int i;
  REAL count,sum;

  count=sum=0.0;
  for(i=0;i<NR_BLOCKS;i++)
  {
    if(BlockWeightedCount[CurrentSystem][i]>0.0)
    {
      sum+=VolumeSquaredAccumulated[CurrentSystem][i];
      count+=BlockWeightedCount[CurrentSystem][i];
    }
  }
  return sum/count;
}


REAL GetAverageIsothermalExpansionCoefficient(void)
{
  int i;
  REAL HV,V,H,T;
  REAL count,sum;

  count=sum=0.0;
  HV=V=H=0.0;
  for(i=0;i<NR_BLOCKS;i++)
  {
    if(BlockWeightedCount[CurrentSystem][i]>0.0)
    {
      HV+=EnthalpyTimesVolumeAccumulated[CurrentSystem][i];
      V+=VolumeAccumulated[CurrentSystem][i];
      H+=EnthalpyAccumulated[CurrentSystem][i];
      count+=BlockWeightedCount[CurrentSystem][i];
    }
  }
  if(count>0.0)
  {
    HV/=count;
    H/=count;
    V/=count;
    T=therm_baro_stats.ExternalTemperature[CurrentSystem];
    return 1e6*VOLUMETRIC_EXPANSION_COEFFICIENT_CONVERSION_FACTOR*(HV-H*V)/(K_B*V*SQR(T));
  }
  else return 0.0;
}

REAL GetAverageIsothermalCompressibilityCoefficient(void)
{
  int i;
  REAL T,V,V2;
  REAL count,sum;

  sum=count=0.0;
  V=V2=0.0;
  for(i=0;i<NR_BLOCKS;i++)
  {
    if(BlockWeightedCount[CurrentSystem][i]>0.0)
    {
      V+=VolumeAccumulated[CurrentSystem][i];
      V2+=VolumeSquaredAccumulated[CurrentSystem][i];
      count+=BlockWeightedCount[CurrentSystem][i];
    }
  }
  if(count>0.0)
  {
    V/=count;
    V2/=count;
    T=therm_baro_stats.ExternalTemperature[CurrentSystem];
    return 1e12*ISOTHERMAL_COMPRESSIBILITY_CONVERSION_FACTOR*(V2-SQR(V))/(K_B*T*V);
  }
  else return 0.0;
}

REAL GetAverageHeatCapacityConstantPressure(void)
{
  int i;
  REAL N,T;
  REAL H,VH,V,U,UH,p;
  REAL count,sum;

  N=NumberOfAdsorbateMolecules[CurrentSystem]+NumberOfCationMolecules[CurrentSystem];
  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
    N+=Framework[CurrentSystem].TotalNumberOfAtoms;
  count=sum=0.0;
  H=VH=V=U=UH=0.0;
  for(i=0;i<NR_BLOCKS;i++)
  {
    if(BlockWeightedCount[CurrentSystem][i]>0.0)
    {
      H+=EnthalpyAccumulated[CurrentSystem][i];
      VH+=EnthalpyTimesVolumeAccumulated[CurrentSystem][i];
      V+=VolumeAccumulated[CurrentSystem][i];
      U+=TotalEnergyAccumulated[CurrentSystem][i];
      UH+=EnthalpyTimesEnergyAccumulated[CurrentSystem][i];
      count+=BlockWeightedCount[CurrentSystem][i];
    }
  }
  if(count>0.0)
  {
    p=therm_baro_stats.ExternalPressure[CurrentSystem][0];
    H/=count;
    VH/=count;
    V/=count;
    U/=count;
    UH/=count;
    T=therm_baro_stats.ExternalTemperature[CurrentSystem];
    return HEAT_CAPACITY_CONVERSION_FACTOR*(((UH-U*H)+p*(VH-V*H))/(N*K_B*SQR(T))+(5.0/2.0)*K_B);
  }
  else return 0.0;
}


REAL GetAverageHeatCapacity(void)
{
  int i;
  REAL N,H2,H,T;
  REAL count,sum;
  REAL sum1,sum2;
  int NumberOfBlocks;

  N=NumberOfAdsorbateMolecules[CurrentSystem]+NumberOfCationMolecules[CurrentSystem];
  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
    N+=Framework[CurrentSystem].TotalNumberOfAtoms;
  count=sum=0.0;
  sum1=sum2=0.0;
  H2=H=0.0;

  T=therm_baro_stats.ExternalTemperature[CurrentSystem];

  NumberOfBlocks=0;
  for(i=0;i<NR_BLOCKS;i++)
  {
    if(BlockWeightedCount[CurrentSystem][i]>0.0)
    {
      H2=EnthalpySquaredAccumulated[CurrentSystem][i]/BlockWeightedCount[CurrentSystem][i];
      H=EnthalpyAccumulated[CurrentSystem][i]/BlockWeightedCount[CurrentSystem][i];
      sum+=HEAT_CAPACITY_CONVERSION_FACTOR*(H2-SQR(H));
      NumberOfBlocks++;
    }
  }
  return sum/NumberOfBlocks/(N*K_B*SQR(T))+HEAT_CAPACITY_CONVERSION_FACTOR*(5.0/2.0)*K_B;
}

/*
REAL GetAverageHeatCapacityConstantPressure(void)
{
  int i;
  REAL N,H2,H,T;
  REAL count,sum;

  N=NumberOfAdsorbateMolecules[CurrentSystem]+NumberOfCationMolecules[CurrentSystem];
  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
    N+=Framework[CurrentSystem].TotalNumberOfAtoms;
  count=sum=0.0;
  H2=H=0.0;
  for(i=0;i<NR_BLOCKS;i++)
  {
    if(BlockWeightedCount[CurrentSystem][i]>0.0)
    {
      H2+=EnthalpySquaredAccumulated[CurrentSystem][i];
      H+=EnthalpyAccumulated[CurrentSystem][i];
      count+=BlockWeightedCount[CurrentSystem][i];
    }
  }
  if(count>0.0)
  {
    H2/=count;
    H/=count;
    T=therm_baro_stats.ExternalTemperature[CurrentSystem];
    return HEAT_CAPACITY_CONVERSION_FACTOR*((H2-SQR(H))/(N*K_B*SQR(T))+(3.0/2.0)*K_B);
  }
  else return 0.0;
}
*/

// return the average Henry coefficient for Component 'comp' in mol/kg/Pa.
REAL GetAverageHenryCoefficient(int comp)
{
  int i;
  REAL FrameworkDensity,Temperature;
  REAL count,sum;

  // get framework density in kg/m^3
  FrameworkDensity=1e-3*Framework[CurrentSystem].FrameworkDensity;
  Temperature=therm_baro_stats.ExternalTemperature[CurrentSystem];
  count=sum=0.0;
  for(i=0;i<NR_BLOCKS;i++)
  {
    if(WidomRosenbluthFactorCount[CurrentSystem][comp][i]>0.0)
    {
      sum+=WidomRosenbluthFactorAccumulated[CurrentSystem][comp][i]/
            Components[comp].IdealGasRosenbluthWeight[CurrentSystem];
      count+=WidomRosenbluthFactorCount[CurrentSystem][comp][i];
    }
  }
  if(count>0.0)
    return (1.0/(MOLAR_GAS_CONSTANT*Temperature*FrameworkDensity))*(sum/count);
  else return 0.0;
}

REAL GetAverageRosenbluthWeight(int comp)
{
  int i;
  REAL count,sum;

  count=sum=0.0;
  for(i=0;i<NR_BLOCKS;i++)
  {
    if(WidomRosenbluthFactorCount[CurrentSystem][comp][i]>0.0)
    {
      sum+=WidomRosenbluthFactorAccumulated[CurrentSystem][comp][i];
      count+=WidomRosenbluthFactorCount[CurrentSystem][comp][i];
    }
  }
  if(count>0.0)
    return sum/count;
  else return 0.0;
}

REAL GetAverageGibbsWidom(int comp)
{
  int i;
  REAL count,sum;

  count=sum=0.0;
  for(i=0;i<NR_BLOCKS;i++)
  {
    if(GibbsWidomRosenbluthFactorCount[CurrentSystem][comp][i]>0.0)
    {
      sum+=GibbsWidomRosenbluthFactorAccumulated[CurrentSystem][comp][i];
      count+=GibbsWidomRosenbluthFactorCount[CurrentSystem][comp][i];
    }
  }
  if(count>0.0)
    return sum/count;
  else return 0.0;
}

REAL GetAverageWidomExcess(int comp)
{
  int i;
  REAL count,sum;

  count=sum=0.0;
  for(i=0;i<NR_BLOCKS;i++)
  {
    if(WidomRosenbluthFactorCount[CurrentSystem][comp][i]>0.0)
    {
      sum+=WidomRosenbluthFactorAccumulated[CurrentSystem][comp][i];
      count+=WidomRosenbluthFactorCount[CurrentSystem][comp][i];
    }
  }
  if(count>0.0)
    return sum/count;
  else return 0.0;
}


REAL GetAverageGibbsWidomIdealGas(int comp)
{
  int i;
  REAL count,sum;

  count=sum=0.0;
  for(i=0;i<NR_BLOCKS;i++)
  {
    if(GibbsWidomRosenbluthFactorCount[CurrentSystem][comp][i]>0.0)
    {
      sum+=GibbsWidomIdealGasAccumulated[CurrentSystem][comp][i];
      count+=GibbsWidomRosenbluthFactorCount[CurrentSystem][comp][i];
    }
  }
  if(count>0.0)
    return sum/count;
  else return 0.0;
}

REAL GetAverageWidomIdealGas(int comp)
{
  int i;
  REAL count,sum;

  count=sum=0.0;
  for(i=0;i<NR_BLOCKS;i++)
  {
    if(WidomRosenbluthFactorCount[CurrentSystem][comp][i]>0.0)
    {
      sum+=WidomIdealGasAccumulated[CurrentSystem][comp][i];
      count+=WidomRosenbluthFactorCount[CurrentSystem][comp][i];
    }
  }
  if(count>0.0)
    return sum/count;
  else return 0.0;
}

REAL GetAverageGibbsInverseDensity()
{
  int i;
  REAL count,sum;

  count=sum=0.0;
  for(i=0;i<NR_BLOCKS;i++)
  {
    sum+=GibbsInverseDensityAccumulated[CurrentSystem][i];
    count+=BlockWeightedCount[CurrentSystem][i];
  }
  if(count>0.0)
    return sum/count;
  else return 0.0;
}

REAL GetAverageGibbsInverseDensityForComponent(int comp)
{
  int i;
  REAL count,sum;

  count=sum=0.0;
  for(i=0;i<NR_BLOCKS;i++)
  {
    sum+=GibbsInverseDensityPerComponentAccumulated[CurrentSystem][comp][i];
    count+=BlockWeightedCount[CurrentSystem][i];
  }
  if(count>0.0)
    return sum/count;
  else return 0.0;
}


REAL GetAverageInverseDensity()
{
  int i;
  REAL count,sum;

  count=sum=0.0;
  for(i=0;i<NR_BLOCKS;i++)
  {
    sum+=InverseDensityAccumulated[CurrentSystem][i];
    count+=BlockWeightedCount[CurrentSystem][i];
  }
  if(count>0.0)
    return sum/count;
  else return 0.0;
}

REAL GetWidomHeat(int comp)
{
  int i;
  REAL count,sum;

  count=sum=0.0;
  for(i=0;i<NR_BLOCKS;i++)
  {
    sum+=WidomEnergyDifferenceAccumulated[CurrentSystem][comp][i];
    count+=WidomRosenbluthFactorAccumulated[CurrentSystem][comp][i];
  }
  if(count>0.0)
    return sum/count;
  else return 0.0;
}

REAL GetFrameworkHeat(int comp)
{
  int i;
  REAL count,sum;

  count=sum=0.0;
  for(i=0;i<NR_BLOCKS;i++)
  {
    sum+=WidomEnergyFrameworkAccumulated[CurrentSystem][comp][i];
    count+=WidomEnergyFrameworkCount[CurrentSystem][comp][i];
  }
  if(count>0.0)
    return sum/count;
  else return 0.0;
}

REAL GetAverageFrameworkSurfaceArea(void)
{
  int i;
  REAL count,sum;

  // get surface area in Angstrom^2
  count=sum=0.0;
  for(i=0;i<NR_BLOCKS;i++)
  {
    if(SurfaceAreaCount[CurrentSystem][i]>0.0)
    {
      sum+=SurfaceAreaFrameworkAccumulated[CurrentSystem][i];
      count+=SurfaceAreaCount[CurrentSystem][i];
    }
  }
  if(count>0.0)
    return sum/count;
  else return 0.0;
}

REAL GetAverageFrameworksSurfaceArea(int f)
{
  int i;
  REAL count,sum;

  // get surface area in Angstrom^2
  count=sum=0.0;
  for(i=0;i<NR_BLOCKS;i++)
  {
    if(SurfaceAreaCount[CurrentSystem][i]>0.0)
    {
      sum+=SurfaceAreaFrameworksAccumulated[CurrentSystem][f][i];
      count+=SurfaceAreaCount[CurrentSystem][i];
    }
  }
  if(count>0.0)
    return sum/count;
  else return 0.0;
}

REAL GetAverageCationSurfaceArea(void)
{
  int i;
  REAL count,sum;

  // get surface area in Angstrom^2
  count=sum=0.0;
  for(i=0;i<NR_BLOCKS;i++)
  {
    if(SurfaceAreaCount[CurrentSystem][i]>0.0)
    {
      sum+=SurfaceAreaCationsAccumulated[CurrentSystem][i];
      count+=SurfaceAreaCount[CurrentSystem][i];
    }
  }
  if(count>0.0)
    return sum/count;
  else return 0.0;
}

void PrintPropertyStatus(long long CurrentCycle,long long NumberOfCycles, FILE *FilePtr)
{
  int i,k;
  REAL_MATRIX6x6 ElasticConstants;
  REAL_MATRIX6x6 BornTermVoigt,StressFluctuationTermVoigt,IdealGasTermVoigt;
  VECTOR InertiaAverage;
  REAL count;

  fprintf(FilePtr,"Average Properties at Current cycle: %lld out of %lld\n",CurrentCycle,NumberOfCycles);
  fprintf(FilePtr,"========================================================================================\n\n");
  fprintf(FilePtr,"Framework surface area: %18.10lf [m^2/g]  %18.10lf [m^2/cm^3] %18.10lf [A^2]\n",
        (double)(GetAverageFrameworkSurfaceArea()*SQR(ANGSTROM)*AVOGADRO_CONSTANT/(Framework[CurrentSystem].FrameworkMass)),
        (double)(1.0e4*GetAverageFrameworkSurfaceArea()/Volume[CurrentSystem]),(double)GetAverageFrameworkSurfaceArea());
  for(i=0;i<Framework[CurrentSystem].NumberOfFrameworks;i++)
    fprintf(FilePtr,"\tFramework %d individual surface area: %18.10lf [m^2/g]  %18.10lf [m^2/cm^3] %18.10f [A^2]\n",i,
        (double)(GetAverageFrameworksSurfaceArea(i)*SQR(ANGSTROM)*AVOGADRO_CONSTANT/(Framework[CurrentSystem].FrameworkMass)),
        (double)(1.0e4*GetAverageFrameworksSurfaceArea(i)/Volume[CurrentSystem]),(double)GetAverageFrameworksSurfaceArea(i));
  fprintf(FilePtr,"\tCation surface area:                 %18.10lf [m^2/g]  %18.10lf [m^2/cm^3] %18.10f [A^2]\n",
        (double)(GetAverageCationSurfaceArea()*SQR(ANGSTROM)*AVOGADRO_CONSTANT/(Framework[CurrentSystem].FrameworkMass)),
        (double)(1.0e4*GetAverageCationSurfaceArea()/Volume[CurrentSystem]),(double)GetAverageCationSurfaceArea());

  fprintf(FilePtr,"Compressibility: %18.10lf [-]\n",(double)GetAverageWeightedProperty(CompressibilityAccumulated));

  fprintf(FilePtr,"Henry coefficients\n");
  for(i=0;i<NumberOfComponents;i++)
    fprintf(FilePtr,"\tComponent %d: %lg [mol/kg/Pa] (Rosenbluth factor new: %lg [-])\n",i,
        (double)GetAverageHenryCoefficient(i),(double)GetAverageRosenbluthWeight(i));

  fprintf(FilePtr,"Energy <U_gh>_1-<U_h>_0 from Widom\n");
  for(i=0;i<NumberOfComponents;i++)
  {
    if(Components[i].Widom)
    {
      fprintf(FilePtr,"\tComponent %d: %18.10lf [K]  (%18.10lf [kJ/mol])\n",i,
        (double)((GetWidomHeat(i)-GetFrameworkHeat(i))*ENERGY_TO_KELVIN),
        (double)((GetWidomHeat(i)-GetFrameworkHeat(i))*ENERGY_TO_KELVIN*KELVIN_TO_KJ_PER_MOL));
    }
    if(Components[i].FractionOfWidomMove>0.0)
    {
      fprintf(FilePtr,"Component [%s] average Widom: %18.10f, average chemical potential: %18.10f [K]\n",
        Components[i].Name,GetAverageWidomExcess(i),
        ((-log(GetAverageWidomExcess(i))-log(GetAverageWidomIdealGas(i)))/Beta[CurrentSystem])*ENERGY_TO_KELVIN);
      fprintf(FilePtr,"\t\t(average excess chemical potential: %18.10f [K], ideal-gas contribution: %18.10f [K])\n",
            (-log(GetAverageWidomExcess(i))/Beta[CurrentSystem])*ENERGY_TO_KELVIN,
            (-log(GetAverageWidomIdealGas(i))/Beta[CurrentSystem])*ENERGY_TO_KELVIN);
            
      /*fprintf(FilePtr,"ideal-gas contribution: %18.10f [K])\n ideal-gas value: %18.10f [K])\n volume: %18.10f [A^3])\n inverse density: %18.10f [])\n Beta: %18.10f [])\n toto: %18.10f [])\n",
            (-log(GetAverageWidomIdealGas(i))/Beta[CurrentSystem])*ENERGY_TO_KELVIN,
            (-log(GetAverageInverseDensity())/Beta[CurrentSystem])*ENERGY_TO_KELVIN,
            Volume[CurrentSystem],
            GetAverageInverseDensity(),
            Beta[CurrentSystem],
            -log(GetAverageInverseDensity())/Beta[CurrentSystem],
            );
       */
    }

    if(Components[i].FractionOfGibbsWidomMove>0.0)
    {
      fprintf(FilePtr,"Component [%s] average Gibbs Widom: %18.10f, average chemical potential: %18.10f [K]\n",
        Components[i].Name,GetAverageGibbsWidom(i),
        (-log(GetAverageGibbsWidom(i))/Beta[CurrentSystem])*ENERGY_TO_KELVIN);
      fprintf(FilePtr,"\t\t(average excess chemical potential: %18.10f [K], ideal-gas contribution: %18.10f [K])\n",
             (-log(GetAverageGibbsWidom(i))/Beta[CurrentSystem])*ENERGY_TO_KELVIN
            -(-log(GetAverageGibbsWidomIdealGas(i))/Beta[CurrentSystem])*ENERGY_TO_KELVIN,
            (-log(GetAverageGibbsWidomIdealGas(i))/Beta[CurrentSystem])*ENERGY_TO_KELVIN);
    }
  }
  fprintf(FilePtr,"\n");

  // output intermediate averages of the elastic constants
  if(ComputeBornTerm)
  {
    fprintf(FilePtr,"Elastic constant [GPa]\n");
    fprintf(FilePtr,"----------------------\n");
    ElasticConstants=ComputeCurrentAverageElasticConstants(&BornTermVoigt,&StressFluctuationTermVoigt,&IdealGasTermVoigt);
    PrintRealMatrix6x6ToFile(&ElasticConstants,FilePtr,1e9/PRESSURE_CONVERSION_FACTOR);
    fprintf(FilePtr,"\n");

    fprintf(FilePtr,"Born-term [GPa]\n");
    fprintf(FilePtr,"---------------\n");
    PrintRealMatrix6x6ToFile(&BornTermVoigt,FilePtr,1e9/PRESSURE_CONVERSION_FACTOR);
    fprintf(FilePtr,"\n");

    fprintf(FilePtr,"Stress-fluctuation-term [GPa]\n");
    fprintf(FilePtr,"-----------------------------\n");
    PrintRealMatrix6x6ToFile(&StressFluctuationTermVoigt,FilePtr,1e9/PRESSURE_CONVERSION_FACTOR);
    fprintf(FilePtr,"\n");

    fprintf(FilePtr,"Ideal-gas-term [GPa]\n");
    fprintf(FilePtr,"-----------------------------\n");
    PrintRealMatrix6x6ToFile(&IdealGasTermVoigt,FilePtr,1e9/PRESSURE_CONVERSION_FACTOR);
    fprintf(FilePtr,"\n");
  }

  // output intermediate averages of the principle moments of inertia per component
  if(ComputePrincipleMomentsOfInertia)
  {
    fprintf(FilePtr,"Principle moments of inertia:\n");
    fprintf(FilePtr,"-----------------------------\n");
    for(k=0;k<NumberOfComponents;k++)
    {
      count=0.0;
      InertiaAverage.x=InertiaAverage.y=InertiaAverage.z=0.0;
      for(i=0;i<NR_BLOCKS;i++)
      {
        if(BlockWeightedCount[CurrentSystem][i]>0.0)
        {
          InertiaAverage.x+=PrincipleMomentsOfInertiaAccumulated[CurrentSystem][k][i].x;
          InertiaAverage.y+=PrincipleMomentsOfInertiaAccumulated[CurrentSystem][k][i].y;
          InertiaAverage.z+=PrincipleMomentsOfInertiaAccumulated[CurrentSystem][k][i].z;
          count+=PrincipleMomentsOfInertiaCount[CurrentSystem][k][i];
        }
      }
      fprintf(FilePtr,"Component [%s]: %g %g %g\n",
         Components[k].Name,InertiaAverage.x/count,InertiaAverage.y/count,InertiaAverage.z/count);
    }
    fprintf(FilePtr,"\n");
  }

  if(ComputeRattleSteps)
  {
    fprintf(FilePtr,"Average amount of RATTLE-steps stage 1: %g\n",NumberOfRattleCyclesStage1[CurrentSystem]/(REAL)CurrentCycle);
    fprintf(FilePtr,"Maximum amount of RATTLE-steps stage 1: %d\n",MaximumNumberOfRattleCyclesStage1[CurrentSystem]);
    fprintf(FilePtr,"Average amount of RATTLE-steps stage 2: %g\n",NumberOfRattleCyclesStage2[CurrentSystem]/(REAL)CurrentCycle);
    fprintf(FilePtr,"Maximum amount of RATTLE-steps stage 2: %d\n",MaximumNumberOfRattleCyclesStage2[CurrentSystem]);
  }

  fprintf(FilePtr,"\n");
}

void PrintIntervalStatusProduction(long long CurrentCycle,long long NumberOfCycles, FILE *FilePtr)
{
  int i,j,k;
  REAL number_of_unit_cells;
  int FractionalMolecule;
  REAL loading,average_loading,fac,fac2,temp;
  REAL Lambda,shift;
  REAL_MATRIX3x3 Stress;
  VECTOR com;

  fprintf(FilePtr,"Current cycle: %lld out of %lld\n",CurrentCycle,NumberOfCycles);
  fprintf(FilePtr,"========================================================================================================\n\n");
  if(SimulationType==MOLECULAR_DYNAMICS)
    fprintf(FilePtr,"Time run: %10.6lf [ps] %10.6lf [ns] \n",(double)(DeltaT*(REAL)(double)CurrentCycle),(double)(1e-3*DeltaT*(REAL)(double)CurrentCycle));

  fprintf(FilePtr,"Net charge: %g (F: %g, A: %g, C: %g)\n",NetChargeSystem[CurrentSystem],NetChargeFramework[CurrentSystem],
                                                           NetChargeAdsorbates[CurrentSystem],NetChargeCations[CurrentSystem]);

  // write out the boxlengths
  fprintf(FilePtr,"Current Box: %9.5lf %9.5lf %9.5lf [A]   Average Box: %9.5lf %9.5lf %9.5lf [A]\n",
          (double)Box[CurrentSystem].ax,(double)Box[CurrentSystem].bx,(double)Box[CurrentSystem].cx,
          (double)GetAverageWeightedProperty(BoxAXAccumulated),(double)GetAverageWeightedProperty(BoxBXAccumulated),(double)GetAverageWeightedProperty(BoxCXAccumulated));
  fprintf(FilePtr,"             %9.5lf %9.5lf %9.5lf [A]                %9.5lf %9.5lf %9.5lf [A]\n",
          (double)Box[CurrentSystem].ay,(double)Box[CurrentSystem].by,(double)Box[CurrentSystem].cy,
          (double)GetAverageWeightedProperty(BoxAYAccumulated),(double)GetAverageWeightedProperty(BoxBYAccumulated),(double)GetAverageWeightedProperty(BoxCYAccumulated));
  fprintf(FilePtr,"             %9.5lf %9.5lf %9.5lf [A]                %9.5lf %9.5lf %9.5lf [A]\n",
          (double)Box[CurrentSystem].az,(double)Box[CurrentSystem].bz,(double)Box[CurrentSystem].cz,
          (double)GetAverageWeightedProperty(BoxAZAccumulated),(double)GetAverageWeightedProperty(BoxBZAccumulated),(double)GetAverageWeightedProperty(BoxCZAccumulated));
  fprintf(FilePtr,"Box-lengths:  %9.5lf %9.5lf %9.5lf [A] Average: %9.5lf %9.5lf %9.5lf [A]\n",
          (double)BoxProperties[CurrentSystem].ax,
          (double)BoxProperties[CurrentSystem].ay,
          (double)BoxProperties[CurrentSystem].az,
          (double)GetAverageVectorProperty(BoxLengthAccumulated).x,
          (double)GetAverageVectorProperty(BoxLengthAccumulated).y,
          (double)GetAverageVectorProperty(BoxLengthAccumulated).z);
  fprintf(FilePtr,"Box-angles:  %9.5lf %9.5lf %9.5lf [degrees] Average: %9.5lf %9.5lf %9.5lf [degrees]\n",
          (double)AlphaAngle[CurrentSystem]*RAD2DEG,
          (double)BetaAngle[CurrentSystem]*RAD2DEG,
          (double)GammaAngle[CurrentSystem]*RAD2DEG,
          (double)GetAverageWeightedProperty(AlphaAngleAccumulated),
          (double)GetAverageWeightedProperty(BetaAngleAccumulated),
          (double)GetAverageWeightedProperty(GammaAngleAccumulated));
  fprintf(FilePtr,"Volume: %9.5lf [A^3] Average Volume: %9.5lf [A^3]\n\n",
          (double)Volume[CurrentSystem],(double)GetAverageWeightedProperty(VolumeAccumulated));
  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
  {
    com=GetFrameworkCenterOfMass();
    fprintf(FilePtr,"Framework center-of-mass drift: %18.10lf [A] %18.10lf [A] %18.10lf [A]\n",
      (double)(com.x-Framework[CurrentSystem].IntialCenterOfMassPosition.x),
      (double)(com.y-Framework[CurrentSystem].IntialCenterOfMassPosition.y),
      (double)(com.z-Framework[CurrentSystem].IntialCenterOfMassPosition.z));
  }

  // write out the StressTensor
  switch(SimulationType)
  {
    case MONTE_CARLO:
      if(ComputeMolecularPressure[CurrentSystem])
      {
        fac=PRESSURE_CONVERSION_FACTOR*1e-6;
        fac2=PRESSURE_CONVERSION_FACTOR*1e-3;
        Stress=GetAverageMolecularStressTensor();
        fprintf(FilePtr,"Current Stress: %9.5lf %9.5lf %9.5lf [MPa]   Average Stress: %9.5lf %9.5lf %9.5lf [kPa]\n",
                (double)fac*MolecularStressTensor[CurrentSystem].ax,(double)fac*MolecularStressTensor[CurrentSystem].bx,(double)fac*MolecularStressTensor[CurrentSystem].cx,
                (double)fac2*Stress.ax,(double)fac2*Stress.bx,(double)fac2*Stress.cx);
        fprintf(FilePtr,"                %9.5lf %9.5lf %9.5lf [MPa]                   %9.5lf %9.5lf %9.5lf [kPa]\n",
                (double)fac*MolecularStressTensor[CurrentSystem].ay,(double)fac*MolecularStressTensor[CurrentSystem].by,(double)fac*MolecularStressTensor[CurrentSystem].cy,
                (double)fac2*Stress.ay,(double)fac2*Stress.by,(double)fac2*Stress.cy);
        fprintf(FilePtr,"                %9.5lf %9.5lf %9.5lf [MPa]                   %9.5lf %9.5lf %9.5lf [kPa]\n",
                (double)fac*MolecularStressTensor[CurrentSystem].az,(double)fac*MolecularStressTensor[CurrentSystem].bz,(double)fac*MolecularStressTensor[CurrentSystem].cz,
                (double)fac2*Stress.az,(double)fac2*Stress.bz,(double)fac2*Stress.cz);
        fprintf(FilePtr,"Average pressure: %18.10lf [kPa] %18.10lf [atm] %18.10lf [bar]\n",
           (double)(fac2*(Stress.ax+Stress.by+Stress.cz)/3.0),
           (double)(fac2*(Stress.ax+Stress.by+Stress.cz)*KPA_TO_ATM/3.0),
           (double)(fac2*(Stress.ax+Stress.by+Stress.cz)*KPA_TO_BAR/3.0));
        temp=GetAverageIdealGasPartPressure();
        fprintf(FilePtr,"\tIdeal gas part:   %18.10lf [kPa] %18.10lf [atm] %18.10lf [bar]\n",
                temp*fac2,temp*fac2*KPA_TO_ATM,temp*fac2*KPA_TO_BAR);
        temp=GetAverageExcessPartPressure();
        fprintf(FilePtr,"\texcess part:      %18.10lf [kPa] %18.10lf [atm] %18.10lf [bar]\n",
                temp*fac2,temp*fac2*KPA_TO_ATM,temp*fac2*KPA_TO_BAR);
        temp=GetAverageTailCorrectionPressure();
        fprintf(FilePtr,"\ttail correction: %18.10lf [kPa] %18.10lf [atm] %18.10lf [bar]\n",
                temp*fac2,temp*fac2*KPA_TO_ATM,temp*fac2*KPA_TO_BAR);
        fprintf(FilePtr,"\n");
      }
      break;
    case MOLECULAR_DYNAMICS:
      com=MeasureVelocityDrift();
      fprintf(FilePtr,"Total linear momentum: %18.10lf [A] %18.10lf [A] %18.10lf [A]\n",(double)com.x,(double)com.y,(double)com.z);
      if(Framework[CurrentSystem].NumberOfCoreShellDefinitions>0)
        fprintf(FilePtr,"Core/Shell temperature: %g\n",GetCoreShellTemperature());
      fprintf(FilePtr,"\n");

      fac=PRESSURE_CONVERSION_FACTOR*1e-6;
      fac2=PRESSURE_CONVERSION_FACTOR*1e-3;
      Stress=ComputeCurrentAverageStressTensor();
      fprintf(FilePtr,"Current Stress: %9.5lf %9.5lf %9.5lf [MPa]   Average Stress: %9.5lf %9.5lf %9.5lf [kPa]\n",
              (double)fac*StressTensor[CurrentSystem].ax,(double)fac*StressTensor[CurrentSystem].bx,(double)fac*StressTensor[CurrentSystem].cx,
              (double)fac2*Stress.ax,(double)fac2*Stress.bx,(double)fac2*Stress.cx);
      fprintf(FilePtr,"                %9.5lf %9.5lf %9.5lf [MPa]                   %9.5lf %9.5lf %9.5lf [kPa]\n",
              (double)fac*StressTensor[CurrentSystem].ay,(double)fac*StressTensor[CurrentSystem].by,(double)fac*StressTensor[CurrentSystem].cy,
              (double)fac2*Stress.ay,(double)fac2*Stress.by,(double)fac2*Stress.cy);
      fprintf(FilePtr,"                %9.5lf %9.5lf %9.5lf [MPa]                   %9.5lf %9.5lf %9.5lf [kPa]\n",
              (double)fac*StressTensor[CurrentSystem].az,(double)fac*StressTensor[CurrentSystem].bz,(double)fac*StressTensor[CurrentSystem].cz,
              (double)fac2*Stress.az,(double)fac2*Stress.bz,(double)fac2*Stress.cz);
      fprintf(FilePtr,"Average pressure: %18.10lf [kPa] %18.10lf [atm] %18.10lf [bar]\n",
         (double)(fac2*(Stress.ax+Stress.by+Stress.cz)/3.0),
         (double)(fac2*(Stress.ax+Stress.by+Stress.cz)*KPA_TO_ATM/3.0),
         (double)(fac2*(Stress.ax+Stress.by+Stress.cz)*KPA_TO_BAR/3.0));
      fprintf(FilePtr,"\n");
      break;
    default:
      break;
  }

  if(NumberOfReactions>0)
  {
    fprintf(FilePtr,"Reactions:\n");
    fprintf(FilePtr,"----------------------------------------------------------------------------------------------------------------------------------------------------\n");
    for(i=0;i<NumberOfReactions;i++)
    {
      fprintf(FilePtr,"Reaction %d, current Lambda: %18.10f, maximum Lambda-change: %18.10f\n",i,CFCRXMCLambda[CurrentSystem][i],MaximumReactionLambdaChange[CurrentSystem][i]);
      fprintf(FilePtr,"\tFractional molecules: ");
      for(j=0;j<NumberOfComponents;j++)
      {
        if(ProductsStoichiometry[i][j]>0)
        {
          for(k=0;k<ProductsStoichiometry[i][j];k++)
            fprintf(FilePtr," %d (%s)",Components[j].ProductFractionalMolecules[CurrentSystem][i][k],Components[j].Name);
        }
      }
      fprintf(FilePtr," <-->");
      for(j=0;j<NumberOfComponents;j++)
      {
        if(ReactantsStoichiometry[i][j]>0)
        {
          for(k=0;k<ReactantsStoichiometry[i][j];k++)
            fprintf(FilePtr," %d (%s)",Components[j].ReactantFractionalMolecules[CurrentSystem][i][k],Components[j].Name);
        }
      }
      fprintf(FilePtr,"\n");

      shift=RXMCBiasingFactors[CurrentSystem][i][0];
      fprintf(FilePtr,"\tBiasing Factors: ");
      for(k=0;k<RXMCLambdaHistogramSize;k++)
      {
        fprintf(FilePtr,"%4f ",RXMCBiasingFactors[CurrentSystem][i][k]-shift);
        if((k+1)%10==0&&(k+1)!=RXMCLambdaHistogramSize)  fprintf(FilePtr,"\n\t                 ");
      }
      fprintf(FilePtr,"\n");
    }
    fprintf(FilePtr,"\n");
  }


  number_of_unit_cells=NumberOfUnitCells[CurrentSystem].x*NumberOfUnitCells[CurrentSystem].y*NumberOfUnitCells[CurrentSystem].z;
  if(Framework[CurrentSystem].FrameworkModel==NONE)
  {
    fprintf(FilePtr,"Amount of molecules per component:\n");
    fprintf(FilePtr,"----------------------------------------------------------------------------------------------------------------------------------------------------\n");
    for(i=0;i<NumberOfComponents;i++)
    {
      fprintf(FilePtr,"Component %d (%s), current number of integer/fractional/reaction molecules: %d/%d/%d (average %9.5lf/%9.5lf), density: %9.5lf (average %9.5lf) kg/m^3]\n",
        i,
        Components[i].Name,
        Components[i].NumberOfMolecules[CurrentSystem]-(Components[i].FractionalMolecule[CurrentSystem]>=0?1:0)-Components[i].NumberOfRXMCMoleculesPresent[CurrentSystem],
        (Components[i].FractionalMolecule[CurrentSystem]>=0?1:0),
        Components[i].NumberOfRXMCMoleculesPresent[CurrentSystem],
        (double)(GetAverageWeightedComponentProperty(NumberOfIntegerMoleculesPerComponentAccumulated,i)/(REAL)number_of_unit_cells),
        (double)GetAverageComponentProperty(NumberOfFractionalMoleculesPerComponentAccumulated,i),
        (double)((Components[i].Mass*(REAL)(Components[i].NumberOfMolecules[CurrentSystem]-(Components[i].FractionalMolecule[CurrentSystem]>=0?1:0)-Components[i].NumberOfRXMCMoleculesPresent[CurrentSystem])/
                 Volume[CurrentSystem])*DENSITY_CONVERSION_FACTOR),
        (double)(GetAverageWeightedComponentProperty(DensityPerComponentAccumulated,i)*DENSITY_CONVERSION_FACTOR));
      if(Components[i].CFMoleculePresent[CurrentSystem])
      {
        char FractionalMoleculeString[6];
        FractionalMolecule=Components[i].FractionalMolecule[CurrentSystem];
        if (FractionalMolecule<0)
          strcpy(FractionalMoleculeString,"-");
        else
          sprintf(FractionalMoleculeString,"%d", FractionalMolecule);
        fprintf(FilePtr,"\tFractional molecule-id: %s, max. Lambda-change: %5lf (CFMC) %5lf (CB/CFMC)\n",FractionalMoleculeString,MaximumCFLambdaChange[CurrentSystem][i],MaximumCBCFLambdaChange[CurrentSystem][i]);
        fprintf(FilePtr,"\tLambda factors: ");
        for(k=0;k<Components[i].NumberOfAtoms;k++)
        {
          Lambda=0.0;
          if(FractionalMolecule>=0)
          {
            if(Components[i].ExtraFrameworkMolecule)
              Lambda=Cations[CurrentSystem][FractionalMolecule].Atoms[k].CFVDWScalingParameter;
            else
              Lambda=Adsorbates[CurrentSystem][FractionalMolecule].Atoms[k].CFVDWScalingParameter;
          }
          fprintf(FilePtr,"%4f ",Lambda);
          if((k+1)%10==0&&(k+1)!=Components[i].NumberOfAtoms)  fprintf(FilePtr,"\n\t                ");
        }
        fprintf(FilePtr,"\n");
        shift=Components[i].CFBiasingFactors[CurrentSystem][0];
        shift=0;
        fprintf(FilePtr,"\tBiasing Factors: ");
        for(k=0;k<Components[i].CFLambdaHistogramSize;k++)
        {
          fprintf(FilePtr,"%4f ",Components[i].CFBiasingFactors[CurrentSystem][k]-shift);
          if((k+1)%10==0&&(k+1)!=Components[i].CFLambdaHistogramSize)  fprintf(FilePtr,"\n\t                 ");
        }
        fprintf(FilePtr,"\n");
      }
    }
    fprintf(FilePtr,"----------------------------------------------------------------------------------------------------------------------------------------------------\n");
  }
  else
  {
    fprintf(FilePtr,"Loadings per component:\n");
    fprintf(FilePtr,"----------------------------------------------------------------------------------------------------------------------------------------------------\n");
    for(i=0;i<NumberOfComponents;i++)
    {
      fprintf(FilePtr,"Component %d (%s), current number of integer/fractional/reaction molecules: %d/%d/%d (avg. %9.5lf), density: %9.5lf (avg. %9.5lf) [kg/m^3]\n",
        i,
        Components[i].Name,
        Components[i].NumberOfMolecules[CurrentSystem]-(Components[i].FractionalMolecule[CurrentSystem]>=0?1:0)-Components[i].NumberOfRXMCMoleculesPresent[CurrentSystem],
        (Components[i].FractionalMolecule[CurrentSystem]>=0?1:0),
        Components[i].NumberOfRXMCMoleculesPresent[CurrentSystem],
        (double)GetAverageWeightedComponentProperty(NumberOfIntegerMoleculesPerComponentAccumulated,i),
        (double)(Components[i].Mass*(REAL)(Components[i].NumberOfMolecules[CurrentSystem]-(Components[i].FractionalMolecule[CurrentSystem]>=0?1:0)-Components[i].NumberOfRXMCMoleculesPresent[CurrentSystem])/
                Volume[CurrentSystem])*DENSITY_CONVERSION_FACTOR,
        (double)GetAverageWeightedComponentProperty(DensityPerComponentAccumulated,i)*DENSITY_CONVERSION_FACTOR);
      if(Components[i].CFMoleculePresent[CurrentSystem])
      {
        FractionalMolecule=Components[i].FractionalMolecule[CurrentSystem];
        fprintf(FilePtr,"\tFractional molecule-id: %d, max. Lambda-change: %5lf (CFMC) %5lf (CB/CFMC)\n",FractionalMolecule,MaximumCFLambdaChange[CurrentSystem][i],MaximumCBCFLambdaChange[CurrentSystem][i]);
        fprintf(FilePtr,"\tLambda factors: ");
        for(k=0;k<Components[i].NumberOfAtoms;k++)
        {
          if(Components[i].ExtraFrameworkMolecule)
            Lambda=Cations[CurrentSystem][FractionalMolecule].Atoms[k].CFVDWScalingParameter;
          else
            Lambda=Adsorbates[CurrentSystem][FractionalMolecule].Atoms[k].CFVDWScalingParameter;
          fprintf(FilePtr,"%4f ",Lambda);
          if((k+1)%10==0&&(k+1)!=Components[i].NumberOfAtoms)  fprintf(FilePtr,"\n\t                ");
        }
        fprintf(FilePtr,"\n");
        shift=Components[i].CFBiasingFactors[CurrentSystem][0];
        fprintf(FilePtr,"\tBiasing Factors: ");
        for(k=0;k<Components[i].CFLambdaHistogramSize;k++)
        {
          fprintf(FilePtr,"%4f ",Components[i].CFBiasingFactors[CurrentSystem][k]-shift);
          if((k+1)%10==0&&(k+1)!=Components[i].CFLambdaHistogramSize)  fprintf(FilePtr,"\n\t                 ");
        }
        fprintf(FilePtr,"\n");
      }

      loading=(REAL)(Components[i].NumberOfMolecules[CurrentSystem]
                   -(Components[i].FractionalMolecule[CurrentSystem]>=0?1:0)
                   -Components[i].NumberOfRXMCMoleculesPresent[CurrentSystem])/(REAL)number_of_unit_cells;
      average_loading=GetAverageWeightedComponentProperty(NumberOfIntegerMoleculesPerComponentAccumulated,i)/(REAL)number_of_unit_cells;
      fprintf(FilePtr,"\tabsolute adsorption: %9.5lf (avg. %9.5lf) [mol/uc], %14.10lf (avg. %14.10lf) [mol/kg], %14.10lf (avg. %14.10lf) [mg/g]\n",
         (double)loading,
         (double)average_loading,
         (double)(Components[i].MOLEC_PER_UC_TO_MOL_PER_KG[CurrentSystem]*loading),
         (double)(Components[i].MOLEC_PER_UC_TO_MOL_PER_KG[CurrentSystem]*average_loading),
         (double)(Components[i].MOLEC_PER_UC_TO_MILLIGRAM_PER_GRAM_OF_FRAMEWORK[CurrentSystem]*loading),
         (double)(Components[i].MOLEC_PER_UC_TO_MILLIGRAM_PER_GRAM_OF_FRAMEWORK[CurrentSystem]*average_loading));
      fprintf(FilePtr,"\t                     %14.10lf (avg. %14.10lf) [cm^3 STP/g],  %14.10lf (avg. %14.10lf) [cm^3 STP/cm^3]\n",
         (double)(Components[i].MOLEC_PER_UC_TO_CC_STP_G[CurrentSystem]*loading),
         (double)(Components[i].MOLEC_PER_UC_TO_CC_STP_G[CurrentSystem]*average_loading),
         (double)(Components[i].MOLEC_PER_UC_TO_CC_STP_CC[CurrentSystem]*loading),
         (double)(Components[i].MOLEC_PER_UC_TO_CC_STP_CC[CurrentSystem]*average_loading));


      loading=((REAL)Components[i].NumberOfMolecules[CurrentSystem]
                    -Components[i].AmountOfExcessMolecules[CurrentSystem]
                    -(Components[i].FractionalMolecule[CurrentSystem]>=0?1:0)
                    -Components[i].NumberOfRXMCMoleculesPresent[CurrentSystem])/(REAL)number_of_unit_cells;
      average_loading=GetAverageWeightedComponentProperty(NumberOfExcessMoleculesPerComponentAccumulated,i)/(REAL)number_of_unit_cells;
      fprintf(FilePtr,"\texcess adsorption:   %14.10lf (avg. %14.10lf) [mol/uc], %14.10lf (avg. %14.10lf) [mol/kg], %14.10lf (avg. %14.10lf) [mg/g]\n",
         (double)loading,
         (double)average_loading,
         (double)(Components[i].MOLEC_PER_UC_TO_MOL_PER_KG[CurrentSystem]*loading),
         (double)(Components[i].MOLEC_PER_UC_TO_MOL_PER_KG[CurrentSystem]*average_loading),
         (double)(Components[i].MOLEC_PER_UC_TO_MILLIGRAM_PER_GRAM_OF_FRAMEWORK[CurrentSystem]*loading),
         (double)(Components[i].MOLEC_PER_UC_TO_MILLIGRAM_PER_GRAM_OF_FRAMEWORK[CurrentSystem]*average_loading));
      fprintf(FilePtr,"\t                     %14.10lf (avg. %14.10lf) [cm^3 STP/g],  %14.10lf (avg. %14.10lf) [cm^3 STP/cm^3]\n",
         (double)(Components[i].MOLEC_PER_UC_TO_CC_STP_G[CurrentSystem]*loading),
         (double)(Components[i].MOLEC_PER_UC_TO_CC_STP_G[CurrentSystem]*average_loading),
         (double)(Components[i].MOLEC_PER_UC_TO_CC_STP_CC[CurrentSystem]*loading),
         (double)(Components[i].MOLEC_PER_UC_TO_CC_STP_CC[CurrentSystem]*average_loading));
    }
    fprintf(FilePtr,"----------------------------------------------------------------------------------------------------------------------------------------------------\n");
  }

  fprintf(FilePtr,"Degrees of freedom: %d %d %d %d\n",DegreesOfFreedom[CurrentSystem],DegreesOfFreedomFramework[CurrentSystem],
          DegreesOfFreedomAdsorbates[CurrentSystem],DegreesOfFreedomCations[CurrentSystem]);
  fprintf(FilePtr,"Number of Framework-atoms: %6d\n",Framework[CurrentSystem].TotalNumberOfAtoms);
  fprintf(FilePtr,"Number of Adsorbates:      %6d (%d integer, %d fractional, %d reaction)\n",
          NumberOfAdsorbateMolecules[CurrentSystem],
          TotalNumberOfIntegerAdsorbates(),
          TotalNumberOfFractionalAdsorbates(),
          NumberOfReactionAdsorbateMolecules[CurrentSystem]);
  fprintf(FilePtr,"Number of Cations:         %6d (%d integer, %d fractional, %d reaction)\n",
          NumberOfCationMolecules[CurrentSystem],
          TotalNumberOfIntegerCations(),
          TotalNumberOfFractionalCations(),
          NumberOfReactionCationMolecules[CurrentSystem]);
  fprintf(FilePtr,"\n");

  if(SimulationType==MOLECULAR_DYNAMICS)
  {
    fprintf(FilePtr,"Conserved energy: % 22.10lf Energy drifts: % 11.10lf % 22.10lf\n",
        (double)(ConservedEnergy[CurrentSystem]*ENERGY_TO_KELVIN),
        (double)(fabs((ConservedEnergy[CurrentSystem]-ReferenceEnergy[CurrentSystem])/ReferenceEnergy[CurrentSystem])*ENERGY_TO_KELVIN),
        (double)(Drift[CurrentSystem]/(CurrentCycle+1.0)));
    if(DegreesOfFreedom[CurrentSystem]>0)
    {
      fprintf(FilePtr,"Temperature:            % 8.3lf (avg. % 8.3lf), Translational (avg. % 8.3lf), Rotational (avg. % 8.3lf)\n",
        (double)(2.0*UKinetic[CurrentSystem]/(K_B*DegreesOfFreedom[CurrentSystem])),
        (double)GetAverageWeightedProperty(TemperatureAccumulated),
        (double)GetAverageWeightedProperty(TemperatureTranslationAccumulated),
        (double)GetAverageWeightedProperty(TemperatureRotationAccumulated));
    }
    if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
      fprintf(FilePtr,"Temperature Framework:  % 8.3lf (avg. % 8.3lf)\n",
        (double)(2.0*UHostKinetic[CurrentSystem]/(K_B*DegreesOfFreedomFramework[CurrentSystem])),
        (double)GetAverageWeightedProperty(TemperatureFrameworkAccumulated));
    if(DegreesOfFreedomAdsorbates[CurrentSystem]>0)
        fprintf(FilePtr,"Temperature Adsorbates: % 8.3lf (avg. % 8.3lf), Translational (avg. % 8.3lf), Rotational (avg. % 8.3lf)\n",
        (double)(2.0*UAdsorbateKinetic[CurrentSystem]/(K_B*DegreesOfFreedomAdsorbates[CurrentSystem])),
        (double)GetAverageWeightedProperty(TemperatureAdsorbatesAccumulated),
        (double)GetAverageWeightedProperty(TemperatureTranslationAdsorbateAccumulated),
        (double)GetAverageWeightedProperty(TemperatureRotationAdsorbateAccumulated));
    if(DegreesOfFreedomCations[CurrentSystem]>0)
      fprintf(FilePtr,"Temperature Cations:    % 8.3lf (avg. % 8.3lf)\n",
        (double)(2.0*UCationKinetic[CurrentSystem]/(K_B*DegreesOfFreedomCations[CurrentSystem])),
        (double)GetAverageWeightedProperty(TemperatureCationsAccumulated));
    fprintf(FilePtr,"Cell temperature: % 8.3lf (avg. % 8.3lf)\n",GetCellTemperature(),GetAverageCellTemperature());
    fprintf(FilePtr,"Current total kinetic energy:       % 22.10lf [K]\n",(double)UKinetic[CurrentSystem]);
    fprintf(FilePtr,"Current total Nose-Hoover energy:   % 22.10lf [K]\n",(double)UNoseHoover[CurrentSystem]);
  }
  fprintf(FilePtr,"Current total potential energy:       % 22.10lf [K]  (avg. % 22.10lf)\n",
      (double)UTotal[CurrentSystem]*ENERGY_TO_KELVIN,(double)GetAverageWeightedProperty(UTotalAccumulated)*ENERGY_TO_KELVIN);

  fprintf(FilePtr,"\tCurrent Host-Host energy:           % 22.10lf [K]  (avg. % 22.10lf)\n",
      (double)UHostHost[CurrentSystem]*ENERGY_TO_KELVIN,(double)GetAverageWeightedProperty(UHostHostAccumulated)*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tCurrent Host-Adsorbate energy:      % 22.10lf [K]  (avg. % 22.10lf)\n",
      (double)UHostAdsorbate[CurrentSystem]*ENERGY_TO_KELVIN,(double)GetAverageWeightedProperty(UHostAdsorbateAccumulated)*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tCurrent Host-Cation energy:         % 22.10lf [K]  (avg. % 22.10lf)\n",
      (double)UHostCation[CurrentSystem]*ENERGY_TO_KELVIN,(double)GetAverageWeightedProperty(UHostCationAccumulated)*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tCurrent Adsorbate-Adsorbate energy: % 22.10lf [K]  (avg. % 22.10lf)\n",
      (double)UAdsorbateAdsorbate[CurrentSystem]*ENERGY_TO_KELVIN,(double)GetAverageWeightedProperty(UAdsorbateAdsorbateAccumulated)*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tCurrent Cation-Cation energy:       % 22.10lf [K]  (avg. % 22.10lf)\n",
      (double)UCationCation[CurrentSystem]*ENERGY_TO_KELVIN,(double)GetAverageWeightedProperty(UCationCationAccumulated)*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tCurrent Adsorbate-Cation energy:    % 22.10lf [K]  (avg. % 22.10lf)\n",
      (double)UAdsorbateCation[CurrentSystem]*ENERGY_TO_KELVIN,(double)GetAverageWeightedProperty(UAdsorbateCationAccumulated)*ENERGY_TO_KELVIN);

  if(ComputePolarization)
  {
    fprintf(FilePtr,"\tCurrent polarization energy:        % 22.10lf [K]  (avg. % 22.10lf)\n",(double)(UHostPolarization[CurrentSystem]+
         UAdsorbatePolarization[CurrentSystem]+UCationPolarization[CurrentSystem])*ENERGY_TO_KELVIN,
         (double)(GetAverageWeightedProperty(UHostPolarizationAccumulated)+GetAverageWeightedProperty(UAdsorbatePolarizationAccumulated)+
                 GetAverageWeightedProperty(UCationPolarizationAccumulated))*ENERGY_TO_KELVIN);
    fprintf(FilePtr,"\tCurrent back-polarization energy:   % 22.10lf [K]  (avg. % 22.10lf)\n",(double)(UHostBackPolarization[CurrentSystem]+
         UAdsorbateBackPolarization[CurrentSystem]+UCationBackPolarization[CurrentSystem])*ENERGY_TO_KELVIN,
         (double)(GetAverageWeightedProperty(UHostBackPolarizationAccumulated)+GetAverageWeightedProperty(UAdsorbateBackPolarizationAccumulated)+
                 GetAverageWeightedProperty(UCationBackPolarizationAccumulated))*ENERGY_TO_KELVIN);
  }

  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
  {
    if(Framework[CurrentSystem].NumberOfBondsDefinitions>0)
      fprintf(FilePtr,"\t\tCurrent Host-Bond energy:             % 22.10lf [K]  (avg. % 22.10lf)\n",
          (double)UHostBond[CurrentSystem]*ENERGY_TO_KELVIN,(double)GetAverageWeightedProperty(UHostBondAccumulated)*ENERGY_TO_KELVIN);
    if(Framework[CurrentSystem].NumberOfUreyBradleyDefinitions>0)
      fprintf(FilePtr,"\t\tCurrent Host-UreyBradly energy:       % 22.10lf [K]  (avg. % 22.10lf)\n",
          (double)UHostUreyBradley[CurrentSystem]*ENERGY_TO_KELVIN,(double)GetAverageWeightedProperty(UHostUreyBradleyAccumulated)*ENERGY_TO_KELVIN);
    if(Framework[CurrentSystem].NumberOfBendDefinitions>0)
      fprintf(FilePtr,"\t\tCurrent Host-Bend energy:             % 22.10lf [K]  (avg. % 22.10lf)\n",
          (double)UHostBend[CurrentSystem]*ENERGY_TO_KELVIN,(double)GetAverageWeightedProperty(UHostBendAccumulated)*ENERGY_TO_KELVIN);
    if(Framework[CurrentSystem].NumberOfInversionBendDefinitions>0)
      fprintf(FilePtr,"\t\tCurrent Host-Inversion Bend energy:   % 22.10lf [K]  (avg. % 22.10lf)\n",
          (double)UHostInversionBend[CurrentSystem]*ENERGY_TO_KELVIN,(double)GetAverageWeightedProperty(UHostInversionBendAccumulated)*ENERGY_TO_KELVIN);
    if(Framework[CurrentSystem].NumberOfTorsionDefinitions>0)
      fprintf(FilePtr,"\t\tCurrent Host-Torsion energy:          % 22.10lf [K]  (avg. % 22.10lf)\n",
          (double)UHostTorsion[CurrentSystem]*ENERGY_TO_KELVIN,(double)GetAverageWeightedProperty(UHostTorsionAccumulated)*ENERGY_TO_KELVIN);
    if(Framework[CurrentSystem].NumberOfImproperTorsionDefinitions>0)
      fprintf(FilePtr,"\t\tCurrent Host-Improper torsion energy: % 22.10lf [K]  (avg. % 22.10lf)\n",
          (double)UHostImproperTorsion[CurrentSystem]*ENERGY_TO_KELVIN,(double)GetAverageWeightedProperty(UHostImproperTorsionAccumulated)*ENERGY_TO_KELVIN);
    if(Framework[CurrentSystem].NumberOfBondBondDefinitions>0)
      fprintf(FilePtr,"\t\tCurrent Host-Bond/Bond energy:        % 22.10lf [K]  (avg. % 22.10lf)\n",
          (double)UHostBondBond[CurrentSystem]*ENERGY_TO_KELVIN,(double)GetAverageWeightedProperty(UHostBondBondAccumulated)*ENERGY_TO_KELVIN);
    if(Framework[CurrentSystem].NumberOfBendBendDefinitions>0)
      fprintf(FilePtr,"\t\tCurrent Host-Bend/Bend energy:        % 22.10lf [K]  (avg. % 22.10lf)\n",
          (double)UHostBendBend[CurrentSystem]*ENERGY_TO_KELVIN,(double)GetAverageWeightedProperty(UHostBendBendAccumulated)*ENERGY_TO_KELVIN);
    if(Framework[CurrentSystem].NumberOfBondBendDefinitions>0)
      fprintf(FilePtr,"\t\tCurrent Host-Bond/Bend energy:        % 22.10lf [K]  (avg. % 22.10lf)\n",
          (double)UHostBondBend[CurrentSystem]*ENERGY_TO_KELVIN,(double)GetAverageWeightedProperty(UHostBondBendAccumulated)*ENERGY_TO_KELVIN);
    if(Framework[CurrentSystem].NumberOfBondTorsionDefinitions>0)
      fprintf(FilePtr,"\t\tCurrent Host-Bond/Torsion energy:     % 22.10lf [K]  (avg. % 22.10lf)\n",
          (double)UHostBondTorsion[CurrentSystem]*ENERGY_TO_KELVIN,(double)GetAverageWeightedProperty(UHostBondTorsionAccumulated)*ENERGY_TO_KELVIN);
    if(Framework[CurrentSystem].NumberOfBendTorsionDefinitions>0)
      fprintf(FilePtr,"\t\tCurrent Host-Bend/Torsion energy:     % 22.10lf [K]  (avg. % 22.10lf)\n\n",
          (double)UHostBendTorsion[CurrentSystem]*ENERGY_TO_KELVIN,(double)GetAverageWeightedProperty(UHostBendTorsionAccumulated)*ENERGY_TO_KELVIN);
  }

  fprintf(FilePtr,"\n");
  PrintWarningStatus();
  fprintf(FilePtr,"\n\n");
  fflush(FilePtr);
}

void PrintProperty(FILE *FilePtr,char *string,char *units,REAL conv_factor,REAL **Property)
{
  int i;
  REAL sum,sum_squared;
  REAL avg,error;

  fprintf(FilePtr,"\n");
  fprintf(FilePtr,"%s",string);
  sum=sum_squared=0.0;
  for(i=0;i<NR_BLOCKS;i++)
  {
    if(BlockWeightedCount[CurrentSystem][i]>0.0)
    {
      REAL tmp=conv_factor*(Property[CurrentSystem][i]/BlockWeightedCount[CurrentSystem][i]);
      sum+=tmp;
      sum_squared+=SQR(tmp);
      fprintf(FilePtr,"\tBlock[%2d] %18.5lf %s\n",i,(double)tmp,units);
    }
    else
      fprintf(FilePtr,"\tBlock[%2d] %18.5lf %s\n",i,(double)0.0,units);
  }
  avg=AVERAGE(sum);
  error=ERROR_CONFIDENCE_INTERVAL_95(sum,sum_squared);
  fprintf(FilePtr,"\t------------------------------------------------------------------------------\n");
  fprintf(FilePtr,"\tAverage   %18.5lf %s +/- %18.5lf %s\n",avg,units,error,units);

}

void PrintEnergies(FILE *FilePtr,char *string,char *units,REAL conv_factor,REAL **Property,
                   REAL **PropertyVDW,REAL **PropertyCoulomb)
{
  int i;
  REAL sum,tmp,sum_vdw,sum_coul,tmp_vdw,tmp_coul;
  REAL sum_squared, sum_vdw_squared, sum_coul_squared;

  fprintf(FilePtr,"\n");
  fprintf(FilePtr,"%s",string);
  sum=sum_vdw=sum_coul=0.0;
  sum_squared=sum_vdw_squared=sum_coul_squared=0.0;
  for(i=0;i<NR_BLOCKS;i++)
  {
    if(BlockWeightedCount[CurrentSystem][i]>0.0)
    {
      tmp=conv_factor*(Property[CurrentSystem][i]/BlockWeightedCount[CurrentSystem][i]);
      tmp_vdw=conv_factor*(PropertyVDW[CurrentSystem][i]/BlockWeightedCount[CurrentSystem][i]);
      tmp_coul=conv_factor*(PropertyCoulomb[CurrentSystem][i]/BlockWeightedCount[CurrentSystem][i]);

      sum+=tmp;
      sum_vdw+=tmp_vdw;
      sum_coul+=tmp_coul;

      sum_squared+=SQR(tmp);
      sum_vdw_squared+=SQR(tmp_vdw);
      sum_coul_squared+=SQR(tmp_coul);

      fprintf(FilePtr,"\tBlock[%2d] %-18.5lf Van der Waals: %-18.5lf Coulomb: %-18.5lf %s\n",
              i,(double)tmp,(double)tmp_vdw,(double)tmp_coul,units);
    }
    else
      fprintf(FilePtr,"\tBlock[%2d] %-18.5lf Van der Waals: %-18.5lf Coulomb: %-18.5lf %s\n",i,(double)0.0,(double)0.0,(double)0.0,units);
  }
  REAL avg=AVERAGE(sum);
  REAL avg_vdw=AVERAGE(sum_vdw);
  REAL avg_coul=AVERAGE(sum_coul);
  REAL error=ERROR_CONFIDENCE_INTERVAL_95(sum,sum_squared);
  REAL error_vdw=ERROR_CONFIDENCE_INTERVAL_95(sum_vdw,sum_vdw_squared);
  REAL error_coul=ERROR_CONFIDENCE_INTERVAL_95(sum_coul,sum_coul_squared);
  fprintf(FilePtr,"\t------------------------------------------------------------------------------\n");
  fprintf(FilePtr,"\tAverage   %-18.5lf Van der Waals: %-18lf Coulomb: %-18.5lf %s\n", avg, avg_vdw, avg_coul, units);
  fprintf(FilePtr,"\t      +/- %-18.5lf            +/- %-18lf      +/- %-18.5lf %s\n", error, error_vdw, error_coul, units);
}


void PrintAverageTotalSystemEnergiesMC(FILE *FilePtr)
{
  REAL sum,sum_squared; 
  REAL avg,error;

  fprintf(FilePtr,"\n\n\n");
  fprintf(FilePtr,"Average properties of the system[%d]:\n",CurrentSystem);
  fprintf(FilePtr,"========================================================================\n");

  // Temperature
  fprintf(FilePtr,"\n");
  fprintf(FilePtr,"Average temperature:\n");
  fprintf(FilePtr,"====================\n");
  sum=sum_squared=0.0;
  for(int i=0;i<NR_BLOCKS;i++)
  {
    if(BlockWeightedCount[CurrentSystem][i]>0.0)
    {
      REAL tmp=TemperatureAccumulated[CurrentSystem][i]/BlockWeightedCount[CurrentSystem][i];
      sum+=tmp;
      sum_squared+=SQR(tmp);
      fprintf(FilePtr,"\tBlock[%2d] %18.5lf [K]\n",i,(double)tmp);
    }
    else
      fprintf(FilePtr,"\tBlock[%2d] %18.5lf [K]\n",i,(double)0.0);
  }
  avg=AVERAGE(sum);
  error=ERROR_CONFIDENCE_INTERVAL_95(sum,sum_squared);
  fprintf(FilePtr,"\t------------------------------------------------------------------------------\n");
  fprintf(FilePtr,"\tAverage   %18.5lf [K] +/- %18.5lf [K]\n", avg, error);

  // Pressure
  fprintf(FilePtr,"\n");
  switch(SimulationType)
  {
    case MONTE_CARLO:
      fprintf(FilePtr,"Average Pressure:\n");
      fprintf(FilePtr,"=================\n");
      sum=sum_squared=0.0;
      for(int i=0;i<NR_BLOCKS;i++)
      {
        if(BlockWeightedCount[CurrentSystem][i]>0.0)
        {
          REAL tmp=((MolecularStressTensorAccumulated[CurrentSystem][i].ax+MolecularStressTensorAccumulated[CurrentSystem][i].by+MolecularStressTensorAccumulated[CurrentSystem][i].cz)/
               (3.0*BlockWeightedCount[CurrentSystem][i]))*PRESSURE_CONVERSION_FACTOR;

          sum+=tmp;
          sum_squared+=SQR(tmp);
          fprintf(FilePtr,"\tBlock[%2d] %18.5lf [Pa]\n",i,(double)tmp);
        }
        else
          fprintf(FilePtr,"\tBlock[%2d] %18.5lf [Pa]\n",i,(double)0.0);
      }
      avg=AVERAGE(sum);
      error=ERROR_CONFIDENCE_INTERVAL_95(sum,sum_squared);
      fprintf(FilePtr,"\t------------------------------------------------------------------------------\n");
      fprintf(FilePtr,"\tAverage   %18.5lf [Pa] +/- %18.5lf [Pa]\n", avg, error);
      fprintf(FilePtr,"\tAverage   %18.5lf [bar] +/- %18.5lf [bar]\n", avg*PA_TO_BAR, error*PA_TO_BAR);
      fprintf(FilePtr,"\tAverage   %18.5lf [atm] +/- %18.5lf [atm]\n", avg*PA_TO_ATM, error*PA_TO_ATM);
      fprintf(FilePtr,"\tAverage   %18.5lf [Torr] +/- %18.5lf [Torr]\n", avg*PA_TO_TORR, error*PA_TO_TORR);
      break;
    case MOLECULAR_DYNAMICS:
      fprintf(FilePtr,"Average Pressure:\n");
      fprintf(FilePtr,"=================\n");
      sum=sum_squared=0.0;
      for(int i=0;i<NR_BLOCKS;i++)
      {
        if(BlockWeightedCount[CurrentSystem][i]>0.0)
        {
          REAL tmp=((StressTensorAccumulated[CurrentSystem][i].ax+StressTensorAccumulated[CurrentSystem][i].by+StressTensorAccumulated[CurrentSystem][i].cz)/
               (3.0*BlockWeightedCount[CurrentSystem][i]))*PRESSURE_CONVERSION_FACTOR;
          sum+=tmp;
          sum_squared+=SQR(tmp);
          fprintf(FilePtr,"\tBlock[%2d] %18.5lf [Pa]\n",i,(double)tmp);
        }
        else
          fprintf(FilePtr,"\tBlock[%2d] %18.5lf [Pa]\n",i,(double)0.0);
      }
      avg=AVERAGE(sum);
      error=ERROR_CONFIDENCE_INTERVAL_95(sum,sum_squared);
      fprintf(FilePtr,"\t------------------------------------------------------------------------------\n");
      fprintf(FilePtr,"\tAverage   %18.5lf [Pa] +/- %18.5lf [Pa]\n", avg, error);
      fprintf(FilePtr,"\tAverage   %18.5lf [bar] +/- %18.5lf [bar]\n", avg*PA_TO_BAR, error*PA_TO_BAR);
      fprintf(FilePtr,"\tAverage   %18.5lf [atm] +/- %18.5lf [atm]\n", avg*PA_TO_ATM, error*PA_TO_ATM);
      fprintf(FilePtr,"\tAverage   %18.5lf [Torr] +/- %18.5lf [Torr]\n", avg*PA_TO_TORR, error*PA_TO_TORR);
      break;
    default:
      break;
  }

  // Volume
  fprintf(FilePtr,"\n");
  fprintf(FilePtr,"Average Volume:\n");
  fprintf(FilePtr,"=================\n");
  sum=sum_squared=0.0;
  for(int i=0;i<NR_BLOCKS;i++)
  {
    if(BlockWeightedCount[CurrentSystem][i]>0.0)
    {
      REAL tmp=(VolumeAccumulated[CurrentSystem][i]/BlockWeightedCount[CurrentSystem][i]);
      sum+=tmp;
      sum_squared+=SQR(tmp);
      fprintf(FilePtr,"\tBlock[%2d] %18.5lf [A^3]\n",i,(double)tmp);
    }
    else
      fprintf(FilePtr,"\tBlock[%2d] %18.5lf [A^3]\n",i,(double)0.0);
  }
  avg=AVERAGE(sum);
  error=ERROR_CONFIDENCE_INTERVAL_95(sum,sum_squared);
  fprintf(FilePtr,"\t------------------------------------------------------------------------------\n");
  fprintf(FilePtr,"\tAverage   %18.5lf [A^3] +/- %18.5lf [A^3]\n", avg, error);

  // Box-lengths
  fprintf(FilePtr,"\n");
  fprintf(FilePtr,"Average Box-lengths:\n");
  fprintf(FilePtr,"====================\n");
  sum=sum_squared=0.0;
  for(int i=0;i<NR_BLOCKS;i++)
  {
    if(BlockWeightedCount[CurrentSystem][i]>0.0)
    {
      REAL tmp=(BoxAccumulated[CurrentSystem][i].x/BlockWeightedCount[CurrentSystem][i]);
      sum+=tmp;
      sum_squared+=SQR(tmp);
      fprintf(FilePtr,"\tBlock[%2d] %18.5lf [A^3]\n",i,(double)tmp);
    }
    else
      fprintf(FilePtr,"\tBlock[%2d] %18.5lf [A^3]\n",i,(double)0.0);
  }
  fprintf(FilePtr,"\t------------------------------------------------------------------------------\n");
  avg=AVERAGE(sum);
  error=ERROR_CONFIDENCE_INTERVAL_95(sum,sum_squared);
  fprintf(FilePtr,"\tAverage Box.ax  %18.5lf [A^3] +/- %18.5lf [A^3]\n\n", avg, error);

  sum=sum_squared=0.0;
  for(int i=0;i<NR_BLOCKS;i++)
  {
    if(BlockWeightedCount[CurrentSystem][i]>0.0)
    {
      REAL tmp=(BoxAccumulated[CurrentSystem][i].y/BlockWeightedCount[CurrentSystem][i]);
      sum+=tmp;
      sum_squared+=SQR(tmp);
      fprintf(FilePtr,"\tBlock[%2d] %18.5lf [A^3]\n",i,(double)tmp);
    }
    else
      fprintf(FilePtr,"\tBlock[%2d] %18.5lf [A^3]\n",i,(double)0.0);
  }
  avg=AVERAGE(sum);
  error=ERROR_CONFIDENCE_INTERVAL_95(sum,sum_squared);
  fprintf(FilePtr,"\t------------------------------------------------------------------------------\n");
  fprintf(FilePtr,"\tAverage Box.by  %18.5lf [A^3] +/- %18.5lf [A^3]\n\n", avg, error);

  sum=sum_squared=0.0;
  for(int i=0;i<NR_BLOCKS;i++)
  {
    if(BlockWeightedCount[CurrentSystem][i]>0.0)
    {
      REAL tmp=(BoxAccumulated[CurrentSystem][i].z/BlockWeightedCount[CurrentSystem][i]);
      sum+=tmp;
      sum_squared+=SQR(tmp);
      fprintf(FilePtr,"\tBlock[%2d] %18.5lf [A^3]\n",i,(double)tmp);
    }
    else
      fprintf(FilePtr,"\tBlock[%2d] %18.5lf [A^3]\n",i,(double)0.0);
  }
  avg=AVERAGE(sum);
  error=ERROR_CONFIDENCE_INTERVAL_95(sum,sum_squared);
  fprintf(FilePtr,"\t------------------------------------------------------------------------------\n");
  fprintf(FilePtr,"\tAverage Box.cz  %18.5lf [A^3] +/- %18.5lf [A^3]\n\n", avg, error);

  sum=sum_squared=0.0;
  for(int i=0;i<NR_BLOCKS;i++)
  {
    if(BlockWeightedCount[CurrentSystem][i]>0.0)
    {
      REAL tmp=(AlphaAngleAccumulated[CurrentSystem][i]/BlockWeightedCount[CurrentSystem][i]);
      sum+=tmp;
      sum_squared+=SQR(tmp);
      fprintf(FilePtr,"\tBlock[%2d] %18.5lf [A^3]\n",i,(double)tmp);
    }
    else
      fprintf(FilePtr,"\tBlock[%2d] %18.5lf [A^3]\n",i,(double)0.0);
  }
  fprintf(FilePtr,"\t------------------------------------------------------------------------------\n");
  avg=AVERAGE(sum);
  error=ERROR_CONFIDENCE_INTERVAL_95(sum,sum_squared);
  fprintf(FilePtr,"\tAverage alpha angle  %18.5lf [degrees] +/- %18.5lf [degrees]\n\n", avg, error);


  sum=sum_squared=0.0;
  for(int i=0;i<NR_BLOCKS;i++)
  {
    if(BlockWeightedCount[CurrentSystem][i]>0.0)
    {
      REAL tmp=(BetaAngleAccumulated[CurrentSystem][i]/BlockWeightedCount[CurrentSystem][i]);
      sum+=tmp;
      sum_squared+=SQR(tmp);
      fprintf(FilePtr,"\tBlock[%2d] %18.5lf [A^3]\n",i,(double)tmp);
    }
    else
      fprintf(FilePtr,"\tBlock[%2d] %18.5lf [A^3]\n",i,(double)0.0);
  }
  avg=AVERAGE(sum);
  error=ERROR_CONFIDENCE_INTERVAL_95(sum,sum_squared);
  fprintf(FilePtr,"\t------------------------------------------------------------------------------\n");
  fprintf(FilePtr,"\tAverage beta angle  %18.5lf [degrees] +/- %18.5lf [degrees]\n\n", avg, error);


  sum=sum_squared=0.0;
  for(int i=0;i<NR_BLOCKS;i++)
  {
    if(BlockWeightedCount[CurrentSystem][i]>0.0)
    {
      REAL tmp=(GammaAngleAccumulated[CurrentSystem][i]/BlockWeightedCount[CurrentSystem][i]);
      sum+=tmp;
      sum_squared+=SQR(tmp);
      fprintf(FilePtr,"\tBlock[%2d] %18.5lf [A^3]\n",i,(double)tmp);
    }
    else
      fprintf(FilePtr,"\tBlock[%2d] %18.5lf [A^3]\n",i,(double)0.0);
  }
  avg=AVERAGE(sum);
  error=ERROR_CONFIDENCE_INTERVAL_95(sum,sum_squared);
  fprintf(FilePtr,"\t------------------------------------------------------------------------------\n");
  fprintf(FilePtr,"\tAverage gamma angle  %18.5lf [degrees] +/- %18.5lf [degrees]\n\n", avg, error);


  // Average Surface area
  fprintf(FilePtr,"\n");
  fprintf(FilePtr,"Average Surface Area:\n");
  fprintf(FilePtr,"=====================\n");
  for(int j=0;j<NumberOfComponents;j++)
  {
    sum=sum_squared=0.0;
    for(int i=0;i<NR_BLOCKS;i++)
    {
      if(SurfaceAreaCount[CurrentSystem][i]>0.0)
      {
        REAL tmp=SurfaceAreaFrameworkAccumulated[CurrentSystem][i]/SurfaceAreaCount[CurrentSystem][i];
        sum+=tmp;
        sum_squared+=SQR(tmp);
        fprintf(FilePtr,"\tBlock[%2d] %-lf [-]\n",i,(double)tmp);
      }
      else
        fprintf(FilePtr,"\tBlock[%2d] %-lf [-]\n",i,(double)0.0);
    }
    avg=AVERAGE(sum);
    error=ERROR_CONFIDENCE_INTERVAL_95(sum,sum_squared);
    fprintf(FilePtr,"\t------------------------------------------------------------------------------\n");
    fprintf(FilePtr,"\tAverage surface area:   %lf +/- %lf [A^2]\n", avg, error);
    REAL convert_to_m2_g = SQR(ANGSTROM)*AVOGADRO_CONSTANT/(Framework[CurrentSystem].FrameworkMass);
    fprintf(FilePtr,"\t                        %lf +/- %lf [m^2/g]\n", avg * convert_to_m2_g, error * convert_to_m2_g);
    REAL convert_to_m2_cm3 = 1.0e4/Volume[CurrentSystem];
    fprintf(FilePtr,"\t                        %lf +/- %lf [m^2/cm^3]\n\n", avg * convert_to_m2_cm3, error * convert_to_m2_cm3);
  }

  // Density
  fprintf(FilePtr,"\n");
  fprintf(FilePtr,"Average density:\n");
  fprintf(FilePtr,"=================\n");
  sum=sum_squared=0.0;
  for(int i=0;i<NR_BLOCKS;i++)
  {
    if(BlockWeightedCount[CurrentSystem][i]>0.0)
    {
      REAL tmp=(DensityAccumulated[CurrentSystem][i]/BlockWeightedCount[CurrentSystem][i])*DENSITY_CONVERSION_FACTOR;
      sum+=tmp;
      sum_squared+=SQR(tmp);
      fprintf(FilePtr,"\tBlock[%2d] %18.5lf [kg/m^3]\n",i,(double)tmp);
    }
    else
      fprintf(FilePtr,"\tBlock[%2d] %18.5lf [kg/m^3]\n",i,(double)0.0);
  }
  avg=AVERAGE(sum);
  error=ERROR_CONFIDENCE_INTERVAL_95(sum,sum_squared);
  fprintf(FilePtr,"\t------------------------------------------------------------------------------\n");
  fprintf(FilePtr,"\tAverage   %18.5lf [kg/m^3] +/- %18.5lf [kg/m^3]\n", avg, error);
  fprintf(FilePtr,"\n");

  for(int j=0;j<NumberOfComponents;j++)
  {
    fprintf(FilePtr,"\tAverage density component %d [%s]\n",j,Components[j].Name);
    fprintf(FilePtr,"\t-------------------------------------------------------------\n");

    sum=sum_squared=0.0;
    for(int i=0;i<NR_BLOCKS;i++)
    {
      if(BlockWeightedCount[CurrentSystem][i]>0.0)
      {
        REAL tmp=(DensityPerComponentAccumulated[CurrentSystem][j][i]/BlockWeightedCount[CurrentSystem][i])*DENSITY_CONVERSION_FACTOR;
        sum+=tmp;
        sum_squared+=SQR(tmp);

        fprintf(FilePtr,"\t\tBlock[%2d] %18.5lf [kg/m^3]\n",i,(double)tmp);
      }
      else
        fprintf(FilePtr,"\t\tBlock[%2d] %-18.5lf [-]\n",i,(double)0.0);
    }
    avg=AVERAGE(sum);
    error=ERROR_CONFIDENCE_INTERVAL_95(sum,sum_squared);
    fprintf(FilePtr,"\t\t------------------------------------------------------------------------------\n");
    fprintf(FilePtr,"\t\tAverage   %18.5lf [kg/m^3] +/- %18.5lf [kg/m^3]\n", avg, error);
  }

  // compressibility
  fprintf(FilePtr,"\n");
  fprintf(FilePtr,"Average compressibility Z:\n");
  fprintf(FilePtr,"=========================\n");
  sum=sum_squared=0.0;
  for(int i=0;i<NR_BLOCKS;i++)
  {
    if(BlockWeightedCount[CurrentSystem][i]>0.0)
    {
      REAL tmp=CompressibilityAccumulated[CurrentSystem][i]/BlockWeightedCount[CurrentSystem][i];
      sum+=tmp;
      sum_squared+=SQR(tmp);
      fprintf(FilePtr,"\tBlock[%2d] %18.5lf [-]\n",i,(double)tmp);
    }
    else
      fprintf(FilePtr,"\tBlock[%2d] %18.5lf [-]\n",i,(double)0.0);
  }
  avg=AVERAGE(sum);
  error=ERROR_CONFIDENCE_INTERVAL_95(sum,sum_squared);
  fprintf(FilePtr,"\t------------------------------------------------------------------------------\n");
  fprintf(FilePtr,"\tAverage   %18.5lf [-] +/- %18.5lf [-]\n", avg, error);

  if(NumberOfSystems==2)
  {
    // Heat of vaporization
    fprintf(FilePtr,"\n");
    fprintf(FilePtr,"Heat of vaporization:\n");
    fprintf(FilePtr,"=====================\n");
    sum=sum_squared=0.0;
    REAL N=NumberOfAdsorbateMolecules[CurrentSystem]+NumberOfCationMolecules[CurrentSystem];
    if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
      N+=Framework[CurrentSystem].TotalNumberOfAtoms;
    for(int i=0;i<NR_BLOCKS;i++)
    {
      if(BlockWeightedCount[CurrentSystem][i]>0.0)
      {
        REAL tmp=ENERGY_TO_KELVIN*((fabs(EnergyPerMolecule[CurrentSystem][i]-EnergyPerMolecule[1-CurrentSystem][i])/BlockWeightedCount[CurrentSystem][i])
                 +((MolecularStressTensorAccumulated[CurrentSystem][i].ax+MolecularStressTensorAccumulated[CurrentSystem][i].by+
                 MolecularStressTensorAccumulated[CurrentSystem][i].cz)/(3.0*BlockWeightedCount[CurrentSystem][i]))*
                 (fabs(VolumePerMolecule[CurrentSystem][i]-VolumePerMolecule[1-CurrentSystem][i])/BlockWeightedCount[CurrentSystem][i]));
        sum+=tmp;
        sum_squared+=SQR(tmp);

        fprintf(FilePtr,"\tBlock[%2d] %-18.5lf [K]\n",i,(double)tmp);
      }
      else
        fprintf(FilePtr,"\tBlock[%2d] %-18.5lf [K]\n",i,(double)0.0);
    }
    avg=AVERAGE(sum);
    error=ERROR_CONFIDENCE_INTERVAL_95(sum,sum_squared);
    fprintf(FilePtr,"\t------------------------------------------------------------------------------\n");
    fprintf(FilePtr,"\tAverage   %18.5lf +/- %18lf [K]\n", avg, error);
    fprintf(FilePtr,"\tAverage   %18.5lf +/- %18lf [kJ/mol]\n", avg * KELVIN_TO_KJ_PER_MOL, error * KELVIN_TO_KJ_PER_MOL);
    fprintf(FilePtr,"\tAverage   %18.5lf +/- %18lf [kcal/mol]\n", avg * KELVIN_TO_KCAL_PER_MOL, error * KELVIN_TO_KCAL_PER_MOL);


    fprintf(FilePtr,"\n");
    fprintf(FilePtr,"Heat of vaporization term (U_v/N_v-U_l/N_l):\n");
    fprintf(FilePtr,"============================================\n");
    sum=sum_squared=0.0;
    for(int i=0;i<NR_BLOCKS;i++)
    {
      if(BlockWeightedCount[CurrentSystem][i]>0.0)
      {
        REAL tmp=ENERGY_TO_KELVIN*fabs(EnergyPerMolecule[CurrentSystem][i]-EnergyPerMolecule[1-CurrentSystem][i])/BlockWeightedCount[CurrentSystem][i];
        sum+=tmp;
        sum_squared+=SQR(tmp);

        fprintf(FilePtr,"\tBlock[%2d] %-18.5lf [K]\n",i,(double)tmp);
      }
      else
        fprintf(FilePtr,"\tBlock[%2d] %-18.5lf [K]\n",i,(double)0.0);
    }
    avg=AVERAGE(sum);
    error=ERROR_CONFIDENCE_INTERVAL_95(sum,sum_squared);
    fprintf(FilePtr,"\t------------------------------------------------------------------------------\n");
    fprintf(FilePtr,"\tAverage   %18.5lf +/- %18lf [K]\n", avg, error);
    fprintf(FilePtr,"\tAverage   %18.5lf +/- %18lf [kJ/mol]\n", avg * KELVIN_TO_KJ_PER_MOL, error * KELVIN_TO_KJ_PER_MOL);
    fprintf(FilePtr,"\tAverage   %18.5lf +/- %18lf [kcal/mol]\n", avg * KELVIN_TO_KCAL_PER_MOL, error * KELVIN_TO_KCAL_PER_MOL);

    fprintf(FilePtr,"\n");
    fprintf(FilePtr,"Heat of vaporization term p (V_v-V_l), [is approximately RT]:\n");
    fprintf(FilePtr,"=============================================================\n");
    sum=sum_squared=0.0;
    N=NumberOfAdsorbateMolecules[CurrentSystem]+NumberOfCationMolecules[CurrentSystem];
    if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
      N+=Framework[CurrentSystem].TotalNumberOfAtoms;
    for(int i=0;i<NR_BLOCKS;i++)
    {
      if(BlockWeightedCount[CurrentSystem][i]>0.0)
      {
        REAL tmp=ENERGY_TO_KELVIN*fabs(((MolecularStressTensorAccumulated[CurrentSystem][i].ax+MolecularStressTensorAccumulated[CurrentSystem][i].by+
                 MolecularStressTensorAccumulated[CurrentSystem][i].cz)/(3.0*BlockWeightedCount[CurrentSystem][i]))*
                 (VolumePerMolecule[CurrentSystem][i]-VolumePerMolecule[1-CurrentSystem][i])/BlockWeightedCount[CurrentSystem][i]);
        sum+=tmp;
        sum_squared+=SQR(tmp);

        fprintf(FilePtr,"\tBlock[%2d] %-18.5lf [K]\n",i,(double)tmp);
      }
      else
        fprintf(FilePtr,"\tBlock[%2d] %-18.5lf [K]\n",i,(double)0.0);
    }
    avg=AVERAGE(sum);
    error=ERROR_CONFIDENCE_INTERVAL_95(sum,sum_squared);
    fprintf(FilePtr,"\t------------------------------------------------------------------------------\n");
    fprintf(FilePtr,"\tAverage   %18.5lf +/- %18lf [K]\n", avg, error);
    fprintf(FilePtr,"\tAverage   %18.5lf +/- %18lf [kJ/mol]\n", avg * KELVIN_TO_KJ_PER_MOL, error * KELVIN_TO_KJ_PER_MOL);
    fprintf(FilePtr,"\tAverage   %18.5lf +/- %18lf [kcal/mol]\n", avg * KELVIN_TO_KCAL_PER_MOL, error * KELVIN_TO_KCAL_PER_MOL);
  }

  // heat capacity
  fprintf(FilePtr,"\n");
  fprintf(FilePtr,"Average Heat Capacity (MC-NPT-ensemble): [1/(kB T^2)]*[<H^2>-<H>^2]\n");
  fprintf(FilePtr,"===================================================================\n");
  sum=sum_squared=0.0;
  for(int i=0;i<NR_BLOCKS;i++)
  {
    if(BlockWeightedCount[CurrentSystem][i]>0.0)
    {
      REAL H2=EnthalpySquaredAccumulated[CurrentSystem][i]/BlockWeightedCount[CurrentSystem][i];
      REAL H=EnthalpyAccumulated[CurrentSystem][i]/BlockWeightedCount[CurrentSystem][i];
      REAL T=therm_baro_stats.ExternalTemperature[CurrentSystem];
      REAL tmp=HEAT_CAPACITY_CONVERSION_FACTOR*((H2-SQR(H))/(K_B*SQR(T)));
      sum+=tmp;
      sum_squared+=SQR(tmp);
      fprintf(FilePtr,"\tBlock[%2d] %lf [J/mol/K]\n",i,(double)tmp);
    }
    else
      fprintf(FilePtr,"\tBlock[%2d] %lf [J/mol/K]\n",i,(double)0.0);
  }
  avg=AVERAGE(sum);
  error=ERROR_CONFIDENCE_INTERVAL_95(sum,sum_squared);
  fprintf(FilePtr,"\t------------------------------------------------------------------------------\n");
  fprintf(FilePtr,"\tAverage   %18.5lf [J/mol/K] +/- %18.5lf [J/mol/K]\n", avg, error);
  fprintf(FilePtr,"\tAverage   %18.5lf [cal/mol/K] +/- %18.5lf [cal/mol/K]\n", avg * J_TO_CAL, error * J_TO_CAL);

  // enthalpy of desorption
  fprintf(FilePtr,"\n");
  fprintf(FilePtr,"Enthalpy of adsorption:\n");
  fprintf(FilePtr,"=======================\n\n");
  sum=sum_squared=0.0;

  if(NumberOfComponents==1)
  {
    fprintf(FilePtr,"\tTotal enthalpy of adsorption\n");
    fprintf(FilePtr,"\t----------------------------\n");
   
    sum=sum_squared=0.0;
    for(int i=0;i<NR_BLOCKS;i++)
    {
      if(BlockWeightedCount[CurrentSystem][i]>0.0)
      {
        REAL UTotalTimesN = TotalEnergyTimesNumberOfMoleculesAccumulated[CurrentSystem][i]/BlockWeightedCount[CurrentSystem][i];
        REAL UTotal = UTotalAccumulated[CurrentSystem][i]/BlockWeightedCount[CurrentSystem][i];
        REAL N = NumberOfIntegerMoleculesAccumulated[CurrentSystem][i]/BlockWeightedCount[CurrentSystem][i];
        REAL NSquared = NumberOfMoleculesSquaredAccumulated[CurrentSystem][i]/BlockWeightedCount[CurrentSystem][i];

        REAL tmp=ENERGY_TO_KELVIN*((UTotalTimesN - UTotal*N)/(NSquared - SQR(N)))-therm_baro_stats.ExternalTemperature[CurrentSystem];
        sum+=tmp;
        sum_squared+=SQR(tmp);
        fprintf(FilePtr,"\tBlock[%2d] %-18.5lf [K]\n",i,(double)tmp);
      }
      else
        fprintf(FilePtr,"\tBlock[%2d] %-18.5lf [K]\n",i,(double)0.0);
    }
    avg=AVERAGE(sum);
    error=ERROR_CONFIDENCE_INTERVAL_95(sum,sum_squared);
    fprintf(FilePtr,"\t------------------------------------------------------------------------------\n");
    fprintf(FilePtr,"\tAverage   %18.5lf +/- %18lf [K]\n", avg, error);
    fprintf(FilePtr,"\t          %18.5lf +/- %18lf [KJ/MOL]\n", avg * KELVIN_TO_KJ_PER_MOL, error * KELVIN_TO_KJ_PER_MOL);
    fprintf(FilePtr,"\tNote: Ug should be subtracted from this value\n");
    fprintf(FilePtr,"\tNote: The heat of adsorption Q=-H\n\n");
  }
  
  if(NumberOfComponents>1)
  {
     REAL (*HeatOfAdsorptionPerComponent)[NR_BLOCKS];
     REAL_MATRIX matrix;

    HeatOfAdsorptionPerComponent=(REAL(*)[NR_BLOCKS])calloc(NumberOfComponents,sizeof(REAL[NR_BLOCKS]));

    for(int i=0;i<NR_BLOCKS;i++)
    {
      matrix=CreateRealMatrix(NumberOfComponents,NumberOfComponents);

      for(int k1=0;k1<NumberOfComponents;k1++)
      {
        for(int k2=0;k2<NumberOfComponents;k2++)
        {
          matrix.element[k1][k2]=NumberOfMoleculesPerComponentSquaredAccumulated[CurrentSystem][k1][k2][i]/BlockWeightedCount[CurrentSystem][i]-
                              ( NumberOfIntegerMoleculesPerComponentAccumulated[CurrentSystem][k1][i]/BlockWeightedCount[CurrentSystem][i])*
                              ( NumberOfIntegerMoleculesPerComponentAccumulated[CurrentSystem][k2][i]/BlockWeightedCount[CurrentSystem][i]);
        }
      }

      InverseRealMatrix(matrix);
	 
      for(int k1=0;k1<NumberOfComponents;k1++)
      {
        HeatOfAdsorptionPerComponent[k1][i]=0.0;
        for(int k2=0;k2<NumberOfComponents;k2++)
        {
          REAL UTotalTimesNcomp = TotalEnergyTimesNumberOfMoleculesPerComponentAccumulated[CurrentSystem][k2][i]/BlockWeightedCount[CurrentSystem][i];
          REAL UTotal = UTotalAccumulated[CurrentSystem][i]/BlockWeightedCount[CurrentSystem][i];
          REAL Ncomp = NumberOfIntegerMoleculesPerComponentAccumulated[CurrentSystem][k2][i]/BlockWeightedCount[CurrentSystem][i];
          HeatOfAdsorptionPerComponent[k1][i]+=ENERGY_TO_KELVIN*((UTotalTimesNcomp-UTotal*Ncomp)*matrix.element[k2][k1]);
        }
        HeatOfAdsorptionPerComponent[k1][i]-=therm_baro_stats.ExternalTemperature[CurrentSystem];
      }

      DeleteRealMatrix(matrix);
    }

    for(int k1=0;k1<NumberOfComponents;k1++)
    {
      sum=sum_squared=0.0;
      fprintf(FilePtr,"\tEnthalpy of adsorption component %d [%s]\n",k1,Components[k1].Name);
      fprintf(FilePtr,"\t-------------------------------------------------------------\n");
      for(int i=0;i<NR_BLOCKS;i++)
      {
        if(BlockWeightedCount[CurrentSystem][i]>0.0)
        {
          REAL tmp=HeatOfAdsorptionPerComponent[k1][i];
          sum+=tmp;
          sum_squared+=SQR(tmp);
          fprintf(FilePtr,"\t\tBlock[%2d] %-18.5lf [-]\n",i,tmp);
        }
        else
          fprintf(FilePtr,"\t\tBlock[%2d] %-18.5lf [K]\n",i,(double)0.0);
      }
      avg=AVERAGE(sum);
      error=ERROR_CONFIDENCE_INTERVAL_95(sum,sum_squared);
      fprintf(FilePtr,"\t\t------------------------------------------------------------------------------\n");
      fprintf(FilePtr,"\t\tAverage   %18.5lf +/- %18lf [K]\n", avg, error);
      fprintf(FilePtr,"\t\t          %18.5lf +/- %18lf [KJ/MOL]\n", avg * KELVIN_TO_KJ_PER_MOL, error * KELVIN_TO_KJ_PER_MOL);
      fprintf(FilePtr,"\t\tNote: Ug should be subtracted to this value\n");
      fprintf(FilePtr,"\t\tNote: The heat of adsorption Q=-H\n\n");
    }

    sum=sum_squared=0.0;
    fprintf(FilePtr,"\tTotal enthalpy of adsorption from components and measured mol-fraction\n");
    fprintf(FilePtr,"\t----------------------------------------------------------------------\n");
    for(int i=0;i<NR_BLOCKS;i++)
    {
      if(BlockWeightedCount[CurrentSystem][i]>0.0)
      {
        double loading_total = 0.0;
        for (int j=0;j<NumberOfComponents;j++)
           loading_total += NumberOfIntegerMoleculesPerComponentAccumulated[CurrentSystem][j][i];
        double sumc=0.0;
        for(int k1=0;k1<NumberOfComponents;k1++)
        {
          double loading_i = NumberOfIntegerMoleculesPerComponentAccumulated[CurrentSystem][k1][i];
          sumc+=(loading_i/loading_total)*HeatOfAdsorptionPerComponent[k1][i];
        }
        sum+=sumc;
        sum_squared+=SQR(sumc);
        fprintf(FilePtr,"\t\tBlock[%2d] %-18.5lf [-]\n",i,(double)sumc);
      }
      else
        fprintf(FilePtr,"\t\tBlock[%2d] %-18.5lf [K]\n",i,(double)0.0);
    }
    avg=AVERAGE(sum);
    error=ERROR_CONFIDENCE_INTERVAL_95(sum,sum_squared);
    fprintf(FilePtr,"\t\t------------------------------------------------------------------------------\n");
    fprintf(FilePtr,"\t\tAverage   %18.5lf +/- %18lf [K]\n", avg, error);
    fprintf(FilePtr,"\t\t          %18.5lf +/- %18lf [KJ/MOL]\n", avg * KELVIN_TO_KJ_PER_MOL, error * KELVIN_TO_KJ_PER_MOL);
    fprintf(FilePtr,"\t\tNote: Ug should be subtracted to this value\n");
    fprintf(FilePtr,"\t\tNote: The heat of adsorption Q=-H\n\n");

    free(HeatOfAdsorptionPerComponent);
  }


  // Derivative
  fprintf(FilePtr,"\n");
  fprintf(FilePtr,"derivative of the chemical potential with respect to density (constant T,V):\n");
  fprintf(FilePtr,"============================================================================\n");
  sum=sum_squared=0.0;
  for(int i=0;i<NR_BLOCKS;i++)
  {
    if(BlockWeightedCount[CurrentSystem][i]>0.0)
    {
      REAL tmp=ENERGY_TO_KELVIN*(Volume[0]*therm_baro_stats.ExternalTemperature[CurrentSystem]/
               (NumberOfMoleculesSquaredAccumulated[CurrentSystem][i]/BlockWeightedCount[CurrentSystem][i]-
               SQR(NumberOfIntegerMoleculesAccumulated[CurrentSystem][i]/BlockWeightedCount[CurrentSystem][i])));
      sum+=tmp;
      sum_squared+=SQR(tmp);
      fprintf(FilePtr,"\tBlock[%2d] %-18.5lf [-]\n",i,(double)tmp);
    }
    else
      fprintf(FilePtr,"\tBlock[%2d] %-18.5lf [-]\n",i,(double)0.0);
  }
  avg=AVERAGE(sum);
  error=ERROR_CONFIDENCE_INTERVAL_95(sum,sum_squared);
  fprintf(FilePtr,"\t------------------------------------------------------------------------------\n");
  fprintf(FilePtr,"\tAverage   %18.5lf +/- %18lf [-]\n", avg, error);

  // output final averages of the principle moments of inertia per component
  if(ComputePrincipleMomentsOfInertia)
  {
    VECTOR sumv1,sumv2,av;

    fprintf(FilePtr,"\n");
    fprintf(FilePtr,"Average Principle Moments of Inertia:\n");
    fprintf(FilePtr,"=====================================\n");
    for(int j=0;j<NumberOfComponents;j++)
    {
      sumv1.x=0.0;
      sumv1.y=0.0;
      sumv1.z=0.0;
      sumv2.x=0.0;
      sumv2.y=0.0;
      sumv2.z=0.0;
      for(int i=0;i<NR_BLOCKS;i++)
      {
        if(PrincipleMomentsOfInertiaCount[CurrentSystem][j][i]>0.0)
        {
          av.x=PrincipleMomentsOfInertiaAccumulated[CurrentSystem][j][i].x/PrincipleMomentsOfInertiaCount[CurrentSystem][j][i];
          av.y=PrincipleMomentsOfInertiaAccumulated[CurrentSystem][j][i].y/PrincipleMomentsOfInertiaCount[CurrentSystem][j][i];
          av.z=PrincipleMomentsOfInertiaAccumulated[CurrentSystem][j][i].z/PrincipleMomentsOfInertiaCount[CurrentSystem][j][i];
          sumv1.x+=av.x;
          sumv1.y+=av.y;
          sumv1.z+=av.z;
          sumv2.x+=SQR(av.x);
          sumv2.y+=SQR(av.y);
          sumv2.z+=SQR(av.z);
          fprintf(FilePtr,"\tBlock[%2d] %-lg %-lg %-lg [-]\n",i,(double)av.x,(double)av.y,(double)av.z);
        }
        else
          fprintf(FilePtr,"\tBlock[%2d] %-lg [-]\n",i,(double)0.0);
      }
      REAL avgx=AVERAGE(sumv1.x);
      REAL avgy=AVERAGE(sumv1.y);
      REAL avgz=AVERAGE(sumv1.z);
      REAL errorx=ERROR_CONFIDENCE_INTERVAL_95(sumv1.x,sumv2.x);
      REAL errory=ERROR_CONFIDENCE_INTERVAL_95(sumv1.y,sumv2.y);
      REAL errorz=ERROR_CONFIDENCE_INTERVAL_95(sumv1.z,sumv2.z);
      fprintf(FilePtr,"\t------------------------------------------------------------------------------\n");
      fprintf(FilePtr,"\t[%s] Average Principle Moment of Inertia:  %lg +/- %lg [-], %lg +/- %lg [-], %lg +/- %lg [-]\n\n",
        Components[j].Name, avgx, errorx, avgy, errory, avgz, errorz);
    }
  }


  fprintf(FilePtr,"\n\n\n\n\n");
  fprintf(FilePtr,"Average energies of the system[%d]:\n",CurrentSystem);
  fprintf(FilePtr,"========================================================================\n");

  // Host intra molecular energies
  PrintProperty(FilePtr,"Average Host Bond stretch energy:\n=================================\n",
                "[K]",ENERGY_TO_KELVIN,UHostBondAccumulated);

  PrintProperty(FilePtr,"Average Host UreyBradley stretch energy:\n==============================================\n",
                "[K]",ENERGY_TO_KELVIN,UHostUreyBradleyAccumulated);

  PrintProperty(FilePtr,"Average Host Bend angle energy:\n===============================\n",
                "[K]",ENERGY_TO_KELVIN,UHostBendAccumulated);

  PrintProperty(FilePtr,"Average Host Bend angle inversion energy:\n=========================================\n",
                "[K]",ENERGY_TO_KELVIN,UHostInversionBendAccumulated);

  PrintProperty(FilePtr,"Average Host Torsion energy:\n============================\n",
                "[K]",ENERGY_TO_KELVIN,UHostTorsionAccumulated);

  PrintProperty(FilePtr,"Average Host Improper Torsion energy:\n=====================================\n",
                "[K]",ENERGY_TO_KELVIN,UHostImproperTorsionAccumulated);

  PrintProperty(FilePtr,"Average Host Bond-Bond cross term energy:\n===============================================\n",
                "[K]",ENERGY_TO_KELVIN,UHostBondBondAccumulated);

  PrintProperty(FilePtr,"Average Host Bend-Bend cross term energy:\n=========================================\n",
                "[K]",ENERGY_TO_KELVIN,UHostBendBendAccumulated);

  PrintProperty(FilePtr,"Average Host Bond-Bend cross term energy:\n============================================\n",
                "[K]",ENERGY_TO_KELVIN,UHostBondBendAccumulated);

  PrintProperty(FilePtr,"Average Host Bond-Torsion cross term energy:\n============================================\n",
                "[K]",ENERGY_TO_KELVIN,UHostBondTorsionAccumulated);

  PrintProperty(FilePtr,"Average Host Bend-Torsion cross term energy:\n============================================\n",
                "[K]",ENERGY_TO_KELVIN,UHostBendTorsionAccumulated);

  // Adsorbate intra molecular energies
  PrintProperty(FilePtr,"Average Adsorbate Bond stretch energy:\n====================================\n",
                "[K]",ENERGY_TO_KELVIN,UAdsorbateBondAccumulated);

  PrintProperty(FilePtr,"Average Adsorbate UreyBradley stretch energy:\n==============================================\n",
                "[K]",ENERGY_TO_KELVIN,UAdsorbateUreyBradleyAccumulated);

  PrintProperty(FilePtr,"Average Adsorbate Bend angle energy:\n====================================\n",
                "[K]",ENERGY_TO_KELVIN,UAdsorbateBendAccumulated);

  PrintProperty(FilePtr,"Average Adsorbate Bend angle inversion energy:\n==============================================\n",
                "[K]",ENERGY_TO_KELVIN,UAdsorbateInversionBendAccumulated);

  PrintProperty(FilePtr,"Average Adsorbate Torsion energy:\n=================================\n",
                "[K]",ENERGY_TO_KELVIN,UAdsorbateTorsionAccumulated);

  PrintProperty(FilePtr,"Average Adsorbate Improper Torsion energy:\n==========================================\n",
                "[K]",ENERGY_TO_KELVIN,UAdsorbateImproperTorsionAccumulated);

  PrintProperty(FilePtr,"Average Adsorbate Bond-Bond cross term energy:\n====================================================\n",
                "[K]",ENERGY_TO_KELVIN,UAdsorbateBondBondAccumulated);

  PrintProperty(FilePtr,"Average Adsorbate Bend-Bend cross term energy:\n==============================================\n",
                "[K]",ENERGY_TO_KELVIN,UAdsorbateBendBendAccumulated);

  PrintProperty(FilePtr,"Average Adsorbate Bond-Bend cross term energy:\n=================================================\n",
                "[K]",ENERGY_TO_KELVIN,UAdsorbateBondBendAccumulated);

  PrintProperty(FilePtr,"Average Adsorbate Bond-Torsion cross term energy:\n=================================================\n",
                "[K]",ENERGY_TO_KELVIN,UAdsorbateBondTorsionAccumulated);

  PrintProperty(FilePtr,"Average Adsorbate Bend-Torsion cross term energy:\n=================================================\n",
                "[K]",ENERGY_TO_KELVIN,UAdsorbateBendTorsionAccumulated);

  PrintProperty(FilePtr,"Average Adsorbate Intra Van der Waals energy:\n=============================================\n",
                "[K]",ENERGY_TO_KELVIN,UAdsorbateIntraVDWAccumulated);

  PrintProperty(FilePtr,"Average Adsorbate Intra charge-charge Coulomb energy:\n=======================================\n",
                "[K]",ENERGY_TO_KELVIN,UAdsorbateIntraChargeChargeAccumulated);

  PrintProperty(FilePtr,"Average Adsorbate Intra charge-bonddipole Coulomb energy:\n=======================================\n",
                "[K]",ENERGY_TO_KELVIN,UAdsorbateIntraChargeBondDipoleAccumulated);

  PrintProperty(FilePtr,"Average Adsorbate Intra bonddipole-bonddipole Coulomb energy:\n=======================================\n",
                "[K]",ENERGY_TO_KELVIN,UAdsorbateIntraBondDipoleBondDipoleAccumulated);


  // Cation intra molecular energies
  PrintProperty(FilePtr,"Average Cation Bond stretch energy:\n=================================\n",
                "[K]",ENERGY_TO_KELVIN,UCationBondAccumulated);

  PrintProperty(FilePtr,"Average Cation UreyBradley stretch energy:\n==========================================\n",
                "[K]",ENERGY_TO_KELVIN,UCationUreyBradleyAccumulated);

  PrintProperty(FilePtr,"Average Cation Bend angle energy:\n=================================\n",
                "[K]",ENERGY_TO_KELVIN,UCationBendAccumulated);

  PrintProperty(FilePtr,"Average Cation Bend angle inversion energy:\n===========================================\n",
                "[K]",ENERGY_TO_KELVIN,UCationInversionBendAccumulated);

  PrintProperty(FilePtr,"Average Cation Torsion energy:\n==============================\n",
                "[K]",ENERGY_TO_KELVIN,UCationTorsionAccumulated);

  PrintProperty(FilePtr,"Average Cation Improper Torsion energy:\n=======================================\n",
                "[K]",ENERGY_TO_KELVIN,UCationImproperTorsionAccumulated);

  PrintProperty(FilePtr,"Average Cation Bond-Bond cross term energy:\n=================================================\n",
                "[K]",ENERGY_TO_KELVIN,UCationBondBondAccumulated);

  PrintProperty(FilePtr,"Average Cation Bend-Bend cross term energy:\n===========================================\n",
                "[K]",ENERGY_TO_KELVIN,UCationBendBendAccumulated);

  PrintProperty(FilePtr,"Average Cation Bond-Bend cross term energy:\n==============================================\n",
                "[K]",ENERGY_TO_KELVIN,UCationBondBendAccumulated);

  PrintProperty(FilePtr,"Average Cation Bond-Torsion cross term energy:\n==============================================\n",
                "[K]",ENERGY_TO_KELVIN,UCationBondTorsionAccumulated);

  PrintProperty(FilePtr,"Average Cation Bend-Torsion cross term energy:\n==============================================\n",
                "[K]",ENERGY_TO_KELVIN,UCationBendTorsionAccumulated);

  PrintProperty(FilePtr,"Average Cation Intra Van der Waals energy:\n==========================================\n",
                "[K]",ENERGY_TO_KELVIN,UCationIntraVDWAccumulated);

  PrintProperty(FilePtr,"Average Cation Intra charge-charge Coulomb energy:\n====================================\n",
                "[K]",ENERGY_TO_KELVIN,UCationIntraChargeChargeAccumulated);

  PrintProperty(FilePtr,"Average Cation Intra charge-bonddipole Coulomb energy:\n====================================\n",
                "[K]",ENERGY_TO_KELVIN,UCationIntraChargeBondDipoleAccumulated);

  PrintProperty(FilePtr,"Average Cation Intra bonddipole-bonddipole Coulomb energy:\n====================================\n",
                "[K]",ENERGY_TO_KELVIN,UCationIntraBondDipoleBondDipoleAccumulated);

  // Host-Host inter molecular energies
  PrintEnergies(FilePtr,"Average Host-Host energy:\n=========================\n","[K]",ENERGY_TO_KELVIN,
                UHostHostAccumulated,UHostHostVDWAccumulated,UHostHostCoulombAccumulated);

  // Adsorbate-Adsorbate inter molecular energies
  PrintEnergies(FilePtr,"Average Adsorbate-Adsorbate energy:\n===================================\n","[K]",ENERGY_TO_KELVIN,
                UAdsorbateAdsorbateAccumulated,UAdsorbateAdsorbateVDWAccumulated,UAdsorbateAdsorbateCoulombAccumulated);

  // Cation-Cation inter molecular energies
  PrintEnergies(FilePtr,"Average Cation-Cation energy:\n=============================\n","[K]",ENERGY_TO_KELVIN,
                UCationCationAccumulated,UCationCationVDWAccumulated,UCationCationCoulombAccumulated);

  // Host-Adsorbate inter molecular energies
  PrintEnergies(FilePtr,"Average Host-Adsorbate energy:\n==============================\n","[K]",ENERGY_TO_KELVIN,
                UHostAdsorbateAccumulated,UHostAdsorbateVDWAccumulated,UHostAdsorbateCoulombAccumulated);

  // Host-Cation inter molecular energies
  PrintEnergies(FilePtr,"Average Host-Cation energy:\n===========================\n","[K]",ENERGY_TO_KELVIN,
                UHostCationAccumulated,UHostCationVDWAccumulated,UHostCationCoulombAccumulated);

  // Adsorbate-Cation inter molecular energies
  PrintEnergies(FilePtr,"Average Adsorbate-Cation energy:\n================================\n","[K]",ENERGY_TO_KELVIN,
                UAdsorbateCationAccumulated,UAdsorbateCationVDWAccumulated,UAdsorbateCationCoulombAccumulated);

  // Host polarization energy
  PrintProperty(FilePtr,"Host polarization energy:\n=======================\n",
                "[K]",ENERGY_TO_KELVIN,UHostPolarizationAccumulated);

  // Adsorbate polarization energy
  PrintProperty(FilePtr,"Adsorbate polarization energy:\n=======================\n",
                "[K]",ENERGY_TO_KELVIN,UAdsorbatePolarizationAccumulated);

  // Cation polarization energy
  PrintProperty(FilePtr,"Cation polarization energy:\n=======================\n",
                "[K]",ENERGY_TO_KELVIN,UCationPolarizationAccumulated);

  // Host back-polarization energy
  PrintProperty(FilePtr,"Host back-polarization energy:\n=======================\n",
                "[K]",ENERGY_TO_KELVIN,UHostBackPolarizationAccumulated);

  // Adsorbate back-polarization energy
  PrintProperty(FilePtr,"Adsorbate back-polarization energy:\n=======================\n",
                "[K]",ENERGY_TO_KELVIN,UAdsorbateBackPolarizationAccumulated);

  // Cation back-polarization energy
  PrintProperty(FilePtr,"Cation back-polarization energy:\n=======================\n",
                "[K]",ENERGY_TO_KELVIN,UCationBackPolarizationAccumulated);


  // Tailcorrection energy
  PrintProperty(FilePtr,"Tail-correction energy:\n=======================\n",
                "[K]",ENERGY_TO_KELVIN,UTailCorrectionAccumulated);

  // Distance constraints energy
  PrintProperty(FilePtr,"Distance-constraints energy:\n=======================\n",
                "[K]",ENERGY_TO_KELVIN,UDistanceConstraintsAccumulated);

  // Angle constraints energy
  PrintProperty(FilePtr,"Angle-constraints energy:\n=======================\n",
                "[K]",ENERGY_TO_KELVIN,UAngleConstraintsAccumulated);

  // Dihedral constraints energy
  PrintProperty(FilePtr,"Dihedral-constraints energy:\n=======================\n",
                "[K]",ENERGY_TO_KELVIN,UDihedralConstraintsAccumulated);

  // Inversion-bend constraints energy
  PrintProperty(FilePtr,"Inversion-bend constraints energy:\n=======================\n",
                "[K]",ENERGY_TO_KELVIN,UInversionBendConstraintsAccumulated);

  // Out-of-plane-distance constraints energy
  PrintProperty(FilePtr,"Out-of-plane-distance constraints energy:\n=======================\n",
                "[K]",ENERGY_TO_KELVIN,UOutOfPlaneDistanceConstraintsAccumulated);

  // Exclusion constraints energy
  PrintProperty(FilePtr,"Exclusion-constraints energy:\n=======================\n",
                "[K]",ENERGY_TO_KELVIN,UExclusionConstraintsAccumulated);

  // Total energy
  PrintProperty(FilePtr,"Total energy:\n=============\n",
                "[K]",ENERGY_TO_KELVIN,UTotalAccumulated);

  // Number of molecules
  fprintf(FilePtr,"\n");
  fprintf(FilePtr,"Number of molecules:\n");
  fprintf(FilePtr,"====================\n\n");

  for(int j=0;j<NumberOfComponents;j++)
  {
    fprintf(FilePtr,"Component %d [%s]\n",j,Components[j].Name);
    fprintf(FilePtr,"-------------------------------------------------------------\n");

    // absolute adsorption
    sum=0.0;
    sum_squared=0.0;
    REAL nr=NumberOfUnitCells[CurrentSystem].x*NumberOfUnitCells[CurrentSystem].y*NumberOfUnitCells[CurrentSystem].z;
    for(int i=0;i<NR_BLOCKS;i++)
    {
      if(BlockWeightedCount[CurrentSystem][i]>0.0)
      {
        REAL tmp=NumberOfIntegerMoleculesPerComponentAccumulated[CurrentSystem][j][i]/(nr * BlockWeightedCount[CurrentSystem][i]);;
        sum+=tmp;
        sum_squared+=SQR(tmp);

        fprintf(FilePtr,"\tBlock[%2d] %-18.5lf [-]\n",i, tmp * nr);
      }
      else
        fprintf(FilePtr,"\tBlock[%2d] %-18.5lf [-]\n",i,(double)0.0);
    }
    avg=AVERAGE(sum);
    error=ERROR_CONFIDENCE_INTERVAL_95(sum,sum_squared);
    fprintf(FilePtr,"\t------------------------------------------------------------------------------\n");
    fprintf(FilePtr,"\tAverage loading absolute                        %18.10lf +/- %18.10lf [-]\n", avg * nr, error * nr);
    fprintf(FilePtr,"\tAverage loading absolute [molecules/unit cell]  %18.10lf +/- %18.10lf [-]\n", avg, error);

    REAL conv_to_mol_per_kg = Components[j].MOLEC_PER_UC_TO_MOL_PER_KG[CurrentSystem];
    fprintf(FilePtr,"\tAverage loading absolute [mol/kg framework]            %18.10lf +/- %18.10lf [-]\n", conv_to_mol_per_kg * avg, conv_to_mol_per_kg * error);

    REAL conv_to_mg_per_g = Components[j].MOLEC_PER_UC_TO_MILLIGRAM_PER_GRAM_OF_FRAMEWORK[CurrentSystem];
    fprintf(FilePtr,"\tAverage loading absolute [milligram/gram framework]    %18.10lf +/- %18.10lf [-]\n", conv_to_mg_per_g * avg, conv_to_mg_per_g * error);

    REAL conv_to_cm3_per_g = Components[j].MOLEC_PER_UC_TO_CC_STP_G[CurrentSystem];
    fprintf(FilePtr,"\tAverage loading absolute [cm^3 (STP)/gr framework]     %18.10lf +/- %18.10lf [-]\n", conv_to_cm3_per_g * avg, conv_to_cm3_per_g * error);

    REAL conv_to_cm3_per_cm3 = Components[j].MOLEC_PER_UC_TO_CC_STP_CC[CurrentSystem];
    fprintf(FilePtr,"\tAverage loading absolute [cm^3 (STP)/cm^3 framework]   %18.10lf +/- %18.10lf [-]\n", conv_to_cm3_per_cm3 * avg, conv_to_cm3_per_cm3 * error);
    fprintf(FilePtr,"\n");

    // excess adsorption
    sum=0.0;
    sum_squared=0.0;
    nr=NumberOfUnitCells[CurrentSystem].x*NumberOfUnitCells[CurrentSystem].y*NumberOfUnitCells[CurrentSystem].z;
    for(int i=0;i<NR_BLOCKS;i++)
    {
      if(BlockWeightedCount[CurrentSystem][i]>0.0)
      {
        REAL tmp=(NumberOfExcessMoleculesPerComponentAccumulated[CurrentSystem][j][i]/(nr * BlockWeightedCount[CurrentSystem][i]));
        sum+=tmp;
        sum_squared+=SQR(tmp);

        fprintf(FilePtr,"\tBlock[%2d] %-18.5lf [-]\n",i, tmp * nr);
      }
      else
        fprintf(FilePtr,"\tBlock[%2d] %-18.5lf [-]\n",i,(double)0.0);
    }
    avg=AVERAGE(sum);
    error=ERROR_CONFIDENCE_INTERVAL_95(sum,sum_squared);
    fprintf(FilePtr,"\t------------------------------------------------------------------------------\n");
    fprintf(FilePtr,"\tAverage loading excess                        %18.10lf +/- %18.10lf [-]\n", avg * nr, error * nr);
    fprintf(FilePtr,"\tAverage loading excess [molecules/unit cell]  %18.10lf +/- %18.10lf [-]\n", avg, error);

    fprintf(FilePtr,"\tAverage loading excess [mol/kg framework]              %18.10lf +/- %18.10lf [-]\n", conv_to_mol_per_kg * avg, conv_to_mol_per_kg * error);
    fprintf(FilePtr,"\tAverage loading excess [milligram/gram framework]      %18.10lf +/- %18.10lf [-]\n", conv_to_mg_per_g * avg, conv_to_mg_per_g * error);
    fprintf(FilePtr,"\tAverage loading excess [cm^3 (STP)/gr framework]       %18.10lf +/- %18.10lf [-]\n", conv_to_cm3_per_g * avg, conv_to_cm3_per_g * error);
    fprintf(FilePtr,"\tAverage loading excess [cm^3 (STP)/cm^3 framework]     %18.10lf +/- %18.10lf [-]\n", conv_to_cm3_per_cm3 * avg, conv_to_cm3_per_cm3 * error);
    fprintf(FilePtr,"\n");
  }

  // Average fixed volume Widom Rosenbluth factor
  fprintf(FilePtr,"\n");
  fprintf(FilePtr,"Average Widom Rosenbluth factor:\n");
  fprintf(FilePtr,"================================\n");
  for(int j=0;j<NumberOfComponents;j++)
  {
    sum=sum_squared=0.0;
    for(int i=0;i<NR_BLOCKS;i++)
    {
      if(WidomRosenbluthFactorCount[CurrentSystem][j][i]>0.0)
      {
        REAL tmp=WidomRosenbluthFactorAccumulated[CurrentSystem][j][i]/
                 WidomRosenbluthFactorCount[CurrentSystem][j][i];
        sum+=tmp;
        sum_squared+=SQR(tmp);
        fprintf(FilePtr,"\tBlock[%2d] %-lg [-]\n",i,(double)tmp);
      }
      else
        fprintf(FilePtr,"\tBlock[%2d] %-lg [-]\n",i,(double)0.0);
    }
    avg=AVERAGE(sum);
    error=ERROR_CONFIDENCE_INTERVAL_95(sum,sum_squared);
    fprintf(FilePtr,"\t------------------------------------------------------------------------------\n");
    fprintf(FilePtr,"\t[%s] Average Widom Rosenbluth-weight:   %lg +/- %lf [-]\n",
      Components[j].Name, avg, error);
  }

  // Average fixed volume chemical potential
  fprintf(FilePtr,"\n");
  fprintf(FilePtr,"Average Widom chemical potential:\n");
  fprintf(FilePtr,"=================================\n");
  for(int j=0;j<NumberOfComponents;j++)
  {
    sum=sum_squared=0.0;
    for(int i=0;i<NR_BLOCKS;i++)
    {
      if(WidomRosenbluthFactorCount[CurrentSystem][j][i]>0.0)
      {
        REAL tmp=(-log(WidomRosenbluthFactorAccumulated[CurrentSystem][j][i]/
                      WidomRosenbluthFactorCount[CurrentSystem][j][i])/Beta[CurrentSystem])*ENERGY_TO_KELVIN+
                 (-log(WidomIdealGasAccumulated[CurrentSystem][j][i]/
                       WidomRosenbluthFactorCount[CurrentSystem][j][i])/Beta[CurrentSystem])*ENERGY_TO_KELVIN;
        sum+=tmp;
        sum_squared+=SQR(tmp);
        fprintf(FilePtr,"\tBlock[%2d] %-lg [-]\n",i,(double)tmp);
      }
      else
        fprintf(FilePtr,"\tBlock[%2d] %-lg [-]\n",i,(double)0.0);
    }
    avg=AVERAGE(sum);
    error=ERROR_CONFIDENCE_INTERVAL_95(sum,sum_squared);
    fprintf(FilePtr,"\t------------------------------------------------------------------------------\n");
    fprintf(FilePtr,"\t[%s] Average chemical potential:   %lg +/- %lf [K]\n",
      Components[j].Name, avg, error);
  }

  // Average fixed volume Widom ideal-gas contribution
  fprintf(FilePtr,"\n");
  fprintf(FilePtr,"Average Widom Ideal-gas contribution:\n");
  fprintf(FilePtr,"=====================================\n");
  for(int j=0;j<NumberOfComponents;j++)
  {
    sum=sum_squared=0.0;
    for(int i=0;i<NR_BLOCKS;i++)
    {
      if(WidomRosenbluthFactorCount[CurrentSystem][j][i]>0.0)
      {
        REAL tmp=(-log(WidomIdealGasAccumulated[CurrentSystem][j][i]/
                     WidomRosenbluthFactorCount[CurrentSystem][j][i])/Beta[CurrentSystem])*ENERGY_TO_KELVIN;
        sum+=tmp;
        sum_squared+=SQR(tmp);
        fprintf(FilePtr,"\tBlock[%2d] %-lg [-]\n",i,(double)tmp);
      }
      else
        fprintf(FilePtr,"\tBlock[%2d] %-lg [-]\n",i,(double)0.0);
    }
    avg=AVERAGE(sum);
    error=ERROR_CONFIDENCE_INTERVAL_95(sum,sum_squared);
    fprintf(FilePtr,"\t------------------------------------------------------------------------------\n");
    fprintf(FilePtr,"\t[%s] Average Widom Ideal-gas chemical potential:   %lg +/- %lf [-]\n",
      Components[j].Name, avg, error);
  }

  // Average fixed volume Widom excess contribution
  fprintf(FilePtr,"\n");
  fprintf(FilePtr,"Average Widom excess contribution:\n");
  fprintf(FilePtr,"==================================\n");
  for(int j=0;j<NumberOfComponents;j++)
  {
    sum=sum_squared=0.0;
    for(int i=0;i<NR_BLOCKS;i++)
    {
      if(WidomRosenbluthFactorCount[CurrentSystem][j][i]>0.0)
      {
        REAL tmp=(-log(WidomRosenbluthFactorAccumulated[CurrentSystem][j][i]/
                  WidomRosenbluthFactorCount[CurrentSystem][j][i])/Beta[CurrentSystem])*ENERGY_TO_KELVIN;
        sum+=tmp;
        sum_squared+=SQR(tmp);
        fprintf(FilePtr,"\tBlock[%2d] %-lg [-]\n",i,(double)tmp);
      }
      else
        fprintf(FilePtr,"\tBlock[%2d] %-lg [-]\n",i,(double)0.0);
    }
    avg=AVERAGE(sum);
    error=ERROR_CONFIDENCE_INTERVAL_95(sum,sum_squared);
    fprintf(FilePtr,"\t------------------------------------------------------------------------------\n");
    fprintf(FilePtr,"\t[%s] Average Widom excess chemical potential:   %lg +/- %lf [-]\n", Components[j].Name, avg, error);
  }



  // Average Gibbs Widom Rosenbluth factor
  fprintf(FilePtr,"\n");
  fprintf(FilePtr,"Average Gibbs Widom Rosenbluth factor:\n");
  fprintf(FilePtr,"======================================\n");
  for(int j=0;j<NumberOfComponents;j++)
  {
    sum=sum_squared=0.0;
    for(int i=0;i<NR_BLOCKS;i++)
    {
      if(GibbsWidomRosenbluthFactorCount[CurrentSystem][j][i]>0.0)
      {
        REAL tmp=GibbsWidomRosenbluthFactorAccumulated[CurrentSystem][j][i]/
                 GibbsWidomRosenbluthFactorCount[CurrentSystem][j][i];
        sum+=tmp;
        sum_squared+=SQR(tmp);
        fprintf(FilePtr,"\tBlock[%2d] %-lg [-]\n",i,(double)tmp);
      }
      else
        fprintf(FilePtr,"\tBlock[%2d] %-lg [-]\n",i,(double)0.0);
    }
    avg=AVERAGE(sum);
    error=ERROR_CONFIDENCE_INTERVAL_95(sum,sum_squared);
    fprintf(FilePtr,"\t------------------------------------------------------------------------------\n");
    fprintf(FilePtr,"\t[%s] Average Gibbs Widom Rosenbluth-weight:   %lg +/- %lf [-]\n",
      Components[j].Name, avg, error);
  }


  // Average Gibbs chemical potential
  fprintf(FilePtr,"\n");
  fprintf(FilePtr,"Average Gibbs Widom chemical potential:\n");
  fprintf(FilePtr,"=======================================\n");
  for(int j=0;j<NumberOfComponents;j++)
  {
    sum=sum_squared=0.0;
    for(int i=0;i<NR_BLOCKS;i++)
    {
      if(GibbsWidomRosenbluthFactorCount[CurrentSystem][j][i]>0.0)
      {
        REAL tmp=(-log(GibbsWidomRosenbluthFactorAccumulated[CurrentSystem][j][i]/
                  GibbsWidomRosenbluthFactorCount[CurrentSystem][j][i])/Beta[CurrentSystem])*ENERGY_TO_KELVIN;
        sum+=tmp;
        sum_squared+=SQR(tmp);
        fprintf(FilePtr,"\tBlock[%2d] %-lg [-]\n",i,(double)tmp);
      }
      else
        fprintf(FilePtr,"\tBlock[%2d] %-lg [-]\n",i,(double)0.0);
    }
    avg=AVERAGE(sum);
    error=ERROR_CONFIDENCE_INTERVAL_95(sum,sum_squared);
    fprintf(FilePtr,"\t------------------------------------------------------------------------------\n");
    fprintf(FilePtr,"\t[%s] Average Gibbs chemical potential:   %lg +/- %lf [K]\n",
      Components[j].Name, avg, error);
  }
  
  
  // Added by Ambroise de Izarra
  //-------------------------------------------------------------------  
  fprintf(FilePtr,"\n");
  fprintf(FilePtr,"Average osmostat calibration for alchemical transformation:\n");
  fprintf(FilePtr,"=======================================\n"); 
   
  if(WidomOsmostat[CurrentSystem][0]>0.0)
  {
	    REAL RatioSpecies=1.00;
  
		for(int i=0;i<=(MultiplicitySalt[CurrentAlchemicalReaction][0]+MultiplicitySalt[CurrentAlchemicalReaction][1]-1);i++)
			RatioSpecies *= ((Components[SolventIndex].NumberOfMolecules[CurrentSystem]-i));

		for(int i=1;i<=MultiplicitySalt[CurrentAlchemicalReaction] [0];i++)	  
			RatioSpecies *= 1.0/(Components[SaltIndex[CurrentAlchemicalReaction][0]].NumberOfMolecules[CurrentSystem]+i);
	
		for(int i=1;i<=MultiplicitySalt[CurrentAlchemicalReaction][1];i++)
			RatioSpecies *= 1.0/(Components[SaltIndex[CurrentAlchemicalReaction][1]].NumberOfMolecules[CurrentSystem]+i);
		
		REAL increase;
		sum=sum_squared=0.0;
		
		int blocksize=floor(SizeAlchemicalWorkStore/NR_BLOCKS);

		for(int i=0;i<NR_BLOCKS;i++)
		{
			increase = 0.0;
			for(int j=0;j<blocksize;j++)
			{	
		        increase+= exp(-AlchemicalWorkStore[i*blocksize+j]*Beta[CurrentSystem]);
			}

			REAL tmp=-(log((RatioSpecies)*(increase/blocksize))/Beta[CurrentSystem])*ENERGY_TO_KELVIN;
			sum+=tmp;
			sum_squared+=SQR(tmp);
			fprintf(FilePtr,"\tBlock[%2d] %-lg [-]\n",i,(double)tmp);
		}
   }
   else
	{
		for(int i=0;i<NR_BLOCKS;i++)
		{
			fprintf(FilePtr,"\tBlock[%2d] %-lg [-]\n",i,(double)0.0);
		}
    } 
      
    avg=AVERAGE(sum);
    error=ERROR_CONFIDENCE_INTERVAL_95(sum,sum_squared);
    fprintf(FilePtr,"\t------------------------------------------------------------------------------\n");
    fprintf(FilePtr,"\tAverage osmostat calibration for alchemical transformation:   %lg +/- %lf [K]\n",
      avg, error);
    fprintf(FilePtr,"\tAverage osmostat calibration for alchemical transformation:   %lg +/- %lf [kJ/mol]\n",
      avg*KELVIN_TO_KJ_PER_MOL, error*KELVIN_TO_KJ_PER_MOL);
  //-------------------------------------------------------------------
  
  
  
  // Average Gibbs Widom ideal-gas contribution
  fprintf(FilePtr,"\n");
  fprintf(FilePtr,"Average Gibbs Widom Ideal-gas contribution:\n");
  fprintf(FilePtr,"===========================================\n");
  for(int j=0;j<NumberOfComponents;j++)
  {
    sum=sum_squared=0.0;
    for(int i=0;i<NR_BLOCKS;i++)
    {
      if(GibbsWidomRosenbluthFactorCount[CurrentSystem][j][i]>0.0)
      {
        REAL tmp=(-log(GibbsWidomIdealGasAccumulated[CurrentSystem][j][i]/
                  GibbsWidomRosenbluthFactorCount[CurrentSystem][j][i])/Beta[CurrentSystem])*ENERGY_TO_KELVIN;
        sum+=tmp;
        sum_squared+=SQR(tmp);
        fprintf(FilePtr,"\tBlock[%2d] %-lg [-]\n",i,(double)tmp);
      }
      else
        fprintf(FilePtr,"\tBlock[%2d] %-lg [-]\n",i,(double)0.0);
    }
    avg=AVERAGE(sum);
    error=ERROR_CONFIDENCE_INTERVAL_95(sum,sum_squared);
    fprintf(FilePtr,"\t------------------------------------------------------------------------------\n");
    fprintf(FilePtr,"\t[%s] Average Gibbs Ideal-gas chemical potential:   %lg +/- %lf [-]\n",
      Components[j].Name, avg, error);
  }

  // Average Gibbs Widom excess contribution
  fprintf(FilePtr,"\n");
  fprintf(FilePtr,"Average Gibbs Widom excess contribution:\n");
  fprintf(FilePtr,"===========================================\n");
  for(int j=0;j<NumberOfComponents;j++)
  {
    sum=sum_squared=0.0;
    for(int i=0;i<NR_BLOCKS;i++)
    {
      if(GibbsWidomRosenbluthFactorCount[CurrentSystem][j][i]>0.0)
      {
        REAL tmp=(-log(GibbsWidomRosenbluthFactorAccumulated[CurrentSystem][j][i]/
                      GibbsWidomRosenbluthFactorCount[CurrentSystem][j][i])/Beta[CurrentSystem])*ENERGY_TO_KELVIN
                 +(log(GibbsWidomIdealGasAccumulated[CurrentSystem][j][i]/
                      GibbsWidomRosenbluthFactorCount[CurrentSystem][j][i])/Beta[CurrentSystem])*ENERGY_TO_KELVIN;
        sum+=tmp;
        sum_squared+=SQR(tmp);
        fprintf(FilePtr,"\tBlock[%2d] %-lg [-]\n",i,(double)tmp);
      }
      else
        fprintf(FilePtr,"\tBlock[%2d] %-lg [-]\n",i,(double)0.0);
    }
    avg=AVERAGE(sum);
    error=ERROR_CONFIDENCE_INTERVAL_95(sum,sum_squared);
    fprintf(FilePtr,"\t------------------------------------------------------------------------------\n");
    fprintf(FilePtr,"\t[%s] Average Gibbs excess chemical potential:   %lg +/- %lf [-]\n",
      Components[j].Name, avg, error);
  }


  // Average Henry Coefficient
  fprintf(FilePtr,"\n");
  fprintf(FilePtr,"Average Henry coefficient:\n");
  fprintf(FilePtr,"==========================\n");
  REAL FrameworkDensity=1e-3*Framework[CurrentSystem].FrameworkMass/(Volume[CurrentSystem]*CUBE(ANGSTROM)*AVOGADRO_CONSTANT);
  REAL Temperature=therm_baro_stats.ExternalTemperature[CurrentSystem];
  for(int j=0;j<NumberOfComponents;j++)
  {
    sum=sum_squared=0.0;
    for(int i=0;i<NR_BLOCKS;i++)
    {
      if(WidomRosenbluthFactorCount[CurrentSystem][j][i]>0.0)
      {
        REAL tmp=(1.0/(MOLAR_GAS_CONSTANT*Temperature*FrameworkDensity))*
                 WidomRosenbluthFactorAccumulated[CurrentSystem][j][i]/
                 (Components[j].IdealGasRosenbluthWeight[CurrentSystem]*WidomRosenbluthFactorCount[CurrentSystem][j][i]);
        sum+=tmp;
        sum_squared+=SQR(tmp);
        fprintf(FilePtr,"\tBlock[%2d] %-lg [mol/kg/Pa]\n",i,(double)tmp);
      }
      else
        fprintf(FilePtr,"\tBlock[%2d] %-lg [mol/kg/Pa]\n",i,(double)0.0);
    }
    avg=AVERAGE(sum);
    error=ERROR_CONFIDENCE_INTERVAL_95(sum,sum_squared);
    fprintf(FilePtr,"\t------------------------------------------------------------------------------\n");
    fprintf(FilePtr,"\t[%s] Average Henry coefficient:  %lg +/- %lg [mol/kg/Pa]\n",
      Components[j].Name, avg, error);
  }

  // Average adsorption energy
  fprintf(FilePtr,"\n");
  fprintf(FilePtr,"Average adsorption energy <U_gh>_1-<U_h>_0 obtained from Widom-insertion:\n");
  fprintf(FilePtr,"(Note: the total heat of adsorption is dH=<U_gh>_1-<U_h>_0 - <U_g> - RT)\n");
  fprintf(FilePtr,"=========================================================================\n");
  Temperature=therm_baro_stats.ExternalTemperature[CurrentSystem];
  for(int j=0;j<NumberOfComponents;j++)
  {
    if(Components[j].Widom)
    {
      sum=sum_squared=0.0;
      for(int i=0;i<NR_BLOCKS;i++)
      {
        if(WidomRosenbluthFactorCount[CurrentSystem][j][i]>0.0)
        {
          REAL tmp=(WidomEnergyDifferenceAccumulated[CurrentSystem][j][i]/WidomRosenbluthFactorAccumulated[CurrentSystem][j][i]-
                    WidomEnergyFrameworkAccumulated[CurrentSystem][j][i]/WidomEnergyFrameworkCount[CurrentSystem][j][i])*ENERGY_TO_KELVIN;
          sum+=tmp;
          sum_squared+=SQR(tmp);
          fprintf(FilePtr,"\tBlock[%2d] %-18.10lf [K]\n",i,(double)tmp);
        }
        else
          fprintf(FilePtr,"\tBlock[%2d] %-18.10lf [K]\n",i,(double)0.0);
      }
      avg=AVERAGE(sum);
      error=ERROR_CONFIDENCE_INTERVAL_95(sum,sum_squared);
      fprintf(FilePtr,"\t------------------------------------------------------------------------------\n");
      fprintf(FilePtr,"\t[%s] Average  <U_gh>_1-<U_h>_0:  %18.10lf +/- %18.10lf [K]       (%18.10lf +/- %18.10lf kJ/mol)\n",
        Components[j].Name, avg, error, avg * KELVIN_TO_KJ_PER_MOL, error * KELVIN_TO_KJ_PER_MOL);
    }
  }
}

static int versionNumber=1;

void WriteRestartStatistics(FILE *FilePtr)
{
  int i,j,k;
  int NumberOfBlocks=NR_BLOCKS;
  REAL Check;

  fwrite(&versionNumber,sizeof(int),1,FilePtr);

  NumberOfBlocks=NR_BLOCKS;

  fwrite(&NumberOfBlocks,sizeof(int),1,FilePtr);
  fwrite(&Block,sizeof(int),1,FilePtr);
  fwrite(BlockCycle,sizeof(long long),NumberOfBlocks,FilePtr);

  for(i=0;i<NumberOfSystems;i++)
  {
    fwrite(BlockWeightedCount[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(BlockCount[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    fwrite(UHostHostAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UAdsorbateAdsorbateAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UCationCationAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UHostAdsorbateAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UHostCationAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UAdsorbateCationAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    fwrite(UHostHostVDWAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UAdsorbateAdsorbateVDWAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UCationCationVDWAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UHostAdsorbateVDWAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UHostCationVDWAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UAdsorbateCationVDWAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    fwrite(UHostHostCoulombAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UAdsorbateAdsorbateCoulombAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UCationCationCoulombAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UHostAdsorbateCoulombAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UHostCationCoulombAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UAdsorbateCationCoulombAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    fwrite(UHostBondAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UHostUreyBradleyAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UHostBendAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UHostInversionBendAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UHostTorsionAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UHostImproperTorsionAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UHostBondBondAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UHostBendBendAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UHostBondBendAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UHostBondTorsionAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UHostBendTorsionAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    fwrite(UAdsorbateBondAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UAdsorbateUreyBradleyAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UAdsorbateBondBondAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UAdsorbateBendAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UAdsorbateInversionBendAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UAdsorbateTorsionAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UAdsorbateImproperTorsionAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UAdsorbateBondBondAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UAdsorbateBendBendAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UAdsorbateBondBendAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UAdsorbateBondTorsionAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UAdsorbateBendTorsionAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UAdsorbateIntraVDWAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UAdsorbateIntraChargeChargeAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UAdsorbateIntraChargeBondDipoleAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UAdsorbateIntraBondDipoleBondDipoleAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    fwrite(UCationBondAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UCationUreyBradleyAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UCationBendAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UCationInversionBendAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UCationTorsionAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UCationImproperTorsionAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UCationBondBondAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UCationBendBendAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UCationBondBendAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UCationBondTorsionAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UCationBendTorsionAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UCationIntraVDWAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UCationIntraChargeChargeAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UCationIntraChargeBondDipoleAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UCationIntraBondDipoleBondDipoleAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    fwrite(UHostPolarizationAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UAdsorbatePolarizationAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UCationPolarizationAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    fwrite(UHostBackPolarizationAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UAdsorbateBackPolarizationAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UCationBackPolarizationAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    fwrite(UTailCorrectionAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UDistanceConstraintsAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UAngleConstraintsAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UDihedralConstraintsAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UInversionBendConstraintsAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UOutOfPlaneDistanceConstraintsAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UExclusionConstraintsAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UTotalAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    fwrite(TotalSystemDipoleAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(TotalSystemDipoleSquaredAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(TotalSystemNormDipoleAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(TotalSystemNormDipoleSquaredAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    fwrite(NumberOfIntegerMoleculesAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(NumberOfMoleculesSquaredAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(DensityAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(InverseDensityAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(GibbsInverseDensityAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    fwrite(BoxAccumulated[i],sizeof(VECTOR),NumberOfBlocks,FilePtr);
    fwrite(BoxAXAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(BoxAYAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(BoxAZAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(BoxBXAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(BoxBYAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(BoxBZAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(BoxCXAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(BoxCYAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(BoxCZAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(BoxLengthAccumulated[i],sizeof(VECTOR),NumberOfBlocks,FilePtr);
    fwrite(AlphaAngleAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(BetaAngleAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(GammaAngleAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(VolumeAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(VolumeSquaredAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    fwrite(TotalEnergyAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(TotalEnergySquaredAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(EnthalpyAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(EnthalpySquaredAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(EnthalpyTimesVolumeAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(EnthalpyTimesEnergyAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    fwrite(TemperatureAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(TemperatureCellAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(TemperatureTranslationAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(TemperatureRotationAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(TemperatureRotationAdsorbateAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(TemperatureTranslationAdsorbateAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    fwrite(TemperatureAdsorbatesAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(TemperatureCationsAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(TemperatureFrameworkAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    fwrite(MolecularPressureAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(MolecularStressTensorAccumulated[i],sizeof(REAL_MATRIX3x3),NumberOfBlocks,FilePtr);
    fwrite(PressureIdealGasPartAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(PressureExcessPartAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(PressureTailCorrectionAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(PressureAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UNoseHooverAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(HeatOfVaporization[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(EnergyPerMolecule[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(VolumePerMolecule[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(CompressibilityAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    fwrite(BornTermAccumulated[i],sizeof(REAL_MATRIX9x9),NumberOfBlocks,FilePtr);
    fwrite(ConfigurationalStressFluctuationTermAccumulated[i],sizeof(REAL_MATRIX9x9),NumberOfBlocks,FilePtr);
    fwrite(ConfigurationalStressTensorAccumulated[i],sizeof(REAL_MATRIX3x3),NumberOfBlocks,FilePtr);
    fwrite(StressTensorAccumulated[i],sizeof(REAL_MATRIX3x3),NumberOfBlocks,FilePtr);

    fwrite(SurfaceAreaFrameworkAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    for(j=0;j<Framework[i].NumberOfFrameworks;j++)
      fwrite(SurfaceAreaFrameworksAccumulated[i][j],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(SurfaceAreaCationsAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(SurfaceAreaCount[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    fwrite(TotalEnergyTimesNumberOfMoleculesAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    for(j=0;j<NumberOfComponents;j++)
    {

      fwrite(TotalEnergyTimesNumberOfMoleculesPerComponentAccumulated[i][j],sizeof(REAL),NumberOfBlocks,FilePtr);
      fwrite(HostAdsorbateEnergyTimesNumberOfMoleculesAccumulated[i][j],sizeof(REAL),NumberOfBlocks,FilePtr);
      fwrite(AdsorbateAdsorbateEnergyTimesNumberOfMoleculesAccumulated[i][j],sizeof(REAL),NumberOfBlocks,FilePtr);

      fwrite(NumberOfIntegerMoleculesPerComponentAccumulated[i][j],sizeof(REAL),NumberOfBlocks,FilePtr);
      fwrite(NumberOfFractionalMoleculesPerComponentAccumulated[i][j],sizeof(REAL),NumberOfBlocks,FilePtr);
      fwrite(NumberOfExcessMoleculesPerComponentAccumulated[i][j],sizeof(REAL),NumberOfBlocks,FilePtr);
      fwrite(DensityPerComponentAccumulated[i][j],sizeof(REAL),NumberOfBlocks,FilePtr);

      fwrite(WidomRosenbluthFactorAccumulated[i][j],sizeof(REAL),NumberOfBlocks,FilePtr);
      fwrite(WidomIdealGasAccumulated[i][j],sizeof(REAL),NumberOfBlocks,FilePtr);
      fwrite(WidomRosenbluthFactorCount[i][j],sizeof(REAL),NumberOfBlocks,FilePtr);
      fwrite(WidomEnergyDifferenceAccumulated[i][j],sizeof(REAL),NumberOfBlocks,FilePtr);
      fwrite(WidomEnergyFrameworkAccumulated[i][j],sizeof(REAL),NumberOfBlocks,FilePtr);
      fwrite(WidomEnergyFrameworkCount[i][j],sizeof(REAL),NumberOfBlocks,FilePtr);

      fwrite(GibbsWidomRosenbluthFactorAccumulated[i][j],sizeof(REAL),NumberOfBlocks,FilePtr);
      fwrite(GibbsWidomIdealGasAccumulated[i][j],sizeof(REAL),NumberOfBlocks,FilePtr);
      fwrite(GibbsWidomRosenbluthFactorCount[i][j],sizeof(REAL),NumberOfBlocks,FilePtr);
      fwrite(GibbsWidomEnergyDifferenceAccumulated[i][j],sizeof(REAL),NumberOfBlocks,FilePtr);
      fwrite(GibbsWidomEnergyFrameworkAccumulated[i][j],sizeof(REAL),NumberOfBlocks,FilePtr);
      fwrite(GibbsWidomEnergyFrameworkCount[i][j],sizeof(REAL),NumberOfBlocks,FilePtr);

      fwrite(PrincipleMomentsOfInertiaAccumulated[i][j],sizeof(REAL),NumberOfBlocks,FilePtr);
      fwrite(PrincipleMomentsOfInertiaCount[i][j],sizeof(REAL),NumberOfBlocks,FilePtr);

      for(k=0;k<NumberOfComponents;k++)
        fwrite(NumberOfMoleculesPerComponentSquaredAccumulated[i][j][k],sizeof(REAL),NumberOfBlocks,FilePtr);
    }
  }

  Check=123456789.0;
  fwrite(&Check,1,sizeof(REAL),FilePtr);
}

void AllocateStatisticsMemory(void)
{
  int i,j,k;
  int NumberOfBlocks;

  NumberOfBlocks=NR_BLOCKS;

  BlockCycle=(long long*)calloc(NumberOfBlocks,sizeof(long long));
  BlockWeightedCount=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  BlockCount=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));

  UHostHostAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UAdsorbateAdsorbateAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UCationCationAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UHostAdsorbateAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UHostCationAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UAdsorbateCationAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));

  UHostHostVDWAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UAdsorbateAdsorbateVDWAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UCationCationVDWAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UHostAdsorbateVDWAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UHostCationVDWAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UAdsorbateCationVDWAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));

  UHostHostCoulombAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UHostAdsorbateCoulombAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UHostCationCoulombAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UAdsorbateAdsorbateCoulombAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UCationCationCoulombAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UAdsorbateCationCoulombAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));

  UHostBondAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UHostUreyBradleyAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UHostBendAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UHostInversionBendAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UHostTorsionAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UHostImproperTorsionAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UHostBondBondAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UHostBendBendAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UHostBondBendAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UHostBondTorsionAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UHostBendTorsionAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));

  UAdsorbateBondAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UAdsorbateUreyBradleyAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UAdsorbateBondBondAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UAdsorbateBendAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UAdsorbateInversionBendAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UAdsorbateTorsionAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UAdsorbateImproperTorsionAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UAdsorbateBondBondAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UAdsorbateBendBendAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UAdsorbateBondBendAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UAdsorbateBondTorsionAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UAdsorbateBendTorsionAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UAdsorbateIntraVDWAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UAdsorbateIntraChargeChargeAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UAdsorbateIntraChargeBondDipoleAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UAdsorbateIntraBondDipoleBondDipoleAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));

  UCationBondAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UCationUreyBradleyAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UCationBendAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UCationInversionBendAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UCationTorsionAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UCationImproperTorsionAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UCationBondBondAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UCationBendBendAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UCationBondBendAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UCationBondTorsionAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UCationBendTorsionAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UCationIntraVDWAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UCationIntraChargeChargeAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UCationIntraChargeBondDipoleAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UCationIntraBondDipoleBondDipoleAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));

  UHostPolarizationAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UAdsorbatePolarizationAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UCationPolarizationAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UHostBackPolarizationAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UAdsorbateBackPolarizationAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UCationBackPolarizationAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));

  UTailCorrectionAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UDistanceConstraintsAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UAngleConstraintsAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UDihedralConstraintsAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UInversionBendConstraintsAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UOutOfPlaneDistanceConstraintsAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UExclusionConstraintsAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UTotalAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));

  TotalSystemDipoleAccumulated=(VECTOR**)calloc(NumberOfSystems,sizeof(VECTOR*));
  TotalSystemDipoleSquaredAccumulated=(VECTOR**)calloc(NumberOfSystems,sizeof(VECTOR*));
  TotalSystemNormDipoleAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  TotalSystemNormDipoleSquaredAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));

  NumberOfIntegerMoleculesAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  NumberOfMoleculesSquaredAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  DensityAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  InverseDensityAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  GibbsInverseDensityAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));

  BoxAccumulated=(VECTOR**)calloc(NumberOfSystems,sizeof(VECTOR*));
  BoxAXAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  BoxAYAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  BoxAZAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  BoxBXAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  BoxBYAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  BoxBZAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  BoxCXAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  BoxCYAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  BoxCZAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  BoxLengthAccumulated=(VECTOR**)calloc(NumberOfSystems,sizeof(VECTOR*));
  AlphaAngleAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  BetaAngleAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  GammaAngleAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  VolumeAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  VolumeSquaredAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));

  TotalEnergyAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  TotalEnergySquaredAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  EnthalpyAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  EnthalpySquaredAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  EnthalpyTimesVolumeAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  EnthalpyTimesEnergyAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));

  TemperatureAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  TemperatureCellAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  TemperatureTranslationAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  TemperatureRotationAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  TemperatureRotationAdsorbateAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  TemperatureTranslationAdsorbateAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));

  TemperatureAdsorbatesAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  TemperatureCationsAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  TemperatureFrameworkAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));

  MolecularPressureAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  MolecularStressTensorAccumulated=(REAL_MATRIX3x3**)calloc(NumberOfSystems,sizeof(REAL_MATRIX3x3*));
  PressureIdealGasPartAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  PressureExcessPartAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  PressureTailCorrectionAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  PressureAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UNoseHooverAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  HeatOfVaporization=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  EnergyPerMolecule=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  VolumePerMolecule=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  CompressibilityAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));

  BornTermAccumulated=(REAL_MATRIX9x9**)calloc(NumberOfSystems,sizeof(REAL_MATRIX9x9*));
  ConfigurationalStressFluctuationTermAccumulated=(REAL_MATRIX9x9**)calloc(NumberOfSystems,sizeof(REAL_MATRIX9x9*));
  ConfigurationalStressTensorAccumulated=(REAL_MATRIX3x3**)calloc(NumberOfSystems,sizeof(REAL_MATRIX3x3*));
  StressTensorAccumulated=(REAL_MATRIX3x3**)calloc(NumberOfSystems,sizeof(REAL_MATRIX3x3*));

  SurfaceAreaFrameworkAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  SurfaceAreaFrameworksAccumulated=(REAL***)calloc(NumberOfSystems,sizeof(REAL*));
  SurfaceAreaCationsAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  SurfaceAreaCount=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));

  TotalEnergyTimesNumberOfMoleculesAccumulated=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  TotalEnergyTimesNumberOfMoleculesPerComponentAccumulated=(REAL***)calloc(NumberOfSystems,sizeof(REAL**));
  HostAdsorbateEnergyTimesNumberOfMoleculesAccumulated=(REAL***)calloc(NumberOfSystems,sizeof(REAL*));
  AdsorbateAdsorbateEnergyTimesNumberOfMoleculesAccumulated=(REAL***)calloc(NumberOfSystems,sizeof(REAL*));

  NumberOfIntegerMoleculesPerComponentAccumulated=(REAL***)calloc(NumberOfSystems,sizeof(REAL**));
  NumberOfFractionalMoleculesPerComponentAccumulated=(REAL***)calloc(NumberOfSystems,sizeof(REAL**));
  NumberOfExcessMoleculesPerComponentAccumulated=(REAL***)calloc(NumberOfSystems,sizeof(REAL**));
  DensityPerComponentAccumulated=(REAL***)calloc(NumberOfSystems,sizeof(REAL**));

  WidomRosenbluthFactorAccumulated=(REAL***)calloc(NumberOfSystems,sizeof(REAL**));
  WidomIdealGasAccumulated=(REAL***)calloc(NumberOfSystems,sizeof(REAL**));
  WidomRosenbluthFactorCount=(REAL***)calloc(NumberOfSystems,sizeof(REAL**));
  WidomEnergyDifferenceAccumulated=(REAL***)calloc(NumberOfSystems,sizeof(REAL**));
  WidomEnergyFrameworkAccumulated=(REAL***)calloc(NumberOfSystems,sizeof(REAL**));
  WidomEnergyFrameworkCount=(REAL***)calloc(NumberOfSystems,sizeof(REAL**));

  GibbsWidomRosenbluthFactorAccumulated=(REAL***)calloc(NumberOfSystems,sizeof(REAL**));
  GibbsWidomIdealGasAccumulated=(REAL***)calloc(NumberOfSystems,sizeof(REAL**));
  GibbsWidomRosenbluthFactorCount=(REAL***)calloc(NumberOfSystems,sizeof(REAL**));
  GibbsWidomEnergyDifferenceAccumulated=(REAL***)calloc(NumberOfSystems,sizeof(REAL**));
  GibbsWidomEnergyFrameworkAccumulated=(REAL***)calloc(NumberOfSystems,sizeof(REAL**));
  GibbsWidomEnergyFrameworkCount=(REAL***)calloc(NumberOfSystems,sizeof(REAL**));
  GibbsInverseDensityPerComponentAccumulated=(REAL***)calloc(NumberOfSystems,sizeof(REAL**));

  PrincipleMomentsOfInertiaAccumulated=(VECTOR***)calloc(NumberOfSystems,sizeof(VECTOR**));
  PrincipleMomentsOfInertiaCount=(REAL***)calloc(NumberOfSystems,sizeof(REAL**));

  NumberOfMoleculesPerComponentSquaredAccumulated=(REAL****)calloc(NumberOfSystems,sizeof(REAL***));

  for(i=0;i<NumberOfSystems;i++)
  {
    BlockWeightedCount[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    BlockCount[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));

    UHostHostAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UAdsorbateAdsorbateAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UCationCationAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UHostAdsorbateAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UHostCationAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UAdsorbateCationAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));

    UHostHostVDWAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UAdsorbateAdsorbateVDWAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UCationCationVDWAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UHostAdsorbateVDWAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UHostCationVDWAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UAdsorbateCationVDWAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));

    UHostHostCoulombAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UAdsorbateAdsorbateCoulombAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UCationCationCoulombAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UHostAdsorbateCoulombAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UHostCationCoulombAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UAdsorbateCationCoulombAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));

    UHostBondAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UHostUreyBradleyAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UHostBendAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UHostInversionBendAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UHostTorsionAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UHostImproperTorsionAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UHostBondBondAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UHostBendBendAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UHostBondBendAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UHostBondTorsionAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UHostBendTorsionAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));

    UAdsorbateBondAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UAdsorbateUreyBradleyAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UAdsorbateBondBondAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UAdsorbateBendAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UAdsorbateInversionBendAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UAdsorbateTorsionAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UAdsorbateImproperTorsionAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UAdsorbateBondBondAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UAdsorbateBendBendAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UAdsorbateBondBendAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UAdsorbateBondTorsionAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UAdsorbateBendTorsionAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UAdsorbateIntraVDWAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UAdsorbateIntraChargeChargeAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UAdsorbateIntraChargeBondDipoleAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UAdsorbateIntraBondDipoleBondDipoleAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));

    UCationBondAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UCationUreyBradleyAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UCationBendAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UCationInversionBendAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UCationTorsionAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UCationImproperTorsionAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UCationBondBondAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UCationBendBendAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UCationBondBendAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UCationBondTorsionAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UCationBendTorsionAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UCationIntraVDWAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UCationIntraChargeChargeAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UCationIntraChargeBondDipoleAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UCationIntraBondDipoleBondDipoleAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));

    UHostPolarizationAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UAdsorbatePolarizationAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UCationPolarizationAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UHostBackPolarizationAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UAdsorbateBackPolarizationAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UCationBackPolarizationAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));

    UTailCorrectionAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UDistanceConstraintsAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UAngleConstraintsAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UDihedralConstraintsAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UInversionBendConstraintsAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UOutOfPlaneDistanceConstraintsAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UExclusionConstraintsAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UTotalAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));

    TotalSystemDipoleAccumulated[i]=(VECTOR*)calloc(NumberOfBlocks,sizeof(VECTOR));
    TotalSystemDipoleSquaredAccumulated[i]=(VECTOR*)calloc(NumberOfBlocks,sizeof(VECTOR));
    TotalSystemNormDipoleAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    TotalSystemNormDipoleSquaredAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));

    NumberOfIntegerMoleculesAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    NumberOfMoleculesSquaredAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    DensityAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    InverseDensityAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    GibbsInverseDensityAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));

    BoxAccumulated[i]=(VECTOR*)calloc(NumberOfBlocks,sizeof(VECTOR));
    BoxAXAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    BoxAYAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    BoxAZAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    BoxBXAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    BoxBYAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    BoxBZAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    BoxCXAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    BoxCYAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    BoxCZAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    BoxLengthAccumulated[i]=(VECTOR*)calloc(NumberOfBlocks,sizeof(VECTOR));
    AlphaAngleAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    BetaAngleAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    GammaAngleAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    VolumeAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    VolumeSquaredAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));

    TotalEnergyAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    TotalEnergySquaredAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    EnthalpyAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    EnthalpySquaredAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    EnthalpyTimesVolumeAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    EnthalpyTimesEnergyAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));

    TemperatureAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    TemperatureCellAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    TemperatureTranslationAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    TemperatureRotationAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    TemperatureRotationAdsorbateAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    TemperatureTranslationAdsorbateAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));

    TemperatureAdsorbatesAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    TemperatureCationsAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    TemperatureFrameworkAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));

    MolecularPressureAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    MolecularStressTensorAccumulated[i]=(REAL_MATRIX3x3*)calloc(NumberOfBlocks,sizeof(REAL_MATRIX3x3));
    PressureIdealGasPartAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    PressureExcessPartAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    PressureTailCorrectionAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    PressureAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UNoseHooverAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    HeatOfVaporization[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    EnergyPerMolecule[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    VolumePerMolecule[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    CompressibilityAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));

    BornTermAccumulated[i]=(REAL_MATRIX9x9*)calloc(NumberOfBlocks,sizeof(REAL_MATRIX9x9));
    ConfigurationalStressFluctuationTermAccumulated[i]=(REAL_MATRIX9x9*)calloc(NumberOfBlocks,sizeof(REAL_MATRIX9x9));
    ConfigurationalStressTensorAccumulated[i]=(REAL_MATRIX3x3*)calloc(NumberOfBlocks,sizeof(REAL_MATRIX3x3));
    StressTensorAccumulated[i]=(REAL_MATRIX3x3*)calloc(NumberOfBlocks,sizeof(REAL_MATRIX3x3));

    SurfaceAreaFrameworkAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    SurfaceAreaFrameworksAccumulated[i]=(REAL**)calloc(Framework[i].NumberOfFrameworks,sizeof(REAL*));
    for(j=0;j<Framework[i].NumberOfFrameworks;j++)
       SurfaceAreaFrameworksAccumulated[i][j]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    SurfaceAreaCationsAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    SurfaceAreaCount[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));

    NumberOfIntegerMoleculesPerComponentAccumulated[i]=(REAL**)calloc(NumberOfComponents,sizeof(REAL*));
    NumberOfFractionalMoleculesPerComponentAccumulated[i]=(REAL**)calloc(NumberOfComponents,sizeof(REAL*));
    NumberOfExcessMoleculesPerComponentAccumulated[i]=(REAL**)calloc(NumberOfComponents,sizeof(REAL*));
    DensityPerComponentAccumulated[i]=(REAL**)calloc(NumberOfComponents,sizeof(REAL*));

    TotalEnergyTimesNumberOfMoleculesAccumulated[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));

    TotalEnergyTimesNumberOfMoleculesPerComponentAccumulated[i]=(REAL**)calloc(NumberOfComponents,sizeof(REAL));
    HostAdsorbateEnergyTimesNumberOfMoleculesAccumulated[i]=(REAL**)calloc(NumberOfComponents,sizeof(REAL));
    AdsorbateAdsorbateEnergyTimesNumberOfMoleculesAccumulated[i]=(REAL**)calloc(NumberOfComponents,sizeof(REAL));

    WidomRosenbluthFactorAccumulated[i]=(REAL**)calloc(NumberOfComponents,sizeof(REAL*));
    WidomIdealGasAccumulated[i]=(REAL**)calloc(NumberOfComponents,sizeof(REAL*));
    WidomRosenbluthFactorCount[i]=(REAL**)calloc(NumberOfComponents,sizeof(REAL*));
    WidomEnergyDifferenceAccumulated[i]=(REAL**)calloc(NumberOfComponents,sizeof(REAL*));
    WidomEnergyFrameworkAccumulated[i]=(REAL**)calloc(NumberOfComponents,sizeof(REAL*));
    WidomEnergyFrameworkCount[i]=(REAL**)calloc(NumberOfComponents,sizeof(REAL*));

    GibbsWidomRosenbluthFactorAccumulated[i]=(REAL**)calloc(NumberOfComponents,sizeof(REAL*));
    GibbsWidomIdealGasAccumulated[i]=(REAL**)calloc(NumberOfComponents,sizeof(REAL*));
    GibbsWidomRosenbluthFactorCount[i]=(REAL**)calloc(NumberOfComponents,sizeof(REAL*));
    GibbsWidomEnergyDifferenceAccumulated[i]=(REAL**)calloc(NumberOfComponents,sizeof(REAL*));
    GibbsWidomEnergyFrameworkAccumulated[i]=(REAL**)calloc(NumberOfComponents,sizeof(REAL*));
    GibbsWidomEnergyFrameworkCount[i]=(REAL**)calloc(NumberOfComponents,sizeof(REAL*));
    GibbsInverseDensityPerComponentAccumulated[i]=(REAL**)calloc(NumberOfComponents,sizeof(REAL*));

    PrincipleMomentsOfInertiaAccumulated[i]=(VECTOR**)calloc(NumberOfComponents,sizeof(VECTOR*));
    PrincipleMomentsOfInertiaCount[i]=(REAL**)calloc(NumberOfComponents,sizeof(REAL*));

    NumberOfMoleculesPerComponentSquaredAccumulated[i]=(REAL***)calloc(NumberOfComponents,sizeof(REAL**));

    for(j=0;j<NumberOfComponents;j++)
    {
      TotalEnergyTimesNumberOfMoleculesPerComponentAccumulated[i][j]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
      HostAdsorbateEnergyTimesNumberOfMoleculesAccumulated[i][j]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
      AdsorbateAdsorbateEnergyTimesNumberOfMoleculesAccumulated[i][j]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));

      NumberOfIntegerMoleculesPerComponentAccumulated[i][j]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
      NumberOfFractionalMoleculesPerComponentAccumulated[i][j]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
      NumberOfExcessMoleculesPerComponentAccumulated[i][j]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
      DensityPerComponentAccumulated[i][j]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));

      WidomRosenbluthFactorAccumulated[i][j]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
      WidomIdealGasAccumulated[i][j]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
      WidomRosenbluthFactorCount[i][j]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
      WidomEnergyDifferenceAccumulated[i][j]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
      WidomEnergyFrameworkAccumulated[i][j]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
      WidomEnergyFrameworkCount[i][j]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));

      GibbsWidomRosenbluthFactorAccumulated[i][j]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
      GibbsWidomIdealGasAccumulated[i][j]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
      GibbsWidomRosenbluthFactorCount[i][j]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
      GibbsWidomEnergyDifferenceAccumulated[i][j]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
      GibbsWidomEnergyFrameworkAccumulated[i][j]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
      GibbsWidomEnergyFrameworkCount[i][j]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
      GibbsInverseDensityPerComponentAccumulated[i][j]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));

      PrincipleMomentsOfInertiaAccumulated[i][j]=(VECTOR*)calloc(NumberOfBlocks,sizeof(VECTOR));
      PrincipleMomentsOfInertiaCount[i][j]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));

      NumberOfMoleculesPerComponentSquaredAccumulated[i][j]=(REAL**)calloc(NumberOfComponents,sizeof(REAL*));

      for(k=0;k<NumberOfComponents;k++)
        NumberOfMoleculesPerComponentSquaredAccumulated[i][j][k]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));

    }
  }
}

void ReadRestartStatistics(FILE *FilePtr)
{
  int i,j,k;
  int NumberOfBlocks;
  REAL Check;
  int readversionNumber=0;

  NumberOfBlocks=NR_BLOCKS;

  AllocateStatisticsMemory();


  fread(&readversionNumber,sizeof(int),1,FilePtr);
  if(readversionNumber > versionNumber)
  {
    fprintf(stderr,"Upgrade to last version of RASPA to read binary restart-file");
    exit(-1);
  }

  fread(&NumberOfBlocks,sizeof(int),1,FilePtr);
  fread(&Block,sizeof(int),1,FilePtr);

  fread(BlockCycle,sizeof(long long),NumberOfBlocks,FilePtr);

  for(i=0;i<NumberOfSystems;i++)
  {

    fread(BlockWeightedCount[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(BlockCount[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    fread(UHostHostAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UAdsorbateAdsorbateAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UCationCationAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UHostAdsorbateAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UHostCationAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UAdsorbateCationAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    fread(UHostHostVDWAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UAdsorbateAdsorbateVDWAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UCationCationVDWAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UHostAdsorbateVDWAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UHostCationVDWAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UAdsorbateCationVDWAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    fread(UHostHostCoulombAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UAdsorbateAdsorbateCoulombAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UCationCationCoulombAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UHostAdsorbateCoulombAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UHostCationCoulombAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UAdsorbateCationCoulombAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    fread(UHostBondAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UHostUreyBradleyAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UHostBendAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UHostInversionBendAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UHostTorsionAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UHostImproperTorsionAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UHostBondBondAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UHostBendBendAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UHostBondBendAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UHostBondTorsionAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UHostBendTorsionAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    fread(UAdsorbateBondAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UAdsorbateUreyBradleyAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UAdsorbateBondBondAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UAdsorbateBendAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UAdsorbateInversionBendAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UAdsorbateTorsionAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UAdsorbateImproperTorsionAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UAdsorbateBondBondAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UAdsorbateBendBendAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UAdsorbateBondBendAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UAdsorbateBondTorsionAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UAdsorbateBendTorsionAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UAdsorbateIntraVDWAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UAdsorbateIntraChargeChargeAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UAdsorbateIntraChargeBondDipoleAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UAdsorbateIntraBondDipoleBondDipoleAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    fread(UCationBondAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UCationUreyBradleyAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UCationBendAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UCationInversionBendAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UCationTorsionAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UCationImproperTorsionAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UCationBondBondAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UCationBendBendAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UCationBondBendAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UCationBondTorsionAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UCationBendTorsionAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UCationIntraVDWAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UCationIntraChargeChargeAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UCationIntraChargeBondDipoleAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UCationIntraBondDipoleBondDipoleAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    fread(UHostPolarizationAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UAdsorbatePolarizationAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UCationPolarizationAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UHostBackPolarizationAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UAdsorbateBackPolarizationAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UCationBackPolarizationAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    fread(UTailCorrectionAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UDistanceConstraintsAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UAngleConstraintsAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UDihedralConstraintsAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UInversionBendConstraintsAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UOutOfPlaneDistanceConstraintsAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UExclusionConstraintsAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UTotalAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    fread(TotalSystemDipoleAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(TotalSystemDipoleSquaredAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(TotalSystemNormDipoleAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(TotalSystemNormDipoleSquaredAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    fread(NumberOfIntegerMoleculesAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(NumberOfMoleculesSquaredAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(DensityAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(InverseDensityAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(GibbsInverseDensityAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    fread(BoxAccumulated[i],sizeof(VECTOR),NumberOfBlocks,FilePtr);
    fread(BoxAXAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(BoxAYAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(BoxAZAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(BoxBXAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(BoxBYAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(BoxBZAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(BoxCXAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(BoxCYAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(BoxCZAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(BoxLengthAccumulated[i],sizeof(VECTOR),NumberOfBlocks,FilePtr);
    fread(AlphaAngleAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(BetaAngleAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(GammaAngleAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(VolumeAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(VolumeSquaredAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    fread(TotalEnergyAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(TotalEnergySquaredAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(EnthalpyAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(EnthalpySquaredAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(EnthalpyTimesVolumeAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(EnthalpyTimesEnergyAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    fread(TemperatureAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(TemperatureCellAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(TemperatureTranslationAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(TemperatureRotationAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(TemperatureRotationAdsorbateAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(TemperatureTranslationAdsorbateAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    fread(TemperatureAdsorbatesAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(TemperatureCationsAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(TemperatureFrameworkAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    fread(MolecularPressureAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(MolecularStressTensorAccumulated[i],NumberOfBlocks,sizeof(REAL_MATRIX3x3),FilePtr);
    fread(PressureIdealGasPartAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(PressureExcessPartAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(PressureTailCorrectionAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(PressureAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UNoseHooverAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(HeatOfVaporization[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(EnergyPerMolecule[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(VolumePerMolecule[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(CompressibilityAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    fread(BornTermAccumulated[i],sizeof(REAL_MATRIX9x9),NumberOfBlocks,FilePtr);
    fread(ConfigurationalStressFluctuationTermAccumulated[i],sizeof(REAL_MATRIX9x9),NumberOfBlocks,FilePtr);
    fread(ConfigurationalStressTensorAccumulated[i],sizeof(REAL_MATRIX3x3),NumberOfBlocks,FilePtr);
    fread(StressTensorAccumulated[i],sizeof(REAL_MATRIX3x3),NumberOfBlocks,FilePtr);

    fread(SurfaceAreaFrameworkAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    for(j=0;j<Framework[i].NumberOfFrameworks;j++)
      fread(SurfaceAreaFrameworksAccumulated[i][j],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(SurfaceAreaCationsAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(SurfaceAreaCount[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    fread(TotalEnergyTimesNumberOfMoleculesAccumulated[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    for(j=0;j<NumberOfComponents;j++)
    {
      fread(TotalEnergyTimesNumberOfMoleculesPerComponentAccumulated[i][j],sizeof(REAL),NumberOfBlocks,FilePtr);
      fread(HostAdsorbateEnergyTimesNumberOfMoleculesAccumulated[i][j],sizeof(REAL),NumberOfBlocks,FilePtr);
      fread(AdsorbateAdsorbateEnergyTimesNumberOfMoleculesAccumulated[i][j],sizeof(REAL),NumberOfBlocks,FilePtr);

      fread(NumberOfIntegerMoleculesPerComponentAccumulated[i][j],sizeof(REAL),NumberOfBlocks,FilePtr);
      fread(NumberOfFractionalMoleculesPerComponentAccumulated[i][j],sizeof(REAL),NumberOfBlocks,FilePtr);
      fread(NumberOfExcessMoleculesPerComponentAccumulated[i][j],sizeof(REAL),NumberOfBlocks,FilePtr);
      fread(DensityPerComponentAccumulated[i][j],sizeof(REAL),NumberOfBlocks,FilePtr);

      fread(WidomRosenbluthFactorAccumulated[i][j],sizeof(REAL),NumberOfBlocks,FilePtr);
      fread(WidomIdealGasAccumulated[i][j],sizeof(REAL),NumberOfBlocks,FilePtr);
      fread(WidomRosenbluthFactorCount[i][j],sizeof(REAL),NumberOfBlocks,FilePtr);
      fread(WidomEnergyDifferenceAccumulated[i][j],sizeof(REAL),NumberOfBlocks,FilePtr);
      fread(WidomEnergyFrameworkAccumulated[i][j],sizeof(REAL),NumberOfBlocks,FilePtr);
      fread(WidomEnergyFrameworkCount[i][j],sizeof(REAL),NumberOfBlocks,FilePtr);

      fread(GibbsWidomRosenbluthFactorAccumulated[i][j],sizeof(REAL),NumberOfBlocks,FilePtr);
      fread(GibbsWidomIdealGasAccumulated[i][j],sizeof(REAL),NumberOfBlocks,FilePtr);
      fread(GibbsWidomRosenbluthFactorCount[i][j],sizeof(REAL),NumberOfBlocks,FilePtr);
      fread(GibbsWidomEnergyDifferenceAccumulated[i][j],sizeof(REAL),NumberOfBlocks,FilePtr);
      fread(GibbsWidomEnergyFrameworkAccumulated[i][j],sizeof(REAL),NumberOfBlocks,FilePtr);
      fread(GibbsWidomEnergyFrameworkCount[i][j],sizeof(REAL),NumberOfBlocks,FilePtr);

      fread(PrincipleMomentsOfInertiaAccumulated[i][j],sizeof(REAL),NumberOfBlocks,FilePtr);
      fread(PrincipleMomentsOfInertiaCount[i][j],sizeof(REAL),NumberOfBlocks,FilePtr);

      for(k=0;k<NumberOfComponents;k++)
        fread(NumberOfMoleculesPerComponentSquaredAccumulated[i][j][k],sizeof(REAL),NumberOfBlocks,FilePtr);
    }
  }

  fread(&Check,1,sizeof(REAL),FilePtr);
  if(fabs(Check-123456789.0)>1e-10)
  {
    fprintf(stderr, "Error in binary restart-file (ReadRestartStatistics)\n");
    ContinueAfterCrash=FALSE;
  }
}


