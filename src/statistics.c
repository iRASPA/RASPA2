/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2015 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'statistics.c' is part of RASPA-2.0

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
#endif
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
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
REAL **BlockCount;
int Block;

// average energies
static REAL **UHostHostAverage;
static REAL **UAdsorbateAdsorbateAverage;
static REAL **UCationCationAverage;
static REAL **UHostAdsorbateAverage;
static REAL **UHostCationAverage;
static REAL **UAdsorbateCationAverage;

static REAL **UHostHostVDWAverage;
static REAL **UAdsorbateAdsorbateVDWAverage;
static REAL **UCationCationVDWAverage;
static REAL **UHostAdsorbateVDWAverage;
static REAL **UHostCationVDWAverage;
static REAL **UAdsorbateCationVDWAverage;

static REAL **UHostHostCoulombAverage;
static REAL **UAdsorbateAdsorbateCoulombAverage;
static REAL **UCationCationCoulombAverage;
static REAL **UHostAdsorbateCoulombAverage;
static REAL **UHostCationCoulombAverage;
static REAL **UAdsorbateCationCoulombAverage;

static REAL **UHostBondAverage;
static REAL **UHostUreyBradleyAverage;
static REAL **UHostBendAverage;
static REAL **UHostInversionBendAverage;
static REAL **UHostTorsionAverage;
static REAL **UHostImproperTorsionAverage;
static REAL **UHostBondBondAverage;
static REAL **UHostBendBendAverage;
static REAL **UHostBondBendAverage;
static REAL **UHostBondTorsionAverage;
static REAL **UHostBendTorsionAverage;

static REAL **UAdsorbateBondAverage;
static REAL **UAdsorbateUreyBradleyAverage;
static REAL **UAdsorbateBondBondAverage;
static REAL **UAdsorbateBendAverage;
static REAL **UAdsorbateInversionBendAverage;
static REAL **UAdsorbateTorsionAverage;
static REAL **UAdsorbateImproperTorsionAverage;
static REAL **UAdsorbateBondBondAverage;
static REAL **UAdsorbateBendBendAverage;
static REAL **UAdsorbateBondBendAverage;
static REAL **UAdsorbateBondTorsionAverage;
static REAL **UAdsorbateBendTorsionAverage;
static REAL **UAdsorbateIntraVDWAverage;
static REAL **UAdsorbateIntraChargeChargeAverage;
static REAL **UAdsorbateIntraChargeBondDipoleAverage;
static REAL **UAdsorbateIntraBondDipoleBondDipoleAverage;

static REAL **UCationBondAverage;
static REAL **UCationUreyBradleyAverage;
static REAL **UCationBendAverage;
static REAL **UCationInversionBendAverage;
static REAL **UCationTorsionAverage;
static REAL **UCationImproperTorsionAverage;
static REAL **UCationBondBondAverage;
static REAL **UCationBendBendAverage;
static REAL **UCationBondBendAverage;
static REAL **UCationBondTorsionAverage;
static REAL **UCationBendTorsionAverage;
static REAL **UCationIntraVDWAverage;
static REAL **UCationIntraChargeChargeAverage;
static REAL **UCationIntraChargeBondDipoleAverage;
static REAL **UCationIntraBondDipoleBondDipoleAverage;

static REAL **UHostPolarizationAverage;
static REAL **UAdsorbatePolarizationAverage;
static REAL **UCationPolarizationAverage;
static REAL **UHostBackPolarizationAverage;
static REAL **UAdsorbateBackPolarizationAverage;
static REAL **UCationBackPolarizationAverage;

static REAL **UTailCorrectionAverage;

static REAL **UDistanceConstraintsAverage;
static REAL **UAngleConstraintsAverage;
static REAL **UDihedralConstraintsAverage;
static REAL **UInversionBendConstraintsAverage;
static REAL **UOutOfPlaneDistanceConstraintsAverage;
static REAL **UExclusionConstraintsAverage;

static REAL **UTotalAverage;

static REAL ***NumberOfMoleculesPerComponentAverage;
static REAL ***NumberOfExcessMoleculesPerComponentAverage;
static REAL ***DensityPerComponentAverage;
static REAL **TotalEnergyTimesNumberOfMoleculesAverage;
static REAL **HostAdsorbateEnergyTimesNumberOfMoleculesAverage;
static REAL **AdsorbateAdsorbateEnergyTimesNumberOfMoleculesAverage;

static VECTOR **TotalSystemDipoleAverage;
static VECTOR **TotalSystemDipoleSquaredAverage;
static REAL **TotalSystemNormDipoleAverage;
static REAL **TotalSystemNormDipoleSquaredAverage;

static REAL **NumberOfMoleculesAverage;
static REAL **NumberOfMoleculesSquaredAverage;
static REAL **DensityAverage;

static VECTOR **BoxAverage;
static REAL **BoxAverageAX;
static REAL **BoxAverageAY;
static REAL **BoxAverageAZ;
static REAL **BoxAverageBX;
static REAL **BoxAverageBY;
static REAL **BoxAverageBZ;
static REAL **BoxAverageCX;
static REAL **BoxAverageCY;
static REAL **BoxAverageCZ;
static VECTOR **BoxLengthAverage;
static REAL **AlphaAngleAverage;
static REAL **BetaAngleAverage;
static REAL **GammaAngleAverage;
static REAL **VolumeAverage;
static REAL **VolumeSquaredAverage;

static REAL **TotalEnergyAverage;
static REAL **TotalEnergySquaredAverage;
static REAL **EnthalpyAverage;
static REAL **EnthalpySquaredAverage;
static REAL **EnthalpyTimesVolumeAverage;
static REAL **EnthalpyTimesEnergyAverage;

static REAL **TemperatureAverage;
static REAL **TemperatureCellAverage;
static REAL **TemperatureTranslationAverage;
static REAL **TemperatureRotationAverage;
static REAL **TemperatureRotationAdsorbateAverage;
static REAL **TemperatureTranslationAdsorbateAverage;

static REAL **TemperatureAdsorbatesAverage;
static REAL **TemperatureCationsAverage;
static REAL **TemperatureFrameworkAverage;

static REAL **MolecularPressureAverage;
static REAL_MATRIX3x3 **MolecularStressTensorAverage;
static REAL **PressureIdealGasPartAverage;
static REAL **PressureExcessPartAverage;
static REAL **PressureTailCorrectionAverage;
static REAL **PressureAverage;
static REAL **UNoseHooverAverage;

static REAL **HeatOfVaporization;
static REAL **EnergyPerMolecule;
static REAL **VolumePerMolecule;
static REAL **CompressibilityAverage;

static REAL_MATRIX9x9 **BornTermAverage;
static REAL_MATRIX9x9 **ConfigurationalStressFluctuationTermAverage;
static REAL_MATRIX3x3 **ConfigurationalStressTensorAverage;
static REAL_MATRIX3x3 **StressTensorAverage;

REAL ***WidomRosenbluthFactorAverage;
REAL ***WidomRosenbluthFactorCount;

REAL ***WidomEnergyDifferenceAverage;

REAL ***WidomEnergyFrameworkAverage;
REAL ***WidomEnergyFrameworkCount;

REAL **SurfaceAreaFrameworkAverage;
REAL ***SurfaceAreaFrameworksAverage;
REAL **SurfaceAreaCationsAverage;
REAL **SurfaceAreaCount;

VECTOR ***PrincipleMomentsOfInertiaAverage;
REAL ***PrincipleMomentsOfInertiaCount;

void AddBornTermToAverages(void)
{

  // BornTerm
  BornTermAverage[CurrentSystem][Block].xxxx+=BornTerm[CurrentSystem].xxxx/Volume[CurrentSystem];
  BornTermAverage[CurrentSystem][Block].yxxx+=BornTerm[CurrentSystem].yxxx/Volume[CurrentSystem];
  BornTermAverage[CurrentSystem][Block].zxxx+=BornTerm[CurrentSystem].zxxx/Volume[CurrentSystem];
  BornTermAverage[CurrentSystem][Block].xxyx+=BornTerm[CurrentSystem].xxyx/Volume[CurrentSystem];
  BornTermAverage[CurrentSystem][Block].yxyx+=BornTerm[CurrentSystem].yxyx/Volume[CurrentSystem];
  BornTermAverage[CurrentSystem][Block].zxyx+=BornTerm[CurrentSystem].zxyx/Volume[CurrentSystem];
  BornTermAverage[CurrentSystem][Block].xxzx+=BornTerm[CurrentSystem].xxzx/Volume[CurrentSystem];
  BornTermAverage[CurrentSystem][Block].yxzx+=BornTerm[CurrentSystem].yxzx/Volume[CurrentSystem];
  BornTermAverage[CurrentSystem][Block].zxzx+=BornTerm[CurrentSystem].zxzx/Volume[CurrentSystem];

  BornTermAverage[CurrentSystem][Block].xyxx+=BornTerm[CurrentSystem].xyxx/Volume[CurrentSystem];
  BornTermAverage[CurrentSystem][Block].yyxx+=BornTerm[CurrentSystem].yyxx/Volume[CurrentSystem];
  BornTermAverage[CurrentSystem][Block].zyxx+=BornTerm[CurrentSystem].zyxx/Volume[CurrentSystem];
  BornTermAverage[CurrentSystem][Block].xyyx+=BornTerm[CurrentSystem].xyyx/Volume[CurrentSystem];
  BornTermAverage[CurrentSystem][Block].yyyx+=BornTerm[CurrentSystem].yyyx/Volume[CurrentSystem];
  BornTermAverage[CurrentSystem][Block].zyyx+=BornTerm[CurrentSystem].zyyx/Volume[CurrentSystem];
  BornTermAverage[CurrentSystem][Block].xyzx+=BornTerm[CurrentSystem].xyzx/Volume[CurrentSystem];
  BornTermAverage[CurrentSystem][Block].yyzx+=BornTerm[CurrentSystem].yyzx/Volume[CurrentSystem];
  BornTermAverage[CurrentSystem][Block].zyzx+=BornTerm[CurrentSystem].zyzx/Volume[CurrentSystem];

  BornTermAverage[CurrentSystem][Block].xzxx+=BornTerm[CurrentSystem].xzxx/Volume[CurrentSystem];
  BornTermAverage[CurrentSystem][Block].yzxx+=BornTerm[CurrentSystem].yzxx/Volume[CurrentSystem];
  BornTermAverage[CurrentSystem][Block].zzxx+=BornTerm[CurrentSystem].zzxx/Volume[CurrentSystem];
  BornTermAverage[CurrentSystem][Block].xzyx+=BornTerm[CurrentSystem].xzyx/Volume[CurrentSystem];
  BornTermAverage[CurrentSystem][Block].yzyx+=BornTerm[CurrentSystem].yzyx/Volume[CurrentSystem];
  BornTermAverage[CurrentSystem][Block].zzyx+=BornTerm[CurrentSystem].zzyx/Volume[CurrentSystem];
  BornTermAverage[CurrentSystem][Block].xzzx+=BornTerm[CurrentSystem].xzzx/Volume[CurrentSystem];
  BornTermAverage[CurrentSystem][Block].yzzx+=BornTerm[CurrentSystem].yzzx/Volume[CurrentSystem];
  BornTermAverage[CurrentSystem][Block].zzzx+=BornTerm[CurrentSystem].zzzx/Volume[CurrentSystem];

  BornTermAverage[CurrentSystem][Block].xxxy+=BornTerm[CurrentSystem].xxxy/Volume[CurrentSystem];
  BornTermAverage[CurrentSystem][Block].yxxy+=BornTerm[CurrentSystem].yxxy/Volume[CurrentSystem];
  BornTermAverage[CurrentSystem][Block].zxxy+=BornTerm[CurrentSystem].zxxy/Volume[CurrentSystem];
  BornTermAverage[CurrentSystem][Block].xxyy+=BornTerm[CurrentSystem].xxyy/Volume[CurrentSystem];
  BornTermAverage[CurrentSystem][Block].yxyy+=BornTerm[CurrentSystem].yxyy/Volume[CurrentSystem];
  BornTermAverage[CurrentSystem][Block].zxyy+=BornTerm[CurrentSystem].zxyy/Volume[CurrentSystem];
  BornTermAverage[CurrentSystem][Block].xxzy+=BornTerm[CurrentSystem].xxzy/Volume[CurrentSystem];
  BornTermAverage[CurrentSystem][Block].yxzy+=BornTerm[CurrentSystem].yxzy/Volume[CurrentSystem];
  BornTermAverage[CurrentSystem][Block].zxzy+=BornTerm[CurrentSystem].zxzy/Volume[CurrentSystem];

  BornTermAverage[CurrentSystem][Block].xyxy+=BornTerm[CurrentSystem].xyxy/Volume[CurrentSystem];
  BornTermAverage[CurrentSystem][Block].yyxy+=BornTerm[CurrentSystem].yyxy/Volume[CurrentSystem];
  BornTermAverage[CurrentSystem][Block].zyxy+=BornTerm[CurrentSystem].zyxy/Volume[CurrentSystem];
  BornTermAverage[CurrentSystem][Block].xyyy+=BornTerm[CurrentSystem].xyyy/Volume[CurrentSystem];
  BornTermAverage[CurrentSystem][Block].yyyy+=BornTerm[CurrentSystem].yyyy/Volume[CurrentSystem];
  BornTermAverage[CurrentSystem][Block].zyyy+=BornTerm[CurrentSystem].zyyy/Volume[CurrentSystem];
  BornTermAverage[CurrentSystem][Block].xyzy+=BornTerm[CurrentSystem].xyzy/Volume[CurrentSystem];
  BornTermAverage[CurrentSystem][Block].yyzy+=BornTerm[CurrentSystem].yyzy/Volume[CurrentSystem];
  BornTermAverage[CurrentSystem][Block].zyzy+=BornTerm[CurrentSystem].zyzy/Volume[CurrentSystem];

  BornTermAverage[CurrentSystem][Block].xzxy+=BornTerm[CurrentSystem].xzxy/Volume[CurrentSystem];
  BornTermAverage[CurrentSystem][Block].yzxy+=BornTerm[CurrentSystem].yzxy/Volume[CurrentSystem];
  BornTermAverage[CurrentSystem][Block].zzxy+=BornTerm[CurrentSystem].zzxy/Volume[CurrentSystem];
  BornTermAverage[CurrentSystem][Block].xzyy+=BornTerm[CurrentSystem].xzyy/Volume[CurrentSystem];
  BornTermAverage[CurrentSystem][Block].yzyy+=BornTerm[CurrentSystem].yzyy/Volume[CurrentSystem];
  BornTermAverage[CurrentSystem][Block].zzyy+=BornTerm[CurrentSystem].zzyy/Volume[CurrentSystem];
  BornTermAverage[CurrentSystem][Block].xzzy+=BornTerm[CurrentSystem].xzzy/Volume[CurrentSystem];
  BornTermAverage[CurrentSystem][Block].yzzy+=BornTerm[CurrentSystem].yzzy/Volume[CurrentSystem];
  BornTermAverage[CurrentSystem][Block].zzzy+=BornTerm[CurrentSystem].zzzy/Volume[CurrentSystem];

  BornTermAverage[CurrentSystem][Block].xxxz+=BornTerm[CurrentSystem].xxxz/Volume[CurrentSystem];
  BornTermAverage[CurrentSystem][Block].yxxz+=BornTerm[CurrentSystem].yxxz/Volume[CurrentSystem];
  BornTermAverage[CurrentSystem][Block].zxxz+=BornTerm[CurrentSystem].zxxz/Volume[CurrentSystem];
  BornTermAverage[CurrentSystem][Block].xxyz+=BornTerm[CurrentSystem].xxyz/Volume[CurrentSystem];
  BornTermAverage[CurrentSystem][Block].yxyz+=BornTerm[CurrentSystem].yxyz/Volume[CurrentSystem];
  BornTermAverage[CurrentSystem][Block].zxyz+=BornTerm[CurrentSystem].zxyz/Volume[CurrentSystem];
  BornTermAverage[CurrentSystem][Block].xxzz+=BornTerm[CurrentSystem].xxzz/Volume[CurrentSystem];
  BornTermAverage[CurrentSystem][Block].yxzz+=BornTerm[CurrentSystem].yxzz/Volume[CurrentSystem];
  BornTermAverage[CurrentSystem][Block].zxzz+=BornTerm[CurrentSystem].zxzz/Volume[CurrentSystem];

  BornTermAverage[CurrentSystem][Block].xyxz+=BornTerm[CurrentSystem].xyxz/Volume[CurrentSystem];
  BornTermAverage[CurrentSystem][Block].yyxz+=BornTerm[CurrentSystem].yyxz/Volume[CurrentSystem];
  BornTermAverage[CurrentSystem][Block].zyxz+=BornTerm[CurrentSystem].zyxz/Volume[CurrentSystem];
  BornTermAverage[CurrentSystem][Block].xyyz+=BornTerm[CurrentSystem].xyyz/Volume[CurrentSystem];
  BornTermAverage[CurrentSystem][Block].yyyz+=BornTerm[CurrentSystem].yyyz/Volume[CurrentSystem];
  BornTermAverage[CurrentSystem][Block].zyyz+=BornTerm[CurrentSystem].zyyz/Volume[CurrentSystem];
  BornTermAverage[CurrentSystem][Block].xyzz+=BornTerm[CurrentSystem].xyzz/Volume[CurrentSystem];
  BornTermAverage[CurrentSystem][Block].yyzz+=BornTerm[CurrentSystem].yyzz/Volume[CurrentSystem];
  BornTermAverage[CurrentSystem][Block].zyzz+=BornTerm[CurrentSystem].zyzz/Volume[CurrentSystem];

  BornTermAverage[CurrentSystem][Block].xzxz+=BornTerm[CurrentSystem].xzxz/Volume[CurrentSystem];
  BornTermAverage[CurrentSystem][Block].yzxz+=BornTerm[CurrentSystem].yzxz/Volume[CurrentSystem];
  BornTermAverage[CurrentSystem][Block].zzxz+=BornTerm[CurrentSystem].zzxz/Volume[CurrentSystem];
  BornTermAverage[CurrentSystem][Block].xzyz+=BornTerm[CurrentSystem].xzyz/Volume[CurrentSystem];
  BornTermAverage[CurrentSystem][Block].yzyz+=BornTerm[CurrentSystem].yzyz/Volume[CurrentSystem];
  BornTermAverage[CurrentSystem][Block].zzyz+=BornTerm[CurrentSystem].zzyz/Volume[CurrentSystem];
  BornTermAverage[CurrentSystem][Block].xzzz+=BornTerm[CurrentSystem].xzzz/Volume[CurrentSystem];
  BornTermAverage[CurrentSystem][Block].yzzz+=BornTerm[CurrentSystem].yzzz/Volume[CurrentSystem];
  BornTermAverage[CurrentSystem][Block].zzzz+=BornTerm[CurrentSystem].zzzz/Volume[CurrentSystem];

  // stress fluctuations
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].xxxx+=ConfigurationalStressTensor[CurrentSystem].ax*ConfigurationalStressTensor[CurrentSystem].ax;
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].yxxx+=ConfigurationalStressTensor[CurrentSystem].bx*ConfigurationalStressTensor[CurrentSystem].ax;
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].zxxx+=ConfigurationalStressTensor[CurrentSystem].cx*ConfigurationalStressTensor[CurrentSystem].ax;
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].xxyx+=ConfigurationalStressTensor[CurrentSystem].ax*ConfigurationalStressTensor[CurrentSystem].bx;
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].yxyx+=ConfigurationalStressTensor[CurrentSystem].bx*ConfigurationalStressTensor[CurrentSystem].bx;
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].zxyx+=ConfigurationalStressTensor[CurrentSystem].cx*ConfigurationalStressTensor[CurrentSystem].bx;
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].xxzx+=ConfigurationalStressTensor[CurrentSystem].ax*ConfigurationalStressTensor[CurrentSystem].cx;
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].yxzx+=ConfigurationalStressTensor[CurrentSystem].bx*ConfigurationalStressTensor[CurrentSystem].cx;
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].zxzx+=ConfigurationalStressTensor[CurrentSystem].cx*ConfigurationalStressTensor[CurrentSystem].cx;

  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].xyxx+=ConfigurationalStressTensor[CurrentSystem].ay*ConfigurationalStressTensor[CurrentSystem].ax;
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].yyxx+=ConfigurationalStressTensor[CurrentSystem].by*ConfigurationalStressTensor[CurrentSystem].ax;
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].zyxx+=ConfigurationalStressTensor[CurrentSystem].cy*ConfigurationalStressTensor[CurrentSystem].ax;
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].xyyx+=ConfigurationalStressTensor[CurrentSystem].ay*ConfigurationalStressTensor[CurrentSystem].bx;
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].yyyx+=ConfigurationalStressTensor[CurrentSystem].by*ConfigurationalStressTensor[CurrentSystem].bx;
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].zyyx+=ConfigurationalStressTensor[CurrentSystem].cy*ConfigurationalStressTensor[CurrentSystem].bx;
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].xyzx+=ConfigurationalStressTensor[CurrentSystem].ay*ConfigurationalStressTensor[CurrentSystem].cx;
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].yyzx+=ConfigurationalStressTensor[CurrentSystem].by*ConfigurationalStressTensor[CurrentSystem].cx;
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].zyzx+=ConfigurationalStressTensor[CurrentSystem].cy*ConfigurationalStressTensor[CurrentSystem].cx;

  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].xzxx+=ConfigurationalStressTensor[CurrentSystem].az*ConfigurationalStressTensor[CurrentSystem].ax;
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].yzxx+=ConfigurationalStressTensor[CurrentSystem].bz*ConfigurationalStressTensor[CurrentSystem].ax;
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].zzxx+=ConfigurationalStressTensor[CurrentSystem].cz*ConfigurationalStressTensor[CurrentSystem].ax;
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].xzyx+=ConfigurationalStressTensor[CurrentSystem].az*ConfigurationalStressTensor[CurrentSystem].bx;
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].yzyx+=ConfigurationalStressTensor[CurrentSystem].bz*ConfigurationalStressTensor[CurrentSystem].bx;
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].zzyx+=ConfigurationalStressTensor[CurrentSystem].cz*ConfigurationalStressTensor[CurrentSystem].bx;
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].xzzx+=ConfigurationalStressTensor[CurrentSystem].az*ConfigurationalStressTensor[CurrentSystem].cx;
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].yzzx+=ConfigurationalStressTensor[CurrentSystem].bz*ConfigurationalStressTensor[CurrentSystem].cx;
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].zzzx+=ConfigurationalStressTensor[CurrentSystem].cz*ConfigurationalStressTensor[CurrentSystem].cx;


  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].xxxy+=ConfigurationalStressTensor[CurrentSystem].ax*ConfigurationalStressTensor[CurrentSystem].ay;
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].yxxy+=ConfigurationalStressTensor[CurrentSystem].bx*ConfigurationalStressTensor[CurrentSystem].ay;
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].zxxy+=ConfigurationalStressTensor[CurrentSystem].cx*ConfigurationalStressTensor[CurrentSystem].ay;
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].xxyy+=ConfigurationalStressTensor[CurrentSystem].ax*ConfigurationalStressTensor[CurrentSystem].by;
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].yxyy+=ConfigurationalStressTensor[CurrentSystem].bx*ConfigurationalStressTensor[CurrentSystem].by;
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].zxyy+=ConfigurationalStressTensor[CurrentSystem].cx*ConfigurationalStressTensor[CurrentSystem].by;
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].xxzy+=ConfigurationalStressTensor[CurrentSystem].ax*ConfigurationalStressTensor[CurrentSystem].cy;
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].yxzy+=ConfigurationalStressTensor[CurrentSystem].bx*ConfigurationalStressTensor[CurrentSystem].cy;
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].zxzy+=ConfigurationalStressTensor[CurrentSystem].cx*ConfigurationalStressTensor[CurrentSystem].cy;

  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].xyxy+=ConfigurationalStressTensor[CurrentSystem].ay*ConfigurationalStressTensor[CurrentSystem].ay;
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].yyxy+=ConfigurationalStressTensor[CurrentSystem].by*ConfigurationalStressTensor[CurrentSystem].ay;
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].zyxy+=ConfigurationalStressTensor[CurrentSystem].cy*ConfigurationalStressTensor[CurrentSystem].ay;
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].xyyy+=ConfigurationalStressTensor[CurrentSystem].ay*ConfigurationalStressTensor[CurrentSystem].by;
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].yyyy+=ConfigurationalStressTensor[CurrentSystem].by*ConfigurationalStressTensor[CurrentSystem].by;
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].zyyy+=ConfigurationalStressTensor[CurrentSystem].cy*ConfigurationalStressTensor[CurrentSystem].by;
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].xyzy+=ConfigurationalStressTensor[CurrentSystem].ay*ConfigurationalStressTensor[CurrentSystem].cy;
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].yyzy+=ConfigurationalStressTensor[CurrentSystem].by*ConfigurationalStressTensor[CurrentSystem].cy;
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].zyzy+=ConfigurationalStressTensor[CurrentSystem].cy*ConfigurationalStressTensor[CurrentSystem].cy;

  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].xzxy+=ConfigurationalStressTensor[CurrentSystem].az*ConfigurationalStressTensor[CurrentSystem].ay;
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].yzxy+=ConfigurationalStressTensor[CurrentSystem].bz*ConfigurationalStressTensor[CurrentSystem].ay;
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].zzxy+=ConfigurationalStressTensor[CurrentSystem].cz*ConfigurationalStressTensor[CurrentSystem].ay;
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].xzyy+=ConfigurationalStressTensor[CurrentSystem].az*ConfigurationalStressTensor[CurrentSystem].by;
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].yzyy+=ConfigurationalStressTensor[CurrentSystem].bz*ConfigurationalStressTensor[CurrentSystem].by;
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].zzyy+=ConfigurationalStressTensor[CurrentSystem].cz*ConfigurationalStressTensor[CurrentSystem].by;
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].xzzy+=ConfigurationalStressTensor[CurrentSystem].az*ConfigurationalStressTensor[CurrentSystem].cy;
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].yzzy+=ConfigurationalStressTensor[CurrentSystem].bz*ConfigurationalStressTensor[CurrentSystem].cy;
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].zzzy+=ConfigurationalStressTensor[CurrentSystem].cz*ConfigurationalStressTensor[CurrentSystem].cy;


  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].xxxz+=ConfigurationalStressTensor[CurrentSystem].ax*ConfigurationalStressTensor[CurrentSystem].az;
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].yxxz+=ConfigurationalStressTensor[CurrentSystem].bx*ConfigurationalStressTensor[CurrentSystem].az;
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].zxxz+=ConfigurationalStressTensor[CurrentSystem].cx*ConfigurationalStressTensor[CurrentSystem].az;
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].xxyz+=ConfigurationalStressTensor[CurrentSystem].ax*ConfigurationalStressTensor[CurrentSystem].bz;
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].yxyz+=ConfigurationalStressTensor[CurrentSystem].bx*ConfigurationalStressTensor[CurrentSystem].bz;
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].zxyz+=ConfigurationalStressTensor[CurrentSystem].cx*ConfigurationalStressTensor[CurrentSystem].bz;
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].xxzz+=ConfigurationalStressTensor[CurrentSystem].ax*ConfigurationalStressTensor[CurrentSystem].cz;
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].yxzz+=ConfigurationalStressTensor[CurrentSystem].bx*ConfigurationalStressTensor[CurrentSystem].cz;
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].zxzz+=ConfigurationalStressTensor[CurrentSystem].cx*ConfigurationalStressTensor[CurrentSystem].cz;

  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].xyxz+=ConfigurationalStressTensor[CurrentSystem].ay*ConfigurationalStressTensor[CurrentSystem].az;
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].yyxz+=ConfigurationalStressTensor[CurrentSystem].by*ConfigurationalStressTensor[CurrentSystem].az;
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].zyxz+=ConfigurationalStressTensor[CurrentSystem].cy*ConfigurationalStressTensor[CurrentSystem].az;
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].xyyz+=ConfigurationalStressTensor[CurrentSystem].ay*ConfigurationalStressTensor[CurrentSystem].bz;
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].yyyz+=ConfigurationalStressTensor[CurrentSystem].by*ConfigurationalStressTensor[CurrentSystem].bz;
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].zyyz+=ConfigurationalStressTensor[CurrentSystem].cy*ConfigurationalStressTensor[CurrentSystem].bz;
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].xyzz+=ConfigurationalStressTensor[CurrentSystem].ay*ConfigurationalStressTensor[CurrentSystem].cz;
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].yyzz+=ConfigurationalStressTensor[CurrentSystem].by*ConfigurationalStressTensor[CurrentSystem].cz;
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].zyzz+=ConfigurationalStressTensor[CurrentSystem].cy*ConfigurationalStressTensor[CurrentSystem].cz;

  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].xzxz+=ConfigurationalStressTensor[CurrentSystem].az*ConfigurationalStressTensor[CurrentSystem].az;
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].yzxz+=ConfigurationalStressTensor[CurrentSystem].bz*ConfigurationalStressTensor[CurrentSystem].az;
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].zzxz+=ConfigurationalStressTensor[CurrentSystem].cz*ConfigurationalStressTensor[CurrentSystem].az;
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].xzyz+=ConfigurationalStressTensor[CurrentSystem].az*ConfigurationalStressTensor[CurrentSystem].bz;
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].yzyz+=ConfigurationalStressTensor[CurrentSystem].bz*ConfigurationalStressTensor[CurrentSystem].bz;
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].zzyz+=ConfigurationalStressTensor[CurrentSystem].cz*ConfigurationalStressTensor[CurrentSystem].bz;
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].xzzz+=ConfigurationalStressTensor[CurrentSystem].az*ConfigurationalStressTensor[CurrentSystem].cz;
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].yzzz+=ConfigurationalStressTensor[CurrentSystem].bz*ConfigurationalStressTensor[CurrentSystem].cz;
  ConfigurationalStressFluctuationTermAverage[CurrentSystem][Block].zzzz+=ConfigurationalStressTensor[CurrentSystem].cz*ConfigurationalStressTensor[CurrentSystem].cz;
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
    if(BlockCount[CurrentSystem][i]>0.0)
    {
      Average.ax+=StressTensorAverage[CurrentSystem][i].ax;
      Average.ay+=StressTensorAverage[CurrentSystem][i].ay;
      Average.az+=StressTensorAverage[CurrentSystem][i].az;

      Average.bx+=StressTensorAverage[CurrentSystem][i].bx;
      Average.by+=StressTensorAverage[CurrentSystem][i].by;
      Average.bz+=StressTensorAverage[CurrentSystem][i].bz;

      Average.cx+=StressTensorAverage[CurrentSystem][i].cx;
      Average.cy+=StressTensorAverage[CurrentSystem][i].cy;
      Average.cz+=StressTensorAverage[CurrentSystem][i].cz;

      count+=BlockCount[CurrentSystem][i];
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
    if(BlockCount[CurrentSystem][i]>0.0)
    {
      Average.ax+=ConfigurationalStressTensorAverage[CurrentSystem][i].ax;
      Average.ay+=ConfigurationalStressTensorAverage[CurrentSystem][i].ay;
      Average.az+=ConfigurationalStressTensorAverage[CurrentSystem][i].az;

      Average.bx+=ConfigurationalStressTensorAverage[CurrentSystem][i].bx;
      Average.by+=ConfigurationalStressTensorAverage[CurrentSystem][i].by;
      Average.bz+=ConfigurationalStressTensorAverage[CurrentSystem][i].bz;

      Average.cx+=ConfigurationalStressTensorAverage[CurrentSystem][i].cx;
      Average.cy+=ConfigurationalStressTensorAverage[CurrentSystem][i].cy;
      Average.cz+=ConfigurationalStressTensorAverage[CurrentSystem][i].cz;

      count+=BlockCount[CurrentSystem][i];
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
    if(BlockCount[CurrentSystem][i]>0.0)
    {
      sum1+=TemperatureAverage[CurrentSystem][i];
      sum2+=BlockCount[CurrentSystem][i];
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
    if(BlockCount[CurrentSystem][i]>0.0)
    {
      sum1+=TemperatureCellAverage[CurrentSystem][i];
      sum2+=BlockCount[CurrentSystem][i];
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
    if(BlockCount[CurrentSystem][i]>0.0)
    {
      sum+=VolumeAverage[CurrentSystem][i];
      count+=BlockCount[CurrentSystem][i];
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
    if(BlockCount[CurrentSystem][i]>0.0)
    {
      AddRealMatrix9x9(&BornTerm,BornTerm,BornTermAverage[CurrentSystem][i]);
      AddRealMatrix9x9(&StressFluctuationTerm,StressFluctuationTerm,ConfigurationalStressFluctuationTermAverage[CurrentSystem][i]);

      Stress.ax+=ConfigurationalStressTensorAverage[CurrentSystem][i].ax;
      Stress.ay+=ConfigurationalStressTensorAverage[CurrentSystem][i].ay;
      Stress.az+=ConfigurationalStressTensorAverage[CurrentSystem][i].az;

      Stress.bx+=ConfigurationalStressTensorAverage[CurrentSystem][i].bx;
      Stress.by+=ConfigurationalStressTensorAverage[CurrentSystem][i].by;
      Stress.bz+=ConfigurationalStressTensorAverage[CurrentSystem][i].bz;

      Stress.cx+=ConfigurationalStressTensorAverage[CurrentSystem][i].cx;
      Stress.cy+=ConfigurationalStressTensorAverage[CurrentSystem][i].cy;
      Stress.cz+=ConfigurationalStressTensorAverage[CurrentSystem][i].cz;

      count+=BlockCount[CurrentSystem][i];
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
  int i,j,k;

  Block=0;

  for(i=0;i<NR_BLOCKS;i++)
  {
    BlockCycle[i]=(long long)((i+1)*NumberOfCycles/(double)NR_BLOCKS);
    for(j=0;j<NumberOfSystems;j++)
      BlockCount[j][i]=0.0;
  }

  for(k=0;k<NumberOfSystems;k++)
  {
    for(i=0;i<NR_BLOCKS;i++)
    {
      ConfigurationalStressTensorAverage[k][i].ax=0.0; ConfigurationalStressTensorAverage[k][i].bx=0.0; ConfigurationalStressTensorAverage[k][i].cx=0.0;
      ConfigurationalStressTensorAverage[k][i].ay=0.0; ConfigurationalStressTensorAverage[k][i].by=0.0; ConfigurationalStressTensorAverage[k][i].cy=0.0;
      ConfigurationalStressTensorAverage[k][i].az=0.0; ConfigurationalStressTensorAverage[k][i].bz=0.0; ConfigurationalStressTensorAverage[k][i].cz=0.0;

      // initialize energies
      UHostBondAverage[k][i]=0.0;
      UHostUreyBradleyAverage[k][i]=0.0;
      UHostBendAverage[k][i]=0.0;
      UHostInversionBendAverage[k][i]=0.0;
      UHostTorsionAverage[k][i]=0.0;
      UHostImproperTorsionAverage[k][i]=0.0;
      UHostBondBondAverage[k][i]=0.0;
      UHostBendBendAverage[k][i]=0.0;
      UHostBondBendAverage[k][i]=0.0;
      UHostBondTorsionAverage[k][i]=0.0;
      UHostBendTorsionAverage[k][i]=0.0;

      UCationBondAverage[k][i]=0.0;
      UCationUreyBradleyAverage[k][i]=0.0;
      UCationBendAverage[k][i]=0.0;
      UCationInversionBendAverage[k][i]=0.0;
      UCationTorsionAverage[k][i]=0.0;
      UCationImproperTorsionAverage[k][i]=0.0;
      UCationBondBondAverage[k][i]=0.0;
      UCationBendBendAverage[k][i]=0.0;
      UCationBondBendAverage[k][i]=0.0;
      UCationBondTorsionAverage[k][i]=0.0;
      UCationBendTorsionAverage[k][i]=0.0;
      UCationIntraVDWAverage[k][i]=0.0;
      UCationIntraChargeChargeAverage[k][i]=0.0;
      UCationIntraChargeBondDipoleAverage[k][i]=0.0;
      UCationIntraBondDipoleBondDipoleAverage[k][i]=0.0;

      UAdsorbateBondAverage[k][i]=0.0;
      UAdsorbateUreyBradleyAverage[k][i]=0.0;
      UAdsorbateBondBondAverage[k][i]=0.0;
      UAdsorbateBendAverage[k][i]=0.0;
      UAdsorbateInversionBendAverage[k][i]=0.0;
      UAdsorbateTorsionAverage[k][i]=0.0;
      UAdsorbateImproperTorsionAverage[k][i]=0.0;
      UAdsorbateBondBondAverage[k][i]=0.0;
      UAdsorbateBendBendAverage[k][i]=0.0;
      UAdsorbateBondBendAverage[k][i]=0.0;
      UAdsorbateBondTorsionAverage[k][i]=0.0;
      UAdsorbateBendTorsionAverage[k][i]=0.0;
      UAdsorbateIntraVDWAverage[k][i]=0.0;
      UAdsorbateIntraChargeChargeAverage[k][i]=0.0;
      UAdsorbateIntraChargeBondDipoleAverage[k][i]=0.0;
      UAdsorbateIntraBondDipoleBondDipoleAverage[k][i]=0.0;

      UHostHostAverage[k][i]=0.0;
      UAdsorbateAdsorbateAverage[k][i]=0.0;
      UCationCationAverage[k][i]=0.0;
      UHostAdsorbateAverage[k][i]=0.0;
      UHostCationAverage[k][i]=0.0;
      UAdsorbateCationAverage[k][i]=0.0;

      UHostHostVDWAverage[k][i]=0.0;
      UAdsorbateAdsorbateVDWAverage[k][i]=0.0;
      UCationCationVDWAverage[k][i]=0.0;
      UHostAdsorbateVDWAverage[k][i]=0.0;
      UHostCationVDWAverage[k][i]=0.0;
      UAdsorbateCationVDWAverage[k][i]=0.0;

      UHostHostCoulombAverage[k][i]=0.0;
      UAdsorbateAdsorbateCoulombAverage[k][i]=0.0;
      UCationCationCoulombAverage[k][i]=0.0;
      UHostAdsorbateCoulombAverage[k][i]=0.0;
      UHostCationCoulombAverage[k][i]=0.0;
      UAdsorbateCationCoulombAverage[k][i]=0.0;

      TotalSystemDipoleAverage[k][i].x=0.0;
      TotalSystemDipoleAverage[k][i].y=0.0;
      TotalSystemDipoleAverage[k][i].z=0.0;
      TotalSystemDipoleSquaredAverage[k][i].x=0.0;
      TotalSystemDipoleSquaredAverage[k][i].y=0.0;
      TotalSystemDipoleSquaredAverage[k][i].z=0.0;
      TotalSystemNormDipoleAverage[k][i]=0.0;
      TotalSystemNormDipoleSquaredAverage[k][i]=0.0;

      UTailCorrectionAverage[k][i]=0.0;
      UDistanceConstraintsAverage[k][i]=0.0;
      UAngleConstraintsAverage[k][i]=0.0;
      UDihedralConstraintsAverage[k][i]=0.0;
      UInversionBendConstraintsAverage[k][i]=0.0;
      UOutOfPlaneDistanceConstraintsAverage[k][i]=0.0;
      UExclusionConstraintsAverage[k][i]=0.0;

      UHostPolarizationAverage[k][i]=0.0;
      UAdsorbatePolarizationAverage[k][i]=0.0;
      UCationPolarizationAverage[k][i]=0.0;
      UHostBackPolarizationAverage[k][i]=0.0;
      UAdsorbateBackPolarizationAverage[k][i]=0.0;
      UCationBackPolarizationAverage[k][i]=0.0;

      UTotalAverage[k][i]=0.0;
      NumberOfMoleculesAverage[k][i]=0.0;
      DensityAverage[k][i]=0.0;
      for(j=0;j<NumberOfComponents;j++)
      {
        NumberOfMoleculesPerComponentAverage[k][j][i]=0.0;
        NumberOfExcessMoleculesPerComponentAverage[k][j][i]=0.0;
        DensityPerComponentAverage[k][j][i]=0.0;

        WidomRosenbluthFactorAverage[k][j][i]=0.0;
        WidomRosenbluthFactorCount[k][j][i]=0.0;

        WidomEnergyDifferenceAverage[k][j][i]=0.0;

        WidomEnergyFrameworkAverage[k][j][i]=0.0;
        WidomEnergyFrameworkCount[k][j][i]=0.0;
      }
      NumberOfMoleculesSquaredAverage[k][i]=0.0;
      TotalEnergyTimesNumberOfMoleculesAverage[k][i]=0.0;
      HostAdsorbateEnergyTimesNumberOfMoleculesAverage[k][i]=0.0;
      AdsorbateAdsorbateEnergyTimesNumberOfMoleculesAverage[k][i]=0.0;

      TemperatureAverage[k][i]=0.0;
      TemperatureCellAverage[k][i]=0.0;
      MolecularPressureAverage[k][i]=0.0;
      CompressibilityAverage[k][i]=0.0;

      MolecularStressTensorAverage[k][i].ax=0.0;
      MolecularStressTensorAverage[k][i].ay=0.0;
      MolecularStressTensorAverage[k][i].az=0.0;
      MolecularStressTensorAverage[k][i].bx=0.0;
      MolecularStressTensorAverage[k][i].by=0.0;
      MolecularStressTensorAverage[k][i].bz=0.0;
      MolecularStressTensorAverage[k][i].cx=0.0;
      MolecularStressTensorAverage[k][i].cy=0.0;
      MolecularStressTensorAverage[k][i].cz=0.0;
      PressureIdealGasPartAverage[k][i]=0.0;
      PressureExcessPartAverage[k][i]=0.0;
      PressureTailCorrectionAverage[k][i]=0.0;
      PressureAverage[k][i]=0.0;

      BoxAverage[k][i].x=0.0;
      BoxAverage[k][i].y=0.0;
      BoxAverage[k][i].z=0.0;

      BoxAverageAX[k][i]=0.0;
      BoxAverageAY[k][i]=0.0;
      BoxAverageAZ[k][i]=0.0;
      BoxAverageBX[k][i]=0.0;
      BoxAverageBY[k][i]=0.0;
      BoxAverageBZ[k][i]=0.0;
      BoxAverageCX[k][i]=0.0;
      BoxAverageCY[k][i]=0.0;
      BoxAverageCZ[k][i]=0.0;

      BoxLengthAverage[k][i].x=0.0;
      BoxLengthAverage[k][i].y=0.0;
      BoxLengthAverage[k][i].z=0.0;
      AlphaAngleAverage[k][i]=0.0;
      BetaAngleAverage[k][i]=0.0;
      GammaAngleAverage[k][i]=0.0;

      VolumeAverage[k][i]=0.0;
      VolumeSquaredAverage[k][i]=0.0;

      TotalEnergyAverage[k][i]=0.0;
      TotalEnergySquaredAverage[k][i]=0.0;

      EnthalpyAverage[k][i]=0.0;
      EnthalpySquaredAverage[k][i]=0.0;

      EnthalpyTimesVolumeAverage[k][i]=0.0;
      EnthalpyTimesEnergyAverage[k][i]=0.0;

      UNoseHooverAverage[k][i]=0.0;

      HeatOfVaporization[k][i]=0.0;
      EnergyPerMolecule[k][i]=0.0;
      VolumePerMolecule[k][i]=0.0;

      InitializeMatrix9x9(&BornTermAverage[k][i]);
    }
  }

}

void UpdateEnergyAveragesCurrentSystem(void)
{
  int i,nr;
  VECTOR dipole_adsorbates,dipole_cations;
  REAL Enthalpy,density;
  REAL NumberOfMolecules;
  REAL PressureIdealGas;
  REAL PressureTail;

  // check for new block
  if(CurrentCycle==BlockCycle[Block])
    Block++;

  BlockCount[CurrentSystem][Block]+=1.0;

  if(ComputeBornTerm)
    AddBornTermToAverages();

  ConfigurationalStressTensorAverage[CurrentSystem][Block].ax+=ConfigurationalStressTensor[CurrentSystem].ax;
  ConfigurationalStressTensorAverage[CurrentSystem][Block].bx+=ConfigurationalStressTensor[CurrentSystem].bx;
  ConfigurationalStressTensorAverage[CurrentSystem][Block].cx+=ConfigurationalStressTensor[CurrentSystem].cx;

  ConfigurationalStressTensorAverage[CurrentSystem][Block].ay+=ConfigurationalStressTensor[CurrentSystem].ay;
  ConfigurationalStressTensorAverage[CurrentSystem][Block].by+=ConfigurationalStressTensor[CurrentSystem].by;
  ConfigurationalStressTensorAverage[CurrentSystem][Block].cy+=ConfigurationalStressTensor[CurrentSystem].cy;

  ConfigurationalStressTensorAverage[CurrentSystem][Block].az+=ConfigurationalStressTensor[CurrentSystem].az;
  ConfigurationalStressTensorAverage[CurrentSystem][Block].bz+=ConfigurationalStressTensor[CurrentSystem].bz;
  ConfigurationalStressTensorAverage[CurrentSystem][Block].cz+=ConfigurationalStressTensor[CurrentSystem].cz;

  StressTensorAverage[CurrentSystem][Block].ax+=StressTensor[CurrentSystem].ax;
  StressTensorAverage[CurrentSystem][Block].bx+=StressTensor[CurrentSystem].bx;
  StressTensorAverage[CurrentSystem][Block].cx+=StressTensor[CurrentSystem].cx;

  StressTensorAverage[CurrentSystem][Block].ay+=StressTensor[CurrentSystem].ay;
  StressTensorAverage[CurrentSystem][Block].by+=StressTensor[CurrentSystem].by;
  StressTensorAverage[CurrentSystem][Block].cy+=StressTensor[CurrentSystem].cy;

  StressTensorAverage[CurrentSystem][Block].az+=StressTensor[CurrentSystem].az;
  StressTensorAverage[CurrentSystem][Block].bz+=StressTensor[CurrentSystem].bz;
  StressTensorAverage[CurrentSystem][Block].cz+=StressTensor[CurrentSystem].cz;



  UHostBondAverage[CurrentSystem][Block]+=UHostBond[CurrentSystem];
  UHostUreyBradleyAverage[CurrentSystem][Block]+=UHostUreyBradley[CurrentSystem];
  UHostBendAverage[CurrentSystem][Block]+=UHostBend[CurrentSystem];
  UHostInversionBendAverage[CurrentSystem][Block]+=UHostInversionBend[CurrentSystem];
  UHostTorsionAverage[CurrentSystem][Block]+=UHostTorsion[CurrentSystem];
  UHostImproperTorsionAverage[CurrentSystem][Block]+=UHostImproperTorsion[CurrentSystem];
  UHostBondBondAverage[CurrentSystem][Block]+=UHostBondBond[CurrentSystem];
  UHostBendBendAverage[CurrentSystem][Block]+=UHostBendBend[CurrentSystem];
  UHostBondBendAverage[CurrentSystem][Block]+=UHostBondBend[CurrentSystem];
  UHostBondTorsionAverage[CurrentSystem][Block]+=UHostBondTorsion[CurrentSystem];
  UHostBendTorsionAverage[CurrentSystem][Block]+=UHostBendTorsion[CurrentSystem];

  UCationBondAverage[CurrentSystem][Block]+=UCationBond[CurrentSystem];
  UCationUreyBradleyAverage[CurrentSystem][Block]+=UCationUreyBradley[CurrentSystem];
  UCationBendAverage[CurrentSystem][Block]+=UCationBend[CurrentSystem];
  UCationInversionBendAverage[CurrentSystem][Block]+=UCationInversionBend[CurrentSystem];
  UCationTorsionAverage[CurrentSystem][Block]+=UCationTorsion[CurrentSystem];
  UCationImproperTorsionAverage[CurrentSystem][Block]+=UCationImproperTorsion[CurrentSystem];
  UCationBondBondAverage[CurrentSystem][Block]+=UCationBondBond[CurrentSystem];
  UCationBendBendAverage[CurrentSystem][Block]+=UCationBendBend[CurrentSystem];
  UCationBondBendAverage[CurrentSystem][Block]+=UCationBondBend[CurrentSystem];
  UCationBondTorsionAverage[CurrentSystem][Block]+=UCationBondTorsion[CurrentSystem];
  UCationBendTorsionAverage[CurrentSystem][Block]+=UCationBendTorsion[CurrentSystem];
  UCationIntraVDWAverage[CurrentSystem][Block]+=UCationIntraVDW[CurrentSystem];
  UCationIntraChargeChargeAverage[CurrentSystem][Block]+=UCationIntraChargeCharge[CurrentSystem];
  UCationIntraChargeBondDipoleAverage[CurrentSystem][Block]+=UCationIntraChargeBondDipole[CurrentSystem];
  UCationIntraBondDipoleBondDipoleAverage[CurrentSystem][Block]+=UCationIntraBondDipoleBondDipole[CurrentSystem];

  UAdsorbateBondAverage[CurrentSystem][Block]+=UAdsorbateBond[CurrentSystem];
  UAdsorbateUreyBradleyAverage[CurrentSystem][Block]+=UAdsorbateUreyBradley[CurrentSystem];
  UAdsorbateBendAverage[CurrentSystem][Block]+=UAdsorbateBend[CurrentSystem];
  UAdsorbateInversionBendAverage[CurrentSystem][Block]+=UAdsorbateInversionBend[CurrentSystem];
  UAdsorbateTorsionAverage[CurrentSystem][Block]+=UAdsorbateTorsion[CurrentSystem];
  UAdsorbateImproperTorsionAverage[CurrentSystem][Block]+=UAdsorbateImproperTorsion[CurrentSystem];
  UAdsorbateBondBondAverage[CurrentSystem][Block]+=UAdsorbateBondBond[CurrentSystem];
  UAdsorbateBendBendAverage[CurrentSystem][Block]+=UAdsorbateBendBend[CurrentSystem];
  UAdsorbateBondBendAverage[CurrentSystem][Block]+=UAdsorbateBondBend[CurrentSystem];
  UAdsorbateBondTorsionAverage[CurrentSystem][Block]+=UAdsorbateBondTorsion[CurrentSystem];
  UAdsorbateBendTorsionAverage[CurrentSystem][Block]+=UAdsorbateBendTorsion[CurrentSystem];
  UAdsorbateIntraVDWAverage[CurrentSystem][Block]+=UAdsorbateIntraVDW[CurrentSystem];
  UAdsorbateIntraChargeChargeAverage[CurrentSystem][Block]+=UAdsorbateIntraChargeCharge[CurrentSystem];
  UAdsorbateIntraChargeBondDipoleAverage[CurrentSystem][Block]+=UAdsorbateIntraChargeBondDipole[CurrentSystem];
  UAdsorbateIntraBondDipoleBondDipoleAverage[CurrentSystem][Block]+=UAdsorbateIntraBondDipoleBondDipole[CurrentSystem];

  UHostHostAverage[CurrentSystem][Block]+=UHostHost[CurrentSystem];
  UAdsorbateAdsorbateAverage[CurrentSystem][Block]+=UAdsorbateAdsorbate[CurrentSystem];
  UCationCationAverage[CurrentSystem][Block]+=UCationCation[CurrentSystem];
  UHostAdsorbateAverage[CurrentSystem][Block]+=UHostAdsorbate[CurrentSystem];
  UHostCationAverage[CurrentSystem][Block]+=UHostCation[CurrentSystem];
  UAdsorbateCationAverage[CurrentSystem][Block]+=UAdsorbateCation[CurrentSystem];

  UHostHostVDWAverage[CurrentSystem][Block]+=UHostHostVDW[CurrentSystem];
  UAdsorbateAdsorbateVDWAverage[CurrentSystem][Block]+=UAdsorbateAdsorbateVDW[CurrentSystem];
  UCationCationVDWAverage[CurrentSystem][Block]+=UCationCationVDW[CurrentSystem];
  UHostAdsorbateVDWAverage[CurrentSystem][Block]+=UHostAdsorbateVDW[CurrentSystem];
  UHostCationVDWAverage[CurrentSystem][Block]+=UHostCationVDW[CurrentSystem];
  UAdsorbateCationVDWAverage[CurrentSystem][Block]+=UAdsorbateCationVDW[CurrentSystem];

  UHostHostCoulombAverage[CurrentSystem][Block]+=UHostHostCoulomb[CurrentSystem];
  UAdsorbateAdsorbateCoulombAverage[CurrentSystem][Block]+=UAdsorbateAdsorbateCoulomb[CurrentSystem];
  UCationCationCoulombAverage[CurrentSystem][Block]+=UCationCationCoulomb[CurrentSystem];
  UHostAdsorbateCoulombAverage[CurrentSystem][Block]+=UHostAdsorbateCoulomb[CurrentSystem];
  UHostCationCoulombAverage[CurrentSystem][Block]+=UHostCationCoulomb[CurrentSystem];
  UAdsorbateCationCoulombAverage[CurrentSystem][Block]+=UAdsorbateCationCoulomb[CurrentSystem];

  dipole_adsorbates=ComputeTotalDipoleMomentSystemAdsorbates();
  dipole_cations=ComputeTotalDipoleMomentSystemCations();
  dipole_cations.x=dipole_cations.y=dipole_cations.z=0.0;

  TotalSystemDipoleAverage[CurrentSystem][Block].x+=dipole_adsorbates.x+dipole_cations.x;
  TotalSystemDipoleAverage[CurrentSystem][Block].y+=dipole_adsorbates.y+dipole_cations.y;
  TotalSystemDipoleAverage[CurrentSystem][Block].z+=dipole_adsorbates.z+dipole_cations.z;
  TotalSystemNormDipoleAverage[CurrentSystem][Block]+=sqrt(SQR(dipole_adsorbates.x+dipole_cations.x)+
        SQR(dipole_adsorbates.y+dipole_cations.y)+SQR(dipole_adsorbates.z+dipole_cations.z));

  TotalSystemDipoleSquaredAverage[CurrentSystem][Block].x+=SQR(dipole_adsorbates.x+dipole_cations.x);
  TotalSystemDipoleSquaredAverage[CurrentSystem][Block].y+=SQR(dipole_adsorbates.y+dipole_cations.y);
  TotalSystemDipoleSquaredAverage[CurrentSystem][Block].z+=SQR(dipole_adsorbates.z+dipole_cations.z);
  TotalSystemNormDipoleSquaredAverage[CurrentSystem][Block]+=SQR(dipole_adsorbates.x+dipole_cations.x)+
        SQR(dipole_adsorbates.y+dipole_cations.y)+SQR(dipole_adsorbates.z+dipole_cations.z);

  UHostPolarizationAverage[CurrentSystem][Block]+=UHostPolarization[CurrentSystem];
  UAdsorbatePolarizationAverage[CurrentSystem][Block]+=UAdsorbatePolarization[CurrentSystem];
  UCationPolarizationAverage[CurrentSystem][Block]+=UCationPolarization[CurrentSystem];
  UHostBackPolarizationAverage[CurrentSystem][Block]+=UHostBackPolarization[CurrentSystem];
  UAdsorbateBackPolarizationAverage[CurrentSystem][Block]+=UAdsorbateBackPolarization[CurrentSystem];
  UCationBackPolarizationAverage[CurrentSystem][Block]+=UCationBackPolarization[CurrentSystem];

  UTailCorrectionAverage[CurrentSystem][Block]+=UTailCorrection[CurrentSystem];
  UDistanceConstraintsAverage[CurrentSystem][Block]+=UDistanceConstraints[CurrentSystem];
  UAngleConstraintsAverage[CurrentSystem][Block]+=UAngleConstraints[CurrentSystem];
  UDihedralConstraintsAverage[CurrentSystem][Block]+=UDihedralConstraints[CurrentSystem];
  UInversionBendConstraintsAverage[CurrentSystem][Block]+=UInversionBendConstraints[CurrentSystem];
  UOutOfPlaneDistanceConstraintsAverage[CurrentSystem][Block]+=UOutOfPlaneDistanceConstraints[CurrentSystem];
  UExclusionConstraintsAverage[CurrentSystem][Block]+=UExclusionConstraints[CurrentSystem];
  UTotalAverage[CurrentSystem][Block]+=UTotal[CurrentSystem];

  NumberOfMoleculesAverage[CurrentSystem][Block]+=NumberOfAdsorbateMolecules[CurrentSystem]
                                                 -NumberOfFractionalAdsorbateMolecules[CurrentSystem];
  NumberOfMoleculesSquaredAverage[CurrentSystem][Block]+=SQR(NumberOfAdsorbateMolecules[CurrentSystem]
                                                        -NumberOfFractionalAdsorbateMolecules[CurrentSystem]);
  TotalEnergyTimesNumberOfMoleculesAverage[CurrentSystem][Block]+=UTotal[CurrentSystem]*(REAL)NumberOfAdsorbateMolecules[CurrentSystem];
  HostAdsorbateEnergyTimesNumberOfMoleculesAverage[CurrentSystem][Block]+=UHostAdsorbate[CurrentSystem]*(REAL)NumberOfAdsorbateMolecules[CurrentSystem];
  AdsorbateAdsorbateEnergyTimesNumberOfMoleculesAverage[CurrentSystem][Block]+=UAdsorbateAdsorbate[CurrentSystem]*(REAL)NumberOfAdsorbateMolecules[CurrentSystem];

  nr=NumberOfUnitCells[0].x*NumberOfUnitCells[0].y*NumberOfUnitCells[0].z;
  for(i=0;i<NumberOfComponents;i++)
  {
    NumberOfMoleculesPerComponentAverage[CurrentSystem][i][Block]+=Components[i].NumberOfMolecules[CurrentSystem]
                         -(Components[i].CFMoleculePresent[CurrentSystem]?1:0)
                         -Components[i].NumberOfRXMCMoleculesPresent[CurrentSystem];
    NumberOfExcessMoleculesPerComponentAverage[CurrentSystem][i][Block]+=
             (REAL)Components[i].NumberOfMolecules[CurrentSystem]
              -Components[i].AmountOfExcessMolecules[CurrentSystem]
              -(Components[i].CFMoleculePresent[CurrentSystem]?1:0)
              -Components[i].NumberOfRXMCMoleculesPresent[CurrentSystem];
    density=Components[i].Mass*(REAL)(Components[i].NumberOfMolecules[CurrentSystem]
                    -(Components[i].CFMoleculePresent[CurrentSystem]?1:0)
                    -Components[i].NumberOfRXMCMoleculesPresent[CurrentSystem])
                    /Volume[CurrentSystem];
    DensityPerComponentAverage[CurrentSystem][i][Block]+=density;
    DensityAverage[CurrentSystem][Block]+=density;
  }

  TemperatureAverage[CurrentSystem][Block]+=2.0*UKinetic[CurrentSystem]/(K_B*DegreesOfFreedom[CurrentSystem]);
  TemperatureCellAverage[CurrentSystem][Block]+=GetCellTemperature();
  TemperatureTranslationAverage[CurrentSystem][Block]+=2.0*(UAdsorbateTranslationalKinetic[CurrentSystem]+
         UCationTranslationalKinetic[CurrentSystem]+UHostKinetic[CurrentSystem])/(K_B*DegreesOfFreedomTranslation[CurrentSystem]);
  TemperatureRotationAverage[CurrentSystem][Block]+=2.0*(UAdsorbateRotationalKinetic[CurrentSystem]+
         UCationRotationalKinetic[CurrentSystem])/(K_B*DegreesOfFreedomRotation[CurrentSystem]);
  TemperatureRotationAdsorbateAverage[CurrentSystem][Block]+=2.0*(UAdsorbateRotationalKinetic[CurrentSystem])/
         (K_B*DegreesOfFreedomRotationalAdsorbates[CurrentSystem]);
  TemperatureTranslationAdsorbateAverage[CurrentSystem][Block]+=2.0*(UAdsorbateTranslationalKinetic[CurrentSystem])/
         (K_B*DegreesOfFreedomTranslationalAdsorbates[CurrentSystem]);

  TemperatureAdsorbatesAverage[CurrentSystem][Block]+=2.0*(UAdsorbateTranslationalKinetic[CurrentSystem]+
              UAdsorbateRotationalKinetic[CurrentSystem])/(K_B*DegreesOfFreedomAdsorbates[CurrentSystem]);
  TemperatureCationsAverage[CurrentSystem][Block]+=2.0*UCationKinetic[CurrentSystem]/
                                             (K_B*DegreesOfFreedomCations[CurrentSystem]);
  TemperatureFrameworkAverage[CurrentSystem][Block]+=2.0*UHostKinetic[CurrentSystem]/
                                             (K_B*DegreesOfFreedomFramework[CurrentSystem]);

  NumberOfMolecules=NumberOfAdsorbateMolecules[CurrentSystem]+NumberOfCationMolecules[CurrentSystem];

  if(ComputeMolecularPressure[CurrentSystem])
  {
    ComputeMolecularPressureTensor(&MolecularStressTensor[CurrentSystem],&PressureIdealGas,&PressureTail);

    PressureIdealGasPartAverage[CurrentSystem][Block]+=PressureIdealGas;
    PressureExcessPartAverage[CurrentSystem][Block]-=(MolecularStressTensor[CurrentSystem].ax+MolecularStressTensor[CurrentSystem].by+
                MolecularStressTensor[CurrentSystem].cz)/(3.0*Volume[CurrentSystem]);
    PressureTailCorrectionAverage[CurrentSystem][Block]+=PressureTail;

    MolecularStressTensor[CurrentSystem].ax=PressureIdealGas+PressureTail-MolecularStressTensor[CurrentSystem].ax/Volume[CurrentSystem];
    MolecularStressTensor[CurrentSystem].ay=-MolecularStressTensor[CurrentSystem].ay/Volume[CurrentSystem];
    MolecularStressTensor[CurrentSystem].az=-MolecularStressTensor[CurrentSystem].az/Volume[CurrentSystem];

    MolecularStressTensor[CurrentSystem].bx=-MolecularStressTensor[CurrentSystem].bx/Volume[CurrentSystem];
    MolecularStressTensor[CurrentSystem].by=PressureIdealGas+PressureTail-MolecularStressTensor[CurrentSystem].by/Volume[CurrentSystem];
    MolecularStressTensor[CurrentSystem].bz=-MolecularStressTensor[CurrentSystem].bz/Volume[CurrentSystem];

    MolecularStressTensor[CurrentSystem].cx=-MolecularStressTensor[CurrentSystem].cx/Volume[CurrentSystem];
    MolecularStressTensor[CurrentSystem].cy=-MolecularStressTensor[CurrentSystem].cy/Volume[CurrentSystem];
    MolecularStressTensor[CurrentSystem].cz=PressureIdealGas+PressureTail-MolecularStressTensor[CurrentSystem].cz/Volume[CurrentSystem];

    MolecularStressTensorAverage[CurrentSystem][Block].ax+=MolecularStressTensor[CurrentSystem].ax;
    MolecularStressTensorAverage[CurrentSystem][Block].ay+=MolecularStressTensor[CurrentSystem].ay;
    MolecularStressTensorAverage[CurrentSystem][Block].az+=MolecularStressTensor[CurrentSystem].az;

    MolecularStressTensorAverage[CurrentSystem][Block].bx+=MolecularStressTensor[CurrentSystem].bx;
    MolecularStressTensorAverage[CurrentSystem][Block].by+=MolecularStressTensor[CurrentSystem].by;
    MolecularStressTensorAverage[CurrentSystem][Block].bz+=MolecularStressTensor[CurrentSystem].bz;

    MolecularStressTensorAverage[CurrentSystem][Block].cx+=MolecularStressTensor[CurrentSystem].cx;
    MolecularStressTensorAverage[CurrentSystem][Block].cy+=MolecularStressTensor[CurrentSystem].cy;
    MolecularStressTensorAverage[CurrentSystem][Block].cz+=MolecularStressTensor[CurrentSystem].cz;
  }

  UNoseHooverAverage[CurrentSystem][Block]+=UNoseHoover[CurrentSystem];

  BoxAverage[CurrentSystem][Block].x+=Box[CurrentSystem].ax;
  BoxAverage[CurrentSystem][Block].y+=Box[CurrentSystem].by;
  BoxAverage[CurrentSystem][Block].z+=Box[CurrentSystem].cz;

  BoxAverageAX[CurrentSystem][Block]+=Box[CurrentSystem].ax;
  BoxAverageAY[CurrentSystem][Block]+=Box[CurrentSystem].ay;
  BoxAverageAZ[CurrentSystem][Block]+=Box[CurrentSystem].az;
  BoxAverageBX[CurrentSystem][Block]+=Box[CurrentSystem].bx;
  BoxAverageBY[CurrentSystem][Block]+=Box[CurrentSystem].by;
  BoxAverageBZ[CurrentSystem][Block]+=Box[CurrentSystem].bz;
  BoxAverageCX[CurrentSystem][Block]+=Box[CurrentSystem].cx;
  BoxAverageCY[CurrentSystem][Block]+=Box[CurrentSystem].cy;
  BoxAverageCZ[CurrentSystem][Block]+=Box[CurrentSystem].cz;

  BoxLengthAverage[CurrentSystem][Block].x+=BoxProperties[CurrentSystem].ax;
  BoxLengthAverage[CurrentSystem][Block].y+=BoxProperties[CurrentSystem].ay;
  BoxLengthAverage[CurrentSystem][Block].z+=BoxProperties[CurrentSystem].az;

  AlphaAngleAverage[CurrentSystem][Block]+=AlphaAngle[CurrentSystem]*RAD2DEG;
  BetaAngleAverage[CurrentSystem][Block]+=BetaAngle[CurrentSystem]*RAD2DEG;
  GammaAngleAverage[CurrentSystem][Block]+=GammaAngle[CurrentSystem]*RAD2DEG;

  VolumeAverage[CurrentSystem][Block]+=Volume[CurrentSystem];

  VolumeSquaredAverage[CurrentSystem][Block]+=SQR(Volume[CurrentSystem]);

  TotalEnergyAverage[CurrentSystem][Block]+=UTotal[CurrentSystem];
  TotalEnergySquaredAverage[CurrentSystem][Block]+=SQR(UTotal[CurrentSystem]);

  Enthalpy=UTotal[CurrentSystem]+Volume[CurrentSystem]*therm_baro_stats.ExternalPressure[CurrentSystem][0];
  EnthalpyAverage[CurrentSystem][Block]+=Enthalpy;
  EnthalpySquaredAverage[CurrentSystem][Block]+=SQR(Enthalpy);

  EnthalpyTimesVolumeAverage[CurrentSystem][Block]+=Enthalpy*Volume[CurrentSystem];
  EnthalpyTimesEnergyAverage[CurrentSystem][Block]+=Enthalpy*UTotal[CurrentSystem];

  HeatOfVaporization[CurrentSystem][Block]+=therm_baro_stats.ExternalTemperature[CurrentSystem]-
                              (UAdsorbateAdsorbate[CurrentSystem]+UCationCation[CurrentSystem])/
                              (NumberOfAdsorbateMolecules[CurrentSystem]+NumberOfCationMolecules[CurrentSystem]);
  if(NumberOfMolecules>0)
  {
    EnergyPerMolecule[CurrentSystem][Block]+=UTotal[CurrentSystem]/NumberOfMolecules;
    VolumePerMolecule[CurrentSystem][Block]+=Volume[CurrentSystem]/NumberOfMolecules;
  }

  if(ComputePrincipleMomentsOfInertia)
    MeasurePrincipleMomentsOfInertia();

  CompressibilityAverage[CurrentSystem][Block]+=((MolecularStressTensor[CurrentSystem].ax+MolecularStressTensor[CurrentSystem].by+MolecularStressTensor[CurrentSystem].cz)/3.0)*
           Volume[CurrentSystem]*Beta[CurrentSystem]/NumberOfMolecules;

  //UpdateCrystallographics();
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

VECTOR GetAverageVectorProperty(VECTOR **Property)
{
  int i;
  VECTOR sum1;
  REAL sum2;

  sum1.x=sum1.y=sum1.z=sum2=0.0;
  for(i=0;i<NR_BLOCKS;i++)
    if(BlockCount[CurrentSystem][i]>0.0)
    {
      sum1.x+=Property[CurrentSystem][i].x;
      sum1.y+=Property[CurrentSystem][i].y;
      sum1.z+=Property[CurrentSystem][i].z;
      sum2+=BlockCount[CurrentSystem][i];
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
    if(BlockCount[CurrentSystem][i]>0.0)
    {
      AddRealMatrix9x9(&sum1,sum1,Property[CurrentSystem][i]);
      sum2+=BlockCount[CurrentSystem][i];
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

VECTOR GetAverageComponentVectorProperty(VECTOR ***Property,int comp)
{
  int i;
  VECTOR sum1;
  REAL sum2;

  sum1.x=sum1.y=sum1.z=0.0;
  sum2=0.0;
  for(i=0;i<NR_BLOCKS;i++)
    if(BlockCount[CurrentSystem][i]>0.0)
    {
      sum1.x+=Property[CurrentSystem][comp][i].x;
      sum1.y+=Property[CurrentSystem][comp][i].y;
      sum1.z+=Property[CurrentSystem][comp][i].z;
      sum2+=BlockCount[CurrentSystem][i];
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
        Components[i].NumberOfMolecules[CurrentSystem]-(Components[i].CFMoleculePresent[CurrentSystem]?1:0)-Components[i].NumberOfRXMCMoleculesPresent[CurrentSystem],
        Components[i].CFMoleculePresent[CurrentSystem]?1:0,
        Components[i].NumberOfRXMCMoleculesPresent[CurrentSystem],
        (double)((Components[i].Mass*(REAL)(Components[i].NumberOfMolecules[CurrentSystem]
                 -(Components[i].CFMoleculePresent[CurrentSystem]?1:0)
                 -Components[i].NumberOfRXMCMoleculesPresent[CurrentSystem])/
                  Volume[CurrentSystem])*DENSITY_CONVERSION_FACTOR));

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
        Components[i].NumberOfMolecules[CurrentSystem]-(Components[i].CFMoleculePresent[CurrentSystem]?1:0)-Components[i].NumberOfRXMCMoleculesPresent[CurrentSystem],
        (Components[i].CFMoleculePresent[CurrentSystem]?1:0),
        Components[i].NumberOfRXMCMoleculesPresent[CurrentSystem],
        (double)(Components[i].Mass*(REAL)(Components[i].NumberOfMolecules[CurrentSystem]
                                           -(Components[i].CFMoleculePresent[CurrentSystem]?1:0)
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
        fprintf(FilePtr,"\tBiasing Factors: ");
        for(k=0;k<Components[i].CFLambdaHistogramSize;k++)
        {
          fprintf(FilePtr,"%4f ",Components[i].CFBiasingFactors[CurrentSystem][k]-shift);
          if((k+1)%10==0&&(k+1)!=Components[i].CFLambdaHistogramSize)  fprintf(FilePtr,"\n\t                 ");
        }
        fprintf(FilePtr,"\n");
      }

      loading=(REAL)(Components[i].NumberOfMolecules[CurrentSystem]
                    -(Components[i].CFMoleculePresent[CurrentSystem]?1:0)
                    -Components[i].NumberOfRXMCMoleculesPresent[CurrentSystem])/(REAL)number_of_unit_cells;
      average_loading=GetAverageComponentProperty(NumberOfMoleculesPerComponentAverage,i)/(REAL)number_of_unit_cells;
      fprintf(FilePtr,"\tabsolute adsorption: %9.5lf [mol/uc], %14.4lf [mol/kg],      %14.4lf [mg/g]\n",
         (double)loading,
         (double)(Components[i].MOLEC_PER_UC_TO_MOL_PER_KG[CurrentSystem]*loading),
         (double)(Components[i].MOLEC_PER_UC_TO_GRAM_PER_GRAM_OF_FRAMEWORK[CurrentSystem]*loading));
      fprintf(FilePtr,"\t                                         %14.4lf [cm^3 STP/g],  %14.4lf [cm^3 STP/cm^3]\n",
         (double)(Components[i].MOLEC_PER_UC_TO_CC_STP_G[CurrentSystem]*loading),
         (double)(Components[i].MOLEC_PER_UC_TO_CC_STP_CC[CurrentSystem]*loading));

      loading=((REAL)Components[i].NumberOfMolecules[CurrentSystem]
                    -Components[i].AmountOfExcessMolecules[CurrentSystem]
                    -(Components[i].CFMoleculePresent[CurrentSystem]?1:0)
                    -Components[i].NumberOfRXMCMoleculesPresent[CurrentSystem])/(REAL)number_of_unit_cells;
      average_loading=GetAverageComponentProperty(NumberOfExcessMoleculesPerComponentAverage,i)/(REAL)number_of_unit_cells;
      fprintf(FilePtr,"\texcess adsorption:   %9.5lf [mol/uc], %14.4lf [mol/kg],      %14.4lf [mg/g]\n",
         (double)loading,
         (double)(Components[i].MOLEC_PER_UC_TO_MOL_PER_KG[CurrentSystem]*loading),
         (double)(Components[i].MOLEC_PER_UC_TO_GRAM_PER_GRAM_OF_FRAMEWORK[CurrentSystem]*loading));
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
          NumberOfAdsorbateMolecules[CurrentSystem]-NumberOfFractionalAdsorbateMolecules[CurrentSystem]-NumberOfReactionAdsorbateMolecules[CurrentSystem],
          NumberOfFractionalAdsorbateMolecules[CurrentSystem],
          NumberOfReactionAdsorbateMolecules[CurrentSystem]);
  fprintf(FilePtr,"Number of Cations:         %6d (%d integer, %d fractional, % reaction\n",
      NumberOfCationMolecules[CurrentSystem],
      NumberOfCationMolecules[CurrentSystem]-NumberOfFractionalCationMolecules[CurrentSystem]-NumberOfReactionCationMolecules[CurrentSystem],
      NumberOfFractionalCationMolecules[CurrentSystem],
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
        Components[i].NumberOfMolecules[CurrentSystem]-(Components[i].CFMoleculePresent[CurrentSystem]?1:0)-Components[i].NumberOfRXMCMoleculesPresent[CurrentSystem],
        Components[i].CFMoleculePresent[CurrentSystem]?1:0,
        Components[i].NumberOfRXMCMoleculesPresent[CurrentSystem],
        (double)((Components[i].Mass*(REAL)(Components[i].NumberOfMolecules[CurrentSystem]
                       -(Components[i].CFMoleculePresent[CurrentSystem]?1:0)
                       -Components[i].NumberOfRXMCMoleculesPresent[CurrentSystem])/
                        Volume[CurrentSystem])*DENSITY_CONVERSION_FACTOR));
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
        Components[i].NumberOfMolecules[CurrentSystem]-(Components[i].CFMoleculePresent[CurrentSystem]?1:0)-Components[i].NumberOfRXMCMoleculesPresent[CurrentSystem],
        (Components[i].CFMoleculePresent[CurrentSystem]?1:0),
        Components[i].NumberOfRXMCMoleculesPresent[CurrentSystem],
        (double)(Components[i].Mass*(REAL)(Components[i].NumberOfMolecules[CurrentSystem]
                                         -(Components[i].CFMoleculePresent[CurrentSystem]?1:0)
                                         -Components[i].NumberOfRXMCMoleculesPresent[CurrentSystem])/
                                          Volume[CurrentSystem])*DENSITY_CONVERSION_FACTOR);
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
                   -(Components[i].CFMoleculePresent[CurrentSystem]?1:0)
                   -Components[i].NumberOfRXMCMoleculesPresent[CurrentSystem])/(REAL)number_of_unit_cells;
      average_loading=GetAverageComponentProperty(NumberOfMoleculesPerComponentAverage,i)/(REAL)number_of_unit_cells;
      fprintf(FilePtr,"\tabsolute adsorption: %9.5lf [mol/uc], %14.4lf [mol/kg],      %14.4lf [mg/g]\n",
         (double)loading,
         (double)(Components[i].MOLEC_PER_UC_TO_MOL_PER_KG[CurrentSystem]*loading),
         (double)(Components[i].MOLEC_PER_UC_TO_GRAM_PER_GRAM_OF_FRAMEWORK[CurrentSystem]*loading));
      fprintf(FilePtr,"\t                                          %14.4lf [cm^3 STP/g],  %14.4lf [cm^3 STP/cm^3]\n",
         (double)(Components[i].MOLEC_PER_UC_TO_CC_STP_G[CurrentSystem]*loading),
         (double)(Components[i].MOLEC_PER_UC_TO_CC_STP_CC[CurrentSystem]*loading));

      loading=((REAL)Components[i].NumberOfMolecules[CurrentSystem]
                    -Components[i].AmountOfExcessMolecules[CurrentSystem]
                    -(Components[i].CFMoleculePresent[CurrentSystem]?1:0)
                    -Components[i].NumberOfRXMCMoleculesPresent[CurrentSystem])/(REAL)number_of_unit_cells;
      average_loading=GetAverageComponentProperty(NumberOfExcessMoleculesPerComponentAverage,i)/(REAL)number_of_unit_cells;
      fprintf(FilePtr,"\texcess adsorption:   %9.5lf [mol/uc], %14.4lf [mol/kg],      %14.4lf [mg/g]\n",
         (double)loading,
         (double)(Components[i].MOLEC_PER_UC_TO_MOL_PER_KG[CurrentSystem]*loading),
         (double)(Components[i].MOLEC_PER_UC_TO_GRAM_PER_GRAM_OF_FRAMEWORK[CurrentSystem]*loading));
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
          NumberOfAdsorbateMolecules[CurrentSystem]-NumberOfFractionalAdsorbateMolecules[CurrentSystem]-NumberOfReactionAdsorbateMolecules[CurrentSystem],
          NumberOfFractionalAdsorbateMolecules[CurrentSystem],
          NumberOfReactionAdsorbateMolecules[CurrentSystem]);
  fprintf(FilePtr,"Number of Cations:         %6d (%d integer, %d fractional, %d reaction)\n",
          NumberOfCationMolecules[CurrentSystem],
          NumberOfCationMolecules[CurrentSystem]-NumberOfFractionalCationMolecules[CurrentSystem]-NumberOfReactionCationMolecules[CurrentSystem],
          NumberOfFractionalCationMolecules[CurrentSystem],
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
    if(BlockCount[CurrentSystem][i]>0.0)
    {
      sum1+=TemperatureAdsorbatesAverage[CurrentSystem][i];
      sum2+=BlockCount[CurrentSystem][i];
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
    if(BlockCount[CurrentSystem][i]>0.0)
    {
      sum1+=TemperatureCationsAverage[CurrentSystem][i];
      sum2+=BlockCount[CurrentSystem][i];
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
    if(BlockCount[CurrentSystem][i]>0.0)
    {
      sum1+=TemperatureFrameworkAverage[CurrentSystem][i];
      sum2+=BlockCount[CurrentSystem][i];
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
    if(BlockCount[CurrentSystem][i]>0.0)
    {
      sum1+=MolecularPressureAverage[CurrentSystem][i];
      sum2+=BlockCount[CurrentSystem][i];
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
    if(BlockCount[CurrentSystem][i]>0.0)
    {
      sum1+=CompressibilityAverage[CurrentSystem][i];
      sum2+=BlockCount[CurrentSystem][i];
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
    if(BlockCount[CurrentSystem][i]>0.0)
    {
      sum1+=MolecularPressureAverage[CurrentSystem][i];
      sum2+=BlockCount[CurrentSystem][i];
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
    if(BlockCount[CurrentSystem][i]>0.0)
    {
      Stress.ax+=MolecularStressTensorAverage[CurrentSystem][i].ax;
      Stress.ay+=MolecularStressTensorAverage[CurrentSystem][i].ay;
      Stress.az+=MolecularStressTensorAverage[CurrentSystem][i].az;
      Stress.bx+=MolecularStressTensorAverage[CurrentSystem][i].bx;
      Stress.by+=MolecularStressTensorAverage[CurrentSystem][i].by;
      Stress.bz+=MolecularStressTensorAverage[CurrentSystem][i].bz;
      Stress.cx+=MolecularStressTensorAverage[CurrentSystem][i].cx;
      Stress.cy+=MolecularStressTensorAverage[CurrentSystem][i].cy;
      Stress.cz+=MolecularStressTensorAverage[CurrentSystem][i].cz;
      sum+=BlockCount[CurrentSystem][i];
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
    if(BlockCount[CurrentSystem][i]>0.0)
    {
      sum+=PressureIdealGasPartAverage[CurrentSystem][i];
      count+=BlockCount[CurrentSystem][i];
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
    if(BlockCount[CurrentSystem][i]>0.0)
    {
      sum+=PressureExcessPartAverage[CurrentSystem][i];
      count+=BlockCount[CurrentSystem][i];
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
    if(BlockCount[CurrentSystem][i]>0.0)
    {
      sum+=PressureTailCorrectionAverage[CurrentSystem][i];
      count+=BlockCount[CurrentSystem][i];
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
    if(BlockCount[CurrentSystem][i]>0.0)
    {
      sum+=VolumeSquaredAverage[CurrentSystem][i];
      count+=BlockCount[CurrentSystem][i];
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
    if(BlockCount[CurrentSystem][i]>0.0)
    {
      HV+=EnthalpyTimesVolumeAverage[CurrentSystem][i];
      V+=VolumeAverage[CurrentSystem][i];
      H+=EnthalpyAverage[CurrentSystem][i];
      count+=BlockCount[CurrentSystem][i];
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
    if(BlockCount[CurrentSystem][i]>0.0)
    {
      V+=VolumeAverage[CurrentSystem][i];
      V2+=VolumeSquaredAverage[CurrentSystem][i];
      count+=BlockCount[CurrentSystem][i];
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
    if(BlockCount[CurrentSystem][i]>0.0)
    {
      H+=EnthalpyAverage[CurrentSystem][i];
      VH+=EnthalpyTimesVolumeAverage[CurrentSystem][i];
      V+=VolumeAverage[CurrentSystem][i];
      U+=TotalEnergyAverage[CurrentSystem][i];
      UH+=EnthalpyTimesEnergyAverage[CurrentSystem][i];
      count+=BlockCount[CurrentSystem][i];
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
    if(BlockCount[CurrentSystem][i]>0.0)
    {
      H2=EnthalpySquaredAverage[CurrentSystem][i]/BlockCount[CurrentSystem][i];
      H=EnthalpyAverage[CurrentSystem][i]/BlockCount[CurrentSystem][i];
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
    if(BlockCount[CurrentSystem][i]>0.0)
    {
      H2+=EnthalpySquaredAverage[CurrentSystem][i];
      H+=EnthalpyAverage[CurrentSystem][i];
      count+=BlockCount[CurrentSystem][i];
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
      sum+=WidomRosenbluthFactorAverage[CurrentSystem][comp][i]/
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
      sum+=WidomRosenbluthFactorAverage[CurrentSystem][comp][i];
      count+=WidomRosenbluthFactorCount[CurrentSystem][comp][i];
    }
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
    sum+=WidomEnergyDifferenceAverage[CurrentSystem][comp][i];
    count+=WidomRosenbluthFactorAverage[CurrentSystem][comp][i];
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
    sum+=WidomEnergyFrameworkAverage[CurrentSystem][comp][i];
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
      sum+=SurfaceAreaFrameworkAverage[CurrentSystem][i];
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
      sum+=SurfaceAreaFrameworksAverage[CurrentSystem][f][i];
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
      sum+=SurfaceAreaCationsAverage[CurrentSystem][i];
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
/*
  fprintf(FilePtr,"Isothermal expansion coefficient: %18.10lf [10^6 K^-1]\n",(double)GetAverageIsothermalExpansionCoefficient());
  fprintf(FilePtr,"Isothermal compressibility coefficient: %18.10lf [10^12 Pa^-1]\n",
          (double)GetAverageIsothermalCompressibilityCoefficient());
  fprintf(FilePtr,"Heat of vaporization: %18.10lf [J/mole/K]\n",(double)GetAverageProperty(HeatOfVaporization));
*/
  fprintf(FilePtr,"Compressibility: %18.10lf [-]\n",(double)GetAverageProperty(CompressibilityAverage));
  //fprintf(FilePtr,"Heat capacity Cp: %18.10lf [J/mole/K]\n",(double)GetAverageHeatCapacity());
  //fprintf(FilePtr,"Heat capacity Cp: %18.10lf [J/mole/K]\n",(double)GetAverageHeatCapacityConstantPressure());

  fprintf(FilePtr,"Henry coefficients\n");
  for(i=0;i<NumberOfComponents;i++)
    fprintf(FilePtr,"\tComponent %d: %lg [mol/kg/Pa] (Rosenbluth factor new: %lg [-])\n",i,
        (double)GetAverageHenryCoefficient(i),(double)GetAverageRosenbluthWeight(i));

  fprintf(FilePtr,"Energy <U_gh>_1-<U_h>_0 from Widom\n");
  for(i=0;i<NumberOfComponents;i++)
  {
    if(Components[i].Widom)
      fprintf(FilePtr,"\tComponent %d: %18.10lf [K]  (%18.10lf [kJ/mol])\n",i,
        (double)((GetWidomHeat(i)-GetFrameworkHeat(i))*ENERGY_TO_KELVIN),
        (double)((GetWidomHeat(i)-GetFrameworkHeat(i))*ENERGY_TO_KELVIN*KELVIN_TO_KJ_PER_MOL));
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
        if(BlockCount[CurrentSystem][i]>0.0)
        {
          InertiaAverage.x+=PrincipleMomentsOfInertiaAverage[CurrentSystem][k][i].x;
          InertiaAverage.y+=PrincipleMomentsOfInertiaAverage[CurrentSystem][k][i].y;
          InertiaAverage.z+=PrincipleMomentsOfInertiaAverage[CurrentSystem][k][i].z;
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

void PrintIntervalStatus(long long CurrentCycle,long long NumberOfCycles, FILE *FilePtr)
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
          (double)GetAverageProperty(BoxAverageAX),(double)GetAverageProperty(BoxAverageBX),(double)GetAverageProperty(BoxAverageCX));
  fprintf(FilePtr,"             %9.5lf %9.5lf %9.5lf [A]                %9.5lf %9.5lf %9.5lf [A]\n",
          (double)Box[CurrentSystem].ay,(double)Box[CurrentSystem].by,(double)Box[CurrentSystem].cy,
          (double)GetAverageProperty(BoxAverageAY),(double)GetAverageProperty(BoxAverageBY),(double)GetAverageProperty(BoxAverageCY));
  fprintf(FilePtr,"             %9.5lf %9.5lf %9.5lf [A]                %9.5lf %9.5lf %9.5lf [A]\n",
          (double)Box[CurrentSystem].az,(double)Box[CurrentSystem].bz,(double)Box[CurrentSystem].cz,
          (double)GetAverageProperty(BoxAverageAZ),(double)GetAverageProperty(BoxAverageBZ),(double)GetAverageProperty(BoxAverageCZ));
  fprintf(FilePtr,"Box-lengths:  %9.5lf %9.5lf %9.5lf [A] Average: %9.5lf %9.5lf %9.5lf [A]\n",
          (double)BoxProperties[CurrentSystem].ax,
          (double)BoxProperties[CurrentSystem].ay,
          (double)BoxProperties[CurrentSystem].az,
          (double)GetAverageVectorProperty(BoxLengthAverage).x,
          (double)GetAverageVectorProperty(BoxLengthAverage).y,
          (double)GetAverageVectorProperty(BoxLengthAverage).z);
  fprintf(FilePtr,"Box-angles:  %9.5lf %9.5lf %9.5lf [degrees] Average: %9.5lf %9.5lf %9.5lf [degrees]\n",
          (double)AlphaAngle[CurrentSystem]*RAD2DEG,
          (double)BetaAngle[CurrentSystem]*RAD2DEG,
          (double)GammaAngle[CurrentSystem]*RAD2DEG,
          (double)GetAverageProperty(AlphaAngleAverage),
          (double)GetAverageProperty(BetaAngleAverage),
          (double)GetAverageProperty(GammaAngleAverage));
  fprintf(FilePtr,"Volume: %9.5lf [A^3] Average Volume: %9.5lf [A^3]\n\n",
          (double)Volume[CurrentSystem],(double)GetAverageProperty(VolumeAverage));
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
      fprintf(FilePtr,"Component %d (%s), current number of integer/fractional/reaction molecules: %d/%d/%d (average %9.5lf), density: %9.5lf (average %9.5lf) kg/m^3]\n",
        i,
        Components[i].Name,
        Components[i].NumberOfMolecules[CurrentSystem]-(Components[i].CFMoleculePresent[CurrentSystem]?1:0)-Components[i].NumberOfRXMCMoleculesPresent[CurrentSystem],
        (Components[i].CFMoleculePresent[CurrentSystem]?1:0),
        Components[i].NumberOfRXMCMoleculesPresent[CurrentSystem],
        (double)(GetAverageComponentProperty(NumberOfMoleculesPerComponentAverage,i)/(REAL)number_of_unit_cells),
        (double)((Components[i].Mass*(REAL)(Components[i].NumberOfMolecules[CurrentSystem]-(Components[i].CFMoleculePresent[CurrentSystem]?1:0)-Components[i].NumberOfRXMCMoleculesPresent[CurrentSystem])/
                 Volume[CurrentSystem])*DENSITY_CONVERSION_FACTOR),
        (double)(GetAverageComponentProperty(DensityPerComponentAverage,i)*DENSITY_CONVERSION_FACTOR));
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
        Components[i].NumberOfMolecules[CurrentSystem]-(Components[i].CFMoleculePresent[CurrentSystem]?1:0)-Components[i].NumberOfRXMCMoleculesPresent[CurrentSystem],
        (Components[i].CFMoleculePresent[CurrentSystem]?1:0),
        Components[i].NumberOfRXMCMoleculesPresent[CurrentSystem],
        (double)GetAverageComponentProperty(NumberOfMoleculesPerComponentAverage,i),
        (double)(Components[i].Mass*(REAL)(Components[i].NumberOfMolecules[CurrentSystem]-(Components[i].CFMoleculePresent[CurrentSystem]?1:0)-Components[i].NumberOfRXMCMoleculesPresent[CurrentSystem])/
                Volume[CurrentSystem])*DENSITY_CONVERSION_FACTOR,
        (double)GetAverageComponentProperty(DensityPerComponentAverage,i)*DENSITY_CONVERSION_FACTOR);
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
                   -(Components[i].CFMoleculePresent[CurrentSystem]?1:0)
                   -Components[i].NumberOfRXMCMoleculesPresent[CurrentSystem])/(REAL)number_of_unit_cells;
      average_loading=GetAverageComponentProperty(NumberOfMoleculesPerComponentAverage,i)/(REAL)number_of_unit_cells;
      fprintf(FilePtr,"\tabsolute adsorption: %9.5lf (avg. %9.5lf) [mol/uc], %14.10lf (avg. %14.10lf) [mol/kg], %14.10lf (avg. %14.10lf) [mg/g]\n",
         (double)loading,
         (double)average_loading,
         (double)(Components[i].MOLEC_PER_UC_TO_MOL_PER_KG[CurrentSystem]*loading),
         (double)(Components[i].MOLEC_PER_UC_TO_MOL_PER_KG[CurrentSystem]*average_loading),
         (double)(Components[i].MOLEC_PER_UC_TO_GRAM_PER_GRAM_OF_FRAMEWORK[CurrentSystem]*loading),
         (double)(Components[i].MOLEC_PER_UC_TO_GRAM_PER_GRAM_OF_FRAMEWORK[CurrentSystem]*average_loading));
      fprintf(FilePtr,"\t                     %14.10lf (avg. %14.10lf) [cm^3 STP/g],  %14.10lf (avg. %14.10lf) [cm^3 STP/cm^3]\n",
         (double)(Components[i].MOLEC_PER_UC_TO_CC_STP_G[CurrentSystem]*loading),
         (double)(Components[i].MOLEC_PER_UC_TO_CC_STP_G[CurrentSystem]*average_loading),
         (double)(Components[i].MOLEC_PER_UC_TO_CC_STP_CC[CurrentSystem]*loading),
         (double)(Components[i].MOLEC_PER_UC_TO_CC_STP_CC[CurrentSystem]*average_loading));


      loading=((REAL)Components[i].NumberOfMolecules[CurrentSystem]
                    -Components[i].AmountOfExcessMolecules[CurrentSystem]
                    -(Components[i].CFMoleculePresent[CurrentSystem]?1:0)
                    -Components[i].NumberOfRXMCMoleculesPresent[CurrentSystem])/(REAL)number_of_unit_cells;
      average_loading=GetAverageComponentProperty(NumberOfExcessMoleculesPerComponentAverage,i)/(REAL)number_of_unit_cells;
      fprintf(FilePtr,"\texcess adsorption:   %14.10lf (avg. %14.10lf) [mol/uc], %14.10lf (avg. %14.10lf) [mol/kg], %14.10lf (avg. %14.10lf) [mg/g]\n",
         (double)loading,
         (double)average_loading,
         (double)(Components[i].MOLEC_PER_UC_TO_MOL_PER_KG[CurrentSystem]*loading),
         (double)(Components[i].MOLEC_PER_UC_TO_MOL_PER_KG[CurrentSystem]*average_loading),
         (double)(Components[i].MOLEC_PER_UC_TO_GRAM_PER_GRAM_OF_FRAMEWORK[CurrentSystem]*loading),
         (double)(Components[i].MOLEC_PER_UC_TO_GRAM_PER_GRAM_OF_FRAMEWORK[CurrentSystem]*average_loading));
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
          NumberOfAdsorbateMolecules[CurrentSystem]-NumberOfFractionalAdsorbateMolecules[CurrentSystem]-NumberOfReactionAdsorbateMolecules[CurrentSystem],
          NumberOfFractionalAdsorbateMolecules[CurrentSystem],
          NumberOfReactionAdsorbateMolecules[CurrentSystem]);
  fprintf(FilePtr,"Number of Cations:         %6d (%d integer, %d fractional, %d reaction)\n",
          NumberOfCationMolecules[CurrentSystem],
          NumberOfCationMolecules[CurrentSystem]-NumberOfFractionalCationMolecules[CurrentSystem]-NumberOfReactionCationMolecules[CurrentSystem],
          NumberOfFractionalCationMolecules[CurrentSystem],
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
        (double)GetAverageProperty(TemperatureAverage),
        (double)GetAverageProperty(TemperatureTranslationAverage),
        (double)GetAverageProperty(TemperatureRotationAverage));
    }
    if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
      fprintf(FilePtr,"Temperature Framework:  % 8.3lf (avg. % 8.3lf)\n",
        (double)(2.0*UHostKinetic[CurrentSystem]/(K_B*DegreesOfFreedomFramework[CurrentSystem])),
        (double)GetAverageProperty(TemperatureFrameworkAverage));
    if(DegreesOfFreedomAdsorbates[CurrentSystem]>0)
        fprintf(FilePtr,"Temperature Adsorbates: % 8.3lf (avg. % 8.3lf), Translational (avg. % 8.3lf), Rotational (avg. % 8.3lf)\n",
        (double)(2.0*UAdsorbateKinetic[CurrentSystem]/(K_B*DegreesOfFreedomAdsorbates[CurrentSystem])),
        (double)GetAverageProperty(TemperatureAdsorbatesAverage),
        (double)GetAverageProperty(TemperatureTranslationAdsorbateAverage),
        (double)GetAverageProperty(TemperatureRotationAdsorbateAverage));
    if(DegreesOfFreedomCations[CurrentSystem]>0)
      fprintf(FilePtr,"Temperature Cations:    % 8.3lf (avg. % 8.3lf)\n",
        (double)(2.0*UCationKinetic[CurrentSystem]/(K_B*DegreesOfFreedomCations[CurrentSystem])),
        (double)GetAverageProperty(TemperatureCationsAverage));
    fprintf(FilePtr,"Cell temperature: % 8.3lf (avg. % 8.3lf)\n",GetCellTemperature(),GetAverageCellTemperature());
    fprintf(FilePtr,"Current total kinetic energy:       % 22.10lf [K]\n",(double)UKinetic[CurrentSystem]);
    fprintf(FilePtr,"Current total Nose-Hoover energy:   % 22.10lf [K]\n",(double)UNoseHoover[CurrentSystem]);
  }
  fprintf(FilePtr,"Current total potential energy:       % 22.10lf [K]  (avg. % 22.10lf)\n",
      (double)UTotal[CurrentSystem]*ENERGY_TO_KELVIN,(double)GetAverageProperty(UTotalAverage)*ENERGY_TO_KELVIN);

  fprintf(FilePtr,"\tCurrent Host-Host energy:           % 22.10lf [K]  (avg. % 22.10lf)\n",
      (double)UHostHost[CurrentSystem]*ENERGY_TO_KELVIN,(double)GetAverageProperty(UHostHostAverage)*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tCurrent Host-Adsorbate energy:      % 22.10lf [K]  (avg. % 22.10lf)\n",
      (double)UHostAdsorbate[CurrentSystem]*ENERGY_TO_KELVIN,(double)GetAverageProperty(UHostAdsorbateAverage)*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tCurrent Host-Cation energy:         % 22.10lf [K]  (avg. % 22.10lf)\n",
      (double)UHostCation[CurrentSystem]*ENERGY_TO_KELVIN,(double)GetAverageProperty(UHostCationAverage)*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tCurrent Adsorbate-Adsorbate energy: % 22.10lf [K]  (avg. % 22.10lf)\n",
      (double)UAdsorbateAdsorbate[CurrentSystem]*ENERGY_TO_KELVIN,(double)GetAverageProperty(UAdsorbateAdsorbateAverage)*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tCurrent Cation-Cation energy:       % 22.10lf [K]  (avg. % 22.10lf)\n",
      (double)UCationCation[CurrentSystem]*ENERGY_TO_KELVIN,(double)GetAverageProperty(UCationCationAverage)*ENERGY_TO_KELVIN);
  fprintf(FilePtr,"\tCurrent Adsorbate-Cation energy:    % 22.10lf [K]  (avg. % 22.10lf)\n",
      (double)UAdsorbateCation[CurrentSystem]*ENERGY_TO_KELVIN,(double)GetAverageProperty(UAdsorbateCationAverage)*ENERGY_TO_KELVIN);

  if(ComputePolarization)
  {
    fprintf(FilePtr,"\tCurrent polarization energy:        % 22.10lf [K]  (avg. % 22.10lf)\n",(double)(UHostPolarization[CurrentSystem]+
         UAdsorbatePolarization[CurrentSystem]+UCationPolarization[CurrentSystem])*ENERGY_TO_KELVIN,
         (double)(GetAverageProperty(UHostPolarizationAverage)+GetAverageProperty(UAdsorbatePolarizationAverage)+
                 GetAverageProperty(UCationPolarizationAverage))*ENERGY_TO_KELVIN);
    fprintf(FilePtr,"\tCurrent back-polarization energy:   % 22.10lf [K]  (avg. % 22.10lf)\n",(double)(UHostBackPolarization[CurrentSystem]+
         UAdsorbateBackPolarization[CurrentSystem]+UCationBackPolarization[CurrentSystem])*ENERGY_TO_KELVIN,
         (double)(GetAverageProperty(UHostBackPolarizationAverage)+GetAverageProperty(UAdsorbateBackPolarizationAverage)+
                 GetAverageProperty(UCationBackPolarizationAverage))*ENERGY_TO_KELVIN);
  }

  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
  {
    if(Framework[CurrentSystem].NumberOfBondsDefinitions>0)
      fprintf(FilePtr,"\t\tCurrent Host-Bond energy:             % 22.10lf [K]  (avg. % 22.10lf)\n",
          (double)UHostBond[CurrentSystem]*ENERGY_TO_KELVIN,(double)GetAverageProperty(UHostBondAverage)*ENERGY_TO_KELVIN);
    if(Framework[CurrentSystem].NumberOfUreyBradleyDefinitions>0)
      fprintf(FilePtr,"\t\tCurrent Host-UreyBradly energy:       % 22.10lf [K]  (avg. % 22.10lf)\n",
          (double)UHostUreyBradley[CurrentSystem]*ENERGY_TO_KELVIN,(double)GetAverageProperty(UHostUreyBradleyAverage)*ENERGY_TO_KELVIN);
    if(Framework[CurrentSystem].NumberOfBendDefinitions>0)
      fprintf(FilePtr,"\t\tCurrent Host-Bend energy:             % 22.10lf [K]  (avg. % 22.10lf)\n",
          (double)UHostBend[CurrentSystem]*ENERGY_TO_KELVIN,(double)GetAverageProperty(UHostBendAverage)*ENERGY_TO_KELVIN);
    if(Framework[CurrentSystem].NumberOfInversionBendDefinitions>0)
      fprintf(FilePtr,"\t\tCurrent Host-Inversion Bend energy:   % 22.10lf [K]  (avg. % 22.10lf)\n",
          (double)UHostInversionBend[CurrentSystem]*ENERGY_TO_KELVIN,(double)GetAverageProperty(UHostInversionBendAverage)*ENERGY_TO_KELVIN);
    if(Framework[CurrentSystem].NumberOfTorsionDefinitions>0)
      fprintf(FilePtr,"\t\tCurrent Host-Torsion energy:          % 22.10lf [K]  (avg. % 22.10lf)\n",
          (double)UHostTorsion[CurrentSystem]*ENERGY_TO_KELVIN,(double)GetAverageProperty(UHostTorsionAverage)*ENERGY_TO_KELVIN);
    if(Framework[CurrentSystem].NumberOfImproperTorsionDefinitions>0)
      fprintf(FilePtr,"\t\tCurrent Host-Improper torsion energy: % 22.10lf [K]  (avg. % 22.10lf)\n",
          (double)UHostImproperTorsion[CurrentSystem]*ENERGY_TO_KELVIN,(double)GetAverageProperty(UHostImproperTorsionAverage)*ENERGY_TO_KELVIN);
    if(Framework[CurrentSystem].NumberOfBondBondDefinitions>0)
      fprintf(FilePtr,"\t\tCurrent Host-Bond/Bond energy:        % 22.10lf [K]  (avg. % 22.10lf)\n",
          (double)UHostBondBond[CurrentSystem]*ENERGY_TO_KELVIN,(double)GetAverageProperty(UHostBondBondAverage)*ENERGY_TO_KELVIN);
    if(Framework[CurrentSystem].NumberOfBendBendDefinitions>0)
      fprintf(FilePtr,"\t\tCurrent Host-Bend/Bend energy:        % 22.10lf [K]  (avg. % 22.10lf)\n",
          (double)UHostBendBend[CurrentSystem]*ENERGY_TO_KELVIN,(double)GetAverageProperty(UHostBendBendAverage)*ENERGY_TO_KELVIN);
    if(Framework[CurrentSystem].NumberOfBondBendDefinitions>0)
      fprintf(FilePtr,"\t\tCurrent Host-Bond/Bend energy:        % 22.10lf [K]  (avg. % 22.10lf)\n",
          (double)UHostBondBend[CurrentSystem]*ENERGY_TO_KELVIN,(double)GetAverageProperty(UHostBondBendAverage)*ENERGY_TO_KELVIN);
    if(Framework[CurrentSystem].NumberOfBondTorsionDefinitions>0)
      fprintf(FilePtr,"\t\tCurrent Host-Bond/Torsion energy:     % 22.10lf [K]  (avg. % 22.10lf)\n",
          (double)UHostBondTorsion[CurrentSystem]*ENERGY_TO_KELVIN,(double)GetAverageProperty(UHostBondTorsionAverage)*ENERGY_TO_KELVIN);
    if(Framework[CurrentSystem].NumberOfBendTorsionDefinitions>0)
      fprintf(FilePtr,"\t\tCurrent Host-Bend/Torsion energy:     % 22.10lf [K]  (avg. % 22.10lf)\n\n",
          (double)UHostBendTorsion[CurrentSystem]*ENERGY_TO_KELVIN,(double)GetAverageProperty(UHostBendTorsionAverage)*ENERGY_TO_KELVIN);
  }

  fprintf(FilePtr,"\n");
  PrintWarningStatus();
  fprintf(FilePtr,"\n\n");
  fflush(FilePtr);
}

void PrintProperty(FILE *FilePtr,char *string,char *units,REAL conv_factor,REAL **Property)
{
  int i;
  REAL sum,sum2,tmp;

  fprintf(FilePtr,"\n");
  fprintf(FilePtr,"%s",string);
  sum=sum2=0.0;
  for(i=0;i<NR_BLOCKS;i++)
  {
    if(BlockCount[CurrentSystem][i]>0.0)
    {
      tmp=conv_factor*(Property[CurrentSystem][i]/BlockCount[CurrentSystem][i]);
      sum+=tmp;
      sum2+=SQR(tmp);
      fprintf(FilePtr,"\tBlock[%2d] %18.5lf %s\n",i,(double)tmp,units);
    }
    else
      fprintf(FilePtr,"\tBlock[%2d] %18.5lf %s\n",i,(double)0.0,units);
  }
  fprintf(FilePtr,"\t------------------------------------------------------------------------------\n");
  tmp=2.0*sqrt(fabs((sum2/(REAL)NR_BLOCKS)-SQR(sum)/(REAL)SQR(NR_BLOCKS)));
  fprintf(FilePtr,"\tAverage   %18.5lf %s +/- %18.5lf %s\n",(double)(sum/(REAL)NR_BLOCKS),units,(double)tmp,units);

}

void PrintEnergies(FILE *FilePtr,char *string,char *units,REAL conv_factor,REAL **Property,
                   REAL **PropertyVDW,REAL **PropertyCoulomb)
{
  int i;
  REAL sum,sum2,tmp,sum_vdw,sum_vdw2,sum_coul,sum_coul2,tmp_vdw,tmp_coul;

  fprintf(FilePtr,"\n");
  fprintf(FilePtr,"%s",string);
  sum=sum_vdw=sum_coul=0.0;
  sum2=sum_vdw2=sum_coul2=0.0;
  for(i=0;i<NR_BLOCKS;i++)
  {
    if(BlockCount[CurrentSystem][i]>0.0)
    {
      tmp=conv_factor*(Property[CurrentSystem][i]/BlockCount[CurrentSystem][i]);
      tmp_vdw=conv_factor*(PropertyVDW[CurrentSystem][i]/BlockCount[CurrentSystem][i]);
      tmp_coul=conv_factor*(PropertyCoulomb[CurrentSystem][i]/BlockCount[CurrentSystem][i]);
      sum+=tmp;
      sum_vdw+=tmp_vdw;
      sum_coul+=tmp_coul;
      sum2+=SQR(tmp);
      sum_vdw2+=SQR(tmp_vdw);
      sum_coul2+=SQR(tmp_coul);

      fprintf(FilePtr,"\tBlock[%2d] %-18.5lf Van der Waals: %-18.5lf Coulomb: %-18.5lf %s\n",
        i,(double)tmp,(double)tmp_vdw,(double)tmp_coul,units);
    }
    else
      fprintf(FilePtr,"\tBlock[%2d] %-18.5lf Van der Waals: %-18.5lf Coulomb: %-18.5lf %s\n",i,(double)0.0,(double)0.0,(double)0.0,units);
  }
  fprintf(FilePtr,"\t------------------------------------------------------------------------------\n");
  tmp=2.0*sqrt(fabs((sum2/(REAL)NR_BLOCKS)-SQR(sum)/(REAL)SQR(NR_BLOCKS)));
  tmp_vdw=2.0*sqrt(fabs((sum_vdw2/(REAL)NR_BLOCKS)-SQR(sum_vdw)/(REAL)SQR(NR_BLOCKS)));
  tmp_coul=2.0*sqrt(fabs((sum_coul2/(REAL)NR_BLOCKS)-SQR(sum_coul)/(REAL)SQR(NR_BLOCKS)));
  fprintf(FilePtr,"\tAverage   %-18.5lf Van der Waals: %-18lf Coulomb: %-18.5lf %s\n",
         (double)(sum/(REAL)NR_BLOCKS),(double)(sum_vdw/(REAL)NR_BLOCKS),(double)(sum_coul/(REAL)NR_BLOCKS),units);
  fprintf(FilePtr,"\t      +/- %-18.5lf            +/- %-18lf      +/- %-18.5lf %s\n",
         (double)tmp,(double)tmp_vdw,(double)tmp_coul,units);
}


void PrintAverageTotalSystemEnergiesMC(FILE *FilePtr)
{
  int i,j,nr;
  REAL sum,sum_vdw,sum_coul,tmp;
  REAL sum2,sum_vdw2,sum_coul2;
  REAL HeatOfDesorption[NR_BLOCKS];
  REAL CationMass,Mass;
  REAL vol,dipole_norm,dipole_norm_squared;
  REAL nr_molecules;
  REAL AverageVolume;
  REAL HV,V,V2,H,H2,N,T;
  REAL FrameworkDensity,Temperature;
  VECTOR sumv1,sumv2,av;

  AverageVolume=GetAverageProperty(VolumeAverage);

  nr_molecules=0.0;
  for(j=0;j<NumberOfComponents;j++)
  {
    fprintf(FilePtr,"Component %d [%s]\n",j,Components[j].Name);
    fprintf(FilePtr,"-------------------------------------------------------------\n");
    sum=sum_vdw=sum_coul=0.0;
    sum2=sum_vdw2=sum_coul2=0.0;
    for(i=0;i<NR_BLOCKS;i++)
    {
      if(BlockCount[CurrentSystem][i]>0.0)
      {
        sum+=NumberOfMoleculesPerComponentAverage[CurrentSystem][j][i]/BlockCount[CurrentSystem][i];
        sum2+=SQR(NumberOfMoleculesPerComponentAverage[CurrentSystem][j][i]/BlockCount[CurrentSystem][i]);

        fprintf(FilePtr,"\tBlock[%2d] %-18.5lf [-]\n",i,
                (double)(NumberOfMoleculesPerComponentAverage[CurrentSystem][j][i]/BlockCount[CurrentSystem][i]));
      }
      else
        fprintf(FilePtr,"\tBlock[%2d] %-18.5lf [-]\n",i,(double)0.0);
    }
    nr_molecules+=sum/(REAL)NR_BLOCKS;
    fprintf(FilePtr,"\n");
  }

  sum=sum2=0.0;
  CationMass=0.0;
  for(j=0;j<NumberOfComponents;j++)
  {
    if(Components[j].ExtraFrameworkMolecule)
    {
      for(i=0;i<NR_BLOCKS;i++)
        if(BlockCount[CurrentSystem][i]>0.0)
        {
          sum+=NumberOfMoleculesPerComponentAverage[CurrentSystem][j][i];
          sum2+=BlockCount[CurrentSystem][i];
        }
      CationMass=(sum/sum2)*Components[j].Mass;
    }
  }

  fprintf(FilePtr,"\n\n\n\n\n");
  fprintf(FilePtr,"Average properties of the system[%d]:\n",CurrentSystem);
  fprintf(FilePtr,"========================================================================\n");

  // Temperature
  fprintf(FilePtr,"\n");
  fprintf(FilePtr,"Average temperature:\n");
  fprintf(FilePtr,"====================\n");
  sum=sum2=0.0;
  for(i=0;i<NR_BLOCKS;i++)
  {
    if(BlockCount[CurrentSystem][i]>0.0)
    {
      tmp=TemperatureAverage[CurrentSystem][i]/BlockCount[CurrentSystem][i];
      sum+=tmp;
      sum2+=SQR(tmp);
      fprintf(FilePtr,"\tBlock[%2d] %18.5lf [K]\n",i,(double)tmp);
    }
    else
      fprintf(FilePtr,"\tBlock[%2d] %18.5lf [K]\n",i,(double)0.0);
  }
  fprintf(FilePtr,"\t------------------------------------------------------------------------------\n");
  tmp=2.0*sqrt(fabs((sum2/(REAL)NR_BLOCKS)-SQR(sum)/(REAL)SQR(NR_BLOCKS)));
  fprintf(FilePtr,"\tAverage   %18.5lf [K] +/- %18.5lf [K]\n",(double)(sum/(REAL)NR_BLOCKS),(double)tmp);

  // Pressure
  fprintf(FilePtr,"\n");
  switch(SimulationType)
  {
    case MONTE_CARLO:
      fprintf(FilePtr,"Average Pressure:\n");
      fprintf(FilePtr,"=================\n");
      sum=sum2=0.0;
      for(i=0;i<NR_BLOCKS;i++)
      {
        if(BlockCount[CurrentSystem][i]>0.0)
        {
          tmp=((MolecularStressTensorAverage[CurrentSystem][i].ax+MolecularStressTensorAverage[CurrentSystem][i].by+MolecularStressTensorAverage[CurrentSystem][i].cz)/
               (3.0*BlockCount[CurrentSystem][i]))*PRESSURE_CONVERSION_FACTOR;

          sum+=tmp;
          sum2+=SQR(tmp);
          fprintf(FilePtr,"\tBlock[%2d] %18.5lf [Pa]\n",i,(double)tmp);
        }
        else
          fprintf(FilePtr,"\tBlock[%2d] %18.5lf [Pa]\n",i,(double)0.0);
      }
      fprintf(FilePtr,"\t------------------------------------------------------------------------------\n");
      tmp=2.0*sqrt(fabs((sum2/(REAL)NR_BLOCKS)-SQR(sum)/(REAL)SQR(NR_BLOCKS)));
      fprintf(FilePtr,"\tAverage   %18.5lf [Pa] +/- %18.5lf [Pa]\n",(double)(sum/(REAL)NR_BLOCKS),(double)tmp);
      fprintf(FilePtr,"\tAverage   %18.5lf [bar] +/- %18.5lf [bar]\n",(double)((sum/(REAL)NR_BLOCKS)*PA_TO_BAR),(double)(tmp*PA_TO_BAR));
      fprintf(FilePtr,"\tAverage   %18.5lf [atm] +/- %18.5lf [atm]\n",(double)((sum/(REAL)NR_BLOCKS)*PA_TO_ATM),(double)(tmp*PA_TO_ATM));
      fprintf(FilePtr,"\tAverage   %18.5lf [Torr] +/- %18.5lf [Torr]\n",(double)((sum/(REAL)NR_BLOCKS)*PA_TO_TORR),(double)(tmp*PA_TO_TORR));
      break;
    case MOLECULAR_DYNAMICS:
      fprintf(FilePtr,"Average Pressure:\n");
      fprintf(FilePtr,"=================\n");
      sum=sum2=0.0;
      for(i=0;i<NR_BLOCKS;i++)
      {
        if(BlockCount[CurrentSystem][i]>0.0)
        {
          // NEW
          //tmp=(MolecularPressureAverage[CurrentSystem][i]/BlockCount[CurrentSystem][i])*PRESSURE_CONVERSION_FACTOR;
          tmp=((StressTensorAverage[CurrentSystem][i].ax+StressTensorAverage[CurrentSystem][i].by+StressTensorAverage[CurrentSystem][i].cz)/
               (3.0*BlockCount[CurrentSystem][i]))*PRESSURE_CONVERSION_FACTOR;
          sum+=tmp;
          sum2+=SQR(tmp);
          fprintf(FilePtr,"\tBlock[%2d] %18.5lf [Pa]\n",i,(double)tmp);
        }
        else
          fprintf(FilePtr,"\tBlock[%2d] %18.5lf [Pa]\n",i,(double)0.0);
      }
      fprintf(FilePtr,"\t------------------------------------------------------------------------------\n");
      tmp=2.0*sqrt(fabs((sum2/(REAL)NR_BLOCKS)-SQR(sum)/(REAL)SQR(NR_BLOCKS)));
      fprintf(FilePtr,"\tAverage   %18.5lf [Pa] +/- %18.5lf [Pa]\n",(double)(sum/(REAL)NR_BLOCKS),(double)tmp);
      fprintf(FilePtr,"\tAverage   %18.5lf [bar] +/- %18.5lf [bar]\n",(double)((sum/(REAL)NR_BLOCKS)*PA_TO_BAR),(double)(tmp*PA_TO_BAR));
      fprintf(FilePtr,"\tAverage   %18.5lf [atm] +/- %18.5lf [atm]\n",(double)((sum/(REAL)NR_BLOCKS)*PA_TO_ATM),(double)(tmp*PA_TO_ATM));
      fprintf(FilePtr,"\tAverage   %18.5lf [Torr] +/- %18.5lf [Torr]\n",(double)((sum/(REAL)NR_BLOCKS)*PA_TO_TORR),(double)(tmp*PA_TO_TORR));
      break;
    default:
      break;
  }

  // Volume
  fprintf(FilePtr,"\n");
  fprintf(FilePtr,"Average Volume:\n");
  fprintf(FilePtr,"=================\n");
  sum=sum2=0.0;
  for(i=0;i<NR_BLOCKS;i++)
  {
    if(BlockCount[CurrentSystem][i]>0.0)
    {
      tmp=(VolumeAverage[CurrentSystem][i]/BlockCount[CurrentSystem][i]);
      sum+=tmp;
      sum2+=SQR(tmp);
      fprintf(FilePtr,"\tBlock[%2d] %18.5lf [A^3]\n",i,(double)tmp);
    }
    else
      fprintf(FilePtr,"\tBlock[%2d] %18.5lf [A^3]\n",i,(double)0.0);
  }
  fprintf(FilePtr,"\t------------------------------------------------------------------------------\n");
  tmp=2.0*sqrt(fabs((sum2/(REAL)NR_BLOCKS)-SQR(sum)/(REAL)SQR(NR_BLOCKS)));
  fprintf(FilePtr,"\tAverage   %18.5lf [A^3] +/- %18.5lf [A^3]\n",(double)(sum/(REAL)NR_BLOCKS),(double)tmp);
  vol=sum/(REAL)NR_BLOCKS;

  // Box-lengths
  fprintf(FilePtr,"\n");
  fprintf(FilePtr,"Average Box-lengths:\n");
  fprintf(FilePtr,"====================\n");
  sum=sum2=0.0;
  for(i=0;i<NR_BLOCKS;i++)
  {
    if(BlockCount[CurrentSystem][i]>0.0)
    {
      tmp=(BoxAverage[CurrentSystem][i].x/BlockCount[CurrentSystem][i]);
      sum+=tmp;
      sum2+=SQR(tmp);
      fprintf(FilePtr,"\tBlock[%2d] %18.5lf [A^3]\n",i,(double)tmp);
    }
    else
      fprintf(FilePtr,"\tBlock[%2d] %18.5lf [A^3]\n",i,(double)0.0);
  }
  fprintf(FilePtr,"\t------------------------------------------------------------------------------\n");
  tmp=2.0*sqrt(fabs((sum2/(REAL)NR_BLOCKS)-SQR(sum)/(REAL)SQR(NR_BLOCKS)));
  fprintf(FilePtr,"\tAverage Box.ax  %18.5lf [A^3] +/- %18.5lf [A^3]\n\n",(double)(sum/(REAL)NR_BLOCKS),(double)tmp);

  sum=sum2=0.0;
  for(i=0;i<NR_BLOCKS;i++)
  {
    if(BlockCount[CurrentSystem][i]>0.0)
    {
      tmp=(BoxAverage[CurrentSystem][i].y/BlockCount[CurrentSystem][i]);
      sum+=tmp;
      sum2+=SQR(tmp);
      fprintf(FilePtr,"\tBlock[%2d] %18.5lf [A^3]\n",i,(double)tmp);
    }
    else
      fprintf(FilePtr,"\tBlock[%2d] %18.5lf [A^3]\n",i,(double)0.0);
  }
  fprintf(FilePtr,"\t------------------------------------------------------------------------------\n");
  tmp=2.0*sqrt(fabs((sum2/(REAL)NR_BLOCKS)-SQR(sum)/(REAL)SQR(NR_BLOCKS)));
  fprintf(FilePtr,"\tAverage Box.by  %18.5lf [A^3] +/- %18.5lf [A^3]\n\n",(double)(sum/(REAL)NR_BLOCKS),(double)tmp);

  sum=sum2=0.0;
  for(i=0;i<NR_BLOCKS;i++)
  {
    if(BlockCount[CurrentSystem][i]>0.0)
    {
      tmp=(BoxAverage[CurrentSystem][i].z/BlockCount[CurrentSystem][i]);
      sum+=tmp;
      sum2+=SQR(tmp);
      fprintf(FilePtr,"\tBlock[%2d] %18.5lf [A^3]\n",i,(double)tmp);
    }
    else
      fprintf(FilePtr,"\tBlock[%2d] %18.5lf [A^3]\n",i,(double)0.0);
  }
  fprintf(FilePtr,"\t------------------------------------------------------------------------------\n");
  tmp=2.0*sqrt(fabs((sum2/(REAL)NR_BLOCKS)-SQR(sum)/(REAL)SQR(NR_BLOCKS)));
  fprintf(FilePtr,"\tAverage Box.cz  %18.5lf [A^3] +/- %18.5lf [A^3]\n",(double)(sum/(REAL)NR_BLOCKS),(double)tmp);

  sum=sum2=0.0;
  for(i=0;i<NR_BLOCKS;i++)
  {
    if(BlockCount[CurrentSystem][i]>0.0)
    {
      tmp=(AlphaAngleAverage[CurrentSystem][i]/BlockCount[CurrentSystem][i]);
      sum+=tmp;
      sum2+=SQR(tmp);
      fprintf(FilePtr,"\tBlock[%2d] %18.5lf [A^3]\n",i,(double)tmp);
    }
    else
      fprintf(FilePtr,"\tBlock[%2d] %18.5lf [A^3]\n",i,(double)0.0);
  }
  fprintf(FilePtr,"\t------------------------------------------------------------------------------\n");
  tmp=2.0*sqrt(fabs((sum2/(REAL)NR_BLOCKS)-SQR(sum)/(REAL)SQR(NR_BLOCKS)));
  fprintf(FilePtr,"\tAverage alpha angle  %18.5lf [degrees] +/- %18.5lf [degrees]\n",(double)(sum/(REAL)NR_BLOCKS),(double)tmp);


  sum=sum2=0.0;
  for(i=0;i<NR_BLOCKS;i++)
  {
    if(BlockCount[CurrentSystem][i]>0.0)
    {
      tmp=(BetaAngleAverage[CurrentSystem][i]/BlockCount[CurrentSystem][i]);
      sum+=tmp;
      sum2+=SQR(tmp);
      fprintf(FilePtr,"\tBlock[%2d] %18.5lf [A^3]\n",i,(double)tmp);
    }
    else
      fprintf(FilePtr,"\tBlock[%2d] %18.5lf [A^3]\n",i,(double)0.0);
  }
  fprintf(FilePtr,"\t------------------------------------------------------------------------------\n");
  tmp=2.0*sqrt(fabs((sum2/(REAL)NR_BLOCKS)-SQR(sum)/(REAL)SQR(NR_BLOCKS)));
  fprintf(FilePtr,"\tAverage beta angle  %18.5lf [degrees] +/- %18.5lf [degrees]\n",(double)(sum/(REAL)NR_BLOCKS),(double)tmp);


  sum=sum2=0.0;
  for(i=0;i<NR_BLOCKS;i++)
  {
    if(BlockCount[CurrentSystem][i]>0.0)
    {
      tmp=(GammaAngleAverage[CurrentSystem][i]/BlockCount[CurrentSystem][i]);
      sum+=tmp;
      sum2+=SQR(tmp);
      fprintf(FilePtr,"\tBlock[%2d] %18.5lf [A^3]\n",i,(double)tmp);
    }
    else
      fprintf(FilePtr,"\tBlock[%2d] %18.5lf [A^3]\n",i,(double)0.0);
  }
  fprintf(FilePtr,"\t------------------------------------------------------------------------------\n");
  tmp=2.0*sqrt(fabs((sum2/(REAL)NR_BLOCKS)-SQR(sum)/(REAL)SQR(NR_BLOCKS)));
  fprintf(FilePtr,"\tAverage gamma angle  %18.5lf [degrees] +/- %18.5lf [degrees]\n",(double)(sum/(REAL)NR_BLOCKS),(double)tmp);


  // Average Surface area
  fprintf(FilePtr,"\n");
  fprintf(FilePtr,"Average Surface Area:\n");
  fprintf(FilePtr,"=====================\n");
  for(j=0;j<NumberOfComponents;j++)
  {
    sum=sum2=0.0;
    for(i=0;i<NR_BLOCKS;i++)
    {
      if(SurfaceAreaCount[CurrentSystem][i]>0.0)
      {
        tmp=SurfaceAreaFrameworkAverage[CurrentSystem][i]/SurfaceAreaCount[CurrentSystem][i];
        sum+=tmp;
        sum2+=SQR(tmp);
        fprintf(FilePtr,"\tBlock[%2d] %-lf [-]\n",i,(double)tmp);
      }
      else
        fprintf(FilePtr,"\tBlock[%2d] %-lf [-]\n",i,(double)0.0);
    }
    fprintf(FilePtr,"\t------------------------------------------------------------------------------\n");
    fprintf(FilePtr,"\tSurface area:   %lf +/- %lf [A^2]\n",
      (double)(sum/(REAL)NR_BLOCKS),(double)(2.0*sqrt(fabs((sum2/(REAL)NR_BLOCKS)-SQR(sum)/(REAL)SQR(NR_BLOCKS)))));
    fprintf(FilePtr,"\tSurface area:   %lf +/- %lf [m^2/g]\n",
      (double)((sum/(REAL)NR_BLOCKS)*SQR(ANGSTROM)*
        AVOGADRO_CONSTANT/(Framework[CurrentSystem].FrameworkMass)),
      (double)(2.0*sqrt(fabs((sum2/(REAL)NR_BLOCKS)-SQR(sum)/(REAL)SQR(NR_BLOCKS)))*SQR(ANGSTROM)*
        AVOGADRO_CONSTANT/(Framework[CurrentSystem].FrameworkMass)));
    fprintf(FilePtr,"\tSurface area:   %lf +/- %g [m^2/cm^3]\n",
      (double)((sum/(REAL)NR_BLOCKS)*1.0e4/Volume[0]),
      (double)(2.0*sqrt(fabs((sum2/(REAL)NR_BLOCKS)-SQR(sum)/(REAL)SQR(NR_BLOCKS)))*1.0e4/Volume[0]));
  }

  // Density
  fprintf(FilePtr,"\n");
  fprintf(FilePtr,"Average Density:\n");
  fprintf(FilePtr,"=================\n");
  sum=sum2=0.0;
  Mass=GetTotalAdsorbateMass()+GetTotalCationMass();
  for(i=0;i<NR_BLOCKS;i++)
  {
    if(BlockCount[CurrentSystem][i]>0.0)
    {
      tmp=(DensityAverage[CurrentSystem][i]/BlockCount[CurrentSystem][i])*DENSITY_CONVERSION_FACTOR;
      sum+=tmp;
      sum2+=SQR(tmp);
      fprintf(FilePtr,"\tBlock[%2d] %18.5lf [kg/m^3]\n",i,(double)tmp);
    }
    else
      fprintf(FilePtr,"\tBlock[%2d] %18.5lf [kg/m^3]\n",i,(double)0.0);
  }
  fprintf(FilePtr,"\t------------------------------------------------------------------------------\n");
  tmp=2.0*sqrt(fabs((sum2/(REAL)NR_BLOCKS)-SQR(sum)/(REAL)SQR(NR_BLOCKS)));
  fprintf(FilePtr,"\tAverage   %18.5lf [kg/m^3] +/- %18.5lf [kg/m^3]\n",(double)(sum/(REAL)NR_BLOCKS),(double)tmp);

  for(j=0;j<NumberOfComponents;j++)
  {
    fprintf(FilePtr,"\tComponent %d [%s]\n",j,Components[j].Name);
    fprintf(FilePtr,"\t-------------------------------------------------------------\n");

    sum=sum2=0.0;
    for(i=0;i<NR_BLOCKS;i++)
    {
      if(BlockCount[CurrentSystem][i]>0.0)
      {
        tmp=(DensityPerComponentAverage[CurrentSystem][j][i]/BlockCount[CurrentSystem][i])*DENSITY_CONVERSION_FACTOR;
        sum+=tmp;
        sum2+=SQR(tmp);

        fprintf(FilePtr,"\t\tBlock[%2d] %18.5lf [kg/m^3]\n",i,(double)tmp);
      }
      else
        fprintf(FilePtr,"\t\tBlock[%2d] %-18.5lf [-]\n",i,(double)0.0);
    }
    fprintf(FilePtr,"\t\t------------------------------------------------------------------------------\n");
    tmp=2.0*sqrt(fabs((sum2/(REAL)NR_BLOCKS)-SQR(sum)/(REAL)SQR(NR_BLOCKS)));
    fprintf(FilePtr,"\t\tAverage   %18.5lf [kg/m^3] +/- %18.5lf [kg/m^3]\n",(double)(sum/(REAL)NR_BLOCKS),(double)tmp);
  }

  // compressibility
  fprintf(FilePtr,"\n");
  fprintf(FilePtr,"Average compressibility Z:\n");
  fprintf(FilePtr,"=========================\n");
  sum=sum2=0.0;
  for(i=0;i<NR_BLOCKS;i++)
  {
    if(BlockCount[CurrentSystem][i]>0.0)
    {
      tmp=CompressibilityAverage[CurrentSystem][i]/BlockCount[CurrentSystem][i];
      sum+=tmp;
      sum2+=SQR(tmp);
      fprintf(FilePtr,"\tBlock[%2d] %18.5lf [-]\n",i,(double)tmp);
    }
    else
      fprintf(FilePtr,"\tBlock[%2d] %18.5lf [-]\n",i,(double)0.0);
  }
  fprintf(FilePtr,"\t------------------------------------------------------------------------------\n");
  tmp=2.0*sqrt(fabs((sum2/(REAL)NR_BLOCKS)-SQR(sum)/(REAL)SQR(NR_BLOCKS)));
  fprintf(FilePtr,"\tAverage   %18.5lf [-] +/- %18.5lf [-]\n",(double)(sum/(REAL)NR_BLOCKS),(double)tmp);

/*
  // Dielectric constant
  fprintf(FilePtr,"\n");
  fprintf(FilePtr,"Average dielectric constant:\n");
  fprintf(FilePtr,"============================\n");
  sum=sum2=0.0;
  dipole_norm=dipole_norm_squared=0.0;
  for(i=0;i<NR_BLOCKS;i++)
  {
    if(BlockCount[CurrentSystem][i]>0.0)
    {
      dipole_norm=TotalSystemNormDipoleAverage[CurrentSystem][i]/BlockCount[CurrentSystem][i];
      dipole_norm_squared=TotalSystemNormDipoleSquaredAverage[CurrentSystem][i]/BlockCount[CurrentSystem][i];
      tmp=1.0+(4.0*M_PI*Beta[CurrentSystem]/(3.0*vol))*(dipole_norm_squared-SQR(dipole_norm))*
          DIELECTRIC_CONSTANT_CONVERSION_FACTOR/DIELECTRIC_CONSTANT_VACUUM;
      sum+=tmp;
      sum2+=SQR(tmp);
      fprintf(FilePtr,"\tBlock[%2d] %18.5lf [-]\n",i,(double)tmp);
    }
    else
      fprintf(FilePtr,"\tBlock[%2d] %18.5lf [-]\n",i,(double)0.0);
  }
  fprintf(FilePtr,"\t------------------------------------------------------------------------------\n");
  tmp=2.0*sqrt(fabs((sum2/(REAL)NR_BLOCKS)-SQR(sum)/(REAL)SQR(NR_BLOCKS)));
  fprintf(FilePtr,"\tAverage   %18.5lf [-] +/- %18.5lf [-]\n",(double)(sum/(REAL)NR_BLOCKS),(double)tmp);

  // isothermal compressibility
  fprintf(FilePtr,"\n");
  fprintf(FilePtr,"Average Isothermal Compressibility:\n");
  fprintf(FilePtr,"===================================\n");
  sum=sum2=0.0;
  for(i=0;i<NR_BLOCKS;i++)
  {
    if(BlockCount[CurrentSystem][i]>0.0)
    {
      T=therm_baro_stats.ExternalTemperature[CurrentSystem];
      V=VolumeAverage[CurrentSystem][i]/BlockCount[CurrentSystem][i];
      V2=VolumeSquaredAverage[CurrentSystem][i]/BlockCount[CurrentSystem][i];
      tmp=1e12*ISOTHERMAL_COMPRESSIBILITY_CONVERSION_FACTOR*(V2-SQR(V))/(K_B*T*V);
      sum+=tmp;
      sum2+=SQR(tmp);
      fprintf(FilePtr,"\tBlock[%2d] %18.5lf [Pa^-1]\n",i,(double)tmp);
    }
    else
      fprintf(FilePtr,"\tBlock[%2d] %18.5lf [Pa^-1]\n",i,(double)0.0);
  }
  fprintf(FilePtr,"\t------------------------------------------------------------------------------\n");
  tmp=2.0*sqrt(fabs((sum2/(REAL)NR_BLOCKS)-SQR(sum)/(REAL)SQR(NR_BLOCKS)));
  fprintf(FilePtr,"\tAverage   %18.5lf [10^12 Pa^-1] +/- %18.5lf [10^12 Pa^-1]\n",(double)(sum/(REAL)NR_BLOCKS),(double)tmp);
  tmp=2.0*sqrt(fabs((sum2/(REAL)NR_BLOCKS)-SQR(sum)/(REAL)SQR(NR_BLOCKS)));
  fprintf(FilePtr,"\tAverage   %18.5lf [10^6 bar^-1] +/- %18.5lf [10^6 bar^-1]\n",
          (double)(sum/(1e6*PA_TO_BAR*(REAL)NR_BLOCKS)),(double)(tmp/(1e6*PA_TO_BAR)));
  tmp=2.0*sqrt(fabs((sum2/(REAL)NR_BLOCKS)-SQR(sum)/(REAL)SQR(NR_BLOCKS)));
  fprintf(FilePtr,"\tAverage   %18.5lf [10^6 atm^-1] +/- %18.5lf [10^6 atm^-1]\n",
          (double)(sum/(1e6*PA_TO_ATM*(REAL)NR_BLOCKS)),(double)(tmp/(1e6*PA_TO_ATM)));

  // isothermal expansion coefficient
  fprintf(FilePtr,"\n");
  fprintf(FilePtr,"Average Volumetric Expansion Coefficient beta:\n");
  fprintf(FilePtr,"===============================================\n");
  sum=sum2=0.0;
  for(i=0;i<NR_BLOCKS;i++)
  {
    if(BlockCount[CurrentSystem][i]>0.0)
    {
      HV=EnthalpyTimesVolumeAverage[CurrentSystem][i]/BlockCount[CurrentSystem][i];
      V=VolumeAverage[CurrentSystem][i]/BlockCount[CurrentSystem][i];
      H=EnthalpyAverage[CurrentSystem][i]/BlockCount[CurrentSystem][i];
      T=therm_baro_stats.ExternalTemperature[CurrentSystem];
      tmp=1e5*VOLUMETRIC_EXPANSION_COEFFICIENT_CONVERSION_FACTOR*(HV-H*V)/(K_B*V*SQR(T));
      sum+=tmp;
      sum2+=SQR(tmp);
      fprintf(FilePtr,"\tBlock[%2d] %18.5lf [10^5 K^-1]\n",i,(double)tmp);
    }
    else
      fprintf(FilePtr,"\tBlock[%2d] %18.5lf [10^5 K^-1]\n",i,(double)0.0);
  }
  fprintf(FilePtr,"\t------------------------------------------------------------------------------\n");
  tmp=2.0*sqrt(fabs((sum2/(REAL)NR_BLOCKS)-SQR(sum)/(REAL)SQR(NR_BLOCKS)));
  fprintf(FilePtr,"\tAverage   %18.5lf [10^5 K^-1] +/- %18.5lf [10^5 K^-1]\n",(double)(sum/(REAL)NR_BLOCKS),(double)tmp);


  // heat capacity
  fprintf(FilePtr,"\n");
  fprintf(FilePtr,"Average Heat Capacity Cp:\n");
  fprintf(FilePtr,"=========================\n");
  N=NumberOfAdsorbateMolecules[CurrentSystem]+NumberOfCationMolecules[CurrentSystem];
  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
    N+=Framework[CurrentSystem].TotalNumberOfAtoms;
  sum=sum2=0.0;
  for(i=0;i<NR_BLOCKS;i++)
  {
    if(BlockCount[CurrentSystem][i]>0.0)
    {
      H2=EnthalpySquaredAverage[CurrentSystem][i]/BlockCount[CurrentSystem][i];
      H=EnthalpyAverage[CurrentSystem][i]/BlockCount[CurrentSystem][i];
      T=therm_baro_stats.ExternalTemperature[CurrentSystem];
      tmp=HEAT_CAPACITY_CONVERSION_FACTOR*((H2-SQR(H))/(N*K_B*SQR(T))+5.0*K_B/2.0);
      sum+=tmp;
      sum2+=SQR(tmp);
      fprintf(FilePtr,"\tBlock[%2d] %lf [J/mol/K]\n",i,(double)tmp);
    }
    else
      fprintf(FilePtr,"\tBlock[%2d] %lf [J/mol/K]\n",i,(double)0.0);
  }
  fprintf(FilePtr,"\t------------------------------------------------------------------------------\n");
  tmp=2.0*sqrt(fabs((sum2/(REAL)NR_BLOCKS)-SQR(sum)/(REAL)SQR(NR_BLOCKS)));
  fprintf(FilePtr,"\tAverage   %18.5lf [J/mol/K] +/- %18.5lf [J/mol/K]\n",
          (double)(sum/(REAL)NR_BLOCKS),(double)tmp);
  fprintf(FilePtr,"\tAverage   %18.5lf [cal/mol/K] +/- %18.5lf [cal/mol/K]\n",
          (double)(sum*J_TO_CAL/(REAL)NR_BLOCKS),(double)tmp*J_TO_CAL);

  // heat capacity
  fprintf(FilePtr,"\n");
  fprintf(FilePtr,"Average Heat Capacity:\n");
  fprintf(FilePtr,"======================\n");
  sum=sum2=0.0;
  N=NumberOfAdsorbateMolecules[CurrentSystem]+NumberOfCationMolecules[CurrentSystem];
  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
    N+=Framework[CurrentSystem].TotalNumberOfAtoms;
  for(i=0;i<NR_BLOCKS;i++)
  {
    if(BlockCount[CurrentSystem][i]>0.0)
    {
      H2=TotalEnergySquaredAverage[CurrentSystem][i]/BlockCount[CurrentSystem][i];
      H=TotalEnergyAverage[CurrentSystem][i]/BlockCount[CurrentSystem][i];
      T=therm_baro_stats.ExternalTemperature[CurrentSystem];
      tmp=HEAT_CAPACITY_CONVERSION_FACTOR*((H2-SQR(H))/(N*K_B*SQR(T))+3.0*K_B/2.0);
      sum+=tmp;
      sum2+=SQR(tmp);
      fprintf(FilePtr,"\tBlock[%2d] %g [J/mol/K]\n",i,(double)tmp);
    }
    else
      fprintf(FilePtr,"\tBlock[%2d] %g [J/mol/K]\n",i,(double)0.0);
  }
  fprintf(FilePtr,"\t------------------------------------------------------------------------------\n");
  tmp=2.0*sqrt(fabs((sum2/(REAL)NR_BLOCKS)-SQR(sum)/(REAL)SQR(NR_BLOCKS)));
  fprintf(FilePtr,"\tAverage   %18.5lf [J/mol/K] +/- %18.5lf [J/mol/K]\n",
          (double)(sum/(REAL)NR_BLOCKS),(double)tmp);
  fprintf(FilePtr,"\tAverage   %18.5lf [cal/mol/K] +/- %18.5lf [cal/mol/K]\n",
          (double)(sum*J_TO_CAL/(REAL)NR_BLOCKS),(double)(tmp*J_TO_CAL));
*/

  if(NumberOfSystems==2)
  {
    // Heat of vaporization
    fprintf(FilePtr,"\n");
    fprintf(FilePtr,"Heat of vaporization:\n");
    fprintf(FilePtr,"=====================\n");
    sum=sum2=0.0;
    N=NumberOfAdsorbateMolecules[CurrentSystem]+NumberOfCationMolecules[CurrentSystem];
    if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
      N+=Framework[CurrentSystem].TotalNumberOfAtoms;
    for(i=0;i<NR_BLOCKS;i++)
    {
      if(BlockCount[CurrentSystem][i]>0.0)
      {
        tmp=ENERGY_TO_KELVIN*((fabs(EnergyPerMolecule[CurrentSystem][i]-EnergyPerMolecule[1-CurrentSystem][i])/BlockCount[CurrentSystem][i])
             +((MolecularStressTensorAverage[CurrentSystem][i].ax+MolecularStressTensorAverage[CurrentSystem][i].by+
              MolecularStressTensorAverage[CurrentSystem][i].cz)/(3.0*BlockCount[CurrentSystem][i]))*
              (fabs(VolumePerMolecule[CurrentSystem][i]-VolumePerMolecule[1-CurrentSystem][i])/BlockCount[CurrentSystem][i]));
        sum+=tmp;
        sum2+=SQR(tmp);

        fprintf(FilePtr,"\tBlock[%2d] %-18.5lf [K]\n",i,(double)tmp);
      }
      else
        fprintf(FilePtr,"\tBlock[%2d] %-18.5lf [K]\n",i,(double)0.0);
    }
    fprintf(FilePtr,"\t------------------------------------------------------------------------------\n");
    tmp=2.0*sqrt(fabs((sum2/(REAL)NR_BLOCKS)-SQR(sum)/(REAL)SQR(NR_BLOCKS)));
    fprintf(FilePtr,"\tAverage   %18.5lf +/- %18lf [K]\n",
      (double)(sum/(REAL)NR_BLOCKS),(double)tmp);
    fprintf(FilePtr,"\tAverage   %18.5lf +/- %18lf [kJ/mol]\n",
      (double)(sum*KELVIN_TO_KJ_PER_MOL/(REAL)NR_BLOCKS),(double)(tmp*KELVIN_TO_KJ_PER_MOL));
    fprintf(FilePtr,"\tAverage   %18.5lf +/- %18lf [kcal/mol]\n",
      (double)(sum*KELVIN_TO_KCAL_PER_MOL/(REAL)NR_BLOCKS),(double)(tmp*KELVIN_TO_KCAL_PER_MOL));


    fprintf(FilePtr,"\n");
    fprintf(FilePtr,"Heat of vaporization term (U_v/N_v-U_l/N_l):\n");
    fprintf(FilePtr,"============================================\n");
    sum=sum2=0.0;
    for(i=0;i<NR_BLOCKS;i++)
    {
      if(BlockCount[CurrentSystem][i]>0.0)
      {
        tmp=ENERGY_TO_KELVIN*fabs(EnergyPerMolecule[CurrentSystem][i]-EnergyPerMolecule[1-CurrentSystem][i])/BlockCount[CurrentSystem][i];
        sum+=tmp;
        sum2+=SQR(tmp);

        fprintf(FilePtr,"\tBlock[%2d] %-18.5lf [K]\n",i,(double)tmp);
      }
      else
        fprintf(FilePtr,"\tBlock[%2d] %-18.5lf [K]\n",i,(double)0.0);
    }
    fprintf(FilePtr,"\t------------------------------------------------------------------------------\n");
    tmp=2.0*sqrt(fabs((sum2/(REAL)NR_BLOCKS)-SQR(sum)/(REAL)SQR(NR_BLOCKS)));
    fprintf(FilePtr,"\tAverage   %18.5lf +/- %18lf [K]\n",
      (double)(sum/(REAL)NR_BLOCKS),(double)tmp);
    fprintf(FilePtr,"\tAverage   %18.5lf +/- %18lf [kJ/mol]\n",
      (double)(sum*KELVIN_TO_KJ_PER_MOL/(REAL)NR_BLOCKS),(double)(tmp*KELVIN_TO_KJ_PER_MOL));
    fprintf(FilePtr,"\tAverage   %18.5lf +/- %18lf [kcal/mol]\n",
      (double)(sum*KELVIN_TO_KCAL_PER_MOL/(REAL)NR_BLOCKS),(double)(tmp*KELVIN_TO_KCAL_PER_MOL));

    fprintf(FilePtr,"\n");
    fprintf(FilePtr,"Heat of vaporization term p (V_v-V_l), [is approximately RT]:\n");
    fprintf(FilePtr,"=============================================================\n");
    sum=sum2=0.0;
    N=NumberOfAdsorbateMolecules[CurrentSystem]+NumberOfCationMolecules[CurrentSystem];
    if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
      N+=Framework[CurrentSystem].TotalNumberOfAtoms;
    for(i=0;i<NR_BLOCKS;i++)
    {
      if(BlockCount[CurrentSystem][i]>0.0)
      {
        tmp=ENERGY_TO_KELVIN*fabs(((MolecularStressTensorAverage[CurrentSystem][i].ax+MolecularStressTensorAverage[CurrentSystem][i].by+
              MolecularStressTensorAverage[CurrentSystem][i].cz)/(3.0*BlockCount[CurrentSystem][i]))*
             (VolumePerMolecule[CurrentSystem][i]-VolumePerMolecule[1-CurrentSystem][i])/BlockCount[CurrentSystem][i]);
        sum+=tmp;
        sum2+=SQR(tmp);

        fprintf(FilePtr,"\tBlock[%2d] %-18.5lf [K]\n",i,(double)tmp);
      }
      else
        fprintf(FilePtr,"\tBlock[%2d] %-18.5lf [K]\n",i,(double)0.0);
    }
    fprintf(FilePtr,"\t------------------------------------------------------------------------------\n");
    tmp=2.0*sqrt(fabs((sum2/(REAL)NR_BLOCKS)-SQR(sum)/(REAL)SQR(NR_BLOCKS)));
    fprintf(FilePtr,"\tAverage   %18.5lf +/- %18lf [K]\n",
      (double)(sum/(REAL)NR_BLOCKS),(double)tmp);
    fprintf(FilePtr,"\tAverage   %18.5lf +/- %18lf [kJ/mol]\n",
      (double)(sum*KELVIN_TO_KJ_PER_MOL/(REAL)NR_BLOCKS),(double)(tmp*KELVIN_TO_KJ_PER_MOL));
    fprintf(FilePtr,"\tAverage   %18.5lf +/- %18lf [kcal/mol]\n",
      (double)(sum*KELVIN_TO_KCAL_PER_MOL/(REAL)NR_BLOCKS),(double)(tmp*KELVIN_TO_KCAL_PER_MOL));
  }

  // Heat of desorption
  fprintf(FilePtr,"\n");
  fprintf(FilePtr,"Heat of desorption:\n");
  fprintf(FilePtr,"===================\n");
  sum=sum2=0.0;
  for(i=0;i<NR_BLOCKS;i++)
  {
    if(BlockCount[CurrentSystem][i]>0.0)
    {
      HeatOfDesorption[i]=ENERGY_TO_KELVIN*(therm_baro_stats.ExternalTemperature[CurrentSystem]-
         (TotalEnergyTimesNumberOfMoleculesAverage[CurrentSystem][i]/BlockCount[CurrentSystem][i]-
         (UTotalAverage[CurrentSystem][i]/BlockCount[CurrentSystem][i])*
         (NumberOfMoleculesAverage[CurrentSystem][i]/BlockCount[CurrentSystem][i]))/
         (NumberOfMoleculesSquaredAverage[CurrentSystem][i]/BlockCount[CurrentSystem][i]-
         SQR(NumberOfMoleculesAverage[CurrentSystem][i]/BlockCount[CurrentSystem][i])));
      sum+=HeatOfDesorption[i];
      sum2+=SQR(HeatOfDesorption[i]);
      fprintf(FilePtr,"\tBlock[%2d] %-18.5lf [K]\n",i,(double)HeatOfDesorption[i]);
    }
    else
      fprintf(FilePtr,"\tBlock[%2d] %-18.5lf [K]\n",i,(double)0.0);
  }
  fprintf(FilePtr,"\t------------------------------------------------------------------------------\n");
  fprintf(FilePtr,"\tAverage   %18.5lf +/- %18lf [K]\n",
    (double)(sum/(REAL)NR_BLOCKS),
    (double)(2.0*sqrt(fabs((sum2/(REAL)NR_BLOCKS)-SQR(sum)/(REAL)SQR(NR_BLOCKS)))));
  fprintf(FilePtr,"\t          %18.5lf +/- %18lf [KJ/MOL]\n",
    (double)(sum/(REAL)NR_BLOCKS)*KELVIN_TO_KJ_PER_MOL,
    (double)(2.0*sqrt(fabs((sum2/(REAL)NR_BLOCKS)-SQR(sum)/(REAL)SQR(NR_BLOCKS)))*KELVIN_TO_KJ_PER_MOL));
  fprintf(FilePtr,"Note: Ug should be added to this value\n");
  fprintf(FilePtr,"Note: The heat of adsorption Q=-H\n");

  // Heat of desorption
  fprintf(FilePtr,"\n");
  fprintf(FilePtr,"Heat of desorption (Host-Adsorbate contribution):\n");
  fprintf(FilePtr,"=================================================\n");
  sum=sum2=0.0;
  for(i=0;i<NR_BLOCKS;i++)
  {
    if(BlockCount[CurrentSystem][i]>0.0)
    {
      HeatOfDesorption[i]=ENERGY_TO_KELVIN*(-
         (HostAdsorbateEnergyTimesNumberOfMoleculesAverage[CurrentSystem][i]/BlockCount[CurrentSystem][i]-
         (UHostAdsorbateAverage[CurrentSystem][i]/BlockCount[CurrentSystem][i])*
         (NumberOfMoleculesAverage[CurrentSystem][i]/BlockCount[CurrentSystem][i]))/
         (NumberOfMoleculesSquaredAverage[CurrentSystem][i]/BlockCount[CurrentSystem][i]-
         SQR(NumberOfMoleculesAverage[CurrentSystem][i]/BlockCount[CurrentSystem][i])));
      sum+=HeatOfDesorption[i];
      sum2+=SQR(HeatOfDesorption[i]);
      fprintf(FilePtr,"\tBlock[%2d] %-18.5lf [K]\n",i,(double)HeatOfDesorption[i]);
    }
    else
      fprintf(FilePtr,"\tBlock[%2d] %-18.5lf [K]\n",i,(double)0.0);
  }
  fprintf(FilePtr,"\t------------------------------------------------------------------------------\n");
  fprintf(FilePtr,"\tAverage   %18.5lf +/- %18lf [K]\n",
    (double)(sum/(REAL)NR_BLOCKS),
    (double)(2.0*sqrt(fabs((sum2/(REAL)NR_BLOCKS)-SQR(sum)/(REAL)SQR(NR_BLOCKS)))));
  fprintf(FilePtr,"\t          %18.5lf +/- %18lf [KJ/MOL]\n",
    (double)(sum/(REAL)NR_BLOCKS)*KELVIN_TO_KJ_PER_MOL,
    (double)(2.0*sqrt(fabs((sum2/(REAL)NR_BLOCKS)-SQR(sum)/(REAL)SQR(NR_BLOCKS)))*KELVIN_TO_KJ_PER_MOL));
  fprintf(FilePtr,"Note: the sum of the parts is *not* exactly equal to the true value (Heat of desorption)\n");
  fprintf(FilePtr,"Note: Ug and K_BT are not included here\n");
  fprintf(FilePtr,"Note: The heat of adsorption Q=-H\n");

  // Heat of desorption
  fprintf(FilePtr,"\n");
  fprintf(FilePtr,"Heat of desorption (Adsorbate-Adsorbate contribution):\n");
  fprintf(FilePtr,"======================================================\n");
  sum=sum2=0.0;
  for(i=0;i<NR_BLOCKS;i++)
  {
    if(BlockCount[CurrentSystem][i]>0.0)
    {
      HeatOfDesorption[i]=ENERGY_TO_KELVIN*(-
         (AdsorbateAdsorbateEnergyTimesNumberOfMoleculesAverage[CurrentSystem][i]/BlockCount[CurrentSystem][i]-
         (UAdsorbateAdsorbateAverage[CurrentSystem][i]/BlockCount[CurrentSystem][i])*
         (NumberOfMoleculesAverage[CurrentSystem][i]/BlockCount[CurrentSystem][i]))/
         (NumberOfMoleculesSquaredAverage[CurrentSystem][i]/BlockCount[CurrentSystem][i]-
         SQR(NumberOfMoleculesAverage[CurrentSystem][i]/BlockCount[CurrentSystem][i])));
      sum+=HeatOfDesorption[i];
      sum2+=SQR(HeatOfDesorption[i]);
      fprintf(FilePtr,"\tBlock[%2d] %-18.5lf [K]\n",i,(double)HeatOfDesorption[i]);
    }
    else
      fprintf(FilePtr,"\tBlock[%2d] %-18.5lf [K]\n",i,(double)0.0);
  }
  fprintf(FilePtr,"\t------------------------------------------------------------------------------\n");
  fprintf(FilePtr,"\tAverage   %18.5lf +/- %18lf [K]\n",
    (double)(sum/(REAL)NR_BLOCKS),
    (double)(2.0*sqrt(fabs((sum2/(REAL)NR_BLOCKS)-SQR(sum)/(REAL)SQR(NR_BLOCKS)))));
  fprintf(FilePtr,"\t          %18.5lf +/- %18lf [KJ/MOL]\n",
    (double)(sum/(REAL)NR_BLOCKS)*KELVIN_TO_KJ_PER_MOL,
    (double)(2.0*sqrt(fabs((sum2/(REAL)NR_BLOCKS)-SQR(sum)/(REAL)SQR(NR_BLOCKS)))*KELVIN_TO_KJ_PER_MOL));
  fprintf(FilePtr,"Note: the sum of the parts is *not* exactly equal to the true value (Heat of desorption)\n");
  fprintf(FilePtr,"Note: Ug and K_BT are not included here\n");
  fprintf(FilePtr,"Note: The heat of adsorption Q=-H\n");



  // Derivative
  fprintf(FilePtr,"\n");
  fprintf(FilePtr,"derivative of the chemical potential with respect to density (constant T,V):\n");
  fprintf(FilePtr,"============================================================================\n");
  sum=sum2=0.0;
  for(i=0;i<NR_BLOCKS;i++)
  {
    if(BlockCount[CurrentSystem][i]>0.0)
    {
      tmp=ENERGY_TO_KELVIN*(Volume[0]*therm_baro_stats.ExternalTemperature[CurrentSystem]/
          (NumberOfMoleculesSquaredAverage[CurrentSystem][i]/BlockCount[CurrentSystem][i]-
          SQR(NumberOfMoleculesAverage[CurrentSystem][i]/BlockCount[CurrentSystem][i])));
      sum+=tmp;
      sum2+=SQR(tmp);
      fprintf(FilePtr,"\tBlock[%2d] %-18.5lf [-]\n",i,(double)tmp);
    }
    else
      fprintf(FilePtr,"\tBlock[%2d] %-18.5lf [-]\n",i,(double)0.0);
  }
  fprintf(FilePtr,"\t------------------------------------------------------------------------------\n");
  fprintf(FilePtr,"\tAverage   %18.5lf +/- %18lf [-]\n",
    (double)(sum/(REAL)NR_BLOCKS),
    (double)(2.0*sqrt(fabs((sum2/(REAL)NR_BLOCKS)-SQR(sum)/(REAL)SQR(NR_BLOCKS)))));


  // output final averages of the principle moments of inertia per component
  if(ComputePrincipleMomentsOfInertia)
  {
    fprintf(FilePtr,"\n");
    fprintf(FilePtr,"Average Principle Moments of Inertia:\n");
    fprintf(FilePtr,"=====================================\n");
    for(j=0;j<NumberOfComponents;j++)
    {
      sumv1.x=sumv1.y=sumv1.z=0.0;
      sumv2.x=sumv2.y=sumv2.z=0.0;
      for(i=0;i<NR_BLOCKS;i++)
      {
        if(PrincipleMomentsOfInertiaCount[CurrentSystem][j][i]>0.0)
        {
          av.x=PrincipleMomentsOfInertiaAverage[CurrentSystem][j][i].x/PrincipleMomentsOfInertiaCount[CurrentSystem][j][i];
          av.y=PrincipleMomentsOfInertiaAverage[CurrentSystem][j][i].y/PrincipleMomentsOfInertiaCount[CurrentSystem][j][i];
          av.z=PrincipleMomentsOfInertiaAverage[CurrentSystem][j][i].z/PrincipleMomentsOfInertiaCount[CurrentSystem][j][i];
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
      fprintf(FilePtr,"\t------------------------------------------------------------------------------\n");
      fprintf(FilePtr,"\t[%s] Average Principle Moment of Inertia:  %lg +/- %lg [-], %lg +/- %lg [-], %lg +/- %lg [-]\n\n",
        Components[j].Name,
        (double)(sumv1.x/(REAL)NR_BLOCKS),
        (double)(2.0*sqrt(fabs((sumv2.x/(REAL)NR_BLOCKS)-SQR(sumv1.x)/(REAL)SQR(NR_BLOCKS)))),
        (double)(sumv1.y/(REAL)NR_BLOCKS),
        (double)(2.0*sqrt(fabs((sumv2.y/(REAL)NR_BLOCKS)-SQR(sumv1.y)/(REAL)SQR(NR_BLOCKS)))),
        (double)(sumv1.z/(REAL)NR_BLOCKS),
        (double)(2.0*sqrt(fabs((sumv2.z/(REAL)NR_BLOCKS)-SQR(sumv1.z)/(REAL)SQR(NR_BLOCKS)))));
    }
  }


  fprintf(FilePtr,"\n\n\n\n\n");
  fprintf(FilePtr,"Average energies of the system[%d]:\n",CurrentSystem);
  fprintf(FilePtr,"========================================================================\n");

  // Host intra molecular energies
  PrintProperty(FilePtr,"Average Host Bond stretch energy:\n=================================\n",
                "[K]",ENERGY_TO_KELVIN,UHostBondAverage);

  PrintProperty(FilePtr,"Average Host UreyBradley stretch energy:\n==============================================\n",
                "[K]",ENERGY_TO_KELVIN,UHostUreyBradleyAverage);

  PrintProperty(FilePtr,"Average Host Bend angle energy:\n===============================\n",
                "[K]",ENERGY_TO_KELVIN,UHostBendAverage);

  PrintProperty(FilePtr,"Average Host Bend angle inversion energy:\n=========================================\n",
                "[K]",ENERGY_TO_KELVIN,UHostInversionBendAverage);

  PrintProperty(FilePtr,"Average Host Torsion energy:\n============================\n",
                "[K]",ENERGY_TO_KELVIN,UHostTorsionAverage);

  PrintProperty(FilePtr,"Average Host Improper Torsion energy:\n=====================================\n",
                "[K]",ENERGY_TO_KELVIN,UHostImproperTorsionAverage);

  PrintProperty(FilePtr,"Average Host Bond-Bond cross term energy:\n===============================================\n",
                "[K]",ENERGY_TO_KELVIN,UHostBondBondAverage);

  PrintProperty(FilePtr,"Average Host Bend-Bend cross term energy:\n=========================================\n",
                "[K]",ENERGY_TO_KELVIN,UHostBendBendAverage);

  PrintProperty(FilePtr,"Average Host Bond-Bend cross term energy:\n============================================\n",
                "[K]",ENERGY_TO_KELVIN,UHostBondBendAverage);

  PrintProperty(FilePtr,"Average Host Bond-Torsion cross term energy:\n============================================\n",
                "[K]",ENERGY_TO_KELVIN,UHostBondTorsionAverage);

  PrintProperty(FilePtr,"Average Host Bend-Torsion cross term energy:\n============================================\n",
                "[K]",ENERGY_TO_KELVIN,UHostBendTorsionAverage);

  // Adsorbate intra molecular energies
  PrintProperty(FilePtr,"Average Adsorbate Bond stretch energy:\n====================================\n",
                "[K]",ENERGY_TO_KELVIN,UAdsorbateBondAverage);

  PrintProperty(FilePtr,"Average Adsorbate UreyBradley stretch energy:\n==============================================\n",
                "[K]",ENERGY_TO_KELVIN,UAdsorbateUreyBradleyAverage);

  PrintProperty(FilePtr,"Average Adsorbate Bend angle energy:\n====================================\n",
                "[K]",ENERGY_TO_KELVIN,UAdsorbateBendAverage);

  PrintProperty(FilePtr,"Average Adsorbate Bend angle inversion energy:\n==============================================\n",
                "[K]",ENERGY_TO_KELVIN,UAdsorbateInversionBendAverage);

  PrintProperty(FilePtr,"Average Adsorbate Torsion energy:\n=================================\n",
                "[K]",ENERGY_TO_KELVIN,UAdsorbateTorsionAverage);

  PrintProperty(FilePtr,"Average Adsorbate Improper Torsion energy:\n==========================================\n",
                "[K]",ENERGY_TO_KELVIN,UAdsorbateImproperTorsionAverage);

  PrintProperty(FilePtr,"Average Adsorbate Bond-Bond cross term energy:\n====================================================\n",
                "[K]",ENERGY_TO_KELVIN,UAdsorbateBondBondAverage);

  PrintProperty(FilePtr,"Average Adsorbate Bend-Bend cross term energy:\n==============================================\n",
                "[K]",ENERGY_TO_KELVIN,UAdsorbateBendBendAverage);

  PrintProperty(FilePtr,"Average Adsorbate Bond-Bend cross term energy:\n=================================================\n",
                "[K]",ENERGY_TO_KELVIN,UAdsorbateBondBendAverage);

  PrintProperty(FilePtr,"Average Adsorbate Bond-Torsion cross term energy:\n=================================================\n",
                "[K]",ENERGY_TO_KELVIN,UAdsorbateBondTorsionAverage);

  PrintProperty(FilePtr,"Average Adsorbate Bend-Torsion cross term energy:\n=================================================\n",
                "[K]",ENERGY_TO_KELVIN,UAdsorbateBendTorsionAverage);

  PrintProperty(FilePtr,"Average Adsorbate Intra Van der Waals energy:\n=============================================\n",
                "[K]",ENERGY_TO_KELVIN,UAdsorbateIntraVDWAverage);

  PrintProperty(FilePtr,"Average Adsorbate Intra charge-charge Coulomb energy:\n=======================================\n",
                "[K]",ENERGY_TO_KELVIN,UAdsorbateIntraChargeChargeAverage);

  PrintProperty(FilePtr,"Average Adsorbate Intra charge-bonddipole Coulomb energy:\n=======================================\n",
                "[K]",ENERGY_TO_KELVIN,UAdsorbateIntraChargeBondDipoleAverage);

  PrintProperty(FilePtr,"Average Adsorbate Intra bonddipole-bonddipole Coulomb energy:\n=======================================\n",
                "[K]",ENERGY_TO_KELVIN,UAdsorbateIntraBondDipoleBondDipoleAverage);


  // Cation intra molecular energies
  PrintProperty(FilePtr,"Average Cation Bond stretch energy:\n=================================\n",
                "[K]",ENERGY_TO_KELVIN,UCationBondAverage);

  PrintProperty(FilePtr,"Average Cation UreyBradley stretch energy:\n==========================================\n",
                "[K]",ENERGY_TO_KELVIN,UCationUreyBradleyAverage);

  PrintProperty(FilePtr,"Average Cation Bend angle energy:\n=================================\n",
                "[K]",ENERGY_TO_KELVIN,UCationBendAverage);

  PrintProperty(FilePtr,"Average Cation Bend angle inversion energy:\n===========================================\n",
                "[K]",ENERGY_TO_KELVIN,UCationInversionBendAverage);

  PrintProperty(FilePtr,"Average Cation Torsion energy:\n==============================\n",
                "[K]",ENERGY_TO_KELVIN,UCationTorsionAverage);

  PrintProperty(FilePtr,"Average Cation Improper Torsion energy:\n=======================================\n",
                "[K]",ENERGY_TO_KELVIN,UCationImproperTorsionAverage);

  PrintProperty(FilePtr,"Average Cation Bond-Bond cross term energy:\n=================================================\n",
                "[K]",ENERGY_TO_KELVIN,UCationBondBondAverage);

  PrintProperty(FilePtr,"Average Cation Bend-Bend cross term energy:\n===========================================\n",
                "[K]",ENERGY_TO_KELVIN,UCationBendBendAverage);

  PrintProperty(FilePtr,"Average Cation Bond-Bend cross term energy:\n==============================================\n",
                "[K]",ENERGY_TO_KELVIN,UCationBondBendAverage);

  PrintProperty(FilePtr,"Average Cation Bond-Torsion cross term energy:\n==============================================\n",
                "[K]",ENERGY_TO_KELVIN,UCationBondTorsionAverage);

  PrintProperty(FilePtr,"Average Cation Bend-Torsion cross term energy:\n==============================================\n",
                "[K]",ENERGY_TO_KELVIN,UCationBendTorsionAverage);

  PrintProperty(FilePtr,"Average Cation Intra Van der Waals energy:\n==========================================\n",
                "[K]",ENERGY_TO_KELVIN,UCationIntraVDWAverage);

  PrintProperty(FilePtr,"Average Cation Intra charge-charge Coulomb energy:\n====================================\n",
                "[K]",ENERGY_TO_KELVIN,UCationIntraChargeChargeAverage);

  PrintProperty(FilePtr,"Average Cation Intra charge-bonddipole Coulomb energy:\n====================================\n",
                "[K]",ENERGY_TO_KELVIN,UCationIntraChargeBondDipoleAverage);

  PrintProperty(FilePtr,"Average Cation Intra bonddipole-bonddipole Coulomb energy:\n====================================\n",
                "[K]",ENERGY_TO_KELVIN,UCationIntraBondDipoleBondDipoleAverage);

  // Host-Host inter molecular energies
  PrintEnergies(FilePtr,"Average Host-Host energy:\n=========================\n","[K]",ENERGY_TO_KELVIN,
                UHostHostAverage,UHostHostVDWAverage,UHostHostCoulombAverage);

  // Adsorbate-Adsorbate inter molecular energies
  PrintEnergies(FilePtr,"Average Adsorbate-Adsorbate energy:\n===================================\n","[K]",ENERGY_TO_KELVIN,
                UAdsorbateAdsorbateAverage,UAdsorbateAdsorbateVDWAverage,UAdsorbateAdsorbateCoulombAverage);

  // Cation-Cation inter molecular energies
  PrintEnergies(FilePtr,"Average Cation-Cation energy:\n=============================\n","[K]",ENERGY_TO_KELVIN,
                UCationCationAverage,UCationCationVDWAverage,UCationCationCoulombAverage);

  // Host-Adsorbate inter molecular energies
  PrintEnergies(FilePtr,"Average Host-Adsorbate energy:\n==============================\n","[K]",ENERGY_TO_KELVIN,
                UHostAdsorbateAverage,UHostAdsorbateVDWAverage,UHostAdsorbateCoulombAverage);

  // Host-Cation inter molecular energies
  PrintEnergies(FilePtr,"Average Host-Cation energy:\n===========================\n","[K]",ENERGY_TO_KELVIN,
                UHostCationAverage,UHostCationVDWAverage,UHostCationCoulombAverage);

  // Adsorbate-Cation inter molecular energies
  PrintEnergies(FilePtr,"Average Adsorbate-Cation energy:\n================================\n","[K]",ENERGY_TO_KELVIN,
                UAdsorbateCationAverage,UAdsorbateCationVDWAverage,UAdsorbateCationCoulombAverage);

  // Host polarization energy
  PrintProperty(FilePtr,"Host polarization energy:\n=======================\n",
                "[K]",ENERGY_TO_KELVIN,UHostPolarizationAverage);

  // Adsorbate polarization energy
  PrintProperty(FilePtr,"Adsorbate polarization energy:\n=======================\n",
                "[K]",ENERGY_TO_KELVIN,UAdsorbatePolarizationAverage);

  // Cation polarization energy
  PrintProperty(FilePtr,"Cation polarization energy:\n=======================\n",
                "[K]",ENERGY_TO_KELVIN,UCationPolarizationAverage);

  // Host back-polarization energy
  PrintProperty(FilePtr,"Host back-polarization energy:\n=======================\n",
                "[K]",ENERGY_TO_KELVIN,UHostBackPolarizationAverage);

  // Adsorbate back-polarization energy
  PrintProperty(FilePtr,"Adsorbate back-polarization energy:\n=======================\n",
                "[K]",ENERGY_TO_KELVIN,UAdsorbateBackPolarizationAverage);

  // Cation back-polarization energy
  PrintProperty(FilePtr,"Cation back-polarization energy:\n=======================\n",
                "[K]",ENERGY_TO_KELVIN,UCationBackPolarizationAverage);


  // Tailcorrection energy
  PrintProperty(FilePtr,"Tail-correction energy:\n=======================\n",
                "[K]",ENERGY_TO_KELVIN,UTailCorrectionAverage);

  // Distance constraints energy
  PrintProperty(FilePtr,"Distance-constraints energy:\n=======================\n",
                "[K]",ENERGY_TO_KELVIN,UDistanceConstraintsAverage);

  // Angle constraints energy
  PrintProperty(FilePtr,"Angle-constraints energy:\n=======================\n",
                "[K]",ENERGY_TO_KELVIN,UAngleConstraintsAverage);

  // Dihedral constraints energy
  PrintProperty(FilePtr,"Dihedral-constraints energy:\n=======================\n",
                "[K]",ENERGY_TO_KELVIN,UDihedralConstraintsAverage);

  // Inversion-bend constraints energy
  PrintProperty(FilePtr,"Inversion-bend constraints energy:\n=======================\n",
                "[K]",ENERGY_TO_KELVIN,UInversionBendConstraintsAverage);

  // Out-of-plane-distance constraints energy
  PrintProperty(FilePtr,"Out-of-plane-distance constraints energy:\n=======================\n",
                "[K]",ENERGY_TO_KELVIN,UOutOfPlaneDistanceConstraintsAverage);

  // Exclusion constraints energy
  PrintProperty(FilePtr,"Exclusion-constraints energy:\n=======================\n",
                "[K]",ENERGY_TO_KELVIN,UExclusionConstraintsAverage);

  // Total energy
  PrintProperty(FilePtr,"Total energy:\n=============\n",
                "[K]",ENERGY_TO_KELVIN,UTotalAverage);

  // Number of molecules
  fprintf(FilePtr,"\n");
  fprintf(FilePtr,"Number of molecules:\n");
  fprintf(FilePtr,"====================\n\n");

  for(j=0;j<NumberOfComponents;j++)
  {
    fprintf(FilePtr,"Component %d [%s]\n",j,Components[j].Name);
    fprintf(FilePtr,"-------------------------------------------------------------\n");

    // absolute adsorption
    sum=sum_vdw=sum_coul=0.0;
    sum2=sum_vdw2=sum_coul2=0.0;
    for(i=0;i<NR_BLOCKS;i++)
    {
      if(BlockCount[CurrentSystem][i]>0.0)
      {
        sum+=NumberOfMoleculesPerComponentAverage[CurrentSystem][j][i]/BlockCount[CurrentSystem][i];
        sum2+=SQR(NumberOfMoleculesPerComponentAverage[CurrentSystem][j][i]/BlockCount[CurrentSystem][i]);

        fprintf(FilePtr,"\tBlock[%2d] %-18.5lf [-]\n",i,
                (double)(NumberOfMoleculesPerComponentAverage[CurrentSystem][j][i]/BlockCount[CurrentSystem][i]));
      }
      else
        fprintf(FilePtr,"\tBlock[%2d] %-18.5lf [-]\n",i,(double)0.0);
    }
    nr_molecules=sum/(REAL)NR_BLOCKS;
    fprintf(FilePtr,"\t------------------------------------------------------------------------------\n");
    nr=NumberOfUnitCells[0].x*NumberOfUnitCells[0].y*NumberOfUnitCells[0].z;
    fprintf(FilePtr,"\tAverage                                %18.10lf +/- %18.10lf [-]\n",
      (double)(sum/(REAL)NR_BLOCKS),(double)(2.0*sqrt(fabs((sum2/(REAL)NR_BLOCKS)-SQR(sum)/(REAL)SQR(NR_BLOCKS)))));
    fprintf(FilePtr,"\tAverage loading absolute [molecules/unit cell]  %18.10lf +/- %18.10lf [-]\n",
      (double)(sum/((REAL)nr*(REAL)NR_BLOCKS)),
      (double)(2.0*sqrt(fabs((sum2/(REAL)NR_BLOCKS)-SQR(sum)/(REAL)SQR(NR_BLOCKS)))/(REAL)nr));

    fprintf(FilePtr,"\tAverage loading absolute [mol/kg framework]    %18.10lf +/- %18.10lf [-]\n",
      (double)(sum*Components[j].MOLEC_PER_UC_TO_MOL_PER_KG[CurrentSystem]/((REAL)nr*(REAL)NR_BLOCKS)),
      (double)(2.0*Components[j].MOLEC_PER_UC_TO_MOL_PER_KG[CurrentSystem]*sqrt(fabs((sum2/(REAL)NR_BLOCKS)-SQR(sum)/(REAL)SQR(NR_BLOCKS)))/(REAL)nr));

    fprintf(FilePtr,"\tAverage loading absolute [milligram/gram framework]    %18.10lf +/- %18.10lf [-]\n",
      (double)(sum*Components[j].MOLEC_PER_UC_TO_GRAM_PER_GRAM_OF_FRAMEWORK[CurrentSystem]/((REAL)nr*(REAL)NR_BLOCKS)),
      (double)(2.0*Components[j].MOLEC_PER_UC_TO_GRAM_PER_GRAM_OF_FRAMEWORK[CurrentSystem]*sqrt(fabs((sum2/(REAL)NR_BLOCKS)-SQR(sum)/(REAL)SQR(NR_BLOCKS)))/(REAL)nr));

    fprintf(FilePtr,"\tAverage loading absolute [cm^3 (STP)/gr framework]    %18.10lf +/- %18.10lf [-]\n",
      (double)(sum*Components[j].MOLEC_PER_UC_TO_CC_STP_G[CurrentSystem]/((REAL)nr*(REAL)NR_BLOCKS)),
      (double)(2.0*Components[j].MOLEC_PER_UC_TO_CC_STP_G[CurrentSystem]*sqrt(fabs((sum2/(REAL)NR_BLOCKS)-SQR(sum)/(REAL)SQR(NR_BLOCKS)))/(REAL)nr));

    fprintf(FilePtr,"\tAverage loading absolute [cm^3 (STP)/cm^3 framework]    %18.10lf +/- %18.10lf [-]\n",
      (double)(sum*Components[j].MOLEC_PER_UC_TO_CC_STP_CC[CurrentSystem]/((REAL)nr*(REAL)NR_BLOCKS)),
      (double)(2.0*Components[j].MOLEC_PER_UC_TO_CC_STP_CC[CurrentSystem]*sqrt(fabs((sum2/(REAL)NR_BLOCKS)-SQR(sum)/(REAL)SQR(NR_BLOCKS)))/(REAL)nr));
    fprintf(FilePtr,"\n");

    // excess adsorption
    sum=sum_vdw=sum_coul=0.0;
    sum2=sum_vdw2=sum_coul2=0.0;
    for(i=0;i<NR_BLOCKS;i++)
    {
      if(BlockCount[CurrentSystem][i]>0.0)
      {
        sum+=(NumberOfExcessMoleculesPerComponentAverage[CurrentSystem][j][i]/BlockCount[CurrentSystem][i]);
        sum2+=SQR(NumberOfExcessMoleculesPerComponentAverage[CurrentSystem][j][i]/BlockCount[CurrentSystem][i]);

        fprintf(FilePtr,"\tBlock[%2d] %-18.5lf [-]\n",i,
                (double)(NumberOfMoleculesPerComponentAverage[CurrentSystem][j][i]/BlockCount[CurrentSystem][i]));
      }
      else
        fprintf(FilePtr,"\tBlock[%2d] %-18.5lf [-]\n",i,(double)0.0);
    }
    nr_molecules=sum/(REAL)NR_BLOCKS;
    fprintf(FilePtr,"\t------------------------------------------------------------------------------\n");
    fprintf(FilePtr,"\tAverage                                %18.10lf +/- %18.10lf [-]\n",
      (double)(sum/(REAL)NR_BLOCKS),(double)(2.0*sqrt(fabs((sum2/(REAL)NR_BLOCKS)-SQR(sum)/(REAL)SQR(NR_BLOCKS)))));
    nr=NumberOfUnitCells[CurrentSystem].x*NumberOfUnitCells[CurrentSystem].y*NumberOfUnitCells[CurrentSystem].z;

    fprintf(FilePtr,"\tAverage loading excess [molecules/unit cell]  %18.10lf +/- %18.10lf [-]\n",
      (double)(sum/((REAL)nr*(REAL)NR_BLOCKS)),
      (double)(2.0*sqrt(fabs((sum2/(REAL)NR_BLOCKS)-SQR(sum)/(REAL)SQR(NR_BLOCKS)))/(REAL)nr));

    fprintf(FilePtr,"\tAverage loading excess [mol/kg framework]    %18.10lf +/- %18.10lf [-]\n",
      (double)(sum*Components[j].MOLEC_PER_UC_TO_MOL_PER_KG[CurrentSystem]/((REAL)nr*(REAL)NR_BLOCKS)),
      (double)(2.0*Components[j].MOLEC_PER_UC_TO_MOL_PER_KG[CurrentSystem]*sqrt(fabs((sum2/(REAL)NR_BLOCKS)-SQR(sum)/(REAL)SQR(NR_BLOCKS)))/(REAL)nr));

    fprintf(FilePtr,"\tAverage loading excess [milligram/gram framework]    %18.10lf +/- %18.10lf [-]\n",
      (double)(sum*Components[j].MOLEC_PER_UC_TO_GRAM_PER_GRAM_OF_FRAMEWORK[CurrentSystem]/((REAL)nr*(REAL)NR_BLOCKS)),
      (double)(2.0*Components[j].MOLEC_PER_UC_TO_GRAM_PER_GRAM_OF_FRAMEWORK[CurrentSystem]*sqrt(fabs((sum2/(REAL)NR_BLOCKS)-SQR(sum)/(REAL)SQR(NR_BLOCKS)))/(REAL)nr));

    fprintf(FilePtr,"\tAverage loading excess [cm^3 (STP)/gr framework]    %18.10lf +/- %18.10lf [-]\n",
      (double)(sum*Components[j].MOLEC_PER_UC_TO_CC_STP_G[CurrentSystem]/((REAL)nr*(REAL)NR_BLOCKS)),
      (double)(2.0*Components[j].MOLEC_PER_UC_TO_CC_STP_G[CurrentSystem]*sqrt(fabs((sum2/(REAL)NR_BLOCKS)-SQR(sum)/(REAL)SQR(NR_BLOCKS)))/(REAL)nr));

    fprintf(FilePtr,"\tAverage loading excess [cm^3 (STP)/cm^3 framework]    %18.10lf +/- %18.10lf [-]\n",
      (double)(sum*Components[j].MOLEC_PER_UC_TO_CC_STP_CC[CurrentSystem]/((REAL)nr*(REAL)NR_BLOCKS)),
      (double)(2.0*Components[j].MOLEC_PER_UC_TO_CC_STP_CC[CurrentSystem]*sqrt(fabs((sum2/(REAL)NR_BLOCKS)-SQR(sum)/(REAL)SQR(NR_BLOCKS)))/(REAL)nr));
    fprintf(FilePtr,"\n");
  }

  // Average Widom Rosenbluth factor
  fprintf(FilePtr,"\n");
  fprintf(FilePtr,"Average Widom Rosenbluth factor:\n");
  fprintf(FilePtr,"================================\n");
  for(j=0;j<NumberOfComponents;j++)
  {
    sum=sum2=0.0;
    for(i=0;i<NR_BLOCKS;i++)
    {
      if(WidomRosenbluthFactorCount[CurrentSystem][j][i]>0.0)
      {
        tmp=WidomRosenbluthFactorAverage[CurrentSystem][j][i]/
            WidomRosenbluthFactorCount[CurrentSystem][j][i];
        sum+=tmp;
        sum2+=SQR(tmp);
        fprintf(FilePtr,"\tBlock[%2d] %-lg [-]\n",i,(double)tmp);
      }
      else
        fprintf(FilePtr,"\tBlock[%2d] %-lg [-]\n",i,(double)0.0);
    }
    fprintf(FilePtr,"\t------------------------------------------------------------------------------\n");
    fprintf(FilePtr,"\t[%s] Average Widom:   %lg +/- %lf [-]\n",
      Components[j].Name,
      (double)(sum/(REAL)NR_BLOCKS),
      (double)(2.0*sqrt(fabs((sum2/(REAL)NR_BLOCKS)-SQR(sum)/(REAL)SQR(NR_BLOCKS)))));
  }

  // Average Henry Coefficient
  fprintf(FilePtr,"\n");
  fprintf(FilePtr,"Average Henry coefficient:\n");
  fprintf(FilePtr,"==========================\n");
  FrameworkDensity=1e-3*Framework[CurrentSystem].FrameworkMass/(Volume[CurrentSystem]*CUBE(ANGSTROM)*AVOGADRO_CONSTANT);
  Temperature=therm_baro_stats.ExternalTemperature[CurrentSystem];
  for(j=0;j<NumberOfComponents;j++)
  {
    sum=sum2=0.0;
    for(i=0;i<NR_BLOCKS;i++)
    {
      if(WidomRosenbluthFactorCount[CurrentSystem][j][i]>0.0)
      {
        tmp=(1.0/(MOLAR_GAS_CONSTANT*Temperature*FrameworkDensity))*
            WidomRosenbluthFactorAverage[CurrentSystem][j][i]/
            (Components[j].IdealGasRosenbluthWeight[CurrentSystem]*WidomRosenbluthFactorCount[CurrentSystem][j][i]);
        sum+=tmp;
        sum2+=SQR(tmp);
        fprintf(FilePtr,"\tBlock[%2d] %-lg [mol/kg/Pa]\n",i,(double)tmp);
      }
      else
        fprintf(FilePtr,"\tBlock[%2d] %-lg [mol/kg/Pa]\n",i,(double)0.0);
    }
    fprintf(FilePtr,"\t------------------------------------------------------------------------------\n");
    fprintf(FilePtr,"\t[%s] Average Henry coefficient:  %lg +/- %lg [mol/kg/Pa]\n",
      Components[j].Name,
      (double)(sum/(REAL)NR_BLOCKS),
      (double)(2.0*sqrt(fabs((sum2/(REAL)NR_BLOCKS)-SQR(sum)/(REAL)SQR(NR_BLOCKS)))));
  }

  // Average adsorption energy
  fprintf(FilePtr,"\n");
  fprintf(FilePtr,"Average adsorption energy <U_gh>_1-<U_h>_0 obtained from Widom-insertion:\n");
  fprintf(FilePtr,"(Note: the total heat of adsorption is dH=<U_gh>_1-<U_h>_0 - <U_g> - RT)\n");
  fprintf(FilePtr,"=========================================================================\n");
  Temperature=therm_baro_stats.ExternalTemperature[CurrentSystem];
  for(j=0;j<NumberOfComponents;j++)
  {
    if(Components[j].Widom)
    {
      sum=sum2=0.0;
      for(i=0;i<NR_BLOCKS;i++)
      {
        if(WidomRosenbluthFactorCount[CurrentSystem][j][i]>0.0)
        {
          tmp=(WidomEnergyDifferenceAverage[CurrentSystem][j][i]/WidomRosenbluthFactorAverage[CurrentSystem][j][i]-
              WidomEnergyFrameworkAverage[CurrentSystem][j][i]/WidomEnergyFrameworkCount[CurrentSystem][j][i])*ENERGY_TO_KELVIN;
          sum+=tmp;
          sum2+=SQR(tmp);
          fprintf(FilePtr,"\tBlock[%2d] %-18.10lf [K]\n",i,(double)tmp);
        }
        else
          fprintf(FilePtr,"\tBlock[%2d] %-18.10lf [K]\n",i,(double)0.0);
      }
      fprintf(FilePtr,"\t------------------------------------------------------------------------------\n");
      fprintf(FilePtr,"\t[%s] Average  <U_gh>_1-<U_h>_0:  %18.10lf +/- %18.10lf [K]       (%18.10lf +/- %18.10lf kJ/mol)\n",
        Components[j].Name,
        (double)(sum/(REAL)NR_BLOCKS),
        (double)(2.0*sqrt(fabs((sum2/(REAL)NR_BLOCKS)-SQR(sum)/(REAL)SQR(NR_BLOCKS)))),
        (double)(sum*KELVIN_TO_KJ_PER_MOL/(REAL)NR_BLOCKS),
        (double)(2.0*KELVIN_TO_KJ_PER_MOL*sqrt(fabs((sum2/(REAL)NR_BLOCKS)-SQR(sum)/(REAL)SQR(NR_BLOCKS)))));
    }
  }
}

void WriteRestartStatistics(FILE *FilePtr)
{
  int i,j;
  int NumberOfBlocks=NR_BLOCKS;
  REAL Check;

  NumberOfBlocks=NR_BLOCKS;

  fwrite(&NumberOfBlocks,sizeof(int),1,FilePtr);
  fwrite(&Block,sizeof(int),1,FilePtr);
  fwrite(BlockCycle,sizeof(long long),NumberOfBlocks,FilePtr);

  for(i=0;i<NumberOfSystems;i++)
  {
    fwrite(BlockCount[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    fwrite(UHostHostAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UAdsorbateAdsorbateAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UCationCationAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UHostAdsorbateAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UHostCationAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UAdsorbateCationAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    fwrite(UHostHostVDWAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UAdsorbateAdsorbateVDWAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UCationCationVDWAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UHostAdsorbateVDWAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UHostCationVDWAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UAdsorbateCationVDWAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    fwrite(UHostHostCoulombAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UAdsorbateAdsorbateCoulombAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UCationCationCoulombAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UHostAdsorbateCoulombAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UHostCationCoulombAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UAdsorbateCationCoulombAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    fwrite(UHostBondAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UHostUreyBradleyAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UHostBendAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UHostInversionBendAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UHostTorsionAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UHostImproperTorsionAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UHostBondBondAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UHostBendBendAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UHostBondBendAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UHostBondTorsionAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UHostBendTorsionAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    fwrite(UAdsorbateBondAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UAdsorbateUreyBradleyAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UAdsorbateBondBondAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UAdsorbateBendAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UAdsorbateInversionBendAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UAdsorbateTorsionAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UAdsorbateImproperTorsionAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UAdsorbateBondBondAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UAdsorbateBendBendAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UAdsorbateBondBendAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UAdsorbateBondTorsionAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UAdsorbateBendTorsionAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UAdsorbateIntraVDWAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UAdsorbateIntraChargeChargeAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UAdsorbateIntraChargeBondDipoleAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UAdsorbateIntraBondDipoleBondDipoleAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    fwrite(UCationBondAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UCationUreyBradleyAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UCationBendAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UCationInversionBendAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UCationTorsionAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UCationImproperTorsionAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UCationBondBondAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UCationBendBendAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UCationBondBendAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UCationBondTorsionAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UCationBendTorsionAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UCationIntraVDWAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UCationIntraChargeChargeAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UCationIntraChargeBondDipoleAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UCationIntraBondDipoleBondDipoleAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    fwrite(UHostPolarizationAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UAdsorbatePolarizationAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UCationPolarizationAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    fwrite(UHostBackPolarizationAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UAdsorbateBackPolarizationAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UCationBackPolarizationAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    fwrite(UTailCorrectionAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UDistanceConstraintsAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UAngleConstraintsAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UDihedralConstraintsAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UInversionBendConstraintsAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UOutOfPlaneDistanceConstraintsAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UExclusionConstraintsAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UTotalAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    fwrite(TotalSystemDipoleAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(TotalSystemDipoleSquaredAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(TotalSystemNormDipoleAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(TotalSystemNormDipoleSquaredAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    fwrite(NumberOfMoleculesAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(NumberOfMoleculesSquaredAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(DensityAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    fwrite(BoxAverage[i],sizeof(VECTOR),NumberOfBlocks,FilePtr);
    fwrite(BoxAverageAX[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(BoxAverageAY[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(BoxAverageAZ[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(BoxAverageBX[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(BoxAverageBY[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(BoxAverageBZ[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(BoxAverageCX[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(BoxAverageCY[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(BoxAverageCZ[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(BoxLengthAverage[i],sizeof(VECTOR),NumberOfBlocks,FilePtr);
    fwrite(AlphaAngleAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(BetaAngleAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(GammaAngleAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(VolumeAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(VolumeSquaredAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    fwrite(TotalEnergyAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(TotalEnergySquaredAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(EnthalpyAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(EnthalpySquaredAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(EnthalpyTimesVolumeAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(EnthalpyTimesEnergyAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    fwrite(TemperatureAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(TemperatureCellAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(TemperatureTranslationAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(TemperatureRotationAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(TemperatureRotationAdsorbateAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(TemperatureTranslationAdsorbateAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    fwrite(TemperatureAdsorbatesAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(TemperatureCationsAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(TemperatureFrameworkAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    fwrite(MolecularPressureAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(MolecularStressTensorAverage[i],sizeof(REAL_MATRIX3x3),NumberOfBlocks,FilePtr);
    fwrite(PressureIdealGasPartAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(PressureExcessPartAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(PressureTailCorrectionAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(PressureAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(UNoseHooverAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(HeatOfVaporization[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(EnergyPerMolecule[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(VolumePerMolecule[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(CompressibilityAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    fwrite(BornTermAverage[i],sizeof(REAL_MATRIX9x9),NumberOfBlocks,FilePtr);
    fwrite(ConfigurationalStressFluctuationTermAverage[i],sizeof(REAL_MATRIX9x9),NumberOfBlocks,FilePtr);
    fwrite(ConfigurationalStressTensorAverage[i],sizeof(REAL_MATRIX3x3),NumberOfBlocks,FilePtr);
    fwrite(StressTensorAverage[i],sizeof(REAL_MATRIX3x3),NumberOfBlocks,FilePtr);

    fwrite(SurfaceAreaFrameworkAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    for(j=0;j<Framework[i].NumberOfFrameworks;j++)
      fwrite(SurfaceAreaFrameworksAverage[i][j],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(SurfaceAreaCationsAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(SurfaceAreaCount[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    fwrite(TotalEnergyTimesNumberOfMoleculesAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(HostAdsorbateEnergyTimesNumberOfMoleculesAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fwrite(AdsorbateAdsorbateEnergyTimesNumberOfMoleculesAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    for(j=0;j<NumberOfComponents;j++)
    {
      fwrite(NumberOfMoleculesPerComponentAverage[i][j],sizeof(REAL),NumberOfBlocks,FilePtr);
      fwrite(NumberOfExcessMoleculesPerComponentAverage[i][j],sizeof(REAL),NumberOfBlocks,FilePtr);
      fwrite(DensityPerComponentAverage[i][j],sizeof(REAL),NumberOfBlocks,FilePtr);

      fwrite(WidomRosenbluthFactorAverage[i][j],sizeof(REAL),NumberOfBlocks,FilePtr);
      fwrite(WidomRosenbluthFactorCount[i][j],sizeof(REAL),NumberOfBlocks,FilePtr);

      fwrite(WidomEnergyDifferenceAverage[i][j],sizeof(REAL),NumberOfBlocks,FilePtr);

      fwrite(WidomEnergyFrameworkAverage[i][j],sizeof(REAL),NumberOfBlocks,FilePtr);
      fwrite(WidomEnergyFrameworkCount[i][j],sizeof(REAL),NumberOfBlocks,FilePtr);

      fwrite(PrincipleMomentsOfInertiaAverage[i][j],sizeof(REAL),NumberOfBlocks,FilePtr);
      fwrite(PrincipleMomentsOfInertiaCount[i][j],sizeof(REAL),NumberOfBlocks,FilePtr);
    }
  }

  Check=123456789.0;
  fwrite(&Check,1,sizeof(REAL),FilePtr);
}

void AllocateStatisticsMemory(void)
{
  int i,j;
  int NumberOfBlocks;

  NumberOfBlocks=NR_BLOCKS;

  BlockCycle=(long long*)calloc(NumberOfBlocks,sizeof(long long));
  BlockCount=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));

  UHostHostAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UAdsorbateAdsorbateAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UCationCationAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UHostAdsorbateAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UHostCationAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UAdsorbateCationAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));

  UHostHostVDWAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UAdsorbateAdsorbateVDWAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UCationCationVDWAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UHostAdsorbateVDWAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UHostCationVDWAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UAdsorbateCationVDWAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));

  UHostHostCoulombAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UHostAdsorbateCoulombAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UHostCationCoulombAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UAdsorbateAdsorbateCoulombAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UCationCationCoulombAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UAdsorbateCationCoulombAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));

  UHostBondAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UHostUreyBradleyAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UHostBendAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UHostInversionBendAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UHostTorsionAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UHostImproperTorsionAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UHostBondBondAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UHostBendBendAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UHostBondBendAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UHostBondTorsionAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UHostBendTorsionAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));

  UAdsorbateBondAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UAdsorbateUreyBradleyAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UAdsorbateBondBondAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UAdsorbateBendAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UAdsorbateInversionBendAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UAdsorbateTorsionAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UAdsorbateImproperTorsionAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UAdsorbateBondBondAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UAdsorbateBendBendAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UAdsorbateBondBendAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UAdsorbateBondTorsionAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UAdsorbateBendTorsionAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UAdsorbateIntraVDWAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UAdsorbateIntraChargeChargeAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UAdsorbateIntraChargeBondDipoleAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UAdsorbateIntraBondDipoleBondDipoleAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));

  UCationBondAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UCationUreyBradleyAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UCationBendAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UCationInversionBendAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UCationTorsionAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UCationImproperTorsionAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UCationBondBondAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UCationBendBendAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UCationBondBendAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UCationBondTorsionAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UCationBendTorsionAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UCationIntraVDWAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UCationIntraChargeChargeAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UCationIntraChargeBondDipoleAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UCationIntraBondDipoleBondDipoleAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));

  UHostPolarizationAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UAdsorbatePolarizationAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UCationPolarizationAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UHostBackPolarizationAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UAdsorbateBackPolarizationAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UCationBackPolarizationAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));

  UTailCorrectionAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UDistanceConstraintsAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UAngleConstraintsAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UDihedralConstraintsAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UInversionBendConstraintsAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UOutOfPlaneDistanceConstraintsAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UExclusionConstraintsAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UTotalAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));

  TotalSystemDipoleAverage=(VECTOR**)calloc(NumberOfSystems,sizeof(VECTOR*));
  TotalSystemDipoleSquaredAverage=(VECTOR**)calloc(NumberOfSystems,sizeof(VECTOR*));
  TotalSystemNormDipoleAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  TotalSystemNormDipoleSquaredAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));

  NumberOfMoleculesAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  NumberOfMoleculesSquaredAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  DensityAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));

  BoxAverage=(VECTOR**)calloc(NumberOfSystems,sizeof(VECTOR*));
  BoxAverageAX=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  BoxAverageAY=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  BoxAverageAZ=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  BoxAverageBX=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  BoxAverageBY=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  BoxAverageBZ=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  BoxAverageCX=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  BoxAverageCY=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  BoxAverageCZ=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  BoxLengthAverage=(VECTOR**)calloc(NumberOfSystems,sizeof(VECTOR*));
  AlphaAngleAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  BetaAngleAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  GammaAngleAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  VolumeAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  VolumeSquaredAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));

  TotalEnergyAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  TotalEnergySquaredAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  EnthalpyAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  EnthalpySquaredAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  EnthalpyTimesVolumeAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  EnthalpyTimesEnergyAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));

  TemperatureAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  TemperatureCellAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  TemperatureTranslationAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  TemperatureRotationAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  TemperatureRotationAdsorbateAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  TemperatureTranslationAdsorbateAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));

  TemperatureAdsorbatesAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  TemperatureCationsAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  TemperatureFrameworkAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));

  MolecularPressureAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  MolecularStressTensorAverage=(REAL_MATRIX3x3**)calloc(NumberOfSystems,sizeof(REAL_MATRIX3x3*));
  PressureIdealGasPartAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  PressureExcessPartAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  PressureTailCorrectionAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  PressureAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  UNoseHooverAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  HeatOfVaporization=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  EnergyPerMolecule=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  VolumePerMolecule=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  CompressibilityAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));

  BornTermAverage=(REAL_MATRIX9x9**)calloc(NumberOfSystems,sizeof(REAL_MATRIX9x9*));
  ConfigurationalStressFluctuationTermAverage=(REAL_MATRIX9x9**)calloc(NumberOfSystems,sizeof(REAL_MATRIX9x9*));
  ConfigurationalStressTensorAverage=(REAL_MATRIX3x3**)calloc(NumberOfSystems,sizeof(REAL_MATRIX3x3*));
  StressTensorAverage=(REAL_MATRIX3x3**)calloc(NumberOfSystems,sizeof(REAL_MATRIX3x3*));

  SurfaceAreaFrameworkAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  SurfaceAreaFrameworksAverage=(REAL***)calloc(NumberOfSystems,sizeof(REAL*));
  SurfaceAreaCationsAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  SurfaceAreaCount=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));

  TotalEnergyTimesNumberOfMoleculesAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  HostAdsorbateEnergyTimesNumberOfMoleculesAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  AdsorbateAdsorbateEnergyTimesNumberOfMoleculesAverage=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));

  NumberOfMoleculesPerComponentAverage=(REAL***)calloc(NumberOfSystems,sizeof(REAL**));
  NumberOfExcessMoleculesPerComponentAverage=(REAL***)calloc(NumberOfSystems,sizeof(REAL**));
  DensityPerComponentAverage=(REAL***)calloc(NumberOfSystems,sizeof(REAL**));

  WidomRosenbluthFactorAverage=(REAL***)calloc(NumberOfSystems,sizeof(REAL**));
  WidomRosenbluthFactorCount=(REAL***)calloc(NumberOfSystems,sizeof(REAL**));

  WidomEnergyDifferenceAverage=(REAL***)calloc(NumberOfSystems,sizeof(REAL**));

  WidomEnergyFrameworkAverage=(REAL***)calloc(NumberOfSystems,sizeof(REAL**));
  WidomEnergyFrameworkCount=(REAL***)calloc(NumberOfSystems,sizeof(REAL**));

  PrincipleMomentsOfInertiaAverage=(VECTOR***)calloc(NumberOfSystems,sizeof(VECTOR**));
  PrincipleMomentsOfInertiaCount=(REAL***)calloc(NumberOfSystems,sizeof(REAL**));

  for(i=0;i<NumberOfSystems;i++)
  {
    BlockCount[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));

    UHostHostAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UAdsorbateAdsorbateAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UCationCationAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UHostAdsorbateAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UHostCationAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UAdsorbateCationAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));

    UHostHostVDWAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UAdsorbateAdsorbateVDWAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UCationCationVDWAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UHostAdsorbateVDWAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UHostCationVDWAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UAdsorbateCationVDWAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));

    UHostHostCoulombAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UAdsorbateAdsorbateCoulombAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UCationCationCoulombAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UHostAdsorbateCoulombAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UHostCationCoulombAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UAdsorbateCationCoulombAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));

    UHostBondAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UHostUreyBradleyAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UHostBendAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UHostInversionBendAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UHostTorsionAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UHostImproperTorsionAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UHostBondBondAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UHostBendBendAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UHostBondBendAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UHostBondTorsionAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UHostBendTorsionAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));

    UAdsorbateBondAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UAdsorbateUreyBradleyAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UAdsorbateBondBondAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UAdsorbateBendAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UAdsorbateInversionBendAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UAdsorbateTorsionAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UAdsorbateImproperTorsionAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UAdsorbateBondBondAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UAdsorbateBendBendAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UAdsorbateBondBendAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UAdsorbateBondTorsionAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UAdsorbateBendTorsionAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UAdsorbateIntraVDWAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UAdsorbateIntraChargeChargeAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UAdsorbateIntraChargeBondDipoleAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UAdsorbateIntraBondDipoleBondDipoleAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));

    UCationBondAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UCationUreyBradleyAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UCationBendAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UCationInversionBendAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UCationTorsionAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UCationImproperTorsionAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UCationBondBondAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UCationBendBendAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UCationBondBendAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UCationBondTorsionAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UCationBendTorsionAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UCationIntraVDWAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UCationIntraChargeChargeAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UCationIntraChargeBondDipoleAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UCationIntraBondDipoleBondDipoleAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));

    UHostPolarizationAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UAdsorbatePolarizationAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UCationPolarizationAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UHostBackPolarizationAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UAdsorbateBackPolarizationAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UCationBackPolarizationAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));

    UTailCorrectionAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UDistanceConstraintsAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UAngleConstraintsAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UDihedralConstraintsAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UInversionBendConstraintsAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UOutOfPlaneDistanceConstraintsAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UExclusionConstraintsAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UTotalAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));

    TotalSystemDipoleAverage[i]=(VECTOR*)calloc(NumberOfBlocks,sizeof(VECTOR));
    TotalSystemDipoleSquaredAverage[i]=(VECTOR*)calloc(NumberOfBlocks,sizeof(VECTOR));
    TotalSystemNormDipoleAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    TotalSystemNormDipoleSquaredAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));

    NumberOfMoleculesAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    NumberOfMoleculesSquaredAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    DensityAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));

    BoxAverage[i]=(VECTOR*)calloc(NumberOfBlocks,sizeof(VECTOR));
    BoxAverageAX[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    BoxAverageAY[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    BoxAverageAZ[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    BoxAverageBX[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    BoxAverageBY[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    BoxAverageBZ[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    BoxAverageCX[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    BoxAverageCY[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    BoxAverageCZ[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    BoxLengthAverage[i]=(VECTOR*)calloc(NumberOfBlocks,sizeof(VECTOR));
    AlphaAngleAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    BetaAngleAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    GammaAngleAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    VolumeAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    VolumeSquaredAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));

    TotalEnergyAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    TotalEnergySquaredAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    EnthalpyAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    EnthalpySquaredAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    EnthalpyTimesVolumeAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    EnthalpyTimesEnergyAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));

    TemperatureAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    TemperatureCellAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    TemperatureTranslationAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    TemperatureRotationAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    TemperatureRotationAdsorbateAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    TemperatureTranslationAdsorbateAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));

    TemperatureAdsorbatesAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    TemperatureCationsAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    TemperatureFrameworkAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));

    MolecularPressureAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    MolecularStressTensorAverage[i]=(REAL_MATRIX3x3*)calloc(NumberOfBlocks,sizeof(REAL_MATRIX3x3));
    PressureIdealGasPartAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    PressureExcessPartAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    PressureTailCorrectionAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    PressureAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    UNoseHooverAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    HeatOfVaporization[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    EnergyPerMolecule[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    VolumePerMolecule[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    CompressibilityAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));

    BornTermAverage[i]=(REAL_MATRIX9x9*)calloc(NumberOfBlocks,sizeof(REAL_MATRIX9x9));
    ConfigurationalStressFluctuationTermAverage[i]=(REAL_MATRIX9x9*)calloc(NumberOfBlocks,sizeof(REAL_MATRIX9x9));
    ConfigurationalStressTensorAverage[i]=(REAL_MATRIX3x3*)calloc(NumberOfBlocks,sizeof(REAL_MATRIX3x3));
    StressTensorAverage[i]=(REAL_MATRIX3x3*)calloc(NumberOfBlocks,sizeof(REAL_MATRIX3x3));

    SurfaceAreaFrameworkAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    SurfaceAreaFrameworksAverage[i]=(REAL**)calloc(Framework[i].NumberOfFrameworks,sizeof(REAL*));
    for(j=0;j<Framework[i].NumberOfFrameworks;j++)
       SurfaceAreaFrameworksAverage[i][j]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    SurfaceAreaCationsAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    SurfaceAreaCount[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));

    NumberOfMoleculesPerComponentAverage[i]=(REAL**)calloc(NumberOfComponents,sizeof(REAL*));
    NumberOfExcessMoleculesPerComponentAverage[i]=(REAL**)calloc(NumberOfComponents,sizeof(REAL*));
    DensityPerComponentAverage[i]=(REAL**)calloc(NumberOfComponents,sizeof(REAL*));

    TotalEnergyTimesNumberOfMoleculesAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    HostAdsorbateEnergyTimesNumberOfMoleculesAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    AdsorbateAdsorbateEnergyTimesNumberOfMoleculesAverage[i]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));

    WidomRosenbluthFactorAverage[i]=(REAL**)calloc(NumberOfComponents,sizeof(REAL*));
    WidomRosenbluthFactorCount[i]=(REAL**)calloc(NumberOfComponents,sizeof(REAL*));

    WidomEnergyDifferenceAverage[i]=(REAL**)calloc(NumberOfComponents,sizeof(REAL*));

    WidomEnergyFrameworkAverage[i]=(REAL**)calloc(NumberOfComponents,sizeof(REAL*));
    WidomEnergyFrameworkCount[i]=(REAL**)calloc(NumberOfComponents,sizeof(REAL*));

    PrincipleMomentsOfInertiaAverage[i]=(VECTOR**)calloc(NumberOfComponents,sizeof(VECTOR*));
    PrincipleMomentsOfInertiaCount[i]=(REAL**)calloc(NumberOfComponents,sizeof(REAL*));

    for(j=0;j<NumberOfComponents;j++)
    {
      NumberOfMoleculesPerComponentAverage[i][j]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
      NumberOfExcessMoleculesPerComponentAverage[i][j]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
      DensityPerComponentAverage[i][j]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));

      WidomRosenbluthFactorAverage[i][j]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
      WidomRosenbluthFactorCount[i][j]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));

      WidomEnergyDifferenceAverage[i][j]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));

      WidomEnergyFrameworkAverage[i][j]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
      WidomEnergyFrameworkCount[i][j]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));

      PrincipleMomentsOfInertiaAverage[i][j]=(VECTOR*)calloc(NumberOfBlocks,sizeof(VECTOR));
      PrincipleMomentsOfInertiaCount[i][j]=(REAL*)calloc(NumberOfBlocks,sizeof(REAL));
    }
  }
}

void ReadRestartStatistics(FILE *FilePtr)
{
  int i,j;
  int NumberOfBlocks;
  REAL Check;

  NumberOfBlocks=NR_BLOCKS;

  AllocateStatisticsMemory();

  fread(&NumberOfBlocks,sizeof(int),1,FilePtr);
  fread(&Block,sizeof(int),1,FilePtr);

  fread(BlockCycle,sizeof(long long),NumberOfBlocks,FilePtr);

  for(i=0;i<NumberOfSystems;i++)
  {

    fread(BlockCount[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    fread(UHostHostAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UAdsorbateAdsorbateAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UCationCationAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UHostAdsorbateAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UHostCationAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UAdsorbateCationAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    fread(UHostHostVDWAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UAdsorbateAdsorbateVDWAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UCationCationVDWAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UHostAdsorbateVDWAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UHostCationVDWAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UAdsorbateCationVDWAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    fread(UHostHostCoulombAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UAdsorbateAdsorbateCoulombAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UCationCationCoulombAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UHostAdsorbateCoulombAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UHostCationCoulombAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UAdsorbateCationCoulombAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    fread(UHostBondAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UHostUreyBradleyAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UHostBendAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UHostInversionBendAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UHostTorsionAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UHostImproperTorsionAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UHostBondBondAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UHostBendBendAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UHostBondBendAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UHostBondTorsionAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UHostBendTorsionAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    fread(UAdsorbateBondAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UAdsorbateUreyBradleyAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UAdsorbateBondBondAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UAdsorbateBendAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UAdsorbateInversionBendAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UAdsorbateTorsionAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UAdsorbateImproperTorsionAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UAdsorbateBondBondAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UAdsorbateBendBendAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UAdsorbateBondBendAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UAdsorbateBondTorsionAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UAdsorbateBendTorsionAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UAdsorbateIntraVDWAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UAdsorbateIntraChargeChargeAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UAdsorbateIntraChargeBondDipoleAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UAdsorbateIntraBondDipoleBondDipoleAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    fread(UCationBondAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UCationUreyBradleyAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UCationBendAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UCationInversionBendAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UCationTorsionAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UCationImproperTorsionAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UCationBondBondAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UCationBendBendAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UCationBondBendAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UCationBondTorsionAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UCationBendTorsionAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UCationIntraVDWAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UCationIntraChargeChargeAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UCationIntraChargeBondDipoleAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UCationIntraBondDipoleBondDipoleAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    fread(UHostPolarizationAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UAdsorbatePolarizationAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UCationPolarizationAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UHostBackPolarizationAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UAdsorbateBackPolarizationAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UCationBackPolarizationAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    fread(UTailCorrectionAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UDistanceConstraintsAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UAngleConstraintsAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UDihedralConstraintsAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UInversionBendConstraintsAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UOutOfPlaneDistanceConstraintsAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UExclusionConstraintsAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UTotalAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    fread(TotalSystemDipoleAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(TotalSystemDipoleSquaredAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(TotalSystemNormDipoleAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(TotalSystemNormDipoleSquaredAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    fread(NumberOfMoleculesAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(NumberOfMoleculesSquaredAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(DensityAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    fread(BoxAverage[i],sizeof(VECTOR),NumberOfBlocks,FilePtr);
    fread(BoxAverageAX[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(BoxAverageAY[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(BoxAverageAZ[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(BoxAverageBX[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(BoxAverageBY[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(BoxAverageBZ[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(BoxAverageCX[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(BoxAverageCY[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(BoxAverageCZ[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(BoxLengthAverage[i],sizeof(VECTOR),NumberOfBlocks,FilePtr);
    fread(AlphaAngleAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(BetaAngleAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(GammaAngleAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(VolumeAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(VolumeSquaredAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    fread(TotalEnergyAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(TotalEnergySquaredAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(EnthalpyAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(EnthalpySquaredAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(EnthalpyTimesVolumeAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(EnthalpyTimesEnergyAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    fread(TemperatureAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(TemperatureCellAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(TemperatureTranslationAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(TemperatureRotationAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(TemperatureRotationAdsorbateAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(TemperatureTranslationAdsorbateAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    fread(TemperatureAdsorbatesAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(TemperatureCationsAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(TemperatureFrameworkAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    fread(MolecularPressureAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(MolecularStressTensorAverage[i],NumberOfBlocks,sizeof(REAL_MATRIX3x3),FilePtr);
    fread(PressureIdealGasPartAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(PressureExcessPartAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(PressureTailCorrectionAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(PressureAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(UNoseHooverAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(HeatOfVaporization[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(EnergyPerMolecule[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(VolumePerMolecule[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(CompressibilityAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    fread(BornTermAverage[i],sizeof(REAL_MATRIX9x9),NumberOfBlocks,FilePtr);
    fread(ConfigurationalStressFluctuationTermAverage[i],sizeof(REAL_MATRIX9x9),NumberOfBlocks,FilePtr);
    fread(ConfigurationalStressTensorAverage[i],sizeof(REAL_MATRIX3x3),NumberOfBlocks,FilePtr);
    fread(StressTensorAverage[i],sizeof(REAL_MATRIX3x3),NumberOfBlocks,FilePtr);

    fread(SurfaceAreaFrameworkAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    for(j=0;j<Framework[i].NumberOfFrameworks;j++)
      fread(SurfaceAreaFrameworksAverage[i][j],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(SurfaceAreaCationsAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(SurfaceAreaCount[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    fread(TotalEnergyTimesNumberOfMoleculesAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(HostAdsorbateEnergyTimesNumberOfMoleculesAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);
    fread(AdsorbateAdsorbateEnergyTimesNumberOfMoleculesAverage[i],sizeof(REAL),NumberOfBlocks,FilePtr);

    for(j=0;j<NumberOfComponents;j++)
    {
      fread(NumberOfMoleculesPerComponentAverage[i][j],sizeof(REAL),NumberOfBlocks,FilePtr);
      fread(NumberOfExcessMoleculesPerComponentAverage[i][j],sizeof(REAL),NumberOfBlocks,FilePtr);
      fread(DensityPerComponentAverage[i][j],sizeof(REAL),NumberOfBlocks,FilePtr);

      fread(WidomRosenbluthFactorAverage[i][j],sizeof(REAL),NumberOfBlocks,FilePtr);
      fread(WidomRosenbluthFactorCount[i][j],sizeof(REAL),NumberOfBlocks,FilePtr);

      fread(WidomEnergyDifferenceAverage[i][j],sizeof(REAL),NumberOfBlocks,FilePtr);

      fread(WidomEnergyFrameworkAverage[i][j],sizeof(REAL),NumberOfBlocks,FilePtr);
      fread(WidomEnergyFrameworkCount[i][j],sizeof(REAL),NumberOfBlocks,FilePtr);

      fread(PrincipleMomentsOfInertiaAverage[i][j],sizeof(REAL),NumberOfBlocks,FilePtr);
      fread(PrincipleMomentsOfInertiaCount[i][j],sizeof(REAL),NumberOfBlocks,FilePtr);
    }
  }

  fread(&Check,1,sizeof(REAL),FilePtr);
  if(fabs(Check-123456789.0)>1e-10)
  {
    fprintf(stderr, "Error in binary restart-file (ReadRestartStatistics)\n");
    ContinueAfterCrash=FALSE;
  }
}
