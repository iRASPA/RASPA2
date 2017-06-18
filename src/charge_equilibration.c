/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2017 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'charge_equilibration.c' is part of RASPA-2.0

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

// July 2011: linear charge equilibration code by Chris Wilmer
// September 2011: Expansion around oxidation state by Lifeng Ding and Ozgur Yazaydin

// Literature:
// ===========
// C.E. Wilmer and R.Q. Snurr,  "Towards rapid computational screening of metal-organic frameworks for carbon dioxide capture:
//    Calculation of framework charges via charge equilibration",Chem. Eng. J., 171(3), 775-781,2011
// A.K. Rappe and Goddard, "charge equilibration for molecular-dynamics simulation", J. Phys. Chem., 95(8), 3358-3363, 1991


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>
#include "molecule.h"
#include "framework.h"
#include "charge_equilibration.h"
#include "potentials.h"
#include "scattering_factors.h"

int ChargeEquilibrationPeriodic;
int ChargeEquilibrationEwald;

static int size;
enum {OKAY,WARN_EA,WARN_IP2,WARN_IP3,WARN_EA_IP3};

void SolveMatrix(REAL_MATRIX A,REAL *b,REAL *x);
int GetChargeEquilibrationIndexFromElementName(char *s);

REAL *J,*X,*Xc,*R;
int *HydrogenList;
int HydroIndex;
VECTOR *Pos;
static REAL Qtot=0.0;

typedef struct
{
  int index;
  const char  *Label;
  int State;
  REAL ElectronAffinity;            // in units of [eV]
  int NumberOfIonizationEnergies;
  REAL IonizationEnergy[8];         // in units of [eV]
} CHARGE_EQUILIBRATION_STRUCT;

#define NUMBER_OF_ELECTRONEGATIVITIES_RG 16

CHARGE_EQUILIBRATION_STRUCT ChargeEquilibrationDataRappeGoddard[NUMBER_OF_ELECTRONEGATIVITIES_RG] =
{
  { 1, "H",OKAY,      -2.4172, 1,{11.4732,0.0,0.0,0.0,0.0,0.0,0.0,0.0}},
  { 3,"Li",OKAY,       0.62, 1,{5.392,0.0,0.0,0.0,0.0,0.0,0.0,0.0}},
  { 6, "C",OKAY,       0.28, 1,{10.406,0.0,0.0,0.0,0.0,0.0,0.0,0.0}},
  { 7, "N",OKAY,       1.019, 1,{12.779,0.0,0.0,0.0,0.0,0.0,0.0,0.0}},
  { 8, "O",OKAY,       2.059, 1,{15.423,0.0,0.0,0.0,0.0,0.0,0.0,0.0}},
  { 9, "F",OKAY,       3.4, 1,{18.348,0.0,0.0,0.0,0.0,0.0,0.0,0.0}},
  {11,"Na",OKAY,       0.547, 1,{5.139,0.0,0.0,0.0,0.0,0.0,0.0,0.0}},
  {14,"Si",OKAY,       0.681, 1,{7.655,0.0,0.0,0.0,0.0,0.0,0.0,0.0}},
  {15, "P",OKAY,       1.463, 1,{9.463,0.0,0.0,0.0,0.0,0.0,0.0,0.0}},
  {16, "S",OKAY,       2.442, 1,{11.414,0.0,0.0,0.0,0.0,0.0,0.0,0.0}},
  {17,"Cl",OKAY,       3.618, 1,{13.51,0.0,0.0,0.0,0.0,0.0,0.0,0.0}},
  {19, "K",OKAY,       0.501, 1,{4.341,0.0,0.0,0.0,0.0,0.0,0.0,0.0}},
  {35,"Br",OKAY,       3.365, 1,{12.215,0.0,0.0,0.0,0.0,0.0,0.0,0.0}},
  {37,"Rb",OKAY,       0.485, 1,{4.177,0.0,0.0,0.0,0.0,0.0,0.0,0.0}},
  {53, "I",OKAY,       3.06, 1,{10.584,0.0,0.0,0.0,0.0,0.0,0.0,0.0}},
  {55,"Cs",OKAY,       0.472, 1,{3.894,0.0,0.0,0.0,0.0,0.0,0.0,0.0}}
};

// data taken from Handbook Chemistry and Physics, edition 92, 2011-2012
// electron affinities: 10-147,10-148,10-149
// ionization energies: 10-196,10-197
#define NUMBER_OF_ELECTRONEGATIVITIES 84
CHARGE_EQUILIBRATION_STRUCT ChargeEquilibrationData[NUMBER_OF_ELECTRONEGATIVITIES] =
{
  { 1, "H",OKAY,       0.7520812, 1,{13.598443,0.0,0.0,0.0,0.0,0.0,0.0,0.0}},
  { 2,"He",OKAY,       0.00000, 2,{24.587387,54.417760,0.0,0.0,0.0,0.0,0.0,0.0}},
  { 3,"Li",OKAY,       0.618049, 3,{5.391719,75.6400,122.45429,0.0,0.0,0.0,0.0,0.0}},
  { 4,"Be",OKAY,       0.00000, 4,{9.32270,18.21114,153.89661,217.71865,0.0,0.0,0.0,0.0}},
  { 5, "B",OKAY,       0.279730,5,{8.29802,25.1548,37.93064,259.37521,340.22580,0.0,0.0,0.0}},
  { 6, "C",OKAY,       1.262119,6,{11.26030,24.3833,47.8878,64.4939,392.087,489.99334,0.0,0.0}},
  { 7, "N",OKAY,       0.000000,7,{14.5341,29.6013,47.44924,77.4735,97.8902,552.0718,667.046,0.0}},
  { 8, "O",OKAY,       1.4611135,8,{13.61805,35.1211,54.9355,77.41353,113.8990,138.1197,739.29,871.4101}},
  { 9, "F",OKAY,       3.4011895,8,{17.4228,34.9708,62.7084,87.1398,114.2428,157.1651,185.186,953.9112}},
  {10,"Ne",OKAY,       0.000000,8,{21.56454,40.96296,63.45,97.12,126.21,157.93,207.2759,239.0989}},
  {11,"Na",OKAY,       0.547926,8,{5.139076,47.2864,71.6200,98.91,138.40,172.18,208.50,264.25}},
  {12,"Mg",OKAY,       0.000000,8,{7.646235,15.03527,80.1437,109.2655,141.27,186.76,225.02,265.96}},
  {13,"Al",OKAY,       0.43283,8,{5.985768,18.82855,28.44765,119.992,153.825,190.49,241.76,284.66}},
  {14,"Si",OKAY,       1.3895213,8,{8.15168,16.34584,33.49302,45.14181,166.767,205.27,246.5,303.54}},
  {15, "P",OKAY,       0.7465,8,{10.48669,19.7695,30.2027,51.4439,65.0251,220.421,263.57,309.60}},
  {16, "S",OKAY,       2.07710418,8,{10.36001,23.33788,34.79,47.222,72.5945,88.0530,280.948,328.75}},
  {17,"Cl",OKAY,       3.612724,8,{12.96763,23.8136,39.61,53.4652,67.8,97.03,114.1958,348.28}},
  {18,"Ar",OKAY,       0.000000,8,{15.759610,27.62966,40.64,59.81,75.02,91.009,124.323,143.460}},
  {19, "K",OKAY,       0.50147,8,{4.3406633,31.63,45.806,60.91,82.66,99.4,117.56,154.88}},
  {20,"Ca",OKAY,       0.02455,8,{6.11316,11.87172,50.9131,67.27,84.50,108.78,127.2,147.24}},
  {21,"Sc",OKAY,       0.188,8,{6.56149,12.79977,24.75666,73.4894,91.65,110.68,138.0,158.1}},
  {22,"Ti",OKAY,       0.079,8,{6.82812,13.5755,27.4917,43.2672,99.30,119.53,140.8,170.4}},
  {23, "V",OKAY,       0.525,8,{6.74619,14.618,29.311,46.709,65.2817,128.13,150.6,173.4}},
  {24,"Cr",OKAY,       0.666,8,{6.76651,16.4857,30.96,49.16,69.46,90.6349,160.18,184.7}},
  {25,"Mn",OKAY,       0.000000,8,{7.43402,15.6400,33.668,51.2,72.4,95.6,119.203,194.5}},
  {26,"Fe",OKAY,       0.151,8,{7.9024,16.1877,30.652,54.8,75.0,99.1,124.98,151.06}},
  {27,"Co",OKAY,       0.662,8,{7.88101,17.084,33.50,51.3,79.5,102.0,128.9,157.8}},
  {28,"Ni",OKAY,       1.156,8,{7.6398,18.16884,35.19,54.9,76.06,108,133,162}},
  {29,"Cu",OKAY,       1.235,8,{7.72638,20.2924,36.841,57.38,79.8,103,139,166}},
  {30,"Zn",OKAY,       0.0,8,{9.394199,17.96439,39.723,59.4,82.6,108,134,174}},
  {31,"Ga",OKAY,       0.43,8,{5.999301,20.51515,30.7258,63.241,86.01,112.7,140.9,169.9}},
  {32,"Ge",OKAY,       1.232712,5,{7.89943,15.93461,34.2241,45.7131,93.5,0.0,0.0,0.0}},
  {33,"As",OKAY,       0.804,6,{9.7886,18.5892,28.351,50.13,62.63,127.6,0.0,0.0}},
  {34,"Se",OKAY,       2.02067,7,{9.75239,21.19,30.8204,42.9450,68.3,81.7,155.4,0.0}},
  {35,"Br",OKAY,       3.363588,8,{11.8138,21.591,36,47.3,59.7,88.6,103.0,192.8}},
  {36,"Kr",OKAY,       0.000000,8,{13.99961,24.35984,36.950,52.5,64.7,78.5,111.0,125.802}},
  {37,"Rb",OKAY,       0.48592,8,{4.177128,27.2895,40,52.6,71.0,84.4,99.2,136}},
  {38,"Sr",OKAY,       0.048,8,{5.69485,11.0301,42.89,57,71.6,90.8,106,122.3}},
  {39, "Y",OKAY,       0.307,8,{6.2173,12.224,20.52,60.597,77.0,93.0,116,129}},
  {40,"Zr",OKAY,       0.426,5,{6.63390,13.1,22.99,34.34,80.348,0.0,0.0,0.0}},
  {41,"Nb",OKAY,       0.916,7,{6.75885,14.0,25.04,38.3,50.55,102.057,125,0.0}},
  {42,"Mo",OKAY,       0.748,8,{7.09243,16.16,27.13,46.4,54.49,68.8276,125.664,143.6}},
  {43,"Tc",OKAY,       0.55,3,{7.28,15.26,29.54,0.0,0.0,0.0,0.0,0.0}},
  {44,"Ru",OKAY,       1.05,3,{7.36050,16.76,28.47,0.0,0.0,0.0,0.0,0.0}},
  {45,"Rh",OKAY,       1.137,3,{7.45890,18.08,31.06,0.0,0.0,0.0,0.0,0.0}},
  {46,"Pd",OKAY,       0.562,3,{8.3369,19.43,32.93,0.0,0.0,0.0,0.0,0.0}},
  {47,"Ag",OKAY,       1.302,3,{7.57623,21.47746,34.83,0.0,0.0,0.0,0.0,0.0}},
  {48,"Cd",OKAY,       0.000000,3,{8.99382,16.90831,37.48,0.0,0.0,0.0,0.0,0.0}},
  {49,"In",OKAY,       0.3,4,{5.78636,18.8703,28.03,54,0.0,0.0,0.0,0.0}},
  {50,"Sn",OKAY,       1.112067,5,{7.34392,14.6322,30.50260,40.73502,72.28,0.0,0.0,0.0}},
  {51,"Sb",OKAY,       1.046,6,{8.60839,16.63,25.3,44.2,56,108,0.0,0.0}},
  {52,"Te",OKAY,       1.970876,7,{9.0096,18.6,27.96,37.41,58.75,70.7,137,0.0}},
  {53, "I",OKAY,       3.0590463,3,{10.45126,19.1313,33,0.0,0.0,0.0,0.0,0.0}},
  {54,"Xe",OKAY,       0.000000,3,{12.12984,20.9750,32.1230,0.0,0.0,0.0,0.0,0.0}},
  {55,"Cs",WARN_IP3,   0.471626,2,{3.893905,23.15744,0.0,0.0,0.0,0.0,0.0,0.0}},
  {56,"Ba",WARN_IP3,   0.14462,2,{5.211664,10.00383,0.0,0.0,0.0,0.0,0.0,0.0}},
  {57,"La",OKAY,       0.47,5,{5.5769,11.059,19.1773,49.950,61.6,0.0,0.0,0.0}},
  {58,"Ce",OKAY,       0.65,6,{5.5387,10.85,20.198,36.758,65.55,77.6,0.0,0.0}},
  {59,"Pr",WARN_EA,    0.962,5,{5.5473,10.55,21.624,38.98,57.53,0.0,0.0,0.0}},
  {60,"Nd",WARN_EA_IP3,1.916,4,{5.5250,10.72,22.1,40.4,0.0,0.0,0.0,0.0}},
  {61,"Pm",WARN_EA_IP3,0.000000,4,{5.582,10.9,22.3,41.1,0.0,0.0,0.0,0.0}},
  {62,"Sm",WARN_EA_IP3,0.000000,4,{5.6437,11.07,23.4,41.4,0.0,0.0,0.0,0.0}},
  {63,"Eu",WARN_EA_IP3,0.864,4,{5.67038,11.25,24.92,42.7,0.0,0.0,0.0,0.0}},
  {64,"Gd",WARN_EA_IP3,0.000000,4,{6.14980,12.09,20.63,44.0,0.0,0.0,0.0,0.0}},
  {65,"Tb",WARN_EA_IP3,1.165,4,{5.8638,11.52,21.91,39.79,0.0,0.0,0.0,0.0}},
  {66,"Dy",WARN_EA_IP3,0.000000,4,{5.9389,11.67,22.8,41.47,0.0,0.0,0.0,0.0}},
  {67,"Ho",WARN_EA_IP3,0.000000,4,{6.0215,11.8,22.84,42.5,0.0,0.0,0.0,0.0}},
  {68,"Er",WARN_EA_IP3,0.000000,4,{6.1077,11.93,22.74,42.7,0.0,0.0,0.0,0.0}},
  {69,"Tm",WARN_EA,    1.029,4,{6.18431,12.05,23.68,42.7,0.0,0.0,0.0,0.0}},
  {70,"Yb",WARN_EA,    -0.020,4,{6.25416,12.176,25.05,43.56,0.0,0.0,0.0,0.0}},
  {71,"Lu",WARN_EA_IP3,0.34,5,{5.42586,13.9,20.9594,45.25,66.8,0.0,0.0,0.0}},
  {72,"Hf",OKAY,       0.017,4,{6.82507,15,23.3,33.33,0.0,0.0,0.0,0.0}},
  {73,"Ta",WARN_IP2,   0.322000,1,{7.54957,0.0,0.0,0.0,0.0,0.0,0.0,0.0}},
  {74, "W",WARN_IP2,   0.81626,1,{7.86403,16.1,0.0,0.0,0.0,0.0,0.0,0.0}},
  {75,"Re",WARN_IP2,   0.15,1,{7.83352,0.0,0.0,0.0,0.0,0.0,0.0,0.0}},
  {76,"Os",WARN_IP2,   1.1,1,{8.43823,0.0,0.0,0.0,0.0,0.0,0.0,0.0}},
  {77,"Ir",WARN_IP2,   1.5638,1,{8.96702,0.0,0.0,0.0,0.0,0.0,0.0,0.0}},
  {78,"Pt",WARN_IP3,   2.128,2,{8.9588,18.563,0.0,0.0,0.0,0.0,0.0,0.0}},
  {79,"Au",WARN_IP3,   2.30863,2,{9.22553,20.2,0.0,0.0,0.0,0.0,0.0,0.0}},
  {80,"Hg",OKAY,       0.000000,3,{10.4375,18.7568,34.2,0.0,0.0,0.0,0.0,0.0}},
  {81,"Tl",OKAY,       0.377,3,{6.108194,20.4283,29.83,0.0,0.0,0.0,0.0,0.0}},
  {82,"Pb",OKAY,       0.364,5,{7.41663,15.03248,31.9373,42.32,68.8,0.0,0.0,0.0}},
  {83,"Bi",OKAY,       0.942362,6,{7.2855,16.703,25.56,45.3,56.0,88.3,0.0,0.0}},
  {84,"Po",WARN_IP2,   1.9,1,{8.414,0.0,0.0,0.0,0.0,0.0,0.0,0.0}}
};

// Values taken from "Absolute Electronegativity and Hardness: Application to Inorganic Chemistry"
// Inorganic Chemistry, 1988, Vol. 27, pg 734 - 740, Author: Ralph G. Pearson
REAL ReferenceTableX(int type,int index)
{
  REAL TempEA,TempIP;

  if(PseudoAtoms[type].OxidationState==0)
  {
    TempEA=ChargeEquilibrationData[index].ElectronAffinity;
    TempIP=ChargeEquilibrationData[index].IonizationEnergy[0];
  }
  else
  {
    TempEA=ChargeEquilibrationData[index].IonizationEnergy[PseudoAtoms[type].OxidationState-1];
    TempIP=ChargeEquilibrationData[index].IonizationEnergy[PseudoAtoms[type].OxidationState];
  }

  return 0.5*(TempEA+TempIP);
}

// Values taken from "Absolute Electronegativity and Hardness: Application to Inorganic Chemistry"
// Inorganic Chemistry, 1988, Vol. 27, pg 734 - 740, Author: Ralph G. Pearson
// *Note*: J is defined as twice that of eta. All values are in eV (electron volts).
REAL ReferenceTableJ(int type,int index)
{
  REAL TempEA, TempIP;

  if(PseudoAtoms[type].OxidationState==0)
  {
    TempEA=ChargeEquilibrationData[index].ElectronAffinity;
    TempIP=ChargeEquilibrationData[index].IonizationEnergy[0];
  }
  else
  {
    TempEA=ChargeEquilibrationData[index].IonizationEnergy[PseudoAtoms[type].OxidationState-1];
    TempIP=ChargeEquilibrationData[index].IonizationEnergy[PseudoAtoms[type].OxidationState];
  }

  return (TempIP-TempEA);
}

//  ReferenceTableXc is for the taylor expansion center Xc = J*a
REAL ReferenceTableXc(int type,int index)
{
  REAL TempEA,TempIP;

  if(PseudoAtoms[type].OxidationState==0)
  {
    TempEA=ChargeEquilibrationData[index].ElectronAffinity;
    TempIP=ChargeEquilibrationData[index].IonizationEnergy[0];
  }
  else
  {
    TempEA=ChargeEquilibrationData[index].IonizationEnergy[PseudoAtoms[type].OxidationState-1];
    TempIP=ChargeEquilibrationData[index].IonizationEnergy[PseudoAtoms[type].OxidationState];
  }

  return (TempIP-TempEA)*PseudoAtoms[type].OxidationState;
}

int GetChargeEquilibrationIndexFromElementName(char *s)
{
  int i;

  for(i=0;i<NUMBER_OF_ELECTRONEGATIVITIES;i++)
    if(strcasecmp(s,ChargeEquilibrationData[i].Label)==0) return i;
  return -1;
}

REAL CoulombIntegral(REAL r, REAL a, REAL b)
{
  REAL orbitalOverlapTerm;

  orbitalOverlapTerm=exp(-(a/b)*r*r)*(1.0/r-a*(1.0-r/b));

  return 1.0/r-orbitalOverlapTerm;
}


REAL getJ(int i, int j)
{
  int u,v,w;
  VECTOR dr,kv;
  REAL r,r2,Jab;
  REAL a,orbitalOverlapTerm;
  REAL sigmaStar,sigma;
  REAL orbital,alphaStar,betaStar;
  int aVnum,bVnum,cVnum;
  REAL minCellLength;
  REAL h,b;
  REAL alpha,beta;

//  aVnum=3;
//  bVnum=3;
//  cVnum=3;
  minCellLength = 250;
  aVnum = ceil(minCellLength / (2.0 * Box[0].ax) ) - 1;
  bVnum = ceil(minCellLength / (2.0 * Box[0].by) ) - 1;
  cVnum = ceil(minCellLength / (2.0 * Box[0].cz) ) - 1;

  REAL k = 14.4; // [Angstroms * electron volts]
  REAL lambda = 1.2; // Global hardness scaling parameter
  REAL eta=50.0;

  if(ChargeEquilibrationPeriodic==FALSE)
  {
    //////////////////////////////////////////////////////////////////////
    //  NonPeriodic                                                     //
    //////////////////////////////////////////////////////////////////////
    if (i == j)
    {
      return J[i]; // Return the hardness/idempotential
    }
    else
    {
      dr.x=Pos[i].x-Pos[j].x;
      dr.y=Pos[i].y-Pos[j].y;
      dr.z=Pos[i].z-Pos[j].z;

      // apply boundary conditions
      //dr=ApplyBoundaryCondition(dr);

      r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
      r=sqrt(r2);

      a=sqrt(J[i]*J[j])/k;
      orbitalOverlapTerm=exp(-(a*a*r2))*(2.0*a-a*a*r-1.0/r);

      Jab=lambda*(k/2.0)*((1.0/r)+orbitalOverlapTerm);

      return Jab;
    }
  }
  else if(ChargeEquilibrationPeriodic==TRUE)
  {
    if(ChargeEquilibrationEwald==FALSE)
    {
      //////////////////////////////////////////////////////////////////////
      // Direct sums                                                      //
      //////////////////////////////////////////////////////////////////////
      if(i==j)
      {
        sigmaStar=0;
        for(u=-aVnum;u<=aVnum;u++)
        {
          for(v=-bVnum;v<=bVnum;v++)
          {
            for(w=-cVnum;w<=cVnum;w++)
            {
              if (!((u==0)&&(v==0)&&(w==0)))
              {
                dr.x=u*Box[0].ax+v*Box[0].bx+w*Box[0].cx;
                dr.y=u*Box[0].ay+v*Box[0].by+w*Box[0].cy;
                dr.z=u*Box[0].az+v*Box[0].bz+w*Box[0].cz;
                r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
                r=sqrt(r2);

                a=sqrt(J[i]*J[j])/k;
                orbitalOverlapTerm=exp(-(a*a*r2))*(2.0*a-a*a*r-1.0/r);

                sigmaStar+=(1.0/r)+orbitalOverlapTerm;
              }
            }
          }
        }
        return J[i]+lambda*(k/2.0)*sigmaStar;
      }
      else
      {
        sigma=0;
        for(u=-aVnum;u<= aVnum;u++)
        {
          for(v=-bVnum;v<= bVnum;v++)
          {
            for(w=-cVnum;w<=cVnum;w++)
            {
              dr.x=Pos[i].x-Pos[j].x+u*Box[0].ax+v*Box[0].bx+w*Box[0].cx;
              dr.y=Pos[i].y-Pos[j].y+u*Box[0].ay+v*Box[0].by+w*Box[0].cy;
              dr.z=Pos[i].z-Pos[j].z+u*Box[0].az+v*Box[0].bz+w*Box[0].cz;
              r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
              r=sqrt(r2);

              a=sqrt(J[i]*J[j])/k;
              orbitalOverlapTerm=exp(-(a*a*r2))*(2.0*a-a*a*r-1.0/r);

              sigma+=(1.0/r)+orbitalOverlapTerm;
            }
          }
        }
        return lambda*(k/2.0)*sigma;
      }
    }
    else
    {
      //////////////////////////////////////////////////////////////////////
      // Ewald sums                                                       //
      //////////////////////////////////////////////////////////////////////
      if (i==j)
      {
        // Orbital energy term
        orbital = 0;
        for(u=-aVnum;u<=aVnum;u++)
        {
          for(v=-bVnum;v<=bVnum;v++)
          {
            for(w=-cVnum;w<=cVnum;w++)
            {
              if((u==0)&&(v==0)&&(w==0))
              {
                // do nothing
              }
              else
              {
                dr.x=u*Box[0].ax+v*Box[0].bx+w*Box[0].cx;
                dr.y=u*Box[0].ay+v*Box[0].by+w*Box[0].cy;
                dr.z=u*Box[0].az+v*Box[0].bz+w*Box[0].cz;
                r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
                r=sqrt(r2);

                a=sqrt(J[i]*J[j])/k;
                orbitalOverlapTerm=exp(-(a*a*r2))*(2.0*a-a*a*r-1.0/r);

                orbital+=orbitalOverlapTerm;
              }
            }
          }
        }

        // Real-space Coulomb component
        alphaStar=0;
        for(u=-aVnum;u<=aVnum;u++)
        {
          for(v=-bVnum;v<=bVnum;v++)
          {
            for(w=-cVnum;w<=cVnum;w++)
            {
              if((u==0)&&(v==0)&&(w==0))
              {
                // do nothing
              }
              else
              {
                dr.x=u*Box[0].ax+v*Box[0].bx+w*Box[0].cx;
                dr.y=u*Box[0].ay+v*Box[0].by+w*Box[0].cy;
                dr.z=u*Box[0].az+v*Box[0].bz+w*Box[0].cz;
                r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
                r=sqrt(r2);

                alphaStar+=erfc(r/eta)/r;
              }
            }
          }
        }

        // K-space component
        betaStar=0;
        h=0.0; b=0.0;
        for(u=-aVnum;u<=aVnum;u++)
        {
          for(v=-bVnum;v<=bVnum;v++)
          {
            for(w=-cVnum;w<=cVnum;w++)
            {
              if((u==0)&&(v==0)&&(w==0))
              {
                // do nothing
              }
              else
              {
                kv.x=u*InverseBox[0].ax+v*InverseBox[0].bx+w*InverseBox[0].cx;
                kv.y=u*InverseBox[0].ay+v*InverseBox[0].by+w*InverseBox[0].cy;
                kv.z=u*InverseBox[0].az+v*InverseBox[0].bz+w*InverseBox[0].cz;

                kv.x *= 2 * M_PI;
                kv.y *= 2 * M_PI;
                kv.z *= 2 * M_PI;

                h=sqrt(SQR(kv.x)+SQR(kv.y)+SQR(kv.z));
                b=0.5*h*eta;

                betaStar+=1.0/(h*h)*exp(-b*b);
              }
            }
          }
        }
        betaStar*=4.0*M_PI/Volume[0];

        return J[i]+lambda*(k/2.0)*(alphaStar+betaStar+orbital-2.0/(eta*sqrt(M_PI)));
      }
      else
      {
        // Orbital energy term
        orbital=0.0;
        for(u=-aVnum;u<=aVnum;u++)
        {
          for(v=-bVnum;v<=bVnum;v++)
          {
            for(w=-cVnum;w<=cVnum;w++)
            {
              dr.x=Pos[i].x-Pos[j].x+u*Box[0].ax+v*Box[0].bx+w*Box[0].cx;
              dr.y=Pos[i].y-Pos[j].y+u*Box[0].ay+v*Box[0].by+w*Box[0].cy;
              dr.z=Pos[i].z-Pos[j].z+u*Box[0].az+v*Box[0].bz+w*Box[0].cz;
              r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
              r=sqrt(r2);

              a=sqrt(J[i]*J[j])/k;
              orbitalOverlapTerm=exp(-(a*a*r2))*(2.0*a-a*a*r-1.0/r);

              orbital+=orbitalOverlapTerm;
            }
          }
        }

        // Real-space Coulomb component
        alpha=0.0;
        for(u=-aVnum;u<=aVnum;u++)
        {
          for(v=-bVnum;v<=bVnum;v++)
          {
            for(w=-cVnum;w<=cVnum;w++)
            {
              dr.x=Pos[i].x-Pos[j].x+u*Box[0].ax+v*Box[0].bx+w*Box[0].cx;
              dr.y=Pos[i].y-Pos[j].y+u*Box[0].ay+v*Box[0].by+w*Box[0].cy;
              dr.z=Pos[i].z-Pos[j].z+u*Box[0].az+v*Box[0].bz+w*Box[0].cz;
              r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
              r=sqrt(r2);

              alpha+=erfc(r/eta)/r;
            }
          }
        }

        // K-space component
        beta=0.0;
        h=0.0; b=0.0;
        for(u=-aVnum;u<=aVnum;u++)
        {
          for(v=-bVnum;v<=bVnum;v++)
          {
            for(w=-cVnum;w<=cVnum;w++)
            {
              if((u==0)&&(v==0)&&(w==0))
              {
                // do nothing
              }
              else
              {
                kv.x= u*InverseBox[0].ax+v*InverseBox[0].bx+w*InverseBox[0].cx;
                kv.y= u*InverseBox[0].ay+v*InverseBox[0].by+w*InverseBox[0].cy;
                kv.z= u*InverseBox[0].az+v*InverseBox[0].bz+w*InverseBox[0].cz;

                kv.x *= 2 * M_PI;
                kv.y *= 2 * M_PI;
                kv.z *= 2 * M_PI;

                h=sqrt(SQR(kv.x)+SQR(kv.y)+SQR(kv.z));
                b=0.5*h*eta;

                dr.x=Pos[i].x-Pos[j].x;
                dr.y=Pos[i].y-Pos[j].y;
                dr.z=Pos[i].z-Pos[j].z;

                beta+=cos(kv.x*dr.x+kv.y*dr.y+kv.z*dr.z )/(h*h)*exp(-b*b);
              }
            }
          }
        }
        beta*=4.0*M_PI/Volume[0];

        return lambda*(k/2.0)*(alpha+beta+orbital);
      }
    }
  }
  else
  {
    fprintf(stderr, "Serious error specifying periodic boundary conditions. Exiting\n");
    exit(1);
  }
}

void Qeq(REAL *Q)
{
  int i,j;
  REAL JQSum;

  REAL_MATRIX A;
  REAL *b, *x0;

  A=CreateRealMatrix(size,size);
  b=(REAL*)calloc(size,sizeof(REAL));
  x0=(REAL*)calloc(size,sizeof(REAL));

  // First row of A is all ones
  for(i=0;i<size;i++)
    A.element[0][i]=1.0;

  // First element in b is the total charge
  b[0]=Qtot;

  // Rest of elements in b are the differences in electronegativity
  // Use lattice potential. If it hasn't been calculated yet, it will just default to zero.
  for(i=1;i<size;i++)
    b[i]=(X[i]-Xc[i])-(X[i-1]-Xc[i-1]);

  // Fill in 2nd to Nth rows of A
  for(i=1;i<size;i++)
  {
    for(j=0;j<size;j++)
      A.element[i][j]=getJ(i-1, j)-getJ(i, j);
  }

  SolveMatrix(A,b,Q);

  // compute electronegativity of the atoms
  for(i=0;i<size;i++)
  {
    JQSum=0.0;
    for(j=0;j<size;j++)
    {
      A.element[i][j]=getJ(i, j);
      JQSum=JQSum+A.element[i][j]*Q[j];
    }
    x0[i]=(X[i]-Xc[i])+JQSum;
  }

  for(i=1;i<size;i++)
    if(fabs(x0[i]-x0[0])>1e-5)
    {
      fprintf(stderr, "ERROR in charge-equilibration: The electronegativity of some atoms is different!\n");
      exit(0);
    }

  free(b);
  DeleteRealMatrix(A);
}

// Assumptions: A x = b, A is a MxN matrix, M = rows, N = cols, x is a vector, b is vector
// matrix has more rows than columns
// number of rows of matrix is equal to size of the vector x
void SolveMatrix(REAL_MATRIX A,REAL *b,REAL *x)
{
  int i, j, k;
  int N,M;
  REAL *d;
  REAL aii,alpha,f,ak;
  REAL max_norm;
  REAL r,norm,sum;

  // Initialize x = b

  for(i=0;i<A.m;i++)
    x[i]=b[i];

  N=A.n;
  M=A.m;

  d=(REAL*)calloc(N,sizeof(REAL));

  /* Perform Householder transformation */
  for(i=0;i<N;i++)
  {
    aii=A.element[i][i];
    max_norm=0.0;
    r=0.0;

    for(k=i;k<M;k++)
    {
      r+=A.element[k][i]*A.element[k][i];
    }

    if(r==0)
    {
      fprintf(stderr, "ERROR in charge-equilibration: Matrix is rank deficient.\n");
      exit(0);
    }

    if(A.element[i][i]<0)
      alpha=(-1)*sqrt(r);
    else
      alpha=sqrt(r);

    ak=1.0/(r+alpha*A.element[i][i]);

    A.element[i][i]+=alpha;

    d[i]=-alpha;

    for(k=i+1;k<N;k++)
    {
      norm=0.0;
      f=0.0;

      for(j=i;j<M;j++)
      {
        norm+=A.element[j][k]*A.element[j][k];
        f+=A.element[j][k]*A.element[j][i];
      }

      max_norm=MAX2(max_norm,norm);

      f*=ak;

      for(j=i;j<M;j++)
      {
        A.element[j][k]-=f*A.element[j][i];
      }
    }

    if(fabs(alpha)<0.00001)
    {
      fprintf(stderr, "ERROR in charge-equilibration: apparent singularity in matrix.\n");
      exit(0);
    }

    f=0.0;
    for(j=i;j<M;j++)
    {
      f+=x[j]*A.element[j][i];
    }

    f*=ak;

    for(j=i;j<M;j++)
    {
      x[j]-=f*A.element[j][i];
    }
  }

  /* Perform back-substitution */

  for(i=N-1;i>=0;i--)
  {
    sum=0.0;

    for(k=i+1;k<N;k++)
    {
      sum+=A.element[i][k]*x[k];
    }

    x[i]=(x[i]-sum)/d[i];
  }
}

void ChargeEquilibration(void)
{
  int i,j,f1;
  int index;
  int element_index;
  int type;
  REAL *Q;
  REAL net;
  REAL k = 14.4; // [Angstroms * electron volts]
  REAL gamma2 = 0.5; // Global atomic radii scaling parameter
  REAL total_charge;
  REAL count;
  REAL covalent_radius;
  int atom_type;
  int NumberOfHydrogens;

  size=Framework[CurrentSystem].TotalNumberOfAtoms;

  J=(REAL*)calloc(size,sizeof(REAL));
  X=(REAL*)calloc(size,sizeof(REAL));
  Xc=(REAL*)calloc(size,sizeof(REAL));
  R=(REAL*)calloc(size,sizeof(REAL));
  Q=(REAL*)calloc(size,sizeof(REAL));
  HydrogenList=(int*)calloc(size,sizeof(int));

  Pos=(VECTOR*)calloc(size,sizeof(VECTOR));

  // overwrite hydrogen electron affinity and ionization potential
  ChargeEquilibrationData[0].ElectronAffinity=-2.4172;
  ChargeEquilibrationData[0].IonizationEnergy[0]=11.4732;


/*
  switch(ChargeEquilibrationMethod)
  {
    case RAPPE_GODDARD:
    case WILMER_SNURR:
      for(i=0;i<NUMBER_OF_ELECTRONEGATIVITIES_RG;i++)
      {
        index=ChargeEquilibrationDataRappeGoddard[i].index;
        ChargeEquilibrationData[index-1]=ChargeEquilibrationDataRappeGoddard[i];
      }
      break;
    case DING_YAZAYDIN_DUBBELDAM:
      break;
  }
*/

  index=0;
  NumberOfHydrogens = 0;
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(j=0;j<Framework[CurrentSystem].NumberOfAtoms[f1];j++)
    {
      type=Framework[CurrentSystem].Atoms[f1][j].Type;
      element_index=GetChargeEquilibrationIndexFromElementName(PseudoAtoms[type].ChemicalElement);
      covalent_radius=GetCovalentRadiusExtended(type,PseudoAtoms[type].ChemicalElement);

      Pos[index]=Framework[CurrentSystem].Atoms[f1][j].Position;

      J[index] = ReferenceTableJ(type,element_index);
      X[index] = ReferenceTableX(type,element_index);
      Xc[index] = ReferenceTableXc(type,element_index);

      switch(ChargeEquilibrationMethod)
      {
        case RAPPE_GODDARD:
          // TODO
          R[index] = gamma2*k*(1.0/J[index]);
          break;
        case WILMER_SNURR:
          R[index] = gamma2*k*(1.0/J[index]);
          break;
        case DING_YAZAYDIN_DUBBELDAM:
          R[index] = covalent_radius;
          break;
      }

      // make list of all the hydrogens
      if(element_index==0)
      {
        HydrogenList[NumberOfHydrogens]=index;
        NumberOfHydrogens++;
      }
      index++;
    }
  }

  // adjust hydrogen
  for(i=0;i<NumberOfHydrogens;i++)
  {
    element_index=HydrogenList[i];
    X[element_index]=0.5*(13.598 -2.0);
    J[element_index]=13.598+2.0;
  }

  if(NumberOfHydrogens>0)
  {
/*
    for(i=0;i<NumberOfHydrogens;i++)
    {
      element_index=HydrogenList[i];
      Q[element_index]=0.0;
    }
*/
/*
    for(j=0;j<10;j++) //iteration for H is 10 times
    {
      for(i=0;i<NumberOfHydrogens;i++)
      {
        element_index=HydrogenList[i];
        J[element_index]=(1.0+Q[element_index]/1.0698)*(ChargeEquilibrationData[0].IonizationEnergy[0]-ChargeEquilibrationData[0].ElectronAffinity);
      }
      Qeq(Q);
    }
*/
    Qeq(Q);
  }
  else
  {
    Qeq(Q);
  }

  //for(i=0;i<index;i++)
  //  fprintf(stderr, "%d %g %g charge: %g\n",i,X[i]-Xc[i],J[i],Q[i]);

  // compute net-charge
  net=0.0;
  for(i=0;i<size;i++)
    net+=Q[i];

  if(SymmetrizeFrameworkCharges)
  {
    for(i=0;i<NumberOfPseudoAtoms;i++)
    {
      index=0;
      // compute the average charge for this pseudo-atom
      for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
      {
        total_charge=0.0;
        count=0.0;
        for(j=0;j<Framework[CurrentSystem].NumberOfAtoms[f1];j++)
        {
          atom_type=Framework[CurrentSystem].Atoms[f1][j].Type;
          if(i==atom_type)
          {
            count+=1.0;
            total_charge+=Q[index];
          }
          index++;
        }
      }
      // set the average charge for this pseudo-atom
      for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
      {
        for(j=0;j<Framework[CurrentSystem].NumberOfAtoms[f1];j++)
        {
          atom_type=Framework[CurrentSystem].Atoms[f1][j].Type;
          if(i==atom_type)
          {
            Framework[CurrentSystem].Atoms[f1][j].Charge=total_charge/count;
            PseudoAtoms[i].Charge1=total_charge/count;
          }
        }
      }
    }
  }
  else
  {
    index=0;
    for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
    {
      for(j=0;j<Framework[CurrentSystem].NumberOfAtoms[f1];j++)
        Framework[CurrentSystem].Atoms[f1][j].Charge=Q[index++];
    }
  }

  free(Q);
  free(R);
  free(X);
  free(Xc);
  free(J);
}
