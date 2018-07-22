/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2017 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'molecule.c' is part of RASPA-2.0

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
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <math.h>
#include <signal.h>
#include "constants.h"
#include "molecule.h"
#include "framework.h"
#include "potentials.h"
#include "rigid.h"
#include "utils.h"
#include "cbmc.h"
#include "mc_moves.h"
#include "ewald.h"
#include "thermo_baro_stats.h"
#include "integration.h"
#include "scattering_factors.h"
#include "input.h"
#include "sample.h"
#include "cubic_spline_1d.h"
#include "output.h"

int CurrentComponent;

int MoleculeId=0;
int AtomId=0;

int WidomParticleInsertionComponent;

int NumberOfPseudoAtoms;
PSEUDO_ATOM *PseudoAtoms;
int **NumberOfPseudoAtomsCount;
int **NumberOfPseudoAtomsType;
int **NumberOfFractionalPseudoAtomsType;
int *NumberOfPseudoAtomsTypeNew;
int *NumberOfPseudoAtomsTypeOld;
int *MapPseudoAtom;

int *MaxNumberOfAdsorbateMolecules;
int *MaxNumberOfCationMolecules;
int *NumberOfFractionalMolecules;
int *NumberOfFractionalAdsorbateMolecules;
int *NumberOfFractionalCationMolecules;
int *NumberOfReactionMolecules;
int *NumberOfReactionAdsorbateMolecules;
int *NumberOfReactionCationMolecules;

int MaxNumberOfCoulombicSites;
int MaxNumberOfBondDipoleSites;
int LargestNumberOfCoulombicSites;
int LargestNumberOfBondDipoleSites;
int *NumberOfAtomsPerSystem;
int *NumberOfChargesPerSystem;
int *NumberOfBondDipolesPerSystem;

int *NumberOfAdsorbateMolecules;
int CurrentAdsorbateMolecule;
ADSORBATE_MOLECULE **Adsorbates;

int *NumberOfCationMolecules;
int CurrentCationMolecule;

CATION_MOLECULE **Cations;

int NumberOfComponents;
int NumberOfAdsorbateComponents;
int NumberOfCationComponents;
COMPONENT *Components;

REAL Derivatives[5];

int ReturnBondDipole(int comp,int A1,int A2)
{
  int i;
  int A,B;

  for(i=0;i<Components[comp].NumberOfBondDipoles;i++)
  {
    A=Components[comp].BondDipoles[i].A;
    B=Components[comp].BondDipoles[i].B;
    if(((A==A1)&&(B==A2))||((A==A2)&&(B==A1))) return i;
  }
  fprintf(stderr, "Dipole-pair is not present\n");
  return -1;
}


int ReturnPseudoAtomNumber(char *buffer)
{
  int j;

  for(j=0;j<NumberOfPseudoAtoms;j++)
    if(strncasecmp(PseudoAtoms[j].Name,buffer,MAX2(strlen(PseudoAtoms[j].Name),strlen(buffer)))==0) return j;

  fprintf(stderr, "ReturnPseudoAtomNumber: Error!!!! :%s\n",buffer);
  exit(0);
  return -1;
}

int ReturnPossiblePseudoAtomNumber(char *buffer)
{
  int j;

  for(j=0;j<NumberOfPseudoAtoms;j++)
    if(strncasecmp(PseudoAtoms[j].Name,buffer,MAX2(strlen(PseudoAtoms[j].Name),strlen(buffer)))==0) return j;

  return -1;
}


int CheckPseudoAtomNumber(char *buffer)
{
  int j;

  for(j=0;j<NumberOfPseudoAtoms;j++)
    if(strncasecmp(PseudoAtoms[j].Name,buffer,MAX2(strlen(PseudoAtoms[j].Name),strlen(buffer)))==0) return TRUE;
  return FALSE;
}

int AddPseudoAtom(PSEUDO_ATOM atom)
{
  int i,j;
  int AlreadyPresent;

  // first check if the pseudo_atom is already present
  AlreadyPresent=FALSE;
  for(j=0;j<NumberOfPseudoAtoms;j++)
  {
    if(strncasecmp(PseudoAtoms[j].Name,atom.Name,MAX2(strlen(PseudoAtoms[j].Name),strlen(atom.Name)))==0)
    {
      AlreadyPresent=TRUE;
      break;
    }
  }

  if(AlreadyPresent)
  {
    // for polarizable atoms which have zero charge, the inclusion is still needed to
    // compute the electric-field at that position
    if(fabs(atom.Charge1)>1e-10)
    {
      PseudoAtoms[j].HasCharges=TRUE;
    }

    return j;
  }
  else
  {
    PseudoAtoms=(PSEUDO_ATOM*)realloc(PseudoAtoms,(NumberOfPseudoAtoms+1)*sizeof(PSEUDO_ATOM));

    PseudoAtoms[NumberOfPseudoAtoms]=atom;

    NumberOfPseudoAtomsTypeNew=(int*)realloc(NumberOfPseudoAtomsTypeNew,(NumberOfPseudoAtoms+1)*sizeof(int));
    NumberOfPseudoAtomsTypeOld=(int*)realloc(NumberOfPseudoAtomsTypeOld,(NumberOfPseudoAtoms+1)*sizeof(int));

    MapPseudoAtom=(int*)realloc(MapPseudoAtom,(NumberOfPseudoAtoms+1)*sizeof(int));

    NumberOfPseudoAtomsTypeNew[NumberOfPseudoAtoms]=0;
    NumberOfPseudoAtomsTypeOld[NumberOfPseudoAtoms]=0;
    MapPseudoAtom[NumberOfPseudoAtoms]=0;

    for(i=0;i<NumberOfSystems;i++)
    {
      NumberOfPseudoAtomsCount[i]=(int*)realloc(NumberOfPseudoAtomsCount[i],(NumberOfPseudoAtoms+1)*sizeof(int));
      NumberOfPseudoAtomsType[i]=(int*)realloc(NumberOfPseudoAtomsType[i],(NumberOfPseudoAtoms+1)*sizeof(int));
      NumberOfFractionalPseudoAtomsType[i]=(int*)realloc(NumberOfFractionalPseudoAtomsType[i],(NumberOfPseudoAtoms+1)*sizeof(int));
      NumberOfPseudoAtomsCount[i][NumberOfPseudoAtoms]=0;
      NumberOfPseudoAtomsType[i][NumberOfPseudoAtoms]=0;
      NumberOfFractionalPseudoAtomsType[i][NumberOfPseudoAtoms]=0;
    }

    NumberOfPseudoAtoms++;
  }
  return NumberOfPseudoAtoms-1;
}


void CheckNumberOfMolecules()
{
  int k,l;
  int total;
  int measuredNumber;

  for(k=0;k<NumberOfSystems;k++)
  {
    total=NumberOfAdsorbateMolecules[k];

    measuredNumber=0;
    for(l=0;l<NumberOfComponents;l++)
    {
      measuredNumber+=Components[l].NumberOfMolecules[k];
    }
    if(measuredNumber!=total)
    {
      fprintf(stdout,"CheckNumberOfMolecules: Number of molecules (%d) not equal to the sum of the components (%d)\n", total, measuredNumber);
      __builtin_trap();
    }
  }
}

void CheckTypeOfMolecules()
{
  int i,k,l;
  int total;
  int measuredNumber;


  for(k=0;k<NumberOfSystems;k++)
  {
    for(l=0;l<NumberOfComponents;l++)
    {
      total=Components[l].NumberOfMolecules[k];

      measuredNumber=0;
      for(i=0;i<NumberOfAdsorbateMolecules[k];i++)
      {
        if(Adsorbates[k][i].Type==l) measuredNumber++;
      }
      if(measuredNumber!=total)
      {
        fprintf(stdout,"CheckTypeOfMolecules: Number of molecules of type %d (%d) not equal to the sum over all molecules (%d)\n",l, total, measuredNumber);
        __builtin_trap();
      }
    }
  }

}




int SelectRandomMoleculeOfType(int comp)
{
  int d,count;
  int CurrentMolecule;

  // choose a random molecule of this component
  d=(int)(RandomNumber()*(REAL)Components[comp].NumberOfMolecules[CurrentSystem]);

  count=-1;
  CurrentMolecule=-1;
  if(Components[comp].ExtraFrameworkMolecule)
  {
    do   // search for d-th molecule of the right type
      if(Cations[CurrentSystem][++CurrentMolecule].Type==comp) count++;
    while(d!=count);
  }
  else
  {
    do   // search for d-th molecule of the right type
      if(Adsorbates[CurrentSystem][++CurrentMolecule].Type==comp) count++;
    while(d!=count);
  }


  return CurrentMolecule;
}

int IsFractionalAdsorbateMolecule(int m)
{
  return ((IsFractionalCFMCAdsorbateMolecule(m))||(IsFractionalReactionAdsorbateMolecule(m)));
}

int IsFractionalCationMolecule(int m)
{
  return ((IsFractionalCFMCCationMolecule(m))||(IsFractionalReactionCationMolecule(m)));
}

int IsFractionalCFMCAdsorbateMolecule(int m)
{
  int i;

  for(i=0;i<NumberOfComponents;i++)
    if((!Components[i].ExtraFrameworkMolecule)&&(Components[i].FractionalMolecule[CurrentSystem]==m)) return TRUE;

  return FALSE;
}

int IsFractionalCFMCCationMolecule(int m)
{
  int i;

  for(i=0;i<NumberOfComponents;i++)
    if((Components[i].ExtraFrameworkMolecule)&&(Components[i].FractionalMolecule[CurrentSystem]==m)) return TRUE;

  return FALSE;
}

int IsFractionalReactionAdsorbateMolecule(int m)
{
  int k,l,type;

  type=Adsorbates[CurrentSystem][m].Type;

  for(k=0;k<NumberOfReactions;k++)
    for(l=0;l<ReactantsStoichiometry[k][type];l++)
      if(Components[type].ReactantFractionalMolecules[CurrentSystem][k][l]==m) return TRUE;

  for(k=0;k<NumberOfReactions;k++)
    for(l=0;l<ProductsStoichiometry[k][type];l++)
      if(Components[type].ProductFractionalMolecules[CurrentSystem][k][l]==m) return TRUE;

  return FALSE;
}

int IsFractionalReactionCationMolecule(int m)
{
  int k,l,type;

  type=Cations[CurrentSystem][m].Type;

  for(k=0;k<NumberOfReactions;k++)
    for(l=0;l<ReactantsStoichiometry[k][type];l++)
      if(Components[type].ReactantFractionalMolecules[CurrentSystem][k][l]==m) return TRUE;

  for(k=0;k<NumberOfReactions;k++)
    for(l=0;l<ProductsStoichiometry[k][type];l++)
      if(Components[type].ProductFractionalMolecules[CurrentSystem][k][l]==m) return TRUE;


  return FALSE;
}

int SelectRandomMoleculeOfTypeExcludingFractionalMolecule(int comp)
{
  int d,count;
  int CurrentMolecule;

  if(Components[comp].ExtraFrameworkMolecule)
  {
    // choose a random molecule of this component
    // the number of possibly selected molecules is reduced by the presence of
    // 1) the CFMC fractional molecule (at most 1)
    // 2) RXMC molecules
    d=(int)(RandomNumber()*(Components[comp].NumberOfMolecules[CurrentSystem]
            -(Components[comp].FractionalMolecule[CurrentSystem]>=0?1:0)
            -Components[comp].NumberOfRXMCMoleculesPresent[CurrentSystem]));

    count=-1;
    CurrentMolecule=-1;
    do   // search for d-th molecule of the right type
    {
      CurrentMolecule++;
      if((Cations[CurrentSystem][CurrentMolecule].Type==comp)&&
         (!IsFractionalCationMolecule(CurrentMolecule))&&
         (!IsFractionalReactionCationMolecule(CurrentMolecule))) count++;
    }
    while(d!=count);
  }
  else
  {
    // choose a random molecule of this component
    // the number of possibly selected molecules is reduced by the presence of
    // 1) the CFMC fractional molecule (at most 1)
    // 2) RXMC molecules
    d=(int)(RandomNumber()*(Components[comp].NumberOfMolecules[CurrentSystem]
            -(Components[comp].FractionalMolecule[CurrentSystem]>=0?1:0)
            -Components[comp].NumberOfRXMCMoleculesPresent[CurrentSystem]));

    count=-1;
    CurrentMolecule=-1;
    do   // search for d-th molecule of the right type
    {
      CurrentMolecule++;
      if((Adsorbates[CurrentSystem][CurrentMolecule].Type==comp)&&
         (!IsFractionalAdsorbateMolecule(CurrentMolecule))&&
         (!IsFractionalReactionAdsorbateMolecule(CurrentMolecule))) count++;
    }
    while(d!=count);
  }

  return CurrentMolecule;
}

int numberOfReactionMoleculesForComponent(int comp)
{
  int i,n;

  n=0;
  for(i=0;i<NumberOfReactions;i++)
  {
    n+=ReactantsStoichiometry[i][comp];
    n+=ProductsStoichiometry[i][comp];
  }
  return n;
}

int numberOfReactantMoleculesForComponent(int comp)
{
  int i,n;

  n=0;
  for(i=0;i<NumberOfReactions;i++)
    n+=ReactantsStoichiometry[i][comp];
  return n;
}

int numberOfProductMoleculesForComponent(int comp)
{
  int i,n;

  n=0;
  for(i=0;i<NumberOfReactions;i++)
    n+=ProductsStoichiometry[i][comp];
  return n;
}

void getListOfMoleculeIdentifiersForReactantsAndProduct(int comp,int *n,int *array)
{
  int i,k,index;

  index=0;
  for(i=0;i<NumberOfReactions;i++)
  {
    for(k=0;k<ReactantsStoichiometry[i][comp];k++)
      array[index++]=Components[comp].ReactantFractionalMolecules[CurrentSystem][i][k];
    for(k=0;k<ProductsStoichiometry[i][comp];k++)
      array[index++]=Components[comp].ProductFractionalMolecules[CurrentSystem][i][k];
  }
  *n=index;
}

void getListOfAllMoleculeIdentifiersForReactantsAndProduct(int *n,int *array)
{
  int i,j,k,index;

  index=0;
  for(i=0;i<NumberOfReactions;i++)
  {
    for(j=0;j<NumberOfComponents;j++)
    {
      for(k=0;k<ReactantsStoichiometry[i][j];k++)
        array[index++]=Components[j].ReactantFractionalMolecules[CurrentSystem][i][k];
      for(k=0;k<ProductsStoichiometry[i][j];k++)
        array[index++]=Components[j].ProductFractionalMolecules[CurrentSystem][i][k];
    }
  }
  *n=index;
}

int SelectRandomMoleculeOfTypeExcludingReactionMolecules(int reaction,int **LambdaRetraceMolecules)
{
  int i,j;
  int d,count;
  int numberOfSelectableMolecules;
  int numberOfSelectedMolecules;
  int numberOfExcludedMolecules;
  int listOfExcludedMolecules[256];
  int CurrentMolecule;

  // get a list of all molecules involved in reactions (they can not be selected)
  getListOfAllMoleculeIdentifiersForReactantsAndProduct(&numberOfExcludedMolecules,listOfExcludedMolecules);

  numberOfSelectedMolecules=0;
  for(i=0;i<NumberOfComponents;i++)
  {
    numberOfSelectableMolecules=Components[i].NumberOfMolecules[CurrentSystem]
                                -(Components[i].FractionalMolecule[CurrentSystem]>=0?1:0)
                                -numberOfReactionMoleculesForComponent(i);
    for(j=0;j<ReactantsStoichiometry[reaction][i];j++)
    {
      if(numberOfSelectableMolecules<=0) return -1;

      // choose a random molecule of this component
      d=(int)(RandomNumber()*numberOfSelectableMolecules);

      count=-1;
      CurrentMolecule=-1;
      do   // search for d-th molecule of the right type
      {
        CurrentMolecule++;
        if((Adsorbates[CurrentSystem][CurrentMolecule].Type==i)&&
           (!IsFractionalCFMCAdsorbateMolecule(CurrentMolecule))&&
           (!isInArrayOfSize(CurrentMolecule,numberOfExcludedMolecules,listOfExcludedMolecules))) count++;
      }
      while(d!=count);

      LambdaRetraceMolecules[i][j]=CurrentMolecule;
      numberOfSelectedMolecules++;

      // add the array with molecules that can not be selected
      listOfExcludedMolecules[numberOfExcludedMolecules++]=CurrentMolecule;
      numberOfSelectableMolecules--;
    }
  }

  return numberOfSelectedMolecules;
}

int SelectRandomMoleculeOfTypeExcludingProductMolecules(int reaction,int **LambdaRetraceMolecules)
{
  int i,j;
  int d,count;
  int numberOfSelectableMolecules;
  int numberOfSelectedMolecules;
  int numberOfExcludedMolecules;
  int listOfExcludedMolecules[256];
  int CurrentMolecule;

  // get a list of all molecules involved in reactions (they can not be selected)
  getListOfAllMoleculeIdentifiersForReactantsAndProduct(&numberOfExcludedMolecules,listOfExcludedMolecules);

  numberOfSelectedMolecules=0;
  for(i=0;i<NumberOfComponents;i++)
  {
    numberOfSelectableMolecules=Components[i].NumberOfMolecules[CurrentSystem]
                                -(Components[i].FractionalMolecule[CurrentSystem]>=0?1:0)
                                -numberOfReactionMoleculesForComponent(i);
    for(j=0;j<ProductsStoichiometry[reaction][i];j++)
    {
      if(numberOfSelectableMolecules<=0) return -1;

      // choose a random molecule of this component
      d=(int)(RandomNumber()*numberOfSelectableMolecules);

      count=-1;
      CurrentMolecule=-1;
      do   // search for d-th molecule of the right type
      {
        CurrentMolecule++;
        if((Adsorbates[CurrentSystem][CurrentMolecule].Type==i)&&
           (!IsFractionalCFMCAdsorbateMolecule(CurrentMolecule))&&
           (!isInArrayOfSize(CurrentMolecule,numberOfExcludedMolecules,listOfExcludedMolecules))) count++;
      }
      while(d!=count);

      LambdaRetraceMolecules[i][j]=CurrentMolecule;
      numberOfSelectedMolecules++;

      // add the array with molecules that can not be selected
      listOfExcludedMolecules[numberOfExcludedMolecules++]=CurrentMolecule;
      numberOfSelectableMolecules--;
    }
  }

  return numberOfSelectedMolecules;
}

// Read the definitions of the pseudo Atoms
void ReadPseudoAtomsDefinitions(void)
{
  int i,n,temp;
  char buffer[256];
  char buffer1[256];
  FILE *FilePtr;
  double temp1,temp2,temp3,temp4,temp5,temp6;
  char line[1024];
  char *arg_pointer;

  NumberOfPseudoAtomsCount=(int**)calloc(NumberOfSystems,sizeof(int*));
  NumberOfPseudoAtomsType=(int**)calloc(NumberOfSystems,sizeof(int*));
  NumberOfFractionalPseudoAtomsType=(int**)calloc(NumberOfSystems,sizeof(int*));

  sprintf(buffer,"./pseudo_atoms.def");
  if(!(FilePtr=fopen(buffer,"r")))
  {
    sprintf(buffer,"%s/share/raspa/forcefield/%s/pseudo_atoms.def",RASPA_DIRECTORY,ForceField);
    if(!(FilePtr=fopen(buffer,"r")))
    {
      fprintf(stderr, "'pseudo_atoms.def' file not found and therefore not used\n");
      return;
    }
  }

  ReadLine(line,1024,FilePtr); // skip line
  ReadLine(line,1024,FilePtr); // skip line
  sscanf(line,"%d\n",&temp);
  NumberOfPseudoAtoms=temp+1;

  PseudoAtoms=(PSEUDO_ATOM*)calloc(NumberOfPseudoAtoms,sizeof(PSEUDO_ATOM));

  NumberOfPseudoAtomsTypeNew=(int*)calloc(NumberOfPseudoAtoms,sizeof(int));
  NumberOfPseudoAtomsTypeOld=(int*)calloc(NumberOfPseudoAtoms,sizeof(int));
  MapPseudoAtom=(int*)calloc(NumberOfPseudoAtoms,sizeof(int));

  for(i=0;i<NumberOfSystems;i++)
  {
    NumberOfPseudoAtomsCount[i]=(int*)calloc(NumberOfPseudoAtoms,sizeof(int));
    NumberOfPseudoAtomsType[i]=(int*)calloc(NumberOfPseudoAtoms,sizeof(int));
    NumberOfFractionalPseudoAtomsType[i]=(int*)calloc(NumberOfPseudoAtoms,sizeof(int));
  }

  for(i=0;i<NumberOfPseudoAtoms;i++)
  {
    PseudoAtoms[i].FrameworkAtom=FALSE;  // set default to false, overwritten when frameworks are read
    PseudoAtoms[i].HasVDWInteraction=TRUE;
  }

  // one pseudo-atom that is always present is 'UNIT'
  strcpy(PseudoAtoms[0].Name,"UNIT");
  strcpy(PseudoAtoms[0].PrintToPDBName,"H");
  strcpy(PseudoAtoms[0].ChemicalElement,"H");
  strcpy(PseudoAtoms[0].OxidationStateString,"");
  PseudoAtoms[0].OxidationState=0.0;
  strcpy(PseudoAtoms[0].ScatteringSource,"H");
  PseudoAtoms[0].Occupancy=1.0;
  PseudoAtoms[0].FrameworkAtom=FALSE;
  PseudoAtoms[0].PrintToPDB=FALSE;
  PseudoAtoms[0].ScatteringType=0;
  PseudoAtoms[0].AnomalousScatteringType=0;
  PseudoAtoms[0].TemperatureFactor=0.0;
  PseudoAtoms[0].AnisotropicTemperatureFactor.ax=0.0;
  PseudoAtoms[0].AnisotropicTemperatureFactor.ay=0.0;
  PseudoAtoms[0].AnisotropicTemperatureFactor.az=0.0;
  PseudoAtoms[0].AnisotropicTemperatureFactor.bx=0.0;
  PseudoAtoms[0].AnisotropicTemperatureFactor.by=0.0;
  PseudoAtoms[0].AnisotropicTemperatureFactor.bz=0.0;
  PseudoAtoms[0].AnisotropicTemperatureFactor.cx=0.0;
  PseudoAtoms[0].AnisotropicTemperatureFactor.cy=0.0;
  PseudoAtoms[0].AnisotropicTemperatureFactor.cz=0.0;
  PseudoAtoms[0].ScatteringDispersionImaginary=0.0;
  PseudoAtoms[0].ScatteringDispersionImaginary=0.0;
  PseudoAtoms[0].Mass=1.0;
  PseudoAtoms[0].Charge1=1.0;
  PseudoAtoms[0].ChargeDefinitionType=CHARGE_ATOM_FROM_PSEUDO_ATOM_DEFINITION;
  PseudoAtoms[0].Polarization.ax=PseudoAtoms[0].Polarization.ay=PseudoAtoms[0].Polarization.az=0.0;
  PseudoAtoms[0].Polarization.bx=PseudoAtoms[0].Polarization.by=PseudoAtoms[0].Polarization.bz=0.0;
  PseudoAtoms[0].Polarization.cx=PseudoAtoms[0].Polarization.cy=PseudoAtoms[0].Polarization.cz=0.0;
  PseudoAtoms[0].HasCharges=TRUE;
  PseudoAtoms[0].IsPolarizable=FALSE;
  PseudoAtoms[0].Interaction=TRUE;
  PseudoAtoms[0].Radius=1.0;
  PseudoAtoms[0].Connectivity=0;
  PseudoAtoms[0].TinkerType=0;
  PseudoAtoms[0].AnisotropicCorrection=FALSE;
  PseudoAtoms[0].AnisotropicDisplacement=0.0;
  PseudoAtoms[0].AnisotropicType=0.0;
  PseudoAtoms[0].HasVDWInteraction=TRUE;
  PseudoAtoms[0].Hybridization=HYBRIDIZATION_UNINITIALIZED;
  PseudoAtoms[0].CF=FALSE;


  ReadLine(line,1024,FilePtr); // skip line

  for(i=1;i<NumberOfPseudoAtoms;i++)
  {
    temp1=temp2=temp3=temp4=temp5=temp6=0.0;
    PseudoAtoms[i].Connectivity=0;
    PseudoAtoms[i].AnisotropicDisplacement=0.0;

    strcpy(buffer1,"");
    ReadLine(line,1024,FilePtr);
    arg_pointer=line;

    sscanf(arg_pointer,"%s %s %s %s %d %lf %lf%n",
       (char*)&PseudoAtoms[i].Name,
       buffer,
       (char*)&PseudoAtoms[i].PrintToPDBName,
       (char*)&PseudoAtoms[i].ChemicalElement,
       &PseudoAtoms[i].OxidationState,
       &temp1,
       &temp2,
       &n);

    sprintf(PseudoAtoms[i].OxidationStateString,"%+d",PseudoAtoms[i].OxidationState);

    PseudoAtoms[i].Mass=(REAL)temp1;
    PseudoAtoms[i].Charge1=(REAL)temp2;
    arg_pointer+=n;

    switch(PolarizationMatrix)
    {
      case ISOTROPIC:
        sscanf(arg_pointer,"%lf%n",&temp1,&n);
        PseudoAtoms[i].Polarization.ax=(REAL)(temp1/COULOMBIC_CONVERSION_FACTOR);
        PseudoAtoms[i].Polarization.by=(REAL)(temp1/COULOMBIC_CONVERSION_FACTOR);
        PseudoAtoms[i].Polarization.cz=(REAL)(temp1/COULOMBIC_CONVERSION_FACTOR);
        break;
      case ANISOTROPIC:
        sscanf(arg_pointer,"%lf %lf %lf%n",&temp1,&temp2,&temp3,&n);
        PseudoAtoms[i].Polarization.ax=(REAL)(temp1/COULOMBIC_CONVERSION_FACTOR);
        PseudoAtoms[i].Polarization.by=(REAL)(temp2/COULOMBIC_CONVERSION_FACTOR);
        PseudoAtoms[i].Polarization.cz=(REAL)(temp3/COULOMBIC_CONVERSION_FACTOR);
        break;
      case REGULAR_UPPER_TRIANGLE:
        sscanf(arg_pointer,"%lf %lf %lf %lf %lf %lf%n",&temp1,&temp2,&temp3,&temp4,&temp5,&temp6,&n);
        PseudoAtoms[i].Polarization.ax=(REAL)(temp1/COULOMBIC_CONVERSION_FACTOR);
        PseudoAtoms[i].Polarization.ay=PseudoAtoms[i].Polarization.bx=(REAL)(temp2/COULOMBIC_CONVERSION_FACTOR);
        PseudoAtoms[i].Polarization.az=PseudoAtoms[i].Polarization.cx=(REAL)(temp3/COULOMBIC_CONVERSION_FACTOR);
        PseudoAtoms[i].Polarization.by=(REAL)(temp4/COULOMBIC_CONVERSION_FACTOR);
        PseudoAtoms[i].Polarization.bz=PseudoAtoms[i].Polarization.cy=(REAL)(temp5/COULOMBIC_CONVERSION_FACTOR);
        PseudoAtoms[i].Polarization.cz=(REAL)(temp6/COULOMBIC_CONVERSION_FACTOR);
        break;
    }
    arg_pointer+=n;
    sscanf(arg_pointer,"%lf %lf %d %lf %s %d",
       &temp4,
       &temp5,
       &PseudoAtoms[i].Connectivity,
       &PseudoAtoms[i].AnisotropicDisplacement,
       buffer1,
       &PseudoAtoms[i].TinkerType);
    PseudoAtoms[i].TemperatureFactor=(REAL)temp4;
    PseudoAtoms[i].Radius=(REAL)temp5;

    if(strcasecmp("Absolute",buffer1)==0) PseudoAtoms[i].AnisotropicType=ABSOLUTE;
    if(strcasecmp("Relative",buffer1)==0) PseudoAtoms[i].AnisotropicType=RELATIVE;

    if((fabs(PseudoAtoms[i].Polarization.ax)<1e-10)&&(fabs(PseudoAtoms[i].Polarization.by)<1e-10)&&(fabs(PseudoAtoms[i].Polarization.cz)<1e-10))
      PseudoAtoms[i].IsPolarizable=FALSE;
    else
      PseudoAtoms[i].IsPolarizable=TRUE;

    if(fabs(PseudoAtoms[i].AnisotropicDisplacement)<1e-10)
      PseudoAtoms[i].AnisotropicCorrection=FALSE;
    else
      PseudoAtoms[i].AnisotropicCorrection=TRUE;

    // for polarizable atoms which have zero charge, the inclusion is still needed to
    // compute the electric-field at that position
    if((fabs(PseudoAtoms[i].Charge1)<1e-10)&&(!PseudoAtoms[i].IsPolarizable))
      PseudoAtoms[i].HasCharges=FALSE;
    else
      PseudoAtoms[i].HasCharges=TRUE;

    if(strncasecmp(buffer,"yes",strlen("yes"))==0)
      PseudoAtoms[i].PrintToPDB=TRUE;
    else
      PseudoAtoms[i].PrintToPDB=FALSE;

    PseudoAtoms[i].ChargeDefinitionType=CHARGE_ATOM_FROM_PSEUDO_ATOM_DEFINITION;

    PseudoAtoms[i].CF=FALSE;
  }
  fclose(FilePtr);
}

void ComputeInertiaTensorGroups(int comp)
{
  int j,k,ill;
  REAL Mass,TotalMass,rotxyz,temp,rotall,rotlim;
  VECTOR com,pos;
  REAL_MATRIX3x3 eigenvectors;
  VECTOR eigenvalues,dr;
  int atom_nr;

  for(j=0;j<Components[comp].NumberOfGroups;j++)
  {
    TotalMass=0.0;
    com.x=com.y=com.z=0.0;

    for(k=0;k<Components[comp].Groups[j].NumberOfGroupAtoms;k++)
    {
      atom_nr=Components[comp].Groups[j].Atoms[k];
      Mass=PseudoAtoms[Components[comp].Type[atom_nr]].Mass;
      com.x+=Mass*Components[comp].Positions[atom_nr].x;
      com.y+=Mass*Components[comp].Positions[atom_nr].y;
      com.z+=Mass*Components[comp].Positions[atom_nr].z;
      TotalMass+=Mass;
    }
    com.x/=TotalMass;
    com.y/=TotalMass;
    com.z/=TotalMass;

    Components[comp].Groups[j].InertiaTensor.ax=0.0;
    Components[comp].Groups[j].InertiaTensor.bx=0.0;
    Components[comp].Groups[j].InertiaTensor.cx=0.0;
    Components[comp].Groups[j].InertiaTensor.ay=0.0;
    Components[comp].Groups[j].InertiaTensor.by=0.0;
    Components[comp].Groups[j].InertiaTensor.cy=0.0;
    Components[comp].Groups[j].InertiaTensor.az=0.0;
    Components[comp].Groups[j].InertiaTensor.bz=0.0;
    Components[comp].Groups[j].InertiaTensor.cz=0.0;
    for(k=0;k<Components[comp].Groups[j].NumberOfGroupAtoms;k++)
    {
      atom_nr=Components[comp].Groups[j].Atoms[k];
      dr.x=Components[comp].Positions[atom_nr].x-com.x;
      dr.y=Components[comp].Positions[atom_nr].y-com.y;
      dr.z=Components[comp].Positions[atom_nr].z-com.z;

      Mass=PseudoAtoms[Components[comp].Type[atom_nr]].Mass;
      Components[comp].Groups[j].InertiaTensor.ax+=Mass*dr.x*dr.x;
      Components[comp].Groups[j].InertiaTensor.bx+=Mass*dr.y*dr.x;
      Components[comp].Groups[j].InertiaTensor.cx+=Mass*dr.z*dr.x;
      Components[comp].Groups[j].InertiaTensor.ay+=Mass*dr.x*dr.y;
      Components[comp].Groups[j].InertiaTensor.by+=Mass*dr.y*dr.y;
      Components[comp].Groups[j].InertiaTensor.cy+=Mass*dr.z*dr.y;
      Components[comp].Groups[j].InertiaTensor.az+=Mass*dr.x*dr.z;
      Components[comp].Groups[j].InertiaTensor.bz+=Mass*dr.y*dr.z;
      Components[comp].Groups[j].InertiaTensor.cz+=Mass*dr.z*dr.z;
    }

    // the local body frame is taken to be that in which the rotational inertia tensor is diagonal
    EigenSystem3x3(Components[comp].Groups[j].InertiaTensor,&eigenvectors,&eigenvalues);

    Components[comp].Groups[j].InertiaVector.x=0.0;
    Components[comp].Groups[j].InertiaVector.y=0.0;
    Components[comp].Groups[j].InertiaVector.z=0.0;
    for(k=0;k<Components[comp].Groups[j].NumberOfGroupAtoms;k++)
    {
      atom_nr=Components[comp].Groups[j].Atoms[k];
      dr.x=Components[comp].Positions[atom_nr].x-com.x;
      dr.y=Components[comp].Positions[atom_nr].y-com.y;
      dr.z=Components[comp].Positions[atom_nr].z-com.z;
      pos.x=eigenvectors.ax*dr.x+eigenvectors.bx*dr.y+eigenvectors.cx*dr.z;
      pos.y=eigenvectors.ay*dr.x+eigenvectors.by*dr.y+eigenvectors.cy*dr.z;
      pos.z=eigenvectors.az*dr.x+eigenvectors.bz*dr.y+eigenvectors.cz*dr.z;
      if(fabs(pos.x)<1e-8) pos.x=0.0;
      if(fabs(pos.y)<1e-8) pos.y=0.0;
      if(fabs(pos.z)<1e-8) pos.z=0.0;
      Components[comp].Positions[atom_nr]=pos;

      atom_nr=Components[comp].Groups[j].Atoms[k];
      Mass=PseudoAtoms[Components[comp].Type[atom_nr]].Mass;
      Components[comp].Groups[j].InertiaVector.x+=Mass*(SQR(pos.y)+SQR(pos.z));
      Components[comp].Groups[j].InertiaVector.y+=Mass*(SQR(pos.x)+SQR(pos.z));
      Components[comp].Groups[j].InertiaVector.z+=Mass*(SQR(pos.x)+SQR(pos.y));
    }
    for(k=0;k<Components[comp].Groups[j].NumberOfPermanentDipoles;k++)
    {
      dr.x=Components[comp].Groups[j].PermanentDipolePositions[k].x-com.x;
      dr.y=Components[comp].Groups[j].PermanentDipolePositions[k].y-com.y;
      dr.z=Components[comp].Groups[j].PermanentDipolePositions[k].z-com.z;
      pos.x=eigenvectors.ax*dr.x+eigenvectors.bx*dr.y+eigenvectors.cx*dr.z;
      pos.y=eigenvectors.ay*dr.x+eigenvectors.by*dr.y+eigenvectors.cy*dr.z;
      pos.z=eigenvectors.az*dr.x+eigenvectors.bz*dr.y+eigenvectors.cz*dr.z;
      if(fabs(pos.x)<1e-8) pos.x=0.0;
      if(fabs(pos.y)<1e-8) pos.y=0.0;
      if(fabs(pos.z)<1e-8) pos.z=0.0;
      Components[comp].Groups[j].PermanentDipolePositions[k]=pos;

      dr.x=Components[comp].Groups[j].PermanentDipoles[k].x;
      dr.y=Components[comp].Groups[j].PermanentDipoles[k].y;
      dr.z=Components[comp].Groups[j].PermanentDipoles[k].z;
      pos.x=eigenvectors.ax*dr.x+eigenvectors.bx*dr.y+eigenvectors.cx*dr.z;
      pos.y=eigenvectors.ay*dr.x+eigenvectors.by*dr.y+eigenvectors.cy*dr.z;
      pos.z=eigenvectors.az*dr.x+eigenvectors.bz*dr.y+eigenvectors.cz*dr.z;
      if(fabs(pos.x)<1e-8) pos.x=0.0;
      if(fabs(pos.y)<1e-8) pos.y=0.0;
      if(fabs(pos.z)<1e-8) pos.z=0.0;
      Components[comp].Groups[j].PermanentDipoles[k].x=pos.x;
      Components[comp].Groups[j].PermanentDipoles[k].y=pos.y;
      Components[comp].Groups[j].PermanentDipoles[k].z=pos.z;
    }
    for(k=0;k<Components[comp].Groups[j].NumberOfPolarizabilities;k++)
    {
      dr.x=Components[comp].Groups[j].PolarizabilityPositions[k].x-com.x;
      dr.y=Components[comp].Groups[j].PolarizabilityPositions[k].y-com.y;
      dr.z=Components[comp].Groups[j].PolarizabilityPositions[k].z-com.z;
      pos.x=eigenvectors.ax*dr.x+eigenvectors.bx*dr.y+eigenvectors.cx*dr.z;
      pos.y=eigenvectors.ay*dr.x+eigenvectors.by*dr.y+eigenvectors.cy*dr.z;
      pos.z=eigenvectors.az*dr.x+eigenvectors.bz*dr.y+eigenvectors.cz*dr.z;
      if(fabs(pos.x)<1e-8) pos.x=0.0;
      if(fabs(pos.y)<1e-8) pos.y=0.0;
      if(fabs(pos.z)<1e-8) pos.z=0.0;
      Components[comp].Groups[j].PolarizabilityPositions[k]=pos;

      dr.x=Components[comp].Groups[j].Polarizabilites[k].ax;
      dr.y=Components[comp].Groups[j].Polarizabilites[k].ay;
      dr.z=Components[comp].Groups[j].Polarizabilites[k].az;
      pos.x=eigenvectors.ax*dr.x+eigenvectors.bx*dr.y+eigenvectors.cx*dr.z;
      pos.y=eigenvectors.ay*dr.x+eigenvectors.by*dr.y+eigenvectors.cy*dr.z;
      pos.z=eigenvectors.az*dr.x+eigenvectors.bz*dr.y+eigenvectors.cz*dr.z;
      if(fabs(pos.x)<1e-8) pos.x=0.0;
      if(fabs(pos.y)<1e-8) pos.y=0.0;
      if(fabs(pos.z)<1e-8) pos.z=0.0;
      Components[comp].Groups[j].Polarizabilites[k].ax=pos.x;
      Components[comp].Groups[j].Polarizabilites[k].ay=pos.y;
      Components[comp].Groups[j].Polarizabilites[k].az=pos.z;

      dr.x=Components[comp].Groups[j].Polarizabilites[k].bx;
      dr.y=Components[comp].Groups[j].Polarizabilites[k].by;
      dr.z=Components[comp].Groups[j].Polarizabilites[k].bz;
      pos.x=eigenvectors.ax*dr.x+eigenvectors.bx*dr.y+eigenvectors.cx*dr.z;
      pos.y=eigenvectors.ay*dr.x+eigenvectors.by*dr.y+eigenvectors.cy*dr.z;
      pos.z=eigenvectors.az*dr.x+eigenvectors.bz*dr.y+eigenvectors.cz*dr.z;
      if(fabs(pos.x)<1e-8) pos.x=0.0;
      if(fabs(pos.y)<1e-8) pos.y=0.0;
      if(fabs(pos.z)<1e-8) pos.z=0.0;
      Components[comp].Groups[j].Polarizabilites[k].bx=pos.x;
      Components[comp].Groups[j].Polarizabilites[k].by=pos.y;
      Components[comp].Groups[j].Polarizabilites[k].bz=pos.z;

      dr.x=Components[comp].Groups[j].Polarizabilites[k].cx;
      dr.y=Components[comp].Groups[j].Polarizabilites[k].cy;
      dr.z=Components[comp].Groups[j].Polarizabilites[k].cz;
      pos.x=eigenvectors.ax*dr.x+eigenvectors.bx*dr.y+eigenvectors.cx*dr.z;
      pos.y=eigenvectors.ay*dr.x+eigenvectors.by*dr.y+eigenvectors.cy*dr.z;
      pos.z=eigenvectors.az*dr.x+eigenvectors.bz*dr.y+eigenvectors.cz*dr.z;
      if(fabs(pos.x)<1e-8) pos.x=0.0;
      if(fabs(pos.y)<1e-8) pos.y=0.0;
      if(fabs(pos.z)<1e-8) pos.z=0.0;
      Components[comp].Groups[j].Polarizabilites[k].cx=pos.x;
      Components[comp].Groups[j].Polarizabilites[k].cy=pos.y;
      Components[comp].Groups[j].Polarizabilites[k].cz=pos.z;
    }

    // set axis system: Ixx >= Iyy >= Izz
    rotxyz=MAX3(Components[comp].Groups[j].InertiaVector.x,Components[comp].Groups[j].InertiaVector.y,
                Components[comp].Groups[j].InertiaVector.z);
    if(rotxyz>=Components[comp].Groups[j].InertiaVector.x)
    {
      if(Components[comp].Groups[j].InertiaVector.y>=rotxyz)
      {
        for(k=0;k<Components[comp].Groups[j].NumberOfGroupAtoms;k++)
        {
          atom_nr=Components[comp].Groups[j].Atoms[k];
          temp=Components[comp].Positions[atom_nr].x;
          Components[comp].Positions[atom_nr].x=Components[comp].Positions[atom_nr].y;
          Components[comp].Positions[atom_nr].y=-temp;
        }
        for(k=0;k<Components[comp].Groups[j].NumberOfPermanentDipoles;k++)
        {
          temp=Components[comp].Groups[j].PermanentDipolePositions[k].x;
          Components[comp].Groups[j].PermanentDipolePositions[k].x=Components[comp].Groups[j].PermanentDipolePositions[k].y;
          Components[comp].Groups[j].PermanentDipolePositions[k].y=-temp;

          temp=Components[comp].Groups[j].PermanentDipoles[k].x;
          Components[comp].Groups[j].PermanentDipoles[k].x=Components[comp].Groups[j].PermanentDipoles[k].y;
          Components[comp].Groups[j].PermanentDipoles[k].y=-temp;
        }
        for(k=0;k<Components[comp].Groups[j].NumberOfPolarizabilities;k++)
        {
          temp=Components[comp].Groups[j].PolarizabilityPositions[k].x;
          Components[comp].Groups[j].PolarizabilityPositions[k].x=Components[comp].Groups[j].PolarizabilityPositions[k].y;
          Components[comp].Groups[j].PolarizabilityPositions[k].y=-temp;

          temp=Components[comp].Groups[j].Polarizabilites[k].ax;
          Components[comp].Groups[j].Polarizabilites[k].ax=Components[comp].Groups[j].Polarizabilites[k].ay;
          Components[comp].Groups[j].Polarizabilites[k].ay=-temp;
          temp=Components[comp].Groups[j].Polarizabilites[k].bx;
          Components[comp].Groups[j].Polarizabilites[k].bx=Components[comp].Groups[j].Polarizabilites[k].by;
          Components[comp].Groups[j].Polarizabilites[k].by=-temp;
          temp=Components[comp].Groups[j].Polarizabilites[k].cx;
          Components[comp].Groups[j].Polarizabilites[k].cx=Components[comp].Groups[j].Polarizabilites[k].cy;
          Components[comp].Groups[j].Polarizabilites[k].cy=-temp;
        }

        Components[comp].Groups[j].InertiaVector.y=Components[comp].Groups[j].InertiaVector.x;
        Components[comp].Groups[j].InertiaVector.x=rotxyz;
      }
      else if(Components[comp].Groups[j].InertiaVector.z>=rotxyz)
      {
        for(k=0;k<Components[comp].Groups[j].NumberOfGroupAtoms;k++)
        {
          atom_nr=Components[comp].Groups[j].Atoms[k];
          temp=Components[comp].Positions[atom_nr].x;
          Components[comp].Positions[atom_nr].x=Components[comp].Positions[atom_nr].z;
          Components[comp].Positions[atom_nr].z=-temp;
        }
        for(k=0;k<Components[comp].Groups[j].NumberOfPermanentDipoles;k++)
        {
          temp=Components[comp].Groups[j].PermanentDipolePositions[k].x;
          Components[comp].Groups[j].PermanentDipolePositions[k].x=Components[comp].Groups[j].PermanentDipolePositions[k].z;
          Components[comp].Groups[j].PermanentDipolePositions[k].z=-temp;

          temp=Components[comp].Groups[j].PermanentDipoles[k].x;
          Components[comp].Groups[j].PermanentDipoles[k].x=Components[comp].Groups[j].PermanentDipoles[k].z;
          Components[comp].Groups[j].PermanentDipoles[k].z=-temp;
        }
        for(k=0;k<Components[comp].Groups[j].NumberOfPolarizabilities;k++)
        {
          temp=Components[comp].Groups[j].PolarizabilityPositions[k].x;
          Components[comp].Groups[j].PolarizabilityPositions[k].x=Components[comp].Groups[j].PolarizabilityPositions[k].z;
          Components[comp].Groups[j].PolarizabilityPositions[k].z=-temp;

          temp=Components[comp].Groups[j].Polarizabilites[k].ax;
          Components[comp].Groups[j].Polarizabilites[k].ax=Components[comp].Groups[j].Polarizabilites[k].az;
          Components[comp].Groups[j].Polarizabilites[k].az=-temp;
          temp=Components[comp].Groups[j].Polarizabilites[k].bx;
          Components[comp].Groups[j].Polarizabilites[k].bx=Components[comp].Groups[j].Polarizabilites[k].bz;
          Components[comp].Groups[j].Polarizabilites[k].bz=-temp;
          temp=Components[comp].Groups[j].Polarizabilites[k].cx;
          Components[comp].Groups[j].Polarizabilites[k].cx=Components[comp].Groups[j].Polarizabilites[k].cz;
          Components[comp].Groups[j].Polarizabilites[k].cz=-temp;
        }

        Components[comp].Groups[j].InertiaVector.z=Components[comp].Groups[j].InertiaVector.x;
        Components[comp].Groups[j].InertiaVector.x=rotxyz;
      }
    }
    if(Components[comp].Groups[j].InertiaVector.z>Components[comp].Groups[j].InertiaVector.y)
    {
      for(k=0;k<Components[comp].Groups[j].NumberOfGroupAtoms;k++)
      {
        atom_nr=Components[comp].Groups[j].Atoms[k];
        temp=Components[comp].Positions[atom_nr].y;
        Components[comp].Positions[atom_nr].y=Components[comp].Positions[atom_nr].z;
        Components[comp].Positions[atom_nr].z=-temp;
      }
      for(k=0;k<Components[comp].Groups[j].NumberOfPermanentDipoles;k++)
      {
        temp=Components[comp].Groups[j].PermanentDipolePositions[k].y;
        Components[comp].Groups[j].PermanentDipolePositions[k].y=Components[comp].Groups[j].PermanentDipolePositions[k].z;
        Components[comp].Groups[j].PermanentDipolePositions[k].z=-temp;

        temp=Components[comp].Groups[j].PermanentDipoles[k].y;
        Components[comp].Groups[j].PermanentDipoles[k].y=Components[comp].Groups[j].PermanentDipoles[k].z;
        Components[comp].Groups[j].PermanentDipoles[k].z=-temp;
      }
      for(k=0;k<Components[comp].Groups[j].NumberOfPolarizabilities;k++)
      {
        temp=Components[comp].Groups[j].PolarizabilityPositions[k].y;
        Components[comp].Groups[j].PolarizabilityPositions[k].y=Components[comp].Groups[j].PolarizabilityPositions[k].z;
        Components[comp].Groups[j].PolarizabilityPositions[k].z=-temp;

        temp=Components[comp].Groups[j].Polarizabilites[k].ay;
        Components[comp].Groups[j].Polarizabilites[k].ay=Components[comp].Groups[j].Polarizabilites[k].az;
        Components[comp].Groups[j].Polarizabilites[k].az=-temp;
        temp=Components[comp].Groups[j].Polarizabilites[k].by;
        Components[comp].Groups[j].Polarizabilites[k].by=Components[comp].Groups[j].Polarizabilites[k].bz;
        Components[comp].Groups[j].Polarizabilites[k].bz=-temp;
        temp=Components[comp].Groups[j].Polarizabilites[k].cy;
        Components[comp].Groups[j].Polarizabilites[k].cy=Components[comp].Groups[j].Polarizabilites[k].cz;
        Components[comp].Groups[j].Polarizabilites[k].cz=-temp;
      }

      temp=Components[comp].Groups[j].InertiaVector.z;
      Components[comp].Groups[j].InertiaVector.z=Components[comp].Groups[j].InertiaVector.y;
      Components[comp].Groups[j].InertiaVector.y=temp;
    }

    // remove singularities
    rotlim=MAX2((REAL)1.0e-2,Components[comp].Groups[j].InertiaVector.x+Components[comp].Groups[j].InertiaVector.y+
            Components[comp].Groups[j].InertiaVector.z)*1.0e-5;

    if(Components[comp].Groups[j].InertiaVector.x<rotlim)
      Components[comp].Groups[j].InverseInertiaVector.x=0.0;
    else
      Components[comp].Groups[j].InverseInertiaVector.x=1.0/Components[comp].Groups[j].InertiaVector.x;

    if(Components[comp].Groups[j].InertiaVector.y<rotlim)
      Components[comp].Groups[j].InverseInertiaVector.y=0.0;
    else
      Components[comp].Groups[j].InverseInertiaVector.y=1.0/Components[comp].Groups[j].InertiaVector.y;

    if(Components[comp].Groups[j].InertiaVector.z<rotlim)
      Components[comp].Groups[j].InverseInertiaVector.z=0.0;
    else
      Components[comp].Groups[j].InverseInertiaVector.z=1.0/Components[comp].Groups[j].InertiaVector.z;


    if((Components[comp].Groups[j].InertiaVector.x+Components[comp].Groups[j].InertiaVector.y+Components[comp].Groups[j].InertiaVector.z)>1.0e-5)
      rotall=Components[comp].Groups[j].InertiaVector.x+Components[comp].Groups[j].InertiaVector.y+Components[comp].Groups[j].InertiaVector.z;
    else
      rotall=1.0;

    // determine point-particle,linear-group or nonlinear-group
    ill=0;
    if(Components[comp].Groups[j].InertiaVector.x/rotall<1.0e-5) ill++;
    if(Components[comp].Groups[j].InertiaVector.y/rotall<1.0e-5) ill++;
    if(Components[comp].Groups[j].InertiaVector.z/rotall<1.0e-5) ill++;

    switch(ill)
    {
      case 0:
        Components[comp].Groups[j].Type=NONLINEAR_MOLECULE;
        break;
      case 1:
        Components[comp].Groups[j].Type=LINEAR_MOLECULE;
        break;
      case 2:
      default:
        Components[comp].Groups[j].Type=POINT_PARTICLE;
        break;
    }
  }

}

// Read the definition of the comp-th Component of Type "Name"
void ReadComponentDefinition(int comp)
{
  int i,j,k,n,nr;
  int A,B,C,D,temp;
  int A1,A2,B1,B2,C1,C2;
  REAL mass,total_mass;
  char LR[256],buffer[256],line[1024];
  char *arg_pointer;
  VECTOR pos;
  FILE *FilePtr;
  VECTOR dr;
  int atom_nr;
  double temp1,temp2,temp3;
  double temp4,temp5,temp6;
  double temp7,temp8,temp9;

  Components[comp].CriticalTemperature=0.0;
  Components[comp].CriticalPressure=0.0;
  Components[comp].AcentricFactor=0.0;

  Components[comp].NumberOfRigidAtoms=0;
  Components[comp].NumberOfFlexibleAtoms=0;

  if(Components[comp].ExtraFrameworkMolecule)
    NumberOfCationComponents++;
  else
    NumberOfAdsorbateComponents++;

  sprintf(buffer,"%s.def",Components[comp].Name);
  if(!(FilePtr=fopen(buffer,"r")))
  {
    sprintf(buffer,"%s/share/raspa/molecules/%s/%s.%s",
       RASPA_DIRECTORY,
       Components[comp].MoleculeDefinition,
       Components[comp].Name,
       "def");

    if(!(FilePtr=fopen(buffer,"r")))
    {
      fprintf(stderr, "Cannot open %s. Error: %s.\n", buffer, strerror(errno));
      exit(1);
    }
  }

  ReadLine(line,1024,FilePtr); // skip line

  ReadLine(line,1024,FilePtr);
  sscanf(line,"%lf",&temp1);
  Components[comp].CriticalTemperature=(REAL)temp1; // read critical temperature

  ReadLine(line,1024,FilePtr);
  sscanf(line,"%lf",&temp1);
  Components[comp].CriticalPressure=(REAL)temp1;    // read critical pressure

  ReadLine(line,1024,FilePtr);
  sscanf(line,"%lf\n",&temp1);
  Components[comp].AcentricFactor=(REAL)temp1;      // read acentric factor

  ReadLine(line,1024,FilePtr);  // skip line

  ReadLine(line,1024,FilePtr);
  sscanf(line,"%d",&Components[comp].NumberOfAtoms);    // read NumberOfAtoms

  if(Components[comp].NumberOfAtoms<=0)
  {
    fprintf(stderr, "Error:  Number of atoms per molecule (%d) smaller or equal than zero \n",
           Components[comp].NumberOfAtoms);
    exit(0);
  }

  // allocate charility-centers
  Components[comp].Chirality=(int*)calloc(Components[comp].NumberOfAtoms,sizeof(int));
  Components[comp].ChiralityType=(int*)calloc(Components[comp].NumberOfAtoms,sizeof(int));
  Components[comp].ChiralA=(int*)calloc(Components[comp].NumberOfAtoms,sizeof(int));
  Components[comp].ChiralB=(int*)calloc(Components[comp].NumberOfAtoms,sizeof(int));
  Components[comp].ChiralC=(int*)calloc(Components[comp].NumberOfAtoms,sizeof(int));
  Components[comp].ChiralD=(int*)calloc(Components[comp].NumberOfAtoms,sizeof(int));

  for(i=0;i<Components[comp].NumberOfAtoms;i++)
  {
    Components[comp].ChiralityType[i]=NO_CHARILITY;
    Components[comp].Chirality[i]=FALSE;
    Components[comp].ChiralA[i]=0;
    Components[comp].ChiralB[i]=0;
    Components[comp].ChiralC[i]=0;
    Components[comp].ChiralD[i]=0;
  }

  Components[comp].Fixed=(int*)calloc(Components[comp].NumberOfAtoms,sizeof(int));
  Components[comp].Type=(int*)calloc(Components[comp].NumberOfAtoms,sizeof(int));
  Components[comp].Charge=(REAL*)calloc(Components[comp].NumberOfAtoms,sizeof(REAL));
  Components[comp].Connectivity=(int*)calloc(Components[comp].NumberOfAtoms,sizeof(int));

  // allocate bond-connectivity matrix
  Components[comp].ConnectivityMatrix=(int**)calloc(Components[comp].NumberOfAtoms,sizeof(int*));
  for(i=0;i<Components[comp].NumberOfAtoms;i++)
    Components[comp].ConnectivityMatrix[i]=(int*)calloc(Components[comp].NumberOfAtoms,sizeof(int));

  // set the bond-connectivity to all FALSE at first
  for(i=0;i<Components[comp].NumberOfAtoms;i++)
    for(j=0;j<Components[comp].NumberOfAtoms;j++)
      Components[comp].ConnectivityMatrix[i][j]=FALSE;

  Components[comp].Positions=(VECTOR*)calloc(Components[comp].NumberOfAtoms,sizeof(VECTOR));
  Components[comp].group=(int*)calloc(Components[comp].NumberOfAtoms,sizeof(int));

  Components[comp].RMCMOL=(VECTOR*)calloc(Components[comp].NumberOfAtoms,sizeof(VECTOR));

  Components[comp].CpuTimeTranslationMove=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  Components[comp].CpuTimeRandomTranslationMove=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  Components[comp].CpuTimeRotationMove=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  Components[comp].CpuTimeRandomRotationMove=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  Components[comp].CpuTimePartialReinsertionMove=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  Components[comp].CpuTimeReinsertionMove=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  Components[comp].CpuTimeReinsertionInPlaceMove=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  Components[comp].CpuTimeReinsertionInPlaneMove=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  Components[comp].CpuTimeIdentityChangeMove=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  Components[comp].CpuTimeSwapMoveInsertion=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  Components[comp].CpuTimeSwapMoveDeletion=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  Components[comp].CpuTimeCFSwapLambdaMove=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  Components[comp].CpuTimeCBCFSwapLambdaMove=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  Components[comp].CpuTimeWidomMove=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  Components[comp].CpuTimeCFWidomLambdaMove=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  Components[comp].CpuTimeGibbsWidomMove=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  Components[comp].CpuTimeSurfaceAreaMove=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  Components[comp].CpuTimeGibbsChangeMove=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  Components[comp].CpuTimeCFGibbsChangeMove=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  Components[comp].CpuTimeCBCFGibbsChangeMove=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  Components[comp].CpuTimeGibbsIdentityChangeMove=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  Components[comp].CpuTimeExchangeFractionalParticleMove=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  Components[comp].CpuTimeCFGibbsSwapFractionalMoleculeToOtherBoxMove=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  Components[comp].CpuTimeCFGibbsLambdaChangeMove=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  Components[comp].CpuTimeCFGibbsFractionalToIntegerMove=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

  Components[comp].MaximumCBMCChangeBondLength=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  Components[comp].MaximumCBMCChangeBendAngle=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  Components[comp].MaximumCBMCRotationOnCone=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  Components[comp].CBMCChangeBondLengthAttempts=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  Components[comp].CBMCChangeBondLengthAccepted=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  Components[comp].CBMCChangeBendAngleAttempts=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  Components[comp].CBMCChangeBendAngleAccepted=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  Components[comp].CBMCRotationOnConeAttempts=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  Components[comp].CBMCRotationOnConeAccepted=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  Components[comp].TotalCBMCChangeBondLengthAttempts=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  Components[comp].TotalCBMCChangeBendAngleAttempts=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  Components[comp].TotalCBMCRotationOnConeAttempts=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  Components[comp].TotalCBMCChangeBondLengthAccepted=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  Components[comp].TotalCBMCChangeBendAngleAccepted=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
  Components[comp].TotalCBMCRotationOnConeAccepted=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));

  Components[comp].ReactantFractionalMolecules=(int***)calloc(NumberOfSystems,sizeof(int**));
  Components[comp].ProductFractionalMolecules=(int***)calloc(NumberOfSystems,sizeof(int**));

  for(i=0;i<NumberOfSystems;i++)
  {
    Components[comp].MaximumCBMCChangeBondLength[i]=(REAL*)calloc(Components[comp].NumberOfAtoms,sizeof(REAL));
    Components[comp].MaximumCBMCChangeBendAngle[i]=(REAL*)calloc(Components[comp].NumberOfAtoms,sizeof(REAL));
    Components[comp].MaximumCBMCRotationOnCone[i]=(REAL*)calloc(Components[comp].NumberOfAtoms,sizeof(REAL));
    Components[comp].CBMCChangeBondLengthAttempts[i]=(REAL*)calloc(Components[comp].NumberOfAtoms,sizeof(REAL));
    Components[comp].CBMCChangeBondLengthAccepted[i]=(REAL*)calloc(Components[comp].NumberOfAtoms,sizeof(REAL));
    Components[comp].CBMCChangeBendAngleAttempts[i]=(REAL*)calloc(Components[comp].NumberOfAtoms,sizeof(REAL));
    Components[comp].CBMCChangeBendAngleAccepted[i]=(REAL*)calloc(Components[comp].NumberOfAtoms,sizeof(REAL));
    Components[comp].CBMCRotationOnConeAttempts[i]=(REAL*)calloc(Components[comp].NumberOfAtoms,sizeof(REAL));
    Components[comp].CBMCRotationOnConeAccepted[i]=(REAL*)calloc(Components[comp].NumberOfAtoms,sizeof(REAL));
    Components[comp].TotalCBMCChangeBondLengthAttempts[i]=(REAL*)calloc(Components[comp].NumberOfAtoms,sizeof(REAL));
    Components[comp].TotalCBMCChangeBendAngleAttempts[i]=(REAL*)calloc(Components[comp].NumberOfAtoms,sizeof(REAL));
    Components[comp].TotalCBMCRotationOnConeAttempts[i]=(REAL*)calloc(Components[comp].NumberOfAtoms,sizeof(REAL));
    Components[comp].TotalCBMCChangeBondLengthAccepted[i]=(REAL*)calloc(Components[comp].NumberOfAtoms,sizeof(REAL));
    Components[comp].TotalCBMCChangeBendAngleAccepted[i]=(REAL*)calloc(Components[comp].NumberOfAtoms,sizeof(REAL));
    Components[comp].TotalCBMCRotationOnConeAccepted[i]=(REAL*)calloc(Components[comp].NumberOfAtoms,sizeof(REAL));

    if(NumberOfReactions>0)
    {
      Components[comp].ReactantFractionalMolecules[i]=(int**)calloc(NumberOfReactions,sizeof(int*));
      Components[comp].ProductFractionalMolecules[i]=(int**)calloc(NumberOfReactions,sizeof(int*));

      for(j=0;j<NumberOfReactions;j++)
      {
        if(ReactantsStoichiometry[j][comp]>0)
          Components[comp].ReactantFractionalMolecules[i][j]=(int*)calloc(ReactantsStoichiometry[j][comp],sizeof(int));
        if(ProductsStoichiometry[j][comp]>0)
          Components[comp].ProductFractionalMolecules[i][j]=(int*)calloc(ProductsStoichiometry[j][comp],sizeof(int));
      }
    }
  }

  ReadLine(line,1024,FilePtr);  // skip line

  ReadLine(line,1024,FilePtr);
  sscanf(line,"%d",&Components[comp].NumberOfGroups);    // read number of groups

  Components[comp].Groups=(GROUP_DEFINITION*)calloc(Components[comp].NumberOfGroups,sizeof(GROUP_DEFINITION));

  if(Components[comp].NumberOfGroups<=0)
  {
    fprintf(stderr, "Error:  Number of groups per molecule (%d) smaller or equal than zero \n",
           Components[comp].NumberOfGroups);
    exit(0);
  }

  total_mass=0.0;
  for(i=0;i<Components[comp].NumberOfGroups;i++)
  {
    ReadLine(line,1024,FilePtr);  // skip line


    ReadLine(line,1024,FilePtr);
    sscanf(line,"%s",buffer);
    if(strcasecmp("rigid",buffer)==0)
      Components[comp].Groups[i].Rigid=TRUE;
    else if(strcasecmp("flexible",buffer)==0)
      Components[comp].Groups[i].Rigid=FALSE;
    else
    {
      fprintf(stderr, "Error:  Unknown group type (should be 'Flexible' or 'Rigid')\n");
      exit(0);
    }

    ReadLine(line,1024,FilePtr);  // skip line

    Components[comp].Groups[i].NumberOfGroupAtoms=0;
    Components[comp].Groups[i].NumberOfPermanentDipoles=0;
    Components[comp].Groups[i].NumberOfPolarizabilities=0;
    ReadLine(line,1024,FilePtr);
    sscanf(line,"%d%d%d",
        &Components[comp].Groups[i].NumberOfGroupAtoms,
        &Components[comp].Groups[i].NumberOfPermanentDipoles,
        &Components[comp].Groups[i].NumberOfPolarizabilities);

    if(Components[comp].Groups[i].Rigid)
      Components[comp].NumberOfRigidAtoms+=Components[comp].Groups[i].NumberOfGroupAtoms;
    else
      Components[comp].NumberOfFlexibleAtoms+=Components[comp].Groups[i].NumberOfGroupAtoms;

    mass=0.0;
    if(Components[comp].Groups[i].NumberOfGroupAtoms>0)
    {
      ReadLine(line,1024,FilePtr);  // skip line
      Components[comp].Groups[i].Atoms=(int*)calloc(Components[comp].Groups[i].NumberOfGroupAtoms,sizeof(int));
      for(j=0;j<Components[comp].Groups[i].NumberOfGroupAtoms;j++)
      {
        temp1=temp2=temp3=0.0;
        if(Components[comp].Groups[i].Rigid)
        {
          ReadLine(line,1024,FilePtr);
          sscanf(line,"%d%s%lf%lf%lf",&temp,buffer,&temp1,&temp2,&temp3);
        }
        else
        {
          ReadLine(line,1024,FilePtr);
          sscanf(line,"%d%s",&temp,buffer);
        }
        pos.x=(REAL)temp1;
        pos.y=(REAL)temp2;
        pos.z=(REAL)temp3;

        k=ReturnPseudoAtomNumber(buffer);
        Components[comp].Type[temp]=k;
        mass+=PseudoAtoms[k].Mass;
        total_mass+=PseudoAtoms[k].Mass;
        atom_nr=Components[comp].Groups[i].Atoms[j];
        Components[comp].Positions[temp]=pos;
        Components[comp].group[temp]=i;
        Components[comp].Groups[i].Atoms[j]=temp;
      }
    }

    if(Components[comp].Groups[i].NumberOfPermanentDipoles>0)
    {
      ReadLine(line,1024,FilePtr);  // skip line
      Components[comp].Groups[i].PermanentDipolePositions=(VECTOR*)calloc(Components[comp].Groups[i].NumberOfPermanentDipoles,sizeof(VECTOR));
      Components[comp].Groups[i].PermanentDipoles=(VECTOR*)calloc(Components[comp].Groups[i].NumberOfPermanentDipoles,sizeof(VECTOR));
      for(j=0;j<Components[comp].Groups[i].NumberOfPermanentDipoles;j++)
      {
        temp1=temp2=temp3=0.0;
        temp4=temp5=temp6=0.0;
        temp7=temp8=temp9=0.0;
        ReadLine(line,1024,FilePtr);
        sscanf(line,"%lf%lf%lf%lf%lf%lf",&temp1,&temp2,&temp3,&temp4,&temp5,&temp6);
        Components[comp].Groups[i].PermanentDipolePositions[j].x=temp1;
        Components[comp].Groups[i].PermanentDipolePositions[j].y=temp2;
        Components[comp].Groups[i].PermanentDipolePositions[j].z=temp3;
        Components[comp].Groups[i].PermanentDipoles[j].x=temp4/DEBYE_CONVERSION_FACTOR;
        Components[comp].Groups[i].PermanentDipoles[j].y=temp5/DEBYE_CONVERSION_FACTOR;
        Components[comp].Groups[i].PermanentDipoles[j].z=temp6/DEBYE_CONVERSION_FACTOR;
      }
    }

    if(Components[comp].Groups[i].NumberOfPolarizabilities>0)
    {
      ReadLine(line,1024,FilePtr);  // skip line
      Components[comp].Groups[i].PolarizabilityPositions=(VECTOR*)calloc(Components[comp].Groups[i].NumberOfPolarizabilities,sizeof(VECTOR));
      Components[comp].Groups[i].Polarizabilites=(REAL_MATRIX3x3*)calloc(Components[comp].Groups[i].NumberOfPolarizabilities,sizeof(REAL_MATRIX3x3));
      for(j=0;j<Components[comp].Groups[i].NumberOfPolarizabilities;j++)
      {
        temp1=temp2=temp3=0.0;
        temp4=temp5=temp6=0.0;
        temp7=temp8=temp9=0.0;
        ReadLine(line,1024,FilePtr);
        sscanf(line,"%lf%lf%lf%lf%lf%lf%lf%lf%lf",&temp1,&temp2,&temp3,&temp4,&temp5,&temp6,&temp7,&temp8,&temp9);
        Components[comp].Groups[i].PolarizabilityPositions[j].x=temp1;
        Components[comp].Groups[i].PolarizabilityPositions[j].y=temp2;
        Components[comp].Groups[i].PolarizabilityPositions[j].z=temp3;
        Components[comp].Groups[i].Polarizabilites[j].ax=temp4/COULOMBIC_CONVERSION_FACTOR;
        Components[comp].Groups[i].Polarizabilites[j].ay=Components[comp].Groups[i].Polarizabilites[j].bx=temp5/COULOMBIC_CONVERSION_FACTOR;
        Components[comp].Groups[i].Polarizabilites[j].az=Components[comp].Groups[i].Polarizabilites[j].cx=temp6/COULOMBIC_CONVERSION_FACTOR;
        Components[comp].Groups[i].Polarizabilites[j].by=temp7/COULOMBIC_CONVERSION_FACTOR;
        Components[comp].Groups[i].Polarizabilites[j].bz=Components[comp].Groups[i].Polarizabilites[j].cy=temp8/COULOMBIC_CONVERSION_FACTOR;
        Components[comp].Groups[i].Polarizabilites[j].cz=temp9/COULOMBIC_CONVERSION_FACTOR;
      }
    }

  }
  Components[comp].Mass=total_mass;

  for(i=0;i<Components[comp].NumberOfGroups;i++)
  {
    Components[comp].Groups[i].Mass=0.0;
    for(j=0;j<Components[comp].Groups[i].NumberOfGroupAtoms;j++)
    {
      A=Components[comp].Groups[i].Atoms[j];
      Components[comp].Groups[i].Mass+=PseudoAtoms[Components[comp].Type[A]].Mass;
    }
  }


  ComputeInertiaTensorGroups(comp);

  // fill in charges
  for(i=0;i<Components[comp].NumberOfAtoms;i++)
    Components[comp].Charge[i]=PseudoAtoms[Components[comp].Type[i]].Charge1;

  // consider molecule uncharged when all atoms are uncharged
  Components[comp].HasCharges=FALSE;
  for(i=0;i<Components[comp].NumberOfAtoms;i++)
    Components[comp].HasCharges|=PseudoAtoms[Components[comp].Type[i]].HasCharges;

  Components[comp].NumberOfCharges=0;
  for(i=0;i<Components[comp].NumberOfAtoms;i++)
    if(PseudoAtoms[Components[comp].Type[i]].HasCharges) Components[comp].NumberOfCharges++;

  Components[comp].IsPolarizable=FALSE;
  for(i=0;i<Components[comp].NumberOfAtoms;i++)
    Components[comp].IsPolarizable|=PseudoAtoms[Components[comp].Type[i]].IsPolarizable;

  ReadLine(line,1024,FilePtr); // skip line

  ReadLine(line,1024,FilePtr);
  sscanf(line,"%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d\n",
       &Components[comp].NumberOfChiralityCenters,
       &Components[comp].NumberOfBonds,
       &Components[comp].NumberOfBondDipoles,
       &Components[comp].NumberOfBends,
       &Components[comp].NumberOfUreyBradleys,
       &Components[comp].NumberOfInversionBends,
       &Components[comp].NumberOfTorsions,
       &Components[comp].NumberOfImproperTorsions,
       //&Components[comp].NumberOfOutOfPlanes,
       &Components[comp].NumberOfBondBonds,
       &Components[comp].NumberOfBondBends,
       &Components[comp].NumberOfBendBends,
       &Components[comp].NumberOfBondTorsions,
       &Components[comp].NumberOfBendTorsions,
       &Components[comp].NumberOfIntraVDW,
       &Components[comp].NumberOfIntraChargeCharge,
       &Components[comp].NumberOfIntraChargeBondDipole,
       &Components[comp].NumberOfIntraBondDipoleBondDipole);

  // update the maximum settings for the components
  // this is used to allocate the proper amount of memory in 'cbmc.c'
  if(Components[comp].NumberOfAtoms>MaxNumberOfBeads) MaxNumberOfBeads=Components[comp].NumberOfAtoms;
  if(Components[comp].NumberOfBonds>MaxNumberOfBonds) MaxNumberOfBonds=Components[comp].NumberOfBonds;
  if(Components[comp].NumberOfBondDipoles>MaxNumberOfBondDipoles) MaxNumberOfBondDipoles=Components[comp].NumberOfBondDipoles;
  if(Components[comp].NumberOfUreyBradleys>MaxNumberOfUreyBradleys) MaxNumberOfUreyBradleys=Components[comp].NumberOfUreyBradleys;
  if(Components[comp].NumberOfBends>MaxNumberOfBends) MaxNumberOfBends=Components[comp].NumberOfBends;
  if(Components[comp].NumberOfBendBends>MaxNumberOfBendBends) MaxNumberOfBendBends=Components[comp].NumberOfBendBends;
  if(Components[comp].NumberOfInversionBends>MaxNumberOfInversionBends)
    MaxNumberOfInversionBends=Components[comp].NumberOfInversionBends;
  if(Components[comp].NumberOfTorsions>MaxNumberOfTorsions) MaxNumberOfTorsions=Components[comp].NumberOfTorsions;
  if(Components[comp].NumberOfImproperTorsions>MaxNumberOfImproperTorsions)
     MaxNumberOfImproperTorsions=Components[comp].NumberOfImproperTorsions;
  if(Components[comp].NumberOfOutOfPlanes>MaxNumberOfOutOfPlanes)
     MaxNumberOfOutOfPlanes=Components[comp].NumberOfOutOfPlanes;

  if(Components[comp].NumberOfBondBonds>MaxNumberOfBondBonds) MaxNumberOfBondBonds=Components[comp].NumberOfBondBonds;
  if(Components[comp].NumberOfBondBends>MaxNumberOfBondBends) MaxNumberOfBondBends=Components[comp].NumberOfBondBends;
  if(Components[comp].NumberOfBondTorsions>MaxNumberOfBondTorsions) MaxNumberOfBondTorsions=Components[comp].NumberOfBondTorsions;
  if(Components[comp].NumberOfBendTorsions>MaxNumberOfBendTorsions) MaxNumberOfBendTorsions=Components[comp].NumberOfBendTorsions;

  if(Components[comp].NumberOfIntraVDW>MaxNumberOfIntraVDW) MaxNumberOfIntraVDW=Components[comp].NumberOfIntraVDW;
  if(Components[comp].NumberOfIntraChargeCharge>MaxNumberOfIntraChargeCharge)
     MaxNumberOfIntraChargeCharge=Components[comp].NumberOfIntraChargeCharge;
  if(Components[comp].NumberOfIntraChargeBondDipole>MaxNumberOfIntraChargeBondDipole)
     MaxNumberOfIntraChargeBondDipole=Components[comp].NumberOfIntraChargeBondDipole;
  if(Components[comp].NumberOfIntraBondDipoleBondDipole>MaxNumberOfIntraBondDipoleBondDipole)
     MaxNumberOfIntraBondDipoleBondDipole=Components[comp].NumberOfIntraBondDipoleBondDipole;

  // allocate bonds
  Components[comp].Bonds=(PAIR*)calloc(Components[comp].NumberOfBonds,sizeof(PAIR));
  Components[comp].BondType=(int*)calloc(Components[comp].NumberOfBonds,sizeof(int));
  Components[comp].BondArguments=(REAL(*)[MAX_BOND_POTENTIAL_ARGUMENTS])calloc(Components[comp].NumberOfBonds,
                                  sizeof(REAL[MAX_BOND_POTENTIAL_ARGUMENTS]));

  // allocate bond-dipoles
  Components[comp].BondDipoles=(PAIR*)calloc(Components[comp].NumberOfBondDipoles,sizeof(PAIR));
  Components[comp].BondDipoleMagnitude=(REAL*)calloc(Components[comp].NumberOfBondDipoles,sizeof(REAL));


  // allocate bends
  Components[comp].Bends=(QUAD*)calloc(Components[comp].NumberOfBends,sizeof(QUAD));
  Components[comp].BendType=(int*)calloc(Components[comp].NumberOfBends,sizeof(int));
  Components[comp].BendArguments=(REAL(*)[MAX_BEND_POTENTIAL_ARGUMENTS])calloc(Components[comp].NumberOfBends,
                                  sizeof(REAL[MAX_BEND_POTENTIAL_ARGUMENTS]));

  // allocate urey-bradleys
  Components[comp].UreyBradleys=(TRIPLE*)calloc(Components[comp].NumberOfUreyBradleys,sizeof(TRIPLE));
  Components[comp].UreyBradleyType=(int*)calloc(Components[comp].NumberOfUreyBradleys,sizeof(int));
  Components[comp].UreyBradleyArguments=(REAL(*)[MAX_UREYBRADLEY_POTENTIAL_ARGUMENTS])calloc(Components[comp].NumberOfUreyBradleys,
                                  sizeof(REAL[MAX_UREYBRADLEY_POTENTIAL_ARGUMENTS]));

  // allocate inversion-bends
  Components[comp].InversionBends=(QUAD*)calloc(Components[comp].NumberOfInversionBends,sizeof(QUAD));
  Components[comp].InversionBendType=(int*)calloc(Components[comp].NumberOfInversionBends,sizeof(int));
  Components[comp].InversionBendArguments=(REAL(*)[MAX_INVERSION_BEND_POTENTIAL_ARGUMENTS])calloc(Components[comp].NumberOfInversionBends,
                                  sizeof(REAL[MAX_INVERSION_BEND_POTENTIAL_ARGUMENTS]));

  // allocate torsions
  Components[comp].Torsions=(QUAD*)calloc(Components[comp].NumberOfTorsions,sizeof(QUAD));
  Components[comp].TorsionType=(int*)calloc(Components[comp].NumberOfTorsions,sizeof(int));
  Components[comp].TorsionArguments=(REAL(*)[MAX_TORSION_POTENTIAL_ARGUMENTS])calloc(Components[comp].NumberOfTorsions,
                                  sizeof(REAL[MAX_TORSION_POTENTIAL_ARGUMENTS]));

  // allocate improper torsions
  Components[comp].ImproperTorsions=(QUAD*)calloc(Components[comp].NumberOfImproperTorsions,sizeof(QUAD));
  Components[comp].ImproperTorsionType=(int*)calloc(Components[comp].NumberOfImproperTorsions,sizeof(int));
  Components[comp].ImproperTorsionArguments=(REAL(*)[MAX_IMPROPER_TORSION_POTENTIAL_ARGUMENTS])calloc(Components[comp].NumberOfImproperTorsions,
                                  sizeof(REAL[MAX_IMPROPER_TORSION_POTENTIAL_ARGUMENTS]));

  // allocate out-of-plane distances
  Components[comp].OutOfPlanes=(QUAD*)calloc(Components[comp].NumberOfOutOfPlanes,sizeof(QUAD));
  Components[comp].OutOfPlaneType=(int*)calloc(Components[comp].NumberOfOutOfPlanes,sizeof(int));
  Components[comp].OutOfPlaneArguments=(REAL(*)[MAX_OUT_OF_PLANE_POTENTIAL_ARGUMENTS])calloc(Components[comp].NumberOfOutOfPlanes,
                                  sizeof(REAL[MAX_OUT_OF_PLANE_POTENTIAL_ARGUMENTS]));

  // allocate bond-bonds
  Components[comp].BondBonds=(TRIPLE*)calloc(Components[comp].NumberOfBondBonds,sizeof(TRIPLE));
  Components[comp].BondBondType=(int*)calloc(Components[comp].NumberOfBondBonds,sizeof(int));
  Components[comp].BondBondArguments=(REAL(*)[MAX_BOND_BOND_POTENTIAL_ARGUMENTS])calloc(Components[comp].NumberOfBondBonds,
                                  sizeof(REAL[MAX_BOND_BOND_POTENTIAL_ARGUMENTS]));

  // allocate bond-bends
  Components[comp].BondBends=(TRIPLE*)calloc(Components[comp].NumberOfBondBends,sizeof(TRIPLE));
  Components[comp].BondBendType=(int*)calloc(Components[comp].NumberOfBondBends,sizeof(int));
  Components[comp].BondBendArguments=(REAL(*)[MAX_BOND_BEND_POTENTIAL_ARGUMENTS])calloc(Components[comp].NumberOfBondBends,
                                  sizeof(REAL[MAX_BOND_BEND_POTENTIAL_ARGUMENTS]));

  // allocate bend-bends
  Components[comp].BendBends=(QUAD*)calloc(Components[comp].NumberOfBendBends,sizeof(QUAD));
  Components[comp].BendBendType=(int*)calloc(Components[comp].NumberOfBendBends,sizeof(int));
  Components[comp].BendBendArguments=(REAL(*)[MAX_BEND_BEND_POTENTIAL_ARGUMENTS])calloc(Components[comp].NumberOfBendBends,
                                  sizeof(REAL[MAX_BEND_BEND_POTENTIAL_ARGUMENTS]));

  // allocate stretch-torsions
  Components[comp].BondTorsions=(QUAD*)calloc(Components[comp].NumberOfBondTorsions,sizeof(QUAD));
  Components[comp].BondTorsionType=(int*)calloc(Components[comp].NumberOfBondTorsions,sizeof(int));
  Components[comp].BondTorsionArguments=(REAL(*)[MAX_BOND_TORSION_POTENTIAL_ARGUMENTS])calloc(Components[comp].NumberOfBondTorsions,
                                  sizeof(REAL[MAX_BOND_TORSION_POTENTIAL_ARGUMENTS]));

  // allocate bend-torsions
  Components[comp].BendTorsions=(QUAD*)calloc(Components[comp].NumberOfBendTorsions,sizeof(QUAD));
  Components[comp].BendTorsionType=(int*)calloc(Components[comp].NumberOfBendTorsions,sizeof(int));
  Components[comp].BendTorsionArguments=(REAL(*)[MAX_BEND_TORSION_POTENTIAL_ARGUMENTS])calloc(Components[comp].NumberOfBendTorsions,
                                  sizeof(REAL[MAX_BEND_TORSION_POTENTIAL_ARGUMENTS]));

  Components[comp].IntraVDW=(PAIR*)calloc(Components[comp].NumberOfIntraVDW,sizeof(PAIR));
  Components[comp].IntraChargeCharge=(PAIR*)calloc(Components[comp].NumberOfIntraChargeCharge,sizeof(PAIR));
  Components[comp].IntraChargeBondDipole=(PAIR*)calloc(Components[comp].NumberOfIntraChargeBondDipole,sizeof(PAIR));
  Components[comp].IntraBondDipoleBondDipole=(PAIR*)calloc(Components[comp].NumberOfIntraBondDipoleBondDipole,sizeof(PAIR));

  Components[comp].IntraVDWScaling=(REAL*)calloc(Components[comp].NumberOfIntraVDW,sizeof(REAL));
  Components[comp].IntraChargeChargeScaling=(REAL*)calloc(Components[comp].NumberOfIntraVDW,sizeof(REAL));

  // allocate excluded pairs
  Components[comp].ExcludedIntraChargeCharge=(PAIR*)calloc(SQR(Components[comp].NumberOfAtoms),sizeof(PAIR));
  Components[comp].ExcludedIntraChargeBondDipole=(PAIR*)calloc(SQR(Components[comp].NumberOfAtoms),sizeof(PAIR));
  Components[comp].ExcludedIntraBondDipoleBondDipole=(PAIR*)calloc(SQR(Components[comp].NumberOfAtoms),sizeof(PAIR));

  // reading Chirality centers
  if(Components[comp].NumberOfChiralityCenters>0)
  {
    ReadLine(line,1024,FilePtr); // skip line
    for(i=0;i<Components[comp].NumberOfChiralityCenters;i++)
    {
      ReadLine(line,1024,FilePtr);
      sscanf(line,"%d%d%d%d%s%s",&A,&B,&C,&D,LR,buffer);
      Components[comp].ChiralA[B]=A;
      Components[comp].ChiralB[B]=B;
      Components[comp].ChiralC[B]=C;
      Components[comp].ChiralD[B]=D;
      Components[comp].Chirality[B]=TRUE;
      if(strncasecmp("L",LR,1)==0)
        Components[comp].ChiralityType[B]=S_CHIRAL;
      else if(strncasecmp("S",LR,1)==0)
        Components[comp].ChiralityType[B]=S_CHIRAL;
      else if(strncasecmp("R",LR,1)==0)
        Components[comp].ChiralityType[B]=R_CHIRAL;
    }
  }

  // reading Bond-data
  if(Components[comp].NumberOfBonds>0)
  {
    ReadLine(line,1024,FilePtr); // skip line
    for(i=0;i<Components[comp].NumberOfBonds;i++)
    {
      ReadLine(line,1024,FilePtr);
      arg_pointer=line;
      sscanf(line,"%d%d%s%n",&A,&B,buffer,&n);

      Components[comp].Bonds[i].A=A;
      Components[comp].Bonds[i].B=B;
      Components[comp].Connectivity[A]++;
      Components[comp].Connectivity[B]++;
      Components[comp].ConnectivityMatrix[A][B]=TRUE;
      Components[comp].ConnectivityMatrix[B][A]=TRUE;

      for(j=0;j<NR_BOND_TYPES;j++)
        if(strncasecmp(BondTypes[j].Name,buffer,MAX2(strlen(BondTypes[j].Name),strlen(buffer)))==0)
          Components[comp].BondType[i]=j;

      for(j=0;j<BondTypes[Components[comp].BondType[i]].nr_args;j++)
      {
        arg_pointer+=n;
        sscanf(arg_pointer,"%lf%n",&temp1,&n);
        Components[comp].BondArguments[i][j]=(REAL)temp1;
      }

      switch(Components[comp].BondType[i])
      {
        case HARMONIC_BOND:
          // 0.5*p0*SQR(r-p1);
          // ===============================================
          // p_0/k_B [K/A^2]   force constant
          // p_1     [A]       reference bond distance
          Components[comp].BondArguments[i][0]/=ENERGY_TO_KELVIN;
          break;
        case CORE_SHELL_SPRING:
          // 0.5*p0*SQR(r);
          // ===============================================
          // p_0/k_B [K/A^2]   force constant
          Components[comp].BondArguments[i][0]/=ENERGY_TO_KELVIN;
          break;
        case MORSE_BOND:
          // p_0*[(1.0-{exp(-p_1*(r-p_2))})^2-1.0]
          // ===============================================
          // p_0/k_B [K]       force constant
          // p_1     [A^-1]    parameter
          // p_2     [A]       reference bond distance
          Components[comp].BondArguments[i][0]/=ENERGY_TO_KELVIN;
          break;
        case LJ_12_6_BOND:
          // A/r_ij^12-B/r_ij^6
          // ===============================================
          // p_0/k_B [K A^12]
          // p_1/k_B [K A^6]
          Components[comp].BondArguments[i][0]/=ENERGY_TO_KELVIN;
          Components[comp].BondArguments[i][1]/=ENERGY_TO_KELVIN;
          break;
        case LENNARD_JONES_BOND:
          // 4*p_0*((p_1/r)^12-(p_1/r)^6)
          // ===============================================
          // p_0/k_B [K]
          // p_1     [A]
          Components[comp].BondArguments[i][0]/=ENERGY_TO_KELVIN;
          break;
        case BUCKINGHAM_BOND:
          // p_0*exp(-p_1 r)-p_2/r^6
          // ===============================================
          // p_0/k_B [K]
          // p_1     [A^-1]
          // p_2/k_B [K A^6]
          Components[comp].BondArguments[i][0]/=ENERGY_TO_KELVIN;
          Components[comp].BondArguments[i][2]/=ENERGY_TO_KELVIN;
          break;
        case RESTRAINED_HARMONIC_BOND:
          // 0.5*p_0*(r-p_1)^2                   |r-p_1|<=p_2
          // 0.5*p_0*p_2^2+p_0*p_2*(|r-p_1|-p_2) |r-p_1|>p_2
          // ===============================================
          // p_0/k_B [K/A^2]
          // p_1     [A]
          // p_2     [A]
          Components[comp].BondArguments[i][0]/=ENERGY_TO_KELVIN;
          break;
        case QUARTIC_BOND:
          // (1/2)*p_0*(r-p_1)^2+(1/3)*p_2*(r-p_1)^3+(1/4)*p_3*(r-p_1)^4
          // ===========================================================
          // p_0/k_B [K/A^2]
          // p_1     [A]
          // p_2/k_B [K/A^3]
          // p_3/k_B [K/A^4]
          Components[comp].BondArguments[i][0]/=ENERGY_TO_KELVIN;
          Components[comp].BondArguments[i][2]/=ENERGY_TO_KELVIN;
          Components[comp].BondArguments[i][3]/=ENERGY_TO_KELVIN;
          break;
        case CFF_QUARTIC_BOND:
          // p_0*(r-p_1)^2+p_2*(r-p_1)^3+p_3*(r-p_1)^4
          // ===============================================
          // p_0/k_B [K/A^2]
          // p_1     [A]
          // p_2/k_B [K/A^3]
          // p_3/k_B [K/A^4]
          Components[comp].BondArguments[i][0]/=ENERGY_TO_KELVIN;
          Components[comp].BondArguments[i][2]/=ENERGY_TO_KELVIN;
          Components[comp].BondArguments[i][3]/=ENERGY_TO_KELVIN;
          break;
        case MM3_BOND:
          // p_0*(r-p_1)^2*(1.0-2.55*(r-p_1)+(7.0/12.0)*2.55^2*(r-p_1)^2)
          // ============================================================
          // p_0     [mdyne/A molecule]
          // p_1     [A]
          Components[comp].BondArguments[i][0]*=71.94*KCAL_PER_MOL_TO_ENERGY;
          break;
        case RIGID_BOND:
          A=Components[comp].Bonds[i].A;
          B=Components[comp].Bonds[i].B;
          dr.x=Components[comp].Positions[A].x-Components[comp].Positions[B].x;
          dr.y=Components[comp].Positions[A].y-Components[comp].Positions[B].y;
          dr.z=Components[comp].Positions[A].z-Components[comp].Positions[B].z;
          Components[comp].BondArguments[i][0]=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
          break;
        case FIXED_BOND:
          Components[comp].NumberOfConstraintBonds++;
          break;
        default:
          fprintf(stderr, "Undefined Bond potential in routine 'ReadComponentDefinition' ('molecule.c')\n");
          exit(0);
          break;
       }
    }
  }

  Components[comp].ConnectivityList=(int**)calloc(Components[comp].NumberOfAtoms,sizeof(int*));
  for(A=0;A<Components[comp].NumberOfAtoms;A++)
  {
    Components[comp].ConnectivityList[A]=(int*)calloc(Components[comp].Connectivity[A],sizeof(int));

    nr=0;
    for(B=0;B<Components[comp].NumberOfAtoms;B++)
      if(Components[comp].ConnectivityMatrix[A][B])
        Components[comp].ConnectivityList[A][nr++]=B;
  }

  // reading Bond-dipole data
  if(Components[comp].NumberOfBondDipoles>0)
  {
    ReadLine(line,1024,FilePtr); // skip line
    for(i=0;i<Components[comp].NumberOfBondDipoles;i++)
    {
      ReadLine(line,1024,FilePtr);
      sscanf(line,"%d%d%lf",&Components[comp].BondDipoles[i].A,&Components[comp].BondDipoles[i].B,
                            &Components[comp].BondDipoleMagnitude[i]);
      Components[comp].BondDipoleMagnitude[i]/=DEBYE_CONVERSION_FACTOR;
    }
  }

  // reading Bend-data
  if(Components[comp].NumberOfBends>0)
  {
    ReadLine(line,1024,FilePtr); // skip line
    for(i=0;i<Components[comp].NumberOfBends;i++)
    {
      ReadLine(line,1024,FilePtr);
      arg_pointer=line;
      sscanf(line,"%d%d%d%s%n",
             &Components[comp].Bends[i].A,
             &Components[comp].Bends[i].B,
             &Components[comp].Bends[i].C,
             buffer,&n);
      for(j=0;j<NR_BEND_TYPES;j++)
        if(strncasecmp(BendTypes[j].Name,buffer,MAX2(strlen(BendTypes[j].Name),strlen(buffer)))==0)
          Components[comp].BendType[i]=j;

      for(j=0;j<BendTypes[Components[comp].BendType[i]].nr_args;j++)
      {
        arg_pointer+=n;
        sscanf(arg_pointer,"%lf%n",&temp1,&n);
        Components[comp].BendArguments[i][j]=(REAL)temp1;
      }

      switch(Components[comp].BendType[i])
      {
        case HARMONIC_BEND:
          // (1/2)p_0*(theta-p_1)^2
          // ===============================================
          // p_0/k_B [K/rad^2]
          // p_1     [degrees]
          Components[comp].BendArguments[i][0]/=ENERGY_TO_KELVIN;
          Components[comp].BendArguments[i][1]*=DEG2RAD;
          break;
        case CORE_SHELL_BEND:
          // (1/2)p_0*(theta-p_1)^2
          // ===============================================
          // p_0/k_B [K/rad^2]
          // p_1     [degrees]
          Components[comp].BendArguments[i][0]/=ENERGY_TO_KELVIN;
          Components[comp].BendArguments[i][1]*=DEG2RAD;
          break;
        case QUARTIC_BEND:
          // (1/2)p_0*(theta-p_1)^2+(1/3)*p_2*(theta-p_1)^3+(1/4)*p_2*(theta-p_1)^4
          // ======================================================================
          // p_0/k_B [K/rad^2]
          // p_1     [degrees]
          // p_2/k_B [K/rad^3]
          // p_3/k_B [K/rad^4]
          Components[comp].BendArguments[i][0]/=ENERGY_TO_KELVIN;
          Components[comp].BendArguments[i][1]*=DEG2RAD;
          Components[comp].BendArguments[i][2]/=ENERGY_TO_KELVIN;
          Components[comp].BendArguments[i][3]/=ENERGY_TO_KELVIN;
          break;
        case CFF_QUARTIC_BEND:
          // p_0*(theta-p_1)^2+p_2*(theta-p_1)^3+p_3*(theta-p_1)^4
          // =====================================================
          // p_0/k_B [K/rad^2]
          // p_1     [degrees]
          // p_2/k_B [K/rad^3]
          // p_3/k_B [K/rad^4]
          Components[comp].BendArguments[i][0]/=ENERGY_TO_KELVIN;
          Components[comp].BendArguments[i][1]*=DEG2RAD;
          Components[comp].BendArguments[i][2]/=ENERGY_TO_KELVIN;
          Components[comp].BendArguments[i][3]/=ENERGY_TO_KELVIN;
          break;
        case HARMONIC_COSINE_BEND:
          // (1/2)*p_0*(cos(theta)-cos(p_1))^2
          // ===============================================
          // p_0/k_B [K]
          // p_1     [degrees]
          Components[comp].BendArguments[i][0]/=ENERGY_TO_KELVIN;
          Components[comp].BendArguments[i][1]=cos(Components[comp].BendArguments[i][1]*DEG2RAD);
          break;
        case COSINE_BEND:
          // p_0*(1+cos(p_1*theta-p_2))
          // ===============================================
          // p_0/k_B [K]
          // p_1     [-]
          // p_2     [degrees]
          Components[comp].BendArguments[i][0]/=ENERGY_TO_KELVIN;
          Components[comp].BendArguments[i][2]*=RAD2DEG;
          break;
        case TAFIPOLSKY_BEND:
          // 0.5*p_0*(1+cos(theta))*(1+cos(2*theta))
          // ===============================================
          // p_0/k_B [K]
          Components[comp].BendArguments[i][0]/=ENERGY_TO_KELVIN;
          break;
        case MM3_BEND:
        case MM3_IN_PLANE_BEND:
          // p_0*(theta-p_1)^2(1-0.014*(theta-p_1)+5.6e-5*(theta-p_1)^2-7e-7*(theta-p_1)^3+2.2e-8(theta-p_1)^4)
          // =================================================================================================
          // p_0/k_B [mdyne A/rad^2]
          // p_1     [degrees]
          Components[comp].BendArguments[i][0]*=0.021914*KCAL_PER_MOL_TO_ENERGY;
          Components[comp].BendArguments[i][1]*=RAD2DEG;
          break;
        case FIXED_BEND:
          Components[comp].BendArguments[i][0]*=DEG2RAD;
          Components[comp].NumberOfConstraintBends++;
          break;
        default:
          fprintf(stderr, "Undefined Bend potential in routine 'ReadComponentDefinition' ('molecule.c')\n");
          exit(0);
          break;
       }
    }
  }

  // reading UreyBradley-data
  if(Components[comp].NumberOfUreyBradleys>0)
  {
    ReadLine(line,1024,FilePtr); // skip line
    for(i=0;i<Components[comp].NumberOfUreyBradleys;i++)
    {
      ReadLine(line,1024,FilePtr);
      arg_pointer=line;
      sscanf(line,"%d%d%d%s%n",&Components[comp].UreyBradleys[i].A,&Components[comp].UreyBradleys[i].B,
                             &Components[comp].UreyBradleys[i].C,buffer,&n);
      for(j=0;j<NR_UREYBRADLEY_TYPES;j++)
        if(strncasecmp(UreyBradleyTypes[j].Name,buffer,MAX2(strlen(UreyBradleyTypes[j].Name),strlen(buffer)))==0)
          Components[comp].UreyBradleyType[i]=j;

      for(j=0;j<UreyBradleyTypes[Components[comp].UreyBradleyType[i]].nr_args;j++)
      {
        arg_pointer+=n;
        sscanf(arg_pointer,"%lf%n",&temp1,&n);
        Components[comp].UreyBradleyArguments[i][j]=(REAL)temp1;
      }

      switch(Components[comp].UreyBradleyType[i])
      {
        case HARMONIC_UREYBRADLEY:
          // 0.5*p0*SQR(r-p1);
          // ===============================================
          // p_0/k_B [K/A^2]   force constant
          // p_1     [A]       reference bond distance
          Components[comp].UreyBradleyArguments[i][0]/=ENERGY_TO_KELVIN;
          break;
        case MORSE_UREYBRADLEY:
          // p_0*[(1.0-{exp(-p_1*(r-p_2))})^2-1.0]
          // ===============================================
          // p_0/k_B [K]       force constant
          // p_1     [A^-1]    parameter
          // p_2     [A]       reference bond distance
          Components[comp].UreyBradleyArguments[i][0]/=ENERGY_TO_KELVIN;
          break;
        case LJ_12_6_UREYBRADLEY:
          // A/r_ij^12-B/r_ij^6
          // ===============================================
          // p_0/k_B [K A^12]
          // p_1/k_B [K A^6]
          Components[comp].UreyBradleyArguments[i][0]/=ENERGY_TO_KELVIN;
          Components[comp].UreyBradleyArguments[i][1]/=ENERGY_TO_KELVIN;
          break;
        case LENNARD_JONES_UREYBRADLEY:
          // 4*p_0*((p_1/r)^12-(p_1/r)^6)
          // ===============================================
          // p_0/k_B [K]
          // p_1     [A]
          Components[comp].UreyBradleyArguments[i][0]/=ENERGY_TO_KELVIN;
          break;
        case BUCKINGHAM_UREYBRADLEY:
          // p_0*exp(-p_1 r)-p_2/r^6
          // ===============================================
          // p_0/k_B [K]
          // p_1     [A^-1]
          // p_2/k_B [K A^6]
          Components[comp].UreyBradleyArguments[i][0]/=ENERGY_TO_KELVIN;
          Components[comp].UreyBradleyArguments[i][2]/=ENERGY_TO_KELVIN;
          break;
        case RESTRAINED_HARMONIC_UREYBRADLEY:
          // 0.5*p_0*(r-p_1)^2                   |r-p_1|<=p_2
          // 0.5*p_0*p_2^2+p_0*p_2*(|r-p_1|-p_2) |r-p_1|>p_2
          // ===============================================
          // p_0/k_B [K/A^2]
          // p_1     [A]
          // p_2     [A]
          Components[comp].UreyBradleyArguments[i][0]/=ENERGY_TO_KELVIN;
          break;
        case QUARTIC_UREYBRADLEY:
          // (1/2)*p_0*(r-p_1)^2+(1/3)*p_2*(r-p_1)^3+(1/4)*p_3*(r-p_1)^4
          // ===========================================================
          // p_0/k_B [K/A^2]
          // p_1     [A]
          // p_2/k_B [K/A^3]
          // p_3/k_B [K/A^4]
          Components[comp].UreyBradleyArguments[i][0]/=ENERGY_TO_KELVIN;
          Components[comp].UreyBradleyArguments[i][2]/=ENERGY_TO_KELVIN;
          Components[comp].UreyBradleyArguments[i][3]/=ENERGY_TO_KELVIN;
          break;
        case CFF_QUARTIC_UREYBRADLEY:
          // p_0*(r-p_1)^2+p_2*(r-p_1)^3+p_3*(r-p_1)^4
          // ===============================================
          // p_0/k_B [K/A^2]
          // p_1     [A]
          // p_2/k_B [K/A^3]
          // p_3/k_B [K/A^4]
          Components[comp].UreyBradleyArguments[i][0]/=ENERGY_TO_KELVIN;
          Components[comp].UreyBradleyArguments[i][2]/=ENERGY_TO_KELVIN;
          Components[comp].UreyBradleyArguments[i][3]/=ENERGY_TO_KELVIN;
          break;
        case MM3_UREYBRADLEY:
          // p_0*(r-p_1)^2*(1.0-2.55*(r-p_1)+(7.0/12.0)*2.55^2*(r-p_1)^2)
          // ============================================================
          // p_0     [mdyne/A molecule]
          // p_1     [A]
          Components[comp].UreyBradleyArguments[i][0]*=71.94*KCAL_PER_MOL_TO_ENERGY;
          break;
        case RIGID_UREYBRADLEY:
          A=Components[comp].UreyBradleys[i].A;
          B=Components[comp].UreyBradleys[i].B;
          dr.x=Components[comp].Positions[A].x-Components[comp].Positions[B].x;
          dr.y=Components[comp].Positions[A].y-Components[comp].Positions[B].y;
          dr.z=Components[comp].Positions[A].z-Components[comp].Positions[B].z;
          Components[comp].BondArguments[i][0]=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
          break;
        case FIXED_UREYBRADLEY:
          break;
        default:
          fprintf(stderr, "Undefined Urey-Bradley potential in routine 'ReadComponentDefinition' ('molecule.c')\n");
          exit(0);
          break;
       }
    }
  }

  // reading inversion Bend-data
  if(Components[comp].NumberOfInversionBends>0)
  {
    ReadLine(line,1024,FilePtr); // skip line
    for(i=0;i<Components[comp].NumberOfInversionBends;i++)
    {
      ReadLine(line,1024,FilePtr);
      arg_pointer=line;
      sscanf(line,"%d%d%d%d%s%n",
        &Components[comp].InversionBends[i].A,
        &Components[comp].InversionBends[i].B,
        &Components[comp].InversionBends[i].C,
        &Components[comp].InversionBends[i].D,
        buffer,&n);
      for(j=0;j<NR_INVERSION_BEND_TYPES;j++)
        if(strncasecmp(InversionBendTypes[j].Name,buffer,MAX2(strlen(InversionBendTypes[j].Name),strlen(buffer)))==0)
          Components[comp].InversionBendType[i]=j;

      for(j=0;j<InversionBendTypes[Components[comp].InversionBendType[i]].nr_args;j++)
      {
        arg_pointer+=n;
        sscanf(arg_pointer,"%lf%n",&temp1,&n);
        Components[comp].InversionBendArguments[i][j]=(REAL)temp1;
      }

      switch(Components[comp].InversionBendType[i])
      {
        case HARMONIC_INVERSION:
        case HARMONIC_INVERSION2:
          // (1/2)*p_0*(chi-p_1)^2
          // ===============================================
          // p_0/k_B [K/rad^2]
          // p_1     [degrees]
          Components[comp].InversionBendArguments[i][0]/=ENERGY_TO_KELVIN;
          Components[comp].InversionBendArguments[i][1]*=DEG2RAD;
          break;
        case HARMONIC_COSINE_INVERSION:
        case HARMONIC_COSINE_INVERSION2:
          // (1/2)*p_0*(cos(phi)-cos(p_1))^2
          // ===============================================
          // p_0/k_B [K]
          // p_1     [degrees]
          Components[comp].InversionBendArguments[i][0]/=ENERGY_TO_KELVIN;
          Components[comp].InversionBendArguments[i][1]=cos(DEG2RAD*Components[comp].InversionBendArguments[i][1]);
          break;
        case PLANAR_INVERSION:
        case PLANAR_INVERSION2:
          // (1/2)*p_0*(1-cos(phi))
          // ===============================================
          // p_0/k_B [K]
          Components[comp].InversionBendArguments[i][0]/=ENERGY_TO_KELVIN;
          break;
        case MM3_INVERSION:
          // p_0*(theta-p_1)^2(1-0.014*(theta-p_1)+5.6e-5*(theta-p_1)^2-7e-7*(theta-p_1)^3+2.2e-8(theta-p_1)^4)
          // =================================================================================================
          // p_0/k_B [mdyne A/rad^2]
          // p_1     [degrees]
          Components[comp].InversionBendArguments[i][0]*=0.02191418*KCAL_PER_MOL_TO_ENERGY;
          break;
        case FIXED_INVERSION_BEND:
          Components[comp].InversionBendArguments[i][0]*=DEG2RAD;
          Components[comp].NumberOfConstraintInversionBends++;
          break;
        default:
          fprintf(stderr, "Undefined Inversion potential in routine 'ReadComponentDefinition' ('molecule.c')\n");
          exit(0);
          break;
       }
    }
  }

  // reading Torsion-data
  if(Components[comp].NumberOfTorsions>0)
  {
    ReadLine(line,1024,FilePtr); // skip line
    for(i=0;i<Components[comp].NumberOfTorsions;i++)
    {
      ReadLine(line,1024,FilePtr);
      arg_pointer=line;
      sscanf(line,"%d%d%d%d%s%n",
             &Components[comp].Torsions[i].A,
             &Components[comp].Torsions[i].B,
             &Components[comp].Torsions[i].C,
             &Components[comp].Torsions[i].D,
             buffer,&n);
      for(j=0;j<NR_TORSION_TYPES;j++)
      {
        if(strncasecmp(TorsionTypes[j].Name,buffer,MAX2(strlen(TorsionTypes[j].Name),strlen(buffer)))==0)
          Components[comp].TorsionType[i]=j;
      }

      for(j=0;j<TorsionTypes[Components[comp].TorsionType[i]].nr_args;j++)
      {
        arg_pointer+=n;
        sscanf(arg_pointer,"%lf%n",&temp1,&n);
        Components[comp].TorsionArguments[i][j]=(REAL)temp1;
      }

      switch(Components[comp].TorsionType[i])
      {
        case HARMONIC_DIHEDRAL:
          // (1/2)*p_0*(phi-p_1)^2
          // ===============================================
          // p_0/k_B [K/rad^2]
          // p_1     [degrees]
          Components[comp].TorsionArguments[i][0]/=ENERGY_TO_KELVIN;
          Components[comp].TorsionArguments[i][1]*=DEG2RAD;
          break;
        case HARMONIC_COSINE_DIHEDRAL:
          // (1/2)*p_0*(cos(phi)-cos(p_1))^2
          // ===============================================
          // p_0/k_B [K]
          // p_1     [degrees]
          Components[comp].TorsionArguments[i][0]/=ENERGY_TO_KELVIN;
          Components[comp].TorsionArguments[i][1]=cos(Components[comp].TorsionArguments[i][1]*DEG2RAD);
          break;
        case THREE_COSINE_DIHEDRAL:
          // (1/2)*p_0*(1+cos(phi))+(1/2)*p_1*(1-cos(2*phi))+(1/2)*p_2*(1+cos(3*phi))
          // ========================================================================
          // p_0/k_B [K]
          // p_1/k_B [K]
          // p_2/k_B [K]
          Components[comp].TorsionArguments[i][0]/=ENERGY_TO_KELVIN;
          Components[comp].TorsionArguments[i][1]/=ENERGY_TO_KELVIN;
          Components[comp].TorsionArguments[i][2]/=ENERGY_TO_KELVIN;
          break;
        case MM3_DIHEDRAL:
          // (1/2)*p_0*(1+cos(phi))+(1/2)*p_1*(1-cos(2*phi))+(1/2)*p_2*(1+cos(3*phi))
          // ========================================================================
          // p_0     [kcal/mol]
          // p_1     [kcal/mol]
          // p_2     [kcal/mol]
          Components[comp].TorsionArguments[i][0]*=KCAL_PER_MOL_TO_ENERGY;
          Components[comp].TorsionArguments[i][1]*=KCAL_PER_MOL_TO_ENERGY;
          Components[comp].TorsionArguments[i][2]*=KCAL_PER_MOL_TO_ENERGY;
        case CVFF_BLOCKED_DIHEDRAL:
          // 
          // ========================================================================
          // p_0     [rad]
          // p_1     [K]
          // p_2     [-]
          // p_3     [rad]
          // p_4     [rad]
          Components[comp].TorsionArguments[i][1]/=ENERGY_TO_KELVIN;
          break;
        case CFF_DIHEDRAL:
          // p_0*(1-cos(phi))+p_1*(1-cos(2*phi))+p_2*(1-cos(3*phi))
          // ======================================================
          // p_0/k_B [K]
          // p_1/k_B [K]
          // p_2/k_B [K]
          Components[comp].TorsionArguments[i][0]/=ENERGY_TO_KELVIN;
          Components[comp].TorsionArguments[i][1]/=ENERGY_TO_KELVIN;
          Components[comp].TorsionArguments[i][2]/=ENERGY_TO_KELVIN;
          break;
        case CFF_DIHEDRAL2:
          // p_0*(1+cos(phi))+p_1*(1+cos(2*phi))+p_2*(1+cos(3*phi))
          // ======================================================
          // p_0/k_B [K]
          // p_1/k_B [K]
          // p_2/k_B [K]
          Components[comp].TorsionArguments[i][0]/=ENERGY_TO_KELVIN;
          Components[comp].TorsionArguments[i][1]/=ENERGY_TO_KELVIN;
          Components[comp].TorsionArguments[i][2]/=ENERGY_TO_KELVIN;
          break;
        case SIX_COSINE_DIHEDRAL:
          // Prod_i=0^5 p_i*cos(phi)^i
          // =========================
          // p_0/k_B [K]
          // p_1/k_B [K]
          // p_2/k_B [K]
          // p_3/k_B [K]
          // p_4/k_B [K]
          // p_5/k_B [K]
          Components[comp].TorsionArguments[i][0]/=ENERGY_TO_KELVIN;
          Components[comp].TorsionArguments[i][1]/=ENERGY_TO_KELVIN;
          Components[comp].TorsionArguments[i][2]/=ENERGY_TO_KELVIN;
          Components[comp].TorsionArguments[i][3]/=ENERGY_TO_KELVIN;
          Components[comp].TorsionArguments[i][4]/=ENERGY_TO_KELVIN;
          Components[comp].TorsionArguments[i][5]/=ENERGY_TO_KELVIN;
          break;
        case TRAPPE_DIHEDRAL:
          // p_0[0]+p_1*(1+cos(phi))+p_2*(1-cos(2*phi))+p_3*(1+cos(3*phi))
          // =============================================================
          // p_0/k_B [K]
          // p_1/k_B [K]
          // p_2/k_B [K]
          // p_3/k_B [K]
          Components[comp].TorsionArguments[i][0]/=ENERGY_TO_KELVIN;
          Components[comp].TorsionArguments[i][1]/=ENERGY_TO_KELVIN;
          Components[comp].TorsionArguments[i][2]/=ENERGY_TO_KELVIN;
          Components[comp].TorsionArguments[i][3]/=ENERGY_TO_KELVIN;
          break;
        case TRAPPE_DIHEDRAL_EXTENDED:
          // p_0[0]+p_1*cos(phi)+p_2*cos(2*phi)+p_3*cos(3*phi)+p_4*cos(4*phi)
          // ================================================================
          // p_0/k_B [K]
          // p_1/k_B [K]
          // p_2/k_B [K]
          // p_3/k_B [K]
          // p_4/k_B [K]
          Components[comp].TorsionArguments[i][0]/=ENERGY_TO_KELVIN;
          Components[comp].TorsionArguments[i][1]/=ENERGY_TO_KELVIN;
          Components[comp].TorsionArguments[i][2]/=ENERGY_TO_KELVIN;
          Components[comp].TorsionArguments[i][3]/=ENERGY_TO_KELVIN;
          Components[comp].TorsionArguments[i][4]/=ENERGY_TO_KELVIN;
          break;
        case MOD_TRAPPE_DIHEDRAL:
          /* Salvador modification: 16/08/2016
           add phase in cos function:
           p_0+p_1*(1+cos(phi-p_4))+p_2*(1-cos(2*(phi-p_4)))+p_3*(1+cos(3*(phi-p_4)))
          */
          Components[comp].TorsionArguments[i][0]/=ENERGY_TO_KELVIN;
          Components[comp].TorsionArguments[i][1]/=ENERGY_TO_KELVIN;
          Components[comp].TorsionArguments[i][2]/=ENERGY_TO_KELVIN;
          Components[comp].TorsionArguments[i][3]/=ENERGY_TO_KELVIN;
          Components[comp].TorsionArguments[i][4]*=DEG2RAD;
          break;
        case CVFF_DIHEDRAL:
          // p_0*(1+cos(p_1*phi-p_2))
          // ========================
          // p_0/k_B [K]
          // p_1     [-]
          // p_2     [degrees]
          Components[comp].TorsionArguments[i][0]/=ENERGY_TO_KELVIN;
          Components[comp].TorsionArguments[i][2]*=DEG2RAD;
          break;
        case OPLS_DIHEDRAL:
          // (1/2)p_0[0]+(1/2)p_1*(1+cos(phi))+(1/2)p_2*(1-cos(2*phi))+(1/2)p_3*(1+cos(3*phi))
          // =================================================================================
          // p_0/k_B [K]
          // p_1/k_B [K]
          // p_2/k_B [K]
          // p_3/k_B [K]
          Components[comp].TorsionArguments[i][0]/=ENERGY_TO_KELVIN;
          Components[comp].TorsionArguments[i][1]/=ENERGY_TO_KELVIN;
          Components[comp].TorsionArguments[i][2]/=ENERGY_TO_KELVIN;
          Components[comp].TorsionArguments[i][3]/=ENERGY_TO_KELVIN;
          break;
        case FOURIER_SERIES_DIHEDRAL:
          // (1/2)p_0*(1+cos(phi))+(1/2)p_1(1-cos(2*phi))+(1/2)*p2_2*(1+cos(3*phi))+
          // (1/2)p_3*(1-cos(4*phi))+(1/2)p_4*(1+cos(5*phi))+(1/2)p_5*(1+cos(6*phi))
          // =======================================================================
          // p_0/k_B [K]
          // p_1/k_B [K]
          // p_2/k_B [K]
          // p_3/k_B [K]
          // p_4/k_B [K]
          // p_5/k_B [K]
          Components[comp].TorsionArguments[i][0]/=ENERGY_TO_KELVIN;
          Components[comp].TorsionArguments[i][1]/=ENERGY_TO_KELVIN;
          Components[comp].TorsionArguments[i][2]/=ENERGY_TO_KELVIN;
          Components[comp].TorsionArguments[i][3]/=ENERGY_TO_KELVIN;
          Components[comp].TorsionArguments[i][4]/=ENERGY_TO_KELVIN;
          Components[comp].TorsionArguments[i][5]/=ENERGY_TO_KELVIN;
          break;
        case FOURIER_SERIES_DIHEDRAL2:
          // (1/2)p_0*(1+cos(phi))+(1/2)p_1(1-cos(2*phi))+(1/2)*p2_2*(1+cos(3*phi))+
          // (1/2)p_3*(1+cos(4*phi))+(1/2)p_4*(1+cos(5*phi))+(1/2)p_5*(1+cos(6*phi))
          // =======================================================================
          // p_0/k_B [K]
          // p_1/k_B [K]
          // p_2/k_B [K]
          // p_3/k_B [K]
          // p_4/k_B [K]
          // p_5/k_B [K]
          Components[comp].TorsionArguments[i][0]/=ENERGY_TO_KELVIN;
          Components[comp].TorsionArguments[i][1]/=ENERGY_TO_KELVIN;
          Components[comp].TorsionArguments[i][2]/=ENERGY_TO_KELVIN;
          Components[comp].TorsionArguments[i][3]/=ENERGY_TO_KELVIN;
          Components[comp].TorsionArguments[i][4]/=ENERGY_TO_KELVIN;
          Components[comp].TorsionArguments[i][5]/=ENERGY_TO_KELVIN;
          break;
        case FIXED_DIHEDRAL:
          Components[comp].TorsionArguments[i][0]*=DEG2RAD;
          Components[comp].NumberOfConstraintTorsions++;
          break;
        default:
          fprintf(stderr, "Undefined Torsion potential in routine 'ReadComponentDefinition' ('molecule.c')\n");
          exit(0);
          break;
       }
    }
  }

  // reading Improper Torsion-data
  if(Components[comp].NumberOfImproperTorsions>0)
  {
    ReadLine(line,1024,FilePtr); // skip line
    for(i=0;i<Components[comp].NumberOfImproperTorsions;i++)
    {
      ReadLine(line,1024,FilePtr);
      arg_pointer=line;
      sscanf(line,"%d%d%d%d%s%n",
             &Components[comp].ImproperTorsions[i].A,
             &Components[comp].ImproperTorsions[i].B,
             &Components[comp].ImproperTorsions[i].C,
             &Components[comp].ImproperTorsions[i].D,
             buffer,&n);
      for(j=0;j<NR_IMPROPER_TORSION_TYPES;j++)
      {
        if(strncasecmp(ImproperTorsionTypes[j].Name,buffer,MAX2(strlen(ImproperTorsionTypes[j].Name),strlen(buffer)))==0)
          Components[comp].ImproperTorsionType[i]=j;
      }

      for(j=0;j<ImproperTorsionTypes[Components[comp].ImproperTorsionType[i]].nr_args;j++)
      {
        arg_pointer+=n;
        sscanf(arg_pointer,"%lf%n",&temp1,&n);
        Components[comp].ImproperTorsionArguments[i][j]=(REAL)temp1;
      }

      switch(Components[comp].ImproperTorsionType[i])
      {
        case HARMONIC_IMPROPER_DIHEDRAL:
          // (1/2)*p_0*(phi-p_1)^2
          // ===============================================
          // p_0/k_B [K/rad^2]
          // p_1     [degrees]
          Components[comp].ImproperTorsionArguments[i][0]/=ENERGY_TO_KELVIN;
          Components[comp].ImproperTorsionArguments[i][1]*=DEG2RAD;
          break;
        case HARMONIC_COSINE_IMPROPER_DIHEDRAL:
          // (1/2)*p_0*(cos(phi)-cos(p_1))^2
          // ===============================================
          // p_0/k_B [K]
          // p_1     [degrees]
          Components[comp].ImproperTorsionArguments[i][0]/=ENERGY_TO_KELVIN;
          Components[comp].ImproperTorsionArguments[i][1]=cos(Components[comp].TorsionArguments[i][1]*DEG2RAD);
          break;
        case THREE_COSINE_IMPROPER_DIHEDRAL:
          // (1/2)*p_0*(1+cos(phi))+(1/2)*p_1*(1-cos(2*phi))+(1/2)*p_2*(1+cos(3*phi))
          // ========================================================================
          // p_0/k_B [K]
          // p_1/k_B [K]
          // p_2/k_B [K]
          Components[comp].ImproperTorsionArguments[i][0]/=ENERGY_TO_KELVIN;
          Components[comp].ImproperTorsionArguments[i][1]/=ENERGY_TO_KELVIN;
          Components[comp].ImproperTorsionArguments[i][2]/=ENERGY_TO_KELVIN;
          break;
        case MM3_IMPROPER_DIHEDRAL:
          // p_0*(1+cos(phi))+p_1*(1+cos(2*phi))+p_2*(1+cos(3*phi))
          // ======================================================
          // p_0/k_B [K]
          // p_1/k_B [K]
          // p_2/k_B [K]
          Components[comp].ImproperTorsionArguments[i][0]*=KCAL_PER_MOL_TO_ENERGY;
          Components[comp].ImproperTorsionArguments[i][1]*=KCAL_PER_MOL_TO_ENERGY;
          Components[comp].ImproperTorsionArguments[i][2]*=KCAL_PER_MOL_TO_ENERGY;
        case CFF_IMPROPER_DIHEDRAL:
          // p_0*(1-cos(phi))+p_1*(1-cos(2*phi))+p_2*(1-cos(3*phi))
          // ======================================================
          // p_0/k_B [K]
          // p_1/k_B [K]
          // p_2/k_B [K]
          Components[comp].ImproperTorsionArguments[i][0]/=ENERGY_TO_KELVIN;
          Components[comp].ImproperTorsionArguments[i][1]/=ENERGY_TO_KELVIN;
          Components[comp].ImproperTorsionArguments[i][2]/=ENERGY_TO_KELVIN;
          break;
        case CFF_IMPROPER_DIHEDRAL2:
          // p_0*(1+cos(phi))+p_1*(1+cos(2*phi))+p_2*(1+cos(3*phi))
          // ======================================================
          // p_0/k_B [K]
          // p_1/k_B [K]
          // p_2/k_B [K]
          Components[comp].ImproperTorsionArguments[i][0]/=ENERGY_TO_KELVIN;
          Components[comp].ImproperTorsionArguments[i][1]/=ENERGY_TO_KELVIN;
          Components[comp].ImproperTorsionArguments[i][2]/=ENERGY_TO_KELVIN;
          break;
        case SIX_COSINE_IMPROPER_DIHEDRAL:
          // Prod_i=0^5 p_i*cos(phi)^i
          // =========================
          // p_0/k_B [K]
          // p_1/k_B [K]
          // p_2/k_B [K]
          // p_3/k_B [K]
          // p_4/k_B [K]
          // p_5/k_B [K]
          Components[comp].ImproperTorsionArguments[i][0]/=ENERGY_TO_KELVIN;
          Components[comp].ImproperTorsionArguments[i][1]/=ENERGY_TO_KELVIN;
          Components[comp].ImproperTorsionArguments[i][2]/=ENERGY_TO_KELVIN;
          Components[comp].ImproperTorsionArguments[i][3]/=ENERGY_TO_KELVIN;
          Components[comp].ImproperTorsionArguments[i][4]/=ENERGY_TO_KELVIN;
          Components[comp].ImproperTorsionArguments[i][5]/=ENERGY_TO_KELVIN;
          break;
        case TRAPPE_IMPROPER_DIHEDRAL:
          // p_0[0]+p_1*(1+cos(phi))+p_2*(1-cos(2*phi))+p_3*(1+cos(3*phi))
          // =============================================================
          // p_0/k_B [K]
          // p_1/k_B [K]
          // p_2/k_B [K]
          // p_3/k_B [K]
          Components[comp].ImproperTorsionArguments[i][0]/=ENERGY_TO_KELVIN;
          Components[comp].ImproperTorsionArguments[i][1]/=ENERGY_TO_KELVIN;
          Components[comp].ImproperTorsionArguments[i][2]/=ENERGY_TO_KELVIN;
          Components[comp].ImproperTorsionArguments[i][3]/=ENERGY_TO_KELVIN;
          break;
        case TRAPPE_IMPROPER_DIHEDRAL_EXTENDED:
          // p_0[0]+p_1*cos(phi)+p_2*cos(2*phi)+p_3*cos(3*phi)+p_4*cos(4*phi)
          // ================================================================
          // p_0/k_B [K]
          // p_1/k_B [K]
          // p_2/k_B [K]
          // p_3/k_B [K]
          // p_4/k_B [K]
          Components[comp].ImproperTorsionArguments[i][0]/=ENERGY_TO_KELVIN;
          Components[comp].ImproperTorsionArguments[i][1]/=ENERGY_TO_KELVIN;
          Components[comp].ImproperTorsionArguments[i][2]/=ENERGY_TO_KELVIN;
          Components[comp].ImproperTorsionArguments[i][3]/=ENERGY_TO_KELVIN;
          Components[comp].ImproperTorsionArguments[i][4]/=ENERGY_TO_KELVIN;
          break;
        case CVFF_IMPROPER_DIHEDRAL:
          // p_0*(1+cos(p_1*phi-p_2))
          // ========================
          // p_0/k_B [K]
          // p_1     [-]
          // p_2     [degrees]
          Components[comp].ImproperTorsionArguments[i][0]/=ENERGY_TO_KELVIN;
          Components[comp].ImproperTorsionArguments[i][2]*=DEG2RAD;
          break;
        case OPLS_IMPROPER_DIHEDRAL:
          // (1/2)p_0[0]+(1/2)p_1*(1+cos(phi))+(1/2)p_2*(1-cos(2*phi))+(1/2)p_3*(1+cos(3*phi))
          // =================================================================================
          // p_0/k_B [K]
          // p_1/k_B [K]
          // p_2/k_B [K]
          // p_3/k_B [K]
          Components[comp].ImproperTorsionArguments[i][0]/=ENERGY_TO_KELVIN;
          Components[comp].ImproperTorsionArguments[i][1]/=ENERGY_TO_KELVIN;
          Components[comp].ImproperTorsionArguments[i][2]/=ENERGY_TO_KELVIN;
          Components[comp].ImproperTorsionArguments[i][3]/=ENERGY_TO_KELVIN;
          break;
        case FOURIER_SERIES_IMPROPER_DIHEDRAL:
          // (1/2)p_0*(1+cos(phi))+(1/2)p_1(1-cos(2*phi))+(1/2)*p2_2*(1+cos(3*phi))+
          // (1/2)p_3*(1-cos(4*phi))+(1/2)p_4*(1+cos(5*phi))+(1/2)p_5*(1+cos(6*phi))
          // =======================================================================
          // p_0/k_B [K]
          // p_1/k_B [K]
          // p_2/k_B [K]
          // p_3/k_B [K]
          // p_4/k_B [K]
          // p_5/k_B [K]
          Components[comp].ImproperTorsionArguments[i][0]/=ENERGY_TO_KELVIN;
          Components[comp].ImproperTorsionArguments[i][1]/=ENERGY_TO_KELVIN;
          Components[comp].ImproperTorsionArguments[i][2]/=ENERGY_TO_KELVIN;
          Components[comp].ImproperTorsionArguments[i][3]/=ENERGY_TO_KELVIN;
          Components[comp].ImproperTorsionArguments[i][4]/=ENERGY_TO_KELVIN;
          Components[comp].ImproperTorsionArguments[i][5]/=ENERGY_TO_KELVIN;
          break;
        case FOURIER_SERIES_IMPROPER_DIHEDRAL2:
          // (1/2)p_0*(1+cos(phi))+(1/2)p_1(1-cos(2*phi))+(1/2)*p2_2*(1+cos(3*phi))+
          // (1/2)p_3*(1+cos(4*phi))+(1/2)p_4*(1+cos(5*phi))+(1/2)p_5*(1+cos(6*phi))
          // =======================================================================
          // p_0/k_B [K]
          // p_1/k_B [K]
          // p_2/k_B [K]
          // p_3/k_B [K]
          // p_4/k_B [K]
          // p_5/k_B [K]
          Components[comp].ImproperTorsionArguments[i][0]/=ENERGY_TO_KELVIN;
          Components[comp].ImproperTorsionArguments[i][1]/=ENERGY_TO_KELVIN;
          Components[comp].ImproperTorsionArguments[i][2]/=ENERGY_TO_KELVIN;
          Components[comp].ImproperTorsionArguments[i][3]/=ENERGY_TO_KELVIN;
          Components[comp].ImproperTorsionArguments[i][4]/=ENERGY_TO_KELVIN;
          Components[comp].ImproperTorsionArguments[i][5]/=ENERGY_TO_KELVIN;
          break;
        case FIXED_IMPROPER_DIHEDRAL:
          Components[comp].ImproperTorsionArguments[i][0]*=DEG2RAD;
          Components[comp].NumberOfConstraintImproperTorsions++;
          break;
        default:
          fprintf(stderr, "Undefined Improper Torsion potential in routine 'ReadComponentDefinition' ('molecule.c')\n");
          exit(0);
          break;
      }
    }
  }

  // READ OUT-OF_PLANE TODO

  // reading Bond/Strech cross term-data
  if(Components[comp].NumberOfBondBonds>0)
  {
    ReadLine(line,1024,FilePtr); // skip line
    for(i=0;i<Components[comp].NumberOfBondBonds;i++)
    {
      ReadLine(line,1024,FilePtr);
      arg_pointer=line;
      sscanf(line,"%d%d%d%s%n",
              &Components[comp].BondBonds[i].A,
              &Components[comp].BondBonds[i].B,
              &Components[comp].BondBonds[i].C,
              buffer,&n);
      for(j=0;j<NR_BOND_BOND_TYPES;j++)
        if(strncasecmp(BondBondTypes[j].Name,buffer,MAX2(strlen(BondBondTypes[j].Name),strlen(buffer)))==0)
          Components[comp].BondBondType[i]=j;

      for(j=0;j<BondBondTypes[Components[comp].BondBondType[i]].nr_args;j++)
      {
        arg_pointer+=n;
        sscanf(arg_pointer,"%lf%n",&temp1,&n);
        Components[comp].BondBondArguments[i][j]=(REAL)temp1;
      }

      switch(Components[comp].BondBondType[i])
      {
        case CVFF_BOND_BOND_CROSS:
        case CFF_BOND_BOND_CROSS:
          // p_0*(rab-p_1)*(rbc-p_2)
          // =======================
          // p_0/k_B [K/A^2]
          // p_1     [A]
          // p_2     [A]
          Components[comp].BondBondArguments[i][0]/=ENERGY_TO_KELVIN;
          break;
        default:
          fprintf(stderr, "Undefined Bond-Bond potential in routine 'ReadComponentDefinition' ('molecule.c')\n");
          exit(0);
          break;
      }
    }
  }

  // reading Bond/Bend cross term-data
  if(Components[comp].NumberOfBondBends>0)
  {
    ReadLine(line,1024,FilePtr); // skip line
    for(i=0;i<Components[comp].NumberOfBondBends;i++)
    {
      ReadLine(line,1024,FilePtr);
      arg_pointer=line;
      sscanf(line,"%d%d%d%s%n",
             &Components[comp].BondBends[i].A,
             &Components[comp].BondBends[i].B,
             &Components[comp].BondBends[i].C,
             buffer,&n);
      for(j=0;j<NR_BOND_BEND_TYPES;j++)
        if(strncasecmp(BondBendTypes[j].Name,buffer,MAX2(strlen(BondBendTypes[j].Name),strlen(buffer)))==0)
          Components[comp].BondBendType[i]=j;

      for(j=0;j<BondBendTypes[Components[comp].BondBendType[i]].nr_args;j++)
      {
        arg_pointer+=n;
        sscanf(arg_pointer,"%lf%n",&temp1,&n);
        Components[comp].BondBendArguments[i][j]=(REAL)temp1;
      }

      switch(Components[comp].BondBendType[i])
      {
        case CVFF_BOND_BEND_CROSS:
        case CFF_BOND_BEND_CROSS:
          // (Theta-p_0)*(p_1*(rab-p_2)+p_3*(rbc-p_4))
          // =========================================
          // p_0     [degrees]
          // p_1/k_B [K/A/rad]
          // p_2     [A]
          // p_3/k_B [K/A/rad]
          // p_4     [A]
          Components[comp].BondBendArguments[i][0]*=DEG2RAD;
          Components[comp].BondBendArguments[i][1]/=ENERGY_TO_KELVIN;
          Components[comp].BondBendArguments[i][3]/=ENERGY_TO_KELVIN;
          break;
        case MM3_BOND_BEND_CROSS:
          // p_0*[(rab-p_1)+(rbc-p_2)]*(Theta-p_3)
          // =====================================
          // p_0     [mdyne/rad]
          // p_1     [A]
          // p_2     [A]
          // p_3     [degrees]
          Components[comp].BondBendArguments[i][0]*=2.51118*KCAL_PER_MOL_TO_ENERGY;
          Components[comp].BondBendArguments[i][3]*=DEG2RAD;
          break;
        case TRUNCATED_HARMONIC:
          // (1/2)*p_0*(Theta-p_1)^2*exp(-(pow(rab,8)+pow(rbc,8))/pow(p_2,8))
          // ================================================================
          // p_0/k_B [K/rad^2]
          // p_1     [degrees]
          // p_2     [A]
          Components[comp].BondBendArguments[i][0]/=ENERGY_TO_KELVIN;
          Components[comp].BondBendArguments[i][1]*=DEG2RAD;
          break;
        case SCREENED_HARMONIC:
          // (1/2)*p_0*(Theta-p_1)^2*exp(-(rab/p_2+rbc/p_3))
          // ===============================================
          // p_0/k_B [K/rad^2]
          // p_1     [degrees]
          // p_2     [A]
          // p_3     [A]
          Components[comp].BondBendArguments[i][0]/=ENERGY_TO_KELVIN;
          Components[comp].BondBendArguments[i][1]*=DEG2RAD;
          break;
        case SCREENED_VESSAL:
          // (p_0/(8.0*(Theta-PI)^2))*((p_1-PI)^2-(Theta-PI)^2)^2*exp(-(rab/p_2+rbc/p_3))
          // ============================================================================
          // p_0/k_B [K/rad^2]
          // p_1     [degrees]
          // p_2     [A]
          // p_3     [A]
          Components[comp].BondBendArguments[i][0]/=ENERGY_TO_KELVIN;
          Components[comp].BondBendArguments[i][1]*=DEG2RAD;
          break;
        case TRUNCATED_VESSAL:
          // p_0*[pow(theta,p_2)*(theta-p_1)^2*(theta+p_1-2.0*PI)^2-0.5*p_2*pow(PI,p_2-1.0)*(theta-p_1)^2*pow(PI-p_1,3)]
          //    *exp(-(pow(rab,8)+pow(rbc,8))/pow(p_3,8))
          // ============================================================================
          // p_0/k_B [K/rad^(4+p_2)]
          // p_1     [degrees]
          // p_2     [-]
          // p_3     [A]
          Components[comp].BondBendArguments[i][0]/=ENERGY_TO_KELVIN;
          Components[comp].BondBendArguments[i][1]*=DEG2RAD;
          break;
        default:
          fprintf(stderr, "Undefined Bond-Bend potential in routine 'ReadComponentDefinition' ('molecule.c')\n");
          exit(0);
          break;
      }
    }
  }

  // reading Bend/Bend cross term-data
  if(Components[comp].NumberOfBendBends>0)
  {
    ReadLine(line,1024,FilePtr); // skip line
    for(i=0;i<Components[comp].NumberOfBendBends;i++)
    {
      ReadLine(line,1024,FilePtr);
      arg_pointer=line;
      sscanf(line,"%d%d%d%d%s%n",
             &Components[comp].BendBends[i].A,
             &Components[comp].BendBends[i].B,
             &Components[comp].BendBends[i].C,
             &Components[comp].BendBends[i].D,
             buffer,&n);
      for(j=0;j<NR_BEND_BEND_TYPES;j++)
        if(strncasecmp(BendBendTypes[j].Name,buffer,MAX2(strlen(BendBendTypes[j].Name),strlen(buffer)))==0)
          Components[comp].BendBendType[i]=j;

      for(j=0;j<BendBendTypes[Components[comp].BendBendType[i]].nr_args;j++)
      {
        arg_pointer+=n;
        sscanf(arg_pointer,"%lf%n",&temp1,&n);
        Components[comp].BendBendArguments[i][j]=(REAL)temp1;
      }

      switch(Components[comp].BendBendType[i])
      {
        case CVFF_BEND_BEND_CROSS:
        case CFF_BEND_BEND_CROSS:
         // p_0*(Theta1-p_1)*(Theta2-p_2)
          // ===================================
          // p_0/k_B [K/rad^2)]
          // p_1     [degrees]
          // p_2     [degrees]
          Components[comp].BendBendArguments[i][0]/=ENERGY_TO_KELVIN;
          Components[comp].BendBendArguments[i][1]*=DEG2RAD;
        case MM3_BEND_BEND_CROSS:
          // -p_0*(Theta1-p_1)*(Theta2-p_2)
          // ===================================
          // p_0     [mdyne A/rad^2]
          // p_1     [degrees]
          // p_2     [degrees]
          Components[comp].BendBendArguments[i][0]*=0.02191418*KCAL_PER_MOL_TO_ENERGY;
          Components[comp].BendBendArguments[i][1]*=DEG2RAD;
          Components[comp].BendBendArguments[i][2]*=DEG2RAD;
          break;
        default:
          fprintf(stderr, "Undefined Bend-Bend potential in routine 'ReadComponentDefinition' ('molecule.c')\n");
          exit(0);
          break;
      }
    }
  }

  // reading Bond/Torsion cross term-data
  if(Components[comp].NumberOfBondTorsions>0)
  {
    ReadLine(line,1024,FilePtr); // skip line
    for(i=0;i<Components[comp].NumberOfBondTorsions;i++)
    {
      ReadLine(line,1024,FilePtr);
      arg_pointer=line;
      sscanf(line,"%d%d%d%d%s%n",
             &Components[comp].BondTorsions[i].A,
             &Components[comp].BondTorsions[i].B,
             &Components[comp].BondTorsions[i].C,
             &Components[comp].BondTorsions[i].D,
             buffer,&n);
      for(j=0;j<NR_BOND_TORSION_TYPES;j++)
        if(strncasecmp(BondTorsionTypes[j].Name,buffer,MAX2(strlen(BondTorsionTypes[j].Name),strlen(buffer)))==0)
          Components[comp].BondTorsionType[i]=j;

      for(j=0;j<BondTorsionTypes[Components[comp].BondTorsionType[i]].nr_args;j++)
      {
        arg_pointer+=n;
        sscanf(arg_pointer,"%lf%n",&temp1,&n);
        Components[comp].BondTorsionArguments[i][j]=(REAL)temp1;
      }

      switch(Components[comp].BondTorsionType[i])
      {
        case MM3_BOND_TORSION_CROSS:
          // (1/2)p_0(r-p_3)(1+cos(phi))+(1/2)p_1(r-p_3)(1+cos(2phi))+(1/2)p_2(r-p_3)(1+cos(3phi))
          // =====================================================================================
          // p_0     [kcal/A mole]
          // p_1     [kcal/A mole]
          // p_2     [kcal/A mole]
          // p_3     [A]
          Components[comp].BondTorsionArguments[i][0]*=KCAL_PER_MOL_TO_ENERGY;
          Components[comp].BondTorsionArguments[i][1]*=KCAL_PER_MOL_TO_ENERGY;
          Components[comp].BondTorsionArguments[i][2]*=KCAL_PER_MOL_TO_ENERGY;
          break;
        default:
          fprintf(stderr, "Undefined Bond-Torsion potential in routine 'ReadComponentDefinition' ('molecule.c')\n");
          exit(0);
          break;
      }
    }
  }

  // reading Bend/Torsion cross term-data
  if(Components[comp].NumberOfBendTorsions>0)
  {
    ReadLine(line,1024,FilePtr); // skip line
    for(i=0;i<Components[comp].NumberOfBendTorsions;i++)
    {
      ReadLine(line,1024,FilePtr);
      arg_pointer=line;
      sscanf(line,"%d%d%d%d%s%n",
             &Components[comp].BendTorsions[i].A,
             &Components[comp].BendTorsions[i].B,
             &Components[comp].BendTorsions[i].C,
             &Components[comp].BendTorsions[i].D,
             buffer,&n);
      for(j=0;j<NR_BEND_TORSION_TYPES;j++)
        if(strncasecmp(BendTorsionTypes[j].Name,buffer,MAX2(strlen(BendTorsionTypes[j].Name),strlen(buffer)))==0)
          Components[comp].BendTorsionType[i]=j;

      for(j=0;j<BendTorsionTypes[Components[comp].BendTorsionType[i]].nr_args;j++)
      {
        arg_pointer+=n;
        sscanf(arg_pointer,"%lf%n",&temp1,&n);
        Components[comp].BendTorsionArguments[i][j]=(REAL)temp1;
      }

      switch(Components[comp].BendTorsionType[i])
      {
        case SMOOTHED_DIHEDRAL:
          // S(Theta1)*[p_0(1+cos(p_1*Phi-p_2)]*S(Theta2)
          // ======================================================================================
          // p_0/k_B [K/rad^2]
          // p_1     [-]
          // p_2     [degrees]
          Components[comp].BendTorsionArguments[i][0]/=ENERGY_TO_KELVIN;
          Components[comp].BendTorsionArguments[i][2]*=DEG2RAD;
          break;
        case SMOOTHED_THREE_COSINE_DIHEDRAL:
          // S(Theta1)*[(1/2)*(1+cos(Phi))+(1/2)*p_1*(1-cos(2*Phi))+(1/2)*(1+cos(3*Phi))]*S(Theta2)
          // ======================================================================================
          // p_0/k_B [K]
          // p_1/k_B [K]
          // p_2/k_B [K]
          Components[comp].BendTorsionArguments[i][0]/=ENERGY_TO_KELVIN;
          Components[comp].BendTorsionArguments[i][1]/=ENERGY_TO_KELVIN;
          Components[comp].BendTorsionArguments[i][2]/=ENERGY_TO_KELVIN;
          break;
        case NICHOLAS_DIHEDRAL:
          // S(Theta1)*[(1/2)*(1+cos(Phi))+(1/2)*p_1*(1-cos(2*Phi))+(1/2)*(1+cos(3*Phi))]*S(Theta2)
          // ======================================================================================
          // p_0/k_B [K]
          // p_1/k_B [K]
          // p_2/k_B [K]
          Components[comp].BendTorsionArguments[i][0]/=ENERGY_TO_KELVIN;
          Components[comp].BendTorsionArguments[i][1]/=ENERGY_TO_KELVIN;
          Components[comp].BendTorsionArguments[i][2]/=ENERGY_TO_KELVIN;
          break;
        case SMOOTHED_CFF_DIHEDRAL:
         // S(Theta1)*[(1-cos(Phi))+p_1*(1-cos(2*Phi))+(1-cos(3*Phi))]*S(Theta2)
          // ======================================================================================
          // p_0/k_B [K]
          // p_1/k_B [K]
          // p_2/k_B [K]
          Components[comp].BendTorsionArguments[i][0]/=ENERGY_TO_KELVIN;
          Components[comp].BendTorsionArguments[i][1]/=ENERGY_TO_KELVIN;
          Components[comp].BendTorsionArguments[i][2]/=ENERGY_TO_KELVIN;
          break;
        case SMOOTHED_CFF_DIHEDRAL2:
          Components[comp].BendTorsionArguments[i][0]/=ENERGY_TO_KELVIN;
          Components[comp].BendTorsionArguments[i][1]/=ENERGY_TO_KELVIN;
          Components[comp].BendTorsionArguments[i][2]/=ENERGY_TO_KELVIN;
          break;
        case CVFF_BEND_TORSION_CROSS:
        case CFF_BEND_TORSION_CROSS:
          // p_0*(Theta1-p_1)*(Theta2-p_2)*cos(Phi)
          // =====================================================================================
          // p_0/k_B [K/rad^3]
          // p_1     [degrees]
          // p_2     [degrees]
          Components[comp].BendTorsionArguments[i][0]/=ENERGY_TO_KELVIN;
          Components[comp].BendTorsionArguments[i][1]*=DEG2RAD;
          Components[comp].BendTorsionArguments[i][2]*=DEG2RAD;
          break;
        case SMOOTHED_CFF_BEND_TORSION_CROSS:
          // S(Theta1)*[p_0*(Theta1-p_1)*(Theta2-p_2)*cos(Phi)]*S(Theta2)
          // ======================================================================================
          // p_0/k_B [K/rad^3]
          // p_1     [degrees]
          // p_2     [degrees]
          Components[comp].BendTorsionArguments[i][0]/=ENERGY_TO_KELVIN;
          Components[comp].BendTorsionArguments[i][1]*=DEG2RAD;
          Components[comp].BendTorsionArguments[i][2]*=DEG2RAD;
          break;
        default:
          fprintf(stderr, "Undefined Bend-Torsion potential in routine 'ReadComponentDefinition' ('molecule.c')\n");
          exit(0);
          break;
      }
    }
  }

  if(Components[comp].NumberOfIntraVDW>0)
  {
    ReadLine(line,1024,FilePtr); // skip line
    for(i=0;i<Components[comp].NumberOfIntraVDW;i++)
    {
      ReadLine(line,1024,FilePtr);
      sscanf(line,"%d%d",
         &Components[comp].IntraVDW[i].A,
         &Components[comp].IntraVDW[i].B);
      Components[comp].IntraVDWScaling[i]=1.0;
    }
  }

  if(Components[comp].NumberOfIntraChargeCharge>0)
  {
    ReadLine(line,1024,FilePtr); // skip line
    for(i=0;i<Components[comp].NumberOfIntraChargeCharge;i++)
    {
      ReadLine(line,1024,FilePtr);
      sscanf(line,"%d%d",
         &Components[comp].IntraChargeCharge[i].A,
         &Components[comp].IntraChargeCharge[i].B);
      Components[comp].IntraChargeChargeScaling[i]=1.0;
    }
  }

  // compute exclusions for the Ewald-summation
  Components[comp].NumberOfExcludedIntraChargeCharge=0;
  for(i=0;i<Components[comp].NumberOfAtoms-1;i++)
  {
    for(j=i+1;j<Components[comp].NumberOfAtoms;j++)
    {
      if((PseudoAtoms[Components[comp].Type[i]].HasCharges)&&(PseudoAtoms[Components[comp].Type[j]].HasCharges))
      {
        Components[comp].ExcludedIntraChargeCharge[Components[comp].NumberOfExcludedIntraChargeCharge].A=i;
        Components[comp].ExcludedIntraChargeCharge[Components[comp].NumberOfExcludedIntraChargeCharge].B=j;
        Components[comp].NumberOfExcludedIntraChargeCharge++;
      }
    }
  }

  if(Components[comp].NumberOfIntraChargeBondDipole>0)
  {
    ReadLine(line,1024,FilePtr); // skip line
    for(i=0;i<Components[comp].NumberOfIntraChargeBondDipole;i++)
    {
      ReadLine(line,1024,FilePtr);
      sscanf(line,"%d%d%d",&A,&B1,&B2);
      Components[comp].IntraChargeBondDipole[i].A=A;
      Components[comp].IntraChargeBondDipole[i].B=ReturnBondDipole(comp,B1,B2);
    }
  }

  // add exclusion based on defined intra charge-bonddipole Coulombic potentials
  Components[comp].NumberOfExcludedIntraChargeBondDipole=0;
  for(i=0;i<Components[comp].NumberOfAtoms;i++)
  {
    if(PseudoAtoms[Components[comp].Type[i]].HasCharges)
    {
      for(j=0;j<Components[comp].NumberOfBondDipoles;j++)
      {
        Components[comp].ExcludedIntraChargeBondDipole[Components[comp].NumberOfExcludedIntraChargeBondDipole].A=i;
        Components[comp].ExcludedIntraChargeBondDipole[Components[comp].NumberOfExcludedIntraChargeBondDipole].B=j;
        Components[comp].NumberOfExcludedIntraChargeBondDipole++;
      }
    }
  }


  if(Components[comp].NumberOfIntraBondDipoleBondDipole>0)
  {
    ReadLine(line,1024,FilePtr); // skip line
    for(i=0;i<Components[comp].NumberOfIntraBondDipoleBondDipole;i++)
    {
      ReadLine(line,1024,FilePtr);
      sscanf(line,"%d%d%d%d",&A1,&A2,&B1,&B2);
      Components[comp].IntraBondDipoleBondDipole[i].A=MIN2(ReturnBondDipole(comp,A1,A2),ReturnBondDipole(comp,B1,B2));
      Components[comp].IntraBondDipoleBondDipole[i].B=MAX2(ReturnBondDipole(comp,A1,A2),ReturnBondDipole(comp,B1,B2));
    }
  }

  // add exclusion based on defined intra bonddipole-bonddipole Coulombic potentials
  Components[comp].NumberOfExcludedIntraBondDipoleBondDipole=0;
  for(i=0;i<Components[comp].NumberOfBondDipoles-1;i++)
  {
    for(j=i+1;j<Components[comp].NumberOfBondDipoles;j++)
    {
      Components[comp].ExcludedIntraBondDipoleBondDipole[Components[comp].NumberOfExcludedIntraBondDipoleBondDipole].A=i;
      Components[comp].ExcludedIntraBondDipoleBondDipole[Components[comp].NumberOfExcludedIntraBondDipoleBondDipole].B=j;
      Components[comp].NumberOfExcludedIntraBondDipoleBondDipole++;
    }
  }

  // read the defined config moves
  ReadLine(line,1024,FilePtr); // skip line
  ReadLine(line,1024,FilePtr);
  sscanf(line,"%d",&Components[comp].NumberOfConfigMoves);
  ReadLine(line,1024,FilePtr); // skip line

  // allocate config-moves
  Components[comp].NumberOfUnchangedAtomsConfig=(int*)calloc(Components[comp].NumberOfConfigMoves,sizeof(int));
  Components[comp].UnchangedAtomsConfig=(int**)calloc(Components[comp].NumberOfConfigMoves,sizeof(int*));

  for(i=0;i<Components[comp].NumberOfConfigMoves;i++)
  {
    ReadLine(line,1024,FilePtr);
    arg_pointer=line;
    sscanf(arg_pointer,"%d%n",&temp,&n);
    Components[comp].UnchangedAtomsConfig[i]=(int*)calloc(temp,sizeof(int));
    Components[comp].NumberOfUnchangedAtomsConfig[i]=temp;
    for(j=0;j<temp;j++)
    {
      arg_pointer+=n;
      sscanf(arg_pointer,"%d%n",&Components[comp].UnchangedAtomsConfig[i][j],&n);
    }
  }

  // read the defined identity moves
  Components[comp].NumberOfIdentityConfigMoves=0;
  strcpy(line,"empty");
  ReadLine(line,1024,FilePtr); // skip line
  ReadLine(line,1024,FilePtr);
  if(sscanf(line,"%d",&temp)==1)
  {
    Components[comp].NumberOfIdentityConfigMoves=temp;
    ReadLine(line,1024,FilePtr); // skip line
    // allocate config-moves
    Components[comp].NumberOfUnchangedAtomsIdentityConfig=(int*)calloc(Components[comp].NumberOfIdentityConfigMoves,sizeof(int));
    Components[comp].UnchangedAtomsIdentityConfig=(int**)calloc(Components[comp].NumberOfIdentityConfigMoves,sizeof(int*));

    for(i=0;i<Components[comp].NumberOfIdentityConfigMoves;i++)
    {
      ReadLine(line,1024,FilePtr);
      arg_pointer=line;
      sscanf(arg_pointer,"%d%n",&temp,&n);
      Components[comp].UnchangedAtomsIdentityConfig[i]=(int*)calloc(temp,sizeof(int));
      Components[comp].NumberOfUnchangedAtomsIdentityConfig[i]=temp;
      for(j=0;j<temp;j++)
      {
        arg_pointer+=n;
        sscanf(arg_pointer,"%d%n",&Components[comp].UnchangedAtomsIdentityConfig[i][j],&n);
      }
    }
  }
  else
  {
    // allocate config-moves
    Components[comp].NumberOfIdentityConfigMoves=1;
    Components[comp].NumberOfUnchangedAtomsIdentityConfig=(int*)calloc(Components[comp].NumberOfIdentityConfigMoves,sizeof(int));
    Components[comp].UnchangedAtomsIdentityConfig=(int**)calloc(Components[comp].NumberOfIdentityConfigMoves,sizeof(int*));
    Components[comp].UnchangedAtomsIdentityConfig[0]=(int*)calloc(1,sizeof(int));
    Components[comp].NumberOfUnchangedAtomsIdentityConfig[0]=1;
    Components[comp].UnchangedAtomsIdentityConfig[0][0]=Components[comp].StartingBead;
  }

  // set defaults
  Components[comp].LMCMOL=FALSE;

  Components[comp].TranslationMatrix.ax=1.0;
  Components[comp].TranslationMatrix.ay=0.0;
  Components[comp].TranslationMatrix.az=0.0;

  Components[comp].TranslationMatrix.bx=0.0;
  Components[comp].TranslationMatrix.by=1.0;
  Components[comp].TranslationMatrix.bz=0.0;

  Components[comp].TranslationMatrix.cx=0.0;
  Components[comp].TranslationMatrix.cy=0.0;
  Components[comp].TranslationMatrix.cz=1.0;

  Components[comp].SwapEvery=3;

  for(i=0;i<Components[comp].NumberOfAtoms;i++)
  {
    for(k=0;k<NumberOfSystems;k++)
    {
      Components[comp].MaximumCBMCChangeBondLength[k][i]=0.3;
      Components[comp].MaximumCBMCChangeBendAngle[k][i]=0.3;
      Components[comp].MaximumCBMCRotationOnCone[k][i]=0.3;
      Components[comp].CBMCChangeBendAngleAttempts[k][i]=0.0;
      Components[comp].CBMCChangeBendAngleAccepted[k][i]=0.0;
      Components[comp].CBMCRotationOnConeAttempts[k][i]=0.0;
      Components[comp].CBMCRotationOnConeAccepted[k][i]=0.0;
    }
  }

  // search for 1-4 pairs and set scaling factors
  // assumption: all bends must be defined
  // there are 4-ways two bend-angle can be combined to define a quad, the non-mathcing pair is 1-4
  for(i=0;i<Components[comp].NumberOfIntraVDW;i++)
  {
    A=Components[comp].IntraVDW[i].A;
    B=Components[comp].IntraVDW[i].B;

    // LOOP
    for(j=0;j<Components[comp].NumberOfBends;j++)
    {
      A1=Components[comp].Bends[j].A;
      B1=Components[comp].Bends[j].B;
      C1=Components[comp].Bends[j].C;
      for(k=0;k<Components[comp].NumberOfBends;k++)
      {
        A2=Components[comp].Bends[k].A;
        B2=Components[comp].Bends[k].B;
        C2=Components[comp].Bends[k].C;

        // 4 cases
        if((B1==A2)&&(C1==B2)&&(((A1==A)&&(C2==B))||((A1==B)&&(C2==A)))) Components[comp].IntraVDWScaling[i]=Components[comp].Intra14VDWScalingValue;
        if((B1==C2)&&(C1==B2)&&(((A1==A)&&(A2==B))||((A1==B)&&(A2==A)))) Components[comp].IntraVDWScaling[i]=Components[comp].Intra14VDWScalingValue;
        if((A1==B2)&&(B1==A2)&&(((C1==A)&&(C2==B))||((C1==B)&&(C2==A)))) Components[comp].IntraVDWScaling[i]=Components[comp].Intra14VDWScalingValue;
        if((A1==B2)&&(B1==C2)&&(((C1==A)&&(A2==B))||((C1==B)&&(A2==A)))) Components[comp].IntraVDWScaling[i]=Components[comp].Intra14VDWScalingValue;
      }
    }
  }

  // search for 1-4 pairs and set scaling factors
  // assumption: all bends must be defined
  // there are 4-ways two bend-angle can be combined to define a quad, the non-mathcing pair is 1-4
  for(i=0;i<Components[comp].NumberOfIntraChargeCharge;i++)
  {
    A=Components[comp].IntraChargeCharge[i].A;
    B=Components[comp].IntraChargeCharge[i].B;

    // LOOP
    for(j=0;j<Components[comp].NumberOfBends;j++)
    {
      A1=Components[comp].Bends[j].A;
      B1=Components[comp].Bends[j].B;
      C1=Components[comp].Bends[j].C;
      for(k=0;k<Components[comp].NumberOfBends;k++)
      {
        A2=Components[comp].Bends[k].A;
        B2=Components[comp].Bends[k].B;
        C2=Components[comp].Bends[k].C;

        // 4 cases
        if((B1==A2)&&(C1==B2)&&(((A1==A)&&(C2==B))||((A1==B)&&(C2==A)))) Components[comp].IntraChargeChargeScaling[i]=Components[comp].Intra14ChargeChargeScalingValue;
        if((B1==C2)&&(C1==B2)&&(((A1==A)&&(A2==B))||((A1==B)&&(A2==A)))) Components[comp].IntraChargeChargeScaling[i]=Components[comp].Intra14ChargeChargeScalingValue;
        if((A1==B2)&&(B1==A2)&&(((C1==A)&&(C2==B))||((C1==B)&&(C2==A)))) Components[comp].IntraChargeChargeScaling[i]=Components[comp].Intra14ChargeChargeScalingValue;
        if((A1==B2)&&(B1==C2)&&(((C1==A)&&(A2==B))||((C1==B)&&(A2==A)))) Components[comp].IntraChargeChargeScaling[i]=Components[comp].Intra14ChargeChargeScalingValue;
      }
    }
  }

  Components[comp].NumberOfHessianIndices=0;
  for(i=0;i<Components[comp].NumberOfGroups;i++)
  {
    if(Components[comp].Groups[i].Rigid)
      Components[comp].NumberOfHessianIndices++;
    else
      Components[comp].NumberOfHessianIndices+=Components[comp].Groups[i].NumberOfGroupAtoms;
  }


  // biasing spline
  if(Components[comp].ReadBiasingFunction)
  {
    ReadBiasingProfile(comp);
  }

  fclose(FilePtr);
}

void InsertAdsorbateMolecule(void)
{
  int i,type,nr_atoms;
  int NewMolecule;

  // add the number of atoms
  nr_atoms=Components[CurrentComponent].NumberOfAtoms;
  NumberOfAtomsPerSystem[CurrentSystem]+=nr_atoms;
  NumberOfChargesPerSystem[CurrentSystem]+=Components[CurrentComponent].NumberOfCharges;
  NumberOfBondDipolesPerSystem[CurrentSystem]+=Components[CurrentComponent].NumberOfBondDipoles;

  // update the maximum amount of atoms (maximum over all systems)
  LargestNumberOfCoulombicSites=NumberOfChargesPerSystem[CurrentSystem];
  LargestNumberOfBondDipoleSites=NumberOfBondDipolesPerSystem[CurrentSystem];
  for(i=0;i<Framework[CurrentSystem].NumberOfFrameworks;i++)
  {
    LargestNumberOfCoulombicSites+=Framework[CurrentSystem].NumberOfCharges[i];
    LargestNumberOfBondDipoleSites+=Framework[CurrentSystem].NumberOfBondDipoles[i];
  }

  // if the number is largest than the currently allocated memory reallocate the memory for Ewald
  // the hard-coded default here is to extend the arrays with 256 atoms
  if(LargestNumberOfCoulombicSites>=MaxNumberOfCoulombicSites)
  {
    MaxNumberOfCoulombicSites+=MAX2(MaxNumberOfBeads,512);
    if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
      ReallocateEwaldChargeMemory();
  }
  if(LargestNumberOfBondDipoleSites>=MaxNumberOfBondDipoleSites)
  {
    MaxNumberOfBondDipoleSites+=MAX2(MaxNumberOfBeads,512);
    if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
      ReallocateEwaldBondDipoleMemory();
  }

  // update the number of adsorbate molecules
  NumberOfAdsorbateMolecules[CurrentSystem]++;
  Components[CurrentComponent].NumberOfMolecules[CurrentSystem]++;

  // if the number of adsorbate molecules is larger than the currently allocated amount, then reallocate the memory
  // the hard-coded default here is to add the arrays with 50 molecules
  if(NumberOfAdsorbateMolecules[CurrentSystem]>=MaxNumberOfAdsorbateMolecules[CurrentSystem])
  {
    MaxNumberOfAdsorbateMolecules[CurrentSystem]+=256;
    Adsorbates[CurrentSystem]=(ADSORBATE_MOLECULE*)realloc(Adsorbates[CurrentSystem],
           MaxNumberOfAdsorbateMolecules[CurrentSystem]*sizeof(ADSORBATE_MOLECULE));
  }

  // add the new data at the last newly created element
  NewMolecule=NumberOfAdsorbateMolecules[CurrentSystem]-1;
  Adsorbates[CurrentSystem][NewMolecule].NumberOfAtoms=nr_atoms;
  Adsorbates[CurrentSystem][NewMolecule].Type=CurrentComponent;

  // allocate the meory for the atoms and groups
  Adsorbates[CurrentSystem][NewMolecule].Atoms=(ATOM*)calloc(nr_atoms,sizeof(ATOM));
  if(Components[CurrentComponent].NumberOfGroups>0)
    Adsorbates[CurrentSystem][NewMolecule].Groups=(GROUP*)calloc(Components[CurrentComponent].NumberOfGroups,sizeof(GROUP));

  // copy the grown positions etc to the current element
  for(i=0;i<nr_atoms;i++)
  {
    Adsorbates[CurrentSystem][NewMolecule].Atoms[i].Position=NewPosition[CurrentSystem][i];
    Adsorbates[CurrentSystem][NewMolecule].Atoms[i].AnisotropicPosition=TrialAnisotropicPosition[CurrentSystem][i];
    Adsorbates[CurrentSystem][NewMolecule].Atoms[i].Velocity=NewVelocity[CurrentSystem][i];
    Adsorbates[CurrentSystem][NewMolecule].Atoms[i].Force=NewForce[CurrentSystem][i];

    // new Continuous-Fraction scaling factors are taken from the component-information
    Adsorbates[CurrentSystem][NewMolecule].Atoms[i].CFVDWScalingParameter=CFVDWScaling[i];
    Adsorbates[CurrentSystem][NewMolecule].Atoms[i].CFChargeScalingParameter=CFChargeScaling[i];

    type=Components[CurrentComponent].Type[i];
    Adsorbates[CurrentSystem][NewMolecule].Atoms[i].Type=type;
    Adsorbates[CurrentSystem][NewMolecule].Atoms[i].Fixed.x=Components[CurrentComponent].Fixed[i];
    Adsorbates[CurrentSystem][NewMolecule].Atoms[i].Fixed.y=Components[CurrentComponent].Fixed[i];
    Adsorbates[CurrentSystem][NewMolecule].Atoms[i].Fixed.z=Components[CurrentComponent].Fixed[i];
    Adsorbates[CurrentSystem][NewMolecule].Atoms[i].Charge=Components[CurrentComponent].Charge[i];
    NumberOfPseudoAtomsType[CurrentSystem][type]++;
  }

  // update the center of mass
  UpdateGroupCenterOfMassAdsorbate(NewMolecule);

  // compute the quaternion (orientation) from the positions
  ComputeQuaternionAdsorbate(NewMolecule);

  // initialize the velocities
  InitializeVelocityAdsorbate(NewMolecule);

  // modify the degrees of freedom
  for(i=0;i<Components[CurrentComponent].NumberOfGroups;i++)
  {
    if(Components[CurrentComponent].Groups[i].Rigid)
    {
      DegreesOfFreedomAdsorbates[CurrentSystem]+=3;
      DegreesOfFreedomTranslation[CurrentSystem]+=3;
      DegreesOfFreedomTranslationalAdsorbates[CurrentSystem]+=3;
      DegreesOfFreedom[CurrentSystem]+=3;

      DegreesOfFreedomRotation[CurrentSystem]+=Components[CurrentComponent].Groups[i].RotationalDegreesOfFreedom;
      DegreesOfFreedomAdsorbates[CurrentSystem]+=Components[CurrentComponent].Groups[i].RotationalDegreesOfFreedom;
      DegreesOfFreedomRotationalAdsorbates[CurrentSystem]+=Components[CurrentComponent].Groups[i].RotationalDegreesOfFreedom;
      DegreesOfFreedom[CurrentSystem]+=Components[CurrentComponent].Groups[i].RotationalDegreesOfFreedom;
    }
    else
    {
      DegreesOfFreedomTranslation[CurrentSystem]+=3*Components[CurrentComponent].Groups[i].NumberOfGroupAtoms;
      DegreesOfFreedomAdsorbates[CurrentSystem]+=3*Components[CurrentComponent].Groups[i].NumberOfGroupAtoms;
      DegreesOfFreedomTranslationalAdsorbates[CurrentSystem]+=3*Components[CurrentComponent].Groups[i].NumberOfGroupAtoms;
      DegreesOfFreedom[CurrentSystem]+=3*Components[CurrentComponent].Groups[i].NumberOfGroupAtoms;
    }
  }
// reinitialize the Nose-Hoover internal variables based on the current number of molecules for MuPT, MuPTPR and MuVT ensembles 
  if((Ensemble[CurrentSystem]==MuPT)||(Ensemble[CurrentSystem]==MuPTPR)||(Ensemble[CurrentSystem]==MuVT))
  {
    InitializeNoseHooverCurrentSystem();
  }
}

void RemoveAdsorbateMolecule(void)
{
  int i,j,k,type,nr_atoms;
  int LastMolecule;

  // remove the the amount of atoms from the total
  nr_atoms=Components[CurrentComponent].NumberOfAtoms;
  NumberOfAtomsPerSystem[CurrentSystem]-=nr_atoms;
  NumberOfChargesPerSystem[CurrentSystem]-=Components[CurrentComponent].NumberOfCharges;
  NumberOfBondDipolesPerSystem[CurrentSystem]-=Components[CurrentComponent].NumberOfBondDipoles;

  // index of the last molecule
  LastMolecule=NumberOfAdsorbateMolecules[CurrentSystem]-1;

  // free the memory of the current molecule
  free(Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms);
  if(Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Groups)
    free(Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Groups);

  // copy the pointer to the last molecule to the current molecule
  Adsorbates[CurrentSystem][CurrentAdsorbateMolecule]=Adsorbates[CurrentSystem][LastMolecule];

  // set the pointers to NULL
  Adsorbates[CurrentSystem][LastMolecule].Atoms=NULL;
  Adsorbates[CurrentSystem][LastMolecule].Groups=NULL;

  // decrease the molecule counters
  NumberOfAdsorbateMolecules[CurrentSystem]--;
  Components[CurrentComponent].NumberOfMolecules[CurrentSystem]--;

  for(i=0;i<nr_atoms;i++)
  {
    type=Components[CurrentComponent].Type[i];
    NumberOfPseudoAtomsType[CurrentSystem][type]--;
  }

  // shift the fraction-molecule references if needed (they are molecule-'numbers' which have to
  // be changed when you delete a molecule with a lower index).
  for(i=0;i<NumberOfComponents;i++)
  {
    if((Components[i].FractionalMolecule[CurrentSystem]==LastMolecule)&&(!Components[i].ExtraFrameworkMolecule))
      Components[i].FractionalMolecule[CurrentSystem]=CurrentAdsorbateMolecule;
  }

  for(i=0;i<NumberOfComponents;i++)
  {
    for(j=0;j<NumberOfReactions;j++)
    {
      for(k=0;k<ReactantsStoichiometry[j][i];k++)
      {
        if(Components[i].ReactantFractionalMolecules[CurrentSystem][j][k]==LastMolecule)
          Components[i].ReactantFractionalMolecules[CurrentSystem][j][k]=CurrentAdsorbateMolecule;
      }
      for(k=0;k<ProductsStoichiometry[j][i];k++)
      {
        if(Components[i].ProductFractionalMolecules[CurrentSystem][j][k]==LastMolecule)
          Components[i].ProductFractionalMolecules[CurrentSystem][j][k]=CurrentAdsorbateMolecule;
      }
    }
  }


  // modify the degrees of freedom
  for(i=0;i<Components[CurrentComponent].NumberOfGroups;i++)
  {
    if(Components[CurrentComponent].Groups[i].Rigid)
    {
      DegreesOfFreedomAdsorbates[CurrentSystem]-=3;
      DegreesOfFreedomTranslation[CurrentSystem]-=3;
      DegreesOfFreedomTranslationalAdsorbates[CurrentSystem]-=3;
      DegreesOfFreedom[CurrentSystem]-=3;

      DegreesOfFreedomRotation[CurrentSystem]-=Components[CurrentComponent].Groups[i].RotationalDegreesOfFreedom;
      DegreesOfFreedomAdsorbates[CurrentSystem]-=Components[CurrentComponent].Groups[i].RotationalDegreesOfFreedom;
      DegreesOfFreedomRotationalAdsorbates[CurrentSystem]-=Components[CurrentComponent].Groups[i].RotationalDegreesOfFreedom;
      DegreesOfFreedom[CurrentSystem]-=Components[CurrentComponent].Groups[i].RotationalDegreesOfFreedom;
    }
    else
    {
      DegreesOfFreedomTranslation[CurrentSystem]-=3*Components[CurrentComponent].Groups[i].NumberOfGroupAtoms;
      DegreesOfFreedomAdsorbates[CurrentSystem]-=3*Components[CurrentComponent].Groups[i].NumberOfGroupAtoms;
      DegreesOfFreedomTranslationalAdsorbates[CurrentSystem]-=3*Components[CurrentComponent].Groups[i].NumberOfGroupAtoms;
      DegreesOfFreedom[CurrentSystem]-=3*Components[CurrentComponent].Groups[i].NumberOfGroupAtoms;
    }
  }
// reinitialize the Nose-Hoover internal variables based on the current number of molecules for MuPT, MuPTPR and MuVT ensembles 
  if((Ensemble[CurrentSystem]==MuPT)||(Ensemble[CurrentSystem]==MuPTPR)||(Ensemble[CurrentSystem]==MuVT))
  {
    InitializeNoseHooverCurrentSystem();
  }
}



void PrintAdsorbateMolecules(void)
{
  int i,j,nr_atoms;

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    fprintf(stderr, "Index: %d\n",i);
    nr_atoms=Adsorbates[CurrentSystem][i].NumberOfAtoms;
    for(j=0;j<nr_atoms;j++)
      fprintf(stderr, "\tAtom: %d Positon: %f %f %f\n",
        j,
        (double)Adsorbates[CurrentSystem][i].Atoms[j].Position.x,
        (double)Adsorbates[CurrentSystem][i].Atoms[j].Position.y,
        (double)Adsorbates[CurrentSystem][i].Atoms[j].Position.z);
  }
}

void InsertCationMolecule(void)
{
  int i,type,nr_atoms;
  int NewMolecule;

  // add the number of atoms
  nr_atoms=Components[CurrentComponent].NumberOfAtoms;
  NumberOfAtomsPerSystem[CurrentSystem]+=nr_atoms;
  NumberOfChargesPerSystem[CurrentSystem]+=Components[CurrentComponent].NumberOfCharges;
  NumberOfBondDipolesPerSystem[CurrentSystem]+=Components[CurrentComponent].NumberOfAtoms;

  // update the maximum amount of atoms (maximum over all systems)
  LargestNumberOfCoulombicSites=NumberOfChargesPerSystem[CurrentSystem];
  LargestNumberOfBondDipoleSites=NumberOfBondDipolesPerSystem[CurrentSystem];
  for(i=0;i<Framework[CurrentSystem].NumberOfFrameworks;i++)
  {
    LargestNumberOfCoulombicSites+=Framework[CurrentSystem].NumberOfCharges[i];
    LargestNumberOfBondDipoleSites+=Framework[CurrentSystem].NumberOfBondDipoles[i];
  }

  // if the number is largest than the currently allocated memory reallocate the memory for Ewald
  // the hard-coded default here is to extend the arrays with 512 atoms
  if(LargestNumberOfCoulombicSites>=MaxNumberOfCoulombicSites)
  {
    MaxNumberOfCoulombicSites+=MAX2(MaxNumberOfBeads,512);
    if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
      ReallocateEwaldChargeMemory();
  }
  if(LargestNumberOfBondDipoleSites>=MaxNumberOfBondDipoleSites)
  {
    MaxNumberOfBondDipoleSites+=MAX2(MaxNumberOfBeads,512);
    if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
      ReallocateEwaldBondDipoleMemory();
  }

  // update the number of adsorbate molecules
  NumberOfCationMolecules[CurrentSystem]++;
  Components[CurrentComponent].NumberOfMolecules[CurrentSystem]++;

  // if the number of adsorbate molecules is larger than the currently allocated amount, then reallocate the memory
  // the hard-coded default here is to add the arrays with 50 molecules
  if(NumberOfCationMolecules[CurrentSystem]>=MaxNumberOfCationMolecules[CurrentSystem])
  {
    MaxNumberOfCationMolecules[CurrentSystem]+=256;
    Cations[CurrentSystem]=(CATION_MOLECULE*)realloc(Cations[CurrentSystem],
           MaxNumberOfCationMolecules[CurrentSystem]*sizeof(CATION_MOLECULE));
  }

  // add the new data at the last newly created element
  NewMolecule=NumberOfCationMolecules[CurrentSystem]-1;
  Cations[CurrentSystem][NewMolecule].NumberOfAtoms=nr_atoms;
  Cations[CurrentSystem][NewMolecule].Type=CurrentComponent;

  // allocate the meory for the atoms and groups
  Cations[CurrentSystem][NewMolecule].Atoms=(ATOM*)calloc(nr_atoms,sizeof(ATOM));
  if(Components[CurrentComponent].NumberOfGroups>0)
    Cations[CurrentSystem][NewMolecule].Groups=(GROUP*)calloc(Components[CurrentComponent].NumberOfGroups,sizeof(GROUP));

  // copy the grown positions etc to the current element
  for(i=0;i<nr_atoms;i++)
  {
    Cations[CurrentSystem][NewMolecule].Atoms[i].Position=NewPosition[CurrentSystem][i];
    Cations[CurrentSystem][NewMolecule].Atoms[i].AnisotropicPosition=TrialAnisotropicPosition[CurrentSystem][i];
    Cations[CurrentSystem][NewMolecule].Atoms[i].Velocity=NewVelocity[CurrentSystem][i];
    Cations[CurrentSystem][NewMolecule].Atoms[i].Force=NewForce[CurrentSystem][i];

    // new Continuous-Fraction scaling factors are taken from the component-information
    Cations[CurrentSystem][NewMolecule].Atoms[i].CFVDWScalingParameter=CFVDWScaling[i];
    Cations[CurrentSystem][NewMolecule].Atoms[i].CFChargeScalingParameter=CFChargeScaling[i];

    type=Components[CurrentComponent].Type[i];
    Cations[CurrentSystem][NewMolecule].Atoms[i].Type=type;
    Cations[CurrentSystem][NewMolecule].Atoms[i].Fixed.x=Components[CurrentComponent].Fixed[i];
    Cations[CurrentSystem][NewMolecule].Atoms[i].Fixed.y=Components[CurrentComponent].Fixed[i];
    Cations[CurrentSystem][NewMolecule].Atoms[i].Fixed.z=Components[CurrentComponent].Fixed[i];
    Cations[CurrentSystem][NewMolecule].Atoms[i].Charge=Components[CurrentComponent].Charge[i];
    NumberOfPseudoAtomsType[CurrentSystem][type]++;
  }

  // update the center of mass
  UpdateGroupCenterOfMassCation(NewMolecule);

  // compute the quaternion (orientation) from the positions
  ComputeQuaternionCation(NewMolecule);

  // initialize the velocities
  InitializeVelocityCation(NewMolecule);

  // modify the degrees of freedom
  for(i=0;i<Components[CurrentComponent].NumberOfGroups;i++)
  {
    if(Components[CurrentComponent].Groups[i].Rigid)
    {
      DegreesOfFreedomCations[CurrentSystem]+=3;
      DegreesOfFreedomTranslation[CurrentSystem]+=3;
      DegreesOfFreedomTranslationalCations[CurrentSystem]+=3;
      DegreesOfFreedom[CurrentSystem]+=3;

      DegreesOfFreedomRotation[CurrentSystem]+=Components[CurrentComponent].Groups[i].RotationalDegreesOfFreedom;
      DegreesOfFreedomCations[CurrentSystem]+=Components[CurrentComponent].Groups[i].RotationalDegreesOfFreedom;
      DegreesOfFreedomRotationalCations[CurrentSystem]+=Components[CurrentComponent].Groups[i].RotationalDegreesOfFreedom;
      DegreesOfFreedom[CurrentSystem]+=Components[CurrentComponent].Groups[i].RotationalDegreesOfFreedom;
    }
    else
    {
      DegreesOfFreedomTranslation[CurrentSystem]+=3*Components[CurrentComponent].Groups[i].NumberOfGroupAtoms;
      DegreesOfFreedomCations[CurrentSystem]+=3*Components[CurrentComponent].Groups[i].NumberOfGroupAtoms;
      DegreesOfFreedomTranslationalCations[CurrentSystem]+=3*Components[CurrentComponent].Groups[i].NumberOfGroupAtoms;
      DegreesOfFreedom[CurrentSystem]+=3*Components[CurrentComponent].Groups[i].NumberOfGroupAtoms;
    }
  }

  //InitializeNoseHooverCurrentSystem();
}

void RemoveCationMolecule(void)
{
  int i,type,nr_atoms;
  int LastMolecule;

  // remove the the amount of atoms from the total
  nr_atoms=Components[CurrentComponent].NumberOfAtoms;
  NumberOfAtomsPerSystem[CurrentSystem]-=nr_atoms;
  NumberOfChargesPerSystem[CurrentSystem]-=Components[CurrentComponent].NumberOfCharges;
  NumberOfBondDipolesPerSystem[CurrentSystem]-=Components[CurrentComponent].NumberOfBondDipoles;

  // index of the last molecule
  LastMolecule=NumberOfCationMolecules[CurrentSystem]-1;

  // free the memory of the current molecule
  free(Cations[CurrentSystem][CurrentCationMolecule].Atoms);
  if(Cations[CurrentSystem][CurrentCationMolecule].Groups)
    free(Cations[CurrentSystem][CurrentCationMolecule].Groups);

  // copy the pointer to the last molecule to the current molecule
  Cations[CurrentSystem][CurrentCationMolecule]=Cations[CurrentSystem][LastMolecule];

  // set the pointers to NULL
  Cations[CurrentSystem][LastMolecule].Atoms=NULL;
  Cations[CurrentSystem][LastMolecule].Groups=NULL;

  // decrease the molecule counters
  NumberOfCationMolecules[CurrentSystem]--;
  Components[CurrentComponent].NumberOfMolecules[CurrentSystem]--;

  for(i=0;i<nr_atoms;i++)
  {
    type=Components[CurrentComponent].Type[i];
    NumberOfPseudoAtomsType[CurrentSystem][type]--;
  }

  // shift the fraction-molecule references if needed (they are molecule-'numbers' which have to
  // be changed when you delete a molecule with a lower index).
  for(i=0;i<NumberOfComponents;i++)
  {
    if((Components[i].FractionalMolecule[CurrentSystem]==LastMolecule)&&(Components[i].ExtraFrameworkMolecule))
      Components[i].FractionalMolecule[CurrentSystem]=CurrentCationMolecule;
  }

  // modify the degrees of freedom
  for(i=0;i<Components[CurrentComponent].NumberOfGroups;i++)
  {
    if(Components[CurrentComponent].Groups[i].Rigid)
    {
      DegreesOfFreedomCations[CurrentSystem]-=3;
      DegreesOfFreedomTranslation[CurrentSystem]-=3;
      DegreesOfFreedomTranslationalCations[CurrentSystem]-=3;
      DegreesOfFreedom[CurrentSystem]-=3;

      DegreesOfFreedomRotation[CurrentSystem]-=Components[CurrentComponent].Groups[i].RotationalDegreesOfFreedom;
      DegreesOfFreedomCations[CurrentSystem]-=Components[CurrentComponent].Groups[i].RotationalDegreesOfFreedom;
      DegreesOfFreedomRotationalCations[CurrentSystem]-=Components[CurrentComponent].Groups[i].RotationalDegreesOfFreedom;
      DegreesOfFreedom[CurrentSystem]-=Components[CurrentComponent].Groups[i].RotationalDegreesOfFreedom;
    }
    else
    {
      DegreesOfFreedomTranslation[CurrentSystem]-=3*Components[CurrentComponent].Groups[i].NumberOfGroupAtoms;
      DegreesOfFreedomCations[CurrentSystem]-=3*Components[CurrentComponent].Groups[i].NumberOfGroupAtoms;
      DegreesOfFreedomTranslationalCations[CurrentSystem]-=3*Components[CurrentComponent].Groups[i].NumberOfGroupAtoms;
      DegreesOfFreedom[CurrentSystem]-=3*Components[CurrentComponent].Groups[i].NumberOfGroupAtoms;
    }
  }

  //InitializeNoseHooverCurrentSystem();
}

void PrintCationMolecules(void)
{
  int i,j,nr_atoms;

  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    fprintf(stderr, "Index: %d\n",i);
    nr_atoms=Cations[CurrentSystem][i].NumberOfAtoms;
    for(j=0;j<nr_atoms;j++)
      fprintf(stderr, "\tAtom: %d Positon: %f %f %f\n",
        j,
        (double)Cations[CurrentSystem][i].Atoms[j].Position.x,
        (double)Cations[CurrentSystem][i].Atoms[j].Position.y,
        (double)Cations[CurrentSystem][i].Atoms[j].Position.z);
  }
}

// Rescale the MC trial move probability
// Note that the trial moves MUST be checked in the same order
void RescaleComponentProbabilities(void)
{
  int i;
  REAL TotProb;

  for(i=0;i<NumberOfComponents;i++)
  {
    TotProb=Components[i].ProbabilityTranslationMove+
            Components[i].ProbabilityRandomTranslationMove+
            Components[i].ProbabilityRotationMove+
            Components[i].ProbabilityRandomRotationMove+
            Components[i].ProbabilityPartialReinsertionMove+
            Components[i].ProbabilityReinsertionMove+
            Components[i].ProbabilityReinsertionInPlaceMove+
            Components[i].ProbabilityReinsertionInPlaneMove+
            Components[i].ProbabilityIdentityChangeMove+
            Components[i].ProbabilitySwapMove+
            Components[i].ProbabilityCFSwapLambdaMove+
            Components[i].ProbabilityCBCFSwapLambdaMove+
            Components[i].ProbabilityWidomMove+
            Components[i].ProbabilityCFWidomLambdaMove+
            Components[i].ProbabilityGibbsWidomMove+
            Components[i].ProbabilitySurfaceAreaMove+
            Components[i].ProbabilityGibbsChangeMove+
            Components[i].ProbabilityGibbsIdentityChangeMove+
            Components[i].ProbabilityCFGibbsChangeMove+
            Components[i].ProbabilityCBCFGibbsChangeMove+
            Components[i].ProbabilityExchangeFractionalParticleMove+
            Components[i].ProbabilityCFGibbsSwapFractionalMoleculeToOtherBoxMove+
            Components[i].ProbabilityCFGibbsLambdaChangeMove+
            Components[i].ProbabilityCFGibbsFractionalToIntegerMove+
            ProbabilityParallelTemperingMove+
            ProbabilityHyperParallelTemperingMove+
            ProbabilityParallelMolFractionMove+
            ProbabilityChiralInversionMove+
            ProbabilityHybridNVEMove+
            ProbabilityHybridNPHMove+
            ProbabilityHybridNPHPRMove+
            ProbabilityVolumeChangeMove+
            ProbabilityBoxShapeChangeMove+
            ProbabilityGibbsVolumeChangeMove+
            ProbabilityFrameworkChangeMove+
            ProbabilityFrameworkShiftMove+
            ProbabilityCFCRXMCLambdaChangeMove;

    Components[i].ProbabilityRandomTranslationMove+=Components[i].ProbabilityTranslationMove;
    Components[i].ProbabilityRotationMove+=Components[i].ProbabilityRandomTranslationMove;
    Components[i].ProbabilityRandomRotationMove+=Components[i].ProbabilityRotationMove;
    Components[i].ProbabilityPartialReinsertionMove+=Components[i].ProbabilityRandomRotationMove;
    Components[i].ProbabilityReinsertionMove+=Components[i].ProbabilityPartialReinsertionMove;
    Components[i].ProbabilityReinsertionInPlaceMove+=Components[i].ProbabilityReinsertionMove;
    Components[i].ProbabilityReinsertionInPlaneMove+=Components[i].ProbabilityReinsertionInPlaceMove;
    Components[i].ProbabilityIdentityChangeMove+=Components[i].ProbabilityReinsertionInPlaneMove;
    Components[i].ProbabilitySwapMove+=Components[i].ProbabilityIdentityChangeMove;
    Components[i].ProbabilityCFSwapLambdaMove+=Components[i].ProbabilitySwapMove;
    Components[i].ProbabilityCBCFSwapLambdaMove+=Components[i].ProbabilityCFSwapLambdaMove;
    Components[i].ProbabilityWidomMove+=Components[i].ProbabilityCBCFSwapLambdaMove;
    Components[i].ProbabilityCFWidomLambdaMove+=Components[i].ProbabilityWidomMove;
    Components[i].ProbabilityGibbsWidomMove+=Components[i].ProbabilityCFWidomLambdaMove;
    Components[i].ProbabilitySurfaceAreaMove+=Components[i].ProbabilityGibbsWidomMove;
    Components[i].ProbabilityGibbsChangeMove+=Components[i].ProbabilitySurfaceAreaMove;
    Components[i].ProbabilityGibbsIdentityChangeMove+=Components[i].ProbabilityGibbsChangeMove;
    Components[i].ProbabilityCFGibbsChangeMove+=Components[i].ProbabilityGibbsIdentityChangeMove;
    Components[i].ProbabilityCBCFGibbsChangeMove+=Components[i].ProbabilityCFGibbsChangeMove;
    Components[i].ProbabilityExchangeFractionalParticleMove+=Components[i].ProbabilityCBCFGibbsChangeMove;

    Components[i].ProbabilityCFGibbsSwapFractionalMoleculeToOtherBoxMove+=Components[i].ProbabilityExchangeFractionalParticleMove;
    Components[i].ProbabilityCFGibbsLambdaChangeMove+=Components[i].ProbabilityCFGibbsSwapFractionalMoleculeToOtherBoxMove;
    Components[i].ProbabilityCFGibbsFractionalToIntegerMove+=Components[i].ProbabilityCFGibbsLambdaChangeMove;

    Components[i].ProbabilityParallelTemperingMove=ProbabilityParallelTemperingMove+Components[i].ProbabilityCFGibbsFractionalToIntegerMove;
    Components[i].ProbabilityHyperParallelTemperingMove=ProbabilityHyperParallelTemperingMove+Components[i].ProbabilityParallelTemperingMove;
    Components[i].ProbabilityParallelMolFractionMove=ProbabilityParallelMolFractionMove+Components[i].ProbabilityHyperParallelTemperingMove;
    Components[i].ProbabilityChiralInversionMove=ProbabilityChiralInversionMove+Components[i].ProbabilityParallelMolFractionMove;
    Components[i].ProbabilityHybridNVEMove=ProbabilityHybridNVEMove+Components[i].ProbabilityChiralInversionMove;
    Components[i].ProbabilityHybridNPHMove=ProbabilityHybridNPHMove+Components[i].ProbabilityHybridNVEMove;
    Components[i].ProbabilityHybridNPHPRMove=ProbabilityHybridNPHPRMove+Components[i].ProbabilityHybridNPHMove;
    Components[i].ProbabilityVolumeChangeMove=ProbabilityVolumeChangeMove+Components[i].ProbabilityHybridNPHPRMove;
    Components[i].ProbabilityBoxShapeChangeMove=ProbabilityBoxShapeChangeMove+Components[i].ProbabilityVolumeChangeMove;
    Components[i].ProbabilityGibbsVolumeChangeMove=ProbabilityGibbsVolumeChangeMove+Components[i].ProbabilityBoxShapeChangeMove;
    Components[i].ProbabilityFrameworkChangeMove=ProbabilityFrameworkChangeMove+Components[i].ProbabilityGibbsVolumeChangeMove;
    Components[i].ProbabilityFrameworkShiftMove=ProbabilityFrameworkShiftMove+Components[i].ProbabilityFrameworkChangeMove;
    Components[i].ProbabilityCFCRXMCLambdaChangeMove=ProbabilityCFCRXMCLambdaChangeMove+Components[i].ProbabilityFrameworkShiftMove;

    if(TotProb>1e-5)
    {
      Components[i].ProbabilityTranslationMove/=TotProb;
      Components[i].ProbabilityRandomTranslationMove/=TotProb;
      Components[i].ProbabilityRotationMove/=TotProb;
      Components[i].ProbabilityRandomRotationMove/=TotProb;
      Components[i].ProbabilityPartialReinsertionMove/=TotProb;
      Components[i].ProbabilityReinsertionMove/=TotProb;
      Components[i].ProbabilityReinsertionInPlaceMove/=TotProb;
      Components[i].ProbabilityReinsertionInPlaneMove/=TotProb;
      Components[i].ProbabilityIdentityChangeMove/=TotProb;
      Components[i].ProbabilitySwapMove/=TotProb;
      Components[i].ProbabilityCFSwapLambdaMove/=TotProb;
      Components[i].ProbabilityCBCFSwapLambdaMove/=TotProb;
      Components[i].ProbabilityWidomMove/=TotProb;
      Components[i].ProbabilityCFWidomLambdaMove/=TotProb;
      Components[i].ProbabilityGibbsWidomMove/=TotProb;
      Components[i].ProbabilitySurfaceAreaMove/=TotProb;
      Components[i].ProbabilityGibbsChangeMove/=TotProb;
      Components[i].ProbabilityGibbsIdentityChangeMove/=TotProb;
      Components[i].ProbabilityCFGibbsChangeMove/=TotProb;
      Components[i].ProbabilityCBCFGibbsChangeMove/=TotProb;
      Components[i].ProbabilityExchangeFractionalParticleMove/=TotProb;
      Components[i].ProbabilityCFGibbsSwapFractionalMoleculeToOtherBoxMove/=TotProb;
      Components[i].ProbabilityCFGibbsLambdaChangeMove/=TotProb;
      Components[i].ProbabilityCFGibbsFractionalToIntegerMove/=TotProb;

      Components[i].ProbabilityParallelTemperingMove/=TotProb;
      Components[i].ProbabilityHyperParallelTemperingMove/=TotProb;
      Components[i].ProbabilityParallelMolFractionMove/=TotProb;
      Components[i].ProbabilityChiralInversionMove/=TotProb;
      Components[i].ProbabilityHybridNVEMove/=TotProb;
      Components[i].ProbabilityHybridNPHMove/=TotProb;
      Components[i].ProbabilityHybridNPHPRMove/=TotProb;
      Components[i].ProbabilityVolumeChangeMove/=TotProb;
      Components[i].ProbabilityBoxShapeChangeMove/=TotProb;
      Components[i].ProbabilityGibbsVolumeChangeMove/=TotProb;
      Components[i].ProbabilityFrameworkChangeMove/=TotProb;
      Components[i].ProbabilityFrameworkShiftMove/=TotProb;
      Components[i].ProbabilityCFCRXMCLambdaChangeMove/=TotProb;
    }

    Components[i].FractionOfTranslationMove=Components[i].ProbabilityTranslationMove;
    Components[i].FractionOfRandomTranslationMove=Components[i].ProbabilityRandomTranslationMove-Components[i].ProbabilityTranslationMove;
    Components[i].FractionOfRotationMove=Components[i].ProbabilityRotationMove-Components[i].ProbabilityRandomTranslationMove;
    Components[i].FractionOfRandomRotationMove=Components[i].ProbabilityRandomRotationMove-Components[i].ProbabilityRotationMove;
    Components[i].FractionOfPartialReinsertionMove=Components[i].ProbabilityPartialReinsertionMove-Components[i].ProbabilityRandomRotationMove;
    Components[i].FractionOfReinsertionMove=Components[i].ProbabilityReinsertionMove-Components[i].ProbabilityPartialReinsertionMove;
    Components[i].FractionOfReinsertionInPlaceMove=Components[i].ProbabilityReinsertionInPlaceMove-Components[i].ProbabilityReinsertionMove;
    Components[i].FractionOfReinsertionInPlaneMove=Components[i].ProbabilityReinsertionInPlaneMove-Components[i].ProbabilityReinsertionInPlaceMove;
    Components[i].FractionOfIdentityChangeMove=Components[i].ProbabilityIdentityChangeMove-Components[i].ProbabilityReinsertionInPlaneMove;
    Components[i].FractionOfSwapMove=Components[i].ProbabilitySwapMove-Components[i].ProbabilityIdentityChangeMove;
    Components[i].FractionOfCFSwapLambdaMove=Components[i].ProbabilityCFSwapLambdaMove-Components[i].ProbabilitySwapMove;
    Components[i].FractionOfCBCFSwapLambdaMove=Components[i].ProbabilityCBCFSwapLambdaMove-Components[i].ProbabilityCFSwapLambdaMove;
    Components[i].FractionOfWidomMove=Components[i].ProbabilityWidomMove-Components[i].ProbabilityCBCFSwapLambdaMove;
    Components[i].FractionOfCFWidomLambdaMove=Components[i].ProbabilityCFWidomLambdaMove-Components[i].ProbabilityWidomMove;
    Components[i].FractionOfGibbsWidomMove=Components[i].ProbabilityGibbsWidomMove-Components[i].ProbabilityCFWidomLambdaMove;
    Components[i].FractionOfSurfaceAreaMove=Components[i].ProbabilitySurfaceAreaMove-Components[i].ProbabilityGibbsWidomMove;
    Components[i].FractionOfGibbsChangeMove=Components[i].ProbabilityGibbsChangeMove-Components[i].ProbabilitySurfaceAreaMove;
    Components[i].FractionOfGibbsIdentityChangeMove=Components[i].ProbabilityGibbsIdentityChangeMove-Components[i].ProbabilityGibbsChangeMove;
    Components[i].FractionOfCFGibbsChangeMove=Components[i].ProbabilityCFGibbsChangeMove-Components[i].ProbabilityGibbsIdentityChangeMove;
    Components[i].FractionOfCBCFGibbsChangeMove=Components[i].ProbabilityCBCFGibbsChangeMove-Components[i].ProbabilityCFGibbsChangeMove;
    Components[i].FractionOfExchangeFractionalParticleMove=Components[i].ProbabilityExchangeFractionalParticleMove-Components[i].ProbabilityCBCFGibbsChangeMove;

    Components[i].FractionOfCFGibbsSwapFractionalMoleculeToOtherBoxMove=Components[i].ProbabilityCFGibbsSwapFractionalMoleculeToOtherBoxMove-Components[i].ProbabilityExchangeFractionalParticleMove;
    Components[i].FractionOfCFGibbsLambdaChangeMove=Components[i].ProbabilityCFGibbsLambdaChangeMove-Components[i].ProbabilityCFGibbsSwapFractionalMoleculeToOtherBoxMove;
    Components[i].FractionOfCFGibbsFractionalToIntegerMove=Components[i].ProbabilityCFGibbsFractionalToIntegerMove-Components[i].ProbabilityCFGibbsLambdaChangeMove;

    Components[i].FractionOfParallelTemperingMove=Components[i].ProbabilityParallelTemperingMove-Components[i].ProbabilityCFGibbsFractionalToIntegerMove;
    Components[i].FractionOfHyperParallelTemperingMove=Components[i].ProbabilityHyperParallelTemperingMove-Components[i].ProbabilityParallelTemperingMove;
    Components[i].FractionOfParallelMolFractionMove=Components[i].ProbabilityParallelMolFractionMove-Components[i].ProbabilityHyperParallelTemperingMove;
    Components[i].FractionOfChiralInversionMove=Components[i].ProbabilityChiralInversionMove-Components[i].ProbabilityParallelMolFractionMove;
    Components[i].FractionOfHybridNVEMove=Components[i].ProbabilityHybridNVEMove-Components[i].ProbabilityChiralInversionMove;
    Components[i].FractionOfHybridNPHMove=Components[i].ProbabilityHybridNPHMove-Components[i].ProbabilityHybridNVEMove;
    Components[i].FractionOfHybridNPHPRMove=Components[i].ProbabilityHybridNPHPRMove-Components[i].ProbabilityHybridNPHMove;
    Components[i].FractionOfVolumeChangeMove=Components[i].ProbabilityVolumeChangeMove-Components[i].ProbabilityHybridNPHPRMove;
    Components[i].FractionOfBoxShapeChangeMove=Components[i].ProbabilityBoxShapeChangeMove-Components[i].ProbabilityVolumeChangeMove;
    Components[i].FractionOfGibbsVolumeChangeMove=Components[i].ProbabilityGibbsVolumeChangeMove-Components[i].ProbabilityBoxShapeChangeMove;
    Components[i].FractionOfFrameworkChangeMove=Components[i].ProbabilityFrameworkChangeMove-Components[i].ProbabilityGibbsVolumeChangeMove;
    Components[i].FractionOfFrameworkShiftMove=Components[i].ProbabilityFrameworkShiftMove-Components[i].ProbabilityFrameworkChangeMove;
    Components[i].FractionOfCFCRXMCLambdaChangeMove=Components[i].ProbabilityCFCRXMCLambdaChangeMove-Components[i].ProbabilityFrameworkShiftMove;
  }
}

void UpdateGroupCenterOfMassAdsorbate(int m)
{
  int i,l,A;
  int Type;
  REAL Mass,TotalMass;
  VECTOR com;

  Type=Adsorbates[CurrentSystem][m].Type;
  for(l=0;l<Components[Type].NumberOfGroups;l++)
  {
    TotalMass=0.0;
    com.x=com.y=com.z=0.0;
    if(Components[Type].Groups[l].Rigid)
    {
      for(i=0;i<Components[Type].Groups[l].NumberOfGroupAtoms;i++)
      {
        A=Components[Type].Groups[l].Atoms[i];

        Mass=PseudoAtoms[Adsorbates[CurrentSystem][m].Atoms[A].Type].Mass;
        TotalMass+=Mass;
        com.x+=Mass*Adsorbates[CurrentSystem][m].Atoms[A].Position.x;
        com.y+=Mass*Adsorbates[CurrentSystem][m].Atoms[A].Position.y;
        com.z+=Mass*Adsorbates[CurrentSystem][m].Atoms[A].Position.z;
      }
      com.x/=TotalMass;
      com.y/=TotalMass;
      com.z/=TotalMass;
      Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassPosition=com;
    }
  }
}

void UpdateGroupCenterOfMassCation(int m)
{
  int i,l,A;
  int Type;
  REAL Mass,TotalMass;
  VECTOR com;

  Type=Cations[CurrentSystem][m].Type;
  for(l=0;l<Components[Type].NumberOfGroups;l++)
  {
    TotalMass=0.0;
    com.x=com.y=com.z=0.0;
    if(Components[Type].Groups[l].Rigid)
    {
      for(i=0;i<Components[Type].Groups[l].NumberOfGroupAtoms;i++)
      {
        A=Components[Type].Groups[l].Atoms[i];

        Mass=PseudoAtoms[Cations[CurrentSystem][m].Atoms[A].Type].Mass;
        TotalMass+=Mass;
        com.x+=Mass*Cations[CurrentSystem][m].Atoms[A].Position.x;
        com.y+=Mass*Cations[CurrentSystem][m].Atoms[A].Position.y;
        com.z+=Mass*Cations[CurrentSystem][m].Atoms[A].Position.z;
      }
      com.x/=TotalMass;
      com.y/=TotalMass;
      com.z/=TotalMass;
      Cations[CurrentSystem][m].Groups[l].CenterOfMassPosition=com;
    }
  }
}

VECTOR GetCenterOfMassCurrentSystem(void)
{
  int i,l,m,f1;
  int Type,A;
  REAL Mass,TotalMass;
  VECTOR com;

  TotalMass=0.0;
  com.x=com.y=com.z=0.0;

  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    Type=Adsorbates[CurrentSystem][m].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid)
      {
        Mass=Components[Type].Mass;
        TotalMass+=Mass;
        com.x+=Mass*Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassPosition.x;
        com.y+=Mass*Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassPosition.y;
        com.z+=Mass*Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassPosition.z;
      }
      else
      {
        for(i=0;i<Components[Type].Groups[l].NumberOfGroupAtoms;i++)
        {
          A=Components[Type].Groups[l].Atoms[i];

          Mass=PseudoAtoms[Adsorbates[CurrentSystem][m].Atoms[A].Type].Mass;
          TotalMass+=Mass;
          com.x+=Mass*Adsorbates[CurrentSystem][m].Atoms[A].Position.x;
          com.y+=Mass*Adsorbates[CurrentSystem][m].Atoms[A].Position.y;
          com.z+=Mass*Adsorbates[CurrentSystem][m].Atoms[A].Position.z;
        }
      }
    }
  }

  for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
  {
    Type=Cations[CurrentSystem][m].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid)
      {
        Mass=Components[Type].Mass;
        TotalMass+=Mass;
        com.x+=Mass*Cations[CurrentSystem][m].Groups[l].CenterOfMassPosition.x;
        com.y+=Mass*Cations[CurrentSystem][m].Groups[l].CenterOfMassPosition.y;
        com.z+=Mass*Cations[CurrentSystem][m].Groups[l].CenterOfMassPosition.z;
      }
      else
      {
        for(i=0;i<Components[Type].Groups[l].NumberOfGroupAtoms;i++)
        {
          A=Components[Type].Groups[l].Atoms[i];

          Mass=PseudoAtoms[Cations[CurrentSystem][m].Atoms[A].Type].Mass;
          TotalMass+=Mass;
          com.x+=Mass*Cations[CurrentSystem][m].Atoms[A].Position.x;
          com.y+=Mass*Cations[CurrentSystem][m].Atoms[A].Position.y;
          com.z+=Mass*Cations[CurrentSystem][m].Atoms[A].Position.z;
        }
      }
    }
  }


  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    if(Framework[CurrentSystem].FrameworkModels[f1]==FLEXIBLE)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
      {
        Mass=PseudoAtoms[Framework[CurrentSystem].Atoms[f1][i].Type].Mass;
        TotalMass+=Mass;
        com.x+=Mass*Framework[CurrentSystem].Atoms[f1][i].Position.x;
        com.y+=Mass*Framework[CurrentSystem].Atoms[f1][i].Position.y;
        com.z+=Mass*Framework[CurrentSystem].Atoms[f1][i].Position.z;
      }
    }
  }
  com.x/=TotalMass;
  com.y/=TotalMass;
  com.z/=TotalMass;
  return com;
}

VECTOR GetCenterOfCurrentSystem(void)
{
  int i,l,m,f1;
  int Type,A;
  REAL TotalMass;
  VECTOR com;

  TotalMass=0.0;
  com.x=com.y=com.z=0.0;

  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    Type=Adsorbates[CurrentSystem][m].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      for(i=0;i<Components[Type].Groups[l].NumberOfGroupAtoms;i++)
      {
        A=Components[Type].Groups[l].Atoms[i];

        TotalMass+=1.0;
        com.x+=Adsorbates[CurrentSystem][m].Atoms[A].Position.x;
        com.y+=Adsorbates[CurrentSystem][m].Atoms[A].Position.y;
        com.z+=Adsorbates[CurrentSystem][m].Atoms[A].Position.z;
      }
    }
  }

  for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
  {
    Type=Cations[CurrentSystem][m].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      for(i=0;i<Components[Type].Groups[l].NumberOfGroupAtoms;i++)
      {
        A=Components[Type].Groups[l].Atoms[i];

        TotalMass+=1.0;
        com.x+=Cations[CurrentSystem][m].Atoms[A].Position.x;
        com.y+=Cations[CurrentSystem][m].Atoms[A].Position.y;
        com.z+=Cations[CurrentSystem][m].Atoms[A].Position.z;
      }
    }
  }

  //if(FrameworkModel==FLEXIBLE)
  {
    for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
      {
        TotalMass+=1.0;
        com.x+=Framework[CurrentSystem].Atoms[f1][i].Position.x;
        com.y+=Framework[CurrentSystem].Atoms[f1][i].Position.y;
        com.z+=Framework[CurrentSystem].Atoms[f1][i].Position.z;
      }
    }
  }
  com.x/=TotalMass;
  com.y/=TotalMass;
  com.z/=TotalMass;
  return com;
}



VECTOR GetCenterOfMassVelocityCurrentSystem(void)
{
  int i,l,m,fr;
  int Type,A;
  REAL Mass,TotalMass;
  VECTOR com;

  TotalMass=0.0;
  com.x=com.y=com.z=0.0;
  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    Type=Adsorbates[CurrentSystem][m].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid)
      {
        Mass=Components[Type].Mass;
        TotalMass+=Mass;
        com.x+=Mass*Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassVelocity.x;
        com.y+=Mass*Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassVelocity.y;
        com.z+=Mass*Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassVelocity.z;
      }
      else
      {
        for(i=0;i<Components[Type].Groups[l].NumberOfGroupAtoms;i++)
        {
          A=Components[Type].Groups[l].Atoms[i];

          Mass=PseudoAtoms[Adsorbates[CurrentSystem][m].Atoms[A].Type].Mass;
          TotalMass+=Mass;
          com.x+=Mass*Adsorbates[CurrentSystem][m].Atoms[A].Velocity.x;
          com.y+=Mass*Adsorbates[CurrentSystem][m].Atoms[A].Velocity.y;
          com.z+=Mass*Adsorbates[CurrentSystem][m].Atoms[A].Velocity.z;
        }
      }
    }
  }

  for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
  {
    Type=Cations[CurrentSystem][m].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid)
      {
        Mass=Components[Type].Mass;
        TotalMass+=Mass;
        com.x+=Mass*Cations[CurrentSystem][m].Groups[l].CenterOfMassVelocity.x;
        com.y+=Mass*Cations[CurrentSystem][m].Groups[l].CenterOfMassVelocity.y;
        com.z+=Mass*Cations[CurrentSystem][m].Groups[l].CenterOfMassVelocity.z;
      }
      else
      {
        for(i=0;i<Components[Type].Groups[l].NumberOfGroupAtoms;i++)
        {
          A=Components[Type].Groups[l].Atoms[i];

          Mass=PseudoAtoms[Cations[CurrentSystem][m].Atoms[A].Type].Mass;
          TotalMass+=Mass;
          com.x+=Mass*Cations[CurrentSystem][m].Atoms[A].Velocity.x;
          com.y+=Mass*Cations[CurrentSystem][m].Atoms[A].Velocity.y;
          com.z+=Mass*Cations[CurrentSystem][m].Atoms[A].Velocity.z;
        }
      }
    }
  }

  for(fr=0;fr<Framework[CurrentSystem].NumberOfFrameworks;fr++)
  {
    if(Framework[CurrentSystem].FrameworkModels[fr]==FLEXIBLE)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[fr];i++)
      {
        Mass=PseudoAtoms[Framework[CurrentSystem].Atoms[fr][i].Type].Mass;
        TotalMass+=Mass;
        com.x+=Mass*Framework[CurrentSystem].Atoms[fr][i].Velocity.x;
        com.y+=Mass*Framework[CurrentSystem].Atoms[fr][i].Velocity.y;
        com.z+=Mass*Framework[CurrentSystem].Atoms[fr][i].Velocity.z;
      }
    }
  }
  com.x/=TotalMass;
  com.y/=TotalMass;
  com.z/=TotalMass;
  return com;
}


VECTOR GetAdsorbateCenterOfMass(int m)
{
  int i,l;
  int A,Type;
  REAL Mass,TotalMass;
  VECTOR com;

  TotalMass=0.0;
  com.x=com.y=com.z=0.0;
  Type=Adsorbates[CurrentSystem][m].Type;
  for(l=0;l<Components[Type].NumberOfGroups;l++)
  {
    if(Components[Type].Groups[l].Rigid)
    {
      Mass=Components[Type].Mass;
      TotalMass+=Mass;
      com.x+=Mass*Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassPosition.x;
      com.y+=Mass*Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassPosition.y;
      com.z+=Mass*Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassPosition.z;
    }
    else
    {
      for(i=0;i<Components[Type].Groups[l].NumberOfGroupAtoms;i++)
      {
        A=Components[Type].Groups[l].Atoms[i];

        Mass=PseudoAtoms[Adsorbates[CurrentSystem][m].Atoms[A].Type].Mass;
        TotalMass+=Mass;
        com.x+=Mass*Adsorbates[CurrentSystem][m].Atoms[A].Position.x;
        com.y+=Mass*Adsorbates[CurrentSystem][m].Atoms[A].Position.y;
        com.z+=Mass*Adsorbates[CurrentSystem][m].Atoms[A].Position.z;
      }
    }
  }

  com.x/=TotalMass;
  com.y/=TotalMass;
  com.z/=TotalMass;

  return com;
}

VECTOR GetAdsorbateCenterOfMassForce(int m)
{
  int i,l;
  int A,Type;
  VECTOR com;

  com.x=com.y=com.z=0.0;
  Type=Adsorbates[CurrentSystem][m].Type;
  for(l=0;l<Components[Type].NumberOfGroups;l++)
  {
    if(Components[Type].Groups[l].Rigid)
    {
      com.x+=Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassForce.x;
      com.y+=Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassForce.y;
      com.z+=Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassForce.z;
    }
    else
    {
      for(i=0;i<Components[Type].Groups[l].NumberOfGroupAtoms;i++)
      {
        A=Components[Type].Groups[l].Atoms[i];

        com.x+=Adsorbates[CurrentSystem][m].Atoms[A].Force.x;
        com.y+=Adsorbates[CurrentSystem][m].Atoms[A].Force.y;
        com.z+=Adsorbates[CurrentSystem][m].Atoms[A].Force.z;
      }
    }
  }
  return com;
}


VECTOR GetAdsorbateCenterOfCharge(int m)
{
  int i;
  REAL Charge,TotalCharge;
  VECTOR com;

  TotalCharge=0.0;
  com.x=com.y=com.z=0.0;
  for(i=0;i<Adsorbates[CurrentSystem][m].NumberOfAtoms;i++)
  {
    Charge=Adsorbates[CurrentSystem][m].Atoms[i].Charge;
    TotalCharge+=fabs(Charge);
    com.x+=Charge*Adsorbates[CurrentSystem][m].Atoms[i].Position.x;
    com.y+=Charge*Adsorbates[CurrentSystem][m].Atoms[i].Position.y;
    com.z+=Charge*Adsorbates[CurrentSystem][m].Atoms[i].Position.z;
  }
  com.x/=TotalCharge;
  com.y/=TotalCharge;
  com.z/=TotalCharge;
  return com;
}


VECTOR GetAdsorbateCenterOfMassVelocity(int m)
{
  int i,l;
  int Type,A;
  REAL Mass,TotalMass;
  VECTOR com;

  Type=Adsorbates[CurrentSystem][m].Type;
  TotalMass=0.0;
  com.x=com.y=com.z=0.0;
  for(l=0;l<Components[Type].NumberOfGroups;l++)
  {
    if(Components[Type].Groups[l].Rigid)
    {
      Mass=Components[Type].Mass;
      TotalMass+=Mass;
      com.x+=Mass*Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassVelocity.x;
      com.y+=Mass*Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassVelocity.y;
      com.z+=Mass*Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassVelocity.z;
    }
    else
    {
      for(i=0;i<Components[Type].Groups[l].NumberOfGroupAtoms;i++)
      {
        A=Components[Type].Groups[l].Atoms[i];
        Mass=PseudoAtoms[Adsorbates[CurrentSystem][m].Atoms[A].Type].Mass;
        TotalMass+=Mass;
        com.x+=Mass*Adsorbates[CurrentSystem][m].Atoms[A].Velocity.x;
        com.y+=Mass*Adsorbates[CurrentSystem][m].Atoms[A].Velocity.y;
        com.z+=Mass*Adsorbates[CurrentSystem][m].Atoms[A].Velocity.z;
      }
    }
  }
  com.x/=TotalMass;
  com.y/=TotalMass;
  com.z/=TotalMass;
  return com;
}

VECTOR GetCationCenterOfMassVelocity(int m)
{
  int i,l;
  int Type,A;
  REAL Mass,TotalMass;
  VECTOR com;

  Type=Cations[CurrentSystem][m].Type;
  TotalMass=0.0;
  com.x=com.y=com.z=0.0;
  for(l=0;l<Components[Type].NumberOfGroups;l++)
  {
    if(Components[Type].Groups[l].Rigid)
    {
      Mass=Components[Type].Mass;
      TotalMass+=Mass;
      com.x+=Mass*Cations[CurrentSystem][m].Groups[l].CenterOfMassVelocity.x;
      com.y+=Mass*Cations[CurrentSystem][m].Groups[l].CenterOfMassVelocity.y;
      com.z+=Mass*Cations[CurrentSystem][m].Groups[l].CenterOfMassVelocity.z;
    }
    else
    {
      for(i=0;i<Components[Type].Groups[l].NumberOfGroupAtoms;i++)
      {
        A=Components[Type].Groups[l].Atoms[i];

        Mass=PseudoAtoms[Cations[CurrentSystem][m].Atoms[A].Type].Mass;
        TotalMass+=Mass;
        com.x+=Mass*Cations[CurrentSystem][m].Atoms[A].Velocity.x;
        com.y+=Mass*Cations[CurrentSystem][m].Atoms[A].Velocity.y;
        com.z+=Mass*Cations[CurrentSystem][m].Atoms[A].Velocity.z;
      }
    }
  }
  com.x/=TotalMass;
  com.y/=TotalMass;
  com.z/=TotalMass;
  return com;
}


VECTOR GetCationCenterOfMass(int m)
{
  int i,l;
  int A,Type;
  REAL Mass,TotalMass;
  VECTOR com;

  TotalMass=0.0;
  com.x=com.y=com.z=0.0;
  Type=Cations[CurrentSystem][m].Type;
  for(l=0;l<Components[Type].NumberOfGroups;l++)
  {
    if(Components[Type].Groups[l].Rigid)
    {
      Mass=Components[Type].Mass;
      TotalMass+=Mass;
      com.x+=Mass*Cations[CurrentSystem][m].Groups[l].CenterOfMassPosition.x;
      com.y+=Mass*Cations[CurrentSystem][m].Groups[l].CenterOfMassPosition.y;
      com.z+=Mass*Cations[CurrentSystem][m].Groups[l].CenterOfMassPosition.z;
    }
    else
    {
      for(i=0;i<Components[Type].Groups[l].NumberOfGroupAtoms;i++)
      {
        A=Components[Type].Groups[l].Atoms[i];

        Mass=PseudoAtoms[Cations[CurrentSystem][m].Atoms[A].Type].Mass;
        TotalMass+=Mass;
        com.x+=Mass*Cations[CurrentSystem][m].Atoms[A].Position.x;
        com.y+=Mass*Cations[CurrentSystem][m].Atoms[A].Position.y;
        com.z+=Mass*Cations[CurrentSystem][m].Atoms[A].Position.z;
      }
    }
  }

  com.x/=TotalMass;
  com.y/=TotalMass;
  com.z/=TotalMass;
  return com;
}

VECTOR GetCationCenterOfCharge(int m)
{
  int i;
  REAL Charge,TotalCharge;
  VECTOR com;

  TotalCharge=0.0;
  com.x=com.y=com.z=0.0;
  for(i=0;i<Cations[CurrentSystem][m].NumberOfAtoms;i++)
  {
    Charge=Cations[CurrentSystem][m].Atoms[i].Charge;
    TotalCharge+=fabs(Charge);
    com.x+=Charge*Cations[CurrentSystem][m].Atoms[i].Position.x;
    com.y+=Charge*Cations[CurrentSystem][m].Atoms[i].Position.y;
    com.z+=Charge*Cations[CurrentSystem][m].Atoms[i].Position.z;
  }
  com.x/=TotalCharge;
  com.y/=TotalCharge;
  com.z/=TotalCharge;
  return com;
}


REAL GetAdsorbateMass(int m)
{
  int i;
  REAL Mass;

  Mass=0.0;
  for(i=0;i<Adsorbates[CurrentSystem][m].NumberOfAtoms;i++)
    Mass+=PseudoAtoms[Adsorbates[CurrentSystem][m].Atoms[i].Type].Mass;
  return Mass;
}

REAL GetTotalAdsorbateMass(void)
{
  int i;
  REAL Mass;

  Mass=0.0;
  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
    Mass+=GetAdsorbateMass(i);
  return Mass;
}

REAL GetCationMass(int m)
{
  int i;
  REAL Mass;

  Mass=0.0;
  for(i=0;i<Cations[CurrentSystem][m].NumberOfAtoms;i++)
    Mass+=PseudoAtoms[Cations[CurrentSystem][m].Atoms[i].Type].Mass;
  return Mass;
}

REAL GetTotalCationMass(void)
{
  int i;
  REAL Mass;

  Mass=0.0;
  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
    Mass+=GetCationMass(i);
  return Mass;
}

VECTOR ComputeDipoleMomentComponent(int m,int g)
{
  int i,A;
  REAL TotalChargeNegative,TotalChargePositive,Charge;
  VECTOR dipole,coc_negative,coc_positive;

  TotalChargePositive=TotalChargeNegative=0.0;
  coc_negative.x=coc_negative.y=coc_negative.z=0.0;
  coc_positive.x=coc_positive.y=coc_positive.z=0.0;
  for(i=0;i<Components[m].Groups[g].NumberOfGroupAtoms;i++)
  {
    A=Components[m].Groups[g].Atoms[i];
    Charge=PseudoAtoms[Components[m].Type[A]].Charge1;
    if(Charge>0.0)
    {
      TotalChargePositive+=Charge;
      coc_positive.x+=Charge*Components[m].Positions[A].x;
      coc_positive.y+=Charge*Components[m].Positions[A].y;
      coc_positive.z+=Charge*Components[m].Positions[A].z;
    }
    else
    {
      TotalChargeNegative-=Charge;
      coc_negative.x-=Charge*Components[m].Positions[A].x;
      coc_negative.y-=Charge*Components[m].Positions[A].y;
      coc_negative.z-=Charge*Components[m].Positions[A].z;
    }
  }

  dipole.x=dipole.y=dipole.z=0.0;
  if((fabs(TotalChargePositive)>1e-8)&&(fabs(TotalChargeNegative)>1e-8))
  {
    dipole.x=coc_positive.x/TotalChargePositive-coc_negative.x/TotalChargeNegative;
    dipole.y=coc_positive.y/TotalChargePositive-coc_negative.y/TotalChargeNegative;
    dipole.z=coc_positive.z/TotalChargePositive-coc_negative.z/TotalChargeNegative;
    dipole.x*=MIN2(TotalChargePositive,TotalChargeNegative);
    dipole.y*=MIN2(TotalChargePositive,TotalChargeNegative);
    dipole.z*=MIN2(TotalChargePositive,TotalChargeNegative);
  }
  return dipole;
}

// The quadrupole moment is a tensor of second rank and is symmetrical and traceless
// The quadrupole moment gives an indication of how much the charge distribution deviates from spherical symmetry
REAL_MATRIX3x3 ComputeQuadrupoleMomentComponent(int m,int g)
{
  int i,A;
  REAL_MATRIX3x3 quadrupole;
  REAL charge,rr;

  quadrupole.ax=quadrupole.bx=quadrupole.cx=0.0;
  quadrupole.ay=quadrupole.by=quadrupole.cy=0.0;
  quadrupole.az=quadrupole.bz=quadrupole.cz=0.0;

  for(i=0;i<Components[m].Groups[g].NumberOfGroupAtoms;i++)
  {
    A=Components[m].Groups[g].Atoms[i];
    charge=PseudoAtoms[Components[m].Type[A]].Charge1;
    rr=SQR(Components[m].Positions[A].x)+SQR(Components[m].Positions[A].y)+SQR(Components[m].Positions[A].z);

    quadrupole.ax+=0.5*charge*(3.0*Components[m].Positions[A].x*Components[m].Positions[A].x-rr);
    quadrupole.ay+=0.5*charge*(3.0*Components[m].Positions[A].x*Components[m].Positions[A].y);
    quadrupole.az+=0.5*charge*(3.0*Components[m].Positions[A].x*Components[m].Positions[A].z);
    quadrupole.bx+=0.5*charge*(3.0*Components[m].Positions[A].y*Components[m].Positions[A].x);
    quadrupole.by+=0.5*charge*(3.0*Components[m].Positions[A].y*Components[m].Positions[A].y-rr);
    quadrupole.bz+=0.5*charge*(3.0*Components[m].Positions[A].y*Components[m].Positions[A].z);
    quadrupole.cx+=0.5*charge*(3.0*Components[m].Positions[A].z*Components[m].Positions[A].x);
    quadrupole.cy+=0.5*charge*(3.0*Components[m].Positions[A].z*Components[m].Positions[A].y);
    quadrupole.cz+=0.5*charge*(3.0*Components[m].Positions[A].z*Components[m].Positions[A].z-rr);
  }
  return quadrupole;
}



VECTOR ComputeDipoleMomentAdsorbate(int m)
{
  int i;
  REAL TotalChargeNegative,TotalChargePositive,Charge;
  VECTOR dipole,coc_negative,coc_positive;

  TotalChargePositive=TotalChargeNegative=0.0;
  coc_negative.x=coc_negative.y=coc_negative.z=0.0;
  coc_positive.x=coc_positive.y=coc_positive.z=0.0;
  for(i=0;i<Adsorbates[CurrentSystem][m].NumberOfAtoms;i++)
  {
    Charge=Adsorbates[CurrentSystem][m].Atoms[i].Charge;
    if(Charge>0.0)
    {
      TotalChargePositive+=Charge;
      coc_positive.x+=Charge*Adsorbates[CurrentSystem][m].Atoms[i].Position.x;
      coc_positive.y+=Charge*Adsorbates[CurrentSystem][m].Atoms[i].Position.y;
      coc_positive.z+=Charge*Adsorbates[CurrentSystem][m].Atoms[i].Position.z;
    }
    else
    {
      TotalChargeNegative-=Charge;
      coc_negative.x-=Charge*Adsorbates[CurrentSystem][m].Atoms[i].Position.x;
      coc_negative.y-=Charge*Adsorbates[CurrentSystem][m].Atoms[i].Position.y;
      coc_negative.z-=Charge*Adsorbates[CurrentSystem][m].Atoms[i].Position.z;
    }
  }
  dipole.x=coc_positive.x/TotalChargePositive-coc_negative.x/TotalChargeNegative;
  dipole.y=coc_positive.y/TotalChargePositive-coc_negative.y/TotalChargeNegative;
  dipole.z=coc_positive.z/TotalChargePositive-coc_negative.z/TotalChargeNegative;
  dipole.x*=MIN2(TotalChargePositive,TotalChargeNegative);
  dipole.y*=MIN2(TotalChargePositive,TotalChargeNegative);
  dipole.z*=MIN2(TotalChargePositive,TotalChargeNegative);
  return dipole;
}

VECTOR ComputeDipoleMomentCation(int m)
{
  int i;
  REAL TotalChargeNegative,TotalChargePositive,Charge;
  VECTOR dipole,coc_negative,coc_positive;

  TotalChargePositive=TotalChargeNegative=0.0;
  coc_negative.x=coc_negative.y=coc_negative.z=0.0;
  coc_positive.x=coc_positive.y=coc_positive.z=0.0;
  for(i=0;i<Cations[CurrentSystem][m].NumberOfAtoms;i++)
  {
    Charge=Cations[CurrentSystem][m].Atoms[i].Charge;
    if(Charge>0.0)
    {
      TotalChargePositive+=Charge;
      coc_positive.x+=Charge*Cations[CurrentSystem][m].Atoms[i].Position.x;
      coc_positive.y+=Charge*Cations[CurrentSystem][m].Atoms[i].Position.y;
      coc_positive.z+=Charge*Cations[CurrentSystem][m].Atoms[i].Position.z;
    }
    else
    {
      TotalChargeNegative-=Charge;
      coc_negative.x-=Charge*Cations[CurrentSystem][m].Atoms[i].Position.x;
      coc_negative.y-=Charge*Cations[CurrentSystem][m].Atoms[i].Position.y;
      coc_negative.z-=Charge*Cations[CurrentSystem][m].Atoms[i].Position.z;
    }
  }
  dipole.x=coc_positive.x/TotalChargePositive-coc_negative.x/TotalChargeNegative;
  dipole.y=coc_positive.y/TotalChargePositive-coc_negative.y/TotalChargeNegative;
  dipole.z=coc_positive.z/TotalChargePositive-coc_negative.z/TotalChargeNegative;
  return dipole;
}




VECTOR ComputeTotalDipoleMomentSystemAdsorbates(void)
{
  int i;
  VECTOR dipole,total_dipole;

  total_dipole.x=total_dipole.y=total_dipole.z=0.0;
  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    dipole=ComputeDipoleMomentAdsorbate(i);
    total_dipole.x+=dipole.x;
    total_dipole.y+=dipole.y;
    total_dipole.z+=dipole.z;
  }
  return total_dipole;
}

VECTOR ComputeTotalDipoleMomentSystemCations(void)
{
  int i;
  VECTOR dipole,total_dipole;

  total_dipole.x=total_dipole.y=total_dipole.z=0.0;
  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    dipole=ComputeDipoleMomentCation(i);
    total_dipole.x+=dipole.x;
    total_dipole.y+=dipole.y;
    total_dipole.z+=dipole.z;
  }
  return total_dipole;
}

void AdjustVelocitiesToTemperature(void)
{
  int i,k,l,m,fr;
  int Type,A;
  REAL FactorTranslational,FactorRotational;
  REAL UKineticTranslational,UKineticRotational;

  UKineticTranslational=GetTranslationKineticEnergy();
  FactorTranslational=sqrt(DegreesOfFreedomTranslation[CurrentSystem]*
    K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]/(2.0*UKineticTranslational));

  UKineticRotational=GetRotationalKineticEnergy();
  FactorRotational=sqrt(DegreesOfFreedomRotation[CurrentSystem]*
    K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]/(2.0*UKineticRotational));

  // scale the velocities of the Framework
  for(fr=0;fr<Framework[CurrentSystem].NumberOfFrameworks;fr++)
  {
    if(Framework[CurrentSystem].FrameworkModels[fr]==FLEXIBLE)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[fr];i++)
      {
        Framework[CurrentSystem].Atoms[fr][i].Velocity.x*=FactorTranslational;
        Framework[CurrentSystem].Atoms[fr][i].Velocity.y*=FactorTranslational;
        Framework[CurrentSystem].Atoms[fr][i].Velocity.z*=FactorTranslational;
      }
    }
  }

  // scale the velocities of the Adsorbates
  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    Type=Adsorbates[CurrentSystem][m].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid)
      {
        Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassVelocity.x*=FactorTranslational;
        Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassVelocity.y*=FactorTranslational;
        Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassVelocity.z*=FactorTranslational;

        Adsorbates[CurrentSystem][m].Groups[l].AngularVelocity.x*=FactorRotational;
        Adsorbates[CurrentSystem][m].Groups[l].AngularVelocity.y*=FactorRotational;
        Adsorbates[CurrentSystem][m].Groups[l].AngularVelocity.z*=FactorRotational;

        Adsorbates[CurrentSystem][m].Groups[l].QuaternionMomentum=AngularVelocityToQuaternionMomentumAdsorbates(m,l);
      }
      else
      {
        for(k=0;k<Components[Type].Groups[l].NumberOfGroupAtoms;k++)
        {
          A=Components[Type].Groups[l].Atoms[k];
          Adsorbates[CurrentSystem][m].Atoms[A].Velocity.x*=FactorTranslational;
          Adsorbates[CurrentSystem][m].Atoms[A].Velocity.y*=FactorTranslational;
          Adsorbates[CurrentSystem][m].Atoms[A].Velocity.z*=FactorTranslational;
        }
      }
    }
  }

  // scale the velocities of the Cations
  for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
  {
    Type=Cations[CurrentSystem][m].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid)
      {
        Cations[CurrentSystem][m].Groups[l].CenterOfMassVelocity.x*=FactorTranslational;
        Cations[CurrentSystem][m].Groups[l].CenterOfMassVelocity.y*=FactorTranslational;
        Cations[CurrentSystem][m].Groups[l].CenterOfMassVelocity.z*=FactorTranslational;

        Cations[CurrentSystem][m].Groups[l].AngularVelocity.x*=FactorRotational;
        Cations[CurrentSystem][m].Groups[l].AngularVelocity.y*=FactorRotational;
        Cations[CurrentSystem][m].Groups[l].AngularVelocity.z*=FactorRotational;

        Cations[CurrentSystem][m].Groups[l].QuaternionMomentum=AngularVelocityToQuaternionMomentumCations(m,l);
      }
      else
      {
        for(k=0;k<Components[Type].Groups[l].NumberOfGroupAtoms;k++)
        {
          A=Components[Type].Groups[l].Atoms[k];
          Cations[CurrentSystem][m].Atoms[A].Velocity.x*=FactorTranslational;
          Cations[CurrentSystem][m].Atoms[A].Velocity.y*=FactorTranslational;
          Cations[CurrentSystem][m].Atoms[A].Velocity.z*=FactorTranslational;
        }
      }
    }
  }
}

void InitializeVelocityAdsorbate(int m)
{
  int k,l,A,B,Type,iter;
  REAL Mass;
  VECTOR I,posA,posB,dr,velA,velB;
  REAL length,wwA,wwB,max_error;
  REAL MassA,MassB,vel;
  const int max_iter=5000;

  Type=Adsorbates[CurrentSystem][m].Type;
  for(l=0;l<Components[Type].NumberOfGroups;l++)
  {
    if(Components[Type].Groups[l].Rigid)
    {
      Mass=Components[Type].Groups[l].Mass;
      Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassVelocity.x=RandomGaussianNumber()*
                sqrt(K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]/Mass);
      Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassVelocity.y=RandomGaussianNumber()*
                sqrt(K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]/Mass);
      Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassVelocity.z=RandomGaussianNumber()*
                sqrt(K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]/Mass);

      I=Components[Type].Groups[l].InverseInertiaVector;
      Adsorbates[CurrentSystem][m].Groups[l].AngularVelocity.x=RandomGaussianNumber()*
                sqrt(K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]*I.x);
      Adsorbates[CurrentSystem][m].Groups[l].AngularVelocity.y=RandomGaussianNumber()*
                sqrt(K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]*I.y);
      Adsorbates[CurrentSystem][m].Groups[l].AngularVelocity.z=RandomGaussianNumber()*
                sqrt(K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]*I.z);

      Adsorbates[CurrentSystem][m].Groups[l].QuaternionMomentum=AngularVelocityToQuaternionMomentumAdsorbates(m,l);
    }
    else
    {
      for(k=0;k<Components[Type].Groups[l].NumberOfGroupAtoms;k++)
      {
        A=Components[Type].Groups[l].Atoms[k];
        Mass=PseudoAtoms[Adsorbates[CurrentSystem][m].Atoms[A].Type].Mass;
        Adsorbates[CurrentSystem][m].Atoms[A].Velocity.x=RandomGaussianNumber()*
                sqrt(K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]/Mass);
        Adsorbates[CurrentSystem][m].Atoms[A].Velocity.y=RandomGaussianNumber()*
                sqrt(K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]/Mass);
        Adsorbates[CurrentSystem][m].Atoms[A].Velocity.z=RandomGaussianNumber()*
                sqrt(K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]/Mass);
      }
    }
  }

  // for constraints the velocities need to be corrected (there are degrees of freedom removed)
  iter=0;
  do
  {
    max_error=0.0;

    for(k=0;k<Components[Type].NumberOfBonds;k++)
    {
      if(Components[Type].BondType[k]==FIXED_BOND)
      {
        A=Components[Type].Bonds[k].A;
        B=Components[Type].Bonds[k].B;

        MassA=PseudoAtoms[Components[Type].Type[A]].Mass;
        MassB=PseudoAtoms[Components[Type].Type[B]].Mass;

        velA=Adsorbates[CurrentSystem][m].Atoms[A].Velocity;
        velB=Adsorbates[CurrentSystem][m].Atoms[B].Velocity;

        posA=Adsorbates[CurrentSystem][m].Atoms[A].Position;
        posB=Adsorbates[CurrentSystem][m].Atoms[B].Position;
        dr.x=posA.x-posB.x;
        dr.y=posA.y-posB.y;
        dr.z=posA.z-posB.z;
        length=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
        dr.x/=length; dr.y/=length; dr.z/=length;

        vel=dr.x*(velA.x-velB.x)+dr.y*(velA.y-velB.y)+dr.z*(velA.z-velB.z);

        max_error=MAX2(max_error,fabs(vel));

        wwA=vel*MassB/(MassA+MassB);
        wwB=vel*MassA/(MassA+MassB);

        Adsorbates[CurrentSystem][m].Atoms[A].Velocity.x-=wwA*dr.x;
        Adsorbates[CurrentSystem][m].Atoms[A].Velocity.y-=wwA*dr.y;
        Adsorbates[CurrentSystem][m].Atoms[A].Velocity.z-=wwA*dr.z;

        Adsorbates[CurrentSystem][m].Atoms[B].Velocity.x+=wwB*dr.x;
        Adsorbates[CurrentSystem][m].Atoms[B].Velocity.y+=wwB*dr.y;
        Adsorbates[CurrentSystem][m].Atoms[B].Velocity.z+=wwB*dr.z;
      }
    }
  }while(max_error>(1e-12)&&(iter++<max_iter));
}

void InitializeAdsorbateVelocities(void)
{
  int i;

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
    InitializeVelocityAdsorbate(i);
}

void InitializeVelocityAdsorbateToZero(int m)
{
  int k,l,A,Type;

  Type=Adsorbates[CurrentSystem][m].Type;
  for(l=0;l<Components[Type].NumberOfGroups;l++)
  {
    if(Components[Type].Groups[l].Rigid)
    {
      Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassVelocity.x=0.0;
      Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassVelocity.y=0.0;
      Adsorbates[CurrentSystem][m].Groups[l].CenterOfMassVelocity.z=0.0;

      Adsorbates[CurrentSystem][m].Groups[l].AngularVelocity.x=0.0;
      Adsorbates[CurrentSystem][m].Groups[l].AngularVelocity.y=0.0;
      Adsorbates[CurrentSystem][m].Groups[l].AngularVelocity.z=0.0;

      Adsorbates[CurrentSystem][m].Groups[l].QuaternionMomentum=AngularVelocityToQuaternionMomentumAdsorbates(m,l);
    }
    else
    {
      for(k=0;k<Components[Type].Groups[l].NumberOfGroupAtoms;k++)
      {
        A=Components[Type].Groups[l].Atoms[k];
        Adsorbates[CurrentSystem][m].Atoms[A].Velocity.x=0.0;
        Adsorbates[CurrentSystem][m].Atoms[A].Velocity.y=0.0;
        Adsorbates[CurrentSystem][m].Atoms[A].Velocity.z=0.0;
      }
    }
  }
}

void InitializeAdsorbateVelocitiesToZero(void)
{
  int i;

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
    InitializeVelocityAdsorbateToZero(i);
}



void InitializeVelocityCation(int m)
{
  int k,l,A,B,Type,iter;
  REAL Mass;
  VECTOR I,posA,posB,dr,velA,velB;
  REAL length,wwA,wwB,max_error;
  REAL MassA,MassB,vel;
  const int max_iter=5000;

  Type=Cations[CurrentSystem][m].Type;
  for(l=0;l<Components[Type].NumberOfGroups;l++)
  {
    if(Components[Type].Groups[l].Rigid)
    {
      Mass=Components[Type].Groups[l].Mass;
      Cations[CurrentSystem][m].Groups[l].CenterOfMassVelocity.x=RandomGaussianNumber()*
                sqrt(K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]/Mass);
      Cations[CurrentSystem][m].Groups[l].CenterOfMassVelocity.y=RandomGaussianNumber()*
                sqrt(K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]/Mass);
      Cations[CurrentSystem][m].Groups[l].CenterOfMassVelocity.z=RandomGaussianNumber()*
                sqrt(K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]/Mass);

      I=Components[Type].Groups[l].InverseInertiaVector;
      Cations[CurrentSystem][m].Groups[l].AngularVelocity.x=RandomGaussianNumber()*
                sqrt(K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]*I.x);
      Cations[CurrentSystem][m].Groups[l].AngularVelocity.y=RandomGaussianNumber()*
                sqrt(K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]*I.y);
      Cations[CurrentSystem][m].Groups[l].AngularVelocity.z=RandomGaussianNumber()*
                sqrt(K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]*I.z);

      Cations[CurrentSystem][m].Groups[l].QuaternionMomentum=AngularVelocityToQuaternionMomentumCations(m,l);
    }
    else
    {
      for(k=0;k<Components[Type].Groups[l].NumberOfGroupAtoms;k++)
      {
        A=Components[Type].Groups[l].Atoms[k];
        Mass=PseudoAtoms[Cations[CurrentSystem][m].Atoms[A].Type].Mass;
        Cations[CurrentSystem][m].Atoms[A].Velocity.x=RandomGaussianNumber()*
                sqrt(K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]/Mass);
        Cations[CurrentSystem][m].Atoms[A].Velocity.y=RandomGaussianNumber()*
                sqrt(K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]/Mass);
        Cations[CurrentSystem][m].Atoms[A].Velocity.z=RandomGaussianNumber()*
                sqrt(K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]/Mass);
      }
    }
  }

  // for constraints the velocities need to be corrected (there are degrees of freedom removed)
  iter=0;
  do
  {
    max_error=0.0;

    for(k=0;k<Components[Type].NumberOfBonds;k++)
    {
      if(Components[Type].BondType[k]==FIXED_BOND)
      {
        A=Components[Type].Bonds[k].A;
        B=Components[Type].Bonds[k].B;

        MassA=PseudoAtoms[Components[Type].Type[A]].Mass;
        MassB=PseudoAtoms[Components[Type].Type[B]].Mass;

        velA=Cations[CurrentSystem][m].Atoms[A].Velocity;
        velB=Cations[CurrentSystem][m].Atoms[B].Velocity;

        posA=Cations[CurrentSystem][m].Atoms[A].Position;
        posB=Cations[CurrentSystem][m].Atoms[B].Position;
        dr.x=posA.x-posB.x;
        dr.y=posA.y-posB.y;
        dr.z=posA.z-posB.z;
        length=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
        dr.x/=length; dr.y/=length; dr.z/=length;

        vel=dr.x*(velA.x-velB.x)+dr.y*(velA.y-velB.y)+dr.z*(velA.z-velB.z);

        max_error=MAX2(max_error,fabs(vel));

        wwA=vel*MassB/(MassA+MassB);
        wwB=vel*MassA/(MassA+MassB);

        Cations[CurrentSystem][m].Atoms[A].Velocity.x-=wwA*dr.x;
        Cations[CurrentSystem][m].Atoms[A].Velocity.y-=wwA*dr.y;
        Cations[CurrentSystem][m].Atoms[A].Velocity.z-=wwA*dr.z;

        Cations[CurrentSystem][m].Atoms[B].Velocity.x+=wwB*dr.x;
        Cations[CurrentSystem][m].Atoms[B].Velocity.y+=wwB*dr.y;
        Cations[CurrentSystem][m].Atoms[B].Velocity.z+=wwB*dr.z;
      }
    }
  }while(max_error>(1e-12)&&(iter++<max_iter));
}

void InitializeCationVelocities(void)
{
  int i;

  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
    InitializeVelocityCation(i);
}

REAL GetAdsorbateKineticEnergy(void)
{
  int i,j;
  REAL Ukin,Mass;

  Ukin=0.0;
  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
    {
      Mass=PseudoAtoms[Adsorbates[CurrentSystem][i].Atoms[j].Type].Mass;
      Ukin+=0.5*Mass*(SQR(Adsorbates[CurrentSystem][i].Atoms[j].Velocity.x)+
                      SQR(Adsorbates[CurrentSystem][i].Atoms[j].Velocity.y)+
                      SQR(Adsorbates[CurrentSystem][i].Atoms[j].Velocity.z));
    }
  }
  return Ukin;
}

REAL GetCationKineticEnergy(void)
{
  int i,j;
  REAL Ukin,Mass;

  Ukin=0.0;
  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
    {
      Mass=PseudoAtoms[Cations[CurrentSystem][i].Atoms[j].Type].Mass;
      Ukin+=0.5*Mass*(SQR(Cations[CurrentSystem][i].Atoms[j].Velocity.x)+
                      SQR(Cations[CurrentSystem][i].Atoms[j].Velocity.y)+
                      SQR(Cations[CurrentSystem][i].Atoms[j].Velocity.z));
    }
  }
  return Ukin;
}

void ScaleVelocitesToTemperature(void)
{
  int i,j;
  REAL Factor;

  UAdsorbateKinetic[CurrentSystem]=GetAdsorbateKineticEnergy();
  UCationKinetic[CurrentSystem]=GetCationKineticEnergy();
  UKinetic[CurrentSystem]=UAdsorbateKinetic[CurrentSystem]+UCationKinetic[CurrentSystem];
  Factor=sqrt(DegreesOfFreedom[CurrentSystem]*K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]/(2.0*UKinetic[CurrentSystem]));
  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
    {
      Adsorbates[CurrentSystem][i].Atoms[j].Velocity.x*=Factor;
      Adsorbates[CurrentSystem][i].Atoms[j].Velocity.y*=Factor;
      Adsorbates[CurrentSystem][i].Atoms[j].Velocity.z*=Factor;
    }
  }
  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
    {
      Cations[CurrentSystem][i].Atoms[j].Velocity.x*=Factor;
      Cations[CurrentSystem][i].Atoms[j].Velocity.y*=Factor;
      Cations[CurrentSystem][i].Atoms[j].Velocity.z*=Factor;
    }
  }
  UAdsorbateKinetic[CurrentSystem]=GetAdsorbateKineticEnergy();
  UCationKinetic[CurrentSystem]=GetCationKineticEnergy();
  UKinetic[CurrentSystem]=UAdsorbateKinetic[CurrentSystem]+UCationKinetic[CurrentSystem];
}

void SetAdsorbateVelocitesToZero(void)
{
  int i,j;
  int MolType;

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
    {
      Adsorbates[CurrentSystem][i].Atoms[j].Velocity.x=0.0;
      Adsorbates[CurrentSystem][i].Atoms[j].Velocity.y=0.0;
      Adsorbates[CurrentSystem][i].Atoms[j].Velocity.z=0.0;
    }

    MolType=Adsorbates[CurrentSystem][i].Type;
    for(j=0;j<Components[MolType].NumberOfGroups;j++)
    {
      if(Components[MolType].Groups[j].Rigid) // rigid unit
      {
        Adsorbates[CurrentSystem][i].Groups[j].CenterOfMassVelocity.x=0.0;
        Adsorbates[CurrentSystem][i].Groups[j].CenterOfMassVelocity.y=0.0;
        Adsorbates[CurrentSystem][i].Groups[j].CenterOfMassVelocity.z=0.0;
        Adsorbates[CurrentSystem][i].Groups[j].QuaternionMomentum.r=0.0;
        Adsorbates[CurrentSystem][i].Groups[j].QuaternionMomentum.i=0.0;
        Adsorbates[CurrentSystem][i].Groups[j].QuaternionMomentum.j=0.0;
        Adsorbates[CurrentSystem][i].Groups[j].QuaternionMomentum.k=0.0;
      }
    }
  }
}

void SetCationVelocitesToZero(void)
{
  int i,j;
  int MolType;

  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
    {
      Cations[CurrentSystem][i].Atoms[j].Velocity.x=0.0;
      Cations[CurrentSystem][i].Atoms[j].Velocity.y=0.0;
      Cations[CurrentSystem][i].Atoms[j].Velocity.z=0.0;
    }

    MolType=Cations[CurrentSystem][i].Type;
    for(j=0;j<Components[MolType].NumberOfGroups;j++)
    {
      if(Components[MolType].Groups[j].Rigid) // rigid unit
      {
        Cations[CurrentSystem][i].Groups[j].CenterOfMassVelocity.x=0.0;
        Cations[CurrentSystem][i].Groups[j].CenterOfMassVelocity.y=0.0;
        Cations[CurrentSystem][i].Groups[j].CenterOfMassVelocity.z=0.0;
        Cations[CurrentSystem][i].Groups[j].QuaternionMomentum.r=0.0;
        Cations[CurrentSystem][i].Groups[j].QuaternionMomentum.i=0.0;
        Cations[CurrentSystem][i].Groups[j].QuaternionMomentum.j=0.0;
        Cations[CurrentSystem][i].Groups[j].QuaternionMomentum.k=0.0;
      }
    }
  }
}

void ScaleAdsorbateVelocitesToTemperature(void)
{
  int i,j;
  REAL Factor;

  UAdsorbateKinetic[CurrentSystem]=GetAdsorbateKineticEnergy();
  Factor=sqrt(DegreesOfFreedomAdsorbates[CurrentSystem]*K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]/(2.0*UAdsorbateKinetic[CurrentSystem]));
  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
    {
      Adsorbates[CurrentSystem][i].Atoms[j].Velocity.x*=Factor;
      Adsorbates[CurrentSystem][i].Atoms[j].Velocity.y*=Factor;
      Adsorbates[CurrentSystem][i].Atoms[j].Velocity.z*=Factor;
    }
  }
  UAdsorbateKinetic[CurrentSystem]=GetAdsorbateKineticEnergy();
}

void ScaleCationVelocitesToTemperature(void)
{
  int i,j;
  REAL Factor;

  UCationKinetic[CurrentSystem]=GetCationKineticEnergy();
  Factor=sqrt(DegreesOfFreedomCations[CurrentSystem]*K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]/(2.0*UCationKinetic[CurrentSystem]));
  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    for(j=0;j<Cations[CurrentSystem][i].NumberOfAtoms;j++)
    {
      Cations[CurrentSystem][i].Atoms[j].Velocity.x*=Factor;
      Cations[CurrentSystem][i].Atoms[j].Velocity.y*=Factor;
      Cations[CurrentSystem][i].Atoms[j].Velocity.z*=Factor;
    }
  }
  UCationKinetic[CurrentSystem]=GetCationKineticEnergy();
}

VECTOR MeasureVelocityDrift(void)
{
  int i,k,l,Type,A,fr;
  REAL Mass,TotalMass;
  VECTOR com;

  TotalMass=0.0;
  com.x=com.y=com.z=0.0;
  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    Type=Adsorbates[CurrentSystem][i].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid)
      {
        Mass=Components[Type].Groups[l].Mass;
        TotalMass+=Mass;
        com.x+=Mass*Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassVelocity.x;
        com.y+=Mass*Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassVelocity.y;
        com.z+=Mass*Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassVelocity.z;
      }
      else
      {
        for(k=0;k<Components[Type].Groups[l].NumberOfGroupAtoms;k++)
        {
          A=Components[Type].Groups[l].Atoms[k];
          Mass=PseudoAtoms[Adsorbates[CurrentSystem][i].Atoms[A].Type].Mass;
          TotalMass+=Mass;
          com.x+=Mass*Adsorbates[CurrentSystem][i].Atoms[A].Velocity.x;
          com.y+=Mass*Adsorbates[CurrentSystem][i].Atoms[A].Velocity.y;
          com.z+=Mass*Adsorbates[CurrentSystem][i].Atoms[A].Velocity.z;
        }
      }
    }
  }

  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    Type=Cations[CurrentSystem][i].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid)
      {
        Mass=Components[Type].Groups[l].Mass;
        TotalMass+=Mass;
        com.x+=Mass*Cations[CurrentSystem][i].Groups[l].CenterOfMassVelocity.x;
        com.y+=Mass*Cations[CurrentSystem][i].Groups[l].CenterOfMassVelocity.y;
        com.z+=Mass*Cations[CurrentSystem][i].Groups[l].CenterOfMassVelocity.z;
      }
      else
      {
        for(k=0;k<Components[Type].Groups[l].NumberOfGroupAtoms;k++)
        {
          A=Components[Type].Groups[l].Atoms[k];
          Mass=PseudoAtoms[Cations[CurrentSystem][i].Atoms[A].Type].Mass;
          TotalMass+=Mass;
          com.x+=Mass*Cations[CurrentSystem][i].Atoms[A].Velocity.x;
          com.y+=Mass*Cations[CurrentSystem][i].Atoms[A].Velocity.y;
          com.z+=Mass*Cations[CurrentSystem][i].Atoms[A].Velocity.z;
        }
      }
    }
  }

  for(fr=0;fr<Framework[CurrentSystem].NumberOfFrameworks;fr++)
  {
    if(Framework[CurrentSystem].FrameworkModels[fr]==FLEXIBLE)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[fr];i++)
      {
        Mass=PseudoAtoms[Framework[CurrentSystem].Atoms[fr][i].Type].Mass;
        TotalMass+=Mass;
        com.x+=Mass*Framework[CurrentSystem].Atoms[fr][i].Velocity.x;
        com.y+=Mass*Framework[CurrentSystem].Atoms[fr][i].Velocity.y;
        com.z+=Mass*Framework[CurrentSystem].Atoms[fr][i].Velocity.z;
      }
    }
  }

  com.x/=TotalMass;
  com.y/=TotalMass;
  com.z/=TotalMass;
  return com;
}

VECTOR MeasureVelocityDriftFramework(void)
{
  int i,fr;
  REAL Mass,TotalMass;
  VECTOR com;

  TotalMass=0.0;
  com.x=com.y=com.z=0.0;
  for(fr=0;fr<Framework[CurrentSystem].NumberOfFrameworks;fr++)
  {
    if(Framework[CurrentSystem].FrameworkModels[fr]==FLEXIBLE)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[fr];i++)
      {
        Mass=PseudoAtoms[Framework[CurrentSystem].Atoms[fr][i].Type].Mass;
        TotalMass+=Mass;
        com.x+=Mass*Framework[CurrentSystem].Atoms[fr][i].Velocity.x;
        com.y+=Mass*Framework[CurrentSystem].Atoms[fr][i].Velocity.y;
        com.z+=Mass*Framework[CurrentSystem].Atoms[fr][i].Velocity.z;
      }
    }
  }

  com.x/=TotalMass;
  com.y/=TotalMass;
  com.z/=TotalMass;
  return com;
}


VECTOR MeasureVelocityDriftAdsorbate(void)
{
  int i,k,l,Type,A;
  REAL Mass,TotalMass;
  VECTOR com;

  TotalMass=0.0;
  com.x=com.y=com.z=0.0;
  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    Type=Adsorbates[CurrentSystem][i].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid)
      {
        Mass=Components[Type].Groups[l].Mass;
        TotalMass+=Mass;
        com.x+=Mass*Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassVelocity.x;
        com.y+=Mass*Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassVelocity.y;
        com.z+=Mass*Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassVelocity.z;
      }
      else
      {
        for(k=0;k<Components[Type].Groups[l].NumberOfGroupAtoms;k++)
        {
          A=Components[Type].Groups[l].Atoms[k];
          Mass=PseudoAtoms[Adsorbates[CurrentSystem][i].Atoms[A].Type].Mass;
          TotalMass+=Mass;
          com.x+=Mass*Adsorbates[CurrentSystem][i].Atoms[A].Velocity.x;
          com.y+=Mass*Adsorbates[CurrentSystem][i].Atoms[A].Velocity.y;
          com.z+=Mass*Adsorbates[CurrentSystem][i].Atoms[A].Velocity.z;
        }
      }
    }
  }
  com.x/=TotalMass;
  com.y/=TotalMass;
  com.z/=TotalMass;
  return com;
}

void RemoveVelocityDrift(void)
{
  int i,l,k,fr;
  int Type,A;
  VECTOR com;

  com=MeasureVelocityDrift();

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    Type=Adsorbates[CurrentSystem][i].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid)
      {
        Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassVelocity.x-=com.x;
        Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassVelocity.y-=com.y;
        Adsorbates[CurrentSystem][i].Groups[l].CenterOfMassVelocity.z-=com.z;
      }
      else
      {
        for(k=0;k<Components[Type].Groups[l].NumberOfGroupAtoms;k++)
        {
          A=Components[Type].Groups[l].Atoms[k];
          Adsorbates[CurrentSystem][i].Atoms[A].Velocity.x-=com.x;
          Adsorbates[CurrentSystem][i].Atoms[A].Velocity.y-=com.y;
          Adsorbates[CurrentSystem][i].Atoms[A].Velocity.z-=com.z;
        }
      }
    }
  }

  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    Type=Cations[CurrentSystem][i].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      if(Components[Type].Groups[l].Rigid)
      {
        Cations[CurrentSystem][i].Groups[l].CenterOfMassVelocity.x-=com.x;
        Cations[CurrentSystem][i].Groups[l].CenterOfMassVelocity.y-=com.y;
        Cations[CurrentSystem][i].Groups[l].CenterOfMassVelocity.z-=com.z;
      }
      else
      {
        for(k=0;k<Components[Type].Groups[l].NumberOfGroupAtoms;k++)
        {
          A=Components[Type].Groups[l].Atoms[k];
          Cations[CurrentSystem][i].Atoms[A].Velocity.x-=com.x;
          Cations[CurrentSystem][i].Atoms[A].Velocity.y-=com.y;
          Cations[CurrentSystem][i].Atoms[A].Velocity.z-=com.z;
        }
      }
    }
  }

  for(fr=0;fr<Framework[CurrentSystem].NumberOfFrameworks;fr++)
  {
    if(Framework[CurrentSystem].FrameworkModels[fr]==FLEXIBLE)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[fr];i++)
      {
        Framework[CurrentSystem].Atoms[fr][i].Velocity.x-=com.x;
        Framework[CurrentSystem].Atoms[fr][i].Velocity.y-=com.y;
        Framework[CurrentSystem].Atoms[fr][i].Velocity.z-=com.z;
      }
    }
  }

  if(SimulationType==MOLECULAR_DYNAMICS)
  {
    switch(Framework[CurrentSystem].FrameworkModel)
    {
      case NONE:
        if(NumberOfAdsorbateMolecules[CurrentSystem]+NumberOfCationMolecules[CurrentSystem]>1)
        {
          DegreesOfFreedom[CurrentSystem]-=3;
          DegreesOfFreedomTranslation[CurrentSystem]-=3;
          if(!NumberOfAdsorbateMolecules[CurrentSystem])
          {
            DegreesOfFreedomCations[CurrentSystem]-=3;
            DegreesOfFreedomTranslationalCations[CurrentSystem]-=3;
          }
          if(!NumberOfCationMolecules[CurrentSystem])
          {
            DegreesOfFreedomAdsorbates[CurrentSystem]-=3;
            DegreesOfFreedomTranslationalAdsorbates[CurrentSystem]-=3;
          }
        }
        break;
      case FLEXIBLE:
        DegreesOfFreedom[CurrentSystem]-=3;
        DegreesOfFreedomTranslation[CurrentSystem]-=3;
        if(NumberOfAdsorbateMolecules[CurrentSystem]+NumberOfAdsorbateMolecules[CurrentSystem]==0)
          DegreesOfFreedomFramework[CurrentSystem]-=3;
        break;
      default:
        break;
    }
  }
}

VECTOR MapToBox(VECTOR pos)
{
  VECTOR s;

  switch(BoundaryCondition[CurrentSystem])
  {
    case RECTANGULAR:
      pos.x-=Box[CurrentSystem].ax*NINT(pos.x/Box[CurrentSystem].ax);
      pos.y-=Box[CurrentSystem].by*NINT(pos.y/Box[CurrentSystem].by);
      pos.z-=Box[CurrentSystem].cz*NINT(pos.z/Box[CurrentSystem].cz);
      if(pos.x<0.0) pos.x+=Box[CurrentSystem].ax;
      if(pos.y<0.0) pos.y+=Box[CurrentSystem].by;
      if(pos.z<0.0) pos.z+=Box[CurrentSystem].cz;
      break;
    case TRICLINIC:
      s.x=InverseBox[CurrentSystem].ax*pos.x+InverseBox[CurrentSystem].bx*pos.y+InverseBox[CurrentSystem].cx*pos.z;
      s.y=InverseBox[CurrentSystem].ay*pos.x+InverseBox[CurrentSystem].by*pos.y+InverseBox[CurrentSystem].cy*pos.z;
      s.z=InverseBox[CurrentSystem].az*pos.x+InverseBox[CurrentSystem].bz*pos.y+InverseBox[CurrentSystem].cz*pos.z;

      // apply boundary condition
      s.x-=NINT(s.x);
      s.y-=NINT(s.y);
      s.z-=NINT(s.z);
      if(s.x<0.0) s.x+=1.0;
      if(s.y<0.0) s.y+=1.0;
      if(s.z<0.0) s.z+=1.0;

      // convert from abc to xyz
      pos.x=Box[CurrentSystem].ax*s.x+Box[CurrentSystem].bx*s.y+Box[CurrentSystem].cx*s.z;
      pos.y=Box[CurrentSystem].ay*s.x+Box[CurrentSystem].by*s.y+Box[CurrentSystem].cy*s.z;
      pos.z=Box[CurrentSystem].az*s.x+Box[CurrentSystem].bz*s.y+Box[CurrentSystem].cz*s.z;
      break;
  }
  return pos;
}


void MoveAdsorbateCenterOfMassBackInBox(int m)
{
  int i;
  VECTOR com,shifted,displacement;

  com=GetAdsorbateCenterOfMass(m);
  shifted=MapToBox(com);
  displacement.x=shifted.x-com.x;
  displacement.y=shifted.y-com.y;
  displacement.z=shifted.z-com.z;

  for(i=0;i<Adsorbates[CurrentSystem][m].NumberOfAtoms;i++)
  {
    Adsorbates[CurrentSystem][m].Atoms[i].Position.x+=displacement.x;
    Adsorbates[CurrentSystem][m].Atoms[i].Position.y+=displacement.y;
    Adsorbates[CurrentSystem][m].Atoms[i].Position.z+=displacement.z;
  }
}

void MoveCationCenterOfMassBackInBox(int m)
{
  int i;
  VECTOR com,shifted,displacement;

  com=GetCationCenterOfMass(m);
  shifted=MapToBox(com);
  displacement.x=shifted.x-com.x;
  displacement.y=shifted.y-com.y;
  displacement.z=shifted.z-com.z;

  for(i=0;i<Cations[CurrentSystem][m].NumberOfAtoms;i++)
  {
    Cations[CurrentSystem][m].Atoms[i].Position.x+=displacement.x;
    Cations[CurrentSystem][m].Atoms[i].Position.y+=displacement.y;
    Cations[CurrentSystem][m].Atoms[i].Position.z+=displacement.z;
  }
}

void ConstructBondDipolesFromBondsAdsorbates(void)
{
  int i,j,Type,A1,A2;
  VECTOR Dipole,posA1,posA2;
  REAL DipoleMagnitudeA,temp,length;


  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    Type=Adsorbates[CurrentSystem][i].Type;
    for(j=0;j<Components[Type].NumberOfBondDipoles;j++)
    {
      DipoleMagnitudeA=Components[Type].BondDipoleMagnitude[j];
      A1=Components[Type].BondDipoles[j].A;
      A2=Components[Type].BondDipoles[j].B;
      posA1=Adsorbates[CurrentSystem][i].Atoms[A1].Position;
      posA2=Adsorbates[CurrentSystem][i].Atoms[A2].Position;
      Dipole.x=posA2.x-posA1.x;
      Dipole.y=posA2.y-posA1.y;
      Dipole.z=posA2.z-posA1.z;
      length=sqrt(SQR(Dipole.x)+SQR(Dipole.y)+SQR(Dipole.z));
      temp=DipoleMagnitudeA/length;
      Dipole.x*=temp; Dipole.y*=temp; Dipole.z*=temp;
      //Adsorbates[CurrentSystem][i].BondDipoleVector[i]=Dipole;
      //Adsorbates[CurrentSystem][i].BondLength[i]=length;
    }
  }
}

int ReturnAtomBondedToHydrogen(int Type,int c)
{
  int i;
  int A,B;

  for(i=0;i<Components[Type].NumberOfBonds;i++)
  {
    A=Components[Type].Bonds[i].A;
    B=Components[Type].Bonds[i].B;
    if(A==c) return B;
    if(B==c) return A;
  }
  fprintf(stderr, "Error in bonding\n");
  return -1;
}


void CalculateAnisotropicSites(void)
{
  int i,k,A,B,f1;
  int typeA,TypeMolA;
  VECTOR posA;
  VECTOR posL,posR,vec;
  REAL length,d;
  VECTOR drAL,drAR;
  REAL ral,rar;

  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
    {
      typeA=Framework[CurrentSystem].Atoms[f1][i].Type;
      posA=Framework[CurrentSystem].Atoms[f1][i].Position;

      if(PseudoAtoms[typeA].AnisotropicCorrection)
      {
        d=PseudoAtoms[typeA].AnisotropicDisplacement;
        switch(Framework[CurrentSystem].Connectivity[f1][i])
        {
          case 0:
            fprintf(stderr, "Error in routine 'CalculateAnisotropicSites'\n");
            exit(0);
            break;
          case 1:
            A=Framework[CurrentSystem].Neighbours[f1][i][0];
            posL=Framework[CurrentSystem].Atoms[f1][A].Position;
            vec.x=posA.x-posL.x;
            vec.y=posA.y-posL.y;
            vec.z=posA.z-posL.z;
            vec=ApplyBoundaryCondition(vec);
            length=sqrt(SQR(vec.x)+SQR(vec.y)+SQR(vec.z));
            if(PseudoAtoms[typeA].AnisotropicType==ABSOLUTE)
            {
              vec.x/=length;
              vec.y/=length;
              vec.z/=length;
            }
            vec.x*=d;
            vec.y*=d;
            vec.z*=d;
            posA.x+=vec.x;
            posA.y+=vec.y;
            posA.z+=vec.z;
            break;
          case 2:
            switch(Framework[CurrentSystem].AnisotropicType)
            {
              case ANISOTROPIC_BISECTION:
                // use the bisection vector: b/(a+b) vec(AB)+a/(a+b) vec (BC)
                A=Framework[CurrentSystem].Neighbours[f1][i][0];
                B=Framework[CurrentSystem].Neighbours[f1][i][1];
                posL=Framework[CurrentSystem].Atoms[f1][A].Position;
                posR=Framework[CurrentSystem].Atoms[f1][B].Position;
                drAL.x=posA.x-posL.x;
                drAL.y=posA.y-posL.y;
                drAL.z=posA.z-posL.z;
                ral=sqrt(SQR(drAL.x)+SQR(drAL.y)+SQR(drAL.z));
                drAR.x=posA.x-posR.x;
                drAR.y=posA.y-posR.y;
                drAR.z=posA.z-posR.z;
                rar=sqrt(SQR(drAR.x)+SQR(drAR.y)+SQR(drAR.z));
                vec.x=(rar*drAL.x+ral*drAR.x)/(ral+rar);
                vec.y=(rar*drAL.y+ral*drAR.y)/(ral+rar);
                vec.z=(rar*drAL.z+ral*drAR.z)/(ral+rar);
                length=sqrt(SQR(vec.x)+SQR(vec.y)+SQR(vec.z));
                vec.x/=length;
                vec.y/=length;
                vec.z/=length;
                vec.x*=d;
                vec.y*=d;
                vec.z*=d;
                posA.x+=vec.x;
                posA.y+=vec.y;
                posA.z+=vec.z;
                break;
              default:
              case ANISOTROPIC_MID_POINT:
                // use vector from (A+C)/2 to B
                A=Framework[CurrentSystem].Neighbours[f1][i][0];
                B=Framework[CurrentSystem].Neighbours[f1][i][1];
                posL=Framework[CurrentSystem].Atoms[f1][A].Position;
                posR=Framework[CurrentSystem].Atoms[f1][B].Position;
                vec.x=posA.x-0.5*(posL.x+posR.x);
                vec.y=posA.y-0.5*(posL.y+posR.y);
                vec.z=posA.z-0.5*(posL.z+posR.z);
                length=sqrt(SQR(vec.x)+SQR(vec.y)+SQR(vec.z));
                vec.x/=length;
                vec.y/=length;
                vec.z/=length;
                vec.x*=d;
                vec.y*=d;
                vec.z*=d;
                posA.x+=vec.x;
                posA.y+=vec.y;
                posA.z+=vec.z;
                break;
            }
            break;
          case 4:
            break;
          default:
            fprintf(stderr, "ERROR: undefined anisotropic atom with connectivity: %d\n",Framework[CurrentSystem].Connectivity[f1][i]);
            exit(0);
            break;
        }
      }

      Framework[CurrentSystem].Atoms[f1][i].AnisotropicPosition=posA;
    }
  }

  for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
  {
    TypeMolA=Adsorbates[CurrentSystem][i].Type;

    for(k=0;k<Adsorbates[CurrentSystem][i].NumberOfAtoms;k++)
    {
      typeA=Adsorbates[CurrentSystem][i].Atoms[k].Type;
      posA=Adsorbates[CurrentSystem][i].Atoms[k].Position;

      if(PseudoAtoms[typeA].AnisotropicCorrection)
      {
        d=PseudoAtoms[typeA].AnisotropicDisplacement;
        switch(Components[TypeMolA].Connectivity[k])
        {
          case 0:
            fprintf(stderr, "Error in routine 'CalculateAnisotropicSites'\n");
            break;
          case 1:
            A=Components[TypeMolA].ConnectivityList[k][0];
            posL=Adsorbates[CurrentSystem][i].Atoms[A].Position;
            vec.x=posA.x-posL.x;
            vec.y=posA.y-posL.y;
            vec.z=posA.z-posL.z;
            length=sqrt(SQR(vec.x)+SQR(vec.y)+SQR(vec.z));
            if(PseudoAtoms[typeA].AnisotropicType==ABSOLUTE)
            {
              vec.x/=length;
              vec.y/=length;
              vec.z/=length;
              }
            vec.x*=d;
            vec.y*=d;
            vec.z*=d;
            posA.x+=vec.x;
            posA.y+=vec.y;
            posA.z+=vec.z;
            break;
          case 2:
            switch(Components[TypeMolA].AnisotropicType)
            {
              case ANISOTROPIC_BISECTION:
                // use the bisection vector: b/(a+b) vec(AB)+a/(a+b) vec (BC)
                A=Components[TypeMolA].ConnectivityList[k][0];
                B=Components[TypeMolA].ConnectivityList[k][1];
                posL=Adsorbates[CurrentSystem][i].Atoms[A].Position;
                posR=Adsorbates[CurrentSystem][i].Atoms[B].Position;
                drAL.x=posA.x-posL.x;
                drAL.y=posA.y-posL.y;
                drAL.z=posA.z-posL.z;
                ral=sqrt(SQR(drAL.x)+SQR(drAL.y)+SQR(drAL.z));
                drAR.x=posA.x-posR.x;
                drAR.y=posA.y-posR.y;
                drAR.z=posA.z-posR.z;
                rar=sqrt(SQR(drAR.x)+SQR(drAR.y)+SQR(drAR.z));
                vec.x=(rar*drAL.x+ral*drAR.x)/(ral+rar);
                vec.y=(rar*drAL.y+ral*drAR.y)/(ral+rar);
                vec.z=(rar*drAL.z+ral*drAR.z)/(ral+rar);
                length=sqrt(SQR(vec.x)+SQR(vec.y)+SQR(vec.z));
                vec.x/=length;
                vec.y/=length;
                vec.z/=length;
                vec.x*=d;
                vec.y*=d;
                vec.z*=d;
                posA.x+=vec.x;
                posA.y+=vec.y;
                posA.z+=vec.z;
                break;
              default:
              case ANISOTROPIC_MID_POINT:
                // use vector from (A+C)/2 to B
                A=Components[TypeMolA].ConnectivityList[k][0];
                B=Components[TypeMolA].ConnectivityList[k][1];
                posL=Adsorbates[CurrentSystem][i].Atoms[A].Position;
                posR=Adsorbates[CurrentSystem][i].Atoms[B].Position;
                vec.x=posA.x-0.5*(posL.x+posR.x);
                vec.y=posA.y-0.5*(posL.y+posR.y);
                vec.z=posA.z-0.5*(posL.z+posR.z);
                length=sqrt(SQR(vec.x)+SQR(vec.y)+SQR(vec.z));
                vec.x/=length;
                vec.y/=length;
                vec.z/=length;
                vec.x*=d;
                vec.y*=d;
                vec.z*=d;
                posA.x+=vec.x;
                posA.y+=vec.y;
                posA.z+=vec.z;
                break;
            }
            break;
          default:
            fprintf(stderr, "ERROR: undefined anisotropic atom with connecvity: %d\n",Components[TypeMolA].Connectivity[k]);
            exit(0);
            break;
        }
      }

      // set the position of the anisotropic VDW site
      Adsorbates[CurrentSystem][i].Atoms[k].AnisotropicPosition=posA;
    }
  }

  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    TypeMolA=Cations[CurrentSystem][i].Type;

    for(k=0;k<Cations[CurrentSystem][i].NumberOfAtoms;k++)
    {
      typeA=Cations[CurrentSystem][i].Atoms[k].Type;
      posA=Cations[CurrentSystem][i].Atoms[k].Position;

      if(PseudoAtoms[typeA].AnisotropicCorrection)
      {
        d=PseudoAtoms[typeA].AnisotropicDisplacement;
        switch(Components[TypeMolA].Connectivity[k])
        {
          case 0:
            fprintf(stderr, "Error in routine 'CalculateAnisotropicSites'\n");
            break;
          case 1:
            A=Components[TypeMolA].ConnectivityList[k][0];
            posL=Cations[CurrentSystem][i].Atoms[A].Position;
            vec.x=posA.x-posL.x;
            vec.y=posA.y-posL.y;
            vec.z=posA.z-posL.z;
            length=sqrt(SQR(vec.x)+SQR(vec.y)+SQR(vec.z));
            if(PseudoAtoms[typeA].AnisotropicType==ABSOLUTE)
            {
              vec.x/=length;
              vec.y/=length;
              vec.z/=length;
              }
            vec.x*=d;
            vec.y*=d;
            vec.z*=d;
            posA.x+=vec.x;
            posA.y+=vec.y;
            posA.z+=vec.z;
            break;
          case 2:
            switch(Components[TypeMolA].AnisotropicType)
            {
              case ANISOTROPIC_BISECTION:
                // use the bisection vector: b/(a+b) vec(AB)+a/(a+b) vec (BC)
                A=Components[TypeMolA].ConnectivityList[k][0];
                B=Components[TypeMolA].ConnectivityList[k][1];
                posL=Cations[CurrentSystem][i].Atoms[A].Position;
                posR=Cations[CurrentSystem][i].Atoms[B].Position;
                drAL.x=posL.x-posA.x;
                drAL.y=posL.y-posA.y;
                drAL.z=posL.z-posA.z;
                ral=sqrt(SQR(drAL.x)+SQR(drAL.y)+SQR(drAL.z));
                drAR.x=posR.x-posA.x;
                drAR.y=posR.y-posA.y;
                drAR.z=posR.z-posA.z;
                rar=sqrt(SQR(drAR.x)+SQR(drAR.y)+SQR(drAR.z));
                vec.x=(rar*drAL.x+ral*drAR.x)/(ral+rar);
                vec.y=(rar*drAL.y+ral*drAR.y)/(ral+rar);
                vec.z=(rar*drAL.z+ral*drAR.z)/(ral+rar);
                length=sqrt(SQR(vec.x)+SQR(vec.y)+SQR(vec.z));
                vec.x/=length;
                vec.y/=length;
                vec.z/=length;
                vec.x*=d;
                vec.y*=d;
                vec.z*=d;
                posA.x+=vec.x;
                posA.y+=vec.y;
                posA.z+=vec.z;
                break;
              default:
              case ANISOTROPIC_MID_POINT:
                // use vector from (A+C)/2 to B
                A=Components[TypeMolA].ConnectivityList[k][0];
                B=Components[TypeMolA].ConnectivityList[k][1];
                posL=Cations[CurrentSystem][i].Atoms[A].Position;
                posR=Cations[CurrentSystem][i].Atoms[B].Position;
                vec.x=posA.x-0.5*(posL.x+posR.x);
                vec.y=posA.y-0.5*(posL.y+posR.y);
                vec.z=posA.z-0.5*(posL.z+posR.z);
                length=sqrt(SQR(vec.x)+SQR(vec.y)+SQR(vec.z));
                vec.x/=length;
                vec.y/=length;
                vec.z/=length;
                vec.x*=d;
                vec.y*=d;
                vec.z*=d;
                posA.x+=vec.x;
                posA.y+=vec.y;
                posA.z+=vec.z;
                break;
            }
            break;
          default:
            fprintf(stderr, "ERROR: undefined anisotropic atom with connecvity: %d\n",Components[TypeMolA].Connectivity[k]);
            exit(0);
            break;
        }
      }

      // set the position of the anisotropic VDW site
      Cations[CurrentSystem][i].Atoms[k].AnisotropicPosition=posA;
    }
  }
}


void ReadBiasingProfile(int comp)
{
  int i,n;
  FILE *FilePtr;
  REAL F,FA,Sum,norm1,norm2;
  char buffer[256];
  REAL Mass;
  int ReactionBead;
  REAL prefactor;
  double temp1,temp2,temp3;
  REAL JumpDistance,QMapping;
  REAL *FreeEnergyYdata;
  REAL *FreeEnergyWeights;

  if((FilePtr=fopen(Components[comp].BiasingFunctionName,"r")))
  {
    fprintf(stderr, "opening Biasing-file for reading\n");

    fscanf(FilePtr,"# %d%*[^\n]",&n);fscanf(FilePtr,"%*c");
    Components[comp].BiasingFunction.n=n;

    fscanf(FilePtr,"# %lf%*[^\n]",&temp1);fscanf(FilePtr,"%*c");
    Components[comp].QStarA=(REAL)temp1;

    fscanf(FilePtr,"# %lf %lf %lf%*[^\n]",&temp1,&temp2,&temp3);fscanf(FilePtr,"%*c");
    Components[comp].gA=(REAL)temp1;
    Components[comp].gB=(REAL)temp2;
    JumpDistance=(REAL)temp3;

    QMapping=JumpDistance/fabs(temp1-temp2);


    fscanf(FilePtr,"# %lf %lf%*[^\n]",&temp1,&temp2);fscanf(FilePtr,"%*c");
    Components[comp].LeftBoundary=(REAL)temp1;
    Components[comp].RightBoundary=(REAL)temp2;

    // read in spline definitions for the zeolite
    Components[comp].FreeEnergyXdata=(REAL*)calloc(n,sizeof(REAL));
    FreeEnergyYdata=(REAL*)calloc(n,sizeof(REAL));
    FreeEnergyWeights=(REAL*)calloc(n,sizeof(REAL));
    for(i=0;i<n;i++)
    {
      fscanf(FilePtr,"%lf %lf %lf\n",&temp1,&temp2,&temp3);
      Components[comp].FreeEnergyXdata[i]=(REAL)temp1;
      if(fabs(temp2)<100.0)
      {
        FreeEnergyYdata[i]=(REAL)temp2;
        FreeEnergyWeights[i]=(REAL)temp3;
        FreeEnergyWeights[i]=150000.0*fabs(FreeEnergyYdata[i]/FreeEnergyWeights[i]);
      }
      else
      {
        FreeEnergyYdata[i]=100.0;
        FreeEnergyWeights[i]=1.0;
      }
    }

    Components[comp].BiasingFunction=
        CreateCubicFittingSpline(n-1,Components[comp].FreeEnergyXdata,FreeEnergyYdata,
                                  FreeEnergyWeights,PERIODIC_SPLINE,0,0);

    free(FreeEnergyWeights);
    free(FreeEnergyYdata);

    Components[comp].Normalization=
              QuadratureSimpson(Components[comp].BiasingFunction,
              Exp,
              Components[comp].gA,
              Components[comp].gB)*
              QuadratureSimpson(Components[comp].BiasingFunction,
              MExp,
              Components[comp].gA,
              Components[comp].gB);

    sprintf(buffer,"BiasingSpline_%s_%d.dat",Components[comp].Name,comp);
    FilePtr=fopen(buffer,"w");

    fprintf(FilePtr,"# Dividing surfaces: %18.12lf [-]\n",
      (double)Components[comp].QStarA);

    fprintf(FilePtr,"# Free energy minima: %18.12lf [-] %18.12lf [-]  lattice distance: %18.12lf [A]\n",
      (double)Components[comp].gA,
      (double)Components[comp].gB,
      (double)JumpDistance);

    fprintf(FilePtr,"# Left and right boundary: %18.12lf [-] %18.12lf [-]\n",
      (double)Components[comp].LeftBoundary,
      (double)Components[comp].RightBoundary);

    FA=EvaluateCubicSpline(Components[comp].BiasingFunction,
                      Components[comp].BiasingFunction.n,
                      Components[comp].QStarA,
                      Components[comp].FreeEnergyXdata,
                      Derivatives);
    fprintf(FilePtr,"# F(QstarA): %18.12lf\n",(double)FA);
    fprintf(FilePtr,"# Exp(-Beta QStarA): %18.12g\n",(double)exp(-FA));

    Sum=QMapping*QuadratureSimpson(Components[comp].BiasingFunction,
                                   MExp,
                                   Components[comp].gA,
                                   Components[comp].gB);
    //fprintf(FilePtr,"# Integral Exp(-Beta q) over region gA (%g) to gB (%g): %18.12g\n",
    //     (double)Components[comp].gA,(double)Components[comp].gB,(double)(Sum*ANGSTROM));

    Sum=QMapping*QuadratureSimpson(Components[comp].BiasingFunction,
                                   MExp,
                                   Components[comp].LeftBoundary,
                                   Components[comp].QStarA);
    fprintf(FilePtr,"# Integral Exp(-Beta q) over region left boundary (%g) to q* (%g): %18.12g\n",
         (double)Components[comp].LeftBoundary,(double)Components[comp].QStarA,(double)(Sum*ANGSTROM));

    ReactionBead=Components[comp].StartingBead;
    Mass=PseudoAtoms[Components[comp].Type[ReactionBead]].Mass;
    prefactor=sqrt((BOLTZMANN_CONSTANT*therm_baro_stats.ExternalTemperature[0]/(2.0*M_PI*Mass*ATOMIC_MASS_UNIT)));
    fprintf(FilePtr,"# Mass reaction bead: %18.12lf [au]\n",(double)Mass);
    fprintf(FilePtr,"# |v|=Sqrt(k_B T/(2.0*PI*Mass)): %18.12lf [m/s]\n",(double)prefactor);
    fprintf(FilePtr,"# P(q*) dq: %18.12g [1/m]\n",(double)((exp(-FA)/(Sum*ANGSTROM))));
    fprintf(FilePtr,"# k^TST= |v| P(q*) dq, i.e. the TST hopping rate: %18.12g [1/s]\n",
       (double)((prefactor*exp(-FA)/(Sum*ANGSTROM))));
    fprintf(FilePtr,"# D^TST: %18.12g [m^2/s]\n\n",
       (double)((prefactor*exp(-FA)/(Sum*ANGSTROM))*SQR(JumpDistance*ANGSTROM)));

    Components[comp].TST_estimate=(prefactor*exp(-FA)/Sum)*DIFFUSION_CONVERSION_FACTOR;

    norm1=QuadratureSimpson(Components[comp].BiasingFunction,
              Exp,
              Components[comp].gA,
              Components[comp].gB);
    fprintf(FilePtr,"# RM Int1, Integral Exp(Beta q) over region gA to gB: %18.12lf\n",(double)norm1);
    norm2=QuadratureSimpson(Components[comp].BiasingFunction,
              MExp,
              Components[comp].LeftBoundary,
              Components[comp].RightBoundary);
    fprintf(FilePtr,"# RM Int2, Integral Exp(-Beta q) over full region: %18.12lf\n",(double)norm2);
    fprintf(FilePtr,"# RM Int1*Int2: %18.12lf\n",(double)(norm1*norm2));
    fprintf(FilePtr,"# RM 1/(Int1*Int2): %18.12lf\n\n",(double)(1.0/(norm1*norm2)));
    fprintf(FilePtr,"# <n_A>: %18.12lf\n",(double)(Sum/norm2));
    fprintf(FilePtr,"# 1/<n_A>: %18.12lf\n\n",(double)(norm2/Sum));

    for(i=0;i<n;i++)
    {
      F=BiasingPotentialDerivatives(comp,Components[comp].FreeEnergyXdata[i],Derivatives);
      fprintf(FilePtr,"%lf %lf %lf %lf %lf %lf\n",(double)Components[comp].FreeEnergyXdata[i],(double)F,(double)exp(-F),
        (double)Derivatives[0],(double)Derivatives[1],(double)Derivatives[2]);
    }
    fclose(FilePtr);
  }
}

/*********************************************************************************************************
 * Name       | TotalNumberOfIntegerMolecules                                                            *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Dynamically computes the number of integer molecule for the current system.              *
 * Parameters | -                                                                                        *
 * Note       | The number of fractional molecules can vary for the new Gibbs method.                    *
 *********************************************************************************************************/
int TotalNumberOfIntegerMolecules()
{
  int i;
  int total=0;

  for(i=0;i<NumberOfComponents;i++)
  {
    total += Components[i].NumberOfMolecules[CurrentSystem]-(Components[i].FractionalMolecule[CurrentSystem]>=0?1:0)-Components[i].NumberOfRXMCMoleculesPresent[CurrentSystem];
  }
  return total;
}

int TotalNumberOfIntegerAdsorbates()
{
  int i;
  int total=0;

  for(i=0;i<NumberOfComponents;i++)
  {
    if(!Components[i].ExtraFrameworkMolecule)
    {
      total += Components[i].NumberOfMolecules[CurrentSystem]-(Components[i].FractionalMolecule[CurrentSystem]>=0?1:0)-Components[i].NumberOfRXMCMoleculesPresent[CurrentSystem];
    }
  }
  return total;
}

int TotalNumberOfFractionalAdsorbates()
{
  int i;
  int total=0;

  for(i=0;i<NumberOfComponents;i++)
  {
    if(!Components[i].ExtraFrameworkMolecule)
    {
      total += Components[i].FractionalMolecule[CurrentSystem]>=0?1:0;
    }
  }
  return total;
}


int TotalNumberOfIntegerCations()
{
  int i;
  int total=0;

  for(i=0;i<NumberOfComponents;i++)
  {
    if(Components[i].ExtraFrameworkMolecule)
    {
      total += Components[i].NumberOfMolecules[CurrentSystem]-(Components[i].FractionalMolecule[CurrentSystem]>=0?1:0)-Components[i].NumberOfRXMCMoleculesPresent[CurrentSystem];
    }
  }
  return total;
}

int TotalNumberOfFractionalCations()
{
  int i;
  int total=0;

  for(i=0;i<NumberOfComponents;i++)
  {
    if(Components[i].ExtraFrameworkMolecule)
    {
      total += Components[i].FractionalMolecule[CurrentSystem]>=0?1:0;
    }
  }
  return total;
}

int TotalNumberOfFractionalMolecules()
{
  int i;
  int total=0;

  for(i=0;i<NumberOfComponents;i++)
  {
    total += (Components[i].FractionalMolecule[CurrentSystem]>=0?1:0);
  }
  return total;
}

int TotalNumberOfIntegerMoleculesForSystem(int k)
{
  int i;
  int total=0;

  for(i=0;i<NumberOfComponents;i++)
  {
    total += Components[i].NumberOfMolecules[k]-
            (Components[i].FractionalMolecule[k]>=0?1:0)-
             Components[i].NumberOfRXMCMoleculesPresent[k];
  }
  return total;
}

int TotalNumberOfFractionalMoleculesForSystem(int k)
{
  int i;
  int total=0;

  for(i=0;i<NumberOfComponents;i++)
  {
    total += (Components[i].FractionalMolecule[k]>=0?1:0);
  }
  return total;
}

// The biasing is the sum of the weights of the components
REAL CFBiasingWeight()
{
  int i,j,index;
  int  FractionalMolecule;
  REAL weight,Lambda;

  weight=0.0;
  for(i=0;i<NumberOfSystems;i++)
  {
    for(j=0;j<NumberOfComponents;j++)
    {
      FractionalMolecule=Components[j].FractionalMolecule[i];
      if(FractionalMolecule>=0)
      {
        if(Components[j].ExtraFrameworkMolecule)
          Lambda=Cations[i][FractionalMolecule].Atoms[0].CFVDWScalingParameter;
        else
          Lambda=Adsorbates[i][FractionalMolecule].Atoms[0].CFVDWScalingParameter;

        index=(int)(Components[j].CFLambdaHistogramSize*Lambda);
        if(index==Components[j].CFLambdaHistogramSize) index--;
        weight+=Components[j].CFBiasingFactors[i][index];
      }
    }
  }

  return exp(-weight);
}

REAL CFBiasingLambda(int system, int comp)
{
  int  FractionalMolecule;
  REAL Lambda;

  Lambda=0.0;
  if(NumberOfComponents>0)
  {
    FractionalMolecule=Components[comp].FractionalMolecule[system];
    if(FractionalMolecule>=0)
    {
      if(Components[comp].ExtraFrameworkMolecule)
        Lambda=Cations[system][FractionalMolecule].Atoms[0].CFVDWScalingParameter;
      else
        Lambda=Adsorbates[system][FractionalMolecule].Atoms[0].CFVDWScalingParameter;
    }
  }

  return Lambda;
}


/*********************************************************************************************************
 * Name       | ValidFractionalPoint                                                                     *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Determines whether a fractional position is allowed or not.                              *
 * Parameters | -                                                                                        *
 * Note       | You can defined up to 4 fractional regions in a,b,c-space that are excluded.             *
 * Use        | dcTST: to keep particles confined to certain regions.                                    *
 *********************************************************************************************************/

int ValidFractionalPoint(int i, POINT s)
{
  VECTOR pos;

  // compute Cartesian Position in box
  pos.x=Box[CurrentSystem].ax*s.x+Box[CurrentSystem].bx*s.y+Box[CurrentSystem].cx*s.z;
  pos.y=Box[CurrentSystem].ay*s.x+Box[CurrentSystem].by*s.y+Box[CurrentSystem].cy*s.z;
  pos.z=Box[CurrentSystem].az*s.x+Box[CurrentSystem].bz*s.y+Box[CurrentSystem].cz*s.z;

  return ValidCartesianPoint(i,pos);
}

/*********************************************************************************************************
 * Name       | ValidCartesianPoint                                                                      *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Determines whether a Cartesian position is allowed or not.                               *
 * Parameters | -                                                                                        *
 * Note       | You can defined up to 4 fractional regions in a,b,c-space that are excluded.             *
 * Use        | dcTST: to keep particles confined to certain regions.                                    *
 *********************************************************************************************************/

int ValidCartesianPoint(int i, POINT pos)
{
  int k;
  VECTOR CartesianCenter,s,dr;
  REAL rr;


  if(!Components[CurrentComponent].RestrictMoves) return TRUE;

  // convert to fractional position between 0 and 1
  s.x=InverseBox[CurrentSystem].ax*pos.x+InverseBox[CurrentSystem].bx*pos.y+InverseBox[CurrentSystem].cx*pos.z;
  s.y=InverseBox[CurrentSystem].ay*pos.x+InverseBox[CurrentSystem].by*pos.y+InverseBox[CurrentSystem].cy*pos.z;
  s.z=InverseBox[CurrentSystem].az*pos.x+InverseBox[CurrentSystem].bz*pos.y+InverseBox[CurrentSystem].cz*pos.z;
  s.x-=(REAL)NINT(s.x);
  s.y-=(REAL)NINT(s.y);
  s.z-=(REAL)NINT(s.z);
  if(s.x<0.0) s.x+=1.0;
  if(s.y<0.0) s.y+=1.0;
  if(s.z<0.0) s.z+=1.0;

  // compute Cartesian Position in box
  pos.x=Box[CurrentSystem].ax*s.x+Box[CurrentSystem].bx*s.y+Box[CurrentSystem].cx*s.z;
  pos.y=Box[CurrentSystem].ay*s.x+Box[CurrentSystem].by*s.y+Box[CurrentSystem].cy*s.z;
  pos.z=Box[CurrentSystem].az*s.x+Box[CurrentSystem].bz*s.y+Box[CurrentSystem].cz*s.z;

  if(Components[i].RestrictMovesToBox)
  {
    if(((s.x>=Components[i].BoxAxisABC_Min.x)&&(s.x<=Components[i].BoxAxisABC_Max.x))&&
       ((s.y>=Components[i].BoxAxisABC_Min.y)&&(s.y<=Components[i].BoxAxisABC_Max.y))&&
       ((s.z>=Components[i].BoxAxisABC_Min.z)&&(s.z<=Components[i].BoxAxisABC_Max.z))) return TRUE;
    if(((s.x>=Components[i].BoxAxisABC_Min2.x)&&(s.x<=Components[i].BoxAxisABC_Max2.x))&&
       ((s.y>=Components[i].BoxAxisABC_Min2.y)&&(s.y<=Components[i].BoxAxisABC_Max2.y))&&
       ((s.z>=Components[i].BoxAxisABC_Min2.z)&&(s.z<=Components[i].BoxAxisABC_Max2.z))) return TRUE;
    if(((s.x>=Components[i].BoxAxisABC_Min3.x)&&(s.x<=Components[i].BoxAxisABC_Max3.x))&&
       ((s.y>=Components[i].BoxAxisABC_Min3.y)&&(s.y<=Components[i].BoxAxisABC_Max3.y))&&
       ((s.z>=Components[i].BoxAxisABC_Min3.z)&&(s.z<=Components[i].BoxAxisABC_Max3.z))) return TRUE;
    if(((s.x>=Components[i].BoxAxisABC_Min4.x)&&(s.x<=Components[i].BoxAxisABC_Max4.x))&&
       ((s.y>=Components[i].BoxAxisABC_Min4.y)&&(s.y<=Components[i].BoxAxisABC_Max4.y))&&
       ((s.z>=Components[i].BoxAxisABC_Min4.z)&&(s.z<=Components[i].BoxAxisABC_Max4.z))) return TRUE;
  }

  if(Components[i].RestrictMovesToPrisms)
  {
    for(k=0;k<MAX_NUMBER_OF_PRISMS;k++)
    {
      if(Components[i].RestrictMovesToPrism[k])
      {
        if(((s.x>=Components[i].RestrictPrismABC_Min[k].x)&&(s.x<=Components[i].RestrictPrismABC_Max[k].x))&&
           ((s.y>=Components[i].RestrictPrismABC_Min[k].y)&&(s.y<=Components[i].RestrictPrismABC_Max[k].y))&&
           ((s.z>=Components[i].RestrictPrismABC_Min[k].z)&&(s.z<=Components[i].RestrictPrismABC_Max[k].z))) return TRUE;
      }
    }
  }

  if(Components[i].RestrictMovesToCylinders)
  {
    for(k=0;k<MAX_NUMBER_OF_CYLINDERS;k++)
    {
      if(Components[i].RestrictMovesToCylinder[k])
      {
        CartesianCenter=ConvertFromABCtoXYZ(Components[i].RestrictCylinderCenter[k]);
        dr.x=CartesianCenter.x-pos.x;
        dr.y=CartesianCenter.y-pos.y;
        dr.z=CartesianCenter.z-pos.z;
        dr=ApplyBoundaryCondition(dr);
        switch(Components[i].RestrictCylinderDirection[k])
        {
          case X_DIR:
            rr=SQR(dr.y)+SQR(dr.z);
            if((rr<SQR(Components[i].RestrictCylinderRadius[k]))&&(s.x>=Components[i].RestrictCylinderABC_Min[k].x)&&(s.x<=Components[i].RestrictCylinderABC_Max[k].x)) return TRUE;
            break;
          case Y_DIR:
            rr=SQR(dr.x)+SQR(dr.z);
            if((rr<SQR(Components[i].RestrictCylinderRadius[k]))&&(s.y>=Components[i].RestrictCylinderABC_Min[k].y)&&(s.y<=Components[i].RestrictCylinderABC_Max[k].y)) return TRUE;
            break;
          case Z_DIR:
            rr=SQR(dr.x)+SQR(dr.y);
            if((rr<SQR(Components[i].RestrictCylinderRadius[k]))&&(s.z>=Components[i].RestrictCylinderABC_Min[k].z)&&(s.z<=Components[i].RestrictCylinderABC_Max[k].z)) return TRUE;
            break;
          default:
            fprintf(stderr, "ERROR: unknown CylinderDirection in ValidCartesianPoint()\n");
            exit(0);
            break;
        }
      }
    }
  }

  if(Components[i].RestrictMovesToSpheres)
  {
    for(k=0;k<MAX_NUMBER_OF_SPHERES;k++)
    {
      if(Components[i].RestrictMovesToSphere[k])
      {
        CartesianCenter=ConvertFromABCtoXYZ(Components[i].RestrictSphereCenter[k]);
        dr.x=CartesianCenter.x-pos.x;
        dr.y=CartesianCenter.y-pos.y;
        dr.z=CartesianCenter.z-pos.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
        if(rr<SQR(Components[i].RestrictSphereRadius[k])) return TRUE;
      }
    }
  }

  // if not in an allowed region, return 'false'
  return FALSE;
}

void PrintCPUStatistics(FILE *FilePtr)
{
  int i,j;
  REAL CpuTimeTranslationMove,CpuTimeRandomTranslationMove,CpuTimeRotationMove,CpuTimeRandomRotationMove,CpuTimePartialReinsertionMove;
  REAL CpuTimeReinsertionMove,CpuTimeReinsertionInPlaceMove,CpuTimeReinsertionInPlaneMove;
  REAL CpuTimeIdentityChangeMove,CpuTimeSwapMoveInsertion,CpuTimeSwapMoveDeletion;
  REAL CpuTimeCFSwapLambdaMove,CpuTimeCBCFSwapLambdaMove,CpuTimeWidomMove,CpuTimeCFWidomLambdaMove,CpuTimeGibbsWidomMove;
  REAL CpuTimeSurfaceAreaMove,CpuTimeGibbsChangeMove,CpuTimeCFGibbsChangeMove;
  REAL CpuTimeCBCFGibbsChangeMove,CpuTimeGibbsIdentityChangeMove,CpuTimeExchangeFractionalParticleMove;
  REAL CpuTimeCFGibbsSwapFractionalMoleculeToOtherBoxMove,CpuTimeCFGibbsLambdaChangeMove,CpuTimeCFGibbsFractionalToIntegerMove;
  REAL CpuTimeParallelTemperingMoveTotal,CpuTimeHyperParallelTemperingMoveTotal,CpuTimeParallelMolFractionMoveTotal;
  REAL CpuTimeChiralInversionMoveTotal,CpuTimeHybridNVEMoveTotal,CpuTimeHybridNPHMoveTotal;
  REAL CpuTimeHybridNPHPRMoveTotal,CpuTimeVolumeChangeMoveTotal,CpuTimeBoxShapeChangeMoveTotal;
  REAL CpuTimeGibbsVolumeChangeMoveTotal,CpuTimeFrameworkChangeMoveTotal,CpuTimeFrameworkShiftMoveTotal;
  REAL CpuTimeCFCRXMCLambdaChangeMoveTotal;

  CpuTimeTranslationMove=0.0;
  CpuTimeRandomTranslationMove=0.0;
  CpuTimeRotationMove=0.0;
  CpuTimeRandomRotationMove=0.0;
  CpuTimePartialReinsertionMove=0.0;
  CpuTimeReinsertionMove=0.0;
  CpuTimeReinsertionInPlaceMove=0.0;
  CpuTimeReinsertionInPlaneMove=0.0;
  CpuTimeIdentityChangeMove=0.0;
  CpuTimeSwapMoveInsertion=0.0;
  CpuTimeSwapMoveDeletion=0.0;
  CpuTimeCFSwapLambdaMove=0.0;
  CpuTimeCBCFSwapLambdaMove=0.0;
  CpuTimeWidomMove=0.0;
  CpuTimeCFWidomLambdaMove=0.0;
  CpuTimeGibbsWidomMove=0.0;
  CpuTimeSurfaceAreaMove=0.0;
  CpuTimeGibbsChangeMove=0.0;
  CpuTimeCFGibbsChangeMove=0.0;
  CpuTimeCBCFGibbsChangeMove=0.0;
  CpuTimeGibbsIdentityChangeMove=0.0;
  CpuTimeExchangeFractionalParticleMove=0.0;
  CpuTimeCFGibbsSwapFractionalMoleculeToOtherBoxMove=0.0;
  CpuTimeCFGibbsLambdaChangeMove=0.0;
  CpuTimeCFGibbsFractionalToIntegerMove=0.0;

  fprintf(FilePtr,"Total CPU timings:\n");
  fprintf(FilePtr,"===========================================\n");
  fprintf(FilePtr,"initialization:    %18.10g [s]\n",CpuTimeInitialization);
  fprintf(FilePtr,"equilibration:     %18.10g [s]\n",CpuTimeEquilibration);
  fprintf(FilePtr,"production run:    %18.10g [s]\n",CpuTimeProductionRun);
  fprintf(FilePtr,"total time:        %18.10g [s]\n\n",CpuTotal);


  fprintf(FilePtr,"Production run CPU timings of the MC moves:\n");
  fprintf(FilePtr,"===========================================\n");
  for(i=0;i<NumberOfComponents;i++)
  {
    fprintf(FilePtr,"Component: %d (%s)\n",i,Components[i].Name);
    fprintf(FilePtr,"\ttranslation:                        %18.10g [s]\n",Components[i].CpuTimeTranslationMove[CurrentSystem]);
    fprintf(FilePtr,"\trandom translation:                 %18.10g [s]\n",Components[i].CpuTimeRandomTranslationMove[CurrentSystem]);
    fprintf(FilePtr,"\trotation:                           %18.10g [s]\n",Components[i].CpuTimeRotationMove[CurrentSystem]);
    fprintf(FilePtr,"\trandom rotation:                    %18.10g [s]\n",Components[i].CpuTimeRandomRotationMove[CurrentSystem]);
    fprintf(FilePtr,"\tpartial reinsertion:                %18.10g [s]\n",Components[i].CpuTimePartialReinsertionMove[CurrentSystem]);
    fprintf(FilePtr,"\treinsertion:                        %18.10g [s]\n",Components[i].CpuTimeReinsertionMove[CurrentSystem]);
    fprintf(FilePtr,"\treinsertion in-place:               %18.10g [s]\n",Components[i].CpuTimeReinsertionInPlaceMove[CurrentSystem]);
    fprintf(FilePtr,"\treinsertion in-plane:               %18.10g [s]\n",Components[i].CpuTimeReinsertionInPlaneMove[CurrentSystem]);
    fprintf(FilePtr,"\tidentity switch:                    %18.10g [s]\n",Components[i].CpuTimeIdentityChangeMove[CurrentSystem]);
    fprintf(FilePtr,"\tswap (insertion):                   %18.10g [s]\n",Components[i].CpuTimeSwapMoveInsertion[CurrentSystem]);
    fprintf(FilePtr,"\tswap (deletion):                    %18.10g [s]\n",Components[i].CpuTimeSwapMoveDeletion[CurrentSystem]);
    fprintf(FilePtr,"\tswap lambda (CFMC):                 %18.10g [s]\n",Components[i].CpuTimeCFSwapLambdaMove[CurrentSystem]);
    fprintf(FilePtr,"\tswap lambda (CB/CFMC):              %18.10g [s]\n",Components[i].CpuTimeCBCFSwapLambdaMove[CurrentSystem]);
    fprintf(FilePtr,"\tWidom:                              %18.10g [s]\n",Components[i].CpuTimeWidomMove[CurrentSystem]);
    fprintf(FilePtr,"\tCF-Widom:                           %18.10g [s]\n",Components[i].CpuTimeCFWidomLambdaMove[CurrentSystem]);
    fprintf(FilePtr,"\tGibbs Widom:                        %18.10g [s]\n",Components[i].CpuTimeGibbsWidomMove[CurrentSystem]);
    fprintf(FilePtr,"\tsurface area:                       %18.10g [s]\n",Components[i].CpuTimeSurfaceAreaMove[CurrentSystem]);
    fprintf(FilePtr,"\tGibbs particle transform:           %18.10g [s]\n",Components[i].CpuTimeGibbsChangeMove[CurrentSystem]);
    fprintf(FilePtr,"\tGibbs particle transform (CFMC):    %18.10g [s]\n",Components[i].CpuTimeCFGibbsChangeMove[CurrentSystem]);
    fprintf(FilePtr,"\tGibbs particle transform (CB/CFMC): %18.10g [s]\n",Components[i].CpuTimeCBCFGibbsChangeMove[CurrentSystem]);
    fprintf(FilePtr,"\tGibbs indentity change:             %18.10g [s]\n",Components[i].CpuTimeGibbsIdentityChangeMove[CurrentSystem]);
    fprintf(FilePtr,"\tExchange fract./int. particle:      %18.10g [s]\n",Components[i].CpuTimeExchangeFractionalParticleMove[CurrentSystem]);
    fprintf(FilePtr,"\tSwap Gibbs-fractional molecules:    %18.10g [s]\n",Components[i].CpuTimeCFGibbsSwapFractionalMoleculeToOtherBoxMove[CurrentSystem]);
    fprintf(FilePtr,"\tChange Gibs-lambda value:           %18.10g [s]\n",Components[i].CpuTimeCFGibbsLambdaChangeMove[CurrentSystem]);
    fprintf(FilePtr,"\tConvert Gibbs fract. to integer:    %18.10g [s]\n",Components[i].CpuTimeCFGibbsFractionalToIntegerMove[CurrentSystem]);

    CpuTimeTranslationMove+=Components[i].CpuTimeTranslationMove[CurrentSystem];
    CpuTimeRandomTranslationMove+=Components[i].CpuTimeRandomTranslationMove[CurrentSystem];
    CpuTimeRotationMove+=Components[i].CpuTimeRotationMove[CurrentSystem];
    CpuTimeRandomRotationMove+=Components[i].CpuTimeRandomRotationMove[CurrentSystem];
    CpuTimePartialReinsertionMove+=Components[i].CpuTimePartialReinsertionMove[CurrentSystem];
    CpuTimeReinsertionMove+=Components[i].CpuTimeReinsertionMove[CurrentSystem];
    CpuTimeReinsertionInPlaceMove+=Components[i].CpuTimeReinsertionInPlaceMove[CurrentSystem];
    CpuTimeReinsertionInPlaneMove+=Components[i].CpuTimeReinsertionInPlaneMove[CurrentSystem];
    CpuTimeIdentityChangeMove+=Components[i].CpuTimeIdentityChangeMove[CurrentSystem];
    CpuTimeSwapMoveInsertion+=Components[i].CpuTimeSwapMoveInsertion[CurrentSystem];
    CpuTimeSwapMoveDeletion+=Components[i].CpuTimeSwapMoveInsertion[CurrentSystem];
    CpuTimeCFSwapLambdaMove+=Components[i].CpuTimeCFSwapLambdaMove[CurrentSystem];
    CpuTimeCBCFSwapLambdaMove+=Components[i].CpuTimeCBCFSwapLambdaMove[CurrentSystem];
    CpuTimeWidomMove+=Components[i].CpuTimeWidomMove[CurrentSystem];
    CpuTimeCFWidomLambdaMove+=Components[i].CpuTimeCFWidomLambdaMove[CurrentSystem];
    CpuTimeGibbsWidomMove+=Components[i].CpuTimeGibbsWidomMove[CurrentSystem];
    CpuTimeSurfaceAreaMove+=Components[i].CpuTimeSurfaceAreaMove[CurrentSystem];
    CpuTimeGibbsChangeMove+=Components[i].CpuTimeGibbsChangeMove[CurrentSystem];
    CpuTimeCFGibbsChangeMove+=Components[i].CpuTimeCFGibbsChangeMove[CurrentSystem];
    CpuTimeCBCFGibbsChangeMove+=Components[i].CpuTimeCBCFGibbsChangeMove[CurrentSystem];
    CpuTimeGibbsIdentityChangeMove+=Components[i].CpuTimeGibbsIdentityChangeMove[CurrentSystem];
    CpuTimeExchangeFractionalParticleMove+=Components[i].CpuTimeExchangeFractionalParticleMove[CurrentSystem];
    CpuTimeCFGibbsSwapFractionalMoleculeToOtherBoxMove+=Components[i].CpuTimeCFGibbsSwapFractionalMoleculeToOtherBoxMove[CurrentSystem];
    CpuTimeCFGibbsLambdaChangeMove+=Components[i].CpuTimeCFGibbsLambdaChangeMove[CurrentSystem];
    CpuTimeCFGibbsFractionalToIntegerMove+=Components[i].CpuTimeCFGibbsFractionalToIntegerMove[CurrentSystem];
  }

  fprintf(FilePtr,"\nTotal all components:\n");
  fprintf(FilePtr,"\ttranslation:                        %18.10g [s]\n",CpuTimeTranslationMove);
  fprintf(FilePtr,"\trandom translation:                 %18.10g [s]\n",CpuTimeRandomTranslationMove);
  fprintf(FilePtr,"\trotation:                           %18.10g [s]\n",CpuTimeRotationMove);
  fprintf(FilePtr,"\trandom rotation:                    %18.10g [s]\n",CpuTimeRandomRotationMove);
  fprintf(FilePtr,"\tpartial reinsertion:                %18.10g [s]\n",CpuTimePartialReinsertionMove);
  fprintf(FilePtr,"\treinsertion:                        %18.10g [s]\n",CpuTimeReinsertionMove);
  fprintf(FilePtr,"\treinsertion in-place:               %18.10g [s]\n",CpuTimeReinsertionInPlaceMove);
  fprintf(FilePtr,"\treinsertion in-plane:               %18.10g [s]\n",CpuTimeReinsertionInPlaneMove);
  fprintf(FilePtr,"\tidentity switch:                    %18.10g [s]\n",CpuTimeIdentityChangeMove);
  fprintf(FilePtr,"\tswap (insertion):                   %18.10g [s]\n",CpuTimeSwapMoveInsertion);
  fprintf(FilePtr,"\tswap (deletion):                    %18.10g [s]\n",CpuTimeSwapMoveDeletion);
  fprintf(FilePtr,"\tswap lambda (CFMC):                 %18.10g [s]\n",CpuTimeCFSwapLambdaMove);
  fprintf(FilePtr,"\tswap lambda (CB/CFMC):              %18.10g [s]\n",CpuTimeCBCFSwapLambdaMove);
  fprintf(FilePtr,"\tWidom:                              %18.10g [s]\n",CpuTimeWidomMove);
  fprintf(FilePtr,"\tCF-Widom:                           %18.10g [s]\n",CpuTimeCFWidomLambdaMove);
  fprintf(FilePtr,"\tGibbs Widom:                        %18.10g [s]\n",CpuTimeGibbsWidomMove);
  fprintf(FilePtr,"\tsurface area:                       %18.10g [s]\n",CpuTimeSurfaceAreaMove);
  fprintf(FilePtr,"\tGibbs particle transform:           %18.10g [s]\n",CpuTimeGibbsChangeMove);
  fprintf(FilePtr,"\tGibbs particle transform (CFMC):    %18.10g [s]\n",CpuTimeCFGibbsChangeMove);
  fprintf(FilePtr,"\tGibbs particle transform (CB/CFMC): %18.10g [s]\n",CpuTimeCBCFGibbsChangeMove);
  fprintf(FilePtr,"\tGibbs identity change:              %18.10g [s]\n",CpuTimeGibbsIdentityChangeMove);
  fprintf(FilePtr,"\tExchange fract./int. particle:      %18.10g [s]\n",CpuTimeExchangeFractionalParticleMove);
  fprintf(FilePtr,"\tSwap Gibbs-fractional molecules:    %18.10g [s]\n",CpuTimeCFGibbsSwapFractionalMoleculeToOtherBoxMove);
  fprintf(FilePtr,"\tChange Gibs-lambda value:           %18.10g [s]\n",CpuTimeCFGibbsLambdaChangeMove);
  fprintf(FilePtr,"\tConvert Gibbs fract. to integer:    %18.10g [s]\n",CpuTimeCFGibbsFractionalToIntegerMove);

  fprintf(FilePtr,"\nSystem moves:\n");
  fprintf(FilePtr,"\tparallel tempering:            %18.10g [s]\n",CpuTimeParallelTemperingMove[CurrentSystem]);
  fprintf(FilePtr,"\thyper parallel tempering:      %18.10g [s]\n",CpuTimeHyperParallelTemperingMove[CurrentSystem]);
  fprintf(FilePtr,"\tmol-fraction replica-exchange: %18.10g [s]\n",CpuTimeParallelMolFractionMove[CurrentSystem]);
  fprintf(FilePtr,"\tchiral inversion:              %18.10g [s]\n",CpuTimeChiralInversionMove[CurrentSystem]);
  fprintf(FilePtr,"\thybrid MC/MD (NVE):            %18.10g [s]\n",CpuTimeHybridNVEMove[CurrentSystem]);
  fprintf(FilePtr,"\thybrid MC/MD (NPH):            %18.10g [s]\n",CpuTimeHybridNPHMove[CurrentSystem]);
  fprintf(FilePtr,"\thybrid MC/MD (NPHPR):          %18.10g [s]\n",CpuTimeHybridNPHPRMove[CurrentSystem]);
  fprintf(FilePtr,"\tvolume change:                 %18.10g [s]\n",CpuTimeVolumeChangeMove[CurrentSystem]);
  fprintf(FilePtr,"\tbox change:                    %18.10g [s]\n",CpuTimeBoxShapeChangeMove[CurrentSystem]);
  fprintf(FilePtr,"\tGibbs volume change:           %18.10g [s]\n",CpuTimeGibbsVolumeChangeMove[CurrentSystem]);
  fprintf(FilePtr,"\tframework change:              %18.10g [s]\n",CpuTimeFrameworkChangeMove[CurrentSystem]);
  fprintf(FilePtr,"\tframework shift:               %18.10g [s]\n",CpuTimeFrameworkShiftMove[CurrentSystem]);
  fprintf(FilePtr,"\treaction MC move:              %18.10g [s]\n",CpuTimeCFCRXMCLambdaChangeMove[CurrentSystem]);
  fprintf(FilePtr,"\n");


  fprintf(FilePtr,"Production run CPU timings of the MC moves summed over all systems and components:\n");
  fprintf(FilePtr,"==================================================================================\n");

  CpuTimeTranslationMove=0.0;
  CpuTimeRandomTranslationMove=0.0;
  CpuTimeRotationMove=0.0;
  CpuTimeRandomRotationMove=0.0;
  CpuTimePartialReinsertionMove=0.0;
  CpuTimeReinsertionMove=0.0;
  CpuTimeReinsertionInPlaceMove=0.0;
  CpuTimeReinsertionInPlaneMove=0.0;
  CpuTimeIdentityChangeMove=0.0;
  CpuTimeSwapMoveInsertion=0.0;
  CpuTimeSwapMoveDeletion=0.0;
  CpuTimeCFSwapLambdaMove=0.0;
  CpuTimeCBCFSwapLambdaMove=0.0;
  CpuTimeWidomMove=0.0;
  CpuTimeCFWidomLambdaMove=0.0;
  CpuTimeGibbsWidomMove=0.0;
  CpuTimeSurfaceAreaMove=0.0;
  CpuTimeGibbsChangeMove=0.0;
  CpuTimeCFGibbsChangeMove=0.0;
  CpuTimeCBCFGibbsChangeMove=0.0;
  CpuTimeGibbsIdentityChangeMove=0.0;
  CpuTimeExchangeFractionalParticleMove=0.0;
  CpuTimeCFGibbsSwapFractionalMoleculeToOtherBoxMove=0.0;
  CpuTimeCFGibbsLambdaChangeMove=0.0;
  CpuTimeCFGibbsFractionalToIntegerMove=0.0;

  for(j=0;j<NumberOfSystems;j++)
  {
    for(i=0;i<NumberOfComponents;i++)
    {
      CpuTimeTranslationMove+=Components[i].CpuTimeTranslationMove[j];
      CpuTimeRandomTranslationMove+=Components[i].CpuTimeRandomTranslationMove[j];
      CpuTimeRotationMove+=Components[i].CpuTimeRotationMove[j];
      CpuTimeRandomRotationMove+=Components[i].CpuTimeRandomRotationMove[j];
      CpuTimePartialReinsertionMove+=Components[i].CpuTimePartialReinsertionMove[j];
      CpuTimeReinsertionMove+=Components[i].CpuTimeReinsertionMove[j];
      CpuTimeReinsertionInPlaceMove+=Components[i].CpuTimeReinsertionInPlaceMove[j];
      CpuTimeReinsertionInPlaneMove+=Components[i].CpuTimeReinsertionInPlaneMove[j];
      CpuTimeIdentityChangeMove+=Components[i].CpuTimeIdentityChangeMove[j];
      CpuTimeSwapMoveInsertion+=Components[i].CpuTimeSwapMoveInsertion[j];
      CpuTimeSwapMoveDeletion+=Components[i].CpuTimeSwapMoveInsertion[j];
      CpuTimeCFSwapLambdaMove+=Components[i].CpuTimeCFSwapLambdaMove[j];
      CpuTimeCBCFSwapLambdaMove+=Components[i].CpuTimeCBCFSwapLambdaMove[j];
      CpuTimeWidomMove+=Components[i].CpuTimeWidomMove[j];
      CpuTimeCFWidomLambdaMove+=Components[i].CpuTimeCFWidomLambdaMove[j];
      CpuTimeGibbsWidomMove+=Components[i].CpuTimeGibbsWidomMove[j];
      CpuTimeSurfaceAreaMove+=Components[i].CpuTimeSurfaceAreaMove[j];
      CpuTimeGibbsChangeMove+=Components[i].CpuTimeGibbsChangeMove[j];
      CpuTimeCFGibbsChangeMove+=Components[i].CpuTimeCFGibbsChangeMove[j];
      CpuTimeCBCFGibbsChangeMove+=Components[i].CpuTimeCBCFGibbsChangeMove[j];
      CpuTimeGibbsIdentityChangeMove+=Components[i].CpuTimeGibbsIdentityChangeMove[j];
      CpuTimeExchangeFractionalParticleMove+=Components[i].CpuTimeExchangeFractionalParticleMove[j];
      CpuTimeCFGibbsSwapFractionalMoleculeToOtherBoxMove+=Components[i].CpuTimeCFGibbsSwapFractionalMoleculeToOtherBoxMove[j];
      CpuTimeCFGibbsLambdaChangeMove+=Components[i].CpuTimeCFGibbsLambdaChangeMove[j];
      CpuTimeCFGibbsFractionalToIntegerMove+=Components[i].CpuTimeCFGibbsFractionalToIntegerMove[j];
    }
  }

  fprintf(FilePtr,"\nParticles moves:\n");
  fprintf(FilePtr,"\ttranslation:                        %18.10g [s]\n",CpuTimeTranslationMove);
  fprintf(FilePtr,"\trandom translation:                 %18.10g [s]\n",CpuTimeRandomTranslationMove);
  fprintf(FilePtr,"\trotation:                           %18.10g [s]\n",CpuTimeRotationMove);
  fprintf(FilePtr,"\trandom rotation:                    %18.10g [s]\n",CpuTimeRandomRotationMove);
  fprintf(FilePtr,"\tpartial reinsertion:                %18.10g [s]\n",CpuTimePartialReinsertionMove);
  fprintf(FilePtr,"\treinsertion:                        %18.10g [s]\n",CpuTimeReinsertionMove);
  fprintf(FilePtr,"\treinsertion in-place:               %18.10g [s]\n",CpuTimeReinsertionInPlaceMove);
  fprintf(FilePtr,"\treinsertion in-plane:               %18.10g [s]\n",CpuTimeReinsertionInPlaneMove);
  fprintf(FilePtr,"\tidentity switch:                    %18.10g [s]\n",CpuTimeIdentityChangeMove);
  fprintf(FilePtr,"\tswap (insertion):                   %18.10g [s]\n",CpuTimeSwapMoveInsertion);
  fprintf(FilePtr,"\tswap (deletion):                    %18.10g [s]\n",CpuTimeSwapMoveDeletion);
  fprintf(FilePtr,"\tswap lambda (CFMC):                 %18.10g [s]\n",CpuTimeCFSwapLambdaMove);
  fprintf(FilePtr,"\tswap lambda (CB/CFMC):              %18.10g [s]\n",CpuTimeCBCFSwapLambdaMove);
  fprintf(FilePtr,"\tWidom:                              %18.10g [s]\n",CpuTimeWidomMove);
  fprintf(FilePtr,"\tCF-Widom:                           %18.10g [s]\n",CpuTimeCFWidomLambdaMove);
  fprintf(FilePtr,"\tGibbs Widom:                        %18.10g [s]\n",CpuTimeGibbsWidomMove);
  fprintf(FilePtr,"\tsurface area:                       %18.10g [s]\n",CpuTimeSurfaceAreaMove);
  fprintf(FilePtr,"\tGibbs particle transform:           %18.10g [s]\n",CpuTimeGibbsChangeMove);
  fprintf(FilePtr,"\tGibbs particle transform (CFMC):    %18.10g [s]\n",CpuTimeCFGibbsChangeMove);
  fprintf(FilePtr,"\tGibbs particle transform (CB/CFMC): %18.10g [s]\n",CpuTimeCBCFGibbsChangeMove);
  fprintf(FilePtr,"\tGibbs indentity change:             %18.10g [s]\n",CpuTimeGibbsIdentityChangeMove);
  fprintf(FilePtr,"\tExchange frac./int. particle:       %18.10g [s]\n",CpuTimeExchangeFractionalParticleMove);
  fprintf(FilePtr,"\tSwap Gibbs-fractional molecules:    %18.10g [s]\n",CpuTimeCFGibbsSwapFractionalMoleculeToOtherBoxMove);
  fprintf(FilePtr,"\tChange Gibs-lambda value:           %18.10g [s]\n",CpuTimeCFGibbsLambdaChangeMove);
  fprintf(FilePtr,"\tConvert Gibbs fract. to integer:    %18.10g [s]\n",CpuTimeCFGibbsFractionalToIntegerMove);

  CpuTimeParallelTemperingMoveTotal=0.0;
  CpuTimeHyperParallelTemperingMoveTotal=0.0;
  CpuTimeParallelMolFractionMoveTotal=0.0;
  CpuTimeChiralInversionMoveTotal=0.0;
  CpuTimeHybridNVEMoveTotal=0.0;
  CpuTimeHybridNPHMoveTotal=0.0;
  CpuTimeHybridNPHPRMoveTotal=0.0;
  CpuTimeVolumeChangeMoveTotal=0.0;
  CpuTimeBoxShapeChangeMoveTotal=0.0;
  CpuTimeGibbsVolumeChangeMoveTotal=0.0;
  CpuTimeFrameworkChangeMoveTotal=0.0;
  CpuTimeFrameworkShiftMoveTotal=0.0;
  CpuTimeCFCRXMCLambdaChangeMoveTotal=0.0;
  for(j=0;j<NumberOfSystems;j++)
  {
    CpuTimeParallelTemperingMoveTotal+=CpuTimeParallelTemperingMove[j];
    CpuTimeHyperParallelTemperingMoveTotal+=CpuTimeHyperParallelTemperingMove[j];
    CpuTimeParallelMolFractionMoveTotal+=CpuTimeParallelMolFractionMove[j];
    CpuTimeChiralInversionMoveTotal+=CpuTimeChiralInversionMove[j];
    CpuTimeHybridNVEMoveTotal+=CpuTimeHybridNVEMove[j];
    CpuTimeHybridNPHMoveTotal+=CpuTimeHybridNPHMove[j];
    CpuTimeHybridNPHPRMoveTotal+=CpuTimeHybridNPHPRMove[j];
    CpuTimeVolumeChangeMoveTotal+=CpuTimeVolumeChangeMove[j];
    CpuTimeBoxShapeChangeMoveTotal+=CpuTimeBoxShapeChangeMove[j];
    CpuTimeGibbsVolumeChangeMoveTotal+=CpuTimeGibbsVolumeChangeMove[j];
    CpuTimeFrameworkChangeMoveTotal+=CpuTimeFrameworkChangeMove[j];
    CpuTimeFrameworkShiftMoveTotal+=CpuTimeFrameworkShiftMove[j];
    CpuTimeCFCRXMCLambdaChangeMoveTotal+=CpuTimeCFCRXMCLambdaChangeMove[j];
  }

  fprintf(FilePtr,"\nSystem moves:\n");
  fprintf(FilePtr,"\tparallel tempering:            %18.10g [s]\n",CpuTimeParallelTemperingMoveTotal);
  fprintf(FilePtr,"\thyper parallel tempering:      %18.10g [s]\n",CpuTimeHyperParallelTemperingMoveTotal);
  fprintf(FilePtr,"\tmol-fraction replica-exchange: %18.10g [s]\n",CpuTimeParallelMolFractionMoveTotal);
  fprintf(FilePtr,"\tchiral inversion:              %18.10g [s]\n",CpuTimeChiralInversionMoveTotal);
  fprintf(FilePtr,"\thybrid MC/MD (NVE):            %18.10g [s]\n",CpuTimeHybridNVEMoveTotal);
  fprintf(FilePtr,"\thybrid MC/MD (NPH):            %18.10g [s]\n",CpuTimeHybridNPHMoveTotal);
  fprintf(FilePtr,"\thybrid MC/MD (NPHPR):          %18.10g [s]\n",CpuTimeHybridNPHPRMoveTotal);
  fprintf(FilePtr,"\tvolume change:                 %18.10g [s]\n",CpuTimeVolumeChangeMoveTotal);
  fprintf(FilePtr,"\tbox change:                    %18.10g [s]\n",CpuTimeBoxShapeChangeMoveTotal);
  fprintf(FilePtr,"\tGibbs volume change:           %18.10g [s]\n",CpuTimeGibbsVolumeChangeMoveTotal);
  fprintf(FilePtr,"\tframework change:              %18.10g [s]\n",CpuTimeFrameworkChangeMoveTotal);
  fprintf(FilePtr,"\tframework shift:               %18.10g [s]\n",CpuTimeFrameworkShiftMoveTotal);
  fprintf(FilePtr,"\treaction MC move:              %18.10g [s]\n",CpuTimeCFCRXMCLambdaChangeMoveTotal);
  fprintf(FilePtr,"\n");
}


static int versionNumberPseudoAtoms=1;

void WriteRestartPseudoAtoms(FILE *FilePtr)
{
  int i;
  REAL Check;

  fwrite(&versionNumberPseudoAtoms,sizeof(int),1,FilePtr);
  fwrite(&ShiftPotentials,sizeof(int),1,FilePtr);
  fwrite(&NumberOfPseudoAtoms,sizeof(NumberOfPseudoAtoms),1,FilePtr);

  fwrite(&IndividualMixingRules,sizeof(int),1,FilePtr);
  fwrite(&IndividualInteractions,sizeof(int),1,FilePtr);
  fwrite(&GeneralMixingRule,sizeof(int),1,FilePtr);

  fwrite(PseudoAtoms,sizeof(PSEUDO_ATOM),NumberOfPseudoAtoms,FilePtr);
  for(i=0;i<NumberOfPseudoAtoms;i++)
  {
    fwrite(PotentialParms[i],sizeof(REAL[MAX_NUMBER_OF_POTENTIAL_ARGUMENTS]),NumberOfPseudoAtoms,FilePtr);
    fwrite(PotentialType[i],sizeof(int),NumberOfPseudoAtoms,FilePtr);
    fwrite(TailCorrection[i],sizeof(int),NumberOfPseudoAtoms,FilePtr);
    fwrite(ShiftPotential[i],sizeof(int),NumberOfPseudoAtoms,FilePtr);
  }
  for(i=0;i<NumberOfSystems;i++)
  {
    fwrite(NumberOfPseudoAtomsCount[i],sizeof(int),NumberOfPseudoAtoms,FilePtr);
    fwrite(NumberOfPseudoAtomsType[i],sizeof(int),NumberOfPseudoAtoms,FilePtr);
    fwrite(NumberOfFractionalPseudoAtomsType[i],sizeof(int),NumberOfPseudoAtoms,FilePtr);
  }
  fwrite(MapPseudoAtom,sizeof(int),NumberOfPseudoAtoms,FilePtr);

  fwrite(SwitchingVDWFactors3,sizeof(REAL),4,FilePtr);
  fwrite(SwitchingVDWFactors5,sizeof(REAL),6,FilePtr);
  fwrite(SwitchingVDWFactors7,sizeof(REAL),8,FilePtr);



  Check=123456789.0;
  fwrite(&Check,1,sizeof(REAL),FilePtr);
}

void ReadRestartPseudoAtoms(FILE *FilePtr)
{
  int i;
  REAL Check;
  int readversionNumber=0;

  fread(&readversionNumber,sizeof(int),1,FilePtr);
  if(readversionNumber > versionNumberPseudoAtoms)
  {
    fprintf(stderr,"Upgrade to last version of RASPA to read binary restart-file");
    exit(-1);
  }

  fread(&ShiftPotentials,sizeof(int),1,FilePtr);
  fread(&NumberOfPseudoAtoms,sizeof(NumberOfPseudoAtoms),1,FilePtr);

  fread(&IndividualMixingRules,sizeof(int),1,FilePtr);
  fread(&IndividualInteractions,sizeof(int),1,FilePtr);
  fread(&GeneralMixingRule,sizeof(int),1,FilePtr);

  PseudoAtoms=(PSEUDO_ATOM*)calloc(NumberOfPseudoAtoms,sizeof(PSEUDO_ATOM));
  PotentialParms=(REAL(**)[MAX_NUMBER_OF_POTENTIAL_ARGUMENTS])calloc(NumberOfPseudoAtoms,sizeof(REAL(*)[MAX_NUMBER_OF_POTENTIAL_ARGUMENTS]));
  PotentialType=(int**)calloc(NumberOfPseudoAtoms,sizeof(int*));
  TailCorrection=(int**)calloc(NumberOfPseudoAtoms,sizeof(int*));
  ShiftPotential=(int**)calloc(NumberOfPseudoAtoms,sizeof(int*));
  NumberOfPseudoAtomsCount=(int**)calloc(NumberOfSystems,sizeof(int*));
  NumberOfPseudoAtomsType=(int**)calloc(NumberOfSystems,sizeof(int*));
  NumberOfFractionalPseudoAtomsType=(int**)calloc(NumberOfSystems,sizeof(int*));
  NumberOfPseudoAtomsTypeNew=(int*)calloc(NumberOfPseudoAtoms,sizeof(int));
  NumberOfPseudoAtomsTypeOld=(int*)calloc(NumberOfPseudoAtoms,sizeof(int));
  MapPseudoAtom=(int*)calloc(NumberOfPseudoAtoms,sizeof(int));
  for(i=0;i<NumberOfPseudoAtoms;i++)
  {
    PotentialParms[i]=(REAL(*)[MAX_NUMBER_OF_POTENTIAL_ARGUMENTS])calloc(NumberOfPseudoAtoms,sizeof(REAL[MAX_NUMBER_OF_POTENTIAL_ARGUMENTS]));
    PotentialType[i]=(int*)calloc(NumberOfPseudoAtoms,sizeof(int));
    TailCorrection[i]=(int*)calloc(NumberOfPseudoAtoms,sizeof(int));
    ShiftPotential[i]=(int*)calloc(NumberOfPseudoAtoms,sizeof(int));
  }
  for(i=0;i<NumberOfSystems;i++)
  {
    NumberOfPseudoAtomsCount[i]=(int*)calloc(NumberOfPseudoAtoms,sizeof(int));
    NumberOfPseudoAtomsType[i]=(int*)calloc(NumberOfPseudoAtoms,sizeof(int));
    NumberOfFractionalPseudoAtomsType[i]=(int*)calloc(NumberOfPseudoAtoms,sizeof(int));
  }

  fread(PseudoAtoms,sizeof(PSEUDO_ATOM),NumberOfPseudoAtoms,FilePtr);
  for(i=0;i<NumberOfPseudoAtoms;i++)
  {
    fread(PotentialParms[i],sizeof(REAL[MAX_NUMBER_OF_POTENTIAL_ARGUMENTS]),NumberOfPseudoAtoms,FilePtr);
    fread(PotentialType[i],sizeof(int),NumberOfPseudoAtoms,FilePtr);
    fread(TailCorrection[i],sizeof(int),NumberOfPseudoAtoms,FilePtr);
    fread(ShiftPotential[i],sizeof(int),NumberOfPseudoAtoms,FilePtr);
  }
  for(i=0;i<NumberOfSystems;i++)
  {
    fread(NumberOfPseudoAtomsCount[i],sizeof(int),NumberOfPseudoAtoms,FilePtr);
    fread(NumberOfPseudoAtomsType[i],sizeof(int),NumberOfPseudoAtoms,FilePtr);
    fread(NumberOfFractionalPseudoAtomsType[i],sizeof(int),NumberOfPseudoAtoms,FilePtr);
  }
  fread(MapPseudoAtom,sizeof(int),NumberOfPseudoAtoms,FilePtr);

  fread(SwitchingVDWFactors3,sizeof(REAL),4,FilePtr);
  fread(SwitchingVDWFactors5,sizeof(REAL),6,FilePtr);
  fread(SwitchingVDWFactors7,sizeof(REAL),8,FilePtr);

  fread(&Check,1,sizeof(REAL),FilePtr);
  if(fabs(Check-123456789.0)>1e-10)
  {
    fprintf(stderr, "Error in binary restart-file (ReadRestartPseudoAtoms)\n");
    ContinueAfterCrash=FALSE;
  }
}

static int versionNumberMolecules=1;
void WriteRestartMolecules(FILE *FilePtr)
{
  int i,j;
  int Type;
  REAL Check;

  fwrite(&versionNumberMolecules,sizeof(int),1,FilePtr);

  fwrite(&NumberOfSystems,sizeof(int),1,FilePtr);
  fwrite(&CurrentSystem,sizeof(int),1,FilePtr);
  fwrite(&CurrentComponent,sizeof(int),1,FilePtr);

  fwrite(MaxNumberOfAdsorbateMolecules,sizeof(int),NumberOfSystems,FilePtr);
  fwrite(NumberOfAdsorbateMolecules,sizeof(int),NumberOfSystems,FilePtr);

  fwrite(MaxNumberOfCationMolecules,NumberOfSystems,sizeof(int),FilePtr);
  fwrite(NumberOfCationMolecules,NumberOfSystems,sizeof(int),FilePtr);

  fwrite(NumberOfFractionalMolecules,sizeof(int),NumberOfSystems,FilePtr);
  fwrite(NumberOfFractionalAdsorbateMolecules,sizeof(int),NumberOfSystems,FilePtr);
  fwrite(NumberOfFractionalCationMolecules,sizeof(int),NumberOfSystems,FilePtr);
  fwrite(NumberOfReactionMolecules,sizeof(int),NumberOfSystems,FilePtr);
  fwrite(NumberOfReactionAdsorbateMolecules,sizeof(int),NumberOfSystems,FilePtr);
  fwrite(NumberOfReactionCationMolecules,sizeof(int),NumberOfSystems,FilePtr);

  fwrite(&MaxNumberOfCoulombicSites,sizeof(int),1,FilePtr);
  fwrite(&LargestNumberOfCoulombicSites,sizeof(int),1,FilePtr);
  fwrite(&LargestNumberOfBondDipoleSites,sizeof(int),1,FilePtr);
  fwrite(NumberOfAtomsPerSystem,sizeof(int),NumberOfSystems,FilePtr);
  fwrite(NumberOfChargesPerSystem,sizeof(int),NumberOfSystems,FilePtr);
  fwrite(NumberOfBondDipolesPerSystem,sizeof(int),NumberOfSystems,FilePtr);

  for(i=0;i<NumberOfSystems;i++)
  {
    fwrite(Adsorbates[i],sizeof(ADSORBATE_MOLECULE),MaxNumberOfAdsorbateMolecules[i],FilePtr);
    fwrite(Cations[i],sizeof(CATION_MOLECULE),MaxNumberOfCationMolecules[i],FilePtr);
  }

  // read the adsorbates and cations atoms and groups
  for(i=0;i<NumberOfSystems;i++)
  {
    for(j=0;j<NumberOfAdsorbateMolecules[i];j++)
    {
      Type=Adsorbates[i][j].Type;
      fwrite(Adsorbates[i][j].Atoms,sizeof(ATOM),Components[Type].NumberOfAtoms,FilePtr);

      if(Components[Type].NumberOfGroups>0)
        fwrite(Adsorbates[i][j].Groups,sizeof(GROUP),Components[Type].NumberOfGroups,FilePtr);
    }

    for(j=0;j<NumberOfCationMolecules[i];j++)
    {
      Type=Cations[i][j].Type;
      fwrite(Cations[i][j].Atoms,sizeof(ATOM),Components[Type].NumberOfAtoms,FilePtr);

      if(Components[Type].NumberOfGroups>0)
        fwrite(Cations[i][j].Groups,sizeof(GROUP),Components[Type].NumberOfGroups,FilePtr);
    }
  }

  Check=123456789.0;
  fwrite(&Check,1,sizeof(REAL),FilePtr);
}

void ReadRestartMolecules(FILE *FilePtr)
{
  int i,j;
  int Type;
  REAL Check;
  int readversionNumber=0;

  fread(&readversionNumber,sizeof(int),1,FilePtr);
  if(readversionNumber > versionNumberMolecules)
  {
    fprintf(stderr,"Upgrade to last version of RASPA to read binary restart-file");
    exit(-1);
  }

  fread(&NumberOfSystems,sizeof(int),1,FilePtr);
  fread(&CurrentSystem,sizeof(int),1,FilePtr);
  fread(&CurrentComponent,sizeof(int),1,FilePtr);

  // allocate memory for the components elements
  MaxNumberOfAdsorbateMolecules=(int*)calloc(NumberOfSystems,sizeof(int));
  NumberOfAdsorbateMolecules=(int*)calloc(NumberOfSystems,sizeof(int));

  MaxNumberOfCationMolecules=(int*)calloc(NumberOfSystems,sizeof(int));
  NumberOfCationMolecules=(int*)calloc(NumberOfSystems,sizeof(int));

  fread(MaxNumberOfAdsorbateMolecules,sizeof(int),NumberOfSystems,FilePtr);
  fread(NumberOfAdsorbateMolecules,sizeof(int),NumberOfSystems,FilePtr);

  fread(MaxNumberOfCationMolecules,sizeof(int),NumberOfSystems,FilePtr);
  fread(NumberOfCationMolecules,sizeof(int),NumberOfSystems,FilePtr);

  NumberOfFractionalMolecules=(int*)calloc(NumberOfSystems,sizeof(int));
  NumberOfFractionalAdsorbateMolecules=(int*)calloc(NumberOfSystems,sizeof(int));
  NumberOfFractionalCationMolecules=(int*)calloc(NumberOfSystems,sizeof(int));
  fread(NumberOfFractionalMolecules,sizeof(int),NumberOfSystems,FilePtr);
  fread(NumberOfFractionalAdsorbateMolecules,sizeof(int),NumberOfSystems,FilePtr);
  fread(NumberOfFractionalCationMolecules,sizeof(int),NumberOfSystems,FilePtr);

  NumberOfReactionMolecules=(int*)calloc(NumberOfSystems,sizeof(int));
  NumberOfReactionAdsorbateMolecules=(int*)calloc(NumberOfSystems,sizeof(int));
  NumberOfReactionCationMolecules=(int*)calloc(NumberOfSystems,sizeof(int));
  fread(NumberOfReactionMolecules,sizeof(int),NumberOfSystems,FilePtr);
  fread(NumberOfReactionAdsorbateMolecules,sizeof(int),NumberOfSystems,FilePtr);
  fread(NumberOfReactionCationMolecules,sizeof(int),NumberOfSystems,FilePtr);

  fread(&MaxNumberOfCoulombicSites,sizeof(int),1,FilePtr);
  fread(&LargestNumberOfCoulombicSites,sizeof(int),1,FilePtr);
  fread(&LargestNumberOfBondDipoleSites,sizeof(int),1,FilePtr);

  NumberOfAtomsPerSystem=(int*)calloc(NumberOfSystems,sizeof(int));
  fread(NumberOfAtomsPerSystem,sizeof(int),NumberOfSystems,FilePtr);

  NumberOfChargesPerSystem=(int*)calloc(NumberOfSystems,sizeof(int));
  fread(NumberOfChargesPerSystem,sizeof(int),NumberOfSystems,FilePtr);

  NumberOfBondDipolesPerSystem=(int*)calloc(NumberOfSystems,sizeof(int));
  fread(NumberOfBondDipolesPerSystem,sizeof(int),NumberOfSystems,FilePtr);

  // allocate memory for the adsorbates and cations
  Adsorbates=(ADSORBATE_MOLECULE**)calloc(NumberOfSystems,sizeof(ADSORBATE_MOLECULE*));
  for(i=0;i<NumberOfSystems;i++)
    Adsorbates[i]=(ADSORBATE_MOLECULE*)calloc(MaxNumberOfAdsorbateMolecules[i],sizeof(ADSORBATE_MOLECULE));

  Cations=(CATION_MOLECULE**)calloc(NumberOfSystems,sizeof(CATION_MOLECULE*));
  for(i=0;i<NumberOfSystems;i++)
    Cations[i]=(CATION_MOLECULE*)calloc(MaxNumberOfCationMolecules[i],sizeof(CATION_MOLECULE));

  for(i=0;i<NumberOfSystems;i++)
  {
    fread(Adsorbates[i],sizeof(ADSORBATE_MOLECULE),MaxNumberOfAdsorbateMolecules[i],FilePtr);
    fread(Cations[i],sizeof(CATION_MOLECULE),MaxNumberOfCationMolecules[i],FilePtr);
  }

  // allocate memory for the adsorbates and cations atoms and groups
  for(i=0;i<NumberOfSystems;i++)
  {
    for(j=0;j<NumberOfAdsorbateMolecules[i];j++)
    {
      Type=Adsorbates[i][j].Type;
      Adsorbates[i][j].Atoms=(ATOM*)calloc(Components[Type].NumberOfAtoms,sizeof(ATOM));
      if(Components[Type].NumberOfGroups>0)
        Adsorbates[i][j].Groups=(GROUP*)calloc(Components[Type].NumberOfGroups,sizeof(GROUP));
    }
    for(j=0;j<NumberOfCationMolecules[i];j++)
    {
      Type=Cations[i][j].Type;
      Cations[i][j].Atoms=(ATOM*)calloc(Components[Type].NumberOfAtoms,sizeof(ATOM));
      if(Components[Type].NumberOfGroups>0)
        Cations[i][j].Groups=(GROUP*)calloc(Components[Type].NumberOfGroups,sizeof(GROUP));
    }

  }

  // read the adsorbates and cations atoms and groups
  for(i=0;i<NumberOfSystems;i++)
  {
    for(j=0;j<NumberOfAdsorbateMolecules[i];j++)
    {
      Type=Adsorbates[i][j].Type;
      fread(Adsorbates[i][j].Atoms,sizeof(ATOM),Components[Type].NumberOfAtoms,FilePtr);

      if(Components[Type].NumberOfGroups>0)
        fread(Adsorbates[i][j].Groups,sizeof(GROUP),Components[Type].NumberOfGroups,FilePtr);
    }
    for(j=0;j<NumberOfCationMolecules[i];j++)
    {
      Type=Cations[i][j].Type;
      fread(Cations[i][j].Atoms,sizeof(ATOM),Components[Type].NumberOfAtoms,FilePtr);

      if(Components[Type].NumberOfGroups>0)
        fread(Cations[i][j].Groups,sizeof(GROUP),Components[Type].NumberOfGroups,FilePtr);
    }
  }

  fread(&Check,1,sizeof(REAL),FilePtr);
  if(fabs(Check-123456789.0)>1e-10)
  {
    fprintf(stderr, "Error in binary restart-file (ReadRestartMolecules)\n");
    ContinueAfterCrash=FALSE;
  }
}

/*********************************************************************************************************
 * Name       | CheckChiralityMolecules                                                                  *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Checks the whether the molecules have the correct chirality.                             *
 * Parameters | -                                                                                        *
 * Used in    | monte-carlo.c                                                                            *
 *********************************************************************************************************/

void CheckChiralityMolecules(void)
{
  int j,k,m;
  int A,B,C,D,type;
  POINT posA,posB,posC,posD;
  VECTOR drA,drC,drD;
  REAL tmp;
  REAL chirality;
  REAL computed_chirality;

  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    type=Adsorbates[CurrentSystem][m].Type;
    for(j=0;j<Adsorbates[CurrentSystem][m].NumberOfAtoms;j++)
    {
      if(Components[type].Chirality[j])
      {
        chirality=Components[type].ChiralityType[j];
        A=Components[type].ChiralA[j];
        B=Components[type].ChiralB[j];
        C=Components[type].ChiralC[j];
        D=Components[type].ChiralD[j];

        posA=Adsorbates[CurrentSystem][m].Atoms[A].Position;
        posB=Adsorbates[CurrentSystem][m].Atoms[B].Position;
        posC=Adsorbates[CurrentSystem][m].Atoms[C].Position;
        posD=Adsorbates[CurrentSystem][m].Atoms[D].Position;

        drA.x=posA.x-posB.x;
        drA.y=posA.y-posB.y;
        drA.z=posA.z-posB.z;

        drC.x=posC.x-posB.x;
        drC.y=posC.y-posB.y;
        drC.z=posC.z-posB.z;

        drD.x=posD.x-posB.x;
        drD.y=posD.y-posB.y;
        drD.z=posD.z-posB.z;

        // calculate (D x C) . A    if > 0 then R else S
        tmp=drA.x*(drD.y*drC.z-drD.z*drC.y)+drA.y*(drD.z*drC.x-drD.x*drC.z)+drA.z*(drD.x*drC.y-drD.y*drC.x);
        computed_chirality=tmp>0.0?R_CHIRAL:S_CHIRAL;

        if(computed_chirality!=chirality)
        {
          for(k=0;k<NumberOfSystems;k++)
          {
            fprintf(OutputFilePtr[k],"ERROR [CheckChiralityMolecules]: system %d adsorbate %d, center %d is %s but should be %s\n",CurrentSystem,
               m,j,computed_chirality==R_CHIRAL?"R":"S",chirality==R_CHIRAL?"R":"S");
            fclose(OutputFilePtr[k]);
          }
          exit(0);
        }
      }
    }
  }

  for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
  {
    type=Cations[CurrentSystem][m].Type;
    for(j=0;j<Cations[CurrentSystem][m].NumberOfAtoms;j++)
    {
      if(Components[type].Chirality[j])
      {
        chirality=Components[type].ChiralityType[j];
        A=Components[type].ChiralA[j];
        B=Components[type].ChiralB[j];
        C=Components[type].ChiralC[j];
        D=Components[type].ChiralD[j];

        posA=Cations[CurrentSystem][m].Atoms[A].Position;
        posB=Cations[CurrentSystem][m].Atoms[B].Position;
        posC=Cations[CurrentSystem][m].Atoms[C].Position;
        posD=Cations[CurrentSystem][m].Atoms[D].Position;

        drA.x=posA.x-posB.x;
        drA.y=posA.y-posB.y;
        drA.z=posA.z-posB.z;

        drC.x=posC.x-posB.x;
        drC.y=posC.y-posB.y;
        drC.z=posC.z-posB.z;

        drD.x=posD.x-posB.x;
        drD.y=posD.y-posB.y;
        drD.z=posD.z-posB.z;

        // calculate (D x C) . A    if > 0 then R else S
        tmp=drA.x*(drD.y*drC.z-drD.z*drC.y)+drA.y*(drD.z*drC.x-drD.x*drC.z)+drA.z*(drD.x*drC.y-drD.y*drC.x);
        computed_chirality=tmp>0.0?R_CHIRAL:S_CHIRAL;

        if(computed_chirality!=chirality)
        {
          for(k=0;k<NumberOfSystems;k++)
          {
            fprintf(OutputFilePtr[k],"ERROR [CheckChiralityMolecules]: system %d cation %d, center %d is %s but should be %s\n",CurrentSystem,
               m,j,computed_chirality==R_CHIRAL?"R":"S",chirality==R_CHIRAL?"R":"S");
            fclose(OutputFilePtr[k]);
          }
          exit(0);
        }
      }
    }
  }
}

static int versionNumberComponent=1;

void WriteRestartComponent(FILE *FilePtr)
{
  int i,j,k,n;
  REAL Check;

  fwrite(&versionNumberComponent,sizeof(int),1,FilePtr);

  fwrite(&NumberOfSystems,sizeof(int),1,FilePtr);
  fwrite(&NumberOfComponents,sizeof(int),1,FilePtr);
  fwrite(&NumberOfAdsorbateComponents,sizeof(int),1,FilePtr);
  fwrite(&NumberOfCationComponents,sizeof(int),1,FilePtr);
  fwrite(Components,sizeof(COMPONENT),NumberOfComponents,FilePtr);

  for(i=0;i<NumberOfReactions;i++)
  {
    fwrite(ReactantsStoichiometry[i],sizeof(int),NumberOfComponents,FilePtr);
    fwrite(ProductsStoichiometry[i],sizeof(int),NumberOfComponents,FilePtr);
  }

  for(i=0;i<NumberOfComponents;i++)
  {
    fwrite(Components[i].IdealGasRosenbluthWeight,sizeof(REAL),NumberOfSystems,FilePtr);
    fwrite(Components[i].IdealGasTotalEnergy,sizeof(REAL),NumberOfSystems,FilePtr);
    fwrite(Components[i].PartialPressure,sizeof(REAL),NumberOfSystems,FilePtr);
    fwrite(Components[i].FugacityCoefficient,sizeof(REAL),NumberOfSystems,FilePtr);
    fwrite(Components[i].BulkFluidDensity,sizeof(REAL),NumberOfSystems,FilePtr);
    fwrite(Components[i].Compressibility,sizeof(REAL),NumberOfSystems,FilePtr);
    fwrite(Components[i].MolFraction,sizeof(REAL),NumberOfSystems,FilePtr);
    fwrite(Components[i].AmountOfExcessMolecules,sizeof(REAL),NumberOfSystems,FilePtr);

    fwrite(Components[i].CreateNumberOfMolecules,sizeof(int),NumberOfSystems,FilePtr);
    fwrite(Components[i].NumberOfMolecules,sizeof(int),NumberOfSystems,FilePtr);

    fwrite(Components[i].MOLEC_PER_UC_TO_MOL_PER_KG,sizeof(REAL),NumberOfSystems,FilePtr);
    fwrite(Components[i].MOLEC_PER_UC_TO_GRAM_PER_GRAM_OF_FRAMEWORK,sizeof(REAL),NumberOfSystems,FilePtr);
    fwrite(Components[i].MOLEC_PER_UC_TO_CC_STP_G,sizeof(REAL),NumberOfSystems,FilePtr);
    fwrite(Components[i].MOLEC_PER_UC_TO_CC_STP_CC,sizeof(REAL),NumberOfSystems,FilePtr);
    fwrite(Components[i].MOL_PER_KG_TO_CC_STP_G,sizeof(REAL),NumberOfSystems,FilePtr);
    fwrite(Components[i].MOL_PER_KG_TO_CC_STP_CC,sizeof(REAL),NumberOfSystems,FilePtr);

    fwrite(Components[i].IdentityChanges,sizeof(int),NumberOfComponents,FilePtr);
    fwrite(Components[i].GibbsIdentityChanges,sizeof(int),NumberOfComponents,FilePtr);

    fwrite(Components[i].BlockPockets,sizeof(int),NumberOfSystems,FilePtr);
    fwrite(Components[i].ComputeFreeEnergyProfile,sizeof(int),NumberOfSystems,FilePtr);
    fwrite(Components[i].BlockPocketsFilename,sizeof(char[256]),NumberOfSystems,FilePtr);
    fwrite(Components[i].NumberOfBlockCenters,sizeof(int),NumberOfSystems,FilePtr);

    if(Components[i].Biased)
    {
      n=Components[i].BiasingFunction.n;
      fwrite(Components[i].BiasingFunction.a,sizeof(REAL),n,FilePtr);
      fwrite(Components[i].BiasingFunction.b,sizeof(REAL),n,FilePtr);
      fwrite(Components[i].BiasingFunction.c,sizeof(REAL),n,FilePtr);
      fwrite(Components[i].BiasingFunction.d,sizeof(REAL),n,FilePtr);
      fwrite(Components[i].BiasingFunction.x,sizeof(REAL),n,FilePtr);
      fwrite(Components[i].BiasingFunction.y,sizeof(REAL),n,FilePtr);
      fwrite(Components[i].FreeEnergyXdata,sizeof(REAL),n+1,FilePtr);
    }

    fwrite(Components[i].Fixed,sizeof(int),Components[i].NumberOfAtoms,FilePtr);
    fwrite(Components[i].Type,sizeof(int),Components[i].NumberOfAtoms,FilePtr);
    fwrite(Components[i].Charge,sizeof(REAL),Components[i].NumberOfAtoms,FilePtr);
    fwrite(Components[i].Connectivity,sizeof(int),Components[i].NumberOfAtoms,FilePtr);
    for(j=0;j<Components[i].NumberOfAtoms;j++)
      fwrite(Components[i].ConnectivityMatrix[j],sizeof(int),Components[i].NumberOfAtoms,FilePtr);
    for(j=0;j<Components[i].NumberOfAtoms;j++)
      fwrite(Components[i].ConnectivityList[j],sizeof(int),Components[i].Connectivity[j],FilePtr);

    fwrite(Components[i].Positions,sizeof(VECTOR),Components[i].NumberOfAtoms,FilePtr);
    if(Components[i].NumberOfGroups>0)
    {
      fwrite(Components[i].Groups,sizeof(GROUP_DEFINITION),Components[i].NumberOfGroups,FilePtr);
      for(j=0;j<Components[i].NumberOfGroups;j++)
        fwrite(Components[i].Groups[j].Atoms,sizeof(int),Components[i].Groups[j].NumberOfGroupAtoms,FilePtr);
    }
    fwrite(Components[i].group,sizeof(int),Components[i].NumberOfAtoms,FilePtr);

    fwrite(Components[i].RMCMOL,sizeof(VECTOR),Components[i].NumberOfAtoms,FilePtr);


    fwrite(Components[i].FractionalMolecule,sizeof(int),NumberOfSystems,FilePtr);
    fwrite(Components[i].CFMoleculePresent,sizeof(int),NumberOfSystems,FilePtr);
    fwrite(Components[i].RXMCMoleculesPresent,sizeof(int),NumberOfSystems,FilePtr);
    fwrite(Components[i].NumberOfRXMCMoleculesPresent,sizeof(int),NumberOfSystems,FilePtr);
    fwrite(Components[i].CFWangLandauScalingFactor,sizeof(REAL),NumberOfSystems,FilePtr);

    for(j=0;j<NumberOfSystems;j++)
    {
      fwrite(Components[i].BlockDistance[j],sizeof(REAL),Components[i].NumberOfBlockCenters[j],FilePtr);
      fwrite(Components[i].BlockCenters[j],sizeof(VECTOR),Components[i].NumberOfBlockCenters[j],FilePtr);

      fwrite(Components[i].MaximumCBMCChangeBondLength[j],sizeof(REAL),Components[i].NumberOfAtoms,FilePtr);
      fwrite(Components[i].MaximumCBMCChangeBendAngle[j],sizeof(REAL),Components[i].NumberOfAtoms,FilePtr);
      fwrite(Components[i].MaximumCBMCRotationOnCone[j],sizeof(REAL),Components[i].NumberOfAtoms,FilePtr);
      fwrite(Components[i].CBMCChangeBondLengthAttempts[j],sizeof(REAL),Components[i].NumberOfAtoms,FilePtr);
      fwrite(Components[i].CBMCChangeBondLengthAccepted[j],sizeof(REAL),Components[i].NumberOfAtoms,FilePtr);
      fwrite(Components[i].CBMCChangeBendAngleAttempts[j],sizeof(REAL),Components[i].NumberOfAtoms,FilePtr);
      fwrite(Components[i].CBMCChangeBendAngleAccepted[j],sizeof(REAL),Components[i].NumberOfAtoms,FilePtr);
      fwrite(Components[i].CBMCRotationOnConeAttempts[j],sizeof(REAL),Components[i].NumberOfAtoms,FilePtr);
      fwrite(Components[i].CBMCRotationOnConeAccepted[j],sizeof(REAL),Components[i].NumberOfAtoms,FilePtr);
      fwrite(Components[i].TotalCBMCChangeBondLengthAttempts[j],sizeof(REAL),Components[i].NumberOfAtoms,FilePtr);
      fwrite(Components[i].TotalCBMCChangeBendAngleAttempts[j],sizeof(REAL),Components[i].NumberOfAtoms,FilePtr);
      fwrite(Components[i].TotalCBMCRotationOnConeAttempts[j],sizeof(REAL),Components[i].NumberOfAtoms,FilePtr);
      fwrite(Components[i].TotalCBMCChangeBondLengthAccepted[j],sizeof(REAL),Components[i].NumberOfAtoms,FilePtr);
      fwrite(Components[i].TotalCBMCChangeBendAngleAccepted[j],sizeof(REAL),Components[i].NumberOfAtoms,FilePtr);
      fwrite(Components[i].TotalCBMCRotationOnConeAccepted[j],sizeof(REAL),Components[i].NumberOfAtoms,FilePtr);

      fwrite(Components[i].CFBiasingFactors[j],sizeof(REAL),Components[i].CFLambdaHistogramSize,FilePtr);

      if(NumberOfReactions>0)
      {
        for(k=0;k<NumberOfReactions;k++)
        {
          if(ReactantsStoichiometry[k][i]>0)
            fwrite(Components[i].ReactantFractionalMolecules[j][k],sizeof(int),ReactantsStoichiometry[k][i],FilePtr);
          if(ProductsStoichiometry[k][i]>0)
            fwrite(Components[i].ProductFractionalMolecules[j][k],sizeof(int),ProductsStoichiometry[k][i],FilePtr);
        }
      }
    }

    // allocate charility-centers
    fwrite(&Components[i].NumberOfChiralityCenters,sizeof(int),1,FilePtr);
    fwrite(Components[i].Chirality,sizeof(int),Components[i].NumberOfAtoms,FilePtr);
    if(Components[i].NumberOfChiralityCenters>0)
    {
      fwrite(Components[i].ChiralityType,sizeof(int),Components[i].NumberOfAtoms,FilePtr);
      fwrite(Components[i].ChiralA,sizeof(int),Components[i].NumberOfAtoms,FilePtr);
      fwrite(Components[i].ChiralB,sizeof(int),Components[i].NumberOfAtoms,FilePtr);
      fwrite(Components[i].ChiralC,sizeof(int),Components[i].NumberOfAtoms,FilePtr);
      fwrite(Components[i].ChiralD,sizeof(int),Components[i].NumberOfAtoms,FilePtr);
    }

    // write bonds
    fwrite(Components[i].Bonds,sizeof(PAIR),Components[i].NumberOfBonds,FilePtr);
    fwrite(Components[i].BondType,sizeof(int),Components[i].NumberOfBonds,FilePtr);
    fwrite(Components[i].BondArguments,sizeof(REAL[MAX_BOND_POTENTIAL_ARGUMENTS]),Components[i].NumberOfBonds,FilePtr);

    // write bond-dipoles
    fwrite(Components[i].BondDipoles,sizeof(PAIR),Components[i].NumberOfBondDipoles,FilePtr);
    fwrite(Components[i].BondDipoleMagnitude,sizeof(REAL),Components[i].NumberOfBondDipoles,FilePtr);

    // write bends
    fwrite(Components[i].Bends,sizeof(QUAD),Components[i].NumberOfBends,FilePtr);
    fwrite(Components[i].BendType,sizeof(int),Components[i].NumberOfBends,FilePtr);
    fwrite(Components[i].BendArguments,sizeof(REAL[MAX_BEND_POTENTIAL_ARGUMENTS]),Components[i].NumberOfBends,FilePtr);

    // write urey-bradleys
    fwrite(Components[i].UreyBradleys,sizeof(TRIPLE),Components[i].NumberOfUreyBradleys,FilePtr);
    fwrite(Components[i].UreyBradleyType,sizeof(int),Components[i].NumberOfUreyBradleys,FilePtr);
    fwrite(Components[i].UreyBradleyArguments,sizeof(REAL[MAX_UREYBRADLEY_POTENTIAL_ARGUMENTS]),Components[i].NumberOfUreyBradleys,FilePtr);

    // write inversion-bends
    fwrite(Components[i].InversionBends,sizeof(QUAD),Components[i].NumberOfInversionBends,FilePtr);
    fwrite(Components[i].InversionBendType,sizeof(int),Components[i].NumberOfInversionBends,FilePtr);
    fwrite(Components[i].InversionBendArguments,sizeof(REAL[MAX_INVERSION_BEND_POTENTIAL_ARGUMENTS]),Components[i].NumberOfInversionBends,FilePtr);

    // write torsions
    fwrite(Components[i].Torsions,sizeof(QUAD),Components[i].NumberOfTorsions,FilePtr);
    fwrite(Components[i].TorsionType,sizeof(int),Components[i].NumberOfTorsions,FilePtr);
    fwrite(Components[i].TorsionArguments,sizeof(REAL[MAX_TORSION_POTENTIAL_ARGUMENTS]),Components[i].NumberOfTorsions,FilePtr);

    // write improper torsions
    fwrite(Components[i].ImproperTorsions,sizeof(QUAD),Components[i].NumberOfImproperTorsions,FilePtr);
    fwrite(Components[i].ImproperTorsionType,sizeof(int),Components[i].NumberOfImproperTorsions,FilePtr);
    fwrite(Components[i].ImproperTorsionArguments,sizeof(REAL[MAX_IMPROPER_TORSION_POTENTIAL_ARGUMENTS]),Components[i].NumberOfImproperTorsions,FilePtr);

    // write out-of-plane
    fwrite(Components[i].OutOfPlanes,sizeof(QUAD),Components[i].NumberOfOutOfPlanes,FilePtr);
    fwrite(Components[i].OutOfPlaneType,sizeof(int),Components[i].NumberOfOutOfPlanes,FilePtr);
    fwrite(Components[i].OutOfPlaneArguments,sizeof(REAL[MAX_OUT_OF_PLANE_POTENTIAL_ARGUMENTS]),Components[i].NumberOfOutOfPlanes,FilePtr);

    // write bond-bonds
    fwrite(Components[i].BondBonds,sizeof(TRIPLE),Components[i].NumberOfBondBonds,FilePtr);
    fwrite(Components[i].BondBondType,sizeof(int),Components[i].NumberOfBondBonds,FilePtr);
    fwrite(Components[i].BondBondArguments,sizeof(REAL[MAX_BOND_BOND_POTENTIAL_ARGUMENTS]),Components[i].NumberOfBondBonds,FilePtr);

    // write bond-bends
    fwrite(Components[i].BondBends,sizeof(TRIPLE),Components[i].NumberOfBondBends,FilePtr);
    fwrite(Components[i].BondBendType,sizeof(int),Components[i].NumberOfBondBends,FilePtr);
    fwrite(Components[i].BondBendArguments,sizeof(REAL[MAX_BOND_BEND_POTENTIAL_ARGUMENTS]),Components[i].NumberOfBondBends,FilePtr);

    // write bend-bends
    fwrite(Components[i].BendBends,sizeof(QUAD),Components[i].NumberOfBendBends,FilePtr);
    fwrite(Components[i].BendBendType,sizeof(int),Components[i].NumberOfBendBends,FilePtr);
    fwrite(Components[i].BendBendArguments,sizeof(REAL[MAX_BEND_BEND_POTENTIAL_ARGUMENTS]),Components[i].NumberOfBendBends,FilePtr);

    // write stretch-torsions
    fwrite(Components[i].BondTorsions,sizeof(QUAD),Components[i].NumberOfBondTorsions,FilePtr);
    fwrite(Components[i].BondTorsionType,sizeof(int),Components[i].NumberOfBondTorsions,FilePtr);
    fwrite(Components[i].BondTorsionArguments,sizeof(REAL[MAX_BOND_TORSION_POTENTIAL_ARGUMENTS]),Components[i].NumberOfBondTorsions,FilePtr);

    // write bend-torsions
    fwrite(Components[i].BendTorsions,sizeof(QUAD),Components[i].NumberOfBendTorsions,FilePtr);
    fwrite(Components[i].BendTorsionType,sizeof(int),Components[i].NumberOfBendTorsions,FilePtr);
    fwrite(Components[i].BendTorsionArguments,sizeof(REAL[MAX_BEND_TORSION_POTENTIAL_ARGUMENTS]),Components[i].NumberOfBendTorsions,FilePtr);

    // write intra
    fwrite(Components[i].IntraVDW,sizeof(PAIR),Components[i].NumberOfIntraVDW,FilePtr);
    fwrite(Components[i].IntraChargeCharge,sizeof(PAIR),Components[i].NumberOfIntraChargeCharge,FilePtr);
    fwrite(Components[i].IntraChargeBondDipole,sizeof(PAIR),Components[i].NumberOfIntraChargeBondDipole,FilePtr);
    fwrite(Components[i].IntraBondDipoleBondDipole,sizeof(PAIR),Components[i].NumberOfIntraBondDipoleBondDipole,FilePtr);

    fwrite(Components[i].IntraVDWScaling,sizeof(REAL),Components[i].NumberOfIntraVDW,FilePtr);
    fwrite(Components[i].IntraChargeChargeScaling,sizeof(REAL),Components[i].NumberOfIntraChargeCharge,FilePtr);

    // write excluded-pairs
    fwrite(Components[i].ExcludedIntraChargeCharge,sizeof(PAIR),Components[i].NumberOfExcludedIntraChargeCharge,FilePtr);
    fwrite(Components[i].ExcludedIntraChargeBondDipole,sizeof(PAIR),Components[i].NumberOfExcludedIntraChargeBondDipole,FilePtr);
    fwrite(Components[i].ExcludedIntraBondDipoleBondDipole,sizeof(PAIR),Components[i].NumberOfExcludedIntraBondDipoleBondDipole,FilePtr);

    // write config moves
    fwrite(&Components[i].NumberOfConfigMoves,sizeof(int),1,FilePtr);
    fwrite(Components[i].NumberOfUnchangedAtomsConfig,sizeof(int),Components[i].NumberOfConfigMoves,FilePtr);
    for(j=0;j<Components[i].NumberOfConfigMoves;j++)
      fwrite(Components[i].UnchangedAtomsConfig[j],sizeof(int),Components[i].NumberOfUnchangedAtomsConfig[j],FilePtr);

    // write identity moves
    fwrite(&Components[i].NumberOfIdentityConfigMoves,sizeof(int),1,FilePtr);
    fwrite(Components[i].NumberOfUnchangedAtomsIdentityConfig,sizeof(int),Components[i].NumberOfIdentityConfigMoves,FilePtr);
    for(j=0;j<Components[i].NumberOfIdentityConfigMoves;j++)
      fwrite(Components[i].UnchangedAtomsIdentityConfig[j],sizeof(int),Components[i].NumberOfUnchangedAtomsIdentityConfig[j],FilePtr);

    fwrite(&Components[i].RestrictMoves,sizeof(int),1,FilePtr);
    fwrite(&Components[i].RestrictMovesToBox,sizeof(int),1,FilePtr);
    fwrite(&Components[i].BoxAxisABC_Min,sizeof(VECTOR),1,FilePtr);
    fwrite(&Components[i].BoxAxisABC_Max,sizeof(VECTOR),1,FilePtr);
    fwrite(&Components[i].BoxAxisABC_Min2,sizeof(VECTOR),1,FilePtr);
    fwrite(&Components[i].BoxAxisABC_Max2,sizeof(VECTOR),1,FilePtr);
    fwrite(&Components[i].BoxAxisABC_Min3,sizeof(VECTOR),1,FilePtr);
    fwrite(&Components[i].BoxAxisABC_Max3,sizeof(VECTOR),1,FilePtr);
    fwrite(&Components[i].BoxAxisABC_Min4,sizeof(VECTOR),1,FilePtr);
    fwrite(&Components[i].BoxAxisABC_Max4,sizeof(VECTOR),1,FilePtr);
    fwrite(&Components[i].RestrictMovesToPrisms,sizeof(int),1,FilePtr);
    fwrite(Components[i].RestrictMovesToPrism,sizeof(int),MAX_NUMBER_OF_PRISMS,FilePtr);
    fwrite(Components[i].RestrictPrismABC_Min,sizeof(VECTOR),MAX_NUMBER_OF_PRISMS,FilePtr);
    fwrite(Components[i].RestrictPrismABC_Max,sizeof(VECTOR),MAX_NUMBER_OF_PRISMS,FilePtr);
    fwrite(&Components[i].RestrictMovesToCylinders,sizeof(int),1,FilePtr);
    fwrite(Components[i].RestrictCylinderABC_Min,sizeof(VECTOR),MAX_NUMBER_OF_CYLINDERS,FilePtr);
    fwrite(Components[i].RestrictCylinderABC_Max,sizeof(VECTOR),MAX_NUMBER_OF_CYLINDERS,FilePtr);
    fwrite(Components[i].RestrictCylinderCenter,sizeof(VECTOR),MAX_NUMBER_OF_CYLINDERS,FilePtr);
    fwrite(Components[i].RestrictCylinderDirection,sizeof(int),MAX_NUMBER_OF_CYLINDERS,FilePtr);
    fwrite(Components[i].RestrictCylinderRadius,sizeof(REAL),MAX_NUMBER_OF_CYLINDERS,FilePtr);
    fwrite(&Components[i].RestrictMovesToSpheres,sizeof(int),1,FilePtr);
    fwrite(Components[i].RestrictMovesToSphere,sizeof(int),MAX_NUMBER_OF_SPHERES,FilePtr);
    fwrite(Components[i].RestrictSphereCenter,sizeof(VECTOR),MAX_NUMBER_OF_SPHERES,FilePtr);
    fwrite(Components[i].RestrictSphereRadius,sizeof(REAL),MAX_NUMBER_OF_SPHERES,FilePtr);

    fwrite(Components[i].CpuTimeTranslationMove,sizeof(REAL),NumberOfSystems,FilePtr);
    fwrite(Components[i].CpuTimeRandomTranslationMove,sizeof(REAL),NumberOfSystems,FilePtr);
    fwrite(Components[i].CpuTimeRotationMove,sizeof(REAL),NumberOfSystems,FilePtr);
    fwrite(Components[i].CpuTimeRandomRotationMove,sizeof(REAL),NumberOfSystems,FilePtr);
    fwrite(Components[i].CpuTimePartialReinsertionMove,sizeof(REAL),NumberOfSystems,FilePtr);
    fwrite(Components[i].CpuTimeReinsertionMove,sizeof(REAL),NumberOfSystems,FilePtr);
    fwrite(Components[i].CpuTimeReinsertionInPlaceMove,sizeof(REAL),NumberOfSystems,FilePtr);
    fwrite(Components[i].CpuTimeReinsertionInPlaneMove,sizeof(REAL),NumberOfSystems,FilePtr);
    fwrite(Components[i].CpuTimeIdentityChangeMove,sizeof(REAL),NumberOfSystems,FilePtr);
    fwrite(Components[i].CpuTimeSwapMoveInsertion,sizeof(REAL),NumberOfSystems,FilePtr);
    fwrite(Components[i].CpuTimeSwapMoveDeletion,sizeof(REAL),NumberOfSystems,FilePtr);
    fwrite(Components[i].CpuTimeCFSwapLambdaMove,sizeof(REAL),NumberOfSystems,FilePtr);
    fwrite(Components[i].CpuTimeCBCFSwapLambdaMove,sizeof(REAL),NumberOfSystems,FilePtr);
    fwrite(Components[i].CpuTimeWidomMove,sizeof(REAL),NumberOfSystems,FilePtr);
    fwrite(Components[i].CpuTimeCFWidomLambdaMove,sizeof(REAL),NumberOfSystems,FilePtr);
    fwrite(Components[i].CpuTimeGibbsWidomMove,sizeof(REAL),NumberOfSystems,FilePtr);
    fwrite(Components[i].CpuTimeSurfaceAreaMove,sizeof(REAL),NumberOfSystems,FilePtr);
    fwrite(Components[i].CpuTimeGibbsChangeMove,sizeof(REAL),NumberOfSystems,FilePtr);
    fwrite(Components[i].CpuTimeCFGibbsChangeMove,sizeof(REAL),NumberOfSystems,FilePtr);
    fwrite(Components[i].CpuTimeCBCFGibbsChangeMove,sizeof(REAL),NumberOfSystems,FilePtr);
    fwrite(Components[i].CpuTimeGibbsIdentityChangeMove,sizeof(REAL),NumberOfSystems,FilePtr);
    fwrite(Components[i].CpuTimeExchangeFractionalParticleMove,sizeof(REAL),NumberOfSystems,FilePtr);
    fwrite(Components[i].CpuTimeCFGibbsSwapFractionalMoleculeToOtherBoxMove,sizeof(REAL),NumberOfSystems,FilePtr);
    fwrite(Components[i].CpuTimeCFGibbsLambdaChangeMove,sizeof(REAL),NumberOfSystems,FilePtr);
    fwrite(Components[i].CpuTimeCFGibbsFractionalToIntegerMove,sizeof(REAL),NumberOfSystems,FilePtr);
  }


  Check=123456789.0;
  fwrite(&Check,1,sizeof(REAL),FilePtr);
}

void ReadRestartComponent(FILE *FilePtr)
{
  int i,j,k,n;
  REAL Check;
  int readversionNumber=0;

  fread(&readversionNumber,sizeof(int),1,FilePtr);
  if(readversionNumber > versionNumberComponent)
  {
    fprintf(stderr,"Upgrade to last version of RASPA to read binary restart-file");
    exit(-1);
  }

  fread(&NumberOfSystems,sizeof(int),1,FilePtr);
  fread(&NumberOfComponents,sizeof(int),1,FilePtr);
  fread(&NumberOfAdsorbateComponents,sizeof(int),1,FilePtr);
  fread(&NumberOfCationComponents,sizeof(int),1,FilePtr);

  // allocate memory for the components
  Components=(COMPONENT*)calloc(NumberOfComponents,sizeof(COMPONENT));
  fread(Components,sizeof(COMPONENT),NumberOfComponents,FilePtr);

  ReactantsStoichiometry=(int**)calloc(NumberOfReactions,sizeof(int*));
  ProductsStoichiometry=(int**)calloc(NumberOfReactions,sizeof(int*));
  for(i=0;i<NumberOfReactions;i++)
  {
    ReactantsStoichiometry[i]=(int*)calloc(NumberOfComponents,sizeof(int));
    fread(ReactantsStoichiometry[i],sizeof(int),NumberOfComponents,FilePtr);

    ProductsStoichiometry[i]=(int*)calloc(NumberOfComponents,sizeof(int));
    fread(ProductsStoichiometry[i],sizeof(int),NumberOfComponents,FilePtr);
  }


  for(i=0;i<NumberOfComponents;i++)
  {
    Components[i].IdealGasRosenbluthWeight=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
    Components[i].IdealGasTotalEnergy=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
    Components[i].PartialPressure=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
    Components[i].FugacityCoefficient=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
    Components[i].BulkFluidDensity=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
    Components[i].Compressibility=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
    Components[i].MolFraction=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
    Components[i].AmountOfExcessMolecules=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

    Components[i].CreateNumberOfMolecules=(int*)calloc(NumberOfSystems,sizeof(int));
    Components[i].NumberOfMolecules=(int*)calloc(NumberOfSystems,sizeof(int));

    Components[i].MOLEC_PER_UC_TO_MOL_PER_KG=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
    Components[i].MOLEC_PER_UC_TO_GRAM_PER_GRAM_OF_FRAMEWORK=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
    Components[i].MOLEC_PER_UC_TO_CC_STP_G=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
    Components[i].MOLEC_PER_UC_TO_CC_STP_CC=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
    Components[i].MOL_PER_KG_TO_CC_STP_G=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
    Components[i].MOL_PER_KG_TO_CC_STP_CC=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

    Components[i].IdentityChanges=(int*)calloc(NumberOfComponents,sizeof(int));
    Components[i].GibbsIdentityChanges=(int*)calloc(NumberOfComponents,sizeof(int));

    Components[i].BlockPockets=(int*)calloc(NumberOfSystems,sizeof(int));
    Components[i].ComputeFreeEnergyProfile=(int*)calloc(NumberOfSystems,sizeof(int));
    Components[i].BlockPocketsFilename=(char(*)[256])calloc(NumberOfSystems,sizeof(char[256]));
    Components[i].NumberOfBlockCenters=(int*)calloc(NumberOfSystems,sizeof(int));
    Components[i].BlockDistance=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
    Components[i].BlockCenters=(VECTOR**)calloc(NumberOfSystems,sizeof(VECTOR*));

    if(Components[i].Biased)
    {
      n=Components[i].BiasingFunction.n;
      Components[i].BiasingFunction.a=(REAL*)calloc(n,sizeof(REAL));
      Components[i].BiasingFunction.b=(REAL*)calloc(n,sizeof(REAL));
      Components[i].BiasingFunction.c=(REAL*)calloc(n,sizeof(REAL));
      Components[i].BiasingFunction.d=(REAL*)calloc(n,sizeof(REAL));
      Components[i].BiasingFunction.x=(REAL*)calloc(n,sizeof(REAL));
      Components[i].BiasingFunction.y=(REAL*)calloc(n,sizeof(REAL));
      Components[i].FreeEnergyXdata=(REAL*)calloc(n+1,sizeof(REAL));
    }

    Components[i].Fixed=(int*)calloc(Components[i].NumberOfAtoms,sizeof(int));
    Components[i].Type=(int*)calloc(Components[i].NumberOfAtoms,sizeof(int));
    Components[i].Charge=(REAL*)calloc(Components[i].NumberOfAtoms,sizeof(REAL));
    Components[i].Connectivity=(int*)calloc(Components[i].NumberOfAtoms,sizeof(int));
    Components[i].ConnectivityMatrix=(int**)calloc(Components[i].NumberOfAtoms,sizeof(int*));
    for(j=0;j<Components[i].NumberOfAtoms;j++)
      Components[i].ConnectivityMatrix[j]=(int*)calloc(Components[i].NumberOfAtoms,sizeof(int));
    Components[i].ConnectivityList=(int**)calloc(Components[i].NumberOfAtoms,sizeof(int*));
    Components[i].Positions=(VECTOR*)calloc(Components[i].NumberOfAtoms,sizeof(VECTOR));
    if(Components[i].NumberOfGroups>0)
      Components[i].Groups=(GROUP_DEFINITION*)calloc(Components[i].NumberOfGroups,sizeof(GROUP_DEFINITION));
    Components[i].group=(int*)calloc(Components[i].NumberOfAtoms,sizeof(int));
    Components[i].RMCMOL=(VECTOR*)calloc(Components[i].NumberOfAtoms,sizeof(VECTOR));

    fread(Components[i].IdealGasRosenbluthWeight,sizeof(REAL),NumberOfSystems,FilePtr);
    fread(Components[i].IdealGasTotalEnergy,sizeof(REAL),NumberOfSystems,FilePtr);
    fread(Components[i].PartialPressure,sizeof(REAL),NumberOfSystems,FilePtr);
    fread(Components[i].FugacityCoefficient,sizeof(REAL),NumberOfSystems,FilePtr);
    fread(Components[i].BulkFluidDensity,sizeof(REAL),NumberOfSystems,FilePtr);
    fread(Components[i].Compressibility,sizeof(REAL),NumberOfSystems,FilePtr);
    fread(Components[i].MolFraction,sizeof(REAL),NumberOfSystems,FilePtr);
    fread(Components[i].AmountOfExcessMolecules,sizeof(REAL),NumberOfSystems,FilePtr);

    fread(Components[i].CreateNumberOfMolecules,sizeof(int),NumberOfSystems,FilePtr);
    fread(Components[i].NumberOfMolecules,sizeof(int),NumberOfSystems,FilePtr);

    fread(Components[i].MOLEC_PER_UC_TO_MOL_PER_KG,sizeof(REAL),NumberOfSystems,FilePtr);
    fread(Components[i].MOLEC_PER_UC_TO_GRAM_PER_GRAM_OF_FRAMEWORK,sizeof(REAL),NumberOfSystems,FilePtr);
    fread(Components[i].MOLEC_PER_UC_TO_CC_STP_G,sizeof(REAL),NumberOfSystems,FilePtr);
    fread(Components[i].MOLEC_PER_UC_TO_CC_STP_CC,sizeof(REAL),NumberOfSystems,FilePtr);
    fread(Components[i].MOL_PER_KG_TO_CC_STP_G,sizeof(REAL),NumberOfSystems,FilePtr);
    fread(Components[i].MOL_PER_KG_TO_CC_STP_CC,sizeof(REAL),NumberOfSystems,FilePtr);

    fread(Components[i].IdentityChanges,sizeof(int),NumberOfComponents,FilePtr);
    fread(Components[i].GibbsIdentityChanges,sizeof(int),NumberOfComponents,FilePtr);

    fread(Components[i].BlockPockets,sizeof(int),NumberOfSystems,FilePtr);
    fread(Components[i].ComputeFreeEnergyProfile,sizeof(int),NumberOfSystems,FilePtr);
    fread(Components[i].BlockPocketsFilename,sizeof(char[256]),NumberOfSystems,FilePtr);
    fread(Components[i].NumberOfBlockCenters,sizeof(int),NumberOfSystems,FilePtr);

    if(Components[i].Biased)
    {
      n=Components[i].BiasingFunction.n;
      fread(Components[i].BiasingFunction.a,sizeof(REAL),n,FilePtr);
      fread(Components[i].BiasingFunction.b,sizeof(REAL),n,FilePtr);
      fread(Components[i].BiasingFunction.c,sizeof(REAL),n,FilePtr);
      fread(Components[i].BiasingFunction.d,sizeof(REAL),n,FilePtr);
      fread(Components[i].BiasingFunction.x,sizeof(REAL),n,FilePtr);
      fread(Components[i].BiasingFunction.y,sizeof(REAL),n,FilePtr);
      fread(Components[i].FreeEnergyXdata,sizeof(REAL),n+1,FilePtr);
    }

    fread(Components[i].Fixed,sizeof(int),Components[i].NumberOfAtoms,FilePtr);
    fread(Components[i].Type,sizeof(int),Components[i].NumberOfAtoms,FilePtr);
    fread(Components[i].Charge,sizeof(REAL),Components[i].NumberOfAtoms,FilePtr);

    fread(Components[i].Connectivity,sizeof(int),Components[i].NumberOfAtoms,FilePtr);
    for(j=0;j<Components[i].NumberOfAtoms;j++)
      fread(Components[i].ConnectivityMatrix[j],sizeof(int),Components[i].NumberOfAtoms,FilePtr);
    for(j=0;j<Components[i].NumberOfAtoms;j++)
      Components[i].ConnectivityList[j]=(int*)calloc(Components[i].Connectivity[j],sizeof(int));
    for(j=0;j<Components[i].NumberOfAtoms;j++)
      fread(Components[i].ConnectivityList[j],sizeof(int),Components[i].Connectivity[j],FilePtr);

    fread(Components[i].Positions,sizeof(VECTOR),Components[i].NumberOfAtoms,FilePtr);
    if(Components[i].NumberOfGroups>0)
    {
      fread(Components[i].Groups,sizeof(GROUP_DEFINITION),Components[i].NumberOfGroups,FilePtr);
      for(j=0;j<Components[i].NumberOfGroups;j++)
      {
        Components[i].Groups[j].Atoms=(int*)calloc(Components[i].Groups[j].NumberOfGroupAtoms,sizeof(int));
        fread(Components[i].Groups[j].Atoms,sizeof(int),Components[i].Groups[j].NumberOfGroupAtoms,FilePtr);
      }
    }
    fread(Components[i].group,sizeof(int),Components[i].NumberOfAtoms,FilePtr);

    fread(Components[i].RMCMOL,sizeof(VECTOR),Components[i].NumberOfAtoms,FilePtr);

    Components[i].MaximumCBMCChangeBondLength=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
    Components[i].MaximumCBMCChangeBendAngle=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
    Components[i].MaximumCBMCRotationOnCone=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
    Components[i].CBMCChangeBondLengthAttempts=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
    Components[i].CBMCChangeBondLengthAccepted=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
    Components[i].CBMCChangeBendAngleAttempts=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
    Components[i].CBMCChangeBendAngleAccepted=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
    Components[i].CBMCRotationOnConeAttempts=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
    Components[i].CBMCRotationOnConeAccepted=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
    Components[i].TotalCBMCChangeBondLengthAttempts=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
    Components[i].TotalCBMCChangeBendAngleAttempts=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
    Components[i].TotalCBMCRotationOnConeAttempts=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
    Components[i].TotalCBMCChangeBondLengthAccepted=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
    Components[i].TotalCBMCChangeBendAngleAccepted=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
    Components[i].TotalCBMCRotationOnConeAccepted=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));

    Components[i].FractionalMolecule=(int*)calloc(NumberOfSystems,sizeof(int));
    Components[i].CFMoleculePresent=(int*)calloc(NumberOfSystems,sizeof(int));
    Components[i].RXMCMoleculesPresent=(int*)calloc(NumberOfSystems,sizeof(int));
    Components[i].NumberOfRXMCMoleculesPresent=(int*)calloc(NumberOfSystems,sizeof(int));

    Components[i].CFWangLandauScalingFactor=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
    Components[i].CFBiasingFactors=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));

    fread(Components[i].FractionalMolecule,sizeof(int),NumberOfSystems,FilePtr);
    fread(Components[i].CFMoleculePresent,sizeof(int),NumberOfSystems,FilePtr);
    fread(Components[i].RXMCMoleculesPresent,sizeof(int),NumberOfSystems,FilePtr);
    fread(Components[i].NumberOfRXMCMoleculesPresent,sizeof(int),NumberOfSystems,FilePtr);
    fread(Components[i].CFWangLandauScalingFactor,sizeof(REAL),NumberOfSystems,FilePtr);

    Components[i].ReactantFractionalMolecules=(int***)calloc(NumberOfSystems,sizeof(int**));
    Components[i].ProductFractionalMolecules=(int***)calloc(NumberOfSystems,sizeof(int**));

    for(j=0;j<NumberOfSystems;j++)
    {
      Components[i].BlockDistance[j]=(REAL*)calloc(Components[i].NumberOfBlockCenters[j],sizeof(REAL));
      Components[i].BlockCenters[j]=(VECTOR*)calloc(Components[i].NumberOfBlockCenters[j],sizeof(VECTOR));
      fread(Components[i].BlockDistance[j],sizeof(REAL),Components[i].NumberOfBlockCenters[j],FilePtr);
      fread(Components[i].BlockCenters[j],sizeof(VECTOR),Components[i].NumberOfBlockCenters[j],FilePtr);

      Components[i].MaximumCBMCChangeBondLength[j]=(REAL*)calloc(Components[i].NumberOfAtoms,sizeof(REAL));
      Components[i].MaximumCBMCChangeBendAngle[j]=(REAL*)calloc(Components[i].NumberOfAtoms,sizeof(REAL));
      Components[i].MaximumCBMCRotationOnCone[j]=(REAL*)calloc(Components[i].NumberOfAtoms,sizeof(REAL));
      Components[i].CBMCChangeBondLengthAttempts[j]=(REAL*)calloc(Components[i].NumberOfAtoms,sizeof(REAL));
      Components[i].CBMCChangeBondLengthAccepted[j]=(REAL*)calloc(Components[i].NumberOfAtoms,sizeof(REAL));
      Components[i].CBMCChangeBendAngleAttempts[j]=(REAL*)calloc(Components[i].NumberOfAtoms,sizeof(REAL));
      Components[i].CBMCChangeBendAngleAccepted[j]=(REAL*)calloc(Components[i].NumberOfAtoms,sizeof(REAL));
      Components[i].CBMCRotationOnConeAttempts[j]=(REAL*)calloc(Components[i].NumberOfAtoms,sizeof(REAL));
      Components[i].CBMCRotationOnConeAccepted[j]=(REAL*)calloc(Components[i].NumberOfAtoms,sizeof(REAL));
      Components[i].TotalCBMCChangeBondLengthAttempts[j]=(REAL*)calloc(Components[i].NumberOfAtoms,sizeof(REAL));
      Components[i].TotalCBMCChangeBendAngleAttempts[j]=(REAL*)calloc(Components[i].NumberOfAtoms,sizeof(REAL));
      Components[i].TotalCBMCRotationOnConeAttempts[j]=(REAL*)calloc(Components[i].NumberOfAtoms,sizeof(REAL));
      Components[i].TotalCBMCChangeBondLengthAccepted[j]=(REAL*)calloc(Components[i].NumberOfAtoms,sizeof(REAL));
      Components[i].TotalCBMCChangeBendAngleAccepted[j]=(REAL*)calloc(Components[i].NumberOfAtoms,sizeof(REAL));
      Components[i].TotalCBMCRotationOnConeAccepted[j]=(REAL*)calloc(Components[i].NumberOfAtoms,sizeof(REAL));

      fread(Components[i].MaximumCBMCChangeBondLength[j],sizeof(REAL),Components[i].NumberOfAtoms,FilePtr);
      fread(Components[i].MaximumCBMCChangeBendAngle[j],sizeof(REAL),Components[i].NumberOfAtoms,FilePtr);
      fread(Components[i].MaximumCBMCRotationOnCone[j],sizeof(REAL),Components[i].NumberOfAtoms,FilePtr);
      fread(Components[i].CBMCChangeBondLengthAttempts[j],sizeof(REAL),Components[i].NumberOfAtoms,FilePtr);
      fread(Components[i].CBMCChangeBondLengthAccepted[j],sizeof(REAL),Components[i].NumberOfAtoms,FilePtr);
      fread(Components[i].CBMCChangeBendAngleAttempts[j],sizeof(REAL),Components[i].NumberOfAtoms,FilePtr);
      fread(Components[i].CBMCChangeBendAngleAccepted[j],sizeof(REAL),Components[i].NumberOfAtoms,FilePtr);
      fread(Components[i].CBMCRotationOnConeAttempts[j],sizeof(REAL),Components[i].NumberOfAtoms,FilePtr);
      fread(Components[i].CBMCRotationOnConeAccepted[j],sizeof(REAL),Components[i].NumberOfAtoms,FilePtr);
      fread(Components[i].TotalCBMCChangeBondLengthAttempts[j],sizeof(REAL),Components[i].NumberOfAtoms,FilePtr);
      fread(Components[i].TotalCBMCChangeBendAngleAttempts[j],sizeof(REAL),Components[i].NumberOfAtoms,FilePtr);
      fread(Components[i].TotalCBMCRotationOnConeAttempts[j],sizeof(REAL),Components[i].NumberOfAtoms,FilePtr);
      fread(Components[i].TotalCBMCChangeBondLengthAccepted[j],sizeof(REAL),Components[i].NumberOfAtoms,FilePtr);
      fread(Components[i].TotalCBMCChangeBendAngleAccepted[j],sizeof(REAL),Components[i].NumberOfAtoms,FilePtr);
      fread(Components[i].TotalCBMCRotationOnConeAccepted[j],sizeof(REAL),Components[i].NumberOfAtoms,FilePtr);

      Components[i].CFBiasingFactors[j]=(REAL*)calloc(Components[i].CFLambdaHistogramSize,sizeof(REAL));
      fread(Components[i].CFBiasingFactors[j],sizeof(REAL),Components[i].CFLambdaHistogramSize,FilePtr);

      if(NumberOfReactions>0)
      {
        Components[i].ReactantFractionalMolecules[j]=(int**)calloc(NumberOfReactions,sizeof(int*));
        Components[i].ProductFractionalMolecules[j]=(int**)calloc(NumberOfReactions,sizeof(int*));

        for(k=0;k<NumberOfReactions;k++)
        {
          if(ReactantsStoichiometry[k][i]>0)
          {
            Components[i].ReactantFractionalMolecules[j][k]=(int*)calloc(ReactantsStoichiometry[j][i],sizeof(int));
            fread(Components[i].ReactantFractionalMolecules[j][k],sizeof(int),ReactantsStoichiometry[k][i],FilePtr);
          }
          if(ProductsStoichiometry[k][i]>0)
          {
            Components[i].ProductFractionalMolecules[j][k]=(int*)calloc(ProductsStoichiometry[j][i],sizeof(int));
            fread(Components[i].ProductFractionalMolecules[j][k],sizeof(int),ProductsStoichiometry[k][i],FilePtr);
          }
        }
      }
    }

    // allocate charility-centers
    fread(&Components[i].NumberOfChiralityCenters,sizeof(int),1,FilePtr);
    Components[i].Chirality=(int*)calloc(Components[i].NumberOfAtoms,sizeof(int));
    fread(Components[i].Chirality,sizeof(int),Components[i].NumberOfAtoms,FilePtr);
    if(Components[i].NumberOfChiralityCenters>0)
    {
      Components[i].ChiralityType=(int*)calloc(Components[i].NumberOfAtoms,sizeof(int));
      Components[i].ChiralA=(int*)calloc(Components[i].NumberOfAtoms,sizeof(int));
      Components[i].ChiralB=(int*)calloc(Components[i].NumberOfAtoms,sizeof(int));
      Components[i].ChiralC=(int*)calloc(Components[i].NumberOfAtoms,sizeof(int));
      Components[i].ChiralD=(int*)calloc(Components[i].NumberOfAtoms,sizeof(int));

      fread(Components[i].ChiralityType,sizeof(int),Components[i].NumberOfAtoms,FilePtr);
      fread(Components[i].ChiralA,sizeof(int),Components[i].NumberOfAtoms,FilePtr);
      fread(Components[i].ChiralB,sizeof(int),Components[i].NumberOfAtoms,FilePtr);
      fread(Components[i].ChiralC,sizeof(int),Components[i].NumberOfAtoms,FilePtr);
      fread(Components[i].ChiralD,sizeof(int),Components[i].NumberOfAtoms,FilePtr);
    }

    // allocate bonds
    Components[i].Bonds=(PAIR*)calloc(Components[i].NumberOfBonds,sizeof(PAIR));
    Components[i].BondType=(int*)calloc(Components[i].NumberOfBonds,sizeof(int));
    Components[i].BondArguments=(REAL(*)[MAX_BOND_POTENTIAL_ARGUMENTS])calloc(Components[i].NumberOfBonds,
                                  sizeof(REAL[MAX_BOND_POTENTIAL_ARGUMENTS]));

    // read bonds
    fread(Components[i].Bonds,sizeof(PAIR),Components[i].NumberOfBonds,FilePtr);
    fread(Components[i].BondType,sizeof(int),Components[i].NumberOfBonds,FilePtr);
    fread(Components[i].BondArguments,sizeof(REAL[MAX_BOND_POTENTIAL_ARGUMENTS]),Components[i].NumberOfBonds,FilePtr);

    // allocate bond-dipoles
    Components[i].BondDipoles=(PAIR*)calloc(Components[i].NumberOfBondDipoles,sizeof(PAIR));
    Components[i].BondDipoleMagnitude=(REAL*)calloc(Components[i].NumberOfBondDipoles,sizeof(REAL));

    // read bond-dipoles
    fread(Components[i].BondDipoles,sizeof(PAIR),Components[i].NumberOfBondDipoles,FilePtr);
    fread(Components[i].BondDipoleMagnitude,sizeof(REAL),Components[i].NumberOfBondDipoles,FilePtr);

    // allocate bends
    Components[i].Bends=(QUAD*)calloc(Components[i].NumberOfBends,sizeof(QUAD));
    Components[i].BendType=(int*)calloc(Components[i].NumberOfBends,sizeof(int));
    Components[i].BendArguments=(REAL(*)[MAX_BEND_POTENTIAL_ARGUMENTS])calloc(Components[i].NumberOfBends,
                                    sizeof(REAL[MAX_BEND_POTENTIAL_ARGUMENTS]));

    // read bends
    fread(Components[i].Bends,sizeof(QUAD),Components[i].NumberOfBends,FilePtr);
    fread(Components[i].BendType,sizeof(int),Components[i].NumberOfBends,FilePtr);
    fread(Components[i].BendArguments,sizeof(REAL[MAX_BEND_POTENTIAL_ARGUMENTS]),Components[i].NumberOfBends,FilePtr);

    // allocate urey-bradleys
    Components[i].UreyBradleys=(TRIPLE*)calloc(Components[i].NumberOfUreyBradleys,sizeof(TRIPLE));
    Components[i].UreyBradleyType=(int*)calloc(Components[i].NumberOfUreyBradleys,sizeof(int));
    Components[i].UreyBradleyArguments=(REAL(*)[MAX_UREYBRADLEY_POTENTIAL_ARGUMENTS])calloc(Components[i].NumberOfUreyBradleys,
                                    sizeof(REAL[MAX_UREYBRADLEY_POTENTIAL_ARGUMENTS]));

    // read urey-bradleys
    fread(Components[i].UreyBradleys,sizeof(TRIPLE),Components[i].NumberOfUreyBradleys,FilePtr);
    fread(Components[i].UreyBradleyType,sizeof(int),Components[i].NumberOfUreyBradleys,FilePtr);
    fread(Components[i].UreyBradleyArguments,sizeof(REAL[MAX_UREYBRADLEY_POTENTIAL_ARGUMENTS]),Components[i].NumberOfUreyBradleys,FilePtr);

    // allocate inversion-bends
    Components[i].InversionBends=(QUAD*)calloc(Components[i].NumberOfInversionBends,sizeof(QUAD));
    Components[i].InversionBendType=(int*)calloc(Components[i].NumberOfInversionBends,sizeof(int));
    Components[i].InversionBendArguments=(REAL(*)[MAX_INVERSION_BEND_POTENTIAL_ARGUMENTS])calloc(Components[i].NumberOfInversionBends,
                                    sizeof(REAL[MAX_INVERSION_BEND_POTENTIAL_ARGUMENTS]));

    // read inversion-bends
    fread(Components[i].InversionBends,sizeof(QUAD),Components[i].NumberOfInversionBends,FilePtr);
    fread(Components[i].InversionBendType,sizeof(int),Components[i].NumberOfInversionBends,FilePtr);
    fread(Components[i].InversionBendArguments,sizeof(REAL[MAX_INVERSION_BEND_POTENTIAL_ARGUMENTS]),Components[i].NumberOfInversionBends,FilePtr);

    // allocate torsions
    Components[i].Torsions=(QUAD*)calloc(Components[i].NumberOfTorsions,sizeof(QUAD));
    Components[i].TorsionType=(int*)calloc(Components[i].NumberOfTorsions,sizeof(int));
    Components[i].TorsionArguments=(REAL(*)[MAX_TORSION_POTENTIAL_ARGUMENTS])calloc(Components[i].NumberOfTorsions,
                                    sizeof(REAL[MAX_TORSION_POTENTIAL_ARGUMENTS]));

    // read torsions
    fread(Components[i].Torsions,sizeof(QUAD),Components[i].NumberOfTorsions,FilePtr);
    fread(Components[i].TorsionType,sizeof(int),Components[i].NumberOfTorsions,FilePtr);
    fread(Components[i].TorsionArguments,sizeof(REAL[MAX_TORSION_POTENTIAL_ARGUMENTS]),Components[i].NumberOfTorsions,FilePtr);

    // allocate improper torsions
    Components[i].ImproperTorsions=(QUAD*)calloc(Components[i].NumberOfImproperTorsions,sizeof(QUAD));
    Components[i].ImproperTorsionType=(int*)calloc(Components[i].NumberOfImproperTorsions,sizeof(int));
    Components[i].ImproperTorsionArguments=(REAL(*)[MAX_IMPROPER_TORSION_POTENTIAL_ARGUMENTS])calloc(Components[i].NumberOfImproperTorsions,
                                    sizeof(REAL[MAX_IMPROPER_TORSION_POTENTIAL_ARGUMENTS]));

    // read improper torsions
    fread(Components[i].ImproperTorsions,sizeof(QUAD),Components[i].NumberOfImproperTorsions,FilePtr);
    fread(Components[i].ImproperTorsionType,sizeof(int),Components[i].NumberOfImproperTorsions,FilePtr);
    fread(Components[i].ImproperTorsionArguments,sizeof(REAL[MAX_IMPROPER_TORSION_POTENTIAL_ARGUMENTS]),Components[i].NumberOfImproperTorsions,FilePtr);

    // allocate out-of-plane
    Components[i].OutOfPlanes=(QUAD*)calloc(Components[i].NumberOfOutOfPlanes,sizeof(QUAD));
    Components[i].OutOfPlaneType=(int*)calloc(Components[i].NumberOfOutOfPlanes,sizeof(int));
    Components[i].OutOfPlaneArguments=(REAL(*)[MAX_OUT_OF_PLANE_POTENTIAL_ARGUMENTS])calloc(Components[i].NumberOfOutOfPlanes,
                                    sizeof(REAL[MAX_OUT_OF_PLANE_POTENTIAL_ARGUMENTS]));

    // read out-of-plane
    fread(Components[i].OutOfPlanes,sizeof(QUAD),Components[i].NumberOfOutOfPlanes,FilePtr);
    fread(Components[i].OutOfPlaneType,sizeof(int),Components[i].NumberOfOutOfPlanes,FilePtr);
    fread(Components[i].OutOfPlaneArguments,sizeof(REAL[MAX_OUT_OF_PLANE_POTENTIAL_ARGUMENTS]),Components[i].NumberOfOutOfPlanes,FilePtr);

    // allocate bond-bonds
    Components[i].BondBonds=(TRIPLE*)calloc(Components[i].NumberOfBondBonds,sizeof(TRIPLE));
    Components[i].BondBondType=(int*)calloc(Components[i].NumberOfBondBonds,sizeof(int));
    Components[i].BondBondArguments=(REAL(*)[MAX_BOND_BOND_POTENTIAL_ARGUMENTS])calloc(Components[i].NumberOfBondBonds,
                                    sizeof(REAL[MAX_BOND_BOND_POTENTIAL_ARGUMENTS]));

    // read bond-bonds
    fread(Components[i].BondBonds,sizeof(TRIPLE),Components[i].NumberOfBondBonds,FilePtr);
    fread(Components[i].BondBondType,sizeof(int),Components[i].NumberOfBondBonds,FilePtr);
    fread(Components[i].BondBondArguments,sizeof(REAL[MAX_BOND_BOND_POTENTIAL_ARGUMENTS]),Components[i].NumberOfBondBonds,FilePtr);

    // allocate bond-bends
    Components[i].BondBends=(TRIPLE*)calloc(Components[i].NumberOfBondBends,sizeof(TRIPLE));
    Components[i].BondBendType=(int*)calloc(Components[i].NumberOfBondBends,sizeof(int));
    Components[i].BondBendArguments=(REAL(*)[MAX_BOND_BEND_POTENTIAL_ARGUMENTS])calloc(Components[i].NumberOfBondBends,
                                    sizeof(REAL[MAX_BOND_BEND_POTENTIAL_ARGUMENTS]));

    // read bond-bends
    fread(Components[i].BondBends,sizeof(TRIPLE),Components[i].NumberOfBondBends,FilePtr);
    fread(Components[i].BondBendType,sizeof(int),Components[i].NumberOfBondBends,FilePtr);
    fread(Components[i].BondBendArguments,sizeof(REAL[MAX_BOND_BEND_POTENTIAL_ARGUMENTS]),Components[i].NumberOfBondBends,FilePtr);

    // allocate bend-bends
    Components[i].BendBends=(QUAD*)calloc(Components[i].NumberOfBendBends,sizeof(QUAD));
    Components[i].BendBendType=(int*)calloc(Components[i].NumberOfBendBends,sizeof(int));
    Components[i].BendBendArguments=(REAL(*)[MAX_BEND_BEND_POTENTIAL_ARGUMENTS])calloc(Components[i].NumberOfBendBends,
                                    sizeof(REAL[MAX_BEND_BEND_POTENTIAL_ARGUMENTS]));

    // read bend-bends
    fread(Components[i].BendBends,sizeof(QUAD),Components[i].NumberOfBendBends,FilePtr);
    fread(Components[i].BendBendType,sizeof(int),Components[i].NumberOfBendBends,FilePtr);
    fread(Components[i].BendBendArguments,sizeof(REAL[MAX_BEND_BEND_POTENTIAL_ARGUMENTS]),Components[i].NumberOfBendBends,FilePtr);

    // allocate stretch-torsions
    Components[i].BondTorsions=(QUAD*)calloc(Components[i].NumberOfBondTorsions,sizeof(QUAD));
    Components[i].BondTorsionType=(int*)calloc(Components[i].NumberOfBondTorsions,sizeof(int));
    Components[i].BondTorsionArguments=(REAL(*)[MAX_BOND_TORSION_POTENTIAL_ARGUMENTS])calloc(Components[i].NumberOfBondTorsions,
                                    sizeof(REAL[MAX_BOND_TORSION_POTENTIAL_ARGUMENTS]));

    // read stretch-torsions
    fread(Components[i].BondTorsions,sizeof(QUAD),Components[i].NumberOfBondTorsions,FilePtr);
    fread(Components[i].BondTorsionType,sizeof(int),Components[i].NumberOfBondTorsions,FilePtr);
    fread(Components[i].BondTorsionArguments,sizeof(REAL[MAX_BOND_TORSION_POTENTIAL_ARGUMENTS]),Components[i].NumberOfBondTorsions,FilePtr);

    // allocate bend-torsions
    Components[i].BendTorsions=(QUAD*)calloc(Components[i].NumberOfBendTorsions,sizeof(QUAD));
    Components[i].BendTorsionType=(int*)calloc(Components[i].NumberOfBendTorsions,sizeof(int));
    Components[i].BendTorsionArguments=(REAL(*)[MAX_BEND_TORSION_POTENTIAL_ARGUMENTS])calloc(Components[i].NumberOfBendTorsions,
                                    sizeof(REAL[MAX_BEND_TORSION_POTENTIAL_ARGUMENTS]));

    // read bend-torsions
    fread(Components[i].BendTorsions,sizeof(QUAD),Components[i].NumberOfBendTorsions,FilePtr);
    fread(Components[i].BendTorsionType,sizeof(int),Components[i].NumberOfBendTorsions,FilePtr);
    fread(Components[i].BendTorsionArguments,sizeof(REAL[MAX_BEND_TORSION_POTENTIAL_ARGUMENTS]),Components[i].NumberOfBendTorsions,FilePtr);

    Components[i].IntraVDW=(PAIR*)calloc(Components[i].NumberOfIntraVDW,sizeof(PAIR));
    Components[i].IntraChargeCharge=(PAIR*)calloc(Components[i].NumberOfIntraChargeCharge,sizeof(PAIR));
    Components[i].IntraChargeBondDipole=(PAIR*)calloc(Components[i].NumberOfIntraChargeBondDipole,sizeof(PAIR));
    Components[i].IntraBondDipoleBondDipole=(PAIR*)calloc(Components[i].NumberOfIntraBondDipoleBondDipole,sizeof(PAIR));

    Components[i].IntraVDWScaling=(REAL*)calloc(Components[i].NumberOfIntraVDW,sizeof(REAL));
    Components[i].IntraChargeChargeScaling=(REAL*)calloc(Components[i].NumberOfIntraChargeCharge,sizeof(REAL));

    // read intra-pairs
    fread(Components[i].IntraVDW,sizeof(PAIR),Components[i].NumberOfIntraVDW,FilePtr);
    fread(Components[i].IntraChargeCharge,sizeof(PAIR),Components[i].NumberOfIntraChargeCharge,FilePtr);
    fread(Components[i].IntraChargeBondDipole,sizeof(PAIR),Components[i].NumberOfIntraChargeBondDipole,FilePtr);
    fread(Components[i].IntraBondDipoleBondDipole,sizeof(PAIR),Components[i].NumberOfIntraBondDipoleBondDipole,FilePtr);

    fread(Components[i].IntraVDWScaling,sizeof(REAL),Components[i].NumberOfIntraVDW,FilePtr);
    fread(Components[i].IntraChargeChargeScaling,sizeof(REAL),Components[i].NumberOfIntraChargeCharge,FilePtr);

    // allocate excluded pairs
    Components[i].ExcludedIntraChargeCharge=(PAIR*)calloc(Components[i].NumberOfExcludedIntraChargeCharge,sizeof(PAIR));
    Components[i].ExcludedIntraChargeBondDipole=(PAIR*)calloc(Components[i].NumberOfExcludedIntraChargeBondDipole,sizeof(PAIR));
    Components[i].ExcludedIntraBondDipoleBondDipole=(PAIR*)calloc(Components[i].NumberOfExcludedIntraBondDipoleBondDipole,sizeof(PAIR));

    // read excluded-pairs
    fread(Components[i].ExcludedIntraChargeCharge,sizeof(PAIR),Components[i].NumberOfExcludedIntraChargeCharge,FilePtr);
    fread(Components[i].ExcludedIntraChargeBondDipole,sizeof(PAIR),Components[i].NumberOfExcludedIntraChargeBondDipole,FilePtr);
    fread(Components[i].ExcludedIntraBondDipoleBondDipole,sizeof(PAIR),Components[i].NumberOfExcludedIntraBondDipoleBondDipole,FilePtr);

    // allocate and read config moves
    fread(&Components[i].NumberOfConfigMoves,sizeof(int),1,FilePtr);
    Components[i].NumberOfUnchangedAtomsConfig=(int*)calloc(Components[i].NumberOfConfigMoves,sizeof(int));
    Components[i].UnchangedAtomsConfig=(int**)calloc(Components[i].NumberOfConfigMoves,sizeof(int*));
    fread(Components[i].NumberOfUnchangedAtomsConfig,sizeof(int),Components[i].NumberOfConfigMoves,FilePtr);
    for(j=0;j<Components[i].NumberOfConfigMoves;j++)
    {
      Components[i].UnchangedAtomsConfig[j]=(int*)calloc(Components[i].NumberOfUnchangedAtomsConfig[j],sizeof(int));
      fread(Components[i].UnchangedAtomsConfig[j],sizeof(int),Components[i].NumberOfUnchangedAtomsConfig[j],FilePtr);
    }

    // allocate and read identity moves
    fread(&Components[i].NumberOfIdentityConfigMoves,sizeof(int),1,FilePtr);
    Components[i].NumberOfUnchangedAtomsIdentityConfig=(int*)calloc(Components[i].NumberOfIdentityConfigMoves,sizeof(int));
    Components[i].UnchangedAtomsIdentityConfig=(int**)calloc(Components[i].NumberOfIdentityConfigMoves,sizeof(int*));
    fread(Components[i].NumberOfUnchangedAtomsIdentityConfig,sizeof(int),Components[i].NumberOfIdentityConfigMoves,FilePtr);
    for(j=0;j<Components[i].NumberOfIdentityConfigMoves;j++)
    {
      Components[i].UnchangedAtomsIdentityConfig[j]=(int*)calloc(Components[i].NumberOfUnchangedAtomsIdentityConfig[j],sizeof(int));
      fread(Components[i].UnchangedAtomsIdentityConfig[j],sizeof(int),Components[i].NumberOfUnchangedAtomsIdentityConfig[j],FilePtr);
    }

    fread(&Components[i].RestrictMoves,sizeof(int),1,FilePtr);
    fread(&Components[i].RestrictMovesToBox,sizeof(int),1,FilePtr);
    fread(&Components[i].BoxAxisABC_Min,sizeof(VECTOR),1,FilePtr);
    fread(&Components[i].BoxAxisABC_Max,sizeof(VECTOR),1,FilePtr);
    fread(&Components[i].BoxAxisABC_Min2,sizeof(VECTOR),1,FilePtr);
    fread(&Components[i].BoxAxisABC_Max2,sizeof(VECTOR),1,FilePtr);
    fread(&Components[i].BoxAxisABC_Min3,sizeof(VECTOR),1,FilePtr);
    fread(&Components[i].BoxAxisABC_Max3,sizeof(VECTOR),1,FilePtr);
    fread(&Components[i].BoxAxisABC_Min4,sizeof(VECTOR),1,FilePtr);
    fread(&Components[i].BoxAxisABC_Max4,sizeof(VECTOR),1,FilePtr);
    fread(&Components[i].RestrictMovesToPrisms,sizeof(int),1,FilePtr);
    fread(Components[i].RestrictMovesToPrism,sizeof(int),MAX_NUMBER_OF_PRISMS,FilePtr);
    fread(Components[i].RestrictPrismABC_Min,sizeof(VECTOR),MAX_NUMBER_OF_PRISMS,FilePtr);
    fread(Components[i].RestrictPrismABC_Max,sizeof(VECTOR),MAX_NUMBER_OF_PRISMS,FilePtr);
    fread(&Components[i].RestrictMovesToCylinders,sizeof(int),1,FilePtr);
    fread(Components[i].RestrictCylinderABC_Min,sizeof(VECTOR),MAX_NUMBER_OF_CYLINDERS,FilePtr);
    fread(Components[i].RestrictCylinderABC_Max,sizeof(VECTOR),MAX_NUMBER_OF_CYLINDERS,FilePtr);
    fread(Components[i].RestrictCylinderCenter,sizeof(VECTOR),MAX_NUMBER_OF_CYLINDERS,FilePtr);
    fread(Components[i].RestrictCylinderDirection,sizeof(int),MAX_NUMBER_OF_CYLINDERS,FilePtr);
    fread(Components[i].RestrictCylinderRadius,sizeof(REAL),MAX_NUMBER_OF_CYLINDERS,FilePtr);
    fread(&Components[i].RestrictMovesToSpheres,sizeof(int),1,FilePtr);
    fread(Components[i].RestrictMovesToSphere,sizeof(int),MAX_NUMBER_OF_SPHERES,FilePtr);
    fread(Components[i].RestrictSphereCenter,sizeof(VECTOR),MAX_NUMBER_OF_SPHERES,FilePtr);
    fread(Components[i].RestrictSphereRadius,sizeof(REAL),MAX_NUMBER_OF_SPHERES,FilePtr);

    Components[i].CpuTimeTranslationMove=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
    Components[i].CpuTimeRandomTranslationMove=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
    Components[i].CpuTimeRotationMove=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
    Components[i].CpuTimeRandomRotationMove=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
    Components[i].CpuTimePartialReinsertionMove=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
    Components[i].CpuTimeReinsertionMove=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
    Components[i].CpuTimeReinsertionInPlaceMove=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
    Components[i].CpuTimeReinsertionInPlaneMove=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
    Components[i].CpuTimeIdentityChangeMove=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
    Components[i].CpuTimeSwapMoveInsertion=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
    Components[i].CpuTimeSwapMoveDeletion=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
    Components[i].CpuTimeCFSwapLambdaMove=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
    Components[i].CpuTimeCBCFSwapLambdaMove=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
    Components[i].CpuTimeWidomMove=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
    Components[i].CpuTimeCFWidomLambdaMove=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
    Components[i].CpuTimeGibbsWidomMove=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
    Components[i].CpuTimeSurfaceAreaMove=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
    Components[i].CpuTimeGibbsChangeMove=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
    Components[i].CpuTimeCFGibbsChangeMove=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
    Components[i].CpuTimeCBCFGibbsChangeMove=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
    Components[i].CpuTimeGibbsIdentityChangeMove=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
    Components[i].CpuTimeExchangeFractionalParticleMove=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
    Components[i].CpuTimeCFGibbsSwapFractionalMoleculeToOtherBoxMove=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
    Components[i].CpuTimeCFGibbsLambdaChangeMove=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
    Components[i].CpuTimeCFGibbsFractionalToIntegerMove=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

    fread(Components[i].CpuTimeTranslationMove,sizeof(REAL),NumberOfSystems,FilePtr);
    fread(Components[i].CpuTimeRandomTranslationMove,sizeof(REAL),NumberOfSystems,FilePtr);
    fread(Components[i].CpuTimeRotationMove,sizeof(REAL),NumberOfSystems,FilePtr);
    fread(Components[i].CpuTimeRandomRotationMove,sizeof(REAL),NumberOfSystems,FilePtr);
    fread(Components[i].CpuTimePartialReinsertionMove,sizeof(REAL),NumberOfSystems,FilePtr);
    fread(Components[i].CpuTimeReinsertionMove,sizeof(REAL),NumberOfSystems,FilePtr);
    fread(Components[i].CpuTimeReinsertionInPlaceMove,sizeof(REAL),NumberOfSystems,FilePtr);
    fread(Components[i].CpuTimeReinsertionInPlaneMove,sizeof(REAL),NumberOfSystems,FilePtr);
    fread(Components[i].CpuTimeIdentityChangeMove,sizeof(REAL),NumberOfSystems,FilePtr);
    fread(Components[i].CpuTimeSwapMoveInsertion,sizeof(REAL),NumberOfSystems,FilePtr);
    fread(Components[i].CpuTimeSwapMoveDeletion,sizeof(REAL),NumberOfSystems,FilePtr);
    fread(Components[i].CpuTimeCFSwapLambdaMove,sizeof(REAL),NumberOfSystems,FilePtr);
    fread(Components[i].CpuTimeCBCFSwapLambdaMove,sizeof(REAL),NumberOfSystems,FilePtr);
    fread(Components[i].CpuTimeWidomMove,sizeof(REAL),NumberOfSystems,FilePtr);
    fread(Components[i].CpuTimeCFWidomLambdaMove,sizeof(REAL),NumberOfSystems,FilePtr);
    fread(Components[i].CpuTimeGibbsWidomMove,sizeof(REAL),NumberOfSystems,FilePtr);
    fread(Components[i].CpuTimeSurfaceAreaMove,sizeof(REAL),NumberOfSystems,FilePtr);
    fread(Components[i].CpuTimeGibbsChangeMove,sizeof(REAL),NumberOfSystems,FilePtr);
    fread(Components[i].CpuTimeCFGibbsChangeMove,sizeof(REAL),NumberOfSystems,FilePtr);
    fread(Components[i].CpuTimeCBCFGibbsChangeMove,sizeof(REAL),NumberOfSystems,FilePtr);
    fread(Components[i].CpuTimeGibbsIdentityChangeMove,sizeof(REAL),NumberOfSystems,FilePtr);
    fread(Components[i].CpuTimeExchangeFractionalParticleMove,sizeof(REAL),NumberOfSystems,FilePtr);
    fread(Components[i].CpuTimeCFGibbsSwapFractionalMoleculeToOtherBoxMove,sizeof(REAL),NumberOfSystems,FilePtr);
    fread(Components[i].CpuTimeCFGibbsLambdaChangeMove,sizeof(REAL),NumberOfSystems,FilePtr);
    fread(Components[i].CpuTimeCFGibbsFractionalToIntegerMove,sizeof(REAL),NumberOfSystems,FilePtr);
  }


  fread(&Check,1,sizeof(REAL),FilePtr);
  if(fabs(Check-123456789.0)>1e-10)
  {
    fprintf(stderr, "Error in binary restart-file (ReadRestartComponent)\n");
    ContinueAfterCrash=FALSE;
  }
  printf("DONE!!!\n");
}
