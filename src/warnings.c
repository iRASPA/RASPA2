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
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "warnings.h"
#include "framework.h"
#include "molecule.h"
#include "potentials.h"
#include "simulation.h"
#include "output.h"
#include "movies.h"
#include "utils.h"

int *NumberOfWarnings;
int (*Warnings)[MAX_NUMBER_OF_WARNINGS];
int (*NumberOfWarningValues)[MAX_NUMBER_OF_WARNINGS];
char (*WarningValues)[MAX_NUMBER_OF_WARNINGS][MAX_NUMBER_OF_WARNING_ARGUMENTS][32];


void CheckForErrors(void)
{
  int i,j;
  // check reactions

  for(i=0;i<NumberOfReactions;i++)
  {
    for(j=0;j<NumberOfComponents;j++)
    {
      if((ReactantsStoichiometry[i][j]>0)&&((Components[i].ExtraFrameworkMolecule)||(Components[j].ExtraFrameworkMolecule)))
      {
        printf("ERROR: rxmc can only be defined for adsorbates and not for cation-components.\n");
        printf("       change your ExtraFrameworkMolecule for reaction-components from YES to NO.\n");
        exit(0);
      }
      if((ProductsStoichiometry[i][j]>0)&&((Components[i].ExtraFrameworkMolecule)||(Components[j].ExtraFrameworkMolecule)))
      {
        printf("ERROR: rxmc can only be defined for adsorbates and not for cation-components.\n");
        printf("       change your ExtraFrameworkMolecule for reaction-components from YES to NO.\n");
        exit(0);
      }
    }
  }
  
  // Added by Ambroise de Izarra
  //-------------------------------------------------------------------
  // Check the alchemical reactions
  int k;
 
  // Read over the chemical reactions
  for(i=0;i<NumberAlchemicalReactions;i++)
  {
	// Read over the index of salt.
	for(j=0;j<2;j++)
	{
		
		// We check if the input index for the chemical transformation is present in the components.
		if((SaltIndex[i][j]<0)||SaltIndex[i][j]>(NumberOfComponents-1))
		{
		        printf("ERROR: alchemical index %d is not present in the component.\n",SaltIndex[i][j]);
				exit(0);	
		}
		
		// We check if the ions are monoatomic.
		if(Components[SaltIndex[i][j]].NumberOfAtoms>1)
		{
		        printf("ERROR: a salt ion is not monoatomic.\n");
				exit(0);					
		}	
	} 
  }

  //-------------------------------------------------------------------
}


void CreateWarnings(void)
{
  int i,j,k,l,m,f1;
  int Type;
  int StartingBead;
  REAL NetSystemCharge,NetCharge;
  int already_present;

  // set running warning
  for(k=0;k<NumberOfSystems;k++)
  {
    // check for the Lowenstein rule for frameworks
    if(!Lowenstein[k])
    {
      if(NumberOfWarnings[k]<MAX_NUMBER_OF_WARNINGS)
        Warnings[k][NumberOfWarnings[k]++]=LOWENSTEIN_RULE_NOT_OBEYED;
    }

    // check for appropriate use of unit cells
    if((double)CutOffVDW*2 > NumberOfReplicaCells[CurrentSystem].x*(double)BoxProperties[k].cx || 
       (double)CutOffVDW*2 > NumberOfReplicaCells[CurrentSystem].y*(double)BoxProperties[k].cy || 
       (double)CutOffVDW*2 > NumberOfReplicaCells[CurrentSystem].z*(double)BoxProperties[k].cz)
    {
      if(NumberOfWarnings[k]<MAX_NUMBER_OF_WARNINGS)
      {
        Warnings[k][NumberOfWarnings[k]]=UNIT_CELL;
        NumberOfWarnings[k]++;
      }
    }

    // check for net-charge of the system
    NetSystemCharge=0.0;
    for(f1=0;f1<Framework[k].NumberOfFrameworks;f1++)
    {
      for(i=0;i<Framework[k].NumberOfAtoms[f1];i++)
      {
        Type=Framework[k].Atoms[f1][i].Type;
        NetSystemCharge+=Framework[k].Atoms[f1][i].Charge;
      }
    }

    for(i=0;i<NumberOfAdsorbateMolecules[k];i++)
      for(j=0;j<Adsorbates[k][i].NumberOfAtoms;j++)
        NetSystemCharge+=Adsorbates[k][i].Atoms[j].Charge;

    for(i=0;i<NumberOfCationMolecules[k];i++)
      for(j=0;j<Cations[k][i].NumberOfAtoms;j++)
        NetSystemCharge+=Cations[k][i].Atoms[j].Charge;

    if(fabs(NetSystemCharge)>1e-4)
    {
      if(NumberOfWarnings[k]<MAX_NUMBER_OF_WARNINGS)
      {
        Warnings[k][NumberOfWarnings[k]]=NET_SYSTEM_CHARGE;
        sprintf(WarningValues[k][NumberOfWarnings[k]][NumberOfWarningValues[k][NumberOfWarnings[k]]],"%lg",NetSystemCharge);
        NumberOfWarnings[k]++;
      }
    }




    for(i=0;i<NumberOfComponents;i++)
    {
      StartingBead=Components[i].StartingBead;
      Type=Components[i].Type[StartingBead];
      if(PotentialType[Type][Type]==MM3_HYDROGEN_VDW)
      {
        fprintf(stderr, "Fatal error: Can not use starting bead (%d) for component (%s) using the MM3-hydrogen potential.\n",StartingBead,Components[Type].Name);
        fprintf(stderr, "             This potential shifts the interaction site along the carbon-hydrogen bond. When starting the CBMC from hydrogen\n");
        fprintf(stderr, "             this interaction site is not yet known. Growing hydrogen from the carbon present no problem.\n");
        exit(0);
      }

      NetCharge=0.0;
      for(j=0;j<Components[i].NumberOfAtoms;j++)
        NetCharge+=Components[i].Charge[j];

      if((fabs(NetCharge)>1e-4)&&((Components[i].FractionOfReinsertionMove>0.0)||
         (Components[i].FractionOfReinsertionInPlaceMove>0.0)||(Components[i].FractionOfReinsertionInPlaneMove>0.0)))
      {
        if(NumberOfWarnings[k]<MAX_NUMBER_OF_WARNINGS)
        {
          Warnings[k][NumberOfWarnings[k]]=REINSERTION_MOVE_IONS;
          NumberOfWarnings[k]++;
        }
      }
    }

    for(i=0;i<NumberOfPseudoAtoms;i++)
      for(j=0;j<NumberOfPseudoAtoms;j++)
      {
        if(((NumberOfPseudoAtomsType[k][i]>0||PseudoAtoms[i].Swapable)&&(NumberOfPseudoAtomsType[k][j]>0||PseudoAtoms[j].Swapable))
           &&(PotentialType[i][j]==UNDEFINED_POTENTIAL)&&(PseudoAtoms[i].HasVDWInteraction)&&(PseudoAtoms[j].HasVDWInteraction))
        {
          already_present=FALSE;
          for(l=0;l<NumberOfWarnings[k];l++)
          {
            if(Warnings[k][l]==MISSING_INTERACTION)
            {
              already_present=TRUE;
              break;
            }
          }

          if(!already_present)
          {
            if(NumberOfWarnings[k]<MAX_NUMBER_OF_WARNINGS)
            {
              strcpy(WarningValues[k][NumberOfWarnings[k]][0],PseudoAtoms[i].Name);
              strcpy(WarningValues[k][NumberOfWarnings[k]][1],PseudoAtoms[j].Name);
              NumberOfWarningValues[k][NumberOfWarnings[k]]=2;
              Warnings[k][NumberOfWarnings[k]]=MISSING_INTERACTION;
              NumberOfWarnings[k]++;
            }
          }
          else // found at index 'l'
          {
            already_present=FALSE;
            for(m=0;m<NumberOfWarningValues[k][l];m+=2)
              if((strcasecmp(WarningValues[k][l][m],PseudoAtoms[i].Name)==0)&&
                 (strcasecmp(WarningValues[k][l][m],PseudoAtoms[j].Name)==0)) already_present=TRUE;
            if(!already_present)
            {
              if(NumberOfWarningValues[k][l]+2<=MAX_NUMBER_OF_WARNING_ARGUMENTS)
              {
                strcpy(WarningValues[k][l][NumberOfWarningValues[k][l]],PseudoAtoms[i].Name);
                strcpy(WarningValues[k][l][NumberOfWarningValues[k][l]+1],PseudoAtoms[j].Name);
                NumberOfWarningValues[k][l]+=2;
              }
            }
          }
        }
      }
  }
}

void PrintWarningStatus(void)
{
  int i,j;
  FILE *FilePtr;

  FilePtr=OutputFilePtr[CurrentSystem];
  for(i=0;i<NumberOfWarnings[CurrentSystem];i++)
  {
    switch(Warnings[CurrentSystem][i])
    {
      case NET_SYSTEM_CHARGE:
        fprintf(FilePtr,"WARNING: THE SYSTEM HAS A NET CHARGE ");
        for(j=0;j<NumberOfWarningValues[CurrentSystem][i];j++)
          fprintf(FilePtr,"%s ",WarningValues[CurrentSystem][i][j]);
        fprintf(FilePtr,"\n");
        break;
      case LOWENSTEIN_RULE_NOT_OBEYED:
        fprintf(FilePtr,"WARNING: THE LOWENSTEIN RULE IS NOT OBEYED\n");
        break;
      case UNIT_CELL:
        fprintf(FilePtr,"WARNING: INAPPROPRIATE NUMBER OF UNIT CELLS USED\n");
        break;
      case OMITTED_HOST_ADSORBATE_VDW_INTERACTIONS:
        fprintf(FilePtr,"WARNING: OMITTED VDW INTERACTIONS HOST-ADSORBATE\n");
        break;
      case OMITTED_HOST_CATION_VDW_INTERACTIONS:
        fprintf(FilePtr,"WARNING: OMITTED VDW INTERACTIONS HOST-CATION\n");
        break;
      case OMITTED_ADSORBATE_ADSORBATE_VDW_INTERACTIONS:
        fprintf(FilePtr,"WARNING: OMITTED VDW INTERACTIONS ADSORBATE-ADSORBATE\n");
        break;
      case OMITTED_ADSORBATE_CATION_VDW_INTERACTIONS:
        fprintf(FilePtr,"WARNING: OMITTED VDW INTERACTIONS ADSORBATE-CATION\n");
        break;
      case OMITTED_CATION_CATION_VDW_INTERACTIONS:
        fprintf(FilePtr,"WARNING: OMITTED VDW INTERACTIONS CATION-CATION\n");
        break;
      case REINSERTION_MOVE_IONS:
        fprintf(FilePtr,"WARNING: REINSERTION MOVE USED ON CHARGED IONS (IF POSSIBLE, CHANGE TO RANDOM TRANSLATION MOVE TO AVOID NUMERICAL PROBLEMS)\n");
        break;
      case GRID_ERROR_ENERGY:
        fprintf(FilePtr,"WARNING: GRID ENERGY INTERPOLATION PROBABLY NOT ACCURATE ENOUGH\n");
        break;
      case GRID_ERROR_FORCE:
        fprintf(FilePtr,"WARNING: GRID FORCE INTERPOLATION PROBABLY NOT ACCURATE ENOUGH\n");
        break;
      case MISSING_INTERACTION:
        fprintf(FilePtr,"WARNING: THERE ARE ATOM-PAIRS WITH NO VDW INTERACTION ");
        for(j=0;j<NumberOfWarningValues[CurrentSystem][i];j+=2)
          fprintf(FilePtr,"%s-%s ",WarningValues[CurrentSystem][i][j],WarningValues[CurrentSystem][i][j+1]);
        fprintf(FilePtr," (maximum %d interactions shown)\n",MAX_NUMBER_OF_WARNING_ARGUMENTS/2);
        break;
      case ENERGY_DRIFT:
        fprintf(FilePtr,"WARNING: ENERGY DRIFT (INTERNAL CONSISTENCY ERROR IN THE CODE), THE SIMULATION RESULTS ARE WRONG!!\n");
        break;
    }
  }
}

static int versionNumber=1;

void WriteRestartWarnings(FILE *FilePtr)
{
  REAL Check;

  fwrite(&versionNumber,sizeof(int),1,FilePtr);

  fwrite(NumberOfWarnings,sizeof(int),NumberOfSystems,FilePtr);
  fwrite(Warnings,sizeof(int[MAX_NUMBER_OF_WARNINGS]),NumberOfSystems,FilePtr);
  fwrite(NumberOfWarningValues,sizeof(int[MAX_NUMBER_OF_WARNINGS]),NumberOfSystems,FilePtr);
  fwrite(WarningValues,sizeof(char[MAX_NUMBER_OF_WARNINGS][MAX_NUMBER_OF_WARNING_ARGUMENTS][32]),NumberOfSystems,FilePtr);

  Check=123456789.0;
  fwrite(&Check,1,sizeof(REAL),FilePtr);
}

void AllocateWarningsMemory(void)
{
  // allocate memory for warnings
  NumberOfWarnings=(int*)calloc(NumberOfSystems,sizeof(int));
  Warnings=(int(*)[MAX_NUMBER_OF_WARNINGS])calloc(NumberOfSystems,sizeof(int[MAX_NUMBER_OF_WARNINGS]));
  NumberOfWarningValues=(int(*)[MAX_NUMBER_OF_WARNINGS])calloc(NumberOfSystems,sizeof(int[MAX_NUMBER_OF_WARNINGS]));
  WarningValues=(char(*)[MAX_NUMBER_OF_WARNINGS][MAX_NUMBER_OF_WARNING_ARGUMENTS][32])calloc(NumberOfSystems,sizeof(char[MAX_NUMBER_OF_WARNINGS][MAX_NUMBER_OF_WARNING_ARGUMENTS][32]));
}

void ReadRestartWarnings(FILE *FilePtr)
{
  REAL Check;
  int readversionNumber=0;

  AllocateWarningsMemory();

  fread(&readversionNumber,sizeof(int),1,FilePtr);
  if(readversionNumber > versionNumber)
  {
    fprintf(stderr,"Upgrade to last version of RASPA to read binary restart-file");
    exit(-1);
  }

  fread(NumberOfWarnings,sizeof(int),NumberOfSystems,FilePtr);
  fread(Warnings,sizeof(int[MAX_NUMBER_OF_WARNINGS]),NumberOfSystems,FilePtr);
  fread(NumberOfWarningValues,sizeof(int[MAX_NUMBER_OF_WARNINGS]),NumberOfSystems,FilePtr);
  fread(WarningValues,sizeof(char[MAX_NUMBER_OF_WARNINGS][MAX_NUMBER_OF_WARNING_ARGUMENTS][32]),NumberOfSystems,FilePtr);

  fread(&Check,1,sizeof(REAL),FilePtr);
  if(fabs(Check-123456789.0)>1e-10)
  {
    fprintf(stderr, "Error in binary restart-file (ReadRestartWarnings)\n");
    ContinueAfterCrash=FALSE;
  }
}
