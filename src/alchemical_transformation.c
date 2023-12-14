// Alchemical_transformation.c Added by Ambroise de Izarra
/*  
 * Alchemical transformation is added to build functions that peform
 * an alchemical transformation at a fixed number of particule in the system (osmotic ensemble)
 */ 
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <stdbool.h>
#include <time.h>
#include <string.h>
#include <ctype.h>
#include <sys/stat.h>
#include "alchemical_transformation.h"
#include "framework_energy.h"
#include "molecule.h"
#include "framework.h"
#include "simulation.h"
#include "ewald.h"
#include "cbmc.h"
#include "potentials.h"
#include "utils.h"
#include "output.h"
#include "input.h"
#include "cbmc.h"
#include "mc_moves.h"
#include "ewald.h"
#include "grids.h"
#include "spacegroup.h"
#include "recrossing.h"
#include "movies.h"
#include "utils.h"
#include "sample.h"
#include "spectra.h"
#include "rigid.h"
#include "integration.h"
#include "statistics.h"
#include "thermo_baro_stats.h"
#include "equations_of_state.h"
#include "inter_energy.h"
#include "minimization.h"
#include "scattering_factors.h"
#include "framework_force.h"
#include "warnings.h"
#include "internal_force.h"
#include "status.h"
#include "charge_equilibration.h"

static REAL (**VectorforMixingRuleAlchemicalTransformation)[2];
static REAL (**VectorChargeAlchemicalTransformation);

int InitialPseudoAtoms;
int NumberExtraPseudoAtoms;
int NumberExtraComponents;
int CurrentAlchemicalReaction;
int IndexExtraComponent[2];
int TransientMoleculeNbAtoms;
static REAL RosenbluthNew;             // the new Rosenbluth weight (grow)
static REAL RosenbluthOld; 

// Save in memory initial coordinates of chosen moities 
ADSORBATE_MOLECULE *AdsorbatesReferenceChosen;
VECTOR **SolventBodyfixedPositions; // the positions in the body-fixed frame to be stored.

// Save in memory initial NumberofPseudoType tables.
int **NumberOfPseudoAtomsReferenceType;
int *NumberOfPseudoAtomsReferenceTypeNew;
int *NumberOfPseudoAtomsReferenceTypeOld;

// Save in memory the degree of freedom
int DegreesOfFreedomReferenceAdsorbates;
int DegreesOfFreedomReferenceTranslation;
int DegreesOfFreedomReferenceTranslationalAdsorbates;
int DegreesOfFreedomReference;
                                  
int DegreesOfFreedomReferenceRotation;
int DegreesOfFreedomReferenceAdsorbates;
int DegreesOfFreedomReferenceRotationalAdsorbates;

// Declare the chemical potential os the osmostat.
REAL *AlchemicalWorkStore;
int SizeAlchemicalWorkStore;

// for the NVE move during the Alchemical transformation.
REAL *HybridNVEAlchDrift;
REAL *HybridNVEAlchDriftCount;

REAL *HybridNVEAlchStartTemperature;
REAL *HybridNVEAlchStartTranslationalTemperature;
REAL *HybridNVEAlchStartRotationalTemperature;
REAL *HybridNVEAlchStartTemperatureFramework;
REAL *HybridNVEAlchStartTemperatureAdsorbate;
REAL *HybridNVEAlchStartTemperatureCation;
REAL *HybridNVEAlchStartTemperatureCount;
REAL *HybridNVEAlchStartTemperatureTranslationCount;
REAL *HybridNVEAlchStartTemperatureRotationCount;
REAL *HybridNVEAlchStartTemperatureFrameworkCount;
REAL *HybridNVEAlchStartTemperatureAdsorbateCount;
REAL *HybridNVEAlchStartTemperatureCationCount;

REAL *HybridNVEAlchEndTemperature;
REAL *HybridNVEAlchEndTranslationalTemperature;
REAL *HybridNVEAlchEndRotationalTemperature;
REAL *HybridNVEAlchEndTemperatureFramework;
REAL *HybridNVEAlchEndTemperatureAdsorbate;
REAL *HybridNVEAlchEndTemperatureCation;
REAL *HybridNVEAlchEndTemperatureCount;
REAL *HybridNVEAlchEndTemperatureTranslationCount;
REAL *HybridNVEAlchEndTemperatureRotationCount;
REAL *HybridNVEAlchEndTemperatureFrameworkCount;
REAL *HybridNVEAlchEndTemperatureAdsorbateCount;
REAL *HybridNVEAlchEndTemperatureCationCount;

// We'll need the energy before and after the set up of the transient molecules.
static REAL EnergyHostVDWFirstBead;
static REAL EnergyAdsorbateVDWFirstBead;
static REAL EnergyCationVDWFirstBead;
static REAL EnergyHostChargeChargeFirstBead;
static REAL EnergyAdsorbateChargeChargeFirstBead;
static REAL EnergyCationChargeChargeFirstBead;
static REAL EnergyHostChargeBondDipoleFirstBead;
static REAL EnergyAdsorbateChargeBondDipoleFirstBead;
static REAL EnergyCationChargeBondDipoleFirstBead;
static REAL EnergyHostBondDipoleBondDipoleFirstBead;
static REAL EnergyAdsorbateBondDipoleBondDipoleFirstBead;
static REAL EnergyCationBondDipoleBondDipoleFirstBead;

// We adapt the Handlefirstbead of cbmc.c to Initialize first beads of transient component.
static void ManageFirstBead(int start);
static void ManageRemainingBeads(int NumberOldComponent);

// File to write Alchemical work thorugh MC.
FILE *OutputOsmostatFilePtr;

/*********************************************************************************************************
 * Name       | InitializeStoreAlchemicalWork   (Added by A. de Izarra)                                  *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Initialize vector to store Alchemical work.                      						 *
 * Parameters | No parameters											                                 *
 *********************************************************************************************************/
void InitializeStoreAlchemicalWork(void)
{
	SizeAlchemicalWorkStore=0;
	AlchemicalWorkStore=(REAL*)calloc(SizeAlchemicalWorkStore,sizeof(REAL));
}

/*********************************************************************************************************
 * Name       | IncreaseSizeAlchemicalWork      (Added by A. de Izarra)                                  *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Increase size of vector to store Alchemical work on the fly.                           	 *
 * Parameters | No parameters											                                 *
 *********************************************************************************************************/
void IncreaseSizeAlchemicalWork(void)
{
	AlchemicalWorkStore=(REAL*)realloc(AlchemicalWorkStore,SizeAlchemicalWorkStore*sizeof(REAL));
}

/*********************************************************************************************************
 * Name       | InitializeFileAlchemicalWork      (Added by A. de Izarra)                                *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Initialize vector to store Alchemical work.                        					 	 *
 * Parameters | int Size => Number of file to store Alchemical work (=1)								 *										                                 		 *
 *********************************************************************************************************/
void InitializeFileAlchemicalWork(int Size)
{	
	int i;
		
	if(ContinueAfterCrash==FALSE) 
	{
		OutputOsmostatFilePtr=fopen("AlchemicalWork.txt","w");
		fprintf(OutputOsmostatFilePtr,"Widom_step\t\tAlchemical_work(k)\t\tAlchemical_work(kJ/mol)\n");	
		fflush(OutputOsmostatFilePtr);
	}
	else
	{
		// Read the osmostat file.
		FILE * BufferOutputOsmostatFilePtr = fopen("AlchemicalWorkbuffer.txt","w");
		OutputOsmostatFilePtr=fopen("AlchemicalWork.txt","r");
		
		char line[1024];
		
		for(i=0;i<=Size;i++)
		{
		  fgets(line,1024,OutputOsmostatFilePtr);
		  fprintf(BufferOutputOsmostatFilePtr,"%s",line);
		}	
		
		fclose(BufferOutputOsmostatFilePtr);
		fclose(OutputOsmostatFilePtr);
		
		remove("AlchemicalWork.txt");
		rename("AlchemicalWorkbuffer.txt","AlchemicalWork.txt");
		
		OutputOsmostatFilePtr=fopen("AlchemicalWork.txt","a");
	}		
}

/*********************************************************************************************************
 * Name       | UpdateFileAlchemicalWork      (Added by A. de Izarra)                                    *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Print in file FilePtr the index-th alchemical work and contribution from elec + vdw.     *                     					 
 * Parameters | int index => index of the computed alchemical potential									 *
 * 			  | FILE *FilePtr => Pointer to File (AlchemicalWork.txt)							 		 *
 *********************************************************************************************************/
void UpdateFileAlchemicalWork(int index, FILE *FilePtr)
{
	fprintf(FilePtr,"%d\t\t%lf\t\t%lf\n",index,AlchemicalWorkStore[index]*ENERGY_TO_KELVIN,AlchemicalWorkStore[index]*ENERGY_TO_KELVIN*KELVIN_TO_KJ_PER_MOL);
	fflush(OutputOsmostatFilePtr);
}

/*********************************************************************************************************
 * Name       | InitializeNVEAlchStatistics     (Added by A. de Izarra)                                  *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Initialize Statistics for the NVE MD during NCMC move.								     *                     					 
 * Parameters | No parameters																			 *
 *********************************************************************************************************/
void InitializeNVEAlchStatistics(void)
{
  int i;	
	
  HybridNVEAlchDrift=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNVEAlchDriftCount=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

  HybridNVEAlchStartTemperature=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNVEAlchStartTranslationalTemperature=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNVEAlchStartRotationalTemperature=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNVEAlchStartTemperatureFramework=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNVEAlchStartTemperatureAdsorbate=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNVEAlchStartTemperatureCation=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNVEAlchStartTemperatureCount=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNVEAlchStartTemperatureTranslationCount=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNVEAlchStartTemperatureRotationCount=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNVEAlchStartTemperatureFrameworkCount=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNVEAlchStartTemperatureAdsorbateCount=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNVEAlchStartTemperatureCationCount=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

  HybridNVEAlchEndTemperature=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNVEAlchEndTranslationalTemperature=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNVEAlchEndRotationalTemperature=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNVEAlchEndTemperatureFramework=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNVEAlchEndTemperatureAdsorbate=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNVEAlchEndTemperatureCation=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNVEAlchEndTemperatureCount=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNVEAlchEndTemperatureTranslationCount=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNVEAlchEndTemperatureRotationCount=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNVEAlchEndTemperatureFrameworkCount=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNVEAlchEndTemperatureAdsorbateCount=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  HybridNVEAlchEndTemperatureCationCount=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
  
  for(i=0;i<NumberOfSystems;i++)
  {
		HybridNVEAlchDrift[i]=0.0;
		HybridNVEAlchDriftCount[i]=0.0;

		HybridNVEAlchStartTemperature[i]=0.0;
		HybridNVEAlchStartTranslationalTemperature[i]=0.0;
		HybridNVEAlchStartRotationalTemperature[i]=0.0;
		HybridNVEAlchStartTemperatureFramework[i]=0.0;
		HybridNVEAlchStartTemperatureAdsorbate[i]=0.0;
		HybridNVEAlchStartTemperatureCation[i]=0.0;
		HybridNVEAlchStartTemperatureCount[i]=0.0;
		HybridNVEAlchStartTemperatureTranslationCount[i]=0.0;
		HybridNVEAlchStartTemperatureRotationCount[i]=0.0;
		HybridNVEAlchStartTemperatureFrameworkCount[i]=0.0;
		HybridNVEAlchStartTemperatureAdsorbateCount[i]=0.0;
		HybridNVEAlchStartTemperatureCationCount[i]=0.0;

		HybridNVEAlchEndTemperature[i]=0.0;
		HybridNVEAlchEndTranslationalTemperature[i]=0.0;
		HybridNVEAlchEndRotationalTemperature[i]=0.0;
		HybridNVEAlchEndTemperatureFramework[i]=0.0;
		HybridNVEAlchEndTemperatureAdsorbate[i]=0.0;
		HybridNVEAlchEndTemperatureCation[i]=0.0;
		HybridNVEAlchEndTemperatureCount[i]=0.0;
		HybridNVEAlchEndTemperatureTranslationCount[i]=0.0;
		HybridNVEAlchEndTemperatureRotationCount[i]=0.0;
		HybridNVEAlchEndTemperatureFrameworkCount[i]=0.0;
		HybridNVEAlchEndTemperatureAdsorbateCount[i]=0.0;
		HybridNVEAlchEndTemperatureCationCount[i]=0.0;
  }
}

/*********************************************************************************************************
 * Name       | InitializeIndexManagingTransientMoities      (Added by A. de Izarra)                     *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Initialize index of the molecules that undergo alchemical transformation.				 *                     					 
 * Parameters | No parameters																			 *
 *********************************************************************************************************/
void InitializeIndexManagingTransientMoities(void)
{
	// The number of extra components.
	NumberExtraComponents = 2;
	// The number of extra pseudoatoms.
	NumberExtraPseudoAtoms = 2*Components[SolventIndex].NumberOfAtoms;
	
	// The number of initial pseudo atoms.
	InitialPseudoAtoms = NumberOfPseudoAtoms;
	
	// Save new transient component index
    IndexExtraComponent[0]=NumberOfComponents;
    IndexExtraComponent[1]=NumberOfComponents+1;
    
    // Store the number of atoms of a transient molecule.
    TransientMoleculeNbAtoms = Components[SolventIndex].NumberOfAtoms; 
    
    // Chemical potential conversion in RASPA units.
    ChemicalPotentialAlchemical *= KJ_PER_MOL_TO_ENERGY;
}

/*********************************************************************************************************
 * Name       | AddExtraPseudoAtoms      (Added by A. de Izarra)                                   		 *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Add the extra pseudo atoms belonging to the molecules that 								 *
 * 			  | undergo alchemical transformation.													     *                     					 
 * Parameters | No parameters																			 *
 *********************************************************************************************************/
void AddExtraPseudoAtoms(void)
{       
	int i,j;

	// Allocate memory to store extra pseudoatoms of transient molecules.
	PseudoAtoms=(PSEUDO_ATOM*)realloc(PseudoAtoms,(InitialPseudoAtoms+NumberExtraPseudoAtoms)*sizeof(PSEUDO_ATOM));
  
	if(!PseudoAtoms) 
    {
      printf("Memory reallocation error of 'PseudoAtoms' in 'molecule.c'\n");
      exit(-1);
    }

    NumberOfPseudoAtomsTypeNew=(int*)realloc(NumberOfPseudoAtomsTypeNew,(InitialPseudoAtoms+NumberExtraPseudoAtoms)*sizeof(int));
    if(!NumberOfPseudoAtomsTypeNew)
    {
      printf("Memory reallocation error of 'NumberOfPseudoAtomsTypeNew' in file %s line %d\n", __FILE__, __LINE__);
      exit(-1);
    }

    NumberOfPseudoAtomsTypeOld=(int*)realloc(NumberOfPseudoAtomsTypeOld,(InitialPseudoAtoms+NumberExtraPseudoAtoms)*sizeof(int));
    if(!NumberOfPseudoAtomsTypeOld)
    {
      printf("Memory reallocation error of 'NumberOfPseudoAtomsTypeOld' in file %s line %d\n", __FILE__, __LINE__);
      exit(-1);
    }


    MapPseudoAtom=(int*)realloc(MapPseudoAtom,(InitialPseudoAtoms+NumberExtraPseudoAtoms)*sizeof(int));
    if(!MapPseudoAtom)
    {
      printf("Memory reallocation error of 'MapPseudoAtom' in file %s line %d\n", __FILE__, __LINE__);
      exit(-1);
    }

	for(j=InitialPseudoAtoms; j<InitialPseudoAtoms+NumberExtraPseudoAtoms; j++)
	{
		NumberOfPseudoAtomsTypeNew[j]=0;
		NumberOfPseudoAtomsTypeOld[j]=0;
		MapPseudoAtom[j]=0;		
	}
	
    for(i=0;i<NumberOfSystems;i++)
	{
		 NumberOfPseudoAtomsCount[i]=(int*)realloc(NumberOfPseudoAtomsCount[i],(InitialPseudoAtoms+NumberExtraPseudoAtoms)*sizeof(int));
		 if(!NumberOfPseudoAtomsCount[i])
		 {
		   printf("Memory reallocation error of 'NumberOfPseudoAtomsCount[i]'\n");
		   exit(-1);
		 }

		 NumberOfPseudoAtomsType[i]=(int*)realloc(NumberOfPseudoAtomsType[i],(InitialPseudoAtoms+NumberExtraPseudoAtoms)*sizeof(int));
		 if(!NumberOfPseudoAtomsType[i])
		 {
		   printf("Memory reallocation error of 'NumberOfPseudoAtomsType[i]'\n");
		   exit(-1);
		 }

		 NumberOfFractionalPseudoAtomsType[i]=(int*)realloc(NumberOfFractionalPseudoAtomsType[i],(InitialPseudoAtoms+NumberExtraPseudoAtoms)*sizeof(int));
		 if(!NumberOfFractionalPseudoAtomsType[i])
		 {
		   printf("Memory reallocation error of 'NumberOfFractionalPseudoAtomsType[i]'\n");
		   exit(-1);
		 }
		 
		 for(j=InitialPseudoAtoms; j<InitialPseudoAtoms+NumberExtraPseudoAtoms; j++)
		 {
			 NumberOfPseudoAtomsCount[i][j]=0;
			 NumberOfPseudoAtomsType[i][j]=0;
			 NumberOfFractionalPseudoAtomsType[i][j]=0;
		 }  
	}
	
 
    // For the new pseudoatom of transient molecules, assign its characteristics
    int Index;
    int atom_count;
    char index_salt_store[10];
    int comp;
   
    // The number of pseudo component to be present in the system
    for(comp=0;comp<2;comp++)
	{
		atom_count = -1;
		Index = (InitialPseudoAtoms) + comp*(Components[SolventIndex].NumberOfAtoms);
		for(i=Index;i<Index+Components[SolventIndex].NumberOfAtoms;i++)
		{
			atom_count++; 
			strcpy(PseudoAtoms[i].Name,"TransientAtom");
			sprintf(index_salt_store,"%d",IndexExtraComponent[comp]);
			strcat(PseudoAtoms[i].Name,index_salt_store);
			strcat(PseudoAtoms[i].Name,"-");
			sprintf(index_salt_store,"%d",atom_count);
			strcat(PseudoAtoms[i].Name,index_salt_store);
			strcpy(PseudoAtoms[i].ChemicalElement,PseudoAtoms[i].Name);
			
			strcpy(PseudoAtoms[i].OxidationStateString,PseudoAtoms[Components[SolventIndex].Type[atom_count]].OxidationStateString);
			PseudoAtoms[i].OxidationState=PseudoAtoms[Components[SolventIndex].Type[atom_count]].OxidationState;
			strcpy(PseudoAtoms[i].ScatteringSource,PseudoAtoms[Components[SolventIndex].Type[atom_count]].ScatteringSource);
		    PseudoAtoms[i].Occupancy=PseudoAtoms[Components[SolventIndex].Type[atom_count]].Occupancy;
		    PseudoAtoms[i].FrameworkAtom=PseudoAtoms[Components[SolventIndex].Type[atom_count]].FrameworkAtom; // This is not framework atom
		    PseudoAtoms[i].PrintToPDB=FALSE;
		    PseudoAtoms[i].ScatteringType=PseudoAtoms[Components[SolventIndex].Type[atom_count]].ScatteringType;
		    PseudoAtoms[i].AnomalousScatteringType=PseudoAtoms[Components[SolventIndex].Type[atom_count]].AnomalousScatteringType;
		    PseudoAtoms[i].TemperatureFactor=PseudoAtoms[Components[SolventIndex].Type[atom_count]].TemperatureFactor;
		    PseudoAtoms[i].AnisotropicTemperatureFactor.ax=PseudoAtoms[Components[SolventIndex].Type[atom_count]].AnisotropicTemperatureFactor.ax;
		    PseudoAtoms[i].AnisotropicTemperatureFactor.ay=PseudoAtoms[Components[SolventIndex].Type[atom_count]].AnisotropicTemperatureFactor.ay;
		    PseudoAtoms[i].AnisotropicTemperatureFactor.az=PseudoAtoms[Components[SolventIndex].Type[atom_count]].AnisotropicTemperatureFactor.az;
		    PseudoAtoms[i].AnisotropicTemperatureFactor.bx=PseudoAtoms[Components[SolventIndex].Type[atom_count]].AnisotropicTemperatureFactor.bx;
		    PseudoAtoms[i].AnisotropicTemperatureFactor.by=PseudoAtoms[Components[SolventIndex].Type[atom_count]].AnisotropicTemperatureFactor.by;
		    PseudoAtoms[i].AnisotropicTemperatureFactor.bz=PseudoAtoms[Components[SolventIndex].Type[atom_count]].AnisotropicTemperatureFactor.bz;
		    PseudoAtoms[i].AnisotropicTemperatureFactor.cx=PseudoAtoms[Components[SolventIndex].Type[atom_count]].AnisotropicTemperatureFactor.cx;
		    PseudoAtoms[i].AnisotropicTemperatureFactor.cy=PseudoAtoms[Components[SolventIndex].Type[atom_count]].AnisotropicTemperatureFactor.cy;
		    PseudoAtoms[i].AnisotropicTemperatureFactor.cz=PseudoAtoms[Components[SolventIndex].Type[atom_count]].AnisotropicTemperatureFactor.cz;
		    PseudoAtoms[i].ScatteringDispersionImaginary=PseudoAtoms[Components[SolventIndex].Type[atom_count]].ScatteringDispersionImaginary;
		    PseudoAtoms[i].Mass=0; // To be assigned in another function.
		    PseudoAtoms[i].Charge1=0; // To be assigned in another function.
		    PseudoAtoms[i].ChargeDefinitionType=PseudoAtoms[Components[SolventIndex].Type[atom_count]].ChargeDefinitionType;
		    PseudoAtoms[i].Polarization.ax=PseudoAtoms[Components[SolventIndex].Type[atom_count]].Polarization.ax;
		    PseudoAtoms[i].Polarization.ay=PseudoAtoms[Components[SolventIndex].Type[atom_count]].Polarization.ay;
		    PseudoAtoms[i].Polarization.az=PseudoAtoms[Components[SolventIndex].Type[atom_count]].Polarization.az;
		    PseudoAtoms[i].Polarization.bx=PseudoAtoms[Components[SolventIndex].Type[atom_count]].Polarization.bx;
		    PseudoAtoms[i].Polarization.by=PseudoAtoms[Components[SolventIndex].Type[atom_count]].Polarization.by;
		    PseudoAtoms[i].Polarization.bz=PseudoAtoms[Components[SolventIndex].Type[atom_count]].Polarization.bz;
		    PseudoAtoms[i].Polarization.cx=PseudoAtoms[Components[SolventIndex].Type[atom_count]].Polarization.cx;
		    PseudoAtoms[i].Polarization.cy=PseudoAtoms[Components[SolventIndex].Type[atom_count]].Polarization.cy;
		    PseudoAtoms[i].Polarization.cz=PseudoAtoms[Components[SolventIndex].Type[atom_count]].Polarization.cz;
		    PseudoAtoms[i].HasCharges=TRUE; // Always has charge.
		    PseudoAtoms[i].IsPolarizable=PseudoAtoms[Components[SolventIndex].Type[atom_count]].IsPolarizable;
		    PseudoAtoms[i].Interaction=PseudoAtoms[Components[SolventIndex].Type[atom_count]].Interaction;
		    PseudoAtoms[i].Radius=1.0; // TO be assigned in another function.
		    PseudoAtoms[i].Connectivity=PseudoAtoms[Components[SolventIndex].Type[atom_count]].Connectivity;
		    PseudoAtoms[i].TinkerType=PseudoAtoms[Components[SolventIndex].Type[atom_count]].TinkerType;
		    PseudoAtoms[i].AnisotropicCorrection=PseudoAtoms[Components[SolventIndex].Type[atom_count]].AnisotropicCorrection;
		    PseudoAtoms[i].AnisotropicDisplacement=PseudoAtoms[Components[SolventIndex].Type[atom_count]].AnisotropicDisplacement;
		    PseudoAtoms[i].AnisotropicType=PseudoAtoms[Components[SolventIndex].Type[atom_count]].AnisotropicType;
		    PseudoAtoms[i].HasVDWInteraction=PseudoAtoms[Components[SolventIndex].Type[atom_count]].HasVDWInteraction;
		    PseudoAtoms[i].Hybridization=PseudoAtoms[Components[SolventIndex].Type[atom_count]].Hybridization;
		    PseudoAtoms[i].CF=PseudoAtoms[Components[SolventIndex].Type[atom_count]].CF;
		    
		}
    }	
    
   NumberOfPseudoAtoms=InitialPseudoAtoms+NumberExtraPseudoAtoms;
}


void GrowTransient(int NumberOldComponent)
{
  int j,start;
  REAL UVDWCorrectionAdsorbate,UVDWCorrectionCation;
  REAL UVDWCorrectionFramework,UVDWCorrectionReplicasNew;
  REAL UChargeChargeCorrectionReplicasNew;

  UBondNew[CurrentSystem]=0.0;
  UBendNew[CurrentSystem]=0.0;
  UBendBendNew[CurrentSystem]=0.0;
  UInversionBendNew[CurrentSystem]=0.0;
  UUreyBradleyNew[CurrentSystem]=0.0;
  UTorsionNew[CurrentSystem]=0.0;
  UImproperTorsionNew[CurrentSystem]=0.0;
  UBondBondNew[CurrentSystem]=0.0;
  UBondBendNew[CurrentSystem]=0.0;
  UBondTorsionNew[CurrentSystem]=0.0;
  UBendTorsionNew[CurrentSystem]=0.0;
  UIntraVDWNew[CurrentSystem]=0.0;
  UIntraChargeChargeNew[CurrentSystem]=0.0;
  UIntraChargeBondDipoleNew[CurrentSystem]=0.0;
  UIntraBondDipoleBondDipoleNew[CurrentSystem]=0.0;

  UHostVDWNew[CurrentSystem]=0.0;
  UAdsorbateVDWNew[CurrentSystem]=0.0;
  UCationVDWNew[CurrentSystem]=0.0;
  UHostChargeChargeNew[CurrentSystem]=0.0;
  UAdsorbateChargeChargeNew[CurrentSystem]=0.0;
  UCationChargeChargeNew[CurrentSystem]=0.0;
  UHostChargeBondDipoleNew[CurrentSystem]=0.0;
  UAdsorbateChargeBondDipoleNew[CurrentSystem]=0.0;
  UCationChargeBondDipoleNew[CurrentSystem]=0.0;
  UHostBondDipoleBondDipoleNew[CurrentSystem]=0.0;
  UAdsorbateBondDipoleBondDipoleNew[CurrentSystem]=0.0;
  UCationBondDipoleBondDipoleNew[CurrentSystem]=0.0;

  RosenbluthNew=1.0;
  OVERLAP=FALSE;

  if(NumberOfBeadsAlreadyPlaced==0)
  {
	start=Components[CurrentComponent].StartingBead;
	FirstBeadPosition=NewPosition[CurrentSystem][start];
	
    ManageFirstBead(start);

    UCationVDWNew[CurrentSystem]=EnergyCationVDWFirstBead;
    UAdsorbateVDWNew[CurrentSystem]=EnergyAdsorbateVDWFirstBead;
    UHostVDWNew[CurrentSystem]=EnergyHostVDWFirstBead;

    UCationChargeChargeNew[CurrentSystem]=EnergyCationChargeChargeFirstBead;
    UAdsorbateChargeChargeNew[CurrentSystem]=EnergyAdsorbateChargeChargeFirstBead;
    UHostChargeChargeNew[CurrentSystem]=EnergyHostChargeChargeFirstBead;

    UCationChargeBondDipoleNew[CurrentSystem]=EnergyCationChargeBondDipoleFirstBead;
    UAdsorbateChargeBondDipoleNew[CurrentSystem]=EnergyAdsorbateChargeBondDipoleFirstBead;
    UHostChargeBondDipoleNew[CurrentSystem]=EnergyHostChargeBondDipoleFirstBead;

    UCationBondDipoleBondDipoleNew[CurrentSystem]=0.0;
    UAdsorbateBondDipoleBondDipoleNew[CurrentSystem]=0.0;
    UHostBondDipoleBondDipoleNew[CurrentSystem]=0.0;

    NumberOfBeadsAlreadyPlaced=1;
  }

  if(Components[CurrentComponent].NumberOfAtoms>1)
  {
    ManageRemainingBeads(NumberOldComponent);
  }

  // copy coordinates for small MC-scheme
  for(j=0;j<Components[CurrentComponent].NumberOfAtoms;j++)
  {
    Components[CurrentComponent].RMCMOL[j]=NewPosition[CurrentSystem][j];
    TrialPosition[CurrentSystem][j]=NewPosition[CurrentSystem][j];
  }

  // calculate anisotropic sites
  CalculateAnisotropicTrialPositions(CurrentComponent,TrialPosition[CurrentSystem],TrialAnisotropicPosition[CurrentSystem]);
  UVDWCorrectionFramework=CalculateFrameworkVDWEnergyCorrection(TrialPosition[CurrentSystem],TrialAnisotropicPosition[CurrentSystem],CFVDWScaling);
  UVDWCorrectionAdsorbate=CalculateInterVDWEnergyCorrectionAdsorbate(TrialPosition[CurrentSystem],TrialAnisotropicPosition[CurrentSystem],CurrentAdsorbateMolecule);
  UVDWCorrectionCation=CalculateInterVDWEnergyCorrectionCation(TrialPosition[CurrentSystem],TrialAnisotropicPosition[CurrentSystem],CurrentCationMolecule);
  RosenbluthNew*=exp(-Beta[CurrentSystem]*(UVDWCorrectionFramework+UVDWCorrectionAdsorbate+UVDWCorrectionCation));

  UHostVDWNew[CurrentSystem]+=UVDWCorrectionFramework;
  UAdsorbateVDWNew[CurrentSystem]+=UVDWCorrectionAdsorbate;
  UCationVDWNew[CurrentSystem]+=UVDWCorrectionCation;

  // correct for self-energy when using replica unit cells
  UVDWCorrectionReplicasNew=CalculateInterVDWSelfEnergyCorrectionNew();
  UChargeChargeCorrectionReplicasNew=CalculateInterChargeChargeSelfEnergyCorrectionNew();
  if(Components[CurrentComponent].ExtraFrameworkMolecule)
  {
    UCationVDWNew[CurrentSystem]+=UVDWCorrectionReplicasNew;
    UCationChargeChargeNew[CurrentSystem]+=UChargeChargeCorrectionReplicasNew;
  }
  else
  {
    UAdsorbateVDWNew[CurrentSystem]+=UVDWCorrectionReplicasNew;
    UAdsorbateChargeChargeNew[CurrentSystem]+=UChargeChargeCorrectionReplicasNew;
  }

  // the old config can be used as a starting point
  Components[CurrentComponent].LMCMOL=TRUE;
}


/*********************************************************************************************************
 * Name       | MakeInitialTransient      (Added by A. de Izarra)                                    	 *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Construct the molecules that undergo alchemical transformation.				     		 *                     					 
 * Parameters | int	NumberOldComponent => determines if old component is water (NumberOldComponent=1)    *
 * 			  |	to be transformed into ions or ions (NumberOldComponent=2) to be transformed into water. *												
 *********************************************************************************************************/
void MakeInitialTransient(int NumberOldComponent)
{	
  int Comp;
  int NumberTransientMoities=-1;
  int i,j;

  REAL UTailNew;
  REAL RosenbluthIdealNew;
  REAL DeltaU;
  
  int StoredNumberOfTrialPositions;
  int StoredNumberOfTrialPositionsFirstBead;

  StoredNumberOfTrialPositions=NumberOfTrialPositions;
  StoredNumberOfTrialPositionsFirstBead=NumberOfTrialPositionsForTheFirstBead;

  for(Comp=0; Comp<NumberExtraComponents; Comp++)
  {
	    UpdateInertiaTensorGroups(IndexExtraComponent[Comp]);
	    
		for(i=0;i<MultiplicitySalt[CurrentAlchemicalReaction][Comp];i++)
	    {  
			NumberTransientMoities++;

			// The current transient moities is added at the top of the adsorbates pile;
			CurrentAdsorbateMolecule=NumberOfAdsorbateMolecules[CurrentSystem];
			
		    CurrentComponent=IndexExtraComponent[Comp];
		    
		    // set Continuous Fraction (CF) atomic scaling-factors to  unity (only integer molecules can be added here)
			for(j=0;j<Components[CurrentComponent].NumberOfAtoms;j++)
			{
				CFVDWScaling[j]=1.0;
				CFChargeScaling[j]=1.0;
			}

			// There is no need to regrow the molecule: past water coordinates into the new transient component.
			// water->ions
			
			NumberOfBeadsAlreadyPlaced=0;
			NumberOfTrialPositions=NumberOfTrialPositionsSwap;
			NumberOfTrialPositionsForTheFirstBead=NumberOfTrialPositionsForTheFirstBeadSwap;
			
			if(NumberOldComponent == 1)
			{
				// Paste in new position the coordinate of water.
				for(j=0;j<Components[CurrentComponent].NumberOfAtoms;j++)
				{
					NewPosition[CurrentSystem][j]=AdsorbatesReferenceChosen[NumberTransientMoities].Atoms[j].Position;
				}
			}
			else
			{
				NewPosition[CurrentSystem][0]=AdsorbatesReferenceChosen[NumberTransientMoities].Atoms[0].Position;
				for(j=1;j<Components[CurrentComponent].NumberOfAtoms;j++)
				{
					NewPosition[CurrentSystem][j].x=0.0;
					NewPosition[CurrentSystem][j].y=0.0;
					NewPosition[CurrentSystem][j].z=0.0;
				}
			}

			GrowTransient(NumberOldComponent);						

		    for(j=0;j<Components[CurrentComponent].NumberOfAtoms;j++)
		    {
			  if(BlockedPocket(TrialPosition[CurrentSystem][j]))
			  return;
		    }

		    UTailNew=TailMolecularEnergyDifferenceAdd();

		    if(ComputePolarization)
		    {
			    ComputeNewPolarizationEnergy(TRUE,CurrentAdsorbateMolecule,-1);
		    }

		    if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
		    {
			    CalculateEwaldFourierAdsorbate(TRUE,FALSE,NumberOfAdsorbateMolecules[CurrentSystem],0);
			}

			UAdsorbateBond[CurrentSystem]+=UBondNew[CurrentSystem];
			UAdsorbateUreyBradley[CurrentSystem]+=UUreyBradleyNew[CurrentSystem];
			UAdsorbateBend[CurrentSystem]+=UBendNew[CurrentSystem];
			UAdsorbateBendBend[CurrentSystem]+=UBendBendNew[CurrentSystem];
			UAdsorbateInversionBend[CurrentSystem]+=UInversionBendNew[CurrentSystem];
			UAdsorbateTorsion[CurrentSystem]+=UTorsionNew[CurrentSystem];
			UAdsorbateImproperTorsion[CurrentSystem]+=UImproperTorsionNew[CurrentSystem];
			UAdsorbateBondBond[CurrentSystem]+=UBondBondNew[CurrentSystem];
			UAdsorbateBondBend[CurrentSystem]+=UBondBendNew[CurrentSystem];
			UAdsorbateBondTorsion[CurrentSystem]+=UBondTorsionNew[CurrentSystem];
			UAdsorbateBendTorsion[CurrentSystem]+=UBendTorsionNew[CurrentSystem];
			UAdsorbateIntraVDW[CurrentSystem]+=UIntraVDWNew[CurrentSystem];

			UAdsorbateAdsorbate[CurrentSystem]+=UAdsorbateVDWNew[CurrentSystem];
			UAdsorbateAdsorbateVDW[CurrentSystem]+=UAdsorbateVDWNew[CurrentSystem];
			UAdsorbateCation[CurrentSystem]+=UCationVDWNew[CurrentSystem];
			UAdsorbateCationVDW[CurrentSystem]+=UCationVDWNew[CurrentSystem];
			UHostAdsorbate[CurrentSystem]+=UHostVDWNew[CurrentSystem];
			UHostAdsorbateVDW[CurrentSystem]+=UHostVDWNew[CurrentSystem];

			UTailCorrection[CurrentSystem]+=UTailNew;

			UHostPolarization[CurrentSystem]=UHostPolarizationNew[CurrentSystem];
			UAdsorbatePolarization[CurrentSystem]=UAdsorbatePolarizationNew[CurrentSystem]+UPolarizationNew[CurrentSystem];
			UCationPolarization[CurrentSystem]=UCationPolarizationNew[CurrentSystem];

			UHostBackPolarization[CurrentSystem]=UHostBackPolarizationNew[CurrentSystem];
			UAdsorbateBackPolarization[CurrentSystem]=UAdsorbateBackPolarizationNew[CurrentSystem]+UBackPolarizationNew[CurrentSystem];
			UCationBackPolarization[CurrentSystem]=UCationBackPolarizationNew[CurrentSystem];

			if(ChargeMethod!=NONE)
			{
			  UAdsorbateIntraChargeCharge[CurrentSystem]+=UIntraChargeChargeNew[CurrentSystem];
			  UAdsorbateIntraChargeBondDipole[CurrentSystem]+=UIntraChargeBondDipoleNew[CurrentSystem];
			  UAdsorbateIntraBondDipoleBondDipole[CurrentSystem]+=UIntraBondDipoleBondDipoleNew[CurrentSystem];

			  UAdsorbateAdsorbateChargeChargeReal[CurrentSystem]+=UAdsorbateChargeChargeNew[CurrentSystem];
			  UAdsorbateAdsorbateChargeBondDipoleReal[CurrentSystem]+=UAdsorbateChargeBondDipoleNew[CurrentSystem];
			  UAdsorbateAdsorbateBondDipoleBondDipoleReal[CurrentSystem]+=UAdsorbateBondDipoleBondDipoleNew[CurrentSystem];
			  UAdsorbateAdsorbateChargeChargeFourier[CurrentSystem]+=UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem];
			  UAdsorbateAdsorbateChargeBondDipoleFourier[CurrentSystem]+=UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem];
			  UAdsorbateAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]+=UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem];
			  UAdsorbateAdsorbateCoulomb[CurrentSystem]+=UAdsorbateChargeChargeNew[CurrentSystem]+
														 UAdsorbateChargeBondDipoleNew[CurrentSystem]+
														 UAdsorbateBondDipoleBondDipoleNew[CurrentSystem]+
														 UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem]+
														 UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+
														 UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem];
			  UAdsorbateAdsorbate[CurrentSystem]+=UAdsorbateChargeChargeNew[CurrentSystem]+
												  UAdsorbateChargeBondDipoleNew[CurrentSystem]+
												  UAdsorbateBondDipoleBondDipoleNew[CurrentSystem]+
												  UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem]+
												  UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+
												  UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem];

			  UHostAdsorbateChargeChargeReal[CurrentSystem]+=UHostChargeChargeNew[CurrentSystem];
			  UHostAdsorbateChargeBondDipoleReal[CurrentSystem]+=UHostChargeBondDipoleNew[CurrentSystem];
			  UHostAdsorbateBondDipoleBondDipoleReal[CurrentSystem]+=UHostBondDipoleBondDipoleNew[CurrentSystem];
			  UHostAdsorbateChargeChargeFourier[CurrentSystem]+=UHostAdsorbateChargeChargeFourierDelta[CurrentSystem];
			  UHostAdsorbateChargeBondDipoleFourier[CurrentSystem]+=UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem];
			  UHostAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]+=UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem];
			  UHostAdsorbateCoulomb[CurrentSystem]+=UHostChargeChargeNew[CurrentSystem]+
													UHostChargeBondDipoleNew[CurrentSystem]+
													UHostBondDipoleBondDipoleNew[CurrentSystem]+
													UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]+
													UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+
													UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem];
			  UHostAdsorbate[CurrentSystem]+=UHostChargeChargeNew[CurrentSystem]+
											 UHostChargeBondDipoleNew[CurrentSystem]+
											 UHostBondDipoleBondDipoleNew[CurrentSystem]+
											 UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]+
											 UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+
											 UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem];

			  UAdsorbateCationChargeChargeReal[CurrentSystem]+=UCationChargeChargeNew[CurrentSystem];
			  UAdsorbateCationChargeBondDipoleReal[CurrentSystem]+=UCationChargeBondDipoleNew[CurrentSystem];
			  UAdsorbateCationBondDipoleBondDipoleReal[CurrentSystem]+=UCationBondDipoleBondDipoleNew[CurrentSystem];
			  UAdsorbateCationChargeChargeFourier[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem];
			  UAdsorbateCationChargeBondDipoleFourier[CurrentSystem]+=UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem];
			  UAdsorbateCationBondDipoleBondDipoleFourier[CurrentSystem]+=UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem];
			  UAdsorbateCationCoulomb[CurrentSystem]+=UCationChargeChargeNew[CurrentSystem]+
													  UCationChargeBondDipoleNew[CurrentSystem]+
													  UCationBondDipoleBondDipoleNew[CurrentSystem]+
													  UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
													  UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
													  UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem];
			  UAdsorbateCation[CurrentSystem]+=UCationChargeChargeNew[CurrentSystem]+
											   UCationChargeBondDipoleNew[CurrentSystem]+
											   UCationBondDipoleBondDipoleNew[CurrentSystem]+
											   UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
											   UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
											   UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem];

			  NetChargeAdsorbates[CurrentSystem]+=NetChargeAdsorbateDelta;
			  NetChargeSystem[CurrentSystem]+=NetChargeAdsorbateDelta;

			  if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
				AcceptEwaldAdsorbateMove(0);
			}
			
			InsertAdsorbateAlchMolecule();

			DeltaU=UBondNew[CurrentSystem]+UUreyBradleyNew[CurrentSystem]+UBendNew[CurrentSystem]+UBendBendNew[CurrentSystem]+UInversionBendNew[CurrentSystem]+UTorsionNew[CurrentSystem]+
				   UImproperTorsionNew[CurrentSystem]+UBondBondNew[CurrentSystem]+UBondBendNew[CurrentSystem]+UBondTorsionNew[CurrentSystem]+UBendTorsionNew[CurrentSystem]+
				   UIntraVDWNew[CurrentSystem]+UIntraChargeChargeNew[CurrentSystem]+
				   UIntraChargeBondDipoleNew[CurrentSystem]+UIntraBondDipoleBondDipoleNew[CurrentSystem]+
				   UAdsorbateVDWNew[CurrentSystem]+UCationVDWNew[CurrentSystem]+UHostVDWNew[CurrentSystem]+
				   UAdsorbateChargeChargeNew[CurrentSystem]+UCationChargeChargeNew[CurrentSystem]+UHostChargeChargeNew[CurrentSystem]+
				   UAdsorbateChargeBondDipoleNew[CurrentSystem]+UCationChargeBondDipoleNew[CurrentSystem]+UHostChargeBondDipoleNew[CurrentSystem]+
				   UAdsorbateBondDipoleBondDipoleNew[CurrentSystem]+UCationBondDipoleBondDipoleNew[CurrentSystem]+UHostBondDipoleBondDipoleNew[CurrentSystem]+
				   UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]+UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem]+UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
				   UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
				   UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
				   UDeltaPolarization+UTailNew;
					
				   UTotal[CurrentSystem]+=DeltaU;
  		}
  }

  return;
}

/*********************************************************************************************************
 * Name       | ReallocateMemoryParameterTab      (Added by A. de Izarra)                                *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Reallocate memory to store potential parameters of peudo atoms for the					 * 
 *            | molecules that undergo the alchemical transformation.						     		 *                     					 
 * Parameters | No parameters																			 *
 *********************************************************************************************************/
void ReallocateMemoryParameterTab(void)
{ 
	 int k,i,j;
	 //-----------------------------------------------
	 // allocate extra memory in PotentialParms, updated at each alchemical transformation step.
	 PotentialParms=(REAL(**)[MAX_NUMBER_OF_POTENTIAL_ARGUMENTS])realloc(PotentialParms,
	 (InitialPseudoAtoms+NumberExtraPseudoAtoms)*sizeof(REAL(*)[MAX_NUMBER_OF_POTENTIAL_ARGUMENTS]));
	
	 for(i = 0; i < InitialPseudoAtoms; i++)
	 {
		  PotentialParms[i] = (REAL(*)[MAX_NUMBER_OF_POTENTIAL_ARGUMENTS])realloc(PotentialParms[i],(InitialPseudoAtoms+NumberExtraPseudoAtoms)*sizeof(REAL[MAX_NUMBER_OF_POTENTIAL_ARGUMENTS]));
	 }

	 for(i = InitialPseudoAtoms; i < InitialPseudoAtoms+NumberExtraPseudoAtoms; i++)
	 {
		  PotentialParms[i] = (REAL(*)[MAX_NUMBER_OF_POTENTIAL_ARGUMENTS])calloc(InitialPseudoAtoms+NumberExtraPseudoAtoms,sizeof(REAL[MAX_NUMBER_OF_POTENTIAL_ARGUMENTS]));
	 }
	 
	 // Fill new initialized PotentialParms to zero.
	 for(i=0;i<InitialPseudoAtoms+NumberExtraPseudoAtoms;i++)
	 {
		  for(j=i+1;j<InitialPseudoAtoms+NumberExtraPseudoAtoms;j++)
		  {
			if((i<InitialPseudoAtoms)&&(j<InitialPseudoAtoms))
			{
				continue;
			}
			
			PotentialParms[i][j][0]=0.0;
			PotentialParms[i][j][1]=0.0;
			PotentialParms[i][j][2]=0.0;
		  }
	 }

	 //-----------------------------------------------
	 // allocate extra memory in Potentialtyoe
	 if(PotentialType[0][0] != LENNARD_JONES)
	 {
		 fprintf(stderr,"Error: Alchemical transformation is set up for LENNARD-JONES interaction.\n");
		 exit(0);
	 }
	 else
	 {
		 // Step 3: allocate extra memory in PotentialType.
		 PotentialType=(int**)realloc(PotentialType,(InitialPseudoAtoms+NumberExtraPseudoAtoms)*(sizeof(int*)));
		 
		 for(i = 0; i < InitialPseudoAtoms; i++)
		 {
			  PotentialType[i] = (int*)realloc(PotentialType[i],(InitialPseudoAtoms+NumberExtraPseudoAtoms)*sizeof(int));
		 }
		 
		 for(i = InitialPseudoAtoms; i < InitialPseudoAtoms+NumberExtraPseudoAtoms; i++)
		 {
			  PotentialType[i] = (int*)calloc(InitialPseudoAtoms+NumberExtraPseudoAtoms,sizeof(int));
		 }

		// Fill the extended potentialTypeTab.
		int index;
		int atom_nr;
		int atom_type;
		
		for(k=0;k<2;k++)
		{
			atom_nr=-1;
			index = (InitialPseudoAtoms) + k*(Components[SolventIndex].NumberOfAtoms);
			for(i=index;i<index+Components[SolventIndex].NumberOfAtoms;i++)
			{
				atom_nr++;
				atom_type = Components[SolventIndex].Type[atom_nr];
				PotentialType[i][i] = PotentialType[atom_type][atom_type];
			}
		}
		
		// PotentialTab updated according to mixing rules.
		for(i=0;i<InitialPseudoAtoms+NumberExtraPseudoAtoms;i++)
		{
		  for(j=i+1;j<InitialPseudoAtoms+NumberExtraPseudoAtoms;j++)
		  {
			if((i<InitialPseudoAtoms)&&(j<InitialPseudoAtoms))
			{
				continue;
			}
			
			if((PotentialType[i][i]==LENNARD_JONES)&&(PotentialType[j][j]==LENNARD_JONES))
			{	
			  PotentialType[i][j]=LENNARD_JONES;
			  PotentialType[j][i]=LENNARD_JONES;
			}
			else
			{
		      PotentialType[i][j]=ZERO_POTENTIAL;
			  PotentialType[j][i]=ZERO_POTENTIAL;
			}
		  }
		}
	}
	
	//-----------------------------------------------
	// Allocate extra memory for tail correction.
	TailCorrection=(int**)realloc(TailCorrection,(InitialPseudoAtoms+NumberExtraPseudoAtoms)*(sizeof(int*)));
	 
	 for(i = 0; i < InitialPseudoAtoms; i++)
	 {
		  TailCorrection[i] = (int*)realloc(TailCorrection[i],(InitialPseudoAtoms+NumberExtraPseudoAtoms)*sizeof(int));
	 }
	 
	 for(i = InitialPseudoAtoms; i < InitialPseudoAtoms+NumberExtraPseudoAtoms; i++)
	 {
		  TailCorrection[i] = (int*)calloc(InitialPseudoAtoms+NumberExtraPseudoAtoms,sizeof(int));
	 }

	// Fill the extended TailCorrection.
	int index;
	
	for(k=0;k<2;k++)
	{
		index = (InitialPseudoAtoms) + k*(Components[SolventIndex].NumberOfAtoms);
		for(i=index;i<index+Components[SolventIndex].NumberOfAtoms;i++)
		{
			TailCorrection[i][i] = TailCorrection[0][0];
		}
	}
	
	// PotentialTab updated according to mixing rules.
	for(i=0;i<InitialPseudoAtoms+NumberExtraPseudoAtoms;i++)
	{
	  for(j=i+1;j<InitialPseudoAtoms+NumberExtraPseudoAtoms;j++)
	  {
		if((i<InitialPseudoAtoms)&&(j<InitialPseudoAtoms))
		{
			continue;
		}
		
		  TailCorrection[i][j]=TailCorrection[0][0];
		  TailCorrection[j][i]=TailCorrection[0][0];
	  }
	}
	
	
	//-----------------------------------------------
	// Allocate extra memory for shift potential
	ShiftPotential=(int**)realloc(ShiftPotential,(InitialPseudoAtoms+NumberExtraPseudoAtoms)*(sizeof(int*)));
	 
	 for(i = 0; i < InitialPseudoAtoms; i++)
	 {
		  ShiftPotential[i] = (int*)realloc(ShiftPotential[i],(InitialPseudoAtoms+NumberExtraPseudoAtoms)*sizeof(int));
	 }
	 
	 for(i = InitialPseudoAtoms; i < InitialPseudoAtoms+NumberExtraPseudoAtoms; i++)
	 {
		  ShiftPotential[i] = (int*)calloc(InitialPseudoAtoms+NumberExtraPseudoAtoms,sizeof(int));
	 }

	// Fill the extended TailCorrection.
	
	for(k=0;k<2;k++)
	{
		index = (InitialPseudoAtoms) + k*(Components[SolventIndex].NumberOfAtoms);
		for(i=index;i<index+Components[SolventIndex].NumberOfAtoms;i++)
		{
			ShiftPotential[i][i] = ShiftPotential[0][0];
		}
	}
	
	// PotentialTab updated according to mixing rules.
	for(i=0;i<InitialPseudoAtoms+NumberExtraPseudoAtoms;i++)
	{
	  for(j=i+1;j<InitialPseudoAtoms+NumberExtraPseudoAtoms;j++)
	  {
		if((i<InitialPseudoAtoms)&&(j<InitialPseudoAtoms))
		{
			continue;
		}
		
		  ShiftPotential[i][j]=ShiftPotential[0][0];
		  ShiftPotential[j][i]=ShiftPotential[0][0];
	  }
    }	
}

/*********************************************************************************************************
 * Name       | InitializeVectorCharge      (Added by A. de Izarra)                               		 *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Initialize vectors to store charges of the pseudo atoms of the molecules that 			 *
 * 			  | undergo the alchemical transformation.										     		 *                     					 
 * Parameters | No parameters																			 *
 *********************************************************************************************************/
void InitializeVectorCharge(void)
{
	
	 // Step 1: preparation of vectors.
	 int Species,i;
	 
	 int type;
											
	 VectorChargeAlchemicalTransformation = (REAL(**))calloc(3,sizeof(REAL(*)));

	 for(Species=0; Species<3; Species++)
	 {	 
		 VectorChargeAlchemicalTransformation[Species] = (REAL(*))calloc(TransientMoleculeNbAtoms,sizeof(REAL));
		 
		 for(i=0;i<TransientMoleculeNbAtoms;i++)
		 {
			 // Manage the solvent.
			 if(Species == SolventIndex)
			 {
				 type = Components[SolventIndex].Type[i];
				 
				 VectorChargeAlchemicalTransformation[Species][i] = PseudoAtoms[type].Charge1;		 
			 }
			 // Manage the ions.
			 else
			 {
				 // for i==0, it is the first atom corresponding to the atomic ions.
				 if(i==0)
				 {
					type = Components[Species].Type[i];
					
					VectorChargeAlchemicalTransformation[Species][i] = PseudoAtoms[type].Charge1;	
				 }
				 else
				 {
					VectorChargeAlchemicalTransformation[Species][i]=0.0;
				 }
			 }
		 }
	 } 	 
}

/*********************************************************************************************************
 * Name       | InitializeChargeTransientMoities    (Added by A. de Izarra)                              *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Initialize charges of the pseudo atoms of the molecules that undergo the alchemical		 *
 * 			  | transformation.										     		 						 *
 * Parameters | int	NumberOldComponent => determines if old component is water (NumberOldComponent=1)    *
 * 			  |	to be transformed into ions or ions (NumberOldComponent=2) to be transformed into water. *
 *********************************************************************************************************/
void InitializeChargeTransientMoities(int NumberOldComponent)
{
	int k,i,j;

	// Index to put the proper parameter in the right place in PotentialParms
	int index;
	int atom_count;
	
	// Initialize pseudoatom charge
	
	// We know water has been chosen as old component
	if(NumberOldComponent == 1)
	{
		for(k=0;k<NumberExtraComponents;k++)
		{
			atom_count = -1;
			index = (InitialPseudoAtoms) + k*(TransientMoleculeNbAtoms);
			for(i=index;i<index+TransientMoleculeNbAtoms;i++)
			{
				atom_count++;
				PseudoAtoms[i].Charge1 = VectorChargeAlchemicalTransformation[SolventIndex][atom_count];
			}
		}
	}
	// We know ions has been chosen as old component
	else
	{
		for(k=0;k<NumberExtraComponents;k++)
		{
			atom_count = -1;
			index = (InitialPseudoAtoms) + k*(TransientMoleculeNbAtoms);
			for(i=index;i<index+TransientMoleculeNbAtoms;i++)
			{
				atom_count++;
				PseudoAtoms[i].Charge1 = VectorChargeAlchemicalTransformation[SaltIndex[CurrentAlchemicalReaction][k]][atom_count];
			}
		}
	}
	
	// Initialize component charge.
	for(i=0;i<NumberExtraComponents;i++)
	{
	  for(j=0;j<Components[IndexExtraComponent[i]].NumberOfAtoms;j++)
	  {
		Components[IndexExtraComponent[i]].Charge[j]=PseudoAtoms[Components[IndexExtraComponent[i]].Type[j]].Charge1;
	  }
	}
	/*
	// Initialize the charge of adsorbate.
	int IndexExtraTransientAdsorbate = NumberOfAdsorbateMolecules[CurrentSystem]-NumberTransientMoities[CurrentAlchemicalReaction];
	
	for(i=0;i<NumberExtraComponents;i++)
	{	  
	  for(j=0;j<Components[IndexExtraComponent[i]].NumberOfMolecules[CurrentSystem];j++)
	  {
		 for(k=0;k<Components[IndexExtraComponent[i]].NumberOfAtoms;k++)
		 {
			 Adsorbates[CurrentSystem][IndexExtraTransientAdsorbate].Atoms[k].Charge=Components[IndexExtraComponent[i]].Charge[k];
		 }
		 IndexExtraTransientAdsorbate++;
	  }
	}
	*/
}

/*********************************************************************************************************
 * Name       | UpdateChargeInterpolationAlchemicalTransformation    (Added by A. de Izarra)             *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Update the charges of the pseudo atoms of the molecules that undergo the alchemical		 *
 * 			  | transformation, during the alchemical transformation.									 *																			           					 
 * Parameters | int	NumberOldComponent => determines if old component is water (NumberOldComponent=1)    *
 * 			  |	to be transformed into ions or ions (NumberOldComponent=2) to be transformed into water. *	
 * 			  | REAL Lambda => Interpolation parameter (between 0 and 1) of the alchemical transformation*											
 *********************************************************************************************************/
void UpdateChargeInterpolationAlchemicalTransformation(int NumberOldComponent, REAL Lambda)
{
	int k,i,j;

	// Index to put the proper parameter in the right place in PotentialParms
	int index;
	int atom_count;
	
	// Update the charge during the alchemical change.
	
	// We know water has been chosen as old component
	if(NumberOldComponent == 1)
	{
		for(k=0;k<NumberExtraComponents;k++)
		{
			atom_count = -1;
			index = (InitialPseudoAtoms) + k*(TransientMoleculeNbAtoms);
			for(i=index;i<index+TransientMoleculeNbAtoms;i++)
			{
				atom_count++;						   																			 
				PseudoAtoms[i].Charge1 = (1.0-Lambda)*(VectorChargeAlchemicalTransformation[SolventIndex][atom_count]) + Lambda*(VectorChargeAlchemicalTransformation[SaltIndex[CurrentAlchemicalReaction][k]][atom_count]);				
			}
		}
	}
	// We know ions has been chosen as old component
	else
	{
		for(k=0;k<NumberExtraComponents;k++)
		{
			atom_count = -1;
			index = (InitialPseudoAtoms) + k*(TransientMoleculeNbAtoms);
			
			for(i=index;i<index+TransientMoleculeNbAtoms;i++)
			{
				atom_count++;			 
				PseudoAtoms[i].Charge1 = (1.0-Lambda)*(VectorChargeAlchemicalTransformation[SaltIndex[CurrentAlchemicalReaction][k]][atom_count]) + Lambda*(VectorChargeAlchemicalTransformation[SolventIndex][atom_count]);
			}
		}
	} 
	
	// Update the charge of the component.
	for(i=0;i<NumberExtraComponents;i++)
	{
	  for(j=0;j<Components[IndexExtraComponent[i]].NumberOfAtoms;j++)
	  {
		Components[IndexExtraComponent[i]].Charge[j]=PseudoAtoms[Components[IndexExtraComponent[i]].Type[j]].Charge1;
	  }
	}	
	
	// Update the charge of adsorbate.
	int IndexExtraTransientAdsorbate = NumberOfAdsorbateMolecules[CurrentSystem]-NumberTransientMoities[CurrentAlchemicalReaction];
	
	for(i=0;i<NumberExtraComponents;i++)
	{	  
	  for(j=0;j<Components[IndexExtraComponent[i]].NumberOfMolecules[CurrentSystem];j++)
	  {
		 for(k=0;k<Components[IndexExtraComponent[i]].NumberOfAtoms;k++)
		 {
			 Adsorbates[CurrentSystem][IndexExtraTransientAdsorbate].Atoms[k].Charge=Components[IndexExtraComponent[i]].Charge[k];
		 }
		 IndexExtraTransientAdsorbate++;
	  }
	}	
	
}

/*********************************************************************************************************
 * Name       | InitializeVectorforMixingRule    (Added by A. de Izarra)            					 *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Initialize the mixing rules for the pseudo atoms of the molecules that undergo 			 *
 * 		      | the alchemical transformation.															 *																			           					 
 * Parameters | No parameters																			 *											
 *********************************************************************************************************/
void InitializeVectorforMixingRule(void)
{
	
	 // Step 1: preparation of vectors.
	 int Species,i;
	 
	 int type;
											
	 VectorforMixingRuleAlchemicalTransformation = (REAL(**)[2])calloc(3,sizeof(REAL(*)[2]));

	 for(Species=0; Species<3; Species++)
	 {
		 VectorforMixingRuleAlchemicalTransformation[Species] = (REAL(*)[2])calloc(TransientMoleculeNbAtoms,sizeof(REAL[2]));
		 
		 for(i=0;i<TransientMoleculeNbAtoms;i++)
		 {
			 // Manage the solvent.
			 if(Species == SolventIndex)
			 {
				 type = Components[SolventIndex].Type[i];
				  
				 VectorforMixingRuleAlchemicalTransformation[Species][i][0]=  PotentialParms[type][type][0];
				 VectorforMixingRuleAlchemicalTransformation[Species][i][1]=  PotentialParms[type][type][1];		 
			 }
			 // Manage the ions.
			 else
			 {
				 // for i==0, it is the first atom corresponding to the atomic ions.
				 if(i==0)
				 {
					type = Components[Species].Type[i];
					
					VectorforMixingRuleAlchemicalTransformation[Species][i][0]=  PotentialParms[type][type][0];
					VectorforMixingRuleAlchemicalTransformation[Species][i][1]=  PotentialParms[type][type][1];
				 }
				 else
				 {
					VectorforMixingRuleAlchemicalTransformation[Species][i][0]=0.0;
					VectorforMixingRuleAlchemicalTransformation[Species][i][1]=0.0;
				 }
			 }
		 }
	 } 	 
}

/*********************************************************************************************************
 * Name       | InitializeVDWTransientMoities    (Added by A. de Izarra)                              	 *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Initialize vdW radius of the pseudo atoms of the molecules that undergo the alchemical	 *
 * 			  | transformation.										     		 						 *
 * Parameters | int	NumberOldComponent => determines if old component is water (NumberOldComponent=1)    *
 * 			  |	to be transformed into ions or ions (NumberOldComponent=2) to be transformed into water. *
 *********************************************************************************************************/
void InitializeVDWTransientMoities(int NumberOldComponent)
{	
	int k,i,j;
	

	// Index to put the proper parameter in the right place in PotentialParms
	int index;
	int atom_count;
	
	// Prepare the interaction between the transformed molecules (aka: PotentialParms[i][i] for alchemical transformation)
	
	// We know water has been chosen as old component
	if(NumberOldComponent == 1)
	{
		for(k=0;k<NumberExtraComponents;k++)
		{
			atom_count = -1;
			index = (InitialPseudoAtoms) + k*(TransientMoleculeNbAtoms);
			for(i=index;i<index+TransientMoleculeNbAtoms;i++)
			{
				atom_count++;
				
				PotentialParms[i][i][0] = VectorforMixingRuleAlchemicalTransformation[SolventIndex][atom_count][0];
				PotentialParms[i][i][1] = VectorforMixingRuleAlchemicalTransformation[SolventIndex][atom_count][1];
			}
		}
	}
	// We know ions has been chosen as old component
	else
	{
		for(k=0;k<NumberExtraComponents;k++)
		{
			atom_count = -1;
			index = (InitialPseudoAtoms) + k*(TransientMoleculeNbAtoms);
			
			for(i=index;i<index+TransientMoleculeNbAtoms;i++)
			{
				atom_count++;
				PotentialParms[i][i][0] = VectorforMixingRuleAlchemicalTransformation[SaltIndex[CurrentAlchemicalReaction][k]][atom_count][0];
				PotentialParms[i][i][1] = VectorforMixingRuleAlchemicalTransformation[SaltIndex[CurrentAlchemicalReaction][k]][atom_count][1];
			}
		}
	}
    // Compute the mixing rule with interpolated non-bonded vdw parameters.
    
    if(GeneralMixingRule=LORENTZ_BERTHELOT)
    {
		for(i=0;i<InitialPseudoAtoms+2*TransientMoleculeNbAtoms;i++)
		{
		  for(j=i+1;j<InitialPseudoAtoms+2*TransientMoleculeNbAtoms;j++)
		  {
			if((i<InitialPseudoAtoms)&&(j<InitialPseudoAtoms))
			{
				continue;
			}
			if((PotentialType[i][i]==LENNARD_JONES)&&(PotentialType[j][j]==LENNARD_JONES))
			{	
			  PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
			  PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
			  PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=(p_1^i+p_1^j)/2
			  PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=(p_1^j+p_1^i)/2
			}
		  }
		}
      }
	  else if(GeneralMixingRule=JORGENSEN)
	  {
		for(i=0;i<InitialPseudoAtoms+2*TransientMoleculeNbAtoms;i++)
		  for(j=i+1;j<InitialPseudoAtoms+2*TransientMoleculeNbAtoms;j++)
		  {
			if((i<InitialPseudoAtoms)&&(j<InitialPseudoAtoms))
			{
				continue;
			}  
			  
			if((PotentialType[i][i]==LENNARD_JONES)&&(PotentialType[j][j]==LENNARD_JONES))
			{
			  PotentialType[i][j]=LENNARD_JONES;
			  PotentialType[j][i]=LENNARD_JONES;
			  PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
			  PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
			  PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
			  PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
			}
		  }
	  }	  
}

/*********************************************************************************************************
 * Name       | UpdateMixingRuleVDWInterpolationAlchemicalTransformation    (Added by A. de Izarra)      *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Update the mixing rules and the vdw parameters during the alchemical transformation      *																		           					 
 * Parameters | int	NumberOldComponent => determines if old component is water (NumberOldComponent=1)    *
 * 			  |	to be transformed into ions or ions (NumberOldComponent=2) to be transformed into water. *	
 * 			  | REAL Lambda => Interpolation parameter (between 0 and 1) of the alchemical transformation*											
 *********************************************************************************************************/
void UpdateMixingRuleVDWInterpolationAlchemicalTransformation(int NumberOldComponent, REAL Lambda)
{
	int k,i,j;

	// Index to put the proper parameter in the right place in PotentialParms
	int index;
	int atom_count;
	
	// Prepare the interaction between the transformed molecules (aka: PotentialParms[i][i] for alchemical transformation)
	
	// We know water has been chosen as old component
	if(NumberOldComponent == 1)
	{
		for(k=0;k<NumberExtraComponents;k++)
		{
			atom_count = -1;
			index = (InitialPseudoAtoms) + k*(TransientMoleculeNbAtoms);
			for(i=index;i<index+TransientMoleculeNbAtoms;i++)
			{
				atom_count++;
				PotentialParms[i][i][0] = (1.0-Lambda)*(VectorforMixingRuleAlchemicalTransformation[SolventIndex][atom_count][0]) + Lambda*(VectorforMixingRuleAlchemicalTransformation[SaltIndex[CurrentAlchemicalReaction][k]][atom_count][0]);
				PotentialParms[i][i][1] = (1.0-Lambda)*(VectorforMixingRuleAlchemicalTransformation[SolventIndex][atom_count][1]) + Lambda*(VectorforMixingRuleAlchemicalTransformation[SaltIndex[CurrentAlchemicalReaction][k]][atom_count][1]);
			}
		}
	}
	// We know ions has been chosen as old component
	else
	{
		for(k=0;k<NumberExtraComponents;k++)
		{
			atom_count = -1;
			index = (InitialPseudoAtoms) + k*(TransientMoleculeNbAtoms);
			
			for(i=index;i<index+TransientMoleculeNbAtoms;i++)
			{
				atom_count++;
				PotentialParms[i][i][0] = (1.0-Lambda)*(VectorforMixingRuleAlchemicalTransformation[SaltIndex[CurrentAlchemicalReaction][k]][atom_count][0]) + Lambda*(VectorforMixingRuleAlchemicalTransformation[SolventIndex][atom_count][0]);
				PotentialParms[i][i][1] = (1.0-Lambda)*(VectorforMixingRuleAlchemicalTransformation[SaltIndex[CurrentAlchemicalReaction][k]][atom_count][1]) + Lambda*(VectorforMixingRuleAlchemicalTransformation[SolventIndex][atom_count][1]);
			}
		}
	}
    // Compute the mixing rule with interpolated non-bonded vdw parameters.
    if(GeneralMixingRule=LORENTZ_BERTHELOT)
    {
		for(i=0;i<InitialPseudoAtoms+2*TransientMoleculeNbAtoms;i++)
		{
		  for(j=i+1;j<InitialPseudoAtoms+2*TransientMoleculeNbAtoms;j++)
		  {
			if((i<InitialPseudoAtoms)&&(j<InitialPseudoAtoms))
			{
				continue;
			}
			if((PotentialType[i][i]==LENNARD_JONES)&&(PotentialType[j][j]==LENNARD_JONES))
			{	
			  PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
			  PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
			  PotentialParms[i][j][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ij=(p_1^i+p_1^j)/2
			  PotentialParms[j][i][1]=0.5*(PotentialParms[i][i][1]+PotentialParms[j][j][1]); // p_1^ji=(p_1^j+p_1^i)/2
			}
		  }
		}
      }
	  else if(GeneralMixingRule=JORGENSEN)
	  {
		for(i=0;i<InitialPseudoAtoms+2*TransientMoleculeNbAtoms;i++)
		  for(j=i+1;j<InitialPseudoAtoms+2*TransientMoleculeNbAtoms;j++)
		  {
			if((i<InitialPseudoAtoms)&&(j<InitialPseudoAtoms))
			{
				continue;
			}  
			  
			if((PotentialType[i][i]==LENNARD_JONES)&&(PotentialType[j][j]==LENNARD_JONES))
			{
			  PotentialType[i][j]=LENNARD_JONES;
			  PotentialType[j][i]=LENNARD_JONES;
			  PotentialParms[i][j][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ij=sqrt(p_0^i*p_0^j)
			  PotentialParms[j][i][0]=sqrt(PotentialParms[i][i][0]*PotentialParms[j][j][0]); // p_0^ji=sqrt(p_0^j*p_0^i)
			  PotentialParms[i][j][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ij=sqrt(p_1^i*p_1^j)
			  PotentialParms[j][i][1]=sqrt(PotentialParms[i][i][1]*PotentialParms[j][j][1]); // p_1^ji=sqrt(p_1^j*p_1^i)
			}
		  }
	  }
}

/*********************************************************************************************************
 * Name       | StoreChosenMoitiesCoordinates    (Added by A. de Izarra)      							 *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Store the coordinates of the chosen molecules that will be alchemically tranformed.      *																		           					 
 * Parameters | int	NumberOldComponent => determines if old component is water (NumberOldComponent=1)    *
 * 			  |	to be transformed into ions or ions (NumberOldComponent=2) to be transformed into water. *									
 *********************************************************************************************************/
void StoreChosenMoitiesCoordinates(int NumberOldComponent)
{	
	int i,j;
	int Type;
	int Index;

	// Copy coordinates and properties of chosen adsorbates.

	AdsorbatesReferenceChosen=(ADSORBATE_MOLECULE*)calloc(NumberTransientMoities[CurrentAlchemicalReaction],sizeof(ADSORBATE_MOLECULE));
	
	for(i=0;i<NumberTransientMoities[CurrentAlchemicalReaction];i++)
	{
	  Index=ChosenMoleculeAlchemicalTransformation[CurrentAlchemicalReaction][i];
		
	  Type = Adsorbates[CurrentSystem][Index].Type;
	  AdsorbatesReferenceChosen[i].Type=Type;

	  AdsorbatesReferenceChosen[i].Atoms=(ATOM*)calloc(Components[Type].NumberOfAtoms,sizeof(ATOM));
	  AdsorbatesReferenceChosen[i].NumberOfAtoms=Components[Type].NumberOfAtoms;

	  if(Components[Type].NumberOfGroups>0)
		AdsorbatesReferenceChosen[i].Groups=(GROUP*)calloc(Components[Type].NumberOfGroups,sizeof(GROUP));
	    	
	  // Copy the pointer by value of positiona nd groups in the new vector.
	  for(j=0;j<Components[Type].NumberOfAtoms;j++)
	  { 
		  AdsorbatesReferenceChosen[i].Atoms[j] = Adsorbates[CurrentSystem][Index].Atoms[j];
	  }

	  for(j=0;j<Components[Type].NumberOfGroups;j++)
	  {
		  AdsorbatesReferenceChosen[i].Groups[j] = Adsorbates[CurrentSystem][Index].Groups[j];
	  }
	}
}

/*********************************************************************************************************
 * Name       | DeleteChosenMoities    (Added by A. de Izarra)      							 		 *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Delete the initial molecules to be alchemically transformed. These molecules are 		 *
 *            | redefined as "transient molecules" in MakeInitialTransient function, because they are 	 *
 * 			  | candidates for the alchemical transformation.											 *						           					 
 * Parameters | No parameters																			 *									
 *********************************************************************************************************/
void DeleteChosenMoities(void)
{
  int i,j;

  REAL DeltaU,UTailOld,RosenbluthOld;
  int LastMolecule;
  int ComponentLastMolecule;
 
  for(i=0;i<NumberTransientMoities[CurrentAlchemicalReaction];i++)
  {
		// The index of the chosen adsorbate molecule.
		CurrentAdsorbateMolecule = ChosenMoleculeAlchemicalTransformation[CurrentAlchemicalReaction][i];
		
		// The index of the chosen component.
		CurrentComponent=Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Type;

		for(j=0;j<Components[CurrentComponent].NumberOfAtoms;j++)
		{
		CFVDWScaling[j]=1.0;
		CFChargeScaling[j]=1.0;
		}
		
		// Retrace the energy of the chosen componenent
		NumberOfBeadsAlreadyPlaced=0;
		
		for(j=0;j<Components[CurrentComponent].NumberOfAtoms;j++)
		{
			OldPosition[j]=Adsorbates[CurrentSystem][CurrentAdsorbateMolecule].Atoms[j].Position;
		}

		NumberOfTrialPositions=NumberOfTrialPositionsSwap;
		NumberOfTrialPositionsForTheFirstBead=NumberOfTrialPositionsForTheFirstBeadSwap;

		RosenbluthOld=RetraceMolecule(CBMC_DELETION);

		UTailOld=TailMolecularEnergyDifferenceRemove();

		if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
			CalculateEwaldFourierAdsorbate(FALSE,TRUE,CurrentAdsorbateMolecule,0);
		

		if(ComputePolarization)
			ComputeNewPolarizationEnergy(FALSE,CurrentAdsorbateMolecule,-1);
		
		UAdsorbateBond[CurrentSystem]-=UBondOld[CurrentSystem];
		UAdsorbateUreyBradley[CurrentSystem]-=UUreyBradleyOld[CurrentSystem];
		UAdsorbateBend[CurrentSystem]-=UBendOld[CurrentSystem];
		UAdsorbateBendBend[CurrentSystem]-=UBendBendOld[CurrentSystem];
		UAdsorbateInversionBend[CurrentSystem]-=UInversionBendOld[CurrentSystem];
		UAdsorbateTorsion[CurrentSystem]-=UTorsionOld[CurrentSystem];
		UAdsorbateImproperTorsion[CurrentSystem]-=UImproperTorsionOld[CurrentSystem];
		UAdsorbateBondBond[CurrentSystem]-=UBondBondOld[CurrentSystem];
		UAdsorbateBondBend[CurrentSystem]-=UBondBendOld[CurrentSystem];
		UAdsorbateBondTorsion[CurrentSystem]-=UBondTorsionOld[CurrentSystem];
		UAdsorbateBendTorsion[CurrentSystem]-=UBendTorsionOld[CurrentSystem];
		UAdsorbateIntraVDW[CurrentSystem]-=UIntraVDWOld[CurrentSystem];
	    
		UAdsorbateAdsorbate[CurrentSystem]-=UAdsorbateVDWOld[CurrentSystem];
		UAdsorbateAdsorbateVDW[CurrentSystem]-=UAdsorbateVDWOld[CurrentSystem];
		UAdsorbateCation[CurrentSystem]-=UCationVDWOld[CurrentSystem];
		UAdsorbateCationVDW[CurrentSystem]-=UCationVDWOld[CurrentSystem];
		UHostAdsorbate[CurrentSystem]-=UHostVDWOld[CurrentSystem];
		UHostAdsorbateVDW[CurrentSystem]-=UHostVDWOld[CurrentSystem];
		
		
		UTailCorrection[CurrentSystem]-=UTailOld;
	
		UHostPolarization[CurrentSystem]=UHostPolarizationNew[CurrentSystem];
	
		UAdsorbatePolarization[CurrentSystem]=UAdsorbatePolarizationNew[CurrentSystem]+UPolarizationNew[CurrentSystem];

		UCationPolarization[CurrentSystem]=UCationPolarizationNew[CurrentSystem];
		
		UHostBackPolarization[CurrentSystem]=UHostBackPolarizationNew[CurrentSystem];
		UAdsorbateBackPolarization[CurrentSystem]=UAdsorbateBackPolarizationNew[CurrentSystem]+UBackPolarizationNew[CurrentSystem];
		UCationBackPolarization[CurrentSystem]=UCationBackPolarizationNew[CurrentSystem];

		if(ChargeMethod!=NONE)
		{
		UAdsorbateIntraChargeCharge[CurrentSystem]-=UIntraChargeChargeOld[CurrentSystem];
		UAdsorbateIntraChargeBondDipole[CurrentSystem]-=UIntraChargeBondDipoleOld[CurrentSystem];
		UAdsorbateIntraBondDipoleBondDipole[CurrentSystem]-=UIntraBondDipoleBondDipoleOld[CurrentSystem];
	
		UAdsorbateAdsorbateChargeChargeReal[CurrentSystem]-=UAdsorbateChargeChargeOld[CurrentSystem];
		UAdsorbateAdsorbateChargeBondDipoleReal[CurrentSystem]-=UAdsorbateChargeBondDipoleOld[CurrentSystem];
		UAdsorbateAdsorbateBondDipoleBondDipoleReal[CurrentSystem]-=UAdsorbateBondDipoleBondDipoleOld[CurrentSystem];
		UAdsorbateAdsorbateChargeChargeFourier[CurrentSystem]+=UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem];
		UAdsorbateAdsorbateChargeBondDipoleFourier[CurrentSystem]+=UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem];
		UAdsorbateAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]+=UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem];
		UAdsorbateAdsorbateCoulomb[CurrentSystem]+=UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem]
													+UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]
													+UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]
													-UAdsorbateChargeChargeOld[CurrentSystem]
													-UAdsorbateChargeBondDipoleOld[CurrentSystem]
													-UAdsorbateBondDipoleBondDipoleOld[CurrentSystem];
		UAdsorbateAdsorbate[CurrentSystem]+=UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem]
											+UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]
											+UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]
											-UAdsorbateChargeChargeOld[CurrentSystem]
											-UAdsorbateChargeBondDipoleOld[CurrentSystem]
											-UAdsorbateBondDipoleBondDipoleOld[CurrentSystem];

		UHostAdsorbateChargeChargeReal[CurrentSystem]-=UHostChargeChargeOld[CurrentSystem];
		UHostAdsorbateChargeBondDipoleReal[CurrentSystem]-=UHostChargeBondDipoleOld[CurrentSystem];
		UHostAdsorbateBondDipoleBondDipoleReal[CurrentSystem]-=UHostBondDipoleBondDipoleOld[CurrentSystem];
		UHostAdsorbateChargeChargeFourier[CurrentSystem]+=UHostAdsorbateChargeChargeFourierDelta[CurrentSystem];
		UHostAdsorbateChargeBondDipoleFourier[CurrentSystem]+=UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem];
		UHostAdsorbateBondDipoleBondDipoleFourier[CurrentSystem]+=UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem];
		UHostAdsorbateCoulomb[CurrentSystem]+=UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]
											+UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]
											+UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]
											-UHostChargeChargeOld[CurrentSystem]
											-UHostChargeBondDipoleOld[CurrentSystem]
											-UHostBondDipoleBondDipoleOld[CurrentSystem];
		UHostAdsorbate[CurrentSystem]+=UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]
										+UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]
										+UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]
										-UHostChargeChargeOld[CurrentSystem]
										-UHostChargeBondDipoleOld[CurrentSystem]
										-UHostBondDipoleBondDipoleOld[CurrentSystem];
	
		UAdsorbateCationChargeChargeReal[CurrentSystem]-=UCationChargeChargeOld[CurrentSystem];
		UAdsorbateCationChargeBondDipoleReal[CurrentSystem]-=UCationChargeBondDipoleOld[CurrentSystem];
		UAdsorbateCationBondDipoleBondDipoleReal[CurrentSystem]-=UCationBondDipoleBondDipoleOld[CurrentSystem];
		UAdsorbateCationChargeChargeFourier[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem];
		UAdsorbateCationChargeBondDipoleFourier[CurrentSystem]+=UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem];
		UAdsorbateCationBondDipoleBondDipoleFourier[CurrentSystem]+=UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem];
		UAdsorbateCationCoulomb[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]
												+UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]
												+UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]
												-UCationChargeChargeOld[CurrentSystem]
												-UCationChargeBondDipoleOld[CurrentSystem]
												-UCationBondDipoleBondDipoleOld[CurrentSystem];
		UAdsorbateCation[CurrentSystem]+=UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]
										+UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]
										+UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]
										-UCationChargeChargeOld[CurrentSystem]
										-UCationChargeBondDipoleOld[CurrentSystem]
										-UCationBondDipoleBondDipoleOld[CurrentSystem];
	
		NetChargeAdsorbates[CurrentSystem]+=NetChargeAdsorbateDelta;
		NetChargeSystem[CurrentSystem]+=NetChargeAdsorbateDelta;
	
		if((ChargeMethod==EWALD)&&(!OmitEwaldFourier))
			AcceptEwaldAdsorbateMove(0);
		}
	
		// Update the list of chosen moities to take account of deleted index.
		LastMolecule=NumberOfAdsorbateMolecules[CurrentSystem]-1;
		ComponentLastMolecule=Adsorbates[CurrentSystem][LastMolecule].Type;
		
		for(j=i+1;j<NumberTransientMoities[CurrentAlchemicalReaction];j++)
		{
			// If a particule in the Chosen index is at the top of index, it takes the place of the current Chosenindex i
			if(ChosenMoleculeAlchemicalTransformation[CurrentAlchemicalReaction][j]==LastMolecule)
			{
				ChosenMoleculeAlchemicalTransformation[CurrentAlchemicalReaction][j]=CurrentAdsorbateMolecule;
			}
		}	
	
		RemoveAdsorbateMolecule();	
	
		DeltaU=-UBondOld[CurrentSystem]-UUreyBradleyOld[CurrentSystem]-UBendOld[CurrentSystem]-UBendBendOld[CurrentSystem]-UInversionBendOld[CurrentSystem]-UTorsionOld[CurrentSystem]
			-UImproperTorsionOld[CurrentSystem]-UBondBondOld[CurrentSystem]-UBondBendOld[CurrentSystem]-UBondTorsionOld[CurrentSystem]-UBendTorsionOld[CurrentSystem]
			-UIntraVDWOld[CurrentSystem]-UIntraChargeChargeOld[CurrentSystem]-UIntraChargeBondDipoleOld[CurrentSystem]-UIntraBondDipoleBondDipoleOld[CurrentSystem]
			-UAdsorbateVDWOld[CurrentSystem]-UCationVDWOld[CurrentSystem]-UHostVDWOld[CurrentSystem]
			-UAdsorbateChargeChargeOld[CurrentSystem]-UCationChargeChargeOld[CurrentSystem]-UHostChargeChargeOld[CurrentSystem]
			-UAdsorbateChargeBondDipoleOld[CurrentSystem]-UCationChargeBondDipoleOld[CurrentSystem]-UHostChargeBondDipoleOld[CurrentSystem]
			-UAdsorbateBondDipoleBondDipoleOld[CurrentSystem]-UCationBondDipoleBondDipoleOld[CurrentSystem]-UHostBondDipoleBondDipoleOld[CurrentSystem]+
				UHostAdsorbateChargeChargeFourierDelta[CurrentSystem]+UAdsorbateAdsorbateChargeChargeFourierDelta[CurrentSystem]+UAdsorbateCationChargeChargeFourierDelta[CurrentSystem]+
				UHostAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateAdsorbateChargeBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationChargeBondDipoleFourierDelta[CurrentSystem]+
				UHostAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateAdsorbateBondDipoleBondDipoleFourierDelta[CurrentSystem]+UAdsorbateCationBondDipoleBondDipoleFourierDelta[CurrentSystem]+
				UDeltaPolarization-UTailOld;
	
		UTotal[CurrentSystem]+=DeltaU;
  }	
}

/*********************************************************************************************************
 * Name       | SwitchMoietiestoRegularComponents    (Added by A. de Izarra)      					 	 *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | If the alchemical transformation is accepted, the transient molecules totally transformed*
 * 			  | become normal molecules.																 *
 * Parameters | int	NumberNewComponent => determines if new component is water (NumberNewComponent=1)    *
 * 			  |	to be transformed into ions or ions (NumberNewComponent=2) to be transformed into water. *	
 *********************************************************************************************************/
void SwitchMoietiestoRegularComponents(int NumberNewComponent)
{		
	int Comp;
	int i,j;
	int m;
	int count = -1;
	int type;
    int nr_atoms;
	int UntouchedAdsorbates=NumberOfAdsorbateMolecules[CurrentSystem]-NumberTransientMoities[CurrentAlchemicalReaction];
  
    // Copy the transient component in new adsorbates.
	ADSORBATE_MOLECULE*AdsorbatesTransientMoities=(ADSORBATE_MOLECULE*)calloc(NumberTransientMoities[CurrentAlchemicalReaction],sizeof(ADSORBATE_MOLECULE));
	
	for(Comp=0; Comp<NumberExtraComponents; Comp++)
	{
		for(i=0;i<MultiplicitySalt[CurrentAlchemicalReaction][Comp];i++)
	    {  
			count++;
			m = count+UntouchedAdsorbates;

			if(NumberNewComponent == 1)
			{
				CurrentComponent=SolventIndex;
			}
			else
			{
				CurrentComponent=SaltIndex[CurrentAlchemicalReaction][Comp];
			}
			
			AdsorbatesTransientMoities[count].Type = CurrentComponent;
			Components[CurrentComponent].NumberOfMolecules[CurrentSystem]++;
			AdsorbatesTransientMoities[count].NumberOfAtoms=Components[CurrentComponent].NumberOfAtoms;
			AdsorbatesTransientMoities[count].Atoms=(ATOM*)calloc(Components[CurrentComponent].NumberOfAtoms,sizeof(ATOM));
			
			for(j=0;j<Components[CurrentComponent].NumberOfAtoms;j++)
			{	
				type=Components[CurrentComponent].Type[j];		
				AdsorbatesTransientMoities[count].Atoms[j].Position=Adsorbates[CurrentSystem][m].Atoms[j].Position;
				AdsorbatesTransientMoities[count].Atoms[j].AnisotropicPosition=Adsorbates[CurrentSystem][m].Atoms[j].AnisotropicPosition;
				AdsorbatesTransientMoities[count].Atoms[j].Velocity=Adsorbates[CurrentSystem][m].Atoms[j].Velocity;
				AdsorbatesTransientMoities[count].Atoms[j].Force=Adsorbates[CurrentSystem][m].Atoms[j].Force;

				// new Continuous-Fraction scaling factors are taken from the component-information
				AdsorbatesTransientMoities[count].Atoms[j].CFVDWScalingParameter=CFVDWScaling[j];
				AdsorbatesTransientMoities[count].Atoms[j].CFChargeScalingParameter=CFChargeScaling[j];

				AdsorbatesTransientMoities[count].Atoms[j].Type=type;
				AdsorbatesTransientMoities[count].Atoms[j].Fixed.x=Components[CurrentComponent].Fixed[j];
				AdsorbatesTransientMoities[count].Atoms[j].Fixed.y=Components[CurrentComponent].Fixed[j];
				AdsorbatesTransientMoities[count].Atoms[j].Fixed.z=Components[CurrentComponent].Fixed[j];
				AdsorbatesTransientMoities[count].Atoms[j].Charge=Adsorbates[CurrentSystem][m].Atoms[j].Charge;
				
			}
			
			if(Components[CurrentComponent].NumberOfGroups>0)
				AdsorbatesTransientMoities[count].Groups=(GROUP*)calloc(Components[CurrentComponent].NumberOfGroups,sizeof(GROUP));
			
			nr_atoms=Components[CurrentComponent].NumberOfAtoms;
			NumberOfAtomsPerSystem[CurrentSystem]+=nr_atoms;
			NumberOfChargesPerSystem[CurrentSystem]+=Components[CurrentComponent].NumberOfCharges;
			NumberOfBondDipolesPerSystem[CurrentSystem]+=Components[CurrentComponent].NumberOfBondDipoles;
			
			// Reinitialize adsorbate pointers.
			Adsorbates[CurrentSystem][m].Atoms=NULL;
			Adsorbates[CurrentSystem][m].Atoms=(ATOM*)calloc(Components[CurrentComponent].NumberOfAtoms,sizeof(ATOM));
			
			Adsorbates[CurrentSystem][m].Groups=NULL;
			Adsorbates[CurrentSystem][m].Groups=(GROUP*)calloc(Components[CurrentComponent].NumberOfGroups,sizeof(GROUP));
			
			Adsorbates[CurrentSystem][m].NumberOfAtoms=AdsorbatesTransientMoities[count].NumberOfAtoms;
			Adsorbates[CurrentSystem][m].Type=AdsorbatesTransientMoities[count].Type;
		
			// Update new adsorbate.
			for(j=0;j<Components[CurrentComponent].NumberOfAtoms;j++)
			{
				type=Components[CurrentComponent].Type[j];
				Adsorbates[CurrentSystem][m].Atoms[j].Position=AdsorbatesTransientMoities[count].Atoms[j].Position;
				Adsorbates[CurrentSystem][m].Atoms[j].AnisotropicPosition=AdsorbatesTransientMoities[count].Atoms[j].AnisotropicPosition;
				Adsorbates[CurrentSystem][m].Atoms[j].Velocity=AdsorbatesTransientMoities[count].Atoms[j].Velocity;
				Adsorbates[CurrentSystem][m].Atoms[j].Force=AdsorbatesTransientMoities[count].Atoms[j].Force;

				// new Continuous-Fraction scaling factors are taken from the component-information
				Adsorbates[CurrentSystem][m].Atoms[j].CFVDWScalingParameter=CFVDWScaling[j];
				Adsorbates[CurrentSystem][m].Atoms[j].CFChargeScalingParameter=CFChargeScaling[j];
	
				Adsorbates[CurrentSystem][m].Atoms[j].Type=AdsorbatesTransientMoities[count].Atoms[j].Type;
				Adsorbates[CurrentSystem][m].Atoms[j].Fixed.x=Components[CurrentComponent].Fixed[j];
				Adsorbates[CurrentSystem][m].Atoms[j].Fixed.y=Components[CurrentComponent].Fixed[j];
				Adsorbates[CurrentSystem][m].Atoms[j].Fixed.z=Components[CurrentComponent].Fixed[j];
				Adsorbates[CurrentSystem][m].Atoms[j].Charge=AdsorbatesTransientMoities[count].Atoms[j].Charge;
				NumberOfPseudoAtomsType[CurrentSystem][type]++;
			}
			
			// modify the degrees of freedom
			 for(j=0;j<Components[CurrentComponent].NumberOfGroups;j++)
			  {
				// Add new degrees of freedom accounting from transient->moieties
				if(Components[CurrentComponent].Groups[j].Rigid)
				{
				  DegreesOfFreedomAdsorbates[CurrentSystem]+=3;
				  DegreesOfFreedomTranslation[CurrentSystem]+=3;
				  DegreesOfFreedomTranslationalAdsorbates[CurrentSystem]+=3;
				  DegreesOfFreedom[CurrentSystem]+=3;
														   
				  DegreesOfFreedomRotation[CurrentSystem]+=Components[CurrentComponent].Groups[j].RotationalDegreesOfFreedom;
				  DegreesOfFreedomAdsorbates[CurrentSystem]+=Components[CurrentComponent].Groups[j].RotationalDegreesOfFreedom;
				  DegreesOfFreedomRotationalAdsorbates[CurrentSystem]+=Components[CurrentComponent].Groups[j].RotationalDegreesOfFreedom;
				  DegreesOfFreedom[CurrentSystem]+=Components[CurrentComponent].Groups[j].RotationalDegreesOfFreedom;
				}
				else
				{
				  DegreesOfFreedomTranslation[CurrentSystem]+=3*Components[CurrentComponent].Groups[j].NumberOfGroupAtoms;
				  DegreesOfFreedomAdsorbates[CurrentSystem]+=3*Components[CurrentComponent].Groups[j].NumberOfGroupAtoms;
				  DegreesOfFreedomTranslationalAdsorbates[CurrentSystem]+=3*Components[CurrentComponent].Groups[j].NumberOfGroupAtoms;
				  DegreesOfFreedom[CurrentSystem]+=3*Components[CurrentComponent].Groups[j].NumberOfGroupAtoms;
				}
				
				// remove degrees of freedom accounting from former transient
				if(Components[IndexExtraComponent[Comp]].Groups[j].Rigid)
				{
				  DegreesOfFreedomAdsorbates[CurrentSystem]-=3;
				  DegreesOfFreedomTranslation[CurrentSystem]-=3;
				  DegreesOfFreedomTranslationalAdsorbates[CurrentSystem]-=3;
				  DegreesOfFreedom[CurrentSystem]-=3;
														   
				  DegreesOfFreedomRotation[CurrentSystem]-=Components[IndexExtraComponent[Comp]].Groups[j].RotationalDegreesOfFreedom;
				  DegreesOfFreedomAdsorbates[CurrentSystem]-=Components[IndexExtraComponent[Comp]].Groups[j].RotationalDegreesOfFreedom;
				  DegreesOfFreedomRotationalAdsorbates[CurrentSystem]-=Components[IndexExtraComponent[Comp]].Groups[j].RotationalDegreesOfFreedom;
				  DegreesOfFreedom[CurrentSystem]-=Components[IndexExtraComponent[Comp]].Groups[j].RotationalDegreesOfFreedom;
				}
				else
				{
				  DegreesOfFreedomTranslation[CurrentSystem]-=3*Components[IndexExtraComponent[Comp]].Groups[j].NumberOfGroupAtoms;
				  DegreesOfFreedomAdsorbates[CurrentSystem]-=3*Components[IndexExtraComponent[Comp]].Groups[j].NumberOfGroupAtoms;
				  DegreesOfFreedomTranslationalAdsorbates[CurrentSystem]-=3*Components[IndexExtraComponent[Comp]].Groups[j].NumberOfGroupAtoms;
				  DegreesOfFreedom[CurrentSystem]-=3*Components[IndexExtraComponent[Comp]].Groups[j].NumberOfGroupAtoms;
				}
				
			  }	
			
			// update the center of mass
			UpdateGroupCenterOfMassAdsorbate(m);

			// compute the quaternion (orientation) from the positions
			ComputeQuaternionAdsorbate(m);
							
			// Update number of atoms and so one (remove atoms and moities of transient component).				
			nr_atoms=Components[IndexExtraComponent[Comp]].NumberOfAtoms;
			NumberOfAtomsPerSystem[CurrentSystem]-=nr_atoms;
			NumberOfChargesPerSystem[CurrentSystem]-=Components[CurrentComponent].NumberOfCharges;
			NumberOfBondDipolesPerSystem[CurrentSystem]-=Components[CurrentComponent].NumberOfBondDipoles;
			
			Components[IndexExtraComponent[Comp]].NumberOfMolecules[CurrentSystem]--;
			
			for(j=0;j<Components[IndexExtraComponent[Comp]].NumberOfAtoms;j++)
			{
				type=Components[IndexExtraComponent[Comp]].Type[j];
				NumberOfPseudoAtomsType[CurrentSystem][type]--;
			}
		}
				
	}
	 
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
	 
	 
	 if((Ensemble[CurrentSystem]==MuPT)||(Ensemble[CurrentSystem]==MuPTPR)||(Ensemble[CurrentSystem]==MuVT))
	 {
		InitializeNoseHooverCurrentSystem();
	 }
 
	 for(i=0;i<NumberTransientMoities[0];i++)
	 {
			AdsorbatesTransientMoities[i].Atoms = NULL;
			AdsorbatesTransientMoities[i].Groups = NULL;
	 }
	
	 free(AdsorbatesTransientMoities);
}

/*********************************************************************************************************
 * Name       | InitializeMassTransientMoities    (Added by A. de Izarra)      					 	 	 *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Initialize the mass of transient moities, to be the mass of water molecule				 *
 * Parameters | int	NumberOldComponent => determines if old component is water (NumberOldComponent=1)    *
 * 			  |	to be transformed into ions or ions (NumberOldComponent=2) to be transformed into water. *	
 *********************************************************************************************************/
void InitializeMassTransientMoities(int NumberOldComponent)
{	
	// The mass of the transient molecule has the mass of water
	
	// for water->ions:  keep water mass during alchemical process and switch to mass of ions at the end.
	// for ions->water:  switch at the beginning of the process the mass of ions to water.
	
	int i,j,k;
	REAL mass,total_mass;
	total_mass=0.0;
	int index;
	int atom_count;
	int A;
	
	// Set up the mass of the pseudo atom. (corresponding to mass of water)
	for(i=0;i<NumberExtraComponents;i++)
	{
		atom_count = -1;
		index = (InitialPseudoAtoms) + i*(TransientMoleculeNbAtoms);
		for(j=index;j<index+TransientMoleculeNbAtoms;j++)
		{
				atom_count++;
				PseudoAtoms[j].Mass=PseudoAtoms[Components[SolventIndex].Type[atom_count]].Mass;
		}
	}
	// Set up the mass of the transient component.
	for(i=0;i<NumberExtraComponents;i++)
	{
		CurrentComponent = IndexExtraComponent[i];
			
		// In order to go from water -> ions, the transient components have the mass of ions.
		Components[CurrentComponent].Mass=Components[SolventIndex].Mass;
						
		for(j=0;j<Components[CurrentComponent].NumberOfGroups;j++)
		{
			Components[CurrentComponent].Groups[j].Mass=0.0;
			for(k=0;k<Components[CurrentComponent].Groups[j].NumberOfGroupAtoms;k++)
			{
				A=Components[CurrentComponent].Groups[j].Atoms[k];
				Components[CurrentComponent].Groups[j].Mass+=PseudoAtoms[Components[CurrentComponent].Type[A]].Mass;
			}
		}
	}
}

/*********************************************************************************************************
 * Name       | EndMassTransientMoities    (Added by A. de Izarra)      					 	 	 	 *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | If the final transient moities are ions -> update the mass to ions.				 		 *
 * 			  | If the final transient moities are water -> keep the mass of water.						 *	
 * Parameters | int	NumberOldComponent => determines if old component is water (NumberOldComponent=1)    *
 * 			  |	to be transformed into ions or ions (NumberOldComponent=2) to be transformed into water. *	
 *********************************************************************************************************/
void EndMassTransientMoities(int NumberOldComponent)
{	
	// The mass of the transient molecule are the mass of the new component.
	int i,j,k;
	REAL mass,total_mass;
	total_mass=0.0;
	int index;
	int atom_count;
	int A;
	
	if(NumberOldComponent == 1)
	{
		// Set up the mass of the pseudo atom. (corresponding to mass of ions)
		for(i=0;i<NumberExtraComponents;i++)
		{
			index = (InitialPseudoAtoms) + i*(TransientMoleculeNbAtoms);
			for(j=index;j<index+TransientMoleculeNbAtoms;j++)
			{
				
				//The mass of the first atom is the mass of the ion.
				if(j==index)
				{   
					PseudoAtoms[j].Mass=PseudoAtoms[Components[SaltIndex[CurrentAlchemicalReaction][i]].Type[0]].Mass;
				}
				else
				{
					PseudoAtoms[j].Mass=0.0;
				}
			}
		}
		// Set up the mass of the transient component.
		for(i=0;i<NumberExtraComponents;i++)
		{
			CurrentComponent = IndexExtraComponent[i];
	
			// In order to go from water -> ions, the transient components have the mass of ions.
			Components[CurrentComponent].Mass=Components[SaltIndex[CurrentAlchemicalReaction][i]].Mass;

			for(j=0;j<Components[CurrentComponent].NumberOfGroups;j++)
			{
				Components[CurrentComponent].Groups[j].Mass=0.0;
				
				for(k=0;k<Components[CurrentComponent].Groups[j].NumberOfGroupAtoms;k++)
				{
					A=Components[CurrentComponent].Groups[j].Atoms[k];
					Components[CurrentComponent].Groups[j].Mass+=PseudoAtoms[Components[CurrentComponent].Type[A]].Mass;
				}
			}
		}
	}
}

void UpdateInertiaTensorGroups(int comp)
{		
  int j,k,ill;
  REAL Mass,TotalMass,rotxyz,temp,rotall,rotlim;
  VECTOR com,pos;
  REAL_MATRIX3x3 eigenvectors;
  VECTOR eigenvalues,dr;
  int atom_nr;


  for(j=0;j<Components[comp].NumberOfGroups;j++)
  {
    for(k=0;k<Components[comp].Groups[j].NumberOfGroupAtoms;k++)
    {
		atom_nr=Components[comp].Groups[j].Atoms[k];
		Components[comp].Positions[atom_nr]=SolventBodyfixedPositions[j][atom_nr];
	}
  }
  
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

/*********************************************************************************************************
 * Name       | AllocateTransientComponentMemory    (Added by A. de Izarra)      					 	 *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Allocate memory to store information component for the molecules that undergo			 *
 * 			  | the alchemical transformation.															 *	
 * Parameters | No parameters																			 *
 *********************************************************************************************************/
void AllocateTransientComponentMemory(void)
{
   int index_salt = -1;
   int index_type;
   int count = -1;
   int comp,i,j,k,n,nr;
   int A,B,C,D,temp;
   int A1,A2,B1,B2,C1,C2;
   char index_salt_store[2];
   VECTOR dr;
	
   // Reallocate memory to store the transient molecule: two type of transient molecules corresponding to cation/water and anion/water.  
   Components=(COMPONENT*)realloc(Components,(NumberOfComponents+NumberExtraComponents)*sizeof(COMPONENT));
 
   // We use only the arguments of the structure "Components" useful for the alchemical calculation.
   for(comp=NumberOfComponents;comp<NumberOfComponents+NumberExtraComponents;comp++)
   {
	Components[comp].NumberOfAtoms=Components[SolventIndex].NumberOfAtoms;   
	   
    Components[comp].IdealGasRosenbluthWeight=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
    Components[comp].IdealGasTotalEnergy=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
    Components[comp].PartialPressure=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
    Components[comp].FugacityCoefficient=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
    Components[comp].BulkFluidDensity=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
    Components[comp].Compressibility=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
    Components[comp].MolFraction=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
    Components[comp].AmountOfExcessMolecules=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

    Components[comp].CreateNumberOfMolecules=(int*)calloc(NumberOfSystems,sizeof(int));
    Components[comp].NumberOfMolecules=(int*)calloc(NumberOfSystems,sizeof(int));

    Components[comp].FractionalMolecule=(int*)calloc(NumberOfSystems,sizeof(int));
    Components[comp].CFMoleculePresent=(int*)calloc(NumberOfSystems,sizeof(int));
    Components[comp].CFWangLandauScalingFactor=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
    Components[comp].CFBiasingFactors=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));

    Components[comp].RXMCMoleculesPresent=(int*)calloc(NumberOfSystems,sizeof(int));
    Components[comp].NumberOfRXMCMoleculesPresent=(int*)calloc(NumberOfSystems,sizeof(int));

    Components[comp].MOLEC_PER_UC_TO_MOL_PER_KG=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
    Components[comp].MOLEC_PER_UC_TO_MILLIGRAM_PER_GRAM_OF_FRAMEWORK=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
    Components[comp].MOLEC_PER_UC_TO_CC_STP_G=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
    Components[comp].MOLEC_PER_UC_TO_CC_STP_CC=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
    Components[comp].MOL_PER_KG_TO_CC_STP_G=(REAL*)calloc(NumberOfSystems,sizeof(REAL));
    Components[comp].MOL_PER_KG_TO_CC_STP_CC=(REAL*)calloc(NumberOfSystems,sizeof(REAL));

    Components[comp].BlockPockets=(int*)calloc(NumberOfSystems,sizeof(int));
    Components[comp].ComputeFreeEnergyProfile=(int*)calloc(NumberOfSystems,sizeof(int));
    Components[comp].BlockPocketsFilename=(char(*)[256])calloc(NumberOfSystems,sizeof(char[256]));
    Components[comp].NumberOfBlockCenters=(int*)calloc(NumberOfSystems,sizeof(int));
    Components[comp].BlockDistance=(REAL**)calloc(NumberOfSystems,sizeof(REAL*));
    Components[comp].BlockCenters=(VECTOR**)calloc(NumberOfSystems,sizeof(VECTOR*));
    
    Components[comp].CFLambdaHistogramSize=Components[SolventIndex].CFLambdaHistogramSize;

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

    for(i=0;i<NumberOfSystems;i++)
    {
      Components[comp].BlockPockets[i]=Components[SolventIndex].BlockPockets[i];
      
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
      
      // start with no defined fractional molecule
      Components[comp].FractionalMolecule[i]=Components[SolventIndex].FractionalMolecule[i];
      Components[comp].CFMoleculePresent[i]=Components[SolventIndex].CFMoleculePresent[i];
      Components[comp].CFWangLandauScalingFactor[i]=Components[SolventIndex].CFWangLandauScalingFactor[i];
    }
	 
    strcpy(Components[comp].MoleculeDefinition,"Transient");
    strcpy(Components[comp].Name,"");
    Components[comp].Swapable=FALSE;
    Components[comp].Widom=FALSE;
	Components[comp].ExtraFrameworkMolecule = TRUE;
	
    // initialize box restriction to false (the full box is selected)
    Components[comp].RestrictMoves=FALSE;
    
    // initialize the first box to all, the other 3 boxes to empty
    Components[comp].RestrictMovesToBox=FALSE;
    Components[comp].BoxAxisABC_Min.x=Components[comp].BoxAxisABC_Min.y=Components[comp].BoxAxisABC_Min.z=0.0;
    Components[comp].BoxAxisABC_Max.x=Components[comp].BoxAxisABC_Max.y=Components[comp].BoxAxisABC_Max.z=1.0;
    Components[comp].BoxAxisABC_Min2.x=Components[comp].BoxAxisABC_Min2.y=Components[comp].BoxAxisABC_Min2.z=1.0;
    Components[comp].BoxAxisABC_Max2.x=Components[comp].BoxAxisABC_Max2.y=Components[comp].BoxAxisABC_Max2.z=0.0;
    Components[comp].BoxAxisABC_Min3.x=Components[comp].BoxAxisABC_Min3.y=Components[comp].BoxAxisABC_Min3.z=1.0;
    Components[comp].BoxAxisABC_Max3.x=Components[comp].BoxAxisABC_Max3.y=Components[comp].BoxAxisABC_Max3.z=0.0;
    Components[comp].BoxAxisABC_Min4.x=Components[comp].BoxAxisABC_Min4.y=Components[comp].BoxAxisABC_Min4.z=1.0;
    Components[comp].BoxAxisABC_Max4.x=Components[comp].BoxAxisABC_Max4.y=Components[comp].BoxAxisABC_Max4.z=0.0;

    Components[comp].RestrictMovesToPrisms=FALSE;
    for(j=0;j<MAX_NUMBER_OF_PRISMS;j++)
    {
      Components[comp].RestrictMovesToPrism[j]  =Components[SolventIndex].RestrictMovesToPrism[j];
      Components[comp].RestrictPrismABC_Min[j].x=Components[SolventIndex].RestrictPrismABC_Min[j].x;
      Components[comp].RestrictPrismABC_Max[j].x=Components[SolventIndex].RestrictPrismABC_Max[j].x;
    }

    Components[comp].RestrictMovesToCylinders=Components[SolventIndex].RestrictMovesToCylinders;
    for(j=0;j<MAX_NUMBER_OF_CYLINDERS;j++)
    {
      Components[comp].RestrictMovesToCylinder[j]  =Components[SolventIndex].RestrictMovesToCylinder[j];
      Components[comp].RestrictCylinderABC_Min[j].x=Components[SolventIndex].RestrictCylinderABC_Min[j].x;
      Components[comp].RestrictCylinderABC_Max[j].x=Components[SolventIndex].RestrictCylinderABC_Max[j].x;
      Components[comp].RestrictCylinderABC_Min[j].y=Components[SolventIndex].RestrictCylinderABC_Min[j].y;
      Components[comp].RestrictCylinderABC_Max[j].y=Components[SolventIndex].RestrictCylinderABC_Max[j].y;
      Components[comp].RestrictCylinderABC_Min[j].z=Components[SolventIndex].RestrictCylinderABC_Min[j].z;
      Components[comp].RestrictCylinderABC_Max[j].z=Components[SolventIndex].RestrictCylinderABC_Max[j].z;
      Components[comp].RestrictCylinderCenter[j].x =Components[SolventIndex].RestrictCylinderCenter[j].x; 
      Components[comp].RestrictCylinderCenter[j].y =Components[SolventIndex].RestrictCylinderCenter[j].y;
      Components[comp].RestrictCylinderCenter[j].z =Components[SolventIndex].RestrictCylinderCenter[j].z; 
      Components[comp].RestrictCylinderDirection[j]=Components[SolventIndex].RestrictCylinderDirection[j];
      Components[comp].RestrictCylinderRadius[j]   =Components[SolventIndex].RestrictCylinderRadius[j];   
    }

    Components[comp].RestrictMovesToSpheres=FALSE;
    for(j=0;j<MAX_NUMBER_OF_SPHERES;j++)
    {
      Components[comp].RestrictMovesToSphere[j] =Components[SolventIndex].RestrictMovesToSphere[j];
      Components[comp].RestrictSphereCenter[j].x=Components[SolventIndex].RestrictSphereCenter[j].x;
      Components[comp].RestrictSphereCenter[j].y=Components[SolventIndex].RestrictSphereCenter[j].y;
      Components[comp].RestrictSphereCenter[j].z=Components[SolventIndex].RestrictSphereCenter[j].z;
      Components[comp].RestrictSphereRadius[j]  =Components[SolventIndex].RestrictSphereRadius[j];  
    }
	
    Components[comp].RuizMonteroFactor=Components[SolventIndex].RuizMonteroFactor;
    Components[comp].UmbrellaFactor=Components[SolventIndex].UmbrellaFactor;
    Components[comp].InvertBlockPockets=Components[SolventIndex].InvertBlockPockets;
	
    Components[comp].Intra14VDWScalingValue=Components[SolventIndex].Intra14VDWScalingValue;
    Components[comp].Intra14ChargeChargeScalingValue=Components[SolventIndex].Intra14ChargeChargeScalingValue;
	
    Components[comp].AnisotropicType=Components[SolventIndex].AnisotropicType;
    Components[comp].BiasingDirection=Components[SolventIndex].BiasingDirection;

  }

  // Update components definition: the transient molecule have the same definition as solvent, by definition.
  for(comp=NumberOfComponents;comp<NumberOfComponents+NumberExtraComponents;comp++)
  {
	   index_salt++;
	   Components[comp].CriticalTemperature=Components[SolventIndex].CriticalTemperature;
	   Components[comp].CriticalPressure=Components[SolventIndex].CriticalPressure;
	   Components[comp].AcentricFactor= Components[SolventIndex].AcentricFactor;
				  
	   Components[comp].NumberOfRigidAtoms=Components[SolventIndex].NumberOfRigidAtoms;
	   Components[comp].NumberOfFlexibleAtoms=Components[SolventIndex].NumberOfFlexibleAtoms;
	   
	   if(Components[comp].ExtraFrameworkMolecule)
		NumberOfCationComponents++;
	   else
		NumberOfAdsorbateComponents++;
		
	   strcpy(Components[comp].Name, "TransientMolecule");
	   sprintf(index_salt_store,"%d",SaltIndex[CurrentAlchemicalReaction][index_salt]);
	   strcat(Components[comp].Name,index_salt_store);
	    
	   // allocate charility-centers
	   Components[comp].Chirality=(int*)calloc(Components[comp].NumberOfAtoms,sizeof(int));
	   Components[comp].ChiralityType=(int*)calloc(Components[comp].NumberOfAtoms,sizeof(int));
	   Components[comp].ChiralA=(int*)calloc(Components[comp].NumberOfAtoms,sizeof(int));
	   Components[comp].ChiralB=(int*)calloc(Components[comp].NumberOfAtoms,sizeof(int));
	   Components[comp].ChiralC=(int*)calloc(Components[comp].NumberOfAtoms,sizeof(int));
	   Components[comp].ChiralD=(int*)calloc(Components[comp].NumberOfAtoms,sizeof(int));
	   
	   for(j=0;j<Components[comp].NumberOfAtoms;j++)
	   {
		 Components[comp].ChiralityType[j]=NO_CHARILITY;
		 Components[comp].Chirality[j]=FALSE;
		 Components[comp].ChiralA[j]=0;
		 Components[comp].ChiralB[j]=0;
		 Components[comp].ChiralC[j]=0;
		 Components[comp].ChiralD[j]=0;
	   }
	  
	  Components[comp].Fixed=(int*)calloc(Components[comp].NumberOfAtoms,sizeof(int));
	  Components[comp].Type=(int*)calloc(Components[comp].NumberOfAtoms,sizeof(int));
      Components[comp].Charge=(REAL*)calloc(Components[comp].NumberOfAtoms,sizeof(REAL));
      Components[comp].Connectivity=(int*)calloc(Components[comp].NumberOfAtoms,sizeof(int));
      
		  // allocate bond-connectivity matrix
		Components[comp].ConnectivityMatrix=(int**)calloc(Components[comp].NumberOfAtoms,sizeof(int*));
	  for(j=0;j<Components[comp].NumberOfAtoms;j++)
		Components[comp].ConnectivityMatrix[j]=(int*)calloc(Components[comp].NumberOfAtoms,sizeof(int));

	  // set the bond-connectivity to all FALSE at first
	  for(j=0;j<Components[comp].NumberOfAtoms;j++)
		for(k=0;k<Components[comp].NumberOfAtoms;k++)
		  Components[comp].ConnectivityMatrix[j][k]=FALSE;
		
	  Components[comp].Positions=(VECTOR*)calloc(Components[comp].NumberOfAtoms,sizeof(VECTOR));
      Components[comp].group=(int*)calloc(Components[comp].NumberOfAtoms,sizeof(int));
				 
      Components[comp].RMCMOL=(VECTOR*)calloc(Components[comp].NumberOfAtoms,sizeof(VECTOR));
      for(i=0;i<Components[comp].NumberOfAtoms;i++)
		Components[comp].RMCMOL[i]=Components[SolventIndex].RMCMOL[i];
      
      // Attribute number of groups
      Components[comp].NumberOfGroups=Components[SolventIndex].NumberOfGroups;
      
      // Define groups
      Components[comp].Groups=(GROUP_DEFINITION*)calloc(Components[comp].NumberOfGroups,sizeof(GROUP_DEFINITION));
	  
	  for(j=0;j<Components[comp].NumberOfGroups;j++)
	  {
		  Components[comp].Groups[j]=Components[SolventIndex].Groups[j];
		  //Components[comp].Groups[j].Rigid = Components[SolventIndex].Groups[j].Rigid;
		  //Components[comp].Groups[j].NumberOfGroupAtoms=Components[SolventIndex].Groups[j].NumberOfGroupAtoms;
		  //Components[comp].Groups[j].NumberOfPermanentDipoles=Components[SolventIndex].Groups[j].NumberOfPermanentDipoles;
		  //Components[comp].Groups[j].NumberOfPolarizabilities=Components[SolventIndex].Groups[j].NumberOfPolarizabilities;
		  
		  if(Components[comp].Groups[j].Rigid)
			Components[comp].NumberOfRigidAtoms+=Components[comp].Groups[j].NumberOfGroupAtoms;
		  else
			Components[comp].NumberOfFlexibleAtoms+=Components[comp].Groups[j].NumberOfGroupAtoms;
		  
		  if(Components[comp].Groups[j].NumberOfGroupAtoms>0)
          {
			  Components[comp].Groups[j].Atoms=(int*)calloc(Components[comp].Groups[j].NumberOfGroupAtoms,sizeof(int));
			  
			  for(k=0;k<Components[comp].Groups[j].NumberOfGroupAtoms;k++)
			  {
				  
				Components[comp].Groups[j].Atoms[k] = Components[SolventIndex].Groups[j].Atoms[k];
				// !!!! TRES IMPORTANT: TYPE D'ATOME
				index_type = (InitialPseudoAtoms) + index_salt*(Components[SolventIndex].NumberOfAtoms) + k;
				Components[comp].Type[k]= index_type;
				// ----------------------
				//Components[comp].Positions[k]=Components[SolventIndex].Positions[k];
				Components[comp].group[k]=Components[SolventIndex].group[k];
			  }
		  }
		 
		  if(Components[comp].Groups[j].NumberOfPermanentDipoles>0)
		  {
		    Components[comp].Groups[j].PermanentDipolePositions=(VECTOR*)calloc(Components[comp].Groups[j].NumberOfPermanentDipoles,sizeof(VECTOR));
		    Components[comp].Groups[j].PermanentDipoles=(VECTOR*)calloc(Components[comp].Groups[j].NumberOfPermanentDipoles,sizeof(VECTOR));
		    for(k=0;k<Components[comp].Groups[j].NumberOfPermanentDipoles;k++)
		    {
				Components[comp].Groups[j].PermanentDipolePositions[k].x=Components[SolventIndex].Groups[j].PermanentDipolePositions[k].x;
				Components[comp].Groups[j].PermanentDipolePositions[k].y=Components[SolventIndex].Groups[j].PermanentDipolePositions[k].y;
				Components[comp].Groups[j].PermanentDipolePositions[k].z=Components[SolventIndex].Groups[j].PermanentDipolePositions[k].z;
				Components[comp].Groups[j].PermanentDipoles[k].x=Components[comp].Groups[j].PermanentDipoles[k].x;
				Components[comp].Groups[j].PermanentDipoles[k].y=Components[comp].Groups[j].PermanentDipoles[k].y;
				Components[comp].Groups[j].PermanentDipoles[k].z=Components[comp].Groups[j].PermanentDipoles[k].z;
		    }
		  }
		  
		  if(Components[comp].Groups[j].NumberOfPolarizabilities>0)
		  {
		    Components[comp].Groups[j].PolarizabilityPositions=(VECTOR*)calloc(Components[comp].Groups[j].NumberOfPolarizabilities,sizeof(VECTOR));
		    Components[comp].Groups[j].Polarizabilites=(REAL_MATRIX3x3*)calloc(Components[comp].Groups[j].NumberOfPolarizabilities,sizeof(REAL_MATRIX3x3));
		    for(k=0;k<Components[comp].Groups[j].NumberOfPolarizabilities;k++)
		    {
				Components[comp].Groups[j].PolarizabilityPositions[k].x=Components[SolventIndex].Groups[j].PolarizabilityPositions[k].x;
				Components[comp].Groups[j].PolarizabilityPositions[k].y=Components[SolventIndex].Groups[j].PolarizabilityPositions[k].y;
				Components[comp].Groups[j].PolarizabilityPositions[k].z=Components[SolventIndex].Groups[j].PolarizabilityPositions[k].z;
				Components[comp].Groups[j].Polarizabilites[k].ax=Components[SolventIndex].Groups[j].Polarizabilites[k].ax;
				Components[comp].Groups[j].Polarizabilites[k].ay=Components[SolventIndex].Groups[j].Polarizabilites[k].ay;
				Components[comp].Groups[j].Polarizabilites[k].az=Components[SolventIndex].Groups[j].Polarizabilites[k].az;
				Components[comp].Groups[j].Polarizabilites[k].by=Components[SolventIndex].Groups[j].Polarizabilites[k].by;
				Components[comp].Groups[j].Polarizabilites[k].bz=Components[SolventIndex].Groups[j].Polarizabilites[k].bz;
				Components[comp].Groups[j].Polarizabilites[k].cz=Components[SolventIndex].Groups[j].Polarizabilites[k].cz;
		    }
		  }
	  }
	   
	  // fill in charges
      for(j=0;j<Components[comp].NumberOfAtoms;j++)
        Components[comp].Charge[j]=PseudoAtoms[Components[SolventIndex].Type[j]].Charge1;
      
      Components[comp].HasCharges=Components[SolventIndex].HasCharges;
      for(j=0;j<Components[comp].NumberOfAtoms;j++)
        Components[comp].HasCharges|=PseudoAtoms[Components[SolventIndex].Type[j]].HasCharges;
      
      Components[comp].NumberOfCharges=0;
      for(j=0;j<Components[comp].NumberOfAtoms;j++)
        if(PseudoAtoms[Components[comp].Type[j]].HasCharges) Components[comp].NumberOfCharges++;
      
      Components[comp].IsPolarizable=FALSE;
      for(j=0;j<Components[comp].NumberOfAtoms;j++)
        Components[comp].IsPolarizable|=PseudoAtoms[Components[SolventIndex].Type[j]].IsPolarizable;
	  
	  // Allocate the number of interaction
	  Components[comp].NumberOfChiralityCenters=Components[SolventIndex].NumberOfChiralityCenters;
	  Components[comp].NumberOfBonds=Components[SolventIndex].NumberOfBonds;
	  Components[comp].NumberOfBondDipoles=Components[SolventIndex].NumberOfBondDipoles;
	  Components[comp].NumberOfBends=Components[SolventIndex].NumberOfBends;
	  Components[comp].NumberOfUreyBradleys=Components[SolventIndex].NumberOfUreyBradleys;
	  Components[comp].NumberOfInversionBends=Components[SolventIndex].NumberOfInversionBends;
	  Components[comp].NumberOfTorsions=Components[SolventIndex].NumberOfTorsions;
	  Components[comp].NumberOfImproperTorsions=Components[SolventIndex].NumberOfImproperTorsions;
	  Components[comp].NumberOfBondBonds=Components[SolventIndex].NumberOfBondBonds;
	  Components[comp].NumberOfBondBends=Components[SolventIndex].NumberOfBondBends;
	  Components[comp].NumberOfBendBends=Components[SolventIndex].NumberOfBendBends;
	  Components[comp].NumberOfBondTorsions=Components[SolventIndex].NumberOfBondTorsions;
	  Components[comp].NumberOfBendTorsions=Components[SolventIndex].NumberOfBendTorsions;
	  Components[comp].NumberOfIntraVDW=Components[SolventIndex].NumberOfIntraVDW;
	  Components[comp].NumberOfIntraChargeCharge=Components[SolventIndex].NumberOfIntraChargeCharge;
	  Components[comp].NumberOfIntraChargeBondDipole=Components[SolventIndex].NumberOfIntraChargeBondDipole;
	  Components[comp].NumberOfIntraBondDipoleBondDipole=Components[SolventIndex].NumberOfIntraBondDipoleBondDipole;
	  
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
	  
	  // Attribute Bond-data
	  if(Components[comp].NumberOfBonds>0)
	  {
		for(j=0;j<Components[comp].NumberOfBonds;j++)
		{
			// Attribute connectivity
			A=Components[SolventIndex].Bonds[j].A;
			B=Components[SolventIndex].Bonds[j].B;
		    Components[comp].Bonds[j].A=A;
		    Components[comp].Bonds[j].B=B;
		    Components[comp].Connectivity[A]++;
		    Components[comp].Connectivity[B]++;
		    Components[comp].ConnectivityMatrix[A][B]=TRUE;
		    Components[comp].ConnectivityMatrix[B][A]=TRUE;

			// Attribute bondtype	
			Components[comp].BondType[j]=Components[SolventIndex].BondType[j];
	
			for(k=0;k<BondTypes[Components[comp].BondType[j]].nr_args;k++)
			{
				Components[comp].BondArguments[j][k]=Components[SolventIndex].BondArguments[j][k];
			}
			
			switch(Components[comp].BondType[j])
			{
			case HARMONIC_BOND:
			  // 0.5*p0*SQR(r-p1);
			  // ===============================================
			  // p_0/k_B [K/A^2]   force constant
			  // p_1     [A]       reference bond distance
			  Components[comp].BondArguments[j][0]=Components[SolventIndex].BondArguments[j][0];
			  break;
			case CORE_SHELL_SPRING:
			  // 0.5*p0*SQR(r);
			  // ===============================================
			  // p_0/k_B [K/A^2]   force constant
			  Components[comp].BondArguments[j][0]=Components[SolventIndex].BondArguments[j][0];
			  break;
			case MORSE_BOND:
			  // p_0*[(1.0-{exp(-p_1*(r-p_2))})^2-1.0]
			  // ===============================================
			  // p_0/k_B [K]       force constant
			  // p_1     [A^-1]    parameter
			  // p_2     [A]       reference bond distance
			  Components[comp].BondArguments[j][0]=Components[SolventIndex].BondArguments[j][0];
			  break;
			case LJ_12_6_BOND:
			  // A/r_ij^12-B/r_ij^6
			  // ===============================================
			  // p_0/k_B [K A^12]
			  // p_1/k_B [K A^6]
			  Components[comp].BondArguments[j][0]=Components[SolventIndex].BondArguments[j][0];
			  Components[comp].BondArguments[j][1]=Components[SolventIndex].BondArguments[j][1];
			  break;
			case LENNARD_JONES_BOND:
			  // 4*p_0*((p_1/r)^12-(p_1/r)^6)
			  // ===============================================
			  // p_0/k_B [K]
			  // p_1     [A]
			  Components[comp].BondArguments[j][0]= Components[SolventIndex].BondArguments[j][0];
			  break;
			case BUCKINGHAM_BOND:
			  // p_0*exp(-p_1 r)-p_2/r^6
			  // ===============================================
			  // p_0/k_B [K]
			  // p_1     [A^-1]
			  // p_2/k_B [K A^6]
			  Components[comp].BondArguments[j][0]=Components[SolventIndex].BondArguments[j][0];
			  Components[comp].BondArguments[j][2]=Components[SolventIndex].BondArguments[j][2];
			  break;
			case RESTRAINED_HARMONIC_BOND:
			  // 0.5*p_0*(r-p_1)^2                   |r-p_1|<=p_2
			  // 0.5*p_0*p_2^2+p_0*p_2*(|r-p_1|-p_2) |r-p_1|>p_2
			  // ===============================================
			  // p_0/k_B [K/A^2]
			  // p_1     [A]
			  // p_2     [A]
			  Components[comp].BondArguments[j][0]=Components[SolventIndex].BondArguments[j][0];
			  break;
			case QUARTIC_BOND:
			  // (1/2)*p_0*(r-p_1)^2+(1/3)*p_2*(r-p_1)^3+(1/4)*p_3*(r-p_1)^4
			  // ===========================================================
			  // p_0/k_B [K/A^2]
			  // p_1     [A]
			  // p_2/k_B [K/A^3]
			  // p_3/k_B [K/A^4]
			  Components[comp].BondArguments[j][0]=Components[SolventIndex].BondArguments[j][0];
			  Components[comp].BondArguments[j][2]=Components[SolventIndex].BondArguments[j][2];
			  Components[comp].BondArguments[j][3]=Components[SolventIndex].BondArguments[j][3];
			  break;
			case CFF_QUARTIC_BOND:
			  // p_0*(r-p_1)^2+p_2*(r-p_1)^3+p_3*(r-p_1)^4
			  // ===============================================
			  // p_0/k_B [K/A^2]
			  // p_1     [A]
			  // p_2/k_B [K/A^3]
			  // p_3/k_B [K/A^4]
			  Components[comp].BondArguments[j][0]=Components[SolventIndex].BondArguments[j][0];
			  Components[comp].BondArguments[j][2]=Components[SolventIndex].BondArguments[j][2];
			  Components[comp].BondArguments[j][3]=Components[SolventIndex].BondArguments[j][3];
			  break;
			case MM3_BOND:
			  // p_0*(r-p_1)^2*(1.0-2.55*(r-p_1)+(7.0/12.0)*2.55^2*(r-p_1)^2)
			  // ============================================================
			  // p_0     [mdyne/A molecule]
			  // p_1     [A]
			  Components[comp].BondArguments[j][0]=Components[SolventIndex].BondArguments[j][0];
			  break;
			case RIGID_BOND:
			  A=Components[SolventIndex].Bonds[j].A;
			  B=Components[SolventIndex].Bonds[j].B;
			  dr.x=Components[SolventIndex].Positions[A].x-Components[SolventIndex].Positions[B].x;
			  dr.y=Components[SolventIndex].Positions[A].y-Components[SolventIndex].Positions[B].y;
			  dr.z=Components[SolventIndex].Positions[A].z-Components[SolventIndex].Positions[B].z;
			  Components[comp].BondArguments[j][0]=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
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

	  // Attribute bond-dipole date.
	  if(Components[comp].NumberOfBondDipoles>0)
	  {
		for(j=0;j<Components[comp].NumberOfBondDipoles;j++)
		{
		  Components[comp].BondDipoles[j].A = Components[SolventIndex].BondDipoles[j].A;
		  Components[comp].BondDipoles[j].B = Components[SolventIndex].BondDipoles[j].B;
		  Components[comp].BondDipoleMagnitude[j]=Components[SolventIndex].BondDipoleMagnitude[j];
		}
	  }
	  
	  
	  // Bend-data
	  if(Components[comp].NumberOfBends>0)
	  {
		for(j=0;j<Components[comp].NumberOfBends;j++)
		{
		 
			Components[comp].Bends[j].A=Components[SolventIndex].Bends[j].A;
			Components[comp].Bends[j].B=Components[SolventIndex].Bends[j].B;
			Components[comp].Bends[j].C=Components[SolventIndex].Bends[j].C;
			
			// Attribute bend-type	
			Components[comp].BendType[j]=Components[SolventIndex].BendType[j];

			for(k=0;j<BendTypes[Components[comp].BendType[j]].nr_args;j++)
			{
				Components[comp].BendArguments[j][k]=Components[SolventIndex].BendArguments[j][k];
			}

			switch(Components[comp].BendType[j])
			{
				case HARMONIC_BEND:
				  // (1/2)p_0*(theta-p_1)^2
				  // ===============================================
				  // p_0/k_B [K/rad^2]
				  // p_1     [degrees]
				  Components[comp].BendArguments[j][0]=Components[SolventIndex].BendArguments[j][0];
				  Components[comp].BendArguments[j][1]=Components[SolventIndex].BendArguments[j][1];
				  break;
				case CORE_SHELL_BEND:
				  // (1/2)p_0*(theta-p_1)^2
				  // ===============================================
				  // p_0/k_B [K/rad^2]
				  // p_1     [degrees]
				  Components[comp].BendArguments[j][0]=Components[SolventIndex].BendArguments[j][0];
				  Components[comp].BendArguments[j][1]=Components[SolventIndex].BendArguments[j][1];
				  break;
				case QUARTIC_BEND:
				  // (1/2)p_0*(theta-p_1)^2+(1/3)*p_2*(theta-p_1)^3+(1/4)*p_2*(theta-p_1)^4
				  // ======================================================================
				  // p_0/k_B [K/rad^2]
				  // p_1     [degrees]
				  // p_2/k_B [K/rad^3]
				  // p_3/k_B [K/rad^4]
				  Components[comp].BendArguments[j][0]=Components[SolventIndex].BendArguments[j][0];
				  Components[comp].BendArguments[j][1]=Components[SolventIndex].BendArguments[j][1];
				  Components[comp].BendArguments[j][2]=Components[SolventIndex].BendArguments[j][2];
				  Components[comp].BendArguments[j][3]=Components[SolventIndex].BendArguments[j][3];
				  break;
				case CFF_QUARTIC_BEND:
				  // p_0*(theta-p_1)^2+p_2*(theta-p_1)^3+p_3*(theta-p_1)^4
				  // =====================================================
				  // p_0/k_B [K/rad^2]
				  // p_1     [degrees]
				  // p_2/k_B [K/rad^3]
				  // p_3/k_B [K/rad^4]
				  Components[comp].BendArguments[j][0]=Components[SolventIndex].BendArguments[j][0];
				  Components[comp].BendArguments[j][1]=Components[SolventIndex].BendArguments[j][1];
				  Components[comp].BendArguments[j][2]=Components[SolventIndex].BendArguments[j][2];
				  Components[comp].BendArguments[j][3]=Components[SolventIndex].BendArguments[j][3];
				  break;
				case HARMONIC_COSINE_BEND:
				  // (1/2)*p_0*(cos(theta)-cos(p_1))^2
				  // ===============================================
				  // p_0/k_B [K]
				  // p_1     [degrees]
				  Components[comp].BendArguments[j][0]=Components[SolventIndex].BendArguments[j][0];
				  Components[comp].BendArguments[j][1]=Components[SolventIndex].BendArguments[j][1];
				  break;
				case COSINE_BEND:
				  // p_0*(1+cos(p_1*theta-p_2))
				  // ===============================================
				  // p_0/k_B [K]
				  // p_1     [-]
				  // p_2     [degrees]
				  Components[comp].BendArguments[j][0]=Components[SolventIndex].BendArguments[j][0];
				  Components[comp].BendArguments[j][2]=Components[SolventIndex].BendArguments[j][2];
				  break;
				case TAFIPOLSKY_BEND:
				  // 0.5*p_0*(1+cos(theta))*(1+cos(2*theta))
				  // ===============================================
				  // p_0/k_B [K]
				  Components[comp].BendArguments[j][0]=Components[SolventIndex].BendArguments[j][0];
				  break;
				case MM3_BEND:
				case MM3_IN_PLANE_BEND:
				  // p_0*(theta-p_1)^2(1-0.014*(theta-p_1)+5.6e-5*(theta-p_1)^2-7e-7*(theta-p_1)^3+2.2e-8(theta-p_1)^4)
				  // =================================================================================================
				  // p_0/k_B [mdyne A/rad^2]
				  // p_1     [degrees]
				  Components[comp].BendArguments[j][0]=Components[SolventIndex].BendArguments[j][0];
				  Components[comp].BendArguments[j][1]=Components[SolventIndex].BendArguments[j][1];
				  break;
				case FIXED_BEND:
				  Components[comp].BendArguments[j][0]=Components[SolventIndex].BendArguments[j][0];
				  Components[comp].NumberOfConstraintBends++;
				  break;
				default:
				  fprintf(stderr, "Undefined Bend potential in routine 'ReadComponentDefinition' ('molecule.c')\n");
				  exit(0);
				  break;
			 }
		}
	  }
	  
	  // UreyBradley-data
	  if(Components[comp].NumberOfUreyBradleys>0)
	  {
		for(j=0;j<Components[comp].NumberOfUreyBradleys;j++)
		{
		  Components[comp].UreyBradleys[j].A=Components[SolventIndex].UreyBradleys[j].A;
		  Components[comp].UreyBradleys[j].B=Components[SolventIndex].UreyBradleys[j].B;
		  Components[comp].UreyBradleys[j].C=Components[SolventIndex].UreyBradleys[j].C;

		  // Attribute Urey-Bradley.
	      Components[comp].UreyBradleyType[j]=Components[SolventIndex].UreyBradleyType[j];

		  for(k=0;k<UreyBradleyTypes[Components[comp].UreyBradleyType[j]].nr_args;k++)
		  {
			Components[comp].UreyBradleyArguments[j][k]=Components[SolventIndex].UreyBradleyArguments[j][k];
		  }

		  switch(Components[comp].UreyBradleyType[j])
		  {
			case HARMONIC_UREYBRADLEY:
			  // 0.5*p0*SQR(r-p1);
			  // ===============================================
			  // p_0/k_B [K/A^2]   force constant
			  // p_1     [A]       reference bond distance
			  Components[comp].UreyBradleyArguments[j][0]=Components[SolventIndex].UreyBradleyArguments[j][0];
			  break;
			case MORSE_UREYBRADLEY:
			  // p_0*[(1.0-{exp(-p_1*(r-p_2))})^2-1.0]
			  // ===============================================
			  // p_0/k_B [K]       force constant
			  // p_1     [A^-1]    parameter
			  // p_2     [A]       reference bond distance
			  Components[comp].UreyBradleyArguments[j][0]=Components[SolventIndex].UreyBradleyArguments[j][0];
			  break;
			case LJ_12_6_UREYBRADLEY:
			  // A/r_ij^12-B/r_ij^6
			  // ===============================================
			  // p_0/k_B [K A^12]
			  // p_1/k_B [K A^6]
			  Components[comp].UreyBradleyArguments[j][0]=Components[SolventIndex].UreyBradleyArguments[j][0];
			  Components[comp].UreyBradleyArguments[j][1]=Components[SolventIndex].UreyBradleyArguments[j][1];
			  break;
			case LENNARD_JONES_UREYBRADLEY:
			  // 4*p_0*((p_1/r)^12-(p_1/r)^6)
			  // ===============================================
			  // p_0/k_B [K]
			  // p_1     [A]
			  Components[comp].UreyBradleyArguments[j][0]=Components[SolventIndex].UreyBradleyArguments[j][0];
			  break;
			case BUCKINGHAM_UREYBRADLEY:
			  // p_0*exp(-p_1 r)-p_2/r^6
			  // ===============================================
			  // p_0/k_B [K]
			  // p_1     [A^-1]
			  // p_2/k_B [K A^6]
			  Components[comp].UreyBradleyArguments[j][0]=Components[SolventIndex].UreyBradleyArguments[j][0];
			  Components[comp].UreyBradleyArguments[j][2]=Components[SolventIndex].UreyBradleyArguments[j][2];
			  break;
			case RESTRAINED_HARMONIC_UREYBRADLEY:
			  // 0.5*p_0*(r-p_1)^2                   |r-p_1|<=p_2
			  // 0.5*p_0*p_2^2+p_0*p_2*(|r-p_1|-p_2) |r-p_1|>p_2
			  // ===============================================
			  // p_0/k_B [K/A^2]
			  // p_1     [A]
			  // p_2     [A]
			  Components[comp].UreyBradleyArguments[j][0]=Components[SolventIndex].UreyBradleyArguments[j][0];
			  break;
			case QUARTIC_UREYBRADLEY:
			  // (1/2)*p_0*(r-p_1)^2+(1/3)*p_2*(r-p_1)^3+(1/4)*p_3*(r-p_1)^4
			  // ===========================================================
			  // p_0/k_B [K/A^2]
			  // p_1     [A]
			  // p_2/k_B [K/A^3]
			  // p_3/k_B [K/A^4]
			  Components[comp].UreyBradleyArguments[j][0]=Components[SolventIndex].UreyBradleyArguments[j][0];
			  Components[comp].UreyBradleyArguments[j][2]=Components[SolventIndex].UreyBradleyArguments[j][2];
			  Components[comp].UreyBradleyArguments[j][3]=Components[SolventIndex].UreyBradleyArguments[j][3];
			  break;
			case CFF_QUARTIC_UREYBRADLEY:
			  // p_0*(r-p_1)^2+p_2*(r-p_1)^3+p_3*(r-p_1)^4
			  // ===============================================
			  // p_0/k_B [K/A^2]
			  // p_1     [A]
			  // p_2/k_B [K/A^3]
			  // p_3/k_B [K/A^4]
			  Components[comp].UreyBradleyArguments[j][0]=Components[SolventIndex].UreyBradleyArguments[j][0];
			  Components[comp].UreyBradleyArguments[j][2]=Components[SolventIndex].UreyBradleyArguments[j][2];
			  Components[comp].UreyBradleyArguments[j][3]=Components[SolventIndex].UreyBradleyArguments[j][3];
			  break;
			case MM3_UREYBRADLEY:
			  // p_0*(r-p_1)^2*(1.0-2.55*(r-p_1)+(7.0/12.0)*2.55^2*(r-p_1)^2)
			  // ============================================================
			  // p_0     [mdyne/A molecule]
			  // p_1     [A]
			  Components[comp].UreyBradleyArguments[j][0]=Components[SolventIndex].UreyBradleyArguments[j][0];
			  break;
			case RIGID_UREYBRADLEY:
			  A=Components[SolventIndex].UreyBradleys[j].A;
			  B=Components[SolventIndex].UreyBradleys[j].B;
			  dr.x=Components[SolventIndex].Positions[A].x-Components[SolventIndex].Positions[B].x;
			  dr.y=Components[SolventIndex].Positions[A].y-Components[SolventIndex].Positions[B].y;
			  dr.z=Components[SolventIndex].Positions[A].z-Components[SolventIndex].Positions[B].z;
			  Components[comp].BondArguments[j][0]=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
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
	  
	    // inversion Bend-data
	  if(Components[comp].NumberOfInversionBends>0)
	  {
		for(j=0;j<Components[comp].NumberOfInversionBends;j++)
		{
			Components[comp].InversionBends[j].A=Components[SolventIndex].InversionBends[j].A;
			Components[comp].InversionBends[j].B=Components[SolventIndex].InversionBends[j].B;
			Components[comp].InversionBends[j].C=Components[SolventIndex].InversionBends[j].C;
			Components[comp].InversionBends[j].D=Components[SolventIndex].InversionBends[j].D;
			
			// Attribute Inversion Bend
			Components[comp].InversionBendType[j]=Components[SolventIndex].InversionBendType[j];
			

		    for(k=0;k<InversionBendTypes[Components[comp].InversionBendType[j]].nr_args;k++)
			{
				Components[comp].InversionBendArguments[j][k]=Components[SolventIndex].InversionBendArguments[j][k];
			}

		    switch(Components[comp].InversionBendType[j])
		    {
				case HARMONIC_INVERSION:
				case HARMONIC_INVERSION2:
				  // (1/2)*p_0*(chi-p_1)^2
				  // ===============================================
				  // p_0/k_B [K/rad^2]
				  // p_1     [degrees]
				  Components[comp].InversionBendArguments[j][0]=Components[SolventIndex].InversionBendArguments[j][0];
				  Components[comp].InversionBendArguments[j][1]=Components[SolventIndex].InversionBendArguments[j][1];
				  break;
				case HARMONIC_COSINE_INVERSION:
				case HARMONIC_COSINE_INVERSION2:
				  // (1/2)*p_0*(cos(phi)-cos(p_1))^2
				  // ===============================================
				  // p_0/k_B [K]
				  // p_1     [degrees]
				  Components[comp].InversionBendArguments[j][0]=Components[SolventIndex].InversionBendArguments[j][0];
				  Components[comp].InversionBendArguments[j][1]=Components[SolventIndex].InversionBendArguments[j][1];
				  break;
				case PLANAR_INVERSION:
				case PLANAR_INVERSION2:
				  // (1/2)*p_0*(1-cos(phi))
				  // ===============================================
				  // p_0/k_B [K]
				  Components[comp].InversionBendArguments[j][0]=Components[SolventIndex].InversionBendArguments[j][0];
				  break;
				case MM3_INVERSION:
				  // p_0*(theta-p_1)^2(1-0.014*(theta-p_1)+5.6e-5*(theta-p_1)^2-7e-7*(theta-p_1)^3+2.2e-8(theta-p_1)^4)
				  // =================================================================================================
				  // p_0/k_B [mdyne A/rad^2]
				  // p_1     [degrees]
				  Components[comp].InversionBendArguments[j][0]=Components[SolventIndex].InversionBendArguments[j][0];
				  break;
				case FIXED_INVERSION_BEND:
				  Components[comp].InversionBendArguments[j][0]=Components[SolventIndex].InversionBendArguments[j][0];
				  Components[comp].NumberOfConstraintInversionBends++;
				  break;
				default:
				  fprintf(stderr, "Undefined Inversion potential in routine 'ReadComponentDefinition' ('molecule.c')\n");
				  exit(0);
				  break;
		   }
		}
	  }
	  
	  // Torsion-data
	  if(Components[comp].NumberOfTorsions>0)
	  {
		for(j=0;j<Components[comp].NumberOfTorsions;j++)
		{
			Components[comp].Torsions[j].A=Components[SolventIndex].Torsions[j].A;
			Components[comp].Torsions[j].B=Components[SolventIndex].Torsions[j].B;
			Components[comp].Torsions[j].C=Components[SolventIndex].Torsions[j].C;
			Components[comp].Torsions[j].D=Components[SolventIndex].Torsions[j].D;
      
			// Attribute Torsion data
			Components[comp].TorsionType[j]=Components[SolventIndex].TorsionType[j];

			for(k=0;k<TorsionTypes[Components[comp].TorsionType[j]].nr_args;k++)
			{
				Components[comp].TorsionArguments[j][k]=Components[SolventIndex].TorsionArguments[j][k];
			}

			switch(Components[comp].TorsionType[j])
			{
				case HARMONIC_DIHEDRAL:
				  // (1/2)*p_0*(phi-p_1)^2
				  // ===============================================
				  // p_0/k_B [K/rad^2]
				  // p_1     [degrees]
				  Components[comp].TorsionArguments[j][0]=Components[SolventIndex].TorsionArguments[j][0];
				  Components[comp].TorsionArguments[j][1]=Components[SolventIndex].TorsionArguments[j][1];
				  break;
				case HARMONIC_COSINE_DIHEDRAL:
				  // (1/2)*p_0*(cos(phi)-cos(p_1))^2
				  // ===============================================
				  // p_0/k_B [K]
				  // p_1     [degrees]
				  Components[comp].TorsionArguments[j][0]=Components[SolventIndex].TorsionArguments[j][0];
				  Components[comp].TorsionArguments[j][1]=Components[SolventIndex].TorsionArguments[j][1];
				  break;
				case THREE_COSINE_DIHEDRAL:
				  // (1/2)*p_0*(1+cos(phi))+(1/2)*p_1*(1-cos(2*phi))+(1/2)*p_2*(1+cos(3*phi))
				  // ========================================================================
				  // p_0/k_B [K]
				  // p_1/k_B [K]
				  // p_2/k_B [K]
				  Components[comp].TorsionArguments[j][0]=Components[SolventIndex].TorsionArguments[j][0];
				  Components[comp].TorsionArguments[j][1]=Components[SolventIndex].TorsionArguments[j][1];
				  Components[comp].TorsionArguments[j][2]=Components[SolventIndex].TorsionArguments[j][2];
				  break;
				case MM3_DIHEDRAL:
				  // (1/2)*p_0*(1+cos(phi))+(1/2)*p_1*(1-cos(2*phi))+(1/2)*p_2*(1+cos(3*phi))
				  // ========================================================================
				  // p_0     [kcal/mol]
				  // p_1     [kcal/mol]
				  // p_2     [kcal/mol]
				  Components[comp].TorsionArguments[j][0]*=Components[SolventIndex].TorsionArguments[j][0];
				  Components[comp].TorsionArguments[j][1]*=Components[SolventIndex].TorsionArguments[j][1];
				  Components[comp].TorsionArguments[j][2]*=Components[SolventIndex].TorsionArguments[j][2];
				case CVFF_BLOCKED_DIHEDRAL:
				  // 
				  // ========================================================================
				  // p_0     [rad]
				  // p_1     [K]
				  // p_2     [-]
				  // p_3     [rad]
				  // p_4     [rad]
				  Components[comp].TorsionArguments[j][1]=Components[SolventIndex].TorsionArguments[j][1];
				  break;
				case CFF_DIHEDRAL:
				  // p_0*(1-cos(phi))+p_1*(1-cos(2*phi))+p_2*(1-cos(3*phi))
				  // ======================================================
				  // p_0/k_B [K]
				  // p_1/k_B [K]
				  // p_2/k_B [K]
				  Components[comp].TorsionArguments[j][0]=Components[SolventIndex].TorsionArguments[j][0];
				  Components[comp].TorsionArguments[j][1]=Components[SolventIndex].TorsionArguments[j][1];
				  Components[comp].TorsionArguments[j][2]=Components[SolventIndex].TorsionArguments[j][2];
				  break;
				case CFF_DIHEDRAL2:
				  // p_0*(1+cos(phi))+p_1*(1+cos(2*phi))+p_2*(1+cos(3*phi))
				  // ======================================================
				  // p_0/k_B [K]
				  // p_1/k_B [K]
				  // p_2/k_B [K]
				  Components[comp].TorsionArguments[j][0]=Components[SolventIndex].TorsionArguments[j][0];
				  Components[comp].TorsionArguments[j][1]=Components[SolventIndex].TorsionArguments[j][1];
				  Components[comp].TorsionArguments[j][2]=Components[SolventIndex].TorsionArguments[j][2];
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
				  Components[comp].TorsionArguments[j][0]=Components[SolventIndex].TorsionArguments[j][0];
				  Components[comp].TorsionArguments[j][1]=Components[SolventIndex].TorsionArguments[j][1];
				  Components[comp].TorsionArguments[j][2]=Components[SolventIndex].TorsionArguments[j][2];
				  Components[comp].TorsionArguments[j][3]=Components[SolventIndex].TorsionArguments[j][3];
				  Components[comp].TorsionArguments[j][4]=Components[SolventIndex].TorsionArguments[j][4];
				  Components[comp].TorsionArguments[j][5]=Components[SolventIndex].TorsionArguments[j][5];
				  break;
				case TRAPPE_DIHEDRAL:
				  // p_0[0]+p_1*(1+cos(phi))+p_2*(1-cos(2*phi))+p_3*(1+cos(3*phi))
				  // =============================================================
				  // p_0/k_B [K]
				  // p_1/k_B [K]
				  // p_2/k_B [K]
				  // p_3/k_B [K]
				  Components[comp].TorsionArguments[j][0]=Components[SolventIndex].TorsionArguments[j][0];
				  Components[comp].TorsionArguments[j][1]=Components[SolventIndex].TorsionArguments[j][1];
				  Components[comp].TorsionArguments[j][2]=Components[SolventIndex].TorsionArguments[j][2];
				  Components[comp].TorsionArguments[j][3]=Components[SolventIndex].TorsionArguments[j][3];
				  break;
				case TRAPPE_DIHEDRAL_EXTENDED:
				  // p_0[0]+p_1*cos(phi)+p_2*cos(2*phi)+p_3*cos(3*phi)+p_4*cos(4*phi)
				  // ================================================================
				  // p_0/k_B [K]
				  // p_1/k_B [K]
				  // p_2/k_B [K]
				  // p_3/k_B [K]
				  // p_4/k_B [K]
				  Components[comp].TorsionArguments[j][0]=Components[SolventIndex].TorsionArguments[j][0];
				  Components[comp].TorsionArguments[j][1]=Components[SolventIndex].TorsionArguments[j][1];
				  Components[comp].TorsionArguments[j][2]=Components[SolventIndex].TorsionArguments[j][2];
				  Components[comp].TorsionArguments[j][3]=Components[SolventIndex].TorsionArguments[j][3];
				  Components[comp].TorsionArguments[j][4]=Components[SolventIndex].TorsionArguments[j][4];
				  break;
				case MOD_TRAPPE_DIHEDRAL:
				  /* Salvador modification: 16/08/2016
				   add phase in cos function:
				   p_0+p_1*(1+cos(phi-p_4))+p_2*(1-cos(2*(phi-p_4)))+p_3*(1+cos(3*(phi-p_4)))
				  */
				  Components[comp].TorsionArguments[j][0]=Components[SolventIndex].TorsionArguments[j][0];
				  Components[comp].TorsionArguments[j][1]=Components[SolventIndex].TorsionArguments[j][1];
				  Components[comp].TorsionArguments[j][2]=Components[SolventIndex].TorsionArguments[j][2];
				  Components[comp].TorsionArguments[j][3]=Components[SolventIndex].TorsionArguments[j][3];
				  Components[comp].TorsionArguments[j][4]=Components[SolventIndex].TorsionArguments[j][4];
				  break;
				case CVFF_DIHEDRAL:
				  // p_0*(1+cos(p_1*phi-p_2))
				  // ========================
				  // p_0/k_B [K]
				  // p_1     [-]
				  // p_2     [degrees]
				  Components[comp].TorsionArguments[j][0]=Components[SolventIndex].TorsionArguments[j][0];
				  Components[comp].TorsionArguments[j][2]=Components[SolventIndex].TorsionArguments[j][2];
				  break;
				case OPLS_DIHEDRAL:
				  // (1/2)p_0[0]+(1/2)p_1*(1+cos(phi))+(1/2)p_2*(1-cos(2*phi))+(1/2)p_3*(1+cos(3*phi))
				  // =================================================================================
				  // p_0/k_B [K]
				  // p_1/k_B [K]
				  // p_2/k_B [K]
				  // p_3/k_B [K]
				  Components[comp].TorsionArguments[j][0]=Components[SolventIndex].TorsionArguments[j][0];
				  Components[comp].TorsionArguments[j][1]=Components[SolventIndex].TorsionArguments[j][1];
				  Components[comp].TorsionArguments[j][2]=Components[SolventIndex].TorsionArguments[j][2];
				  Components[comp].TorsionArguments[j][3]=Components[SolventIndex].TorsionArguments[j][3];
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
				  Components[comp].TorsionArguments[j][0]=Components[SolventIndex].TorsionArguments[j][0];
				  Components[comp].TorsionArguments[j][1]=Components[SolventIndex].TorsionArguments[j][1];
				  Components[comp].TorsionArguments[j][2]=Components[SolventIndex].TorsionArguments[j][2];
				  Components[comp].TorsionArguments[j][3]=Components[SolventIndex].TorsionArguments[j][3];
				  Components[comp].TorsionArguments[j][4]=Components[SolventIndex].TorsionArguments[j][4];
				  Components[comp].TorsionArguments[j][5]=Components[SolventIndex].TorsionArguments[j][5];
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
				  Components[comp].TorsionArguments[j][0]=Components[SolventIndex].TorsionArguments[j][0];
				  Components[comp].TorsionArguments[j][1]=Components[SolventIndex].TorsionArguments[j][1];
				  Components[comp].TorsionArguments[j][2]=Components[SolventIndex].TorsionArguments[j][2];
				  Components[comp].TorsionArguments[j][3]=Components[SolventIndex].TorsionArguments[j][3];
				  Components[comp].TorsionArguments[j][4]=Components[SolventIndex].TorsionArguments[j][4];
				  Components[comp].TorsionArguments[j][5]=Components[SolventIndex].TorsionArguments[j][5];
				  break;
				case FIXED_DIHEDRAL:
				  Components[comp].TorsionArguments[i][0]=Components[SolventIndex].TorsionArguments[i][0];
				  Components[comp].NumberOfConstraintTorsions++;
				  break;
				default:
				  fprintf(stderr, "Undefined Torsion potential in routine 'ReadComponentDefinition' ('molecule.c')\n");
				  exit(0);
				  break;
			}
		}
	  }

	  //  Improper Torsion-data
	  if(Components[comp].NumberOfImproperTorsions>0)
	  {
		for(j=0;j<Components[comp].NumberOfImproperTorsions;j++)
		{
            Components[comp].ImproperTorsions[j].A=Components[SolventIndex].ImproperTorsions[j].A;
            Components[comp].ImproperTorsions[j].B=Components[SolventIndex].ImproperTorsions[j].B;
            Components[comp].ImproperTorsions[j].C=Components[SolventIndex].ImproperTorsions[j].C;
            Components[comp].ImproperTorsions[j].D=Components[SolventIndex].ImproperTorsions[j].D;
             
            // Attribute Torsion data
			Components[comp].ImproperTorsionType[j]=Components[SolventIndex].ImproperTorsionType[j];
			 
			for(k=0;k<ImproperTorsionTypes[Components[comp].ImproperTorsionType[j]].nr_args;k++)
			{
				Components[comp].ImproperTorsionArguments[j][k]=Components[SolventIndex].ImproperTorsionArguments[j][k];
			}

			switch(Components[comp].ImproperTorsionType[j])
			{
			case HARMONIC_IMPROPER_DIHEDRAL:
			  // (1/2)*p_0*(phi-p_1)^2
			  // ===============================================
			  // p_0/k_B [K/rad^2]
			  // p_1     [degrees]
			  Components[comp].ImproperTorsionArguments[j][0]=Components[SolventIndex].ImproperTorsionArguments[j][0];
			  Components[comp].ImproperTorsionArguments[j][1]=Components[SolventIndex].ImproperTorsionArguments[j][1];
			  break;
			case HARMONIC_COSINE_IMPROPER_DIHEDRAL:
			  // (1/2)*p_0*(cos(phi)-cos(p_1))^2
			  // ===============================================
			  // p_0/k_B [K]
			  // p_1     [degrees]
			  Components[comp].ImproperTorsionArguments[j][0]=Components[SolventIndex].ImproperTorsionArguments[j][0];
			  Components[comp].ImproperTorsionArguments[j][1]=Components[SolventIndex].ImproperTorsionArguments[j][1];
			  break;
			case THREE_COSINE_IMPROPER_DIHEDRAL:
			  // (1/2)*p_0*(1+cos(phi))+(1/2)*p_1*(1-cos(2*phi))+(1/2)*p_2*(1+cos(3*phi))
			  // ========================================================================
			  // p_0/k_B [K]
			  // p_1/k_B [K]
			  // p_2/k_B [K]
			  Components[comp].ImproperTorsionArguments[j][0]=Components[SolventIndex].ImproperTorsionArguments[j][0];
			  Components[comp].ImproperTorsionArguments[j][1]=Components[SolventIndex].ImproperTorsionArguments[j][1];
			  Components[comp].ImproperTorsionArguments[j][2]=Components[SolventIndex].ImproperTorsionArguments[j][2];
			  break;
			case MM3_IMPROPER_DIHEDRAL:
			  // p_0*(1+cos(phi))+p_1*(1+cos(2*phi))+p_2*(1+cos(3*phi))
			  // ======================================================
			  // p_0/k_B [K]
			  // p_1/k_B [K]
			  // p_2/k_B [K]
			  Components[comp].ImproperTorsionArguments[j][0]=Components[SolventIndex].ImproperTorsionArguments[j][0];
			  Components[comp].ImproperTorsionArguments[j][1]=Components[SolventIndex].ImproperTorsionArguments[j][1];
			  Components[comp].ImproperTorsionArguments[j][2]=Components[SolventIndex].ImproperTorsionArguments[j][2];
			case CFF_IMPROPER_DIHEDRAL:
			  // p_0*(1-cos(phi))+p_1*(1-cos(2*phi))+p_2*(1-cos(3*phi))
			  // ======================================================
			  // p_0/k_B [K]
			  // p_1/k_B [K]
			  // p_2/k_B [K]
			  Components[comp].ImproperTorsionArguments[j][0]=Components[SolventIndex].ImproperTorsionArguments[j][0];
			  Components[comp].ImproperTorsionArguments[j][1]=Components[SolventIndex].ImproperTorsionArguments[j][1];
			  Components[comp].ImproperTorsionArguments[j][2]=Components[SolventIndex].ImproperTorsionArguments[j][2];
			  break;
			case CFF_IMPROPER_DIHEDRAL2:
			  // p_0*(1+cos(phi))+p_1*(1+cos(2*phi))+p_2*(1+cos(3*phi))
			  // ======================================================
			  // p_0/k_B [K]
			  // p_1/k_B [K]
			  // p_2/k_B [K]
			  Components[comp].ImproperTorsionArguments[j][0]=Components[SolventIndex].ImproperTorsionArguments[j][0];
			  Components[comp].ImproperTorsionArguments[j][1]=Components[SolventIndex].ImproperTorsionArguments[j][1];
			  Components[comp].ImproperTorsionArguments[j][2]=Components[SolventIndex].ImproperTorsionArguments[j][2];
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
			  Components[comp].ImproperTorsionArguments[j][0]=Components[SolventIndex].ImproperTorsionArguments[j][0];
			  Components[comp].ImproperTorsionArguments[j][1]=Components[SolventIndex].ImproperTorsionArguments[j][1];
			  Components[comp].ImproperTorsionArguments[j][2]=Components[SolventIndex].ImproperTorsionArguments[j][2];
			  Components[comp].ImproperTorsionArguments[j][3]=Components[SolventIndex].ImproperTorsionArguments[j][3];
			  Components[comp].ImproperTorsionArguments[j][4]=Components[SolventIndex].ImproperTorsionArguments[j][4];
			  Components[comp].ImproperTorsionArguments[j][5]=Components[SolventIndex].ImproperTorsionArguments[j][5];
			  break;
			case TRAPPE_IMPROPER_DIHEDRAL:
			  // p_0[0]+p_1*(1+cos(phi))+p_2*(1-cos(2*phi))+p_3*(1+cos(3*phi))
			  // =============================================================
			  // p_0/k_B [K]
			  // p_1/k_B [K]
			  // p_2/k_B [K]
			  // p_3/k_B [K]
			  Components[comp].ImproperTorsionArguments[j][0]=Components[SolventIndex].ImproperTorsionArguments[j][0];
			  Components[comp].ImproperTorsionArguments[j][1]=Components[SolventIndex].ImproperTorsionArguments[j][1];
			  Components[comp].ImproperTorsionArguments[j][2]=Components[SolventIndex].ImproperTorsionArguments[j][2];
			  Components[comp].ImproperTorsionArguments[j][3]=Components[SolventIndex].ImproperTorsionArguments[j][3];
			  break;
			case TRAPPE_IMPROPER_DIHEDRAL_EXTENDED:
			  // p_0[0]+p_1*cos(phi)+p_2*cos(2*phi)+p_3*cos(3*phi)+p_4*cos(4*phi)
			  // ================================================================
			  // p_0/k_B [K]
			  // p_1/k_B [K]
			  // p_2/k_B [K]
			  // p_3/k_B [K]
			  // p_4/k_B [K]
			  Components[comp].ImproperTorsionArguments[j][0]=Components[SolventIndex].ImproperTorsionArguments[j][0];
			  Components[comp].ImproperTorsionArguments[j][1]=Components[SolventIndex].ImproperTorsionArguments[j][1];
			  Components[comp].ImproperTorsionArguments[j][2]=Components[SolventIndex].ImproperTorsionArguments[j][2];
			  Components[comp].ImproperTorsionArguments[j][3]=Components[SolventIndex].ImproperTorsionArguments[j][3];
			  Components[comp].ImproperTorsionArguments[j][4]=Components[SolventIndex].ImproperTorsionArguments[j][4];
			  break;
			case CVFF_IMPROPER_DIHEDRAL:
			  // p_0*(1+cos(p_1*phi-p_2))
			  // ========================
			  // p_0/k_B [K]
			  // p_1     [-]
			  // p_2     [degrees]
			  Components[comp].ImproperTorsionArguments[j][0]=Components[SolventIndex].ImproperTorsionArguments[j][0];
			  Components[comp].ImproperTorsionArguments[j][2]=Components[SolventIndex].ImproperTorsionArguments[j][2];
			  break;
			case OPLS_IMPROPER_DIHEDRAL:
			  // (1/2)p_0[0]+(1/2)p_1*(1+cos(phi))+(1/2)p_2*(1-cos(2*phi))+(1/2)p_3*(1+cos(3*phi))
			  // =================================================================================
			  // p_0/k_B [K]
			  // p_1/k_B [K]
			  // p_2/k_B [K]
			  // p_3/k_B [K]
			  Components[comp].ImproperTorsionArguments[j][0]=Components[SolventIndex].ImproperTorsionArguments[j][0];
			  Components[comp].ImproperTorsionArguments[j][1]=Components[SolventIndex].ImproperTorsionArguments[j][1];
			  Components[comp].ImproperTorsionArguments[j][2]=Components[SolventIndex].ImproperTorsionArguments[j][2];
			  Components[comp].ImproperTorsionArguments[j][3]=Components[SolventIndex].ImproperTorsionArguments[j][3];
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
			  Components[comp].ImproperTorsionArguments[j][0]=Components[SolventIndex].ImproperTorsionArguments[j][0];
			  Components[comp].ImproperTorsionArguments[j][1]=Components[SolventIndex].ImproperTorsionArguments[j][1];
			  Components[comp].ImproperTorsionArguments[j][2]=Components[SolventIndex].ImproperTorsionArguments[j][2];
			  Components[comp].ImproperTorsionArguments[j][3]=Components[SolventIndex].ImproperTorsionArguments[j][3];
			  Components[comp].ImproperTorsionArguments[j][4]=Components[SolventIndex].ImproperTorsionArguments[j][4];
			  Components[comp].ImproperTorsionArguments[j][5]=Components[SolventIndex].ImproperTorsionArguments[j][5];
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
			  Components[comp].ImproperTorsionArguments[j][0]=Components[SolventIndex].ImproperTorsionArguments[j][0];
			  Components[comp].ImproperTorsionArguments[j][1]=Components[SolventIndex].ImproperTorsionArguments[j][1];
			  Components[comp].ImproperTorsionArguments[j][2]=Components[SolventIndex].ImproperTorsionArguments[j][2];
			  Components[comp].ImproperTorsionArguments[j][3]=Components[SolventIndex].ImproperTorsionArguments[j][3];
			  Components[comp].ImproperTorsionArguments[j][4]=Components[SolventIndex].ImproperTorsionArguments[j][4];
			  Components[comp].ImproperTorsionArguments[j][5]=Components[SolventIndex].ImproperTorsionArguments[j][5];
			  break;
			case FIXED_IMPROPER_DIHEDRAL:
			  Components[comp].ImproperTorsionArguments[j][0]=Components[SolventIndex].ImproperTorsionArguments[j][0];
			  Components[comp].NumberOfConstraintImproperTorsions++;
			  break;
			default:
			  fprintf(stderr, "Undefined Improper Torsion potential in routine 'ReadComponentDefinition' ('molecule.c')\n");
			  exit(0);
			  break;
      }
    }
  }

  // Bond/Strech cross term-data
  if(Components[comp].NumberOfBondBonds>0)
  {

    for(j=0;j<Components[comp].NumberOfBondBonds;j++)
    {
         Components[comp].BondBonds[j].A=Components[SolventIndex].BondBonds[j].A;
         Components[comp].BondBonds[j].B=Components[SolventIndex].BondBonds[j].B;
         Components[comp].BondBonds[j].C=Components[SolventIndex].BondBonds[j].C;
		
		 // Attribute Bond/stretch cross term
		 Components[comp].BondBondType[j]=Components[SolventIndex].BondBondType[j];

		 for(k=0;k<BondBondTypes[Components[comp].BondBondType[j]].nr_args;k++)
		 {
			Components[comp].BondBondArguments[j][k]=Components[SolventIndex].BondBondArguments[j][k];
		 }

		 switch(Components[comp].BondBondType[j])
		 {
			case CVFF_BOND_BOND_CROSS:
			case CFF_BOND_BOND_CROSS:
			  // p_0*(rab-p_1)*(rbc-p_2)
			  // =======================
			  // p_0/k_B [K/A^2]
			  // p_1     [A]
			  // p_2     [A]
			  Components[comp].BondBondArguments[j][0]=Components[SolventIndex].BondBondArguments[j][0];
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
    for(j=0;j<Components[comp].NumberOfBondBends;j++)
    {

         Components[comp].BondBends[j].A=Components[SolventIndex].BondBends[j].A;
         Components[comp].BondBends[j].B=Components[SolventIndex].BondBends[j].B;
         Components[comp].BondBends[j].C=Components[SolventIndex].BondBends[j].C;
          
         // Attribute Bond/Bend cross term
         Components[comp].BondBendType[j]=Components[SolventIndex].BondBendType[j];
          
		 for(k=0;k<BondBendTypes[Components[comp].BondBendType[j]].nr_args;k++)
		 {
			Components[comp].BondBendArguments[j][k]=Components[SolventIndex].BondBendArguments[j][k];
		 }

		 switch(Components[comp].BondBendType[j])
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
			  Components[comp].BondBendArguments[j][0]=Components[SolventIndex].BondBendArguments[j][0];
			  Components[comp].BondBendArguments[j][1]=Components[SolventIndex].BondBendArguments[j][1];
			  Components[comp].BondBendArguments[j][3]=Components[SolventIndex].BondBendArguments[j][3];
			  break;
			case MM3_BOND_BEND_CROSS:
			  // p_0*[(rab-p_1)+(rbc-p_2)]*(Theta-p_3)
			  // =====================================
			  // p_0     [mdyne/rad]
			  // p_1     [A]
			  // p_2     [A]
			  // p_3     [degrees]
			  Components[comp].BondBendArguments[j][0]=Components[SolventIndex].BondBendArguments[j][0];
			  Components[comp].BondBendArguments[j][3]=Components[SolventIndex].BondBendArguments[j][3];
			  break;
			case TRUNCATED_HARMONIC:
			  // (1/2)*p_0*(Theta-p_1)^2*exp(-(pow(rab,8)+pow(rbc,8))/pow(p_2,8))
			  // ================================================================
			  // p_0/k_B [K/rad^2]
			  // p_1     [degrees]
			  // p_2     [A]
			  Components[comp].BondBendArguments[j][0]=Components[SolventIndex].BondBendArguments[j][0];
			  Components[comp].BondBendArguments[j][1]=Components[SolventIndex].BondBendArguments[j][1];
			  break;
			case SCREENED_HARMONIC:
			  // (1/2)*p_0*(Theta-p_1)^2*exp(-(rab/p_2+rbc/p_3))
			  // ===============================================
			  // p_0/k_B [K/rad^2]
			  // p_1     [degrees]
			  // p_2     [A]
			  // p_3     [A]
			  Components[comp].BondBendArguments[j][0]=Components[SolventIndex].BondBendArguments[j][0];
			  Components[comp].BondBendArguments[j][1]=Components[SolventIndex].BondBendArguments[j][1];
			  break;
			case SCREENED_VESSAL:
			  // (p_0/(8.0*(Theta-PI)^2))*((p_1-PI)^2-(Theta-PI)^2)^2*exp(-(rab/p_2+rbc/p_3))
			  // ============================================================================
			  // p_0/k_B [K/rad^2]
			  // p_1     [degrees]
			  // p_2     [A]
			  // p_3     [A]
			  Components[comp].BondBendArguments[j][0]=Components[SolventIndex].BondBendArguments[j][0];
			  Components[comp].BondBendArguments[j][1]=Components[SolventIndex].BondBendArguments[j][1];
			  break;
			case TRUNCATED_VESSAL:
			  // p_0*[pow(theta,p_2)*(theta-p_1)^2*(theta+p_1-2.0*PI)^2-0.5*p_2*pow(PI,p_2-1.0)*(theta-p_1)^2*pow(PI-p_1,3)]
			  //    *exp(-(pow(rab,8)+pow(rbc,8))/pow(p_3,8))
			  // ============================================================================
			  // p_0/k_B [K/rad^(4+p_2)]
			  // p_1     [degrees]
			  // p_2     [-]
			  // p_3     [A]
			  Components[comp].BondBendArguments[j][0]=Components[SolventIndex].BondBendArguments[j][0];
			  Components[comp].BondBendArguments[j][1]=Components[SolventIndex].BondBendArguments[j][1];
			  break;
			default:
			  fprintf(stderr, "Undefined Bond-Bend potential in routine 'ReadComponentDefinition' ('molecule.c')\n");
			  exit(0);
			  break;
      }
    }
  }
	  
  // Bend/Bend cross term-data
  if(Components[comp].NumberOfBendBends>0)
  {
    for(j=0;j<Components[comp].NumberOfBendBends;j++)
    {
		 Components[comp].BendBends[j].A=Components[SolventIndex].BendBends[j].A;
		 Components[comp].BendBends[j].B=Components[SolventIndex].BendBends[j].B;
		 Components[comp].BendBends[j].C=Components[SolventIndex].BendBends[j].C;
		 Components[comp].BendBends[j].D=Components[SolventIndex].BendBends[j].D;

		 // Attribute Bend/Bend cross term.
		 Components[comp].BendBendType[j]=Components[SolventIndex].BendBendType[j];

      for(k=0;k<BendBendTypes[Components[comp].BendBendType[j]].nr_args;k++)
      {
         Components[comp].BendBendArguments[j][k]=Components[SolventIndex].BendBendArguments[j][k];
      }

      switch(Components[comp].BendBendType[j])
      {
        case CVFF_BEND_BEND_CROSS:
        case CFF_BEND_BEND_CROSS:
         // p_0*(Theta1-p_1)*(Theta2-p_2)
          // ===================================
          // p_0/k_B [K/rad^2)]
          // p_1     [degrees]
          // p_2     [degrees]
          Components[comp].BendBendArguments[j][0]=Components[SolventIndex].BendBendArguments[j][0];
          Components[comp].BendBendArguments[j][1]=Components[SolventIndex].BendBendArguments[j][1];
        case MM3_BEND_BEND_CROSS:
          // -p_0*(Theta1-p_1)*(Theta2-p_2)
          // ===================================
          // p_0     [mdyne A/rad^2]
          // p_1     [degrees]
          // p_2     [degrees]
          Components[comp].BendBendArguments[j][0]=Components[SolventIndex].BendBendArguments[j][0];
          Components[comp].BendBendArguments[j][1]=Components[SolventIndex].BendBendArguments[j][1];
          Components[comp].BendBendArguments[j][2]=Components[SolventIndex].BendBendArguments[j][2];
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
    for(j=0;j<Components[comp].NumberOfBondTorsions;j++)
    {
           Components[comp].BondTorsions[j].A=Components[SolventIndex].BondTorsions[j].A;
           Components[comp].BondTorsions[j].B=Components[SolventIndex].BondTorsions[j].B;
           Components[comp].BondTorsions[j].C=Components[SolventIndex].BondTorsions[j].C;
           Components[comp].BondTorsions[j].D=Components[SolventIndex].BondTorsions[j].D;
           
           // Attribute Bond/Torsion cross term
           Components[comp].BondTorsionType[j]=Components[SolventIndex].BondTorsionType[j];

		   for(k=0;k<BondTorsionTypes[Components[comp].BondTorsionType[j]].nr_args;k++)
		   {
				Components[comp].BondTorsionArguments[j][k]=Components[SolventIndex].BondTorsionArguments[j][k];
		   }

      switch(Components[comp].BondTorsionType[j])
      {
        case MM3_BOND_TORSION_CROSS:
          // (1/2)p_0(r-p_3)(1+cos(phi))+(1/2)p_1(r-p_3)(1+cos(2phi))+(1/2)p_2(r-p_3)(1+cos(3phi))
          // =====================================================================================
          // p_0     [kcal/A mole]
          // p_1     [kcal/A mole]
          // p_2     [kcal/A mole]
          // p_3     [A]
          Components[comp].BondTorsionArguments[j][0]=Components[SolventIndex].BondTorsionArguments[j][0];
          Components[comp].BondTorsionArguments[j][1]=Components[SolventIndex].BondTorsionArguments[j][1];
          Components[comp].BondTorsionArguments[j][2]=Components[SolventIndex].BondTorsionArguments[j][2];
          break;
        default:
          fprintf(stderr, "Undefined Bond-Torsion potential in routine 'ReadComponentDefinition' ('molecule.c')\n");
          exit(0);
          break;
      }
    }
  }

  // Bend/Torsion cross term-data
  if(Components[comp].NumberOfBendTorsions>0)
  {
    for(j=0;j<Components[comp].NumberOfBendTorsions;j++)
    {
           Components[comp].BendTorsions[j].A=Components[SolventIndex].BendTorsions[j].A;
           Components[comp].BendTorsions[j].B=Components[SolventIndex].BendTorsions[j].B;
           Components[comp].BendTorsions[j].C=Components[SolventIndex].BendTorsions[j].C;
           Components[comp].BendTorsions[j].D=Components[SolventIndex].BendTorsions[j].D;
           
           // Attribute Bend/Torsion cross term
           Components[comp].BendTorsionType[j]=Components[SolventIndex].BendTorsionType[j];

      for(k=0;k<BendTorsionTypes[Components[comp].BendTorsionType[j]].nr_args;k++)
      {
        Components[comp].BendTorsionArguments[j][k]=Components[SolventIndex].BendTorsionArguments[j][k];
      }

      switch(Components[comp].BendTorsionType[j])
      {
        case SMOOTHED_DIHEDRAL:
          // S(Theta1)*[p_0(1+cos(p_1*Phi-p_2)]*S(Theta2)
          // ======================================================================================
          // p_0/k_B [K/rad^2]
          // p_1     [-]
          // p_2     [degrees]
          Components[comp].BendTorsionArguments[j][0]=Components[SolventIndex].BendTorsionArguments[j][0];
          Components[comp].BendTorsionArguments[j][2]=Components[SolventIndex].BendTorsionArguments[j][2];
          break;
        case SMOOTHED_THREE_COSINE_DIHEDRAL:
          // S(Theta1)*[(1/2)*(1+cos(Phi))+(1/2)*p_1*(1-cos(2*Phi))+(1/2)*(1+cos(3*Phi))]*S(Theta2)
          // ======================================================================================
          // p_0/k_B [K]
          // p_1/k_B [K]
          // p_2/k_B [K]
          Components[comp].BendTorsionArguments[j][0]=Components[SolventIndex].BendTorsionArguments[j][0];
          Components[comp].BendTorsionArguments[j][1]=Components[SolventIndex].BendTorsionArguments[j][1];
          Components[comp].BendTorsionArguments[j][2]=Components[SolventIndex].BendTorsionArguments[j][2];
          break;
        case NICHOLAS_DIHEDRAL:
          // S(Theta1)*[(1/2)*(1+cos(Phi))+(1/2)*p_1*(1-cos(2*Phi))+(1/2)*(1+cos(3*Phi))]*S(Theta2)
          // ======================================================================================
          // p_0/k_B [K]
          // p_1/k_B [K]
          // p_2/k_B [K]
          Components[comp].BendTorsionArguments[j][0]=Components[SolventIndex].BendTorsionArguments[j][0];
          Components[comp].BendTorsionArguments[j][1]=Components[SolventIndex].BendTorsionArguments[j][1];
          Components[comp].BendTorsionArguments[j][2]=Components[SolventIndex].BendTorsionArguments[j][2];
          break;
        case SMOOTHED_CFF_DIHEDRAL:
         // S(Theta1)*[(1-cos(Phi))+p_1*(1-cos(2*Phi))+(1-cos(3*Phi))]*S(Theta2)
          // ======================================================================================
          // p_0/k_B [K]
          // p_1/k_B [K]
          // p_2/k_B [K]
          Components[comp].BendTorsionArguments[j][0]=Components[SolventIndex].BendTorsionArguments[j][0];
          Components[comp].BendTorsionArguments[j][1]=Components[SolventIndex].BendTorsionArguments[j][1];
          Components[comp].BendTorsionArguments[j][2]=Components[SolventIndex].BendTorsionArguments[j][2];
          break;
        case SMOOTHED_CFF_DIHEDRAL2:
          Components[comp].BendTorsionArguments[j][0]=Components[SolventIndex].BendTorsionArguments[j][0];
          Components[comp].BendTorsionArguments[j][1]=Components[SolventIndex].BendTorsionArguments[j][1];
          Components[comp].BendTorsionArguments[j][2]=Components[SolventIndex].BendTorsionArguments[j][2];
          break;
        case CVFF_BEND_TORSION_CROSS:
        case CFF_BEND_TORSION_CROSS:
          // p_0*(Theta1-p_1)*(Theta2-p_2)*cos(Phi)
          // =====================================================================================
          // p_0/k_B [K/rad^3]
          // p_1     [degrees]
          // p_2     [degrees]
          Components[comp].BendTorsionArguments[j][0]=Components[SolventIndex].BendTorsionArguments[j][0];
          Components[comp].BendTorsionArguments[j][1]=Components[SolventIndex].BendTorsionArguments[j][1];
          Components[comp].BendTorsionArguments[j][2]=Components[SolventIndex].BendTorsionArguments[j][2];
          break;
        case SMOOTHED_CFF_BEND_TORSION_CROSS:
          // S(Theta1)*[p_0*(Theta1-p_1)*(Theta2-p_2)*cos(Phi)]*S(Theta2)
          // ======================================================================================
          // p_0/k_B [K/rad^3]
          // p_1     [degrees]
          // p_2     [degrees]
          Components[comp].BendTorsionArguments[j][0]=Components[SolventIndex].BendTorsionArguments[j][0];
          Components[comp].BendTorsionArguments[j][1]=Components[SolventIndex].BendTorsionArguments[j][1];
          Components[comp].BendTorsionArguments[j][2]=Components[SolventIndex].BendTorsionArguments[j][2];
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
    for(j=0;j<Components[comp].NumberOfIntraVDW;j++)
    {
        Components[comp].IntraVDW[j].A=Components[SolventIndex].IntraVDW[j].A;
        Components[comp].IntraVDW[j].B=Components[SolventIndex].IntraVDW[j].B;
		Components[comp].IntraVDWScaling[j]=Components[SolventIndex].IntraVDWScaling[j];
    }
  }  
	 
  if(Components[comp].NumberOfIntraChargeCharge>0)
  {
    for(j=0;j<Components[comp].NumberOfIntraChargeCharge;j++)
    {

        Components[comp].IntraChargeCharge[j].A=Components[SolventIndex].IntraChargeCharge[j].A;
        Components[comp].IntraChargeCharge[j].B=Components[SolventIndex].IntraChargeCharge[j].B;
        Components[comp].IntraChargeChargeScaling[j]=Components[SolventIndex].IntraChargeChargeScaling[j];
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
    for(i=0;i<Components[comp].NumberOfIntraChargeBondDipole;i++)
    {
      Components[comp].IntraChargeBondDipole[i].A=Components[SolventIndex].IntraChargeBondDipole[i].A;
      Components[comp].IntraChargeBondDipole[i].B=Components[SolventIndex].IntraChargeBondDipole[i].B;
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
    for(i=0;i<Components[comp].NumberOfIntraBondDipoleBondDipole;i++)
    {

      Components[comp].IntraBondDipoleBondDipole[i].A=Components[SolventIndex].IntraBondDipoleBondDipole[i].A;
      Components[comp].IntraBondDipoleBondDipole[i].B=Components[SolventIndex].IntraBondDipoleBondDipole[i].B;
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
  Components[comp].NumberOfConfigMoves=Components[SolventIndex].NumberOfConfigMoves;
	  
  // allocate config-moves
  Components[comp].NumberOfUnchangedAtomsConfig=(int*)calloc(Components[comp].NumberOfConfigMoves,sizeof(int));
  Components[comp].UnchangedAtomsConfig=(int**)calloc(Components[comp].NumberOfConfigMoves,sizeof(int*));
  
  for(i=0;i<Components[comp].NumberOfConfigMoves;i++)
  {
    Components[comp].UnchangedAtomsConfig[i]=(int*)calloc(temp,sizeof(int));
    Components[comp].NumberOfUnchangedAtomsConfig[i]=Components[SolventIndex].NumberOfUnchangedAtomsConfig[i];
    for(j=0;j<temp;j++)
    {
		Components[comp].UnchangedAtomsConfig[i][j]=Components[SolventIndex].UnchangedAtomsConfig[i][j];
    }
  }
  
   // read the defined identity moves
  temp=Components[SolventIndex].NumberOfIdentityConfigMoves;
  Components[comp].NumberOfIdentityConfigMoves=0;
  if(temp>1)
  {
    Components[comp].NumberOfIdentityConfigMoves=temp;
    // allocate config-moves
    Components[comp].NumberOfUnchangedAtomsIdentityConfig=(int*)calloc(Components[comp].NumberOfIdentityConfigMoves,sizeof(int));
    Components[comp].UnchangedAtomsIdentityConfig=(int**)calloc(Components[comp].NumberOfIdentityConfigMoves,sizeof(int*));

    for(i=0;i<Components[comp].NumberOfIdentityConfigMoves;i++)
    {
      Components[comp].UnchangedAtomsIdentityConfig[i]=(int*)calloc(temp,sizeof(int));
      Components[comp].NumberOfUnchangedAtomsIdentityConfig[i]=temp;
      for(j=0;j<temp;j++)
      {
        Components[comp].UnchangedAtomsIdentityConfig[i][j]=Components[SolventIndex].UnchangedAtomsIdentityConfig[i][j];
      }
    }
  }
  else
  {
    // allocate config-moves
    Components[comp].StartingBead=Components[SolventIndex].StartingBead;
    Components[comp].NumberOfIdentityConfigMoves=1;
    Components[comp].NumberOfUnchangedAtomsIdentityConfig=(int*)calloc(Components[comp].NumberOfIdentityConfigMoves,sizeof(int));
    Components[comp].UnchangedAtomsIdentityConfig=(int**)calloc(Components[comp].NumberOfIdentityConfigMoves,sizeof(int*));
    Components[comp].UnchangedAtomsIdentityConfig[0]=(int*)calloc(1,sizeof(int));
    Components[comp].NumberOfUnchangedAtomsIdentityConfig[0]=1;
    Components[comp].UnchangedAtomsIdentityConfig[0][0]=Components[SolventIndex].UnchangedAtomsIdentityConfig[0][0]; // A MODIFIER
  }

  // Attribute defaults
  Components[comp].LMCMOL=Components[SolventIndex].LMCMOL;
 
  Components[comp].TranslationMatrix.ax=Components[SolventIndex].TranslationMatrix.ax;
  Components[comp].TranslationMatrix.ay=Components[SolventIndex].TranslationMatrix.ay;
  Components[comp].TranslationMatrix.az=Components[SolventIndex].TranslationMatrix.az;
											
  Components[comp].TranslationMatrix.bx=Components[SolventIndex].TranslationMatrix.bx;
  Components[comp].TranslationMatrix.by=Components[SolventIndex].TranslationMatrix.by;
  Components[comp].TranslationMatrix.bz=Components[SolventIndex].TranslationMatrix.bz;
												
  Components[comp].TranslationMatrix.cx=Components[SolventIndex].TranslationMatrix.cx;
  Components[comp].TranslationMatrix.cy=Components[SolventIndex].TranslationMatrix.cy;
  Components[comp].TranslationMatrix.cz=Components[SolventIndex].TranslationMatrix.cz;

  Components[comp].SwapEvery=Components[SolventIndex].SwapEvery;

  for(i=0;i<Components[comp].NumberOfAtoms;i++)
  {
    for(k=0;k<NumberOfSystems;k++)
    {
      Components[comp].MaximumCBMCChangeBondLength[k][i]=Components[SolventIndex].MaximumCBMCChangeBondLength[k][i];
      Components[comp].MaximumCBMCChangeBendAngle[k][i] =Components[SolventIndex].MaximumCBMCChangeBendAngle[k][i];
      Components[comp].MaximumCBMCRotationOnCone[k][i]  =Components[SolventIndex].MaximumCBMCRotationOnCone[k][i];
      Components[comp].CBMCChangeBendAngleAttempts[k][i]=Components[SolventIndex].CBMCChangeBendAngleAttempts[k][i];
      Components[comp].CBMCChangeBendAngleAccepted[k][i]=Components[SolventIndex].CBMCChangeBendAngleAccepted[k][i];
      Components[comp].CBMCRotationOnConeAttempts[k][i] =Components[SolventIndex].CBMCRotationOnConeAttempts[k][i];
      Components[comp].CBMCRotationOnConeAccepted[k][i] =Components[SolventIndex].CBMCRotationOnConeAccepted[k][i];
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

  }// end of loop over comp
}

/*********************************************************************************************************
 * Name       | InitializeLambda    (Added by A. de Izarra)      					 	 				 *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Initialization of the interpolation parameter of the alchemical transformation.			 *	
 * Parameters | No parameters																			 *
 *********************************************************************************************************/
void InitializeLambda(void)
{
  int i;
  
  Lambda    = (REAL*)calloc(AlchReacLambda+1,sizeof(REAL));
  
  for(i=0;i<=AlchReacLambda;i++)
  {
	Lambda[i]    = ((REAL)i)/AlchReacLambda;
  }
}


void ManageFirstBead(int start)
{
  int i,type;
  REAL EnergyHostVDW,EnergyAdsorbateVDW,EnergyCationVDW;
  REAL EnergyHostChargeCharge,EnergyAdsorbateChargeCharge,EnergyCationChargeCharge;
  REAL EnergyHostChargeBondDipole,EnergyAdsorbateChargeBondDipole,EnergyCationChargeBondDipole;
  POINT posA;

  // Remark: i=0 here because it is the first bead of the transient molecule.

  // There is already a given position for the first bead.
  posA=FirstBeadPosition;
  type=Components[CurrentComponent].Type[start];

  EnergyHostVDW=EnergyAdsorbateVDW=EnergyCationVDW=0.0;
  EnergyHostChargeCharge=EnergyAdsorbateChargeCharge=EnergyCationChargeCharge=0.0;
  EnergyHostChargeBondDipole=EnergyAdsorbateChargeBondDipole=EnergyCationChargeBondDipole=0.0;

  // calculate energies
  EnergyHostVDW=CalculateFrameworkVDWEnergyAtPosition(posA,type,CFVDWScaling[start]);
  CalculateFrameworkChargeEnergyAtPosition(posA,type,&EnergyHostChargeCharge,&EnergyHostChargeBondDipole,CFChargeScaling[start]);

  // compute VDW energy with adsorbates if no omit of adsorbate-adsorbate or the current molecule is a cation
  if(Components[CurrentComponent].ExtraFrameworkMolecule||(!OmitAdsorbateAdsorbateVDWInteractions))
    EnergyAdsorbateVDW=CalculateInterVDWEnergyAdsorbateAtPosition(posA,type,CurrentAdsorbateMolecule,CFVDWScaling[start]);

  // compute Coulomb energy with adsorbates if no omit of adsorbate-adsorbate or the current molecule is a cation
  if(Components[CurrentComponent].ExtraFrameworkMolecule||(!OmitAdsorbateAdsorbateCoulombInteractions))
    CalculateInterChargeEnergyAdsorbateAtPosition(posA,type,&EnergyAdsorbateChargeCharge,&EnergyAdsorbateChargeBondDipole,CurrentAdsorbateMolecule,CFChargeScaling[start] * PseudoAtoms[type].Charge1);

  // compute VDW energy with cations if no omit of cation-cation or the current molecule is an adsorbate
  if((!Components[CurrentComponent].ExtraFrameworkMolecule)||(!OmitCationCationVDWInteractions))
    EnergyCationVDW=CalculateInterVDWEnergyCationAtPosition(posA,type,CurrentCationMolecule,CFVDWScaling[start]);

  // compute Coulomb energy with cations if no omit of cation-cation or the current molecule is an adsorbate
  if((!Components[CurrentComponent].ExtraFrameworkMolecule)||(!OmitCationCationCoulombInteractions))
    CalculateInterChargeEnergyCationAtPosition(posA,type,&EnergyCationChargeCharge,&EnergyCationChargeBondDipole,CurrentCationMolecule,CFChargeScaling[start] * PseudoAtoms[type].Charge1);

    EnergyHostVDWFirstBead=EnergyHostVDW;
    EnergyAdsorbateVDWFirstBead=EnergyAdsorbateVDW;

    EnergyCationVDWFirstBead=EnergyCationVDW;

    EnergyHostChargeChargeFirstBead=EnergyHostChargeCharge;
    EnergyAdsorbateChargeChargeFirstBead=EnergyAdsorbateChargeCharge;
    EnergyCationChargeChargeFirstBead=EnergyCationChargeCharge;

    EnergyHostChargeBondDipoleFirstBead=EnergyHostChargeBondDipole;
    EnergyAdsorbateChargeBondDipoleFirstBead=EnergyAdsorbateChargeBondDipole;
    EnergyCationChargeBondDipoleFirstBead=EnergyCationChargeBondDipole;

    EnergyHostBondDipoleBondDipoleFirstBead=0.0;
    EnergyAdsorbateBondDipoleBondDipoleFirstBead=0.0;
    EnergyCationBondDipoleBondDipoleFirstBead=0.0;

}


void ManageRemainingBeads(int NumberOldComponent)
{
  int j,k,ip,iu;

  SetConnectivityMatrix();

	if(NumberOldComponent==1) // water->ions
	{
		SetGrowingStatus();

		Interactions();

		for(k=0;k<Components[CurrentComponent].NumberOfAtoms;k++)
		{
		  // The new position to be inserted are in fact the old one (we do this to match raspa function.
		  TrialPositions[0][k]=NewPosition[CurrentSystem][k];
		}

		ComputeExternalEnergies();
			
		UBondNew[CurrentSystem]+=UBondTrial[0];
		UUreyBradleyNew[CurrentSystem]+=UUreyBradleyTrial[0];
		UBendNew[CurrentSystem]+=UBendTrial[0];
		UBendBendNew[CurrentSystem]+=UBendBendTrial[0];
		UInversionBendNew[CurrentSystem]+=UInversionBendTrial[0];
		UTorsionNew[CurrentSystem]+=UTorsionTrial[0];
		UImproperTorsionNew[CurrentSystem]+=UImproperTorsionTrial[0];
		UBondBondNew[CurrentSystem]+=UBondBondTrial[0];
		UBondBendNew[CurrentSystem]+=UBondBendTrial[0];
		UBondTorsionNew[CurrentSystem]+=UBondTorsionTrial[0];
		UBendTorsionNew[CurrentSystem]+=UBendTorsionTrial[0];
		UIntraVDWNew[CurrentSystem]+=UIntraVDWTrial[0];
		UIntraChargeChargeNew[CurrentSystem]+=UIntraChargeChargeTrial[0];
		UIntraChargeBondDipoleNew[CurrentSystem]+=UIntraChargeBondDipoleTrial[0];
		UIntraBondDipoleBondDipoleNew[CurrentSystem]+=UIntraBondDipoleBondDipoleTrial[0];

		UHostVDWNew[CurrentSystem]+=UHostVDWTrial[0];
		UAdsorbateVDWNew[CurrentSystem]+=UAdsorbateVDWTrial[0];
		UCationVDWNew[CurrentSystem]+=UCationVDWTrial[0];
		UHostChargeChargeNew[CurrentSystem]+=UHostChargeChargeTrial[0];
		UAdsorbateChargeChargeNew[CurrentSystem]+=UAdsorbateChargeChargeTrial[0];
		UCationChargeChargeNew[CurrentSystem]+=UCationChargeChargeTrial[0];
		UHostChargeBondDipoleNew[CurrentSystem]+=UHostChargeBondDipoleTrial[0];
		UAdsorbateChargeBondDipoleNew[CurrentSystem]+=UAdsorbateChargeBondDipoleTrial[0];
		UCationChargeBondDipoleNew[CurrentSystem]+=UCationChargeBondDipoleTrial[0];
		UHostBondDipoleBondDipoleNew[CurrentSystem]+=UHostBondDipoleBondDipoleTrial[0];
		UAdsorbateBondDipoleBondDipoleNew[CurrentSystem]+=UAdsorbateBondDipoleBondDipoleTrial[0];
		UCationBondDipoleBondDipoleNew[CurrentSystem]+=UCationBondDipoleBondDipoleTrial[0];
	}
	else // ions->water
	{
		SetGrowingStatus();

		Interactions();

		// The position of the first bead is already known. (it is the position of the former ion).

		TrialPositions[0][0]=NewPosition[CurrentSystem][0];
		  int comp = CurrentComponent;

		GenerateTrialOrientationsSimpleSphere(FALSE);

		ComputeExternalEnergies();

		UBondNew[CurrentSystem]+=UBondTrial[0];
		UBendNew[CurrentSystem]+=UBendTrial[0];
		UBendBendNew[CurrentSystem]+=UBendBendTrial[0];
		UInversionBendNew[CurrentSystem]+=UInversionBendTrial[0];
		UUreyBradleyNew[CurrentSystem]+=UUreyBradleyTrial[0];
		UTorsionNew[CurrentSystem]+=UTorsionTrial[0];
		UImproperTorsionNew[CurrentSystem]+=UImproperTorsionTrial[0];
		UBondBondNew[CurrentSystem]+=UBondBondTrial[0];
		UBondBendNew[CurrentSystem]+=UBondBendTrial[0];
		UBondTorsionNew[CurrentSystem]+=UBondTorsionTrial[0];
		UBendTorsionNew[CurrentSystem]+=UBendTorsionTrial[0];
		UIntraVDWNew[CurrentSystem]+=UIntraVDWTrial[0];
		UIntraChargeChargeNew[CurrentSystem]+=UIntraChargeChargeTrial[0];
		UIntraChargeBondDipoleNew[CurrentSystem]+=UIntraChargeBondDipoleTrial[0];
		UIntraBondDipoleBondDipoleNew[CurrentSystem]+=UIntraBondDipoleBondDipoleTrial[0];

		UCationVDWNew[CurrentSystem]+=UCationVDWTrial[0];
		UAdsorbateVDWNew[CurrentSystem]+=UAdsorbateVDWTrial[0];
		UHostVDWNew[CurrentSystem]+=UHostVDWTrial[0];
		UHostChargeChargeNew[CurrentSystem]+=UHostChargeChargeTrial[0];
		UAdsorbateChargeChargeNew[CurrentSystem]+=UAdsorbateChargeChargeTrial[0];
		UCationChargeChargeNew[CurrentSystem]+=UCationChargeChargeTrial[0];
		UHostChargeBondDipoleNew[CurrentSystem]+=UHostChargeBondDipoleTrial[0];
		UAdsorbateChargeBondDipoleNew[CurrentSystem]+=UAdsorbateChargeBondDipoleTrial[0];
		UCationChargeBondDipoleNew[CurrentSystem]+=UCationChargeBondDipoleTrial[0];
		UHostBondDipoleBondDipoleNew[CurrentSystem]+=UHostBondDipoleBondDipoleTrial[0];
		UAdsorbateBondDipoleBondDipoleNew[CurrentSystem]+=UAdsorbateBondDipoleBondDipoleTrial[0];
		UCationBondDipoleBondDipoleNew[CurrentSystem]+=UCationBondDipoleBondDipoleTrial[0];
	
		// update the selected trial-position to the 'NewPosition'
		for(j=NumberOfBeadsAlreadyPlaced;j<Components[CurrentComponent].NumberOfAtoms;j++)
			NewPosition[CurrentSystem][j]=TrialPositions[0][j];
	}
	for(j=0;j<Components[CurrentComponent].Groups[CurrentGroup].NumberOfGroupAtoms;j++)
    {
          int atom_nr=Components[CurrentComponent].Groups[CurrentGroup].Atoms[j];
    } 
}

/*********************************************************************************************************
 * Name       | RemoveExtraPseudoAtoms    (Added by A. de Izarra)      					 	 		     *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Once alchemical transformation is finished, remove extra pseudo atoms of the molecules	 *
 * 			  | that undergo the alchemical transformation.												 *	
 * Parameters | No parameters																			 *
 *********************************************************************************************************/
void RemoveExtraPseudoAtoms(void)
{
	
	int i;
	
	PseudoAtoms=(PSEUDO_ATOM*)realloc(PseudoAtoms,(InitialPseudoAtoms)*sizeof(PSEUDO_ATOM));

	NumberOfPseudoAtomsTypeNew=(int*)realloc(NumberOfPseudoAtomsTypeNew,(InitialPseudoAtoms)*sizeof(int));
    NumberOfPseudoAtomsTypeOld=(int*)realloc(NumberOfPseudoAtomsTypeOld,(InitialPseudoAtoms)*sizeof(int));
    MapPseudoAtom=(int*)realloc(MapPseudoAtom,(InitialPseudoAtoms)*sizeof(int));
    

    for(i=0;i<NumberOfSystems;i++)
	{
		 NumberOfPseudoAtomsCount[i]=(int*)realloc(NumberOfPseudoAtomsCount[i],(InitialPseudoAtoms)*sizeof(int));
		 NumberOfPseudoAtomsType[i]=(int*)realloc(NumberOfPseudoAtomsType[i],(InitialPseudoAtoms)*sizeof(int));
		 NumberOfFractionalPseudoAtomsType[i]=(int*)realloc(NumberOfFractionalPseudoAtomsType[i],(InitialPseudoAtoms)*sizeof(int));	
	}
	
	// Reinitialize the size of Numberof pseudo atoms.
	NumberOfPseudoAtoms=InitialPseudoAtoms;
	
	
}

/*********************************************************************************************************
 * Name       | DeallocateTransientComponentMemory    (Added by A. de Izarra)      					     *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Deallocate memory for the component of the molecules that undergo the alchemical         *
 * 			  | transformation.																			 *
 * Parameters | No parameters																			 *
 *********************************************************************************************************/
void DeallocateTransientComponentMemory(void)
{
	Components=(COMPONENT*)realloc(Components,(NumberOfComponents)*sizeof(COMPONENT));
}

/*********************************************************************************************************
 * Name       | DeallocateMemoryParameterTab    (Added by A. de Izarra)      					     	 *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Deallocate memory for the paramter tab that store the vdw interaction for the 			 *
 * 			  | pseudoatoms of the molecules that undergo the alchemical transformation.				 *
 * Parameters | No parameters																			 *
 *********************************************************************************************************/
void DeallocateMemoryParameterTab(void)
{
	int i;
	 //-----------------------------------------------
	 // reallocate extra memory in PotentialParms, updated at each alchemical transformation step.
	 PotentialParms=(REAL(**)[MAX_NUMBER_OF_POTENTIAL_ARGUMENTS])realloc(PotentialParms,
	 (InitialPseudoAtoms)*sizeof(REAL(*)[MAX_NUMBER_OF_POTENTIAL_ARGUMENTS]));

	 for(i = 0; i < InitialPseudoAtoms; i++)
	 {
		  PotentialParms[i] = (REAL(*)[MAX_NUMBER_OF_POTENTIAL_ARGUMENTS])realloc(PotentialParms[i],(InitialPseudoAtoms)*sizeof(REAL[MAX_NUMBER_OF_POTENTIAL_ARGUMENTS]));
	 }

	 for(i = InitialPseudoAtoms; i < InitialPseudoAtoms+NumberExtraPseudoAtoms; i++)
	 {
		  //free(PotentialParms[i]); 
	 }

     
	 //-----------------------------------------------
	 // reallocate extra memory in Potentialtyoe
	 PotentialType=(int**)realloc(PotentialType,(InitialPseudoAtoms)*(sizeof(int*)));
		 
	 for(i = 0; i < InitialPseudoAtoms; i++)
	 {	 
	
		PotentialType[i] = (int*)realloc(PotentialType[i],(InitialPseudoAtoms)*sizeof(int));
	 }
		 
	 for(i = InitialPseudoAtoms; i < InitialPseudoAtoms+NumberExtraPseudoAtoms; i++)
	 {
		//free(PotentialType[i]);
	 }

		
	
	//-----------------------------------------------
	// Allocate extra memory for tail correction.
	TailCorrection=(int**)realloc(TailCorrection,(InitialPseudoAtoms)*(sizeof(int*)));
	 
	 for(i = 0; i < InitialPseudoAtoms; i++)
	 {
		  TailCorrection[i] = (int*)realloc(TailCorrection[i],(InitialPseudoAtoms)*sizeof(int));
	 }
	 
	 for(i = InitialPseudoAtoms; i < InitialPseudoAtoms+NumberExtraPseudoAtoms; i++)
	 {
		 // free(TailCorrection[i]);
	 }

	
	
	//-----------------------------------------------
	// Allocate extra memory for shift potential
	ShiftPotential=(int**)realloc(ShiftPotential,(InitialPseudoAtoms)*(sizeof(int*)));
	 
	 for(i = 0; i < InitialPseudoAtoms; i++)
	 {
		  ShiftPotential[i] = (int*)realloc(ShiftPotential[i],(InitialPseudoAtoms)*sizeof(int));
	 }
	 
	 for(i = InitialPseudoAtoms; i < InitialPseudoAtoms+NumberExtraPseudoAtoms; i++)
	 {
		  //free(ShiftPotential[i]);
	 }	 
}

/*********************************************************************************************************
 * Name       | DeallocateChosenMoitiesCoordinates    (Added by A. de Izarra)      					     *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Deallocate the vectors that store the initial position of the molecules that			 *
 * 			  | undergo the alchemical transformation.													 *
 * Parameters | No parameters																			 *
 *********************************************************************************************************/
void DeallocateChosenMoitiesCoordinates()
{
	int i,j;
	
	for(i=0;i<NumberTransientMoities[CurrentAlchemicalReaction];i++)
	{
		free(AdsorbatesReferenceChosen[i].Atoms);
		free(AdsorbatesReferenceChosen[i].Groups);
	}
	
	free(AdsorbatesReferenceChosen);
}
