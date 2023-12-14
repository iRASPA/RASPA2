// Alchemical_transformation.h Added by Ambroise de Izarra
/*  
 * Alchemical transformation is added to build functions that perform
 * an alchemical transformation at a fixed number of particule in the system (osmotic ensemble)
 */ 
#ifndef ALCHEMICAL_H
#define ALCHEMICAL_H

#include "molecule.h"

extern int IsAlchemicalorOsmoticStep;

extern int InitialPseudoAtoms;
extern int NumberExtraPseudoAtoms;
extern int NumberExtraComponents;
extern int CurrentAlchemicalReaction;
extern int IndexExtraComponent[2];
extern int TransientMoleculeNbAtoms;

// Save in memory initial coordinates of chosen moities 
extern ADSORBATE_MOLECULE *AdsorbatesReferenceChosen;
extern VECTOR **SolventBodyfixedPositions;

// Save in memory initial NumberofPseudoType tables.
extern int **NumberOfPseudoAtomsReferenceType;
extern int *NumberOfPseudoAtomsReferenceTypeNew;
extern int *NumberOfPseudoAtomsReferenceTypeOld;

// Save in memory the degree of freedom
extern int DegreesOfFreedomReferenceAdsorbates;
extern int DegreesOfFreedomReferenceTranslation;
extern int DegreesOfFreedomReferenceTranslationalAdsorbates;
extern int DegreesOfFreedomReference;
                              
extern int DegreesOfFreedomReferenceRotation;
extern int DegreesOfFreedomReferenceAdsorbates;
extern int DegreesOfFreedomReferenceRotationalAdsorbates;

// Declare the chemical potential os the osmostat.
extern REAL *AlchemicalWorkStore;
extern int SizeAlchemicalWorkStore;

// for the NVE move during the Alchemical transformation.
extern REAL *HybridNVEAlchDrift;
extern REAL *HybridNVEAlchDriftCount;

extern REAL *HybridNVEAlchStartTemperature;
extern REAL *HybridNVEAlchStartTranslationalTemperature;
extern REAL *HybridNVEAlchStartRotationalTemperature;
extern REAL *HybridNVEAlchStartTemperatureFramework;
extern REAL *HybridNVEAlchStartTemperatureAdsorbate;
extern REAL *HybridNVEAlchStartTemperatureCation;
extern REAL *HybridNVEAlchStartTemperatureCount;
extern REAL *HybridNVEAlchStartTemperatureTranslationCount;
extern REAL *HybridNVEAlchStartTemperatureRotationCount;
extern REAL *HybridNVEAlchStartTemperatureFrameworkCount;
extern REAL *HybridNVEAlchStartTemperatureAdsorbateCount;
extern REAL *HybridNVEAlchStartTemperatureCationCount;

extern REAL *HybridNVEAlchEndTemperature;
extern REAL *HybridNVEAlchEndTranslationalTemperature;
extern REAL *HybridNVEAlchEndRotationalTemperature;
extern REAL *HybridNVEAlchEndTemperatureFramework;
extern REAL *HybridNVEAlchEndTemperatureAdsorbate;
extern REAL *HybridNVEAlchEndTemperatureCation;
extern REAL *HybridNVEAlchEndTemperatureCount;
extern REAL *HybridNVEAlchEndTemperatureTranslationCount;
extern REAL *HybridNVEAlchEndTemperatureRotationCount;
extern REAL *HybridNVEAlchEndTemperatureFrameworkCount;
extern REAL *HybridNVEAlchEndTemperatureAdsorbateCount;
extern REAL *HybridNVEAlchEndTemperatureCationCount;

// File to write Alchemical work thorugh MC.
extern FILE *OutputOsmostatFilePtr;

// Manage storage of the Alchemicalwork and write Alchemical work
void InitializeStoreAlchemicalWork(void);
void IncreaseSizeAlchemicalWork(void);
void InitializeFileAlchemicalWork(int Size);
void UpdateFileAlchemicalWork(int index, FILE *FilePtr);

// Initialize statistics for NVE-MD during alchemical transformation.
void InitializeNVEAlchStatistics(void);

// Index and constant for managing transient components.
void InitializeIndexManagingTransientMoities(void);

// Prepare addition of transient
void MakeInitialTransient(int NumberOldComponent);
void GrowTransient(int Iicode);
void ReallocateMemoryParameterTab(void);

// Management of charge during alchemical transformation.
void InitializeVectorCharge(void);
void InitializeChargeTransientMoities(int NumberOldComponent);
void UpdateChargeInterpolationAlchemicalTransformation(int NumberOldComponent, REAL Lambda);

// Management of vdw para during alchemical transformation.
void InitializeVectorforMixingRule(void);
void InitializeVDWTransientMoities(int NumberOldComponent);
void UpdateMixingRuleVDWInterpolationAlchemicalTransformation(int NumberOldComponent, REAL Lambda);

// Preparation of the initial system: store coordinates of Chosen moities and delete chosen moities from adsorbates + Update the energy. 
void StoreChosenMoitiesCoordinates(int NumberOldComponent);
void DeleteChosenMoities(void);

// If step accepted, transfert the moities into the "regular" components.
void SwitchMoietiestoRegularComponents(int NumberNewComponent);

// Initialize vectors, memory
void InitializeMassTransientMoities(int NumberOldComponent);
void EndMassTransientMoities(int NumberOldComponent);
void UpdateInertiaTensorGroups(int comp);
void AllocateTransientComponentMemory(void);
void AddExtraPseudoAtoms(void);

// Initialize lambda for interpolation of non-bonded parameter.
void InitializeLambda(void);

// Deallocate memory at the end of each step.
void RemoveExtraPseudoAtoms();
void DeallocateTransientComponentMemory();
void DeallocateMemoryParameterTab();
void DeallocateChosenMoitiesCoordinates();


#endif

