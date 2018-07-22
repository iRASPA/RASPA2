/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2017 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'framework.c' is part of RASPA-2.0

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
#include <stdbool.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <errno.h>
#include <pwd.h>
#include "constants.h"
#include "molecule.h"
#include "framework.h"
#include "framework_energy.h"
#include "simulation.h"
#include "ewald.h"
#include "potentials.h"
#include "utils.h"
#include "input.h"
#include "output.h"
#include "cbmc.h"
#include "mc_moves.h"
#include "grids.h"
#include "spacegroup.h"
#include "thermo_baro_stats.h"
#include "scattering_factors.h"
#include "internal_force.h"
#include "minimization.h"
#include "rigid.h"

extern bool STREAM;
extern char *INPUT_CRYSTAL;

CRYSTALLOGRAPHIC_STATISTICS *crystallographic_stats;

int CorrectNetChargeOnPseudoAtom;
REAL CorrectNetChargeOnPseudoAtomValue;

VECTOR *UnitCellSize;
INT_VECTOR3 *NumberOfUnitCells;
REAL_MATRIX3x3 *UnitCellBox;
REAL_MATRIX3x3 *InverseUnitCellBox;
int *Lowenstein;

int AsymmetricIons;
int SpaceGroupIons;
REAL CutOffIons;

FRAMEWORK_COMPONENT *Framework;                            // list of frameworks
int CurrentFramework;                                      // global variable: the current framework

int RemoveBondNeighboursFromLongRangeInteraction;          // remove 1-2 interactions only between atoms in a defined bond
int RemoveBendNeighboursFromLongRangeInteraction;          // remove 1-2 interactions only between atoms in a defined bend
int RemoveTorsionNeighboursFromLongRangeInteraction;       // remove 1-3 interactions only between atoms in a defined torsion

int Remove12NeighboursFromVDWInteraction;                  // remove 1-2 neighbor interaction from VDW interaction
int Remove13NeighboursFromVDWInteraction;                  // remove 1-3 neighbor interaction from VDW interaction
int Remove14NeighboursFromVDWInteraction;                  // remove 1-4 neighbor interaction from VDW interaction

int Remove12NeighboursFromChargeChargeInteraction;         // remove 1-2 neighbor interaction from charge-charge interaction
int Remove13NeighboursFromChargeChargeInteraction;         // remove 1-3 neighbor interaction from charge-charge interaction
int Remove14NeighboursFromChargeChargeInteraction;         // remove 1-4 neighbor interaction from charge-charge interaction

int Remove11NeighboursFromChargeBondDipoleInteraction;     // remove 1-1 neighbor interaction from charge-bonddipole interaction
int Remove12NeighboursFromChargeBondDipoleInteraction;     // remove 1-2 neighbor interaction from charge-bonddipole interaction
int Remove13NeighboursFromChargeBondDipoleInteraction;     // remove 1-3 neighbor interaction from charge-bonddipole interaction
int Remove14NeighboursFromChargeBondDipoleInteraction;     // remove 1-4 neighbor interaction from charge-bonddipole interaction

int Remove12NeighboursFromBondDipoleBondDipoleInteraction; // remove 1-2 neighbor interaction from bonddipole-bonddipole interaction
int Remove13NeighboursFromBondDipoleBondDipoleInteraction; // remove 1-3 neighbor interaction from bonddipole-bonddipole interaction
int Remove14NeighboursFromBondDipoleBondDipoleInteraction; // remove 1-4 neighbor interaction from bonddipole-bonddipole interaction

int InternalFrameworkLennardJonesInteractions;
int ImproperTorsionScanType;

int NumberOfSingleSubstitutionRules;
char (*SubstitutionSingleFrameworkAtom)[3][256];
int (*SubstitutionSingleFrameworkAtomTypes)[3];

int NumberOfRandomSubstitutions;
int NumberOfSubstitutionRules;
char (*SubstitutionFrameworkAtoms)[3][256];
int (*SubstitutionFrameworkAtomTypes)[3];

int NumberOfSubstitutions;
int (*ListOfAtomSubstitutions)[3];

int NumberOfModificationRules;
int *ModificationRuleType;
char (*ModifyFrameworkAtoms)[10][256];
int (*ModifyFrameworkAtomTypes)[10];

int NumberOfForbiddenConnectivityRules;
char (*ForbiddenConnectivityAtoms)[3][256];
int (*ForbiddenConnectivityTypes)[3];

void AllocateFrameworkCellList(void)
{
  int i;

  if((Framework[CurrentSystem].FrameworkModel!=NONE)&&(UseCellLists[CurrentSystem]))
  {
    Framework[CurrentSystem].CellListHead=(int**)calloc(Framework[CurrentSystem].NumberOfFrameworks,sizeof(int*));
    Framework[CurrentSystem].CellList=(int**)calloc(Framework[CurrentSystem].NumberOfFrameworks,sizeof(int*));
    for(i=0;i<Framework[CurrentSystem].NumberOfFrameworks;i++)
    {
      Framework[CurrentSystem].CellListHead[i]=(int*)calloc(NumberOfCellLists[CurrentSystem],sizeof(int));
      Framework[CurrentSystem].CellList[i]=(int*)calloc(Framework[CurrentSystem].NumberOfAtoms[i],sizeof(int));
    }
  }
}

void MakeFrameworkCellList(void)
{
  int i,f1,icell,M;
  VECTOR pos,s;

  if((Framework[CurrentSystem].FrameworkModel!=NONE)&&(UseCellLists[CurrentSystem]))
  {
    M=NumberOfCellLists[CurrentSystem];

    for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
    {
      for(i=0;i<NumberOfCellLists[CurrentSystem];i++)
        Framework[CurrentSystem].CellListHead[f1][i]=-1;

      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
        Framework[CurrentSystem].CellList[f1][i]=-1;

      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
      {
        pos=Framework[CurrentSystem].Atoms[f1][i].Position;

        // convert from xyz to abc
        s.x=InverseBox[CurrentSystem].ax*pos.x+InverseBox[CurrentSystem].bx*pos.y+InverseBox[CurrentSystem].cx*pos.z;
        s.y=InverseBox[CurrentSystem].ay*pos.x+InverseBox[CurrentSystem].by*pos.y+InverseBox[CurrentSystem].cy*pos.z;
        s.z=InverseBox[CurrentSystem].az*pos.x+InverseBox[CurrentSystem].bz*pos.y+InverseBox[CurrentSystem].cz*pos.z;

        // apply boundary condition
        s.x-=(REAL)NINT(s.x);
        s.y-=(REAL)NINT(s.y);
        s.z-=(REAL)NINT(s.z);

        // s between 0 and 1
        s.x+=0.5;
        s.y+=0.5;
        s.z+=0.5;

        // compute the corresponding cell-id
        icell=(int)(s.x*NumberOfCellListCells[CurrentSystem].x)+
              ((int)(s.y*NumberOfCellListCells[CurrentSystem].y))*NumberOfCellListCells[CurrentSystem].x+
              ((int)(s.z*NumberOfCellListCells[CurrentSystem].z))*NumberOfCellListCells[CurrentSystem].x*NumberOfCellListCells[CurrentSystem].y;

        Framework[CurrentSystem].CellList[f1][i]=Framework[CurrentSystem].CellListHead[f1][icell];
        Framework[CurrentSystem].CellListHead[f1][icell]=i;
      }
    }
  }
}



void CheckFrameworkCharges(void)
{
  int i,j,f1;
  int type;
  REAL charge,total_charge,total_count;

  for(i=0;i<NumberOfPseudoAtoms;i++)
  {
    if(PseudoAtoms[i].IsPolarizable)
      PseudoAtoms[i].HasCharges=TRUE;
    else
      PseudoAtoms[i].HasCharges=FALSE;

    // if using charges from CIF-files, then average the pseudo-atom charge over all types
    // for computation, the charge per atom is always used (not the charge of the type)
    if(PseudoAtoms[i].FrameworkAtom)
    {
      total_count=0.0;
      total_charge=0.0;
      for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
      {
        for(j=0;j<Framework[CurrentSystem].NumberOfAtoms[f1];j++)
        {
          charge=Framework[CurrentSystem].Atoms[f1][j].Charge;

          // if one atom of that type has a charge then the pseudo-atom type has charge
          if(fabs(charge)>1e-10)
            PseudoAtoms[i].HasCharges=TRUE;

          type=Framework[CurrentSystem].Atoms[f1][j].Type;
          if(type==i)
          {
            total_charge+=charge;
            total_count+=1.0;
          }
        }
      }
      PseudoAtoms[i].Charge1=total_charge/total_count;
    }
    else
    {
      if(fabs(PseudoAtoms[i].Charge1)>1e-10)
         PseudoAtoms[i].HasCharges=TRUE;
    }
  }
}

void AddAsymmetricAtom(FRAMEWORK_ASYMMETRIC_ATOM atom)
{
  int nr_atoms;

  nr_atoms=Framework[CurrentSystem].NumberOfAsymmetricAtoms[CurrentFramework];
  Framework[CurrentSystem].AtomsAsymmetric[CurrentFramework]=(FRAMEWORK_ASYMMETRIC_ATOM*)realloc(Framework[CurrentSystem].AtomsAsymmetric[CurrentFramework],
              (nr_atoms+1)*sizeof(FRAMEWORK_ASYMMETRIC_ATOM));
  Framework[CurrentSystem].AtomsAsymmetric[CurrentFramework][nr_atoms]=atom;
  Framework[CurrentSystem].NumberOfAsymmetricAtoms[CurrentFramework]++;
}

void RemoveQuotesAroundString(char *string)
{
  if(string[0]=='\'')
  {
    memmove(&string[0],&string[1],strlen(string));
    if(string[strlen(string)-1]=='\'') string[strlen(string)-1]='\0';
  }

  if(string[0]=='\"')
  {
    memmove(&string[0],&string[1],strlen(string));
    if(string[strlen(string)-1]=='\"') string[strlen(string)-1]='\0';
  }
}


/*********************************************************************************************************
 * Name       | ReadFrameworkDefinitionCIF                                                               *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | reads a CIF-file and expands the structure to 'P1' (no symmetry)                         *
 * Parameters | -                                                                                        *
 * Note       | The CIF specification can be found online: http://www.iucr.org/resources/cif/spec        *
 *            | and S.R. Hall, F.H. Allen, and I.D. Brown, Acta Cryst. A47, 655-685, 1991.               *
 *            | The definition is based on the free format STAR procedure. Not many CIF-file closely     *
 *            | follow this format striclty. Therefore the reading of CIF-files can be challenging.      *
 *********************************************************************************************************/

// 'keyword': first string a line with spaces in front removed and separated to the second string by a space.
// 'arguments': the part following a 'keyword',

// The routine implemented here only parses a subset of the CIF-file (only the needed data is processed).
// There are two groups of data, 1) single lines, and 2) looped data.
//
// The single lines are e.g. '_cell_length_a', '_cell_length_b', '_cell_length_c', '_cell_angle_alpha',
// '_cell_angle_beta', '_cell_angle_gamma', '_symmetry_space_group_name_Hall', '_symmetry_space_group_name_H-M',
// and '_symmetry_Int_Tables_number'. These can be parsed without difficulty. This input routine requires that
// the argument follow these 'data names' is on the same line.
//
// The looped data are blocks of data items that are repeated in a list. The first command of such a list is "loop_", eg.
//
// loop_
// _atom_site_label
// _atom_site_type_symbol
// _atom_site_fract_x
// _atom_site_fract_y
// _atom_site_fract_z
// C1 C 0.3671(3) 0.1045(3) 0.2279(5)
// H1 H 0.3756 0.0732 0.2106 0.131
// C2 C 0.3282(3) 0.1718(3) 0.2500
//
// However, we never parse the loop-command here. Instead we check for whether the first keyword starts with
// '_atom_type_' or '_atom_site_' or '_symmetry_equiv_pos'. This is all the data we require for simulations.
// The block is parsed as follows. Take for example '_atom_site_'. We keep on reading keywords until it does not
// start with '_atom_site_' anymore. For each keyword we keep track of whether it is the first, second, third, etc.
// This is because the order of the keyword is also the order in which the data is later provided. After this
// part we know how many elements there should be read in and what each element means. The routine here assumes the
// elements per line should correspond to the amount of "_atom_site" commands. Also unknown commands are allowed
// starting with '_atom_site' (they are simply not used). The elements are read in string by string seperated by
// white spaces. White spaces in front are always removed. If the first character is a ' or " symbol, then we need
// to find the corresponding closing symbol and treat the found string as one element (with white spaces inside).
// e.g.
//
// loop_
// _atom_type_symbol
// _atom_type_scat_source
// C 'International Tables Vol C Tables 4.2.6.8 and 6.1.1.4'
//
// The listed items are parsed until a 'loop_' statement of the next block is encountered or a command (starting
// with '_' is found.
// Note that after a full line of elements is read, the data is stored in the appropriate RASPA structures.
// Note that comment lines starting with '#' are skipped.
// Note that empty lines are skipped.

// Restrictions: 1) a single (quoted) term needs to be a single line
//               2) a new list of data items needs to start at a new line

void ReadFrameworkDefinitionCIF(void)
{
#define MAX_CIF_LINE_SIZE 1024
#define MAX_NUMBER_OF_CIF_TERMS 55
  enum{ATOM_SITE_LABEL,ATOM_SITE_TYPE_SYMBOL,ATOM_SITE_FRACT_X,ATOM_SITE_FRACT_Y,ATOM_SITE_FRACT_Z,
    ATOM_SITE_U_ISO_OR_EQUIV,ATOM_SITE_ADP_TYPE,ATOM_SITE_OCCUPANCY,ATOM_SITE_SYMMETRY_MULTIPLICITY,
    ATOM_SITE_CALC_FLAG,ATOM_SITE_REFINEMENT_FLAG,ATOM_SITE_DISORDER_ASSEMBLY,ATOM_SITE_DISORDER_GROUP,
    ATOM_SITE_ANISO_LABEL,ATOM_SITE_ANISO_U_11,ATOM_SITE_ANISO_U_22,ATOM_SITE_ANISO_U_33,
    ATOM_SITE_ANISO_U_23,ATOM_SITE_ANISO_U_13,ATOM_SITE_ANISO_U_12,ATOM_SITE_THERMAL_DISPLACE_TYPE,
    ATOM_SITE_CHARGE,ATOM_SITE_POLARIZATION,ATOM_SITE_ANISOTROPIC_TYPE,ATOM_SITE_ANISOTROPIC_DISPLACEMENT,
    ATOM_SITE_PRINT_TO_PDB,ATOM_SITE_HYBRIDIZATION,
    CITATION_ID,CITATION_AUTHOR_NAME,CITATION_COORDINATE_LINKAGE,CITATION_TITLE,CITATION_COUNTRY,CITATION_PAGE_FIRST,CITATION_PAGE_LAST,
    CITATION_YEAR,CITATION_JOURNAL_ABBREV,CITATION_JOURNAL_VOLUME,CITATION_JOURNAL_ISSUE,CITATION_JOURNAL_ID_ASTM,
    CITATION_JOURNAL_ID_ISSN,CITATION_BOOKTITLE,CITATION_BOOK_PUBLISHER,CITATION_BOOK_ID_ISBN,CITATION_SPECIAL_DETAILS,
    ATOM_TYPE_SYMBOL,ATOM_TYPE_DESCRIPTION,ATOM_TYPE_OXIDATION_NUMBER,ATOM_TYPE_NUMBER_IN_CELL,
    ATOM_TYPE_SCAT_DISPERSION_REAL,ATOM_TYPE_SCAT_DISPERSION_IMAG,ATOM_TYPE_SCAT_SOURCE,
    ATOM_TYPE_RADIUS_BOND,ATOM_TYPE_RADIUS_CONTACT,
    SYMMETRY_EQUIV_POS_SITE_ID,SYMMETRY_EQUIV_POS_AS_XYZ};

  int i,n;
  char buffer[256];
  char line[MAX_CIF_LINE_SIZE],keyword[MAX_CIF_LINE_SIZE],arguments[MAX_CIF_LINE_SIZE];
  FILE *FilePtr;
  int NumberOfAtomTypeElementsInBlock;
  int NumberOfAtomSiteElementsInBlock;
  int CifBlockData[MAX_NUMBER_OF_CIF_TERMS];
  char *arg_pointer;
  char SpaceGroup[192][32];
  int nr_sg_elements;
  REAL A,B,C;
  REAL alpha,beta,gamma;
  REAL tempd,det;
  PSEUDO_ATOM CurrentPseudoAtom;
  FRAMEWORK_ASYMMETRIC_ATOM CurrentAsymmetricAtom;
  CITATION_INFORMATION CitationInformation;
  char *LoopCondition;
  int FoundSpaceGroupSymmetryElements;
  int FoundSpaceGroupHall;
  int FoundSpaceGroupHermannMauguinString;
  int FoundSpaceGroupIntTables;
  char FoundSpaceGroupOption[32];
  char string1[32];
  char *p;
  int DefinedInLoop,PositionDefined;
  int OriginalNumberOfPseudoAtoms;
  VECTOR pos;
  int BoolChargeDefinedInCIFFile;

  // unless later '_atom_site_charge' is used, the charges are assumed not present in the CIF-file
  BoolChargeDefinedInCIFFile=FALSE;

  // set the number of pseudo-atoms before reading in the CIF-file
  OriginalNumberOfPseudoAtoms=NumberOfPseudoAtoms;

  FoundSpaceGroupSymmetryElements=0;
  FoundSpaceGroupHall=0;
  FoundSpaceGroupHermannMauguinString=0;
  FoundSpaceGroupIntTables=0;
  Framework[CurrentSystem].SpaceGroupIdentifier[CurrentFramework]=0;
  Framework[CurrentSystem].NumberOfAsymmetricAtoms[CurrentFramework]=0;
  Framework[CurrentSystem].NumberOfAtoms[CurrentFramework]=0;


  NumberOfAtomTypeElementsInBlock=0;
  NumberOfAtomSiteElementsInBlock=0;

  strcpy(FoundSpaceGroupOption,"");

  if (STREAM)
  {
#ifdef __unix__
    if (!(FilePtr=fmemopen((void *)INPUT_CRYSTAL, strlen(INPUT_CRYSTAL), "r")))
    {
      printf("Error reading streamed CIF molecule.");
      exit(1);
    }
#else
    fprintf(stderr, "Streaming only allowed on POSIX systems (for now)\n.");
    exit(1);
#endif
  }
  else
  {
    // first try to open the framework-file in the current directory,
    // and next from the repository
    sprintf(buffer,"%s.%s",
            Framework[CurrentSystem].Name[CurrentFramework],
            "cif");
    if(!(FilePtr=fopen(buffer,"r")))
    {
      sprintf(buffer,"%s/share/raspa/structures/cif/%s.%s",
              RASPA_DIRECTORY,
              Framework[CurrentSystem].Name[CurrentFramework],
              "cif");

      if(!(FilePtr=fopen(buffer,"r")))
      {
        fprintf(stderr, "Error:  file %s does not exist.\n",buffer);
        exit(1);
      }
    }
  }

  // read one line of data, and parse in 'keyword' and 'arguments'
  // LoopCondition is TRUE when data can be read and FALSE when end-of-file is reached
  LoopCondition=fgets(line,MAX_CIF_LINE_SIZE,FilePtr);
  strcpy(keyword,"");
  strcpy(arguments,"");
  sscanf(line,"%s%[^\n]",keyword,arguments);

  strcpy(CitationInformation.CitationId,"");
  strcpy(CitationInformation.CitationAuthorName,"");
  strcpy(CitationInformation.CitationCoordinateLinkage,"");
  strcpy(CitationInformation.CitationTitle,"");
  strcpy(CitationInformation.CitationCountry,"");
  strcpy(CitationInformation.CitationPageFirst,"");
  strcpy(CitationInformation.CitationPageFirst,"");
  strcpy(CitationInformation.CitationYear,"");
  strcpy(CitationInformation.CitationJournalAbbrev,"");
  strcpy(CitationInformation.CitationJournalVolume,"");
  strcpy(CitationInformation.CitationJournalIssue,"");
  strcpy(CitationInformation.CitationJournalID_ASTM,"");
  strcpy(CitationInformation.CitationJournalID_ISSN,"");
  strcpy(CitationInformation.CitationBookTitle,"");
  strcpy(CitationInformation.CitationBookPublisher,"");
  strcpy(CitationInformation.CitationBookID_ISBN,"");
  strcpy(CitationInformation.CitationSpecialDetails,"");


  do
  {
    // read loop citation info
    // ====================================
    if(strncasecmp(keyword,"_citation_",strlen("_citation_"))==0)
    {
      NumberOfAtomSiteElementsInBlock=0;
      DefinedInLoop=TRUE; // assume as default: citations are in loop

      for(i=0;i<MAX_NUMBER_OF_CIF_TERMS;i++)
        CifBlockData[i]=-1;
      do
      {
        TrimStringInPlace(arguments);

        // parse the '_citation_' data names
        if(strncasecmp(keyword,"_citation_id",strlen("_citation_id"))==0)
        {
          CifBlockData[CITATION_ID]=NumberOfAtomSiteElementsInBlock;
          if((strlen(arguments)>0)&&(arguments[0]!='_'))
          {
            DefinedInLoop=FALSE;
            RemoveQuotesAroundString(arguments);
            strcpy(CitationInformation.CitationId,arguments);
          }
          else
            NumberOfAtomSiteElementsInBlock++;
        }
        else if(strncasecmp(keyword,"_citation_author_name",strlen("_citation_author_name"))==0)
        {
          CifBlockData[CITATION_AUTHOR_NAME]=NumberOfAtomSiteElementsInBlock;
          if((strlen(arguments)>0)&&(arguments[0]!='_'))
          {
            DefinedInLoop=FALSE;
            RemoveQuotesAroundString(arguments);
            strcpy(CitationInformation.CitationAuthorName,arguments);
          }
          else
            NumberOfAtomSiteElementsInBlock++;
        }
        else if(strncasecmp(keyword,"_citation_coordinate_linkage",strlen("_citation_coordinate_linkage"))==0)
        {
          CifBlockData[CITATION_COORDINATE_LINKAGE]=NumberOfAtomSiteElementsInBlock;
          if((strlen(arguments)>0)&&(arguments[0]!='_'))
          {
            DefinedInLoop=FALSE;
            RemoveQuotesAroundString(arguments);
            strcpy(CitationInformation.CitationCoordinateLinkage,arguments);
          }
          else
            NumberOfAtomSiteElementsInBlock++;
        }
        else if(strncasecmp(keyword,"_citation_title",strlen("_citation_title"))==0)
        {
          CifBlockData[CITATION_TITLE]=NumberOfAtomSiteElementsInBlock;
          if((strlen(arguments)>0)&&(arguments[0]!='_'))
          {
            DefinedInLoop=FALSE;
            RemoveQuotesAroundString(arguments);
            strcpy(CitationInformation.CitationTitle,arguments);
          }
          else
            NumberOfAtomSiteElementsInBlock++;
        }
        else if(strncasecmp(keyword,"_citation_country",strlen("_citation_country"))==0)
        {
          CifBlockData[CITATION_COUNTRY]=NumberOfAtomSiteElementsInBlock;
          if((strlen(arguments)>0)&&(arguments[0]!='_'))
          {
            DefinedInLoop=FALSE;
            RemoveQuotesAroundString(arguments);
            strcpy(CitationInformation.CitationCountry,arguments);
          }
          else
            NumberOfAtomSiteElementsInBlock++;
        }
        else if(strncasecmp(keyword,"_citation_page_first",strlen("_citation_page_first"))==0)
        {
          CifBlockData[CITATION_PAGE_FIRST]=NumberOfAtomSiteElementsInBlock;
          if((strlen(arguments)>0)&&(arguments[0]!='_'))
          {
            DefinedInLoop=FALSE;
            RemoveQuotesAroundString(arguments);
            strcpy(CitationInformation.CitationPageFirst,arguments);
          }
          else
            NumberOfAtomSiteElementsInBlock++;
        }
        else if(strncasecmp(keyword,"_citation_page_last",strlen("_citation_page_last"))==0)
        {
          CifBlockData[CITATION_PAGE_LAST]=NumberOfAtomSiteElementsInBlock;
          if((strlen(arguments)>0)&&(arguments[0]!='_'))
          {
            DefinedInLoop=FALSE;
            RemoveQuotesAroundString(arguments);
            strcpy(CitationInformation.CitationPageLast,arguments);
          }
          else
            NumberOfAtomSiteElementsInBlock++;
        }
        else if(strncasecmp(keyword,"_citation_year",strlen("_citation_year"))==0)
        {
          CifBlockData[CITATION_YEAR]=NumberOfAtomSiteElementsInBlock;
          if((strlen(arguments)>0)&&(arguments[0]!='_'))
          {
            DefinedInLoop=FALSE;
            RemoveQuotesAroundString(arguments);
            strcpy(CitationInformation.CitationYear,arguments);
          }
          else
            NumberOfAtomSiteElementsInBlock++;
        }
        else if(strncasecmp(keyword,"_citation_journal_abbrev",strlen("_citation_journal_abbrev"))==0)
        {
          CifBlockData[CITATION_JOURNAL_ABBREV]=NumberOfAtomSiteElementsInBlock;
          if((strlen(arguments)>0)&&(arguments[0]!='_'))
          {
            DefinedInLoop=FALSE;
            RemoveQuotesAroundString(arguments);
            strcpy(CitationInformation.CitationJournalAbbrev,arguments);
          }
          else
            NumberOfAtomSiteElementsInBlock++;
        }
        else if(strncasecmp(keyword,"_citation_journal_volume",strlen("_citation_journal_volume"))==0)
        {
          CifBlockData[CITATION_JOURNAL_VOLUME]=NumberOfAtomSiteElementsInBlock;
          if((strlen(arguments)>0)&&(arguments[0]!='_'))
          {
            DefinedInLoop=FALSE;
            RemoveQuotesAroundString(arguments);
            strcpy(CitationInformation.CitationJournalVolume,arguments);
          }
          else
            NumberOfAtomSiteElementsInBlock++;
        }
        else if(strncasecmp(keyword,"_citation_journal_issue",strlen("_citation_journal_issue"))==0)
        {
          CifBlockData[CITATION_JOURNAL_ISSUE]=NumberOfAtomSiteElementsInBlock;
          if((strlen(arguments)>0)&&(arguments[0]!='_'))
          {
            DefinedInLoop=FALSE;
            RemoveQuotesAroundString(arguments);
            strcpy(CitationInformation.CitationJournalIssue,arguments);
          }
          else
            NumberOfAtomSiteElementsInBlock++;
        }
        else if(strncasecmp(keyword,"_citation_journal_id_ASTM",strlen("_citation_journal_id_ASTM"))==0)
        {
          CifBlockData[CITATION_JOURNAL_ID_ASTM]=NumberOfAtomSiteElementsInBlock;
          if((strlen(arguments)>0)&&(arguments[0]!='_'))
          {
            DefinedInLoop=FALSE;
            RemoveQuotesAroundString(arguments);
            strcpy(CitationInformation.CitationJournalID_ASTM,arguments);
          }
          else
            NumberOfAtomSiteElementsInBlock++;
        }
        else if(strncasecmp(keyword,"_citation_journal_id_ISSN",strlen("_citation_journal_id_ISSN"))==0)
        {
          CifBlockData[CITATION_JOURNAL_ID_ISSN]=NumberOfAtomSiteElementsInBlock;
          if((strlen(arguments)>0)&&(arguments[0]!='_'))
          {
            DefinedInLoop=FALSE;
            RemoveQuotesAroundString(arguments);
            strcpy(CitationInformation.CitationJournalID_ISSN,arguments);
          }
          else
            NumberOfAtomSiteElementsInBlock++;
        }
        else if(strncasecmp(keyword,"_citation_book_title",strlen("_citation_book_title"))==0)
        {
          CifBlockData[CITATION_BOOKTITLE]=NumberOfAtomSiteElementsInBlock;
          if((strlen(arguments)>0)&&(arguments[0]!='_'))
          {
            DefinedInLoop=FALSE;
            RemoveQuotesAroundString(arguments);
            strcpy(CitationInformation.CitationBookTitle,arguments);
          }
          else
            NumberOfAtomSiteElementsInBlock++;
        }
        else if(strncasecmp(keyword,"_citation_book_publisher",strlen("_citation_book_publisher"))==0)
        {
          CifBlockData[CITATION_BOOK_PUBLISHER]=NumberOfAtomSiteElementsInBlock;
          if((strlen(arguments)>0)&&(arguments[0]!='_'))
          {
            DefinedInLoop=FALSE;
            strcpy(CitationInformation.CitationBookPublisher,arguments);
          }
          else
            NumberOfAtomSiteElementsInBlock++;
        }
        else if(strncasecmp(keyword,"_citation_book_id_ISBN",strlen("_citation_book_id_ISBN"))==0)
        {
          CifBlockData[CITATION_BOOK_ID_ISBN]=NumberOfAtomSiteElementsInBlock;
          if((strlen(arguments)>0)&&(arguments[0]!='_'))
          {
            DefinedInLoop=FALSE;
            strcpy(CitationInformation.CitationBookID_ISBN,arguments);
          }
          else
            NumberOfAtomSiteElementsInBlock++;
        }
        else if(strncasecmp(keyword,"_citation_special_details",strlen("_citation_special_details"))==0)
        {
          CifBlockData[CITATION_SPECIAL_DETAILS]=NumberOfAtomSiteElementsInBlock;
          if((strlen(arguments)>0)&&(arguments[0]!='_'))
          {
            DefinedInLoop=FALSE;
            RemoveQuotesAroundString(arguments);
            strcpy(CitationInformation.CitationSpecialDetails,arguments);
          }
          else
            NumberOfAtomSiteElementsInBlock++;
        }

        // check if another data name follows
        // if so, copy the 'arguments' to 'line'
        // otherwise read 'line' from the file
        if((strlen(arguments)>0)&&(arguments[0]=='_'))
          strcpy(line,arguments);
        else
          LoopCondition=fgets(line,MAX_CIF_LINE_SIZE,FilePtr);

        // split the new line into 'keyword' and 'arguments'
        TrimStringInPlace(line);
        strcpy(keyword,"");
        strcpy(arguments,"");
        sscanf(line,"%s%[^\n]",keyword,arguments);

        // skip comments
        if(keyword[0]=='#')
        {
          strcpy(keyword,"");
          strcpy(arguments,"");
        }
      } // continue when new keyword starts with '_atom_type_', also continue when empty line
      while(((strncasecmp(keyword,"_citation_",strlen("_citation_"))==0)||(strlen(keyword)==0))&&LoopCondition);

      if(!DefinedInLoop)
      {
        Framework[CurrentSystem].CitationInformation[CurrentFramework]=(CITATION_INFORMATION*)calloc(
                                  (Framework[CurrentSystem].NumberOfCitations[CurrentFramework]+1),sizeof(CITATION_INFORMATION));
        Framework[CurrentSystem].CitationInformation[CurrentFramework][Framework[CurrentSystem].NumberOfCitations[CurrentFramework]]=CitationInformation;
        Framework[CurrentSystem].NumberOfCitations[CurrentFramework]++;
      }
      else do
      {
        // loop over data items
        if(strlen(keyword)>0)
        {
          arg_pointer=line;

          // remove spaces
          TrimStringInPlace(arg_pointer);

          strcpy(CitationInformation.CitationId,"");
          strcpy(CitationInformation.CitationAuthorName,"");
          strcpy(CitationInformation.CitationCoordinateLinkage,"");
          strcpy(CitationInformation.CitationTitle,"");
          strcpy(CitationInformation.CitationCountry,"");
          strcpy(CitationInformation.CitationPageFirst,"");
          strcpy(CitationInformation.CitationPageFirst,"");
          strcpy(CitationInformation.CitationYear,"");
          strcpy(CitationInformation.CitationJournalAbbrev,"");
          strcpy(CitationInformation.CitationJournalVolume,"");
          strcpy(CitationInformation.CitationJournalIssue,"");
          strcpy(CitationInformation.CitationJournalID_ASTM,"");
          strcpy(CitationInformation.CitationJournalID_ISSN,"");
          strcpy(CitationInformation.CitationBookTitle,"");
          strcpy(CitationInformation.CitationBookPublisher,"");
          strcpy(CitationInformation.CitationBookID_ISBN,"");
          strcpy(CitationInformation.CitationSpecialDetails,"");

          for(i=0;i<NumberOfAtomSiteElementsInBlock;i++)
          {
            strcpy(keyword,"");
            strcpy(arguments,"");

            if(sscanf(arg_pointer,"%s%n",keyword,&n)==1)
            {
              arg_pointer+=n;
              if(keyword[0]=='\'')
              {
                memmove(&keyword[0],&keyword[1],strlen(keyword));
                if(keyword[strlen(keyword)-1]!='\'')
                {
                  sscanf(arg_pointer,"%[^\']%n",arguments,&n);
                  arg_pointer+=n+1;
                  strcat(keyword, arguments);
                }
                else keyword[strlen(keyword)-1]='\0';
              }
            }
            else
            {
              // too few data items, next item might be on next line
              LoopCondition=fgets(line,MAX_CIF_LINE_SIZE,FilePtr);
              if(!LoopCondition) // if nothing could be read we have too few data items
              {
                fprintf(stderr, "Too few data items found in loop_ _citation_\n");
                exit(0);
              }
              else
              {
                strcpy(keyword,"");
                strcpy(arguments,"");
                sscanf(line,"%s%[^\n]",keyword,arguments);

                if(keyword[0]=='_') // it should not be a data name (but a data item)
                {
                  fprintf(stderr, "Too few data items found in loop_ _citation_\n");
                  exit(0);
                }

                arg_pointer=line;
                if(sscanf(arg_pointer,"%s%n",keyword,&n)==1)
                {
                  arg_pointer+=n;
                  if(keyword[0]=='\'')
                  {
                    memmove(&keyword[0],&keyword[1],strlen(keyword));
                    if(keyword[strlen(keyword)-1]!='\'')
                    {
                      sscanf(arg_pointer,"%[^\']%n",arguments,&n);
                      arg_pointer+=n+1;
                      strcat(keyword, arguments);
                    }
                    else keyword[strlen(keyword)-1]='\0';
                  }
                }
              }
            }

            if(i==CifBlockData[CITATION_ID])
            {
              if(!((strcmp(keyword,"?")==0)||(strcmp(keyword,"'?'")==0)||(strcmp(keyword,"\"?\"")==0)))
                strcpy(CitationInformation.CitationId,keyword);
            }
            else if(i==CifBlockData[CITATION_AUTHOR_NAME])
            {
              if(!((strcmp(keyword,"?")==0)||(strcmp(keyword,"'?'")==0)||(strcmp(keyword,"\"?\"")==0)))
                strcpy(CitationInformation.CitationAuthorName,keyword);
            }
            else if(i==CifBlockData[CITATION_COORDINATE_LINKAGE])
            {
              if(!((strcmp(keyword,"?")==0)||(strcmp(keyword,"'?'")==0)||(strcmp(keyword,"\"?\"")==0)))
                strcpy(CitationInformation.CitationCoordinateLinkage,keyword);
            }
            else if(i==CifBlockData[CITATION_TITLE])
            {
              if(!((strcmp(keyword,"?")==0)||(strcmp(keyword,"'?'")==0)||(strcmp(keyword,"\"?\"")==0)))
                strcpy(CitationInformation.CitationTitle,keyword);
            }
            else if(i==CifBlockData[CITATION_COUNTRY])
            {
              if(!((strcmp(keyword,"?")==0)||(strcmp(keyword,"'?'")==0)||(strcmp(keyword,"\"?\"")==0)))
                strcpy(CitationInformation.CitationCountry,keyword);
            }
            else if(i==CifBlockData[CITATION_PAGE_FIRST])
            {
              if(!((strcmp(keyword,"?")==0)||(strcmp(keyword,"'?'")==0)||(strcmp(keyword,"\"?\"")==0)))
                strcpy(CitationInformation.CitationPageFirst,keyword);
            }
            else if(i==CifBlockData[CITATION_PAGE_LAST])
            {
              if(!((strcmp(keyword,"?")==0)||(strcmp(keyword,"'?'")==0)||(strcmp(keyword,"\"?\"")==0)))
                strcpy(CitationInformation.CitationPageLast,keyword);
            }
            else if(i==CifBlockData[CITATION_YEAR])
            {
              if(!((strcmp(keyword,"?")==0)||(strcmp(keyword,"'?'")==0)||(strcmp(keyword,"\"?\"")==0)))
                strcpy(CitationInformation.CitationYear,keyword);
            }
            else if(i==CifBlockData[CITATION_JOURNAL_ABBREV])
            {
              if(!((strcmp(keyword,"?")==0)||(strcmp(keyword,"'?'")==0)||(strcmp(keyword,"\"?\"")==0)))
                strcpy(CitationInformation.CitationJournalAbbrev,keyword);
            }
            else if(i==CifBlockData[CITATION_JOURNAL_VOLUME])
            {
              if(!((strcmp(keyword,"?")==0)||(strcmp(keyword,"'?'")==0)||(strcmp(keyword,"\"?\"")==0)))
                strcpy(CitationInformation.CitationJournalVolume,keyword);
            }
            else if(i==CifBlockData[CITATION_JOURNAL_ISSUE])
            {
              if(!((strcmp(keyword,"?")==0)||(strcmp(keyword,"'?'")==0)||(strcmp(keyword,"\"?\"")==0)))
                strcpy(CitationInformation.CitationJournalIssue,keyword);
            }
            else if(i==CifBlockData[CITATION_JOURNAL_ID_ASTM])
            {
              if(!((strcmp(keyword,"?")==0)||(strcmp(keyword,"'?'")==0)||(strcmp(keyword,"\"?\"")==0)))
                strcpy(CitationInformation.CitationJournalID_ASTM,keyword);
            }
            else if(i==CifBlockData[CITATION_JOURNAL_ID_ISSN])
            {
              if(!((strcmp(keyword,"?")==0)||(strcmp(keyword,"'?'")==0)||(strcmp(keyword,"\"?\"")==0)))
                strcpy(CitationInformation.CitationJournalID_ISSN,keyword);
            }
            else if(i==CifBlockData[CITATION_BOOKTITLE])
            {
              if(!((strcmp(keyword,"?")==0)||(strcmp(keyword,"'?'")==0)||(strcmp(keyword,"\"?\"")==0)))
                strcpy(CitationInformation.CitationBookTitle,keyword);
            }
            else if(i==CifBlockData[CITATION_BOOK_PUBLISHER])
            {
              if(!((strcmp(keyword,"?")==0)||(strcmp(keyword,"'?'")==0)||(strcmp(keyword,"\"?\"")==0)))
                strcpy(CitationInformation.CitationBookPublisher,keyword);
            }
            else if(i==CifBlockData[CITATION_BOOK_ID_ISBN])
            {
              if(!((strcmp(keyword,"?")==0)||(strcmp(keyword,"'?'")==0)||(strcmp(keyword,"\"?\"")==0)))
                strcpy(CitationInformation.CitationBookID_ISBN,keyword);
            }
            else if(i==CifBlockData[CITATION_SPECIAL_DETAILS])
            {
              if(!((strcmp(keyword,"?")==0)||(strcmp(keyword,"'?'")==0)||(strcmp(keyword,"\"?\"")==0)))
                strcpy(CitationInformation.CitationSpecialDetails,keyword);
            }
            else
              fprintf(stderr, "unknown: %s\n",keyword);
          }

          if(sscanf(arg_pointer,"%s%n",keyword,&n)==1)
            fprintf(stderr, "Warning Too many data items found in loop_ _citation_, excess data ignored\n");


          Framework[CurrentSystem].CitationInformation[CurrentFramework]=(CITATION_INFORMATION*)realloc(
                                 Framework[CurrentSystem].CitationInformation[CurrentFramework],
                                (Framework[CurrentSystem].NumberOfCitations[CurrentFramework]+1)*sizeof(CITATION_INFORMATION));
          Framework[CurrentSystem].CitationInformation[CurrentFramework][Framework[CurrentSystem].NumberOfCitations[CurrentFramework]]=CitationInformation;
          Framework[CurrentSystem].NumberOfCitations[CurrentFramework]++;
        }

        LoopCondition=fgets(line,MAX_CIF_LINE_SIZE,FilePtr);
        strcpy(keyword,"");
        strcpy(arguments,"");
        sscanf(line,"%s%[^\n]",keyword,arguments);

        // skip comments
        if(keyword[0]=='#')
        {
          strcpy(keyword,"");
          strcpy(arguments,"");
        }
      } // continue as long as we do not encounter 'loop_' or a new data name and as long as we did not reach EOF
      while(((strncasecmp(keyword,"loop_",strlen("loop_"))!=0)&&(keyword[0]!='_'))&&LoopCondition);
    }

    // read loop atom-types
    // ====================================

    if(strncasecmp(keyword,"_atom_type_",strlen("_atom_type_"))==0)
    {
      NumberOfAtomSiteElementsInBlock=0;

      for(i=0;i<MAX_NUMBER_OF_CIF_TERMS;i++)
        CifBlockData[i]=-1;
      do
      {
        // parse the '_atom_type_' data names
        if(strncasecmp(keyword,"_atom_type_symbol",strlen("_atom_type_symbol"))==0) CifBlockData[ATOM_TYPE_SYMBOL]=NumberOfAtomSiteElementsInBlock;
        else if(strncasecmp(keyword,"_atom_type_description",strlen("_atom_type_description"))==0) CifBlockData[ATOM_TYPE_DESCRIPTION]=NumberOfAtomSiteElementsInBlock;
        else if(strncasecmp(keyword,"_atom_type_oxidation_number",strlen("_atom_type_oxidation_number"))==0) CifBlockData[ATOM_TYPE_OXIDATION_NUMBER]=NumberOfAtomSiteElementsInBlock;
        else if(strncasecmp(keyword,"_atom_type_number_in_cell",strlen("_atom_type_number_in_cell"))==0) CifBlockData[ATOM_TYPE_NUMBER_IN_CELL]=NumberOfAtomSiteElementsInBlock;
        else if(strncasecmp(keyword,"_atom_type_radius_bond",strlen("_atom_type_radius_bond"))==0) CifBlockData[ATOM_TYPE_RADIUS_BOND]=NumberOfAtomSiteElementsInBlock;
        else if(strncasecmp(keyword,"_atom_type_radius_contact",strlen("_atom_type_radius_contact"))==0) CifBlockData[ATOM_TYPE_RADIUS_CONTACT]=NumberOfAtomSiteElementsInBlock;
        else if(strncasecmp(keyword,"_atom_type_scat_dispersion_real",strlen("_atom_type_scat_dispersion_real"))==0) CifBlockData[ATOM_TYPE_SCAT_DISPERSION_REAL]=NumberOfAtomSiteElementsInBlock;
        else if(strncasecmp(keyword,"_atom_type_scat_dispersion_imag",strlen("_atom_type_scat_dispersion_imag"))==0) CifBlockData[ATOM_TYPE_SCAT_DISPERSION_IMAG]=NumberOfAtomSiteElementsInBlock;
        else if(strncasecmp(keyword,"_atom_type_scat_source",strlen("_atom_type_scat_source"))==0) CifBlockData[ATOM_TYPE_SCAT_SOURCE]=NumberOfAtomSiteElementsInBlock;

        // count the number of elements
        if(strlen(keyword)>0)
          NumberOfAtomSiteElementsInBlock++;

        // check if another data name follows
        // if so, copy the 'arguments' to 'line'
        // otherwise read 'line' from the file
        TrimStringInPlace(arguments);
        if((strlen(arguments)>0)&&(arguments[0]=='_'))
          strcpy(line,arguments);
        else
          LoopCondition=fgets(line,MAX_CIF_LINE_SIZE,FilePtr);

        // split the new line into 'keyword' and 'arguments'
        TrimStringInPlace(line);
        strcpy(keyword,"");
        strcpy(arguments,"");
        sscanf(line,"%s%[^\n]",keyword,arguments);

        // skip comments
        if(keyword[0]=='#')
        {
          strcpy(keyword,"");
          strcpy(arguments,"");
        }
      } // continue when new keyword starts with '_atom_type_', also continue when empty line
      while(((strncasecmp(keyword,"_atom_type_",strlen("_atom_type_"))==0)||(strlen(keyword)==0))&&LoopCondition);

      do
      {
        // loop over data items
        if(strlen(keyword)>0)
        {
          arg_pointer=line;

          // remove spaces
          TrimStringInPlace(arg_pointer);

          strcpy(CurrentPseudoAtom.Name,"");
          strcpy(CurrentPseudoAtom.PrintToPDBName,"");
          strcpy(CurrentPseudoAtom.ChemicalElement,"");
          strcpy(CurrentPseudoAtom.ScatteringSource,"");
          CurrentPseudoAtom.Occupancy=1.0;
          CurrentPseudoAtom.FrameworkAtom=TRUE;
          CurrentPseudoAtom.PrintToPDB=TRUE;
          CurrentPseudoAtom.ScatteringType=0;
          CurrentPseudoAtom.AnomalousScatteringType=0;
          CurrentPseudoAtom.TemperatureFactor=0.0;
          CurrentPseudoAtom.AnisotropicTemperatureFactor.ax=0.0;
          CurrentPseudoAtom.AnisotropicTemperatureFactor.ay=0.0;
          CurrentPseudoAtom.AnisotropicTemperatureFactor.az=0.0;
          CurrentPseudoAtom.AnisotropicTemperatureFactor.bx=0.0;
          CurrentPseudoAtom.AnisotropicTemperatureFactor.by=0.0;
          CurrentPseudoAtom.AnisotropicTemperatureFactor.bz=0.0;
          CurrentPseudoAtom.AnisotropicTemperatureFactor.cx=0.0;
          CurrentPseudoAtom.AnisotropicTemperatureFactor.cy=0.0;
          CurrentPseudoAtom.AnisotropicTemperatureFactor.cz=0.0;
          CurrentPseudoAtom.ScatteringDispersionImaginary=0.0;
          CurrentPseudoAtom.ScatteringDispersionImaginary=0.0;
          CurrentPseudoAtom.Mass=0.0;
          CurrentPseudoAtom.Charge1=0.0;
          CurrentPseudoAtom.Polarization.ax=CurrentPseudoAtom.Polarization.ay=CurrentPseudoAtom.Polarization.az=0.0;
          CurrentPseudoAtom.Polarization.bx=CurrentPseudoAtom.Polarization.by=CurrentPseudoAtom.Polarization.bz=0.0;
          CurrentPseudoAtom.Polarization.cx=CurrentPseudoAtom.Polarization.cy=CurrentPseudoAtom.Polarization.cz=0.0;
          CurrentPseudoAtom.HasCharges=TRUE;
          CurrentPseudoAtom.IsPolarizable=FALSE;
          CurrentPseudoAtom.Interaction=TRUE;
          CurrentPseudoAtom.Radius=0.0;
          CurrentPseudoAtom.Connectivity=0;
          CurrentPseudoAtom.TinkerType=0;
          CurrentPseudoAtom.AnisotropicCorrection=FALSE;
          CurrentPseudoAtom.AnisotropicDisplacement=0.0;
          CurrentPseudoAtom.AnisotropicType=0.0;

          for(i=0;i<NumberOfAtomSiteElementsInBlock;i++)
          {
            strcpy(keyword,"");
            strcpy(arguments,"");

            if(sscanf(arg_pointer,"%s%n",keyword,&n)==1)
            {
              arg_pointer+=n;
              if(keyword[0]=='\'')
              {
                memmove(&keyword[0],&keyword[1],strlen(keyword));
                if(keyword[strlen(keyword)-1]!='\'')
                {
                  sscanf(arg_pointer,"%[^\']%n",arguments,&n);
                  arg_pointer+=n+1;
                  strcat(keyword, arguments);
                }
                else keyword[strlen(keyword)-1]='\0';
              }
            }
            else
            {
              LoopCondition=fgets(line,MAX_CIF_LINE_SIZE,FilePtr);
              if(!LoopCondition) // if nothing could be read we have too few data items
              {
                fprintf(stderr, "Too few data items found in loop_ _atom_site_\n");
                exit(0);
              }
              else
              {
                strcpy(keyword,"");
                strcpy(arguments,"");
                sscanf(line,"%s%[^\n]",keyword,arguments);

                if(keyword[0]=='_') // it should not be a data name (but a data item)
                {
                  fprintf(stderr, "Too few data items found in loop_ _atom_site_\n");
                  exit(0);
                }

                arg_pointer=line;
                TrimStringInPlace(arg_pointer);

                if(sscanf(arg_pointer,"%s%n",keyword,&n)==1)
                {
                  arg_pointer+=n;
                  if(keyword[0]=='\'')
                  {
                    memmove(&keyword[0],&keyword[1],strlen(keyword));
                    if(keyword[strlen(keyword)-1]!='\'')
                    {
                      sscanf(arg_pointer,"%[^\']%n",arguments,&n);
                      arg_pointer+=n+1;
                      strcat(keyword, arguments);
                    }
                    else keyword[strlen(keyword)-1]='\0';
                  }
                }
              }
            }

            if(i==CifBlockData[ATOM_TYPE_SYMBOL])
            {
              // must match '_atom_site_type_symbol'
              // e.g. Cu2+
            }
            else if(i==CifBlockData[ATOM_TYPE_DESCRIPTION])
            {
              // not used
            }
            else if(i==CifBlockData[ATOM_TYPE_NUMBER_IN_CELL])
            {
              // not used (recomputed after expanding from space group to P1)
            }
            else if(i==CifBlockData[ATOM_TYPE_OXIDATION_NUMBER])
            {
              // not used
            }
            else if(i==CifBlockData[ATOM_TYPE_RADIUS_BOND])
            {
              // not used
            }
            else if(i==CifBlockData[ATOM_TYPE_RADIUS_CONTACT])
            {
              // not used
            }
            else if(i==CifBlockData[ATOM_TYPE_SCAT_DISPERSION_REAL])
            {
              // not read, valuea are used from tables
            }
            else if(i==CifBlockData[ATOM_TYPE_SCAT_DISPERSION_IMAG])
            {
              // not read, valuea are used from tables
            }
            else if(i==CifBlockData[ATOM_TYPE_SCAT_SOURCE])
            {
              // not used
            }
            else
              fprintf(stderr, "unknown: %s\n",keyword);
          }

          if(sscanf(arg_pointer,"%s%n",keyword,&n)==1)
            fprintf(stderr, "Warning Too many data items found in loop_ _atom_type_, excess data ignored\n");
        }

        LoopCondition=fgets(line,MAX_CIF_LINE_SIZE,FilePtr);
        strcpy(keyword,"");
        strcpy(arguments,"");
        sscanf(line,"%s%[^\n]",keyword,arguments);

        // skip comments
        if(keyword[0]=='#')
        {
          strcpy(keyword,"");
          strcpy(arguments,"");
        }
      } // continue as long as we do not encounter 'loop_' or a new data name and as long as we did not reach EOF
      while(((strncasecmp(keyword,"loop_",strlen("loop_"))!=0)&&(keyword[0]!='_'))&&LoopCondition);
    }

    // read loop atom-sites
    // ====================================

    if((strncasecmp(keyword,"_atom_site_",strlen("_atom_site_"))==0)&&(strncasecmp(keyword,"_atom_site_aniso",strlen("_atom_site_aniso"))!=0))
    {
      NumberOfAtomSiteElementsInBlock=0;

      for(i=0;i<MAX_NUMBER_OF_CIF_TERMS;i++)
        CifBlockData[i]=-1;
      do
      {
        // parse the '_atom_site_' data names
        if(strncasecmp(keyword,"_atom_site_label",strlen("_atom_site_label"))==0) CifBlockData[ATOM_SITE_LABEL]=NumberOfAtomSiteElementsInBlock;
        else if(strncasecmp(keyword,"_atom_site_type_symbol",strlen("_atom_site_type_symbol"))==0) CifBlockData[ATOM_SITE_TYPE_SYMBOL]=NumberOfAtomSiteElementsInBlock;
        else if(strncasecmp(keyword,"_atom_site_fract_x",strlen("_atom_site_fract_x"))==0) CifBlockData[ATOM_SITE_FRACT_X]=NumberOfAtomSiteElementsInBlock;
        else if(strncasecmp(keyword,"_atom_site_fract_y",strlen("_atom_site_fract_y"))==0) CifBlockData[ATOM_SITE_FRACT_Y]=NumberOfAtomSiteElementsInBlock;
        else if(strncasecmp(keyword,"_atom_site_fract_z",strlen("_atom_site_fract_z"))==0) CifBlockData[ATOM_SITE_FRACT_Z]=NumberOfAtomSiteElementsInBlock;
        else if(strncasecmp(keyword,"_atom_site_U_iso_or_equiv",strlen("_atom_site_U_iso_or_equiv"))==0) CifBlockData[ATOM_SITE_U_ISO_OR_EQUIV]=NumberOfAtomSiteElementsInBlock;
        else if(strncasecmp(keyword,"_atom_site_thermal_displace_type",strlen("_atom_site_thermal_displace_type"))==0) CifBlockData[ATOM_SITE_THERMAL_DISPLACE_TYPE]=NumberOfAtomSiteElementsInBlock;
        else if(strncasecmp(keyword,"_atom_site_adp_type",strlen("_atom_site_adp_type"))==0) CifBlockData[ATOM_SITE_ADP_TYPE]=NumberOfAtomSiteElementsInBlock;
        else if(strncasecmp(keyword,"_atom_site_occupancy",strlen("_atom_site_occupancy"))==0) CifBlockData[ATOM_SITE_OCCUPANCY]=NumberOfAtomSiteElementsInBlock;
        else if(strncasecmp(keyword,"_atom_site_symmetry_multiplicity",strlen("_atom_site_symmetry_multiplicity"))==0) CifBlockData[ATOM_SITE_SYMMETRY_MULTIPLICITY]=NumberOfAtomSiteElementsInBlock;
        else if(strncasecmp(keyword,"_atom_site_calc_flag",strlen("_atom_site_calc_flag"))==0) CifBlockData[ATOM_SITE_CALC_FLAG]=NumberOfAtomSiteElementsInBlock;
        else if(strncasecmp(keyword,"_atom_site_refinement_flags",strlen("_atom_site_refinement_flags"))==0) CifBlockData[ATOM_SITE_REFINEMENT_FLAG]=NumberOfAtomSiteElementsInBlock;
        else if(strncasecmp(keyword,"_atom_site_disorder_assembly",strlen("_atom_site_disorder_assembly"))==0) CifBlockData[ATOM_SITE_DISORDER_ASSEMBLY]=NumberOfAtomSiteElementsInBlock;
        else if(strncasecmp(keyword,"_atom_site_disorder_group",strlen("_atom_site_disorder_group"))==0) CifBlockData[ATOM_SITE_DISORDER_GROUP]=NumberOfAtomSiteElementsInBlock;
        else if(strncasecmp(keyword,"_atom_site_aniso_label",strlen("_atom_site_aniso_label"))==0) CifBlockData[ATOM_SITE_ANISO_LABEL]=NumberOfAtomSiteElementsInBlock;
        else if(strncasecmp(keyword,"_atom_site_aniso_U_11",strlen("_atom_site_aniso_U_11"))==0) CifBlockData[ATOM_SITE_ANISO_U_11]=NumberOfAtomSiteElementsInBlock;
        else if(strncasecmp(keyword,"_atom_site_aniso_U_22",strlen("_atom_site_aniso_U_22"))==0) CifBlockData[ATOM_SITE_ANISO_U_22]=NumberOfAtomSiteElementsInBlock;
        else if(strncasecmp(keyword,"_atom_site_aniso_U_33",strlen("_atom_site_aniso_U_33"))==0) CifBlockData[ATOM_SITE_ANISO_U_33]=NumberOfAtomSiteElementsInBlock;
        else if(strncasecmp(keyword,"_atom_site_aniso_U_23",strlen("_atom_site_aniso_U_23"))==0) CifBlockData[ATOM_SITE_ANISO_U_23]=NumberOfAtomSiteElementsInBlock;
        else if(strncasecmp(keyword,"_atom_site_aniso_U_13",strlen("_atom_site_aniso_U_13"))==0) CifBlockData[ATOM_SITE_ANISO_U_13]=NumberOfAtomSiteElementsInBlock;
        else if(strncasecmp(keyword,"_atom_site_aniso_U_12",strlen("_atom_site_aniso_U_12"))==0) CifBlockData[ATOM_SITE_ANISO_U_12]=NumberOfAtomSiteElementsInBlock;
        else if(strncasecmp(keyword,"_atom_site_charge",strlen("_atom_site_charge"))==0)
        {
          BoolChargeDefinedInCIFFile=TRUE;
          CifBlockData[ATOM_SITE_CHARGE]=NumberOfAtomSiteElementsInBlock;
        }
        else if(strncasecmp(keyword,"_atom_site_polarization",strlen("_atom_site_polarization"))==0) CifBlockData[ATOM_SITE_POLARIZATION]=NumberOfAtomSiteElementsInBlock;
        else if(strncasecmp(keyword,"_atom_site_anisotropic_type",strlen("_atom_site_anisotropic_type"))==0) CifBlockData[ATOM_SITE_ANISOTROPIC_TYPE]=NumberOfAtomSiteElementsInBlock;
        else if(strncasecmp(keyword,"_atom_site_anisotropic_displacement",strlen("_atom_site_anisotropic_displacement"))==0) CifBlockData[ATOM_SITE_ANISOTROPIC_DISPLACEMENT]=NumberOfAtomSiteElementsInBlock;
        else if(strncasecmp(keyword,"_atom_site_print_to_pdb",strlen("_atom_site_print_to_pdb"))==0) CifBlockData[ATOM_SITE_PRINT_TO_PDB]=NumberOfAtomSiteElementsInBlock;
        else if(strncasecmp(keyword,"_atom_site_hybridization",strlen("_atom_site_hybridization"))==0) CifBlockData[ATOM_SITE_HYBRIDIZATION]=NumberOfAtomSiteElementsInBlock;


        // count the number of elements
        if(strlen(keyword)>0)
          NumberOfAtomSiteElementsInBlock++;

        // check if another data name follows
        // if so, copy the 'arguments' to 'line'
        // otherwise read 'line' from the file
        TrimStringInPlace(arguments);
        if((strlen(arguments)>0)&&(arguments[0]=='_'))
          strcpy(line,arguments);
        else
          LoopCondition=fgets(line,MAX_CIF_LINE_SIZE,FilePtr);

        // split the new line into 'keyword' and 'arguments'
        strcpy(keyword,"");
        strcpy(arguments,"");
        sscanf(line,"%s%[^\n]",keyword,arguments);

        // skip comments
        if(keyword[0]=='#')
        {
          strcpy(keyword,"");
          strcpy(arguments,"");
        }
      }
      while(((strncasecmp(keyword,"_atom_site_",strlen("_atom_site_"))==0)||(strlen(keyword)==0))&&LoopCondition);

      do
      {
        // loop over data items
        if(strlen(keyword)>0)
        {
          arg_pointer=line;

          // remove spaces
          TrimStringInPlace(arg_pointer);

          strcpy(CurrentPseudoAtom.Name,"");
          strcpy(CurrentPseudoAtom.PrintToPDBName,"");
          strcpy(CurrentPseudoAtom.ChemicalElement,"");
          strcpy(CurrentPseudoAtom.OxidationStateString,"");
          CurrentPseudoAtom.OxidationState=0.0;
          strcpy(CurrentPseudoAtom.ScatteringSource,"na");
          CurrentPseudoAtom.Occupancy=1.0;
          CurrentPseudoAtom.FrameworkAtom=TRUE;
          CurrentPseudoAtom.PrintToPDB=TRUE;
          CurrentPseudoAtom.ScatteringType=0;
          CurrentPseudoAtom.AnomalousScatteringType=0;
          CurrentPseudoAtom.TemperatureFactor=0.0;
          CurrentPseudoAtom.AnisotropicTemperatureFactor.ax=0.0;
          CurrentPseudoAtom.AnisotropicTemperatureFactor.ay=0.0;
          CurrentPseudoAtom.AnisotropicTemperatureFactor.az=0.0;
          CurrentPseudoAtom.AnisotropicTemperatureFactor.bx=0.0;
          CurrentPseudoAtom.AnisotropicTemperatureFactor.by=0.0;
          CurrentPseudoAtom.AnisotropicTemperatureFactor.bz=0.0;
          CurrentPseudoAtom.AnisotropicTemperatureFactor.cx=0.0;
          CurrentPseudoAtom.AnisotropicTemperatureFactor.cy=0.0;
          CurrentPseudoAtom.AnisotropicTemperatureFactor.cz=0.0;
          CurrentPseudoAtom.ScatteringDispersionImaginary=0.0;
          CurrentPseudoAtom.ScatteringDispersionImaginary=0.0;
          CurrentPseudoAtom.Mass=0.0;
          CurrentPseudoAtom.Charge1=0.0;
          CurrentPseudoAtom.ChargeDefinitionType=CHARGE_NOT_DEFINED;   // charge not yet known
          CurrentPseudoAtom.Polarization.ax=CurrentPseudoAtom.Polarization.ay=CurrentPseudoAtom.Polarization.az=0.0;
          CurrentPseudoAtom.Polarization.bx=CurrentPseudoAtom.Polarization.by=CurrentPseudoAtom.Polarization.bz=0.0;
          CurrentPseudoAtom.Polarization.cx=CurrentPseudoAtom.Polarization.cy=CurrentPseudoAtom.Polarization.cz=0.0;
          CurrentPseudoAtom.HasCharges=TRUE;
          CurrentPseudoAtom.IsPolarizable=FALSE;
          CurrentPseudoAtom.Interaction=TRUE;
          CurrentPseudoAtom.Radius=0.0;
          CurrentPseudoAtom.Connectivity=0;
          CurrentPseudoAtom.TinkerType=0;
          CurrentPseudoAtom.AnisotropicCorrection=FALSE;
          CurrentPseudoAtom.AnisotropicDisplacement=0.0;
          CurrentPseudoAtom.AnisotropicType=0.0;
          CurrentPseudoAtom.HasVDWInteraction=TRUE;
          CurrentPseudoAtom.Hybridization=HYBRIDIZATION_UNINITIALIZED;

          CurrentAsymmetricAtom.Position.x=0.0;
          CurrentAsymmetricAtom.Position.y=0.0;
          CurrentAsymmetricAtom.Position.z=0.0;
          CurrentAsymmetricAtom.Charge=0.0;

          PositionDefined=TRUE;

          for(i=0;i<NumberOfAtomSiteElementsInBlock;i++)
          {
            strcpy(keyword,"");
            strcpy(arguments,"");

            if(sscanf(arg_pointer,"%s%n",keyword,&n)==1)
            {
              // update the pointer to the next element
              arg_pointer+=n;

              // if starting with a single quote then look for next single quote
              // consider everything in between as the 'argument'
              // update the pointer to element after the terminal single quote
              if(keyword[0]=='\'')
              {
                memmove(&keyword[0],&keyword[1],strlen(keyword));
                if(keyword[strlen(keyword)-1]!='\'')
                {
                  sscanf(arg_pointer,"%[^\']%n",arguments,&n);
                  arg_pointer+=n+1;
                  strcat(keyword, arguments);
                }
                else keyword[strlen(keyword)-1]='\0';
              }

              // if starting with a double quote then look for next double quote
              // consider everything in between as the 'argument'
              // update the pointer to element after the terminal double quote
              if(keyword[0]=='\"')
              {
                memmove(&keyword[0],&keyword[1],strlen(keyword));
                if(keyword[strlen(keyword)-1]!='\"')
                {
                  sscanf(arg_pointer,"%[^\"]%n",arguments,&n);
                  arg_pointer+=n+1;
                  strcat(keyword, arguments);
                }
                else keyword[strlen(keyword)-1]='\0';
              }
            }
            else
            {
              LoopCondition=fgets(line,MAX_CIF_LINE_SIZE,FilePtr);
              if(!LoopCondition) // if nothing could be read we have too few data items
              {
                fprintf(stderr, "Too few data items found in loop_ _atom_site_\n");
                exit(0);
              }
              else
              {
                strcpy(keyword,"");
                strcpy(arguments,"");
                sscanf(line,"%s%[^\n]",keyword,arguments);

                if(keyword[0]=='_') // it should not be a data name (but a data item)
                {
                  fprintf(stderr, "Too few data items found in loop_ _atom_site_\n");
                  exit(0);
                }

                arg_pointer=line;
                TrimStringInPlace(arg_pointer);

                if(sscanf(arg_pointer,"%s%n",keyword,&n)==1)
                {
                  // update the pointer to the next element
                  arg_pointer+=n;

                  // if starting with a single quote then look for next single quote
                  // consider everything in between as the 'argument'
                  // update the pointer to element after the terminal single quote
                  if(keyword[0]=='\'')
                  {
                    memmove(&keyword[0],&keyword[1],strlen(keyword));
                    if(keyword[strlen(keyword)-1]!='\'')
                    {
                      sscanf(arg_pointer,"%[^\']%n",arguments,&n);
                      arg_pointer+=n+1;
                      strcat(keyword, arguments);
                    }
                    else keyword[strlen(keyword)-1]='\0';
                  }

                  // if starting with a double quote then look for next double quote
                  // consider everything in between as the 'argument'
                  // update the pointer to element after the terminal double quote
                  if(keyword[0]=='\"')
                  {
                    memmove(&keyword[0],&keyword[1],strlen(keyword));
                    if(keyword[strlen(keyword)-1]!='\"')
                    {
                      sscanf(arg_pointer,"%[^\"]%n",arguments,&n);
                      arg_pointer+=n+1;
                      strcat(keyword, arguments);
                    }
                    else keyword[strlen(keyword)-1]='\0';
                  }
                }
              }
            }

            if(i==CifBlockData[ATOM_SITE_LABEL])
            {
              strcpy(CurrentPseudoAtom.Name,keyword);
              if(Framework[CurrentSystem].RemoveAtomNumberCodeFromLabel[CurrentFramework])
              {
                // remove everything after the first '_'
                p=strchr(keyword,'_');
                if(p) *p='\0';

                // get CIF component-0
                sscanf(keyword,"%[a-zA-Z]%n",CurrentPseudoAtom.ChemicalElement,&n);
                p=keyword+n;
                n=0;
                strcpy(string1,"");
                string1[2]='\0';
                while(sscanf(p,"%[0-9]%[+-]%n",&string1[0],&string1[1],&n)>=2)
                {
                  strcat(CurrentPseudoAtom.OxidationStateString,string1);
                  p+=n;
                }
                strcpy(string1,"");
                sscanf(p,"%s",string1);

                // 'Name' consists of the chemical-element part concatenated with the oxydation state
                strcpy(CurrentPseudoAtom.Name,CurrentPseudoAtom.ChemicalElement);
                strcat(CurrentPseudoAtom.Name,CurrentPseudoAtom.OxidationStateString);
              }
            }
            else if(i==CifBlockData[ATOM_SITE_TYPE_SYMBOL])
            {
              // must match '_atom_type_symbol'
              sscanf(keyword,"%[a-zA-Z]%n",CurrentPseudoAtom.ChemicalElement,&n);
              p=keyword+n;
              n=0;
              strcpy(string1,"");
              string1[2]='\0';
              while(sscanf(p,"%[0-9]%[+-]%n",&string1[0],&string1[1],&n)>=2)
              {
                strcat(CurrentPseudoAtom.OxidationStateString,string1);
                p+=n;
              }

              // 'Name' consists of the chemical-element part concatenated with the oxydation state
              strcpy(CurrentPseudoAtom.PrintToPDBName,CurrentPseudoAtom.ChemicalElement);
              CurrentPseudoAtom.OxidationState=atoi(CurrentPseudoAtom.OxidationStateString);
            }
            else if(i==CifBlockData[ATOM_SITE_FRACT_X])
            {
              if(sscanf(keyword,"%lf",&pos.x)==0)
                PositionDefined=FALSE;
              if(Framework[CurrentSystem].RestrictFrameworkAtomsToBox)
              {
                pos.x-=(REAL)NINT(pos.x);
                if(pos.x<0) pos.x+=1.0;
              }
              CurrentAsymmetricAtom.Position.x=pos.x;

            }
            else if(i==CifBlockData[ATOM_SITE_FRACT_Y])
            {
              if(sscanf(keyword,"%lf",&pos.y)==0)
                PositionDefined=FALSE;
              if(Framework[CurrentSystem].RestrictFrameworkAtomsToBox)
              {
                pos.y-=(REAL)NINT(pos.y);
                if(pos.y<0) pos.y+=1.0;
              }
              CurrentAsymmetricAtom.Position.y=pos.y;
            }
            else if(i==CifBlockData[ATOM_SITE_FRACT_Z])
            {
              if(sscanf(keyword,"%lf",&pos.z)==0)
                PositionDefined=FALSE;
              if(Framework[CurrentSystem].RestrictFrameworkAtomsToBox)
              {
                pos.z-=(REAL)NINT(pos.z);
                if(pos.z<0) pos.z+=1.0;
              }
              CurrentAsymmetricAtom.Position.z=pos.z;

            }
            else if(i==CifBlockData[ATOM_SITE_CHARGE])
            {
              if(UseChargesFromCIFFile)
              {
                // set the charge to the asymmetric atom and to the pseudo-atom
                sscanf(keyword,"%lf",&CurrentAsymmetricAtom.Charge);

    CurrentPseudoAtom.Charge1=CurrentAsymmetricAtom.Charge;
                CurrentPseudoAtom.ChargeDefinitionType=CHARGE_ATOM_FROM_STRUCTURE_FILE;
              }
            }
            else if(i==CifBlockData[ATOM_SITE_POLARIZATION])
            {
              sscanf(keyword,"%lf",&tempd);
              CurrentPseudoAtom.Polarization.ax=tempd;
              CurrentPseudoAtom.Polarization.by=tempd;
              CurrentPseudoAtom.Polarization.cz=tempd;
            }
            else if(i==CifBlockData[ATOM_SITE_ANISOTROPIC_TYPE])
            {
              if(strncmp(keyword,"relative",MAX2(strlen(keyword),strlen("relative")))==0) CurrentPseudoAtom.AnisotropicType=RELATIVE;
              if(strncmp(keyword,"absolute",MAX2(strlen(keyword),strlen("absolute")))==0) CurrentPseudoAtom.AnisotropicType=ABSOLUTE;
            }
            else if(i==CifBlockData[ATOM_SITE_ANISOTROPIC_DISPLACEMENT])
            {
              sscanf(keyword,"%lf",&CurrentPseudoAtom.AnisotropicDisplacement);
            }
            else if(i==CifBlockData[ATOM_SITE_PRINT_TO_PDB])
            {
              if(strncmp(keyword,"yes",MAX2(strlen(keyword),strlen("yes")))==0) CurrentPseudoAtom.PrintToPDB=TRUE;
              if(strncmp(keyword,"no",MAX2(strlen(keyword),strlen("no")))==0) CurrentPseudoAtom.PrintToPDB=FALSE;
            }
            else if(i==CifBlockData[ATOM_SITE_U_ISO_OR_EQUIV])
            {
            }
            else if(i==CifBlockData[ATOM_SITE_THERMAL_DISPLACE_TYPE])
            {
            }
            else if(i==CifBlockData[ATOM_SITE_ADP_TYPE])
            {
            }
            else if(i==CifBlockData[ATOM_SITE_OCCUPANCY])
            {
            }
            else if(i==CifBlockData[ATOM_SITE_SYMMETRY_MULTIPLICITY])
            {
            }
            else if(i==CifBlockData[ATOM_SITE_CALC_FLAG])
            {
            }
            else if(i==CifBlockData[ATOM_SITE_REFINEMENT_FLAG])
            {
            }
            else if(i==CifBlockData[ATOM_SITE_DISORDER_ASSEMBLY])
            {
            }
            else if(i==CifBlockData[ATOM_SITE_DISORDER_GROUP])
            {
            }
            else if(i==CifBlockData[ATOM_SITE_ANISO_LABEL])
            {
            }
            else if(i==CifBlockData[ATOM_SITE_ANISO_U_11])
            {
            }
            else if(i==CifBlockData[ATOM_SITE_ANISO_U_22])
            {
            }
            else if(i==CifBlockData[ATOM_SITE_ANISO_U_33])
            {
            }
            else if(i==CifBlockData[ATOM_SITE_ANISO_U_23])
            {
            }
            else if(i==CifBlockData[ATOM_SITE_ANISO_U_13])
            {
            }
            else if(i==CifBlockData[ATOM_SITE_ANISO_U_12])
            {
            }
            else if(i==CifBlockData[ATOM_SITE_HYBRIDIZATION])
            {
              if(strncmp(keyword,"sp3",MAX2(strlen(keyword),strlen("sp3")))==0) CurrentPseudoAtom.Hybridization=SP3;
              if(strncmp(keyword,"sp2",MAX2(strlen(keyword),strlen("sp2")))==0) CurrentPseudoAtom.Hybridization=SP2;
              if(strncmp(keyword,"sp1",MAX2(strlen(keyword),strlen("sp1")))==0) CurrentPseudoAtom.Hybridization=SP;
              if(strncmp(keyword,"sp",MAX2(strlen(keyword),strlen("sp")))==0) CurrentPseudoAtom.Hybridization=SP;
            }
            else
              fprintf(stderr, "unknown: %s\n",keyword);

          }
          if(sscanf(arg_pointer,"%s%n",keyword,&n)==1)
            fprintf(stderr, "Warning Too many data items found in loop_ _atom_site_, excess data ignored\n");

          // _atom_site_type_symbol not set, get symbol from _atom_site_label component 0
          if(strlen(CurrentPseudoAtom.ChemicalElement)==0)
          {
            sscanf(CurrentPseudoAtom.Name,"%[a-zA-Z]",CurrentPseudoAtom.ChemicalElement);
            strcpy(CurrentPseudoAtom.PrintToPDBName,CurrentPseudoAtom.ChemicalElement);
          }

          // no label-name set for the atom
          if(strlen(CurrentPseudoAtom.Name)==0)
          {
            sscanf(CurrentPseudoAtom.ChemicalElement,"%[a-zA-Z]",CurrentPseudoAtom.Name);
            strcpy(CurrentPseudoAtom.PrintToPDBName,CurrentPseudoAtom.Name);
          }

         

          // chemical element set, remove oxidation state to get element for use in pdb
          sscanf(CurrentPseudoAtom.ChemicalElement,"%[a-zA-Z]",CurrentPseudoAtom.PrintToPDBName);

          CurrentPseudoAtom.Mass=GetAtomicMass(CurrentPseudoAtom.PrintToPDBName);
          CurrentPseudoAtom.Radius=GetCovalentRadius(CurrentPseudoAtom.PrintToPDBName);
          if(CurrentPseudoAtom.Polarization.ax<1e-8) // if not yet set, look up in table
            CurrentPseudoAtom.Polarization.ax=GetAtomPolarization(CurrentPseudoAtom.ChemicalElement)/COULOMBIC_CONVERSION_FACTOR;
          if(CurrentPseudoAtom.Polarization.by<1e-8) // if not yet set, look up in table
            CurrentPseudoAtom.Polarization.by=GetAtomPolarization(CurrentPseudoAtom.ChemicalElement)/COULOMBIC_CONVERSION_FACTOR;
          if(CurrentPseudoAtom.Polarization.cz<1e-8) // if not yet set, look up in table
            CurrentPseudoAtom.Polarization.cz=GetAtomPolarization(CurrentPseudoAtom.ChemicalElement)/COULOMBIC_CONVERSION_FACTOR;

          if((!ComputePolarization)||((fabs(CurrentPseudoAtom.Polarization.ax)<1e-10)&&(fabs(CurrentPseudoAtom.Polarization.by)<1e-10)&&(fabs(CurrentPseudoAtom.Polarization.cz)<1e-10)))
            CurrentPseudoAtom.IsPolarizable=FALSE;
          else
            CurrentPseudoAtom.IsPolarizable=TRUE;

          // for polarizable atoms which have zero charge, the inclusion is still needed to
          // compute the electric-field at that position
          if((fabs(CurrentPseudoAtom.Charge1)<1e-10)&&(!CurrentPseudoAtom.IsPolarizable))
            CurrentPseudoAtom.HasCharges=FALSE;
          else
            CurrentPseudoAtom.HasCharges=TRUE;

          if(fabs(CurrentPseudoAtom.AnisotropicDisplacement)<1e-10)
            CurrentPseudoAtom.AnisotropicCorrection=FALSE;
          else
            CurrentPseudoAtom.AnisotropicCorrection=TRUE;

          // only add when the pseudo-atom type is not yet present
          // this means that the 'pseudo_atoms.def' file has precedence over the cif-file settings
          CurrentAsymmetricAtom.Type=AddPseudoAtom(CurrentPseudoAtom);

          // 1) pseudo-atom already exist -> charge from the pseudo-atom to the asymmetric atom
          // 2) pseudo-atom from CIF is new -> the new pseudo-atom has been filled in with the charge from the cif-file
          // 3) if e.g. 'O2' is new, and then another 'O2' is read, then the charge is taken from CIF-file

          // MODIFIED !!
          //if((CurrentAsymmetricAtom.Type<OriginalNumberOfPseudoAtoms)||(!UseChargesFromCIFFile))
          //if((CurrentAsymmetricAtom.Type<OriginalNumberOfPseudoAtoms)&&(!UseChargesFromCIFFile)&&BoolChargeDefinedInCIFFile)
          //  CurrentAsymmetricAtom.Charge=PseudoAtoms[CurrentAsymmetricAtom.Type].Charge1;

          if(CurrentAsymmetricAtom.Type<OriginalNumberOfPseudoAtoms) // pseudo-atom already defined
          {
            // modify the existing pseudo-atom type as being a 'framework-atom'
            PseudoAtoms[CurrentAsymmetricAtom.Type].FrameworkAtom=TRUE;

            if(UseChargesFromCIFFile)  // if 'UseChargesFromCIFFile' then use the charges from the CIF-file
            {
              PseudoAtoms[CurrentAsymmetricAtom.Type].ChargeDefinitionType=CurrentPseudoAtom.ChargeDefinitionType;
              if(CurrentPseudoAtom.ChargeDefinitionType!=CHARGE_NOT_DEFINED)
                PseudoAtoms[CurrentAsymmetricAtom.Type].ChargeDefinitionType=CHARGE_ATOM_FROM_STRUCTURE_FILE;
              else
              {
                PseudoAtoms[CurrentAsymmetricAtom.Type].ChargeDefinitionType=CHARGE_ATOM_FROM_PSEUDO_ATOM_DEFINITION;
                CurrentAsymmetricAtom.Charge=PseudoAtoms[CurrentAsymmetricAtom.Type].Charge1;
              }
            }
            else
            {
              CurrentAsymmetricAtom.Charge=PseudoAtoms[CurrentAsymmetricAtom.Type].Charge1;
              PseudoAtoms[CurrentAsymmetricAtom.Type].ChargeDefinitionType=CHARGE_ATOM_FROM_PSEUDO_ATOM_DEFINITION;
            }
          }

          if(PositionDefined)
            AddAsymmetricAtom(CurrentAsymmetricAtom);
        }

        // split the new line into 'keyword' and 'arguments'
        LoopCondition=fgets(line,MAX_CIF_LINE_SIZE,FilePtr);
        strcpy(keyword,"");
        strcpy(arguments,"");
        sscanf(line,"%s%[^\n]",keyword,arguments);

        // skip comments
        if(keyword[0]=='#')
        {
          strcpy(keyword,"");
          strcpy(arguments,"");
        }
      } // continue as long as we do not encounter 'loop_' or a new data name and as long as we did not reach EOF
      while(((strncasecmp(keyword,"loop_",strlen("loop_"))!=0)&&(keyword[0]!='_'))&&LoopCondition);
    }

    // There are 3 ways of setting the space group
    // a) using '_symmetry_space_group_name_Hall'
    // b) using '_symmetry_equiv_pos_as_xyz'
    // c) using '_symmetry_space_group_name_H-M'
    // d) using '_symmetry_Int_Tables_number'
    //
    // Option (a) is prefered because it uniquely defines the space group (standard or non-standard as well as alternative origins)
    // Second best is (b) which defines the space group also uniquely. Note that the proper definition is that the list conforms to
    // 'International Tables for Crystallography, Volume A: Space Group Symmetry: Space Group Symmetry v.A'. However, common cif-files
    // use a different order of the elements.
    // Option (c) and (d) are considered last resorts, because they can only set space groups with a default origin setting.
    //
    // The input methods here use this preference. Also note that we allow for more flexible reading then the strict CIF-definition
    // Space groups should be listed with single spaces between the symmetry elements. But we also allow for more spaces or no spaces,
    // as well as space group names with underscores (common in older CIF-files).

    // read loop symmetry elements
    // ====================================
    if(strncasecmp(keyword,"_symmetry_equiv_pos",strlen("_symmetry_equiv_pos"))==0)
    {
      DefinedInLoop=TRUE;
      NumberOfAtomSiteElementsInBlock=0;
      nr_sg_elements=0;

      for(i=0;i<MAX_NUMBER_OF_CIF_TERMS;i++)
        CifBlockData[i]=-1;
      do
      {
        // parse the '_symmetry_equiv_pos' data names
        if(strncasecmp(keyword,"_symmetry_equiv_pos_site_id",strlen("_symmetry_equiv_pos_site_id"))==0)
          CifBlockData[SYMMETRY_EQUIV_POS_SITE_ID]=NumberOfAtomSiteElementsInBlock;
        else if(strncasecmp(keyword,"_symmetry_equiv_pos_as_xyz",strlen("_symmetry_equiv_pos_as_xyz"))==0)
        {
          CifBlockData[SYMMETRY_EQUIV_POS_AS_XYZ]=NumberOfAtomSiteElementsInBlock;
          TrimStringInPlace(arguments);
          if((strlen(arguments)>0)&&(arguments[0]!='_'))
          {
            DefinedInLoop=FALSE;
            RemoveQuotesAroundString(arguments);
            strcpy(SpaceGroup[nr_sg_elements++],arguments);
          }
          else
            NumberOfAtomSiteElementsInBlock++;
        }


        // check if another data name follows
        // if so, copy the 'arguments' to 'line'
        // otherwise read 'line' from the file
        TrimStringInPlace(arguments);
        if((strlen(arguments)>0)&&(arguments[0]=='_'))
          strcpy(line,arguments);
        else
          LoopCondition=fgets(line,MAX_CIF_LINE_SIZE,FilePtr);

        strcpy(keyword,"");
        strcpy(arguments,"");
        sscanf(line,"%s%[^\n]",keyword,arguments);

        // skip comments
        if(keyword[0]=='#')
        {
          strcpy(keyword,"");
          strcpy(arguments,"");
        }
      } // also continue when empty line
      while(((strncasecmp(keyword,"_symmetry_equiv_pos",strlen("_symmetry_equiv_pos"))==0)||(strlen(keyword)==0))&&LoopCondition);

      if(DefinedInLoop)
      {
        nr_sg_elements=0;
        do
        {
          // loop over data items
          if(strlen(keyword)>0)
          {
            arg_pointer=line;

            // remove spaces
            TrimStringInPlace(arg_pointer);

            for(i=0;i<NumberOfAtomSiteElementsInBlock;i++)
            {
              strcpy(keyword,"");
              strcpy(arguments,"");

              if(sscanf(arg_pointer,"%s%n",keyword,&n)==1)
              {
                arg_pointer+=n;
                if(keyword[0]=='\'')
                {
                  memmove(&keyword[0],&keyword[1],strlen(keyword));
                  if(keyword[strlen(keyword)-1]!='\'')
                  {
                    sscanf(arg_pointer,"%[^\']%n",arguments,&n);
                    arg_pointer+=n+1;
                    strcat(keyword, arguments);
                  }
                  else keyword[strlen(keyword)-1]='\0';
                }
              }
              else
              {
                fprintf(stderr, "Too few data items found in loop_ _symmetry_equiv_pos\n");
                exit(0);
              }

              if(i==CifBlockData[SYMMETRY_EQUIV_POS_SITE_ID])
              {
              }
              else if(i==CifBlockData[SYMMETRY_EQUIV_POS_AS_XYZ])
              {
                TrimStringInPlace(keyword);
                StripLeadingAndTrailingQuotesInPlace(keyword);
                StripWhiteSpacesInPlace(keyword);
                strcpy(SpaceGroup[nr_sg_elements++],keyword);
              }
              else
                fprintf(stderr, "unknown: %s\n",keyword);

            }
            if(sscanf(arg_pointer,"%s%n",keyword,&n)==1)
              fprintf(stderr, "Warning Too many data items found in loop_ _symmetry_equiv_pos, excess data ignored\n");
          }

          // split the new line into 'keyword' and 'arguments'
          LoopCondition=fgets(line,MAX_CIF_LINE_SIZE,FilePtr);
          strcpy(keyword,"");
          strcpy(arguments,"");
          sscanf(line,"%s%[^\n]",keyword,arguments);

          // skip comments
          if(keyword[0]=='#')
          {
            strcpy(keyword,"");
            strcpy(arguments,"");
          }
        } // continue as long as we do not encounter 'loop_' or a new data name and as long as we did not reach EOF
        while(((strncasecmp(keyword,"loop_",strlen("loop_"))!=0)&&(keyword[0]!='_'))&&LoopCondition);
      }
      FoundSpaceGroupSymmetryElements=GetSpacegroupFromSymmetryElements(nr_sg_elements,SpaceGroup);

      if (!STREAM)
      {
        fprintf(stderr, "space group found from symmetry elements: %d (nr elements: %d)\n",FoundSpaceGroupSymmetryElements,nr_sg_elements);
      }
    }

    // the Hall symbol is always read and overrules '_symmetry_space_group_name_H-M' and '_symmetry_equiv_pos_as_xyz'
    if(strncasecmp(keyword,"_symmetry_space_group_name_Hall",strlen("_symmetry_space_group_name_Hall"))==0)
    {
      TrimStringInPlace(arguments);                    // remove leading and trailing spaces
      StripLeadingAndTrailingQuotesInPlace(arguments); // strip possible leading and trailing quotes
      TrimStringInPlace(arguments);                    // remove leading and trailing spaces
      CompressSpacesInString(arguments);               // reduces multiple spaces to single spaces

      FoundSpaceGroupHall=GetSpaceGroupFromHallString(arguments);

      if(FoundSpaceGroupHall<1)  // if not found try with all spaces removed
        FoundSpaceGroupHall=GetSpaceGroupFromHallStringWithSpacesRemoved(arguments);

      if(FoundSpaceGroupHall<1)  // if not found try if spaces are '_' (common in older cif-files)
      {
        ReplaceCharacterInString(arguments,'_',' ');
        FoundSpaceGroupHall=GetSpaceGroupFromHallString(arguments);
      }

      fprintf(stderr, "_symmetry_space_group_name_Hall: %s found space group: %d\n",arguments,FoundSpaceGroupHall);
    }

    if((strncasecmp(keyword,"_symmetry_space_group_name_H-M",strlen("_symmetry_space_group_name_H-M"))==0)&&
       (Framework[CurrentSystem].SpaceGroupIdentifier[CurrentFramework]<1))
    {
      TrimStringInPlace(arguments);                    // remove leading and trailing spaces
      StripLeadingAndTrailingQuotesInPlace(arguments); // strip possible leading and trailing quotes
      TrimStringInPlace(arguments);                    // remove leading and trailing spaces
      CompressSpacesInString(arguments);               // reduces multiple spaces to single spaces

      FoundSpaceGroupHermannMauguinString=GetSpaceGroupFromHermannMauguinString(arguments);

      if(FoundSpaceGroupHermannMauguinString<1)  // if not found try with all spaces removed
        FoundSpaceGroupHermannMauguinString=GetSpaceGroupFromHermannMauguinStringWithSpacesRemoved(arguments);

      if(FoundSpaceGroupHermannMauguinString<1)  // if not found try if spaces are '_' (common in older cif-files)
      {
        ReplaceCharacterInString(arguments,'_',' ');
        FoundSpaceGroupHermannMauguinString=GetSpaceGroupFromHermannMauguinString(arguments);
      }

      fprintf(stderr, "_symmetry_space_group_name_H-M: %s found space group: %d\n",arguments,FoundSpaceGroupHermannMauguinString);
    }

    if((strncasecmp(keyword,"_symmetry_Int_Tables_number",strlen("_symmetry_Int_Tables_number"))==0)&&
       (Framework[CurrentSystem].SpaceGroupIdentifier[CurrentFramework]<1))
    {
      TrimStringInPlace(arguments);                    // remove leading and trailing spaces
      StripLeadingAndTrailingQuotesInPlace(arguments); // strip possible leading and trailing quotes
      TrimStringInPlace(arguments);                    // remove leading and trailing spaces
      sscanf(arguments,"%d",&FoundSpaceGroupIntTables);
      FoundSpaceGroupIntTables=GetSpaceGroupFromITSpaceGroupNumber(FoundSpaceGroupIntTables);

      fprintf(stderr, "_symmetry_Int_Tables_number: %d\n",FoundSpaceGroupIntTables);
    }

    if(strncasecmp(keyword,"_symmetry_cell_setting",strlen("_symmetry_cell_setting"))==0)
    {
      // not necessary to read
      // fprintf(stderr, "_symmetry_cell_setting: %s\n",arguments);
    }

    if(strncasecmp(keyword,"_space_group.IT_coordinate_system_code",strlen("_space_group.IT_coordinate_system_code"))==0)
    {
      TrimStringInPlace(arguments);                    // remove leading and trailing spaces
      StripLeadingAndTrailingQuotesInPlace(arguments); // strip possible leading and trailing quotes
      TrimStringInPlace(arguments);                    // remove leading and trailing spaces
      sscanf(arguments,"%s",FoundSpaceGroupOption);
      fprintf(stderr, "_space_group.IT_coordinate_system_code: %s\n",FoundSpaceGroupOption);
    }

    if(strncasecmp(keyword,"_cell_length_a",strlen("_cell_length_a"))==0)
    {
      if(sscanf(arguments,"%lf",&A)!=1)
      {
        fprintf(stderr, "ERROR in cif-file: cell length A not recognized\n");
        exit(0);
      }
      if (!STREAM)
      {
        fprintf(stderr, "_cell_length_a: %lf\n",A);
      }
    }

    if(strncasecmp(keyword,"_cell_length_b",strlen("_cell_length_b"))==0)
    {
      if(sscanf(arguments,"%lf",&B)!=1)
      {
        fprintf(stderr, "ERROR in cif-file: cell length B not recognized\n");
        exit(0);
      }
      if (!STREAM)
      {
        fprintf(stderr, "_cell_length_b: %lf\n",B);
      }
    }

    if(strncasecmp(keyword,"_cell_length_c",strlen("_cell_length_c"))==0)
    {
      if(sscanf(arguments,"%lf",&C)!=1)
      {
        fprintf(stderr, "ERROR in cif-file: cell length C not recognized\n");
        exit(0);
      }
      if (!STREAM)
      {
        fprintf(stderr, "_cell_length_c: %lf\n",C);
      }
    }

    if(strncasecmp(keyword,"_cell_angle_alpha",strlen("_cell_angle_alpha"))==0)
    {
      if(sscanf(arguments,"%lf",&alpha)!=1)
      {
        fprintf(stderr, "ERROR in cif-file: cell angle alpha not recognized\n");
        exit(0);
      }
      if (!STREAM)
      {
        fprintf(stderr, "_cell_length_alpha: %lf\n",alpha);
      }
    }

    if(strncasecmp(keyword,"_cell_angle_beta",strlen("_cell_angle_beta"))==0)
    {
      if(sscanf(arguments,"%lf",&beta)!=1)
      {
        fprintf(stderr, "ERROR in cif-file: cell angle beta not recognized\n");
        exit(0);
      }
      if (!STREAM)
      {
        fprintf(stderr, "_cell_length_beta: %lf\n",beta);
      }
    }

    if(strncasecmp(keyword,"_cell_angle_gamma",strlen("_cell_angle_gamma"))==0)
    {
      if(sscanf(arguments,"%lf",&gamma)!=1)
      {
        fprintf(stderr, "ERROR in cif-file: cell angle gamma not recognized\n");
        exit(0);
      }
      if (!STREAM)
      {
        fprintf(stderr, "_cell_length_gamma: %lf\n",gamma);
      }
    }

    TrimStringInPlace(arguments);
    if((strlen(arguments)>0)&&(arguments[0]=='_'))
      strcpy(line,arguments);
    else
      LoopCondition=fgets(line,MAX_CIF_LINE_SIZE,FilePtr);

    strcpy(keyword,"");
    strcpy(arguments,"");
    sscanf(line,"%s%[^\n]",keyword,arguments);

    // skip comments
    if(keyword[0]=='#')
    {
      strcpy(keyword,"");
      strcpy(arguments,"");
    }
  }
  while(LoopCondition);

  fclose(FilePtr);

  // Preference of finding spacegroups: Hall > Symmetry elements > HM symbol > IT number
  if(FoundSpaceGroupHall>=1)
    Framework[CurrentSystem].SpaceGroupIdentifier[CurrentFramework]=FoundSpaceGroupHall;
  else if(FoundSpaceGroupSymmetryElements>=1)
    Framework[CurrentSystem].SpaceGroupIdentifier[CurrentFramework]=FoundSpaceGroupSymmetryElements;
  else
  {
    if(FoundSpaceGroupHermannMauguinString>=1)
      Framework[CurrentSystem].SpaceGroupIdentifier[CurrentFramework]=FoundSpaceGroupHermannMauguinString;
    else if(FoundSpaceGroupIntTables>=1)
      Framework[CurrentSystem].SpaceGroupIdentifier[CurrentFramework]=FoundSpaceGroupIntTables;
    else
    {
      fprintf(stderr, "ERROR in cif-file: no proper space group definition found\n");
      exit(0);
    }

    if(strlen(FoundSpaceGroupOption)>0)
    {
      Framework[CurrentSystem].SpaceGroupIdentifier[CurrentFramework]=AdjustForNonStandardCIFOption(Framework[CurrentSystem].SpaceGroupIdentifier[CurrentFramework],FoundSpaceGroupOption);
      fprintf(stderr, "After adjustement: %d\n",Framework[CurrentSystem].SpaceGroupIdentifier[CurrentFramework]);
    }
  }

  // convert angles from degrees to radians
  // determine boundary conditions from angles
  if(BoundaryCondition[CurrentSystem]==UNINITIALIZED_BOUNDARY_CONDITION)
  {
    if((fabs(alpha-90.0)>0.001)||(fabs(beta-90.0)>0.001)||(fabs(gamma-90.0)>0.001))
      BoundaryCondition[CurrentSystem]=TRICLINIC;
    else
    {
      if((fabs(A-B)>0.00001)||(fabs(A-C)>0.00001)||(fabs(B-C)>0.00001))
        BoundaryCondition[CurrentSystem]=RECTANGULAR;
      else BoundaryCondition[CurrentSystem]=TRICLINIC;
    }
  }
  AlphaAngle[CurrentSystem]=alpha*DEG2RAD;
  BetaAngle[CurrentSystem]=beta*DEG2RAD;
  GammaAngle[CurrentSystem]=gamma*DEG2RAD;

  UnitCellSize[CurrentSystem].x=A;
  UnitCellSize[CurrentSystem].y=B;
  UnitCellSize[CurrentSystem].z=C;

  // first vector
  UnitCellBox[CurrentSystem].ax=A;
  UnitCellBox[CurrentSystem].ay=0.0;
  UnitCellBox[CurrentSystem].az=0.0;

  // second vector
  UnitCellBox[CurrentSystem].bx=B*cos(GammaAngle[CurrentSystem]);
  UnitCellBox[CurrentSystem].by=B*sin(GammaAngle[CurrentSystem]);
  UnitCellBox[CurrentSystem].bz=0.0;

  // third vector
  tempd=(cos(AlphaAngle[CurrentSystem])-cos(GammaAngle[CurrentSystem])*cos(BetaAngle[CurrentSystem]))/sin(GammaAngle[CurrentSystem]);
  UnitCellBox[CurrentSystem].cx=C*cos(BetaAngle[CurrentSystem]);
  UnitCellBox[CurrentSystem].cy=C*tempd;
  UnitCellBox[CurrentSystem].cz=C*sqrt(1.0-SQR(cos(BetaAngle[CurrentSystem]))-SQR(tempd));

  Invert3x3Matrix(&UnitCellBox[CurrentSystem],&InverseUnitCellBox[CurrentSystem],&det);

  A=(REAL)NumberOfUnitCells[CurrentSystem].x*UnitCellSize[CurrentSystem].x;
  B=(REAL)NumberOfUnitCells[CurrentSystem].y*UnitCellSize[CurrentSystem].y;
  C=(REAL)NumberOfUnitCells[CurrentSystem].z*UnitCellSize[CurrentSystem].z;
  tempd=(cos(AlphaAngle[CurrentSystem])-cos(GammaAngle[CurrentSystem])*cos(BetaAngle[CurrentSystem]))/sin(GammaAngle[CurrentSystem]);

  // first vector
  Box[CurrentSystem].ax=A;
  Box[CurrentSystem].ay=0.0;
  Box[CurrentSystem].az=0.0;

  // second vector
  Box[CurrentSystem].bx=B*cos(GammaAngle[CurrentSystem]);
  Box[CurrentSystem].by=B*sin(GammaAngle[CurrentSystem]);
  Box[CurrentSystem].bz=0.0;

  // third vector
  Box[CurrentSystem].cx=C*cos(BetaAngle[CurrentSystem]);
  Box[CurrentSystem].cy=C*tempd;
  Box[CurrentSystem].cz=C*sqrt(1.0-SQR(cos(BetaAngle[CurrentSystem]))-SQR(tempd));

  Invert3x3Matrix(&Box[CurrentSystem],&InverseBox[CurrentSystem],&det);

  // calculate box-properties
  CellProperties(&InverseBox[CurrentSystem],&InverseBoxProperties[CurrentSystem],&det);
  CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);

  // no charge definition found in the CIF-file, so UseChargesFromCIFFile must be set to false
  if(!BoolChargeDefinedInCIFFile)
    UseChargesFromCIFFile=FALSE;

  // expand asymmetric to full framework
  ExpandAsymmetricFrameworkToFullFramework();

  if (!STREAM)
  {
    fprintf(stderr, "End reading cif-file\n");
  }
}

void WriteFrameworkDefinitionShell(char * string)
{
  int i,j;
  int Type;
  VECTOR pos;
  char symbol[256];
  char buffer[256],name[256];
  FILE *FilePtr;

  mkdir("Movies",S_IRWXU);
  for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
  {
    sprintf(buffer,"Movies/System_%d",CurrentSystem);
    mkdir(buffer,S_IRWXU);

    for(CurrentFramework=0;CurrentFramework<Framework[CurrentSystem].NumberOfFrameworks;CurrentFramework++)
    {

      sprintf(buffer,"Movies/System_%d/Framework_%d_%s_%d_%d_%d%s.cell",
              CurrentSystem,CurrentFramework,string,
              NumberOfUnitCells[CurrentSystem].x,NumberOfUnitCells[CurrentSystem].y,NumberOfUnitCells[CurrentSystem].z,
              FileNameAppend);

      FilePtr=fopen(buffer,"w");

      fprintf(FilePtr,"%%BLOCK LATTICE_CART\n");
      fprintf(FilePtr,"ang    # angstrom units\n");
      fprintf(FilePtr," % -18.12f % -18.12f % -18.12f\n",UnitCellBox[CurrentSystem].ax,UnitCellBox[CurrentSystem].ay,UnitCellBox[CurrentSystem].az);
      fprintf(FilePtr," % -18.12f % -18.12f % -18.12f\n",UnitCellBox[CurrentSystem].bx,UnitCellBox[CurrentSystem].by,UnitCellBox[CurrentSystem].bz);
      fprintf(FilePtr," % -18.12f % -18.12f % -18.12f\n",UnitCellBox[CurrentSystem].cx,UnitCellBox[CurrentSystem].cy,UnitCellBox[CurrentSystem].cz);
      fprintf(FilePtr,"%%ENDBLOCK LATTICE_CART\n");
      fprintf(FilePtr,"\n");
      fprintf(FilePtr,"%%BLOCK POSITIONS_FRAC\n");

      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];i++)
      {
        // convert position from Cartesian to fractional positions
        pos.x=Framework[CurrentSystem].Atoms[CurrentFramework][i].Position.x;
        pos.y=Framework[CurrentSystem].Atoms[CurrentFramework][i].Position.y;
        pos.z=Framework[CurrentSystem].Atoms[CurrentFramework][i].Position.z;
        pos=ConvertFromXYZtoABC(pos);

        Type=Framework[CurrentSystem].Atoms[CurrentFramework][i].Type;

        strcpy(symbol,PseudoAtoms[Type].ChemicalElement);

        fprintf(FilePtr,"%-5s % -18.12f % -18.12f % -18.12f\n",
                symbol,
                fabs(pos.x)<1e-12?0.0:pos.x,
                fabs(pos.y)<1e-12?0.0:pos.y,
                fabs(pos.z)<1e-12?0.0:pos.z);
      }

      fprintf(FilePtr,"%%ENDBLOCK POSITIONS_FRAC\n");
      fclose(FilePtr);


      sprintf(buffer,"Movies/System_%d/Framework_%d_%s_%d_%d_%d%s.txt",
              CurrentSystem,CurrentFramework,string,
              NumberOfUnitCells[CurrentSystem].x,NumberOfUnitCells[CurrentSystem].y,NumberOfUnitCells[CurrentSystem].z,
              FileNameAppend);

      FilePtr=fopen(buffer,"w");
      fprintf(FilePtr,"[double3(%14.12f,%14.12f,%14.12f),double3(%14.12f,%14.12f,%14.12f),double3(%14.12f,%14.12f,%14.12f)]\n",
           UnitCellBox[CurrentSystem].ax,UnitCellBox[CurrentSystem].ay,UnitCellBox[CurrentSystem].az,
           UnitCellBox[CurrentSystem].bx,UnitCellBox[CurrentSystem].by,UnitCellBox[CurrentSystem].bz,
           UnitCellBox[CurrentSystem].cx,UnitCellBox[CurrentSystem].cy,UnitCellBox[CurrentSystem].cz);
      fprintf(FilePtr,"\n");

      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];i++)
      {
        // convert position from Cartesian to fractional positions
        pos.x=Framework[CurrentSystem].Atoms[CurrentFramework][i].Position.x;
        pos.y=Framework[CurrentSystem].Atoms[CurrentFramework][i].Position.y;
        pos.z=Framework[CurrentSystem].Atoms[CurrentFramework][i].Position.z;
        pos=ConvertFromXYZtoABC(pos);

        fprintf(FilePtr,"double3(%13.12f,%13.12f,%13.12f),\n",
                fabs(pos.x)<1e-12?0.0:pos.x,
                fabs(pos.y)<1e-12?0.0:pos.y,
                fabs(pos.z)<1e-12?0.0:pos.z);
      }

      fclose(FilePtr);
    }
  }
}

/*********************************************************************************************************
 * Name       | WriteFrameworkDefinitionCIF                                                              *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | writes a CIF-file                                                                        *
 * Parameters | -                                                                                        *
 * Note       | The CIF specification can be found online: http://www.iucr.org/resources/cif/spec        *
 *            | and S.R. Hall, F.H. Allen, and I.D. Brown, Acta Cryst. A47, 655-685, 1991.               *
 *            | The definition is based on the free format STAR procedure.                               *
 *            | The written CIF-files is readable by e.g. Materials Studio, jmol, and crystalmaker.      *
 *********************************************************************************************************/

void WriteFrameworkDefinitionCIF(char * string)
{
  int i,j,n;
  REAL A,B,C;
  VECTOR pos;
  FILE *FilePtr;
  char buffer[256],name[256];
  int Type;
  time_t curtime;
  struct tm *loctime;
  int *AtomIdentifier;
  char fullname[256];
  char symbol[256];
  struct passwd *p;
  VECTOR flexible_drift;
  VECTOR com;

  if (STREAM)
  {
    fprintf(stderr, "File writing not allowed in streaming mode!");
    return;
  }

  AtomIdentifier=(int*)calloc(NumberOfPseudoAtoms,sizeof(int));

  curtime = time (NULL);
  loctime = localtime (&curtime);

  for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
  {
    A=(REAL)UnitCellSize[CurrentSystem].x;
    B=(REAL)UnitCellSize[CurrentSystem].y;
    C=(REAL)UnitCellSize[CurrentSystem].z;

    flexible_drift.x=flexible_drift.y=flexible_drift.z=0.0;
    if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
    {
      com=GetFrameworkCenterOfMass();
      flexible_drift.x=com.x-Framework[CurrentSystem].IntialCenterOfMassPosition.x;
      flexible_drift.y=com.y-Framework[CurrentSystem].IntialCenterOfMassPosition.y;
      flexible_drift.z=com.z-Framework[CurrentSystem].IntialCenterOfMassPosition.z;
    }

    for(CurrentFramework=0;CurrentFramework<Framework[CurrentSystem].NumberOfFrameworks;CurrentFramework++)
    {
      sprintf(buffer,"Movies/System_%d",CurrentSystem);
      mkdir(buffer,S_IRWXU);

      sprintf(buffer,"Movies/System_%d/Framework_%d_%s_%d_%d_%d_VASP%s.cif",
              CurrentSystem,CurrentFramework,string,
              NumberOfUnitCells[CurrentSystem].x,NumberOfUnitCells[CurrentSystem].y,NumberOfUnitCells[CurrentSystem].z,
              FileNameAppend);

      FilePtr=fopen(buffer,"w");

      fprintf(FilePtr,"data_%s\n\n",Framework[CurrentSystem].Name[CurrentFramework]);

      fprintf(FilePtr,"_audit_creation_method RASPA-1.0\n");
      fprintf(FilePtr,"_audit_creation_date %d-%d-%d\n",loctime->tm_year+1900,loctime->tm_mon+1,loctime->tm_mday);


      #if defined(__INTEL_COMPILER)
        getlogin_r(fullname,256);
        fprintf(FilePtr,"_audit_author_name '%s'\n",fullname);
      #else
        p=getpwuid(geteuid());
        if(p)
        {
          n=strcspn(p->pw_gecos,",");
          memcpy(fullname,p->pw_gecos,n);
          fullname[n] = '\0';
          fprintf(FilePtr,"_audit_author_name '%s'\n",fullname);
        }
      #endif
      fprintf(FilePtr,"\n");

      if(Framework[CurrentSystem].NumberOfCitations[CurrentFramework]>0)
      {
        //for(i=0;i<Framework[CurrentSystem].NumberOfCitations[CurrentFramework];i++)
        if(Framework[CurrentSystem].NumberOfCitations[CurrentFramework]==1)
        {
          if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationId)>0)
            fprintf(FilePtr,"_citation_id                 %s\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationId);
          if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationAuthorName)>0)
            fprintf(FilePtr,"_citation_author_name        '%s'\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationAuthorName);
          if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationCoordinateLinkage)>0)
            fprintf(FilePtr,"_citation_coordinate_linkage %s\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationCoordinateLinkage);
          if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationTitle)>0)
            fprintf(FilePtr,"_citation_title              '%s'\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationTitle);
          if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationCountry)>0)
            fprintf(FilePtr,"_citation_country            '%s'\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationCountry);
          if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationJournalAbbrev)>0)
            fprintf(FilePtr,"_citation_journal_abbrev     '%s'\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationJournalAbbrev);
          if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationJournalVolume)>0)
            fprintf(FilePtr,"_citation_journal_volume     %s\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationJournalVolume);
          if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationJournalIssue)>0)
            fprintf(FilePtr,"_citation_journal_issue      %s\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationJournalIssue);
          if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationPageFirst)>0)
            fprintf(FilePtr,"_citation_page_first         %s\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationPageFirst);
          if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationPageLast)>0)
            fprintf(FilePtr,"_citation_page_last          %s\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationPageLast);
          if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationYear)>0)
            fprintf(FilePtr,"_citation_year               %s\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationYear);
          if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationJournalID_ASTM)>0)
            fprintf(FilePtr,"_citation_journal_ID_ASTM    %s\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationJournalID_ASTM);
          if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationJournalID_ISSN)>0)
            fprintf(FilePtr,"_citation_journal_ID_ISSN    %s\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationJournalID_ISSN);
          if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationBookTitle)>0)
            fprintf(FilePtr,"_citation_book_title         '%s'\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationBookTitle);
          if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationBookPublisher)>0)
            fprintf(FilePtr,"_citation_book_publisher     '%s'\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationBookPublisher);
          if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationBookID_ISBN)>0)
            fprintf(FilePtr,"_citation_book_ID_ISBN       %s\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationBookID_ISBN);
          if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationSpecialDetails)>0)
            fprintf(FilePtr,"_citation_special_details    '%s'\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationSpecialDetails);
          fprintf(FilePtr,"\n");
        }
        else
        {
          fprintf(FilePtr,"loop_\n");

          if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationId)>0) fprintf(FilePtr,"_citation_id\n");
          if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationAuthorName)>0) fprintf(FilePtr,"_citation_author_name\n");
          if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationCoordinateLinkage)>0) fprintf(FilePtr,"_citation_coordinate_linkage\n");
          if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationTitle)>0) fprintf(FilePtr,"_citation_title\n");
          if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationCountry)>0) fprintf(FilePtr,"_citation_country\n");
          if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationJournalAbbrev)>0) fprintf(FilePtr,"_citation_journal_abbrev\n");
          if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationJournalVolume)>0) fprintf(FilePtr,"_citation_journal_volume\n");
          if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationJournalIssue)>0) fprintf(FilePtr,"_citation_journal_issue\n");
          if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationPageFirst)>0) fprintf(FilePtr,"_citation_page_first\n");
          if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationPageLast)>0) fprintf(FilePtr,"_citation_page_last\n");
          if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationYear)>0) fprintf(FilePtr,"_citation_year\n");
          if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationJournalID_ASTM)>0) fprintf(FilePtr,"_citation_journal_ID_ASTM\n");
          if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationJournalID_ISSN)>0) fprintf(FilePtr,"_citation_journal_ID_ISSN\n");
          if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationBookTitle)>0) fprintf(FilePtr,"_citation_book_title\n");
          if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationBookPublisher)>0) fprintf(FilePtr,"_citation_book_publisher\n");
          if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationBookID_ISBN)>0) fprintf(FilePtr,"_citation_book_ID_ISBN\n");
          if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationSpecialDetails)>0) fprintf(FilePtr,"_citation_special_details\n");

          for(i=0;i<Framework[CurrentSystem].NumberOfCitations[CurrentFramework];i++)
          {
            if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationId)>0)
              fprintf(FilePtr,"%s\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationId);
            if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationAuthorName)>0)
              fprintf(FilePtr,"'%s'\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationAuthorName);
            if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationCoordinateLinkage)>0)
              fprintf(FilePtr,"%s\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationCoordinateLinkage);
            if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationTitle)>0)
              fprintf(FilePtr,"'%s'\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationTitle);
            if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationCountry)>0)
              fprintf(FilePtr,"'%s'\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationCountry);
            if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationJournalAbbrev)>0)
              fprintf(FilePtr,"'%s'\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationJournalAbbrev);
            if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationJournalVolume)>0)
              fprintf(FilePtr,"%s\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationJournalVolume);
            if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationJournalIssue)>0)
              fprintf(FilePtr,"%s\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationJournalIssue);
            if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationPageFirst)>0)
              fprintf(FilePtr,"%s\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationPageFirst);
            if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationPageLast)>0)
              fprintf(FilePtr,"%s\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationPageLast);
            if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationYear)>0)
              fprintf(FilePtr,"%s\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationYear);
            if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationJournalID_ASTM)>0)
              fprintf(FilePtr,"%s\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationJournalID_ASTM);
            if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationJournalID_ISSN)>0)
              fprintf(FilePtr,"%s\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationJournalID_ISSN);
            if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationBookTitle)>0)
              fprintf(FilePtr,"'%s'\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationBookTitle);
            if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationBookPublisher)>0)
              fprintf(FilePtr,"'%s'\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationBookPublisher);
            if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationBookID_ISBN)>0)
              fprintf(FilePtr,"%s\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationBookID_ISBN);
            if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationSpecialDetails)>0)
              fprintf(FilePtr,"'%s'\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationSpecialDetails);
          }
          fprintf(FilePtr,"\n");
        }
      }


      fprintf(FilePtr,"_cell_length_a    %g\n",BoxProperties[CurrentSystem].ax/(REAL)NumberOfUnitCells[CurrentSystem].x);
      fprintf(FilePtr,"_cell_length_b    %g\n",BoxProperties[CurrentSystem].ay/(REAL)NumberOfUnitCells[CurrentSystem].y);
      fprintf(FilePtr,"_cell_length_c    %g\n",BoxProperties[CurrentSystem].az/(REAL)NumberOfUnitCells[CurrentSystem].z);

      fprintf(FilePtr,"_cell_angle_alpha %g\n",(REAL)AlphaAngle[CurrentSystem]*RAD2DEG);
      fprintf(FilePtr,"_cell_angle_beta  %g\n",(REAL)BetaAngle[CurrentSystem]*RAD2DEG);
      fprintf(FilePtr,"_cell_angle_gamma %g\n",(REAL)GammaAngle[CurrentSystem]*RAD2DEG);
      fprintf(FilePtr,"_cell_volume      %g\n\n",(REAL)Volume[CurrentSystem]);

      switch(SpaceGroupData[Framework[CurrentSystem].SpaceGroupIdentifier[CurrentFramework]].CrystalSystem)
      {
        case TRICLINIC_SPACEGROUP:
          fprintf(FilePtr,"_symmetry_cell_setting          triclinic\n");
          break;
        case MONOCLINIC_SPACEGROUP:
          fprintf(FilePtr,"_symmetry_cell_setting          monoclinic\n");
          break;
        case ORTHORHOMBIC_SPACEGROUP:
          fprintf(FilePtr,"_symmetry_cell_setting          orthorhombic\n");
          break;
        case TETRAGONAL_SPACEGROUP:
          fprintf(FilePtr,"_symmetry_cell_setting          tetragonal\n");
          break;
        case TRIGONAL_SPACEGROUP:
          fprintf(FilePtr,"_symmetry_cell_setting          trigonal\n");
          break;
        case HEXAGONAL_SPACEGROUP:
          fprintf(FilePtr,"_symmetry_cell_setting          hexagonal\n");
          break;
        case CUBIC_SPACEGROUP:
          fprintf(FilePtr,"_symmetry_cell_setting          cubic\n");
          break;
      }
      fprintf(FilePtr,"_symmetry_space_group_name_Hall '%s'\n",SpaceGroupData[1].HallSpaceGroupSymbol);
      fprintf(FilePtr,"_symmetry_space_group_name_H-M  '%s'\n",SpaceGroupData[1].ShortInternationalHermannMauguinSpaceGroupSymbol);
      fprintf(FilePtr,"_symmetry_Int_Tables_number     %d\n\n",SpaceGroupData[1].Number);


      fprintf(FilePtr,"_symmetry_equiv_pos_as_xyz 'x,y,z'\n");
      fprintf(FilePtr,"\n");


      fprintf(FilePtr,"loop_\n");
      fprintf(FilePtr,"_atom_site_label\n");
      fprintf(FilePtr,"_atom_site_type_symbol\n");
      fprintf(FilePtr,"_atom_site_fract_x\n");
      fprintf(FilePtr,"_atom_site_fract_y\n");
      fprintf(FilePtr,"_atom_site_fract_z\n");
      fprintf(FilePtr,"_atom_site_charge\n");

      for(i=0;i<NumberOfPseudoAtoms;i++)
        AtomIdentifier[i]=0;

      for(i=0;i<NumberOfPseudoAtoms;i++)
      {
        for(j=0;j<Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];j++)
        {
          if((Framework[CurrentSystem].Atoms[CurrentFramework][j].Type==i)&&(PseudoAtoms[Framework[CurrentSystem].Atoms[CurrentFramework][j].Type].PrintToPDB))
          {
            pos=ConvertFromXYZtoABC(Framework[CurrentSystem].Atoms[CurrentFramework][j].Position);
            pos.x=Framework[CurrentSystem].Atoms[CurrentFramework][j].Position.x-flexible_drift.x;
            pos.y=Framework[CurrentSystem].Atoms[CurrentFramework][j].Position.y-flexible_drift.y;
            pos.z=Framework[CurrentSystem].Atoms[CurrentFramework][j].Position.z-flexible_drift.z;
            pos=ConvertFromXYZtoABC(pos);

            if(Framework[CurrentSystem].AddAtomNumberCodeToLabel[CurrentFramework])
            {
              AtomIdentifier[i]++;
              sprintf(name,"%s%d",PseudoAtoms[i].Name,AtomIdentifier[i]);
            }
            else
              sprintf(name,"%s",PseudoAtoms[i].Name);

            strcpy(symbol,PseudoAtoms[i].ChemicalElement);
            if(PseudoAtoms[i].OxidationState!=0)
              strcat(symbol,PseudoAtoms[i].OxidationStateString);

            if(fabs(pos.x)<1e-10) pos.x=fabs(pos.x);
            if(fabs(pos.y)<1e-10) pos.y=fabs(pos.y);
            if(fabs(pos.z)<1e-10) pos.z=fabs(pos.z);

            fprintf(FilePtr,"%-8s %-5s % -18.12f % -18.12f % -18.12f % -12g\n",
                    name,
                    symbol,
                    fabs(pos.x)<1e-12?0.0:pos.x,
                    fabs(pos.y)<1e-12?0.0:pos.y,
                    fabs(pos.z)<1e-12?0.0:pos.z,
                    Framework[CurrentSystem].Atoms[CurrentFramework][i].Charge);
          }
        }
      }
      fprintf(FilePtr,"\n\n");
      fclose(FilePtr);

      sprintf(buffer,"Movies/System_%d/Framework_%d_%s_%d_%d_%d_P1%s.cif",
              CurrentSystem,CurrentFramework,string,
              NumberOfUnitCells[CurrentSystem].x,NumberOfUnitCells[CurrentSystem].y,NumberOfUnitCells[CurrentSystem].z,
              FileNameAppend);

      FilePtr=fopen(buffer,"w");

      fprintf(FilePtr,"data_%s\n\n",Framework[CurrentSystem].Name[CurrentFramework]);

      fprintf(FilePtr,"_audit_creation_method RASPA-1.0\n");
      fprintf(FilePtr,"_audit_creation_date %d-%d-%d\n",loctime->tm_year+1900,loctime->tm_mon+1,loctime->tm_mday);
      #if defined(__INTEL_COMPILER)
        getlogin_r(fullname,256);
        fprintf(FilePtr,"_audit_author_name '%s'\n",fullname);
      #else
        p=getpwuid(geteuid());
        if(p)
        {
          n=strcspn(p->pw_gecos,",");
          memcpy(fullname,p->pw_gecos,n);
          fullname[n] = '\0';
          fprintf(FilePtr,"_audit_author_name '%s'\n",fullname);
        }
      #endif
      fprintf(FilePtr,"\n");

      if(Framework[CurrentSystem].NumberOfCitations[CurrentFramework]>0)
      {
        //for(i=0;i<Framework[CurrentSystem].NumberOfCitations[CurrentFramework];i++)
        if(Framework[CurrentSystem].NumberOfCitations[CurrentFramework]==1)
        {
          if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationId)>0)
            fprintf(FilePtr,"_citation_id                 %s\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationId);
          if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationAuthorName)>0)
            fprintf(FilePtr,"_citation_author_name        '%s'\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationAuthorName);
          if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationCoordinateLinkage)>0)
            fprintf(FilePtr,"_citation_coordinate_linkage %s\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationCoordinateLinkage);
          if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationTitle)>0)
            fprintf(FilePtr,"_citation_title              '%s'\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationTitle);
          if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationCountry)>0)
            fprintf(FilePtr,"_citation_country            '%s'\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationCountry);
          if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationJournalAbbrev)>0)
            fprintf(FilePtr,"_citation_journal_abbrev     '%s'\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationJournalAbbrev);
          if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationJournalVolume)>0)
            fprintf(FilePtr,"_citation_journal_volume     %s\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationJournalVolume);
          if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationJournalIssue)>0)
            fprintf(FilePtr,"_citation_journal_issue      %s\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationJournalIssue);
          if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationPageFirst)>0)
            fprintf(FilePtr,"_citation_page_first         %s\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationPageFirst);
          if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationPageLast)>0)
            fprintf(FilePtr,"_citation_page_last          %s\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationPageLast);
          if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationYear)>0)
            fprintf(FilePtr,"_citation_year               %s\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationYear);
          if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationJournalID_ASTM)>0)
            fprintf(FilePtr,"_citation_journal_ID_ASTM    %s\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationJournalID_ASTM);
          if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationJournalID_ISSN)>0)
            fprintf(FilePtr,"_citation_journal_ID_ISSN    %s\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationJournalID_ISSN);
          if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationBookTitle)>0)
            fprintf(FilePtr,"_citation_book_title         '%s'\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationBookTitle);
          if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationBookPublisher)>0)
            fprintf(FilePtr,"_citation_book_publisher     '%s'\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationBookPublisher);
          if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationBookID_ISBN)>0)
            fprintf(FilePtr,"_citation_book_ID_ISBN       %s\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationBookID_ISBN);
          if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationSpecialDetails)>0)
            fprintf(FilePtr,"_citation_special_details    '%s'\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationSpecialDetails);
          fprintf(FilePtr,"\n");
        }
        else
        {
          fprintf(FilePtr,"loop_\n");

          if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationId)>0) fprintf(FilePtr,"_citation_id\n");
          if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationAuthorName)>0) fprintf(FilePtr,"_citation_author_name\n");
          if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationCoordinateLinkage)>0) fprintf(FilePtr,"_citation_coordinate_linkage\n");
          if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationTitle)>0) fprintf(FilePtr,"_citation_title\n");
          if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationCountry)>0) fprintf(FilePtr,"_citation_country\n");
          if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationJournalAbbrev)>0) fprintf(FilePtr,"_citation_journal_abbrev\n");
          if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationJournalVolume)>0) fprintf(FilePtr,"_citation_journal_volume\n");
          if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationJournalIssue)>0) fprintf(FilePtr,"_citation_journal_issue\n");
          if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationPageFirst)>0) fprintf(FilePtr,"_citation_page_first\n");
          if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationPageLast)>0) fprintf(FilePtr,"_citation_page_last\n");
          if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationYear)>0) fprintf(FilePtr,"_citation_year\n");
          if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationJournalID_ASTM)>0) fprintf(FilePtr,"_citation_journal_ID_ASTM\n");
          if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationJournalID_ISSN)>0) fprintf(FilePtr,"_citation_journal_ID_ISSN\n");
          if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationBookTitle)>0) fprintf(FilePtr,"_citation_book_title\n");
          if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationBookPublisher)>0) fprintf(FilePtr,"_citation_book_publisher\n");
          if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationBookID_ISBN)>0) fprintf(FilePtr,"_citation_book_ID_ISBN\n");
          if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationSpecialDetails)>0) fprintf(FilePtr,"_citation_special_details\n");

          for(i=0;i<Framework[CurrentSystem].NumberOfCitations[CurrentFramework];i++)
          {
            if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationId)>0)
              fprintf(FilePtr,"%s\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationId);
            if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationAuthorName)>0)
              fprintf(FilePtr,"'%s'\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationAuthorName);
            if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationCoordinateLinkage)>0)
              fprintf(FilePtr,"%s\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationCoordinateLinkage);
            if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationTitle)>0)
              fprintf(FilePtr,"'%s'\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationTitle);
            if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationCountry)>0)
              fprintf(FilePtr,"'%s'\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationCountry);
            if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationJournalAbbrev)>0)
              fprintf(FilePtr,"'%s'\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationJournalAbbrev);
            if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationJournalVolume)>0)
              fprintf(FilePtr,"%s\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationJournalVolume);
            if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationJournalIssue)>0)
              fprintf(FilePtr,"%s\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationJournalIssue);
            if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationPageFirst)>0)
              fprintf(FilePtr,"%s\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationPageFirst);
            if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationPageLast)>0)
              fprintf(FilePtr,"%s\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationPageLast);
            if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationYear)>0)
              fprintf(FilePtr,"%s\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationYear);
            if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationJournalID_ASTM)>0)
              fprintf(FilePtr,"%s\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationJournalID_ASTM);
            if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationJournalID_ISSN)>0)
              fprintf(FilePtr,"%s\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationJournalID_ISSN);
            if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationBookTitle)>0)
              fprintf(FilePtr,"'%s'\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationBookTitle);
            if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationBookPublisher)>0)
              fprintf(FilePtr,"'%s'\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationBookPublisher);
            if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationBookID_ISBN)>0)
              fprintf(FilePtr,"%s\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationBookID_ISBN);
            if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationSpecialDetails)>0)
              fprintf(FilePtr,"'%s'\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationSpecialDetails);
          }
          fprintf(FilePtr,"\n");
        }
      }


      fprintf(FilePtr,"_cell_length_a    %g\n",BoxProperties[CurrentSystem].ax);
      fprintf(FilePtr,"_cell_length_b    %g\n",BoxProperties[CurrentSystem].ay);
      fprintf(FilePtr,"_cell_length_c    %g\n",BoxProperties[CurrentSystem].az);
      fprintf(FilePtr,"_cell_angle_alpha %g\n",(REAL)AlphaAngle[CurrentSystem]*RAD2DEG);
      fprintf(FilePtr,"_cell_angle_beta  %g\n",(REAL)BetaAngle[CurrentSystem]*RAD2DEG);
      fprintf(FilePtr,"_cell_angle_gamma %g\n",(REAL)GammaAngle[CurrentSystem]*RAD2DEG);
      fprintf(FilePtr,"_cell_volume      %g\n\n",(REAL)Volume[CurrentSystem]);

      switch(SpaceGroupData[Framework[CurrentSystem].SpaceGroupIdentifier[CurrentFramework]].CrystalSystem)
      {
        case TRICLINIC_SPACEGROUP:
          fprintf(FilePtr,"_symmetry_cell_setting          triclinic\n");
          break;
        case MONOCLINIC_SPACEGROUP:
          fprintf(FilePtr,"_symmetry_cell_setting          monoclinic\n");
          break;
        case ORTHORHOMBIC_SPACEGROUP:
          fprintf(FilePtr,"_symmetry_cell_setting          orthorhombic\n");
          break;
        case TETRAGONAL_SPACEGROUP:
          fprintf(FilePtr,"_symmetry_cell_setting          tetragonal\n");
          break;
        case TRIGONAL_SPACEGROUP:
          fprintf(FilePtr,"_symmetry_cell_setting          trigonal\n");
          break;
        case HEXAGONAL_SPACEGROUP:
          fprintf(FilePtr,"_symmetry_cell_setting          hexagonal\n");
          break;
        case CUBIC_SPACEGROUP:
          fprintf(FilePtr,"_symmetry_cell_setting          cubic\n");
          break;
      }
      fprintf(FilePtr,"_symmetry_space_group_name_Hall '%s'\n",SpaceGroupData[1].HallSpaceGroupSymbol);
      fprintf(FilePtr,"_symmetry_space_group_name_H-M  '%s'\n",SpaceGroupData[1].ShortInternationalHermannMauguinSpaceGroupSymbol);
      fprintf(FilePtr,"_symmetry_Int_Tables_number     %d\n\n",SpaceGroupData[1].Number);


      fprintf(FilePtr,"_symmetry_equiv_pos_as_xyz 'x,y,z'\n");
      fprintf(FilePtr,"\n");


      fprintf(FilePtr,"loop_\n");
      fprintf(FilePtr,"_atom_site_label\n");
      fprintf(FilePtr,"_atom_site_type_symbol\n");
      fprintf(FilePtr,"_atom_site_fract_x\n");
      fprintf(FilePtr,"_atom_site_fract_y\n");
      fprintf(FilePtr,"_atom_site_fract_z\n");
      fprintf(FilePtr,"_atom_site_charge\n");

      for(i=0;i<NumberOfPseudoAtoms;i++)
        AtomIdentifier[i]=0;

      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];i++)
      {
        // convert position from Cartesian to fractional positions
        pos.x=Framework[CurrentSystem].Atoms[CurrentFramework][i].Position.x-flexible_drift.x;
        pos.y=Framework[CurrentSystem].Atoms[CurrentFramework][i].Position.y-flexible_drift.y;
        pos.z=Framework[CurrentSystem].Atoms[CurrentFramework][i].Position.z-flexible_drift.z;
        pos=ConvertFromXYZtoABC(pos);

        Type=Framework[CurrentSystem].Atoms[CurrentFramework][i].Type;

        if(Framework[CurrentSystem].AddAtomNumberCodeToLabel[CurrentFramework])
        {
          AtomIdentifier[Type]++;
          sprintf(name,"%s%d",PseudoAtoms[Type].Name,AtomIdentifier[Type]);
        }
        else
          sprintf(name,"%s",PseudoAtoms[Type].Name);

        strcpy(symbol,PseudoAtoms[Type].ChemicalElement);
        if(PseudoAtoms[Type].OxidationState!=0)
          strcat(symbol,PseudoAtoms[Type].OxidationStateString);

        fprintf(FilePtr,"%-8s %-5s % -18.12f % -18.12f % -18.12f % -12g\n",
                name,
                symbol,
                fabs(pos.x)<1e-12?0.0:pos.x,
                fabs(pos.y)<1e-12?0.0:pos.y,
                fabs(pos.z)<1e-12?0.0:pos.z,
                Framework[CurrentSystem].Atoms[CurrentFramework][i].Charge);
      }

      fprintf(FilePtr,"\n\n");
      fclose(FilePtr);

      if(Framework[CurrentSystem].SpaceGroupIdentifier[CurrentFramework]>1)
      {
        sprintf(buffer,"Movies/System_%d/Framework_%d_%s_%d_%d_%d%s.cif",
                CurrentSystem,CurrentFramework,string,1,1,1,FileNameAppend);

        FilePtr=fopen(buffer,"w");

        fprintf(FilePtr,"data_%s\n\n",Framework[CurrentSystem].Name[CurrentFramework]);

        fprintf(FilePtr,"_audit_creation_method RASPA-1.0\n");
        fprintf(FilePtr,"_audit_creation_date %d-%d-%d\n",loctime->tm_year+1900,loctime->tm_mon+1,loctime->tm_mday);
        #if defined(__INTEL_COMPILER)
          getlogin_r(fullname,256);
          fprintf(FilePtr,"_audit_author_name '%s'\n",fullname);
        #else
          p=getpwuid(geteuid());
          if(p)
          {
            n=strcspn(p->pw_gecos,",");
            memcpy(fullname,p->pw_gecos,n);
            fullname[n] = '\0';
            fprintf(FilePtr,"_audit_author_name '%s'\n",fullname);
          }
        #endif
        fprintf(FilePtr,"\n");

        if(Framework[CurrentSystem].NumberOfCitations[CurrentFramework]>0)
        {
          if(Framework[CurrentSystem].NumberOfCitations[CurrentFramework]==1)
          {
            if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationId)>0)
              fprintf(FilePtr,"_citation_id                 %s\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationId);
            if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationAuthorName)>0)
              fprintf(FilePtr,"_citation_author_name        '%s'\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationAuthorName);
            if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationCoordinateLinkage)>0)
              fprintf(FilePtr,"_citation_coordinate_linkage %s\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationCoordinateLinkage);
            if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationTitle)>0)
              fprintf(FilePtr,"_citation_title              '%s'\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationTitle);
            if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationCountry)>0)
              fprintf(FilePtr,"_citation_country            '%s'\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationCountry);
            if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationJournalAbbrev)>0)
              fprintf(FilePtr,"_citation_journal_abbrev     '%s'\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationJournalAbbrev);
            if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationJournalVolume)>0)
              fprintf(FilePtr,"_citation_journal_volume     %s\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationJournalVolume);
            if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationJournalIssue)>0)
              fprintf(FilePtr,"_citation_journal_issue      %s\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationJournalIssue);
            if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationPageFirst)>0)
              fprintf(FilePtr,"_citation_page_first         %s\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationPageFirst);
            if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationPageLast)>0)
              fprintf(FilePtr,"_citation_page_last          %s\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationPageLast);
            if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationYear)>0)
              fprintf(FilePtr,"_citation_year               %s\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationYear);
            if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationJournalID_ASTM)>0)
              fprintf(FilePtr,"_citation_journal_ID_ASTM    %s\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationJournalID_ASTM);
            if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationJournalID_ISSN)>0)
              fprintf(FilePtr,"_citation_journal_ID_ISSN    %s\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationJournalID_ISSN);
            if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationBookTitle)>0)
              fprintf(FilePtr,"_citation_book_title         '%s'\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationBookTitle);
            if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationBookPublisher)>0)
              fprintf(FilePtr,"_citation_book_publisher     '%s'\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationBookPublisher);
            if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationBookID_ISBN)>0)
              fprintf(FilePtr,"_citation_book_ID_ISBN       %s\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationBookID_ISBN);
            if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationSpecialDetails)>0)
              fprintf(FilePtr,"_citation_special_details    '%s'\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationSpecialDetails);
            fprintf(FilePtr,"\n");
          }
          else
          {
            fprintf(FilePtr,"loop_\n");

            if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationId)>0) fprintf(FilePtr,"_citation_id\n");
            if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationAuthorName)>0) fprintf(FilePtr,"_citation_author_name\n");
            if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationCoordinateLinkage)>0) fprintf(FilePtr,"_citation_coordinate_linkage\n");
            if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationTitle)>0) fprintf(FilePtr,"_citation_title\n");
            if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationCountry)>0) fprintf(FilePtr,"_citation_country\n");
            if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationJournalAbbrev)>0) fprintf(FilePtr,"_citation_journal_abbrev\n");
            if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationJournalVolume)>0) fprintf(FilePtr,"_citation_journal_volume\n");
            if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationJournalIssue)>0) fprintf(FilePtr,"_citation_journal_issue\n");
            if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationPageFirst)>0) fprintf(FilePtr,"_citation_page_first\n");
            if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationPageLast)>0) fprintf(FilePtr,"_citation_page_last\n");
            if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationYear)>0) fprintf(FilePtr,"_citation_year\n");
            if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationJournalID_ASTM)>0) fprintf(FilePtr,"_citation_journal_ID_ASTM\n");
            if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationJournalID_ISSN)>0) fprintf(FilePtr,"_citation_journal_ID_ISSN\n");
            if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationBookTitle)>0) fprintf(FilePtr,"_citation_book_title\n");
            if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationBookPublisher)>0) fprintf(FilePtr,"_citation_book_publisher\n");
            if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationBookID_ISBN)>0) fprintf(FilePtr,"_citation_book_ID_ISBN\n");
            if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][0].CitationSpecialDetails)>0) fprintf(FilePtr,"_citation_special_details\n");

            for(i=0;i<Framework[CurrentSystem].NumberOfCitations[CurrentFramework];i++)
            {
              if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationId)>0)
                fprintf(FilePtr,"%s\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationId);
              if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationAuthorName)>0)
                fprintf(FilePtr,"'%s'\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationAuthorName);
              if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationCoordinateLinkage)>0)
                fprintf(FilePtr,"%s\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationCoordinateLinkage);
              if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationTitle)>0)
                fprintf(FilePtr,"'%s'\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationTitle);
              if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationCountry)>0)
                fprintf(FilePtr,"'%s'\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationCountry);
              if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationJournalAbbrev)>0)
                fprintf(FilePtr,"'%s'\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationJournalAbbrev);
              if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationJournalVolume)>0)
                fprintf(FilePtr,"%s\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationJournalVolume);
              if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationJournalIssue)>0)
                fprintf(FilePtr,"%s\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationJournalIssue);
              if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationPageFirst)>0)
                fprintf(FilePtr,"%s\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationPageFirst);
              if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationPageLast)>0)
                fprintf(FilePtr,"%s\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationPageLast);
              if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationYear)>0)
                fprintf(FilePtr,"%s\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationYear);
              if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationJournalID_ASTM)>0)
                fprintf(FilePtr,"%s\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationJournalID_ASTM);
              if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationJournalID_ISSN)>0)
                fprintf(FilePtr,"%s\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationJournalID_ISSN);
              if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationBookTitle)>0)
                fprintf(FilePtr,"'%s'\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationBookTitle);
              if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationBookPublisher)>0)
                fprintf(FilePtr,"'%s'\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationBookPublisher);
              if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationBookID_ISBN)>0)
                fprintf(FilePtr,"%s\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationBookID_ISBN);
              if(strlen(Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationSpecialDetails)>0)
                fprintf(FilePtr,"'%s'\n",Framework[CurrentSystem].CitationInformation[CurrentFramework][i].CitationSpecialDetails);
            }
            fprintf(FilePtr,"\n");
          }
        }


        fprintf(FilePtr,"_cell_length_a    %g\n",UnitCellSize[CurrentSystem].x);
        fprintf(FilePtr,"_cell_length_b    %g\n",UnitCellSize[CurrentSystem].y);
        fprintf(FilePtr,"_cell_length_c    %g\n",UnitCellSize[CurrentSystem].z);
        fprintf(FilePtr,"_cell_angle_alpha %g\n",(REAL)AlphaAngle[CurrentSystem]*RAD2DEG);
        fprintf(FilePtr,"_cell_angle_beta  %g\n",(REAL)BetaAngle[CurrentSystem]*RAD2DEG);
        fprintf(FilePtr,"_cell_angle_gamma %g\n",(REAL)GammaAngle[CurrentSystem]*RAD2DEG);
        fprintf(FilePtr,"_cell_volume      %g\n\n",(REAL)Volume[CurrentSystem]);

        switch(SpaceGroupData[Framework[CurrentSystem].SpaceGroupIdentifier[CurrentFramework]].CrystalSystem)
        {
          case TRICLINIC_SPACEGROUP:
            fprintf(FilePtr,"_symmetry_cell_setting          triclinic\n");
            break;
          case MONOCLINIC_SPACEGROUP:
            fprintf(FilePtr,"_symmetry_cell_setting          monoclinic\n");
            break;
          case ORTHORHOMBIC_SPACEGROUP:
            fprintf(FilePtr,"_symmetry_cell_setting          orthorhombic\n");
            break;
          case TETRAGONAL_SPACEGROUP:
            fprintf(FilePtr,"_symmetry_cell_setting          tetragonal\n");
            break;
          case TRIGONAL_SPACEGROUP:
            fprintf(FilePtr,"_symmetry_cell_setting          trigonal\n");
            break;
          case HEXAGONAL_SPACEGROUP:
            fprintf(FilePtr,"_symmetry_cell_setting          hexagonal\n");
            break;
          case CUBIC_SPACEGROUP:
            fprintf(FilePtr,"_symmetry_cell_setting          cubic\n");
            break;
        }
        fprintf(FilePtr,"_symmetry_space_group_name_Hall '%s'\n",SpaceGroupData[Framework[CurrentSystem].SpaceGroupIdentifier[CurrentFramework]].HallSpaceGroupSymbol);
        // remove ":1", ":2", ":h" or ":r" from the short HM space group symbol
        sscanf(SpaceGroupData[Framework[CurrentSystem].SpaceGroupIdentifier[CurrentFramework]].ShortInternationalHermannMauguinSpaceGroupSymbol,"%[^:]%*[^\n]",buffer);
        fprintf(FilePtr,"_symmetry_space_group_name_H-M  '%s'\n",buffer);
        fprintf(FilePtr,"_symmetry_Int_Tables_number     %d\n\n",SpaceGroupData[Framework[CurrentSystem].SpaceGroupIdentifier[CurrentFramework]].Number);

        // printf _symmetry_equiv_pos_as_xyz for all space groups that have different origins
        // for these space groups, the H-M symbol is not sufficient to uniquely set the spacegroup
        fprintf(FilePtr,"loop_\n");
        fprintf(FilePtr,"_symmetry_equiv_pos_as_xyz\n");
        for(j=0;j<SpaceGroupData[Framework[CurrentSystem].SpaceGroupIdentifier[CurrentFramework]].NumberOfOperators;j++)
          fprintf(FilePtr," '%s'\n",StringSpaceGroupElements[Framework[CurrentSystem].SpaceGroupIdentifier[CurrentFramework]][j]);
        fprintf(FilePtr,"\n");

        fprintf(FilePtr,"loop_\n");
        fprintf(FilePtr,"_atom_site_label\n");
        fprintf(FilePtr,"_atom_site_type_symbol\n");
        fprintf(FilePtr,"_atom_site_fract_x\n");
        fprintf(FilePtr,"_atom_site_fract_y\n");
        fprintf(FilePtr,"_atom_site_fract_z\n");
        fprintf(FilePtr,"_atom_site_charge\n");

        for(i=0;i<NumberOfPseudoAtoms;i++)
          AtomIdentifier[i]=0;

        for(i=0;i<Framework[CurrentSystem].NumberOfAsymmetricAtoms[CurrentFramework];i++)
        {
          Type=Framework[CurrentSystem].AtomsAsymmetric[CurrentFramework][i].Type;

          if(Framework[CurrentSystem].AddAtomNumberCodeToLabel[CurrentFramework])
          {
            AtomIdentifier[Type]++;
            sprintf(name,"%s%d",PseudoAtoms[Type].Name,AtomIdentifier[Type]);
          }
          else
            sprintf(name,"%s",PseudoAtoms[Type].Name);

          // for asymmetric atoms the position is already in fractional coordinates
          pos=Framework[CurrentSystem].AtomsAsymmetric[CurrentFramework][i].Position;

          //fprintf(FilePtr,"%-8s %-5s %-4s % -12g % -12g % -12g % -7g % -7g\n",
          fprintf(FilePtr,"%-8s %-5s % -18.12f % -18.12f % -18.12f % -7g\n",
                  name,
                  //PseudoAtoms[Type].PrintToPDBName,
                  PseudoAtoms[Type].ChemicalElement,
                  //PseudoAtoms[Type].Hybridization==SP3?"sp3":(PseudoAtoms[Type].Hybridization==SP2?"sp2":"sp"),
                  fabs(pos.x)<1e-12?0.0:pos.x,
                  fabs(pos.y)<1e-12?0.0:pos.y,
                  fabs(pos.z)<1e-12?0.0:pos.z,
                  Framework[CurrentSystem].Atoms[CurrentFramework][i].Charge);
                  //PseudoAtoms[Type].Charge1);
                  //ComputePolarization?PseudoAtoms[Type].Polarization*COULOMBIC_CONVERSION_FACTOR:0);
        }

        fprintf(FilePtr,"\n\n");
        fclose(FilePtr);
      }
    }
  }
  free(AtomIdentifier);
}

/*********************************************************************************************************
 * Name       | ReadFrameworkDefinitionMOL                                                               *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | reads a MUSIC CMOL-file                                                                  *
 * Parameters | -                                                                                        *
 *********************************************************************************************************/

void ReadFrameworkDefinitionMOL(void)
{
  int i,j,k,l,m,temp,index;
  FILE* FilePtr;
  char buffer[256],AtomType[10];
  REAL tempd,det;
  REAL tempx,tempy,tempz;
  REAL charge;
  VECTOR shift;
  double A,B,C;
  int Type=0;

  Framework[CurrentSystem].FrameworkDensityPerComponent[CurrentFramework]=0.0;
  Framework[CurrentSystem].FrameworkMassPerComponent[CurrentFramework]=0.0;

  if (STREAM)
  {
#ifdef __unix__
    if (!(FilePtr=fmemopen((void *)INPUT_CRYSTAL, strlen(INPUT_CRYSTAL), "r")))
    {
      printf("Error reading streamed MOL molecule.");
      exit(1);
    }
#else
    fprintf(stderr, "Streaming only allowed on POSIX systems (for now)\n.");
    exit(1);
#endif
  }
  else
  {
    // first try to open the framework-file in the current directory,
    // and next from the repository
    sprintf(buffer,"%s.%s",
            Framework[CurrentSystem].Name[CurrentFramework],
            "mol");

    if(!(FilePtr=fopen(buffer,"r")))
    {
      sprintf(buffer,"%s/share/raspa/structures/xyz/%s.%s",
              RASPA_DIRECTORY,
              Framework[CurrentSystem].Name[CurrentFramework],
              "mol");

      if(!(FilePtr=fopen(buffer,"r")))
      {
        fprintf(stderr, "Error:  file %s does not exist.\n",buffer);
        exit(1);
      }
    }
  }

  // skip the first 3 lines
  fscanf(FilePtr,"%[^\n]\n",buffer);
  fscanf(FilePtr,"%[^\n]\n",buffer);

  // get number of atoms of one unitcell
  fscanf(FilePtr,"%d%*[^\n]]",&temp);
  Framework[CurrentSystem].NumberOfUnitCellAtoms[CurrentFramework]=temp;

  temp=temp*NumberOfUnitCells[CurrentSystem].x*NumberOfUnitCells[CurrentSystem].y*NumberOfUnitCells[CurrentSystem].z;
  Framework[CurrentSystem].NumberOfAtoms[CurrentFramework]=temp;
  Framework[CurrentSystem].MaxNumberOfAtoms[CurrentFramework]=temp;

  Framework[CurrentSystem].Atoms[CurrentFramework]=(ATOM*)calloc(temp,sizeof(ATOM));

  for(i=0;i<Framework[CurrentSystem].NumberOfUnitCellAtoms[CurrentFramework];i++)
  {
    fscanf(FilePtr,"%d %lf %lf %lf %s %lf %lf %lf\n",&index,&tempx,&tempy,&tempz,(char*)AtomType,&charge,&B,&C);

    // set type of the atom
    Type=-1;
    for(m=0;m<NumberOfPseudoAtoms;m++)
    {
      if(strcasecmp(PseudoAtoms[m].Name,AtomType)==0)
        Type=m;
    }
    if(Type<0)
    {
      fprintf(stderr, "Unknown PseudoAtom-type %s !!\n",AtomType);
      exit(0);
    }

    if(UseChargesFromMOLFile)
    {
      Framework[CurrentSystem].Atoms[CurrentFramework][i].Charge=charge;
      PseudoAtoms[Type].Charge1=charge;
    }
    else
      Framework[CurrentSystem].Atoms[CurrentFramework][i].Charge=PseudoAtoms[Type].Charge1;

    Framework[CurrentSystem].Atoms[CurrentFramework][i].Type=Type;
    PseudoAtoms[Type].FrameworkAtom=TRUE;
    NumberOfPseudoAtomsCount[CurrentSystem][Type]++;
    NumberOfPseudoAtomsType[CurrentSystem][Type]++;
    Framework[CurrentSystem].FrameworkMassPerComponent[CurrentFramework]+=PseudoAtoms[Type].Mass;
    Framework[CurrentSystem].Atoms[CurrentFramework][i].Position.x=tempx;
    Framework[CurrentSystem].Atoms[CurrentFramework][i].Position.y=tempy;
    Framework[CurrentSystem].Atoms[CurrentFramework][i].Position.z=tempz;
    if(fabs(charge)>1e-10)
    {
      Framework[CurrentSystem].NumberOfCharges[CurrentFramework]++;
      Framework[CurrentSystem].TotalNumberOfCharges++;
    }
  }

  // skip the 4 lines
  fscanf(FilePtr,"%[^\n]",buffer);
  fscanf(FilePtr,"%[^\n]",buffer);
  fscanf(FilePtr,"%[^\n]",buffer);
  fscanf(FilePtr,"%[^\n]",buffer);
  // read unit cell size
  fscanf(FilePtr,"%lf %lf %lf%*[^\n]\n",&tempx,&tempy,&tempz);
  UnitCellSize[CurrentSystem].x=(REAL)tempx;
  UnitCellSize[CurrentSystem].y=(REAL)tempy;
  UnitCellSize[CurrentSystem].z=(REAL)tempz;
  A=(REAL)NumberOfUnitCells[CurrentSystem].x*UnitCellSize[CurrentSystem].x;
  B=(REAL)NumberOfUnitCells[CurrentSystem].y*UnitCellSize[CurrentSystem].y;
  C=(REAL)NumberOfUnitCells[CurrentSystem].z*UnitCellSize[CurrentSystem].z;

  // read cell angles and convert to radians
  fscanf(FilePtr,"%lf %lf %lf%*[^\n]",&tempx,&tempy,&tempz);
  AlphaAngle[CurrentSystem]=tempx;
  BetaAngle[CurrentSystem]=tempy;
  GammaAngle[CurrentSystem]=tempz;

  if(BoundaryCondition[CurrentSystem]==UNINITIALIZED_BOUNDARY_CONDITION)
  {
    // determine boundary conditions from angles
    if((fabs(AlphaAngle[CurrentSystem]-90.0)>0.01)||
       (fabs(BetaAngle[CurrentSystem]-90.0)>0.01)||
       (fabs(GammaAngle[CurrentSystem]-90.0)>0.01))
      BoundaryCondition[CurrentSystem]=TRICLINIC;
    else
    {
      if((fabs(A-B)>0.01)||(fabs(A-C)>0.01)||(fabs(B-C)>0.01))
        BoundaryCondition[CurrentSystem]=RECTANGULAR;
      else BoundaryCondition[CurrentSystem]=RECTANGULAR;
    }
  }
  AlphaAngle[CurrentSystem]*=M_PI/180.0;
  BetaAngle[CurrentSystem]*=M_PI/180.0;
  GammaAngle[CurrentSystem]*=M_PI/180.0;

  // construct transformation matrix Box to go from abc-coordinates to xyz-coordinates
  A=UnitCellSize[CurrentSystem].x;
  B=UnitCellSize[CurrentSystem].y;
  C=UnitCellSize[CurrentSystem].z;
  tempd=(cos(AlphaAngle[CurrentSystem])-cos(GammaAngle[CurrentSystem])*cos(BetaAngle[CurrentSystem]))/sin(GammaAngle[CurrentSystem]);

  // first vector
  UnitCellBox[CurrentSystem].ax=A;
  UnitCellBox[CurrentSystem].ay=0.0;
  UnitCellBox[CurrentSystem].az=0.0;

  // second vector
  UnitCellBox[CurrentSystem].bx=B*cos(GammaAngle[CurrentSystem]);
  UnitCellBox[CurrentSystem].by=B*sin(GammaAngle[CurrentSystem]);
  UnitCellBox[CurrentSystem].bz=0.0;

  // third vector
  UnitCellBox[CurrentSystem].cx=C*cos(BetaAngle[CurrentSystem]);
  UnitCellBox[CurrentSystem].cy=C*tempd;
  UnitCellBox[CurrentSystem].cz=C*sqrt(1.0-SQR(cos(BetaAngle[CurrentSystem]))-SQR(tempd));

  Invert3x3Matrix(&UnitCellBox[CurrentSystem],&InverseUnitCellBox[CurrentSystem],&det);

  A=(REAL)NumberOfUnitCells[CurrentSystem].x*UnitCellSize[CurrentSystem].x;
  B=(REAL)NumberOfUnitCells[CurrentSystem].y*UnitCellSize[CurrentSystem].y;
  C=(REAL)NumberOfUnitCells[CurrentSystem].z*UnitCellSize[CurrentSystem].z;
  tempd=(cos(AlphaAngle[CurrentSystem])-cos(GammaAngle[CurrentSystem])*cos(BetaAngle[CurrentSystem]))/sin(GammaAngle[CurrentSystem]);

  // first vector
  Box[CurrentSystem].ax=A;
  Box[CurrentSystem].ay=0.0;
  Box[CurrentSystem].az=0.0;

  // second vector
  Box[CurrentSystem].bx=B*cos(GammaAngle[CurrentSystem]);
  Box[CurrentSystem].by=B*sin(GammaAngle[CurrentSystem]);
  Box[CurrentSystem].bz=0.0;

  // third vector
  Box[CurrentSystem].cx=C*cos(BetaAngle[CurrentSystem]);
  Box[CurrentSystem].cy=C*tempd;
  Box[CurrentSystem].cz=C*sqrt(1.0-SQR(cos(BetaAngle[CurrentSystem]))-SQR(tempd));

  // calculate box-properties
  CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);

  // calculate inverse box
  Invert3x3Matrix(&Box[CurrentSystem],&InverseBox[CurrentSystem],&det);

  // calculate box-properties
  CellProperties(&InverseBox[CurrentSystem],&InverseBoxProperties[CurrentSystem],&det);

  index=Framework[CurrentSystem].NumberOfUnitCellAtoms[CurrentFramework];
  for(i=0;i<Framework[CurrentSystem].NumberOfUnitCellAtoms[CurrentFramework];i++)
  {
    //fscanf(FilePtr,"%d %s %lf %lf %lf %*[^\n]",&temp,(char*)&AtomType,&tempx,&tempy,&tempz);
    tempx=Framework[CurrentSystem].Atoms[CurrentFramework][i].Position.x;
    tempy=Framework[CurrentSystem].Atoms[CurrentFramework][i].Position.y;
    tempz=Framework[CurrentSystem].Atoms[CurrentFramework][i].Position.z;
    Type=Framework[CurrentSystem].Atoms[CurrentFramework][i].Type;
    charge=Framework[CurrentSystem].Atoms[CurrentFramework][i].Charge;

    for(j=0;j<NumberOfUnitCells[CurrentSystem].x;j++)
      for(k=0;k<NumberOfUnitCells[CurrentSystem].y;k++)
        for(l=0;l<NumberOfUnitCells[CurrentSystem].z;l++)
        {
          if(!((j==0)&&(k==0)&&(l==0)))
          {
            Framework[CurrentSystem].Atoms[CurrentFramework][index].Type=Type;
            NumberOfPseudoAtomsCount[CurrentSystem][Type]++;
            NumberOfPseudoAtomsType[CurrentSystem][Type]++;
            Framework[CurrentSystem].FrameworkMassPerComponent[CurrentFramework]+=PseudoAtoms[Type].Mass;
            shift.x=(REAL)j/NumberOfUnitCells[CurrentSystem].x;
            shift.y=(REAL)k/NumberOfUnitCells[CurrentSystem].y;
            shift.z=(REAL)l/NumberOfUnitCells[CurrentSystem].z;
            shift=ConvertFromABCtoXYZ(shift);
            Framework[CurrentSystem].Atoms[CurrentFramework][index].Position.x=tempx+shift.x;
            Framework[CurrentSystem].Atoms[CurrentFramework][index].Position.y=tempy+shift.y;
            Framework[CurrentSystem].Atoms[CurrentFramework][index].Position.z=tempz+shift.z;
            Framework[CurrentSystem].Atoms[CurrentFramework][index].Charge=charge;
            if(fabs(charge)>1e-10)
            {
              Framework[CurrentSystem].NumberOfCharges[CurrentFramework]++;
              Framework[CurrentSystem].TotalNumberOfCharges++;
            }
            index++;
          }
        }
  }
  fclose(FilePtr);

  Framework[CurrentSystem].FrameworkDensityPerComponent[CurrentFramework]=
  Framework[CurrentSystem].FrameworkMassPerComponent[CurrentFramework]/
  (Volume[CurrentSystem]*CUBE(ANGSTROM)*AVOGADRO_CONSTANT);
  Framework[CurrentSystem].NumberOfAtoms[CurrentFramework]=index;

  Framework[CurrentSystem].TotalNumberOfAtoms+=Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];
  Framework[CurrentSystem].FrameworkDensity+=Framework[CurrentSystem].FrameworkDensityPerComponent[CurrentFramework];
  Framework[CurrentSystem].FrameworkMass+=Framework[CurrentSystem].FrameworkMassPerComponent[CurrentFramework];
}

void WriteFrameworkDefinitionMOL(char *string)
{
  int j,Type;
  char buffer[256];
  FILE *FilePtr;
  int *AtomId;
  VECTOR r;
  REAL charge;
  char Name[256];

  if (STREAM)
  {
    fprintf(stderr, "File writing not allowed in streaming mode!");
    return;
  }

  AtomId=(int*)calloc(NumberOfPseudoAtoms,sizeof(int));
  mkdir("Movies",S_IRWXU);

  for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
  {
    if(Framework[CurrentSystem].FrameworkModel!=NONE)
    {
      for(CurrentFramework=0;CurrentFramework<Framework[CurrentSystem].NumberOfFrameworks;CurrentFramework++)
      {
        sprintf(buffer,"Movies/System_%d",CurrentSystem);
        mkdir(buffer,S_IRWXU);

        // write mol-format
        sprintf(buffer,"Movies/System_%d/Framework_%d_%s%s.mol",CurrentSystem,CurrentFramework,string,FileNameAppend);
        FilePtr=fopen(buffer,"w");
        fprintf(FilePtr," Molecule_name: %s\n",Framework[CurrentSystem].Name[CurrentFramework]);
        fprintf(FilePtr,"\n");
        fprintf(FilePtr,"  Coord_Info: Listed Cartesian None\n");
        fprintf(FilePtr,"%12d\n",Framework[CurrentSystem].NumberOfAtoms[CurrentFramework]);
        for(j=0;j<Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];j++)
        {
          Type=Framework[CurrentSystem].Atoms[CurrentFramework][j].Type;
          AtomId[Type]++;

          sprintf(Name,"%s",PseudoAtoms[Type].PrintToPDBName);

          r=Framework[CurrentSystem].Atoms[CurrentFramework][j].Position;
          charge=Framework[CurrentSystem].Atoms[CurrentFramework][j].Charge;

          fprintf(FilePtr,"%5d %9.5lf %9.5lf %9.5lf   %-10s % lf  0  0\n",
                  j+1,(double)r.x,(double)r.y,(double)r.z,Name,(double)charge);
        }
        fprintf(FilePtr,"\n\n\n");
        fprintf(FilePtr,"  Fundcell_Info: Listed\n");
        fprintf(FilePtr," %14.5lf %14.5lf%14.5lf\n",(double)BoxProperties[CurrentSystem].ax,(double)BoxProperties[CurrentSystem].ay,(double)BoxProperties[CurrentSystem].az);
        fprintf(FilePtr," %14.5lf %14.5lf%14.5lf\n",(double)AlphaAngle[CurrentSystem]*RAD2DEG,
                (double)BetaAngle[CurrentSystem]*RAD2DEG,(double)GammaAngle[CurrentSystem]*RAD2DEG);
        fprintf(FilePtr," %14.5lf %14.5lf%14.5lf\n",0.0,0.0,0.0);
        fprintf(FilePtr," %14.5lf %14.5lf%14.5lf\n",(double)BoxProperties[CurrentSystem].ax,(double)BoxProperties[CurrentSystem].ay,                    (double)BoxProperties[CurrentSystem].az);
        fclose(FilePtr);

        // write asymetric unit-cell in pdb-format
        /*
         SerialNumber=1;
         sprintf(buffer,"Movies/System_%d/AsymetricUnitCell.pdb",CurrentSystem);
         FilePtr=fopen(buffer,"w");
         fprintf(FilePtr,"REMARK   Raspa-1.0 PDB file\n");
         for(j=0;j<Framework[CurrentSystem].NumberOfAsymmetricAtoms[CurrentFramework];j++)
         {
         Type=Framework[CurrentSystem].AtomsAsymmetric[CurrentFramework][j].Type;
         if(PseudoAtoms[Type].PrintToPDB)
         {
         pos=Framework[CurrentSystem].AtomsAsymmetric[CurrentFramework][j].Position;

         r.x=Box[CurrentSystem].ax*pos.x+Box[CurrentSystem].ay*pos.y+Box[CurrentSystem].az*pos.z;
         r.y=Box[CurrentSystem].bx*pos.x+Box[CurrentSystem].by*pos.y+Box[CurrentSystem].bz*pos.z;
         r.z=Box[CurrentSystem].cx*pos.x+Box[CurrentSystem].cy*pos.y+Box[CurrentSystem].cz*pos.z;

         sprintf(AtomName,"%2s",PseudoAtoms[Type].PrintToPDBName);
         sprintf(Element,"%2s",PseudoAtoms[Type].PrintToPDBName);
         fprintf(FilePtr,"%s%5d %2s%c%c%c%3s %c%4s%c   %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf%6s%-4s%2s%2s\n",
         RecordName,SerialNumber++,AtomName,RemotenessIndicator,BranchDesignator,AltLoc,ResIdueName,
         ChainId,ResSeq,iCode,(double)r.x,(double)r.y,(double)r.z,(double)Occupancy,(double)Temp,"      ",SegID,Element,charge);
         }
         }
         fclose(FilePtr);
         */
      }
    }
  }
  free(AtomId);
}


/*********************************************************************************************************
 * Name       | ReadFrameworkDefinitionDLPOLY                                                            *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | reads a dlpoly CONFIG-file (no symmetry)                                                 *
 * Parameters | -                                                                                        *
 *********************************************************************************************************/

void ReadFrameworkDefinitionDLPOLY(void)
{
  int i;
  FILE* FilePtr;
  char buffer[256],AtomType[10];
  REAL det;
  int Type=0;
  int int_temp1,int_temp2,int_temp3,int_temp4;
  REAL real_temp1,real_temp2,real_temp3;
  REAL_MATRIX3x3 UnitCellBoxProperties;
  VECTOR pos;

  fprintf(stderr, "CurrentSystem: %d CurrentFramework: %d\n",CurrentSystem,CurrentFramework);
  Framework[CurrentSystem].SpaceGroupIdentifier[CurrentFramework]=1;

  if (STREAM)
  {
#ifdef __unix__
    if (!(FilePtr=fmemopen((void *)INPUT_CRYSTAL, strlen(INPUT_CRYSTAL), "r")))
    {
      printf("Error reading streamed CSSR molecule.");
      exit(1);
    }
#else
    fprintf(stderr, "Streaming only allowed on POSIX systems (for now)\n.");
    exit(1);
#endif
  }
  else
  {
    // first try to open the framework-file in the current directory,
    // and next from the repository
    sprintf(buffer,"./%s.%s",
          Framework[CurrentSystem].Name[CurrentFramework],
          "dlpoly");
    if(!(FilePtr=fopen(buffer,"r")))
    {
      fprintf(stderr, "Error:  file %s does not exist.\n",buffer);
      exit(1);
    }
  }

  fscanf(FilePtr,"%*[^\n]"); // skip first line
  fscanf(FilePtr,"%d %d %d %d %lf %lf\n",&int_temp1,&int_temp2,&int_temp3,&int_temp4,&real_temp1,&real_temp2);

  // NOTE: transposed for DLPOLY input (vectors in DLPOLY are given at line 1,2,3.
  fscanf(FilePtr,"%lf %lf %lf\n",&UnitCellBox[CurrentSystem].ax,&UnitCellBox[CurrentSystem].ay,&UnitCellBox[CurrentSystem].az);
  fscanf(FilePtr,"%lf %lf %lf\n",&UnitCellBox[CurrentSystem].bx,&UnitCellBox[CurrentSystem].by,&UnitCellBox[CurrentSystem].bz);
  fscanf(FilePtr,"%lf %lf %lf\n",&UnitCellBox[CurrentSystem].cx,&UnitCellBox[CurrentSystem].cy,&UnitCellBox[CurrentSystem].cz);

  Invert3x3Matrix(&UnitCellBox[CurrentSystem],&InverseUnitCellBox[CurrentSystem],&det);

  // calculate box-properties
  CellProperties(&UnitCellBox[CurrentSystem],&UnitCellBoxProperties,&det);

  UnitCellSize[CurrentSystem].x=UnitCellBoxProperties.ax;
  UnitCellSize[CurrentSystem].y=UnitCellBoxProperties.ay;
  UnitCellSize[CurrentSystem].z=UnitCellBoxProperties.az;

  Box[CurrentSystem].ax=NumberOfUnitCells[CurrentSystem].x*UnitCellBox[CurrentSystem].ax;
  Box[CurrentSystem].ay=NumberOfUnitCells[CurrentSystem].x*UnitCellBox[CurrentSystem].ay;
  Box[CurrentSystem].az=NumberOfUnitCells[CurrentSystem].x*UnitCellBox[CurrentSystem].az;

  Box[CurrentSystem].bx=NumberOfUnitCells[CurrentSystem].y*UnitCellBox[CurrentSystem].bx;
  Box[CurrentSystem].by=NumberOfUnitCells[CurrentSystem].y*UnitCellBox[CurrentSystem].by;
  Box[CurrentSystem].bz=NumberOfUnitCells[CurrentSystem].y*UnitCellBox[CurrentSystem].bz;

  Box[CurrentSystem].cx=NumberOfUnitCells[CurrentSystem].z*UnitCellBox[CurrentSystem].cx;
  Box[CurrentSystem].cy=NumberOfUnitCells[CurrentSystem].z*UnitCellBox[CurrentSystem].cy;
  Box[CurrentSystem].cz=NumberOfUnitCells[CurrentSystem].z*UnitCellBox[CurrentSystem].cz;

  Invert3x3Matrix(&Box[CurrentSystem],&InverseBox[CurrentSystem],&det);

  // calculate box-properties
  CellProperties(&InverseBox[CurrentSystem],&InverseBoxProperties[CurrentSystem],&det);
  CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);

  AlphaAngle[CurrentSystem]=acos(BoxProperties[CurrentSystem].bx);
  BetaAngle[CurrentSystem]=acos(BoxProperties[CurrentSystem].by);
  GammaAngle[CurrentSystem]=acos(BoxProperties[CurrentSystem].bz);

  // determine boundary conditions from angles
  if((fabs(AlphaAngle[CurrentSystem]-90.0*DEG2RAD)>0.001)||
     (fabs(BetaAngle[CurrentSystem]-90.0*DEG2RAD)>0.001)||
     (fabs(GammaAngle[CurrentSystem]-90.0*DEG2RAD)>0.001))
    BoundaryCondition[CurrentSystem]=TRICLINIC;
  else BoundaryCondition[CurrentSystem]=RECTANGULAR;

  BoundaryCondition[CurrentSystem]=TRICLINIC;

  Framework[CurrentSystem].NumberOfAsymmetricAtoms[CurrentFramework]=int_temp3;
  Framework[CurrentSystem].AtomsAsymmetric[CurrentFramework]=(FRAMEWORK_ASYMMETRIC_ATOM*)calloc(int_temp3,sizeof(FRAMEWORK_ASYMMETRIC_ATOM));

  fprintf(stderr, "int_temp3: %d\n",int_temp3);

  for(i=0;i<int_temp3;i++)
  {
    fscanf(FilePtr,"%s%*[^\n]\n",AtomType);

    if(Framework[CurrentSystem].RemoveAtomNumberCodeFromLabel[CurrentFramework])
      sscanf(AtomType,"%[a-zA-Z]",AtomType);

    // set type of the atom
    Type=ReturnPseudoAtomNumber(AtomType);
    if(Type<0)
    {
      fprintf(stderr, "Unknown PseudoAtom-type %s !! (in'ReadFrameworkDefinitionDLPOLY')\n",AtomType);
      exit(0);
    }

    Framework[CurrentSystem].AtomsAsymmetric[CurrentFramework][i].Type=Type;
    Framework[CurrentSystem].AtomsAsymmetric[CurrentFramework][i].Charge=PseudoAtoms[Type].Charge1;
    PseudoAtoms[Type].FrameworkAtom=TRUE;
    fscanf(FilePtr,"%lf %lf %lf\n",&pos.x,&pos.y,&pos.z);
    Framework[CurrentSystem].AtomsAsymmetric[CurrentFramework][i].Position=ConvertFromXYZtoABC(pos);

    fscanf(FilePtr,"%lf %lf %lf\n",&real_temp1,&real_temp2,&real_temp3);
    fscanf(FilePtr,"%lf %lf %lf\n",&real_temp1,&real_temp2,&real_temp3);
  }
  fclose(FilePtr);

  // create the full framework from the asymmetric definition
  ExpandAsymmetricFrameworkToFullFramework();
}



/*********************************************************************************************************
 * Name       | ReadFrameworkDefinitionCSSR                                                              *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | reads a CSSR-file and expands the structure to 'P1' (no symmetry)                        *
 * Parameters | -                                                                                        *
 *********************************************************************************************************/

void ReadFrameworkDefinitionCSSR(void)
{
  int i,temp,index;
  FILE* FilePtr;
  char buffer[256],AtomType[10];
  REAL tempd,det;
  REAL A,B,C;
  int Type=0;
  VECTOR pos;
  char SpaceGroupName[32];
  int SpaceGroup;
  int Option;

  // first try to open the framework-file in the current directory,
  // and next from the repository
  sprintf(buffer,"./%s.%s",
        Framework[CurrentSystem].Name[CurrentFramework],
        "cssr");
  if(!(FilePtr=fopen(buffer,"r")))
  {
    sprintf(buffer,"%s/share/raspa/structures/cssr/%s.%s",
          RASPA_DIRECTORY,
          Framework[CurrentSystem].Name[CurrentFramework],
          "cssr");

    if(!(FilePtr=fopen(buffer,"r")))
    {
      fprintf(stderr, "Error:  file %s does not exist.\n",buffer);
      exit(1);
    }
  }

  // read unit cell size
  fscanf(FilePtr,"%lf %lf %lf\n",&UnitCellSize[CurrentSystem].x,&UnitCellSize[CurrentSystem].y,&UnitCellSize[CurrentSystem].z);
  A=(REAL)NumberOfUnitCells[CurrentSystem].x*UnitCellSize[CurrentSystem].x;
  B=(REAL)NumberOfUnitCells[CurrentSystem].y*UnitCellSize[CurrentSystem].y;
  C=(REAL)NumberOfUnitCells[CurrentSystem].z*UnitCellSize[CurrentSystem].z;

  // read cell angles and convert to radians
  fscanf(FilePtr,"%lf %lf %lf %[^\n]",&AlphaAngle[CurrentSystem],&BetaAngle[CurrentSystem],&GammaAngle[CurrentSystem],buffer);
  Option=0;
  strcpy(SpaceGroupName,"");
  if(sscanf(buffer,"SPGR =%3d %11[a-zA-Z0-9-/ ] OPT =%2d",&SpaceGroup,SpaceGroupName,&Option)<1)
  {
    fprintf(stderr, "Error reading cssr file \n");
  }
  fprintf(stderr, "%d %s %d\n",SpaceGroup,SpaceGroupName,Option);

  if(strlen(SpaceGroupName)==0)
  {
    if((SpaceGroup<1)||(SpaceGroup>230))
    {
      fprintf(stderr, "Wrong space group number: %d\n",SpaceGroup);
      exit(0);
    }
    Framework[CurrentSystem].SpaceGroupIdentifier[CurrentFramework]=GetSpaceGroupFromITSpaceGroupNumber(SpaceGroup);
  }
  else
  {
    TrimStringInPlace(SpaceGroupName);
    CompressSpacesInString(SpaceGroupName);
    Framework[CurrentSystem].SpaceGroupIdentifier[CurrentFramework]=GetSpaceGroupFromHermannMauguinString(SpaceGroupName);
    fprintf(stderr, "Framework[CurrentSystem].SpaceGroupIdentifier[CurrentFramework]: %d %s\n",Framework[CurrentSystem].SpaceGroupIdentifier[CurrentFramework],SpaceGroupName);
  }
  Framework[CurrentSystem].SpaceGroupIdentifier[CurrentFramework]=AdjustForNonStandardCSSROption(Framework[CurrentSystem].SpaceGroupIdentifier[CurrentFramework],Option);

  if(BoundaryCondition[CurrentSystem]==UNINITIALIZED_BOUNDARY_CONDITION)
  {
    // determine boundary conditions from angles
    if((fabs(AlphaAngle[CurrentSystem]-90.0)>0.001)||(fabs(BetaAngle[CurrentSystem]-90.0)>0.001)||(fabs(GammaAngle[CurrentSystem]-90.0)>0.001))
      BoundaryCondition[CurrentSystem]=TRICLINIC;
    else
    {
      if((fabs(A-B)>0.00001)||(fabs(A-C)>0.00001)||(fabs(B-C)>0.00001))
        BoundaryCondition[CurrentSystem]=RECTANGULAR;
      else BoundaryCondition[CurrentSystem]=TRICLINIC;
    }
  }
  AlphaAngle[CurrentSystem]*=M_PI/180.0;
  BetaAngle[CurrentSystem]*=M_PI/180.0;
  GammaAngle[CurrentSystem]*=M_PI/180.0;

  // construct transformation matrix Box to go from abc-coordinates to xyz-coordinates
  A=UnitCellSize[CurrentSystem].x;
  B=UnitCellSize[CurrentSystem].y;
  C=UnitCellSize[CurrentSystem].z;
  tempd=(cos(AlphaAngle[CurrentSystem])-cos(GammaAngle[CurrentSystem])*cos(BetaAngle[CurrentSystem]))/sin(GammaAngle[CurrentSystem]);

  // first vector
  UnitCellBox[CurrentSystem].ax=A;
  UnitCellBox[CurrentSystem].ay=0.0;
  UnitCellBox[CurrentSystem].az=0.0;

  // second vector
  UnitCellBox[CurrentSystem].bx=B*cos(GammaAngle[CurrentSystem]);
  UnitCellBox[CurrentSystem].by=B*sin(GammaAngle[CurrentSystem]);
  UnitCellBox[CurrentSystem].bz=0.0;

  // third vector
  UnitCellBox[CurrentSystem].cx=C*cos(BetaAngle[CurrentSystem]);
  UnitCellBox[CurrentSystem].cy=C*tempd;
  UnitCellBox[CurrentSystem].cz=C*sqrt(1.0-SQR(cos(BetaAngle[CurrentSystem]))-SQR(tempd));

  Invert3x3Matrix(&UnitCellBox[CurrentSystem],&InverseUnitCellBox[CurrentSystem],&det);

  A=(REAL)NumberOfUnitCells[CurrentSystem].x*UnitCellSize[CurrentSystem].x;
  B=(REAL)NumberOfUnitCells[CurrentSystem].y*UnitCellSize[CurrentSystem].y;
  C=(REAL)NumberOfUnitCells[CurrentSystem].z*UnitCellSize[CurrentSystem].z;
  tempd=(cos(AlphaAngle[CurrentSystem])-cos(GammaAngle[CurrentSystem])*cos(BetaAngle[CurrentSystem]))/sin(GammaAngle[CurrentSystem]);

  // first vector
  Box[CurrentSystem].ax=A;
  Box[CurrentSystem].ay=0.0;
  Box[CurrentSystem].az=0.0;

  // second vector
  Box[CurrentSystem].bx=B*cos(GammaAngle[CurrentSystem]);
  Box[CurrentSystem].by=B*sin(GammaAngle[CurrentSystem]);
  Box[CurrentSystem].bz=0.0;

  // third vector
  Box[CurrentSystem].cx=C*cos(BetaAngle[CurrentSystem]);
  Box[CurrentSystem].cy=C*tempd;
  Box[CurrentSystem].cz=C*sqrt(1.0-SQR(cos(BetaAngle[CurrentSystem]))-SQR(tempd));

  Invert3x3Matrix(&Box[CurrentSystem],&InverseBox[CurrentSystem],&det);

  // calculate box-properties
  CellProperties(&InverseBox[CurrentSystem],&InverseBoxProperties[CurrentSystem],&det);
  CellProperties(&Box[CurrentSystem],&BoxProperties[CurrentSystem],&Volume[CurrentSystem]);

  // get number of atoms of one unitcell
  fscanf(FilePtr,"%d %[^\n]",&temp,buffer);

  Framework[CurrentSystem].NumberOfAsymmetricAtoms[CurrentFramework]=temp;
  Framework[CurrentSystem].AtomsAsymmetric[CurrentFramework]=(FRAMEWORK_ASYMMETRIC_ATOM*)calloc(temp,sizeof(FRAMEWORK_ASYMMETRIC_ATOM));

  index=0;
  fscanf(FilePtr,"%d %[^\n]",&temp,buffer);
  for(i=0;i<Framework[CurrentSystem].NumberOfAsymmetricAtoms[CurrentFramework];i++)
  {
    pos.x=pos.y=pos.z=0.0;
    fscanf(FilePtr,"%d %s %lf %lf %lf %[^\n]",&temp,AtomType,&pos.x,&pos.y,&pos.z,buffer);

    if(Framework[CurrentSystem].RestrictFrameworkAtomsToBox)
    {
      pos.x-=(REAL)NINT(pos.x);
      pos.y-=(REAL)NINT(pos.y);
      pos.z-=(REAL)NINT(pos.z);
      if(pos.x<0.0) pos.x+=1.0;
      if(pos.y<0.0) pos.y+=1.0;
      if(pos.z<0.0) pos.z+=1.0;
    }

    if(Framework[CurrentSystem].RemoveAtomNumberCodeFromLabel[CurrentFramework])
      sscanf(AtomType,"%[a-zA-Z]",AtomType);

    // set type of the atom
    Type=ReturnPseudoAtomNumber(AtomType);
    if(Type<0)
    {
      fprintf(stderr, "Unknown PseudoAtom-type %s !!\n",AtomType);
      exit(0);
    }

    Framework[CurrentSystem].AtomsAsymmetric[CurrentFramework][i].Position=pos;
    Framework[CurrentSystem].AtomsAsymmetric[CurrentFramework][i].Type=Type;
    Framework[CurrentSystem].AtomsAsymmetric[CurrentFramework][i].Charge=PseudoAtoms[Type].Charge1;
    PseudoAtoms[Type].FrameworkAtom=TRUE;
  }
  fclose(FilePtr);

  // create the full framework from the asymmetric definition
  ExpandAsymmetricFrameworkToFullFramework();
}

/*********************************************************************************************************
 * Name       | WriteFrameworkDefinitionCSSR                                                             *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | write a CSSR-file using the current space group and 'P1' (no symmetry)                   *
 * Parameters | -                                                                                        *
 *********************************************************************************************************/

void WriteFrameworkDefinitionCSSR(char *string)
{
  int j;
  REAL A,B,C;
  VECTOR pos;
  char buffer[256];
  FILE *FilePtr;
  char AtomName[10]="      C";                         // Atom Name
  int  SerialNumber;                             // Atom serial number
  int *AtomId;

  if (STREAM)
  {
    fprintf(stderr, "File writing not allowed in streaming mode!");
    return;
  }

  AtomId=(int*)calloc(NumberOfPseudoAtoms,sizeof(int));
  mkdir("Movies",S_IRWXU);
  for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
  {
    if(Framework[CurrentSystem].FrameworkModel!=NONE)
    {
      for(CurrentFramework=0;CurrentFramework<Framework[CurrentSystem].NumberOfFrameworks;CurrentFramework++)
      {
        sprintf(buffer,"Movies/System_%d",CurrentSystem);
        mkdir(buffer,S_IRWXU);

        sprintf(buffer,"Movies/System_%d/Framework_%d_%s%s.cssr",CurrentSystem,CurrentFramework,string,FileNameAppend);
        FilePtr=fopen(buffer,"w");

        A=(REAL)NumberOfUnitCells[CurrentSystem].x*UnitCellSize[CurrentSystem].x;
        B=(REAL)NumberOfUnitCells[CurrentSystem].y*UnitCellSize[CurrentSystem].y;
        C=(REAL)NumberOfUnitCells[CurrentSystem].z*UnitCellSize[CurrentSystem].z;

        fprintf(FilePtr,"%38c %7.3lf %7.3lf %7.3lf\n",' ',
                (double)A,
                (double)B,
                (double)C);

        sscanf(SpaceGroupData[Framework[CurrentSystem].SpaceGroupIdentifier[CurrentFramework]].ShortInternationalHermannMauguinSpaceGroupSymbol,"%[^:]%*[^\n]",buffer);
        fprintf(FilePtr,"%21c %7.3lf %7.3lf %7.3lf    SPGR =%3d %-11s",
                ' ',
                (double)AlphaAngle[CurrentSystem]*RAD2DEG,
                (double)BetaAngle[CurrentSystem]*RAD2DEG,
                (double)GammaAngle[CurrentSystem]*RAD2DEG,
                SpaceGroupData[Framework[CurrentSystem].SpaceGroupIdentifier[CurrentFramework]].Number,buffer);
        if(SpaceGroupData[Framework[CurrentSystem].SpaceGroupIdentifier[CurrentFramework]].SpaceGroupOptionCSSR>1)
          fprintf(FilePtr," OPT = %d",SpaceGroupData[Framework[CurrentSystem].SpaceGroupIdentifier[CurrentFramework]].SpaceGroupOptionCSSR);
        fprintf(FilePtr,"\n");

        if(SpaceGroupData[Framework[CurrentSystem].SpaceGroupIdentifier[CurrentFramework]].Number==1)
        {
          fprintf(FilePtr,"%4d   %d %s\n",Framework[CurrentSystem].TotalNumberOfAtoms,0,"Created by Raspa-1.0");
          fprintf(FilePtr,"     0 %s         : %s\n",Framework[CurrentSystem].Name[0],
                  Framework[CurrentSystem].Name[0]);
          SerialNumber=1;
          for(j=0;j<Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];j++)
          {
            if(PseudoAtoms[Framework[CurrentSystem].Atoms[CurrentFramework][j].Type].PrintToPDB)
            {
              pos=ConvertFromXYZtoABC(Framework[CurrentSystem].Atoms[CurrentFramework][j].Position);

              if(fabs(pos.x)<1e-10) pos.x=fabs(pos.x);
              if(fabs(pos.y)<1e-10) pos.y=fabs(pos.y);
              if(fabs(pos.z)<1e-10) pos.z=fabs(pos.z);

              snprintf(AtomName,10,"%-6s",PseudoAtoms[Framework[CurrentSystem].Atoms[CurrentFramework][j].Type].PrintToPDBName);
              fprintf(FilePtr,"%4d %-4s  %9.5f %9.5f %9.5f %4d%4d%4d%4d%4d%4d%4d%4d %7.3f\n",
                      SerialNumber++,
                      AtomName,
                      (double)pos.x,
                      (double)pos.y,
                      (double)pos.z,
                      0,0,0,0,0,0,0,0,
                      (double)0.0);
            }
          }
        }
        else
        {
          fprintf(FilePtr,"%4d   %d %s\n",Framework[CurrentSystem].NumberOfAsymmetricAtoms[CurrentFramework],0,"Created by Raspa-1.0");
          fprintf(FilePtr,"     0 %s         : %s\n",Framework[CurrentSystem].Name[0],
                  Framework[CurrentSystem].Name[0]);
          SerialNumber=1;
          for(j=0;j<Framework[CurrentSystem].NumberOfAsymmetricAtoms[CurrentFramework];j++)
          {
            if(PseudoAtoms[Framework[CurrentSystem].AtomsAsymmetric[CurrentFramework][j].Type].PrintToPDB)
            {
              pos=Framework[CurrentSystem].AtomsAsymmetric[CurrentFramework][j].Position;

              if(fabs(pos.x)<1e-10) pos.x=fabs(pos.x);
              if(fabs(pos.y)<1e-10) pos.y=fabs(pos.y);
              if(fabs(pos.z)<1e-10) pos.z=fabs(pos.z);

              snprintf(AtomName,10,"%-6s",PseudoAtoms[Framework[CurrentSystem].AtomsAsymmetric[CurrentFramework][j].Type].PrintToPDBName);
              fprintf(FilePtr,"%4d %-4s  %9.5f %9.5f %9.5f %4d%4d%4d%4d%4d%4d%4d%4d %7.3f\n",
                      SerialNumber++,
                      AtomName,
                      (double)pos.x,
                      (double)pos.y,
                      (double)pos.z,
                      0,0,0,0,0,0,0,0,
                      (double)0.0);
            }
          }
        }
        fclose(FilePtr);
      }
    }
  }
  free(AtomId);
}

/*********************************************************************************************************
 * Name       | WriteFrameworkDefinitionPDB                                                              *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | writes a PDB-file in 'P1' (no symmetry)                                                  *
 * Parameters | -                                                                                        *
 *********************************************************************************************************/

void WriteFrameworkDefinitionPDB(char *string)
{
  int j,Type;
  char buffer[256];
  FILE *FilePtr;
  FILE *FilePtrAll;
  char RecordName[7]="ATOM  ";  // ATOM record
  int  SerialNumber;            // Atom serial number
  int  SerialNumberAll;         // Atom serial number
  char AtomName[10]="      C";  // Atom Name
  char RemotenessIndicator=' '; // Remoteness indicator
  char BranchDesignator=' ';    // Branch designator
  char AltLoc=' ';              // Alternate location indicator
  char ResIdueName[4]="MOL";    // ResIdue Name
  char ChainId=' ';             // Chain Identifier
  char ResSeq[5]="    ";        // ResIdue sequence number
  char iCode=' ';               // code for insertion of resIdues
  REAL Occupancy=1.0;           // Occupancy
  REAL Temp=0.0;                // Temperature factor
  char SegID[5]="    ";         // Segment Identifier, left-justified
  char Element[10]="         "; // Element symbol, right-justified
  char charge[3]="  ";          // Charge
  VECTOR r;
  REAL MovieScale=1.0;

  if (STREAM)
  {
    fprintf(stderr, "File writing not allowed in streaming mode!");
    return;
  }

  mkdir("Movies",S_IRWXU);

  for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
  {
    if(Framework[CurrentSystem].FrameworkModel!=NONE)
    {
      sprintf(buffer,"Movies/System_%d",CurrentSystem);
      mkdir(buffer,S_IRWXU);

      sprintf(buffer,"Movies/System_%d/Framework_%s%s.pdb",CurrentSystem,string,FileNameAppend);
      FilePtrAll=fopen(buffer,"w");
      SerialNumberAll=1;

      for(CurrentFramework=0;CurrentFramework<Framework[CurrentSystem].NumberOfFrameworks;CurrentFramework++)
      {
        // write pdb-format
        sprintf(buffer,"Movies/System_%d/Framework_%d_%s%s.pdb",CurrentSystem,CurrentFramework,string,FileNameAppend);
        FilePtr=fopen(buffer,"w");
        SerialNumber=1;
        fprintf(FilePtrAll,"REMARK   Raspa-1.0 PDB file\n");
        fprintf(FilePtrAll,"CRYST1%9.3lf%9.3lf%9.3lf%7.2lf%7.2lf%7.2lf\n",
                (double)(MovieScale*BoxProperties[CurrentSystem].ax),
                (double)(MovieScale*BoxProperties[CurrentSystem].ay),
                (double)(MovieScale*BoxProperties[CurrentSystem].az),
                (double)AlphaAngle[CurrentSystem]*RAD2DEG,(double)BetaAngle[CurrentSystem]*RAD2DEG,(double)GammaAngle[CurrentSystem]*RAD2DEG);
        fprintf(FilePtr,"REMARK   Raspa-1.0 PDB file\n");
        fprintf(FilePtr,"CRYST1%9.3lf%9.3lf%9.3lf%7.2lf%7.2lf%7.2lf\n",
                (double)(MovieScale*BoxProperties[CurrentSystem].ax),
                (double)(MovieScale*BoxProperties[CurrentSystem].ay),
                (double)(MovieScale*BoxProperties[CurrentSystem].az),
                (double)AlphaAngle[CurrentSystem]*RAD2DEG,(double)BetaAngle[CurrentSystem]*RAD2DEG,(double)GammaAngle[CurrentSystem]*RAD2DEG);
        for(j=0;j<Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];j++)
        {
          Type=Framework[CurrentSystem].Atoms[CurrentFramework][j].Type;
          if(PseudoAtoms[Framework[CurrentSystem].Atoms[CurrentFramework][j].Type].PrintToPDB)
          {
            snprintf(AtomName,10,"%2s",PseudoAtoms[Framework[CurrentSystem].Atoms[CurrentFramework][j].Type].PrintToPDBName);
            snprintf(Element,10,"%2s",PseudoAtoms[Framework[CurrentSystem].Atoms[CurrentFramework][j].Type].PrintToPDBName);
            r=Framework[CurrentSystem].Atoms[CurrentFramework][j].Position;
            fprintf(FilePtr,"%s%5d %2s%c%c%c%3s %c%4s%c   %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf%6s%-4s%2s%2s\n",
                    RecordName,SerialNumber++,AtomName,RemotenessIndicator,BranchDesignator,AltLoc,ResIdueName,
                    ChainId,ResSeq,iCode,(double)(MovieScale*r.x),(double)(MovieScale*r.y),(double)(MovieScale*r.z),(double)Occupancy,(double)Temp,"      ",SegID,Element,charge);
            fprintf(FilePtrAll,"%s%5d %2s%c%c%c%3s %c%4s%c   %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf%6s%-4s%2s%2s\n",
                    RecordName,SerialNumberAll++,AtomName,RemotenessIndicator,BranchDesignator,AltLoc,ResIdueName,
                    ChainId,ResSeq,iCode,(double)(MovieScale*r.x),(double)(MovieScale*r.y),(double)(MovieScale*r.z),(double)Occupancy,(double)Temp,"      ",SegID,Element,charge);
          }
        }
        fclose(FilePtr);
      }
      fclose(FilePtrAll);
    }
  }
}

/*********************************************************************************************************
 * Name       | WriteFrameworkDefinitionGulp                                                             *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | writes a GULP-file in 'P1' (no symmetry)                                                 *
 * Parameters | -                                                                                        *
 *********************************************************************************************************/

void WriteFrameworkDefinitionGulp(char *string)
{
  int i,j;
  VECTOR pos;
  char buffer[256];
  FILE *FilePtr;
  int  SerialNumber;                             // Atom serial number
  int *AtomId;
  int nr_cores;
  int count;
  int TypeA,TypeB,TypeC;

  if (STREAM)
  {
    fprintf(stderr, "File writing not allowed in streaming mode!");
    return;
  }

  AtomId=(int*)calloc(NumberOfPseudoAtoms,sizeof(int));
  mkdir("Movies",S_IRWXU);

  for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
  {
    if(Framework[CurrentSystem].FrameworkModel!=NONE)
    {
      for(CurrentFramework=0;CurrentFramework<Framework[CurrentSystem].NumberOfFrameworks;CurrentFramework++)
      {
        sprintf(buffer,"Movies/System_%d",CurrentSystem);
        mkdir(buffer,S_IRWXU);

        sprintf(buffer,"Movies/System_%d/Framework_%d_%s%s.gulp",CurrentSystem,CurrentFramework,string,FileNameAppend);
        FilePtr=fopen(buffer,"w");
        fprintf(FilePtr,"opti prop conp full nosymmetry phonon\n");
        fprintf(FilePtr,"gmax % g\n",MaxGradientTolerance);
        fprintf(FilePtr,"accuracy %g\n",-log10(EwaldPrecision));
        fprintf(FilePtr,"pressure % g\n",1e9*therm_baro_stats.ExternalPressure[CurrentSystem][0]);
        fprintf(FilePtr,"switch_min rfo gnorm 50.0\n\n");
        fprintf(FilePtr,"cell\n");
        fprintf(FilePtr,"%g %g %g %g %g %g\n",
                (double)(BoxProperties[CurrentSystem].ax),
                (double)(BoxProperties[CurrentSystem].ay),
                (double)(BoxProperties[CurrentSystem].az),
                (double)AlphaAngle[CurrentSystem]*RAD2DEG,(double)BetaAngle[CurrentSystem]*RAD2DEG,(double)GammaAngle[CurrentSystem]*RAD2DEG);
        fprintf(FilePtr,"space\n1\n\n");

        fprintf(FilePtr,"species\n");
        for(i=0;i<NumberOfPseudoAtoms;i++)
          if(NumberOfPseudoAtomsType[CurrentSystem][i]>0)
            fprintf(FilePtr,"%s %s %g\n",PseudoAtoms[i].PrintToPDBName,PseudoAtoms[i].CoreShell==CORE?"core":"shel",PseudoAtoms[i].Charge1);

        count=0;
        for(i=0;i<NumberOfPseudoAtoms;i++)
          for(j=i;j<NumberOfPseudoAtoms;j++)
            if(PotentialType[i][j]==BUCKINGHAM) count++;

        if(count>0)
        {
          fprintf(FilePtr,"buckingham\n");
          for(i=0;i<NumberOfPseudoAtoms;i++)
            for(j=i;j<NumberOfPseudoAtoms;j++)
              if((PotentialType[i][j]==BUCKINGHAM)&&(NumberOfPseudoAtomsType[CurrentSystem][i]>0)&&(NumberOfPseudoAtomsType[CurrentSystem][j]>0))
              {
                 fprintf(FilePtr,"%s %s %s %s %g %g %g %g %g\n",PseudoAtoms[i].PrintToPDBName,PseudoAtoms[i].CoreShell==CORE?"core":"shel",
                       PseudoAtoms[j].PrintToPDBName,PseudoAtoms[j].CoreShell==CORE?"core":"shel",
                       (double)PotentialParms[i][j][0]*ENERGY_TO_EV,
                       (double)(1.0/PotentialParms[i][j][1]),
                       (double)PotentialParms[i][j][2]*ENERGY_TO_EV,
                       0.0,
                       CutOffVDW);
              }
        }

        if(Framework[CurrentSystem].NumberOfCoreShellDefinitions>0)
        {
          fprintf(FilePtr,"spring\n");
          for(i=0;i<Framework[CurrentSystem].NumberOfCoreShellDefinitions;i++)
          {
            TypeA=Framework[CurrentSystem].CoreShellDefinitions[i].A;
            if((Framework[CurrentSystem].BondDefinitionType[i]==CORE_SHELL_SPRING)&&(NumberOfPseudoAtomsType[CurrentSystem][TypeA]>0))
              fprintf(FilePtr,"%s %g\n",
                   PseudoAtoms[TypeA].PrintToPDBName,
                   Framework[CurrentSystem].BondArgumentDefinitions[i][0]*KELVIN_TO_ENERGY*ENERGY_TO_EV);
          }
        }

        if(Framework[CurrentSystem].NumberOfBendDefinitions>0)
        {
          fprintf(FilePtr,"three\n");
          for(i=0;i<Framework[CurrentSystem].NumberOfBendDefinitions;i++)
          {
            TypeA=Framework[CurrentSystem].BendDefinitions[i].A;
            TypeB=Framework[CurrentSystem].BendDefinitions[i].B;
            TypeC=Framework[CurrentSystem].BendDefinitions[i].C;
            if((Framework[CurrentSystem].BendDefinitionType[i]==CORE_SHELL_BEND)&&(NumberOfPseudoAtomsType[CurrentSystem][TypeA]>0)&&
               (NumberOfPseudoAtomsType[CurrentSystem][TypeB]>0)&&(NumberOfPseudoAtomsType[CurrentSystem][TypeC]>0))
              fprintf(FilePtr,"%s %s %s %s %s %s %g %g 1.9 1.9 3.5\n",
                   PseudoAtoms[TypeB].PrintToPDBName,
                   PseudoAtoms[TypeB].CoreShell==CORE?"core":"shel",
                   PseudoAtoms[TypeA].PrintToPDBName,
                   PseudoAtoms[TypeA].CoreShell==CORE?"core":"shel",
                   PseudoAtoms[TypeC].PrintToPDBName,
                   PseudoAtoms[TypeC].CoreShell==CORE?"core":"shel",
                   Framework[CurrentSystem].BendArgumentDefinitions[i][0]*KELVIN_TO_ENERGY*ENERGY_TO_EV,
                   Framework[CurrentSystem].BendArgumentDefinitions[i][1]);
          }
        }
        fprintf(FilePtr,"\n");
        fprintf(FilePtr,"frac\n");
        SerialNumber=1;
        nr_cores=Framework[CurrentSystem].NumberOfAtoms[CurrentFramework]-Framework[CurrentSystem].NumberOfCoreShells[CurrentFramework];
        for(j=0;j<Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];j++)
        {
          pos=ConvertFromXYZtoABC(Framework[CurrentSystem].Atoms[CurrentFramework][j].Position);

          if(fabs(pos.x)<1e-10) pos.x=fabs(pos.x);
          if(fabs(pos.y)<1e-10) pos.y=fabs(pos.y);
          if(fabs(pos.z)<1e-10) pos.z=fabs(pos.z);

          fprintf(FilePtr,"%-6s %s %9.5lf %9.5lf %9.5lf\n",
                  PseudoAtoms[Framework[CurrentSystem].Atoms[CurrentFramework][j].Type].PrintToPDBName,
                  j<nr_cores?"core":"shel",
                  (double)pos.x,
                  (double)pos.y,
                  (double)pos.z);
        }
        fclose(FilePtr);
      }
    }
  }
  free(AtomId);
}

/*********************************************************************************************************
 * Name       | WriteFrameworkDefinitionVASP                                                             *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | writes a VASP-file in 'P1' (no symmetry)                                                 *
 * Parameters | -                                                                                        *
 *********************************************************************************************************/

void WriteFrameworkDefinitionVASP(char *string)
{
  int i,j;
  VECTOR pos;
  char buffer[256];
  FILE *FilePtr;
  int *AtomId;
  int count;

  if (STREAM)
  {
    fprintf(stderr, "File writing not allowed in streaming mode!");
    return;
  }

  AtomId=(int*)calloc(NumberOfPseudoAtoms,sizeof(int));
  mkdir("Movies",S_IRWXU);

  for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
  {
    if(Framework[CurrentSystem].FrameworkModel!=NONE)
    {
      for(CurrentFramework=0;CurrentFramework<Framework[CurrentSystem].NumberOfFrameworks;CurrentFramework++)
      {
        sprintf(buffer,"Movies/System_%d",CurrentSystem);
        mkdir(buffer,S_IRWXU);

        sprintf(buffer,"Movies/System_%d/Framework_%d_%s%s.vasp",CurrentSystem,CurrentFramework,string,FileNameAppend);
        FilePtr=fopen(buffer,"w");

        fprintf(FilePtr,"# %s\n",Framework[CurrentSystem].Name[0]);
/*
        for(i=0;i<NumberOfPseudoAtoms;i++)
        {
          count=0;
          for(j=0;j<Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];j++)
            if((Framework[CurrentSystem].Atoms[CurrentFramework][j].Type==i)&&(PseudoAtoms[Framework[CurrentSystem].Atoms[CurrentFramework][j].Type].PrintToPDB))
              count++;
          if(count>0) fprintf(FilePtr,"%s (%d) ",PseudoAtoms[i].ChemicalElement,count);
        }
        fprintf(FilePtr,"\n");
*/

        fprintf(FilePtr,"%f\n",1.0);
        fprintf(FilePtr,"%20.16f %20.16f %20.10f\n",Box[CurrentSystem].ax,Box[CurrentSystem].ay,Box[CurrentSystem].az);
        fprintf(FilePtr,"%20.16f %20.16f %20.10f\n",Box[CurrentSystem].bx,Box[CurrentSystem].by,Box[CurrentSystem].bz);
        fprintf(FilePtr,"%20.16f %20.16f %20.10f\n",Box[CurrentSystem].cx,Box[CurrentSystem].cy,Box[CurrentSystem].cz);

        for(i=0;i<NumberOfPseudoAtoms;i++)
          for(j=0;j<Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];j++)
          {
            if((Framework[CurrentSystem].Atoms[CurrentFramework][j].Type==i)&&(PseudoAtoms[Framework[CurrentSystem].Atoms[CurrentFramework][j].Type].PrintToPDB))
            {
              fprintf(FilePtr,"%s ",PseudoAtoms[i].ChemicalElement);
              break;
            }
          }
        fprintf(FilePtr,"\n");


        for(i=0;i<NumberOfPseudoAtoms;i++)
        {
          count=0;
          for(j=0;j<Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];j++)
            if((Framework[CurrentSystem].Atoms[CurrentFramework][j].Type==i)&&(PseudoAtoms[Framework[CurrentSystem].Atoms[CurrentFramework][j].Type].PrintToPDB))
              count++;
          if(count>0) fprintf(FilePtr,"%d ",count);
        }
        fprintf(FilePtr,"\n");
        fprintf(FilePtr,"Direct\n");
        for(i=0;i<NumberOfPseudoAtoms;i++)
        {
          for(j=0;j<Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];j++)
          {
            if((Framework[CurrentSystem].Atoms[CurrentFramework][j].Type==i)&&(PseudoAtoms[Framework[CurrentSystem].Atoms[CurrentFramework][j].Type].PrintToPDB))
            {
              pos=ConvertFromXYZtoABC(Framework[CurrentSystem].Atoms[CurrentFramework][j].Position);

              if(fabs(pos.x)<1e-10) pos.x=fabs(pos.x);
              if(fabs(pos.y)<1e-10) pos.y=fabs(pos.y);
              if(fabs(pos.z)<1e-10) pos.z=fabs(pos.z);

              fprintf(FilePtr,"%4.16lf %4.16lf %4.16lf\n",
                      (double)pos.x,
                      (double)pos.y,
                      (double)pos.z);
            }
          }
        }
        fclose(FilePtr);
      }
    }
  }
  free(AtomId);
}

/*********************************************************************************************************
 * Name       | WriteFrameworkDefinitionTinker                                                           *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | writes a Tinker-file in 'P1' (no symmetry)                                               *
 * Parameters | -                                                                                        *
 *********************************************************************************************************/

void WriteFrameworkDefinitionTinker(char *string)
{
  int i,j,Type;
  char buffer[256];
  FILE *FilePtr;
  int  SerialNumber;           // Atom serial number
  char AtomName[10]="      C"; // Atom Name
  int *AtomId;
  VECTOR r;

  if (STREAM)
  {
    fprintf(stderr, "File writing not allowed in streaming mode!");
    return;
  }

  AtomId=(int*)calloc(NumberOfPseudoAtoms,sizeof(int));
  mkdir("Movies",S_IRWXU);

  for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
  {
    if(Framework[CurrentSystem].FrameworkModel!=NONE)
    {
      for(CurrentFramework=0;CurrentFramework<Framework[CurrentSystem].NumberOfFrameworks;CurrentFramework++)
      {
        sprintf(buffer,"Movies/System_%d",CurrentSystem);
        mkdir(buffer,S_IRWXU);

        // write tinker xyz-format
        for(j=0;j<NumberOfPseudoAtoms;j++)
          AtomId[j]=0;
        sprintf(buffer,"Movies/System_%d/Framework_%d_%s%s.txyz",CurrentSystem,CurrentFramework,string,FileNameAppend);
        FilePtr=fopen(buffer,"w");
        fprintf(FilePtr,"%8d %s\n",Framework[CurrentSystem].NumberOfAtoms[CurrentFramework],Framework[CurrentSystem].Name[CurrentFramework]);
        SerialNumber=1;
        for(j=0;j<Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];j++)
        {
          Type=Framework[CurrentSystem].Atoms[CurrentFramework][j].Type;
          if(PseudoAtoms[Type].PrintToPDB)
          {
            r=Framework[CurrentSystem].Atoms[CurrentFramework][j].Position;
            sprintf(AtomName,"%3s",PseudoAtoms[Type].PrintToPDBName);
            fprintf(FilePtr,"%6d%3s %12.6lf %12.6lf %12.6lf %6d",
                    SerialNumber++,AtomName,(double)r.x,(double)r.y,(double)r.z,PseudoAtoms[Type].TinkerType);
            for(i=0;i<Framework[CurrentSystem].Connectivity[CurrentFramework][j];i++)
              fprintf(FilePtr," %4d",
                      Framework[CurrentSystem].Neighbours[CurrentFramework][j][i]+1);
            fprintf(FilePtr,"\n");
          }
        }
        fclose(FilePtr);
      }
    }
  }
  free(AtomId);
}

/*********************************************************************************************************
 * Name       | ExpandAsymmetricFrameworkToFullFramework                                                 *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Expands the asymmetric framework atoms to the full framework.                            *
 * Parameters | -                                                                                        *
 *********************************************************************************************************/

void ExpandAsymmetricFrameworkToFullFramework(void)
{
  int i,j,k;
  int k1,k2,k3;
  int index,insert_index;
  int Type;
  VECTOR dr,ds,pos;
  REAL rr;
  REAL charge;

  // allocate memory for an initial maximum of 1024 atoms, reallocate more when needed
  // initialize the mass of the framework
  Framework[CurrentSystem].NumberOfAtoms[CurrentFramework]=0;
  Framework[CurrentSystem].NumberOfCharges[CurrentFramework]=0;
  Framework[CurrentSystem].MaxNumberOfAtoms[CurrentFramework]=1024;
  Framework[CurrentSystem].Atoms[CurrentFramework]=(ATOM*)calloc(Framework[CurrentSystem].MaxNumberOfAtoms[CurrentFramework],sizeof(ATOM));
  Framework[CurrentSystem].FrameworkMassPerComponent[CurrentFramework]=0.0;

  // get the size of the spacegroup
  pos.x=pos.y=pos.z=0.0;
  SpaceGroupSymmetry(Framework[CurrentSystem].SpaceGroupIdentifier[CurrentFramework],pos);

  // loop over the unit cells, then space group elements, and then the asymmetric atoms
  // this will list the asymmetric atoms together
  insert_index=0;
  for(k1=0;k1<NumberOfUnitCells[CurrentSystem].x;k1++)
    for(k2=0;k2<NumberOfUnitCells[CurrentSystem].y;k2++)
      for(k3=0;k3<NumberOfUnitCells[CurrentSystem].z;k3++)
      {
        for(j=0;j<SpaceGroupSize;j++)
        {
          for(i=0;i<Framework[CurrentSystem].NumberOfAsymmetricAtoms[CurrentFramework];i++)
          {
            Type=Framework[CurrentSystem].AtomsAsymmetric[CurrentFramework][i].Type;
            pos=Framework[CurrentSystem].AtomsAsymmetric[CurrentFramework][i].Position;
            charge=Framework[CurrentSystem].AtomsAsymmetric[CurrentFramework][i].Charge;

            // get all the symmetric positions for the current asymmetric atom position
            SpaceGroupSymmetry(Framework[CurrentSystem].SpaceGroupIdentifier[CurrentFramework],pos);

            if(!Framework[CurrentSystem].ReadCIFAsCartesian)
            {
              if(Framework[CurrentSystem].RestrictFrameworkAtomsToBox)
              {
                // apply boundary condition
                SpaceGroupElement[j].x-=NINT(SpaceGroupElement[j].x);
                SpaceGroupElement[j].y-=NINT(SpaceGroupElement[j].y);
                SpaceGroupElement[j].z-=NINT(SpaceGroupElement[j].z);

                if(SpaceGroupElement[j].x<0.0) SpaceGroupElement[j].x+=1.0;
                if(SpaceGroupElement[j].y<0.0) SpaceGroupElement[j].y+=1.0;
                if(SpaceGroupElement[j].z<0.0) SpaceGroupElement[j].z+=1.0;
              }
            }

            // adjust position taking the number of unit cells into account
            pos.x=(SpaceGroupElement[j].x+(REAL)k1)/(REAL)NumberOfUnitCells[CurrentSystem].x;
            pos.y=(SpaceGroupElement[j].y+(REAL)k2)/(REAL)NumberOfUnitCells[CurrentSystem].y;
            pos.z=(SpaceGroupElement[j].z+(REAL)k3)/(REAL)NumberOfUnitCells[CurrentSystem].z;

            // check if the symmetry position is new or if it is already generated
            index=-1;
            for(k=0;k<insert_index;k++)
            {
              ds.x=Framework[CurrentSystem].Atoms[CurrentFramework][k].Position.x-pos.x;
              ds.y=Framework[CurrentSystem].Atoms[CurrentFramework][k].Position.y-pos.y;
              ds.z=Framework[CurrentSystem].Atoms[CurrentFramework][k].Position.z-pos.z;
              ds.x-=NINT(ds.x);
              ds.y-=NINT(ds.y);
              ds.z-=NINT(ds.z);

              dr=ConvertFromABCtoXYZ(ds); 
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
              if(rr<1e-3) index=k;
            }

            // the symmetry position is new -> insert it as a framework atom
            if(index<0)
            {
               // reallocate memory if needed
               if(insert_index>=Framework[CurrentSystem].MaxNumberOfAtoms[CurrentFramework])
               {
                 Framework[CurrentSystem].MaxNumberOfAtoms[CurrentFramework]+=1024;
                 Framework[CurrentSystem].Atoms[CurrentFramework]=(ATOM*)realloc(Framework[CurrentSystem].Atoms[CurrentFramework],
                   Framework[CurrentSystem].MaxNumberOfAtoms[CurrentFramework]*sizeof(ATOM));
               }

              // insert atom into framework
              Framework[CurrentSystem].Atoms[CurrentFramework][insert_index].Position=pos;
              Framework[CurrentSystem].Atoms[CurrentFramework][insert_index].Type=Type;

              // Fill in charge from the pseudo-atom definition
              Framework[CurrentSystem].Atoms[CurrentFramework][insert_index].Charge=charge;

              // adjust the counter for the type of pseudo atoms
              NumberOfPseudoAtomsCount[CurrentSystem][Type]++;
              NumberOfPseudoAtomsType[CurrentSystem][Type]++;

              if(PseudoAtoms[Type].HasCharges)
                Framework[CurrentSystem].NumberOfCharges[CurrentFramework]++;

              // add the mass of the atom to the total framework mass
              Framework[CurrentSystem].FrameworkMassPerComponent[CurrentFramework]+=PseudoAtoms[Type].Mass;

              // increase the number of atoms by one
              insert_index++;
            }
          }
    }
  }

  // update bookkeeping for framework properties
  Framework[CurrentSystem].NumberOfAtoms[CurrentFramework]=insert_index;
  Framework[CurrentSystem].NumberOfUnitCellAtoms[CurrentFramework]=insert_index/(NumberOfUnitCells[CurrentSystem].x*
                                              NumberOfUnitCells[CurrentSystem].y*NumberOfUnitCells[CurrentSystem].z);

  Framework[CurrentSystem].FrameworkDensityPerComponent[CurrentFramework]=
          Framework[CurrentSystem].FrameworkMassPerComponent[CurrentFramework]/
         (Volume[CurrentSystem]*CUBE(ANGSTROM)*AVOGADRO_CONSTANT);

  Framework[CurrentSystem].TotalNumberOfAtoms+=Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];
  Framework[CurrentSystem].TotalNumberOfCharges+=Framework[CurrentSystem].NumberOfCharges[CurrentFramework];
  Framework[CurrentSystem].TotalNumberOfUnitCellAtoms+=Framework[CurrentSystem].NumberOfUnitCellAtoms[CurrentFramework];
  Framework[CurrentSystem].FrameworkDensity+=Framework[CurrentSystem].FrameworkDensityPerComponent[CurrentFramework];
  Framework[CurrentSystem].FrameworkMass+=Framework[CurrentSystem].FrameworkMassPerComponent[CurrentFramework];


  // transform from fractional coordinates (abc) to Cartesian coordinates (xyz)
  for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];i++)
  {
    pos=Framework[CurrentSystem].Atoms[CurrentFramework][i].Position;
    pos.x+=Framework[CurrentSystem].ShiftUnitCell[CurrentFramework].x;
    pos.y+=Framework[CurrentSystem].ShiftUnitCell[CurrentFramework].y;
    pos.z+=Framework[CurrentSystem].ShiftUnitCell[CurrentFramework].z;
    if(!Framework[CurrentSystem].ReadCIFAsCartesian)
      Framework[CurrentSystem].Atoms[CurrentFramework][i].Position=ConvertFromABCtoXYZ(pos);
  }
}

/*********************************************************************************************************
 * Name       | ExpandAsymmetricIonsToFullIons                                                           *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Expands the asymmetric ion sitings to the full simulation cell.                          *
 * Parameters | -                                                                                        *
 *********************************************************************************************************/

void ExpandAsymmetricIonsToFullCell(void)
{
  int i,j,k;
  int k1,k2,k3;
  int index,insert_index;
  int Type;
  VECTOR dr,ds,pos;
  REAL rr;

  // allocate memory for an initial maximum of 1024 atoms, reallocate more when needed
  // initialize the mass of the framework
  Framework[CurrentSystem].NumberOfIons=0;
  Framework[CurrentSystem].MaxNumberOfIons=1024;
  Framework[CurrentSystem].Ions=(ATOM*)calloc(Framework[CurrentSystem].MaxNumberOfIons,sizeof(ATOM));

  // get the size of the spacegroup
  pos.x=pos.y=pos.z=0.0;
  SpaceGroupSymmetry(Framework[CurrentSystem].SpaceGroupIdentifier[CurrentFramework],pos);

  // loop over the unit cells, then space group elements, and then the asymmetric atoms
  // this will list the asymmetric atoms together
  insert_index=0;
  for(k1=0;k1<NumberOfUnitCells[CurrentSystem].x;k1++)
    for(k2=0;k2<NumberOfUnitCells[CurrentSystem].y;k2++)
      for(k3=0;k3<NumberOfUnitCells[CurrentSystem].z;k3++)
      {
        for(j=0;j<SpaceGroupSize;j++)
        {
          for(i=0;i<Framework[CurrentSystem].NumberOfAsymmetricIons;i++)
          {
            Type=Framework[CurrentSystem].IonsAsymmetric[i].Type;
            pos=Framework[CurrentSystem].IonsAsymmetric[i].Position;

            // get all the symmetric positions for the current asymmetric atom position
            SpaceGroupSymmetry(Framework[CurrentSystem].SpaceGroupIdentifier[CurrentFramework],pos);

            // apply boundary condition
            SpaceGroupElement[j].x-=NINT(SpaceGroupElement[j].x);
            SpaceGroupElement[j].y-=NINT(SpaceGroupElement[j].y);
            SpaceGroupElement[j].z-=NINT(SpaceGroupElement[j].z);

            if(SpaceGroupElement[j].x<0.0) SpaceGroupElement[j].x+=1.0;
            if(SpaceGroupElement[j].y<0.0) SpaceGroupElement[j].y+=1.0;
            if(SpaceGroupElement[j].z<0.0) SpaceGroupElement[j].z+=1.0;

            // adjust position taking the number of unit cells into account
            pos.x=(SpaceGroupElement[j].x+(REAL)k1)/(REAL)NumberOfUnitCells[CurrentSystem].x;
            pos.y=(SpaceGroupElement[j].y+(REAL)k2)/(REAL)NumberOfUnitCells[CurrentSystem].y;
            pos.z=(SpaceGroupElement[j].z+(REAL)k3)/(REAL)NumberOfUnitCells[CurrentSystem].z;

            // check if the symmetry position is new or if it is already generated
            index=-1;
            for(k=0;k<insert_index;k++)
            {
              ds.x=Framework[CurrentSystem].Ions[k].Position.x-pos.x;
              ds.y=Framework[CurrentSystem].Ions[k].Position.y-pos.y;
              ds.z=Framework[CurrentSystem].Ions[k].Position.z-pos.z;
              ds.x-=NINT(ds.x);
              ds.y-=NINT(ds.y);
              ds.z-=NINT(ds.z);

              dr=ConvertFromABCtoXYZ(ds);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
              if(rr<1e-3) index=k;
            }

            // the symmetry position is new -> insert it as a framework atom
            if(index<0)
            {
               // reallocate memory if needed
               if(insert_index>=Framework[CurrentSystem].MaxNumberOfIons)
               {
                 Framework[CurrentSystem].MaxNumberOfIons+=1024;
                 Framework[CurrentSystem].Ions=(ATOM*)realloc(Framework[CurrentSystem].Ions,
                   Framework[CurrentSystem].MaxNumberOfIons*sizeof(ATOM));
               }

              // insert atom into framework
              Framework[CurrentSystem].Ions[insert_index].Position=pos;
              Framework[CurrentSystem].Ions[insert_index].Type=Type;
              Framework[CurrentSystem].Ions[insert_index].AssymetricType=i;

              crystallographic_stats[CurrentSystem].NumberOfCationSites[i]++;

              // increase the number of atoms by one
              insert_index++;
            }
          }
    }
  }

  // update bookkeeping for framework properties
  Framework[CurrentSystem].NumberOfIons=insert_index;

  // transform from fractional coordinates (abc) to Cartesian coordinates (xyz)
  for(i=0;i<Framework[CurrentSystem].NumberOfIons;i++)
  {
    pos=Framework[CurrentSystem].Ions[i].Position;
    pos.x+=Framework[CurrentSystem].ShiftUnitCell[CurrentFramework].x;
    pos.y+=Framework[CurrentSystem].ShiftUnitCell[CurrentFramework].y;
    pos.z+=Framework[CurrentSystem].ShiftUnitCell[CurrentFramework].z;
    Framework[CurrentSystem].Ions[i].Position=ConvertFromABCtoXYZ(pos);
  }
}

int **Stack;
int *StackSize;
int *MaxStackSize;

void AllocateStacks(int nr_stacks)
{
  Stack=(int**)calloc(nr_stacks,sizeof(int*));
  StackSize=(int*)calloc(nr_stacks,sizeof(int));
  MaxStackSize=(int*)calloc(nr_stacks,sizeof(int));
}

void AllocateStack(int stack)
{
  MaxStackSize[stack]=1024;
  Stack[stack]=(int*)calloc(MaxStackSize[stack],sizeof(int));
}

void InitializeStack(int stack)
{
  StackSize[stack]=0;
}

void Push(int stack,int A)
{
  Stack[stack][StackSize[stack]++]=A;
}

int Pop(int stack)
{
  if(StackSize[stack]>0) return Stack[stack][StackSize[stack]--];
  return -1;
}

int IsElementInStack(int stack,int element)
{
  int i;

  for(i=0;i<StackSize[stack];i++)
    if(Stack[stack][i]==element) return TRUE;
  return FALSE;
}

void SetDifference(int stack,int stack1,int stack2)
{
  int i;
  int element;

  StackSize[stack]=0;
  for(i=0;i<StackSize[stack1];i++)
  {
    element=Stack[stack1][i];
    if(!IsElementInStack(stack2,element))
     Push(stack,element);
  }
}

void PrintStack(int stack)
{
  int i;

  for(i=0;i<StackSize[stack];i++)
    fprintf(stderr, "%d ",Stack[stack][i]);
  fprintf(stderr, "\n");
}

void DeallocateStack(int stack)
{
  free(Stack[stack]);
}

void DeallocateStacks(int stack)
{
  free(Stack);
  free(StackSize);
  free(MaxStackSize);
}


enum{DISORDER,NO_DISORDER,SELECTED,DELETED};

void FollowHydrogens(int A)
{
  int j,B;

  // stop criteria for recursion
  if(IsElementInStack(0,A)) return;

  Push(0,A);

  // Recursively handle all neighboring atoms
  for(j=0;j<Framework[CurrentSystem].Connectivity[0][A];j++)
  {
    B=Framework[CurrentSystem].Neighbours[0][A][j];
    if(strcmp(PseudoAtoms[Framework[CurrentSystem].Atoms[0][B].Type].PrintToPDBName,"H")==0)
      FollowHydrogens(B);
  }
}

void FollowFrameworkGetSmallSet(int A)
{
  int j,B;
  int min,fmin;

  // stop criteria for recursion
  //if(IsElementInStack(0,A)||Framework[CurrentSystem].Atoms[0][A].AssymetricType<0) return;
  if(IsElementInStack(0,A)||(Framework[CurrentSystem].Atoms[0][A].CreationState==NO_DISORDER)) return;

  Push(0,A);

  // Find the framework with the lowest connectivity
  // Sometimes a certain choice can lead to two configurations being too close together, i.e. they appear bonded
  // By taking the lowest connectivity we follow the configurations as if they were not-bonded
  min=100;
  fmin=100;
  for(j=0;j<Framework[CurrentSystem].NumberOfFrameworks;j++)
  {
    if(Framework[CurrentSystem].Connectivity[j][A]<min)
    {
      min=Framework[CurrentSystem].Connectivity[j][A];
      fmin=j;
    }
  }

  // Recursively handle all neighboring atoms
  for(j=0;j<Framework[CurrentSystem].Connectivity[fmin][A];j++)
  {
    B=Framework[CurrentSystem].Neighbours[fmin][A][j];
    FollowFrameworkGetSmallSet(B);
  }
}


void FollowFrameworkGetLargestSet(int A)
{
  int j,B;
  int max,fmax;

  // stop criteria for recursion
  //if(IsElementInStack(0,A)||Framework[CurrentSystem].Atoms[0][A].AssymetricType<0) return;
  if(IsElementInStack(1,A)||(Framework[CurrentSystem].Atoms[0][A].CreationState==NO_DISORDER)) return;


  Push(1,A);

  // Find the framework with the lowest connectivity
  // Sometimes a certain choice can lead to two configurations being too close together, i.e. they appear bonded
  // By taking the lowest conenctivity we follow the configurations as if they were not-bonded
  max=0;
  fmax=0;
  for(j=0;j<Framework[CurrentSystem].NumberOfFrameworks;j++)
  {
    if(Framework[CurrentSystem].Connectivity[j][A]>max)
    {
      max=Framework[CurrentSystem].Connectivity[j][A];
      fmax=j;
    }
  }

  // Recursively handle all neighboring atoms
  for(j=0;j<Framework[CurrentSystem].Connectivity[fmax][A];j++)
  {
    B=Framework[CurrentSystem].Neighbours[fmax][A][j];
    FollowFrameworkGetLargestSet(B);
  }
}

void FollowFrameworkDisorderedSet(int fr,int stack,int A)
{
  int j,B;
  int TypeA,TypeB;

  // stop criteria for recursion: if element is already in set or no disorder for this atom
  if(IsElementInStack(0,A)||(Framework[CurrentSystem].Atoms[0][A].CreationState==NO_DISORDER)) return;

  fprintf(stderr, "add %d (%s) to stack %d\n",A,PseudoAtoms[Framework[CurrentSystem].Atoms[fr][A].Type].Name,stack);
  Push(0,A);
  Push((stack%2)+1,A);


  // Recursively handle all neighboring atoms
  for(j=0;j<Framework[CurrentSystem].Connectivity[fr][A];j++)
  {
    B=Framework[CurrentSystem].Neighbours[fr][A][j];

    TypeA=Framework[CurrentSystem].Atoms[fr][A].Type;
    TypeB=Framework[CurrentSystem].Atoms[fr][B].Type;

    // keep track of crossing intra-hydrogen-bonds
    if((strcasecmp(PseudoAtoms[TypeA].PrintToPDBName,"H")==0)&&
       (strcasecmp(PseudoAtoms[TypeB].PrintToPDBName,"H")==0)&&
       (!IsElementInStack(0,B))&&
       (Framework[CurrentSystem].Atoms[0][B].CreationState!=NO_DISORDER))
    {
      fprintf(stderr, "hydrogen-hydrogen bond between %d and %d\n",A,B);
      FollowFrameworkDisorderedSet(fr,stack+1,B);
    }
    else
      FollowFrameworkDisorderedSet(fr,stack,B);
  }
}


void TypeAtomsInStack(int stack,int fr)
{
  int i,j,A;

  for(i=0;i<StackSize[stack];i++)
  {
    A=Stack[stack][i];

    // set the type to the type of the selected framework
    Framework[CurrentSystem].Atoms[0][A].Type=Framework[CurrentSystem].Atoms[fr][A].Type;

    for(j=0;j<Framework[CurrentSystem].NumberOfFrameworks;j++)
      Framework[CurrentSystem].Atoms[j][A].CreationState=NO_DISORDER;

    Framework[CurrentSystem].Atoms[fr][A].AssymetricType=-1;
    fprintf(stderr, "Type: fr: %d atom: %d\n",fr,A);
  }
}

int SelectFramework(void)
{
  REAL ws,cumw;
  int selected;

  ws=0.0;
  cumw=Framework[CurrentSystem].FrameworkProbability[0];
  selected=0;
  ws=RandomNumber();
  while(cumw<ws)
    cumw+=Framework[CurrentSystem].FrameworkProbability[++selected];

  return selected;
}

/*********************************************************************************************************
 * Name       | GenerateFramework                                                                        *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Generates a P1 framework without disorder from a collection of possible frameworks       *
 * Parameters | -                                                                                        *
 * Note       |                                                                                          *
 *********************************************************************************************************/


void GenerateFramework(void)
{
  int i,j,fr;
  VECTOR posA,posB,dr;
  REAL rr;
  int count;
  int count2;
  int disorder;
  int A;
  FILE *FilePtr;
  int index;
  REAL TotalProbability;
  int selected;
  int selected2;
  REAL SelectionProbability[100];
  int Type;

  for(i=0;i<100;i++)
    SelectionProbability[i]=0.0;

  // create connectivity-list for the frameworks
  MakeConnectivityList();

  // set all atoms with no disorder to type '-1'
  // atoms with disorder are labeled accoording to their frameworknumber
  for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[0];i++)
  {
    // get the position of an atom in the first framework
    posA=Framework[CurrentSystem].Atoms[0][i].Position;
    for(fr=0;fr<Framework[CurrentSystem].NumberOfFrameworks;fr++)
    {
      Framework[CurrentSystem].Atoms[fr][i].AssymetricType=-1;
      Framework[CurrentSystem].Atoms[fr][i].CreationState=NO_DISORDER;
    }

    for(fr=1;fr<Framework[CurrentSystem].NumberOfFrameworks;fr++)
    {
      // get position in other framework
      posB=Framework[CurrentSystem].Atoms[fr][i].Position;
      dr.x=posA.x-posB.x;
      dr.y=posA.y-posB.y;
      dr.z=posA.z-posB.z;
      dr=ApplyBoundaryCondition(dr);
      rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

      // detect disorder (the distance between atoms in the different frameworks is non-zero)
      if(rr>1e-4)
      {
        for(j=0;j<Framework[CurrentSystem].NumberOfFrameworks;j++)
        {
          Framework[CurrentSystem].Atoms[j][i].AssymetricType=j;
          Framework[CurrentSystem].Atoms[j][i].CreationState=DISORDER;
        }
      }
    }
  }

  count=0;
  count2=0;
  for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[0];i++)
    if(Framework[CurrentSystem].Atoms[0][i].CreationState==DISORDER)
      count++;
    else
      count2++;
  fprintf(stderr, "Disorder: %d no-disorder: %d\n",count,count2);


  AllocateStacks(9); // allocate space for 3 stacks
  AllocateStack(0);  // allocate space for stack 0
  AllocateStack(1);  // allocate space for stack 1
  AllocateStack(2);  // allocate space for stack 3
  AllocateStack(3);  // allocate space for stack 4
  AllocateStack(4);  // allocate space for stack 5
  AllocateStack(5);  // allocate space for stack 6
  AllocateStack(6);  // allocate space for stack 6
  AllocateStack(7);  // allocate space for stack 6
  AllocateStack(8);  // allocate space for stack 6

  // get the total probability of the frameworks
  TotalProbability=0.0;
  for(fr=0;fr<Framework[CurrentSystem].NumberOfFrameworks;fr++)
    TotalProbability+=Framework[CurrentSystem].FrameworkProbability[fr];

  // normalize the probability
  if(TotalProbability>0.0)
    for(fr=0;fr<Framework[CurrentSystem].NumberOfFrameworks;fr++)
      Framework[CurrentSystem].FrameworkProbability[fr]/=TotalProbability;

  // cehck whether there is disorder and the first disordered atom as a seed
  disorder=FALSE;
  for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[0];i++)
  {
    if(Framework[CurrentSystem].Atoms[0][i].CreationState==DISORDER)
    {
      disorder=TRUE;
      A=i;
      break;
    }
  }

//  for(A=0;A<Framework[CurrentSystem].NumberOfAtoms[0];A++)
//  {
//    if(Framework[CurrentSystem].Atoms[0][A].CreationState==DISORDER)
//    {
//      // initialize the 3 stacks
//      InitializeStack(0);
//      InitializeStack(1);
//      InitializeStack(2);
//      InitializeStack(3);
//      InitializeStack(4);
//      InitializeStack(5);
//
//      fprintf(stderr, "A: %d\n",A);
//      fprintf(stderr, "====================\n");
//
//      FollowFrameworkDisorderedSet(0,1,A);
//
//      fprintf(stderr, "(size: %d): ",StackSize[0]);
//      PrintStack(0);
//      fprintf(stderr, "(size: %d): ",StackSize[1]);
//      PrintStack(1);
//      fprintf(stderr, "(size: %d): ",StackSize[2]);
//      PrintStack(2);
//      fprintf(stderr, "(size: %d): \n",StackSize[3]);
//      fprintf(stderr, "(size: %d): \n",StackSize[4]);
//      fprintf(stderr, "(size: %d): \n",StackSize[5]);
//      fprintf(stderr, "\n\n\n");
//    }
//  }


  while(disorder==TRUE)
  {
    // initialize the 3 stacks
    InitializeStack(0);
    InitializeStack(1);
    InitializeStack(2);


//    fprintf(stderr, "i %d (size: %d): \n",A,StackSize[0]);
//    PrintStack(0);
//    fprintf(stderr, "i %d (size: %d): \n",A,StackSize[1]);
//    PrintStack(1);
//    fprintf(stderr, "i %d (size: %d): \n",A,StackSize[2]);
//    PrintStack(2);
//    fprintf(stderr, "i %d (size: %d): \n",A,StackSize[3]);
//    PrintStack(3);
//    fprintf(stderr, "i %d (size: %d): \n",A,StackSize[4]);
//    PrintStack(4);
//    fprintf(stderr, "i %d (size: %d): \n",A,StackSize[5]);
//    PrintStack(5);




//    FollowFrameworkGetSmallSet(A);   // stack 0 is the smallest set starting with A
//    FollowFrameworkGetLargestSet(A); // stack 1 is the largest set starting with A
//    SetDifference(2,1,0);            // stack 2 is the difference between A and B

    //printf("i %d (size: %d): \n",A,StackSize[0]);
    //PrintStack(0);

    //printf("i %d (size: %d): \n",A,StackSize[1]);
    //PrintStack(1);

    //printf("i %d (size: %d): \n",A,StackSize[2]);
    //PrintStack(2);

    fprintf(stderr, "size [%d] 0: %d, 1: %d, 2: %d\n",A,StackSize[0],StackSize[1],StackSize[2]);

    if(Framework[CurrentSystem].FrameworkExclusion) // exlusion between two groups of disordered atoms
    {
      selected=SelectFramework();

      FollowFrameworkDisorderedSet(1,1,A);

      SelectionProbability[selected]+=1.0;

      // type the atoms in framework 'fr' as "disorder resolved"
      TypeAtomsInStack(0,selected);

      if(selected==0) selected2=1;
      else selected2=0;

      TypeAtomsInStack(1,selected);
      TypeAtomsInStack(2,selected2);

      SelectionProbability[selected2]+=1.0;
    }
    else
    {
      // select a framework with the correct probability
      selected=SelectFramework();

      SelectionProbability[selected]+=1.0;

      // type the atoms in framework 'fr' as "disorder resolved"
      TypeAtomsInStack(0,selected);
    }

    disorder=FALSE;
    for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[0];i++)
    {
      if(Framework[CurrentSystem].Atoms[0][i].CreationState==DISORDER)
      {
        disorder=TRUE;
        A=i;
        break;
      }
    }
  }

  fprintf(stderr, "Framework[CurrentSystem].NumberOfAtoms[0]: %d\n",Framework[CurrentSystem].NumberOfAtoms[0]);
  fprintf(stderr, "Framework[CurrentSystem].NumberOfAtoms[1]: %d\n",Framework[CurrentSystem].NumberOfAtoms[1]);
  fprintf(stderr, "All disorder removed!\n");

  TotalProbability=0.0;
  for(fr=0;fr<Framework[CurrentSystem].NumberOfFrameworks;fr++)
    TotalProbability+=SelectionProbability[fr];
  for(fr=0;fr<Framework[CurrentSystem].NumberOfFrameworks;fr++)
    fprintf(stderr, "%d selected: %g (%g percent)\n",fr,SelectionProbability[fr],100.0*SelectionProbability[fr]/TotalProbability);

  //for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[0];i++)
  //    fprintf(stderr, "%d %d %d\n",i,Framework[CurrentSystem].Atoms[0][i].AssymetricType,Framework[CurrentSystem].Atoms[1][i].AssymetricType);
  //printf("\n");



  // put selected framework in framework '0'
  for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[0];i++)
  {
    for(fr=0;fr<Framework[CurrentSystem].NumberOfFrameworks;fr++)
    {
      if(Framework[CurrentSystem].Atoms[fr][i].AssymetricType==-1)
      {
        Framework[CurrentSystem].Atoms[0][i].Position=Framework[CurrentSystem].Atoms[fr][i].Position;
        Framework[CurrentSystem].Atoms[0][i].Type=Framework[CurrentSystem].Atoms[fr][i].Type;
        Framework[CurrentSystem].Atoms[0][i].AssymetricType=-1;
        break;
      }
    }
  }

  // make sure all frameworks are now the same
  for(fr=1;fr<Framework[CurrentSystem].NumberOfFrameworks;fr++)
    for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[fr];i++)
    {
      Framework[CurrentSystem].Atoms[fr][i].Position=Framework[CurrentSystem].Atoms[0][i].Position;
      Framework[CurrentSystem].Atoms[fr][i].Type=Framework[CurrentSystem].Atoms[0][i].Type;
      Framework[CurrentSystem].Atoms[fr][i].AssymetricType=-1;
    }

  //Framework[CurrentSystem].NumberOfFrameworks=1;
  MakeConnectivityList();
  //PrintConnectivityList();


  // remove disorder of hydrogens
  // ZIF-8: 6 hydrogens connected to a single carbon, 3 hydrogens needs to be removed
  //        strategy is to find a 6 ring of hydrogens and either delete the odd or the even hdyrogens.
  //        the remaining hydrogens have the proper 120 degree angles.
  if(Framework[CurrentSystem].RemoveHydrogenDisorder)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[0];i++)
      Framework[CurrentSystem].Atoms[0][i].CreationState=DISORDER;

    MakeConnectivityList();

    for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[0];i++)
    {
      if(Framework[CurrentSystem].Atoms[0][i].CreationState==DISORDER)
      {
        Type=Framework[CurrentSystem].Atoms[0][i].Type;
        if(strcmp(PseudoAtoms[Type].PrintToPDBName,"H")==0) // select hydrogen
        {
          InitializeStack(0);

          // follow the chain of hydrogens
          FollowHydrogens(i);

          fprintf(stderr, "hydro i %d  %s (size: %d): ",i,PseudoAtoms[Type].Name,StackSize[0]);
          PrintStack(0);

          switch(StackSize[0])
          {
            case 2:
              fprintf(stderr, "Deleted hydrogen\n");
              if(RandomNumber()<0.5)
              {
                Framework[CurrentSystem].Atoms[0][Stack[0][0]].CreationState=DELETED;
                Framework[CurrentSystem].Atoms[0][Stack[0][1]].CreationState=NO_DISORDER;
              }
              else
              {
                Framework[CurrentSystem].Atoms[0][Stack[0][1]].CreationState=DELETED;
                Framework[CurrentSystem].Atoms[0][Stack[0][0]].CreationState=NO_DISORDER;
              }
              break;
            case 6:
              if(RandomNumber()<0.5)
              {
                Framework[CurrentSystem].Atoms[0][Stack[0][0]].CreationState=DELETED;
                Framework[CurrentSystem].Atoms[0][Stack[0][2]].CreationState=DELETED;
                Framework[CurrentSystem].Atoms[0][Stack[0][4]].CreationState=DELETED;
                Framework[CurrentSystem].Atoms[0][Stack[0][1]].CreationState=NO_DISORDER;
                Framework[CurrentSystem].Atoms[0][Stack[0][3]].CreationState=NO_DISORDER;
                Framework[CurrentSystem].Atoms[0][Stack[0][5]].CreationState=NO_DISORDER;
              }
              else
              {
                Framework[CurrentSystem].Atoms[0][Stack[0][1]].CreationState=DELETED;
                Framework[CurrentSystem].Atoms[0][Stack[0][3]].CreationState=DELETED;
                Framework[CurrentSystem].Atoms[0][Stack[0][5]].CreationState=DELETED;
                Framework[CurrentSystem].Atoms[0][Stack[0][0]].CreationState=NO_DISORDER;
                Framework[CurrentSystem].Atoms[0][Stack[0][2]].CreationState=NO_DISORDER;
                Framework[CurrentSystem].Atoms[0][Stack[0][4]].CreationState=NO_DISORDER;
              }
              break;
            default:
              fprintf(stderr, "SIZE: %d\n",StackSize[0]);
              break;
          }

        }
      }
    }
    count=0;
    count2=0;
    for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[0];i++)
      if(Framework[CurrentSystem].Atoms[0][i].CreationState==DELETED) count++; else count2++;

    fprintf(stderr, "Deleting %d hydrogens, number of atoms: %d\n",count,count2);

    for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[0];i++)
      if(Framework[CurrentSystem].Atoms[0][i].CreationState==DELETED)
      {
        for(j=i;j<Framework[CurrentSystem].NumberOfAtoms[0]-1;j++)
          Framework[CurrentSystem].Atoms[0][j]=Framework[CurrentSystem].Atoms[0][j+1];
        Framework[CurrentSystem].NumberOfAtoms[0]--;
        i--;
      }

    count=0;
    count2=0;
    for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[0];i++)
      if(Framework[CurrentSystem].Atoms[0][i].CreationState==DELETED) count++; else count2++;

    fprintf(stderr, "Deleting %d hydrogens, number of atoms: %d\n",count,count2);
    fprintf(stderr, "Number of atoms: %d\n",Framework[CurrentSystem].NumberOfAtoms[0]);

  }

  //Framework[CurrentSystem].NumberOfFrameworks=1;

  //MakeConnectivityList();
  //PrintConnectivityList();



  char buffer[1024];
  sprintf(buffer,"framework_generated_%d_%d_%d_%d_%ld%s.cssr",(int)SelectionProbability[0],(int)SelectionProbability[1],(int)SelectionProbability[2],(int)SelectionProbability[3],seed,FileNameAppend);

  FilePtr=fopen(buffer,"w");
  fprintf(FilePtr,"%38c%7.4lf %7.4lf% 7.4lf\n",' ',
          (double)Box[CurrentSystem].ax,
          (double)Box[CurrentSystem].by,
          (double)Box[CurrentSystem].cz);
  fprintf(FilePtr,"%21c%8.3lf%8.3lf%8.3lf%4cSPGR =  1 P 1         OPT = 1\n",
          ' ',
          (double)AlphaAngle[CurrentSystem]*RAD2DEG,
          (double)BetaAngle[CurrentSystem]*RAD2DEG,
          (double)GammaAngle[CurrentSystem]*RAD2DEG,
          ' ');
  fprintf(FilePtr,"%4d%4d %s\n",Framework[CurrentSystem].NumberOfAtoms[0],0,"Created by Raspa-1.0");
  fprintf(FilePtr,"     0 %s         : %s\n",Framework[CurrentSystem].Name[0],Framework[CurrentSystem].Name[0]);
  index=0;
  for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[0];i++)
  {
    //if(Framework[CurrentSystem].Atoms[0][i].AssymetricType<=0)
    {
      posA=ConvertFromXYZtoABC(Framework[CurrentSystem].Atoms[0][i].Position);
      fprintf(FilePtr,"%4d %-4s  %9.5lf %9.5lf %9.5lf %4d%4d%4d%4d%4d%4d%4d%4d %7.3lf\n",
              index++,
              PseudoAtoms[Framework[CurrentSystem].Atoms[0][i].Type].Name,
              (double)posA.x,
              (double)posA.y,
              (double)posA.z,
              0,0,0,0,0,0,0,0,
              (double)0.0);
    }
  }
  fclose(FilePtr);

  WriteFrameworkDefinitionCIF("generated");
}

void ReadIonSitingDefinition(void)
{
  int i,temp,index;
  FILE* FilePtr;
  char buffer[256];
  POINT pos;

  sprintf(buffer,"%s.%s",Framework[CurrentSystem].NameIons,"ions");
  if(!(FilePtr=fopen(buffer,"r")))
  {
    sprintf(buffer,"%s/share/raspa/structures/ions/%s.%s",
            RASPA_DIRECTORY,
            Framework[CurrentSystem].NameIons,
            "ions");
    if(!(FilePtr=fopen(buffer,"r")))
    {
      fprintf(stderr, "Error:  file %s does not exist.\n",buffer);
      exit(1);
    }
  }

  fscanf(FilePtr,"%d%*[^\n]",&temp);fscanf(FilePtr,"%*c");
  Framework[CurrentSystem].NumberOfAsymmetricIons=temp;

  crystallographic_stats[CurrentSystem].Position=(VECTOR*)calloc(temp,sizeof(VECTOR));
  crystallographic_stats[CurrentSystem].PositionSquared=(VECTOR*)calloc(temp,sizeof(VECTOR));
  crystallographic_stats[CurrentSystem].Distance=(VECTOR*)calloc(temp,sizeof(VECTOR));
  crystallographic_stats[CurrentSystem].Occupation=(REAL*)calloc(temp,sizeof(REAL));
  crystallographic_stats[CurrentSystem].AverageDistance=(REAL*)calloc(temp,sizeof(REAL));
  crystallographic_stats[CurrentSystem].RelativeOccupation=(REAL*)calloc(temp,sizeof(REAL));
  crystallographic_stats[CurrentSystem].NumberOfCationSites=(int*)calloc(temp,sizeof(int));
  crystallographic_stats[CurrentSystem].Count=(REAL*)calloc(temp,sizeof(REAL));
  crystallographic_stats[CurrentSystem].TemperatureFactor=(REAL_MATRIX3x3*)calloc(temp,sizeof(REAL_MATRIX3x3));

  Framework[CurrentSystem].IonsAsymmetric=(FRAMEWORK_ASYMMETRIC_ATOM*)calloc(temp,sizeof(FRAMEWORK_ASYMMETRIC_ATOM));

  // generate all other atoms from the ones read from the asymmetric-cell
  index=0;
  for(i=0;i<Framework[CurrentSystem].NumberOfAsymmetricIons;i++)
  {
    pos.x=pos.y=pos.z=0.0;
    fscanf(FilePtr,"%lf %lf %lf%*[^\n]",&pos.x,&pos.y,&pos.z);

    Framework[CurrentSystem].IonsAsymmetric[i].Position=pos;
    Framework[CurrentSystem].IonsAsymmetric[i].Type=i;
  }

  ExpandAsymmetricIonsToFullCell();
}

int ClosestCrystallographicPosition(VECTOR pos)
{
  int i,closest;
  REAL r,minimum_distance;
  VECTOR dr;

  closest=0;
  minimum_distance=100.0;
  for(i=0;i<Framework[CurrentSystem].NumberOfIons;i++)
  {
    dr.x=Framework[CurrentSystem].Ions[i].Position.x-pos.x;
    dr.y=Framework[CurrentSystem].Ions[i].Position.y-pos.y;
    dr.z=Framework[CurrentSystem].Ions[i].Position.z-pos.z;
    dr=ApplyBoundaryCondition(dr);
    r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
    if(r<minimum_distance)
    {
      minimum_distance=r;
      closest=i;
    }
  }
  return closest;
}

void ClosestCrystallographicPosition2(VECTOR pos,int *closest,REAL *minimum_distance)
{
  int i;
  REAL r;
  VECTOR dr;

  *closest=0;
  *minimum_distance=100.0;
  for(i=0;i<Framework[CurrentSystem].NumberOfIons;i++)
  {
    dr.x=Framework[CurrentSystem].Ions[i].Position.x-pos.x;
    dr.y=Framework[CurrentSystem].Ions[i].Position.y-pos.y;
    dr.z=Framework[CurrentSystem].Ions[i].Position.z-pos.z;
    dr=ApplyBoundaryCondition(dr);
    r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
    if(r<*minimum_distance)
    {
      *minimum_distance=r;
      *closest=i;
    }
  }
}


void UpdateCrystallographics(void)
{
  int i,ion_type,type;
  int closest,nr_sites,starting_bead;
  VECTOR pos,dr,s;
  REAL distance;

  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    type=Cations[CurrentSystem][i].Type;
    starting_bead=Components[type].StartingBead;
    starting_bead=1;

    // determine position (using the starting bead)
    pos=Cations[CurrentSystem][i].Atoms[starting_bead].Position;

    // look for the closest crystallographic adsorption site
    ClosestCrystallographicPosition2(pos,&closest,&distance);

    // if the distance is less than a certain cutoff consider them 'classified'
    if(distance<CutOffIons)
    {
      ion_type=Framework[CurrentSystem].Ions[closest].AssymetricType;

      // compute distance to closest crystallogarphic site
      dr.x=Framework[CurrentSystem].Ions[closest].Position.x-pos.x;
      dr.y=Framework[CurrentSystem].Ions[closest].Position.y-pos.y;
      dr.z=Framework[CurrentSystem].Ions[closest].Position.z-pos.z;
      dr=ApplyBoundaryCondition(dr);

      s=ConvertToAsymetricUnitCell(pos);

      crystallographic_stats[CurrentSystem].Position[ion_type].x+=s.x;
      crystallographic_stats[CurrentSystem].Position[ion_type].y+=s.y;
      crystallographic_stats[CurrentSystem].Position[ion_type].z+=s.z;

      crystallographic_stats[CurrentSystem].Distance[ion_type].x+=fabs(dr.x);
      crystallographic_stats[CurrentSystem].Distance[ion_type].y+=fabs(dr.y);
      crystallographic_stats[CurrentSystem].Distance[ion_type].z+=fabs(dr.z);
      crystallographic_stats[CurrentSystem].AverageDistance[ion_type]+=distance;

      crystallographic_stats[CurrentSystem].RelativeOccupation[ion_type]+=1.0;
      nr_sites=crystallographic_stats[CurrentSystem].NumberOfCationSites[ion_type];
      crystallographic_stats[CurrentSystem].Occupation[ion_type]+=1.0/nr_sites;
      crystallographic_stats[CurrentSystem].Count[ion_type]+=1.0;

      crystallographic_stats[CurrentSystem].TemperatureFactor[ion_type].ax+=dr.x;
      crystallographic_stats[CurrentSystem].TemperatureFactor[ion_type].ay+=dr.y;
      crystallographic_stats[CurrentSystem].TemperatureFactor[ion_type].az+=dr.z;

      crystallographic_stats[CurrentSystem].TemperatureFactor[ion_type].bx+=SQR(dr.x);
      crystallographic_stats[CurrentSystem].TemperatureFactor[ion_type].by+=SQR(dr.y);
      crystallographic_stats[CurrentSystem].TemperatureFactor[ion_type].bz+=SQR(dr.z);

      crystallographic_stats[CurrentSystem].TemperatureFactor[ion_type].cx+=dr.x*dr.y;
      crystallographic_stats[CurrentSystem].TemperatureFactor[ion_type].cy+=dr.x*dr.z;
      crystallographic_stats[CurrentSystem].TemperatureFactor[ion_type].cz+=dr.y*dr.z;
    }
    else
      crystallographic_stats[CurrentSystem].Unclassified+=1.0;

    crystallographic_stats[CurrentSystem].count+=1.0;
  }
  crystallographic_stats[CurrentSystem].count2+=1.0;
}

void WriteSymmetricFrameworkCssr(void)
{
  int i;
  REAL A,B,C;
  VECTOR pos;
  FILE *FilePtr;

  A=(REAL)NumberOfUnitCells[0].x*UnitCellSize[0].x;
  B=(REAL)NumberOfUnitCells[0].y*UnitCellSize[0].y;
  C=(REAL)NumberOfUnitCells[0].z*UnitCellSize[0].z;

  FilePtr=fopen("framework_asymmetric.cssr","w");
  fprintf(FilePtr,"%38c%8.3lf%8.3lf%8.3lf\n",' ',
     (double)A,
     (double)B,
     (double)C);
  fprintf(FilePtr,"%21c%8.3lf%8.3lf%8.3lf%4cSPGR =  1 P 1         OPT = 1\n",
     ' ',
     (double)AlphaAngle[0]*RAD2DEG,
     (double)BetaAngle[0]*RAD2DEG,
     (double)GammaAngle[0]*RAD2DEG,
     ' ');
  fprintf(FilePtr,"%4d%4d %s\n",Framework[CurrentSystem].NumberOfAsymmetricAtoms[CurrentFramework],0,"Created by Raspa-0.1");
  fprintf(FilePtr,"     0 %s         : %s\n",Framework[CurrentSystem].Name[CurrentFramework],
                Framework[CurrentSystem].Name[CurrentFramework]);
  for(i=0;i<Framework[CurrentSystem].NumberOfAsymmetricAtoms[CurrentFramework];i++)
  {
    pos=Framework[CurrentSystem].AtomsAsymmetric[CurrentFramework][i].Position;
    fprintf(FilePtr,"%4d %-4s  %9.5lf %9.5lf %9.5lf %4d%4d%4d%4d%4d%4d%4d%4d %7.3lf\n",
      i+1,
      PseudoAtoms[Framework[CurrentSystem].Atoms[CurrentFramework][i].Type].Name,
      (double)pos.x,
      (double)pos.y,
      (double)pos.z,
      0,0,0,0,0,0,0,0,
      (double)0.0);
  }
  fclose(FilePtr);
}


int ReturnNumberOfAsymmetricAtoms(int sg)
{
  int i1,i2,i3,k,nr;
  VECTOR pos,s,t;
  int asym;

  nr=0;
  for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];k++)
  {
    pos=Framework[CurrentSystem].Atoms[CurrentFramework][k].Position;

    s=ConvertFromXYZtoABC(pos);
    if(s.x<0.0) s.x+=1.0;
    if(s.y<0.0) s.y+=1.0;
    if(s.z<0.0) s.z+=1.0;

    asym=FALSE;
    for(i1=-2;i1<2;i1++)
    for(i2=-2;i2<2;i2++)
    for(i3=-2;i3<2;i3++)
    {
      t.x=s.x+(REAL)i1;
      t.y=s.y+(REAL)i2;
      t.z=s.z+(REAL)i3;
      if(AsymmetricUnit(sg,s,0.0)) asym=TRUE;
    }
    if(asym) nr++;
  }
  return nr;
}

void CreateAsymetricFrameworkAtoms(void)
{
  int i,j,k;
  int index,Type;
  REAL rr;
  VECTOR pos,dr,s;

  Framework[CurrentSystem].NumberOfAsymmetricAtoms[CurrentFramework]=0;
  for(k=0;k<Framework[CurrentSystem].NumberOfUnitCellAtoms[CurrentFramework];k++)
  {
    Type=Framework[CurrentSystem].Atoms[CurrentFramework][k].Type;
    pos=Framework[CurrentSystem].Atoms[CurrentFramework][k].Position;

    s=ConvertToAsymetricUnitCell(pos);

    SpaceGroupInverseSymmetry(Framework[CurrentSystem].SpaceGroupIdentifier[CurrentFramework],s);
    for(j=0;j<SpaceGroupSize;j++)
    {
      if(AsymmetricUnit(Framework[CurrentSystem].SpaceGroupIdentifier[CurrentFramework],SpaceGroupElement[j],0.0))
      {
        index=-1;
        for(i=0;i<Framework[CurrentSystem].NumberOfAsymmetricAtoms[CurrentFramework];i++)
        {
          dr.x=Framework[CurrentSystem].AtomsAsymmetric[CurrentFramework][i].Position.x-SpaceGroupElement[j].x;
          dr.y=Framework[CurrentSystem].AtomsAsymmetric[CurrentFramework][i].Position.y-SpaceGroupElement[j].y;
          dr.z=Framework[CurrentSystem].AtomsAsymmetric[CurrentFramework][i].Position.z-SpaceGroupElement[j].z;
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
          if((rr<1e-6)&&(Type==Framework[CurrentSystem].AtomsAsymmetric[CurrentFramework][i].Type)) index=i;
        }
        if(index<0)
        {
          Framework[CurrentSystem].AtomsAsymmetric[CurrentFramework][Framework[CurrentSystem].NumberOfAsymmetricAtoms[CurrentFramework]].Position=SpaceGroupElement[j];
          Framework[CurrentSystem].AtomsAsymmetric[CurrentFramework][Framework[CurrentSystem].NumberOfAsymmetricAtoms[CurrentFramework]].Type=Type;
          Framework[CurrentSystem].Atoms[CurrentFramework][k].AssymetricType=Framework[CurrentSystem].NumberOfAsymmetricAtoms[CurrentFramework];
          Framework[CurrentSystem].NumberOfAsymmetricAtoms[CurrentFramework]++;
        }
        else
         Framework[CurrentSystem].Atoms[0][k].AssymetricType=index;
      }
    }
  }
}


void CreateAsymetricFrameworkAtoms2(void)
{
  int i,j,k;
  int index,Type;
  REAL rr;
  VECTOR pos,dr;
  VECTOR s,t;
  REAL tr[27][3]={{-1,-1,-1},{-1,-1,0},{-1,-1,1},{-1,0,-1},{-1,0,0},{-1,0,1},{-1,1,-1},{-1,1,0},{-1,1,1},{0,-1,-1},{0,-1,0},{0,-1,1},{0,0,-1},{0,0,0},{0,0,1},{0,1,-1},{0,1,0},{0,1,1},{1,-1,-1},{1,-1,0},   {1,-1,1},{1,0,-1},{1,0,0},{1,0,1},{1,1,-1},{1,1,0},{1,1,1}};
  int InAsymmetricExists;
  int nr_atoms;


  CurrentFramework=0;
  CurrentSystem=0;
  Framework[CurrentSystem].NumberOfAsymmetricAtoms[CurrentFramework]=0;

  for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];k++)
    Framework[CurrentSystem].Atoms[CurrentFramework][k].ReferencePosition=Framework[CurrentSystem].Atoms[CurrentFramework][k].Position;

  for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];k++)
  {
    Type=Framework[CurrentSystem].Atoms[CurrentFramework][k].Type;
    pos=Framework[CurrentSystem].Atoms[CurrentFramework][k].Position;

    s=ConvertFromXYZtoABC(pos);

    InAsymmetricExists=-1;
    SpaceGroupInverseSymmetry(Framework[CurrentSystem].SpaceGroupIdentifier[CurrentFramework],s);
    for(i=0;i<SpaceGroupSize;i++)
    {
      for(j=0;j<27;j++)
      {
        t.x=SpaceGroupElement[i].x-tr[j][0];
        t.y=SpaceGroupElement[i].y-tr[j][1];
        t.z=SpaceGroupElement[i].z-tr[j][2];
        if(AsymmetricUnit(Framework[CurrentSystem].SpaceGroupIdentifier[CurrentFramework],t,1e-8))
        {
          InAsymmetricExists=i;
          goto l1;
        }
      }
    }
    return;

l1:;
    if(InAsymmetricExists>=0)
    {
      index=-1;
      for(i=0;i<Framework[CurrentSystem].NumberOfAsymmetricAtoms[CurrentFramework];i++)
      {
        dr.x=Framework[CurrentSystem].AtomsAsymmetric[CurrentFramework][i].Position.x-t.x;
        dr.y=Framework[CurrentSystem].AtomsAsymmetric[CurrentFramework][i].Position.y-t.y;
        dr.z=Framework[CurrentSystem].AtomsAsymmetric[CurrentFramework][i].Position.z-t.z;
        dr.x-=NINT(dr.x);
        dr.y-=NINT(dr.y);
        dr.z-=NINT(dr.z);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
        if((rr<1e-4)&&(Type==Framework[CurrentSystem].AtomsAsymmetric[CurrentFramework][i].Type)) index=i;
      }
      if(index<0)
      {
        Framework[CurrentSystem].AtomsAsymmetric[CurrentFramework][Framework[CurrentSystem].NumberOfAsymmetricAtoms[CurrentFramework]].Position=t;
        Framework[CurrentSystem].AtomsAsymmetric[CurrentFramework][Framework[CurrentSystem].NumberOfAsymmetricAtoms[CurrentFramework]].Type=Type;
        Framework[CurrentSystem].Atoms[CurrentFramework][k].AssymetricType=Framework[CurrentSystem].NumberOfAsymmetricAtoms[CurrentFramework];
        Framework[CurrentSystem].NumberOfAsymmetricAtoms[CurrentFramework]++;
      }
      else
        Framework[CurrentSystem].Atoms[0][k].AssymetricType=index;
    }
  }

  // let's generate the structure from the asymmetric atoms

  nr_atoms=0;
  for(j=0;j<SpaceGroupSize;j++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfAsymmetricAtoms[CurrentFramework];i++)
    {
      pos=Framework[CurrentSystem].AtomsAsymmetric[CurrentFramework][i].Position;
      SpaceGroupSymmetry(Framework[CurrentSystem].SpaceGroupIdentifier[CurrentFramework],pos);

      Type=Framework[CurrentSystem].AtomsAsymmetric[CurrentFramework][i].Type;

      // apply boundary condition
      SpaceGroupElement[j].x-=NINT(SpaceGroupElement[j].x);
      SpaceGroupElement[j].y-=NINT(SpaceGroupElement[j].y);
      SpaceGroupElement[j].z-=NINT(SpaceGroupElement[j].z);

      if(SpaceGroupElement[j].x<0.0) SpaceGroupElement[j].x+=1.0;
      if(SpaceGroupElement[j].y<0.0) SpaceGroupElement[j].y+=1.0;
      if(SpaceGroupElement[j].z<0.0) SpaceGroupElement[j].z+=1.0;


      index=-1;
      for(k=0;k<nr_atoms;k++)
      {
        dr.x=Framework[CurrentSystem].Atoms[CurrentFramework][k].Position.x-SpaceGroupElement[j].x;
        dr.y=Framework[CurrentSystem].Atoms[CurrentFramework][k].Position.y-SpaceGroupElement[j].y;
        dr.z=Framework[CurrentSystem].Atoms[CurrentFramework][k].Position.z-SpaceGroupElement[j].z;
        dr.x-=NINT(dr.x);
        dr.y-=NINT(dr.y);
        dr.z-=NINT(dr.z);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
        if(rr<1e-5) index=k;
      }
      if(index<0)
      {
        Framework[CurrentSystem].Atoms[CurrentFramework][nr_atoms].Position=SpaceGroupElement[j];
        Framework[CurrentSystem].Atoms[CurrentFramework][nr_atoms].Type=Type;

        nr_atoms++;
      }
    }
  }
  fprintf(stderr, "nr of atoms: %d\n",nr_atoms);
  //if(nr_atoms!=Framework[CurrentSystem].NumberOfAtoms[CurrentFramework])
  //  Framework[CurrentSystem].NumberOfAsymmetricAtoms[CurrentFramework]=0;
}


void AllocateConnectivityList(void)
{
  int k,l,f;
  for(f=0;f<Framework[CurrentSystem].NumberOfFrameworks;f++)
  {
    Framework[CurrentSystem].Connectivity[f]=(int*)calloc(
         TotalNumberOfReplicaCells[CurrentSystem]*Framework[CurrentSystem].NumberOfAtoms[f],sizeof(int));

    Framework[CurrentSystem].Neighbours[f]=(int**)calloc(
       TotalNumberOfReplicaCells[CurrentSystem]*Framework[CurrentSystem].NumberOfAtoms[f],sizeof(int*));
    for(k=0;k<TotalNumberOfReplicaCells[CurrentSystem]*Framework[CurrentSystem].NumberOfAtoms[f];k++)
    {
      Framework[CurrentSystem].Neighbours[f][k]=(int*)calloc(64,sizeof(int));
      for(l=0;l<60;l++)
        Framework[CurrentSystem].Neighbours[f][k][l]=-1;
    }
  }
}

void FreeAllocateConnectivityList(void)
{
  int k,f;

  for(f=0;f<Framework[CurrentSystem].NumberOfFrameworks;f++)
  {

    for(k=0;k<TotalNumberOfReplicaCells[CurrentSystem]*(Framework[CurrentSystem].NumberOfAtoms[f]-Framework[CurrentSystem].NumberOfCoreShells[f]);k++)
      free(Framework[CurrentSystem].Neighbours[f][k]);

    free(Framework[CurrentSystem].Neighbours[f]);
    free(Framework[CurrentSystem].Connectivity[f]);
  }
}

/*********************************************************************************************************
 * Name       | GetNeighbour                                                                             *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Returns the neighbor in the unit-cell                                                    *
 * Note       | The Neighbours-list is made for the replica-cell and of size                             *
 *            | 'TotalNumberOfReplicaCells[CurrentSystem]*Framework[CurrentSystem].NumberOfAtoms[f]'     *
 *            | The first 'Framework[CurrentSystem].NumberOfAtoms[f]'-elements are the unit cell ones.   *
 *            | The order is cell by cell and the atoms are ordered identically.                         *
 *            | Therefore to get the neighbor in the unit cell the modules-operator is used.             *
 *********************************************************************************************************/

int GetNeighbour(int system,int f1,int A,int k)
{
  int B;

  A%=Framework[system].NumberOfAtoms[f1];
  B=Framework[system].Neighbours[f1][A][k]%Framework[system].NumberOfAtoms[f1];
  return B;
}

/*********************************************************************************************************
 * Name       | GetReplicaNeighbour                                                                      *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Returns the neighbor in the replica-cell                                                 *
 * Note       | The Neighbours-list is made for the replica-cell and of size                             *
 *            | 'TotalNumberOfReplicaCells[CurrentSystem]*Framework[CurrentSystem].NumberOfAtoms[f]'     *
 *            | The first 'Framework[CurrentSystem].NumberOfAtoms[f]'-elements are the unit cell ones.   *
 *********************************************************************************************************/

int GetReplicaNeighbour(int system,int f1,int A,int k)
{
  int B;

  B=Framework[system].Neighbours[f1][A][k];
  return B;
}


/*********************************************************************************************************
 * Name       | MakeConnectivityList                                                                     *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Creates a connectivity-list of the replica-cell                                          *
 * Note       | The Neighbours-list is of size                                                           *
 *            | 'TotalNumberOfReplicaCells[CurrentSystem]*Framework[CurrentSystem].NumberOfAtoms[f]'     *
 *            | The first 'Framework[CurrentSystem].NumberOfAtoms[f]'-elements are the unit cell ones.   *
 *            | The order is cell by cell and the atoms are ordered identically.                         *
 *            | Two atoms are considered bonded when the distance between them is smaller than the sum   *
 *            | of their covalent radii plus 0.56 Angstroms.                                             *
 *********************************************************************************************************/

void MakeConnectivityList(void)
{
  int i,j,k,l,f;
  int typeA,typeB,ncell1,ncell2;
  REAL r,Bond_length;
  VECTOR dr,posA,posB;
  int index_i,index_j;

  for(f=0;f<Framework[CurrentSystem].NumberOfFrameworks;f++)
  {
    for(k=0;k<TotalNumberOfReplicaCells[CurrentSystem]*Framework[CurrentSystem].NumberOfAtoms[f];k++)
    {
      for(l=0;l<60;l++)
        Framework[CurrentSystem].Neighbours[f][k][l]=-1;
    }
  }

  for(f=0;f<Framework[CurrentSystem].NumberOfFrameworks;f++)
  {
    for(ncell1=0;ncell1<TotalNumberOfReplicaCells[CurrentSystem];ncell1++)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f];i++)
      {
        index_i=i+ncell1*Framework[CurrentSystem].NumberOfAtoms[f];
        Framework[CurrentSystem].Connectivity[f][index_i]=0;

        posA=Framework[CurrentSystem].Atoms[f][i].Position;
        typeA=Framework[CurrentSystem].Atoms[f][i].Type;
        posA.x+=ReplicaShift[ncell1].x;
        posA.y+=ReplicaShift[ncell1].y;
        posA.z+=ReplicaShift[ncell1].z;

        for(ncell2=0;ncell2<TotalNumberOfReplicaCells[CurrentSystem];ncell2++)
        {
          for(j=0;j<Framework[CurrentSystem].NumberOfAtoms[f];j++)
          {
            index_j=j+ncell2*Framework[CurrentSystem].NumberOfAtoms[f];

            posB=Framework[CurrentSystem].Atoms[f][j].Position;
            typeB=Framework[CurrentSystem].Atoms[f][j].Type;
            posB.x+=ReplicaShift[ncell2].x;
            posB.y+=ReplicaShift[ncell2].y;
            posB.z+=ReplicaShift[ncell2].z;

            Bond_length=0.56+PseudoAtoms[typeA].Radius+PseudoAtoms[typeB].Radius;
            dr.x=posA.x-posB.x;
            dr.y=posA.y-posB.y;
            dr.z=posA.z-posB.z;
            dr=ApplyReplicaBoundaryCondition(dr);
            r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));

            if((index_i!=index_j)&&(r<Bond_length))
              Framework[CurrentSystem].Neighbours[f][index_i][Framework[CurrentSystem].Connectivity[f][index_i]++]=index_j;
          }
        }
      }
    }
  }
}


void PrintConnectivityList(void)
{
  int i;
  int f,ncell,index_i;

  if(UseReplicas[CurrentSystem])
  {
    for(f=0;f<Framework[CurrentSystem].NumberOfFrameworks;f++)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f];i++)
      {
        for(ncell=0;ncell<TotalNumberOfReplicaCells[CurrentSystem];ncell++)
        {
          index_i=i+ncell*Framework[CurrentSystem].NumberOfAtoms[f];

          fprintf(stderr, "%d connectivity: %d type: %s %d %d %d %d %d %d %d %d %lf %lf %lf\n",
            index_i,
            Framework[CurrentSystem].Connectivity[f][index_i],
            PseudoAtoms[Framework[CurrentSystem].Atoms[f][i].Type].Name,
            Framework[CurrentSystem].Neighbours[f][index_i][0],
            Framework[CurrentSystem].Neighbours[f][index_i][1],
            Framework[CurrentSystem].Neighbours[f][index_i][2],
            Framework[CurrentSystem].Neighbours[f][index_i][3],
            Framework[CurrentSystem].Neighbours[f][index_i][4],
            Framework[CurrentSystem].Neighbours[f][index_i][5],
            Framework[CurrentSystem].Neighbours[f][index_i][6],
            Framework[CurrentSystem].Neighbours[f][index_i][7],
            (double)Framework[CurrentSystem].Atoms[f][i].Position.x+ReplicaShift[ncell].x,
            (double)Framework[CurrentSystem].Atoms[f][i].Position.y+ReplicaShift[ncell].y,
            (double)Framework[CurrentSystem].Atoms[f][i].Position.z+ReplicaShift[ncell].z);
        }
      }
    }
  }
  else
  {
    for(CurrentFramework=0;CurrentFramework<Framework[CurrentSystem].NumberOfFrameworks;CurrentFramework++)
    {
      for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];i++)
      {
        fprintf(stderr, "%d connectivity: %d type: %s %d %d %d %d %d %d %d %d %lf %lf %lf\n",
           i,
           Framework[CurrentSystem].Connectivity[CurrentFramework][i],
           PseudoAtoms[Framework[CurrentSystem].Atoms[CurrentFramework][i].Type].Name,
           Framework[CurrentSystem].Neighbours[CurrentFramework][i][0],
           Framework[CurrentSystem].Neighbours[CurrentFramework][i][1],
           Framework[CurrentSystem].Neighbours[CurrentFramework][i][2],
           Framework[CurrentSystem].Neighbours[CurrentFramework][i][3],
           Framework[CurrentSystem].Neighbours[CurrentFramework][i][4],
           Framework[CurrentSystem].Neighbours[CurrentFramework][i][5],
           Framework[CurrentSystem].Neighbours[CurrentFramework][i][6],
           Framework[CurrentSystem].Neighbours[CurrentFramework][i][7],
           (double)Framework[CurrentSystem].Atoms[CurrentFramework][i].Position.x,
           (double)Framework[CurrentSystem].Atoms[CurrentFramework][i].Position.y,
           (double)Framework[CurrentSystem].Atoms[CurrentFramework][i].Position.z);
      }
    }
  }
}

void SubstituteAtoms(void)
{
  int i,j,k,nr,f1;
  int TypeSource,TypeDest;
  int index,index2,count;

  count=0;
  NumberOfRandomSubstitutions=0;
  for(i=0;i<NumberOfSubstitutionRules;i++)
    NumberOfRandomSubstitutions+=atoi(SubstitutionFrameworkAtoms[i][0]);
  NumberOfSubstitutions=NumberOfSingleSubstitutionRules+NumberOfRandomSubstitutions;

  ListOfAtomSubstitutions=(int(*)[3])calloc(NumberOfSubstitutions,sizeof(int[3]));

  for(i=0;i<NumberOfSingleSubstitutionRules;i++)
  {
    nr=atoi(SubstitutionSingleFrameworkAtom[i][0]);
    TypeSource=ReturnPseudoAtomNumber(SubstitutionSingleFrameworkAtom[i][1]);
    TypeDest=ReturnPseudoAtomNumber(SubstitutionSingleFrameworkAtom[i][2]);

    index=0;
    index2=0;
    for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
    {
      for(j=0;j<Framework[CurrentSystem].NumberOfAtoms[f1];j++)
      {
        if(Framework[CurrentSystem].Atoms[f1][j].Type==TypeSource)
        {
          if(index2==nr)
          {
            Framework[CurrentSystem].Atoms[f1][j].Type=TypeDest;
            Framework[CurrentSystem].Atoms[f1][j].Modified=TRUE;
            Framework[CurrentSystem].Atoms[f1][j].OriginalType=TypeSource;
            Framework[CurrentSystem].Atoms[f1][j].Charge=PseudoAtoms[TypeDest].Charge1;
            NumberOfPseudoAtomsType[CurrentSystem][TypeDest]++;
            NumberOfPseudoAtomsType[CurrentSystem][TypeSource]--;

            ListOfAtomSubstitutions[count][0]=index2;
            ListOfAtomSubstitutions[count][1]=TypeSource;
            ListOfAtomSubstitutions[count][2]=TypeDest;
            count++;
            goto done1;
          }
          index++;
          index2++;
        }
        if(Framework[CurrentSystem].Atoms[f1][j].Modified&&(Framework[CurrentSystem].Atoms[f1][j].OriginalType==TypeSource))
          index2++;
      }
    }
    done1:;
  }

  // the remaining atoms can be randomly substituted
  for(i=0;i<NumberOfSubstitutionRules;i++)
  {
    nr=atoi(SubstitutionFrameworkAtoms[i][0]);
    TypeSource=ReturnPseudoAtomNumber(SubstitutionFrameworkAtoms[i][1]);
    TypeDest=ReturnPseudoAtomNumber(SubstitutionFrameworkAtoms[i][2]);

    if(nr>NumberOfPseudoAtomsType[CurrentSystem][TypeSource])
      nr=NumberOfPseudoAtomsType[CurrentSystem][TypeSource];

    do
    {
      // substitute the 'k-th' atom of type 'source' by type 'dest'
      k=RandomNumber()*NumberOfPseudoAtomsType[CurrentSystem][TypeSource];;

      index=0;
      index2=0;
      for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
      {
        for(j=0;j<Framework[CurrentSystem].NumberOfAtoms[f1];j++)
        {
          if(Framework[CurrentSystem].Atoms[f1][j].Type==TypeSource)
          {
            if(index==k)
            {
              Framework[CurrentSystem].Atoms[f1][j].Type=TypeDest;
              Framework[CurrentSystem].Atoms[f1][j].Modified=TRUE;
              Framework[CurrentSystem].Atoms[f1][j].OriginalType=TypeSource;
              Framework[CurrentSystem].Atoms[f1][j].Charge=PseudoAtoms[TypeDest].Charge1;
              NumberOfPseudoAtomsType[CurrentSystem][TypeDest]++;
              NumberOfPseudoAtomsType[CurrentSystem][TypeSource]--;
              ListOfAtomSubstitutions[count][0]=index2;
              ListOfAtomSubstitutions[count][1]=TypeSource;
              ListOfAtomSubstitutions[count][2]=TypeDest;
              count++;
              goto done;
            }
            index++;
            index2++;
          }
          if(Framework[CurrentSystem].Atoms[f1][j].Modified&&(Framework[CurrentSystem].Atoms[f1][j].OriginalType==TypeSource))
            index2++;
        }
      }
      done: nr--;
    }
    while(nr>0);
  }
}


void ModifyAtomsConnectedToDefinedNeighbours(void)
{
  int i,j,k,l,f1;
  int A,B,C,D,E,F;
  int found[12];
  int nr_found;
  REAL temp_real;
  REAL Angle1,Angle2;

  Lowenstein[CurrentSystem]=TRUE;

  for(l=0;l<NumberOfModificationRules;l++)
  {
    switch(ModificationRuleType[l])
    {
      case MODIFY_FRAMEWORKATOM_TRIPLE:
        for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
        {
          for(A=0;A<Framework[CurrentSystem].NumberOfAtoms[f1];A++)
          {
            for(j=0;j<Framework[CurrentSystem].Connectivity[f1][A];j++)
            {
              B=GetNeighbour(CurrentSystem,f1,A,j);
              for(k=0;k<Framework[CurrentSystem].Connectivity[f1][B];k++)
              {
                C=GetNeighbour(CurrentSystem,f1,B,k);

                if((A!=C)&&
                   (Framework[CurrentSystem].Atoms[f1][A].Type==ModifyFrameworkAtomTypes[l][0])&&
                   (Framework[CurrentSystem].Atoms[f1][B].Type==ModifyFrameworkAtomTypes[l][1])&&
                   (Framework[CurrentSystem].Atoms[f1][C].Type==ModifyFrameworkAtomTypes[l][2]))
                {
                  Framework[CurrentSystem].Atoms[f1][A].Type=ModifyFrameworkAtomTypes[l][3];
                  if(!UseChargesFromCIFFile)
                    Framework[CurrentSystem].Atoms[f1][A].Charge=PseudoAtoms[ModifyFrameworkAtomTypes[l][3]].Charge1;
                  NumberOfPseudoAtomsCount[CurrentSystem][ModifyFrameworkAtomTypes[l][3]]++;
                  NumberOfPseudoAtomsType[CurrentSystem][ModifyFrameworkAtomTypes[l][3]]++;
                  NumberOfPseudoAtomsCount[CurrentSystem][ModifyFrameworkAtomTypes[l][0]]--;
                  NumberOfPseudoAtomsType[CurrentSystem][ModifyFrameworkAtomTypes[l][0]]--;

                  Framework[CurrentSystem].Atoms[f1][B].Type=ModifyFrameworkAtomTypes[l][4];
                  if(!UseChargesFromCIFFile)
                    Framework[CurrentSystem].Atoms[f1][B].Charge=PseudoAtoms[ModifyFrameworkAtomTypes[l][4]].Charge1;
                  NumberOfPseudoAtomsCount[CurrentSystem][ModifyFrameworkAtomTypes[l][4]]++;
                  NumberOfPseudoAtomsType[CurrentSystem][ModifyFrameworkAtomTypes[l][4]]++;
                  NumberOfPseudoAtomsCount[CurrentSystem][ModifyFrameworkAtomTypes[l][1]]--;
                  NumberOfPseudoAtomsType[CurrentSystem][ModifyFrameworkAtomTypes[l][1]]--;

                  Framework[CurrentSystem].Atoms[f1][C].Type=ModifyFrameworkAtomTypes[l][5];
                  if(!UseChargesFromCIFFile)
                    Framework[CurrentSystem].Atoms[f1][C].Charge=PseudoAtoms[ModifyFrameworkAtomTypes[l][5]].Charge1;
                  NumberOfPseudoAtomsCount[CurrentSystem][ModifyFrameworkAtomTypes[l][5]]++;
                  NumberOfPseudoAtomsType[CurrentSystem][ModifyFrameworkAtomTypes[l][5]]++;
                  NumberOfPseudoAtomsCount[CurrentSystem][ModifyFrameworkAtomTypes[l][2]]--;
                  NumberOfPseudoAtomsType[CurrentSystem][ModifyFrameworkAtomTypes[l][2]]--;
                }
              }
            }
          }
        }
        break;
      case MODIFY_FRAMEWORKATOM_PLANAR:
        for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
        {
          for(A=0;A<Framework[CurrentSystem].NumberOfAtoms[f1];A++)
          {
            if((Framework[CurrentSystem].Connectivity[f1][A]>=4)&&(Framework[CurrentSystem].Atoms[f1][A].Type==ModifyFrameworkAtomTypes[l][0]))
            {
              nr_found=0;
              for(j=0;j<Framework[CurrentSystem].Connectivity[f1][A];j++)
              {
                for(k=nr_found;k<4;k++)
                {
                  found[k]=-1;
                  B=GetNeighbour(CurrentSystem,f1,A,j);
                  if(Framework[CurrentSystem].Atoms[f1][B].Type==ModifyFrameworkAtomTypes[l][k+1])
                  {
                    found[k]=B;
                    nr_found++;
                    break;
                  }
                }
              }
              fprintf(stderr, "found: %d\n",nr_found);

              Angle1=ReturnConstraintBendAngle(Framework[CurrentSystem].Atoms[f1][found[0]].Position,Framework[CurrentSystem].Atoms[f1][A].Position,
                     Framework[CurrentSystem].Atoms[f1][found[1]].Position);
              Angle2=ReturnConstraintBendAngle(Framework[CurrentSystem].Atoms[f1][found[0]].Position,Framework[CurrentSystem].Atoms[f1][A].Position,
                     Framework[CurrentSystem].Atoms[f1][found[2]].Position);
              if(Angle1>Angle2)
                SWAP(found[1],found[2],temp_real);


              Angle1=ReturnConstraintBendAngle(Framework[CurrentSystem].Atoms[f1][found[0]].Position,Framework[CurrentSystem].Atoms[f1][A].Position,
                     Framework[CurrentSystem].Atoms[f1][found[3]].Position);
              Angle2=ReturnConstraintBendAngle(Framework[CurrentSystem].Atoms[f1][found[0]].Position,Framework[CurrentSystem].Atoms[f1][A].Position,
                     Framework[CurrentSystem].Atoms[f1][found[2]].Position);
              if(Angle1>Angle2)
                SWAP(found[2],found[3],temp_real);

              Framework[CurrentSystem].Atoms[f1][A].Type=ModifyFrameworkAtomTypes[l][5];
              if(!UseChargesFromCIFFile)
                Framework[CurrentSystem].Atoms[f1][A].Charge=PseudoAtoms[ModifyFrameworkAtomTypes[l][5]].Charge1;
              NumberOfPseudoAtomsCount[CurrentSystem][ModifyFrameworkAtomTypes[l][0]]--;
              NumberOfPseudoAtomsType[CurrentSystem][ModifyFrameworkAtomTypes[l][0]]--;
              NumberOfPseudoAtomsCount[CurrentSystem][ModifyFrameworkAtomTypes[l][5]]++;
              NumberOfPseudoAtomsType[CurrentSystem][ModifyFrameworkAtomTypes[l][5]]++;

              Framework[CurrentSystem].Atoms[f1][found[0]].Type=ModifyFrameworkAtomTypes[l][6];
              if(!UseChargesFromCIFFile)
                Framework[CurrentSystem].Atoms[f1][found[0]].Charge=PseudoAtoms[ModifyFrameworkAtomTypes[l][6]].Charge1;
              NumberOfPseudoAtomsCount[CurrentSystem][ModifyFrameworkAtomTypes[l][1]]--;
              NumberOfPseudoAtomsType[CurrentSystem][ModifyFrameworkAtomTypes[l][1]]--;
              NumberOfPseudoAtomsCount[CurrentSystem][ModifyFrameworkAtomTypes[l][6]]++;
              NumberOfPseudoAtomsType[CurrentSystem][ModifyFrameworkAtomTypes[l][6]]++;

              Framework[CurrentSystem].Atoms[f1][found[1]].Type=ModifyFrameworkAtomTypes[l][7];
              if(!UseChargesFromCIFFile)
                Framework[CurrentSystem].Atoms[f1][found[1]].Charge=PseudoAtoms[ModifyFrameworkAtomTypes[l][7]].Charge1;
              NumberOfPseudoAtomsCount[CurrentSystem][ModifyFrameworkAtomTypes[l][2]]--;
              NumberOfPseudoAtomsType[CurrentSystem][ModifyFrameworkAtomTypes[l][2]]--;
              NumberOfPseudoAtomsCount[CurrentSystem][ModifyFrameworkAtomTypes[l][7]]++;
              NumberOfPseudoAtomsType[CurrentSystem][ModifyFrameworkAtomTypes[l][7]]++;

              Framework[CurrentSystem].Atoms[f1][found[2]].Type=ModifyFrameworkAtomTypes[l][8];
              if(!UseChargesFromCIFFile)
                Framework[CurrentSystem].Atoms[f1][found[2]].Charge=PseudoAtoms[ModifyFrameworkAtomTypes[l][8]].Charge1;
              NumberOfPseudoAtomsCount[CurrentSystem][ModifyFrameworkAtomTypes[l][3]]--;
              NumberOfPseudoAtomsType[CurrentSystem][ModifyFrameworkAtomTypes[l][3]]--;
              NumberOfPseudoAtomsCount[CurrentSystem][ModifyFrameworkAtomTypes[l][8]]++;
              NumberOfPseudoAtomsType[CurrentSystem][ModifyFrameworkAtomTypes[l][8]]++;

              Framework[CurrentSystem].Atoms[f1][found[3]].Type=ModifyFrameworkAtomTypes[l][9];
              if(!UseChargesFromCIFFile)
                Framework[CurrentSystem].Atoms[f1][found[3]].Charge=PseudoAtoms[ModifyFrameworkAtomTypes[l][9]].Charge1;
              NumberOfPseudoAtomsCount[CurrentSystem][ModifyFrameworkAtomTypes[l][4]]--;
              NumberOfPseudoAtomsType[CurrentSystem][ModifyFrameworkAtomTypes[l][4]]--;
              NumberOfPseudoAtomsCount[CurrentSystem][ModifyFrameworkAtomTypes[l][9]]++;
              NumberOfPseudoAtomsType[CurrentSystem][ModifyFrameworkAtomTypes[l][9]]++;
            }
          }
        }
        break;
      case MODIFY_FRAMEWORKATOM_DIMER:
        for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
        {
          for(A=0;A<Framework[CurrentSystem].NumberOfAtoms[f1];A++)
          {
            for(j=0;j<Framework[CurrentSystem].Connectivity[f1][A];j++)
            {
              B=GetNeighbour(CurrentSystem,f1,A,j);

              if((Framework[CurrentSystem].Atoms[f1][A].Type==ModifyFrameworkAtomTypes[l][0])&&
                 (Framework[CurrentSystem].Atoms[f1][B].Type==ModifyFrameworkAtomTypes[l][1]))
              {
                Framework[CurrentSystem].Atoms[f1][A].Type=ModifyFrameworkAtomTypes[l][2];
                if(!UseChargesFromCIFFile)
                  Framework[CurrentSystem].Atoms[f1][A].Charge=PseudoAtoms[ModifyFrameworkAtomTypes[l][2]].Charge1;
                NumberOfPseudoAtomsCount[CurrentSystem][ModifyFrameworkAtomTypes[l][2]]++;
                NumberOfPseudoAtomsType[CurrentSystem][ModifyFrameworkAtomTypes[l][2]]++;
                NumberOfPseudoAtomsCount[CurrentSystem][ModifyFrameworkAtomTypes[l][0]]--;
                NumberOfPseudoAtomsType[CurrentSystem][ModifyFrameworkAtomTypes[l][0]]--;

                Framework[CurrentSystem].Atoms[f1][B].Type=ModifyFrameworkAtomTypes[l][3];
                if(!UseChargesFromCIFFile)
                  Framework[CurrentSystem].Atoms[f1][B].Charge=PseudoAtoms[ModifyFrameworkAtomTypes[l][3]].Charge1;
                NumberOfPseudoAtomsCount[CurrentSystem][ModifyFrameworkAtomTypes[l][3]]++;
                NumberOfPseudoAtomsType[CurrentSystem][ModifyFrameworkAtomTypes[l][3]]++;
                NumberOfPseudoAtomsCount[CurrentSystem][ModifyFrameworkAtomTypes[l][1]]--;
                NumberOfPseudoAtomsType[CurrentSystem][ModifyFrameworkAtomTypes[l][1]]--;
              }
            }
          }
        }
        break;
      case MODIFY_FRAMEWORKATOM_CONNECTED_TO:
        // case that there are no connected atoms defined
        if(ModifyFrameworkAtomTypes[l][2]<0)
        {
          for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
          {
            for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
            {
              for(j=0;j<Framework[CurrentSystem].Connectivity[f1][i];j++)
              {
                if(Framework[CurrentSystem].Atoms[f1][i].Type==ModifyFrameworkAtomTypes[l][0])
                {
                  Framework[CurrentSystem].Atoms[f1][i].Type=ModifyFrameworkAtomTypes[l][1];
                  if(!UseChargesFromCIFFile)
                    Framework[CurrentSystem].Atoms[f1][i].Charge=PseudoAtoms[ModifyFrameworkAtomTypes[l][1]].Charge1;

                  NumberOfPseudoAtomsCount[CurrentSystem][ModifyFrameworkAtomTypes[l][1]]++;
                  NumberOfPseudoAtomsType[CurrentSystem][ModifyFrameworkAtomTypes[l][1]]++;
                  NumberOfPseudoAtomsCount[CurrentSystem][ModifyFrameworkAtomTypes[l][0]]--;
                  NumberOfPseudoAtomsType[CurrentSystem][ModifyFrameworkAtomTypes[l][0]]--;
                }
              }
            }
          }
        }
        // case that there is only one connected atom defined
        else if(ModifyFrameworkAtomTypes[l][3]<0)
        {
          for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
          {
            for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
            {
              for(j=0;j<Framework[CurrentSystem].Connectivity[f1][i];j++)
              {
                if((Framework[CurrentSystem].Atoms[f1][i].Type==ModifyFrameworkAtomTypes[l][0])&&
                   (Framework[CurrentSystem].Atoms[f1][GetNeighbour(CurrentSystem,f1,i,j)].Type==ModifyFrameworkAtomTypes[l][2]))
                {
                  Framework[CurrentSystem].Atoms[f1][i].Type=ModifyFrameworkAtomTypes[l][1];
                  if(!UseChargesFromCIFFile)
                    Framework[CurrentSystem].Atoms[f1][i].Charge=PseudoAtoms[ModifyFrameworkAtomTypes[l][1]].Charge1;

                  NumberOfPseudoAtomsCount[CurrentSystem][ModifyFrameworkAtomTypes[l][1]]++;
                  NumberOfPseudoAtomsType[CurrentSystem][ModifyFrameworkAtomTypes[l][1]]++;
                  NumberOfPseudoAtomsCount[CurrentSystem][ModifyFrameworkAtomTypes[l][0]]--;
                  NumberOfPseudoAtomsType[CurrentSystem][ModifyFrameworkAtomTypes[l][0]]--;
                }
              }
            }
          }
        }
        else
        {
          for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
          {
            for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
            {
              switch(Framework[CurrentSystem].Connectivity[f1][i])
              {
                case 0:
                  break;
                case 1:
                  break;
                case 2:
                  A=GetNeighbour(CurrentSystem,f1,i,0);
                  B=GetNeighbour(CurrentSystem,f1,i,1);
                  if((Framework[CurrentSystem].Atoms[f1][i].Type==ModifyFrameworkAtomTypes[l][0])&&
                     (((Framework[CurrentSystem].Atoms[f1][A].Type==ModifyFrameworkAtomTypes[l][2])&&
                      (Framework[CurrentSystem].Atoms[f1][B].Type==ModifyFrameworkAtomTypes[l][3]))||
                     ((Framework[CurrentSystem].Atoms[f1][B].Type==ModifyFrameworkAtomTypes[l][2])&&
                      (Framework[CurrentSystem].Atoms[f1][A].Type==ModifyFrameworkAtomTypes[l][3]))))
                  {
                    Framework[CurrentSystem].Atoms[f1][i].Type=ModifyFrameworkAtomTypes[l][1];
                    if(!UseChargesFromCIFFile)
                      Framework[CurrentSystem].Atoms[f1][i].Charge=PseudoAtoms[ModifyFrameworkAtomTypes[l][1]].Charge1;

                    NumberOfPseudoAtomsCount[CurrentSystem][ModifyFrameworkAtomTypes[l][1]]++;
                    NumberOfPseudoAtomsType[CurrentSystem][ModifyFrameworkAtomTypes[l][1]]++;
                    NumberOfPseudoAtomsCount[CurrentSystem][ModifyFrameworkAtomTypes[l][0]]--;
                    NumberOfPseudoAtomsType[CurrentSystem][ModifyFrameworkAtomTypes[l][0]]--;
                  }
                  break;
                case 3:
                  A=GetNeighbour(CurrentSystem,f1,i,0);
                  B=GetNeighbour(CurrentSystem,f1,i,1);
                  C=GetNeighbour(CurrentSystem,f1,i,2);
                  if((Framework[CurrentSystem].Atoms[f1][i].Type==ModifyFrameworkAtomTypes[l][0])&&
                    ((((Framework[CurrentSystem].Atoms[f1][A].Type==ModifyFrameworkAtomTypes[l][2])&&
                       (Framework[CurrentSystem].Atoms[f1][B].Type==ModifyFrameworkAtomTypes[l][3]))||
                      ((Framework[CurrentSystem].Atoms[f1][B].Type==ModifyFrameworkAtomTypes[l][2])&&
                       (Framework[CurrentSystem].Atoms[f1][A].Type==ModifyFrameworkAtomTypes[l][3])))||
                     (((Framework[CurrentSystem].Atoms[f1][A].Type==ModifyFrameworkAtomTypes[l][2])&&
                       (Framework[CurrentSystem].Atoms[f1][C].Type==ModifyFrameworkAtomTypes[l][3]))||
                      ((Framework[CurrentSystem].Atoms[f1][C].Type==ModifyFrameworkAtomTypes[l][2])&&
                       (Framework[CurrentSystem].Atoms[f1][A].Type==ModifyFrameworkAtomTypes[l][3])))||
                     (((Framework[CurrentSystem].Atoms[f1][B].Type==ModifyFrameworkAtomTypes[l][2])&&
                       (Framework[CurrentSystem].Atoms[f1][C].Type==ModifyFrameworkAtomTypes[l][3]))||
                      ((Framework[CurrentSystem].Atoms[f1][C].Type==ModifyFrameworkAtomTypes[l][2])&&
                       (Framework[CurrentSystem].Atoms[f1][B].Type==ModifyFrameworkAtomTypes[l][3])))))
                  {
                    Framework[CurrentSystem].Atoms[f1][i].Type=ModifyFrameworkAtomTypes[l][1];
                    if(!UseChargesFromCIFFile)
                      Framework[CurrentSystem].Atoms[f1][i].Charge=PseudoAtoms[ModifyFrameworkAtomTypes[l][1]].Charge1;

                    NumberOfPseudoAtomsCount[CurrentSystem][ModifyFrameworkAtomTypes[l][1]]++;
                    NumberOfPseudoAtomsType[CurrentSystem][ModifyFrameworkAtomTypes[l][1]]++;
                    NumberOfPseudoAtomsCount[CurrentSystem][ModifyFrameworkAtomTypes[l][0]]--;
                    NumberOfPseudoAtomsType[CurrentSystem][ModifyFrameworkAtomTypes[l][0]]--;
                  }
                  break;
                case 4:
                  A=GetNeighbour(CurrentSystem,f1,i,0);
                  B=GetNeighbour(CurrentSystem,f1,i,1);
                  C=GetNeighbour(CurrentSystem,f1,i,2);
                  D=GetNeighbour(CurrentSystem,f1,i,3);
                  if((Framework[CurrentSystem].Atoms[f1][i].Type==ModifyFrameworkAtomTypes[l][0])&&
                    ((((Framework[CurrentSystem].Atoms[f1][A].Type==ModifyFrameworkAtomTypes[l][2])&&
                       (Framework[CurrentSystem].Atoms[f1][B].Type==ModifyFrameworkAtomTypes[l][3]))||
                      ((Framework[CurrentSystem].Atoms[f1][B].Type==ModifyFrameworkAtomTypes[l][2])&&
                       (Framework[CurrentSystem].Atoms[f1][A].Type==ModifyFrameworkAtomTypes[l][3])))||
                     (((Framework[CurrentSystem].Atoms[f1][A].Type==ModifyFrameworkAtomTypes[l][2])&&
                       (Framework[CurrentSystem].Atoms[f1][C].Type==ModifyFrameworkAtomTypes[l][3]))||
                      ((Framework[CurrentSystem].Atoms[f1][C].Type==ModifyFrameworkAtomTypes[l][2])&&
                       (Framework[CurrentSystem].Atoms[f1][A].Type==ModifyFrameworkAtomTypes[l][3])))||
                     (((Framework[CurrentSystem].Atoms[f1][A].Type==ModifyFrameworkAtomTypes[l][2])&&
                       (Framework[CurrentSystem].Atoms[f1][D].Type==ModifyFrameworkAtomTypes[l][3]))||
                      ((Framework[CurrentSystem].Atoms[f1][D].Type==ModifyFrameworkAtomTypes[l][2])&&
                       (Framework[CurrentSystem].Atoms[f1][A].Type==ModifyFrameworkAtomTypes[l][3])))||
                     (((Framework[CurrentSystem].Atoms[f1][B].Type==ModifyFrameworkAtomTypes[l][2])&&
                       (Framework[CurrentSystem].Atoms[f1][C].Type==ModifyFrameworkAtomTypes[l][3]))||
                      ((Framework[CurrentSystem].Atoms[f1][C].Type==ModifyFrameworkAtomTypes[l][2])&&
                       (Framework[CurrentSystem].Atoms[f1][B].Type==ModifyFrameworkAtomTypes[l][3])))||
                     (((Framework[CurrentSystem].Atoms[f1][B].Type==ModifyFrameworkAtomTypes[l][2])&&
                       (Framework[CurrentSystem].Atoms[f1][D].Type==ModifyFrameworkAtomTypes[l][3]))||
                      ((Framework[CurrentSystem].Atoms[f1][D].Type==ModifyFrameworkAtomTypes[l][2])&&
                       (Framework[CurrentSystem].Atoms[f1][B].Type==ModifyFrameworkAtomTypes[l][3])))||
                     (((Framework[CurrentSystem].Atoms[f1][C].Type==ModifyFrameworkAtomTypes[l][2])&&
                       (Framework[CurrentSystem].Atoms[f1][D].Type==ModifyFrameworkAtomTypes[l][3]))||
                      ((Framework[CurrentSystem].Atoms[f1][D].Type==ModifyFrameworkAtomTypes[l][2])&&
                       (Framework[CurrentSystem].Atoms[f1][C].Type==ModifyFrameworkAtomTypes[l][3])))))
                  {
                    Framework[CurrentSystem].Atoms[f1][i].Type=ModifyFrameworkAtomTypes[l][1];
                    if(!UseChargesFromCIFFile)
                      Framework[CurrentSystem].Atoms[f1][i].Charge=PseudoAtoms[ModifyFrameworkAtomTypes[l][1]].Charge1;

                    NumberOfPseudoAtomsCount[CurrentSystem][ModifyFrameworkAtomTypes[l][1]]++;
                    NumberOfPseudoAtomsType[CurrentSystem][ModifyFrameworkAtomTypes[l][1]]++;
                    NumberOfPseudoAtomsCount[CurrentSystem][ModifyFrameworkAtomTypes[l][0]]--;
                    NumberOfPseudoAtomsType[CurrentSystem][ModifyFrameworkAtomTypes[l][0]]--;
                  }
                  break;
                case 5:
                  A=GetNeighbour(CurrentSystem,f1,i,0);
                  B=GetNeighbour(CurrentSystem,f1,i,1);
                  C=GetNeighbour(CurrentSystem,f1,i,2);
                  D=GetNeighbour(CurrentSystem,f1,i,3);
                  E=GetNeighbour(CurrentSystem,f1,i,4);
                  if((Framework[CurrentSystem].Atoms[f1][i].Type==ModifyFrameworkAtomTypes[l][0])&&
                    ((((Framework[CurrentSystem].Atoms[f1][A].Type==ModifyFrameworkAtomTypes[l][2])&&
                       (Framework[CurrentSystem].Atoms[f1][B].Type==ModifyFrameworkAtomTypes[l][3]))||
                      ((Framework[CurrentSystem].Atoms[f1][B].Type==ModifyFrameworkAtomTypes[l][2])&&
                       (Framework[CurrentSystem].Atoms[f1][A].Type==ModifyFrameworkAtomTypes[l][3])))||
                     (((Framework[CurrentSystem].Atoms[f1][A].Type==ModifyFrameworkAtomTypes[l][2])&&
                       (Framework[CurrentSystem].Atoms[f1][C].Type==ModifyFrameworkAtomTypes[l][3]))||
                      ((Framework[CurrentSystem].Atoms[f1][C].Type==ModifyFrameworkAtomTypes[l][2])&&
                       (Framework[CurrentSystem].Atoms[f1][A].Type==ModifyFrameworkAtomTypes[l][3])))||
                     (((Framework[CurrentSystem].Atoms[f1][A].Type==ModifyFrameworkAtomTypes[l][2])&&
                       (Framework[CurrentSystem].Atoms[f1][D].Type==ModifyFrameworkAtomTypes[l][3]))||
                      ((Framework[CurrentSystem].Atoms[f1][D].Type==ModifyFrameworkAtomTypes[l][2])&&
                       (Framework[CurrentSystem].Atoms[f1][A].Type==ModifyFrameworkAtomTypes[l][3])))||
                     (((Framework[CurrentSystem].Atoms[f1][A].Type==ModifyFrameworkAtomTypes[l][2])&&
                       (Framework[CurrentSystem].Atoms[f1][E].Type==ModifyFrameworkAtomTypes[l][3]))||
                      ((Framework[CurrentSystem].Atoms[f1][E].Type==ModifyFrameworkAtomTypes[l][2])&&
                       (Framework[CurrentSystem].Atoms[f1][A].Type==ModifyFrameworkAtomTypes[l][3])))||
                     (((Framework[CurrentSystem].Atoms[f1][B].Type==ModifyFrameworkAtomTypes[l][2])&&
                       (Framework[CurrentSystem].Atoms[f1][C].Type==ModifyFrameworkAtomTypes[l][3]))||
                      ((Framework[CurrentSystem].Atoms[f1][C].Type==ModifyFrameworkAtomTypes[l][2])&&
                       (Framework[CurrentSystem].Atoms[f1][B].Type==ModifyFrameworkAtomTypes[l][3])))||
                     (((Framework[CurrentSystem].Atoms[f1][B].Type==ModifyFrameworkAtomTypes[l][2])&&
                       (Framework[CurrentSystem].Atoms[f1][D].Type==ModifyFrameworkAtomTypes[l][3]))||
                      ((Framework[CurrentSystem].Atoms[f1][D].Type==ModifyFrameworkAtomTypes[l][2])&&
                       (Framework[CurrentSystem].Atoms[f1][B].Type==ModifyFrameworkAtomTypes[l][3])))||
                     (((Framework[CurrentSystem].Atoms[f1][B].Type==ModifyFrameworkAtomTypes[l][2])&&
                       (Framework[CurrentSystem].Atoms[f1][E].Type==ModifyFrameworkAtomTypes[l][3]))||
                      ((Framework[CurrentSystem].Atoms[f1][E].Type==ModifyFrameworkAtomTypes[l][2])&&
                       (Framework[CurrentSystem].Atoms[f1][B].Type==ModifyFrameworkAtomTypes[l][3])))||
                     (((Framework[CurrentSystem].Atoms[f1][C].Type==ModifyFrameworkAtomTypes[l][2])&&
                       (Framework[CurrentSystem].Atoms[f1][D].Type==ModifyFrameworkAtomTypes[l][3]))||
                      ((Framework[CurrentSystem].Atoms[f1][D].Type==ModifyFrameworkAtomTypes[l][2])&&
                       (Framework[CurrentSystem].Atoms[f1][C].Type==ModifyFrameworkAtomTypes[l][3])))||
                     (((Framework[CurrentSystem].Atoms[f1][C].Type==ModifyFrameworkAtomTypes[l][2])&&
                       (Framework[CurrentSystem].Atoms[f1][E].Type==ModifyFrameworkAtomTypes[l][3]))||
                      ((Framework[CurrentSystem].Atoms[f1][E].Type==ModifyFrameworkAtomTypes[l][2])&&
                       (Framework[CurrentSystem].Atoms[f1][C].Type==ModifyFrameworkAtomTypes[l][3])))||
                     (((Framework[CurrentSystem].Atoms[f1][D].Type==ModifyFrameworkAtomTypes[l][2])&&
                       (Framework[CurrentSystem].Atoms[f1][E].Type==ModifyFrameworkAtomTypes[l][3]))||
                      ((Framework[CurrentSystem].Atoms[f1][E].Type==ModifyFrameworkAtomTypes[l][2])&&
                       (Framework[CurrentSystem].Atoms[f1][D].Type==ModifyFrameworkAtomTypes[l][3])))))
                  {
                    Framework[CurrentSystem].Atoms[f1][i].Type=ModifyFrameworkAtomTypes[l][1];
                    if(!UseChargesFromCIFFile)
                      Framework[CurrentSystem].Atoms[f1][i].Charge=PseudoAtoms[ModifyFrameworkAtomTypes[l][1]].Charge1;

                    NumberOfPseudoAtomsCount[CurrentSystem][ModifyFrameworkAtomTypes[l][1]]++;
                    NumberOfPseudoAtomsType[CurrentSystem][ModifyFrameworkAtomTypes[l][1]]++;
                    NumberOfPseudoAtomsCount[CurrentSystem][ModifyFrameworkAtomTypes[l][0]]--;
                    NumberOfPseudoAtomsType[CurrentSystem][ModifyFrameworkAtomTypes[l][0]]--;
                  }

                  break;
                case 6:
                  A=GetNeighbour(CurrentSystem,f1,i,0);
                  B=GetNeighbour(CurrentSystem,f1,i,1);
                  C=GetNeighbour(CurrentSystem,f1,i,2);
                  D=GetNeighbour(CurrentSystem,f1,i,3);
                  E=GetNeighbour(CurrentSystem,f1,i,4);
                  F=GetNeighbour(CurrentSystem,f1,i,5);
                  if((Framework[CurrentSystem].Atoms[f1][i].Type==ModifyFrameworkAtomTypes[l][0])&&
                    ((((Framework[CurrentSystem].Atoms[f1][A].Type==ModifyFrameworkAtomTypes[l][2])&&
                       (Framework[CurrentSystem].Atoms[f1][B].Type==ModifyFrameworkAtomTypes[l][3]))||
                      ((Framework[CurrentSystem].Atoms[f1][B].Type==ModifyFrameworkAtomTypes[l][2])&&
                       (Framework[CurrentSystem].Atoms[f1][A].Type==ModifyFrameworkAtomTypes[l][3])))||
                     (((Framework[CurrentSystem].Atoms[f1][A].Type==ModifyFrameworkAtomTypes[l][2])&&
                       (Framework[CurrentSystem].Atoms[f1][C].Type==ModifyFrameworkAtomTypes[l][3]))||
                      ((Framework[CurrentSystem].Atoms[f1][C].Type==ModifyFrameworkAtomTypes[l][2])&&
                       (Framework[CurrentSystem].Atoms[f1][A].Type==ModifyFrameworkAtomTypes[l][3])))||
                     (((Framework[CurrentSystem].Atoms[f1][A].Type==ModifyFrameworkAtomTypes[l][2])&&
                       (Framework[CurrentSystem].Atoms[f1][D].Type==ModifyFrameworkAtomTypes[l][3]))||
                      ((Framework[CurrentSystem].Atoms[f1][D].Type==ModifyFrameworkAtomTypes[l][2])&&
                       (Framework[CurrentSystem].Atoms[f1][A].Type==ModifyFrameworkAtomTypes[l][3])))||
                     (((Framework[CurrentSystem].Atoms[f1][A].Type==ModifyFrameworkAtomTypes[l][2])&&
                       (Framework[CurrentSystem].Atoms[f1][E].Type==ModifyFrameworkAtomTypes[l][3]))||
                      ((Framework[CurrentSystem].Atoms[f1][E].Type==ModifyFrameworkAtomTypes[l][2])&&
                       (Framework[CurrentSystem].Atoms[f1][A].Type==ModifyFrameworkAtomTypes[l][3])))||
                     (((Framework[CurrentSystem].Atoms[f1][A].Type==ModifyFrameworkAtomTypes[l][2])&&
                       (Framework[CurrentSystem].Atoms[f1][F].Type==ModifyFrameworkAtomTypes[l][3]))||
                      ((Framework[CurrentSystem].Atoms[f1][F].Type==ModifyFrameworkAtomTypes[l][2])&&
                       (Framework[CurrentSystem].Atoms[f1][A].Type==ModifyFrameworkAtomTypes[l][3])))||
                     (((Framework[CurrentSystem].Atoms[f1][B].Type==ModifyFrameworkAtomTypes[l][2])&&
                       (Framework[CurrentSystem].Atoms[f1][C].Type==ModifyFrameworkAtomTypes[l][3]))||
                      ((Framework[CurrentSystem].Atoms[f1][C].Type==ModifyFrameworkAtomTypes[l][2])&&
                       (Framework[CurrentSystem].Atoms[f1][B].Type==ModifyFrameworkAtomTypes[l][3])))||
                     (((Framework[CurrentSystem].Atoms[f1][B].Type==ModifyFrameworkAtomTypes[l][2])&&
                       (Framework[CurrentSystem].Atoms[f1][D].Type==ModifyFrameworkAtomTypes[l][3]))||
                      ((Framework[CurrentSystem].Atoms[f1][D].Type==ModifyFrameworkAtomTypes[l][2])&&
                       (Framework[CurrentSystem].Atoms[f1][B].Type==ModifyFrameworkAtomTypes[l][3])))||
                     (((Framework[CurrentSystem].Atoms[f1][B].Type==ModifyFrameworkAtomTypes[l][2])&&
                       (Framework[CurrentSystem].Atoms[f1][E].Type==ModifyFrameworkAtomTypes[l][3]))||
                      ((Framework[CurrentSystem].Atoms[f1][E].Type==ModifyFrameworkAtomTypes[l][2])&&
                       (Framework[CurrentSystem].Atoms[f1][B].Type==ModifyFrameworkAtomTypes[l][3])))||
                     (((Framework[CurrentSystem].Atoms[f1][B].Type==ModifyFrameworkAtomTypes[l][2])&&
                       (Framework[CurrentSystem].Atoms[f1][F].Type==ModifyFrameworkAtomTypes[l][3]))||
                      ((Framework[CurrentSystem].Atoms[f1][F].Type==ModifyFrameworkAtomTypes[l][2])&&
                       (Framework[CurrentSystem].Atoms[f1][B].Type==ModifyFrameworkAtomTypes[l][3])))||
                     (((Framework[CurrentSystem].Atoms[f1][C].Type==ModifyFrameworkAtomTypes[l][2])&&
                       (Framework[CurrentSystem].Atoms[f1][D].Type==ModifyFrameworkAtomTypes[l][3]))||
                      ((Framework[CurrentSystem].Atoms[f1][D].Type==ModifyFrameworkAtomTypes[l][2])&&
                       (Framework[CurrentSystem].Atoms[f1][C].Type==ModifyFrameworkAtomTypes[l][3])))||
                     (((Framework[CurrentSystem].Atoms[f1][C].Type==ModifyFrameworkAtomTypes[l][2])&&
                       (Framework[CurrentSystem].Atoms[f1][E].Type==ModifyFrameworkAtomTypes[l][3]))||
                      ((Framework[CurrentSystem].Atoms[f1][E].Type==ModifyFrameworkAtomTypes[l][2])&&
                       (Framework[CurrentSystem].Atoms[f1][C].Type==ModifyFrameworkAtomTypes[l][3])))||
                     (((Framework[CurrentSystem].Atoms[f1][C].Type==ModifyFrameworkAtomTypes[l][2])&&
                       (Framework[CurrentSystem].Atoms[f1][F].Type==ModifyFrameworkAtomTypes[l][3]))||
                      ((Framework[CurrentSystem].Atoms[f1][F].Type==ModifyFrameworkAtomTypes[l][2])&&
                       (Framework[CurrentSystem].Atoms[f1][C].Type==ModifyFrameworkAtomTypes[l][3])))||
                     (((Framework[CurrentSystem].Atoms[f1][D].Type==ModifyFrameworkAtomTypes[l][2])&&
                       (Framework[CurrentSystem].Atoms[f1][E].Type==ModifyFrameworkAtomTypes[l][3]))||
                      ((Framework[CurrentSystem].Atoms[f1][E].Type==ModifyFrameworkAtomTypes[l][2])&&
                       (Framework[CurrentSystem].Atoms[f1][D].Type==ModifyFrameworkAtomTypes[l][3])))||
                     (((Framework[CurrentSystem].Atoms[f1][D].Type==ModifyFrameworkAtomTypes[l][2])&&
                       (Framework[CurrentSystem].Atoms[f1][F].Type==ModifyFrameworkAtomTypes[l][3]))||
                      ((Framework[CurrentSystem].Atoms[f1][F].Type==ModifyFrameworkAtomTypes[l][2])&&
                       (Framework[CurrentSystem].Atoms[f1][D].Type==ModifyFrameworkAtomTypes[l][3])))||
                     (((Framework[CurrentSystem].Atoms[f1][E].Type==ModifyFrameworkAtomTypes[l][2])&&
                       (Framework[CurrentSystem].Atoms[f1][F].Type==ModifyFrameworkAtomTypes[l][3]))||
                      ((Framework[CurrentSystem].Atoms[f1][F].Type==ModifyFrameworkAtomTypes[l][2])&&
                       (Framework[CurrentSystem].Atoms[f1][E].Type==ModifyFrameworkAtomTypes[l][3])))))
                  {
                    Framework[CurrentSystem].Atoms[f1][i].Type=ModifyFrameworkAtomTypes[l][1];
                    if(!UseChargesFromCIFFile)
                      Framework[CurrentSystem].Atoms[f1][i].Charge=PseudoAtoms[ModifyFrameworkAtomTypes[l][1]].Charge1;

                    NumberOfPseudoAtomsCount[CurrentSystem][ModifyFrameworkAtomTypes[l][1]]++;
                    NumberOfPseudoAtomsType[CurrentSystem][ModifyFrameworkAtomTypes[l][1]]++;
                    NumberOfPseudoAtomsCount[CurrentSystem][ModifyFrameworkAtomTypes[l][0]]--;
                    NumberOfPseudoAtomsType[CurrentSystem][ModifyFrameworkAtomTypes[l][0]]--;
                  }
                  break;
                default:
                  fprintf(stderr, "Routine ModifyAtomsConnectedToDefinedNeighbours: case of >6 neighbors not implemented\n");
                  exit(0);
                  break;
              }
            }
          }
        }
        break;
    }
  }

  for(l=0;l<NumberOfForbiddenConnectivityRules;l++)
  {
    for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
    {
      for(A=0;A<Framework[CurrentSystem].NumberOfAtoms[f1];A++)
      {
        if(Framework[CurrentSystem].Atoms[f1][A].Type==ForbiddenConnectivityTypes[l][0])
        {
          for(j=0;j<Framework[CurrentSystem].Connectivity[f1][A];j++)
          {
            B=GetNeighbour(CurrentSystem,f1,A,j);
            if(Framework[CurrentSystem].Atoms[f1][B].Type==ForbiddenConnectivityTypes[l][1])
            {
              for(k=0;k<Framework[CurrentSystem].Connectivity[f1][B];k++)
              {
                C=GetNeighbour(CurrentSystem,f1,B,k);
                if((Framework[CurrentSystem].Atoms[f1][C].Type==ForbiddenConnectivityTypes[l][2])&&(A!=C))
                {
                  fprintf(stderr, "Illegal framework atom sequence (%s-%s-%s) found! (rule: %d)\n",
                    ForbiddenConnectivityAtoms[l][0],ForbiddenConnectivityAtoms[l][1],ForbiddenConnectivityAtoms[l][2],
                    l);
                  exit(0);
                }
              }
            }
          }
        }
      }
    }
  }
}

/*********************************************************************************************************
 * Name       | ReadBlockingPockets                                                                      *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Read the "blocking-pockets" from a file.                                                 *
 * Parameters | -                                                                                        *
 *********************************************************************************************************/

void ReadBlockingPockets(void)
{
  int i,j,k,l,index,nr_pockets;
  double temp;
  VECTOR tempr;
  char buffer[1024];
  FILE *FilePtr;
  VECTOR vec;

  for(CurrentComponent=0;CurrentComponent<NumberOfComponents;CurrentComponent++)
  {
    if(Components[CurrentComponent].BlockPockets[CurrentSystem])
    {
      sprintf(buffer,"./%s.%s",Components[CurrentComponent].BlockPocketsFilename[CurrentSystem],"block");
      if(!(FilePtr=fopen(buffer,"r")))
      {
        sprintf(buffer,"%s/share/raspa/structures/block/%s.%s",
                RASPA_DIRECTORY,
                Components[CurrentComponent].BlockPocketsFilename[CurrentSystem],
                "block");
        if(!(FilePtr=fopen(buffer,"r")))
        {
          fprintf(stderr, "'Blocking-pocket' file not found and therefore not used\n");
          return;
        }
      }

      fscanf(FilePtr,"%d\n",&Components[CurrentComponent].NumberOfBlockCenters[CurrentSystem]);

      nr_pockets=Components[CurrentComponent].NumberOfBlockCenters[CurrentSystem]*NumberOfUnitCells[CurrentSystem].x*
                 NumberOfUnitCells[CurrentSystem].y*NumberOfUnitCells[CurrentSystem].z;

      Components[CurrentComponent].BlockDistance[CurrentSystem]=(REAL*)calloc(nr_pockets,sizeof(REAL));
      Components[CurrentComponent].BlockCenters[CurrentSystem]=(VECTOR*)calloc(nr_pockets,sizeof(VECTOR));

      index=0;
      for(i=0;i<Components[CurrentComponent].NumberOfBlockCenters[CurrentSystem];i++)
      {
        fscanf(FilePtr,"%lf %lf %lf %lf\n",
              &tempr.x,&tempr.y,&tempr.z,&temp);
        tempr.x+=Framework[CurrentSystem].ShiftUnitCell[0].x;
        tempr.y+=Framework[CurrentSystem].ShiftUnitCell[0].y;
        tempr.z+=Framework[CurrentSystem].ShiftUnitCell[0].z;
        for(j=0;j<NumberOfUnitCells[CurrentSystem].x;j++)
          for(k=0;k<NumberOfUnitCells[CurrentSystem].y;k++)
            for(l=0;l<NumberOfUnitCells[CurrentSystem].z;l++)
            {
              // convert to xyz (zeolite atoms are stored in xyz)

              vec.x=(tempr.x+j)/NumberOfUnitCells[CurrentSystem].x;
              vec.y=(tempr.y+k)/NumberOfUnitCells[CurrentSystem].y;
              vec.z=(tempr.z+l)/NumberOfUnitCells[CurrentSystem].z;

              Components[CurrentComponent].BlockCenters[CurrentSystem][index]=ConvertFromABCtoXYZ(vec);
              Components[CurrentComponent].BlockDistance[CurrentSystem][index]=temp;

              index++;
            }
      }
      Components[CurrentComponent].NumberOfBlockCenters[CurrentSystem]*=
         NumberOfUnitCells[CurrentSystem].x*NumberOfUnitCells[CurrentSystem].y*NumberOfUnitCells[CurrentSystem].z;
      fclose(FilePtr);
    }
  }
}

/*********************************************************************************************************
 * Name       | BlockedPocket                                                                            *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Check whether a position is blocked or not.                                              *
 * Parameters | VECTOR pos: the position of the 'probe-point'.                                           *
 *********************************************************************************************************/

int BlockedPocket(VECTOR pos)
{
  int i;
  VECTOR dr;
  REAL r;
  int blocked;

  if(NumberOfComponents==0) return FALSE;

  if(Components[CurrentComponent].BlockPockets[CurrentSystem])
  {
    // molecules must be inside the block-pockets for 'InvertBlockPockets'
    if(Components[CurrentComponent].InvertBlockPockets)
    {
      for(i=0;i<Components[CurrentComponent].NumberOfBlockCenters[CurrentSystem];i++)
      {
        dr.x=Components[CurrentComponent].BlockCenters[CurrentSystem][i].x-pos.x;
        dr.y=Components[CurrentComponent].BlockCenters[CurrentSystem][i].y-pos.y;
        dr.z=Components[CurrentComponent].BlockCenters[CurrentSystem][i].z-pos.z;
        dr=ApplyBoundaryConditionUnitCell(dr);
        r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));

        // if inside blockPockets, then it is allowed (return 'false' for blocking)
        if(r<Components[CurrentComponent].BlockDistance[CurrentSystem][i])
        {
          return FALSE;
        }
      }

      // molecule is not in any of the blocking-pockets, so block it
      return TRUE; 
    }
    else
    {
      // molecules are blocked when inside any of the block-pockets
      for(i=0;i<Components[CurrentComponent].NumberOfBlockCenters[CurrentSystem];i++)
      {
        dr.x=Components[CurrentComponent].BlockCenters[CurrentSystem][i].x-pos.x;
        dr.y=Components[CurrentComponent].BlockCenters[CurrentSystem][i].y-pos.y;
        dr.z=Components[CurrentComponent].BlockCenters[CurrentSystem][i].z-pos.z;
        dr=ApplyBoundaryConditionUnitCell(dr);
        r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));

        // if inside block-pocket, then block (return 'true')
        if(r<Components[CurrentComponent].BlockDistance[CurrentSystem][i])
        {
          return TRUE;
        }
      }

      // molecule is not in any of the blocking-pockets, is is allowed (return 'false' for blocking)
      return FALSE;
    }
  }
  return FALSE;
}

int BlockedPocketGridTest(VECTOR pos)
{
  int i;
  VECTOR dr;
  REAL r;

  for(i=0;i<Components[CurrentComponent].NumberOfBlockCenters[CurrentSystem];i++)
  {
    dr.x=Components[CurrentComponent].BlockCenters[CurrentSystem][i].x-pos.x;
    dr.y=Components[CurrentComponent].BlockCenters[CurrentSystem][i].y-pos.y;
    dr.z=Components[CurrentComponent].BlockCenters[CurrentSystem][i].z-pos.z;
    dr=ApplyBoundaryConditionUnitCell(dr);
    r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
    if(r<Components[CurrentComponent].BlockDistance[CurrentSystem][i])
    {
      return TRUE;
    }
  }
  return FALSE;
}


/*********************************************************************************************************
 * Name       | CellProperties                                                                           *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Calculates the dimensional properties of a simulation cell.                              *
 * Parameters | REAL_MATRIX3x3 *in : the simulation cell in matrix form.                                 *
 *            | REAL_MATRIX3x3 *out: the computed properties of the cell.                                *
 *            | REAL *Volume       : the computed volume of the cell.                                    *
 * Note       | out[x][0] lengths of cell vectors.                                                       *
 *            | out[x][1] cosines of cell angles.                                                        *
 *            | out[x][2] perpendicular cell widths.                                                     *
 *            | The perpendicular width is the shortest distance between opposing cell faces.            *
 *********************************************************************************************************/

void CellProperties(REAL_MATRIX3x3 *in,REAL_MATRIX3x3 *out,REAL *Volume)
{
  REAL axb1,axb2,axb3,bxc1,bxc2,bxc3,cxa1,cxa2,cxa3;

  switch(Dimension)
  {
    case 2:
      // calculate lengths of cell vectors
      out->ax=sqrt(SQR(in->ax)+SQR(in->ay)+SQR(in->az));
      out->ay=sqrt(SQR(in->bx)+SQR(in->by)+SQR(in->bz));
      out->az=0.0;

      // calculate cosines of cell angles
      out->bx=0.0;
      out->by=0.0;
      out->bz=(in->ax*in->bx+in->ay*in->by+in->az*in->bz)/(out->ax*out->ay);

      // calculate volume of cell
      *Volume=in->ax*in->by-in->ay*in->bx;

      // calculate cell perpendicular widths
      out->cx=*Volume/out->ay;
      out->cy=*Volume/out->ax;
      out->cz=0.0;
      break;
    case 3:
      // calculate lengths of cell vectors
      out->ax=sqrt(SQR(in->ax)+SQR(in->ay)+SQR(in->az));
      out->ay=sqrt(SQR(in->bx)+SQR(in->by)+SQR(in->bz));
      out->az=sqrt(SQR(in->cx)+SQR(in->cy)+SQR(in->cz));

      // calculate cosines of cell angles
      out->bx=(in->bx*in->cx+in->by*in->cy+in->bz*in->cz)/(out->ay*out->az);
      out->by=(in->ax*in->cx+in->ay*in->cy+in->az*in->cz)/(out->ax*out->az);
      out->bz=(in->ax*in->bx+in->ay*in->by+in->az*in->bz)/(out->ax*out->ay);

      // calculate vector products of cell vectors
      axb1=in->ay*in->bz-in->az*in->by;
      axb2=in->az*in->bx-in->ax*in->bz;
      axb3=in->ax*in->by-in->ay*in->bx;

      bxc1=in->by*in->cz-in->bz*in->cy;
      bxc2=in->bz*in->cx-in->bx*in->cz;
      bxc3=in->bx*in->cy-in->by*in->cx;

      cxa1=in->cy*in->az-in->ay*in->cz;
      cxa2=in->ax*in->cz-in->az*in->cx;
      cxa3=in->ay*in->cx-in->ax*in->cy;

      // calculate volume of cell
      *Volume=fabs(in->ax*bxc1+in->ay*bxc2+in->az*bxc3);

      // calculate cell perpendicular widths
      out->cx=*Volume/sqrt(bxc1*bxc1+bxc2*bxc2+bxc3*bxc3);
      out->cy=*Volume/sqrt(cxa1*cxa1+cxa2*cxa2+cxa3*cxa3);
      out->cz=*Volume/sqrt(axb1*axb1+axb2*axb2+axb3*axb3);
      break;
  }
}

void PrintFrameworkStatus(void)
{
  int i;
  VECTOR pos;

  fprintf(stderr, "Framework contains %d atoms\n",Framework[CurrentSystem].TotalNumberOfAtoms);
  fprintf(stderr, "===================================================\n");
  for(i=0;i<Framework[CurrentSystem].TotalNumberOfAtoms;i++)
  {
    pos=Framework[CurrentSystem].Atoms[0][i].Position;
    fprintf(stderr, "%d %lf %lf %lf\n",i,(double)pos.x,(double)pos.y,(double)pos.z);
  }
}

VECTOR GetFrameworkCenterOfMassPosition(void)
{
  int i,f1;
  REAL Mass,TotalMass;
  VECTOR com;

  TotalMass=0.0;
  com.x=com.y=com.z=0.0;
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
    {
      Mass=PseudoAtoms[Framework[CurrentSystem].Atoms[f1][i].Type].Mass;
      TotalMass+=Mass;
      com.x=com.x+Mass*Framework[CurrentSystem].Atoms[f1][i].Position.x;
      com.y=com.y+Mass*Framework[CurrentSystem].Atoms[f1][i].Position.y;
      com.z=com.z+Mass*Framework[CurrentSystem].Atoms[f1][i].Position.z;
    }
  }
  com.x=com.x/TotalMass;
  com.y=com.y/TotalMass;
  com.z=com.z/TotalMass;
  return com;
}

VECTOR GetIndividualFrameworkCenterOfMass(int f1)
{
  int i;
  REAL Mass,TotalMass;
  VECTOR com;

  TotalMass=0.0;
  com.x=com.y=com.z=0.0;
  for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
  {
    Mass=PseudoAtoms[Framework[CurrentSystem].Atoms[f1][i].Type].Mass;
    TotalMass+=Mass;
    com.x=com.x+Mass*Framework[CurrentSystem].Atoms[f1][i].Position.x;
    com.y=com.y+Mass*Framework[CurrentSystem].Atoms[f1][i].Position.y;
    com.z=com.z+Mass*Framework[CurrentSystem].Atoms[f1][i].Position.z;
  }
  com.x=com.x/TotalMass;
  com.y=com.y/TotalMass;
  com.z=com.z/TotalMass;
  return com;
}

VECTOR GetFrameworkCenterOfMass(void)
{
  int i,f1;
  REAL Mass,TotalMass;
  VECTOR com;

  TotalMass=0.0;
  com.x=com.y=com.z=0.0;
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
    {
      Mass=PseudoAtoms[Framework[CurrentSystem].Atoms[f1][i].Type].Mass;
      TotalMass+=Mass;
      com.x=com.x+Mass*Framework[CurrentSystem].Atoms[f1][i].Position.x;
      com.y=com.y+Mass*Framework[CurrentSystem].Atoms[f1][i].Position.y;
      com.z=com.z+Mass*Framework[CurrentSystem].Atoms[f1][i].Position.z;
    }
  }
  com.x=com.x/TotalMass;
  com.y=com.y/TotalMass;
  com.z=com.z/TotalMass;
  return com;
}

VECTOR GetFrameworkCenterOfMassVelocity(void)
{
  int i,f1;
  REAL Mass,TotalMass;
  VECTOR com;

  TotalMass=0.0;
  com.x=com.y=com.z=0.0;
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
    {
      Mass=PseudoAtoms[Framework[CurrentSystem].Atoms[f1][i].Type].Mass;
      TotalMass+=Mass;
      com.x+=Mass*Framework[CurrentSystem].Atoms[f1][i].Velocity.x;
      com.y+=Mass*Framework[CurrentSystem].Atoms[f1][i].Velocity.y;
      com.z+=Mass*Framework[CurrentSystem].Atoms[f1][i].Velocity.z;
    }
  }
  com.x/=TotalMass;
  com.y/=TotalMass;
  com.z/=TotalMass;
  return com;
}

void InitializeFrameworkVelocities(void)
{
  int i,f1;
  int A,B;
  int nr_core_shells,nr_atoms;
  REAL Mass;

  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    nr_core_shells=Framework[CurrentSystem].NumberOfCoreShells[f1];
    nr_atoms=Framework[CurrentSystem].NumberOfAtoms[f1];

    for(i=0;i<nr_atoms-nr_core_shells;i++)
    {
      Mass=PseudoAtoms[Framework[CurrentSystem].Atoms[f1][i].Type].Mass;
      Framework[CurrentSystem].Atoms[f1][i].Velocity.x=0.0;
      Framework[CurrentSystem].Atoms[f1][i].Velocity.y=0.0;
      Framework[CurrentSystem].Atoms[f1][i].Velocity.z=0.0;
      if(!Framework[CurrentSystem].Atoms[f1][i].Fixed.x)
        Framework[CurrentSystem].Atoms[f1][i].Velocity.x=RandomGaussianNumber()*
          sqrt(K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]/Mass);
      if(!Framework[CurrentSystem].Atoms[f1][i].Fixed.y)
        Framework[CurrentSystem].Atoms[f1][i].Velocity.y=RandomGaussianNumber()*
          sqrt(K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]/Mass);
      if(!Framework[CurrentSystem].Atoms[f1][i].Fixed.z)
        Framework[CurrentSystem].Atoms[f1][i].Velocity.z=RandomGaussianNumber()*
          sqrt(K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]/Mass);
    }
    if(Framework[CurrentSystem].NumberOfCoreShells[f1]>0)
    {
      // loop over shells
      for(A=0;A<nr_atoms-nr_core_shells;A++)
      {
        // get the core for this shell
        B=Framework[CurrentSystem].CoreShellConnectivity[f1][A];
        if(B>0)
          Framework[CurrentSystem].Atoms[f1][B].Velocity=Framework[CurrentSystem].Atoms[f1][A].Velocity;
      }
    }
  }
}

VECTOR MeasureFrameworkVelocityDrift(void)
{
  int i,f1;
  REAL Mass,TotalMass;
  VECTOR com;

  TotalMass=0.0;
  com.x=com.y=com.z=0.0;
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
    {
      if(!Framework[CurrentSystem].Atoms[f1][i].Fixed.x)
      {
        TotalMass+=(Mass=PseudoAtoms[Framework[CurrentSystem].Atoms[f1][i].Type].Mass);
        com.x=com.x+Mass*Framework[CurrentSystem].Atoms[f1][i].Velocity.x;
        com.y=com.y+Mass*Framework[CurrentSystem].Atoms[f1][i].Velocity.y;
        com.z=com.z+Mass*Framework[CurrentSystem].Atoms[f1][i].Velocity.z;
      }
    }
  }
  com.x=com.x/TotalMass;
  com.y=com.y/TotalMass;
  com.z=com.z/TotalMass;
  return com;
}

void RemoveFrameworkVelocityDrift(void)
{
  int i,f1;
  VECTOR com;

  com=MeasureFrameworkVelocityDrift();
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
    {
      if(!Framework[CurrentSystem].Atoms[f1][i].Fixed.x)
      {
        Framework[CurrentSystem].Atoms[f1][i].Velocity.x-=com.x;
        Framework[CurrentSystem].Atoms[f1][i].Velocity.y-=com.y;
        Framework[CurrentSystem].Atoms[f1][i].Velocity.z-=com.z;
      }
    }
  }
  //DegreesOfFreedomFramework[0]-=3;
}



REAL GetFrameworkKineticEnergy(void)
{
  int i,f1;
  REAL Ukin,Mass;

  Ukin=0.0;

  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
    {
      if(!Framework[CurrentSystem].Atoms[f1][i].Fixed.x)
      {
        Mass=PseudoAtoms[Framework[CurrentSystem].Atoms[f1][i].Type].Mass;
        Ukin+=0.5*Mass*(SQR(Framework[CurrentSystem].Atoms[f1][i].Velocity.x)+
                        SQR(Framework[CurrentSystem].Atoms[f1][i].Velocity.y)+
                        SQR(Framework[CurrentSystem].Atoms[f1][i].Velocity.z));
      }
    }
  }
  return Ukin;
}


void SetFrameworkVelocitesToZero(void)
{
  int i,f1;

  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
    {
     Framework[CurrentSystem].Atoms[f1][i].Velocity.x=0.0;
     Framework[CurrentSystem].Atoms[f1][i].Velocity.y=0.0;
     Framework[CurrentSystem].Atoms[f1][i].Velocity.z=0.0;
    }
  }
}


void ScaleFrameworkVelocitesToTemperature(void)
{
  int i,f1;
  REAL Factor;

  UHostKinetic[CurrentSystem]=GetFrameworkKineticEnergy();

  Factor=sqrt(DegreesOfFreedomFramework[CurrentSystem]*K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]/(2.0*UHostKinetic[CurrentSystem]));
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
    {
     Framework[CurrentSystem].Atoms[f1][i].Velocity.x*=Factor;
     Framework[CurrentSystem].Atoms[f1][i].Velocity.y*=Factor;
     Framework[CurrentSystem].Atoms[f1][i].Velocity.z*=Factor;
    }
  }
  UHostKinetic[CurrentSystem]=GetFrameworkKineticEnergy();
}

int BondNeighbours(int A,int B)
{
  int i;

  for(i=0;i<Framework[CurrentSystem].Connectivity[CurrentFramework][A];i++)
    if(Framework[CurrentSystem].Neighbours[CurrentFramework][A][i]==B) return TRUE;
  return FALSE;
}

int BendNeighbours(int A,int C)
{
  int i,j,B;

  for(i=0;i<Framework[CurrentSystem].Connectivity[CurrentFramework][A];i++)
  {
    B=Framework[CurrentSystem].Neighbours[CurrentFramework][A][i];
    for(j=0;j<Framework[CurrentSystem].Connectivity[CurrentFramework][B];j++)
      if(Framework[CurrentSystem].Neighbours[CurrentFramework][B][j]==C) return TRUE;
  }
  return FALSE;
}

int TorsionNeighbours(int A,int D)
{
  int i,j,k,B,C;

  for(i=0;i<Framework[CurrentSystem].Connectivity[CurrentFramework][A];i++)
  {
    B=Framework[CurrentSystem].Neighbours[CurrentFramework][A][i];
    for(j=0;j<Framework[CurrentSystem].Connectivity[CurrentFramework][B];j++)
    {
      C=Framework[CurrentSystem].Neighbours[CurrentFramework][B][j];
      for(k=0;k<Framework[CurrentSystem].Connectivity[CurrentFramework][C];k++)
        if(Framework[CurrentSystem].Neighbours[CurrentFramework][C][k]==D) return TRUE;
    }
  }
  return FALSE;
}

int ReturnDipoleIndex(int f1,int A,int B)
{
  int i;

  for(i=0;i<Framework[CurrentSystem].NumberOfBondDipoles[f1];i++)
  {
    if(((Framework[CurrentSystem].BondDipoles[f1][i].A==A)&&(Framework[CurrentSystem].BondDipoles[f1][i].B==B))||
      ((Framework[CurrentSystem].BondDipoles[f1][i].A==B)&&(Framework[CurrentSystem].BondDipoles[f1][i].B==A)))
      return i;
  }
  return -1;
}

void AddBondTypeToDefinitions(int TypeA,int TypeB,int BondType,REAL *parms)
{
  int i,j,index;
  int AlreadyPresent;
  int NumberOfArguments;
  REAL parameters[10];

  if(strcasecmp(BondTypes[BondType].Name,"MORSE_BOND")==0)
  {
    NumberOfArguments=3;
    parameters[0]=parms[0]*ENERGY_TO_KELVIN;
    parameters[1]=parms[1];
    parameters[2]=parms[2];
  }
  else if(strcasecmp(BondTypes[BondType].Name,"HARMONIC_BOND")==0)
  {
    NumberOfArguments=2;
    parameters[0]=parms[0]*ENERGY_TO_KELVIN;
    parameters[1]=parms[1];
  }

  if(Framework[CurrentSystem].NumberOfBondsDefinitions==0)
  {
    index=0;
    Framework[CurrentSystem].NumberOfBondsDefinitions=1;
    Framework[CurrentSystem].BondDefinitionType=(int*)calloc(Framework[CurrentSystem].NumberOfBondsDefinitions,sizeof(int));
    Framework[CurrentSystem].BondDefinitions=(PAIR*)calloc(Framework[CurrentSystem].NumberOfBondsDefinitions,sizeof(PAIR));
    Framework[CurrentSystem].NumberOfBondsPerType=(int*)calloc(Framework[CurrentSystem].NumberOfBondsDefinitions,sizeof(int));
    Framework[CurrentSystem].BondArgumentDefinitions=(REAL(*)[MAX_BOND_POTENTIAL_ARGUMENTS])
       calloc(Framework[CurrentSystem].NumberOfBondsDefinitions,sizeof(REAL[MAX_BOND_POTENTIAL_ARGUMENTS]));

    Framework[CurrentSystem].BondDefinitionType[index]=BondType;
    Framework[CurrentSystem].BondDefinitions[index].A=MIN2(TypeA,TypeB);
    Framework[CurrentSystem].BondDefinitions[index].B=MAX2(TypeA,TypeB);
    Framework[CurrentSystem].NumberOfBondsPerType[index]++;
    BondTypes[Framework[CurrentFramework].BondDefinitionType[index]].nr_args=NumberOfArguments;
    for(j=0;j<NumberOfArguments;j++)
      Framework[CurrentSystem].BondArgumentDefinitions[index][j]=parameters[j];
  }
  else
  {
    AlreadyPresent=FALSE;
    for(i=0;i<Framework[CurrentSystem].NumberOfBondsDefinitions;i++)
    {
      if(((Framework[CurrentSystem].BondDefinitions[i].A==TypeA)&&(Framework[CurrentSystem].BondDefinitions[i].B==TypeB))||
         ((Framework[CurrentSystem].BondDefinitions[i].A==TypeB)&&(Framework[CurrentSystem].BondDefinitions[i].B==TypeA)))
      {
        AlreadyPresent=TRUE;
        for(j=0;j<NumberOfArguments;j++)
          AlreadyPresent=AlreadyPresent&&(fabs(parameters[j]-Framework[CurrentSystem].BondArgumentDefinitions[i][j])<1e-5);

        if(AlreadyPresent)
        {
          Framework[CurrentSystem].NumberOfBondsPerType[i]++;
          break;
        }
      }
    }
    if(!AlreadyPresent)
    {
      index=Framework[CurrentSystem].NumberOfBondsDefinitions;
      Framework[CurrentSystem].NumberOfBondsDefinitions++;
      Framework[CurrentSystem].BondDefinitionType=(int*)realloc(Framework[CurrentSystem].BondDefinitionType,
          Framework[CurrentSystem].NumberOfBondsDefinitions*sizeof(int));
      Framework[CurrentSystem].BondDefinitions=(PAIR*)realloc(Framework[CurrentSystem].BondDefinitions,
          Framework[CurrentSystem].NumberOfBondsDefinitions*sizeof(PAIR));
      Framework[CurrentSystem].NumberOfBondsPerType=(int*)realloc(Framework[CurrentSystem].NumberOfBondsPerType,
          Framework[CurrentSystem].NumberOfBondsDefinitions*sizeof(int));
      Framework[CurrentSystem].BondArgumentDefinitions=(REAL(*)[MAX_BOND_POTENTIAL_ARGUMENTS])
          realloc(Framework[CurrentSystem].BondArgumentDefinitions,
          Framework[CurrentSystem].NumberOfBondsDefinitions*sizeof(REAL[MAX_BOND_POTENTIAL_ARGUMENTS]));

      Framework[CurrentSystem].BondDefinitionType[index]=BondType;
      Framework[CurrentSystem].BondDefinitions[index].A=MIN2(TypeA,TypeB);
      Framework[CurrentSystem].BondDefinitions[index].B=MAX2(TypeA,TypeB);
      Framework[CurrentSystem].NumberOfBondsPerType[index]=1;
      BondTypes[Framework[CurrentFramework].BondDefinitionType[index]].nr_args=NumberOfArguments;
      for(j=0;j<NumberOfArguments;j++)
        Framework[CurrentSystem].BondArgumentDefinitions[index][j]=parameters[j];
    }

  }
}

void AddBendTypeToDefinitions(int TypeA,int TypeB,int TypeC,int BendType,REAL *parms)
{
  int i,j,index;
  int AlreadyPresent;
  int NumberOfArguments;
  REAL parameters[10];

  if(strcasecmp(BendTypes[BendType].Name,"HARMONIC_BEND")==0)
  {
    NumberOfArguments=2;
    parameters[0]=parms[0]*ENERGY_TO_KELVIN;
    parameters[1]=parms[1]*RAD2DEG;
  }

  if(Framework[CurrentSystem].NumberOfBendDefinitions==0)
  {
    index=0;
    Framework[CurrentSystem].NumberOfBendDefinitions=1;
    Framework[CurrentSystem].BendDefinitionType=(int*)calloc(Framework[CurrentSystem].NumberOfBendDefinitions,sizeof(int));
    Framework[CurrentSystem].BendDefinitions=(QUAD*)calloc(Framework[CurrentSystem].NumberOfBendDefinitions,sizeof(QUAD));
    Framework[CurrentSystem].NumberOfBendsPerType=(int*)calloc(Framework[CurrentSystem].NumberOfBendDefinitions,sizeof(int));
    Framework[CurrentSystem].BendArgumentDefinitions=(REAL(*)[MAX_BEND_POTENTIAL_ARGUMENTS])
       calloc(Framework[CurrentSystem].NumberOfBendDefinitions,sizeof(REAL[MAX_BEND_POTENTIAL_ARGUMENTS]));


    Framework[CurrentSystem].BendDefinitionType[index]=BendType;
    Framework[CurrentSystem].BendDefinitions[index].A=TypeA;
    Framework[CurrentSystem].BendDefinitions[index].B=TypeB;
    Framework[CurrentSystem].BendDefinitions[index].C=TypeC;
    Framework[CurrentSystem].NumberOfBendsPerType[index]++;
    BendTypes[Framework[CurrentFramework].BendDefinitionType[index]].nr_args=NumberOfArguments;
    for(j=0;j<NumberOfArguments;j++)
      Framework[CurrentSystem].BendArgumentDefinitions[index][j]=parameters[j];
  }
  else
  {
    AlreadyPresent=FALSE;
    for(i=0;i<Framework[CurrentSystem].NumberOfBendDefinitions;i++)
    {
      if(((Framework[CurrentSystem].BendDefinitions[i].A==TypeA)&&(Framework[CurrentSystem].BendDefinitions[i].B==TypeB)&&
         (Framework[CurrentSystem].BendDefinitions[i].C==TypeC))||
         ((Framework[CurrentSystem].BendDefinitions[i].A==TypeC)&&(Framework[CurrentSystem].BendDefinitions[i].B==TypeB)&&
         (Framework[CurrentSystem].BendDefinitions[i].C==TypeA)))
      {
        AlreadyPresent=TRUE;
        for(j=0;j<NumberOfArguments;j++)
          AlreadyPresent=AlreadyPresent&&(fabs(parameters[j]-Framework[CurrentSystem].BendArgumentDefinitions[i][j])<1e-5);;

        if(AlreadyPresent)
        {
          Framework[CurrentSystem].NumberOfBendsPerType[i]++;
          break;
        }
      }
    }
    if(!AlreadyPresent)
    {
      index=Framework[CurrentSystem].NumberOfBendDefinitions;
      Framework[CurrentSystem].NumberOfBendDefinitions++;
      Framework[CurrentSystem].BendDefinitionType=(int*)realloc(Framework[CurrentSystem].BendDefinitionType,
          Framework[CurrentSystem].NumberOfBendDefinitions*sizeof(int));
      Framework[CurrentSystem].BendDefinitions=(QUAD*)realloc(Framework[CurrentSystem].BendDefinitions,
          Framework[CurrentSystem].NumberOfBendDefinitions*sizeof(QUAD));
      Framework[CurrentSystem].NumberOfBendsPerType=(int*)realloc(Framework[CurrentSystem].NumberOfBendsPerType,
          Framework[CurrentSystem].NumberOfBendDefinitions*sizeof(int));
      Framework[CurrentSystem].BendArgumentDefinitions=(REAL(*)[MAX_BEND_POTENTIAL_ARGUMENTS])
          realloc(Framework[CurrentSystem].BendArgumentDefinitions,
          Framework[CurrentSystem].NumberOfBendDefinitions*sizeof(REAL[MAX_BEND_POTENTIAL_ARGUMENTS]));

      Framework[CurrentSystem].BendDefinitionType[index]=BendType;
      Framework[CurrentSystem].BendDefinitions[index].A=TypeA;
      Framework[CurrentSystem].BendDefinitions[index].B=TypeB;
      Framework[CurrentSystem].BendDefinitions[index].C=TypeC;
      Framework[CurrentSystem].NumberOfBendsPerType[index]=1;
      BendTypes[Framework[CurrentFramework].BendDefinitionType[index]].nr_args=NumberOfArguments;
      for(j=0;j<NumberOfArguments;j++)
        Framework[CurrentSystem].BendArgumentDefinitions[index][j]=parameters[j];
    }

  }
}

void AddTorsionTypeToDefinitions(int TypeA,int TypeB,int TypeC,int TypeD,int TorsionType,REAL *parms)
{
  int i,j,index;
  int AlreadyPresent;
  int CompareParameters;
  int NumberOfArguments;
  REAL parameters[10];

  if(strcasecmp(TorsionTypes[TorsionType].Name,"TRAPPE_DIHEDRAL")==0)
  {
    NumberOfArguments=8;
    parameters[0]=parms[0]*ENERGY_TO_KELVIN;
    parameters[1]=parms[1]*ENERGY_TO_KELVIN;
    parameters[2]=parms[2]*ENERGY_TO_KELVIN;
    parameters[3]=parms[3]*ENERGY_TO_KELVIN;
    parameters[4]=parms[4];
    parameters[5]=parms[5];
    parameters[6]=parms[6];
    parameters[7]=parms[7];
  }

  if(Framework[CurrentSystem].NumberOfTorsionDefinitions==0)
  {
    index=0;
    Framework[CurrentSystem].NumberOfTorsionDefinitions=1;
    Framework[CurrentSystem].TorsionDefinitionType=(int*)calloc(Framework[CurrentSystem].NumberOfTorsionDefinitions,sizeof(int));
    Framework[CurrentSystem].TorsionDefinitions=(QUAD*)calloc(Framework[CurrentSystem].NumberOfTorsionDefinitions,sizeof(QUAD));
    Framework[CurrentSystem].NumberOfTorsionsPerType=(int*)calloc(Framework[CurrentSystem].NumberOfTorsionDefinitions,sizeof(int));
    Framework[CurrentSystem].TorsionArgumentDefinitions=(REAL(*)[MAX_TORSION_POTENTIAL_ARGUMENTS])
       calloc(Framework[CurrentSystem].NumberOfTorsionDefinitions,sizeof(REAL[MAX_TORSION_POTENTIAL_ARGUMENTS]));


    Framework[CurrentSystem].TorsionDefinitionType[index]=TorsionType;
    Framework[CurrentSystem].TorsionDefinitions[index].A=TypeA;
    Framework[CurrentSystem].TorsionDefinitions[index].B=TypeB;
    Framework[CurrentSystem].TorsionDefinitions[index].C=TypeC;
    Framework[CurrentSystem].TorsionDefinitions[index].D=TypeD;
    Framework[CurrentSystem].NumberOfTorsionsPerType[index]++;
    TorsionTypes[Framework[CurrentFramework].TorsionDefinitionType[index]].nr_args=NumberOfArguments;
    for(j=0;j<NumberOfArguments;j++)
      Framework[CurrentSystem].TorsionArgumentDefinitions[index][j]=parameters[j];
  }
  else
  {
    AlreadyPresent=FALSE;
    for(i=0;i<Framework[CurrentSystem].NumberOfTorsionDefinitions;i++)
    {
      if(((Framework[CurrentSystem].TorsionDefinitions[i].A==TypeA)&&(Framework[CurrentSystem].TorsionDefinitions[i].B==TypeB)&&
         (Framework[CurrentSystem].TorsionDefinitions[i].C==TypeC)&&(Framework[CurrentSystem].TorsionDefinitions[i].D==TypeD))||
         ((Framework[CurrentSystem].TorsionDefinitions[i].A==TypeD)&&(Framework[CurrentSystem].TorsionDefinitions[i].B==TypeC)&&
         (Framework[CurrentSystem].TorsionDefinitions[i].C==TypeB)&&(Framework[CurrentSystem].TorsionDefinitions[i].D==TypeA)))
      {
        AlreadyPresent=TRUE;
        for(j=0;j<NumberOfArguments;j++)
          AlreadyPresent=AlreadyPresent&&(fabs(parameters[j]-Framework[CurrentSystem].TorsionArgumentDefinitions[i][j])<1e-5);;

        if(AlreadyPresent)
        {
          Framework[CurrentSystem].NumberOfTorsionsPerType[i]++;
          break;
        }
      }
    }
    if(!AlreadyPresent)
    {
      index=Framework[CurrentSystem].NumberOfTorsionDefinitions;
      Framework[CurrentSystem].NumberOfTorsionDefinitions++;
      Framework[CurrentSystem].TorsionDefinitionType=(int*)realloc(Framework[CurrentSystem].TorsionDefinitionType,
          Framework[CurrentSystem].NumberOfTorsionDefinitions*sizeof(int));
      Framework[CurrentSystem].TorsionDefinitions=(QUAD*)realloc(Framework[CurrentSystem].TorsionDefinitions,
          Framework[CurrentSystem].NumberOfTorsionDefinitions*sizeof(QUAD));
      Framework[CurrentSystem].NumberOfTorsionsPerType=(int*)realloc(Framework[CurrentSystem].NumberOfTorsionsPerType,
          Framework[CurrentSystem].NumberOfTorsionDefinitions*sizeof(int));
      Framework[CurrentSystem].TorsionArgumentDefinitions=(REAL(*)[MAX_TORSION_POTENTIAL_ARGUMENTS])
          realloc(Framework[CurrentSystem].TorsionArgumentDefinitions,
          Framework[CurrentSystem].NumberOfTorsionDefinitions*sizeof(REAL[MAX_TORSION_POTENTIAL_ARGUMENTS]));

      Framework[CurrentSystem].TorsionDefinitionType[index]=TorsionType;
      Framework[CurrentSystem].TorsionDefinitions[index].A=TypeA;
      Framework[CurrentSystem].TorsionDefinitions[index].B=TypeB;
      Framework[CurrentSystem].TorsionDefinitions[index].C=TypeC;
      Framework[CurrentSystem].TorsionDefinitions[index].D=TypeD;
      Framework[CurrentSystem].NumberOfTorsionsPerType[index]=1;
      TorsionTypes[Framework[CurrentFramework].TorsionDefinitionType[index]].nr_args=NumberOfArguments;
      for(j=0;j<NumberOfArguments;j++)
        Framework[CurrentSystem].TorsionArgumentDefinitions[index][j]=parameters[j];
    }

  }
}

void AddImproperTorsionTypeToDefinitions(int TypeA,int TypeB,int TypeC,int TypeD,int ImproperTorsionType,REAL *parms)
{
  int i,j,index;
  int AlreadyPresent;
  int CompareParameters;
  int NumberOfArguments;
  REAL parameters[10];

  if(strcasecmp(ImproperTorsionTypes[ImproperTorsionType].Name,"TRAPPE_IMPROPER_DIHEDRAL")==0)
  {
    NumberOfArguments=8;
    parameters[0]=parms[0]*ENERGY_TO_KELVIN;
    parameters[1]=parms[1]*ENERGY_TO_KELVIN;
    parameters[2]=parms[2]*ENERGY_TO_KELVIN;
    parameters[3]=parms[3]*ENERGY_TO_KELVIN;
    parameters[4]=parms[4];
    parameters[5]=parms[5];
    parameters[6]=parms[6];
    parameters[7]=parms[7];
  }

  if(Framework[CurrentSystem].NumberOfImproperTorsionDefinitions==0)
  {
    index=0;
    Framework[CurrentSystem].NumberOfImproperTorsionDefinitions=1;
    Framework[CurrentSystem].ImproperTorsionDefinitionType=(int*)calloc(Framework[CurrentSystem].NumberOfImproperTorsionDefinitions,sizeof(int));
    Framework[CurrentSystem].ImproperTorsionDefinitions=(QUAD*)calloc(Framework[CurrentSystem].NumberOfImproperTorsionDefinitions,sizeof(QUAD));
    Framework[CurrentSystem].NumberOfImproperTorsionsPerType=(int*)calloc(Framework[CurrentSystem].NumberOfImproperTorsionDefinitions,sizeof(int));
    Framework[CurrentSystem].ImproperTorsionArgumentDefinitions=(REAL(*)[MAX_IMPROPER_TORSION_POTENTIAL_ARGUMENTS])
       calloc(Framework[CurrentSystem].NumberOfImproperTorsionDefinitions,sizeof(REAL[MAX_IMPROPER_TORSION_POTENTIAL_ARGUMENTS]));


    Framework[CurrentSystem].ImproperTorsionDefinitionType[index]=ImproperTorsionType;
    Framework[CurrentSystem].ImproperTorsionDefinitions[index].A=TypeA;
    Framework[CurrentSystem].ImproperTorsionDefinitions[index].B=TypeB;
    Framework[CurrentSystem].ImproperTorsionDefinitions[index].C=TypeC;
    Framework[CurrentSystem].ImproperTorsionDefinitions[index].D=TypeD;
    Framework[CurrentSystem].NumberOfImproperTorsionsPerType[index]++;
    ImproperTorsionTypes[Framework[CurrentFramework].ImproperTorsionDefinitionType[index]].nr_args=NumberOfArguments;
    for(j=0;j<NumberOfArguments;j++)
      Framework[CurrentSystem].ImproperTorsionArgumentDefinitions[index][j]=parameters[j];
  }
  else
  {
    AlreadyPresent=FALSE;
    for(i=0;i<Framework[CurrentSystem].NumberOfImproperTorsionDefinitions;i++)
    {
      if(((Framework[CurrentSystem].ImproperTorsionDefinitions[i].A==TypeA)&&(Framework[CurrentSystem].ImproperTorsionDefinitions[i].B==TypeB)&&
         (Framework[CurrentSystem].ImproperTorsionDefinitions[i].C==TypeC)&&(Framework[CurrentSystem].ImproperTorsionDefinitions[i].D==TypeD))||
         ((Framework[CurrentSystem].ImproperTorsionDefinitions[i].A==TypeD)&&(Framework[CurrentSystem].ImproperTorsionDefinitions[i].B==TypeC)&&
         (Framework[CurrentSystem].ImproperTorsionDefinitions[i].C==TypeB)&&(Framework[CurrentSystem].ImproperTorsionDefinitions[i].D==TypeA)))
      {
        AlreadyPresent=TRUE;
        for(j=0;j<NumberOfArguments;j++)
          AlreadyPresent=AlreadyPresent&&(fabs(parameters[j]-Framework[CurrentSystem].ImproperTorsionArgumentDefinitions[i][j])<1e-5);;

        if(AlreadyPresent)
        {
          Framework[CurrentSystem].NumberOfImproperTorsionsPerType[i]++;
          break;
        }
      }
    }
    if(!AlreadyPresent)
    {
      index=Framework[CurrentSystem].NumberOfImproperTorsionDefinitions;
      Framework[CurrentSystem].NumberOfImproperTorsionDefinitions++;
      Framework[CurrentSystem].ImproperTorsionDefinitionType=(int*)realloc(Framework[CurrentSystem].ImproperTorsionDefinitionType,
          Framework[CurrentSystem].NumberOfImproperTorsionDefinitions*sizeof(int));
      Framework[CurrentSystem].ImproperTorsionDefinitions=(QUAD*)realloc(Framework[CurrentSystem].ImproperTorsionDefinitions,
          Framework[CurrentSystem].NumberOfImproperTorsionDefinitions*sizeof(QUAD));
      Framework[CurrentSystem].NumberOfImproperTorsionsPerType=(int*)realloc(Framework[CurrentSystem].NumberOfImproperTorsionsPerType,
          Framework[CurrentSystem].NumberOfImproperTorsionDefinitions*sizeof(int));
      Framework[CurrentSystem].ImproperTorsionArgumentDefinitions=(REAL(*)[MAX_IMPROPER_TORSION_POTENTIAL_ARGUMENTS])
          realloc(Framework[CurrentSystem].ImproperTorsionArgumentDefinitions,
          Framework[CurrentSystem].NumberOfImproperTorsionDefinitions*sizeof(REAL[MAX_IMPROPER_TORSION_POTENTIAL_ARGUMENTS]));

      Framework[CurrentSystem].ImproperTorsionDefinitionType[index]=ImproperTorsionType;
      Framework[CurrentSystem].ImproperTorsionDefinitions[index].A=TypeA;
      Framework[CurrentSystem].ImproperTorsionDefinitions[index].B=TypeB;
      Framework[CurrentSystem].ImproperTorsionDefinitions[index].C=TypeC;
      Framework[CurrentSystem].ImproperTorsionDefinitions[index].D=TypeD;
      Framework[CurrentSystem].NumberOfImproperTorsionsPerType[index]=1;
      ImproperTorsionTypes[Framework[CurrentFramework].ImproperTorsionDefinitionType[index]].nr_args=NumberOfArguments;
      for(j=0;j<NumberOfArguments;j++)
        Framework[CurrentSystem].ImproperTorsionArgumentDefinitions[index][j]=parameters[j];
    }
  }
}

/*********************************************************************************************************
 * Name       | ReadFrameworkSpecificDefinitions                                                         *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Read the definition (bond, bend, torsion, etc) of a framework component.                 *
 * Note       | The excluded Van der Waals and Coulombic interactions are automatically computed from    *
 *            | framework exclusion rules (1-2, 1-3, 1-4 omission, bond/bend/torsion omission).          *
 *********************************************************************************************************/

void ReadFrameworkSpecificDefinition(void)
{
  int i,j,k,l,m,n,f1,index_excluded,present,nr_atoms;
  char line[16384],string[1024],buffer[256];
  char keyword[1024],arguments[16384],firstargument[1024],*arg_pointer;
  int A,B,C,D,A1,B1,C1,TypeA,TypeB,TypeC,TypeD,CurrentTypeA,CurrentTypeB,CurrentTypeC,CurrentTypeD,index,index2,ListTypeC,ListTypeD;
  char TypeNameA[256],TypeNameB[256],TypeNameC[256],TypeNameD[256];
  REAL potential_arguments[10];
  int NumberOfAtoms,NumberOfBonds,NumberOfBends,NumberOfTorsions,NumberOfVDW;
  int NumberOfArguments;
  int BondType=0,BendType=0,TorsionType=0,InversionBendType=0,ImproperTorsionType=0;
  int BondBondType=0,BondBendType=0,BendBendType=0,BondTorsionType=0,BendTorsionType=0;
  REAL r,rab,rbc,theta;
  double temp;
  VECTOR Rab,Rbc,dr;
  FILE *FilePtr;
  int int_temp1,int_temp2;
  REAL UnitConverter;

  UnitConverter=KELVIN_TO_ENERGY;

  CurrentFramework=0;

  if(strlen(Framework[CurrentSystem].FrameworkDefinitions)<1) return;

  switch(Framework[CurrentSystem].FlexibleModelInputType)
  {
    case FLEXIBLE_FILE_TYPE_RASPA:
      break;
    case FLEXIBLE_FILE_TYPE_DLPOLY:

      // try "framework.dlpoly", then "Framework.dlpoly"
      sprintf(buffer,"%s.%s","framework","dlpoly");
      if(!(FilePtr=fopen(buffer,"r")))
      {
        sprintf(buffer,"%s.%s","Framework","dlpoly");
        if(!(FilePtr=fopen(buffer,"r")))
        {
          fprintf(stderr, "Error:  file %s does not exist.\n",buffer);
          exit(1);
        }
      }

      while(ReadLine(line,16384,FilePtr))
      {
        // extract first word
        strcpy(keyword,"keyword");
        sscanf(line,"%s%[^\n]",keyword,arguments);
        sscanf(arguments,"%s",firstargument);

        if(strcasecmp("UNITS",keyword)==0)
        {
          if(strcasecmp("kcal",firstargument)==0) UnitConverter=KCAL_PER_MOL_TO_KELVIN*KELVIN_TO_ENERGY;
          if(strcasecmp("eV",firstargument)==0) UnitConverter=EV_TO_KELVIN*KELVIN_TO_ENERGY;
          if(strcasecmp("kJ",firstargument)==0) UnitConverter=KJ_PER_MOL_TO_KELVIN*KELVIN_TO_ENERGY;
          if(strcasecmp("Kelvin",firstargument)==0) UnitConverter=KELVIN_TO_ENERGY;
          if(strcasecmp("internal",firstargument)==0) UnitConverter=1.0;
        }
        if(strcasecmp("ATOMS",keyword)==0)
        {
          sscanf(arguments,"%d",&NumberOfAtoms);
          fprintf(stderr, "Reading atoms: %d\n",NumberOfAtoms);

          for(i=0;i<NumberOfAtoms;i++)
          {
            ReadLine(line,16384,FilePtr);
            sscanf(line,"%*s %*f %lf %*d",&temp);
            if(UseChargesFromCIFFile)
              Framework[CurrentSystem].Atoms[CurrentFramework][i].Charge=temp;
          }
        }
        if(strcasecmp("BONDS",keyword)==0)
        {
          sscanf(arguments,"%d",&NumberOfBonds);
          fprintf(stderr, "Reading bonds: %d\n",NumberOfBonds);

          // realloc memory

          if(Framework[CurrentSystem].MaxNumberOfBonds[CurrentFramework]>0)
          {
            Framework[CurrentSystem].MaxNumberOfBonds[CurrentFramework]+=NumberOfBonds;
            Framework[CurrentSystem].Bonds[CurrentFramework]=(PAIR*)realloc(Framework[CurrentSystem].Bonds[CurrentFramework],
                     Framework[CurrentSystem].MaxNumberOfBonds[CurrentFramework]*sizeof(PAIR));
            Framework[CurrentSystem].BondType[CurrentFramework]=(int*)realloc(Framework[CurrentSystem].BondType[CurrentFramework],
                     Framework[CurrentSystem].MaxNumberOfBonds[CurrentFramework]*sizeof(int));
            Framework[CurrentSystem].BondArguments[CurrentFramework]=(REAL(*)[MAX_BOND_POTENTIAL_ARGUMENTS])
                    realloc(Framework[CurrentSystem].BondArguments[CurrentFramework],
                            Framework[CurrentSystem].MaxNumberOfBonds[CurrentFramework]*sizeof(REAL[MAX_BOND_POTENTIAL_ARGUMENTS]));
          }
          else
          {
            Framework[CurrentSystem].MaxNumberOfBonds[CurrentFramework]=NumberOfBonds;
            Framework[CurrentSystem].Bonds[CurrentFramework]=(PAIR*)calloc(Framework[CurrentSystem].MaxNumberOfBonds[CurrentFramework],sizeof(PAIR));
            Framework[CurrentSystem].BondType[CurrentFramework]=(int*)calloc(Framework[CurrentSystem].MaxNumberOfBonds[CurrentFramework],sizeof(int));
            Framework[CurrentSystem].BondArguments[CurrentFramework]=(REAL(*)[MAX_BOND_POTENTIAL_ARGUMENTS])
                    calloc(Framework[CurrentSystem].MaxNumberOfBonds[CurrentFramework],sizeof(REAL[MAX_BOND_POTENTIAL_ARGUMENTS]));
          }


          index=Framework[CurrentSystem].NumberOfBonds[CurrentFramework];
          for(i=index;i<index+NumberOfBonds;i++)
          {
            ReadLine(line,16384,FilePtr);
            sscanf(line,"%s %d %d %[^\n]",buffer,&A,&B,arguments);
            Framework[CurrentSystem].Bonds[CurrentFramework][i].A=A-1;
            Framework[CurrentSystem].Bonds[CurrentFramework][i].B=B-1;
            if(strcasecmp(buffer,"mors")==0)
            {
              for(j=0;j<NR_BOND_TYPES;j++)
                if(strcasecmp(BondTypes[j].Name,"MORSE_BOND")==0)
                  BondType=j;
              Framework[CurrentSystem].BondType[CurrentFramework][i]=BondType;
              sscanf(arguments,"%lf %lf %lf\n",&potential_arguments[0],&potential_arguments[1],&potential_arguments[2]);
              Framework[CurrentSystem].BondArguments[CurrentFramework][i][0]=potential_arguments[0]*UnitConverter;
              Framework[CurrentSystem].BondArguments[CurrentFramework][i][1]=potential_arguments[2];
              Framework[CurrentSystem].BondArguments[CurrentFramework][i][2]=potential_arguments[1];
              NumberOfArguments=3;
            }
            else if(strcasecmp(buffer,"harm")==0)
            {
              for(j=0;j<NR_BOND_TYPES;j++)
                if(strcasecmp(BondTypes[j].Name,"HARMONIC_BOND")==0)
                  BondType=j;
              Framework[CurrentSystem].BondType[CurrentFramework][i]=BondType;
              sscanf(arguments,"%lf %lf %lf\n",&potential_arguments[0],&potential_arguments[1],&potential_arguments[2]);
              Framework[CurrentSystem].BondArguments[CurrentFramework][i][0]=potential_arguments[0]*UnitConverter;
              Framework[CurrentSystem].BondArguments[CurrentFramework][i][1]=potential_arguments[1];
              NumberOfArguments=2;
            }

            A=Framework[CurrentSystem].Bonds[CurrentFramework][i].A;
            B=Framework[CurrentSystem].Bonds[CurrentFramework][i].B;
            TypeA=Framework[CurrentSystem].Atoms[CurrentFramework][A].Type;
            TypeB=Framework[CurrentSystem].Atoms[CurrentFramework][B].Type;
            AddBondTypeToDefinitions(TypeA,TypeB,BondType,Framework[CurrentSystem].BondArguments[CurrentFramework][i]);
          }

          Framework[CurrentSystem].NumberOfBonds[CurrentFramework]+=NumberOfBonds;
        }

        if(strcasecmp("ANGLES",keyword)==0)
        {
          sscanf(arguments,"%d",&NumberOfBends);
          fprintf(stderr, "Reading bends: %d\n",NumberOfBends);

          if(Framework[CurrentSystem].MaxNumberOfBonds[CurrentFramework]>0)
          {
            Framework[CurrentSystem].MaxNumberOfBends[CurrentFramework]+=NumberOfBends;
            Framework[CurrentSystem].Bends[CurrentFramework]=(QUAD*)realloc(Framework[CurrentSystem].Bends[CurrentFramework],
                             Framework[CurrentSystem].MaxNumberOfBends[CurrentFramework]*sizeof(QUAD));
            Framework[CurrentSystem].BendType[CurrentFramework]=(int*)realloc(Framework[CurrentSystem].BendType[CurrentFramework],
                             Framework[CurrentSystem].MaxNumberOfBends[CurrentFramework]*sizeof(int));
            Framework[CurrentSystem].BendArguments[CurrentFramework]=(REAL(*)[MAX_BEND_POTENTIAL_ARGUMENTS])
                    realloc(Framework[CurrentSystem].BendArguments[CurrentFramework],
                            Framework[CurrentSystem].MaxNumberOfBends[CurrentFramework]*sizeof(REAL[MAX_BEND_POTENTIAL_ARGUMENTS]));
          }
          else
          {
            Framework[CurrentSystem].MaxNumberOfBends[CurrentFramework]=NumberOfBends;
            Framework[CurrentSystem].Bends[CurrentFramework]=(QUAD*)calloc(Framework[CurrentSystem].MaxNumberOfBends[CurrentFramework],sizeof(QUAD));
            Framework[CurrentSystem].BendType[CurrentFramework]=(int*)calloc(Framework[CurrentSystem].MaxNumberOfBends[CurrentFramework],sizeof(int));
            Framework[CurrentSystem].BendArguments[CurrentFramework]=(REAL(*)[MAX_BEND_POTENTIAL_ARGUMENTS])
                    calloc(Framework[CurrentSystem].MaxNumberOfBends[CurrentFramework],sizeof(REAL[MAX_BEND_POTENTIAL_ARGUMENTS]));
          }

          index=Framework[CurrentSystem].NumberOfBends[CurrentFramework];

          for(i=index;i<index+NumberOfBends;i++)
          {
            ReadLine(line,16384,FilePtr);
            sscanf(line,"%s %d %d %d %[^\n]",buffer,&A,&B,&C,arguments);
            Framework[CurrentSystem].Bends[CurrentFramework][i].A=A-1;
            Framework[CurrentSystem].Bends[CurrentFramework][i].B=B-1;
            Framework[CurrentSystem].Bends[CurrentFramework][i].C=C-1;
            Framework[CurrentSystem].Bends[CurrentFramework][i].D=0;
            if(strcasecmp(buffer,"harm")==0)
            {
              for(j=0;j<NR_BEND_TYPES;j++)
                if(strcasecmp(BendTypes[j].Name,"HARMONIC_BEND")==0)
                  BendType=j;
              Framework[CurrentSystem].BendType[CurrentFramework][i]=BendType;
              sscanf(arguments,"%lf %lf %lf\n",&potential_arguments[0],&potential_arguments[1],&potential_arguments[2]);
              Framework[CurrentSystem].BendArguments[CurrentFramework][i][0]=potential_arguments[0]*UnitConverter;
              Framework[CurrentSystem].BendArguments[CurrentFramework][i][1]=potential_arguments[1]*DEG2RAD;
            }
            A=Framework[CurrentSystem].Bends[CurrentFramework][i].A;
            B=Framework[CurrentSystem].Bends[CurrentFramework][i].B;
            C=Framework[CurrentSystem].Bends[CurrentFramework][i].C;
            TypeA=Framework[CurrentSystem].Atoms[CurrentFramework][A].Type;
            TypeB=Framework[CurrentSystem].Atoms[CurrentFramework][B].Type;
            TypeC=Framework[CurrentSystem].Atoms[CurrentFramework][C].Type;
            AddBendTypeToDefinitions(TypeA,TypeB,TypeC,BendType,Framework[CurrentSystem].BendArguments[CurrentFramework][i]);
          }

          Framework[CurrentSystem].NumberOfBends[CurrentFramework]+=NumberOfBends;
        }

        if(strcasecmp("DIHEDRALS",keyword)==0)
        {
          sscanf(arguments,"%d",&NumberOfTorsions);
          fprintf(stderr, "Reading torsions: %d\n",NumberOfTorsions);

          if(Framework[CurrentSystem].MaxNumberOfTorsions[CurrentFramework]>0)
          {
            Framework[CurrentSystem].MaxNumberOfTorsions[CurrentFramework]+=NumberOfTorsions;
            Framework[CurrentSystem].Torsions[CurrentFramework]=(QUAD*)realloc(Framework[CurrentSystem].Torsions[CurrentFramework],
                   Framework[CurrentSystem].MaxNumberOfTorsions[CurrentFramework]*sizeof(QUAD));
            Framework[CurrentSystem].TorsionType[CurrentFramework]=(int*)realloc(Framework[CurrentSystem].TorsionType[CurrentFramework],
                   Framework[CurrentSystem].MaxNumberOfTorsions[CurrentFramework]*sizeof(int));
            Framework[CurrentSystem].TorsionArguments[CurrentFramework]=(REAL(*)[MAX_TORSION_POTENTIAL_ARGUMENTS])
                    realloc(Framework[CurrentSystem].TorsionArguments[CurrentFramework],
                    Framework[CurrentSystem].MaxNumberOfTorsions[CurrentFramework]*sizeof(REAL[MAX_TORSION_POTENTIAL_ARGUMENTS]));
          }
          else
          {
            Framework[CurrentSystem].MaxNumberOfTorsions[CurrentFramework]=NumberOfTorsions;
            Framework[CurrentSystem].Torsions[CurrentFramework]=(QUAD*)calloc(Framework[CurrentSystem].MaxNumberOfTorsions[CurrentFramework],sizeof(QUAD));
            Framework[CurrentSystem].TorsionType[CurrentFramework]=(int*)calloc(Framework[CurrentSystem].MaxNumberOfTorsions[CurrentFramework],sizeof(int));
            Framework[CurrentSystem].TorsionArguments[CurrentFramework]=(REAL(*)[MAX_TORSION_POTENTIAL_ARGUMENTS])
                    calloc(Framework[CurrentSystem].MaxNumberOfTorsions[CurrentFramework],sizeof(REAL[MAX_TORSION_POTENTIAL_ARGUMENTS]));
          }

          if(Framework[CurrentSystem].MaxNumberOfImproperTorsions[CurrentFramework]>0)
          {
            Framework[CurrentSystem].MaxNumberOfImproperTorsions[CurrentFramework]+=NumberOfTorsions;
            Framework[CurrentSystem].ImproperTorsions[CurrentFramework]=(QUAD*)realloc(Framework[CurrentSystem].ImproperTorsions[CurrentFramework],
                   Framework[CurrentSystem].MaxNumberOfImproperTorsions[CurrentFramework]*sizeof(QUAD));
            Framework[CurrentSystem].ImproperTorsionType[CurrentFramework]=(int*)realloc(Framework[CurrentSystem].ImproperTorsionType[CurrentFramework],
                   Framework[CurrentSystem].MaxNumberOfImproperTorsions[CurrentFramework]*sizeof(int));
            Framework[CurrentSystem].ImproperTorsionArguments[CurrentFramework]=(REAL(*)[MAX_TORSION_POTENTIAL_ARGUMENTS])
                    realloc(Framework[CurrentSystem].ImproperTorsionArguments[CurrentFramework],
                    Framework[CurrentSystem].MaxNumberOfImproperTorsions[CurrentFramework]*sizeof(REAL[MAX_TORSION_POTENTIAL_ARGUMENTS]));
          }
          else
          {
            Framework[CurrentSystem].MaxNumberOfImproperTorsions[CurrentFramework]=NumberOfTorsions;
            Framework[CurrentSystem].ImproperTorsions[CurrentFramework]=(QUAD*)calloc(Framework[CurrentSystem].MaxNumberOfImproperTorsions[CurrentFramework],sizeof(QUAD));
            Framework[CurrentSystem].ImproperTorsionType[CurrentFramework]=(int*)calloc(Framework[CurrentSystem].MaxNumberOfImproperTorsions[CurrentFramework],sizeof(int));
            Framework[CurrentSystem].ImproperTorsionArguments[CurrentFramework]=(REAL(*)[MAX_TORSION_POTENTIAL_ARGUMENTS])
                    calloc(Framework[CurrentSystem].MaxNumberOfImproperTorsions[CurrentFramework],sizeof(REAL[MAX_TORSION_POTENTIAL_ARGUMENTS]));
          }

          index=Framework[CurrentSystem].NumberOfTorsions[CurrentFramework];
          index2=Framework[CurrentSystem].NumberOfImproperTorsions[CurrentFramework];

          for(i=0;i<NumberOfTorsions;i++)
          {
            ReadLine(line,16384,FilePtr);
            sscanf(line,"%s %d %d %d %d %[^\n]",buffer,&A,&B,&C,&D,arguments);

            if(((BondNeighbours(A-1,B-1)&&BondNeighbours(B-1,C-1)&&BondNeighbours(C-1,D-1)&&(!BondNeighbours(A-1,D-1)))))
            {
              Framework[CurrentSystem].Torsions[CurrentFramework][index].A=A-1;
              Framework[CurrentSystem].Torsions[CurrentFramework][index].B=B-1;
              Framework[CurrentSystem].Torsions[CurrentFramework][index].C=C-1;
              Framework[CurrentSystem].Torsions[CurrentFramework][index].D=D-1;

              if(strcasecmp(buffer,"cos")==0)
              {
                for(j=0;j<NR_TORSION_TYPES;j++)
                  if(strcasecmp(TorsionTypes[j].Name,"TRAPPE_DIHEDRAL")==0)
                    TorsionType=j;
                Framework[CurrentSystem].TorsionType[CurrentFramework][index]=TorsionType;
                for(j=0;j<MAX_TORSION_POTENTIAL_ARGUMENTS;j++)
                  potential_arguments[j]=0.0;
                sscanf(arguments,"%lf %lf %d %lf %lf\n",&potential_arguments[0],&potential_arguments[1],&int_temp2,&potential_arguments[2],&potential_arguments[3]);
                Framework[CurrentSystem].TorsionArguments[CurrentFramework][index][0]=0.0;
                Framework[CurrentSystem].TorsionArguments[CurrentFramework][index][1]=0.0;
                Framework[CurrentSystem].TorsionArguments[CurrentFramework][index][2]=0.0;
                Framework[CurrentSystem].TorsionArguments[CurrentFramework][index][3]=0.0;

                // settings scaling parameters for VDW (second parameter in DLPOLY)
                Framework[CurrentSystem].TorsionArguments[CurrentFramework][index][6]=potential_arguments[3];

                // settings scaling parameters for electrostatics (first parameter in DLPOLY)
                Framework[CurrentSystem].TorsionArguments[CurrentFramework][index][7]=potential_arguments[2];

                switch(int_temp2)
                {
                  case 2:
                    Framework[CurrentSystem].TorsionArguments[CurrentFramework][index][2]=potential_arguments[0]*UnitConverter;
                    break;
                  case 3:
                    Framework[CurrentSystem].TorsionArguments[CurrentFramework][index][3]=potential_arguments[0]*UnitConverter;
                    break;
                }
              }
              A=Framework[CurrentSystem].Torsions[CurrentFramework][index].A;
              B=Framework[CurrentSystem].Torsions[CurrentFramework][index].B;
              C=Framework[CurrentSystem].Torsions[CurrentFramework][index].C;
              D=Framework[CurrentSystem].Torsions[CurrentFramework][index].D;
              TypeA=Framework[CurrentSystem].Atoms[CurrentFramework][A].Type;
              TypeB=Framework[CurrentSystem].Atoms[CurrentFramework][B].Type;
              TypeC=Framework[CurrentSystem].Atoms[CurrentFramework][C].Type;
              TypeD=Framework[CurrentSystem].Atoms[CurrentFramework][D].Type;
              AddTorsionTypeToDefinitions(TypeA,TypeB,TypeC,TypeD,TorsionType,Framework[CurrentSystem].TorsionArguments[CurrentFramework][index]);
              index++;
            }
            else if(((BondNeighbours(A-1,C-1)&&BondNeighbours(B-1,C-1)&&BondNeighbours(D-1,C-1))))
            {
              Framework[CurrentSystem].ImproperTorsions[CurrentFramework][index2].A=A-1;
              Framework[CurrentSystem].ImproperTorsions[CurrentFramework][index2].B=B-1;
              Framework[CurrentSystem].ImproperTorsions[CurrentFramework][index2].C=C-1;
              Framework[CurrentSystem].ImproperTorsions[CurrentFramework][index2].D=D-1;

              if(strcasecmp(buffer,"cos")==0)
              {
                for(j=0;j<NR_IMPROPER_TORSION_TYPES;j++)
                  if(strcasecmp(ImproperTorsionTypes[j].Name,"TRAPPE_IMPROPER_DIHEDRAL")==0)
                    ImproperTorsionType=j;
                Framework[CurrentSystem].ImproperTorsionType[CurrentFramework][index2]=ImproperTorsionType;
                sscanf(arguments,"%lf %lf %d %lf %lf\n",&potential_arguments[0],&potential_arguments[1],&int_temp2,&potential_arguments[2],&potential_arguments[3]);
                Framework[CurrentSystem].ImproperTorsionArguments[CurrentFramework][index2][0]=0.0;
                Framework[CurrentSystem].ImproperTorsionArguments[CurrentFramework][index2][1]=0.0;
                Framework[CurrentSystem].ImproperTorsionArguments[CurrentFramework][index2][2]=0.0;
                Framework[CurrentSystem].ImproperTorsionArguments[CurrentFramework][index2][3]=0.0;

                // settings scaling parameters for VDW (second parameter in DLPOLY)
                Framework[CurrentSystem].ImproperTorsionArguments[CurrentFramework][index2][6]=potential_arguments[3];

                // settings scaling parameters for electrostatics (first parameter in DLPOLY)
                Framework[CurrentSystem].ImproperTorsionArguments[CurrentFramework][index2][7]=potential_arguments[2];

                switch(int_temp2)
                {
                  case 2:
                    Framework[CurrentSystem].ImproperTorsionArguments[CurrentFramework][index2][2]=potential_arguments[0]*UnitConverter;
                    break;
                  case 3:
                    Framework[CurrentSystem].ImproperTorsionArguments[CurrentFramework][index2][3]=potential_arguments[0]*UnitConverter;
                    break;
                }
              }
              A=Framework[CurrentSystem].ImproperTorsions[CurrentFramework][index2].A;
              B=Framework[CurrentSystem].ImproperTorsions[CurrentFramework][index2].B;
              C=Framework[CurrentSystem].ImproperTorsions[CurrentFramework][index2].C;
              D=Framework[CurrentSystem].ImproperTorsions[CurrentFramework][index2].D;
              TypeA=Framework[CurrentSystem].Atoms[CurrentFramework][A].Type;
              TypeB=Framework[CurrentSystem].Atoms[CurrentFramework][B].Type;
              TypeC=Framework[CurrentSystem].Atoms[CurrentFramework][C].Type;
              TypeD=Framework[CurrentSystem].Atoms[CurrentFramework][D].Type;
              AddImproperTorsionTypeToDefinitions(TypeA,TypeB,TypeC,TypeD,ImproperTorsionType,Framework[CurrentSystem].ImproperTorsionArguments[CurrentFramework][index2]);
              index2++;
            }
          }
          Framework[CurrentSystem].NumberOfTorsions[CurrentFramework]=index;
          Framework[CurrentSystem].NumberOfImproperTorsions[CurrentFramework]=index2;
        }

        if(strcasecmp("VDW",keyword)==0)
        {
          sscanf(arguments,"%d",&NumberOfVDW);
          fprintf(stderr, "Reading VDW: %d\n",NumberOfVDW);
          for(i=0;i<NumberOfVDW;i++)
          {
            ReadLine(line,16384,FilePtr);
            sscanf(line,"%s %s %s %lf %lf %[^\n]",TypeNameA,TypeNameB,buffer,&potential_arguments[0],&potential_arguments[1],arguments);
            if(strcasecmp(buffer,"lj")==0)
            {
            }
          }
        }
      }

      break;
  }
}

/*********************************************************************************************************
 * Name       | WriteFrameworkDefinitionDLPOLY                                                           *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Write the definition (bond, bend, torsion, etc) of a framework component.                *
 * Note       | The excluded Van der Waals and Coulombic interactions are automatically computed from    *
 *            | framework exclusion rules (1-2, 1-3, 1-4 omission, bond/bend/torsion omission).          *
 *********************************************************************************************************/

void WriteFrameworkDefinitionDLPOLY(void)
{
  int i,j,f1;
  FILE *FilePtr;
  int A,B,C,D;
  int typeA,typeB,typeC,typeD;
  int total;
  REAL charge,mass;
  REAL *parms;
  VECTOR pos,vel,force;
  int index;
  INT_VECTOR3 fixed;
  int m,k,l,Type;

  // Don't make files and folders when streaming, dammit!
  if (STREAM)
    return;

  mkdir("DLPOLY",S_IRWXU);

  FilePtr=fopen("DLPOLY/CONTROL","w");
  fprintf(FilePtr,"%s\n",Framework[CurrentSystem].Name[0]);
  fprintf(FilePtr,"\n");

  fprintf(FilePtr,"temperature %12.6lf\n",(double)therm_baro_stats.ExternalTemperature[CurrentSystem]);
  fprintf(FilePtr,"pressure    %12.6lf\n",(double)therm_baro_stats.ExternalPressure[CurrentSystem][0]*PRESSURE_CONVERSION_FACTOR*PA_TO_ATM*0.001);  // k atm
  fprintf(FilePtr,"\n");

  switch(Ensemble[CurrentSystem])
  {
    case NPT:
      // contribution is zero (i.e. cancels out)
      break;
    case NPTPR:
    case NPHPR:
      fprintf(FilePtr,"ensemble nst lang    1.0  5.0\n");
      break;
  }

  fprintf(FilePtr,"vdw direct\n");
  fprintf(FilePtr,"spme precision  1.0E-20\n");
  fprintf(FilePtr,"\n");

  fprintf(FilePtr,"steps       %12lld\n",NumberOfCycles);
  fprintf(FilePtr,"timestep    %12.6f ps\n",DeltaT);
  fprintf(FilePtr,"multiple    1 steps\n");
  fprintf(FilePtr,"\n");

  fprintf(FilePtr,"cutoff   %12.6lf  angstrom\n",CutOffVDW);
  fprintf(FilePtr,"rvdw     %12.6lf  angstrom\n",CutOffVDW);
  fprintf(FilePtr,"\n");

  fprintf(FilePtr,"job time     7200000000 seconds\n");
  fprintf(FilePtr,"close time   100000 seconds\n");
  fprintf(FilePtr,"\n");

  fprintf(FilePtr,"finish\n");
  fprintf(FilePtr,"\n");

  fclose(FilePtr);

  FilePtr=fopen("DLPOLY/CONFIG","w");
  fprintf(FilePtr,"Raspa-1.0: input files for dlpoly\n");
  fprintf(FilePtr,"%10d%10d\n",2,3);
  fprintf(FilePtr,"%20.10f%20.10f%20.10f\n",Box[0].ax,Box[0].ay,Box[0].az);
  fprintf(FilePtr,"%20.10f%20.10f%20.10f\n",Box[0].bx,Box[0].by,Box[0].bz);
  fprintf(FilePtr,"%20.10f%20.10f%20.10f\n",Box[0].cx,Box[0].cy,Box[0].cz);

  index=1;
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
     for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
     {
       typeA=Framework[CurrentSystem].Atoms[f1][i].Type;
       pos=Framework[CurrentSystem].Atoms[f1][i].Position;
       vel=Framework[CurrentSystem].Atoms[f1][i].Velocity;
       force=Framework[CurrentSystem].Atoms[f1][i].Force;
       fprintf(FilePtr,"%-8s%10d\n",PseudoAtoms[typeA].Name,index);
       fprintf(FilePtr,"%20.10f%20.10f%20.10f\n",pos.x,pos.y,pos.z);
       fprintf(FilePtr,"%20.10f%20.10f%20.10f\n",vel.x,vel.y,vel.z);
       fprintf(FilePtr,"%20.10f%20.10f%20.10f\n",force.x,force.y,force.z);
       index++;
     }
  }
  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    Type=Adsorbates[CurrentSystem][m].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      for(k=0;k<Components[Type].Groups[l].NumberOfGroupAtoms;k++)
      {
        A=Components[Type].Groups[l].Atoms[k];
        typeA=Adsorbates[CurrentSystem][m].Atoms[A].Type;
        pos=Adsorbates[CurrentSystem][m].Atoms[A].Position;
        vel=Adsorbates[CurrentSystem][m].Atoms[A].Velocity;
        force=Adsorbates[CurrentSystem][m].Atoms[A].Force;
        fprintf(FilePtr,"%8s%10d\n",PseudoAtoms[typeA].Name,index);
        fprintf(FilePtr,"%20.10lf%20.10lf%20.10lf\n",pos.x,pos.y,pos.z);
        fprintf(FilePtr,"%20.10lf%20.10lf%20.10lf\n",vel.x,vel.y,vel.z);
        fprintf(FilePtr,"%20.10lf%20.10lf%20.10lf\n",force.x,force.y,force.z);
        index++;
      }
    }
  }
  for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
  {
    Type=Cations[CurrentSystem][m].Type;
    for(l=0;l<Components[Type].NumberOfGroups;l++)
    {
      for(k=0;k<Components[Type].Groups[l].NumberOfGroupAtoms;k++)
      {
        A=Components[Type].Groups[l].Atoms[k];
        typeA=Cations[CurrentSystem][m].Atoms[A].Type;
        pos=Cations[CurrentSystem][m].Atoms[A].Position;
        vel=Cations[CurrentSystem][m].Atoms[A].Velocity;
        force=Cations[CurrentSystem][m].Atoms[A].Force;
        fprintf(FilePtr,"%8s%10d\n",PseudoAtoms[typeA].Name,index);
        fprintf(FilePtr,"%20.10lf%20.10lf%20.10lf\n",pos.x,pos.y,pos.z);
        fprintf(FilePtr,"%20.10lf%20.10lf%20.10lf\n",vel.x,vel.y,vel.z);
        fprintf(FilePtr,"%20.10lf%20.10lf%20.10lf\n",force.x,force.y,force.z);
        index++;
      }
    }
  }
  fclose(FilePtr);



  FilePtr=fopen("DLPOLY/FIELD","w");
  fprintf(FilePtr,"%s\n",Framework[CurrentSystem].Name[0]);
  fprintf(FilePtr,"UNITS Kelvin\n");
  fprintf(FilePtr,"MOLECULES %d\n",1+NumberOfComponents);
  fprintf(FilePtr,"%s\n",Framework[CurrentSystem].Name[0]);
  fprintf(FilePtr,"NUMMOL 1\n");


  fprintf(FilePtr,"ATOMS %d\n",Framework[CurrentSystem].TotalNumberOfAtoms);
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
     for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[f1];i++)
     {
       typeA=Framework[CurrentSystem].Atoms[f1][i].Type;
       mass=PseudoAtoms[typeA].Mass;
       charge=Framework[CurrentSystem].Atoms[f1][i].Charge;
       fixed=Framework[CurrentSystem].Atoms[f1][i].Fixed;
       fprintf(FilePtr,"%-6s %12.6f %12.6f   1    %5d\n",PseudoAtoms[typeA].Name,mass,charge,(fixed.x&&fixed.y&&fixed.z)?1:0);
     }
  }

  total=0;
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
    total+=Framework[CurrentSystem].NumberOfBonds[f1];

  fprintf(FilePtr,"BONDS %d\n",total);
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfBonds[f1];i++)
    {
      A=Framework[CurrentSystem].Bonds[f1][i].A;
      B=Framework[CurrentSystem].Bonds[f1][i].B;
      typeA=Framework[CurrentSystem].Atoms[f1][A].Type;
      typeB=Framework[CurrentSystem].Atoms[f1][B].Type;

      parms=(REAL*)&Framework[CurrentSystem].BondArguments[f1][i];

      switch(Framework[CurrentSystem].BondType[f1][i])
      {
        case HARMONIC_BOND:
          // 0.5*p0*SQR(r-p1);
          // ===============================================
          // p_0/k_B [K/A^2]   force constant
          // p_1     [A]       reference bond distance
          fprintf(FilePtr,"harm        %5d %5d  %12.6f %12.6f\n",
            A+1,B+1,parms[0]*ENERGY_TO_KELVIN,parms[1]);
          break;
        case MORSE_BOND:
          // p_0*[(1.0-{exp(-p_1*(r-p_2))})^2-1.0]
          // ===============================================
          // p_0/k_B [K]       force constant
          // p_1     [A^-1]    parameter
          // p_2     [A]       reference bond distance
          fprintf(FilePtr,"morse       %5d %5d  %12.6f %12.6f %12.6f\n",
             A+1,B+1,parms[0]*ENERGY_TO_KELVIN,parms[2],parms[1]);
          break;
      }
    }
  }

  total=0;
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
    total+=Framework[CurrentSystem].NumberOfBends[f1];

  fprintf(FilePtr,"ANGLES %d\n",total);
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfBends[f1];i++)
    {
      A=Framework[CurrentSystem].Bends[f1][i].A;
      B=Framework[CurrentSystem].Bends[f1][i].B;
      C=Framework[CurrentSystem].Bends[f1][i].C;
      typeA=Framework[CurrentSystem].Atoms[f1][A].Type;
      typeB=Framework[CurrentSystem].Atoms[f1][B].Type;
      typeC=Framework[CurrentSystem].Atoms[f1][C].Type;

      parms=(REAL*)&Framework[CurrentSystem].BendArguments[f1][i];

      switch(Framework[CurrentSystem].BendType[f1][i])
      {
         case HARMONIC_BEND:
            // (1/2)p_0*(theta-p_1)^2
            // ===============================================
            // p_0/k_B [K/rad^2]
            // p_1     [degrees]
          fprintf(FilePtr,"harm        %5d %5d %5d  %12.6f %12.6f\n",
            A+1,B+1,C+1,parms[0]*ENERGY_TO_KELVIN,parms[1]*RAD2DEG);
          break;

      }
    }
  }

  total=0;
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    total+=Framework[CurrentSystem].NumberOfTorsions[f1];
    total+=Framework[CurrentSystem].NumberOfImproperTorsions[f1];
  }

  fprintf(FilePtr,"DIHEDRALS %d\n",total);
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfTorsions[f1];i++)
    {
      A=Framework[CurrentSystem].Torsions[f1][i].A;
      B=Framework[CurrentSystem].Torsions[f1][i].B;
      C=Framework[CurrentSystem].Torsions[f1][i].C;
      D=Framework[CurrentSystem].Torsions[f1][i].D;
      typeA=Framework[CurrentSystem].Atoms[f1][A].Type;
      typeB=Framework[CurrentSystem].Atoms[f1][B].Type;
      typeC=Framework[CurrentSystem].Atoms[f1][C].Type;
      typeD=Framework[CurrentSystem].Atoms[f1][D].Type;

      parms=(REAL*)&Framework[CurrentSystem].TorsionArguments[f1][i];

      switch(Framework[CurrentSystem].TorsionType[f1][i])
      {
        case TRAPPE_DIHEDRAL:
          // p_0[0]+p_1*(1+cos(phi))+p_2*(1-cos(2*phi))+p_3*(1+cos(3*phi))
          // =============================================================
          // p_0/k_B [K]
          // p_1/k_B [K]
          // p_2/k_B [K]
          // p_3/k_B [K]
          fprintf(FilePtr,"cos3        %5d %5d %5d %5d  %12.6f %12.6f %12.6f %12.6f %12.6f\n",
            A+1,B+1,C+1,D+1,
            2.0*parms[1]*ENERGY_TO_KELVIN,
            2.0*parms[2]*ENERGY_TO_KELVIN,
            2.0*parms[3]*ENERGY_TO_KELVIN,
            parms[7],
            parms[6]);
          break;

      }
    }
  }
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfImproperTorsions[f1];i++)
    {
      A=Framework[CurrentSystem].ImproperTorsions[f1][i].A;
      B=Framework[CurrentSystem].ImproperTorsions[f1][i].B;
      C=Framework[CurrentSystem].ImproperTorsions[f1][i].C;
      D=Framework[CurrentSystem].ImproperTorsions[f1][i].D;
      typeA=Framework[CurrentSystem].Atoms[f1][A].Type;
      typeB=Framework[CurrentSystem].Atoms[f1][B].Type;
      typeC=Framework[CurrentSystem].Atoms[f1][C].Type;
      typeD=Framework[CurrentSystem].Atoms[f1][D].Type;

      parms=(REAL*)&Framework[CurrentSystem].ImproperTorsionArguments[f1][i];

      switch(Framework[CurrentSystem].ImproperTorsionType[f1][i])
      {
        case TRAPPE_IMPROPER_DIHEDRAL:
          // p_0[0]+p_1*(1+cos(phi))+p_2*(1-cos(2*phi))+p_3*(1+cos(3*phi))
          // =============================================================
          // p_0/k_B [K]
          // p_1/k_B [K]
          // p_2/k_B [K]
          // p_3/k_B [K]
          fprintf(FilePtr,"cos3        %5d %5d %5d %5d  %12.6f %12.6f %12.6f %12.6f %12.6f\n",
            A+1,B+1,C+1,D+1,
            2.0*parms[1]*ENERGY_TO_KELVIN,
            2.0*parms[2]*ENERGY_TO_KELVIN,
            2.0*parms[3]*ENERGY_TO_KELVIN,
            parms[7],
            parms[6]);
          break;

      }
    }
  }
  fprintf(FilePtr,"FINISH\n");

  for(i=0;i<NumberOfComponents;i++)
  {
    if(Components[i].NumberOfMolecules[CurrentSystem]>0)
    {

      total=0;
      for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
      {
        Type=Adsorbates[CurrentSystem][m].Type;
        if(Type==i)
        {
          for(l=0;l<Components[Type].NumberOfGroups;l++)
          {
            if(Components[Type].Groups[l].Rigid) // rigid unit
            {
              total+=Components[Type].Groups[l].NumberOfGroupAtoms;
            }
          }
        }
      }

      fprintf(FilePtr,"%s\n",Components[i].Name);
      fprintf(FilePtr,"NUMMOL %d\n",Components[i].NumberOfMolecules[CurrentSystem]);
      fprintf(FilePtr,"ATOMS %d\n",Components[i].NumberOfAtoms);
      for(l=0;l<Components[i].NumberOfGroups;l++)
      {
        for(k=0;k<Components[i].Groups[l].NumberOfGroupAtoms;k++)
        {
          A=Components[i].Groups[l].Atoms[k];
          typeA=Components[i].Type[A];
          mass=PseudoAtoms[typeA].Mass;
          charge=PseudoAtoms[typeA].Charge1;
          fprintf(FilePtr,"%-6s %12.6f %12.6f 1   0\n",PseudoAtoms[typeA].Name,mass,charge);
        }
      }

      // print internal bonds
      fprintf(FilePtr,"BONDS %d\n",Components[i].NumberOfBonds);
      for(k=0;k<Components[i].NumberOfBonds;k++)
      {
        A=Components[i].Bonds[k].A;
        B=Components[i].Bonds[k].B;

        parms=Components[i].BondArguments[k];
        switch(Components[i].BondType[k])
        {
          case RIGID_BOND:
            break;
          case HARMONIC_BOND:
            // 0.5*p0*SQR(r-p1);
            // ===============================================
            // p_0/k_B [K/A^2]   force constant
            // p_1     [A]       reference bond distance
            fprintf(FilePtr,"harm        %5d %5d  %12.6f %12.6f\n",
              A+1,B+1,parms[0]*ENERGY_TO_KELVIN,parms[1]);
            break;
          case MORSE_BOND:
            // p_0*[(1.0-{exp(-p_1*(r-p_2))})^2-1.0]
            // ===============================================
            // p_0/k_B [K]       force constant
            // p_1     [A^-1]    parameter
            // p_2     [A]       reference bond distance
            fprintf(FilePtr,"morse       %5d %5d  %12.6f %12.6f %12.6f\n",
               A+1,B+1,parms[0]*ENERGY_TO_KELVIN,parms[2],parms[1]);
            break;
        }
      }

      // print internal bends
      fprintf(FilePtr,"ANGLES %d\n",Components[i].NumberOfBends);
      for(k=0;k<Components[i].NumberOfBends;k++)
      {
        A=Components[i].Bends[k].A;
        B=Components[i].Bends[k].B;
        C=Components[i].Bends[k].C;

        parms=Components[i].BendArguments[k];
        switch(Components[i].BendType[k])
        {
          case HARMONIC_BEND:
             // (1/2)p_0*(theta-p_1)^2
             // ===============================================
             // p_0/k_B [K/rad^2]
             // p_1     [degrees]
           fprintf(FilePtr,"harm        %5d %5d %5d  %12.6f %12.6f\n",
             A+1,B+1,C+1,parms[0]*ENERGY_TO_KELVIN,parms[1]*RAD2DEG);
           break;
        }
      }

      // print internal dihedral
      fprintf(FilePtr,"DIHEDRALS %d\n",Components[i].NumberOfTorsions);
      for(k=0;k<Components[i].NumberOfTorsions;k++)
      {
        A=Components[i].Torsions[k].A;
        B=Components[i].Torsions[k].B;
        C=Components[i].Torsions[k].C;
        D=Components[i].Torsions[k].D;

        parms=Components[i].TorsionArguments[k];
        switch(Components[i].TorsionType[k])
        {
          case TRAPPE_DIHEDRAL:
            // p_0[0]+p_1*(1+cos(phi))+p_2*(1-cos(2*phi))+p_3*(1+cos(3*phi))
            // =============================================================
            // p_0/k_B [K]
            // p_1/k_B [K]
            // p_2/k_B [K]
            // p_3/k_B [K]
            fprintf(FilePtr,"cos3        %5d %5d %5d %5d  %12.6f %12.6f %12.6f %12.6f %12.6f\n",
              A+1,B+1,C+1,D+1,
              2.0*parms[1]*ENERGY_TO_KELVIN,
              2.0*parms[2]*ENERGY_TO_KELVIN,
              2.0*parms[3]*ENERGY_TO_KELVIN,
              parms[7],
              parms[6]);
            break;
        }
      }

      total=0;
      for(l=0;l<Components[i].NumberOfGroups;l++)
        if(Components[i].Groups[l].Rigid) total++;

      if(total>0)
      {
        fprintf(FilePtr,"RIGID %d\n",total);
        for(l=0;l<Components[i].NumberOfGroups;l++)
        {
          if(Components[i].Groups[l].Rigid) // rigid unit
          {
            fprintf(FilePtr,"%5d ",Components[i].Groups[l].NumberOfGroupAtoms);
            for(k=0;k<Components[i].Groups[l].NumberOfGroupAtoms;k++)
            {
              A=Components[i].Groups[l].Atoms[k];
              fprintf(FilePtr,"%5d",A+1);
            }
          }
        }
        fprintf(FilePtr,"\n");
      }
      fprintf(FilePtr,"FINISH\n");
    }
  }

  total=0;
  for(i=0;i<NumberOfPseudoAtoms;i++)
  {
    if(NumberOfPseudoAtomsType[CurrentSystem][i]>0)
    {
      for(j=i;j<NumberOfPseudoAtoms;j++)
      {
        if(NumberOfPseudoAtomsType[CurrentSystem][j]>0)
        {
          if(PotentialType[i][j]==LENNARD_JONES) total++;
        }
      }
    }
  }

  fprintf(FilePtr,"VDW %d\n",total);
  for(i=0;i<NumberOfPseudoAtoms;i++)
  {
    if(NumberOfPseudoAtomsType[CurrentSystem][i]>0)
    {
      for(j=i;j<NumberOfPseudoAtoms;j++)
      {
        if(NumberOfPseudoAtomsType[CurrentSystem][j]>0)
        {
          switch(PotentialType[i][j])
          {
            case ZERO_POTENTIAL:
            case ZERO_POTENTIAL_CONTINUOUS_FRACTIONAL:
              break;
            case LENNARD_JONES:
              // 4*p_0*((p_1/r)^12-(p_1/r)^6)
              // ======================================================================================
              // p_0/k_B [K]    strength parameter epsilon
              // p_1     [A]    size parameter sigma
              // p_2/k_B [K]    (non-zero for a shifted potential)
              fprintf(FilePtr,"%6s %6s lj %12.6lf %12.6lf\n",
                PseudoAtoms[i].Name,
                PseudoAtoms[j].Name,
                (double)PotentialParms[i][j][0]*ENERGY_TO_KELVIN,
                (double)PotentialParms[i][j][1]);
              break;

          }
        }
      }
    }
  }

  fprintf(FilePtr,"CLOSE\n");
  fclose(FilePtr);
}


/*********************************************************************************************************
 * Name       | ReadFrameworkDefinition                                                                  *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Read the definition (bond, bend, torsion, etc) of a framework component.                 *
 * Note       | The excluded Van der Waals and Coulombic interactions are automatically computed from    *
 *            | framework exclusion rules (1-2, 1-3, 1-4 omission, bond/bend/torsion omission).          *
 *********************************************************************************************************/

int ReadFrameworkDefinition(void)
{
  int i,j,k,l,m,n,f1,index_excluded,present,nr_atoms;
  int A,B,C,D,A1,B1,C1,TypeA,TypeB,TypeC,TypeD,CurrentTypeA,CurrentTypeB,CurrentTypeC,CurrentTypeD,index,ListTypeC,ListTypeD;
  char buffer[256],TypeNameA[256],TypeNameB[256],TypeNameC[256],TypeNameD[256];
  char line[1024];
  char *arg_pointer;
  REAL arguments[10];
  int BondType=0,BendType=0,TorsionType=0,InversionBendType=0,ImproperTorsionType=0;
  int BondBondType=0,BondBendType=0,BendBendType=0,BondTorsionType=0,BendTorsionType=0;
  REAL r,rab,rbc,theta;
  double temp;
  VECTOR Rab,Rbc,dr;
  FILE *FilePtr;
  int int_temp1,int_temp2;
  char string[256];

  // Create an initial connectivity list
  // Based on the connectivity and defined rules the type of an atom can be modified (e.g. oxygen connected to an aluminum)
  // The list has to be created first because the shells of cores is possibly based on the modified atoms
  AllocateConnectivityList();
  MakeConnectivityList();
  SubstituteAtoms();
  ModifyAtomsConnectedToDefinedNeighbours();

  if(strlen(Framework[CurrentSystem].FrameworkDefinitions)<1) return 0;

  // try "framework.def", then "Framework,def", and then the repository
  sprintf(buffer,"%s.%s","framework","def");
  if(!(FilePtr=fopen(buffer,"r")))
  {
    sprintf(buffer,"%s.%s","Framework","def");
    if(!(FilePtr=fopen(buffer,"r")))
    {
      sprintf(buffer,"%s/share/raspa/framework/%s/%s.%s",RASPA_DIRECTORY,
              Framework[CurrentSystem].FrameworkDefinitions,"framework","def");

      if(!(FilePtr=fopen(buffer,"r")))
      {
        sprintf(buffer,"%s/share/raspa/framework/%s/%s.%s",RASPA_DIRECTORY,
                Framework[CurrentSystem].FrameworkDefinitions,"Framework","def");

        if(!(FilePtr=fopen(buffer,"r")))
        {
          fprintf(stderr, "Error:  file %s does not exist.\n",buffer);
          exit(1);
        }
      }
    }
  }

  ReadLine(line,1024,FilePtr); // skip line

  ReadLine(line,1024,FilePtr);
  sscanf(line,"%d%d%d%d%d%d%d%d%d%d%d%d%d",
       &Framework[CurrentSystem].NumberOfCoreShellDefinitions,
       &Framework[CurrentSystem].NumberOfBondsDefinitions,
       &Framework[CurrentSystem].NumberOfBondDipoleDefinitions,
       &Framework[CurrentSystem].NumberOfUreyBradleyDefinitions,
       &Framework[CurrentSystem].NumberOfBendDefinitions,
       &Framework[CurrentSystem].NumberOfInversionBendDefinitions,
       &Framework[CurrentSystem].NumberOfTorsionDefinitions,
       &Framework[CurrentSystem].NumberOfImproperTorsionDefinitions,
       //&Framework[CurrentSystem].NumberOfOutOfPlaneDefinitions,
       &Framework[CurrentSystem].NumberOfBondBondDefinitions,
       &Framework[CurrentSystem].NumberOfBondBendDefinitions,
       &Framework[CurrentSystem].NumberOfBendBendDefinitions,
       &Framework[CurrentSystem].NumberOfBondTorsionDefinitions,
       &Framework[CurrentSystem].NumberOfBendTorsionDefinitions);


  // read the core and the shells
  index_excluded=0;
  if(Framework[CurrentSystem].NumberOfCoreShellDefinitions>0)
  {
    Framework[CurrentSystem].CoreShellConnectivity=(int**)calloc(Framework[CurrentSystem].NumberOfFrameworks,sizeof(int*));

    for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
      Framework[CurrentSystem].NumberOfCoreShells[f1]=0;

    Framework[CurrentSystem].CoreShellDefinitions=(PAIR*)calloc(Framework[CurrentSystem].NumberOfCoreShellDefinitions,sizeof(PAIR));
    Framework[CurrentSystem].NumberOfCoreShellsPerType=(int*)calloc(Framework[CurrentSystem].NumberOfCoreShellDefinitions,sizeof(int));

    ReadLine(line,1024,FilePtr); // skip line
    for(i=0;i<Framework[CurrentSystem].NumberOfCoreShellDefinitions;i++)
    {
      ReadLine(line,1024,FilePtr);
      sscanf(line,"%s%s",TypeNameA,TypeNameB);
      TypeA=ReturnPseudoAtomNumber(TypeNameA);
      TypeB=ReturnPseudoAtomNumber(TypeNameB);
      PseudoAtoms[TypeA].CoreShell=CORE;
      PseudoAtoms[TypeB].CoreShell=SHELL;

      Framework[CurrentSystem].CoreShellDefinitions[i].A=TypeA;
      Framework[CurrentSystem].CoreShellDefinitions[i].B=TypeB;

      // count shells
      for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
      {
        for(A=0;A<Framework[CurrentSystem].NumberOfAtoms[f1];A++)
          if(Framework[CurrentSystem].Atoms[f1][A].Type==TypeA)
            Framework[CurrentSystem].NumberOfCoreShells[f1]++;
      }
    }
  }

  if(Framework[CurrentSystem].NumberOfCoreShellDefinitions>0)
  {
    // allocate the shells
    for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
    {
      if(Framework[CurrentSystem].NumberOfCoreShells[f1]>0)
      {
        nr_atoms=Framework[CurrentSystem].NumberOfAtoms[f1]+Framework[CurrentSystem].NumberOfCoreShells[f1];

        // create memory for the core-shell connectivity array
        Framework[CurrentSystem].CoreShellConnectivity[f1]=(int*)calloc(nr_atoms,sizeof(int));
        for(i=0;i<nr_atoms;i++)
          Framework[CurrentSystem].CoreShellConnectivity[f1][i]=-1;

        // create additional memory for the shells
        Framework[CurrentSystem].Atoms[f1]=(ATOM*)realloc(Framework[CurrentSystem].Atoms[f1],
              nr_atoms*sizeof(ATOM));
      }
    }
  }


  // create shells around the cores
  for(i=0;i<Framework[CurrentSystem].NumberOfCoreShellDefinitions;i++)
  {
    TypeA=Framework[CurrentSystem].CoreShellDefinitions[i].A;
    TypeB=Framework[CurrentSystem].CoreShellDefinitions[i].B;
    for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
    {
      index=Framework[CurrentSystem].NumberOfAtoms[f1];
      for(A=0;A<Framework[CurrentSystem].NumberOfAtoms[f1];A++)
        if(Framework[CurrentSystem].Atoms[f1][A].Type==TypeA)
        {
          Framework[CurrentSystem].Atoms[f1][index].Type=TypeB;

          NumberOfPseudoAtomsCount[CurrentSystem][TypeB]++;
          NumberOfPseudoAtomsType[CurrentSystem][TypeB]++;
          Framework[CurrentSystem].FrameworkMass+=PseudoAtoms[TypeB].Mass;
          Framework[CurrentSystem].Atoms[f1][index].Charge+=PseudoAtoms[TypeB].Charge1;

          // set position of the shell close to the core
          Framework[CurrentSystem].Atoms[f1][index].Position.x=Framework[CurrentSystem].Atoms[f1][A].Position.x+(RandomNumber()*0.1);
          Framework[CurrentSystem].Atoms[f1][index].Position.y=Framework[CurrentSystem].Atoms[f1][A].Position.y+(RandomNumber()*0.1);
          Framework[CurrentSystem].Atoms[f1][index].Position.z=Framework[CurrentSystem].Atoms[f1][A].Position.z+(RandomNumber()*0.1);

          Framework[CurrentSystem].CoreShellConnectivity[f1][index]=A;
          Framework[CurrentSystem].CoreShellConnectivity[f1][A]=index;

          Framework[CurrentSystem].NumberOfCoreShellsPerType[i]++;

          index++;
          Framework[CurrentSystem].TotalNumberOfAtoms++;
        }
      Framework[CurrentSystem].NumberOfAtoms[f1]=index;
    }
  }

  // Now that shells have been created the list has to be recreated including the shells around the cores
  FreeAllocateConnectivityList();
  AllocateConnectivityList();
  MakeConnectivityList();
  //PrintConnectivityList();
  ModifyAtomsConnectedToDefinedNeighbours();

  // reading Bond-data
  if(Framework[CurrentSystem].NumberOfBondsDefinitions>0)
  {
    // alloc memory
    Framework[CurrentSystem].BondDefinitionType=(int*)calloc(Framework[CurrentSystem].NumberOfBondsDefinitions,sizeof(int));
    Framework[CurrentSystem].BondDefinitions=(PAIR*)calloc(Framework[CurrentSystem].NumberOfBondsDefinitions,sizeof(PAIR));
    Framework[CurrentSystem].BondArgumentDefinitions=(REAL(*)[MAX_BOND_POTENTIAL_ARGUMENTS])
      calloc(Framework[CurrentSystem].NumberOfBondsDefinitions,sizeof(REAL[MAX_BOND_POTENTIAL_ARGUMENTS]));
    Framework[CurrentSystem].NumberOfBondsPerType=(int*)calloc(Framework[CurrentSystem].NumberOfBondsDefinitions,sizeof(int));

    for(CurrentFramework=0;CurrentFramework<Framework[CurrentSystem].NumberOfFrameworks;CurrentFramework++)
    {
      Framework[CurrentSystem].MaxNumberOfBonds[CurrentFramework]=4096;
      // allocate enough memory, it will be resized when needed
      Framework[CurrentSystem].Bonds[CurrentFramework]=(PAIR*)calloc(Framework[CurrentSystem].MaxNumberOfBonds[CurrentFramework],sizeof(PAIR));
      Framework[CurrentSystem].BondType[CurrentFramework]=(int*)calloc(Framework[CurrentSystem].MaxNumberOfBonds[CurrentFramework],sizeof(int));
      Framework[CurrentSystem].BondArguments[CurrentFramework]=(REAL(*)[MAX_BOND_POTENTIAL_ARGUMENTS])
              calloc(Framework[CurrentSystem].MaxNumberOfBonds[CurrentFramework],sizeof(REAL[MAX_BOND_POTENTIAL_ARGUMENTS]));

      Framework[CurrentSystem].NumberOfBonds[CurrentFramework]=0;
    }

    ReadLine(line,1024,FilePtr); // skip line
    for(i=0;i<Framework[CurrentSystem].NumberOfBondsDefinitions;i++)
    {
      ReadLine(line,1024,FilePtr);
      arg_pointer=line;
      sscanf(line,"%s%s%s%n",TypeNameA,TypeNameB,buffer,&n);
      TypeA=ReturnPseudoAtomNumber(TypeNameA);
      TypeB=ReturnPseudoAtomNumber(TypeNameB);

      // determine bond-type
      for(j=0;j<NR_BOND_TYPES;j++)
        if(strcasecmp(BondTypes[j].Name,buffer)==0)
          BondType=j;

      // read arguments
      for(j=0;j<BondTypes[BondType].nr_args;j++)
      {
        arg_pointer+=n;
        sscanf(arg_pointer,"%lf%n",&temp,&n);
        arguments[j]=(REAL)temp;
      }

      Framework[CurrentSystem].BondDefinitionType[i]=BondType;
      Framework[CurrentSystem].BondDefinitions[i].A=TypeA;
      Framework[CurrentSystem].BondDefinitions[i].B=TypeB;
      for(j=0;j<BondTypes[BondType].nr_args;j++)
        Framework[CurrentSystem].BondArgumentDefinitions[i][j]=arguments[j];

      // fill the appropriate bond-data
      for(CurrentFramework=0;CurrentFramework<Framework[CurrentSystem].NumberOfFrameworks;CurrentFramework++)
      {
        index=Framework[CurrentSystem].NumberOfBonds[CurrentFramework];
        for(A=0;A<Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];A++)
        {
          CurrentTypeA=Framework[CurrentSystem].Atoms[CurrentFramework][A].Type;
          for(k=0;k<Framework[CurrentSystem].Connectivity[CurrentFramework][A];k++)
          {
            B=GetNeighbour(CurrentSystem,CurrentFramework,A,k);
            CurrentTypeB=Framework[CurrentSystem].Atoms[CurrentFramework][B].Type;
            if((A<B)&&((TypeA==CurrentTypeA&&TypeB==CurrentTypeB)||
                     ((TypeA==CurrentTypeB&&TypeB==CurrentTypeA))))
            {
              Framework[CurrentSystem].Bonds[CurrentFramework][index].A=A;
              Framework[CurrentSystem].Bonds[CurrentFramework][index].B=B;

              Framework[CurrentSystem].BondType[CurrentFramework][index]=BondType;

              Framework[CurrentSystem].NumberOfBondsPerType[i]++;

              for(j=0;j<BondTypes[BondType].nr_args;j++)
                Framework[CurrentSystem].BondArguments[CurrentFramework][index][j]=arguments[j];

              // set to appropriate bond-distance
              switch(BondType)
              {
                case HARMONIC_BOND:
                  // 0.5*p0*SQR(r-p1);
                  // ===============================================
                  // p_0/k_B [K/A^2]   force constant
                  // p_1     [A]       reference bond distance
                  r=Framework[CurrentSystem].BondArguments[CurrentFramework][index][1];
                  if(r<1e-5)
                  {
                    dr.x=Framework[CurrentSystem].Atoms[CurrentFramework][A].Position.x-
                         Framework[CurrentSystem].Atoms[CurrentFramework][B].Position.x;
                    dr.y=Framework[CurrentSystem].Atoms[CurrentFramework][A].Position.y-
                         Framework[CurrentSystem].Atoms[CurrentFramework][B].Position.y;
                    dr.z=Framework[CurrentSystem].Atoms[CurrentFramework][A].Position.z-
                         Framework[CurrentSystem].Atoms[CurrentFramework][B].Position.z;
                    dr=ApplyBoundaryCondition(dr);
                    r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
                    Framework[CurrentSystem].BondArguments[CurrentFramework][index][1]=r;
                  }
                  Framework[CurrentSystem].BondArguments[CurrentFramework][index][0]*=KELVIN_TO_ENERGY;
                  break;
                case CORE_SHELL_SPRING:
                 // 0.5*p0*SQR(r);
                  // ===============================================
                  // p_0/k_B [K/A^2]   force constant
                  Framework[CurrentSystem].BondArguments[CurrentFramework][index][0]*=KELVIN_TO_ENERGY;
                  break;
                case MORSE_BOND:
                  // p_0*[(1.0-{exp(-p_1*(r-p_2))})^2-1.0]
                  // ===============================================
                  // p_0/k_B [K]       force constant
                  // p_1     [A^-1]    parameter
                  // p_2     [A]       reference bond distance
                  r=Framework[CurrentSystem].BondArguments[CurrentFramework][index][2];
                  if(r<1e-5)
                  {
                    dr.x=Framework[CurrentSystem].Atoms[CurrentFramework][A].Position.x-
                         Framework[CurrentSystem].Atoms[CurrentFramework][B].Position.x;
                    dr.y=Framework[CurrentSystem].Atoms[CurrentFramework][A].Position.y-
                         Framework[CurrentSystem].Atoms[CurrentFramework][B].Position.y;
                    dr.z=Framework[CurrentSystem].Atoms[CurrentFramework][A].Position.z-
                         Framework[CurrentSystem].Atoms[CurrentFramework][B].Position.z;
                    dr=ApplyBoundaryCondition(dr);
                    r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
                    Framework[CurrentSystem].BondArguments[CurrentFramework][index][2]=r;
                  }
                  Framework[CurrentSystem].BondArguments[CurrentFramework][index][0]*=KELVIN_TO_ENERGY;
                  break;
                case LJ_12_6_BOND:
                  // A/r_ij^12-B/r_ij^6
                  // ===============================================
                  // p_0/k_B [K A^12]
                  // p_1/k_B [K A^6]
                  Framework[CurrentSystem].BondArguments[CurrentFramework][index][0]*=KELVIN_TO_ENERGY;
                  Framework[CurrentSystem].BondArguments[CurrentFramework][index][1]*=KELVIN_TO_ENERGY;
                  break;
                case LENNARD_JONES_BOND:
                  // 4*p_0*((p_1/r)^12-(p_1/r)^6)
                  // ===============================================
                  // p_0/k_B [K]
                  // p_1     [A]
                  Framework[CurrentSystem].BondArguments[CurrentFramework][index][0]*=KELVIN_TO_ENERGY;
                  break;
                case BUCKINGHAM_BOND:
                  // p_0*exp(-p_1 r)-p_2/r^6
                  // ===============================================
                  // p_0/k_B [K]
                  // p_1     [A^-1]
                  // p_2/k_B [K A^6]
                  Framework[CurrentSystem].BondArguments[CurrentFramework][index][0]*=KELVIN_TO_ENERGY;
                  Framework[CurrentSystem].BondArguments[CurrentFramework][index][2]*=KELVIN_TO_ENERGY;
                  break;
                case RESTRAINED_HARMONIC_BOND:
                  // 0.5*p_0*(r-p_1)^2                   |r-p_1|<=p_2
                  // 0.5*p_0*p_2^2+p_0*p_2*(|r-p_1|-p_2) |r-p_1|>p_2
                  // ===============================================
                  // p_0/k_B [K/A^2]
                  // p_1     [A]
                  // p_2     [A]
                  Framework[CurrentSystem].BondArguments[CurrentFramework][index][0]*=KELVIN_TO_ENERGY;
                  break;
                case QUARTIC_BOND:
                  // (1/2)*p_0*(r-p_1)^2+(1/3)*p_2*(r-p_1)^3+(1/4)*p_3*(r-p_1)^4
                  // ===========================================================
                  // p_0/k_B [K/A^2]
                  // p_1     [A]
                  // p_2/k_B [K/A^3]
                  // p_3/k_B [K/A^4]
                case CFF_QUARTIC_BOND:
                  // p_0*(r-p_1)^2+p_2*(r-p_1)^3+p_3*(r-p_1)^4
                  // ===============================================
                  // p_0/k_B [K/A^2]
                  // p_1     [A]
                  // p_2/k_B [K/A^3]
                  // p_3/k_B [K/A^4]
                  r=Framework[CurrentSystem].BondArguments[CurrentFramework][index][1];
                  if(r<1e-5)
                  {
                    dr.x=Framework[CurrentSystem].Atoms[CurrentFramework][A].Position.x-
                         Framework[CurrentSystem].Atoms[CurrentFramework][B].Position.x;
                    dr.y=Framework[CurrentSystem].Atoms[CurrentFramework][A].Position.y-
                         Framework[CurrentSystem].Atoms[CurrentFramework][B].Position.y;
                    dr.z=Framework[CurrentSystem].Atoms[CurrentFramework][A].Position.z-
                         Framework[CurrentSystem].Atoms[CurrentFramework][B].Position.z;
                    dr=ApplyBoundaryCondition(dr);
                    r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
                    Framework[CurrentSystem].BondArguments[CurrentFramework][index][1]=r;
                  }
                  Framework[CurrentSystem].BondArguments[CurrentFramework][index][0]*=KELVIN_TO_ENERGY;
                  Framework[CurrentSystem].BondArguments[CurrentFramework][index][2]*=KELVIN_TO_ENERGY;
                  Framework[CurrentSystem].BondArguments[CurrentFramework][index][3]*=KELVIN_TO_ENERGY;
                  break;
                case MM3_BOND:
                  // p_0*(r-p_1)^2*(1.0-2.55*(r-p_1)+(7.0/12.0)*2.55^2*(r-p_1)^2)
                  // =================================================================
                  // p_0     [mdyne/A molecule]
                  // p_1     [A]
                  r=Framework[CurrentSystem].BondArguments[CurrentFramework][index][1];
                  if(r<1e-5)
                  {
                    dr.x=Framework[CurrentSystem].Atoms[CurrentFramework][A].Position.x-
                         Framework[CurrentSystem].Atoms[CurrentFramework][B].Position.x;
                    dr.y=Framework[CurrentSystem].Atoms[CurrentFramework][A].Position.y-
                         Framework[CurrentSystem].Atoms[CurrentFramework][B].Position.y;
                    dr.z=Framework[CurrentSystem].Atoms[CurrentFramework][A].Position.z-
                         Framework[CurrentSystem].Atoms[CurrentFramework][B].Position.z;
                    dr=ApplyBoundaryCondition(dr);
                    r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
                    Framework[CurrentSystem].BondArguments[CurrentFramework][index][1]=r;
                  }
                  // MM3 input is in mdyne/A-molecule
                  // the MM3 force field uses the factor 71.94 to convert to kcal/mole
                  Framework[CurrentSystem].BondArguments[CurrentFramework][index][0]*=71.94*KCAL_PER_MOL_TO_ENERGY;
                  break;
                case MEASURE_BOND:
                  break;
                default:
                  fprintf(stderr, "Undefined Bond potential in routine 'ReadFrameworkDefinition' ('framework.c')\n");
                  exit(0);
                  break;
              }
              Framework[CurrentSystem].NumberOfBonds[CurrentFramework]++;
              index=Framework[CurrentSystem].NumberOfBonds[CurrentFramework];
              if(index>=Framework[CurrentSystem].MaxNumberOfBonds[CurrentFramework])
              {
                Framework[CurrentSystem].MaxNumberOfBonds[CurrentFramework]+=4096;
                Framework[CurrentSystem].Bonds[CurrentFramework]=(PAIR*)realloc(Framework[CurrentSystem].Bonds[CurrentFramework],
                             Framework[CurrentSystem].MaxNumberOfBonds[CurrentFramework]*sizeof(PAIR));
                Framework[CurrentSystem].BondType[CurrentFramework]=(int*)realloc(Framework[CurrentSystem].BondType[CurrentFramework],
                             Framework[CurrentSystem].MaxNumberOfBonds[CurrentFramework]*sizeof(int));
                Framework[CurrentSystem].BondArguments[CurrentFramework]=(REAL(*)[MAX_BOND_POTENTIAL_ARGUMENTS])
                    realloc(Framework[CurrentSystem].BondArguments[CurrentFramework],
                    Framework[CurrentSystem].MaxNumberOfBonds[CurrentFramework]*sizeof(REAL[MAX_BOND_POTENTIAL_ARGUMENTS]));
              }
            }
          }
        }
      }
    }
  }

  // reading Bond-Dipole-data
  index=0;
  if(Framework[CurrentSystem].NumberOfBondDipoleDefinitions>0)
  {
    Framework[CurrentSystem].BondDipoleDefinitions=(PAIR*)calloc(Framework[CurrentSystem].NumberOfBondDipoleDefinitions,sizeof(PAIR));
    Framework[CurrentSystem].BondDipoleArgumentDefinition=(REAL*)calloc(Framework[CurrentSystem].NumberOfBondDipoleDefinitions,sizeof(REAL));
    Framework[CurrentSystem].NumberOfBondDipolesPerType=(int*)calloc(Framework[CurrentSystem].NumberOfBondDipoleDefinitions,sizeof(int));

      // allocate first dimension of bond-dipoles
    for(CurrentFramework=0;CurrentFramework<Framework[CurrentSystem].NumberOfFrameworks;CurrentFramework++)
    {
      Framework[CurrentSystem].MaxNumberOfBondDipoles[CurrentFramework]=4096;

      Framework[CurrentSystem].BondDipoles[CurrentFramework]=(PAIR*)calloc(Framework[CurrentSystem].MaxNumberOfBondDipoles[CurrentFramework],sizeof(PAIR));
      Framework[CurrentSystem].BondDipoleMagnitude[CurrentFramework]=(REAL*)calloc(Framework[CurrentSystem].MaxNumberOfBondDipoles[CurrentFramework],sizeof(REAL));
    }

    ReadLine(line,1024,FilePtr); // skip line
    for(i=0;i<Framework[CurrentSystem].NumberOfBondDipoleDefinitions;i++)
    {
      ReadLine(line,1024,FilePtr);
      sscanf(line,"%s%s%lf",TypeNameA,TypeNameB,&arguments[0]);
      TypeA=ReturnPseudoAtomNumber(TypeNameA);
      TypeB=ReturnPseudoAtomNumber(TypeNameB);

      Framework[CurrentSystem].BondDipoleDefinitions[i].A=TypeA;
      Framework[CurrentSystem].BondDipoleDefinitions[i].B=TypeB;
      Framework[CurrentSystem].BondDipoleArgumentDefinition[i]=arguments[0]/DEBYE_CONVERSION_FACTOR;

      // fill the appropriate bonddipole-data
      for(CurrentFramework=0;CurrentFramework<Framework[CurrentSystem].NumberOfFrameworks;CurrentFramework++)
      {
        index=Framework[CurrentSystem].NumberOfBondDipoles[CurrentFramework];
        for(A=0;A<Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];A++)
        {
          CurrentTypeA=Framework[CurrentSystem].Atoms[CurrentFramework][A].Type;
          for(k=0;k<Framework[CurrentSystem].Connectivity[CurrentFramework][A];k++)
          {
            B=GetNeighbour(CurrentSystem,CurrentFramework,A,k);
            CurrentTypeB=Framework[CurrentSystem].Atoms[CurrentFramework][B].Type;
            // only allow for pairs in the correct order
            if((TypeA==CurrentTypeA)&&(TypeB==CurrentTypeB))
            {
              Framework[CurrentSystem].BondDipoles[CurrentFramework][index].A=A;
              Framework[CurrentSystem].BondDipoles[CurrentFramework][index].B=B;

              Framework[CurrentSystem].NumberOfBondDipolesPerType[i]++;

              Framework[CurrentSystem].BondDipoleMagnitude[CurrentFramework][index]=arguments[0]/DEBYE_CONVERSION_FACTOR;

              Framework[CurrentSystem].NumberOfBondDipoles[CurrentFramework]++;
              index=Framework[CurrentSystem].NumberOfBondDipoles[CurrentFramework];
              if(index>=Framework[CurrentSystem].MaxNumberOfBondDipoles[CurrentFramework])
              {
                 Framework[CurrentSystem].MaxNumberOfBondDipoles[CurrentFramework]+=4096;

                 Framework[CurrentSystem].BondDipoles[CurrentFramework]=(PAIR*)realloc(Framework[CurrentSystem].BondDipoles[CurrentFramework],
                   Framework[CurrentSystem].MaxNumberOfBondDipoles[CurrentFramework]*sizeof(PAIR));
                 Framework[CurrentSystem].BondDipoleMagnitude[CurrentFramework]=(REAL*)realloc(Framework[CurrentSystem].BondDipoleMagnitude[CurrentFramework],
                   Framework[CurrentSystem].MaxNumberOfBondDipoles[CurrentFramework]*sizeof(REAL));
              }
            }
          }
        }
      }
    }
  }

  // reading Urey-Bradley-data
  if(Framework[CurrentSystem].NumberOfUreyBradleyDefinitions>0)
  {
    // alloc memory
    Framework[CurrentSystem].UreyBradleyDefinitionType=(int*)calloc(Framework[CurrentSystem].NumberOfUreyBradleyDefinitions,sizeof(int));
    Framework[CurrentSystem].UreyBradleyDefinitions=(TRIPLE*)calloc(Framework[CurrentSystem].NumberOfUreyBradleyDefinitions,sizeof(TRIPLE));
    Framework[CurrentSystem].UreyBradleyArgumentDefinitions=(REAL(*)[MAX_UREYBRADLEY_POTENTIAL_ARGUMENTS])
      calloc(Framework[CurrentSystem].NumberOfUreyBradleyDefinitions,sizeof(REAL[MAX_UREYBRADLEY_POTENTIAL_ARGUMENTS]));
    Framework[CurrentSystem].NumberOfUreyBradleysPerType=(int*)calloc(Framework[CurrentSystem].NumberOfUreyBradleyDefinitions,sizeof(int));

    for(CurrentFramework=0;CurrentFramework<Framework[CurrentSystem].NumberOfFrameworks;CurrentFramework++)
    {
      Framework[CurrentSystem].MaxNumberOfUreyBradleys[CurrentFramework]=4096;
      // allocate enough memory, it will be resized when needed
      Framework[CurrentSystem].UreyBradleys[CurrentFramework]=(TRIPLE*)calloc(Framework[CurrentSystem].MaxNumberOfUreyBradleys[CurrentFramework],sizeof(TRIPLE));
      Framework[CurrentSystem].UreyBradleyType[CurrentFramework]=(int*)calloc(Framework[CurrentSystem].MaxNumberOfUreyBradleys[CurrentFramework],sizeof(int));
      Framework[CurrentSystem].UreyBradleyArguments[CurrentFramework]=(REAL(*)[MAX_UREYBRADLEY_POTENTIAL_ARGUMENTS])
              calloc(Framework[CurrentSystem].MaxNumberOfUreyBradleys[CurrentFramework],sizeof(REAL[MAX_UREYBRADLEY_POTENTIAL_ARGUMENTS]));
    }

    ReadLine(line,1024,FilePtr); // skip line
    for(i=0;i<Framework[CurrentSystem].NumberOfUreyBradleyDefinitions;i++)
    {
      ReadLine(line,1024,FilePtr);
      arg_pointer=line;
      sscanf(line,"%s%s%s%s%n",TypeNameA,TypeNameB,TypeNameC,buffer,&n);
      TypeA=ReturnPseudoAtomNumber(TypeNameA);
      TypeB=ReturnPseudoAtomNumber(TypeNameB);
      TypeC=ReturnPseudoAtomNumber(TypeNameC);

      // determine bond-type
      for(j=0;j<NR_BOND_TYPES;j++)
        if(strcasecmp(BondTypes[j].Name,buffer)==0)
          BondType=j;

      // read arguments
      for(j=0;j<BondTypes[BondType].nr_args;j++)
      {
        arg_pointer+=n;
        sscanf(arg_pointer,"%lf%n",&temp,&n);
        arguments[j]=(REAL)temp;
      }

      Framework[CurrentSystem].UreyBradleyDefinitionType[i]=BondType;
      Framework[CurrentSystem].UreyBradleyDefinitions[i].A=TypeA;
      Framework[CurrentSystem].UreyBradleyDefinitions[i].B=TypeB;
      Framework[CurrentSystem].UreyBradleyDefinitions[i].C=TypeC;
      for(j=0;j<UreyBradleyTypes[BondType].nr_args;j++)
        Framework[CurrentSystem].UreyBradleyArgumentDefinitions[i][j]=arguments[j];

      // fill the appropriate bond-data
      for(CurrentFramework=0;CurrentFramework<Framework[CurrentSystem].NumberOfFrameworks;CurrentFramework++)
      {
        index=Framework[CurrentSystem].NumberOfUreyBradleys[CurrentFramework];
        for(A=0;A<Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];A++)
        {
          CurrentTypeA=Framework[CurrentSystem].Atoms[CurrentFramework][A].Type;
          for(k=0;k<Framework[CurrentSystem].Connectivity[CurrentFramework][A];k++)
          {
            B=GetNeighbour(CurrentSystem,CurrentFramework,A,k);
            CurrentTypeB=Framework[CurrentSystem].Atoms[CurrentFramework][B].Type;

            if(CurrentTypeB==TypeB)
            {
              for(l=0;l<Framework[CurrentSystem].Connectivity[CurrentFramework][B];l++)
              {
                C=GetNeighbour(CurrentSystem,CurrentFramework,B,l);
                CurrentTypeC=Framework[CurrentSystem].Atoms[CurrentFramework][C].Type;

                if((A<C)&&((CurrentTypeA==TypeA&&CurrentTypeC==TypeC)||(CurrentTypeA==TypeC&&CurrentTypeC==TypeA)))
                {
                  Framework[CurrentSystem].UreyBradleys[CurrentFramework][index].A=A;
                  Framework[CurrentSystem].UreyBradleys[CurrentFramework][index].B=B;
                  Framework[CurrentSystem].UreyBradleys[CurrentFramework][index].C=C;

                  Framework[CurrentSystem].UreyBradleyType[CurrentFramework][index]=BondType;
                  Framework[CurrentSystem].NumberOfUreyBradleysPerType[i]++;

                  for(j=0;j<BondTypes[BondType].nr_args;j++)
                    Framework[CurrentSystem].UreyBradleyArguments[CurrentFramework][index][j]=arguments[j];

                  // set to appropriate bond-distance
                  switch(BondType)
                  {
                    case HARMONIC_UREYBRADLEY:
                      // 0.5*p0*SQR(r-p1);
                      // ===============================================
                      // p_0/k_B [K/A^2]   force constant
                      // p_1     [A]       reference bond distance
                      r=Framework[CurrentSystem].UreyBradleyArguments[CurrentFramework][index][1];
                      if(r<1e-5)
                      {
                        dr.x=Framework[CurrentSystem].Atoms[CurrentFramework][A].Position.x-
                             Framework[CurrentSystem].Atoms[CurrentFramework][C].Position.x;
                        dr.y=Framework[CurrentSystem].Atoms[CurrentFramework][A].Position.y-
                             Framework[CurrentSystem].Atoms[CurrentFramework][C].Position.y;
                        dr.z=Framework[CurrentSystem].Atoms[CurrentFramework][A].Position.z-
                             Framework[CurrentSystem].Atoms[CurrentFramework][C].Position.z;
                        dr=ApplyBoundaryCondition(dr);
                        r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
                        Framework[CurrentSystem].UreyBradleyArguments[CurrentFramework][index][1]=r;
                      }
                      Framework[CurrentSystem].UreyBradleyArguments[CurrentFramework][index][0]*=KELVIN_TO_ENERGY;
                      break;
                    case MORSE_UREYBRADLEY:
                      // p_0*[(1.0-{exp(-p_1*(r-p_2))})^2-1.0]
                      // ===============================================
                      // p_0/k_B [K]       force constant
                      // p_1     [A^-1]    parameter
                      // p_2     [A]       reference bond distance
                      r=Framework[CurrentSystem].UreyBradleyArguments[CurrentFramework][index][2];
                      if(r<1e-5)
                      {
                        dr.x=Framework[CurrentSystem].Atoms[CurrentFramework][A].Position.x-
                             Framework[CurrentSystem].Atoms[CurrentFramework][C].Position.x;
                        dr.y=Framework[CurrentSystem].Atoms[CurrentFramework][A].Position.y-
                             Framework[CurrentSystem].Atoms[CurrentFramework][C].Position.y;
                        dr.z=Framework[CurrentSystem].Atoms[CurrentFramework][A].Position.z-
                             Framework[CurrentSystem].Atoms[CurrentFramework][C].Position.z;
                        dr=ApplyBoundaryCondition(dr);
                        r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
                        Framework[CurrentSystem].UreyBradleyArguments[CurrentFramework][index][2]=r;
                      }
                      Framework[CurrentSystem].UreyBradleyArguments[CurrentFramework][index][0]*=KELVIN_TO_ENERGY;
                      break;
                    case LJ_12_6_UREYBRADLEY:
                      // A/r_ij^12-B/r_ij^6
                      // ===============================================
                      // p_0/k_B [K A^12]
                      // p_1/k_B [K A^6]
                      Framework[CurrentSystem].UreyBradleyArguments[CurrentFramework][index][0]*=KELVIN_TO_ENERGY;
                      Framework[CurrentSystem].UreyBradleyArguments[CurrentFramework][index][1]*=KELVIN_TO_ENERGY;
                      break;
                    case LENNARD_JONES_UREYBRADLEY:
                      // 4*p_0*((p_1/r)^12-(p_1/r)^6)
                      // ===============================================
                      // p_0/k_B [K]
                      // p_1     [A]
                      Framework[CurrentSystem].UreyBradleyArguments[CurrentFramework][index][0]*=KELVIN_TO_ENERGY;
                      break;
                    case BUCKINGHAM_UREYBRADLEY:
                      // p_0*exp(-p_1 r)-p_2/r^6
                      // ===============================================
                      // p_0/k_B [K]
                      // p_1     [A^-1]
                      // p_2/k_B [K A^6]
                      Framework[CurrentSystem].UreyBradleyArguments[CurrentFramework][index][0]*=KELVIN_TO_ENERGY;
                      Framework[CurrentSystem].UreyBradleyArguments[CurrentFramework][index][2]*=KELVIN_TO_ENERGY;
                      break;
                    case RESTRAINED_HARMONIC_UREYBRADLEY:
                      // 0.5*p_0*(r-p_1)^2                   |r-p_1|<=p_2
                      // 0.5*p_0*p_2^2+p_0*p_2*(|r-p_1|-p_2) |r-p_1|>p_2
                      // ===============================================
                      // p_0/k_B [K/A^2]
                      // p_1     [A]
                      // p_2     [A]
                      Framework[CurrentSystem].UreyBradleyArguments[CurrentFramework][index][0]*=KELVIN_TO_ENERGY;
                      break;
                    case QUARTIC_UREYBRADLEY:
                      // p_0*(r-p_1)^2+p_2*(r-p_1)^3+p_3*(r-p_1)^4
                      // ===============================================
                      // p_0/k_B [K/A^2]
                      // p_1     [A]
                      // p_2/k_B [K/A^3]
                      // p_3/k_B [K/A^4]
                    case CFF_QUARTIC_UREYBRADLEY:
                      // p_0*(r-p_1)^2+p_2*(r-p_1)^3+p_3*(r-p_1)^4
                      // ===============================================
                      // p_0/k_B [K/A^2]
                      // p_1     [A]
                      // p_2/k_B [K/A^3]
                      // p_3/k_B [K/A^4]
                      r=Framework[CurrentSystem].UreyBradleyArguments[CurrentFramework][index][1];
                      if(r<1e-5)
                      {
                        dr.x=Framework[CurrentSystem].Atoms[CurrentFramework][A].Position.x-
                             Framework[CurrentSystem].Atoms[CurrentFramework][C].Position.x;
                        dr.y=Framework[CurrentSystem].Atoms[CurrentFramework][A].Position.y-
                             Framework[CurrentSystem].Atoms[CurrentFramework][C].Position.y;
                        dr.z=Framework[CurrentSystem].Atoms[CurrentFramework][A].Position.z-
                             Framework[CurrentSystem].Atoms[CurrentFramework][C].Position.z;
                        dr=ApplyBoundaryCondition(dr);
                        r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
                        Framework[CurrentSystem].UreyBradleyArguments[CurrentFramework][index][1]=r;
                      }
                      Framework[CurrentSystem].UreyBradleyArguments[CurrentFramework][index][0]*=KELVIN_TO_ENERGY;
                      Framework[CurrentSystem].UreyBradleyArguments[CurrentFramework][index][2]*=KELVIN_TO_ENERGY;
                      Framework[CurrentSystem].UreyBradleyArguments[CurrentFramework][index][3]*=KELVIN_TO_ENERGY;
                      break;
                    case MM3_UREYBRADLEY:
                      // p_0*(r-p_1)^2*(1.0-2.55*(r-p_1)+(7.0/12.0)*2.55^2*(r-p_1)^2)
                      // =================================================================
                      // p_0     [mdyne/A molecule]
                      // p_1     [A]
                      r=Framework[CurrentSystem].UreyBradleyArguments[CurrentFramework][index][1];
                      if(r<1e-5)
                      {
                        dr.x=Framework[CurrentSystem].Atoms[CurrentFramework][A].Position.x-
                             Framework[CurrentSystem].Atoms[CurrentFramework][C].Position.x;
                        dr.y=Framework[CurrentSystem].Atoms[CurrentFramework][A].Position.y-
                             Framework[CurrentSystem].Atoms[CurrentFramework][C].Position.y;
                        dr.z=Framework[CurrentSystem].Atoms[CurrentFramework][A].Position.z-
                             Framework[CurrentSystem].Atoms[CurrentFramework][C].Position.z;
                        dr=ApplyBoundaryCondition(dr);
                        r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
                        Framework[CurrentSystem].UreyBradleyArguments[CurrentFramework][index][1]=r;
                      }
                      // MM3 input is in mdyne/A
                      // the MM3 force field uses the factor 71.94 to convert to kcal/mole
                      Framework[CurrentSystem].UreyBradleyArguments[CurrentFramework][index][0]*=71.94*KCAL_PER_MOL_TO_ENERGY;
                      break;
                    default:
                      fprintf(stderr, "Undefined Urey-Bradley potential in routine 'ReadFrameworkDefinition' ('framework.c')\n");
                      exit(0);
                      break;
                  }
                  Framework[CurrentSystem].NumberOfUreyBradleys[CurrentFramework]++;
                  index=Framework[CurrentSystem].NumberOfUreyBradleys[CurrentFramework];

                  if(index>=Framework[CurrentSystem].MaxNumberOfUreyBradleys[CurrentFramework])
                  {
                    Framework[CurrentSystem].MaxNumberOfUreyBradleys[CurrentFramework]+=4096;
                    Framework[CurrentSystem].UreyBradleys[CurrentFramework]=(TRIPLE*)realloc(Framework[CurrentSystem].UreyBradleys[CurrentFramework],
                                 Framework[CurrentSystem].MaxNumberOfUreyBradleys[CurrentFramework]*sizeof(TRIPLE));
                    Framework[CurrentSystem].UreyBradleyType[CurrentFramework]=(int*)realloc(Framework[CurrentSystem].UreyBradleyType[CurrentFramework],
                                 Framework[CurrentSystem].MaxNumberOfUreyBradleys[CurrentFramework]*sizeof(int));
                    Framework[CurrentSystem].UreyBradleyArguments[CurrentFramework]=(REAL(*)[MAX_UREYBRADLEY_POTENTIAL_ARGUMENTS])
                         realloc(Framework[CurrentSystem].UreyBradleyArguments[CurrentFramework],
                         Framework[CurrentSystem].MaxNumberOfUreyBradleys[CurrentFramework]*sizeof(REAL[MAX_UREYBRADLEY_POTENTIAL_ARGUMENTS]));
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  // reading Bend-data
  if(Framework[CurrentSystem].NumberOfBendDefinitions>0)
  {
    // alloc memory
    Framework[CurrentSystem].BendDefinitionType=(int*)calloc(Framework[CurrentSystem].NumberOfBendDefinitions,sizeof(int));
    Framework[CurrentSystem].BendDefinitions=(QUAD*)calloc(Framework[CurrentSystem].NumberOfBendDefinitions,sizeof(QUAD));
    Framework[CurrentSystem].BendArgumentDefinitions=(REAL(*)[MAX_BEND_POTENTIAL_ARGUMENTS])
      calloc(Framework[CurrentSystem].NumberOfBendDefinitions,sizeof(REAL[MAX_BEND_POTENTIAL_ARGUMENTS]));
    Framework[CurrentSystem].NumberOfBendsPerType=(int*)calloc(Framework[CurrentSystem].NumberOfBendDefinitions,sizeof(int));

    for(CurrentFramework=0;CurrentFramework<Framework[CurrentSystem].NumberOfFrameworks;CurrentFramework++)
    {
      Framework[CurrentSystem].MaxNumberOfBends[CurrentFramework]=4096;
      // allocate enough memory, it will be resized when needed
      Framework[CurrentSystem].Bends[CurrentFramework]=(QUAD*)calloc(Framework[CurrentSystem].MaxNumberOfBends[CurrentFramework],sizeof(QUAD));
      Framework[CurrentSystem].BendType[CurrentFramework]=(int*)calloc(Framework[CurrentSystem].MaxNumberOfBends[CurrentFramework],sizeof(int));
      Framework[CurrentSystem].BendArguments[CurrentFramework]=(REAL(*)[MAX_BEND_POTENTIAL_ARGUMENTS])
              calloc(Framework[CurrentSystem].MaxNumberOfBends[CurrentFramework],sizeof(REAL[MAX_BEND_POTENTIAL_ARGUMENTS]));
    }

    ReadLine(line,1024,FilePtr); // skip line
    for(i=0;i<Framework[CurrentSystem].NumberOfBendDefinitions;i++)
    {
      ReadLine(line,1024,FilePtr);
      arg_pointer=line;
      sscanf(line,"%s%s%s%s%n",TypeNameA,TypeNameB,TypeNameC,buffer,&n);
      TypeA=ReturnPseudoAtomNumber(TypeNameA);
      TypeB=ReturnPseudoAtomNumber(TypeNameB);
      TypeC=ReturnPseudoAtomNumber(TypeNameC);

      // determine bend-type
      for(j=0;j<NR_BEND_TYPES;j++)
        if(strncasecmp(BendTypes[j].Name,buffer,MAX2(strlen(BendTypes[j].Name),strlen(buffer)))==0)
          BendType=j;

      // read arguments
      for(j=0;j<BendTypes[BendType].nr_args;j++)
      {
        arg_pointer+=n;
        sscanf(arg_pointer,"%lf%n",&temp,&n);
        arguments[j]=(REAL)temp;
      }

      Framework[CurrentSystem].BendDefinitionType[i]=BendType;
      Framework[CurrentSystem].BendDefinitions[i].A=TypeA;
      Framework[CurrentSystem].BendDefinitions[i].B=TypeB;
      Framework[CurrentSystem].BendDefinitions[i].C=TypeC;
      for(j=0;j<BendTypes[BendType].nr_args;j++)
        Framework[CurrentSystem].BendArgumentDefinitions[i][j]=arguments[j];

      // fill the appropriate bend-data
      for(CurrentFramework=0;CurrentFramework<Framework[CurrentSystem].NumberOfFrameworks;CurrentFramework++)
      {
        index=Framework[CurrentSystem].NumberOfBends[CurrentFramework];
        for(B=0;B<Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];B++)
        {
          CurrentTypeB=Framework[CurrentSystem].Atoms[CurrentFramework][B].Type;
          if(CurrentTypeB==TypeB)
          {
            // look for core-shell bends
            if(Framework[CurrentSystem].NumberOfCoreShells[CurrentFramework]>0)
            {
              for(k=0;k<Framework[CurrentSystem].Connectivity[CurrentFramework][B];k++)
              {
                A=GetNeighbour(CurrentSystem,CurrentFramework,B,k);
                if(Framework[CurrentSystem].CoreShellConnectivity[CurrentFramework][A]>0)
                {
                  // loop not over the core, but its shell
                  A=Framework[CurrentSystem].CoreShellConnectivity[CurrentFramework][A];
                  CurrentTypeA=Framework[CurrentSystem].Atoms[CurrentFramework][A].Type;

                  for(l=0;l<Framework[CurrentSystem].Connectivity[CurrentFramework][B];l++)
                  {
                    C=GetNeighbour(CurrentSystem,CurrentFramework,B,l);
                    if(Framework[CurrentSystem].CoreShellConnectivity[CurrentFramework][C]>0)
                    {
                      // loop not over the core, but its shell
                      C=Framework[CurrentSystem].CoreShellConnectivity[CurrentFramework][C];
                      CurrentTypeC=Framework[CurrentSystem].Atoms[CurrentFramework][C].Type;

                      // A<C removes duplicates
                      if((A!=C)&&((CurrentTypeA==TypeA&&CurrentTypeC==TypeC)||(CurrentTypeA==TypeC&&CurrentTypeC==TypeA)))
                      {
                        present=FALSE;
                        for(n=0;n<index;n++)
                        {
                          if(((Framework[CurrentSystem].Bends[CurrentFramework][n].A==A)&&(Framework[CurrentSystem].Bends[CurrentFramework][n].B==B)&&
                             (Framework[CurrentSystem].Bends[CurrentFramework][n].C==C))||((Framework[CurrentSystem].Bends[CurrentFramework][n].A==C)&&
                             (Framework[CurrentSystem].Bends[CurrentFramework][n].B==B)&&(Framework[CurrentSystem].Bends[CurrentFramework][n].C==A)))
                          {
                            present=TRUE;
                            if(present) break;
                          }
                        }
                        if(!present)
                        {
                          Framework[CurrentSystem].Bends[CurrentFramework][index].A=A;
                          Framework[CurrentSystem].Bends[CurrentFramework][index].B=B;
                          Framework[CurrentSystem].Bends[CurrentFramework][index].C=C;
                          Framework[CurrentSystem].Bends[CurrentFramework][index].D=-1;

                          Framework[CurrentSystem].BendType[CurrentFramework][index]=BendType;
                          Framework[CurrentSystem].NumberOfBendsPerType[i]++;

                          for(j=0;j<BendTypes[BendType].nr_args;j++)
                            Framework[CurrentSystem].BendArguments[CurrentFramework][index][j]=arguments[j];

                          // set to appropriate bend-angle
                          switch(BendType)
                          {
                            case CORE_SHELL_BEND:
                              theta=DEG2RAD*Framework[CurrentSystem].BendArguments[CurrentFramework][index][1];
                              Framework[CurrentSystem].BendArguments[CurrentFramework][index][0]*=KELVIN_TO_ENERGY;
                              Framework[CurrentSystem].BendArguments[CurrentFramework][index][1]=theta;
                              break;
                            default:
                              fprintf(stderr, "Undefined Core-Shell Bend potential in routine 'ReadFrameworkDefinition' ('framework.c')\n");
                              exit(0);
                              break;
                          }
                          Framework[CurrentSystem].NumberOfBends[CurrentFramework]++;
                          index=Framework[CurrentSystem].NumberOfBends[CurrentFramework];

                          if(index>=Framework[CurrentSystem].MaxNumberOfBends[CurrentFramework])
                          {
                            Framework[CurrentSystem].MaxNumberOfBends[CurrentFramework]+=4096;
                            Framework[CurrentSystem].Bends[CurrentFramework]=(QUAD*)realloc(Framework[CurrentSystem].Bends[CurrentFramework],
                                         Framework[CurrentSystem].MaxNumberOfBends[CurrentFramework]*sizeof(QUAD));
                            Framework[CurrentSystem].BendType[CurrentFramework]=(int*)realloc(Framework[CurrentSystem].BendType[CurrentFramework],
                                         Framework[CurrentSystem].MaxNumberOfBends[CurrentFramework]*sizeof(int));
                            Framework[CurrentSystem].BendArguments[CurrentFramework]=(REAL(*)[MAX_BEND_POTENTIAL_ARGUMENTS])
                                realloc(Framework[CurrentSystem].BendArguments[CurrentFramework],
                                Framework[CurrentSystem].MaxNumberOfBends[CurrentFramework]*sizeof(REAL[MAX_BEND_POTENTIAL_ARGUMENTS]));
                          }
                        }
                      }
                    }
                  }
                }
              }
            }


            for(k=0;k<Framework[CurrentSystem].Connectivity[CurrentFramework][B];k++)
            {
              A=GetNeighbour(CurrentSystem,CurrentFramework,B,k);
              CurrentTypeA=Framework[CurrentSystem].Atoms[CurrentFramework][A].Type;

              for(l=0;l<Framework[CurrentSystem].Connectivity[CurrentFramework][B];l++)
              {
                C=GetNeighbour(CurrentSystem,CurrentFramework,B,l);
                CurrentTypeC=Framework[CurrentSystem].Atoms[CurrentFramework][C].Type;

                // A<C removes duplicates
                if((A!=C)&&((CurrentTypeA==TypeA&&CurrentTypeC==TypeC)||(CurrentTypeA==TypeC&&CurrentTypeC==TypeA)))
                {
                  present=FALSE;
                  for(n=0;n<index;n++)
                  {
                    if(((Framework[CurrentSystem].Bends[CurrentFramework][n].A==A)&&(Framework[CurrentSystem].Bends[CurrentFramework][n].B==B)&&
                       (Framework[CurrentSystem].Bends[CurrentFramework][n].C==C))||((Framework[CurrentSystem].Bends[CurrentFramework][n].A==C)&&
                       (Framework[CurrentSystem].Bends[CurrentFramework][n].B==B)&&(Framework[CurrentSystem].Bends[CurrentFramework][n].C==A)))
                    {
                      present=TRUE;
                      if(present) break;
                    }
                  }

                  if(!present)
                  {
                    Framework[CurrentSystem].Bends[CurrentFramework][index].A=A;
                    Framework[CurrentSystem].Bends[CurrentFramework][index].B=B;
                    Framework[CurrentSystem].Bends[CurrentFramework][index].C=C;
                    Framework[CurrentSystem].Bends[CurrentFramework][index].D=-1;

                    Framework[CurrentSystem].BendType[CurrentFramework][index]=BendType;
                    Framework[CurrentSystem].NumberOfBendsPerType[i]++;

                    for(j=0;j<BendTypes[BendType].nr_args;j++)
                      Framework[CurrentSystem].BendArguments[CurrentFramework][index][j]=arguments[j];

                    // set to appropriate bend-angle
                    switch(BendType)
                    {
                      case HARMONIC_BEND:
                        // (1/2)p_0*(theta-p_1)^2
                        // ===============================================
                        // p_0/k_B [K/rad^2]
                        // p_1     [degrees]
                        theta=DEG2RAD*Framework[CurrentSystem].BendArguments[CurrentFramework][index][1];
                        if(fabs(theta)<1e-5)
                        {
                          Rab.x=Framework[CurrentSystem].Atoms[CurrentFramework][A].Position.x-
                                Framework[CurrentSystem].Atoms[CurrentFramework][B].Position.x;
                          Rab.y=Framework[CurrentSystem].Atoms[CurrentFramework][A].Position.y-
                                Framework[CurrentSystem].Atoms[CurrentFramework][B].Position.y;
                          Rab.z=Framework[CurrentSystem].Atoms[CurrentFramework][A].Position.z-
                                Framework[CurrentSystem].Atoms[CurrentFramework][B].Position.z;
                          Rab=ApplyBoundaryCondition(Rab);
                          rab=sqrt(SQR(Rab.x)+SQR(Rab.y)+SQR(Rab.z));
                          Rab.x/=rab;
                          Rab.y/=rab;
                          Rab.z/=rab;

                          Rbc.x=Framework[CurrentSystem].Atoms[CurrentFramework][C].Position.x-
                                Framework[CurrentSystem].Atoms[CurrentFramework][B].Position.x;
                          Rbc.y=Framework[CurrentSystem].Atoms[CurrentFramework][C].Position.y-
                                Framework[CurrentSystem].Atoms[CurrentFramework][B].Position.y;
                          Rbc.z=Framework[CurrentSystem].Atoms[CurrentFramework][C].Position.z-
                                Framework[CurrentSystem].Atoms[CurrentFramework][B].Position.z;
                          Rbc=ApplyBoundaryCondition(Rbc);
                          rbc=sqrt(SQR(Rbc.x)+SQR(Rbc.y)+SQR(Rbc.z));
                          Rbc.x/=rbc;
                          Rbc.y/=rbc;
                          Rbc.z/=rbc;

                          theta=acos(Rab.x*Rbc.x+Rab.y*Rbc.y+Rab.z*Rbc.z);
                        }
                        Framework[CurrentSystem].BendArguments[CurrentFramework][index][0]*=KELVIN_TO_ENERGY;
                        Framework[CurrentSystem].BendArguments[CurrentFramework][index][1]=theta;
                        break;
                      case QUARTIC_BEND:
                        // (1/2)p_0*(theta-p_1)^2+(1/3)*p_2*(theta-p_1)^3+(1/4)*p_2*(theta-p_1)^4
                        // ======================================================================
                        // p_0/k_B [K/rad^2]
                        // p_1     [degrees]
                        // p_2/k_B [K/rad^3]
                        // p_3/k_B [K/rad^4]
                      case CFF_QUARTIC_BEND:
                        // p_0*(theta-p_1)^2+p_2*(theta-p_1)^3+p_3*(theta-p_1)^4
                        // =====================================================
                        // p_0/k_B [K/rad^2]
                        // p_1     [degrees]
                        // p_2/k_B [K/rad^3]
                        // p_3/k_B [K/rad^4]
                        theta=DEG2RAD*Framework[CurrentSystem].BendArguments[CurrentFramework][index][1];
                        if(fabs(theta)<1e-5)
                        {
                          Rab.x=Framework[CurrentSystem].Atoms[CurrentFramework][A].Position.x-
                                Framework[CurrentSystem].Atoms[CurrentFramework][B].Position.x;
                          Rab.y=Framework[CurrentSystem].Atoms[CurrentFramework][A].Position.y-
                                Framework[CurrentSystem].Atoms[CurrentFramework][B].Position.y;
                          Rab.z=Framework[CurrentSystem].Atoms[CurrentFramework][A].Position.z-
                                Framework[CurrentSystem].Atoms[CurrentFramework][B].Position.z;
                          Rab=ApplyBoundaryCondition(Rab);
                          rab=sqrt(SQR(Rab.x)+SQR(Rab.y)+SQR(Rab.z));
                          Rab.x/=rab;
                          Rab.y/=rab;
                          Rab.z/=rab;

                          Rbc.x=Framework[CurrentSystem].Atoms[CurrentFramework][C].Position.x-
                                Framework[CurrentSystem].Atoms[CurrentFramework][B].Position.x;
                          Rbc.y=Framework[CurrentSystem].Atoms[CurrentFramework][C].Position.y-
                                Framework[CurrentSystem].Atoms[CurrentFramework][B].Position.y;
                          Rbc.z=Framework[CurrentSystem].Atoms[CurrentFramework][C].Position.z-
                                Framework[CurrentSystem].Atoms[CurrentFramework][B].Position.z;
                          Rbc=ApplyBoundaryCondition(Rbc);
                          rbc=sqrt(SQR(Rbc.x)+SQR(Rbc.y)+SQR(Rbc.z));
                          Rbc.x=Rbc.x/rbc;
                          Rbc.y=Rbc.y/rbc;
                          Rbc.z=Rbc.z/rbc;

                          theta=acos(Rab.x*Rbc.x+Rab.y*Rbc.y+Rab.z*Rbc.z);
                        }
                        Framework[CurrentSystem].BendArguments[CurrentFramework][index][0]*=KELVIN_TO_ENERGY;
                        Framework[CurrentSystem].BendArguments[CurrentFramework][index][1]=theta;
                        Framework[CurrentSystem].BendArguments[CurrentFramework][index][2]*=KELVIN_TO_ENERGY;
                        Framework[CurrentSystem].BendArguments[CurrentFramework][index][3]*=KELVIN_TO_ENERGY;
                        break;
                     case HARMONIC_COSINE_BEND:
                        // (1/2)*p_0*(cos(theta)-cos(p_1))^2
                        // ===============================================
                        // p_0/k_B [K]
                        // p_1     [degrees]
                        theta=cos(DEG2RAD*Framework[CurrentSystem].BendArguments[CurrentFramework][index][1]);
                        if(fabs(theta)<1e-5)
                        {
                          Rab.x=Framework[CurrentSystem].Atoms[CurrentFramework][A].Position.x-
                                Framework[CurrentSystem].Atoms[CurrentFramework][B].Position.x;
                          Rab.y=Framework[CurrentSystem].Atoms[CurrentFramework][A].Position.y-
                                Framework[CurrentSystem].Atoms[CurrentFramework][B].Position.y;
                          Rab.z=Framework[CurrentSystem].Atoms[CurrentFramework][A].Position.z-
                                Framework[CurrentSystem].Atoms[CurrentFramework][B].Position.z;
                          Rab=ApplyBoundaryCondition(Rab);
                          rab=sqrt(SQR(Rab.x)+SQR(Rab.y)+SQR(Rab.z));
                          Rab.x/=rab;
                          Rab.y/=rab;
                          Rab.z/=rab;

                          Rbc.x=Framework[CurrentSystem].Atoms[CurrentFramework][C].Position.x-
                                Framework[CurrentSystem].Atoms[CurrentFramework][B].Position.x;
                          Rbc.y=Framework[CurrentSystem].Atoms[CurrentFramework][C].Position.y-
                                Framework[CurrentSystem].Atoms[CurrentFramework][B].Position.y;
                          Rbc.z=Framework[CurrentSystem].Atoms[CurrentFramework][C].Position.z-
                                Framework[CurrentSystem].Atoms[CurrentFramework][B].Position.z;
                          Rbc=ApplyBoundaryCondition(Rbc);
                          rbc=sqrt(SQR(Rbc.x)+SQR(Rbc.y)+SQR(Rbc.z));
                          Rbc.x/=rbc;
                          Rbc.y/=rbc;
                          Rbc.z/=rbc;

                          theta=Rab.x*Rbc.x+Rab.y*Rbc.y+Rab.z*Rbc.z;
                        }
                        Framework[CurrentSystem].BendArguments[CurrentFramework][index][0]*=KELVIN_TO_ENERGY;
                        Framework[CurrentSystem].BendArguments[CurrentFramework][index][1]=theta;
                        break;
                      case COSINE_BEND:
                        // p_0*(1+cos(p_1*theta-p_2))
                        // ===============================================
                        // p_0/k_B [K]
                        // p_1     [-]
                        // p_2     [degrees]
                        theta=DEG2RAD*Framework[CurrentSystem].BendArguments[CurrentFramework][index][2];
                        if(fabs(theta)<1e-5)
                        {
                          Rab.x=Framework[CurrentSystem].Atoms[CurrentFramework][A].Position.x-
                                Framework[CurrentSystem].Atoms[CurrentFramework][B].Position.x;
                          Rab.y=Framework[CurrentSystem].Atoms[CurrentFramework][A].Position.y-
                                Framework[CurrentSystem].Atoms[CurrentFramework][B].Position.y;
                          Rab.z=Framework[CurrentSystem].Atoms[CurrentFramework][A].Position.z-
                                Framework[CurrentSystem].Atoms[CurrentFramework][B].Position.z;
                          Rab=ApplyBoundaryCondition(Rab);
                          rab=sqrt(SQR(Rab.x)+SQR(Rab.y)+SQR(Rab.z));
                          Rab.x/=rab;
                          Rab.y/=rab;
                          Rab.z/=rab;

                          Rbc.x=Framework[CurrentSystem].Atoms[CurrentFramework][C].Position.x-
                                Framework[CurrentSystem].Atoms[CurrentFramework][B].Position.x;
                          Rbc.y=Framework[CurrentSystem].Atoms[CurrentFramework][C].Position.y-
                                Framework[CurrentSystem].Atoms[CurrentFramework][B].Position.y;
                          Rbc.z=Framework[CurrentSystem].Atoms[CurrentFramework][C].Position.z-
                                Framework[CurrentSystem].Atoms[CurrentFramework][B].Position.z;
                          Rbc=ApplyBoundaryCondition(Rbc);
                          rbc=sqrt(SQR(Rbc.x)+SQR(Rbc.y)+SQR(Rbc.z));
                          Rbc.x/=rbc;
                          Rbc.y/=rbc;
                          Rbc.z/=rbc;

                          theta=Rab.x*Rbc.x+Rab.y*Rbc.y+Rab.z*Rbc.z;
                        }
                        Framework[CurrentSystem].BendArguments[CurrentFramework][index][0]*=KELVIN_TO_ENERGY;
                        Framework[CurrentSystem].BendArguments[CurrentFramework][index][2]=theta;
                        break;
                      case TAFIPOLSKY_BEND:
                        // 0.5*p_0*(1+cos(theta))*(1+cos(2*theta))
                        // ===============================================
                        // p_0/k_B [K]
                        Framework[CurrentSystem].BendArguments[CurrentFramework][index][0]*=KELVIN_TO_ENERGY;
                        break;
                      case MM3_BEND:
                      case MM3_IN_PLANE_BEND:
                        // p_0*(theta-p_1)^2(1-0.014*(theta-p_1)+5.6e-5*(theta-p_1)^2-7e-7*(theta-p_1)^3+2.2e-8(theta-p_1)^4)
                        // =================================================================================================
                        // p_0/k_B [mdyne A/rad^2]
                        // p_1     [degrees]
                        theta=DEG2RAD*Framework[CurrentSystem].BendArguments[CurrentFramework][index][1];
                        if(fabs(theta)<1e-5)
                        {
                          Rab.x=Framework[CurrentSystem].Atoms[CurrentFramework][A].Position.x-
                                Framework[CurrentSystem].Atoms[CurrentFramework][B].Position.x;
                          Rab.y=Framework[CurrentSystem].Atoms[CurrentFramework][A].Position.y-
                                Framework[CurrentSystem].Atoms[CurrentFramework][B].Position.y;
                          Rab.z=Framework[CurrentSystem].Atoms[CurrentFramework][A].Position.z-
                                Framework[CurrentSystem].Atoms[CurrentFramework][B].Position.z;
                          Rab=ApplyBoundaryCondition(Rab);
                          rab=sqrt(SQR(Rab.x)+SQR(Rab.y)+SQR(Rab.z));
                          Rab.x/=rab;
                          Rab.y/=rab;
                          Rab.z/=rab;

                          Rbc.x=Framework[CurrentSystem].Atoms[CurrentFramework][C].Position.x-
                                Framework[CurrentSystem].Atoms[CurrentFramework][B].Position.x;
                          Rbc.y=Framework[CurrentSystem].Atoms[CurrentFramework][C].Position.y-
                                Framework[CurrentSystem].Atoms[CurrentFramework][B].Position.y;
                          Rbc.z=Framework[CurrentSystem].Atoms[CurrentFramework][C].Position.z-
                                Framework[CurrentSystem].Atoms[CurrentFramework][B].Position.z;
                          Rbc=ApplyBoundaryCondition(Rbc);
                          rbc=sqrt(SQR(Rbc.x)+SQR(Rbc.y)+SQR(Rbc.z));
                          Rbc.x=Rbc.x/rbc;
                          Rbc.y=Rbc.y/rbc;
                          Rbc.z=Rbc.z/rbc;

                          theta=acos(Rab.x*Rbc.x+Rab.y*Rbc.y+Rab.z*Rbc.z);
                        }
                        // MM3 input is in mdyne A/rad^2
                        // the MM3 force field uses the factor 0.021914 to convert to kcal/mole
                        Framework[CurrentSystem].BendArguments[CurrentFramework][index][0]*=0.021914*KCAL_PER_MOL_TO_ENERGY;
                        Framework[CurrentSystem].BendArguments[CurrentFramework][index][1]=theta;
                        break;
                      case MEASURE_BEND:
                        break;
                      default:
                        fprintf(stderr, "Undefined Bend potential in routine 'ReadFrameworkDefinition' ('framework.c')\n");
                        exit(0);
                        break;
                    }
                    Framework[CurrentSystem].NumberOfBends[CurrentFramework]++;
                    index=Framework[CurrentSystem].NumberOfBends[CurrentFramework];

                    if(index>=Framework[CurrentSystem].MaxNumberOfBends[CurrentFramework])
                    {
                      Framework[CurrentSystem].MaxNumberOfBends[CurrentFramework]+=4096;
                      Framework[CurrentSystem].Bends[CurrentFramework]=(QUAD*)realloc(Framework[CurrentSystem].Bends[CurrentFramework],
                                  Framework[CurrentSystem].MaxNumberOfBends[CurrentFramework]*sizeof(QUAD));
                      Framework[CurrentSystem].BendType[CurrentFramework]=(int*)realloc(Framework[CurrentSystem].BendType[CurrentFramework],
                                  Framework[CurrentSystem].MaxNumberOfBends[CurrentFramework]*sizeof(int));
                      Framework[CurrentSystem].BendArguments[CurrentFramework]=(REAL(*)[MAX_BEND_POTENTIAL_ARGUMENTS])
                          realloc(Framework[CurrentSystem].BendArguments[CurrentFramework],
                          Framework[CurrentSystem].MaxNumberOfBends[CurrentFramework]*sizeof(REAL[MAX_BEND_POTENTIAL_ARGUMENTS]));
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  // reading InversionBend-data
  index=0;
  if(Framework[CurrentSystem].NumberOfInversionBendDefinitions>0)
  {
    // alloc memory
    Framework[CurrentSystem].InversionBendDefinitionType=(int*)calloc(Framework[CurrentSystem].NumberOfInversionBendDefinitions,sizeof(int));
    Framework[CurrentSystem].InversionBendDefinitions=(QUAD*)calloc(Framework[CurrentSystem].NumberOfInversionBendDefinitions,sizeof(QUAD));
    Framework[CurrentSystem].InversionBendArgumentDefinitions=(REAL(*)[MAX_INVERSION_BEND_POTENTIAL_ARGUMENTS])
      calloc(Framework[CurrentSystem].NumberOfInversionBendDefinitions,sizeof(REAL[MAX_INVERSION_BEND_POTENTIAL_ARGUMENTS]));
    Framework[CurrentSystem].NumberOfInversionBendsPerType=(int*)calloc(Framework[CurrentSystem].NumberOfInversionBendDefinitions,sizeof(int));

    for(CurrentFramework=0;CurrentFramework<Framework[CurrentSystem].NumberOfFrameworks;CurrentFramework++)
    {
      Framework[CurrentSystem].MaxNumberOfInversionBends[CurrentFramework]=4096;
      // allocate enough memory, it will be resized when needed
      Framework[CurrentSystem].InversionBends[CurrentFramework]=(QUAD*)calloc(Framework[CurrentSystem].MaxNumberOfInversionBends[CurrentFramework],sizeof(QUAD));
      Framework[CurrentSystem].InversionBendType[CurrentFramework]=(int*)calloc(Framework[CurrentSystem].MaxNumberOfInversionBends[CurrentFramework],sizeof(int));
      Framework[CurrentSystem].InversionBendArguments[CurrentFramework]=(REAL(*)[MAX_INVERSION_BEND_POTENTIAL_ARGUMENTS])
              calloc(Framework[CurrentSystem].MaxNumberOfInversionBends[CurrentFramework],sizeof(REAL[MAX_INVERSION_BEND_POTENTIAL_ARGUMENTS]));
    }

    ReadLine(line,1024,FilePtr); // skip line
    for(i=0;i<Framework[CurrentSystem].NumberOfInversionBendDefinitions;i++)
    {
      ReadLine(line,1024,FilePtr);
      arg_pointer=line;
      sscanf(line,"%s%s%s%s%s%n",TypeNameA,TypeNameB,TypeNameC,TypeNameD,buffer,&n);
      TypeA=ReturnPseudoAtomNumber(TypeNameA);
      TypeB=ReturnPseudoAtomNumber(TypeNameB);
      TypeC=ReturnPseudoAtomNumber(TypeNameC);
      TypeD=ReturnPseudoAtomNumber(TypeNameD);

      // determine inversionbend-type
      for(j=0;j<NR_INVERSION_BEND_TYPES;j++)
        if(strcasecmp(InversionBendTypes[j].Name,buffer)==0)
          InversionBendType=j;

      // read arguments
      for(j=0;j<InversionBendTypes[InversionBendType].nr_args;j++)
      {
        arg_pointer+=n;
        sscanf(arg_pointer,"%lf%n",&temp,&n);
        arguments[j]=(REAL)temp;
      }

      Framework[CurrentSystem].InversionBendDefinitionType[i]=InversionBendType;
      Framework[CurrentSystem].InversionBendDefinitions[i].A=TypeA;
      Framework[CurrentSystem].InversionBendDefinitions[i].B=TypeB;
      Framework[CurrentSystem].InversionBendDefinitions[i].C=TypeC;
      Framework[CurrentSystem].InversionBendDefinitions[i].D=TypeD;

      for(j=0;j<InversionBendTypes[InversionBendType].nr_args;j++)
        Framework[CurrentSystem].InversionBendArgumentDefinitions[i][j]=arguments[j];

      // fill the appropriate inversion bend-data
      for(CurrentFramework=0;CurrentFramework<Framework[CurrentSystem].NumberOfFrameworks;CurrentFramework++)
      {
        index=Framework[CurrentSystem].NumberOfInversionBends[CurrentFramework];
        for(A=0;A<Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];A++)
        {
          CurrentTypeA=Framework[CurrentSystem].Atoms[CurrentFramework][A].Type;

          if(CurrentTypeA==TypeA)
          {
            // atom A is of type 'TypeA'
            for(k=0;k<Framework[CurrentSystem].Connectivity[CurrentFramework][A];k++)
            {
              B=GetNeighbour(CurrentSystem,CurrentFramework,A,k);
              CurrentTypeB=Framework[CurrentSystem].Atoms[CurrentFramework][B].Type;

              if(CurrentTypeB==TypeB)
              {
                // atom B is of type 'TypeB'
                for(l=0;l<Framework[CurrentSystem].Connectivity[CurrentFramework][B];l++)
                {
                  C=GetNeighbour(CurrentSystem,CurrentFramework,B,l);
                  CurrentTypeC=Framework[CurrentSystem].Atoms[CurrentFramework][C].Type;
                  if((A!=C)&&(CurrentTypeC==TypeC))
                  {
                    // atom C is of type 'TypeC' and not equal to atom A
                    for(m=0;m<Framework[CurrentSystem].Connectivity[CurrentFramework][B];m++)
                    {
                      D=GetNeighbour(CurrentSystem,CurrentFramework,B,m);
                      CurrentTypeD=Framework[CurrentSystem].Atoms[CurrentFramework][D].Type;
                      if((D!=C)&&(D!=A)&&(CurrentTypeD==TypeD))
                      {
                        // atom D is of type 'TypeD' and not equal to atom A nor to atom C
                        present=FALSE;
                        for(n=0;n<index;n++)
                        {
                          ListTypeC=Framework[CurrentSystem].Atoms[CurrentFramework][Framework[CurrentSystem].InversionBends[CurrentFramework][n].C].Type;
                          ListTypeD=Framework[CurrentSystem].Atoms[CurrentFramework][Framework[CurrentSystem].InversionBends[CurrentFramework][n].D].Type;
                          if((Framework[CurrentSystem].InversionBends[CurrentFramework][n].A==A)&&(Framework[CurrentSystem].InversionBends[CurrentFramework][n].B==B)&&
                             (((ListTypeC==CurrentTypeC)&&(ListTypeD==CurrentTypeD))||((ListTypeC==CurrentTypeD)&&(ListTypeD==CurrentTypeC))))
                          {
                            present=TRUE;
                            if(present) break;
                          }
                        }

                        if(!present)
                        {
                          Framework[CurrentSystem].InversionBends[CurrentFramework][index].A=A;
                          Framework[CurrentSystem].InversionBends[CurrentFramework][index].B=B;
                          Framework[CurrentSystem].InversionBends[CurrentFramework][index].C=C;
                          Framework[CurrentSystem].InversionBends[CurrentFramework][index].D=D;

                          Framework[CurrentSystem].InversionBendType[CurrentFramework][index]=InversionBendType;
                          Framework[CurrentSystem].NumberOfInversionBendsPerType[i]++;

                          for(j=0;j<InversionBendTypes[InversionBendType].nr_args;j++)
                            Framework[CurrentSystem].InversionBendArguments[CurrentFramework][index][j]=arguments[j];
                          switch(InversionBendType)
                          {
                            case HARMONIC_INVERSION:
                            case HARMONIC_INVERSION2:
                              // (1/2)*p_0*(chi-p_1)^2
                              // ===============================================
                              // p_0/k_B [K]
                              // p_1     [degrees]
                              Framework[CurrentSystem].InversionBendArguments[CurrentFramework][index][0]*=KELVIN_TO_ENERGY;
                              Framework[CurrentSystem].InversionBendArguments[CurrentFramework][index][1]*=DEG2RAD;
                              break;
                            case HARMONIC_COSINE_INVERSION:
                            case HARMONIC_COSINE_INVERSION2:
                              // (1/2)*p_0*(cos(phi)-cos(p_1))^2
                              // ===============================================
                              // p_0/k_B [K]
                              // p_1     [degrees]
                              Framework[CurrentSystem].InversionBendArguments[CurrentFramework][index][0]*=KELVIN_TO_ENERGY;
                              Framework[CurrentSystem].InversionBendArguments[CurrentFramework][index][1]=cos(DEG2RAD*Framework[CurrentSystem].InversionBendArguments[CurrentFramework][index][1]);
                              break;
                            case PLANAR_INVERSION:
                            case PLANAR_INVERSION2:
                              // (1/2)*p_0*(1-cos(phi))
                              // ===============================================
                              // p_0/k_B [K]
                              Framework[CurrentSystem].InversionBendArguments[CurrentFramework][index][0]*=KELVIN_TO_ENERGY;
                              break;
                            case MM3_INVERSION:
                              // p_0*(chi-p_1)^2(1-0.014*(chi-p_1)+5.6e-5*(chi-p_1)^2-7e-7*(chi-p_1)^3+2.2e-8(chi-p_1)^4)
                              // =================================================================================================
                              // p_0/k_B [mdyne A/rad^2]
                              // p_1     [degrees]
                              Framework[CurrentSystem].InversionBendArguments[CurrentFramework][index][0]*=0.02191418*KCAL_PER_MOL_TO_ENERGY;
                              Framework[CurrentSystem].InversionBendArguments[CurrentFramework][index][1]*=DEG2RAD;
                              break;
                            default:
                              fprintf(stderr, "Undefined Inversion Bend potential in routine 'ReadFrameworkDefinition' ('framework.c')\n");
                              exit(0);
                              break;
                          }
                          Framework[CurrentSystem].NumberOfInversionBends[CurrentFramework]++;
                          index=Framework[CurrentSystem].NumberOfInversionBends[CurrentFramework];

                          if(index>=Framework[CurrentSystem].MaxNumberOfInversionBends[CurrentFramework])
                          {
                            Framework[CurrentSystem].MaxNumberOfInversionBends[CurrentFramework]+=4096;
                            Framework[CurrentSystem].InversionBends[CurrentFramework]=(QUAD*)realloc(Framework[CurrentSystem].InversionBends[CurrentFramework],
                                        Framework[CurrentSystem].MaxNumberOfInversionBends[CurrentFramework]*sizeof(QUAD));
                            Framework[CurrentSystem].InversionBendType[CurrentFramework]=(int*)realloc(Framework[CurrentSystem].InversionBendType[CurrentFramework],
                                        Framework[CurrentSystem].MaxNumberOfInversionBends[CurrentFramework]*sizeof(int));
                            Framework[CurrentSystem].InversionBendArguments[CurrentFramework]=(REAL(*)[MAX_INVERSION_BEND_POTENTIAL_ARGUMENTS])
                                 realloc(Framework[CurrentSystem].InversionBendArguments[CurrentFramework],
                                 Framework[CurrentSystem].MaxNumberOfInversionBends[CurrentFramework]*sizeof(REAL[MAX_INVERSION_BEND_POTENTIAL_ARGUMENTS]));
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }


  // determine the fourth atom for in-plane bends by searching the bends in the list of inversion-bends
  for(CurrentFramework=0;CurrentFramework<Framework[CurrentSystem].NumberOfFrameworks;CurrentFramework++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfInversionBends[CurrentFramework];i++)
    {
      A=Framework[CurrentSystem].InversionBends[CurrentFramework][i].A;
      B=Framework[CurrentSystem].InversionBends[CurrentFramework][i].B;
      C=Framework[CurrentSystem].InversionBends[CurrentFramework][i].C;
      D=Framework[CurrentSystem].InversionBends[CurrentFramework][i].D;

      for(j=0;j<Framework[CurrentSystem].NumberOfBends[CurrentFramework];j++)
      {
        if(Framework[CurrentSystem].BendType[CurrentFramework][j]==MM3_IN_PLANE_BEND)
        {
          A1=Framework[CurrentSystem].Bends[CurrentFramework][j].A;
          B1=Framework[CurrentSystem].Bends[CurrentFramework][j].B;
          C1=Framework[CurrentSystem].Bends[CurrentFramework][j].C;

          if(B1==B)
          {
            if(((A1==A)&&(C1==C))||((A1==C)&&(C1==A))) Framework[CurrentSystem].Bends[CurrentFramework][j].D=D;
            if(((A1==A)&&(C1==D))||((A1==D)&&(C1==A))) Framework[CurrentSystem].Bends[CurrentFramework][j].D=C;
            if(((A1==C)&&(C1==D))||((A1==D)&&(C1==C))) Framework[CurrentSystem].Bends[CurrentFramework][j].D=A;
          }
        }
      }
    }
  }


  // reading Torsion-data
  // Note: Torsions on the same four atoms but different potential and/or arguments are consider different
  if(Framework[CurrentSystem].NumberOfTorsionDefinitions>0)
  {
    // alloc memory
    Framework[CurrentSystem].TorsionDefinitionType=(int*)calloc(Framework[CurrentSystem].NumberOfTorsionDefinitions,sizeof(int));
    Framework[CurrentSystem].TorsionDefinitions=(QUAD*)calloc(Framework[CurrentSystem].NumberOfTorsionDefinitions,sizeof(QUAD));
    Framework[CurrentSystem].TorsionArgumentDefinitions=(REAL(*)[MAX_TORSION_POTENTIAL_ARGUMENTS])
      calloc(Framework[CurrentSystem].NumberOfTorsionDefinitions,sizeof(REAL[MAX_TORSION_POTENTIAL_ARGUMENTS]));
    Framework[CurrentSystem].NumberOfTorsionsPerType=(int*)calloc(Framework[CurrentSystem].NumberOfTorsionDefinitions,sizeof(int));

    for(CurrentFramework=0;CurrentFramework<Framework[CurrentSystem].NumberOfFrameworks;CurrentFramework++)
    {
      Framework[CurrentSystem].MaxNumberOfTorsions[CurrentFramework]=4096;
      // allocate enough memory, it will be resized when needed
      Framework[CurrentSystem].Torsions[CurrentFramework]=(QUAD*)calloc(Framework[CurrentSystem].MaxNumberOfTorsions[CurrentFramework],sizeof(QUAD));
      Framework[CurrentSystem].TorsionType[CurrentFramework]=(int*)calloc(Framework[CurrentSystem].MaxNumberOfTorsions[CurrentFramework],sizeof(int));
      Framework[CurrentSystem].TorsionArguments[CurrentFramework]=(REAL(*)[MAX_TORSION_POTENTIAL_ARGUMENTS])
              calloc(Framework[CurrentSystem].MaxNumberOfTorsions[CurrentFramework],sizeof(REAL[MAX_TORSION_POTENTIAL_ARGUMENTS]));
    }

    ReadLine(line,1024,FilePtr); // skip line
    for(i=0;i<Framework[CurrentSystem].NumberOfTorsionDefinitions;i++)
    {
      ReadLine(line,1024,FilePtr);
      arg_pointer=line;
      sscanf(line,"%s%s%s%s%s%n",TypeNameA,TypeNameB,TypeNameC,TypeNameD,buffer,&n);
      TypeA=ReturnPseudoAtomNumber(TypeNameA);
      TypeB=ReturnPseudoAtomNumber(TypeNameB);
      TypeC=ReturnPseudoAtomNumber(TypeNameC);
      TypeD=ReturnPseudoAtomNumber(TypeNameD);

      // determine torsion-type
      for(j=0;j<NR_TORSION_TYPES;j++)
      {
        if(strcasecmp(TorsionTypes[j].Name,buffer)==0)
          TorsionType=j;
      }

      // read arguments
      arg_pointer+=n;
      for(j=0;j<TorsionTypes[TorsionType].nr_args+2;j++)
        arguments[j]=0;

      for(j=0;j<TorsionTypes[TorsionType].nr_args+2;j++)
      {
        arguments[j]=0;
        temp=0.0;
        if(sscanf(arg_pointer,"%lf%n",&temp,&n)>=1)
        {
          arguments[j]=(REAL)temp;
          arg_pointer+=n;
        }
      }

      Framework[CurrentSystem].TorsionDefinitionType[i]=TorsionType;
      Framework[CurrentSystem].TorsionDefinitions[i].A=TypeA;
      Framework[CurrentSystem].TorsionDefinitions[i].B=TypeB;
      Framework[CurrentSystem].TorsionDefinitions[i].C=TypeC;
      Framework[CurrentSystem].TorsionDefinitions[i].D=TypeD;

      for(j=0;j<TorsionTypes[TorsionType].nr_args;j++)
        Framework[CurrentSystem].TorsionArgumentDefinitions[i][j]=arguments[j];
      Framework[CurrentSystem].TorsionArgumentDefinitions[i][6]=arguments[TorsionTypes[TorsionType].nr_args];
      Framework[CurrentSystem].TorsionArgumentDefinitions[i][7]=arguments[TorsionTypes[TorsionType].nr_args+1];

      // fill the appropriate torsion-data
      for(CurrentFramework=0;CurrentFramework<Framework[CurrentSystem].NumberOfFrameworks;CurrentFramework++)
      {
        index=Framework[CurrentSystem].NumberOfTorsions[CurrentFramework];
        for(B=0;B<Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];B++)
        {
          CurrentTypeB=Framework[CurrentSystem].Atoms[CurrentFramework][B].Type;

          if(CurrentTypeB==TypeB)
          {
            for(k=0;k<Framework[CurrentSystem].Connectivity[CurrentFramework][B];k++)
            {
              C=GetNeighbour(CurrentSystem,CurrentFramework,B,k);
              CurrentTypeC=Framework[CurrentSystem].Atoms[CurrentFramework][C].Type;

              if(CurrentTypeC==TypeC)
              {
                for(l=0;l<Framework[CurrentSystem].Connectivity[CurrentFramework][B];l++)
                {
                  A=GetNeighbour(CurrentSystem,CurrentFramework,B,l);
                  CurrentTypeA=Framework[CurrentSystem].Atoms[CurrentFramework][A].Type;
                  if((A!=C)&&(CurrentTypeA==TypeA))
                  {
                    for(m=0;m<Framework[CurrentSystem].Connectivity[CurrentFramework][C];m++)
                    {
                      D=GetNeighbour(CurrentSystem,CurrentFramework,C,m);
                      CurrentTypeD=Framework[CurrentSystem].Atoms[CurrentFramework][D].Type;
                      if((B!=D)&&(A!=D)&&(C!=D)&&(CurrentTypeD==TypeD))
                      {
                        switch(TorsionType)
                        {
                          case HARMONIC_DIHEDRAL:
                            // (1/2)*p_0*(phi-p_1)^2
                            // ===============================================
                            // p_0/k_B [K/rad^2]
                            // p_1     [degrees]
                            arguments[0]=Framework[CurrentSystem].TorsionArgumentDefinitions[i][0]*KELVIN_TO_ENERGY;
                            arguments[1]=Framework[CurrentSystem].TorsionArgumentDefinitions[i][1];
                            break;
                          case HARMONIC_COSINE_DIHEDRAL:
                            // (1/2)*p_0*(cos(phi)-cos(p_1))^2
                            // ===============================================
                            // p_0/k_B [K]
                            // p_1     [degrees]
                            arguments[0]=Framework[CurrentSystem].TorsionArgumentDefinitions[i][0]*KELVIN_TO_ENERGY;
                            arguments[1]=cos(RAD2DEG*Framework[CurrentSystem].TorsionArgumentDefinitions[i][1]);
                            break;
                          case THREE_COSINE_DIHEDRAL:
                            // (1/2)*p_0*(1+cos(phi))+(1/2)*p_1*(1-cos(2*phi))+(1/2)*p_2*(1+cos(3*phi))
                            // ========================================================================
                            // p_0/k_B [K]
                            // p_1/k_B [K]
                            // p_2/k_B [K]
                            arguments[0]=Framework[CurrentSystem].TorsionArgumentDefinitions[i][0]*KELVIN_TO_ENERGY;
                            arguments[1]=Framework[CurrentSystem].TorsionArgumentDefinitions[i][1]*KELVIN_TO_ENERGY;
                            arguments[2]=Framework[CurrentSystem].TorsionArgumentDefinitions[i][2]*KELVIN_TO_ENERGY;
                            break;
                          case MM3_DIHEDRAL:
                            // (1/2)*p_0*(1+cos(phi))+(1/2)*p_1*(1-cos(2*phi))+(1/2)*p_2*(1+cos(3*phi))
                            // ========================================================================
                            // p_0     [kcal/mol]
                            // p_1     [kcal/mol]
                            // p_2     [kcal/mol]
                            arguments[0]=Framework[CurrentSystem].TorsionArgumentDefinitions[i][0]*KCAL_PER_MOL_TO_ENERGY;
                            arguments[1]=Framework[CurrentSystem].TorsionArgumentDefinitions[i][1]*KCAL_PER_MOL_TO_ENERGY;
                            arguments[2]=Framework[CurrentSystem].TorsionArgumentDefinitions[i][2]*KCAL_PER_MOL_TO_ENERGY;
                            break;
                          case CFF_DIHEDRAL:
                            // p_0*(1-cos(phi))+p_1*(1-cos(2*phi))+p_2*(1-cos(3*phi))
                            // ======================================================
                            // p_0/k_B [K]
                            // p_1/k_B [K]
                            // p_2/k_B [K]
                            arguments[0]=Framework[CurrentSystem].TorsionArgumentDefinitions[i][0]*KELVIN_TO_ENERGY;
                            arguments[1]=Framework[CurrentSystem].TorsionArgumentDefinitions[i][1]*KELVIN_TO_ENERGY;
                            arguments[2]=Framework[CurrentSystem].TorsionArgumentDefinitions[i][2]*KELVIN_TO_ENERGY;
                            break;
                          case CFF_DIHEDRAL2:
                            // p_0*(1+cos(phi))+p_1*(1+cos(2*phi))+p_2*(1+cos(3*phi))
                            // ======================================================
                            // p_0/k_B [K]
                            // p_1/k_B [K]
                            // p_2/k_B [K]
                            arguments[0]=Framework[CurrentSystem].TorsionArgumentDefinitions[i][0]*KELVIN_TO_ENERGY;
                            arguments[1]=Framework[CurrentSystem].TorsionArgumentDefinitions[i][1]*KELVIN_TO_ENERGY;
                            arguments[2]=Framework[CurrentSystem].TorsionArgumentDefinitions[i][2]*KELVIN_TO_ENERGY;
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
                            arguments[0]=Framework[CurrentSystem].TorsionArgumentDefinitions[i][0]*KELVIN_TO_ENERGY;
                            arguments[1]=Framework[CurrentSystem].TorsionArgumentDefinitions[i][1]*KELVIN_TO_ENERGY;
                            arguments[2]=Framework[CurrentSystem].TorsionArgumentDefinitions[i][2]*KELVIN_TO_ENERGY;
                            arguments[3]=Framework[CurrentSystem].TorsionArgumentDefinitions[i][3]*KELVIN_TO_ENERGY;
                            arguments[4]=Framework[CurrentSystem].TorsionArgumentDefinitions[i][4]*KELVIN_TO_ENERGY;
                            arguments[5]=Framework[CurrentSystem].TorsionArgumentDefinitions[i][5]*KELVIN_TO_ENERGY;
                            break;
                          case TRAPPE_DIHEDRAL:
                            // p_0[0]+p_1*(1+cos(phi))+p_2*(1-cos(2*phi))+p_3*(1+cos(3*phi))
                            // =============================================================
                            // p_0/k_B [K]
                            // p_1/k_B [K]
                            // p_2/k_B [K]
                            // p_3/k_B [K]
                            arguments[0]=Framework[CurrentSystem].TorsionArgumentDefinitions[i][0]*KELVIN_TO_ENERGY;
                            arguments[1]=Framework[CurrentSystem].TorsionArgumentDefinitions[i][1]*KELVIN_TO_ENERGY;
                            arguments[2]=Framework[CurrentSystem].TorsionArgumentDefinitions[i][2]*KELVIN_TO_ENERGY;
                            arguments[3]=Framework[CurrentSystem].TorsionArgumentDefinitions[i][3]*KELVIN_TO_ENERGY;
                            break;
                          case CVFF_DIHEDRAL:
                            // p_0*(1+cos(p_1*phi-p_2))
                            // ========================
                            // p_0/k_B [K]
                            // p_1     [-]
                            // p_2     [degrees]
                            arguments[0]=Framework[CurrentSystem].TorsionArgumentDefinitions[i][0]*KELVIN_TO_ENERGY;
                            arguments[1]=Framework[CurrentSystem].TorsionArgumentDefinitions[i][1];
                            arguments[2]=Framework[CurrentSystem].TorsionArgumentDefinitions[i][2]*DEG2RAD;
                            break;
                          case OPLS_DIHEDRAL:
                            // (1/2)p_0[0]+(1/2)p_1*(1+cos(phi))+(1/2)p_2*(1-cos(2*phi))+(1/2)p_3*(1+cos(3*phi))
                            // =================================================================================
                            // p_0/k_B [K]
                            // p_1/k_B [K]
                            // p_2/k_B [K]
                            // p_3/k_B [K]
                            arguments[0]=Framework[CurrentSystem].TorsionArgumentDefinitions[i][0]*KELVIN_TO_ENERGY;
                            arguments[1]=Framework[CurrentSystem].TorsionArgumentDefinitions[i][1]*KELVIN_TO_ENERGY;
                            arguments[2]=Framework[CurrentSystem].TorsionArgumentDefinitions[i][2]*KELVIN_TO_ENERGY;
                            arguments[3]=Framework[CurrentSystem].TorsionArgumentDefinitions[i][3]*KELVIN_TO_ENERGY;
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
                            arguments[0]=Framework[CurrentSystem].TorsionArgumentDefinitions[i][0]*KELVIN_TO_ENERGY;
                            arguments[1]=Framework[CurrentSystem].TorsionArgumentDefinitions[i][1]*KELVIN_TO_ENERGY;
                            arguments[2]=Framework[CurrentSystem].TorsionArgumentDefinitions[i][2]*KELVIN_TO_ENERGY;
                            arguments[3]=Framework[CurrentSystem].TorsionArgumentDefinitions[i][3]*KELVIN_TO_ENERGY;
                            arguments[4]=Framework[CurrentSystem].TorsionArgumentDefinitions[i][4]*KELVIN_TO_ENERGY;
                            arguments[5]=Framework[CurrentSystem].TorsionArgumentDefinitions[i][5]*KELVIN_TO_ENERGY;
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
                            arguments[0]=Framework[CurrentSystem].TorsionArgumentDefinitions[i][0]*KELVIN_TO_ENERGY;
                            arguments[1]=Framework[CurrentSystem].TorsionArgumentDefinitions[i][1]*KELVIN_TO_ENERGY;
                            arguments[2]=Framework[CurrentSystem].TorsionArgumentDefinitions[i][2]*KELVIN_TO_ENERGY;
                            arguments[3]=Framework[CurrentSystem].TorsionArgumentDefinitions[i][3]*KELVIN_TO_ENERGY;
                            arguments[4]=Framework[CurrentSystem].TorsionArgumentDefinitions[i][4]*KELVIN_TO_ENERGY;
                            arguments[5]=Framework[CurrentSystem].TorsionArgumentDefinitions[i][5]*KELVIN_TO_ENERGY;
                            break;
                          default:
                            fprintf(stderr, "Undefined Torsion potential in routine 'ReadFrameworkDefinition' ('framework.c')\n");
                            exit(0);
                            break;
                        }

                        // we now have a quad-triple
                        present=FALSE;
                        for(n=0;n<index;n++)
                        {
                          if((((Framework[CurrentSystem].Torsions[CurrentFramework][n].A==A)&&(Framework[CurrentSystem].Torsions[CurrentFramework][n].B==B)&&
                            (Framework[CurrentSystem].Torsions[CurrentFramework][n].C==C)&&(Framework[CurrentSystem].Torsions[CurrentFramework][n].D==D))||
                            ((Framework[CurrentSystem].Torsions[CurrentFramework][n].A==D)&&(Framework[CurrentSystem].Torsions[CurrentFramework][n].B==C)&&
                            (Framework[CurrentSystem].Torsions[CurrentFramework][n].C==B)&&(Framework[CurrentSystem].Torsions[CurrentFramework][n].D==A)))&&
                            (Framework[CurrentSystem].TorsionType[CurrentFramework][n]==TorsionType))
                            {
                               present=TRUE;
                               for(j=0;j<TorsionTypes[TorsionType].nr_args;j++)
                                 present=present&&(fabs(Framework[CurrentSystem].TorsionArguments[CurrentFramework][n][j]-arguments[j])<1e-5);

                               if(present) break;
                            }
                        }

                        if(!present)
                        {
                          Framework[CurrentSystem].Torsions[CurrentFramework][index].A=A;
                          Framework[CurrentSystem].Torsions[CurrentFramework][index].B=B;
                          Framework[CurrentSystem].Torsions[CurrentFramework][index].C=C;
                          Framework[CurrentSystem].Torsions[CurrentFramework][index].D=D;

                          Framework[CurrentSystem].TorsionType[CurrentFramework][index]=TorsionType;
                          Framework[CurrentSystem].NumberOfTorsionsPerType[i]++;

                          // copy parsed and converted arguments to the torsion-potential parameters
                          for(j=0;j<TorsionTypes[TorsionType].nr_args;j++)
                            Framework[CurrentSystem].TorsionArguments[CurrentFramework][index][j]=arguments[j];

                          Framework[CurrentSystem].TorsionArguments[CurrentFramework][index][6]=Framework[CurrentSystem].TorsionArgumentDefinitions[i][6];
                          Framework[CurrentSystem].TorsionArguments[CurrentFramework][index][7]=Framework[CurrentSystem].TorsionArgumentDefinitions[i][7];

                          Framework[CurrentSystem].NumberOfTorsions[CurrentFramework]++;
                          index=Framework[CurrentSystem].NumberOfTorsions[CurrentFramework];

                          if(index>=Framework[CurrentSystem].MaxNumberOfTorsions[CurrentFramework])
                          {
                            Framework[CurrentSystem].MaxNumberOfTorsions[CurrentFramework]+=4096;
                            Framework[CurrentSystem].Torsions[CurrentFramework]=(QUAD*)realloc(Framework[CurrentSystem].Torsions[CurrentFramework],
                                        Framework[CurrentSystem].MaxNumberOfTorsions[CurrentFramework]*sizeof(QUAD));
                            Framework[CurrentSystem].TorsionType[CurrentFramework]=(int*)realloc(Framework[CurrentSystem].TorsionType[CurrentFramework],
                                        Framework[CurrentSystem].MaxNumberOfTorsions[CurrentFramework]*sizeof(int));
                            Framework[CurrentSystem].TorsionArguments[CurrentFramework]=(REAL(*)[MAX_TORSION_POTENTIAL_ARGUMENTS])
                                 realloc(Framework[CurrentSystem].TorsionArguments[CurrentFramework],
                                 Framework[CurrentSystem].MaxNumberOfTorsions[CurrentFramework]*sizeof(REAL[MAX_TORSION_POTENTIAL_ARGUMENTS]));
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  // reading Improper torsion-data
  if(Framework[CurrentSystem].NumberOfImproperTorsionDefinitions>0)
  {
    // alloc memory
    Framework[CurrentSystem].ImproperTorsionDefinitionType=(int*)calloc(Framework[CurrentSystem].NumberOfImproperTorsionDefinitions,sizeof(int));
    Framework[CurrentSystem].ImproperTorsionDefinitions=(QUAD*)calloc(Framework[CurrentSystem].NumberOfImproperTorsionDefinitions,sizeof(QUAD));
    Framework[CurrentSystem].ImproperTorsionArgumentDefinitions=(REAL(*)[MAX_IMPROPER_TORSION_POTENTIAL_ARGUMENTS])
      calloc(Framework[CurrentSystem].NumberOfImproperTorsionDefinitions,sizeof(REAL[MAX_IMPROPER_TORSION_POTENTIAL_ARGUMENTS]));
    Framework[CurrentSystem].NumberOfImproperTorsionsPerType=(int*)calloc(Framework[CurrentSystem].NumberOfImproperTorsionDefinitions,sizeof(int));

    for(CurrentFramework=0;CurrentFramework<Framework[CurrentSystem].NumberOfFrameworks;CurrentFramework++)
    {
      Framework[CurrentSystem].MaxNumberOfImproperTorsions[CurrentFramework]=4096;
      // allocate enough memory, it will be resized when needed
      Framework[CurrentSystem].ImproperTorsions[CurrentFramework]=(QUAD*)calloc(Framework[CurrentSystem].MaxNumberOfImproperTorsions[CurrentFramework],sizeof(QUAD));
      Framework[CurrentSystem].ImproperTorsionType[CurrentFramework]=(int*)calloc(Framework[CurrentSystem].MaxNumberOfImproperTorsions[CurrentFramework],sizeof(int));
      Framework[CurrentSystem].ImproperTorsionArguments[CurrentFramework]=(REAL(*)[MAX_IMPROPER_TORSION_POTENTIAL_ARGUMENTS])
              calloc(Framework[CurrentSystem].MaxNumberOfImproperTorsions[CurrentFramework],sizeof(REAL[MAX_IMPROPER_TORSION_POTENTIAL_ARGUMENTS]));
    }

    ReadLine(line,1024,FilePtr); // skip line
    for(i=0;i<Framework[CurrentSystem].NumberOfImproperTorsionDefinitions;i++)
    {
      ReadLine(line,1024,FilePtr);
      arg_pointer=line;
      sscanf(line,"%s%s%s%s%s%n",TypeNameA,TypeNameB,TypeNameC,TypeNameD,buffer,&n);
      TypeA=ReturnPseudoAtomNumber(TypeNameA);
      TypeB=ReturnPseudoAtomNumber(TypeNameB);
      TypeC=ReturnPseudoAtomNumber(TypeNameC);
      TypeD=ReturnPseudoAtomNumber(TypeNameD);

      // determine inversionbend-type
      for(j=0;j<NR_IMPROPER_TORSION_TYPES;j++)
        if(strncasecmp(ImproperTorsionTypes[j].Name,buffer,MAX2(strlen(ImproperTorsionTypes[j].Name),strlen(buffer)))==0)
          ImproperTorsionType=j;

      // read arguments
      for(j=0;j<ImproperTorsionTypes[ImproperTorsionType].nr_args;j++)
      {
        arg_pointer+=n;
        sscanf(arg_pointer,"%lf%n",&temp,&n);
        arguments[j]=(REAL)temp;
      }

      Framework[CurrentSystem].ImproperTorsionDefinitionType[i]=ImproperTorsionType;
      Framework[CurrentSystem].ImproperTorsionDefinitions[i].A=TypeA;
      Framework[CurrentSystem].ImproperTorsionDefinitions[i].B=TypeB;
      Framework[CurrentSystem].ImproperTorsionDefinitions[i].C=TypeC;
      Framework[CurrentSystem].ImproperTorsionDefinitions[i].D=TypeD;

      for(j=0;j<ImproperTorsionTypes[ImproperTorsionType].nr_args;j++)
        Framework[CurrentSystem].ImproperTorsionArgumentDefinitions[i][j]=arguments[j];

      // fill the appropriate improper torsion-data
      for(CurrentFramework=0;CurrentFramework<Framework[CurrentSystem].NumberOfFrameworks;CurrentFramework++)
      {
        index=Framework[CurrentSystem].NumberOfImproperTorsions[CurrentFramework];
        for(A=0;A<Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];A++)
        {
          CurrentTypeA=Framework[CurrentSystem].Atoms[CurrentFramework][A].Type;

          if(CurrentTypeA==TypeA)
          {
            for(k=0;k<Framework[CurrentSystem].Connectivity[CurrentFramework][A];k++)
            {
              B=GetNeighbour(CurrentSystem,CurrentFramework,A,k);
              CurrentTypeB=Framework[CurrentSystem].Atoms[CurrentFramework][B].Type;

              if(CurrentTypeB==TypeB)
              {
                for(l=0;l<Framework[CurrentSystem].Connectivity[CurrentFramework][B];l++)
                {
                  C=GetNeighbour(CurrentSystem,CurrentFramework,B,l);
                  CurrentTypeC=Framework[CurrentSystem].Atoms[CurrentFramework][C].Type;
                  if((A!=C)&&(CurrentTypeC==TypeC))
                  {
                    for(m=0;m<Framework[CurrentSystem].Connectivity[CurrentFramework][B];m++)
                    {
                      D=GetNeighbour(CurrentSystem,CurrentFramework,B,m);
                      CurrentTypeD=Framework[CurrentSystem].Atoms[CurrentFramework][D].Type;
                      if((D!=A)&&(D!=C)&&(CurrentTypeD==TypeD))
                      {
                        present=FALSE;

                        switch(ImproperTorsionScanType)
                        {
                          case IMPROPER_TORSION_SCAN_GENERAL:
                            for(n=0;n<index;n++)
                              if(((Framework[CurrentSystem].ImproperTorsions[CurrentFramework][n].A==A)&&(Framework[CurrentSystem].ImproperTorsions[CurrentFramework][n].B==B)&&
                                 (Framework[CurrentSystem].ImproperTorsions[CurrentFramework][n].C==C)&&(Framework[CurrentSystem].ImproperTorsions[CurrentFramework][n].D==D))||
                                 ((Framework[CurrentSystem].ImproperTorsions[CurrentFramework][n].A==A)&&(Framework[CurrentSystem].ImproperTorsions[CurrentFramework][n].B==B)&&
                                 (Framework[CurrentSystem].ImproperTorsions[CurrentFramework][n].C==D)&&(Framework[CurrentSystem].ImproperTorsions[CurrentFramework][n].D==C))||
                                 ((Framework[CurrentSystem].ImproperTorsions[CurrentFramework][n].A==C)&&(Framework[CurrentSystem].ImproperTorsions[CurrentFramework][n].B==B)&&
                                 (Framework[CurrentSystem].ImproperTorsions[CurrentFramework][n].C==A)&&(Framework[CurrentSystem].ImproperTorsions[CurrentFramework][n].D==D))||
                                 ((Framework[CurrentSystem].ImproperTorsions[CurrentFramework][n].A==C)&&(Framework[CurrentSystem].ImproperTorsions[CurrentFramework][n].B==B)&&
                                 (Framework[CurrentSystem].ImproperTorsions[CurrentFramework][n].C==D)&&(Framework[CurrentSystem].ImproperTorsions[CurrentFramework][n].D==A))||
                                 ((Framework[CurrentSystem].ImproperTorsions[CurrentFramework][n].A==D)&&(Framework[CurrentSystem].ImproperTorsions[CurrentFramework][n].B==B)&&
                                 (Framework[CurrentSystem].ImproperTorsions[CurrentFramework][n].C==A)&&(Framework[CurrentSystem].ImproperTorsions[CurrentFramework][n].D==C))||
                                 ((Framework[CurrentSystem].ImproperTorsions[CurrentFramework][n].A==D)&&(Framework[CurrentSystem].ImproperTorsions[CurrentFramework][n].B==B)&&
                                 (Framework[CurrentSystem].ImproperTorsions[CurrentFramework][n].C==C)&&(Framework[CurrentSystem].ImproperTorsions[CurrentFramework][n].D==A)))
                                   present=TRUE;
                            break;
                          case IMPROPER_TORSION_SCAN_UNIQUE:
                            if(((Framework[CurrentSystem].ImproperTorsions[CurrentFramework][n].A==A)&&(Framework[CurrentSystem].ImproperTorsions[CurrentFramework][n].B==B)&&
                               (Framework[CurrentSystem].ImproperTorsions[CurrentFramework][n].C==C)&&(Framework[CurrentSystem].ImproperTorsions[CurrentFramework][n].D==D))&&
                               (Framework[CurrentSystem].ImproperTorsionType[CurrentFramework][n]==ImproperTorsionType))
                            {
                               present=TRUE;
                               for(j=0;j<ImproperTorsionTypes[ImproperTorsionType].nr_args;j++)
                                 present=present&&(fabs(Framework[CurrentSystem].ImproperTorsionArguments[CurrentFramework][n][j]-arguments[j])<1e-5);
                            }
                            break;
                        }

                        if(!present)
                        {
                          Framework[CurrentSystem].ImproperTorsions[CurrentFramework][index].A=A;
                          Framework[CurrentSystem].ImproperTorsions[CurrentFramework][index].B=B;
                          Framework[CurrentSystem].ImproperTorsions[CurrentFramework][index].C=C;
                          Framework[CurrentSystem].ImproperTorsions[CurrentFramework][index].D=D;

                          Framework[CurrentSystem].ImproperTorsionType[CurrentFramework][index]=ImproperTorsionType;
                          Framework[CurrentSystem].NumberOfImproperTorsionsPerType[i]++;

                          for(j=0;j<ImproperTorsionTypes[ImproperTorsionType].nr_args;j++)
                            Framework[CurrentSystem].ImproperTorsionArguments[CurrentFramework][index][j]=arguments[j];
                          switch(ImproperTorsionType)
                          {
                            case HARMONIC_IMPROPER_DIHEDRAL:
                              // (1/2)*p_0*(phi-p_1)^2
                              // ===============================================
                              // p_0/k_B [K/rad^2]
                              // p_1     [degrees]
                              Framework[CurrentSystem].ImproperTorsionArguments[CurrentFramework][index][0]*=KELVIN_TO_ENERGY;
                              break;
                            case HARMONIC_COSINE_IMPROPER_DIHEDRAL:
                              // (1/2)*p_0*(cos(phi)-cos(p_1))^2
                              // ===============================================
                              // p_0/k_B [K]
                              // p_1     [degrees]
                              Framework[CurrentSystem].ImproperTorsionArguments[CurrentFramework][index][0]*=KELVIN_TO_ENERGY;
                              Framework[CurrentSystem].ImproperTorsionArguments[CurrentFramework][index][1]=acos(RAD2DEG*Framework[CurrentSystem].ImproperTorsionArguments[CurrentFramework][index][1]);
                              break;
                            case THREE_COSINE_IMPROPER_DIHEDRAL:
                              // (1/2)*p_0*(1+cos(phi))+(1/2)*p_1*(1-cos(2*phi))+(1/2)*p_2*(1+cos(3*phi))
                              // ========================================================================
                              // p_0/k_B [K]
                              // p_1/k_B [K]
                              // p_2/k_B [K]
                              Framework[CurrentSystem].ImproperTorsionArguments[CurrentFramework][index][0]*=KELVIN_TO_ENERGY;
                              Framework[CurrentSystem].ImproperTorsionArguments[CurrentFramework][index][1]*=KELVIN_TO_ENERGY;
                              Framework[CurrentSystem].ImproperTorsionArguments[CurrentFramework][index][2]*=KELVIN_TO_ENERGY;
                              break;
                            case MM3_IMPROPER_DIHEDRAL:
                              // (1/2)*p_0*(1+cos(phi))+(1/2)*p_1*(1-cos(2*phi))+(1/2)*p_2*(1+cos(3*phi))
                              // ========================================================================
                              // p_0     [kcal/mol]
                              // p_1     [kcal/mol]
                              // p_2     [kcal/mol]
                              Framework[CurrentSystem].ImproperTorsionArguments[CurrentFramework][index][0]*=KCAL_PER_MOL_TO_ENERGY;
                              Framework[CurrentSystem].ImproperTorsionArguments[CurrentFramework][index][1]*=KCAL_PER_MOL_TO_ENERGY;
                              Framework[CurrentSystem].ImproperTorsionArguments[CurrentFramework][index][2]*=KCAL_PER_MOL_TO_ENERGY;
                              break;
                            case CFF_IMPROPER_DIHEDRAL:
                              // p_0*(1-cos(phi))+p_1*(1-cos(2*phi))+p_2*(1-cos(3*phi))
                              // ======================================================
                              // p_0/k_B [K]
                              // p_1/k_B [K]
                              // p_2/k_B [K]
                              Framework[CurrentSystem].ImproperTorsionArguments[CurrentFramework][index][0]*=KELVIN_TO_ENERGY;
                              Framework[CurrentSystem].ImproperTorsionArguments[CurrentFramework][index][1]*=KELVIN_TO_ENERGY;
                              Framework[CurrentSystem].ImproperTorsionArguments[CurrentFramework][index][2]*=KELVIN_TO_ENERGY;
                              break;
                            case CFF_IMPROPER_DIHEDRAL2:
                              // p_0*(1+cos(phi))+p_1*(1+cos(2*phi))+p_2*(1+cos(3*phi))
                              // ======================================================
                              // p_0/k_B [K]
                              // p_1/k_B [K]
                              // p_2/k_B [K]
                              Framework[CurrentSystem].ImproperTorsionArguments[CurrentFramework][index][0]*=KELVIN_TO_ENERGY;
                              Framework[CurrentSystem].ImproperTorsionArguments[CurrentFramework][index][1]*=KELVIN_TO_ENERGY;
                              Framework[CurrentSystem].ImproperTorsionArguments[CurrentFramework][index][2]*=KELVIN_TO_ENERGY;
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
                              Framework[CurrentSystem].ImproperTorsionArguments[CurrentFramework][index][0]*=KELVIN_TO_ENERGY;
                              Framework[CurrentSystem].ImproperTorsionArguments[CurrentFramework][index][1]*=KELVIN_TO_ENERGY;
                              Framework[CurrentSystem].ImproperTorsionArguments[CurrentFramework][index][2]*=KELVIN_TO_ENERGY;
                              Framework[CurrentSystem].ImproperTorsionArguments[CurrentFramework][index][3]*=KELVIN_TO_ENERGY;
                              Framework[CurrentSystem].ImproperTorsionArguments[CurrentFramework][index][4]*=KELVIN_TO_ENERGY;
                              Framework[CurrentSystem].ImproperTorsionArguments[CurrentFramework][index][5]*=KELVIN_TO_ENERGY;
                              break;
                            case TRAPPE_IMPROPER_DIHEDRAL:
                              // p_0[0]+p_1*(1+cos(phi))+p_2*(1-cos(2*phi))+p_3*(1+cos(3*phi))
                              // =============================================================
                              // p_0/k_B [K]
                              // p_1/k_B [K]
                              // p_2/k_B [K]
                              // p_3/k_B [K]
                              Framework[CurrentSystem].ImproperTorsionArguments[CurrentFramework][index][0]*=KELVIN_TO_ENERGY;
                              Framework[CurrentSystem].ImproperTorsionArguments[CurrentFramework][index][1]*=KELVIN_TO_ENERGY;
                              Framework[CurrentSystem].ImproperTorsionArguments[CurrentFramework][index][2]*=KELVIN_TO_ENERGY;
                              Framework[CurrentSystem].ImproperTorsionArguments[CurrentFramework][index][3]*=KELVIN_TO_ENERGY;
                              break;
                            case CVFF_IMPROPER_DIHEDRAL:
                              // p_0*(1+cos(p_1*phi-p_2))
                              // ========================
                              // p_0/k_B [K]
                              // p_1     [-]
                              // p_2     [degrees]
                              Framework[CurrentSystem].ImproperTorsionArguments[CurrentFramework][index][0]*=KELVIN_TO_ENERGY;
                              Framework[CurrentSystem].ImproperTorsionArguments[CurrentFramework][index][2]*=DEG2RAD;
                              break;
                            case OPLS_IMPROPER_DIHEDRAL:
                              // (1/2)p_0[0]+(1/2)p_1*(1+cos(phi))+(1/2)p_2*(1-cos(2*phi))+(1/2)p_3*(1+cos(3*phi))
                              // =================================================================================
                              // p_0/k_B [K]
                              // p_1/k_B [K]
                              // p_2/k_B [K]
                              // p_3/k_B [K]
                              Framework[CurrentSystem].ImproperTorsionArguments[CurrentFramework][index][0]*=KELVIN_TO_ENERGY;
                              Framework[CurrentSystem].ImproperTorsionArguments[CurrentFramework][index][1]*=KELVIN_TO_ENERGY;
                              Framework[CurrentSystem].ImproperTorsionArguments[CurrentFramework][index][2]*=KELVIN_TO_ENERGY;
                              Framework[CurrentSystem].ImproperTorsionArguments[CurrentFramework][index][3]*=KELVIN_TO_ENERGY;
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
                              Framework[CurrentSystem].ImproperTorsionArguments[CurrentFramework][index][0]*=KELVIN_TO_ENERGY;
                              Framework[CurrentSystem].ImproperTorsionArguments[CurrentFramework][index][1]*=KELVIN_TO_ENERGY;
                              Framework[CurrentSystem].ImproperTorsionArguments[CurrentFramework][index][2]*=KELVIN_TO_ENERGY;
                              Framework[CurrentSystem].ImproperTorsionArguments[CurrentFramework][index][3]*=KELVIN_TO_ENERGY;
                              Framework[CurrentSystem].ImproperTorsionArguments[CurrentFramework][index][4]*=KELVIN_TO_ENERGY;
                              Framework[CurrentSystem].ImproperTorsionArguments[CurrentFramework][index][5]*=KELVIN_TO_ENERGY;
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
                              Framework[CurrentSystem].ImproperTorsionArguments[CurrentFramework][index][0]*=KELVIN_TO_ENERGY;
                              Framework[CurrentSystem].ImproperTorsionArguments[CurrentFramework][index][1]*=KELVIN_TO_ENERGY;
                              Framework[CurrentSystem].ImproperTorsionArguments[CurrentFramework][index][2]*=KELVIN_TO_ENERGY;
                              Framework[CurrentSystem].ImproperTorsionArguments[CurrentFramework][index][3]*=KELVIN_TO_ENERGY;
                              Framework[CurrentSystem].ImproperTorsionArguments[CurrentFramework][index][4]*=KELVIN_TO_ENERGY;
                              Framework[CurrentSystem].ImproperTorsionArguments[CurrentFramework][index][5]*=KELVIN_TO_ENERGY;
                              break;
                            default:
                              fprintf(stderr, "Undefined Improper Torsion potential in routine 'ReadFrameworkDefinition' ('framework.c')\n");
                              exit(0);
                              break;
                          }
                          Framework[CurrentSystem].NumberOfImproperTorsions[CurrentFramework]++;
                          index=Framework[CurrentSystem].NumberOfImproperTorsions[CurrentFramework];

                          if(index>=Framework[CurrentSystem].MaxNumberOfImproperTorsions[CurrentFramework])
                          {
                            Framework[CurrentSystem].MaxNumberOfImproperTorsions[CurrentFramework]+=4096;
                            Framework[CurrentSystem].ImproperTorsions[CurrentFramework]=(QUAD*)realloc(Framework[CurrentSystem].ImproperTorsions[CurrentFramework],
                                        Framework[CurrentSystem].MaxNumberOfImproperTorsions[CurrentFramework]*sizeof(QUAD));
                            Framework[CurrentSystem].ImproperTorsionType[CurrentFramework]=(int*)realloc(Framework[CurrentSystem].ImproperTorsionType[CurrentFramework],
                                        Framework[CurrentSystem].MaxNumberOfImproperTorsions[CurrentFramework]*sizeof(int));
                            Framework[CurrentSystem].ImproperTorsionArguments[CurrentFramework]=(REAL(*)[MAX_IMPROPER_TORSION_POTENTIAL_ARGUMENTS])
                                 realloc(Framework[CurrentSystem].ImproperTorsionArguments[CurrentFramework],
                                 Framework[CurrentSystem].MaxNumberOfImproperTorsions[CurrentFramework]*sizeof(REAL[MAX_IMPROPER_TORSION_POTENTIAL_ARGUMENTS]));
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  // READ OUT_OF_PLANE TODO


  // reading Bond-Bond cross-term data
  index=0;
  if(Framework[CurrentSystem].NumberOfBondBondDefinitions>0)
  {
    // alloc memory
    Framework[CurrentSystem].BondBondDefinitionType=(int*)calloc(Framework[CurrentSystem].NumberOfBondBondDefinitions,sizeof(int));
    Framework[CurrentSystem].BondBondDefinitions=(TRIPLE*)calloc(Framework[CurrentSystem].NumberOfBondBondDefinitions,sizeof(TRIPLE));
    Framework[CurrentSystem].BondBondArgumentDefinitions=(REAL(*)[MAX_BOND_BOND_POTENTIAL_ARGUMENTS])
      calloc(Framework[CurrentSystem].NumberOfBondBondDefinitions,sizeof(REAL[MAX_BOND_BOND_POTENTIAL_ARGUMENTS]));
    Framework[CurrentSystem].NumberOfBondBondsPerType=(int*)calloc(Framework[CurrentSystem].NumberOfBondBondDefinitions,sizeof(int));

    for(CurrentFramework=0;CurrentFramework<Framework[CurrentSystem].NumberOfFrameworks;CurrentFramework++)
    {
      Framework[CurrentSystem].MaxNumberOfBondBonds[CurrentFramework]=4096;
      // allocate enough memory, it will be resized when needed
      Framework[CurrentSystem].BondBonds[CurrentFramework]=(TRIPLE*)calloc(Framework[CurrentSystem].MaxNumberOfBondBonds[CurrentFramework],sizeof(TRIPLE));
      Framework[CurrentSystem].BondBondType[CurrentFramework]=(int*)calloc(Framework[CurrentSystem].MaxNumberOfBondBonds[CurrentFramework],sizeof(int));
      Framework[CurrentSystem].BondBondArguments[CurrentFramework]=(REAL(*)[MAX_BOND_BOND_POTENTIAL_ARGUMENTS])
              calloc(Framework[CurrentSystem].MaxNumberOfBondBonds[CurrentFramework],sizeof(REAL[MAX_BOND_BOND_POTENTIAL_ARGUMENTS]));
    }

    ReadLine(line,1024,FilePtr); // skip line
    for(i=0;i<Framework[CurrentSystem].NumberOfBondBondDefinitions;i++)
    {
      ReadLine(line,1024,FilePtr);
      arg_pointer=line;
      sscanf(line,"%s%s%s%s%n",TypeNameA,TypeNameB,TypeNameC,buffer,&n);
      TypeA=ReturnPseudoAtomNumber(TypeNameA);
      TypeB=ReturnPseudoAtomNumber(TypeNameB);
      TypeC=ReturnPseudoAtomNumber(TypeNameC);

      // determine bond/bond-type
      for(j=0;j<NR_BOND_BOND_TYPES;j++)
        if(strncasecmp(BondBondTypes[j].Name,buffer,MAX2(strlen(BondBondTypes[j].Name),strlen(buffer)))==0)
          BondBondType=j;

      // read arguments
      for(j=0;j<BondBondTypes[BondBondType].nr_args;j++)
      {
        arg_pointer+=n;
        sscanf(arg_pointer,"%lf%n",&temp,&n);
        arguments[j]=(REAL)temp;
      }

      Framework[CurrentSystem].BondBondDefinitionType[i]=BondBondType;
      Framework[CurrentSystem].BondBondDefinitions[i].A=TypeA;
      Framework[CurrentSystem].BondBondDefinitions[i].B=TypeB;
      Framework[CurrentSystem].BondBondDefinitions[i].C=TypeC;
      for(j=0;j<BondBondTypes[BondBondType].nr_args;j++)
        Framework[CurrentSystem].BondBondArgumentDefinitions[i][j]=arguments[j];

      // fill the appropriate bond-bond data
      for(CurrentFramework=0;CurrentFramework<Framework[CurrentSystem].NumberOfFrameworks;CurrentFramework++)
      {
        index=Framework[CurrentSystem].NumberOfBondBonds[CurrentSystem];
        for(B=0;B<Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];B++)
        {
          CurrentTypeB=Framework[CurrentSystem].Atoms[CurrentFramework][B].Type;
          if(CurrentTypeB==TypeB)
          {
            for(k=0;k<Framework[CurrentSystem].Connectivity[CurrentFramework][B];k++)
            {
              A=GetNeighbour(CurrentSystem,CurrentFramework,B,k);
              CurrentTypeA=Framework[CurrentSystem].Atoms[CurrentFramework][A].Type;

              if(CurrentTypeA==TypeA)
              {
                for(l=0;l<Framework[CurrentSystem].Connectivity[CurrentFramework][B];l++)
                {
                  C=GetNeighbour(CurrentSystem,CurrentFramework,B,l);
                  CurrentTypeC=Framework[CurrentSystem].Atoms[CurrentFramework][C].Type;

                  if((A!=C)&&(CurrentTypeC==TypeC))
                  {
                    // check whether triplet A-B-C is already present in the list
                    present=FALSE;
                    for(n=0;n<index;n++)
                      if(((Framework[CurrentSystem].BondBonds[CurrentFramework][n].A==A)&&(Framework[CurrentSystem].BondBonds[CurrentFramework][n].B==B)&&
                        (Framework[CurrentSystem].BondBonds[CurrentFramework][n].C==C))||((Framework[CurrentSystem].BondBonds[CurrentFramework][n].A==C)&&
                        (Framework[CurrentSystem].BondBonds[CurrentFramework][n].B==B)&&(Framework[CurrentSystem].BondBonds[CurrentFramework][n].C==A)))
                          present=TRUE;
                    if(!present)
                    {
                      Framework[CurrentSystem].BondBonds[CurrentFramework][index].A=A;
                      Framework[CurrentSystem].BondBonds[CurrentFramework][index].B=B;
                      Framework[CurrentSystem].BondBonds[CurrentFramework][index].C=C;

                      Framework[CurrentSystem].BondBondType[CurrentFramework][index]=BondBondType;
                      Framework[CurrentSystem].NumberOfBondBondsPerType[i]++;

                      for(j=0;j<BondBondTypes[BondBondType].nr_args;j++)
                        Framework[CurrentSystem].BondBondArguments[CurrentFramework][index][j]=arguments[j];

                      // set to appropriate parameters
                      // the order of the A-B-C triplet does matter for bond-bond cross terms
                      // if C-B-A was found reverse the appropriate arguments
                    switch(BondBondType)
                      {
                        case CVFF_BOND_BOND_CROSS:
                        case CFF_BOND_BOND_CROSS:
                          // p_0*(rab-p_1)*(rbc-p_2)
                          // =======================
                          // p_0/k_B [K/A^2]
                          // p_1     [A]
                          // p_2     [A]
                          Framework[CurrentSystem].BondBondArguments[CurrentFramework][index][0]*=KELVIN_TO_ENERGY;
                          break;
                        default:
                          fprintf(stderr, "Undefined Bond-Bond potential in routine 'ReadFrameworkDefinition' ('framework.c')\n");
                          exit(0);
                          break;
                      }
                      Framework[CurrentSystem].NumberOfBondBonds[CurrentSystem]++;
                      index=Framework[CurrentSystem].NumberOfBondBonds[CurrentSystem];

                      if(index>=Framework[CurrentSystem].MaxNumberOfBondBonds[CurrentFramework])
                      {
                        Framework[CurrentSystem].MaxNumberOfBondBonds[CurrentFramework]+=4096;
                        Framework[CurrentSystem].BondBonds[CurrentFramework]=(TRIPLE*)realloc(Framework[CurrentSystem].BondBonds[CurrentFramework],
                                    Framework[CurrentSystem].MaxNumberOfBondBonds[CurrentFramework]*sizeof(TRIPLE));
                        Framework[CurrentSystem].BondBondType[CurrentFramework]=(int*)realloc(Framework[CurrentSystem].BondBondType[CurrentFramework],
                                    Framework[CurrentSystem].MaxNumberOfBondBonds[CurrentFramework]*sizeof(int));
                        Framework[CurrentSystem].BondBondArguments[CurrentFramework]=(REAL(*)[MAX_BOND_BOND_POTENTIAL_ARGUMENTS])
                              realloc(Framework[CurrentSystem].BondBondArguments[CurrentFramework],
                        Framework[CurrentSystem].MaxNumberOfBondBonds[CurrentFramework]*sizeof(REAL[MAX_BOND_BOND_POTENTIAL_ARGUMENTS]));
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  // reading Bond-Bend cross-term data
  index=0;
  if(Framework[CurrentSystem].NumberOfBondBendDefinitions>0)
  {
    // alloc memory
    Framework[CurrentSystem].BondBendDefinitionType=(int*)calloc(Framework[CurrentSystem].NumberOfBondBendDefinitions,sizeof(int));
    Framework[CurrentSystem].BondBendDefinitions=(TRIPLE*)calloc(Framework[CurrentSystem].NumberOfBondBendDefinitions,sizeof(TRIPLE));
    Framework[CurrentSystem].BondBendArgumentDefinitions=(REAL(*)[MAX_BOND_BEND_POTENTIAL_ARGUMENTS])
      calloc(Framework[CurrentSystem].NumberOfBondBendDefinitions,sizeof(REAL[MAX_BOND_BEND_POTENTIAL_ARGUMENTS]));
    Framework[CurrentSystem].NumberOfBondBendsPerType=(int*)calloc(Framework[CurrentSystem].NumberOfBondBendDefinitions,sizeof(int));

    for(CurrentFramework=0;CurrentFramework<Framework[CurrentSystem].NumberOfFrameworks;CurrentFramework++)
    {
      Framework[CurrentSystem].MaxNumberOfBondBends[CurrentFramework]=4096;
      // allocate enough memory, it will be resized when needed
      Framework[CurrentSystem].BondBends[CurrentFramework]=(TRIPLE*)calloc(Framework[CurrentSystem].MaxNumberOfBondBends[CurrentFramework],sizeof(TRIPLE));
      Framework[CurrentSystem].BondBendType[CurrentFramework]=(int*)calloc(Framework[CurrentSystem].MaxNumberOfBondBends[CurrentFramework],sizeof(int));
      Framework[CurrentSystem].BondBendArguments[CurrentFramework]=(REAL(*)[MAX_BOND_BEND_POTENTIAL_ARGUMENTS])
              calloc(Framework[CurrentSystem].MaxNumberOfBondBends[CurrentFramework],sizeof(REAL[MAX_BOND_BEND_POTENTIAL_ARGUMENTS]));
    }

    ReadLine(line,1024,FilePtr); // skip line
    for(i=0;i<Framework[CurrentSystem].NumberOfBondBendDefinitions;i++)
    {
      ReadLine(line,1024,FilePtr);
      arg_pointer=line;
      sscanf(line,"%s%s%s%s%n",TypeNameA,TypeNameB,TypeNameC,buffer,&n);
      TypeA=ReturnPseudoAtomNumber(TypeNameA);
      TypeB=ReturnPseudoAtomNumber(TypeNameB);
      TypeC=ReturnPseudoAtomNumber(TypeNameC);

      // determine bond/bond-type
      for(j=0;j<NR_BOND_BEND_TYPES;j++)
        if(strcasecmp(BondBendTypes[j].Name,buffer)==0)
          BondBendType=j;

      // read arguments
      for(j=0;j<BondBendTypes[BondBendType].nr_args;j++)
      {
        arg_pointer+=n;
        sscanf(arg_pointer,"%lf%n",&temp,&n);
        arguments[j]=(REAL)temp;
      }

      Framework[CurrentSystem].BondBendDefinitionType[i]=BondBendType;
      Framework[CurrentSystem].BondBendDefinitions[i].A=TypeA;
      Framework[CurrentSystem].BondBendDefinitions[i].B=TypeB;
      Framework[CurrentSystem].BondBendDefinitions[i].C=TypeC;
      for(j=0;j<BondBendTypes[BondBendType].nr_args;j++)
        Framework[CurrentSystem].BondBendArgumentDefinitions[i][j]=arguments[j];

      // fill the appropriate stretch bend-data
      for(CurrentFramework=0;CurrentFramework<Framework[CurrentSystem].NumberOfFrameworks;CurrentFramework++)
      {
        index=Framework[CurrentSystem].NumberOfBondBends[CurrentFramework];
        for(A=0;A<Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];A++)
        {
          CurrentTypeA=Framework[CurrentSystem].Atoms[CurrentFramework][A].Type;
          if(CurrentTypeA==TypeA)
          {
            for(k=0;k<Framework[CurrentSystem].Connectivity[CurrentFramework][A];k++)
            {
              B=GetNeighbour(CurrentSystem,CurrentFramework,A,k);
              CurrentTypeB=Framework[CurrentSystem].Atoms[CurrentFramework][B].Type;

              if(CurrentTypeB==TypeB)
              for(l=0;l<Framework[CurrentSystem].Connectivity[CurrentFramework][B];l++)
              {
                C=GetNeighbour(CurrentSystem,CurrentFramework,B,l);
                CurrentTypeC=Framework[CurrentSystem].Atoms[CurrentFramework][C].Type;

                if((A!=C)&&(CurrentTypeC==TypeC))
                {
                  // check whether triplet A-B-C is already present in the list
                  present=FALSE;
                  for(n=0;n<index;n++)
                  {
                    if(((Framework[CurrentSystem].BondBends[CurrentFramework][n].A==A)&&
                        (Framework[CurrentSystem].BondBends[CurrentFramework][n].B==B)&&
                        (Framework[CurrentSystem].BondBends[CurrentFramework][n].C==C))||
                       ((Framework[CurrentSystem].BondBends[CurrentFramework][n].A==C)&&
                        (Framework[CurrentSystem].BondBends[CurrentFramework][n].B==B)&&
                        (Framework[CurrentSystem].BondBends[CurrentFramework][n].C==A)))
                      present=TRUE;
                  }
                  if(!present)
                  {
                    Framework[CurrentSystem].BondBends[CurrentFramework][index].A=A;
                    Framework[CurrentSystem].BondBends[CurrentFramework][index].B=B;
                    Framework[CurrentSystem].BondBends[CurrentFramework][index].C=C;

                    Framework[CurrentSystem].BondBendType[CurrentFramework][index]=BondBendType;
                    Framework[CurrentSystem].NumberOfBondBendsPerType[i]++;

                    for(j=0;j<BondBendTypes[BondBendType].nr_args;j++)
                      Framework[CurrentSystem].BondBendArguments[CurrentFramework][index][j]=arguments[j];

                    // set to appropriate parameters
                    // the order of the A-B-C triplet does matter for bond-bond cross terms
                    // if C-B-A was found reverse the appropriate arguments
                    switch(BondBendType)
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
                        Framework[CurrentSystem].BondBendArguments[CurrentFramework][index][0]*=DEG2RAD;
                        Framework[CurrentSystem].BondBendArguments[CurrentFramework][index][1]*=KELVIN_TO_ENERGY;
                        Framework[CurrentSystem].BondBendArguments[CurrentFramework][index][3]*=KELVIN_TO_ENERGY;
                        break;
                      case MM3_BOND_BEND_CROSS:
                        // p_0*[(rab-p_1)+(rbc-p_2)]*(Theta-p_3)
                        // =====================================
                        // p_0     [mdyne/rad]
                        // p_1     [A]
                        // p_2     [A]
                        // p_3     [degrees]
                        Framework[CurrentSystem].BondBendArguments[CurrentFramework][index][0]*=2.51118*KCAL_PER_MOL_TO_ENERGY;
                        Framework[CurrentSystem].BondBendArguments[CurrentFramework][index][3]*=DEG2RAD;
                        break;
                      case TRUNCATED_HARMONIC:
                        // (1/2)*p_0*(Theta-p_1)^2*exp(-(pow(rab,8)+pow(rbc,8))/pow(p_2,8))
                        // ================================================================
                        // p_0/k_B [K/rad^2]
                        // p_1     [degrees]
                        // p_2     [A]
                        Framework[CurrentSystem].BondBendArguments[CurrentFramework][index][0]*=KELVIN_TO_ENERGY;
                        Framework[CurrentSystem].BondBendArguments[CurrentFramework][index][1]*=DEG2RAD;
                        break;
                      case SCREENED_HARMONIC:
                        // (1/2)*p_0*(Theta-p_1)^2*exp(-(rab/p_2+rbc/p_3))
                        // ===============================================
                        // p_0/k_B [K/rad^2]
                        // p_1     [degrees]
                        // p_2     [A]
                        // p_3     [A]
                        Framework[CurrentSystem].BondBendArguments[CurrentFramework][index][0]*=KELVIN_TO_ENERGY;
                        Framework[CurrentSystem].BondBendArguments[CurrentFramework][index][1]*=DEG2RAD;
                        break;
                      case SCREENED_VESSAL:
                        // (p_0/(8.0*(Theta-PI)^2))*((p_1-PI)^2-(Theta-PI)^2)^2*exp(-(rab/p_2+rbc/p_3))
                        // ============================================================================
                        // p_0/k_B [K/rad^2]
                        // p_1     [degrees]
                        // p_2     [A]
                        // p_3     [A]
                        Framework[CurrentSystem].BondBendArguments[CurrentFramework][index][0]*=KELVIN_TO_ENERGY;
                        Framework[CurrentSystem].BondBendArguments[CurrentFramework][index][1]*=DEG2RAD;
                        break;
                      case TRUNCATED_VESSAL:
                        // p_0*[pow(theta,p_2)*(theta-p_1)^2*(theta+p_1-2.0*PI)^2-0.5*p_2*pow(PI,p_2-1.0)*(theta-p_1)^2*pow(PI-p_1,3)]
                        //    *exp(-(pow(rab,8)+pow(rbc,8))/pow(p_3,8))
                        // ============================================================================
                        // p_0/k_B [K/rad^(4+p_2)]
                        // p_1     [degrees]
                        // p_2     [-]
                        // p_3     [A]
                        Framework[CurrentSystem].BondBendArguments[CurrentFramework][index][0]*=KELVIN_TO_ENERGY;
                        Framework[CurrentSystem].BondBendArguments[CurrentFramework][index][1]*=DEG2RAD;
                        break;
                      default:
                        fprintf(stderr, "Undefined Bond-Bend potential in routine 'ReadFrameworkDefinition' ('framework.c')\n");
                        exit(0);
                        break;
                    }
                    Framework[CurrentSystem].NumberOfBondBends[CurrentFramework]++;
                    index=Framework[CurrentSystem].NumberOfBondBends[CurrentFramework];

                    if(index>=Framework[CurrentSystem].MaxNumberOfBondBends[CurrentFramework])
                    {
                      Framework[CurrentSystem].MaxNumberOfBondBends[CurrentFramework]+=4096;
                      Framework[CurrentSystem].BondBends[CurrentFramework]=(TRIPLE*)realloc(Framework[CurrentSystem].BondBends[CurrentFramework],
                                  Framework[CurrentSystem].MaxNumberOfBondBends[CurrentFramework]*sizeof(TRIPLE));
                      Framework[CurrentSystem].BondBendType[CurrentFramework]=(int*)realloc(Framework[CurrentSystem].BondBendType[CurrentFramework],
                                  Framework[CurrentSystem].MaxNumberOfBondBends[CurrentFramework]*sizeof(int));
                      Framework[CurrentSystem].BondBendArguments[CurrentFramework]=(REAL(*)[MAX_BOND_BEND_POTENTIAL_ARGUMENTS])
                            realloc(Framework[CurrentSystem].BondBendArguments[CurrentFramework],
                      Framework[CurrentSystem].MaxNumberOfBondBends[CurrentFramework]*sizeof(REAL[MAX_BOND_BEND_POTENTIAL_ARGUMENTS]));
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  // reading Bend-Bend cross-term data
  index=0;
  if(Framework[CurrentSystem].NumberOfBendBendDefinitions>0)
  {
    // alloc memory
    Framework[CurrentSystem].BendBendDefinitionType=(int*)calloc(Framework[CurrentSystem].NumberOfBendBendDefinitions,sizeof(int));
    Framework[CurrentSystem].BendBendDefinitions=(QUAD*)calloc(Framework[CurrentSystem].NumberOfBendBendDefinitions,sizeof(QUAD));
    Framework[CurrentSystem].BendBendArgumentDefinitions=(REAL(*)[MAX_BEND_BEND_POTENTIAL_ARGUMENTS])
      calloc(Framework[CurrentSystem].NumberOfBendBendDefinitions,sizeof(REAL[MAX_BEND_BEND_POTENTIAL_ARGUMENTS]));
    Framework[CurrentSystem].NumberOfBendBendsPerType=(int*)calloc(Framework[CurrentSystem].NumberOfBendBendDefinitions,sizeof(int));

    for(CurrentFramework=0;CurrentFramework<Framework[CurrentSystem].NumberOfFrameworks;CurrentFramework++)
    {
      Framework[CurrentSystem].MaxNumberOfBendBends[CurrentFramework]=4096;
      // allocate enough memory, it will be resized when needed
      Framework[CurrentSystem].BendBends[CurrentFramework]=(QUAD*)calloc(Framework[CurrentSystem].MaxNumberOfBendBends[CurrentFramework],sizeof(QUAD));
      Framework[CurrentSystem].BendBendType[CurrentFramework]=(int*)calloc(Framework[CurrentSystem].MaxNumberOfBendBends[CurrentFramework],sizeof(int));
      Framework[CurrentSystem].BendBendArguments[CurrentFramework]=(REAL(*)[MAX_BEND_BEND_POTENTIAL_ARGUMENTS])
              calloc(Framework[CurrentSystem].MaxNumberOfBendBends[CurrentFramework],sizeof(REAL[MAX_BEND_BEND_POTENTIAL_ARGUMENTS]));
    }

    ReadLine(line,1024,FilePtr); // skip line
    for(i=0;i<Framework[CurrentSystem].NumberOfBendBendDefinitions;i++)
    {
      ReadLine(line,1024,FilePtr);
      arg_pointer=line;
      sscanf(line,"%s%s%s%s%s%n",TypeNameA,TypeNameB,TypeNameC,TypeNameD,buffer,&n);
      TypeA=ReturnPseudoAtomNumber(TypeNameA);
      TypeB=ReturnPseudoAtomNumber(TypeNameB);
      TypeC=ReturnPseudoAtomNumber(TypeNameC);
      TypeD=ReturnPseudoAtomNumber(TypeNameD);

      // determine bend/bend-type
      for(j=0;j<NR_BEND_BEND_TYPES;j++)
        if(strncasecmp(BendBendTypes[j].Name,buffer,MAX2(strlen(BendBendTypes[j].Name),strlen(buffer)))==0)
          BendBendType=j;

      // read arguments
      for(j=0;j<BendBendTypes[BendBendType].nr_args;j++)
      {
        arg_pointer+=n;
        sscanf(arg_pointer,"%lf%n",&temp,&n);
        arguments[j]=(REAL)temp;
      }

      Framework[CurrentSystem].BendBendDefinitionType[i]=BendBendType;
      Framework[CurrentSystem].BendBendDefinitions[i].A=TypeA;
      Framework[CurrentSystem].BendBendDefinitions[i].B=TypeB;
      Framework[CurrentSystem].BendBendDefinitions[i].C=TypeC;
      Framework[CurrentSystem].BendBendDefinitions[i].D=TypeD;
      for(j=0;j<BendBendTypes[BendBendType].nr_args;j++)
        Framework[CurrentSystem].BendBendArgumentDefinitions[i][j]=arguments[j];

      // fill the appropriate bend/bend-data
      for(CurrentFramework=0;CurrentFramework<Framework[CurrentSystem].NumberOfFrameworks;CurrentFramework++)
      {
        index=Framework[CurrentSystem].NumberOfBendBends[CurrentFramework];
        for(A=0;A<Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];A++)
        {
          CurrentTypeA=Framework[CurrentSystem].Atoms[CurrentFramework][A].Type;

          if(CurrentTypeA==TypeA)
          {
            for(k=0;k<Framework[CurrentSystem].Connectivity[CurrentFramework][A];k++)
            {
              B=GetNeighbour(CurrentSystem,CurrentFramework,A,k);
              CurrentTypeB=Framework[CurrentSystem].Atoms[CurrentFramework][B].Type;

              if(CurrentTypeB==TypeB)
              {
                for(l=0;l<Framework[CurrentSystem].Connectivity[CurrentFramework][B];l++)
                {
                  C=GetNeighbour(CurrentSystem,CurrentFramework,B,l);
                  CurrentTypeC=Framework[CurrentSystem].Atoms[CurrentFramework][C].Type;
                  if((A!=C)&&(CurrentTypeC==TypeC))
                  {
                    for(m=0;m<Framework[CurrentSystem].Connectivity[CurrentFramework][B];m++)
                    {
                      D=GetNeighbour(CurrentSystem,CurrentFramework,B,m);
                      CurrentTypeD=Framework[CurrentSystem].Atoms[CurrentFramework][D].Type;
                      if((D!=A)&&(D!=B)&&(D!=C)&&(CurrentTypeD==TypeD))
                      {
                        present=FALSE;
                        for(n=0;n<index;n++)
                          if(((Framework[CurrentSystem].BendBends[CurrentFramework][n].A==A)&&
                              (Framework[CurrentSystem].BendBends[CurrentFramework][n].B==B)&&
                              (Framework[CurrentSystem].BendBends[CurrentFramework][n].C==C)&&
                              (Framework[CurrentSystem].BendBends[CurrentFramework][n].D==D))||
                             ((Framework[CurrentSystem].BendBends[CurrentFramework][n].A==A)&&
                              (Framework[CurrentSystem].BendBends[CurrentFramework][n].B==B)&&
                              (Framework[CurrentSystem].BendBends[CurrentFramework][n].C==D)&&
                              (Framework[CurrentSystem].BendBends[CurrentFramework][n].D==C)))
                             present=TRUE;
                        if(!present)
                        {
                          Framework[CurrentSystem].BendBends[CurrentFramework][index].A=A;
                          Framework[CurrentSystem].BendBends[CurrentFramework][index].B=B;
                          Framework[CurrentSystem].BendBends[CurrentFramework][index].C=C;
                          Framework[CurrentSystem].BendBends[CurrentFramework][index].D=D;

                          Framework[CurrentSystem].BendBendType[CurrentFramework][index]=BendBendType;
                          Framework[CurrentSystem].NumberOfBendBendsPerType[i]++;

                          for(j=0;j<BendBendTypes[BendBendType].nr_args;j++)
                            Framework[CurrentSystem].BendBendArguments[CurrentFramework][index][j]=arguments[j];
                          switch(BendBendType)
                          {
                            case CVFF_BEND_BEND_CROSS:
                            case CFF_BEND_BEND_CROSS:
                              // p_0*(Theta1-p_1)*(Theta2-p_2)
                              // ===================================
                              // p_0/k_B [K/rad^2)]
                              // p_1     [degrees]
                              // p_2     [degrees]
                              Framework[CurrentSystem].BendBendArguments[CurrentFramework][index][0]*=KELVIN_TO_ENERGY;
                              Framework[CurrentSystem].BendBendArguments[CurrentFramework][index][1]*=DEG2RAD;
                              Framework[CurrentSystem].BendBendArguments[CurrentFramework][index][2]*=DEG2RAD;
                              break;
                            case MM3_BEND_BEND_CROSS:
                              // -p_0*(Theta1-p_1)*(Theta2-p_2)
                              // ===================================
                              // p_0     [mdyne A/rad^2]
                              // p_1     [degrees]
                              // p_2     [degrees]
                              Framework[CurrentSystem].BendBendArguments[CurrentFramework][index][0]*=0.02191418*KCAL_PER_MOL_TO_ENERGY;
                              Framework[CurrentSystem].BendBendArguments[CurrentFramework][index][1]*=DEG2RAD;
                              Framework[CurrentSystem].BendBendArguments[CurrentFramework][index][2]*=DEG2RAD;
                              break;
                            default:
                              fprintf(stderr, "Undefined Bend-Bend potential in routine 'ReadFrameworkDefinition' ('framework.c')\n");
                              exit(0);
                              break;
                          }
                          Framework[CurrentSystem].NumberOfBendBends[CurrentFramework]++;
                          index=Framework[CurrentSystem].NumberOfBendBends[CurrentFramework];

                          if(index>=Framework[CurrentSystem].MaxNumberOfBendBends[CurrentFramework])
                          {
                            Framework[CurrentSystem].MaxNumberOfBendBends[CurrentFramework]+=4096;
                            Framework[CurrentSystem].BendBends[CurrentFramework]=(QUAD*)realloc(Framework[CurrentSystem].BendBends[CurrentFramework],
                                      Framework[CurrentSystem].MaxNumberOfBendBends[CurrentFramework]*sizeof(QUAD));
                            Framework[CurrentSystem].BendBendType[CurrentFramework]=(int*)realloc(Framework[CurrentSystem].BendBendType[CurrentFramework],
                                      Framework[CurrentSystem].MaxNumberOfBendBends[CurrentFramework]*sizeof(int));
                            Framework[CurrentSystem].BendBendArguments[CurrentFramework]=(REAL(*)[MAX_BEND_BEND_POTENTIAL_ARGUMENTS])
                                 realloc(Framework[CurrentSystem].BendBendArguments[CurrentFramework],
                            Framework[CurrentSystem].MaxNumberOfBendBends[CurrentFramework]*sizeof(REAL[MAX_BEND_BEND_POTENTIAL_ARGUMENTS]));
                          }

                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  // reading Bond/Torsion-data
  if(Framework[CurrentSystem].NumberOfBondTorsionDefinitions>0)
  {
    // alloc memory
    Framework[CurrentSystem].BondTorsionDefinitionType=(int*)calloc(Framework[CurrentSystem].NumberOfBondTorsionDefinitions,sizeof(int));
    Framework[CurrentSystem].BondTorsionDefinitions=(QUAD*)calloc(Framework[CurrentSystem].NumberOfBondTorsionDefinitions,sizeof(QUAD));
    Framework[CurrentSystem].BondTorsionArgumentDefinitions=(REAL(*)[MAX_BOND_TORSION_POTENTIAL_ARGUMENTS])
      calloc(Framework[CurrentSystem].NumberOfBondTorsionDefinitions,sizeof(REAL[MAX_BOND_TORSION_POTENTIAL_ARGUMENTS]));
    Framework[CurrentSystem].NumberOfBondTorsionsPerType=(int*)calloc(Framework[CurrentSystem].NumberOfBondTorsionDefinitions,sizeof(int));

    for(CurrentFramework=0;CurrentFramework<Framework[CurrentSystem].NumberOfFrameworks;CurrentFramework++)
    {
      Framework[CurrentSystem].MaxNumberOfBondTorsions[CurrentFramework]=4096;
      // allocate enough memory, it will be resized when needed
      Framework[CurrentSystem].BondTorsions[CurrentFramework]=(QUAD*)calloc(Framework[CurrentSystem].MaxNumberOfBondTorsions[CurrentFramework],sizeof(QUAD));
      Framework[CurrentSystem].BondTorsionType[CurrentFramework]=(int*)calloc(Framework[CurrentSystem].MaxNumberOfBondTorsions[CurrentFramework],sizeof(int));
      Framework[CurrentSystem].BondTorsionArguments[CurrentFramework]=(REAL(*)[MAX_BOND_TORSION_POTENTIAL_ARGUMENTS])
              calloc(Framework[CurrentSystem].MaxNumberOfBondTorsions[CurrentFramework],sizeof(REAL[MAX_BOND_TORSION_POTENTIAL_ARGUMENTS]));
    }

    ReadLine(line,1024,FilePtr); // skip line
    for(i=0;i<Framework[CurrentSystem].NumberOfBondTorsionDefinitions;i++)
    {
      ReadLine(line,1024,FilePtr);
      arg_pointer=line;
      sscanf(line,"%s%s%s%s%s%n",TypeNameA,TypeNameB,TypeNameC,TypeNameD,buffer,&n);
      TypeA=ReturnPseudoAtomNumber(TypeNameA);
      TypeB=ReturnPseudoAtomNumber(TypeNameB);
      TypeC=ReturnPseudoAtomNumber(TypeNameC);
      TypeD=ReturnPseudoAtomNumber(TypeNameD);

      // determine torsion-type
      for(j=0;j<NR_BOND_TORSION_TYPES;j++)
        if(strncasecmp(BondTorsionTypes[j].Name,buffer,MAX2(strlen(BondTorsionTypes[j].Name),strlen(buffer)))==0)
          BondTorsionType=j;

      // read arguments
      for(j=0;j<BondTorsionTypes[BondTorsionType].nr_args;j++)
      {
        arg_pointer+=n;
        sscanf(arg_pointer,"%lf%n",&temp,&n);
        arguments[j]=(REAL)temp;
      }

      Framework[CurrentSystem].BondTorsionDefinitionType[i]=BondTorsionType;
      Framework[CurrentSystem].BondTorsionDefinitions[i].A=TypeA;
      Framework[CurrentSystem].BondTorsionDefinitions[i].B=TypeB;
      Framework[CurrentSystem].BondTorsionDefinitions[i].C=TypeC;
      Framework[CurrentSystem].BondTorsionDefinitions[i].D=TypeD;
      for(j=0;j<BondTorsionTypes[BondTorsionType].nr_args;j++)
        Framework[CurrentSystem].BondTorsionArgumentDefinitions[i][j]=arguments[j];

      // fill the appropriate stretch/torsion-data
      for(CurrentFramework=0;CurrentFramework<Framework[CurrentSystem].NumberOfFrameworks;CurrentFramework++)
      {
        index=Framework[CurrentSystem].NumberOfBondTorsions[CurrentFramework];

        for(B=0;B<Framework[CurrentSystem].TotalNumberOfAtoms;B++)
        {
          CurrentTypeB=Framework[CurrentSystem].Atoms[CurrentFramework][B].Type;

          if(CurrentTypeB==TypeB)
          {
            for(k=0;k<Framework[CurrentSystem].Connectivity[CurrentFramework][B];k++)
            {
              C=GetNeighbour(CurrentSystem,CurrentFramework,B,k);
              CurrentTypeC=Framework[CurrentSystem].Atoms[CurrentFramework][C].Type;

              if((CurrentTypeC==TypeC)&&(B<C))
              {
                // we now have pairs B-C with B<C
                for(l=0;l<Framework[CurrentSystem].Connectivity[CurrentFramework][B];l++)
                {
                  A=GetNeighbour(CurrentSystem,CurrentFramework,B,l);
                  CurrentTypeA=Framework[CurrentSystem].Atoms[CurrentFramework][A].Type;
                  if((A!=C)&&(CurrentTypeA==TypeA))
                  {
                    for(m=0;m<Framework[CurrentSystem].Connectivity[CurrentFramework][C];m++)
                    {
                      D=GetNeighbour(CurrentSystem,CurrentFramework,C,m);
                      CurrentTypeD=Framework[CurrentSystem].Atoms[CurrentFramework][D].Type;
                      if((B!=D)&&(CurrentTypeD==TypeD))
                      {
                        Framework[CurrentSystem].BondTorsions[CurrentFramework][index].A=A;
                        Framework[CurrentSystem].BondTorsions[CurrentFramework][index].B=B;
                        Framework[CurrentSystem].BondTorsions[CurrentFramework][index].C=C;
                        Framework[CurrentSystem].BondTorsions[CurrentFramework][index].D=D;

                        Framework[CurrentSystem].BondTorsionType[CurrentFramework][index]=BondTorsionType;
                        Framework[CurrentSystem].NumberOfBondTorsionsPerType[i]++;

                        for(j=0;j<BondTorsionTypes[BondTorsionType].nr_args;j++)
                          Framework[CurrentSystem].BondTorsionArguments[CurrentFramework][index][j]=arguments[j];
                        switch(BondTorsionType)
                        {
                          case MM3_BOND_TORSION_CROSS:
                            // (1/2)p_0(r-p_3)(1+cos(phi))+(1/2)p_1(r-p_3)(1+cos(2phi))+(1/2)p_2(r-p_3)(1+cos(3phi))
                            // =====================================================================================
                            // p_0     [kcal/A mole]
                            // p_1     [kcal/A mole]
                            // p_2     [kcal/A mole]
                            // p_3     [A]
                            Framework[CurrentSystem].BondTorsionArguments[CurrentFramework][index][0]*=KCAL_PER_MOL_TO_ENERGY;
                            Framework[CurrentSystem].BondTorsionArguments[CurrentFramework][index][1]*=KCAL_PER_MOL_TO_ENERGY;
                            Framework[CurrentSystem].BondTorsionArguments[CurrentFramework][index][2]*=KCAL_PER_MOL_TO_ENERGY;
                            break;
                          default:
                            fprintf(stderr, "Undefined Bond-Torsion potential in routine 'ReadFrameworkDefinition' ('framework.c')\n");
                            exit(0);
                            break;
                        }
                        Framework[CurrentSystem].NumberOfBondTorsions[CurrentFramework]++;
                        index=Framework[CurrentSystem].NumberOfBondTorsions[CurrentFramework];

                        if(index>=Framework[CurrentSystem].MaxNumberOfBondTorsions[CurrentFramework])
                        {
                          Framework[CurrentSystem].MaxNumberOfBondTorsions[CurrentFramework]+=4096;
                          Framework[CurrentSystem].BondTorsions[CurrentFramework]=(QUAD*)realloc(Framework[CurrentSystem].BondTorsions[CurrentFramework],
                                    Framework[CurrentSystem].MaxNumberOfBondTorsions[CurrentFramework]*sizeof(QUAD));
                          Framework[CurrentSystem].BondTorsionType[CurrentFramework]=(int*)realloc(Framework[CurrentSystem].BondTorsionType[CurrentFramework],
                                    Framework[CurrentSystem].MaxNumberOfBondTorsions[CurrentFramework]*sizeof(int));
                          Framework[CurrentSystem].BondTorsionArguments[CurrentFramework]=(REAL(*)[MAX_BOND_TORSION_POTENTIAL_ARGUMENTS])
                               realloc(Framework[CurrentSystem].BondTorsionArguments[CurrentFramework],
                          Framework[CurrentSystem].MaxNumberOfBondTorsions[CurrentFramework]*sizeof(REAL[MAX_BOND_TORSION_POTENTIAL_ARGUMENTS]));
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  // reading Bend/Torsion-data
  index=0;
  if(Framework[CurrentSystem].NumberOfBendTorsionDefinitions>0)
  {
    // alloc memory
    Framework[CurrentSystem].BendTorsionDefinitionType=(int*)calloc(Framework[CurrentSystem].NumberOfBendTorsionDefinitions,sizeof(int));
    Framework[CurrentSystem].BendTorsionDefinitions=(QUAD*)calloc(Framework[CurrentSystem].NumberOfBendTorsionDefinitions,sizeof(QUAD));
    Framework[CurrentSystem].BendTorsionArgumentDefinitions=(REAL(*)[MAX_BEND_TORSION_POTENTIAL_ARGUMENTS])
      calloc(Framework[CurrentSystem].NumberOfBendTorsionDefinitions,sizeof(REAL[MAX_BEND_TORSION_POTENTIAL_ARGUMENTS]));
    Framework[CurrentSystem].NumberOfBendTorsionsPerType=(int*)calloc(Framework[CurrentSystem].NumberOfBendTorsionDefinitions,sizeof(int));

    for(CurrentFramework=0;CurrentFramework<Framework[CurrentSystem].NumberOfFrameworks;CurrentFramework++)
    {
      Framework[CurrentSystem].MaxNumberOfBendTorsions[CurrentFramework]=4096;
      // allocate enough memory, it will be resized when needed
      Framework[CurrentSystem].BendTorsions[CurrentFramework]=(QUAD*)calloc(Framework[CurrentSystem].MaxNumberOfBendTorsions[CurrentFramework],sizeof(QUAD));
      Framework[CurrentSystem].BendTorsionType[CurrentFramework]=(int*)calloc(Framework[CurrentSystem].MaxNumberOfBendTorsions[CurrentFramework],sizeof(int));
      Framework[CurrentSystem].BendTorsionArguments[CurrentFramework]=(REAL(*)[MAX_BEND_TORSION_POTENTIAL_ARGUMENTS])
              calloc(Framework[CurrentSystem].MaxNumberOfBendTorsions[CurrentFramework],sizeof(REAL[MAX_BEND_TORSION_POTENTIAL_ARGUMENTS]));
    }

    ReadLine(line,1024,FilePtr); // skip line
    for(i=0;i<Framework[CurrentSystem].NumberOfBendTorsionDefinitions;i++)
    {
      ReadLine(line,1024,FilePtr);
      arg_pointer=line;
      sscanf(line,"%s%s%s%s%s%n",TypeNameA,TypeNameB,TypeNameC,TypeNameD,buffer,&n);
      TypeA=ReturnPseudoAtomNumber(TypeNameA);
      TypeB=ReturnPseudoAtomNumber(TypeNameB);
      TypeC=ReturnPseudoAtomNumber(TypeNameC);
      TypeD=ReturnPseudoAtomNumber(TypeNameD);

      // determine torsion-type
      for(j=0;j<NR_BEND_TORSION_TYPES;j++)
      {
        if(strcasecmp(BendTorsionTypes[j].Name,buffer)==0)
          BendTorsionType=j;
      }

      // read arguments
      for(j=0;j<BendTorsionTypes[BendTorsionType].nr_args;j++)
      {
        arg_pointer+=n;
        sscanf(arg_pointer,"%lf%n",&temp,&n);
        arguments[j]=(REAL)temp;
      }

      Framework[CurrentSystem].BendTorsionDefinitionType[i]=BendTorsionType;
      Framework[CurrentSystem].BendTorsionDefinitions[i].A=TypeA;
      Framework[CurrentSystem].BendTorsionDefinitions[i].B=TypeB;
      Framework[CurrentSystem].BendTorsionDefinitions[i].C=TypeC;
      Framework[CurrentSystem].BendTorsionDefinitions[i].D=TypeD;
      for(j=0;j<BendTorsionTypes[BendTorsionType].nr_args;j++)
        Framework[CurrentSystem].BendTorsionArgumentDefinitions[i][j]=arguments[j];

      // fill the appropriate bend/torsion-data
      for(CurrentFramework=0;CurrentFramework<Framework[CurrentSystem].NumberOfFrameworks;CurrentFramework++)
      {
        index=Framework[CurrentSystem].NumberOfBendTorsions[CurrentFramework];
        for(B=0;B<Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];B++)
        {
          CurrentTypeB=Framework[CurrentSystem].Atoms[CurrentFramework][B].Type;

          if(CurrentTypeB==TypeB)
          {
            for(k=0;k<Framework[CurrentSystem].Connectivity[CurrentFramework][B];k++)
            {
              C=GetNeighbour(CurrentSystem,CurrentFramework,B,k);
              CurrentTypeC=Framework[CurrentSystem].Atoms[CurrentFramework][C].Type;

              if(CurrentTypeC==TypeC)
              {
                // we now have pairs A-B
                for(l=0;l<Framework[CurrentSystem].Connectivity[CurrentFramework][B];l++)
                {
                  A=GetNeighbour(CurrentSystem,CurrentFramework,B,l);
                  CurrentTypeA=Framework[CurrentSystem].Atoms[CurrentFramework][A].Type;
                  if((A!=C)&&(CurrentTypeA==TypeA))
                  {
                    for(m=0;m<Framework[CurrentSystem].Connectivity[CurrentFramework][C];m++)
                    {
                      D=GetNeighbour(CurrentSystem,CurrentFramework,C,m);
                      CurrentTypeD=Framework[CurrentSystem].Atoms[CurrentFramework][D].Type;
                      if((D!=A)&&(D!=B)&&(CurrentTypeD==TypeD))
                      {
                        present=FALSE;
                        for(n=0;n<index;n++)
                        {
                          if((((Framework[CurrentSystem].BendTorsions[CurrentFramework][n].A==A)&&
                              (Framework[CurrentSystem].BendTorsions[CurrentFramework][n].B==B)&&
                              (Framework[CurrentSystem].BendTorsions[CurrentFramework][n].C==C)&&
                              (Framework[CurrentSystem].BendTorsions[CurrentFramework][n].D==D))||
                             ((Framework[CurrentSystem].BendTorsions[CurrentFramework][n].A==D)&&
                              (Framework[CurrentSystem].BendTorsions[CurrentFramework][n].B==C)&&
                              (Framework[CurrentSystem].BendTorsions[CurrentFramework][n].C==B)&&
                              (Framework[CurrentSystem].BendTorsions[CurrentFramework][n].D==A)))&&
                              (Framework[CurrentSystem].BendTorsionType[CurrentFramework][n]==BendTorsionType))
                             present=TRUE;
                        }
                        if(!present)
                        {
                          Framework[CurrentSystem].BendTorsions[CurrentFramework][index].A=A;
                          Framework[CurrentSystem].BendTorsions[CurrentFramework][index].B=B;
                          Framework[CurrentSystem].BendTorsions[CurrentFramework][index].C=C;
                          Framework[CurrentSystem].BendTorsions[CurrentFramework][index].D=D;

                          Framework[CurrentSystem].BendTorsionType[CurrentFramework][index]=BendTorsionType;
                          Framework[CurrentSystem].NumberOfBendTorsionsPerType[i]++;

                          for(j=0;j<BendTorsionTypes[BendTorsionType].nr_args;j++)
                            Framework[CurrentSystem].BendTorsionArguments[CurrentFramework][index][j]=arguments[j];

                          switch(BendTorsionType)
                          {
                            case SMOOTHED_DIHEDRAL:
                              // S(Theta1)*[p_0(1+cos(p_1*Phi-p_2)]*S(Theta2)
                              // ======================================================================================
                              // p_0/k_B [K/rad^2]
                              // p_1     [-]
                              // p_2     [degrees]
                              Framework[CurrentSystem].BendTorsionArguments[CurrentFramework][index][0]*=KELVIN_TO_ENERGY;
                              Framework[CurrentSystem].BendTorsionArguments[CurrentFramework][index][2]*=DEG2RAD;
                              break;
                            case SMOOTHED_THREE_COSINE_DIHEDRAL:
                              // S(Theta1)*[(1/2)*(1+cos(Phi))+(1/2)*p_1*(1-cos(2*Phi))+(1/2)*(1+cos(3*Phi))]*S(Theta2)
                              // ======================================================================================
                              // p_0/k_B [K]
                              // p_1/k_B [K]
                              // p_2/k_B [K]
                              Framework[CurrentSystem].BendTorsionArguments[CurrentFramework][index][0]*=KELVIN_TO_ENERGY;
                              Framework[CurrentSystem].BendTorsionArguments[CurrentFramework][index][1]*=KELVIN_TO_ENERGY;
                              Framework[CurrentSystem].BendTorsionArguments[CurrentFramework][index][2]*=KELVIN_TO_ENERGY;
                              break;
                            case NICHOLAS_DIHEDRAL:
                              // S(Theta1)*[(1/2)*(1+cos(Phi))+(1/2)*p_1*(1-cos(2*Phi))+(1/2)*(1+cos(3*Phi))]*S(Theta2)
                              // ======================================================================================
                              // p_0/k_B [K]
                              // p_1/k_B [K]
                              // p_2/k_B [K]
                              Framework[CurrentSystem].BendTorsionArguments[CurrentFramework][index][0]*=KELVIN_TO_ENERGY;
                              Framework[CurrentSystem].BendTorsionArguments[CurrentFramework][index][1]*=KELVIN_TO_ENERGY;
                              Framework[CurrentSystem].BendTorsionArguments[CurrentFramework][index][2]*=KELVIN_TO_ENERGY;
                              break;
                            case SMOOTHED_CFF_DIHEDRAL:
                              // S(Theta1)*[(1-cos(Phi))+p_1*(1-cos(2*Phi))+(1-cos(3*Phi))]*S(Theta2)
                              // ======================================================================================
                              // p_0/k_B [K]
                              // p_1/k_B [K]
                              // p_2/k_B [K]
                              Framework[CurrentSystem].BendTorsionArguments[CurrentFramework][index][0]*=KELVIN_TO_ENERGY;
                              Framework[CurrentSystem].BendTorsionArguments[CurrentFramework][index][1]*=KELVIN_TO_ENERGY;
                              Framework[CurrentSystem].BendTorsionArguments[CurrentFramework][index][2]*=KELVIN_TO_ENERGY;
                              break;
                            case SMOOTHED_CFF_DIHEDRAL2:
                              // S(Theta1)*[(1+cos(Phi))+p_1*(1+cos(2*Phi))+(1+cos(3*Phi))]*S(Theta2)
                              // ======================================================================================
                              // p_0/k_B [K]
                              // p_1/k_B [K]
                              // p_2/k_B [K]
                              Framework[CurrentSystem].BendTorsionArguments[CurrentFramework][index][0]*=KELVIN_TO_ENERGY;
                              Framework[CurrentSystem].BendTorsionArguments[CurrentFramework][index][1]*=KELVIN_TO_ENERGY;
                              Framework[CurrentSystem].BendTorsionArguments[CurrentFramework][index][2]*=KELVIN_TO_ENERGY;
                              break;
                            case CVFF_BEND_TORSION_CROSS:
                            case CFF_BEND_TORSION_CROSS:
                              // p_0*(Theta1-p_1)*(Theta2-p_2)*cos(Phi)
                              // =====================================================================================
                              // p_0/k_B [K/rad^3]
                              // p_1     [degrees]
                              // p_2     [degrees]
                              Framework[CurrentSystem].BendTorsionArguments[CurrentFramework][index][0]*=KELVIN_TO_ENERGY;
                              Framework[CurrentSystem].BendTorsionArguments[CurrentFramework][index][1]*=DEG2RAD;
                              Framework[CurrentSystem].BendTorsionArguments[CurrentFramework][index][2]*=DEG2RAD;
                              break;
                            case SMOOTHED_CFF_BEND_TORSION_CROSS:
                              // S(Theta1)*[p_0*(Theta1-p_1)*(Theta2-p_2)*cos(Phi)]*S(Theta2)
                              // ======================================================================================
                              // p_0/k_B [K/rad^3]
                              // p_1     [degrees]
                              // p_2     [degrees]
                              Framework[CurrentSystem].BendTorsionArguments[CurrentFramework][index][0]*=KELVIN_TO_ENERGY;
                              Framework[CurrentSystem].BendTorsionArguments[CurrentFramework][index][1]*=DEG2RAD;
                              Framework[CurrentSystem].BendTorsionArguments[CurrentFramework][index][2]*=DEG2RAD;
                              break;
                            default:
                              fprintf(stderr, "Undefined Bend-Torsion potential in routine 'ReadFrameworkDefinition' ('framework.c')\n");
                              exit(0);
                              break;
                          }
                          Framework[CurrentSystem].NumberOfBendTorsions[CurrentFramework]++;
                          index=Framework[CurrentSystem].NumberOfBendTorsions[CurrentFramework];

                          if(index>=Framework[CurrentSystem].MaxNumberOfBendTorsions[CurrentFramework])
                          {
                            Framework[CurrentSystem].MaxNumberOfBendTorsions[CurrentFramework]+=4096;
                            Framework[CurrentSystem].BendTorsions[CurrentFramework]=(QUAD*)realloc(Framework[CurrentSystem].BendTorsions[CurrentFramework],
                                      Framework[CurrentSystem].MaxNumberOfBendTorsions[CurrentFramework]*sizeof(QUAD));
                            Framework[CurrentSystem].BendTorsionType[CurrentFramework]=(int*)realloc(Framework[CurrentSystem].BendTorsionType[CurrentFramework],
                                      Framework[CurrentSystem].MaxNumberOfBendTorsions[CurrentFramework]*sizeof(int));
                            Framework[CurrentSystem].BendTorsionArguments[CurrentFramework]=(REAL(*)[MAX_BEND_TORSION_POTENTIAL_ARGUMENTS])
                                 realloc(Framework[CurrentSystem].BendTorsionArguments[CurrentFramework],
                            Framework[CurrentSystem].MaxNumberOfBendTorsions[CurrentFramework]*sizeof(REAL[MAX_BEND_TORSION_POTENTIAL_ARGUMENTS]));
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  fprintf(stderr, "Number of bonds: %d %d %d\n",Framework[0].NumberOfBonds[0],CurrentSystem,CurrentFramework);
  return 0;
}

int IsDefinedBondType(int system,int f1,int A,int B)
{
  int i;
  int typeA,typeB;

  A%=Framework[system].NumberOfAtoms[f1];
  B%=Framework[system].NumberOfAtoms[f1];
  typeA=Framework[system].Atoms[f1][A].Type;
  typeB=Framework[system].Atoms[f1][B].Type;
  for(i=0;i<Framework[system].NumberOfBondsDefinitions;i++)
  {
    if(((Framework[system].BondDefinitions[i].A==typeA)&&(Framework[system].BondDefinitions[i].B==typeB))||
       ((Framework[system].BondDefinitions[i].A==typeB)&&(Framework[system].BondDefinitions[i].B==typeA)))
       return TRUE;
  }
  return FALSE;
}

int IsDefinedBendType(int system,int f1,int A,int B,int C)
{
  int i;
  int typeA,typeB,typeC;

  A%=Framework[system].NumberOfAtoms[f1];
  B%=Framework[system].NumberOfAtoms[f1];
  C%=Framework[system].NumberOfAtoms[f1];
  typeA=Framework[system].Atoms[f1][A].Type;
  typeB=Framework[system].Atoms[f1][B].Type;
  typeC=Framework[system].Atoms[f1][C].Type;
  for(i=0;i<Framework[system].NumberOfBendDefinitions;i++)
  {
    if(((Framework[system].BendDefinitions[i].A==typeA)&&(Framework[system].BendDefinitions[i].B==typeB)&&
       (Framework[system].BendDefinitions[i].C==typeC))||((Framework[system].BendDefinitions[i].A==typeC)&&
       (Framework[system].BendDefinitions[i].B==typeB)&&(Framework[system].BendDefinitions[i].C==typeA)))
       return TRUE;
  }
  return FALSE;
}

int IsDefinedTorsionType(int system,int f1,int A,int B,int C,int D)
{
  int i;
  int typeA,typeB,typeC,typeD;

  typeA=Framework[system].Atoms[f1][A].Type;
  typeB=Framework[system].Atoms[f1][B].Type;
  typeC=Framework[system].Atoms[f1][C].Type;
  typeD=Framework[system].Atoms[f1][D].Type;
  for(i=0;i<Framework[system].NumberOfTorsionDefinitions;i++)
  {
    if(((Framework[system].TorsionDefinitions[i].A==typeA)&&(Framework[system].TorsionDefinitions[i].B==typeB)&&
       (Framework[system].TorsionDefinitions[i].C==typeC)&&(Framework[system].TorsionDefinitions[i].D==typeD))||
       ((Framework[system].TorsionDefinitions[i].A==typeD)&&(Framework[system].TorsionDefinitions[i].B==typeC)&&
       (Framework[system].TorsionDefinitions[i].C==typeB)&&(Framework[system].TorsionDefinitions[i].D==typeA)))
       return TRUE;
  }
  return FALSE;
}

int IsDefinedBond(int system,int f1,int A,int B)
{
  int i;

  A%=Framework[system].NumberOfAtoms[f1];
  B%=Framework[system].NumberOfAtoms[f1];

  for(i=0;i<Framework[system].NumberOfBonds[f1];i++)
  {
    if(((Framework[system].Bonds[f1][i].A==A)&&(Framework[system].Bonds[f1][i].B==B))||
       ((Framework[system].Bonds[f1][i].A==B)&&(Framework[system].Bonds[f1][i].B==A)))
       return TRUE;
  }

  return FALSE;
}

int IsDefinedBend(int system,int f1,int A,int B,int C)
{
  int i;

  A%=Framework[system].NumberOfAtoms[f1];
  B%=Framework[system].NumberOfAtoms[f1];
  C%=Framework[system].NumberOfAtoms[f1];

  for(i=0;i<Framework[system].NumberOfBends[f1];i++)
  {
    if(((Framework[system].Bends[f1][i].A==A)&&(Framework[system].Bends[f1][i].B==B)&&(Framework[system].Bends[f1][i].C==C))||
       ((Framework[system].Bends[f1][i].A==C)&&(Framework[system].Bends[f1][i].B==B)&&(Framework[system].Bends[f1][i].C==A)))
       return TRUE;
  }

  return FALSE;
}

int IsDefinedTorsion(int system,int f1,int A,int B,int C,int D)
{
  int i;

  A%=Framework[system].NumberOfAtoms[f1];
  B%=Framework[system].NumberOfAtoms[f1];
  C%=Framework[system].NumberOfAtoms[f1];
  D%=Framework[system].NumberOfAtoms[f1];
  for(i=0;i<Framework[system].NumberOfTorsions[f1];i++)
  {
    if(((Framework[system].Torsions[f1][i].A==A)&&(Framework[system].Torsions[f1][i].B==B)&&(Framework[system].Torsions[f1][i].C==C)&&(Framework[system].Torsions[f1][i].D==D))||
       ((Framework[system].Torsions[f1][i].A==D)&&(Framework[system].Torsions[f1][i].B==C)&&(Framework[system].Torsions[f1][i].C==B)&&(Framework[system].Torsions[f1][i].D==A)))
       return TRUE;
  }
  return FALSE;
}

int IsDefinedImproperTorsion(int system,int f1,int A,int B,int C,int D)
{
  int i;

  A%=Framework[system].NumberOfAtoms[f1];
  B%=Framework[system].NumberOfAtoms[f1];
  C%=Framework[system].NumberOfAtoms[f1];
  D%=Framework[system].NumberOfAtoms[f1];
  for(i=0;i<Framework[system].NumberOfImproperTorsions[f1];i++)
  {
    if(((Framework[system].ImproperTorsions[f1][i].A==A)&&(Framework[system].ImproperTorsions[f1][i].B==B)&&(Framework[system].ImproperTorsions[f1][i].C==C)&&(Framework[system].ImproperTorsions[f1][i].D==D))||
       ((Framework[system].ImproperTorsions[f1][i].A==D)&&(Framework[system].ImproperTorsions[f1][i].B==C)&&(Framework[system].ImproperTorsions[f1][i].C==B)&&(Framework[system].ImproperTorsions[f1][i].D==A)))
       return TRUE;
  }
  return FALSE;
}

/*********************************************************************************************************
 * Name       | MakeExclusionMatrix                                                                      *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Creates a binary exclusion-matrix for VDW, charge, charge-bonddipole, and                *
 *            | bonddipole-bondipole.                                                                    *
 * Note       | Bit 0: 1 omit VDW interaction, 0 pair has VDW interaction                                *
 *            | Bit 1: 1 omit charge-charge interaction, 0 pair has charge-charge interaction            *
 *            | Bit 2: 1 omit charge-bonddipole interaction, 0 pair has charge-bonddipole interaction    *
 *            | Bit 3: 1 omit bondipole-bonddipole interaction, 0 pair has bonddipole-bonddipole inter.  *
 *            | Bit 4,5,6: same as 1,2,3 but for the unit cell (not the replica cell)                    *
 *********************************************************************************************************/

void MakeExclusionMatrix(int system)
{
  int i,j,k,l,m,f1;
  int A,B,C,D;
  int index1,index2;
  int largest_size;

  // allocate memory for the exclusion-matrix
  Framework[system].ExclusionMatrix=(char***)calloc(Framework[system].NumberOfFrameworks,sizeof(char**));
  for(f1=0;f1<Framework[system].NumberOfFrameworks;f1++)
  {
    largest_size=MAX2(Framework[system].NumberOfAtoms[f1],Framework[system].NumberOfBonds[f1]);

    Framework[system].ExclusionMatrix[f1]=(char**)calloc(largest_size,sizeof(char*));
    for(i=0;i<largest_size;i++)
    {
      Framework[system].ExclusionMatrix[f1][i]=(char*)calloc(TotalNumberOfReplicaCells[system]*largest_size,sizeof(char));
      for(j=0;j<TotalNumberOfReplicaCells[system]*largest_size;j++)
        Framework[system].ExclusionMatrix[f1][i][j]=0;
    }
  }

  Framework[system].NumberOfIntra12Interactions=0;
  Framework[system].NumberOfIntra13Interactions=0;
  Framework[system].NumberOfIntra14Interactions=0;
  Framework[system].NumberOfIntra123Interactions=0;
  Framework[system].NumberOfIntra1234Interactions=0;

  // General: count all possible 1-2, 1-3, and 1-4 interactions in the unit cell
  // ==============================================================================================
  for(f1=0;f1<Framework[system].NumberOfFrameworks;f1++)
  {

    for(A=0;A<Framework[system].NumberOfAtoms[f1];A++)
    {
      for(k=0;k<Framework[system].Connectivity[f1][A];k++)
      {
        B=GetNeighbour(system,f1,A,k);

        if(A<B) 
        {
          Framework[system].NumberOfIntra12Interactions++;
          Framework[system].NumberOfIntra123Interactions++;
          Framework[system].NumberOfIntra1234Interactions++;
        }

        for(l=0;l<Framework[system].Connectivity[f1][B];l++)
        {
          C=GetNeighbour(system,f1,B,l);
          if((A!=C)&&(C!=B))
          {
            // always mark as a 1-3 interaction
            if(A<C)
            {
              Framework[system].NumberOfIntra13Interactions++;
              Framework[system].NumberOfIntra123Interactions++;
              Framework[system].NumberOfIntra1234Interactions++;
            }

            for(m=0;m<Framework[system].Connectivity[f1][C];m++)
            {
              D=GetNeighbour(system,f1,C,m);
              if((D!=B)&&(D!=A)&&(D!=C))
              {
                if(A<D)
                {
                  Framework[system].NumberOfIntra14Interactions++;
                  Framework[system].NumberOfIntra1234Interactions++;
                }
              }
            }
          }
        }
      }
    }
  }


  // VDW-VDW
  // ==============================================================================================
  // set general 1-2, 1-3, 1-4 exclusions for VDW
  // VDW is based on replica-method using 'GetReplicaNeighbour'
  for(f1=0;f1<Framework[system].NumberOfFrameworks;f1++)
  {
    for(A=0;A<Framework[system].NumberOfAtoms[f1];A++)
    {
      SETBIT(Framework[system].ExclusionMatrix[f1][A][A],0);
      for(k=0;k<Framework[system].Connectivity[f1][A];k++)
      {
        B=GetReplicaNeighbour(system,f1,A,k);

        // we now have 2 connected atoms: A and B
        if(Remove12NeighboursFromVDWInteraction||(RemoveBondNeighboursFromLongRangeInteraction&&IsDefinedBond(system,f1,A,B)))
        {
          SETBIT(Framework[system].ExclusionMatrix[f1][A][B],0);
          if(B<Framework[system].NumberOfAtoms[f1])
            SETBIT(Framework[system].ExclusionMatrix[f1][B][A],0);
        }

        for(l=0;l<Framework[system].Connectivity[f1][B];l++)
        {
          C=GetReplicaNeighbour(system,f1,B,l);
          if(A!=C)
          {
            // we now have 3 connected atoms: A, B and C
            if(Remove13NeighboursFromVDWInteraction||(RemoveBendNeighboursFromLongRangeInteraction&&IsDefinedBend(system,f1,A,B,C)))
            {
              SETBIT(Framework[system].ExclusionMatrix[f1][A][C],0);
              if(C<Framework[system].NumberOfAtoms[f1])
                SETBIT(Framework[system].ExclusionMatrix[f1][C][A],0);
            }
          }

          for(m=0;m<Framework[system].Connectivity[f1][C];m++)
          {
            D=GetReplicaNeighbour(system,f1,C,m);
            if((D!=B)&&(D!=A))
            {
              // we now have 4 connected atoms: A, B, C, D
              if(Remove14NeighboursFromVDWInteraction||(RemoveTorsionNeighboursFromLongRangeInteraction&&IsDefinedTorsion(system,f1,A,B,C,D)))
              {
                SETBIT(Framework[system].ExclusionMatrix[f1][A][D],0);
                if(D<Framework[system].NumberOfAtoms[f1])
                  SETBIT(Framework[system].ExclusionMatrix[f1][D][A],0);
              }
            }
          }
        }
      }
    }
  }

  // Charge-charge
  // ==============================================================================================

  // set general 1-2, 1-3, 1-4 exclusions for charge-charge interactions
  // charge-charge is based on replica-method using 'GetReplicaNeighbour'
  for(f1=0;f1<Framework[system].NumberOfFrameworks;f1++)
  {
    for(A=0;A<Framework[system].NumberOfAtoms[f1];A++)
    {
      SETBIT(Framework[system].ExclusionMatrix[f1][A][A],1);
      for(k=0;k<Framework[system].Connectivity[f1][A];k++)
      {
        B=GetReplicaNeighbour(system,f1,A,k);

        // we now have 2 connected atoms: A and B
        if(Remove12NeighboursFromChargeChargeInteraction||(RemoveBondNeighboursFromLongRangeInteraction&&IsDefinedBond(system,f1,A,B)))
        {
          SETBIT(Framework[system].ExclusionMatrix[f1][A][B],1);
          if(B<Framework[system].NumberOfAtoms[f1])
            SETBIT(Framework[system].ExclusionMatrix[f1][B][A],1);
        }
        for(l=0;l<Framework[system].Connectivity[f1][B];l++)
        {
          C=GetReplicaNeighbour(system,f1,B,l);
          if(A!=C)
          {
            // we now have 3 connected atoms: A, B and C
            if(Remove13NeighboursFromChargeChargeInteraction||(RemoveBendNeighboursFromLongRangeInteraction&&IsDefinedBend(system,f1,A,B,C)))
            {
              SETBIT(Framework[system].ExclusionMatrix[f1][A][C],1);
              if(C<Framework[system].NumberOfAtoms[f1])
                SETBIT(Framework[system].ExclusionMatrix[f1][C][A],1);
            }
          }

          for(m=0;m<Framework[system].Connectivity[f1][C];m++)
          {
            D=GetReplicaNeighbour(system,f1,C,m);
            if((D!=B)&&(D!=A))
            {
              // we now have 4 connected atoms: A, B, C, D
              if(Remove14NeighboursFromChargeChargeInteraction||(RemoveTorsionNeighboursFromLongRangeInteraction&&IsDefinedTorsion(system,f1,A,B,C,D)))
              {
                SETBIT(Framework[system].ExclusionMatrix[f1][A][D],1);
                if(D<Framework[system].NumberOfAtoms[f1])
                  SETBIT(Framework[system].ExclusionMatrix[f1][D][A],1);
              }
            }
          }
        }
      }
    }
  }

  // set general 1-2, 1-3, 1-4 exclusions for charge-charge exclusion-list
  // charge exclusion in Fourier-space is based on the unit-cell using 'GetNeighbour'
  for(f1=0;f1<Framework[system].NumberOfFrameworks;f1++)
  {
    for(A=0;A<Framework[system].NumberOfAtoms[f1];A++)
    {
      SETBIT(Framework[system].ExclusionMatrix[f1][A][A],4);
      for(k=0;k<Framework[system].Connectivity[f1][A];k++)
      {
        B=GetNeighbour(system,f1,A,k);

        // we now have 2 connected atoms: A and B
        if(Remove12NeighboursFromChargeChargeInteraction||(RemoveBondNeighboursFromLongRangeInteraction&&IsDefinedBond(system,f1,A,B)))
        {
          SETBIT(Framework[system].ExclusionMatrix[f1][A][B],4);
          if(B<Framework[system].NumberOfAtoms[f1])
            SETBIT(Framework[system].ExclusionMatrix[f1][B][A],4);
        }
        for(l=0;l<Framework[system].Connectivity[f1][B];l++)
        {
          C=GetNeighbour(system,f1,B,l);
          if(A!=C)
          {
            // we now have 3 connected atoms: A, B and C
            if(Remove13NeighboursFromChargeChargeInteraction||(RemoveBendNeighboursFromLongRangeInteraction&&IsDefinedBend(system,f1,A,B,C)))
            {
              SETBIT(Framework[system].ExclusionMatrix[f1][A][C],4);
              if(C<Framework[system].NumberOfAtoms[f1])
                SETBIT(Framework[system].ExclusionMatrix[f1][C][A],4);
            }
          }

          for(m=0;m<Framework[system].Connectivity[f1][C];m++)
          {
            D=GetNeighbour(system,f1,C,m);
            if((D!=B)&&(D!=A))
            {
              // we now have 4 connected atoms: A, B, C, D
              if(Remove14NeighboursFromChargeChargeInteraction||(RemoveTorsionNeighboursFromLongRangeInteraction&&IsDefinedTorsion(system,f1,A,B,C,D)))
              {
                SETBIT(Framework[system].ExclusionMatrix[f1][A][D],4);
                if(D<Framework[system].NumberOfAtoms[f1])
                  SETBIT(Framework[system].ExclusionMatrix[f1][D][A],4);
              }
            }
          }
        }
      }
    }
  }

  // Charge-bonddipole
  // ==============================================================================================

  // set general 1-2, 1-3, 1-4 exclusions for charge-bonddipole interactions
  // charge-bonddipole is based on replica-method using 'GetReplicaNeighbour'
  for(f1=0;f1<Framework[system].NumberOfFrameworks;f1++)
  {
    for(A=0;A<Framework[system].NumberOfAtoms[f1];A++)
    {
      for(k=0;k<Framework[system].Connectivity[f1][A];k++)
      {
        B=GetReplicaNeighbour(system,f1,A,k);

        // we now have 2 connected atoms: A and B
        index1=ReturnDipoleIndex(f1,A,B);
        if(Remove11NeighboursFromChargeBondDipoleInteraction&&(index1>=0))
        {
          SETBIT(Framework[system].ExclusionMatrix[f1][A][index1],2);
          if(B<Framework[system].NumberOfAtoms[f1])
            SETBIT(Framework[system].ExclusionMatrix[f1][B][index1],2);
        }
        for(l=0;l<Framework[system].Connectivity[f1][B];l++)
        {
          C=GetReplicaNeighbour(system,f1,B,l);
          if(A!=C)
          {
            // we now have 3 connected atoms: A, B and C
            index1=ReturnDipoleIndex(f1,A,B);
            if(Remove12NeighboursFromChargeBondDipoleInteraction&&(index1>=0))
            {
              if(C<Framework[system].NumberOfAtoms[f1])
                SETBIT(Framework[system].ExclusionMatrix[f1][C][index1],2);
            }
            index1=ReturnDipoleIndex(f1,B,C);
            if(Remove12NeighboursFromChargeBondDipoleInteraction&&(index1>=0))
              SETBIT(Framework[system].ExclusionMatrix[f1][A][index1],2);
          }

          for(m=0;m<Framework[system].Connectivity[f1][C];m++)
          {
            D=GetReplicaNeighbour(system,f1,C,m);
            if((D!=B)&&(D!=A))
            {
              // we now have 4 connected atoms: A, B, C, D
              index1=ReturnDipoleIndex(f1,A,B);
              if(Remove13NeighboursFromChargeBondDipoleInteraction&&(index1>=0))
                SETBIT(Framework[system].ExclusionMatrix[f1][D][index1],2);
              index1=ReturnDipoleIndex(f1,C,D);
              if(Remove13NeighboursFromChargeBondDipoleInteraction&&(index1>=0))
                SETBIT(Framework[system].ExclusionMatrix[f1][A][index1],2);
            }
          }
        }
      }
    }
  }

  // set general 1-2, 1-3, 1-4 exclusions for charge-bonddipole exclusion-list
  // charge-bonddipole exclusion in Fourier-space is based on the unit-cell using 'GetNeighbour'
  for(f1=0;f1<Framework[system].NumberOfFrameworks;f1++)
  {
    for(A=0;A<Framework[system].NumberOfAtoms[f1];A++)
    {
      for(k=0;k<Framework[system].Connectivity[f1][A];k++)
      {
        B=GetNeighbour(system,f1,A,k);

        // we now have 2 connected atoms: A and B
        index1=ReturnDipoleIndex(f1,A,B);
        if(Remove11NeighboursFromChargeBondDipoleInteraction&&(index1>=0))
        {
          SETBIT(Framework[system].ExclusionMatrix[f1][A][index1],5);
          if(B<Framework[system].NumberOfAtoms[f1])
            SETBIT(Framework[system].ExclusionMatrix[f1][B][index1],5);
        }
        for(l=0;l<Framework[system].Connectivity[f1][B];l++)
        {
          C=GetNeighbour(system,f1,B,l);
          if(A!=C)
          {
            // we now have 3 connected atoms: A, B and C
            index1=ReturnDipoleIndex(f1,A,B);
            if(Remove12NeighboursFromChargeBondDipoleInteraction&&(index1>=0))
            {
              if(C<Framework[system].NumberOfAtoms[f1])
                SETBIT(Framework[system].ExclusionMatrix[f1][C][index1],5);
            }
            index1=ReturnDipoleIndex(f1,B,C);
            if(Remove12NeighboursFromChargeBondDipoleInteraction&&(index1>=0))
              SETBIT(Framework[system].ExclusionMatrix[f1][A][index1],5);
          }

          for(m=0;m<Framework[system].Connectivity[f1][C];m++)
          {
            D=GetNeighbour(system,f1,C,m);
            if((D!=B)&&(D!=A))
            {
              // we now have 4 connected atoms: A, B, C, D
              index1=ReturnDipoleIndex(f1,A,B);
              if(Remove13NeighboursFromChargeBondDipoleInteraction&&(index1>=0))
                SETBIT(Framework[system].ExclusionMatrix[f1][D][index1],5);
              index1=ReturnDipoleIndex(f1,C,D);
              if(Remove13NeighboursFromChargeBondDipoleInteraction&&(index1>=0))
                SETBIT(Framework[system].ExclusionMatrix[f1][A][index1],5);
            }
          }
        }
      }
    }
  }


  // Bonddipole-bonddipole
  // ==============================================================================================

  // set general 1-2, 1-3, 1-4 exclusions for bonddipole-bonddipole interactions
  // bonddipole-bonddipole is based on replica-method using 'GetReplicaNeighbour'
  for(f1=0;f1<Framework[system].NumberOfFrameworks;f1++)
  {
    for(A=0;A<Framework[system].NumberOfAtoms[f1];A++)
    {
      SETBIT(Framework[system].ExclusionMatrix[f1][A][A],3);
      for(k=0;k<Framework[system].Connectivity[f1][A];k++)
      {
        B=GetReplicaNeighbour(system,f1,A,k);

        for(l=0;l<Framework[system].Connectivity[f1][B];l++)
        {
          C=GetReplicaNeighbour(system,f1,B,l);
          if(A!=C)
          {
            // we now have 3 connected atoms: A, B and C
            index1=ReturnDipoleIndex(f1,A,B);
            index2=ReturnDipoleIndex(f1,B,C);
            if(Remove12NeighboursFromBondDipoleBondDipoleInteraction&&(index1>=0)&&(index2>=0))
            {
              SETBIT(Framework[system].ExclusionMatrix[f1][index1][index2],3);
              SETBIT(Framework[system].ExclusionMatrix[f1][index2][index1],3);
            }
          }

          for(m=0;m<Framework[system].Connectivity[f1][C];m++)
          {
            D=GetReplicaNeighbour(system,f1,C,m);
            if((D!=B)&&(D!=A))
            {
              // we now have 4 connected atoms: A, B, C, D
              index1=ReturnDipoleIndex(f1,A,B);
              index2=ReturnDipoleIndex(f1,C,D);
              if(Remove13NeighboursFromBondDipoleBondDipoleInteraction&&(index1>=0)&&(index2>=0))
              {
                SETBIT(Framework[system].ExclusionMatrix[f1][index1][index2],3);
                SETBIT(Framework[system].ExclusionMatrix[f1][index2][index1],3);
              }
            }
          }
        }
      }
    }
  }

  // set general 1-2, 1-3, 1-4 exclusions for bonddipole-bonddipole exclusion-list
  // bonddipole-bonddipole exclusion in Fourier-space is based on the unit-cell using 'GetNeighbour'
  for(f1=0;f1<Framework[system].NumberOfFrameworks;f1++)
  {
    for(A=0;A<Framework[system].NumberOfAtoms[f1];A++)
    {
      SETBIT(Framework[system].ExclusionMatrix[f1][A][A],6);
      for(k=0;k<Framework[system].Connectivity[f1][A];k++)
      {
        B=GetNeighbour(system,f1,A,k);

        for(l=0;l<Framework[system].Connectivity[f1][B];l++)
        {
          C=GetNeighbour(system,f1,B,l);
          if(A!=C)
          {
            // we now have 3 connected atoms: A, B and C
            index1=ReturnDipoleIndex(f1,A,B);
            index2=ReturnDipoleIndex(f1,B,C);
            if(Remove12NeighboursFromBondDipoleBondDipoleInteraction&&(index1>=0)&&(index2>=0))
            {
              SETBIT(Framework[system].ExclusionMatrix[f1][index1][index2],6);
              SETBIT(Framework[system].ExclusionMatrix[f1][index2][index1],6);
            }
          }

          for(m=0;m<Framework[system].Connectivity[f1][C];m++)
          {
            D=GetNeighbour(system,f1,C,m);
            if((D!=B)&&(D!=A))
            {
              // we now have 4 connected atoms: A, B, C, D
              index1=ReturnDipoleIndex(f1,A,B);
              index2=ReturnDipoleIndex(f1,C,D);
              if(Remove13NeighboursFromBondDipoleBondDipoleInteraction&&(index1>=0)&&(index2>=0))
              {
                SETBIT(Framework[system].ExclusionMatrix[f1][index1][index2],6);
                SETBIT(Framework[system].ExclusionMatrix[f1][index2][index1],6);
              }
            }
          }
        }
      }
    }
  }
}

void MakeExcludedInteractionLists(int system)
{
  int i,j,f1;
  int index_excluded;
  int TypeA,TypeB;

  // make lists for exclusion in the Ewald summation (subtraction using 'erf')
  // count the amount of interactions and exclusions

  for(f1=0;f1<Framework[system].NumberOfFrameworks;f1++)
  {
    // VDW-VDW
    // ==============================================================================================

    // count the number of intra VDW pairs and excluded VDW pairs
    Framework[system].NumberOfIntraVDW[f1]=0;
    Framework[system].NumberOfExcludedIntraVDW[f1]=0.0;
    for(i=0;i<Framework[system].NumberOfAtoms[f1];i++)
      for(j=i+1;j<TotalNumberOfReplicaCells[system]*Framework[system].NumberOfAtoms[f1];j++)
      {
        if(!BITVAL(Framework[system].ExclusionMatrix[f1][i][j],0))
          Framework[system].NumberOfIntraVDW[f1]++;
        else
          Framework[system].NumberOfExcludedIntraVDW[f1]++;
      }

    // Charge-charge
    // ==============================================================================================

    // count the number of intra charge-charge pairs
    Framework[system].NumberOfIntraCharges[f1]=0;
    for(i=0;i<Framework[system].NumberOfAtoms[f1];i++)
      for(j=i+1;j<TotalNumberOfReplicaCells[system]*Framework[system].NumberOfAtoms[f1];j++)
      {
        if(!BITVAL(Framework[system].ExclusionMatrix[f1][i][j],1))
          Framework[system].NumberOfIntraCharges[f1]++;
      }

    // count the number of Fourier excluded intra charge-charge pairs
    Framework[system].NumberOfExcludedIntraChargeCharge[f1]=0;
    for(i=0;i<Framework[system].NumberOfAtoms[f1];i++)
    {
      TypeA=Framework[system].Atoms[f1][i].Type;
      for(j=i+1;j<Framework[system].NumberOfAtoms[f1];j++)
      {
        TypeB=Framework[system].Atoms[f1][j].Type;
        if(BITVAL(Framework[system].ExclusionMatrix[f1][i][j],4))
          Framework[system].NumberOfExcludedIntraChargeCharge[f1]++;
      }
    }

    // allocate memory for Fourier excluded intra charge-charge pairs
    Framework[system].ExcludedIntraChargeCharge[f1]=(PAIR*)calloc(
        Framework[system].NumberOfExcludedIntraChargeCharge[f1],sizeof(PAIR));


    // create pair-list for Fourier excluded intra charge-charge pairs
    index_excluded=0;
    Framework[system].NumberOfExcludedIntraChargeCharge[f1]=0.0;
    for(i=0;i<Framework[system].NumberOfAtoms[f1];i++)
    {
      TypeA=Framework[system].Atoms[f1][i].Type;
      for(j=i+1;j<Framework[system].NumberOfAtoms[f1];j++)
      {
        TypeB=Framework[system].Atoms[f1][j].Type;
        if(BITVAL(Framework[system].ExclusionMatrix[f1][i][j],4))
        {
          // store A and B with (A<B)
          Framework[system].ExcludedIntraChargeCharge[f1][index_excluded].A=i;
          Framework[system].ExcludedIntraChargeCharge[f1][index_excluded].B=j;
          Framework[system].NumberOfExcludedIntraChargeCharge[f1]++;
          index_excluded++;
        }
      }
    }

    // Charge-bonddipole
    // ==============================================================================================

    // count the number of intra charge-bondipole pairs
    Framework[system].NumberOfIntraChargeBondDipoles[f1]=0;
    for(i=0;i<Framework[system].NumberOfAtoms[f1];i++)
    {
      for(j=i+1;j<Framework[system].NumberOfBondDipoles[f1];j++)
      {
        if(!BITVAL(Framework[system].ExclusionMatrix[f1][i][j],2))
          Framework[system].NumberOfIntraChargeBondDipoles[f1]++;
      }
    }

    // count the number of Fourier excluded intra charge-bonddipole pairs
    Framework[system].NumberOfExcludedIntraChargeBondDipole[f1]=0;
    for(i=0;i<Framework[system].NumberOfAtoms[f1];i++)
    {
      TypeA=Framework[system].Atoms[f1][i].Type;
      for(j=0;j<Framework[system].NumberOfBondDipoles[f1];j++)
      {
        if(BITVAL(Framework[system].ExclusionMatrix[f1][i][j],5))
          Framework[system].NumberOfExcludedIntraChargeBondDipole[f1]++;
      }
    }

    // allocate memory for Fourier excluded intra charge-bonddipole pairs
    Framework[system].ExcludedIntraChargeBondDipole[f1]=(PAIR*)calloc(
        Framework[system].NumberOfExcludedIntraChargeBondDipole[f1],sizeof(PAIR));

    // create pair-list for Fourier excluded intra charge-bonddipole pairs
    index_excluded=0;
    for(i=0;i<Framework[system].NumberOfAtoms[f1];i++)
    {
      TypeA=Framework[system].Atoms[f1][i].Type;
      for(j=0;j<Framework[system].NumberOfBondDipoles[f1];j++)
      {
        if(BITVAL(Framework[system].ExclusionMatrix[f1][i][j],5))
        {
          // store A and B with (A<B)
          Framework[system].ExcludedIntraChargeBondDipole[f1][index_excluded].A=i;
          Framework[system].ExcludedIntraChargeBondDipole[f1][index_excluded].B=j;
          index_excluded++;
        }
      }
    }

    // Bonddipole-bonddipole
    // ==============================================================================================

    // count the number of intra bonddipole-bondipole pairs
    Framework[system].NumberOfIntraBondDipoles[f1]=0;
    for(i=0;i<Framework[system].NumberOfBondDipoles[f1];i++)
    {
      for(j=i+1;j<Framework[system].NumberOfBondDipoles[f1];j++)
      {
        if(!BITVAL(Framework[system].ExclusionMatrix[f1][i][j],3))
          Framework[system].NumberOfIntraBondDipoles[f1]++;
      }
    }

    // count the number of Fourier excluded intra bonddipole-bonddipole pairs
    Framework[system].NumberOfExcludedIntraBondDipoleBondDipole[f1]=0;
    for(i=0;i<Framework[system].NumberOfBondDipoles[f1];i++)
    {
      for(j=i+1;j<Framework[system].NumberOfBondDipoles[f1];j++)
      {
        if(BITVAL(Framework[system].ExclusionMatrix[f1][i][j],6))
          Framework[system].NumberOfExcludedIntraBondDipoleBondDipole[f1]++;
      }
    }

    // allocate memory for Fourier excluded intra bonddipole-bonddipole pairs
    Framework[system].ExcludedIntraBondDipoleBondDipole[f1]=(PAIR*)calloc(
        Framework[system].NumberOfExcludedIntraBondDipoleBondDipole[f1],sizeof(PAIR));


    // create pair-list for Fourier excluded intra bonddipole-bonddipole pairs
    index_excluded=0;
    for(i=0;i<Framework[system].NumberOfBondDipoles[f1];i++)
    {
      for(j=i+1;j<Framework[system].NumberOfBondDipoles[f1];j++)
      {
        if(BITVAL(Framework[system].ExclusionMatrix[f1][i][j],6))
        {
          // store A and B with (A<B)
          Framework[system].ExcludedIntraBondDipoleBondDipole[f1][index_excluded].A=i;
          Framework[system].ExcludedIntraBondDipoleBondDipole[f1][index_excluded].B=j;
          index_excluded++;
        }
      }
    }
  }
}


void PrintFlexibleFrameworkModel(void)
{
  int i,j;
  int A,B,C,D,BondType,BendType,TorsionType;

  fprintf(stderr, "Number of bonds: %d\n",Framework[CurrentSystem].NumberOfBonds[CurrentFramework]);
  for(i=0;i<Framework[CurrentSystem].NumberOfBonds[CurrentFramework];i++)
  {
    A=Framework[CurrentSystem].Bonds[CurrentFramework][i].A;
    B=Framework[CurrentSystem].Bonds[CurrentFramework][i].B;
    BondType=Framework[CurrentSystem].BondType[CurrentFramework][i];
    fprintf(stderr, "Bond %s--%s [%d-%d]  Type: %d ",
            PseudoAtoms[Framework[CurrentSystem].Atoms[CurrentFramework][A].Type].Name,
            PseudoAtoms[Framework[CurrentSystem].Atoms[CurrentFramework][B].Type].Name,
            A,
            B,
            BondType);
    for(j=0;j<BondTypes[BondType].nr_args;j++)
      fprintf(stderr, " %lf",(double)Framework[CurrentSystem].BondArguments[CurrentFramework][i][j]);
    fprintf(stderr, "\n");
  }

  fprintf(stderr, "Number of Urey-Bradley: %d\n",Framework[CurrentSystem].NumberOfUreyBradleys[CurrentFramework]);
  for(i=0;i<Framework[CurrentSystem].NumberOfUreyBradleys[CurrentFramework];i++)
  {
    A=Framework[CurrentSystem].UreyBradleys[CurrentFramework][i].A;
    B=Framework[CurrentSystem].UreyBradleys[CurrentFramework][i].B;
    BondType=Framework[CurrentSystem].UreyBradleyType[CurrentFramework][i];
    fprintf(stderr, "Urey-Bradley %s--%s [%d-%d]  Type: %d ",
            PseudoAtoms[Framework[CurrentSystem].Atoms[CurrentFramework][A].Type].Name,
            PseudoAtoms[Framework[CurrentSystem].Atoms[CurrentFramework][B].Type].Name,
            A,
            B,
            BondType);
    for(j=0;j<BondTypes[BondType].nr_args;j++)
      fprintf(stderr, " %lf",(double)Framework[CurrentSystem].UreyBradleyArguments[CurrentFramework][i][j]);
    fprintf(stderr, "\n");
  }

  for(i=0;i<Framework[CurrentSystem].NumberOfBends[CurrentFramework];i++)
  {
    A=Framework[CurrentSystem].Bends[CurrentFramework][i].A;
    B=Framework[CurrentSystem].Bends[CurrentFramework][i].B;
    C=Framework[CurrentSystem].Bends[CurrentFramework][i].C;
    BendType=Framework[CurrentSystem].BendType[CurrentFramework][i];
    fprintf(stderr, "Bend %s--%s--%s [%d-%d-%d]  Type: %d ",
            PseudoAtoms[Framework[CurrentSystem].Atoms[CurrentFramework][A].Type].Name,
            PseudoAtoms[Framework[CurrentSystem].Atoms[CurrentFramework][B].Type].Name,
            PseudoAtoms[Framework[CurrentSystem].Atoms[CurrentFramework][C].Type].Name,
            A,
            B,
            C,
            BendType);
    for(j=0;j<BendTypes[BendType].nr_args;j++)
      fprintf(stderr, " %lf",(double)(RAD2DEG*(Framework[CurrentSystem].BendArguments[CurrentFramework][i][j])));
    fprintf(stderr, "\n");
  }

  fprintf(stderr, "Number of Torsions: %d\n",Framework[CurrentSystem].NumberOfTorsions[CurrentFramework]);
  for(i=0;i<Framework[CurrentSystem].NumberOfTorsions[CurrentFramework];i++)
  {
    A=Framework[CurrentSystem].Torsions[CurrentFramework][i].A;
    B=Framework[CurrentSystem].Torsions[CurrentFramework][i].B;
    C=Framework[CurrentSystem].Torsions[CurrentFramework][i].C;
    D=Framework[CurrentSystem].Torsions[CurrentFramework][i].D;
    TorsionType=Framework[CurrentSystem].TorsionType[CurrentFramework][i];
    fprintf(stderr, "Torsion %s--%s--%s--%s [%d-%d-%d-%d]  Type: %d ",
            PseudoAtoms[Framework[CurrentSystem].Atoms[CurrentFramework][A].Type].Name,
            PseudoAtoms[Framework[CurrentSystem].Atoms[CurrentFramework][B].Type].Name,
            PseudoAtoms[Framework[CurrentSystem].Atoms[CurrentFramework][C].Type].Name,
            PseudoAtoms[Framework[CurrentSystem].Atoms[CurrentFramework][D].Type].Name,
            A,
            B,
            C,
            D,
            TorsionType);
    for(j=0;j<TorsionTypes[TorsionType].nr_args;j++)
      fprintf(stderr, " %lf",(double)Framework[CurrentSystem].TorsionArguments[CurrentFramework][i][j]);
    fprintf(stderr, "\n");
  }
}

void QuenchCoreSHellVelocities(void)
{
  int f1;
  int A,B;
  int nr_core_shells,nr_atoms;
  REAL pke,rmu;
  REAL MassA,MassB,scale;
  VECTOR dv,tm;

  // permitted core-shell internal kinetic energy
  pke=K_B*therm_baro_stats.ExternalTemperature[CurrentSystem]*1e-4;

  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    nr_core_shells=Framework[CurrentSystem].NumberOfCoreShells[f1];
    nr_atoms=Framework[CurrentSystem].NumberOfAtoms[f1];
    if(nr_core_shells>0)
    {
      // loop over shells
      for(A=0;A<nr_atoms-nr_core_shells;A++)
      {
        // get the core for this shell
        B=Framework[CurrentSystem].CoreShellConnectivity[f1][A];
        if(B>0)
        {
          MassA=PseudoAtoms[Framework[CurrentSystem].Atoms[f1][A].Type].Mass;
          MassB=PseudoAtoms[Framework[CurrentSystem].Atoms[f1][B].Type].Mass;

          rmu=MassA*MassB/(MassA+MassB);

          dv.x=Framework[CurrentSystem].Atoms[f1][B].Velocity.x-Framework[CurrentSystem].Atoms[f1][A].Velocity.x;
          dv.y=Framework[CurrentSystem].Atoms[f1][B].Velocity.y-Framework[CurrentSystem].Atoms[f1][A].Velocity.y;
          dv.z=Framework[CurrentSystem].Atoms[f1][B].Velocity.z-Framework[CurrentSystem].Atoms[f1][A].Velocity.z;

          scale=sqrt(pke/(rmu*(SQR(dv.x)+SQR(dv.y)+SQR(dv.z))));

          tm.x=MassA*Framework[CurrentSystem].Atoms[f1][A].Velocity.x+
               MassB*Framework[CurrentSystem].Atoms[f1][B].Velocity.x;
          tm.y=MassA*Framework[CurrentSystem].Atoms[f1][A].Velocity.y+
               MassB*Framework[CurrentSystem].Atoms[f1][B].Velocity.y;
          tm.z=MassA*Framework[CurrentSystem].Atoms[f1][A].Velocity.z+
               MassB*Framework[CurrentSystem].Atoms[f1][B].Velocity.z;

          Framework[CurrentSystem].Atoms[f1][A].Velocity.x=tm.x/(MassA+MassB)-scale*rmu*dv.x/MassA;
          Framework[CurrentSystem].Atoms[f1][B].Velocity.x=tm.x/(MassA+MassB)+scale*rmu*dv.x/MassB;
          Framework[CurrentSystem].Atoms[f1][A].Velocity.y=tm.y/(MassA+MassB)-scale*rmu*dv.y/MassA;
          Framework[CurrentSystem].Atoms[f1][B].Velocity.y=tm.y/(MassA+MassB)+scale*rmu*dv.y/MassB;
          Framework[CurrentSystem].Atoms[f1][A].Velocity.z=tm.z/(MassA+MassB)-scale*rmu*dv.z/MassA;
          Framework[CurrentSystem].Atoms[f1][B].Velocity.z=tm.z/(MassA+MassB)+scale*rmu*dv.z/MassB;
        }
      }
    }
  }
}



int CheckSurfaceAreaOverlap(int typeA,VECTOR posA,int skipf,int skipfa,int skipcm,int skipca)
{
  int j,m;
  int typeB,f1;
  VECTOR dr,posB;
  REAL rr,well_depth_factor;

  if(BlockedPocket(posA)||(!ValidCartesianPoint(CurrentComponent,posA)))
    return TRUE;

  well_depth_factor=Framework[CurrentSystem].SurfaceAreaProbeDistance;

  // check overlap with framework
  for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
  {
    for(j=0;j<Framework[CurrentSystem].NumberOfAtoms[f1];j++)
    {
      if((f1!=skipf)||(j!=skipfa))
      {
        posB=Framework[CurrentSystem].Atoms[f1][j].Position;
        typeB=Framework[CurrentSystem].Atoms[f1][j].Type;
        dr.x=posA.x-posB.x;
        dr.y=posA.y-posB.y;
        dr.z=posA.z-posB.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

        if(rr<SQR(well_depth_factor*PotentialParms[typeA][typeB][1]))
          return TRUE;
      }
    }
  }

  // check overlap with possible cations
  for(m=0;m<NumberOfCationMolecules[CurrentSystem];m++)
  {
    for(j=0;j<Cations[CurrentSystem][m].NumberOfAtoms;j++)
    {
      if((m!=skipcm)||(j!=skipca))
      {
        posB=Cations[CurrentSystem][m].Atoms[j].Position;
        typeB=Cations[CurrentSystem][m].Atoms[j].Type;
        dr.x=posA.x-posB.x;
        dr.y=posA.y-posB.y;
        dr.z=posA.z-posB.z;
        dr=ApplyBoundaryCondition(dr);
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

        if(rr<SQR(well_depth_factor*PotentialParms[typeA][typeB][1]))
          return TRUE;
      }
    }
  }

  // check overlap with possible adsorbates
  for(m=0;m<NumberOfAdsorbateMolecules[CurrentSystem];m++)
  {
    for(j=0;j<Adsorbates[CurrentSystem][m].NumberOfAtoms;j++)
    {
      posB=Adsorbates[CurrentSystem][m].Atoms[j].Position;
      typeB=Adsorbates[CurrentSystem][m].Atoms[j].Type;
      dr.x=posA.x-posB.x;
      dr.y=posA.y-posB.y;
      dr.z=posA.z-posB.z;
      dr=ApplyBoundaryCondition(dr);
      rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

      if(rr<SQR(well_depth_factor*PotentialParms[typeA][typeB][1]))
        return TRUE;
    }
  }

  return FALSE;
}


void ConstructBondDipolesFromBondsFramework(void)
{
  //int i,A1,A2;
  //VECTOR Dipole,posA1,posA2;
  //REAL DipoleMagnitudeA,temp,length;
/*
  for(CurrentFramework=0;CurrentFramework<Framework[CurrentSystem].NumberOfFrameworks;CurrentFramework++)
  {
    for(i=0;i<Framework[CurrentSystem].NumberOfBondDipoles[CurrentFramework];i++)
    {
      DipoleMagnitudeA=Framework[CurrentSystem].BondDipoleMagnitude[CurrentFramework][i];
      A1=Framework[CurrentSystem].BondDipoles[CurrentFramework][i].A;
      A2=Framework[CurrentSystem].BondDipoles[CurrentFramework][i].B;
      posA1=Framework[CurrentSystem].Atoms[A1][CurrentFramework].Position;
      posA2=Framework[CurrentSystem].Atoms[A2][CurrentFramework].Position;
      Dipole.x=posA2.x-posA1.x;
      Dipole.y=posA2.y-posA1.y;
      Dipole.z=posA2.z-posA1.z;
      Dipole=ApplyBoundaryCondition(Dipole);
      length=sqrt(SQR(Dipole.x)+SQR(Dipole.y)+SQR(Dipole.z));
      if(fabs(length)>1e-10)
        temp=DipoleMagnitudeA/length;
      else
        temp=0.0;
      Dipole.x*=temp; Dipole.y*=temp; Dipole.z*=temp;
      Framework[CurrentSystem].BondDipoleVector[CurrentFramework][i]=Dipole;
      Framework[CurrentSystem].BondLength[CurrentFramework][i]=length;
    }
  }
*/
}

// returns TRUE when an asymmetric atom is find (it is stored in 'r'),
// otherwise FALSE is returned
int FindAsymmetricAtom(int sg,VECTOR s,VECTOR *r)
{
  int i;
  int k1,k2,k3;
  VECTOR t;

  s.x-=NINT(s.x);
  s.y-=NINT(s.y);
  s.z-=NINT(s.z);

  if(s.x<0.0) s.x+=1.0;
  if(s.y<0.0) s.y+=1.0;
  if(s.z<0.0) s.z+=1.0;

  // try all inverse transformation and find one that lies in the asymmetric unit cell
  SpaceGroupInverseSymmetry(sg,s);
  for(i=0;i<SpaceGroupSize;i++)
  {
    for(k1=-2;k1<=2;k1++)
    for(k2=-2;k2<=2;k2++)
    for(k3=-2;k3<=2;k3++)
    {
      t.x=SpaceGroupElement[i].x+(REAL)k1;
      t.y=SpaceGroupElement[i].y+(REAL)k2;
      t.z=SpaceGroupElement[i].z+(REAL)k3;
      if(AsymmetricUnit(sg,t,1e-6))
      {
        *r=t;
        return TRUE;
      }
    }
  }
  return FALSE;
}


/*********************************************************************************************************
 * Name       | DetermineSpaceGroup                                                                      *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | determines the space group of a P1-structure                                             *
 * Parameters | -                                                                                        *
 * Note       | Called when the option 'CalculateSpaceGroup' is set to 'yes'                             *
 *            | In case the structure was generated using a space group and one would still like to      *
 *            | redetermine the space group, use: 'ForceSpaceGroupDetection yes'.
 *********************************************************************************************************/

void DetermineSpaceGroup(void)
{
  int i,j,k;
  int Type,index;
  int sg;
  VECTOR pos,s,dr,t;
  REAL rr;
  int nr_atoms;
  int sg_possible;
  int FoundSpaceGroup;
  int LowestNumberOfAsymmetricUnits;

  // return if no framework present
  if(Framework[CurrentSystem].FrameworkModel==NONE) return;

  // return if the space group was already set in the structure file and if the spacegroup is not 'forced' to be determined
  // in this case the asymmetric atoms are already determined as read in from the structure file
  if((Framework[CurrentSystem].SpaceGroupIdentifier[CurrentFramework]!=1)&&(!Framework[CurrentSystem].ForceSpaceGroupDetection)) return;

  // space group P1 is always possible
  FoundSpaceGroup=1;

  // determine asymmetric atoms if "forced' to, or if space group is P1 and the framework is rigid
  // a flexible structure will always be set to P1
  if(((Framework[CurrentSystem].SpaceGroupIdentifier[CurrentFramework]==1)&&(Framework[CurrentSystem].FrameworkModels[CurrentFramework]!=FLEXIBLE))||
       Framework[CurrentSystem].ForceSpaceGroupDetection)
  {
    // space group P1 is always possible
    LowestNumberOfAsymmetricUnits=Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];

    Framework[CurrentSystem].AtomsAsymmetric[CurrentFramework]=(FRAMEWORK_ASYMMETRIC_ATOM*)realloc(Framework[CurrentSystem].AtomsAsymmetric[CurrentFramework],
                       Framework[CurrentSystem].NumberOfAtoms[CurrentFramework]*sizeof(FRAMEWORK_ASYMMETRIC_ATOM));

    // loop over all space-groups and try them one by one
    for(sg=1;sg<=NUMBER_OF_SPACEGROUPS;sg++)
    {
      // assume a certain space-group
      Framework[CurrentSystem].SpaceGroupIdentifier[CurrentFramework]=sg;

      Framework[CurrentSystem].NumberOfAsymmetricAtoms[CurrentFramework]=0;
      for(k=0;k<Framework[CurrentSystem].NumberOfUnitCellAtoms[CurrentFramework];k++)
      {
        Type=Framework[CurrentSystem].Atoms[CurrentFramework][k].Type;
        pos=Framework[CurrentSystem].Atoms[CurrentFramework][k].Position;

        s=ConvertFromXYZtoABC(pos);

        s.x*=NumberOfUnitCells[CurrentSystem].x;
        s.y*=NumberOfUnitCells[CurrentSystem].y;
        s.z*=NumberOfUnitCells[CurrentSystem].z;

        // if no asymmetric atom is found the loop over atoms is stopped, otherwise the found asymmetric atom is 't'
        if(!FindAsymmetricAtom(sg,s,&t))
        {
          Framework[CurrentSystem].NumberOfAsymmetricAtoms[CurrentFramework]=0;
          break;
        }

        index=-1;
        for(i=0;i<Framework[CurrentSystem].NumberOfAsymmetricAtoms[CurrentFramework];i++)
        {
          dr.x=Framework[CurrentSystem].AtomsAsymmetric[CurrentFramework][i].Position.x-t.x;
          dr.y=Framework[CurrentSystem].AtomsAsymmetric[CurrentFramework][i].Position.y-t.y;
          dr.z=Framework[CurrentSystem].AtomsAsymmetric[CurrentFramework][i].Position.z-t.z;
          rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
          if((rr<1e-8)&&(Type==Framework[CurrentSystem].AtomsAsymmetric[CurrentFramework][i].Type))
          {
            index=i;
            break;
          }
        }
        if(index<0)
        {
          Framework[CurrentSystem].AtomsAsymmetric[CurrentFramework][Framework[CurrentSystem].NumberOfAsymmetricAtoms[CurrentFramework]].Position=t;
          Framework[CurrentSystem].AtomsAsymmetric[CurrentFramework][Framework[CurrentSystem].NumberOfAsymmetricAtoms[CurrentFramework]].Type=Type;
          Framework[CurrentSystem].NumberOfAsymmetricAtoms[CurrentFramework]++;
        }
      }

      // we have the asymmetric atoms, let's form the framework atoms from these and check whether we get the same number of atoms
      nr_atoms=0;
      sg_possible=TRUE;
      for(i=0;sg_possible&&(i<Framework[CurrentSystem].NumberOfAsymmetricAtoms[CurrentFramework]);i++)
      {
        Type=Framework[CurrentSystem].AtomsAsymmetric[CurrentFramework][i].Type;
        s=Framework[CurrentSystem].AtomsAsymmetric[CurrentFramework][i].Position;
        SpaceGroupSymmetry(sg,s);

        for(j=0;sg_possible&&(j<SpaceGroupSize);j++)
        {
          // apply boundary condition
          SpaceGroupElement[j].x-=NINT(SpaceGroupElement[j].x);
          SpaceGroupElement[j].y-=NINT(SpaceGroupElement[j].y);
          SpaceGroupElement[j].z-=NINT(SpaceGroupElement[j].z);

          if(SpaceGroupElement[j].x<0.0) SpaceGroupElement[j].x+=1.0;
          if(SpaceGroupElement[j].y<0.0) SpaceGroupElement[j].y+=1.0;
          if(SpaceGroupElement[j].z<0.0) SpaceGroupElement[j].z+=1.0;

          index=-1;
          for(k=0;k<nr_atoms;k++)
          {
            dr.x=Framework[CurrentSystem].Atoms[CurrentFramework][k].ReferencePosition.x-SpaceGroupElement[j].x;
            dr.y=Framework[CurrentSystem].Atoms[CurrentFramework][k].ReferencePosition.y-SpaceGroupElement[j].y;
            dr.z=Framework[CurrentSystem].Atoms[CurrentFramework][k].ReferencePosition.z-SpaceGroupElement[j].z;
            dr.x-=NINT(dr.x);
            dr.y-=NINT(dr.y);
            dr.z-=NINT(dr.z);
            rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
            if(rr<1e-8) index=k;
          }
          if(index<0)
          {
            Framework[CurrentSystem].Atoms[CurrentFramework][nr_atoms].ReferencePosition=SpaceGroupElement[j];
            nr_atoms++;
          }
          if(nr_atoms>Framework[CurrentSystem].NumberOfUnitCellAtoms[CurrentFramework])
          {
            sg_possible=FALSE;
            break;
          }
        }
      }

      // if the asymmetric unit does not generate the proper amount of atoms then it is certainly the wrong space group
      if(nr_atoms!=Framework[CurrentSystem].NumberOfUnitCellAtoms[CurrentFramework]) sg_possible=FALSE;

      // choice the space group with the lowest amount of atoms and the highest symmetry (Note the '<=')
      if(sg_possible&&(Framework[CurrentSystem].NumberOfAsymmetricAtoms[CurrentFramework]<=LowestNumberOfAsymmetricUnits))
      {
        LowestNumberOfAsymmetricUnits=Framework[CurrentSystem].NumberOfAsymmetricAtoms[CurrentFramework];
        FoundSpaceGroup=sg;
        //WriteAsymmetricUnitCell(FoundSpaceGroup);
      }
      if(sg_possible) fprintf(stderr, "Possible space group: %d\n",FoundSpaceGroup);
    }

    // assume a certain space-group
    Framework[CurrentSystem].SpaceGroupIdentifier[CurrentFramework]=FoundSpaceGroup;
  }

  Framework[CurrentSystem].NumberOfAsymmetricAtoms[CurrentFramework]=0;
  for(k=0;k<Framework[CurrentSystem].NumberOfUnitCellAtoms[CurrentFramework];k++)
  {
    Type=Framework[CurrentSystem].Atoms[CurrentFramework][k].Type;
    pos=Framework[CurrentSystem].Atoms[CurrentFramework][k].Position;

    s=ConvertFromXYZtoABC(pos);

    s.x*=NumberOfUnitCells[CurrentSystem].x;
    s.y*=NumberOfUnitCells[CurrentSystem].y;
    s.z*=NumberOfUnitCells[CurrentSystem].z;

    FindAsymmetricAtom(Framework[CurrentSystem].SpaceGroupIdentifier[CurrentFramework],s,&t);

    index=-1;
    for(i=0;i<Framework[CurrentSystem].NumberOfAsymmetricAtoms[CurrentFramework];i++)
    {
      dr.x=Framework[CurrentSystem].AtomsAsymmetric[CurrentFramework][i].Position.x-t.x;
      dr.y=Framework[CurrentSystem].AtomsAsymmetric[CurrentFramework][i].Position.y-t.y;
      dr.z=Framework[CurrentSystem].AtomsAsymmetric[CurrentFramework][i].Position.z-t.z;
      rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
      if((rr<1e-8)&&(Type==Framework[CurrentSystem].AtomsAsymmetric[CurrentFramework][i].Type))
      {
        index=i;
        break;
      }
    }
    if(index<0)
    {
      Framework[CurrentSystem].AtomsAsymmetric[CurrentFramework][Framework[CurrentSystem].NumberOfAsymmetricAtoms[CurrentFramework]].Position=t;
      Framework[CurrentSystem].AtomsAsymmetric[CurrentFramework][Framework[CurrentSystem].NumberOfAsymmetricAtoms[CurrentFramework]].Type=Type;
      Framework[CurrentSystem].NumberOfAsymmetricAtoms[CurrentFramework]++;
    }
  }

  //WriteAsymmetricUnitCell(FoundSpaceGroup);
}


void WriteAsymmetricUnitCell(int sp)
{
  int i;
  REAL A,B,C;
  VECTOR pos;
  FILE *FilePtr;
  char buffer[256];

  A=(REAL)UnitCellSize[CurrentSystem].x;
  B=(REAL)UnitCellSize[CurrentSystem].y;
  C=(REAL)UnitCellSize[CurrentSystem].z;


  sprintf(buffer,"framework_asymmetric_%d%s.cssr",sp,FileNameAppend);

  FilePtr=fopen(buffer,"w");
  fprintf(FilePtr,"%38c%8.3lf%8.3lf%8.3lf\n",' ',
     (double)A,
     (double)B,
     (double)C);
  fprintf(FilePtr,"%21c%8.3lf%8.3lf%8.3lf%4cSPGR =  %s  OPT = 1\n",
     ' ',
     (double)AlphaAngle[CurrentSystem]*RAD2DEG,
     (double)BetaAngle[CurrentSystem]*RAD2DEG,
     (double)GammaAngle[CurrentSystem]*RAD2DEG,
     ' ',SpaceGroupData[sp].ShortInternationalHermannMauguinSpaceGroupSymbol);

  fprintf(FilePtr,"%4d%4d %s\n",Framework[CurrentSystem].NumberOfAsymmetricAtoms[CurrentFramework],0,"Created by Raspa-1.0");
  fprintf(FilePtr,"     0 %s         : %s\n",Framework[CurrentSystem].Name[CurrentFramework],
                Framework[CurrentSystem].Name[CurrentFramework]);
  for(i=0;i<Framework[CurrentSystem].NumberOfAsymmetricAtoms[CurrentFramework];i++)
  {
    pos=Framework[CurrentSystem].AtomsAsymmetric[CurrentFramework][i].Position;
    fprintf(FilePtr,"%4d %-4s  %9.5lf %9.5lf %9.5lf %4d%4d%4d%4d%4d%4d%4d%4d %7.3lf\n",
      i+1,
      PseudoAtoms[Framework[CurrentSystem].AtomsAsymmetric[CurrentFramework][i].Type].Name,
      (double)pos.x,
      (double)pos.y,
      (double)pos.z,
      0,0,0,0,0,0,0,0,
      (double)0.0);
  }
  fclose(FilePtr);
}



void FindSpatialGroup2(void)
{
  int i,j,k,PossibleSpaceGroup;
  int spaceGroupNumber,nAtoms,type;
  int *atomInList,nAtomsInList;
  int inverseAtomsFound,newAtomFound;

  VECTOR pos,dr,*positionAtomAsymmetric;
  int *typeAsymmetric,sizeLastSpaceGroup;
  REAL rr;

  pos.x=pos.y=pos.z=0.0;

  //Reserve enough memory for the asymmetric atoms
  nAtoms=Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];
  atomInList=(int*)calloc(nAtoms,sizeof(int));
  positionAtomAsymmetric=(VECTOR *)calloc(nAtoms,sizeof(VECTOR));
  typeAsymmetric=(int *)calloc(nAtoms,sizeof(int));

  //The space group 1 is always possible, so we will start with space group 2
  spaceGroupNumber=2;
  SpaceGroupSymmetry(spaceGroupNumber,pos);
  //Condition that will tell us when we have found an appropriate space group or not
  Framework[CurrentSystem].NumberOfAsymmetricAtoms[CurrentFramework]=nAtoms;
  Framework[CurrentSystem].SpaceGroupIdentifier[CurrentFramework]=1; //Default space group
  sizeLastSpaceGroup=1;

  //Convert everything to unitary units
  for(i=0;i<nAtoms;i++)
  {
    pos=Framework[CurrentSystem].Atoms[CurrentFramework][i].Position;
    pos=ConvertFromXYZtoABC(pos);
    Framework[CurrentSystem].Atoms[CurrentFramework][i].Position=pos;
  }

  while(SpaceGroupSize!=0)
  {
    PossibleSpaceGroup=TRUE;
    for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];i++)
    {
      atomInList[i]=FALSE;
    }
    nAtomsInList=0;
    nAtoms=0; //In this loop it will store the number of asymmetric atoms found

    //Apply the spatial group and check that the position generated belongs to the framework
    for(i=0;PossibleSpaceGroup&&(i<Framework[CurrentSystem].NumberOfAtoms[CurrentFramework]);i++)
    {
      //Check only the atoms that have not been generated with the inverse group yet
      if(!atomInList[i])
      {
        //Get the framework atom position in unitary units
        pos=Framework[CurrentSystem].Atoms[CurrentFramework][i].Position;
        type=Framework[CurrentSystem].Atoms[CurrentFramework][i].Type;

        //Get the symmetry positions and check against the original framework
        SpaceGroupSymmetry(spaceGroupNumber,pos);

        inverseAtomsFound=0;
        newAtomFound=FALSE;
        for(j=0;j<SpaceGroupSize;j++)
        {
          // apply boundary condition
          SpaceGroupElement[j].x-=NINT(SpaceGroupElement[j].x);
          SpaceGroupElement[j].y-=NINT(SpaceGroupElement[j].y);
          SpaceGroupElement[j].z-=NINT(SpaceGroupElement[j].z);

          if(SpaceGroupElement[j].x<0.0) SpaceGroupElement[j].x+=1.0;
          if(SpaceGroupElement[j].y<0.0) SpaceGroupElement[j].y+=1.0;
          if(SpaceGroupElement[j].z<0.0) SpaceGroupElement[j].z+=1.0;


          //Check the position against the original framework
          for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];k++)
          {
            //Check only the atoms of the same type
            if(type==Framework[CurrentSystem].Atoms[CurrentFramework][k].Type)
            {
              dr.x=Framework[CurrentSystem].Atoms[CurrentFramework][k].Position.x-SpaceGroupElement[j].x;
              dr.y=Framework[CurrentSystem].Atoms[CurrentFramework][k].Position.y-SpaceGroupElement[j].y;
              dr.z=Framework[CurrentSystem].Atoms[CurrentFramework][k].Position.z-SpaceGroupElement[j].z;
              dr.x-=NINT(dr.x);
              dr.y-=NINT(dr.y);
              dr.z-=NINT(dr.z);
              rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
              if(rr<1e-5)
              {
                inverseAtomsFound++;
                if(!atomInList[k])
                {
                  newAtomFound=TRUE;
                  atomInList[k]=TRUE;
                  nAtomsInList++;
                }
                break;
              }
            }
          }

        } //for SpaceGroupSize

        //If only one atom has been found it is due to the identity element and the group is discarded
        if(inverseAtomsFound==1)
        {
          PossibleSpaceGroup=FALSE;
          break;
        }
        //If this atom has generated new positions, include it in the list of asymmetric
        if(newAtomFound)
        {
          nAtoms++;
          positionAtomAsymmetric[nAtoms-1]=Framework[CurrentSystem].Atoms[CurrentFramework][i].Position;
          typeAsymmetric[nAtoms-1]=Framework[CurrentSystem].Atoms[CurrentFramework][i].Type;
        }
      } //if !atomInList
    } //for Atoms in framework

    //The generated symmetric structure must have the same number of atoms as the original framework
    //If the current space group is a candidate, check the best space group
    if(PossibleSpaceGroup&&nAtomsInList==Framework[CurrentSystem].NumberOfAtoms[CurrentFramework])
    {
      //If we successully generate the original framework with this space group, check the best successful space group
      //(larger size, then the asymmetric cell is smaller, and larger symmetry) and copy asymmetric cell to the original framework

      if(nAtoms<Framework[CurrentSystem].NumberOfAsymmetricAtoms[CurrentFramework]||
        ((nAtoms==Framework[CurrentSystem].NumberOfAsymmetricAtoms[CurrentFramework])&&(SpaceGroupSize>sizeLastSpaceGroup)))
      //if(nAtoms<Framework[CurrentSystem].NumberOfAsymmetricAtoms[CurrentFramework])
      {
        //Prepare the memory needed to store the asymmetric atoms
        if(Framework[CurrentSystem].NumberOfAsymmetricAtoms[CurrentFramework]==Framework[CurrentSystem].NumberOfAtoms[CurrentFramework])
        {
          Framework[CurrentSystem].AtomsAsymmetric[CurrentFramework]=(FRAMEWORK_ASYMMETRIC_ATOM*)calloc(nAtoms,sizeof(FRAMEWORK_ASYMMETRIC_ATOM));
        }
        else
        {
          Framework[CurrentSystem].AtomsAsymmetric[CurrentFramework]=(FRAMEWORK_ASYMMETRIC_ATOM*)realloc(
                                             Framework[CurrentSystem].AtomsAsymmetric[CurrentFramework],nAtoms*sizeof(FRAMEWORK_ASYMMETRIC_ATOM));
        }

        Framework[CurrentSystem].SpaceGroupIdentifier[CurrentFramework]=spaceGroupNumber;
        Framework[CurrentSystem].NumberOfAsymmetricAtoms[CurrentFramework]=nAtoms;

        fprintf(stderr, "Possible: %d %s %d nr. asym: %d %s\n",spaceGroupNumber,SpaceGroupData[spaceGroupNumber].ShortInternationalHermannMauguinSpaceGroupSymbol,
                         SpaceGroupData[spaceGroupNumber].Number,Framework[CurrentSystem].NumberOfAsymmetricAtoms[CurrentFramework],
                         SpaceGroupData[spaceGroupNumber].Standard?"standard":"non-standard");

        for(i=0;i<nAtoms;i++)
        {
          Framework[CurrentSystem].AtomsAsymmetric[CurrentFramework][i].Position=positionAtomAsymmetric[i];
          Framework[CurrentSystem].AtomsAsymmetric[CurrentFramework][i].Type=typeAsymmetric[i];
        }
        sizeLastSpaceGroup=SpaceGroupSize;

        WriteAsymmetricUnitCell(spaceGroupNumber);
      }
    }

    spaceGroupNumber++;
    SpaceGroupSymmetry(spaceGroupNumber,pos);

  }//while space group size != 0

  //Free memory that is not needed anymore
  free(atomInList);
  free(positionAtomAsymmetric);
  free(typeAsymmetric);

  fprintf(stderr, "Framework[CurrentSystem].NumberOfAsymmetricAtoms[CurrentFramework]: %d\n",Framework[CurrentSystem].NumberOfAsymmetricAtoms[CurrentFramework]);

  //If we have not found the space group, it is set to 1 by default
  if(Framework[CurrentSystem].NumberOfAsymmetricAtoms[CurrentFramework]==Framework[CurrentSystem].NumberOfAtoms[CurrentFramework])
  {
    nAtoms=Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];
    Framework[CurrentSystem].AtomsAsymmetric[CurrentFramework]=(FRAMEWORK_ASYMMETRIC_ATOM*)calloc(nAtoms,sizeof(FRAMEWORK_ASYMMETRIC_ATOM));
    for(i=0;i<nAtoms;i++)
    {
      Framework[CurrentSystem].AtomsAsymmetric[CurrentFramework][i].Position=Framework[CurrentSystem].Atoms[CurrentFramework][i].Position;
      Framework[CurrentSystem].AtomsAsymmetric[CurrentFramework][i].Type=Framework[CurrentSystem].Atoms[CurrentFramework][i].Type;
    }
  }

  //Convert back to non unitary units
  for(i=0;i<Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];i++)
  {
    pos=Framework[CurrentSystem].Atoms[CurrentFramework][i].Position;
    pos=ConvertFromABCtoXYZ(pos);
    Framework[CurrentSystem].Atoms[CurrentFramework][i].Position=pos;
  }
  for(i=0;i<Framework[CurrentSystem].NumberOfAsymmetricAtoms[CurrentFramework];i++)
  {
    pos=Framework[CurrentSystem].AtomsAsymmetric[CurrentFramework][i].Position;
    pos=ConvertFromABCtoXYZ(pos);
    Framework[CurrentSystem].AtomsAsymmetric[CurrentFramework][i].Position=pos;
  }
}

void PutNoiseOnFrameworkAtomicPositions(void)
{
  int f1,i;
  int system;

  for(system=0;system<NumberOfSystems;system++)
  {
    for(f1=0;f1<Framework[system].NumberOfFrameworks;f1++)
    {
      for(i=0;i<Framework[system].NumberOfAtoms[f1];i++)
      {
        Framework[system].Atoms[f1][i].Position.x+=0.01*(2.0*RandomNumber()-1.0);
        Framework[system].Atoms[f1][i].Position.y+=0.01*(2.0*RandomNumber()-1.0);
        Framework[system].Atoms[f1][i].Position.z+=0.01*(2.0*RandomNumber()-1.0);
      }
    }
  }
}

static int versionNumber=1;

void WriteRestartFramework(FILE *FilePtr)
{
  int i,j,k;
  REAL Check;

  fwrite(&versionNumber,sizeof(int),1,FilePtr);

  fwrite(&NumberOfSystems,sizeof(NumberOfSystems),1,FilePtr);
  fwrite(crystallographic_stats,sizeof(CRYSTALLOGRAPHIC_STATISTICS),NumberOfSystems,FilePtr);

  fwrite(Framework,sizeof(FRAMEWORK_COMPONENT),NumberOfSystems,FilePtr);
  fwrite(&CurrentFramework,sizeof(CurrentFramework),1,FilePtr);

  fwrite(UnitCellSize,sizeof(VECTOR),NumberOfSystems,FilePtr);
  fwrite(NumberOfUnitCells,sizeof(INT_VECTOR3),NumberOfSystems,FilePtr);
  fwrite(UnitCellBox,sizeof(REAL_MATRIX3x3),NumberOfSystems,FilePtr);
  fwrite(InverseUnitCellBox,sizeof(REAL_MATRIX3x3),NumberOfSystems,FilePtr);
  fwrite(Lowenstein,sizeof(int),NumberOfSystems,FilePtr);

  fwrite(&RemoveBondNeighboursFromLongRangeInteraction,sizeof(int),1,FilePtr);
  fwrite(&RemoveBendNeighboursFromLongRangeInteraction,sizeof(int),1,FilePtr);
  fwrite(&RemoveTorsionNeighboursFromLongRangeInteraction,sizeof(int),1,FilePtr);

  fwrite(&Remove12NeighboursFromVDWInteraction,sizeof(int),1,FilePtr);
  fwrite(&Remove13NeighboursFromVDWInteraction,sizeof(int),1,FilePtr);
  fwrite(&Remove14NeighboursFromVDWInteraction,sizeof(int),1,FilePtr);

  fwrite(&Remove12NeighboursFromChargeChargeInteraction,sizeof(int),1,FilePtr);
  fwrite(&Remove13NeighboursFromChargeChargeInteraction,sizeof(int),1,FilePtr);
  fwrite(&Remove14NeighboursFromChargeChargeInteraction,sizeof(int),1,FilePtr);

  fwrite(&Remove11NeighboursFromChargeBondDipoleInteraction,sizeof(int),1,FilePtr);
  fwrite(&Remove12NeighboursFromChargeBondDipoleInteraction,sizeof(int),1,FilePtr);
  fwrite(&Remove13NeighboursFromChargeBondDipoleInteraction,sizeof(int),1,FilePtr);
  fwrite(&Remove14NeighboursFromChargeBondDipoleInteraction,sizeof(int),1,FilePtr);

  fwrite(&Remove12NeighboursFromBondDipoleBondDipoleInteraction,sizeof(int),1,FilePtr);
  fwrite(&Remove13NeighboursFromBondDipoleBondDipoleInteraction,sizeof(int),1,FilePtr);
  fwrite(&Remove14NeighboursFromBondDipoleBondDipoleInteraction,sizeof(int),1,FilePtr);

  fwrite(&ImproperTorsionScanType,sizeof(int),1,FilePtr);

  fwrite(&InternalFrameworkLennardJonesInteractions,sizeof(int),1,FilePtr);

  for(i=0;i<NumberOfSystems;i++)
  {
    fwrite(&Framework[i].NumberOfAsymmetricIons,sizeof(int),1,FilePtr);

    if(Framework[i].NumberOfAsymmetricIons>0)
    {
      fwrite(&crystallographic_stats[i].Position,sizeof(VECTOR),Framework[i].NumberOfAsymmetricIons,FilePtr);
      fwrite(&crystallographic_stats[i].PositionSquared,sizeof(VECTOR),Framework[i].NumberOfAsymmetricIons,FilePtr);
      fwrite(&crystallographic_stats[i].Distance,sizeof(VECTOR),Framework[i].NumberOfAsymmetricIons,FilePtr);
      fwrite(&crystallographic_stats[i].Occupation,sizeof(REAL),Framework[i].NumberOfAsymmetricIons,FilePtr);
      fwrite(&crystallographic_stats[i].AverageDistance,sizeof(REAL),Framework[i].NumberOfAsymmetricIons,FilePtr);
      fwrite(&crystallographic_stats[i].RelativeOccupation,sizeof(REAL),Framework[i].NumberOfAsymmetricIons,FilePtr);
      fwrite(&crystallographic_stats[i].NumberOfCationSites,sizeof(int),Framework[i].NumberOfAsymmetricIons,FilePtr);
      fwrite(&crystallographic_stats[i].Count,sizeof(REAL),Framework[i].NumberOfAsymmetricIons,FilePtr);
      fwrite(&crystallographic_stats[i].TemperatureFactor,sizeof(REAL_MATRIX3x3),Framework[i].NumberOfAsymmetricIons,FilePtr);
    }
    fwrite(&Framework[i].NumberOfIons,sizeof(int),1,FilePtr);
    fwrite(&Framework[i].MaxNumberOfIons,sizeof(int),1,FilePtr);
    if(Framework[i].NumberOfIons>0)
      fwrite(Framework[i].Ions,sizeof(ATOM),Framework[i].NumberOfIons,FilePtr);

    fwrite(&Framework[i].NumberOfFrameworks,sizeof(int),1,FilePtr);

    if(Framework[i].NumberOfFrameworks==0)
    {
      fwrite(Framework[i].Name,sizeof(char[256]),1,FilePtr);
    }
    else
    {
      fwrite(Framework[i].NumberOfAtoms,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);
      fwrite(Framework[i].NumberOfFixedAtoms,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);
      fwrite(Framework[i].NumberOfFreeAtoms,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);
      fwrite(Framework[i].NumberOfCharges,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);

      fwrite(Framework[i].NumberOfUnitCellAtoms,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);
      fwrite(Framework[i].MaxNumberOfAtoms,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);
      fwrite(Framework[i].FrameworkProbability,sizeof(REAL),Framework[i].NumberOfFrameworks,FilePtr);
      for(j=0;j<Framework[i].NumberOfFrameworks;j++)
        fwrite(Framework[i].Atoms[j],sizeof(ATOM),Framework[i].NumberOfAtoms[j],FilePtr);

      fwrite(Framework[i].NumberOfAsymmetricAtoms,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);
      for(j=0;j<Framework[i].NumberOfFrameworks;j++)
        fwrite(Framework[i].AtomsAsymmetric[j],sizeof(FRAMEWORK_ASYMMETRIC_ATOM),Framework[i].NumberOfAsymmetricAtoms[j],FilePtr);

      fwrite(Framework[i].Name,sizeof(char[256]),Framework[i].NumberOfFrameworks,FilePtr);
      fwrite(Framework[i].FrameworkModels,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);
      fwrite(Framework[i].ShiftUnitCell,sizeof(VECTOR),Framework[i].NumberOfFrameworks,FilePtr);
      fwrite(Framework[i].Asymmetric,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);
      fwrite(Framework[i].SpaceGroupIdentifier,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);
      fwrite(Framework[i].CalculateSpaceGroup,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);
      fwrite(Framework[i].RemoveAtomNumberCodeFromLabel,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);
      fwrite(Framework[i].AddAtomNumberCodeToLabel,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);
      fwrite(Framework[i].InputFileType,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);

      fwrite(Framework[i].SurfaceAreas,sizeof(REAL),Framework[i].NumberOfFrameworks,FilePtr);
      fwrite(Framework[i].FrameworkDensityPerComponent,sizeof(REAL),Framework[i].NumberOfFrameworks,FilePtr);
      fwrite(Framework[i].FrameworkMassPerComponent,sizeof(REAL),Framework[i].NumberOfFrameworks,FilePtr);

      fwrite(Framework[i].NumberOfCitations,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);
      for(j=0;j<Framework[i].NumberOfFrameworks;j++)
        fwrite(Framework[i].CitationInformation[j],sizeof(CITATION_INFORMATION),Framework[i].NumberOfCitations[j],FilePtr);


      // write the connectity and neighbor-list of the frameworks
      for(j=0;j<Framework[i].NumberOfFrameworks;j++)
      {
        fwrite(Framework[i].Connectivity[j],sizeof(int),Framework[i].NumberOfAtoms[j],FilePtr);
        for(k=0;k<Framework[i].NumberOfAtoms[j];k++)
          fwrite(Framework[i].Neighbours[j][k],sizeof(int),Framework[i].Connectivity[j][k],FilePtr);
      }

      // write the core/shell information
      fwrite(&Framework[i].NumberOfCoreShellDefinitions,sizeof(int),1,FilePtr);
      fwrite(Framework[i].NumberOfCoreShells,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);
      if(Framework[i].NumberOfCoreShellDefinitions>0)
      {
        fwrite(Framework[i].CoreShellDefinitions,sizeof(PAIR),Framework[i].NumberOfCoreShellDefinitions,FilePtr);
        fwrite(Framework[i].NumberOfCoreShellsPerType,sizeof(int),Framework[i].NumberOfCoreShellDefinitions,FilePtr);

        for(j=0;j<Framework[i].NumberOfFrameworks;j++)
        {
          if(Framework[i].NumberOfCoreShells[j]>0)
            fwrite(Framework[i].CoreShellConnectivity[j],sizeof(int),Framework[i].NumberOfAtoms[j],FilePtr);
        }
      }

      fwrite(&Framework[i].NumberOfBondsDefinitions,sizeof(int),1,FilePtr);
      if(Framework[i].NumberOfBondsDefinitions>0)
      {
        fwrite(Framework[i].BondDefinitionType,sizeof(int),Framework[i].NumberOfBondsDefinitions,FilePtr);
        fwrite(Framework[i].BondDefinitions,sizeof(PAIR),Framework[i].NumberOfBondsDefinitions,FilePtr);
        fwrite(Framework[i].BondArgumentDefinitions,sizeof(REAL[MAX_BOND_POTENTIAL_ARGUMENTS]),Framework[i].NumberOfBondsDefinitions,FilePtr);
        fwrite(Framework[i].NumberOfBondsPerType,sizeof(int),Framework[i].NumberOfBondsDefinitions,FilePtr);

        fwrite(Framework[i].NumberOfBonds,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);
        fwrite(Framework[i].MaxNumberOfBonds,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);
        for(j=0;j<Framework[i].NumberOfFrameworks;j++)
        {
          if(Framework[i].NumberOfBonds[j]>0)
          {
            fwrite(Framework[i].Bonds[j],sizeof(PAIR),Framework[i].NumberOfBonds[j],FilePtr);
            fwrite(Framework[i].BondType[j],sizeof(int),Framework[i].NumberOfBonds[j],FilePtr);
            fwrite(Framework[i].BondArguments[j],sizeof(REAL[MAX_BOND_POTENTIAL_ARGUMENTS]),Framework[i].NumberOfBonds[j],FilePtr);
          }
        }
      }

      fwrite(&Framework[i].NumberOfBondDipoleDefinitions,sizeof(int),1,FilePtr);
      if(Framework[i].NumberOfBondDipoleDefinitions>0)
      {
        fwrite(Framework[i].NumberOfBondDipoles,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);
        fwrite(Framework[i].MaxNumberOfBondDipoles,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);

        fwrite(Framework[i].BondDipoleDefinitions,sizeof(PAIR),Framework[i].NumberOfBondDipoleDefinitions,FilePtr);
        fwrite(Framework[i].BondDipoleArgumentDefinition,sizeof(REAL),Framework[i].NumberOfBondDipoleDefinitions,FilePtr);
        fwrite(Framework[i].NumberOfBondDipolesPerType,sizeof(int),Framework[i].NumberOfBondDipoleDefinitions,FilePtr);

        // allocate first dimension of bond-dipoles
        for(j=0;j<Framework[i].NumberOfFrameworks;j++)
        {
          if(Framework[i].NumberOfBondDipoles[j]>0)
          {
            fwrite(Framework[i].BondDipoles[j],sizeof(PAIR),Framework[i].NumberOfBondDipoles[j],FilePtr);
            fwrite(Framework[i].BondDipoleMagnitude[j],sizeof(REAL),Framework[i].NumberOfBondDipoles[j],FilePtr);
          }
        }
      }

      fwrite(&Framework[i].NumberOfUreyBradleyDefinitions,sizeof(int),1,FilePtr);
      if(Framework[i].NumberOfUreyBradleyDefinitions>0)
      {
        fwrite(Framework[i].UreyBradleyDefinitionType,sizeof(int),Framework[i].NumberOfUreyBradleyDefinitions,FilePtr);
        fwrite(Framework[i].UreyBradleyDefinitions,sizeof(TRIPLE),Framework[i].NumberOfUreyBradleyDefinitions,FilePtr);
        fwrite(Framework[i].UreyBradleyArgumentDefinitions,sizeof(REAL[MAX_UREYBRADLEY_POTENTIAL_ARGUMENTS]),Framework[i].NumberOfUreyBradleyDefinitions,FilePtr);
        fwrite(Framework[i].NumberOfUreyBradleysPerType,sizeof(int),Framework[i].NumberOfUreyBradleyDefinitions,FilePtr);

        fwrite(Framework[i].NumberOfUreyBradleys,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);
        fwrite(Framework[i].MaxNumberOfUreyBradleys,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);

        for(j=0;j<Framework[i].NumberOfFrameworks;j++)
        {
          if(Framework[i].NumberOfUreyBradleys[j]>0)
          {
            fwrite(Framework[i].UreyBradleys[j],sizeof(TRIPLE),Framework[i].NumberOfUreyBradleys[j],FilePtr);
            fwrite(Framework[i].UreyBradleyType[j],sizeof(int),Framework[i].NumberOfUreyBradleys[j],FilePtr);
            fwrite(Framework[i].UreyBradleyArguments[j],sizeof(REAL[MAX_UREYBRADLEY_POTENTIAL_ARGUMENTS]),Framework[i].NumberOfUreyBradleys[j],FilePtr);
          }
        }
      }

      fwrite(&Framework[i].NumberOfBendDefinitions,sizeof(int),1,FilePtr);
      if(Framework[i].NumberOfBendDefinitions>0)
      {
        fwrite(Framework[i].BendDefinitionType,sizeof(int),Framework[i].NumberOfBendDefinitions,FilePtr);
        fwrite(Framework[i].BendDefinitions,sizeof(QUAD),Framework[i].NumberOfBendDefinitions,FilePtr);
        fwrite(Framework[i].BendArgumentDefinitions,sizeof(REAL[MAX_BEND_POTENTIAL_ARGUMENTS]),Framework[i].NumberOfBendDefinitions,FilePtr);
        fwrite(Framework[i].NumberOfBendsPerType,sizeof(int),Framework[i].NumberOfBendDefinitions,FilePtr);

        fwrite(Framework[i].NumberOfBends,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);
        fwrite(Framework[i].MaxNumberOfBends,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);

        for(j=0;j<Framework[i].NumberOfFrameworks;j++)
        {
          if(Framework[i].NumberOfBends[j]>0)
          {
            fwrite(Framework[i].Bends[j],sizeof(QUAD),Framework[i].NumberOfBends[j],FilePtr);
            fwrite(Framework[i].BendType[j],sizeof(int),Framework[i].NumberOfBends[j],FilePtr);
            fwrite(Framework[i].BendArguments[j],sizeof(REAL[MAX_BEND_POTENTIAL_ARGUMENTS]),Framework[i].NumberOfBends[j],FilePtr);
          }
        }
      }

      fwrite(&Framework[i].NumberOfInversionBendDefinitions,sizeof(int),1,FilePtr);
      if(Framework[i].NumberOfInversionBendDefinitions>0)
      {
        fwrite(Framework[i].InversionBendDefinitionType,sizeof(int),Framework[i].NumberOfInversionBendDefinitions,FilePtr);
        fwrite(Framework[i].InversionBendDefinitions,sizeof(QUAD),Framework[i].NumberOfInversionBendDefinitions,FilePtr);
        fwrite(Framework[i].InversionBendArgumentDefinitions,sizeof(REAL[MAX_INVERSION_BEND_POTENTIAL_ARGUMENTS]),Framework[i].NumberOfInversionBendDefinitions,FilePtr);
        fwrite(Framework[i].NumberOfInversionBendsPerType,sizeof(int),Framework[i].NumberOfInversionBendDefinitions,FilePtr);

        fwrite(Framework[i].NumberOfInversionBends,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);
        fwrite(Framework[i].MaxNumberOfInversionBends,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);

        for(j=0;j<Framework[i].NumberOfFrameworks;j++)
        {
          if(Framework[i].NumberOfInversionBends[j]>0)
          {
            fwrite(Framework[i].InversionBends[j],sizeof(QUAD),Framework[i].NumberOfInversionBends[j],FilePtr);
            fwrite(Framework[i].InversionBendType[j],sizeof(int),Framework[i].NumberOfInversionBends[j],FilePtr);
            fwrite(Framework[i].InversionBendArguments[j],sizeof(REAL[MAX_INVERSION_BEND_POTENTIAL_ARGUMENTS]),Framework[i].NumberOfInversionBends[j],FilePtr);
          }
        }
      }

      fwrite(&Framework[i].NumberOfTorsionDefinitions,sizeof(int),1,FilePtr);
      if(Framework[i].NumberOfTorsionDefinitions>0)
      {
        fwrite(Framework[i].TorsionDefinitionType,sizeof(int),Framework[i].NumberOfTorsionDefinitions,FilePtr);
        fwrite(Framework[i].TorsionDefinitions,sizeof(QUAD),Framework[i].NumberOfTorsionDefinitions,FilePtr);
        fwrite(Framework[i].TorsionArgumentDefinitions,sizeof(REAL[MAX_TORSION_POTENTIAL_ARGUMENTS]),Framework[i].NumberOfTorsionDefinitions,FilePtr);
        fwrite(Framework[i].NumberOfTorsionsPerType,sizeof(int),Framework[i].NumberOfTorsionDefinitions,FilePtr);

        fwrite(Framework[i].NumberOfTorsions,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);
        fwrite(Framework[i].MaxNumberOfTorsions,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);

        for(j=0;j<Framework[i].NumberOfFrameworks;j++)
        {
          if(Framework[i].NumberOfTorsions[j]>0)
          {
            fwrite(Framework[i].Torsions[j],sizeof(QUAD),Framework[i].NumberOfTorsions[j],FilePtr);
            fwrite(Framework[i].TorsionType[j],sizeof(int),Framework[i].NumberOfTorsions[j],FilePtr);
            fwrite(Framework[i].TorsionArguments[j],sizeof(REAL[MAX_TORSION_POTENTIAL_ARGUMENTS]),Framework[i].NumberOfTorsions[j],FilePtr);
          }
        }
      }

      fwrite(&Framework[i].NumberOfImproperTorsionDefinitions,sizeof(int),1,FilePtr);
      if(Framework[i].NumberOfImproperTorsionDefinitions>0)
      {
        fwrite(Framework[i].ImproperTorsionDefinitionType,sizeof(int),Framework[i].NumberOfImproperTorsionDefinitions,FilePtr);
        fwrite(Framework[i].ImproperTorsionDefinitions,sizeof(QUAD),Framework[i].NumberOfImproperTorsionDefinitions,FilePtr);
        fwrite(Framework[i].ImproperTorsionArgumentDefinitions,sizeof(REAL[MAX_IMPROPER_TORSION_POTENTIAL_ARGUMENTS]),Framework[i].NumberOfImproperTorsionDefinitions,FilePtr);
        fwrite(Framework[i].NumberOfImproperTorsionsPerType,sizeof(int),Framework[i].NumberOfImproperTorsionDefinitions,FilePtr);

        fwrite(Framework[i].NumberOfImproperTorsions,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);
        fwrite(Framework[i].MaxNumberOfImproperTorsions,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);

        for(j=0;j<Framework[i].NumberOfFrameworks;j++)
        {
          if(Framework[i].NumberOfImproperTorsions[j]>0)
          {
            fwrite(Framework[i].ImproperTorsions[j],sizeof(QUAD),Framework[i].NumberOfImproperTorsions[j],FilePtr);
            fwrite(Framework[i].ImproperTorsionType[j],sizeof(int),Framework[i].NumberOfImproperTorsions[j],FilePtr);
            fwrite(Framework[i].ImproperTorsionArguments[j],sizeof(REAL[MAX_IMPROPER_TORSION_POTENTIAL_ARGUMENTS]),Framework[i].NumberOfImproperTorsions[j],FilePtr);
          }
        }
      }

      fwrite(&Framework[i].NumberOfOutOfPlaneDefinitions,sizeof(int),1,FilePtr);
      if(Framework[i].NumberOfOutOfPlaneDefinitions>0)
      {
        fwrite(Framework[i].OutOfPlaneDefinitionType,sizeof(int),Framework[i].NumberOfOutOfPlaneDefinitions,FilePtr);
        fwrite(Framework[i].OutOfPlaneDefinitions,sizeof(QUAD),Framework[i].NumberOfOutOfPlaneDefinitions,FilePtr);
        fwrite(Framework[i].OutOfPlaneArgumentDefinitions,sizeof(REAL[MAX_OUT_OF_PLANE_POTENTIAL_ARGUMENTS]),Framework[i].NumberOfOutOfPlaneDefinitions,FilePtr);
        fwrite(Framework[i].NumberOfOutOfPlanesPerType,sizeof(int),Framework[i].NumberOfOutOfPlaneDefinitions,FilePtr);

        fwrite(Framework[i].NumberOfOutOfPlanes,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);
        fwrite(Framework[i].MaxNumberOfOutOfPlanes,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);

        for(j=0;j<Framework[i].NumberOfFrameworks;j++)
        {
          if(Framework[i].NumberOfOutOfPlanes[j]>0)
          {
            fwrite(Framework[i].OutOfPlanes[j],sizeof(QUAD),Framework[i].NumberOfOutOfPlanes[j],FilePtr);
            fwrite(Framework[i].OutOfPlaneType[j],sizeof(int),Framework[i].NumberOfOutOfPlanes[j],FilePtr);
            fwrite(Framework[i].OutOfPlaneArguments[j],sizeof(REAL[MAX_OUT_OF_PLANE_POTENTIAL_ARGUMENTS]),Framework[i].NumberOfOutOfPlanes[j],FilePtr);
          }
        }
      }

      fwrite(&Framework[i].NumberOfBondBondDefinitions,sizeof(int),1,FilePtr);
      if(Framework[i].NumberOfBondBondDefinitions>0)
      {
        fwrite(Framework[i].BondBondDefinitionType,sizeof(int),Framework[i].NumberOfBondBondDefinitions,FilePtr);
        fwrite(Framework[i].BondBondDefinitions,sizeof(TRIPLE),Framework[i].NumberOfBondBondDefinitions,FilePtr);
        fwrite(Framework[i].BondBondArgumentDefinitions,sizeof(REAL[MAX_BOND_BOND_POTENTIAL_ARGUMENTS]),Framework[i].NumberOfBondBondDefinitions,FilePtr);
        fwrite(Framework[i].NumberOfBondBondsPerType,sizeof(int),Framework[i].NumberOfBondBondDefinitions,FilePtr);

        fwrite(Framework[i].NumberOfBondBonds,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);
        fwrite(Framework[i].MaxNumberOfBondBonds,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);

        for(j=0;j<Framework[i].NumberOfFrameworks;j++)
        {
          if(Framework[i].NumberOfBondBonds[j]>0)
          {
            fwrite(Framework[i].BondBonds[j],sizeof(TRIPLE),Framework[i].NumberOfBondBonds[j],FilePtr);
            fwrite(Framework[i].BondBondType[j],sizeof(int),Framework[i].NumberOfBondBonds[j],FilePtr);
            fwrite(Framework[i].BondBondArguments[j],sizeof(REAL[MAX_BOND_BOND_POTENTIAL_ARGUMENTS]),Framework[i].NumberOfBondBonds[j],FilePtr);
          }
        }
      }

      fwrite(&Framework[i].NumberOfBondBendDefinitions,sizeof(int),1,FilePtr);
      if(Framework[i].NumberOfBondBendDefinitions>0)
      {
        fwrite(Framework[i].BondBendDefinitionType,sizeof(int),Framework[i].NumberOfBondBendDefinitions,FilePtr);
        fwrite(Framework[i].BondBendDefinitions,sizeof(TRIPLE),Framework[i].NumberOfBondBendDefinitions,FilePtr);
        fwrite(Framework[i].BondBendArgumentDefinitions,sizeof(REAL[MAX_BOND_BEND_POTENTIAL_ARGUMENTS]),Framework[i].NumberOfBondBendDefinitions,FilePtr);
        fwrite(Framework[i].NumberOfBondBendsPerType,sizeof(int),Framework[i].NumberOfBondBendDefinitions,FilePtr);

        fwrite(Framework[i].NumberOfBondBends,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);
        fwrite(Framework[i].MaxNumberOfBondBends,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);

        for(j=0;j<Framework[i].NumberOfFrameworks;j++)
        {
          if(Framework[i].NumberOfBondBends[j]>0)
          {
            fwrite(Framework[i].BondBends[j],sizeof(TRIPLE),Framework[i].NumberOfBondBends[j],FilePtr);
            fwrite(Framework[i].BondBendType[j],sizeof(int),Framework[i].NumberOfBondBends[j],FilePtr);
            fwrite(Framework[i].BondBendArguments[j],sizeof(REAL[MAX_BOND_BEND_POTENTIAL_ARGUMENTS]),Framework[i].NumberOfBondBends[j],FilePtr);
          }
        }
      }

      fwrite(&Framework[i].NumberOfBendBendDefinitions,sizeof(int),1,FilePtr);
      if(Framework[i].NumberOfBendBendDefinitions>0)
      {
        fwrite(Framework[i].BendBendDefinitionType,sizeof(int),Framework[i].NumberOfBendBendDefinitions,FilePtr);
        fwrite(Framework[i].BendBendDefinitions,sizeof(QUAD),Framework[i].NumberOfBendBendDefinitions,FilePtr);
        fwrite(Framework[i].BendBendArgumentDefinitions,sizeof(REAL[MAX_BEND_BEND_POTENTIAL_ARGUMENTS]),Framework[i].NumberOfBendBendDefinitions,FilePtr);
        fwrite(Framework[i].NumberOfBendBendsPerType,sizeof(int),Framework[i].NumberOfBendBendDefinitions,FilePtr);

        fwrite(Framework[i].NumberOfBendBends,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);
        fwrite(Framework[i].MaxNumberOfBendBends,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);

        for(j=0;j<Framework[i].NumberOfFrameworks;j++)
        {
          if(Framework[i].NumberOfBendBends[j]>0)
          {
            fwrite(Framework[i].BendBends[j],sizeof(QUAD),Framework[i].NumberOfBendBends[j],FilePtr);
            fwrite(Framework[i].BendBendType[j],sizeof(int),Framework[i].NumberOfBendBends[j],FilePtr);
            fwrite(Framework[i].BendBendArguments[j],sizeof(REAL[MAX_BEND_BEND_POTENTIAL_ARGUMENTS]),Framework[i].NumberOfBendBends[j],FilePtr);
          }
        }
      }

      fwrite(&Framework[i].NumberOfBondTorsionDefinitions,sizeof(int),1,FilePtr);
      if(Framework[i].NumberOfBondTorsionDefinitions>0)
      {
        fwrite(Framework[i].BondTorsionDefinitionType,sizeof(int),Framework[i].NumberOfBondTorsionDefinitions,FilePtr);
        fwrite(Framework[i].BondTorsionDefinitions,sizeof(QUAD),Framework[i].NumberOfBondTorsionDefinitions,FilePtr);
        fwrite(Framework[i].BondTorsionArgumentDefinitions,sizeof(REAL[MAX_BOND_TORSION_POTENTIAL_ARGUMENTS]),Framework[i].NumberOfBondTorsionDefinitions,FilePtr);
        fwrite(Framework[i].NumberOfBondTorsionsPerType,sizeof(int),Framework[i].NumberOfBondTorsionDefinitions,FilePtr);

        fwrite(Framework[i].NumberOfBondTorsions,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);
        fwrite(Framework[i].MaxNumberOfBondTorsions,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);

        for(j=0;j<Framework[i].NumberOfFrameworks;j++)
        {
          if(Framework[i].NumberOfBondTorsions[j]>0)
          {
            fwrite(Framework[i].BondTorsions[j],sizeof(QUAD),Framework[i].NumberOfBondTorsions[j],FilePtr);
            fwrite(Framework[i].BondTorsionType[j],sizeof(int),Framework[i].NumberOfBondTorsions[j],FilePtr);
            fwrite(Framework[i].BondTorsionArguments[j],sizeof(REAL[MAX_BOND_TORSION_POTENTIAL_ARGUMENTS]),Framework[i].NumberOfBondTorsions[j],FilePtr);
          }
        }
      }

      fwrite(&Framework[i].NumberOfBendTorsionDefinitions,sizeof(int),1,FilePtr);
      if(Framework[i].NumberOfBendTorsionDefinitions>0)
      {
        fwrite(Framework[i].BendTorsionDefinitionType,sizeof(int),Framework[i].NumberOfBendTorsionDefinitions,FilePtr);
        fwrite(Framework[i].BendTorsionDefinitions,sizeof(QUAD),Framework[i].NumberOfBendTorsionDefinitions,FilePtr);
        fwrite(Framework[i].BendTorsionArgumentDefinitions,sizeof(REAL[MAX_BEND_TORSION_POTENTIAL_ARGUMENTS]),Framework[i].NumberOfBendTorsionDefinitions,FilePtr);
        fwrite(Framework[i].NumberOfBendTorsionsPerType,sizeof(int),Framework[i].NumberOfBendTorsionDefinitions,FilePtr);

        fwrite(Framework[i].NumberOfBendTorsions,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);
        fwrite(Framework[i].MaxNumberOfBendTorsions,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);

        for(j=0;j<Framework[i].NumberOfFrameworks;j++)
        {
          if(Framework[i].NumberOfBendTorsions[j]>0)
          {
            fwrite(Framework[i].BendTorsions[j],sizeof(QUAD),Framework[i].NumberOfBendTorsions[j],FilePtr);
            fwrite(Framework[i].BendTorsionType[j],sizeof(int),Framework[i].NumberOfBendTorsions[j],FilePtr);
            fwrite(Framework[i].BendTorsionArguments[j],sizeof(REAL[MAX_BEND_TORSION_POTENTIAL_ARGUMENTS]),Framework[i].NumberOfBendTorsions[j],FilePtr);
          }
        }
      }

      fwrite(Framework[i].NumberOfIntraVDW,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);
      fwrite(Framework[i].NumberOfExcludedIntraVDW,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);

      fwrite(Framework[i].NumberOfIntraCharges,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);
      fwrite(Framework[i].NumberOfIntraChargeBondDipoles,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);
      fwrite(Framework[i].NumberOfIntraBondDipoles,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);

      fwrite(Framework[i].NumberOfExcludedIntraChargeCharge,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);
      fwrite(Framework[i].MaxNumberOfExcludedIntraChargeCharge,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);

      fwrite(Framework[i].NumberOfExcludedIntraChargeBondDipole,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);
      fwrite(Framework[i].MaxNumberOfExcludedIntraChargeBondDipole,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);

      fwrite(Framework[i].NumberOfExcludedIntraBondDipoleBondDipole,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);
      fwrite(Framework[i].MaxNumberOfExcludedIntraBondDipoleBondDipole,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);

      if(Framework[i].FrameworkModel==FLEXIBLE)
      {
        for(j=0;j<Framework[i].NumberOfFrameworks;j++)
        {
          if(Framework[i].FrameworkModels[j]==FLEXIBLE)
          {
            if(Framework[i].NumberOfExcludedIntraChargeCharge[j]>0)
              fwrite(Framework[i].ExcludedIntraChargeCharge[j],sizeof(PAIR),Framework[i].NumberOfExcludedIntraChargeCharge[j],FilePtr);

            if(Framework[i].NumberOfExcludedIntraChargeBondDipole[j])
              fwrite(Framework[i].ExcludedIntraChargeBondDipole[j],sizeof(PAIR),Framework[i].NumberOfExcludedIntraChargeBondDipole[j],FilePtr);

            if(Framework[i].NumberOfExcludedIntraBondDipoleBondDipole[j]>0)
                fwrite(Framework[i].ExcludedIntraBondDipoleBondDipole[j],sizeof(PAIR),Framework[i].NumberOfExcludedIntraBondDipoleBondDipole[j],FilePtr);
          }
        }
      }
    }
  }

  fwrite(&NumberOfSingleSubstitutionRules,sizeof(int),1,FilePtr);
  if(NumberOfSingleSubstitutionRules>0)
  {
    fwrite(SubstitutionSingleFrameworkAtom,sizeof(char[3][256]),NumberOfSingleSubstitutionRules,FilePtr);
    fwrite(SubstitutionSingleFrameworkAtomTypes,sizeof(int[3]),NumberOfSingleSubstitutionRules,FilePtr);
  }

  fwrite(&NumberOfSubstitutionRules,sizeof(int),1,FilePtr);
  fwrite(&NumberOfRandomSubstitutions,sizeof(int),1,FilePtr);
  if(NumberOfSubstitutionRules>0)
  {
    fwrite(SubstitutionFrameworkAtoms,sizeof(char[3][256]),NumberOfSubstitutionRules,FilePtr);
    fwrite(SubstitutionFrameworkAtomTypes,sizeof(int[3]),NumberOfSubstitutionRules,FilePtr);
  }

  fwrite(&NumberOfSubstitutions,sizeof(int),1,FilePtr);
  if(NumberOfSubstitutions>0)
    fwrite(ListOfAtomSubstitutions,sizeof(int[3]),NumberOfSubstitutions,FilePtr);

  // write the framework atom modification rules
  fwrite(&NumberOfModificationRules,sizeof(int),1,FilePtr);
  if(NumberOfModificationRules>0)
  {
    fwrite(ModificationRuleType,sizeof(int),NumberOfModificationRules,FilePtr);
    fwrite(ModifyFrameworkAtoms,sizeof(char[10][256]),NumberOfModificationRules,FilePtr);
    fwrite(ModifyFrameworkAtomTypes,sizeof(int[10]),NumberOfModificationRules,FilePtr);
  }
  // write the framework atom forbidden connectivity rules
  fwrite(&NumberOfForbiddenConnectivityRules,sizeof(int),1,FilePtr);
  if(NumberOfForbiddenConnectivityRules>0)
  {
    fwrite(ForbiddenConnectivityAtoms,sizeof(char[3][256]),NumberOfForbiddenConnectivityRules,FilePtr);
    fwrite(ForbiddenConnectivityTypes,sizeof(int[3]),NumberOfForbiddenConnectivityRules,FilePtr);
  }

  Check=123456789.0;
  fwrite(&Check,1,sizeof(REAL),FilePtr);
}


void ReadRestartFramework(FILE *FilePtr)
{
  int i,j,k;
  REAL Check;
  int readversionNumber=0;

  fread(&readversionNumber,sizeof(int),1,FilePtr);
  if(readversionNumber > versionNumber)
  {
    fprintf(stderr,"Upgrade to last version of RASPA to read binary restart-file");
    exit(-1);
  }

  fread(&NumberOfSystems,sizeof(NumberOfSystems),1,FilePtr);
  Framework=(FRAMEWORK_COMPONENT*)calloc(NumberOfSystems,sizeof(FRAMEWORK_COMPONENT));

  crystallographic_stats=(CRYSTALLOGRAPHIC_STATISTICS*)calloc(NumberOfSystems,sizeof(CRYSTALLOGRAPHIC_STATISTICS));
  fread(crystallographic_stats,sizeof(CRYSTALLOGRAPHIC_STATISTICS),NumberOfSystems,FilePtr);

  fread(Framework,NumberOfSystems,sizeof(FRAMEWORK_COMPONENT),FilePtr);
  fread(&CurrentFramework,sizeof(CurrentFramework),1,FilePtr);

  UnitCellSize=(VECTOR*)calloc(NumberOfSystems,sizeof(VECTOR));
  NumberOfUnitCells=(INT_VECTOR3*)calloc(NumberOfSystems,sizeof(INT_VECTOR3));
  UnitCellBox=(REAL_MATRIX3x3*)calloc(NumberOfSystems,sizeof(REAL_MATRIX3x3));
  InverseUnitCellBox=(REAL_MATRIX3x3*)calloc(NumberOfSystems,sizeof(REAL_MATRIX3x3));
  Lowenstein=(int*)calloc(NumberOfSystems,sizeof(int));

  fread(UnitCellSize,sizeof(VECTOR),NumberOfSystems,FilePtr);
  fread(NumberOfUnitCells,sizeof(INT_VECTOR3),NumberOfSystems,FilePtr);
  fread(UnitCellBox,sizeof(REAL_MATRIX3x3),NumberOfSystems,FilePtr);
  fread(InverseUnitCellBox,sizeof(REAL_MATRIX3x3),NumberOfSystems,FilePtr);
  fread(Lowenstein,sizeof(int),NumberOfSystems,FilePtr);

  fread(&RemoveBondNeighboursFromLongRangeInteraction,sizeof(int),1,FilePtr);
  fread(&RemoveBendNeighboursFromLongRangeInteraction,sizeof(int),1,FilePtr);
  fread(&RemoveTorsionNeighboursFromLongRangeInteraction,sizeof(int),1,FilePtr);

  fread(&Remove12NeighboursFromVDWInteraction,sizeof(int),1,FilePtr);
  fread(&Remove13NeighboursFromVDWInteraction,sizeof(int),1,FilePtr);
  fread(&Remove14NeighboursFromVDWInteraction,sizeof(int),1,FilePtr);

  fread(&Remove12NeighboursFromChargeChargeInteraction,sizeof(int),1,FilePtr);
  fread(&Remove13NeighboursFromChargeChargeInteraction,sizeof(int),1,FilePtr);
  fread(&Remove14NeighboursFromChargeChargeInteraction,sizeof(int),1,FilePtr);

  fread(&Remove11NeighboursFromChargeBondDipoleInteraction,sizeof(int),1,FilePtr);
  fread(&Remove12NeighboursFromChargeBondDipoleInteraction,sizeof(int),1,FilePtr);
  fread(&Remove13NeighboursFromChargeBondDipoleInteraction,sizeof(int),1,FilePtr);
  fread(&Remove14NeighboursFromChargeBondDipoleInteraction,sizeof(int),1,FilePtr);

  fread(&Remove12NeighboursFromBondDipoleBondDipoleInteraction,sizeof(int),1,FilePtr);
  fread(&Remove13NeighboursFromBondDipoleBondDipoleInteraction,sizeof(int),1,FilePtr);
  fread(&Remove14NeighboursFromBondDipoleBondDipoleInteraction,sizeof(int),1,FilePtr);

  fread(&ImproperTorsionScanType,sizeof(int),1,FilePtr);

  fread(&InternalFrameworkLennardJonesInteractions,sizeof(int),1,FilePtr);

  for(i=0;i<NumberOfSystems;i++)
  {
    fread(&Framework[i].NumberOfAsymmetricIons,sizeof(int),1,FilePtr);

    crystallographic_stats[i].Position=(VECTOR*)calloc(Framework[i].NumberOfAsymmetricIons,sizeof(VECTOR));
    crystallographic_stats[i].PositionSquared=(VECTOR*)calloc(Framework[i].NumberOfAsymmetricIons,sizeof(VECTOR));
    crystallographic_stats[i].Distance=(VECTOR*)calloc(Framework[i].NumberOfAsymmetricIons,sizeof(VECTOR));
    crystallographic_stats[i].Occupation=(REAL*)calloc(Framework[i].NumberOfAsymmetricIons,sizeof(REAL));
    crystallographic_stats[i].AverageDistance=(REAL*)calloc(Framework[i].NumberOfAsymmetricIons,sizeof(REAL));
    crystallographic_stats[i].RelativeOccupation=(REAL*)calloc(Framework[i].NumberOfAsymmetricIons,sizeof(REAL));
    crystallographic_stats[i].NumberOfCationSites=(int*)calloc(Framework[i].NumberOfAsymmetricIons,sizeof(int));
    crystallographic_stats[i].Count=(REAL*)calloc(Framework[i].NumberOfAsymmetricIons,sizeof(REAL));
    crystallographic_stats[i].TemperatureFactor=(REAL_MATRIX3x3*)calloc(Framework[i].NumberOfAsymmetricIons,sizeof(REAL_MATRIX3x3));

    if(Framework[i].NumberOfAsymmetricIons>0)
    {
      fread(&crystallographic_stats[i].Position,sizeof(VECTOR),Framework[i].NumberOfAsymmetricIons,FilePtr);
      fread(&crystallographic_stats[i].PositionSquared,sizeof(VECTOR),Framework[i].NumberOfAsymmetricIons,FilePtr);
      fread(&crystallographic_stats[i].Distance,sizeof(VECTOR),Framework[i].NumberOfAsymmetricIons,FilePtr);
      fread(&crystallographic_stats[i].Occupation,sizeof(REAL),Framework[i].NumberOfAsymmetricIons,FilePtr);
      fread(&crystallographic_stats[i].AverageDistance,sizeof(REAL),Framework[i].NumberOfAsymmetricIons,FilePtr);
      fread(&crystallographic_stats[i].RelativeOccupation,sizeof(REAL),Framework[i].NumberOfAsymmetricIons,FilePtr);
      fread(&crystallographic_stats[i].NumberOfCationSites,sizeof(int),Framework[i].NumberOfAsymmetricIons,FilePtr);
      fread(&crystallographic_stats[i].Count,sizeof(REAL),Framework[i].NumberOfAsymmetricIons,FilePtr);
      fread(&crystallographic_stats[i].TemperatureFactor,sizeof(REAL_MATRIX3x3),Framework[i].NumberOfAsymmetricIons,FilePtr);
    }

    fread(&Framework[i].NumberOfIons,sizeof(int),1,FilePtr);
    fread(&Framework[i].MaxNumberOfIons,sizeof(int),1,FilePtr);
    Framework[i].MaxNumberOfIons=Framework[i].NumberOfIons;
    if(Framework[i].NumberOfIons>0)
    {
      Framework[i].Ions=(ATOM*)calloc(Framework[i].NumberOfIons,sizeof(ATOM));
      fread(Framework[i].Ions,Framework[i].NumberOfIons,sizeof(ATOM),FilePtr);
    }

    fread(&Framework[i].NumberOfFrameworks,sizeof(int),1,FilePtr);

    if(Framework[i].NumberOfFrameworks==0)
    {
      Framework[i].NumberOfAtoms=(int*)calloc(1,sizeof(int));
      Framework[i].NumberOfAtoms[0]=0;
      Framework[i].NumberOfFixedAtoms=(int*)calloc(1,sizeof(int));
      Framework[i].NumberOfFixedAtoms[0]=0;
      Framework[i].NumberOfFreeAtoms=(int*)calloc(1,sizeof(int));
      Framework[i].NumberOfFreeAtoms[0]=0;
      Framework[i].NumberOfCharges=(int*)calloc(1,sizeof(int));
      Framework[i].NumberOfCharges[0]=0;
      Framework[i].NumberOfUnitCellAtoms=(int*)calloc(1,sizeof(int));
      Framework[i].NumberOfUnitCellAtoms[0]=0;
      Framework[i].NumberOfAsymmetricAtoms=(int*)calloc(1,sizeof(int));
      Framework[i].NumberOfAsymmetricAtoms[0]=0;
      Framework[i].Name=(char(*)[256])calloc(1,sizeof(char[256]));
      fread(Framework[i].Name,sizeof(char[256]),1,FilePtr);
    }
    else
    {
      Framework[i].Atoms=(ATOM**)calloc(Framework[i].NumberOfFrameworks,sizeof(ATOM*));
      Framework[i].NumberOfAtoms=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      Framework[i].NumberOfFixedAtoms=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      Framework[i].NumberOfFreeAtoms=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      Framework[i].NumberOfCharges=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      Framework[i].NumberOfUnitCellAtoms=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      Framework[i].MaxNumberOfAtoms=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      Framework[i].FrameworkProbability=(REAL*)calloc(Framework[i].NumberOfFrameworks,sizeof(REAL));
      fread(Framework[i].NumberOfAtoms,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);
      fread(Framework[i].NumberOfFixedAtoms,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);
      fread(Framework[i].NumberOfFreeAtoms,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);
      fread(Framework[i].NumberOfCharges,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);
      fread(Framework[i].NumberOfUnitCellAtoms,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);
      fread(Framework[i].MaxNumberOfAtoms,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);
      fread(Framework[i].FrameworkProbability,sizeof(REAL),Framework[i].NumberOfFrameworks,FilePtr);

      for(j=0;j<Framework[i].NumberOfFrameworks;j++)
      {
        Framework[i].MaxNumberOfAtoms=Framework[i].NumberOfAtoms;
        Framework[i].Atoms[j]=(ATOM*)calloc(Framework[i].NumberOfAtoms[j],sizeof(ATOM));
        fread(Framework[i].Atoms[j],sizeof(ATOM),Framework[i].NumberOfAtoms[j],FilePtr);
      }

      Framework[i].AtomsAsymmetric=(FRAMEWORK_ASYMMETRIC_ATOM**)calloc(Framework[i].NumberOfFrameworks,sizeof(FRAMEWORK_ASYMMETRIC_ATOM*));
      Framework[i].NumberOfAsymmetricAtoms=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      fread(Framework[i].NumberOfAsymmetricAtoms,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);

      for(j=0;j<Framework[i].NumberOfFrameworks;j++)
      {
        Framework[i].AtomsAsymmetric[j]=(FRAMEWORK_ASYMMETRIC_ATOM*)calloc(Framework[i].NumberOfAsymmetricAtoms[j],sizeof(FRAMEWORK_ASYMMETRIC_ATOM));
        fread(Framework[i].AtomsAsymmetric[j],sizeof(FRAMEWORK_ASYMMETRIC_ATOM),Framework[i].NumberOfAsymmetricAtoms[j],FilePtr);
      }


      Framework[i].Name=(char(*)[256])calloc(Framework[i].NumberOfFrameworks,sizeof(char[256]));
      Framework[i].FrameworkModels=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      Framework[i].ShiftUnitCell=(VECTOR*)calloc(Framework[i].NumberOfFrameworks,sizeof(VECTOR));
      Framework[i].Asymmetric=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      Framework[i].SpaceGroupIdentifier=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      Framework[i].CalculateSpaceGroup=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      Framework[i].RemoveAtomNumberCodeFromLabel=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      Framework[i].AddAtomNumberCodeToLabel=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      Framework[i].InputFileType=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));

      Framework[i].SurfaceAreas=(REAL*)calloc(Framework[i].NumberOfFrameworks,sizeof(REAL));
      Framework[i].FrameworkDensityPerComponent=(REAL*)calloc(Framework[i].NumberOfFrameworks,sizeof(REAL));
      Framework[i].FrameworkMassPerComponent=(REAL*)calloc(Framework[i].NumberOfFrameworks,sizeof(REAL));

      Framework[i].NumberOfCitations=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      Framework[i].CitationInformation=(CITATION_INFORMATION**)calloc(Framework[i].NumberOfFrameworks,sizeof(CITATION_INFORMATION*));

      fread(Framework[i].Name,sizeof(char[256]),Framework[i].NumberOfFrameworks,FilePtr);
      fread(Framework[i].FrameworkModels,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);
      fread(Framework[i].ShiftUnitCell,sizeof(VECTOR),Framework[i].NumberOfFrameworks,FilePtr);
      fread(Framework[i].Asymmetric,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);
      fread(Framework[i].SpaceGroupIdentifier,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);
      fread(Framework[i].CalculateSpaceGroup,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);
      fread(Framework[i].RemoveAtomNumberCodeFromLabel,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);
      fread(Framework[i].AddAtomNumberCodeToLabel,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);
      fread(Framework[i].InputFileType,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);

      fread(Framework[i].SurfaceAreas,sizeof(REAL),Framework[i].NumberOfFrameworks,FilePtr);
      fread(Framework[i].FrameworkDensityPerComponent,sizeof(REAL),Framework[i].NumberOfFrameworks,FilePtr);
      fread(Framework[i].FrameworkMassPerComponent,sizeof(REAL),Framework[i].NumberOfFrameworks,FilePtr);

      fread(Framework[i].NumberOfCitations,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);
      for(j=0;j<Framework[i].NumberOfFrameworks;j++)
      {
        Framework[i].CitationInformation[j]=(CITATION_INFORMATION*)calloc(Framework[i].NumberOfCitations[j],sizeof(CITATION_INFORMATION));
        fread(Framework[i].CitationInformation[j],sizeof(CITATION_INFORMATION),Framework[i].NumberOfCitations[j],FilePtr);
      }

      // read the connectity and neighbor-list of the frameworks
      Framework[i].Connectivity=(int**)calloc(Framework[i].NumberOfFrameworks,sizeof(int*));
      Framework[i].Neighbours=(int***)calloc(Framework[i].NumberOfFrameworks,sizeof(int**));
      for(j=0;j<Framework[i].NumberOfFrameworks;j++)
      {
        Framework[i].Connectivity[j]=(int*)calloc(Framework[i].NumberOfAtoms[j],sizeof(int));
        fread(Framework[i].Connectivity[j],sizeof(int),Framework[i].NumberOfAtoms[j],FilePtr);

        Framework[i].Neighbours[j]=(int**)calloc(Framework[i].NumberOfAtoms[j],sizeof(int*));
        for(k=0;k<Framework[i].NumberOfAtoms[j];k++)
        {
          Framework[i].Neighbours[j][k]=(int*)calloc(Framework[i].Connectivity[j][k],sizeof(int));
          fread(Framework[i].Neighbours[j][k],sizeof(int),Framework[i].Connectivity[j][k],FilePtr);
        }
      }

      // read the core/shell information
      fread(&Framework[i].NumberOfCoreShellDefinitions,1,sizeof(int),FilePtr);
      Framework[i].NumberOfCoreShells=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      fread(Framework[i].NumberOfCoreShells,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);
      if(Framework[i].NumberOfCoreShellDefinitions>0)
      {
        Framework[i].CoreShellDefinitions=(PAIR*)calloc(Framework[i].NumberOfCoreShellDefinitions,sizeof(PAIR));
        Framework[i].NumberOfCoreShellsPerType=(int*)calloc(Framework[i].NumberOfCoreShellDefinitions,sizeof(int));
        Framework[i].CoreShellConnectivity=(int**)calloc(Framework[i].NumberOfFrameworks,sizeof(int*));

        fread(Framework[i].CoreShellDefinitions,sizeof(PAIR),Framework[i].NumberOfCoreShellDefinitions,FilePtr);
        fread(Framework[i].NumberOfCoreShellsPerType,sizeof(int),Framework[i].NumberOfCoreShellDefinitions,FilePtr);


        for(j=0;j<Framework[i].NumberOfFrameworks;j++)
        {
          if(Framework[i].NumberOfCoreShells[j]>0)
          {
            Framework[i].CoreShellConnectivity[j]=(int*)calloc(Framework[i].NumberOfAtoms[j],sizeof(int));
            fread(Framework[i].CoreShellConnectivity[j],sizeof(int),Framework[i].NumberOfAtoms[j],FilePtr);
          }
        }
      }


      fread(&Framework[i].NumberOfBondsDefinitions,sizeof(int),1,FilePtr);
      Framework[i].NumberOfBonds=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      if(Framework[i].NumberOfBondsDefinitions>0)
      {
        // alloc memory
        Framework[i].BondDefinitionType=(int*)calloc(Framework[i].NumberOfBondsDefinitions,sizeof(int));
        Framework[i].BondDefinitions=(PAIR*)calloc(Framework[i].NumberOfBondsDefinitions,sizeof(PAIR));
        Framework[i].BondArgumentDefinitions=(REAL(*)[MAX_BOND_POTENTIAL_ARGUMENTS])
          calloc(Framework[i].NumberOfBondsDefinitions,sizeof(REAL[MAX_BOND_POTENTIAL_ARGUMENTS]));
        Framework[i].NumberOfBondsPerType=(int*)calloc(Framework[i].NumberOfBondsDefinitions,sizeof(int));

        fread(Framework[i].BondDefinitionType,sizeof(int),Framework[i].NumberOfBondsDefinitions,FilePtr);
        fread(Framework[i].BondDefinitions,sizeof(PAIR),Framework[i].NumberOfBondsDefinitions,FilePtr);
        fread(Framework[i].BondArgumentDefinitions,sizeof(REAL[MAX_BOND_POTENTIAL_ARGUMENTS]),Framework[i].NumberOfBondsDefinitions,FilePtr);
        fread(Framework[i].NumberOfBondsPerType,sizeof(int),Framework[i].NumberOfBondsDefinitions,FilePtr);

        Framework[i].MaxNumberOfBonds=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
        Framework[i].Bonds=(PAIR**)calloc(Framework[i].NumberOfFrameworks,sizeof(PAIR*));
        Framework[i].BondType=(int**)calloc(Framework[i].NumberOfFrameworks,sizeof(int*));
        Framework[i].BondArguments=(REAL(**)[MAX_BOND_POTENTIAL_ARGUMENTS])calloc(Framework[i].NumberOfFrameworks,
                                  sizeof(REAL(*)[MAX_BOND_POTENTIAL_ARGUMENTS]));

        fread(Framework[i].NumberOfBonds,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);
        fread(Framework[i].MaxNumberOfBonds,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);

        for(j=0;j<Framework[i].NumberOfFrameworks;j++)
        {
          // the first run, the MaxNumberOfBonds is larger than NumberOfBonds, but once the amount of bonds is computed
          // it will never change, so set the maximum number of bonds to the current number
          Framework[i].MaxNumberOfBonds[j]=Framework[i].NumberOfBonds[j];

          if(Framework[i].NumberOfBonds[j]>0)
          {
            Framework[i].Bonds[j]=(PAIR*)calloc(Framework[i].NumberOfBonds[j],sizeof(PAIR));
            Framework[i].BondType[j]=(int*)calloc(Framework[i].NumberOfBonds[j],sizeof(int));
            Framework[i].BondArguments[j]=(REAL(*)[MAX_BOND_POTENTIAL_ARGUMENTS])
                    calloc(Framework[i].NumberOfBonds[j],sizeof(REAL[MAX_BOND_POTENTIAL_ARGUMENTS]));

            fread(Framework[i].Bonds[j],sizeof(PAIR),Framework[i].NumberOfBonds[j],FilePtr);
            fread(Framework[i].BondType[j],sizeof(int),Framework[i].NumberOfBonds[j],FilePtr);
            fread(Framework[i].BondArguments[j],sizeof(REAL[MAX_BOND_POTENTIAL_ARGUMENTS]),Framework[i].NumberOfBonds[j],FilePtr);
          }
        }
      }

      fread(&Framework[i].NumberOfBondDipoleDefinitions,sizeof(int),1,FilePtr);
      Framework[i].NumberOfBondDipoles=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      if(Framework[i].NumberOfBondDipoleDefinitions>0)
      {
        // allocate first dimension of bond-dipoles
        Framework[i].MaxNumberOfBondDipoles=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
        Framework[i].BondDipoles=(PAIR**)calloc(Framework[i].NumberOfFrameworks,sizeof(PAIR*));
        Framework[i].BondDipoleMagnitude=(REAL**)calloc(Framework[i].NumberOfFrameworks,sizeof(REAL*));

        fread(Framework[i].NumberOfBondDipoles,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);
        fread(Framework[i].MaxNumberOfBondDipoles,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);

        Framework[i].BondDipoleDefinitions=(PAIR*)calloc(Framework[i].NumberOfBondDipoleDefinitions,sizeof(PAIR));
        Framework[i].BondDipoleArgumentDefinition=(REAL*)calloc(Framework[i].NumberOfBondDipoleDefinitions,sizeof(REAL));
        Framework[i].NumberOfBondDipolesPerType=(int*)calloc(Framework[i].NumberOfBondDipoleDefinitions,sizeof(int));

        fread(Framework[i].BondDipoleDefinitions,sizeof(PAIR),Framework[i].NumberOfBondDipoleDefinitions,FilePtr);
        fread(Framework[i].BondDipoleArgumentDefinition,sizeof(REAL),Framework[i].NumberOfBondDipoleDefinitions,FilePtr);
        fread(Framework[i].NumberOfBondDipolesPerType,sizeof(int),Framework[i].NumberOfBondDipoleDefinitions,FilePtr);

        // allocate first dimension of bond-dipoles
        for(j=0;j<Framework[i].NumberOfFrameworks;j++)
        {
          // the first run, the MaxNumberOfBonds is larger than NumberOfBonds, but once the amount of bonds is computed
          // it will never change, so set the maximum number of bonds to the current number
          Framework[i].MaxNumberOfBondDipoles[j]=Framework[i].NumberOfBondDipoles[j];

          if(Framework[i].NumberOfBondDipoles[j]>0)
          {
            Framework[i].BondDipoles[j]=(PAIR*)calloc(Framework[i].NumberOfBondDipoles[j],sizeof(PAIR));
            Framework[i].BondDipoleMagnitude[j]=(REAL*)calloc(Framework[i].NumberOfBondDipoles[j],sizeof(REAL));

            fread(Framework[i].BondDipoles[j],sizeof(PAIR),Framework[i].NumberOfBondDipoles[j],FilePtr);
            fread(Framework[i].BondDipoleMagnitude[j],sizeof(REAL),Framework[i].NumberOfBondDipoles[j],FilePtr);
          }
        }
      }

      fread(&Framework[i].NumberOfUreyBradleyDefinitions,sizeof(int),1,FilePtr);
      Framework[i].NumberOfUreyBradleys=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      if(Framework[i].NumberOfUreyBradleyDefinitions>0)
      {
        // alloc memory
        Framework[i].UreyBradleyDefinitionType=(int*)calloc(Framework[i].NumberOfUreyBradleyDefinitions,sizeof(int));
        Framework[i].UreyBradleyDefinitions=(TRIPLE*)calloc(Framework[i].NumberOfUreyBradleyDefinitions,sizeof(TRIPLE));
        Framework[i].UreyBradleyArgumentDefinitions=(REAL(*)[MAX_UREYBRADLEY_POTENTIAL_ARGUMENTS])
          calloc(Framework[i].NumberOfUreyBradleyDefinitions,sizeof(REAL[MAX_UREYBRADLEY_POTENTIAL_ARGUMENTS]));
        Framework[i].NumberOfUreyBradleysPerType=(int*)calloc(Framework[i].NumberOfUreyBradleyDefinitions,sizeof(int));

        fread(Framework[i].UreyBradleyDefinitionType,sizeof(int),Framework[i].NumberOfUreyBradleyDefinitions,FilePtr);
        fread(Framework[i].UreyBradleyDefinitions,sizeof(TRIPLE),Framework[i].NumberOfUreyBradleyDefinitions,FilePtr);
        fread(Framework[i].UreyBradleyArgumentDefinitions,sizeof(REAL[MAX_UREYBRADLEY_POTENTIAL_ARGUMENTS]),Framework[i].NumberOfUreyBradleyDefinitions,FilePtr);
        fread(Framework[i].NumberOfUreyBradleysPerType,sizeof(int),Framework[i].NumberOfUreyBradleyDefinitions,FilePtr);

        Framework[i].MaxNumberOfUreyBradleys=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
        Framework[i].UreyBradleys=(TRIPLE**)calloc(Framework[i].NumberOfFrameworks,sizeof(TRIPLE*));
        Framework[i].UreyBradleyType=(int**)calloc(Framework[i].NumberOfFrameworks,sizeof(int*));
        Framework[i].UreyBradleyArguments=(REAL(**)[MAX_UREYBRADLEY_POTENTIAL_ARGUMENTS])calloc(Framework[i].NumberOfFrameworks,
                                    sizeof(REAL(*)[MAX_UREYBRADLEY_POTENTIAL_ARGUMENTS]));

        fread(Framework[i].NumberOfUreyBradleys,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);
        fread(Framework[i].MaxNumberOfUreyBradleys,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);

        for(j=0;j<Framework[i].NumberOfFrameworks;j++)
        {
          // the first run, the MaxNumberOfUreyBradleys is larger than NumberOfUreyBradley, but once the amount of bonds is computed
          // it will never change, so set the maximum number of bonds to the current number
          Framework[i].MaxNumberOfUreyBradleys[j]=Framework[i].NumberOfUreyBradleys[j];

          if(Framework[i].NumberOfUreyBradleys[j]>0)
          {
            Framework[i].UreyBradleys[j]=(TRIPLE*)calloc(Framework[i].NumberOfUreyBradleys[j],sizeof(TRIPLE));
            Framework[i].UreyBradleyType[j]=(int*)calloc(Framework[i].NumberOfUreyBradleys[j],sizeof(int));
            Framework[i].UreyBradleyArguments[j]=(REAL(*)[MAX_UREYBRADLEY_POTENTIAL_ARGUMENTS])
                    calloc(Framework[i].NumberOfUreyBradleys[j],sizeof(REAL[MAX_UREYBRADLEY_POTENTIAL_ARGUMENTS]));

            fread(Framework[i].UreyBradleys[j],sizeof(TRIPLE),Framework[i].NumberOfUreyBradleys[j],FilePtr);
            fread(Framework[i].UreyBradleyType[j],sizeof(int),Framework[i].NumberOfUreyBradleys[j],FilePtr);
            fread(Framework[i].UreyBradleyArguments[j],sizeof(REAL[MAX_UREYBRADLEY_POTENTIAL_ARGUMENTS]),Framework[i].NumberOfUreyBradleys[j],FilePtr);
          }
        }
      }

      fread(&Framework[i].NumberOfBendDefinitions,sizeof(int),1,FilePtr);
      Framework[i].NumberOfBends=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      if(Framework[i].NumberOfBendDefinitions>0)
      {
        // alloc memory
        Framework[i].BendDefinitionType=(int*)calloc(Framework[i].NumberOfBendDefinitions,sizeof(int));
        Framework[i].BendDefinitions=(QUAD*)calloc(Framework[i].NumberOfBendDefinitions,sizeof(QUAD));
        Framework[i].BendArgumentDefinitions=(REAL(*)[MAX_BEND_POTENTIAL_ARGUMENTS])
          calloc(Framework[i].NumberOfBendDefinitions,sizeof(REAL[MAX_BEND_POTENTIAL_ARGUMENTS]));
        Framework[i].NumberOfBendsPerType=(int*)calloc(Framework[i].NumberOfBendDefinitions,sizeof(int));

        fread(Framework[i].BendDefinitionType,sizeof(int),Framework[i].NumberOfBendDefinitions,FilePtr);
        fread(Framework[i].BendDefinitions,sizeof(QUAD),Framework[i].NumberOfBendDefinitions,FilePtr);
        fread(Framework[i].BendArgumentDefinitions,sizeof(REAL[MAX_BEND_POTENTIAL_ARGUMENTS]),Framework[i].NumberOfBendDefinitions,FilePtr);
        fread(Framework[i].NumberOfBendsPerType,sizeof(int),Framework[i].NumberOfBendDefinitions,FilePtr);

        Framework[i].MaxNumberOfBends=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
        Framework[i].Bends=(QUAD**)calloc(Framework[i].NumberOfFrameworks,sizeof(QUAD*));
        Framework[i].BendType=(int**)calloc(Framework[i].NumberOfFrameworks,sizeof(int*));
        Framework[i].BendArguments=(REAL(**)[MAX_BEND_POTENTIAL_ARGUMENTS])calloc(Framework[i].NumberOfFrameworks,
                                    sizeof(REAL(*)[MAX_BEND_POTENTIAL_ARGUMENTS]));

        fread(Framework[i].NumberOfBends,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);
        fread(Framework[i].MaxNumberOfBends,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);

        for(j=0;j<Framework[i].NumberOfFrameworks;j++)
        {
          // the first run, the MaxNumberOfBends is larger than NumberOfBends, but once the amount of bonds is computed
          // it will never change, so set the maximum number of bonds to the current number
          Framework[i].MaxNumberOfBends[j]=Framework[i].NumberOfBends[j];

          if(Framework[i].NumberOfBends[j]>0)
          {
            Framework[i].Bends[j]=(QUAD*)calloc(Framework[i].NumberOfBends[j],sizeof(QUAD));
            Framework[i].BendType[j]=(int*)calloc(Framework[i].NumberOfBends[j],sizeof(int));
            Framework[i].BendArguments[j]=(REAL(*)[MAX_BEND_POTENTIAL_ARGUMENTS])
                    calloc(Framework[i].NumberOfBends[j],sizeof(REAL[MAX_BEND_POTENTIAL_ARGUMENTS]));

            fread(Framework[i].Bends[j],sizeof(QUAD),Framework[i].NumberOfBends[j],FilePtr);
            fread(Framework[i].BendType[j],sizeof(int),Framework[i].NumberOfBends[j],FilePtr);
            fread(Framework[i].BendArguments[j],sizeof(REAL[MAX_BEND_POTENTIAL_ARGUMENTS]),Framework[i].NumberOfBends[j],FilePtr);
          }
        }
      }

      fread(&Framework[i].NumberOfInversionBendDefinitions,sizeof(int),1,FilePtr);
      Framework[i].NumberOfInversionBends=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      if(Framework[i].NumberOfInversionBendDefinitions>0)
      {
        // alloc memory
        Framework[i].InversionBendDefinitionType=(int*)calloc(Framework[i].NumberOfInversionBendDefinitions,sizeof(int));
        Framework[i].InversionBendDefinitions=(QUAD*)calloc(Framework[i].NumberOfInversionBendDefinitions,sizeof(QUAD));
        Framework[i].InversionBendArgumentDefinitions=(REAL(*)[MAX_INVERSION_BEND_POTENTIAL_ARGUMENTS])
          calloc(Framework[i].NumberOfInversionBendDefinitions,sizeof(REAL[MAX_INVERSION_BEND_POTENTIAL_ARGUMENTS]));
        Framework[i].NumberOfInversionBendsPerType=(int*)calloc(Framework[i].NumberOfInversionBendDefinitions,sizeof(int));

        fread(Framework[i].InversionBendDefinitionType,sizeof(int),Framework[i].NumberOfInversionBendDefinitions,FilePtr);
        fread(Framework[i].InversionBendDefinitions,sizeof(QUAD),Framework[i].NumberOfInversionBendDefinitions,FilePtr);
        fread(Framework[i].InversionBendArgumentDefinitions,sizeof(REAL[MAX_INVERSION_BEND_POTENTIAL_ARGUMENTS]),Framework[i].NumberOfInversionBendDefinitions,FilePtr);
        fread(Framework[i].NumberOfInversionBendsPerType,sizeof(int),Framework[i].NumberOfInversionBendDefinitions,FilePtr);

        Framework[i].MaxNumberOfInversionBends=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
        Framework[i].InversionBends=(QUAD**)calloc(Framework[i].NumberOfFrameworks,sizeof(QUAD*));
        Framework[i].InversionBendType=(int**)calloc(Framework[i].NumberOfFrameworks,sizeof(int*));
        Framework[i].InversionBendArguments=(REAL(**)[MAX_INVERSION_BEND_POTENTIAL_ARGUMENTS])calloc(Framework[i].NumberOfFrameworks,
                                    sizeof(REAL(*)[MAX_INVERSION_BEND_POTENTIAL_ARGUMENTS]));

        fread(Framework[i].NumberOfInversionBends,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);
        fread(Framework[i].MaxNumberOfInversionBends,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);

        for(j=0;j<Framework[i].NumberOfFrameworks;j++)
        {
          // the first run, the MaxNumberOfBends is larger than NumberOfBends, but once the amount of bonds is computed
          // it will never change, so set the maximum number of bonds to the current number
          Framework[i].MaxNumberOfInversionBends[j]=Framework[i].NumberOfInversionBends[j];

          if(Framework[i].NumberOfInversionBends[j]>0)
          {
            Framework[i].InversionBends[j]=(QUAD*)calloc(Framework[i].NumberOfInversionBends[j],sizeof(QUAD));
            Framework[i].InversionBendType[j]=(int*)calloc(Framework[i].NumberOfInversionBends[j],sizeof(int));
            Framework[i].InversionBendArguments[j]=(REAL(*)[MAX_INVERSION_BEND_POTENTIAL_ARGUMENTS])
                    calloc(Framework[i].NumberOfInversionBends[j],sizeof(REAL[MAX_INVERSION_BEND_POTENTIAL_ARGUMENTS]));

            fread(Framework[i].InversionBends[j],sizeof(QUAD),Framework[i].NumberOfInversionBends[j],FilePtr);
            fread(Framework[i].InversionBendType[j],sizeof(int),Framework[i].NumberOfInversionBends[j],FilePtr);
            fread(Framework[i].InversionBendArguments[j],sizeof(REAL[MAX_INVERSION_BEND_POTENTIAL_ARGUMENTS]),Framework[i].NumberOfInversionBends[j],FilePtr);
          }
        }
      }

      fread(&Framework[i].NumberOfTorsionDefinitions,sizeof(int),1,FilePtr);
      Framework[i].NumberOfTorsions=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      if(Framework[i].NumberOfTorsionDefinitions>0)
      {
        // alloc memory
        Framework[i].TorsionDefinitionType=(int*)calloc(Framework[i].NumberOfTorsionDefinitions,sizeof(int));
        Framework[i].TorsionDefinitions=(QUAD*)calloc(Framework[i].NumberOfTorsionDefinitions,sizeof(QUAD));
        Framework[i].TorsionArgumentDefinitions=(REAL(*)[MAX_TORSION_POTENTIAL_ARGUMENTS])
          calloc(Framework[i].NumberOfTorsionDefinitions,sizeof(REAL[MAX_TORSION_POTENTIAL_ARGUMENTS]));
        Framework[i].NumberOfTorsionsPerType=(int*)calloc(Framework[i].NumberOfTorsionDefinitions,sizeof(int));

        fread(Framework[i].TorsionDefinitionType,sizeof(int),Framework[i].NumberOfTorsionDefinitions,FilePtr);
        fread(Framework[i].TorsionDefinitions,sizeof(QUAD),Framework[i].NumberOfTorsionDefinitions,FilePtr);
        fread(Framework[i].TorsionArgumentDefinitions,sizeof(REAL[MAX_TORSION_POTENTIAL_ARGUMENTS]),Framework[i].NumberOfTorsionDefinitions,FilePtr);
        fread(Framework[i].NumberOfTorsionsPerType,sizeof(int),Framework[i].NumberOfTorsionDefinitions,FilePtr);

        Framework[i].MaxNumberOfTorsions=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
        Framework[i].Torsions=(QUAD**)calloc(Framework[i].NumberOfFrameworks,sizeof(QUAD*));
        Framework[i].TorsionType=(int**)calloc(Framework[i].NumberOfFrameworks,sizeof(int*));
        Framework[i].TorsionArguments=(REAL(**)[MAX_TORSION_POTENTIAL_ARGUMENTS])calloc(Framework[i].NumberOfFrameworks,
                                    sizeof(REAL(*)[MAX_TORSION_POTENTIAL_ARGUMENTS]));

        fread(Framework[i].NumberOfTorsions,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);
        fread(Framework[i].MaxNumberOfTorsions,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);

        for(j=0;j<Framework[i].NumberOfFrameworks;j++)
        {
          // the first run, the MaxNumberOfBends is larger than NumberOfBends, but once the amount of bonds is computed
          // it will never change, so set the maximum number of bonds to the current number
          Framework[i].MaxNumberOfTorsions[j]=Framework[i].NumberOfTorsions[j];

          if(Framework[i].NumberOfTorsions[j]>0)
          {
            Framework[i].Torsions[j]=(QUAD*)calloc(Framework[i].NumberOfTorsions[j],sizeof(QUAD));
            Framework[i].TorsionType[j]=(int*)calloc(Framework[i].NumberOfTorsions[j],sizeof(int));
            Framework[i].TorsionArguments[j]=(REAL(*)[MAX_TORSION_POTENTIAL_ARGUMENTS])
                    calloc(Framework[i].NumberOfTorsions[j],sizeof(REAL[MAX_TORSION_POTENTIAL_ARGUMENTS]));

            fread(Framework[i].Torsions[j],sizeof(QUAD),Framework[i].NumberOfTorsions[j],FilePtr);
            fread(Framework[i].TorsionType[j],sizeof(int),Framework[i].NumberOfTorsions[j],FilePtr);
            fread(Framework[i].TorsionArguments[j],sizeof(REAL[MAX_TORSION_POTENTIAL_ARGUMENTS]),Framework[i].NumberOfTorsions[j],FilePtr);
          }
        }
      }

      fread(&Framework[i].NumberOfImproperTorsionDefinitions,sizeof(int),1,FilePtr);
      Framework[i].NumberOfImproperTorsions=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      if(Framework[i].NumberOfImproperTorsionDefinitions>0)
      {
        // alloc memory
        Framework[i].ImproperTorsionDefinitionType=(int*)calloc(Framework[i].NumberOfImproperTorsionDefinitions,sizeof(int));
        Framework[i].ImproperTorsionDefinitions=(QUAD*)calloc(Framework[i].NumberOfImproperTorsionDefinitions,sizeof(QUAD));
        Framework[i].ImproperTorsionArgumentDefinitions=(REAL(*)[MAX_IMPROPER_TORSION_POTENTIAL_ARGUMENTS])
          calloc(Framework[i].NumberOfImproperTorsionDefinitions,sizeof(REAL[MAX_IMPROPER_TORSION_POTENTIAL_ARGUMENTS]));
        Framework[i].NumberOfImproperTorsionsPerType=(int*)calloc(Framework[i].NumberOfImproperTorsionDefinitions,sizeof(int));

        fread(Framework[i].ImproperTorsionDefinitionType,sizeof(int),Framework[i].NumberOfImproperTorsionDefinitions,FilePtr);
        fread(Framework[i].ImproperTorsionDefinitions,sizeof(QUAD),Framework[i].NumberOfImproperTorsionDefinitions,FilePtr);
        fread(Framework[i].ImproperTorsionArgumentDefinitions,sizeof(REAL[MAX_IMPROPER_TORSION_POTENTIAL_ARGUMENTS]),Framework[i].NumberOfImproperTorsionDefinitions,FilePtr);
        fread(Framework[i].NumberOfImproperTorsionsPerType,sizeof(int),Framework[i].NumberOfImproperTorsionDefinitions,FilePtr);

        Framework[i].MaxNumberOfImproperTorsions=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
        Framework[i].ImproperTorsions=(QUAD**)calloc(Framework[i].NumberOfFrameworks,sizeof(QUAD*));
        Framework[i].ImproperTorsionType=(int**)calloc(Framework[i].NumberOfFrameworks,sizeof(int*));
        Framework[i].ImproperTorsionArguments=(REAL(**)[MAX_IMPROPER_TORSION_POTENTIAL_ARGUMENTS])calloc(Framework[i].NumberOfFrameworks,
                                    sizeof(REAL(*)[MAX_IMPROPER_TORSION_POTENTIAL_ARGUMENTS]));

        fread(Framework[i].NumberOfImproperTorsions,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);
        fread(Framework[i].MaxNumberOfImproperTorsions,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);

        for(j=0;j<Framework[i].NumberOfFrameworks;j++)
        {
          // the first run, the MaxNumberOfBends is larger than NumberOfBends, but once the amount of bonds is computed
          // it will never change, so set the maximum number of bonds to the current number
          Framework[i].MaxNumberOfImproperTorsions[j]=Framework[i].NumberOfImproperTorsions[j];

          if(Framework[i].NumberOfImproperTorsions[j]>0)
          {
            Framework[i].ImproperTorsions[j]=(QUAD*)calloc(Framework[i].NumberOfImproperTorsions[j],sizeof(QUAD));
            Framework[i].ImproperTorsionType[j]=(int*)calloc(Framework[i].NumberOfImproperTorsions[j],sizeof(int));
            Framework[i].ImproperTorsionArguments[j]=(REAL(*)[MAX_IMPROPER_TORSION_POTENTIAL_ARGUMENTS])
                    calloc(Framework[i].MaxNumberOfImproperTorsions[j],sizeof(REAL[MAX_IMPROPER_TORSION_POTENTIAL_ARGUMENTS]));

            fread(Framework[i].ImproperTorsions[j],sizeof(QUAD),Framework[i].NumberOfImproperTorsions[j],FilePtr);
            fread(Framework[i].ImproperTorsionType[j],sizeof(int),Framework[i].NumberOfImproperTorsions[j],FilePtr);
            fread(Framework[i].ImproperTorsionArguments[j],sizeof(REAL[MAX_IMPROPER_TORSION_POTENTIAL_ARGUMENTS]),Framework[i].NumberOfImproperTorsions[j],FilePtr);
          }
        }
      }

      fread(&Framework[i].NumberOfOutOfPlaneDefinitions,sizeof(int),1,FilePtr);
      Framework[i].NumberOfOutOfPlanes=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      if(Framework[i].NumberOfOutOfPlaneDefinitions>0)
      {
        // alloc memory
        Framework[i].OutOfPlaneDefinitionType=(int*)calloc(Framework[i].NumberOfOutOfPlaneDefinitions,sizeof(int));
        Framework[i].OutOfPlaneDefinitions=(QUAD*)calloc(Framework[i].NumberOfOutOfPlaneDefinitions,sizeof(QUAD));
        Framework[i].OutOfPlaneArgumentDefinitions=(REAL(*)[MAX_OUT_OF_PLANE_POTENTIAL_ARGUMENTS])
          calloc(Framework[i].NumberOfOutOfPlaneDefinitions,sizeof(REAL[MAX_OUT_OF_PLANE_POTENTIAL_ARGUMENTS]));
        Framework[i].NumberOfOutOfPlanesPerType=(int*)calloc(Framework[i].NumberOfOutOfPlaneDefinitions,sizeof(int));

        fread(Framework[i].OutOfPlaneDefinitionType,sizeof(int),Framework[i].NumberOfOutOfPlaneDefinitions,FilePtr);
        fread(Framework[i].OutOfPlaneDefinitions,sizeof(QUAD),Framework[i].NumberOfOutOfPlaneDefinitions,FilePtr);
        fread(Framework[i].OutOfPlaneArgumentDefinitions,sizeof(REAL[MAX_OUT_OF_PLANE_POTENTIAL_ARGUMENTS]),Framework[i].NumberOfOutOfPlaneDefinitions,FilePtr);
        fread(Framework[i].NumberOfOutOfPlanesPerType,sizeof(int),Framework[i].NumberOfOutOfPlaneDefinitions,FilePtr);

        Framework[i].MaxNumberOfOutOfPlanes=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
        Framework[i].OutOfPlanes=(QUAD**)calloc(Framework[i].NumberOfFrameworks,sizeof(QUAD*));
        Framework[i].OutOfPlaneType=(int**)calloc(Framework[i].NumberOfFrameworks,sizeof(int*));
        Framework[i].OutOfPlaneArguments=(REAL(**)[MAX_OUT_OF_PLANE_POTENTIAL_ARGUMENTS])calloc(Framework[i].NumberOfFrameworks,
                                    sizeof(REAL(*)[MAX_OUT_OF_PLANE_POTENTIAL_ARGUMENTS]));

        fread(Framework[i].NumberOfOutOfPlanes,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);
        fread(Framework[i].MaxNumberOfOutOfPlanes,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);

        for(j=0;j<Framework[i].NumberOfFrameworks;j++)
        {
          Framework[i].MaxNumberOfOutOfPlanes[j]=Framework[i].NumberOfOutOfPlanes[j];

          if(Framework[i].NumberOfOutOfPlanes[j]>0)
          {
            Framework[i].OutOfPlanes[j]=(QUAD*)calloc(Framework[i].NumberOfOutOfPlanes[j],sizeof(QUAD));
            Framework[i].OutOfPlaneType[j]=(int*)calloc(Framework[i].NumberOfOutOfPlanes[j],sizeof(int));
            Framework[i].OutOfPlaneArguments[j]=(REAL(*)[MAX_OUT_OF_PLANE_POTENTIAL_ARGUMENTS])
                    calloc(Framework[i].MaxNumberOfOutOfPlanes[j],sizeof(REAL[MAX_OUT_OF_PLANE_POTENTIAL_ARGUMENTS]));

            fread(Framework[i].OutOfPlanes[j],sizeof(QUAD),Framework[i].NumberOfOutOfPlanes[j],FilePtr);
            fread(Framework[i].OutOfPlaneType[j],sizeof(int),Framework[i].NumberOfOutOfPlanes[j],FilePtr);
            fread(Framework[i].OutOfPlaneArguments[j],sizeof(REAL[MAX_OUT_OF_PLANE_POTENTIAL_ARGUMENTS]),Framework[i].NumberOfOutOfPlanes[j],FilePtr);
          }
        }
      }


      fread(&Framework[i].NumberOfBondBondDefinitions,sizeof(int),1,FilePtr);
      Framework[i].NumberOfBondBonds=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      if(Framework[i].NumberOfBondBondDefinitions>0)
      {
        // alloc memory
        Framework[i].BondBondDefinitionType=(int*)calloc(Framework[i].NumberOfBondBondDefinitions,sizeof(int));
        Framework[i].BondBondDefinitions=(TRIPLE*)calloc(Framework[i].NumberOfBondBondDefinitions,sizeof(TRIPLE));
        Framework[i].BondBondArgumentDefinitions=(REAL(*)[MAX_BOND_BOND_POTENTIAL_ARGUMENTS])
          calloc(Framework[i].NumberOfBondBondDefinitions,sizeof(REAL[MAX_BOND_BOND_POTENTIAL_ARGUMENTS]));
        Framework[i].NumberOfBondBondsPerType=(int*)calloc(Framework[i].NumberOfBondBondDefinitions,sizeof(int));

        fread(Framework[i].BondBondDefinitionType,sizeof(int),Framework[i].NumberOfBondBondDefinitions,FilePtr);
        fread(Framework[i].BondBondDefinitions,sizeof(TRIPLE),Framework[i].NumberOfBondBondDefinitions,FilePtr);
        fread(Framework[i].BondBondArgumentDefinitions,sizeof(REAL[MAX_BOND_BOND_POTENTIAL_ARGUMENTS]),Framework[i].NumberOfBondBondDefinitions,FilePtr);
        fread(Framework[i].NumberOfBondBondsPerType,sizeof(int),Framework[i].NumberOfBondBondDefinitions,FilePtr);

        Framework[i].MaxNumberOfBondBonds=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
        Framework[i].BondBonds=(TRIPLE**)calloc(Framework[i].NumberOfFrameworks,sizeof(TRIPLE*));
        Framework[i].BondBondType=(int**)calloc(Framework[i].NumberOfFrameworks,sizeof(int*));
        Framework[i].BondBondArguments=(REAL(**)[MAX_BOND_BOND_POTENTIAL_ARGUMENTS])calloc(Framework[i].NumberOfFrameworks,
                                    sizeof(REAL(*)[MAX_BOND_BOND_POTENTIAL_ARGUMENTS]));

        fread(Framework[i].NumberOfBondBonds,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);
        fread(Framework[i].MaxNumberOfBondBonds,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);

        for(j=0;j<Framework[i].NumberOfFrameworks;j++)
        {
          // the first run, the MaxNumberOfBends is larger than NumberOfBends, but once the amount of bonds is computed
          // it will never change, so set the maximum number of bonds to the current number
          Framework[i].MaxNumberOfBondBonds[j]=Framework[i].NumberOfBondBonds[j];

          if(Framework[i].NumberOfBondBonds[j]>0)
          {
            Framework[i].BondBonds[j]=(TRIPLE*)calloc(Framework[i].NumberOfBondBonds[j],sizeof(TRIPLE));
            Framework[i].BondBondType[j]=(int*)calloc(Framework[i].NumberOfBondBonds[j],sizeof(int));
            Framework[i].BondBondArguments[j]=(REAL(*)[MAX_IMPROPER_TORSION_POTENTIAL_ARGUMENTS])
                    calloc(Framework[i].MaxNumberOfBondBonds[j],sizeof(REAL[MAX_BOND_BOND_POTENTIAL_ARGUMENTS]));

            fread(Framework[i].BondBonds[j],sizeof(TRIPLE),Framework[i].NumberOfBondBonds[j],FilePtr);
            fread(Framework[i].BondBondType[j],sizeof(int),Framework[i].NumberOfBondBonds[j],FilePtr);
            fread(Framework[i].BondBondArguments[j],sizeof(REAL[MAX_BOND_BOND_POTENTIAL_ARGUMENTS]),Framework[i].NumberOfBondBonds[j],FilePtr);
          }
        }
      }


      fread(&Framework[i].NumberOfBondBendDefinitions,sizeof(int),1,FilePtr);
      Framework[i].NumberOfBondBends=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      if(Framework[i].NumberOfBondBendDefinitions>0)
      {
        // alloc memory
        Framework[i].BondBendDefinitionType=(int*)calloc(Framework[i].NumberOfBondBendDefinitions,sizeof(int));
        Framework[i].BondBendDefinitions=(TRIPLE*)calloc(Framework[i].NumberOfBondBendDefinitions,sizeof(TRIPLE));
        Framework[i].BondBendArgumentDefinitions=(REAL(*)[MAX_BOND_BEND_POTENTIAL_ARGUMENTS])
          calloc(Framework[i].NumberOfBondBendDefinitions,sizeof(REAL[MAX_BOND_BEND_POTENTIAL_ARGUMENTS]));
        Framework[i].NumberOfBondBendsPerType=(int*)calloc(Framework[i].NumberOfBondBendDefinitions,sizeof(int));

        fread(Framework[i].BondBendDefinitionType,sizeof(int),Framework[i].NumberOfBondBendDefinitions,FilePtr);
        fread(Framework[i].BondBendDefinitions,sizeof(TRIPLE),Framework[i].NumberOfBondBendDefinitions,FilePtr);
        fread(Framework[i].BondBendArgumentDefinitions,sizeof(REAL[MAX_BOND_BEND_POTENTIAL_ARGUMENTS]),Framework[i].NumberOfBondBendDefinitions,FilePtr);
        fread(Framework[i].NumberOfBondBendsPerType,sizeof(int),Framework[i].NumberOfBondBendDefinitions,FilePtr);

        Framework[i].MaxNumberOfBondBends=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
        Framework[i].BondBends=(TRIPLE**)calloc(Framework[i].NumberOfFrameworks,sizeof(TRIPLE*));
        Framework[i].BondBendType=(int**)calloc(Framework[i].NumberOfFrameworks,sizeof(int*));
        Framework[i].BondBendArguments=(REAL(**)[MAX_BOND_BEND_POTENTIAL_ARGUMENTS])calloc(Framework[i].NumberOfFrameworks,
                                    sizeof(REAL(*)[MAX_BOND_BEND_POTENTIAL_ARGUMENTS]));

        fread(Framework[i].NumberOfBondBends,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);
        fread(Framework[i].MaxNumberOfBondBends,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);

        for(j=0;j<Framework[i].NumberOfFrameworks;j++)
        {
          // the first run, the MaxNumberOfBends is larger than NumberOfBends, but once the amount of bonds is computed
          // it will never change, so set the maximum number of bonds to the current number
          Framework[i].MaxNumberOfBondBends[j]=Framework[i].NumberOfBondBends[j];

          if(Framework[i].NumberOfBondBends[j]>0)
          {
            Framework[i].BondBends[j]=(TRIPLE*)calloc(Framework[i].NumberOfBondBends[j],sizeof(TRIPLE));
            Framework[i].BondBendType[j]=(int*)calloc(Framework[i].NumberOfBondBends[j],sizeof(int));
            Framework[i].BondBendArguments[j]=(REAL(*)[MAX_IMPROPER_TORSION_POTENTIAL_ARGUMENTS])
                    calloc(Framework[i].MaxNumberOfBondBends[j],sizeof(REAL[MAX_BOND_BEND_POTENTIAL_ARGUMENTS]));

            fread(Framework[i].BondBends[j],sizeof(TRIPLE),Framework[i].NumberOfBondBends[j],FilePtr);
            fread(Framework[i].BondBendType[j],sizeof(int),Framework[i].NumberOfBondBends[j],FilePtr);
            fread(Framework[i].BondBendArguments[j],sizeof(REAL[MAX_BOND_BEND_POTENTIAL_ARGUMENTS]),Framework[i].NumberOfBondBends[j],FilePtr);
          }
        }
      }

      fread(&Framework[i].NumberOfBendBendDefinitions,sizeof(int),1,FilePtr);
      Framework[i].NumberOfBendBends=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      if(Framework[i].NumberOfBendBendDefinitions>0)
      {
        // alloc memory
        Framework[i].BendBendDefinitionType=(int*)calloc(Framework[i].NumberOfBendBendDefinitions,sizeof(int));
        Framework[i].BendBendDefinitions=(QUAD*)calloc(Framework[i].NumberOfBendBendDefinitions,sizeof(QUAD));
        Framework[i].BendBendArgumentDefinitions=(REAL(*)[MAX_BEND_BEND_POTENTIAL_ARGUMENTS])
          calloc(Framework[i].NumberOfBendBendDefinitions,sizeof(REAL[MAX_BEND_BEND_POTENTIAL_ARGUMENTS]));
        Framework[i].NumberOfBendBendsPerType=(int*)calloc(Framework[i].NumberOfBendBendDefinitions,sizeof(int));

        fread(Framework[i].BendBendDefinitionType,sizeof(int),Framework[i].NumberOfBendBendDefinitions,FilePtr);
        fread(Framework[i].BendBendDefinitions,sizeof(QUAD),Framework[i].NumberOfBendBendDefinitions,FilePtr);
        fread(Framework[i].BendBendArgumentDefinitions,sizeof(REAL[MAX_BEND_BEND_POTENTIAL_ARGUMENTS]),Framework[i].NumberOfBendBendDefinitions,FilePtr);
        fread(Framework[i].NumberOfBendBendsPerType,sizeof(int),Framework[i].NumberOfBendBendDefinitions,FilePtr);

        Framework[i].MaxNumberOfBendBends=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
        Framework[i].BendBends=(QUAD**)calloc(Framework[i].NumberOfFrameworks,sizeof(QUAD*));
        Framework[i].BendBendType=(int**)calloc(Framework[i].NumberOfFrameworks,sizeof(int*));
        Framework[i].BendBendArguments=(REAL(**)[MAX_BEND_BEND_POTENTIAL_ARGUMENTS])calloc(Framework[i].NumberOfFrameworks,
                                    sizeof(REAL(*)[MAX_BEND_BEND_POTENTIAL_ARGUMENTS]));

        fread(Framework[i].NumberOfBendBends,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);
        fread(Framework[i].MaxNumberOfBendBends,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);

        for(j=0;j<Framework[i].NumberOfFrameworks;j++)
        {
          // the first run, the MaxNumberOfBends is larger than NumberOfBends, but once the amount of bonds is computed
          // it will never change, so set the maximum number of bonds to the current number
          Framework[i].MaxNumberOfBendBends[j]=Framework[i].NumberOfBendBends[j];

          if(Framework[i].NumberOfBendBends[j]>0)
          {
            Framework[i].BendBends[j]=(QUAD*)calloc(Framework[i].NumberOfBendBends[j],sizeof(QUAD));
            Framework[i].BendBendType[j]=(int*)calloc(Framework[i].NumberOfBendBends[j],sizeof(int));
            Framework[i].BendBendArguments[j]=(REAL(*)[MAX_IMPROPER_TORSION_POTENTIAL_ARGUMENTS])
                    calloc(Framework[i].MaxNumberOfBendBends[j],sizeof(REAL[MAX_BEND_BEND_POTENTIAL_ARGUMENTS]));

            fread(Framework[i].BendBends[j],sizeof(QUAD),Framework[i].NumberOfBendBends[j],FilePtr);
            fread(Framework[i].BendBendType[j],sizeof(int),Framework[i].NumberOfBendBends[j],FilePtr);
            fread(Framework[i].BendBendArguments[j],sizeof(REAL[MAX_BEND_BEND_POTENTIAL_ARGUMENTS]),Framework[i].NumberOfBendBends[j],FilePtr);
          }
        }
      }

      fread(&Framework[i].NumberOfBondTorsionDefinitions,sizeof(int),1,FilePtr);
      Framework[i].NumberOfBondTorsions=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      if(Framework[i].NumberOfBondTorsionDefinitions>0)
      {
        // alloc memory
        Framework[i].BondTorsionDefinitionType=(int*)calloc(Framework[i].NumberOfBondTorsionDefinitions,sizeof(int));
        Framework[i].BondTorsionDefinitions=(QUAD*)calloc(Framework[i].NumberOfBondTorsionDefinitions,sizeof(QUAD));
        Framework[i].BondTorsionArgumentDefinitions=(REAL(*)[MAX_BOND_TORSION_POTENTIAL_ARGUMENTS])
          calloc(Framework[i].NumberOfBondTorsionDefinitions,sizeof(REAL[MAX_BOND_TORSION_POTENTIAL_ARGUMENTS]));
        Framework[i].NumberOfBondTorsionsPerType=(int*)calloc(Framework[i].NumberOfBondTorsionDefinitions,sizeof(int));

        fread(Framework[i].BondTorsionDefinitionType,sizeof(int),Framework[i].NumberOfBondTorsionDefinitions,FilePtr);
        fread(Framework[i].BondTorsionDefinitions,sizeof(QUAD),Framework[i].NumberOfBondTorsionDefinitions,FilePtr);
        fread(Framework[i].BondTorsionArgumentDefinitions,sizeof(REAL[MAX_BOND_TORSION_POTENTIAL_ARGUMENTS]),Framework[i].NumberOfBondTorsionDefinitions,FilePtr);
        fread(Framework[i].NumberOfBondTorsionsPerType,sizeof(int),Framework[i].NumberOfBondTorsionDefinitions,FilePtr);

        Framework[i].MaxNumberOfBondTorsions=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
        Framework[i].BondTorsions=(QUAD**)calloc(Framework[i].NumberOfFrameworks,sizeof(QUAD*));
        Framework[i].BondTorsionType=(int**)calloc(Framework[i].NumberOfFrameworks,sizeof(int*));
        Framework[i].BondTorsionArguments=(REAL(**)[MAX_BOND_TORSION_POTENTIAL_ARGUMENTS])calloc(Framework[i].NumberOfFrameworks,
                                    sizeof(REAL(*)[MAX_BOND_TORSION_POTENTIAL_ARGUMENTS]));

        fread(Framework[i].NumberOfBondTorsions,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);
        fread(Framework[i].MaxNumberOfBondTorsions,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);

        for(j=0;j<Framework[i].NumberOfFrameworks;j++)
        {
          // the first run, the MaxNumberOfBends is larger than NumberOfBends, but once the amount of bonds is computed
          // it will never change, so set the maximum number of bonds to the current number
          Framework[i].MaxNumberOfBondTorsions[j]=Framework[i].NumberOfBondTorsions[j];

          if(Framework[i].NumberOfBondTorsions[j]>0)
          {
            Framework[i].BondTorsions[j]=(QUAD*)calloc(Framework[i].NumberOfBondTorsions[j],sizeof(QUAD));
            Framework[i].BondTorsionType[j]=(int*)calloc(Framework[i].NumberOfBondTorsions[j],sizeof(int));
            Framework[i].BondTorsionArguments[j]=(REAL(*)[MAX_IMPROPER_TORSION_POTENTIAL_ARGUMENTS])
                    calloc(Framework[i].MaxNumberOfBondTorsions[j],sizeof(REAL[MAX_BOND_TORSION_POTENTIAL_ARGUMENTS]));

            fread(Framework[i].BondTorsions[j],sizeof(QUAD),Framework[i].NumberOfBondTorsions[j],FilePtr);
            fread(Framework[i].BondTorsionType[j],sizeof(int),Framework[i].NumberOfBondTorsions[j],FilePtr);
           fread(Framework[i].BondTorsionArguments[j],sizeof(REAL[MAX_BOND_TORSION_POTENTIAL_ARGUMENTS]),Framework[i].NumberOfBondTorsions[j],FilePtr);
          }
        }
      }

      fread(&Framework[i].NumberOfBendTorsionDefinitions,sizeof(int),1,FilePtr);
      Framework[i].NumberOfBendTorsions=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      if(Framework[i].NumberOfBendTorsionDefinitions>0)
      {
        // alloc memory
        Framework[i].BendTorsionDefinitionType=(int*)calloc(Framework[i].NumberOfBendTorsionDefinitions,sizeof(int));
        Framework[i].BendTorsionDefinitions=(QUAD*)calloc(Framework[i].NumberOfBendTorsionDefinitions,sizeof(QUAD));
        Framework[i].BendTorsionArgumentDefinitions=(REAL(*)[MAX_BEND_TORSION_POTENTIAL_ARGUMENTS])
          calloc(Framework[i].NumberOfBendTorsionDefinitions,sizeof(REAL[MAX_BEND_TORSION_POTENTIAL_ARGUMENTS]));
        Framework[i].NumberOfBendTorsionsPerType=(int*)calloc(Framework[i].NumberOfBendTorsionDefinitions,sizeof(int));

        fread(Framework[i].BendTorsionDefinitionType,sizeof(int),Framework[i].NumberOfBendTorsionDefinitions,FilePtr);
        fread(Framework[i].BendTorsionDefinitions,sizeof(QUAD),Framework[i].NumberOfBendTorsionDefinitions,FilePtr);
        fread(Framework[i].BendTorsionArgumentDefinitions,sizeof(REAL[MAX_BEND_TORSION_POTENTIAL_ARGUMENTS]),Framework[i].NumberOfBendTorsionDefinitions,FilePtr);
        fread(Framework[i].NumberOfBendTorsionsPerType,sizeof(int),Framework[i].NumberOfBendTorsionDefinitions,FilePtr);

        Framework[i].MaxNumberOfBendTorsions=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
        Framework[i].BendTorsions=(QUAD**)calloc(Framework[i].NumberOfFrameworks,sizeof(QUAD*));
        Framework[i].BendTorsionType=(int**)calloc(Framework[i].NumberOfFrameworks,sizeof(int*));
        Framework[i].BendTorsionArguments=(REAL(**)[MAX_BEND_TORSION_POTENTIAL_ARGUMENTS])calloc(Framework[i].NumberOfFrameworks,
                                    sizeof(REAL(*)[MAX_BEND_TORSION_POTENTIAL_ARGUMENTS]));

        fread(Framework[i].NumberOfBendTorsions,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);
        fread(Framework[i].MaxNumberOfBendTorsions,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);

        for(j=0;j<Framework[i].NumberOfFrameworks;j++)
        {
          // the first run, the MaxNumberOfBends is larger than NumberOfBends, but once the amount of bonds is computed
          // it will never change, so set the maximum number of bonds to the current number
          Framework[i].MaxNumberOfBendTorsions[j]=Framework[i].NumberOfBendTorsions[j];

          if(Framework[i].NumberOfBendTorsions[j]>0)
          {
            Framework[i].BendTorsions[j]=(QUAD*)calloc(Framework[i].NumberOfBendTorsions[j],sizeof(QUAD));
            Framework[i].BendTorsionType[j]=(int*)calloc(Framework[i].NumberOfBendTorsions[j],sizeof(int));
            Framework[i].BendTorsionArguments[j]=(REAL(*)[MAX_IMPROPER_TORSION_POTENTIAL_ARGUMENTS])
                    calloc(Framework[i].MaxNumberOfBendTorsions[j],sizeof(REAL[MAX_BEND_TORSION_POTENTIAL_ARGUMENTS]));

            fread(Framework[i].BendTorsions[j],sizeof(QUAD),Framework[i].NumberOfBendTorsions[j],FilePtr);
            fread(Framework[i].BendTorsionType[j],sizeof(int),Framework[i].NumberOfBendTorsions[j],FilePtr);
            fread(Framework[i].BendTorsionArguments[j],sizeof(REAL[MAX_BEND_TORSION_POTENTIAL_ARGUMENTS]),Framework[i].NumberOfBendTorsions[j],FilePtr);
          }
        }
      }

      Framework[i].NumberOfIntraVDW=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      Framework[i].NumberOfExcludedIntraVDW=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));

      fread(Framework[i].NumberOfIntraVDW,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);
      fread(Framework[i].NumberOfExcludedIntraVDW,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);

      Framework[i].NumberOfIntraCharges=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      Framework[i].NumberOfIntraChargeBondDipoles=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      Framework[i].NumberOfIntraBondDipoles=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));

      fread(Framework[i].NumberOfIntraCharges,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);
      fread(Framework[i].NumberOfIntraChargeBondDipoles,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);
      fread(Framework[i].NumberOfIntraBondDipoles,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);

      Framework[i].NumberOfExcludedIntraChargeCharge=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      Framework[i].MaxNumberOfExcludedIntraChargeCharge=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      fread(Framework[i].NumberOfExcludedIntraChargeCharge,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);
      fread(Framework[i].MaxNumberOfExcludedIntraChargeCharge,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);
      Framework[i].ExcludedIntraChargeCharge=(PAIR**)calloc(Framework[i].NumberOfFrameworks,sizeof(PAIR*));

      Framework[i].NumberOfExcludedIntraChargeBondDipole=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      Framework[i].MaxNumberOfExcludedIntraChargeBondDipole=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      fread(Framework[i].NumberOfExcludedIntraChargeBondDipole,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);
      fread(Framework[i].MaxNumberOfExcludedIntraChargeBondDipole,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);
      Framework[i].ExcludedIntraChargeBondDipole=(PAIR**)calloc(Framework[i].NumberOfFrameworks,sizeof(PAIR*));

      Framework[i].NumberOfExcludedIntraBondDipoleBondDipole=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      Framework[i].MaxNumberOfExcludedIntraBondDipoleBondDipole=(int*)calloc(Framework[i].NumberOfFrameworks,sizeof(int));
      fread(Framework[i].NumberOfExcludedIntraBondDipoleBondDipole,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);
      fread(Framework[i].MaxNumberOfExcludedIntraBondDipoleBondDipole,sizeof(int),Framework[i].NumberOfFrameworks,FilePtr);
      Framework[i].ExcludedIntraBondDipoleBondDipole=(PAIR**)calloc(Framework[i].NumberOfFrameworks,sizeof(PAIR*));


      if(Framework[i].FrameworkModel==FLEXIBLE)
      {
        for(j=0;j<Framework[i].NumberOfFrameworks;j++)
        {
          if(Framework[i].FrameworkModels[j]==FLEXIBLE)
          {
            Framework[i].MaxNumberOfExcludedIntraChargeCharge[j]=Framework[i].NumberOfExcludedIntraChargeCharge[j];
            if(Framework[i].NumberOfExcludedIntraChargeCharge[j]>0)
            {
              Framework[i].ExcludedIntraChargeCharge[j]=(PAIR*)calloc(Framework[i].NumberOfExcludedIntraChargeCharge[j],sizeof(PAIR));
              fread(Framework[i].ExcludedIntraChargeCharge[j],sizeof(PAIR),Framework[i].NumberOfExcludedIntraChargeCharge[j],FilePtr);
            }

            Framework[i].MaxNumberOfExcludedIntraChargeBondDipole[j]=Framework[i].NumberOfExcludedIntraChargeBondDipole[j];
            if(Framework[i].NumberOfExcludedIntraChargeBondDipole[j])
            {
              Framework[i].ExcludedIntraChargeBondDipole[j]=(PAIR*)calloc(Framework[i].NumberOfExcludedIntraChargeBondDipole[j],sizeof(PAIR));
              fread(Framework[i].ExcludedIntraChargeBondDipole[j],sizeof(PAIR),Framework[i].NumberOfExcludedIntraChargeBondDipole[j],FilePtr);
            }

            Framework[i].MaxNumberOfExcludedIntraBondDipoleBondDipole[j]=Framework[i].NumberOfExcludedIntraBondDipoleBondDipole[j];
            if(Framework[i].NumberOfExcludedIntraBondDipoleBondDipole[j]>0)
            {
              Framework[i].ExcludedIntraBondDipoleBondDipole[j]=(PAIR*)calloc(Framework[i].NumberOfExcludedIntraBondDipoleBondDipole[j],sizeof(PAIR));
              fread(Framework[i].ExcludedIntraBondDipoleBondDipole[j],sizeof(PAIR),Framework[i].NumberOfExcludedIntraBondDipoleBondDipole[j],FilePtr);
            }
          }
        }
      }

      // rather then store the exclusion matrix, recompute it when all required information is read
      CurrentSystem=i;
      MakeExclusionMatrix(i);
    }
  }

  fread(&NumberOfSingleSubstitutionRules,sizeof(int),1,FilePtr);
  if(NumberOfSingleSubstitutionRules>0)
  {
    SubstitutionSingleFrameworkAtom=(char(*)[3][256])calloc(NumberOfSingleSubstitutionRules,sizeof(char[3][256]));
    SubstitutionSingleFrameworkAtomTypes=(int(*)[3])calloc(NumberOfSingleSubstitutionRules,sizeof(int[3]));
    fread(SubstitutionSingleFrameworkAtom,sizeof(char[3][256]),NumberOfSingleSubstitutionRules,FilePtr);
    fread(SubstitutionSingleFrameworkAtomTypes,sizeof(int[3]),NumberOfSingleSubstitutionRules,FilePtr);
  }

  fread(&NumberOfSubstitutionRules,sizeof(int),1,FilePtr);
  fread(&NumberOfRandomSubstitutions,sizeof(int),1,FilePtr);
  if(NumberOfSubstitutionRules>0)
  {
    SubstitutionFrameworkAtoms=(char(*)[3][256])calloc(NumberOfSubstitutionRules,sizeof(char[3][256]));
    SubstitutionFrameworkAtomTypes=(int(*)[3])calloc(NumberOfSubstitutionRules,sizeof(int[3]));
    fread(SubstitutionFrameworkAtoms,sizeof(char[3][256]),NumberOfSubstitutionRules,FilePtr);
    fread(SubstitutionFrameworkAtomTypes,sizeof(int[3]),NumberOfSubstitutionRules,FilePtr);
  }

  fread(&NumberOfSubstitutions,sizeof(int),1,FilePtr);
  if(NumberOfSubstitutions>0)
    fread(ListOfAtomSubstitutions,sizeof(int[3]),NumberOfSubstitutions,FilePtr);

  // read the framework atom modification rules
  fread(&NumberOfModificationRules,sizeof(int),1,FilePtr);
  if(NumberOfModificationRules>0)
  {
    ModificationRuleType=(int*)calloc(NumberOfModificationRules,sizeof(int));
    ModifyFrameworkAtoms=(char(*)[10][256])calloc(NumberOfModificationRules,sizeof(char[10][256]));
    ModifyFrameworkAtomTypes=(int(*)[10])calloc(NumberOfModificationRules,sizeof(int[10]));
    fread(ModificationRuleType,sizeof(int),NumberOfModificationRules,FilePtr);
    fread(ModifyFrameworkAtoms,sizeof(char[10][256]),NumberOfModificationRules,FilePtr);
    fread(ModifyFrameworkAtomTypes,sizeof(int[10]),NumberOfModificationRules,FilePtr);
  }
  fread(&NumberOfForbiddenConnectivityRules,sizeof(int),1,FilePtr);
  if(NumberOfForbiddenConnectivityRules>0)
  {
    ForbiddenConnectivityAtoms=(char(*)[3][256])calloc(NumberOfForbiddenConnectivityRules,sizeof(char[3][256]));
    ForbiddenConnectivityTypes=(int(*)[3])calloc(NumberOfForbiddenConnectivityRules,sizeof(int[3]));
    fread(ForbiddenConnectivityAtoms,sizeof(char[3][256]),NumberOfForbiddenConnectivityRules,FilePtr);
    fread(ForbiddenConnectivityTypes,sizeof(int[3]),NumberOfForbiddenConnectivityRules,FilePtr);
  }

  for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
  {
    AllocateFrameworkCellList();
    MakeFrameworkCellList();
  }

  fread(&Check,1,sizeof(REAL),FilePtr);
  if(fabs(Check-123456789.0)>1e-10)
  {
    fprintf(stderr, "Error in binary restart-file (ReadRestartFramework)\n");
    ContinueAfterCrash=FALSE;
  }

}
