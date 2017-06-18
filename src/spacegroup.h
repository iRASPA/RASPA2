/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2017 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'spacegroup.h' is part of RASPA-2.0

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

#ifndef SPACEGROUP_H
#define SPACEGROUP_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#define NUMBER_OF_SPACEGROUPS 530

enum {MONOCLINIC_SPACEGROUP,ORTHORHOMBIC_SPACEGROUP,TRICLINIC_SPACEGROUP,TETRAGONAL_SPACEGROUP,TRIGONAL_SPACEGROUP,HEXAGONAL_SPACEGROUP,CUBIC_SPACEGROUP};

extern int SpaceGroupSize;
extern VECTOR SpaceGroupElement[200];

typedef struct space_group_data
{
  int Identifier;
  char Name[32];
  int Number;
  int SpaceGroupOption;
  int SpaceGroupOptionCSSR;
  char ShortInternationalHermannMauguinSpaceGroupSymbol[64];
  char LongInternationalHermannMauguinSpaceGroupSymbol[64];
  char SchoenfliesSpaceGroupSymbol[32];
  char HallSpaceGroupSymbol[32];
  //char ExplicitHallSymbols[32];
  int CrystalSystem;
  char LauClass[32];
  int Centered;
  int NumberOfLatticeTranslations;
  VECTOR LatticeTranslation[4];
  int NumberOfOperators;
  int Chiral;
  int Centric;
  int Enantiomorphic;
  int Standard;
} SPACE_GROUP_DATA;

extern SPACE_GROUP_DATA SpaceGroupData[NUMBER_OF_SPACEGROUPS+1];

typedef struct string_space_group_elements
{
  int Identifier;
  char Name[192][32];
} STRING_SPACE_GROUP_ELEMENTS;

extern char StringSpaceGroupElements[NUMBER_OF_SPACEGROUPS+1][192][32];

void TestSpacegroup(void);
void SpaceGroupSymmetry(int spacegroup,VECTOR pos);
void SpaceGroupInverseSymmetry(int spacegroup,VECTOR pos);
int AsymmetricUnit(int Spacegroup,VECTOR p,REAL eps);

VECTOR ConvertToAsymetricUnitCell(VECTOR pos);
VECTOR ConvertIonsToAsymetricUnitCell(VECTOR pos);

int GetSpacegroupFromSymmetryElements(int nr,char elements[192][32]);
int GetSpaceGroupFromITSpaceGroupNumber(int sg);
int AdjustForNonStandardCSSROption(int sg,int Option);
int AdjustForNonStandardCIFOption(int sg,char *string);
int GetSpaceGroupFromHermannMauguinString(char *string);
int GetSpaceGroupFromHallString(char *hall);

int GetSpaceGroupFromHallStringWithSpacesRemoved(char *hall);
int GetSpaceGroupFromHermannMauguinStringWithSpacesRemoved(char *hm);

void PrintSpaceGroupElements(int spacegroup,FILE *Stream);

void ConvertStringElement(char *string);

void PrintSpaceGroupInformationLatex(void);

VECTOR ChiralInversion(int sg,VECTOR pos);
void PrintSpaceGroupChiralInversion(void);

#endif
