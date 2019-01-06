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
