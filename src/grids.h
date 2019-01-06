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

#ifndef GRIDS_H
#define GRIDS_H

extern int UseTabularGrid;
extern int NumberOfGrids;
extern int *GridTypeList;
extern char (*GridTypeListName)[256];

extern float *****VDWGrid;
extern float ****CoulombGrid;

extern INT_VECTOR3 NumberOfVDWGridPoints;

extern REAL SpacingVDWGrid;
extern REAL SpacingCoulombGrid;

extern int BlockEnergyGrids;
extern int BlockGridPockets;
extern int BlockGridPores;
extern REAL BlockEnergyGridOverlapCriteria;
extern int NumberOfGridSeeds;
extern VECTOR *GridSeeds;

void MakeASCIGrid(void);


VECTOR MapToUnitCell(VECTOR pos);
VECTOR MapZToBox(VECTOR pos);
void MakeRigidFrameworkList(void);
void RigidFrameworkGrid(VECTOR pos,int typeA,REAL *Uvdw,REAL *Ucoul);

void MakeGrid(void);
int WriteVDWGrid(int l);
void ReadVDWGrid(void);
int WriteCoulombGrid(void);
void ReadCoulombGrid(void);
REAL InterpolateVDWGrid(int typeA,VECTOR pos);
REAL InterpolateVDWForceGrid(int typeA,VECTOR pos,VECTOR *Force);
REAL InterpolateCoulombGrid(int typeA,VECTOR pos);
REAL InterpolateCoulombForceGrid(int typeA,VECTOR pos,VECTOR *Force);
void TestGrid(FILE *FilePtr);
void TestForceGrid(FILE *FilePtr);
INT_VECTOR3 ConvertXYZPositionToGridIndex(VECTOR pos);
VECTOR ConvertGridIndexToXYZIndex(INT_VECTOR3 GridIndex);
void BlockingVDWGrid(void);

void WriteRestartGrids(FILE *FilePtr);
void AllocateGridMemory(void);
void ReadRestartGrids(FILE *FilePtr);
#endif
