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
#include <stdbool.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>
#include "constants.h"
#include "simulation.h"
#include "molecule.h"
#include "output.h"
#include "molecule.h"
#include "framework.h"
#include "framework_energy.h"
#include "grids.h"
#include "potentials.h"
#include "spacegroup.h"
#include "movies.h"
#include "sample.h"
#include "cbmc.h"
#include "thermo_baro_stats.h"

extern bool STREAM;

REAL MovieScale;

VECTOR VTKFractionalFrameworkAtomsMin;
VECTOR VTKFractionalFrameworkAtomsMax;
VECTOR VTKFractionalFrameworkBondsMin;
VECTOR VTKFractionalFrameworkBondsMax;
VECTOR VTKFractionalAdsorbateComMin;
VECTOR VTKFractionalAdsorbateComMax;
VECTOR VTKFractionalCationComMin;
VECTOR VTKFractionalCationComMax;

int FreeEnergyAveragingTypeVTK;
int DensityAveragingTypeVTK;
int AverageDensityOverUnitCellsVTK;

int WriteVTKGrids;

typedef struct color
{
  double r;
  double g;
  double b;
  double alpha;
} COLOR;

enum {BLUE,RED,GRAY,ORANGE,YELLOW,TAN,SILVER,GREEN,WHITE,PINK,CYAN,PURPLE,LIME,MAUVRE,OCHRE,
      ICEBLUE,BLACK,YELLOW2,YELLOW3,GREEN2,GREEN3,CYAN2,CYAN3,BLUE2,BLUE3,VIOLET,VIOLET2,
      MAGENTA,MAGENTA2,OCHRE2,OCHRE3,ORANGE2,ORANGE3};

#define NUMBER_OF_COLORS 33
COLOR Colortable[NUMBER_OF_COLORS]=
{{0.0, 0.0, 1.0, 1.0}, // blue
 {1.0, 0.0, 0.0, 1.0}, // red
 {0.35,0.35,0.35,1.0}, // gray
 {1.0, 0.5, 0.0, 1.0}, // orange
 {1.0, 1.0, 0.0, 1.0}, // yellow
 {0.5, 0.5, 0.2, 1.0}, // tan
 {0.6, 0.6, 0.6, 1.0}, // silver
 {0.0, 1.0, 0.0, 1.0}, // green
 {1.0, 1.0, 1.0, 1.0}, // white
 {1.0, 0.6, 0.6, 1.0}, // pink
 {0.25,0.75,0.75,1.0}, // cyan
 {0.65,0.0, 0.65,1.0}, // purple
 {0.5, 0.9, 0.4, 1.0}, // lime
 {0.9, 0.4, 0.7, 1.0}, // mauvre
 {0.5, 0.3, 0.0, 1.0}, // ochre
 {0.5, 0.5, 0.75,1.0}, // iceblue
 {0.0, 0.0, 0.0, 1.0}, // black
 {0.88,0.97,0.02,1.0}, // yellow2
 {0.55,0.9 ,0.02,1.0}, // yellow3
 {0.0, 0.9 ,0.04,1.0}, // green2
 {0.0, 0.9 ,0.5 ,1.0}, // green3
 {0.0, 0.88,1.0 ,1.0}, // cyan2
 {0.0, 0.76,1.0 ,1.0}, // cyan3
 {0.02,0.38,0.67,1.0}, // blue2
 {0.01,0.04,0.93,1.0}, // blue3
 {0.27,0.0, 0.98,1.0}, // violet
 {0.45,0.0, 0.9, 1.0}, // violet2
 {0.9 ,0.0, 0.9, 1.0}, // magenta
 {1.0 ,0.0, 0.66,1.0}, // magenta2
 {0.98,0.0, 0.23,1.0}, // red2
 {0.81,0.0, 0.0 ,1.0}, // red3
 {0.89,0.35,0.0 ,1.0}, // orange2
 {0.96,0.72,0.0 ,1.0}};// orange3

typedef struct atom_info
{
  char *Name;
  int Color;
  REAL Mass;
  REAL CovalentRadius;
  REAL VDWRadius;
} ATOM_INFO;

#define NUMBER_OF_ELEMENTS 110
ATOM_INFO ColorElements[NUMBER_OF_ELEMENTS]=
{{"X",WHITE,0.0,0.0,0.0},         // Empty
 {"H",WHITE,1.00,0.23,1.09},      // Hydrogen
 {"He",OCHRE,4.003,1.5,1.4},      // Helium
 {"Li",OCHRE,6.941,0.68,1.82},    // Lithium
 {"Be",OCHRE,9.012,0.35,2.00},    // Beryllium
 {"B",OCHRE,10.811,0.83,2.00},    // Boron
 {"C",CYAN,12.011,0.68,1.70},     // Carbon
 {"N",BLUE,14.007,0.68,1.55},     // Nitrogen
 {"O",RED,15.999,0.68,1.52},      // Oxygen
 {"F",OCHRE,18.998,0.64,1.47},    // Fluorine
 {"Ne",OCHRE,20.180,1.50,1.54},   // Neon
 {"Na",BLUE2,22.991,0.97,2.27},   // Sodium
 {"Mg",OCHRE,24.305,1.10,1.73},   // Magnesium
 {"Al",GREEN,26.982,1.35,2.00},   // Aluminium
 {"Si",YELLOW,28.086,1.20,2.10},   // Silicon
 {"P",TAN,30.974,1.05,1.80},      // Phosphorus
 {"S",YELLOW,32.066,1.02,1.80},   // Sulphur
 {"Cl",OCHRE,35.453,0.99,1.75},   // Chlorine
 {"Ar",OCHRE,39.948,1.51,1.88},   // Argon
 {"K",OCHRE,39.098,1.33,2.75},    // Potassium
 {"Ca",OCHRE,40.078,0.99,2.00},   // Calcium
 {"Sc",OCHRE,44.956,1.44,2.00},   // Scandium
 {"Ti",OCHRE,47.867,1.47,2.00},   // Scandium
 {"V",OCHRE,50.942,1.33,2.00},    // Vanadium
 {"Cr",OCHRE,51.996,1.35,2.00},   // Chromium
 {"Mn",OCHRE,54.938,1.35,2.00},   // Manganese
 {"Fe",OCHRE,55.845,1.34,2.00},   // Iron
 {"Co",OCHRE,58.933,1.33,2.00},   // Cobalt
 {"Ni",OCHRE,58.693,1.50,1.63},   // Nickel
 {"Cu",OCHRE,63.546,1.52,1.40},   // Copper
 {"Zn",OCHRE,65.390,1.45,1.39},   // Zinc
 {"Ga",OCHRE,69.723,1.22,1.87},   // Gallium
 {"Ge",OCHRE,72.610,1.17,2.00},   // Germanium
 {"As",OCHRE,74.922,1.21,1.85},   // Arsenic
 {"Se",OCHRE,78.960,1.22,1.90},   // Selenium
 {"Br",OCHRE,79.904,1.21,1.85},   // Bromine
 {"Kr",OCHRE,83.800,1.50,2.02},   // Krypton
 {"Rb",OCHRE,85.468,1.47,2.00},   // Rubidium
 {"Sr",OCHRE,87.620,1.12,2.00},   // Strontium
 {"Y",OCHRE,88.906,1.78,2.00},    // Yttrium
 {"Zr",OCHRE,91.224,1.56,2.00},   // Zirconium
 {"Nb",OCHRE,92.906,1.48,2.00},   // Niobium
 {"Mo",OCHRE,95.940,1.47,2.00},   // Molybdenum
 {"Tc",OCHRE,98,1.35,2.00},       // Technetium
 {"Ru",OCHRE,101.070,1.40,2.00},  // Ruthenium
 {"Rh",OCHRE,102.906,1.45,2.00},  // Rhodium
 {"Pd",OCHRE,106.420,1.50,1.63},  // Palladium
 {"Ag",OCHRE,107.868,1.59,1.72},  // Silver
 {"Cd",OCHRE,112.411,1.69,1.58},  // Cadmium
 {"In",OCHRE,114.818,1.63,1.93},  // Indium
 {"Sn",OCHRE,118.71,1.46,2.17},   // Tin
 {"Sb",OCHRE,121.760,1.46,2.00},  // Antimony
 {"Te",OCHRE,127.600,1.47,2.06},  // Tellurium
 {"I",OCHRE,126.904,1.40,1.98},   // Iodine
 {"Xe",OCHRE,131.290,1.50,2.16},  // Xenon
 {"Cs",OCHRE,132.905,1.67,2.00},  // Caesium
 {"Ba",OCHRE,137.327,1.34,2.00},  // Barium
 {"Lu",OCHRE,174.967,1.72,2.00},  // Lutetium
 {"Hf",OCHRE,178.490,1.57,2.00},  // Hafnium
 {"Ta",OCHRE,180.948,1.43,2.00},  // Tantalum
 {"W",OCHRE,183.840,1.37,2.00},   // Tungsten
 {"Re",OCHRE,186.207,1.35,2.00},  // Rhenium
 {"Os",OCHRE,190.230,1.37,2.00},  // Osmium
 {"Ir",OCHRE,192.217,1.32,2.00},  // Iridium
 {"Pt",OCHRE,195.078,1.50,1.72},  // Platinum
 {"Au",OCHRE,196.967,1.50,1.66},  // Gold
 {"Hg",OCHRE,200.590,1.70,1.55},  // Mercury
 {"Tl",OCHRE,204.383,1.55,1.96},  // Thallium
 {"Pb",OCHRE,207.200,1.54,2.02},  // Lead
 {"Bi",OCHRE,208.980,1.54,2.00},  // Bismuth
 {"Po",OCHRE,210,1.68,2.00},      // Polonium
 {"At",OCHRE,210,1.21,2.00},      // Astatine
 {"Rn",OCHRE,222,1.50,2.00},      // Radon
 {"Ce",OCHRE,140.116,1.83,2.00},  // Cerium
 {"Dy",OCHRE,162.500,1.75,2.00},  // Dysprosium
 {"Er",OCHRE,167.260,1.73,2.00},  // Erbium
 {"Gd",OCHRE,157.250,1.79,2.00},  // Gadolinium
 {"Ho",OCHRE,164.930,1.74,2.00},  // Holmium
 {"La",OCHRE,138.906,1.87,2.00},  // Lanthanum
 {"Nd",OCHRE,144.240,1.81,2.00},  // Neodymium
 {"Pm",OCHRE,145,1.80,2.00},      // Promethium
 {"Pr",OCHRE,140.908,1.82,2.00},  // Praseodymium
 {"Sm",OCHRE,150.360,1.80,2.00},  // Samarium
 {"Tb",OCHRE,158.925,1.76,2.00},  // Terbium
 {"Tm",OCHRE,168.934,1.72,2.00},  // Thulium
 {"Yb",OCHRE,173.04,1.94,2.00},   // Ytterbium
 {"Fr",OCHRE,223,1.50,2.00},      // Francium
 {"Ra",OCHRE,226,1.90,2.00},      // Radium
 {"Lr",OCHRE,262,1.50,2.00},      // Lawrencium
 {"Rf",OCHRE,261,1.50,2.00},      // Rutherfordium
 {"Db",OCHRE,262,1.50,2.00},      // Dubnium
 {"Sg",OCHRE,266,1.50,2.00},      // Seaborgium
 {"Bh",OCHRE,264,1.50,2.00},      // Bohrium
 {"Hs",OCHRE,269,1.50,2.00},      // Hassium
 {"Mt",OCHRE,268,1.50,2.00},      // Meitnerium
 {"Ds",OCHRE,271,1.50,2.00},      // Darmstadtium
 {"Ac",OCHRE,227,1.88,2.00},      // Actinium
 {"Am",OCHRE,243,1.51,2.00},      // Americium
 {"Bk",OCHRE,247,1.54,2.00},      // Berkelium
 {"Cf",OCHRE,251,1.83,2.00},      // Berkelium
 {"Cm",OCHRE,247,0.99,2.00},      // Curium
 {"Es",OCHRE,252,1.50,2.00},      // Einsteinium
 {"Fm",OCHRE,257,1.50,2.00},      // Fermium
 {"Md",OCHRE,258,1.50,2.00},      // Mendelevium
 {"No",OCHRE,259,1.50,2.00},      // Nobelium
 {"Np",OCHRE,237,1.55,2.00},      // Neptunium
 {"Pa",OCHRE,231.036,1.61,2.00},  // Protactinium
 {"Pu",OCHRE,244,1.53,2.00},      // Plutonium
 {"Th",OCHRE,232.038,1.79,2.00},  // Thorium
 {"U",OCHRE,238.029,1.58,1.86}};  // Uranium


int GetColorIndex(int type)
{
  int i;

  for(i=0;i<NUMBER_OF_ELEMENTS;i++)
    if(strcasecmp(PseudoAtoms[type].PrintToPDBName,ColorElements[i].Name)==0) return ColorElements[i].Color;
  return -1;
}

REAL GetVDWRadius(int type)
{
  int i;

  for(i=0;i<NUMBER_OF_ELEMENTS;i++)
    if(strcasecmp(PseudoAtoms[type].PrintToPDBName,ColorElements[i].Name)==0) return ColorElements[i].VDWRadius;
  return -1;
}



int *Movies;
int *WriteMoviesEvery;

// PDB file-format
//  1 -  6        Record name     "ATOM  "
//  7 - 11        Integer         serial        Atom serial number
// 13 - 14        Atom            name          Chemical symbol (right justified)
//      15        Remoteness indicator
//      16        Branch designator
// 17             Character       altLoc        Alternate location indicator
// 18 - 20        Residue name    resName       Residue name
//      21        Reserved
// 22             Character       chainID       Chain identifier
// 23 - 26        Integer         resSeq        Residue sequence number
// 27             AChar           iCode         Code for insertion of residues
// 31 - 38        Real(8.3)       x             Orthogonal coordinates for X
// 39 - 46        Real(8.3)       y             Orthogonal coordinates for Y
// 47 - 54        Real(8.3)       z             Orthogonal coordinates for Z
// 55 - 60        Real(6.2)       occupancy     Occupancy
// 61 - 66        Real(6.2)       tempFactor    Isotropic B-factor
// 73 - 76        LString(4)      segID         Segment identifier, left-justified, may
//                                              include a space, e.g., CH86, A 1, NASE.
// 77 - 78        LString(2)      element       Element symbol, right-justified
// 79 - 80        LString(2)      charge        Charge on the atom
// Typical Format:  (6A1,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,6X,2A4)

// Cols.  1-6    Record name "CRYST1"
//      7-15    a (Angstrom)
//     16-24    b (Angstrom)
//     25-33    c (Angstrom)
//     34-40    alpha (degrees)
//     41-47    beta  (degrees)
//     48-54    gamma (degrees)
//     56-66    Space group symbol, left justified
//     67-70    Z   Z value is the number of polymeric chains in a unit cell. In the case of heteropolymers,
//                  Z is the number of occurrences of the most populous chain.
// Typical Format:  (6A1,3F9.3,3F7.2,1X,11A1,I4)


// transformation from the orthogonal coordinates contained in the entry to fractional
// crystallographic coordinates
// # Notes:  If the orthogonal Angstroms coordinates are X, Y, Z, and the fractional
// cell coordinates are xfrac, yfrac, zfrac, then:
// xfrac = S11X + S12Y + S13Z + U1
// yfrac = S21X + S22Y + S23Z + U2
// zfrac = S31X + S32Y + S33Z + U3

// Col. 1-6     Record name SCALEn
//  1 -  6       Record name    "SCALEn" (n=1, 2, or 3)
// 11 - 20       Real(10.6)     s[n][1]
// 21 - 30       Real(10.6)     s[n][2]
// 31 - 40       Real(10.6)     s[n][3]
// 46 - 55       Real(10.5)     u[n]

// Record:   MODEL
// Contains:  the model serial number when a single coordinate entry contains multiple structures
// # Notes:  Models are numbered sequentially beginning with 1.
// # If an entry contains more than 99,999 total atoms,
// then it must be divided among multiple models.
// # Each MODEL must have a corresponding ENDMDL record.
// # In the case of an NMR entry the EXPDTA record states the number of model structures
// that are present in the individual entry.
//  1 -  6       Record name    "MODEL "
// 11 - 14       Integer        Model serial number

// Record: ENDMDL
// Contains:  these records are paired with MODEL records to group individual structures found in a coordinate entry
// # Notes:   MODEL/ENDMDL records are used only when more than one structure
// is presented in the entry, or if there are more than 99,999 atoms.
// # Every MODEL record has an associated ENDMDL record.
//  1 -  6         Record name      "ENDMDL"

static FILE ***PDBFilePtr;
static FILE **PDBFilePtrwork;
static FILE **PDBFilePtrAll;

const int SIZE_X=150;
const int SIZE_Y=150;
const int SIZE_Z=150;


VECTOR ConvertPositionToVTKPosition(VECTOR pos3)
{
  int system;
  VECTOR Size,shift,C,pos2,pos;
  REAL max;


  system=0;
  Size.x=Size.y=Size.z=0.0;
  shift.x=shift.y=shift.z=0.0;
  C.x=1.0;
  C.y=0.0;
  C.z=0.0;
  pos.x=UnitCellBox[system].ax*C.x+UnitCellBox[system].bx*C.y+UnitCellBox[system].cx*C.z;
  pos.y=UnitCellBox[system].ay*C.x+UnitCellBox[system].by*C.y+UnitCellBox[system].cy*C.z;
  pos.z=UnitCellBox[system].az*C.x+UnitCellBox[system].bz*C.y+UnitCellBox[system].cz*C.z;
  Size.x+=fabs(pos.x);
  Size.y+=fabs(pos.y);
  Size.z+=fabs(pos.z);
  if(pos.x<0.0) shift.x+=pos.x;
  if(pos.y<0.0) shift.y+=pos.y;
  if(pos.z<0.0) shift.z+=pos.z;

  C.x=0.0;
  C.y=1.0;
  C.z=0.0;
  pos.x=UnitCellBox[system].ax*C.x+UnitCellBox[system].bx*C.y+UnitCellBox[system].cx*C.z;
  pos.y=UnitCellBox[system].ay*C.x+UnitCellBox[system].by*C.y+UnitCellBox[system].cy*C.z;
  pos.z=UnitCellBox[system].az*C.x+UnitCellBox[system].bz*C.y+UnitCellBox[system].cz*C.z;
  Size.x+=fabs(pos.x);
  Size.y+=fabs(pos.y);
  Size.z+=fabs(pos.z);
  if(pos.x<0.0) shift.x+=pos.x;
  if(pos.y<0.0) shift.y+=pos.y;
  if(pos.z<0.0) shift.z+=pos.z;

  C.x=0.0;
  C.y=0.0;
  C.z=1.0;
  pos.x=UnitCellBox[system].ax*C.x+UnitCellBox[system].bx*C.y+UnitCellBox[system].cx*C.z;
  pos.y=UnitCellBox[system].ay*C.x+UnitCellBox[system].by*C.y+UnitCellBox[system].cy*C.z;
  pos.z=UnitCellBox[system].az*C.x+UnitCellBox[system].bz*C.y+UnitCellBox[system].cz*C.z;
  Size.x+=fabs(pos.x);
  Size.y+=fabs(pos.y);
  Size.z+=fabs(pos.z);
  if(pos.x<0.0) shift.x+=pos.x;
  if(pos.y<0.0) shift.y+=pos.y;
  if(pos.z<0.0) shift.z+=pos.z;

  max=MAX2(Size.x,MAX2(Size.y,Size.z));

  pos2.x=(pos3.x-shift.x)*(REAL)SIZE_X/max;
  pos2.y=(pos3.y-shift.y)*(REAL)SIZE_Y/max;
  pos2.z=(pos3.z-shift.z)*(REAL)SIZE_Z/max;
  return pos2;

}

void WriteVTK(int system)
{
  int i,j,k,l,nr_atoms,f1,nr_bonds,count;
  int Type,TypeA,TypeB;
  FILE *FilePtr;
  VECTOR fac,pos;
  char buffer[256];
  VECTOR posA,posB,dr,Size,shift,C,s,t;
  VECTOR flexible_drift,com,com_pdb;
  REAL r,Bond_length,max,value;
  REAL tr[27][3]={{-1,-1,-1},{-1,-1,0},{-1,-1,1},{-1,0,-1},{-1,0,0},{-1,0,1},{-1,1,-1},{-1,1,0},{-1,1,1},{0,-1,-1},{0,-1,0},{0,-1,1},{0,0,-1},{0,0,0},{0,0,1},{0,1,-1},{0,1,0},{0,1,1},{1,-1,-1},{1,-1,0},{1,-1,1},{1,0,-1},{1,0,0},{1,0,1},{1,1,-1},{1,1,0},{1,1,1}};

  int **Connectivity;
  int ***Neighbours;
  VECTOR **Positions;
  int *AtomTypes;
  int *NumberOfListAtoms;
  int MaxNumberOfListAtoms;
  int TotalNumberOfListAtoms;
  int index;

  // Don't make files and folders when streaming, dammit!
  if (STREAM)
    return;

  CurrentSystem=system;

  mkdir("VTK",S_IRWXU);
  sprintf(buffer,"VTK/System_%d",system);
  mkdir(buffer,S_IRWXU);

  Size.x=Size.y=Size.z=0.0;
  shift.x=shift.y=shift.z=0.0;
  C.x=1.0;
  C.y=0.0;
  C.z=0.0;
  pos.x=UnitCellBox[system].ax*C.x+UnitCellBox[system].bx*C.y+UnitCellBox[system].cx*C.z;
  pos.y=UnitCellBox[system].ay*C.x+UnitCellBox[system].by*C.y+UnitCellBox[system].cy*C.z;
  pos.z=UnitCellBox[system].az*C.x+UnitCellBox[system].bz*C.y+UnitCellBox[system].cz*C.z;
  Size.x+=fabs(pos.x);
  Size.y+=fabs(pos.y);
  Size.z+=fabs(pos.z);
  if(pos.x<0.0) shift.x+=pos.x;
  if(pos.y<0.0) shift.y+=pos.y;
  if(pos.z<0.0) shift.z+=pos.z;

  C.x=0.0;
  C.y=1.0;
  C.z=0.0;
  pos.x=UnitCellBox[system].ax*C.x+UnitCellBox[system].bx*C.y+UnitCellBox[system].cx*C.z;
  pos.y=UnitCellBox[system].ay*C.x+UnitCellBox[system].by*C.y+UnitCellBox[system].cy*C.z;
  pos.z=UnitCellBox[system].az*C.x+UnitCellBox[system].bz*C.y+UnitCellBox[system].cz*C.z;
  Size.x+=fabs(pos.x);
  Size.y+=fabs(pos.y);
  Size.z+=fabs(pos.z);
  if(pos.x<0.0) shift.x+=pos.x;
  if(pos.y<0.0) shift.y+=pos.y;
  if(pos.z<0.0) shift.z+=pos.z;

  C.x=0.0;
  C.y=0.0;
  C.z=1.0;
  pos.x=UnitCellBox[system].ax*C.x+UnitCellBox[system].bx*C.y+UnitCellBox[system].cx*C.z;
  pos.y=UnitCellBox[system].ay*C.x+UnitCellBox[system].by*C.y+UnitCellBox[system].cy*C.z;
  pos.z=UnitCellBox[system].az*C.x+UnitCellBox[system].bz*C.y+UnitCellBox[system].cz*C.z;
  Size.x+=fabs(pos.x);
  Size.y+=fabs(pos.y);
  Size.z+=fabs(pos.z);
  if(pos.x<0.0) shift.x+=pos.x;
  if(pos.y<0.0) shift.y+=pos.y;
  if(pos.z<0.0) shift.z+=pos.z;

  max=MAX2(Size.x,MAX2(Size.y,Size.z));

  sprintf(buffer,"VTK/System_%d/Frame%s.vtk",system,FileNameAppend);
  FilePtr=fopen(buffer,"w");
  fprintf(FilePtr,"# vtk DataFile Version 1.0\n");
  fprintf(FilePtr,"Frame\n");
  fprintf(FilePtr,"ASCII\n");
  fprintf(FilePtr,"\n");
  fprintf(FilePtr,"DATASET POLYDATA\n");
  fprintf(FilePtr,"POINTS 8 float\n");

  C.x=0.0;
  C.y=0.0;
  C.z=0.0;
  pos.x=Box[system].ax*C.x+Box[system].bx*C.y+Box[system].cx*C.z;
  pos.y=Box[system].ay*C.x+Box[system].by*C.y+Box[system].cy*C.z;
  pos.z=Box[system].az*C.x+Box[system].bz*C.y+Box[system].cz*C.z;
  fprintf(FilePtr,"%lf %lf %lf\n",(double)pos.x,(double)pos.y,(double)pos.z);

  C.x=1.0;
  C.y=0.0;
  C.z=0.0;
  pos.x=Box[system].ax*C.x+Box[system].bx*C.y+Box[system].cx*C.z;
  pos.y=Box[system].ay*C.x+Box[system].by*C.y+Box[system].cy*C.z;
  pos.z=Box[system].az*C.x+Box[system].bz*C.y+Box[system].cz*C.z;
  fprintf(FilePtr,"%lf %lf %lf\n",(double)pos.x,(double)pos.y,(double)pos.z);

  C.x=1.0;
  C.y=1.0;
  C.z=0.0;
  pos.x=Box[system].ax*C.x+Box[system].bx*C.y+Box[system].cx*C.z;
  pos.y=Box[system].ay*C.x+Box[system].by*C.y+Box[system].cy*C.z;
  pos.z=Box[system].az*C.x+Box[system].bz*C.y+Box[system].cz*C.z;
  fprintf(FilePtr,"%lf %lf %lf\n",(double)pos.x,(double)pos.y,(double)pos.z);

  C.x=0.0;
  C.y=1.0;
  C.z=0.0;
  pos.x=Box[system].ax*C.x+Box[system].bx*C.y+Box[system].cx*C.z;
  pos.y=Box[system].ay*C.x+Box[system].by*C.y+Box[system].cy*C.z;
  pos.z=Box[system].az*C.x+Box[system].bz*C.y+Box[system].cz*C.z;
  fprintf(FilePtr,"%lf %lf %lf\n",(double)pos.x,(double)pos.y,(double)pos.z);

  C.x=0.0;
  C.y=0.0;
  C.z=1.0;
  pos.x=Box[system].ax*C.x+Box[system].bx*C.y+Box[system].cx*C.z;
  pos.y=Box[system].ay*C.x+Box[system].by*C.y+Box[system].cy*C.z;
  pos.z=Box[system].az*C.x+Box[system].bz*C.y+Box[system].cz*C.z;
  fprintf(FilePtr,"%lf %lf %lf\n",(double)pos.x,(double)pos.y,(double)pos.z);

  C.x=1.0;
  C.y=0.0;
  C.z=1.0;
  pos.x=Box[system].ax*C.x+Box[system].bx*C.y+Box[system].cx*C.z;
  pos.y=Box[system].ay*C.x+Box[system].by*C.y+Box[system].cy*C.z;
  pos.z=Box[system].az*C.x+Box[system].bz*C.y+Box[system].cz*C.z;
  fprintf(FilePtr,"%lf %lf %lf\n",(double)pos.x,(double)pos.y,(double)pos.z);

  C.x=1.0;
  C.y=1.0;
  C.z=1.0;
  pos.x=Box[system].ax*C.x+Box[system].bx*C.y+Box[system].cx*C.z;
  pos.y=Box[system].ay*C.x+Box[system].by*C.y+Box[system].cy*C.z;
  pos.z=Box[system].az*C.x+Box[system].bz*C.y+Box[system].cz*C.z;
  fprintf(FilePtr,"%lf %lf %lf\n",(double)pos.x,(double)pos.y,(double)pos.z);

  C.x=0.0;
  C.y=1.0;
  C.z=1.0;
  pos.x=Box[system].ax*C.x+Box[system].bx*C.y+Box[system].cx*C.z;
  pos.y=Box[system].ay*C.x+Box[system].by*C.y+Box[system].cy*C.z;
  pos.z=Box[system].az*C.x+Box[system].bz*C.y+Box[system].cz*C.z;
  fprintf(FilePtr,"%lf %lf %lf\n",(double)pos.x,(double)pos.y,(double)pos.z);

  fprintf(FilePtr,"LINES 6 36\n");
  fprintf(FilePtr,"5 0 1 2 3 0\n");
  fprintf(FilePtr,"5 4 5 6 7 4\n");
  fprintf(FilePtr,"5 0 1 5 4 0\n");
  fprintf(FilePtr,"5 2 3 7 6 2\n");
  fprintf(FilePtr,"5 0 4 7 3 0\n");
  fprintf(FilePtr,"5 1 2 6 5 1\n");
  fprintf(FilePtr,"\n");
  fclose(FilePtr);


  fac.x=Box[system].ax/MAX3(Box[system].ax,Box[system].by,Box[system].cz);
  fac.y=Box[system].by/MAX3(Box[system].ax,Box[system].by,Box[system].cz);
  fac.z=Box[system].cz/MAX3(Box[system].ax,Box[system].by,Box[system].cz);

  // allocate memory for connectivity list
  MaxNumberOfListAtoms=27*Framework[system].TotalNumberOfAtoms;
  NumberOfListAtoms=(int*)calloc(Framework[system].NumberOfFrameworks,sizeof(int));
  Positions=(VECTOR**)calloc(Framework[system].NumberOfFrameworks,sizeof(VECTOR*));
  AtomTypes=(int*)calloc(MaxNumberOfListAtoms,sizeof(int));

  Neighbours=(int***)calloc(MaxNumberOfListAtoms,sizeof(int**));
  Connectivity=(int**)calloc(Framework[system].NumberOfFrameworks,sizeof(int*));

  for(f1=0;f1<Framework[system].NumberOfFrameworks;f1++)
  {
    Positions[f1]=(VECTOR*)calloc(MaxNumberOfListAtoms,sizeof(VECTOR));
    Connectivity[f1]=(int*)calloc(MaxNumberOfListAtoms,sizeof(int));
    Neighbours[f1]=(int**)calloc(MaxNumberOfListAtoms,sizeof(int*));
    for(i=0;i<MaxNumberOfListAtoms;i++)
    {
      Neighbours[f1][i]=(int*)calloc(64,sizeof(int));
      for(l=0;l<60;l++)
        Neighbours[f1][i][l]=-1;
    }
  }

  TotalNumberOfListAtoms=0;
  for(f1=0;f1<Framework[system].NumberOfFrameworks;f1++)
  {
    NumberOfListAtoms[f1]=0;
    // make bond list
    for(i=0;i<Framework[system].NumberOfAtoms[f1];i++)
    {
      posA=Framework[system].Atoms[f1][i].Position;
      Type=Framework[system].Atoms[f1][i].Type;

      // example of a skipped atoms are virtual charge sites (Tip5p model)
      if(PseudoAtoms[Type].PrintToPDB)
      {
        s=ConvertFromXYZtoABC(posA);

        // make sure 's' is between 0 and 1.
        if(s.x<0.0) s.x+=1.0;
        if(s.y<0.0) s.y+=1.0;
        if(s.z<0.0) s.z+=1.0;
        if(s.x>=1.0) s.x-=1.0;
        if(s.y>=1.0) s.y-=1.0;
        if(s.z>=1.0) s.z-=1.0;

        for(j=0;j<27;j++)
        {
          t.x=s.x+tr[j][0];
          t.y=s.y+tr[j][1];
          t.z=s.z+tr[j][2];

          if((t.x>=VTKFractionalFrameworkAtomsMin.x)&&(t.x<=VTKFractionalFrameworkAtomsMax.x)&&
             (t.y>=VTKFractionalFrameworkAtomsMin.y)&&(t.y<=VTKFractionalFrameworkAtomsMax.y)&&
             (t.z>=VTKFractionalFrameworkAtomsMin.z)&&(t.z<=VTKFractionalFrameworkAtomsMax.z))
          {
            Positions[f1][NumberOfListAtoms[f1]]=ConvertFromABCtoXYZ(t);
            AtomTypes[NumberOfListAtoms[f1]]=Type;
            NumberOfListAtoms[f1]++;
            TotalNumberOfListAtoms++;
          }
        }
      }
    }
  }

  // determine how much the framework has drifted
  flexible_drift.x=flexible_drift.y=flexible_drift.z=0.0;
  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
  {
    com=GetFrameworkCenterOfMass();
    flexible_drift.x=com.x-Framework[CurrentSystem].IntialCenterOfMassPosition.x;
    flexible_drift.y=com.y-Framework[CurrentSystem].IntialCenterOfMassPosition.y;
    flexible_drift.z=com.z-Framework[CurrentSystem].IntialCenterOfMassPosition.z;
  }


  sprintf(buffer,"VTK/System_%d/FrameworkAtoms%s.vtk",system,FileNameAppend);
  FilePtr=fopen(buffer,"w");
  fprintf(FilePtr,"# vtk DataFile Version 1.0\n");
  fprintf(FilePtr,"Cube\n");
  fprintf(FilePtr,"ASCII\n");
  fprintf(FilePtr,"\n");
  fprintf(FilePtr,"DATASET POLYDATA\n");
  fprintf(FilePtr,"POINTS %d float\n",TotalNumberOfListAtoms);
  for(f1=0;f1<Framework[system].NumberOfFrameworks;f1++)
  {
    for(i=0;i<NumberOfListAtoms[f1];i++)
    {
      pos.x=Positions[f1][i].x-flexible_drift.x;
      pos.y=Positions[f1][i].y-flexible_drift.y;
      pos.z=Positions[f1][i].z-flexible_drift.z;
      fprintf(FilePtr,"%lf %lf %lf\n",(double)pos.x,(double)pos.y,(double)pos.z);
    }
  }

  // print the scaling factors for the size of atoms
  fprintf(FilePtr,"POINT_DATA %d\n",TotalNumberOfListAtoms);
  fprintf(FilePtr,"SCALARS my_scalars float\n");
  fprintf(FilePtr,"LOOKUP_TABLE default\n");
  for(f1=0;f1<Framework[system].NumberOfFrameworks;f1++)
  {
    for(i=0;i<NumberOfListAtoms[f1];i++)
    {
      Type=AtomTypes[i];
      fprintf(FilePtr,"%g\n",GetVDWRadius(Type));
    }
  }

  // print the colors of the atoms
  fprintf(FilePtr,"VECTORS vectors float\n");
  for(f1=0;f1<Framework[system].NumberOfFrameworks;f1++)
  {
    for(i=0;i<NumberOfListAtoms[f1];i++)
    {
      Type=AtomTypes[i];

      // transform color to a number between 0 and 1 (required by VTK)
      fprintf(FilePtr,"%g 0 0\n",(REAL)GetColorIndex(Type)/(NUMBER_OF_COLORS-1.0));
    }
  }
  fclose(FilePtr);

  for(f1=0;f1<Framework[system].NumberOfFrameworks;f1++)
  {
    for(i=0;i<MaxNumberOfListAtoms;i++)
    {
      Connectivity[f1][i]=0;
      for(l=0;l<60;l++)
        Neighbours[f1][i][l]=-1;
    }
  }

  TotalNumberOfListAtoms=0;
  for(f1=0;f1<Framework[system].NumberOfFrameworks;f1++)
  {
    NumberOfListAtoms[f1]=0;

    // make bond list
    for(i=0;i<Framework[system].NumberOfAtoms[f1];i++)
    {
      posA=Framework[system].Atoms[f1][i].Position;
      Type=Framework[system].Atoms[f1][i].Type;

      // example of a skipped atoms are virtual charge sites (Tip5p model)
      if(PseudoAtoms[Type].PrintToPDB)
      {
        s=ConvertFromXYZtoABC(posA);

        // make sure 's' is between 0 and 1.
        if(s.x<0.0) s.x+=1.0;
        if(s.y<0.0) s.y+=1.0;
        if(s.z<0.0) s.z+=1.0;
        if(s.x>=1.0) s.x-=1.0;
        if(s.y>=1.0) s.y-=1.0;
        if(s.z>=1.0) s.z-=1.0;

        for(j=0;j<27;j++)
        {
          t.x=s.x+tr[j][0];
          t.y=s.y+tr[j][1];
          t.z=s.z+tr[j][2];

          if((t.x>=VTKFractionalFrameworkBondsMin.x)&&(t.x<=VTKFractionalFrameworkBondsMax.x)&&
             (t.y>=VTKFractionalFrameworkBondsMin.y)&&(t.y<=VTKFractionalFrameworkBondsMax.y)&&
             (t.z>=VTKFractionalFrameworkBondsMin.z)&&(t.z<=VTKFractionalFrameworkBondsMax.z))
          {
            Positions[f1][NumberOfListAtoms[f1]]=ConvertFromABCtoXYZ(t);
            AtomTypes[NumberOfListAtoms[f1]]=Type;
            NumberOfListAtoms[f1]++;
            TotalNumberOfListAtoms++;
          }
        }
      }
    }

    // make bond list
    for(i=0;i<NumberOfListAtoms[f1];i++)
    {
      posA=Positions[f1][i];
      Connectivity[f1][i]=0;
      for(j=0;j<NumberOfListAtoms[f1];j++)
      {
        posB=Positions[f1][j];

        Bond_length=0.56+PseudoAtoms[AtomTypes[i]].Radius+PseudoAtoms[AtomTypes[j]].Radius;
        dr.x=posA.x-posB.x;
        dr.y=posA.y-posB.y;
        dr.z=posA.z-posB.z;

        // No boundary condition
        r=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
        if((i!=j)&&(r<Bond_length))
           Neighbours[f1][i][Connectivity[f1][i]++]=j;
      }
    }
  }

  // determine how much the framework has drifted
  flexible_drift.x=flexible_drift.y=flexible_drift.z=0.0;
  if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
  {
    com=GetFrameworkCenterOfMass();
    flexible_drift.x=com.x-Framework[CurrentSystem].IntialCenterOfMassPosition.x;
    flexible_drift.y=com.y-Framework[CurrentSystem].IntialCenterOfMassPosition.y;
    flexible_drift.z=com.z-Framework[CurrentSystem].IntialCenterOfMassPosition.z;
  }


  sprintf(buffer,"VTK/System_%d/FrameworkBonds%s.vtk",system,FileNameAppend);
  FilePtr=fopen(buffer,"w");
  fprintf(FilePtr,"# vtk DataFile Version 1.0\n");
  fprintf(FilePtr,"Cube\n");
  fprintf(FilePtr,"ASCII\n");
  fprintf(FilePtr,"\n");
  fprintf(FilePtr,"DATASET POLYDATA\n");
  fprintf(FilePtr,"POINTS %d float\n",TotalNumberOfListAtoms);
  for(f1=0;f1<Framework[system].NumberOfFrameworks;f1++)
  {
    for(i=0;i<NumberOfListAtoms[f1];i++)
    {
      pos.x=Positions[f1][i].x-flexible_drift.x;
      pos.y=Positions[f1][i].y-flexible_drift.y;
      pos.z=Positions[f1][i].z-flexible_drift.z;
      fprintf(FilePtr,"%lf %lf %lf\n",(double)pos.x,(double)pos.y,(double)pos.z);
    }
  }

  nr_bonds=0;
  for(f1=0;f1<Framework[system].NumberOfFrameworks;f1++)
  {
    // write zeolite connection tubes
    for(i=0;i<NumberOfListAtoms[f1];i++)
      nr_bonds+=Connectivity[f1][i];
  }
  fprintf(FilePtr,"LINES %d %d\n",nr_bonds/2,(int)(nr_bonds*1.5));
  index=0;
  for(f1=0;f1<Framework[system].NumberOfFrameworks;f1++)
  {
    for(i=0;i<NumberOfListAtoms[f1];i++)
    {
      for(j=0;j<Connectivity[f1][i];j++)
      {
        if(i<Neighbours[f1][i][j]) // avoid double counting
        {
          fprintf(FilePtr,"%d %d %d\n",
            2,
            i+index,
            Neighbours[f1][i][j]+index);
        }
      }
    }
    index=NumberOfListAtoms[f1];
  }
  fclose(FilePtr);



  // precompute the number of adsorbate atoms
  // this is needed at the top of the VTK-file
  // example of a skipped atoms are virtual charge sites (Tip5p model)
  nr_atoms=0;
  for(i=0;i<NumberOfAdsorbateMolecules[system];i++)
  {
    com=GetAdsorbateCenterOfMass(i);
    com.x-=flexible_drift.x;
    com.y-=flexible_drift.y;
    com.z-=flexible_drift.z;
    com_pdb=MapToBox(com);
    dr.x=com.x-com_pdb.x;
    dr.y=com.y-com_pdb.y;
    dr.z=com.z-com_pdb.z;
    for(k=0;k<27;k++)
    {
      s=ConvertFromXYZtoABC(com_pdb);

      t.x=s.x+tr[k][0];
      t.y=s.y+tr[k][1];
      t.z=s.z+tr[k][2];

      if((t.x>=VTKFractionalAdsorbateComMin.x)&&(t.x<=VTKFractionalAdsorbateComMax.x)&&
         (t.y>=VTKFractionalAdsorbateComMin.y)&&(t.y<=VTKFractionalAdsorbateComMax.y)&&
         (t.z>=VTKFractionalAdsorbateComMin.z)&&(t.z<=VTKFractionalAdsorbateComMax.z))
      {
        nr_atoms+=Components[Adsorbates[system][i].Type].NumberOfAtoms;
      }
    }
  }

  sprintf(buffer,"VTK/System_%d/AdsorbateAtoms%s.vtk",system,FileNameAppend);
  FilePtr=fopen(buffer,"w");
  fprintf(FilePtr,"# vtk DataFile Version 1.0\n");
  fprintf(FilePtr,"Cube\n");
  fprintf(FilePtr,"ASCII\n");
  fprintf(FilePtr,"\n");
  fprintf(FilePtr,"DATASET POLYDATA\n");
  fprintf(FilePtr,"POINTS %d float\n",nr_atoms);
  nr_bonds=0;
  for(i=0;i<NumberOfAdsorbateMolecules[system];i++)
  {
    com=GetAdsorbateCenterOfMass(i);
    com.x-=flexible_drift.x;
    com.y-=flexible_drift.y;
    com.z-=flexible_drift.z;
    com_pdb=MapToBox(com);
    dr.x=com.x-com_pdb.x;
    dr.y=com.y-com_pdb.y;
    dr.z=com.z-com_pdb.z;

    for(k=0;k<27;k++)
    {
      s=ConvertFromXYZtoABC(com_pdb);

      t.x=s.x+tr[k][0];
      t.y=s.y+tr[k][1];
      t.z=s.z+tr[k][2];

      if((t.x>=VTKFractionalAdsorbateComMin.x)&&(t.x<=VTKFractionalAdsorbateComMax.x)&&
         (t.y>=VTKFractionalAdsorbateComMin.y)&&(t.y<=VTKFractionalAdsorbateComMax.y)&&
         (t.z>=VTKFractionalAdsorbateComMin.z)&&(t.z<=VTKFractionalAdsorbateComMax.z))
      {
        for(j=0;j<Adsorbates[system][i].NumberOfAtoms;j++)
        {
          pos.x=Adsorbates[system][i].Atoms[j].Position.x-dr.x-flexible_drift.x;
          pos.y=Adsorbates[system][i].Atoms[j].Position.y-dr.y-flexible_drift.y;
          pos.z=Adsorbates[system][i].Atoms[j].Position.z-dr.z-flexible_drift.z;

          pos.x+=Box[CurrentSystem].ax*tr[k][0]+Box[CurrentSystem].bx*tr[k][1]+Box[CurrentSystem].cx*tr[k][2];
          pos.y+=Box[CurrentSystem].ay*tr[k][0]+Box[CurrentSystem].by*tr[k][1]+Box[CurrentSystem].cy*tr[k][2];
          pos.z+=Box[CurrentSystem].az*tr[k][0]+Box[CurrentSystem].bz*tr[k][1]+Box[CurrentSystem].cz*tr[k][2];

          fprintf(FilePtr,"%lf %lf %lf\n",(double)pos.x,(double)pos.y,(double)pos.z);
        }

        Type=Adsorbates[system][i].Type;
        for(j=0;j<Components[Type].NumberOfBonds;j++)
        {
          TypeA=Components[Type].Type[Components[Type].Bonds[j].A];
          TypeB=Components[Type].Type[Components[Type].Bonds[j].B];
          if((PseudoAtoms[TypeA].PrintToPDB)&&(PseudoAtoms[TypeB].PrintToPDB))
             nr_bonds++;
        }
      }
    }
  }

  fprintf(FilePtr,"LINES %d %d\n",nr_bonds,nr_bonds*3);
  count=0;
  for(i=0;i<NumberOfAdsorbateMolecules[system];i++)
  {
    Type=Adsorbates[system][i].Type;

    com=GetAdsorbateCenterOfMass(i);
    com.x-=flexible_drift.x;
    com.y-=flexible_drift.y;
    com.z-=flexible_drift.z;
    com_pdb=MapToBox(com);
    dr.x=com.x-com_pdb.x;
    dr.y=com.y-com_pdb.y;
    dr.z=com.z-com_pdb.z;
    for(k=0;k<27;k++)
    {
      s=ConvertFromXYZtoABC(com_pdb);

      t.x=s.x+tr[k][0];
      t.y=s.y+tr[k][1];
      t.z=s.z+tr[k][2];

      if((t.x>=VTKFractionalAdsorbateComMin.x)&&(t.x<=VTKFractionalAdsorbateComMax.x)&&
         (t.y>=VTKFractionalAdsorbateComMin.y)&&(t.y<=VTKFractionalAdsorbateComMax.y)&&
         (t.z>=VTKFractionalAdsorbateComMin.z)&&(t.z<=VTKFractionalAdsorbateComMax.z))
      {
        for(j=0;j<Components[Type].NumberOfBonds;j++)
        {
          TypeA=Components[Type].Type[Components[Type].Bonds[j].A];
          TypeB=Components[Type].Type[Components[Type].Bonds[j].B];
          if((PseudoAtoms[TypeA].PrintToPDB)&&(PseudoAtoms[TypeB].PrintToPDB))
          {
            fprintf(FilePtr,"%d %d %d\n",2,
               Components[Type].Bonds[j].A+count,
               Components[Type].Bonds[j].B+count);
          }
        }
        count+=Components[Adsorbates[system][i].Type].NumberOfAtoms;

      }
    }
  }

  // print the scaling factors for the size of atoms
  fprintf(FilePtr,"POINT_DATA %d\n",nr_atoms);
  fprintf(FilePtr,"SCALARS my_scalars float\n");
  fprintf(FilePtr,"LOOKUP_TABLE default\n");
  for(i=0;i<NumberOfAdsorbateMolecules[system];i++)
  {
    com=GetAdsorbateCenterOfMass(i);
    com.x-=flexible_drift.x;
    com.y-=flexible_drift.y;
    com.z-=flexible_drift.z;
    com_pdb=MapToBox(com);
    dr.x=com.x-com_pdb.x;
    dr.y=com.y-com_pdb.y;
    dr.z=com.z-com_pdb.z;
    for(k=0;k<27;k++)
    {
      s=ConvertFromXYZtoABC(com_pdb);

      t.x=s.x+tr[k][0];
      t.y=s.y+tr[k][1];
      t.z=s.z+tr[k][2];

      if((t.x>=VTKFractionalAdsorbateComMin.x)&&(t.x<=VTKFractionalAdsorbateComMax.x)&&
         (t.y>=VTKFractionalAdsorbateComMin.y)&&(t.y<=VTKFractionalAdsorbateComMax.y)&&
         (t.z>=VTKFractionalAdsorbateComMin.z)&&(t.z<=VTKFractionalAdsorbateComMax.z))
      {
        for(j=0;j<Adsorbates[system][i].NumberOfAtoms;j++)
        {
          Type=Adsorbates[system][i].Atoms[j].Type;
          if(PseudoAtoms[Type].PrintToPDB)
            fprintf(FilePtr,"%g\n",GetVDWRadius(Type));
          else
            fprintf(FilePtr,"%g\n",0.0);
        }
      }
    }
  }

  // print the colors of the atoms
  fprintf(FilePtr,"VECTORS vectors float\n");
  for(i=0;i<NumberOfAdsorbateMolecules[system];i++)
  {
    com=GetAdsorbateCenterOfMass(i);
    com.x-=flexible_drift.x;
    com.y-=flexible_drift.y;
    com.z-=flexible_drift.z;
    com_pdb=MapToBox(com);
    dr.x=com.x-com_pdb.x;
    dr.y=com.y-com_pdb.y;
    dr.z=com.z-com_pdb.z;
    for(k=0;k<27;k++)
    {
      s=ConvertFromXYZtoABC(com_pdb);

      t.x=s.x+tr[k][0];
      t.y=s.y+tr[k][1];
      t.z=s.z+tr[k][2];

      if((t.x>=VTKFractionalAdsorbateComMin.x)&&(t.x<=VTKFractionalAdsorbateComMax.x)&&
         (t.y>=VTKFractionalAdsorbateComMin.y)&&(t.y<=VTKFractionalAdsorbateComMax.y)&&
         (t.z>=VTKFractionalAdsorbateComMin.z)&&(t.z<=VTKFractionalAdsorbateComMax.z))
      {
        for(j=0;j<Adsorbates[system][i].NumberOfAtoms;j++)
        {
          Type=Adsorbates[system][i].Atoms[j].Type;
          if(PseudoAtoms[Type].PrintToPDB)
            fprintf(FilePtr,"%g 0 0\n",(REAL)GetColorIndex(Type)/(NUMBER_OF_COLORS-1.0));
          else
            fprintf(FilePtr,"0 0 0\n");
        }
      }
    }
  }
  fclose(FilePtr);


  sprintf(buffer,"VTK/System_%d/CationsAtoms%s.vtk",system,FileNameAppend);
  FilePtr=fopen(buffer,"w");
  fprintf(FilePtr,"# vtk DataFile Version 1.0\n");
  fprintf(FilePtr,"Cube\n");
  fprintf(FilePtr,"ASCII\n");
  fprintf(FilePtr,"\n");
  fprintf(FilePtr,"DATASET POLYDATA\n");
  fprintf(FilePtr,"POINTS %d float\n",nr_atoms);
  nr_bonds=0;
  for(i=0;i<NumberOfCationMolecules[system];i++)
  {
    com=GetCationCenterOfMass(i);
    com.x-=flexible_drift.x;
    com.y-=flexible_drift.y;
    com.z-=flexible_drift.z;
    com_pdb=MapToBox(com);
    dr.x=com.x-com_pdb.x;
    dr.y=com.y-com_pdb.y;
    dr.z=com.z-com_pdb.z;

    for(k=0;k<27;k++)
    {
      s=ConvertFromXYZtoABC(com_pdb);

      t.x=s.x+tr[k][0];
      t.y=s.y+tr[k][1];
      t.z=s.z+tr[k][2];

      if((t.x>=VTKFractionalCationComMin.x)&&(t.x<=VTKFractionalCationComMax.x)&&
         (t.y>=VTKFractionalCationComMin.y)&&(t.y<=VTKFractionalCationComMax.y)&&
         (t.z>=VTKFractionalCationComMin.z)&&(t.z<=VTKFractionalCationComMax.z))
      {
        for(j=0;j<Cations[system][i].NumberOfAtoms;j++)
        {
          pos.x=Cations[system][i].Atoms[j].Position.x-dr.x-flexible_drift.x;
          pos.y=Cations[system][i].Atoms[j].Position.y-dr.y-flexible_drift.y;
          pos.z=Cations[system][i].Atoms[j].Position.z-dr.z-flexible_drift.z;

          pos.x+=Box[CurrentSystem].ax*tr[k][0]+Box[CurrentSystem].bx*tr[k][1]+Box[CurrentSystem].cx*tr[k][2];
          pos.y+=Box[CurrentSystem].ay*tr[k][0]+Box[CurrentSystem].by*tr[k][1]+Box[CurrentSystem].cy*tr[k][2];
          pos.z+=Box[CurrentSystem].az*tr[k][0]+Box[CurrentSystem].bz*tr[k][1]+Box[CurrentSystem].cz*tr[k][2];

          fprintf(FilePtr,"%lf %lf %lf\n",(double)pos.x,(double)pos.y,(double)pos.z);
        }

        Type=Cations[system][i].Type;
        for(j=0;j<Components[Type].NumberOfBonds;j++)
        {
          TypeA=Components[Type].Type[Components[Type].Bonds[j].A];
          TypeB=Components[Type].Type[Components[Type].Bonds[j].B];
          if((PseudoAtoms[TypeA].PrintToPDB)&&(PseudoAtoms[TypeB].PrintToPDB))
             nr_bonds++;
        }
      }
    }
  }

  fprintf(FilePtr,"LINES %d %d\n",nr_bonds,nr_bonds*3);
  count=0;
  for(i=0;i<NumberOfCationMolecules[system];i++)
  {
    Type=Cations[system][i].Type;

    com=GetCationCenterOfMass(i);
    com.x-=flexible_drift.x;
    com.y-=flexible_drift.y;
    com.z-=flexible_drift.z;
    com_pdb=MapToBox(com);
    dr.x=com.x-com_pdb.x;
    dr.y=com.y-com_pdb.y;
    dr.z=com.z-com_pdb.z;
    for(k=0;k<27;k++)
    {
      s=ConvertFromXYZtoABC(com_pdb);

      t.x=s.x+tr[k][0];
      t.y=s.y+tr[k][1];
      t.z=s.z+tr[k][2];

      if((t.x>=VTKFractionalCationComMin.x)&&(t.x<=VTKFractionalCationComMax.x)&&
         (t.y>=VTKFractionalCationComMin.y)&&(t.y<=VTKFractionalCationComMax.y)&&
         (t.z>=VTKFractionalCationComMin.z)&&(t.z<=VTKFractionalCationComMax.z))
      {
        for(j=0;j<Components[Type].NumberOfBonds;j++)
        {
          TypeA=Components[Type].Type[Components[Type].Bonds[j].A];
          TypeB=Components[Type].Type[Components[Type].Bonds[j].B];
          if((PseudoAtoms[TypeA].PrintToPDB)&&(PseudoAtoms[TypeB].PrintToPDB))
          {
            fprintf(FilePtr,"%d %d %d\n",2,
               Components[Type].Bonds[j].A+count,
               Components[Type].Bonds[j].B+count);
          }
        }
        count+=Components[Cations[system][i].Type].NumberOfAtoms;

      }
    }
  }

  // print the scaling factors for the size of atoms
  fprintf(FilePtr,"POINT_DATA %d\n",nr_atoms);
  fprintf(FilePtr,"SCALARS my_scalars float\n");
  fprintf(FilePtr,"LOOKUP_TABLE default\n");
  for(i=0;i<NumberOfCationMolecules[system];i++)
  {
    com=GetCationCenterOfMass(i);
    com.x-=flexible_drift.x;
    com.y-=flexible_drift.y;
    com.z-=flexible_drift.z;
    com_pdb=MapToBox(com);
    dr.x=com.x-com_pdb.x;
    dr.y=com.y-com_pdb.y;
    dr.z=com.z-com_pdb.z;
    for(k=0;k<27;k++)
    {
      s=ConvertFromXYZtoABC(com_pdb);

      t.x=s.x+tr[k][0];
      t.y=s.y+tr[k][1];
      t.z=s.z+tr[k][2];

      if((t.x>=VTKFractionalCationComMin.x)&&(t.x<=VTKFractionalCationComMax.x)&&
         (t.y>=VTKFractionalCationComMin.y)&&(t.y<=VTKFractionalCationComMax.y)&&
         (t.z>=VTKFractionalCationComMin.z)&&(t.z<=VTKFractionalCationComMax.z))
      {
        for(j=0;j<Cations[system][i].NumberOfAtoms;j++)
        {
          Type=Cations[system][i].Atoms[j].Type;
          if(PseudoAtoms[Type].PrintToPDB)
            fprintf(FilePtr,"%g\n",GetVDWRadius(Type));
          else
            fprintf(FilePtr,"%g\n",0.0);
        }
      }
    }
  }

  // print the colors of the atoms
  fprintf(FilePtr,"VECTORS vectors float\n");
  for(i=0;i<NumberOfCationMolecules[system];i++)
  {
    com=GetCationCenterOfMass(i);
    com.x-=flexible_drift.x;
    com.y-=flexible_drift.y;
    com.z-=flexible_drift.z;
    com_pdb=MapToBox(com);
    dr.x=com.x-com_pdb.x;
    dr.y=com.y-com_pdb.y;
    dr.z=com.z-com_pdb.z;
    for(k=0;k<27;k++)
    {
      s=ConvertFromXYZtoABC(com_pdb);

      t.x=s.x+tr[k][0];
      t.y=s.y+tr[k][1];
      t.z=s.z+tr[k][2];

      if((t.x>=VTKFractionalCationComMin.x)&&(t.x<=VTKFractionalCationComMax.x)&&
         (t.y>=VTKFractionalCationComMin.y)&&(t.y<=VTKFractionalCationComMax.y)&&
         (t.z>=VTKFractionalCationComMin.z)&&(t.z<=VTKFractionalCationComMax.z))
      {
        for(j=0;j<Cations[system][i].NumberOfAtoms;j++)
        {
          Type=Cations[system][i].Atoms[j].Type;
          if(PseudoAtoms[Type].PrintToPDB)
            fprintf(FilePtr,"%g 0 0\n",(REAL)GetColorIndex(Type)/(NUMBER_OF_COLORS-1.0));
          else
            fprintf(FilePtr,"0 0 0\n");
        }
      }
    }
  }
  fclose(FilePtr);



  if(WriteVTKGrids)
  {
    if(Framework[CurrentSystem].FrameworkModel==GRID)
    {
      for(l=0;l<NumberOfGrids;l++)
      {
        sprintf(buffer,"VTK/System_%d/Grid_%d%s.vtk",system,l,FileNameAppend);
        FilePtr=fopen(buffer,"w");
        fprintf(FilePtr,"# vtk DataFile Version 1.0\n");
        fprintf(FilePtr,"energy zeolite\n");
        fprintf(FilePtr,"ASCII\n");
        fprintf(FilePtr,"DATASET STRUCTURED_POINTS\n");
        fprintf(FilePtr,"DIMENSIONS %d %d %d\n",NumberOfVDWGridPoints.x+1,
          NumberOfVDWGridPoints.y+1,NumberOfVDWGridPoints.z+1);

        max=MAX2(Size.x,MAX2(Size.y,Size.z));
        fprintf(FilePtr,"ASPECT_RATIO %lf %lf %lf\n",(double)(Size.x/max),(double)(Size.y/max),(double)(Size.z/max));
        fprintf(FilePtr,"ORIGIN 0.0 0.0 0.0\n");
        fprintf(FilePtr,"\n");
        fprintf(FilePtr,"POINT_DATA %d\n",(NumberOfVDWGridPoints.x+1)*(NumberOfVDWGridPoints.y+1)*
                (NumberOfVDWGridPoints.z+1));
        fprintf(FilePtr,"SCALARS scalars float\n");
        fprintf(FilePtr,"LOOKUP_TABLE default\n");
        for(i=0;i<=NumberOfVDWGridPoints.x;i++)
          for(j=0;j<=NumberOfVDWGridPoints.y;j++)
            for(k=0;k<=NumberOfVDWGridPoints.z;k++)
              fprintf(FilePtr,"%f\n",VDWGrid[GridTypeList[l]][i][j][k][0]);
        fclose(FilePtr);
      }
    }
    else
    {
      CurrentSystem=system;
      // compute the number of grid points
      NumberOfVDWGridPoints.x=(int)(Size.x/SpacingVDWGrid);
      NumberOfVDWGridPoints.y=(int)(Size.y/SpacingVDWGrid);
      NumberOfVDWGridPoints.z=(int)(Size.z/SpacingVDWGrid);

      for(l=0;l<NumberOfGrids;l++)
      {
        sprintf(buffer,"VTK/System_%d/Grid_%d%s.vtk",system,l,FileNameAppend);
        FilePtr=fopen(buffer,"w");
        fprintf(FilePtr,"# vtk DataFile Version 1.0\n");
        fprintf(FilePtr,"energy zeolite\n");
        fprintf(FilePtr,"ASCII\n");
        fprintf(FilePtr,"DATASET STRUCTURED_POINTS\n");
        fprintf(FilePtr,"DIMENSIONS %d %d %d\n",NumberOfVDWGridPoints.x+1,
          NumberOfVDWGridPoints.y+1,NumberOfVDWGridPoints.z+1);

        max=MAX2(Size.x,MAX2(Size.y,Size.z));
        fprintf(FilePtr,"ASPECT_RATIO %lf %lf %lf\n",(double)(Size.x/max),(double)(Size.y/max),(double)(Size.z/max));
        fprintf(FilePtr,"ORIGIN 0.0 0.0 0.0\n");
        fprintf(FilePtr,"\n");
        fprintf(FilePtr,"POINT_DATA %d\n",(NumberOfVDWGridPoints.x+1)*(NumberOfVDWGridPoints.y+1)*
                (NumberOfVDWGridPoints.z+1));
        fprintf(FilePtr,"SCALARS scalars float\n");
        fprintf(FilePtr,"LOOKUP_TABLE default\n");

        for(i=0;i<=NumberOfVDWGridPoints.x;i++)
        {
          for(j=0;j<=NumberOfVDWGridPoints.y;j++)
          {
            for(k=0;k<=NumberOfVDWGridPoints.z;k++)
            {
              posA.x=i*Size.x/NumberOfVDWGridPoints.x+shift.x;
              posA.y=j*Size.y/NumberOfVDWGridPoints.y+shift.y;
              posA.z=k*Size.z/NumberOfVDWGridPoints.z+shift.z;

              value=CalculateFrameworkVDWEnergyAtPosition(posA,GridTypeList[l],1.0);
              fprintf(FilePtr,"%f\n",value);
            }
          }
        }
        fclose(FilePtr);
      }
    }
  }

  // Free up this huge memory structure by working up the pointer tree
  for(f1=0; f1<Framework[system].NumberOfFrameworks; f1++)
  {
    for(i=0; i<MaxNumberOfListAtoms; i++) {
      free(Neighbours[f1][i]);
    }
    free(Positions[f1]);
    free(Connectivity[f1]);
    free(Neighbours[f1]);
  }

  free(Positions);
  free(Neighbours);
  free(Connectivity);
  free(AtomTypes);
  free(NumberOfListAtoms);
}

void WriteMolecule(int mol)
{
  int i;
  static int count=0;
  char buffer[256];
  VECTOR shift,com,fac;
  FILE *FilePtr;

  fac.x=UnitCellSize[0].x/MAX3(UnitCellSize[0].x,UnitCellSize[0].y,UnitCellSize[0].z);
  fac.y=UnitCellSize[0].y/MAX3(UnitCellSize[0].x,UnitCellSize[0].y,UnitCellSize[0].z);
  fac.z=UnitCellSize[0].z/MAX3(UnitCellSize[0].x,UnitCellSize[0].y,UnitCellSize[0].z);


  com.x=0.0;
  com.y=0.0;
  com.z=0.0;
  for(i=0;i<Adsorbates[CurrentSystem][mol].NumberOfAtoms;i++)
  {
    com.x+=Adsorbates[CurrentSystem][mol].Atoms[i].Position.x;
    com.y+=Adsorbates[CurrentSystem][mol].Atoms[i].Position.y;
    com.z+=Adsorbates[CurrentSystem][mol].Atoms[i].Position.z;
  }
  com.x/=Adsorbates[CurrentSystem][mol].NumberOfAtoms;
  com.y/=Adsorbates[CurrentSystem][mol].NumberOfAtoms;
  com.z/=Adsorbates[CurrentSystem][mol].NumberOfAtoms;

  shift=MapToUnitCell(com);
  shift.x-=com.x;
  shift.y-=com.y;
  shift.z-=com.z;

  sprintf(buffer,"molecule_%d%s.vtk",count,FileNameAppend);
  FilePtr=fopen(buffer,"w");
  fprintf(FilePtr,"# vtk DataFile Version 1.0\n");
  fprintf(FilePtr,"Cube\n");
  fprintf(FilePtr,"ASCII\n");
  fprintf(FilePtr,"\n");
  fprintf(FilePtr,"DATASET POLYDATA\n");
  fprintf(FilePtr,"POINTS %d float\n",Adsorbates[CurrentSystem][mol].NumberOfAtoms);
  for(i=0;i<Adsorbates[CurrentSystem][mol].NumberOfAtoms;i++)
    fprintf(FilePtr,"%lf %lf %lf\n",
      (double)(fac.x*(Adsorbates[CurrentSystem][mol].Atoms[i].Position.x+shift.x)*150.0/UnitCellSize[0].x),
      (double)(fac.y*(Adsorbates[CurrentSystem][mol].Atoms[i].Position.y+shift.y)*150.0/UnitCellSize[0].y),
      (double)(fac.z*(Adsorbates[CurrentSystem][mol].Atoms[i].Position.z+shift.z)*150.0/UnitCellSize[0].z));
  fprintf(FilePtr,"LINES 1 %d\n",Adsorbates[CurrentSystem][mol].NumberOfAtoms+1);
  fprintf(FilePtr,"%d",Adsorbates[CurrentSystem][mol].NumberOfAtoms);
  for(i=0;i<Adsorbates[CurrentSystem][mol].NumberOfAtoms;i++)
    fprintf(FilePtr," %d",i);
  fprintf(FilePtr,"\n");
  fclose(FilePtr);
  count++;
}

static int *nr_frame;
int SamplePDBMovies(int Choice,int Subdir)
{
  int i,j,k,f1;
  char buffer[1024];
  int Type;
  int AtomType;
  int index;
  char RecordName[7]="ATOM  ";                   // ATOM record
  static int  *SerialNumber;   // Atom serial number
  int  SerialNumberAll,SerialNumberFramework;
  char AtomName[5]=" C";                         // Atom Name
  char RemotenessIndicator=' ';                  // Remoteness indicator
  char BranchDesignator=' ';                     // Branch designator
  char AltLoc=' ';                               // Alternate location indicator
  char ResIdueName[4]="MOL";                     // ResIdue Name
  char ChainId=' ';                              // Chain Identifier
  char ResSeq[5]="    ";                         // ResIdue sequence number
  char iCode=' ';                                // code for insertion of resIdues
  REAL Occupancy=1.0;                            // Occupancy
  REAL Temp=0.0;                                 // Temperature factor
  char SegID[5]="    ";                          // Segment Identifier, left-justified
  char Element[3]="  ";                          // Element symbol, right-justified
  char charge[3]="  ";                           // Charge
  VECTOR dr,r,com,com_pdb,flexible_drift;

  if (STREAM)
    return 0;

  switch(Choice)
  {
    case ALLOCATE:
      nr_frame=(int*)calloc(NumberOfSystems,sizeof(int));
      SerialNumber=(int*)calloc(NumberOfComponents,sizeof(int));
      break;
    case INITIALIZE:
      mkdir("Movies",S_IRWXU);
      for(i=0;i<NumberOfSystems;i++)
      {
        nr_frame[i]=1;
        if(Movies[i])
        {
          if(Subdir<0)
          {
            sprintf(buffer,"Movies/System_%d",i);
            mkdir(buffer,S_IRWXU);

            for(j=0;j<NumberOfComponents;j++)
            {
              sprintf(buffer,"Movies/System_%d/Movie_%s_%d.%d.%d_%lf_%lf_component_%s_%d%s.pdb",
                    i,
                    Framework[i].Name[0],
                    NumberOfUnitCells[i].x,
                    NumberOfUnitCells[i].y,
                    NumberOfUnitCells[i].z,
                    (double)therm_baro_stats.ExternalTemperature[i],
                    (double)(therm_baro_stats.ExternalPressure[i][CurrentIsothermPressure]*PRESSURE_CONVERSION_FACTOR),
                    Components[j].Name,
                    j,
                    FileNameAppend);
              PDBFilePtr[i][j]=fopen(buffer,"w");
            }

            sprintf(buffer,"Movies/System_%d/Movie_%s_%d.%d.%d_%lf_%lf_allcomponents%s.pdb",
                    i,
                    Framework[i].Name[0],
                    NumberOfUnitCells[i].x,
                    NumberOfUnitCells[i].y,
                    NumberOfUnitCells[i].z,
                    (double)therm_baro_stats.ExternalTemperature[i],
                    (double)(therm_baro_stats.ExternalPressure[i][CurrentIsothermPressure]*PRESSURE_CONVERSION_FACTOR),
                    FileNameAppend);
            PDBFilePtrAll[i]=fopen(buffer,"w");

            sprintf(buffer,"Movies/System_%d/Movie_%s_%d.%d.%d_%lf_%lf_frameworks%s.pdb",
                    i,
                    Framework[i].Name[0],
                    NumberOfUnitCells[i].x,
                    NumberOfUnitCells[i].y,
                    NumberOfUnitCells[i].z,
                    (double)therm_baro_stats.ExternalTemperature[i],
                    (double)(therm_baro_stats.ExternalPressure[i][CurrentIsothermPressure]*PRESSURE_CONVERSION_FACTOR),
                    FileNameAppend);
            PDBFilePtrwork[i]=fopen(buffer,"w");
          }
          else
          {
            sprintf(buffer,"Movies/System_%d",i);
            mkdir(buffer,S_IRWXU);

            sprintf(buffer,"Movies/System_%d/Run_%d",i,Subdir);
            mkdir(buffer,S_IRWXU);

            for(j=0;j<NumberOfComponents;j++)
            {
              sprintf(buffer,"Movies/System_%d/Run_%d/Movie_%s_%d.%d.%d_%lf_%lf_component_%s_%d%s.pdb",
                    i,
                    Subdir,
                    Framework[i].Name[0],
                    NumberOfUnitCells[i].x,
                    NumberOfUnitCells[i].y,
                    NumberOfUnitCells[i].z,
                    (double)therm_baro_stats.ExternalTemperature[i],
                    (double)(therm_baro_stats.ExternalPressure[i][CurrentIsothermPressure]*PRESSURE_CONVERSION_FACTOR),
                    Components[j].Name,
                    j,
                    FileNameAppend);
              PDBFilePtr[i][j]=fopen(buffer,"w");
            }

            sprintf(buffer,"Movies/System_%d/Run_%d/Movie_%s_%d.%d.%d_%lf_%lf_allcomponents%s.pdb",
                    i,
                    Subdir,
                    Framework[i].Name[0],
                    NumberOfUnitCells[i].x,
                    NumberOfUnitCells[i].y,
                    NumberOfUnitCells[i].z,
                    (double)therm_baro_stats.ExternalTemperature[i],
                    (double)(therm_baro_stats.ExternalPressure[i][CurrentIsothermPressure]*PRESSURE_CONVERSION_FACTOR),
                    FileNameAppend);
            PDBFilePtrAll[i]=fopen(buffer,"w");

            sprintf(buffer,"Movies/System_%d/Run_%d/Movie_%s_%d.%d.%d_%lf_%lf_frameworks%s.pdb",
                    i,
                    Subdir,
                    Framework[i].Name[0],
                    NumberOfUnitCells[i].x,
                    NumberOfUnitCells[i].y,
                    NumberOfUnitCells[i].z,
                    (double)therm_baro_stats.ExternalTemperature[i],
                    (double)(therm_baro_stats.ExternalPressure[i][CurrentIsothermPressure]*PRESSURE_CONVERSION_FACTOR),
                    FileNameAppend);
            PDBFilePtrwork[i]=fopen(buffer,"w");
          }
        }
      }
      break;
    case SAMPLE:
      if(!(Movies[CurrentSystem]&&(CurrentCycle%WriteMoviesEvery[CurrentSystem]==0))) return 0;

      flexible_drift.x=flexible_drift.y=flexible_drift.z=0.0;
      if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
      {
        com=GetFrameworkCenterOfMass();
        flexible_drift.x=com.x-Framework[CurrentSystem].IntialCenterOfMassPosition.x;
        flexible_drift.y=com.y-Framework[CurrentSystem].IntialCenterOfMassPosition.y;
        flexible_drift.z=com.z-Framework[CurrentSystem].IntialCenterOfMassPosition.z;
      }

      index=1;

      fprintf(PDBFilePtrwork[CurrentSystem],"MODEL %4d\n",nr_frame[CurrentSystem]);
      fprintf(PDBFilePtrwork[CurrentSystem],"CRYST1%9.3lf%9.3lf%9.3lf%7.2lf%7.2lf%7.2lf\n",
          (double)(MovieScale*BoxProperties[CurrentSystem].ax),
          (double)(MovieScale*BoxProperties[CurrentSystem].ay),
          (double)(MovieScale*BoxProperties[CurrentSystem].az),
          (double)AlphaAngle[CurrentSystem]*RAD2DEG,
          (double)BetaAngle[CurrentSystem]*RAD2DEG,
          (double)GammaAngle[CurrentSystem]*RAD2DEG);

      SerialNumberFramework=1;

      fprintf(PDBFilePtrAll[CurrentSystem],"MODEL %4d\n",nr_frame[CurrentSystem]);
      fprintf(PDBFilePtrAll[CurrentSystem],"CRYST1%9.3lf%9.3lf%9.3lf%7.2lf%7.2lf%7.2lf\n",
          (double)(MovieScale*BoxProperties[CurrentSystem].ax),
          (double)(MovieScale*BoxProperties[CurrentSystem].ay),
          (double)(MovieScale*BoxProperties[CurrentSystem].az),
          (double)AlphaAngle[CurrentSystem]*RAD2DEG,
          (double)BetaAngle[CurrentSystem]*RAD2DEG,
          (double)GammaAngle[CurrentSystem]*RAD2DEG);


      SerialNumberAll=1;

      for(j=0;j<NumberOfComponents;j++)
      {
        fprintf(PDBFilePtr[CurrentSystem][j],"MODEL %4d\n",nr_frame[CurrentSystem]);
        fprintf(PDBFilePtr[CurrentSystem][j],"CRYST1%9.3lf%9.3lf%9.3lf%7.2lf%7.2lf%7.2lf\n",
          (double)(MovieScale*BoxProperties[CurrentSystem].ax),
          (double)(MovieScale*BoxProperties[CurrentSystem].ay),
          (double)(MovieScale*BoxProperties[CurrentSystem].az),
          (double)AlphaAngle[CurrentSystem]*RAD2DEG,
          (double)BetaAngle[CurrentSystem]*RAD2DEG,
          (double)GammaAngle[CurrentSystem]*RAD2DEG);

        SerialNumber[j]=1;
      }
      nr_frame[CurrentSystem]++;

      if(Framework[CurrentSystem].FrameworkModel==FLEXIBLE)
      {
        SerialNumberFramework=1;
        Type=NumberOfComponents;
        for(f1=0;f1<Framework[CurrentSystem].NumberOfFrameworks;f1++)
        {
          for(k=0;k<Framework[CurrentSystem].NumberOfAtoms[f1];k++)
          {
            if(PseudoAtoms[Framework[CurrentSystem].Atoms[f1][k].Type].PrintToPDB)
            {
              AtomType=Framework[CurrentSystem].Atoms[f1][k].Type;
              snprintf(AtomName,5,"%2s",PseudoAtoms[AtomType].PrintToPDBName);
              snprintf(Element,3,"%2s",PseudoAtoms[AtomType].PrintToPDBName);

              // shift positions to remove the drift of the framework
              r.x=MovieScale*(Framework[CurrentSystem].Atoms[f1][k].Position.x-flexible_drift.x);
              r.y=MovieScale*(Framework[CurrentSystem].Atoms[f1][k].Position.y-flexible_drift.y);
              r.z=MovieScale*(Framework[CurrentSystem].Atoms[f1][k].Position.z-flexible_drift.z);
              fprintf(PDBFilePtrwork[CurrentSystem],"%s%5d %2s%c%c%c%3s %c%4s%c   %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf%6s%-4s%2s%2s\n",
                   RecordName,SerialNumberFramework++,AtomName,RemotenessIndicator,BranchDesignator,AltLoc,ResIdueName,
                   ChainId,ResSeq,iCode,(double)r.x,(double)r.y,(double)r.z,(double)Occupancy,(double)Temp,"      ",SegID,Element,charge);
              fprintf(PDBFilePtrAll[CurrentSystem],"%s%5d %2s%c%c%c%3s %c%4s%c   %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf%6s%-4s%2s%2s\n",
                   RecordName,SerialNumberAll++,AtomName,RemotenessIndicator,BranchDesignator,AltLoc,ResIdueName,
                   ChainId,ResSeq,iCode,(double)r.x,(double)r.y,(double)r.z,(double)Occupancy,(double)Temp,"      ",SegID,Element,charge);
            }
          }
        }
      }

      for(j=0;j<NumberOfAdsorbateMolecules[CurrentSystem];j++)
      {
        Type=Adsorbates[CurrentSystem][j].Type;

        com=GetAdsorbateCenterOfMass(j);
        com.x-=flexible_drift.x;
        com.y-=flexible_drift.y;
        com.z-=flexible_drift.z;
        com_pdb=MapToBox(com);
        dr.x=com.x-com_pdb.x;
        dr.y=com.y-com_pdb.y;
        dr.z=com.z-com_pdb.z;

        for(k=0;k<Adsorbates[CurrentSystem][j].NumberOfAtoms;k++)
        {
          if(PseudoAtoms[Adsorbates[CurrentSystem][j].Atoms[k].Type].PrintToPDB)
          {
            r.x=MovieScale*(Adsorbates[CurrentSystem][j].Atoms[k].Position.x-dr.x-flexible_drift.x);
            r.y=MovieScale*(Adsorbates[CurrentSystem][j].Atoms[k].Position.y-dr.y-flexible_drift.y);
            r.z=MovieScale*(Adsorbates[CurrentSystem][j].Atoms[k].Position.z-dr.z-flexible_drift.z);

            Occupancy = 1.0;
            if (IsFractionalReactionAdsorbateMolecule(j))
            {
              Occupancy = 0.5;
            }
	    Occupancy=Adsorbates[CurrentSystem][j].Atoms[k].CFVDWScalingParameter;

            AtomType=Adsorbates[CurrentSystem][j].Atoms[k].Type;
            snprintf(AtomName,5,"%2s",PseudoAtoms[AtomType].PrintToPDBName);
            snprintf(Element,3,"%2s",PseudoAtoms[AtomType].PrintToPDBName);

            fprintf(PDBFilePtr[CurrentSystem][Type],"%s%5d %2s%c%c%c%3s %c%4s%c   %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf%6s%-4s%2s%2s\n",
                 RecordName,SerialNumber[Type]++,AtomName,RemotenessIndicator,BranchDesignator,AltLoc,ResIdueName,
                 ChainId,ResSeq,iCode,(double)r.x,(double)r.y,(double)r.z,(double)Occupancy,(double)Temp,"      ",SegID,Element,charge);

            fprintf(PDBFilePtrAll[CurrentSystem],"%s%5d %2s%c%c%c%3s %c%4s%c   %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf%6s%-4s%2s%2s\n",
                 RecordName,SerialNumberAll++,AtomName,RemotenessIndicator,BranchDesignator,AltLoc,ResIdueName,
                 ChainId,ResSeq,iCode,(double)r.x,(double)r.y,(double)r.z,(double)Occupancy,(double)Temp,"      ",SegID,Element,charge);
          }
        }
      }

      for(j=0;j<NumberOfCationMolecules[CurrentSystem];j++)
      {
        Type=Cations[CurrentSystem][j].Type;

        com=GetCationCenterOfMass(j);
        com.x-=flexible_drift.x;
        com.y-=flexible_drift.y;
        com.z-=flexible_drift.z;
        com_pdb=MapToBox(com);
        dr.x=com.x-com_pdb.x;
        dr.y=com.y-com_pdb.y;
        dr.z=com.z-com_pdb.z;

        for(k=0;k<Cations[CurrentSystem][j].NumberOfAtoms;k++)
        {
          if(PseudoAtoms[Cations[CurrentSystem][j].Atoms[k].Type].PrintToPDB)
          {
            r.x=MovieScale*(Cations[CurrentSystem][j].Atoms[k].Position.x-dr.x-flexible_drift.x);
            r.y=MovieScale*(Cations[CurrentSystem][j].Atoms[k].Position.y-dr.y-flexible_drift.y);
            r.z=MovieScale*(Cations[CurrentSystem][j].Atoms[k].Position.z-dr.z-flexible_drift.z);

            AtomType=Cations[CurrentSystem][j].Atoms[k].Type;
            snprintf(AtomName,5,"%2s",PseudoAtoms[AtomType].PrintToPDBName);
            snprintf(Element,3,"%2s",PseudoAtoms[AtomType].PrintToPDBName);

            Occupancy = 1.0;
            if (IsFractionalReactionCationMolecule(j))
            {
              Occupancy = 0.5;
            }
	    Occupancy=Cations[CurrentSystem][j].Atoms[k].CFVDWScalingParameter;

            fprintf(PDBFilePtr[CurrentSystem][Type],"%s%5d %2s%c%c%c%3s %c%4s%c   %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf%6s%-4s%2s%2s\n",
                 RecordName,SerialNumber[Type]++,AtomName,RemotenessIndicator,BranchDesignator,AltLoc,ResIdueName,
                 ChainId,ResSeq,iCode,(double)r.x,(double)r.y,(double)r.z,(double)Occupancy,(double)Temp,"      ",SegID,Element,charge);

            fprintf(PDBFilePtrAll[CurrentSystem],"%s%5d %2s%c%c%c%3s %c%4s%c   %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf%6s%-4s%2s%2s\n",
                 RecordName,SerialNumberAll++,AtomName,RemotenessIndicator,BranchDesignator,AltLoc,ResIdueName,
                 ChainId,ResSeq,iCode,(double)r.x,(double)r.y,(double)r.z,(double)Occupancy,(double)Temp,"      ",SegID,Element,charge);
          }
        }
      }

      fprintf(PDBFilePtrwork[CurrentSystem],"ENDMDL\n");
      fprintf(PDBFilePtrAll[CurrentSystem],"ENDMDL\n");
      for(j=0;j<NumberOfComponents;j++)
        fprintf(PDBFilePtr[CurrentSystem][j],"ENDMDL\n");
      break;
    case PRINT:
      break;
    case FINALIZE:
      for(i=0;i<NumberOfSystems;i++)
      {
        if(Movies[i])
        {
          fclose(PDBFilePtrAll[i]);
          fclose(PDBFilePtrwork[i]);
          for(j=0;j<NumberOfComponents;j++)
           fclose(PDBFilePtr[i][j]);
        }
      }
      break;
  }
  return 0;
}

void WriteSnapshotIonsCssr(void)
{
  int i;
  REAL A,B,C;
  VECTOR pos;
  FILE *FilePtr;

  if (STREAM)
    return;

  A=(REAL)NumberOfUnitCells[0].x*UnitCellSize[0].x;
  B=(REAL)NumberOfUnitCells[0].y*UnitCellSize[0].y;
  C=(REAL)NumberOfUnitCells[0].z*UnitCellSize[0].z;

  FilePtr=fopen("IonsSnapshot.cssr","w");
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
  fprintf(FilePtr,"%4d%4d %s\n",NumberOfCationMolecules[CurrentSystem],0,"Created by Raspa-1.0");
  fprintf(FilePtr,"     0 %s         : %s\n",Framework[CurrentSystem].Name[0],Framework[CurrentSystem].Name[0]);
  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    pos=ConvertFromXYZtoABC(Cations[CurrentSystem][i].Atoms[0].Position);

    if(fabs(pos.x)<1e-10) pos.x=fabs(pos.x);
    if(fabs(pos.y)<1e-10) pos.y=fabs(pos.y);
    if(fabs(pos.z)<1e-10) pos.z=fabs(pos.z);

    fprintf(FilePtr,"%4d %-4s  %9.5lf %9.5lf %9.5lf %4d%4d%4d%4d%4d%4d%4d%4d %7.3lf\n",
      i+1,
      PseudoAtoms[Cations[CurrentSystem][i].Atoms[0].Type].Name,
      (double)pos.x,
      (double)pos.y,
      (double)pos.z,
      0,0,0,0,0,0,0,0,
      (double)0.0);
  }
  fclose(FilePtr);
}

void WriteSnapshotIonsCssrUsingSymmetry(void)
{
  int i,j,index;
  REAL A,B,C;
  VECTOR pos;
  FILE *FilePtr;
  VECTOR s;

  A=(REAL)NumberOfUnitCells[0].x*UnitCellSize[0].x;
  B=(REAL)NumberOfUnitCells[0].y*UnitCellSize[0].y;
  C=(REAL)NumberOfUnitCells[0].z*UnitCellSize[0].z;

  FilePtr=fopen("IonsSnapshotSymmetry.cssr","w");
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
  fprintf(FilePtr,"%4d%4d %s\n",NumberOfCationMolecules[CurrentSystem]*SpaceGroupSize,0,"Created by Raspa-1.0");
  fprintf(FilePtr,"     0 %s         : %s\n",Framework[CurrentSystem].Name[0],Framework[CurrentSystem].Name[0]);

  index=1;
  for(i=0;i<NumberOfCationMolecules[CurrentSystem];i++)
  {
    pos=Cations[CurrentSystem][i].Atoms[0].Position;
    s=ConvertToAsymetricUnitCell(pos);

    SpaceGroupSymmetry(Framework[CurrentSystem].SpaceGroupIdentifier[CurrentFramework],s);
    for(j=0;j<SpaceGroupSize;j++)
    {
      s.x=SpaceGroupElement[j].x-NINT(SpaceGroupElement[j].x);
      s.y=SpaceGroupElement[j].y-NINT(SpaceGroupElement[j].y);
      s.z=SpaceGroupElement[j].z-NINT(SpaceGroupElement[j].z);

      if(s.x<0.0) s.x+=1.0;
      if(s.y<0.0) s.y+=1.0;
      if(s.z<0.0) s.z+=1.0;

      fprintf(FilePtr,"%4d %-4s  %9.5lf %9.5lf %9.5lf %4d%4d%4d%4d%4d%4d%4d%4d %7.3lf\n",
        index++,
        "Ca",
        (double)s.x,
        (double)s.y,
        (double)s.z,
        0,0,0,0,0,0,0,0,
        (double)0.0);
    }
  }
  fclose(FilePtr);
}

void WriteAsymetricPositions(void)
{
  int i,j,k;
  char buffer[1024];
  FILE *FilePtr;
  VECTOR pos;

  if (STREAM)
    return;

  mkdir("Movies",S_IRWXU);
  for(i=0;i<NumberOfSystems;i++)
  {
    sprintf(buffer,"Movies/System_%d",i);
    mkdir(buffer,S_IRWXU);
  }

  for(k=0;k<NumberOfSystems;k++)
  {
    sprintf(buffer,"Movies/System_%d/AsymetricPositions.dat",k);
    FilePtr=fopen(buffer,"w");

    for(i=0;i<NumberOfAdsorbateMolecules[k];i++)
    {
      for(j=0;j<Adsorbates[k][i].NumberOfAtoms;j++)
      {
        pos=ConvertToAsymetricUnitCell(Adsorbates[k][i].Atoms[j].Position);
        fprintf(FilePtr,"%18.10f %18.10f %18.10f\n",pos.x,pos.y,pos.z);
      }
    }
    fclose(FilePtr);
  }
}

void WriteDlpolyInputFiles(void)
{
  int i,j;
  FILE *FilePtr;
  char buffer[1024];
  int total_atoms;

  mkdir("DlpolyInputFiles",S_IRWXU);

/*
  sprintf(buffer,"DlpolyInputFiles/FIELD");
  FilePtr=fopen(buffer,"w");
  fclose(FilePtr);
*/

  sprintf(buffer,"DlpolyInputFiles/CONFIG");
  FilePtr=fopen(buffer,"w");

  fprintf(FilePtr,"Raspa-1.0: input files for dlpoly\n");

  fprintf(FilePtr,"%10d%10d\n",0,3);
  fprintf(FilePtr,"%20f%20f%20f\n",Box[0].ax,Box[0].ay,Box[0].az);
  fprintf(FilePtr,"%20f%20f%20f\n",Box[0].bx,Box[0].by,Box[0].bz);
  fprintf(FilePtr,"%20f%20f%20f\n",Box[0].cx,Box[0].cy,Box[0].cz);

  total_atoms=1;
  for(i=0;i<NumberOfAdsorbateMolecules[0];i++)
  {
    for(j=0;j<Adsorbates[0][i].NumberOfAtoms;j++)
    {
      fprintf(FilePtr,"%s%10d\n",PseudoAtoms[Adsorbates[0][i].Atoms[j].Type].Name,total_atoms++);
      fprintf(FilePtr,"%20f%20f%20f\n",Adsorbates[0][i].Atoms[j].Position.x,
           Adsorbates[0][i].Atoms[j].Position.y,Adsorbates[0][i].Atoms[j].Position.z);
    }
  }

  fclose(FilePtr);
}

void AllocateMovieMemory(void)
{
  int i;

  Movies=(int*)calloc(NumberOfSystems,sizeof(int));
  WriteMoviesEvery=(int*)calloc(NumberOfSystems,sizeof(int));
  PDBFilePtrwork=(FILE**)calloc(NumberOfSystems,sizeof(FILE*));
  PDBFilePtrAll=(FILE**)calloc(NumberOfSystems,sizeof(FILE*));
  PDBFilePtr=(FILE***)calloc(NumberOfSystems,sizeof(FILE**));
  for(i=0;i<NumberOfSystems;i++)
    PDBFilePtr[i]=(FILE**)calloc(NumberOfComponents,sizeof(FILE*));
}


REAL *Histogram_3D;
REAL *Histogram_Count_3D;

void FreeEnergyProfile3D(void)
{
  int i;
  int x,y,z,typeA,index,temp,start;
  REAL value,min;
  VECTOR pos,s;
  char buffer[256];
  FILE *FilePtr;
  REAL_MATRIX3x3 Cell;
  REAL_MATRIX3x3 InverseCell;
  REAL_MATRIX3x3 cellProperties;
  REAL det;

  // Modify overlap criteria to remove artifacts in the VTK pictures
  EnergyOverlapCriteria=5e4;

  Histogram_3D=(REAL*)calloc(SIZE_X*SIZE_Y*SIZE_Z,sizeof(REAL));
  Histogram_Count_3D=(REAL*)calloc(SIZE_X*SIZE_Y*SIZE_Z,sizeof(REAL));

  CurrentSystem=0;
  CurrentComponent=0;
  typeA=Components[0].Type[0];
  ChargeMethod=NONE;

  switch(FreeEnergyAveragingTypeVTK)
  {
    case VTK_UNIT_CELL:
      Cell=UnitCellBox[0];
      InverseCell=InverseUnitCellBox[0];
      CellProperties(&UnitCellBox[0],&cellProperties,&det);
      break;
    case VTK_FULL_BOX:
    default:
      Cell=Box[CurrentSystem];
      InverseCell=InverseBox[CurrentSystem];
      CellProperties(&Box[0],&cellProperties,&det);
      break;
  }

  for(i=0;i<NumberOfCycles;i++)
  {
    if(i%PrintEvery==0) fprintf(stderr, "iteration: %d\n",i);

    // generate random number in enclosed box
    
    s.x=RandomNumber();
    s.y=RandomNumber();
    s.z=RandomNumber();

    switch(FreeEnergyAveragingTypeVTK)
    {
      case VTK_UNIT_CELL:
        pos=ConvertFromABCtoXYZUnitCell(s);
        break;
      case VTK_FULL_BOX:
      default:
        pos=ConvertFromABCtoXYZ(s);
        break;
    }

    NumberOfBeadsAlreadyPlaced=0;
    NumberOfTrialPositionsForTheFirstBead=1;
    start=Components[CurrentComponent].StartingBead;
    NewPosition[0][start].x=pos.x;
    NewPosition[0][start].y=pos.y;
    NewPosition[0][start].z=pos.z;
    FirstBeadPosition.x=pos.x;
    FirstBeadPosition.y=pos.y;
    FirstBeadPosition.z=pos.z;

    NewPosition[CurrentSystem][start]=FirstBeadPosition;

    if(!BlockedPocket(FirstBeadPosition))
       value=GrowMolecule(CBMC_PARTIAL_INSERTION);
    else
       value=0.0;

    x=s.x*(REAL)SIZE_X;
    y=s.y*(REAL)SIZE_Y;
    z=s.z*(REAL)SIZE_Z;
    index=x+y*SIZE_Y+z*SIZE_X*SIZE_Y;
    if((index>0)&&(index<SIZE_X*SIZE_Y*SIZE_Z))
    {
      Histogram_3D[index]+=value;
      Histogram_Count_3D[index]+=1.0;
    }

    if(i%PrintEvery==0)
    {
      mkdir("VTK",S_IRWXU);
      sprintf(buffer,"VTK/System_%d",CurrentSystem);
      mkdir(buffer,S_IRWXU);
      sprintf(buffer,"VTK/System_%d/FrameworkSurface%s.vtk",CurrentSystem,FileNameAppend);

      FilePtr=fopen(buffer,"w");
      fprintf(FilePtr,"# vtk DataFile Version 1.0\n");
      fprintf(FilePtr,"CELL_PARAMETERS %lf %lf %lf %lf %lf %lf\n",
              cellProperties.ax,cellProperties.ay,cellProperties.az,
              (double)AlphaAngle[CurrentSystem]*180.0/M_PI,(double)BetaAngle[CurrentSystem]*180.0/M_PI,(double)GammaAngle[CurrentSystem]*180.0/M_PI);
      fprintf(FilePtr,"ASCII\n");
      fprintf(FilePtr,"DATASET STRUCTURED_POINTS\n");
      fprintf(FilePtr,"DIMENSIONS %d %d %d\n",SIZE_X,SIZE_Y,SIZE_Z);
      fprintf(FilePtr,"ORIGIN %lf %lf %lf\n",0.0,0.0,0.0);
      fprintf(FilePtr,"SPACING %lf %lf %lf\n",cellProperties.ax/(double)SIZE_X, cellProperties.ay/(double)SIZE_Y,cellProperties.az/(double)SIZE_Z);
      fprintf(FilePtr,"POINT_DATA %d\n",SIZE_X*SIZE_Y*SIZE_Z);
      fprintf(FilePtr,"SCALARS scalars unsigned_short\n");
      fprintf(FilePtr,"LOOKUP_TABLE default\n");


      min=1000.0;
      for(z=0;z<SIZE_Z;z++)
        for(y=0;y<SIZE_Y;y++)
          for(x=0;x<SIZE_X;x++)
          {
            index=x+y*SIZE_Y+z*SIZE_X*SIZE_Y;
            if(Histogram_Count_3D[index]>0.0&&Histogram_3D[index]>0.0)
            {
              value=-log(Histogram_3D[index]/Histogram_Count_3D[index]);
              if(value<min) min=value;
            }
          }
      fprintf(stderr, "min: %lf\n",(double)min);

      for(z=0;z<SIZE_Z;z++)
        for(y=0;y<SIZE_Y;y++)
          for(x=0;x<SIZE_X;x++)
          {
            index=x+y*SIZE_Y+z*SIZE_X*SIZE_Y;
            if((Histogram_Count_3D[index]>0.0)&&(Histogram_3D[index]>0.0))
            {
              temp=1+(int)(1000.0*(-log(Histogram_3D[index]/
                (REAL)Histogram_Count_3D[index])-min));
              if(temp<0)
                fprintf(FilePtr,"0\n");
              else if(temp>64000)
                fprintf(FilePtr,"65535\n");
              else
                fprintf(FilePtr,"%d\n",temp);
            }
            else
              fprintf(FilePtr,"%d\n",65535);

          }
      fclose(FilePtr);
    }
  }
  free(Histogram_3D);
  free(Histogram_Count_3D);
}

static int versionNumber=1;

void WriteRestartMovies(FILE *FilePtr)
{
  int i,j;
  fpos_t pos;
  REAL Check;

  fwrite(&versionNumber,sizeof(int),1,FilePtr);

  fwrite(Movies,NumberOfSystems,sizeof(int),FilePtr);
  fwrite(WriteMoviesEvery,NumberOfSystems,sizeof(int),FilePtr);

  // write the file offset of all the movie-files
  // flush the files to be sure the last part of all movies files have been written (e.g. in case of a computer-crash)
  for(i=0;i<NumberOfSystems;i++)
  {
    if(Movies[i])
    {
      fwrite(&nr_frame[i],1,sizeof(int),FilePtr);

      for(j=0;j<NumberOfComponents;j++)
      {
        if(PDBFilePtr[i][j])
        {
          fflush(PDBFilePtr[i][j]);
          fgetpos(PDBFilePtr[i][j],&pos);
        }
        fwrite(&pos,1,sizeof(fpos_t),FilePtr);
      }
      if(PDBFilePtrAll[i])
      {
        fflush(PDBFilePtrAll[i]);
        fgetpos(PDBFilePtrAll[i],&pos);
      }
      fwrite(&pos,1,sizeof(fpos_t),FilePtr);

      if(PDBFilePtrwork[i])
      {
        fflush(PDBFilePtrwork[i]);
        fgetpos(PDBFilePtrwork[i],&pos);
      }
      fwrite(&pos,1,sizeof(fpos_t),FilePtr);
    }
  }

  Check=123456789.0;
  fwrite(&Check,1,sizeof(REAL),FilePtr);
}

void ReadRestartMovies(FILE *FilePtr)
{
  int i,j;
  fpos_t pos;
  char buffer[1024];
  REAL Check;
  int readversionNumber=0;

  AllocateMovieMemory();

  fread(&readversionNumber,sizeof(int),1,FilePtr);
  if(readversionNumber > versionNumber)
  {
    fprintf(stderr,"Upgrade to last version of RASPA to read binary restart-file");
    exit(-1);
  }

  fread(Movies,NumberOfSystems,sizeof(int),FilePtr);
  fread(WriteMoviesEvery,NumberOfSystems,sizeof(int),FilePtr);

  SamplePDBMovies(ALLOCATE,-1);

  mkdir("Movies",S_IRWXU);
  for(i=0;i<NumberOfSystems;i++)
  {
    if(Movies[i])
    {
      fread(&nr_frame[i],1,sizeof(int),FilePtr);

      sprintf(buffer,"Movies/System_%d",i);
      mkdir(buffer,S_IRWXU);

      for(j=0;j<NumberOfComponents;j++)
      {
        fread(&pos,1,sizeof(fpos_t),FilePtr);
        sprintf(buffer,"Movies/System_%d/Component_%s_%d%s.pdb",
                i,
                Components[j].Name,
                j,
                FileNameAppend);

        // check if the file exist
        if((PDBFilePtr[i][j]=fopen(buffer,"r+")))
        {
          // the file exist, try to reposition to the saved state, if not possible start at the beginning of the file
          if(fsetpos(PDBFilePtr[i][j],&pos)!=0)
            rewind(PDBFilePtr[i][j]);
        }
        else // the file does not exist
          PDBFilePtr[i][j]=fopen(buffer,"a"); // create new file
      }

      sprintf(buffer,"Movies/System_%d/AllComponents%s.pdb",i,FileNameAppend);
      fread(&pos,1,sizeof(fpos_t),FilePtr);
      if((PDBFilePtrAll[i]=fopen(buffer,"r+")))
      {
        // the file exist, try to reposition to the saved state, if not possible start at the beginning of the file
        if(fsetpos(PDBFilePtrAll[i],&pos)!=0)
          rewind(PDBFilePtrAll[i]);
      }
      else
        PDBFilePtrAll[i]=fopen(buffer,"a");

      sprintf(buffer,"Movies/System_%d/Frameworks%s.pdb",i,FileNameAppend);
      fread(&pos,1,sizeof(fpos_t),FilePtr);
      if((PDBFilePtrwork[i]=fopen(buffer,"r+")))
      {
        // the file exist, try to reposition to the saved state, if not possible start at the beginning of the file
        if(fsetpos(PDBFilePtrwork[i],&pos)!=0)
          rewind(PDBFilePtrwork[i]);
      }
      else
        PDBFilePtrwork[i]=fopen(buffer,"a");
    }
  }

  fread(&Check,1,sizeof(REAL),FilePtr);
  if(fabs(Check-123456789.0)>1e-10)
  {
    fprintf(stderr, "Error in binary restart-file (ReadRestartMovies)\n");
    ContinueAfterCrash=FALSE;
  }
}

/*********************************************************************************************************
 * Name       | WriteSnapshotTinker                                                                      *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | writes a Tinker-file in 'P1' (no symmetry)                                               *
 * Parameters | -                                                                                        *
 *********************************************************************************************************/

void WriteSnapshotTinker(char *string)
{
  int i,j,k,Type;
  char buffer[256];
  FILE *FilePtr;
  int  SerialNumber;           // Atom serial number
  char AtomName[10]="      C"; // Atom Name
  int *AtomId;
  VECTOR r;
  int count;
  int A,B,C,D;
  REAL *parms;

  AtomId=(int*)calloc(NumberOfPseudoAtoms,sizeof(int));
  mkdir("Tinker",S_IRWXU);

  for(CurrentSystem=0;CurrentSystem<NumberOfSystems;CurrentSystem++)
  {
    sprintf(buffer,"Tinker/System_%d",CurrentSystem);
    mkdir(buffer,S_IRWXU);

    // write tinker xyz-format
    for(j=0;j<NumberOfPseudoAtoms;j++)
      AtomId[j]=0;
    sprintf(buffer,"Tinker/System_%d/Snapshot_%s%s.txyz",CurrentSystem,string,FileNameAppend);
    FilePtr=fopen(buffer,"w");

    count=0;
    for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
      for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
        if(PseudoAtoms[Adsorbates[CurrentSystem][i].Atoms[j].Type].PrintToPDB) count++;

    SerialNumber=1;

    fprintf(FilePtr,"%8d %s\n",Framework[CurrentSystem].TotalNumberOfAtoms+count,Framework[CurrentSystem].Name[0]);
    if(Framework[CurrentSystem].FrameworkModel!=NONE)
    {
      for(CurrentFramework=0;CurrentFramework<Framework[CurrentSystem].NumberOfFrameworks;CurrentFramework++)
      {
        for(j=0;j<Framework[CurrentSystem].NumberOfAtoms[CurrentFramework];j++)
        {
          Type=Framework[CurrentSystem].Atoms[CurrentFramework][j].Type;
          if(PseudoAtoms[Type].PrintToPDB)
          {
            r=Framework[CurrentSystem].Atoms[CurrentFramework][j].Position;
            sprintf(AtomName,"%3s",PseudoAtoms[Type].PrintToPDBName);
            fprintf(FilePtr,"%6d%3s %12.6lf %12.6lf %12.6lf %6d",
                    SerialNumber++,AtomName,(double)r.x,(double)r.y,(double)r.z,Type+1);
            for(i=0;i<Framework[CurrentSystem].Connectivity[CurrentFramework][j];i++)
              fprintf(FilePtr," %4d",
                      Framework[CurrentSystem].Neighbours[CurrentFramework][j][i]+1);
            fprintf(FilePtr,"\n");
          }
        }
      }
    }

    count=Framework[CurrentSystem].TotalNumberOfAtoms;
    for(i=0;i<NumberOfAdsorbateMolecules[CurrentSystem];i++)
    {
      Type=Adsorbates[CurrentSystem][i].Type;

      for(j=0;j<Adsorbates[CurrentSystem][i].NumberOfAtoms;j++)
      {
        if(PseudoAtoms[Adsorbates[CurrentSystem][i].Atoms[j].Type].PrintToPDB)
        {
          r=Adsorbates[CurrentSystem][i].Atoms[j].Position;
          sprintf(AtomName,"%3s",PseudoAtoms[Adsorbates[CurrentSystem][i].Atoms[j].Type].PrintToPDBName);
          fprintf(FilePtr,"%6d%3s %18.12lf %18.12lf %18.12lf %6d",
                  SerialNumber++,AtomName,(double)r.x,(double)r.y,(double)r.z,Adsorbates[CurrentSystem][i].Atoms[j].Type+1);
          for(k=0;k<Components[Type].Connectivity[j];k++)
            fprintf(FilePtr," %4d",
                    Components[Type].ConnectivityList[j][k]+1+count);
          fprintf(FilePtr,"\n");
        }
      }
      count+=Adsorbates[CurrentSystem][i].NumberOfAtoms;
    }

    fclose(FilePtr);

    sprintf(buffer,"Tinker/System_%d/RASPA%s.key",CurrentSystem,FileNameAppend);
    FilePtr=fopen(buffer,"w");
    fprintf(FilePtr,"bondunit                %18.12lf\n",1.0);
    fprintf(FilePtr,"angleunit               %18.12lf\n",SQR(M_PI/180.0));
    fprintf(FilePtr,"torsionunit             %18.12lf\n",1.0);
    fprintf(FilePtr,"imptorunit              %18.12lf\n",1.0);
    fprintf(FilePtr,"vdw-14-scale            0.0\n");
    fprintf(FilePtr,"chg-14-scale            0.0\n");
    fprintf(FilePtr,"electric                332.06\n");
    fprintf(FilePtr,"dielectric              1.0\n");
    fprintf(FilePtr,"\n");
    fprintf(FilePtr,"ewald\n");
    fprintf(FilePtr,"ewald-alpha 0.3\n");
    fprintf(FilePtr,"ewald-cutoff 9.0\n");
    fprintf(FilePtr,"\n");
    fprintf(FilePtr,"vdw-cutoff     12.0\n");
    fprintf(FilePtr,"vdw-taper      12.0\n");
    fprintf(FilePtr,"radiustype     SIGMA\n");
    fprintf(FilePtr,"\n");
    fprintf(FilePtr,"a-axis         25.0\n");
    fprintf(FilePtr,"b-axis         25.0\n");
    fprintf(FilePtr,"c-axis         25.0\n");
    fprintf(FilePtr,"alpha          90.0\n");
    fprintf(FilePtr,"beta           90.0\n");
    fprintf(FilePtr,"gamma          90.0\n");
    for(i=0;i<NumberOfPseudoAtoms;i++)
    {
      fprintf(FilePtr,"atom         %d    %s    %s    %d    %lf     %d\n",
          i+1,
          PseudoAtoms[i].ChemicalElement,
          PseudoAtoms[i].Name,
          i+1,
          PseudoAtoms[i].Mass,
          PseudoAtoms[i].Connectivity);
    }
    fprintf(FilePtr,"\n");

    for(i=0;i<NumberOfComponents;i++)
    {
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
            fprintf(FilePtr,"bond %d %d %lg %lg\n",
              Components[i].Type[A]+1,Components[i].Type[B]+1,0.5*parms[0]*ENERGY_TO_KELVIN,parms[1]);
            break;
          case MORSE_BOND:
            // p_0*[(1.0-{exp(-p_1*(r-p_2))})^2-1.0]
            // ===============================================
            // p_0/k_B [K]       force constant
            // p_1     [A^-1]    parameter
            // p_2     [A]       reference bond distance
            fprintf(FilePtr,"bond %d %d %lg %lg %lg\n",
               Components[i].Type[A]+1,Components[i].Type[B]+1,parms[0]*ENERGY_TO_KELVIN,parms[2],parms[1]);
            break;
        }
      }
    }

    for(i=0;i<NumberOfComponents;i++)
    {
      for(k=0;k<Components[i].NumberOfBends;k++)
      {
        A=Components[i].Bends[k].A;
        B=Components[i].Bends[k].B;
        C=Components[i].Bends[k].C;

        parms=Components[i].BendArguments[k];
        switch(Components[i].BendType[k])
        {
          case HARMONIC_BEND:
            // 0.5*p0*SQR(r-p1);
            // ===============================================
            // p_0/k_B [K/A^2]   force constant
            // p_1     [A]       reference bond distance
            fprintf(FilePtr,"angle %d %d %d %lg %lg 0.0 0.0\n",
              Components[i].Type[A]+1,Components[i].Type[B]+1,Components[i].Type[C]+1,0.5*parms[0]*ENERGY_TO_KELVIN,parms[1]*RAD2DEG);
            break;
        }
      }
    }

    for(i=0;i<NumberOfComponents;i++)
    {
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
            fprintf(FilePtr,"torsion %d %d %d %d %lg %lg %d %lg %lg %d %lg %lg %d\n",
              Components[i].Type[A]+1,Components[i].Type[B]+1,Components[i].Type[C]+1,Components[i].Type[D]+1,
              parms[1]*ENERGY_TO_KELVIN,
              0.0,
              1,
              parms[2]*ENERGY_TO_KELVIN,
              180.0,
              2,
              parms[3]*ENERGY_TO_KELVIN,
              0.0,
              3);
            break;
        }
      }
    }

    for(i=0;i<NumberOfPseudoAtoms;i++)
    {
      if(NumberOfPseudoAtomsType[CurrentSystem][i]>0)
      {
        for(j=i;j<NumberOfPseudoAtoms;j++)
        {
          if(NumberOfPseudoAtomsType[CurrentSystem][j]>0)
          {
            //if(!((PseudoAtoms[i].FrameworkAtom)&&(PseudoAtoms[j].FrameworkAtom)))
            {
              switch(PotentialType[i][j])
              {
                case ZERO_POTENTIAL:
                  break;
                case LENNARD_JONES:
                  // 4*p_0*((p_1/r)^12-(p_1/r)^6)
                  // ======================================================================================
                  // p_0/k_B [K]    strength parameter epsilon
                  // p_1     [A]    size parameter sigma
                  // p_2/k_B [K]    (non-zero for a shifted potential)
                  fprintf(FilePtr,"vdwpr %d %d %lg %lg\n",
                    i+1,
                    j+1,
                    (double)PotentialParms[i][j][1],
                    (double)PotentialParms[i][j][0]*ENERGY_TO_KELVIN);
                  break;

              }
            }
          }
        }
      }
    }

    for(i=0;i<NumberOfPseudoAtoms;i++)
    {
      if(NumberOfPseudoAtomsType[CurrentSystem][i]>0)
      {
        fprintf(FilePtr,"charge %d %lg\n",i+1,PseudoAtoms[i].Charge1);
      }
    }

    fclose(FilePtr);
  }
  free(AtomId);
}

