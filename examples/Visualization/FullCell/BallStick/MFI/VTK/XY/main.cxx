#include <math.h>
#include "vtkVersion.h"
#include "vtkPlaneCollection.h"
#include "vtkDoubleArray.h"
#include "vtkDataArray.h"
#include "vtkPoints.h"
#include "vtkCutter.h"
#include "vtkClipDataSet.h"
#include "vtkClipConvexPolyData.h"
#include "vtkClipClosedSurface.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkStructuredPointsReader.h"
#include "vtkStructuredPointsWriter.h"
#include "vtkPiecewiseFunction.h"
#include "vtkVolume.h"
#include "vtkVolumeProperty.h"
#include "vtkVolumeRayCastCompositeFunction.h"
#include "vtkVolumeRayCastIsosurfaceFunction.h"
#include <vtkVolumeRayCastMIPFunction.h>
#include <vtkVolumeRayCastMapper.h>
#include "vtkVolumeMapper.h"
#include "vtkVolume.h"
#include "vtkSphereSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkTubeFilter.h"
#include "vtkLineSource.h"
#include "vtkGlyph3D.h"
#include "vtkTIFFWriter.h"
#include "vtkPNGWriter.h"
#include "vtkPostScriptWriter.h"
#include "vtkWindowToImageFilter.h"
#include "vtkRenderLargeImage.h"
#include "vtkAxes.h"
#include "vtkTubeFilter.h"
#include "vtkCubeAxesActor2D.h"
#include "vtkOutlineFilter.h"
#include "vtkColorTransferFunction.h"
#include "vtkMarchingCubes.h"
#include "vtkImageData.h"
#include "vtkProperty.h"
#include "vtkProperty2D.h"
#include "vtkCamera.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkStructuredPoints.h"
#include "vtkContourFilter.h"
#include "vtkVolumeMapper.h"
#include "vtkVolume.h"
#include "vtkVolumeRayCastMapper.h"
#include "vtkVolumeRayCastCompositeFunction.h"
#include "vtkVolume.h"
#include "vtkVolumeProperty.h"
#include "vtkPiecewiseFunction.h"
#include "vtkPolyDataReader.h"
#include "vtkArrowSource.h"
#include "vtkVectorNorm.h"
#include "vtkGlyph3D.h"
#include "vtkPolyData.h"
#include "vtkTextProperty.h"
#include "vtkTextMapper.h"
#include "vtkAxisActor2D.h"
#include "vtkLookupTable.h"
#include "vtkBox.h"
#include "vtkPlane.h"
#include "vtkPlanes.h"
#include "vtkCylinder.h"
#include "vtkSphere.h"
#include "vtkTransform.h"
#include "vtkClipPolyData.h"
#include "vtkJPEGWriter.h"
#include "vtkPNGWriter.h"
#include "vtkImageCast.h"
#include "vtkImageMapToColors.h"
#include "vtkDataSetMapper.h"
#include "vtkPolyDataNormals.h"
#include "vtkSmartPointer.h"
#include "vtkCamera.h"
#include "vtkLight.h"
#include "vtkLightKit.h"
#include "vtkCleanPolyData.h"
#include "vtkExtractEdges.h"
#include "vtkPolyLine.h"
#include "vtkCellArray.h"
#include "vtkPointData.h"
#include "vtkDoubleArray.h"

typedef struct point
{
  double x;
  double y;
  double z;
} VECTOR;

int ReadInputFile(void);

// Copyright © 2009,2010,2011 David Dubbeldam
// David Dubbeldam: version 1.0  31 October 2009 D.Dubbeldam@uva.nl
// David Dubbeldam: version 1.01 23 November 2009 D.Dubbeldam@uva.nl
// David Dubbeldam: version 1.1  17 September 2010 D.Dubbeldam@uva.nl
// David Dubbeldam: version 1.2  11 March 2012 D.Dubbeldam@uva.nl
// David Dubbeldam: version 1.3  23 March 2012 D.Dubbeldam@uva.nl

// Changes:
//
// 23 Nov 2009, version 1.01 (corrected a bug for non-orthorhombic unit cells, added isocontours)
// ==============================================================================================
// the problem is that very high values are at clipped to 65535 (2^16), the grid is the rectabgular box enclosing
// the non-orthorhombic unit cell. The outside value are also put at 65535, but because an open pore has very low value
// the interpolation scheme within VTK generates intermnediate values which are visibly rendered, i.e. there appears a
// very thin layer at the boundaries of the unit cell. Putting the outside values at 0 and also overlap inside the box 
// to 0 leads to "ugly" pictures (the interpolation problem at the walls of the pore). Example of the problem: DDR-zeolite.
// The chosen solution is to generate the whole orthorhombic VTK grid and then to clip the unit cell so that the original
// non-orthorhombic unit cell is recovered. Note that the unit cell properties (length, angles) MUST be defined.
// In addition to volume-rendering, one can now also use isocontouring using various materials (like glass, metal, etc).

// 31 March 2011, version 1.3
// ==============================================================================================
// Modified the coordinates, the input vtk-files are now in Cartesian positions, except for the surface- and density grids.
// Because they can not be scaled properly (the colors and rendering changes) the Cartesian positions of atoms are
// scaled and translated to map in the positional system of the grids (e.g. 150x150x150). Note that in the grid-files
// the ASPECT_RATIO is used. Triclinic cells still have rectangular grids but everything outside the cell is clipped away.
// New:
// 1) Frameworks-bonds are now drawn when they extend out of the cell (they are clipped at the cell-boundary).
// 2) Framework-atoms outside the cell are clipped, atoms at the boundary do have sphere-glyphs that extends out of the cell.
// 3) You can clip adsorbates/cation atoms. Both atom and bonds are clipped and the surfaces filled. Note: the VTK algorithm
//    often fails, but when projected orthogonal it is not noticable.
// 4) You can draw the cell-frame for the full cell, or for all the individual unit-cells. You can specify exactly which
//    unit cells to include.
// 5) You can highlight a specific volume, with a cylinder, a rectangular channel, or a rotated rectangular channel (e.g. diamond-shaped).
//    The high-lighted volume is for example a channel, the low-lighted volume can be made completely transparant.
//    Note: does not work for volume-rendering of the framework (VTK limitation)

enum{VOLUME_RENDERING,ISOSURFACE};
enum{GLASS,BRUSHED_METAL,METALLIC_PASTEL,TRANSPARENT,RASPA};
enum{OFF=0,ON=1};
enum{NONE,X,Y,Z,XY,XZ,YZ,XYZ};
enum{CYLINDER,RECTANGULAR_CYLINDER,SPHERE,CUBE};
enum{NO_FRAME,FULL_CELL,UNIT_CELL};

// sets a name and address in the lower-left corner using a grey color
// replace by your own name and address
int Credits=OFF;
char credits[256]={"David Dubbeldam, University of Amsterdam\nVan 't Hoff Institute for Molecular Sciences\nD.Dubbeldam@uva.nl"};

int Title=OFF;
char title[256]={"MgMOF-74 (single unitcell); XY-view"};
char OutputFileName[256]="./MgMOF-74-XY-view-nC6.jpg";

int CutOutShape=OFF;           // Cut-out a shape (e.g. cylinder) to highlight everything inside the volume
int CutOutDirection=Z;         // The orientation of the cut-out shape
double CutOutRotationAngle=0;  // the angle of rotation (around the axis of the recangular channel, e.g. to make a diamond-shape)
int CutOutType=CYLINDER;       // the type of the shape: CYLINDER, RECTANGULAR_CYLINDER, SPHERE
double CutOutRadius1=8.0;      // the radius of the cylinder, sphere
double CutOutRadius2=8.0;      // second radius, used to make rectangular channels
VECTOR CutOutFractionalShift={0.5,0.5,0.5}; // the center of the shape in fractional positions

// picture quality settings
int Resolution=10;    // the resolution of spheres and tubes, the higher the more smooth (use 10, but 50 for the final picture)
int AA=1;                       // anti-aliasing, use 1, but 8 for final picture
double SampleDistance=0.1;      // sample spacing, use 1, but 0.1 for final picture
double ImageSampleDistance=0.5; // rays per pixel, use 1, but 0.5 for final picture (4 rays per pixel)

// write a serie of pictures while rotating the structure
int Movie=OFF;

// control the properties of framework-atoms
int FrameworkAtoms=ON;                                   // draw the framework atoms
int FrameworkCutOutShape=OFF;                            // whether or not cutting applies to the framework
char FrameworkAtomsFilename[256]="./FrameworkAtoms.vtk"; // framework atoms file
char FrameworkBondsFilename[256]="./FrameworkBonds.vtk"; // framework bonding file
double FrameworkOpacity=1.0;                             // set the opacity of the atoms
double FrameworkBondThickness=1.0;                       // set the thickness of the bonds

char FrameworkSurfaceFilename[256]="./FrameworkSurface.vtk";
int FrameworkSurface=OFF;                 // by default do not render the surface
int FrameworkRenderingMethod=ISOSURFACE;  // the rendering method of the framework-surface: ISOSURFACE, VOLUME_RENDERING
int IsoSurfaceMaterial=GLASS;             // set the type of isosuface-material: GLASS,BRUSHED_METAL,METALLIC_PASTEL,TRANSPARENT,RASPA
double IsoSurfaceOpacity=1.0;             // set the opacity of the usual/low-lighted volume
int HighLightedIsoSurfaceMaterial=GLASS;           // set the type of isosuface-material of the highlighted volume
VECTOR HighLightedIsoSurfaceColor={0.5, 0.75,1.0}; // set the color of the highlighted surface
double HighLightedIsoSurfaceOpacity=0.35;          // set the opacity of the highlighted surface

int Density=OFF;                           // by default do not render the density
char DensityFilename[256]="./Density.vtk";

// control the properties of adsorbates
int Adsorbates=ON;
char AdsorbateAtomsFilename[256]="./AdsorbateAtoms.vtk";
double AdsorbateOpacity=1.0;
double HighLightedAdsorbateOpacity=1.0;
int ClipAdsorbates=NONE;

// control the properties of cations
int Cations=ON;
char CationAtomsFilename[256]="./CationAtoms.vtk";
double CationOpacity=1.0;
double HighLightedCationOpacity=1.0;
int ClipCations=NONE;

// zoom in or out by increasing/decreasing the zoom-factor
double ZoomFactor=1.0;

// scale the size of the atoms and bonds
double ScaleFactor=0.175;

// control the view-point of the oject (input in degrees)
double Azimuth=0.0;
double Elevation=0.0;
double Roll=0.0;

// control the view-point of the oject (input in degrees)
double AdditionalAzimuth=0.0;
double AdditionalElevation=0.0;
double AdditionalRoll=0.0;

VECTOR VisibilityFractionSurfaceMax={1.0,1.0,1.0};
VECTOR VisibilityFractionAdsorbateMax={1.0,1.0,1.0};
VECTOR VisibilityFractionCationMax={1.0,1.0,1.0};

// control of the cell-frame
int Frame=FULL_CELL;
double FrameOpacity=1.0;
double FrameThickness=1.0;
int NrDuplicatesX=1;
int NrDuplicatesY=1;
int NrDuplicatesZ=1;
int FrameMin[3]={0,0,0};
int FrameMax[3]={0,0,0};

// control of the cell-axes (works only for rectangular cells [VTK limitation])
int Axes=OFF;
double AxesFontFactor=1.0;
double AxesLabelFactor=1.0;
VECTOR AxesColor={0.0,0.0,0.0};
VECTOR AxesLabelColor={0.5,0.2,0.2};
int NumberOfAxesLabels=4;
int SetXAxisVisibility=ON;
int SetYAxisVisibility=ON;
int SetZAxisVisibility=ON;

// the size of the image in pixels
int ImageSizeX=500;
int ImageSizeY=500;

// the lengths of the cell-vectors
double A=51.75300;
double B=51.75300;
double C=27.14240;

// the angles between the cell-vectors
double AlphaAngle=90.0;
double BetaAngle=90.0;
double GammaAngle=90.0;

// the size of the framework suface- and density grids.
int SIZE_X=150-1;
int SIZE_Y=150-1;
int SIZE_Z=150-1;

#define SQR(x) ((x)*(x))
#define MAX2(x,y) (((x)>(y))?(x):(y))

double epsilon=1e-3;

double magnification=1;

typedef struct real_matrix3x3
{
  double ax;
  double ay;
  double az;

  double bx;
  double by;
  double bz;

  double cx;
  double cy;
  double cz;

} double_MATRIX3x3;

vtkActor* IsoContourArray;
vtkVolume* VolumeArray;
vtkVolume* DensityArray;
vtkActor* FrameArray[1024];

vtkDataSetMapper *FrameworkSurfaceHighLightMapper;
vtkDataSetMapper *FrameworkSurfaceLowLightMapper;
vtkClipDataSet *CutOutLowLightedShape;
vtkClipDataSet *CutFullCell;
vtkMarchingCubes *FrameworkSurfaceIsocontour;
vtkPolyDataNormals *FrameworkSurfaceNormals;

bool FileExists(const char *filename)
{
  if (FILE * file = fopen(filename, "r"))
  {
    fclose(file);
    return true;
  }
  return false;
}

static double x[8][3]={{0,0,0}, {1,0,0}, {1,1,0}, {0,1,0},
                       {0,0,1}, {1,0,1}, {1,1,1}, {0,1,1}};
static vtkIdType pts[6][5]={{0,1,2,3,0}, {4,5,6,7,4}, {0,1,5,4,0},
                      {1,2,6,5,1}, {2,3,7,6,2}, {3,0,4,7,3}};

int main( int argc, char *argv[] )
{
  int i,j,k,index;
  double tempd;
  double a,b,c;
  double_MATRIX3x3 CellBox,UnitCellBox;
  VECTOR Size,shift,edge_vector,center_vector,pos;
  double max;
  VECTOR v1,v2;
  VECTOR plane_normal1,plane_normal2,plane_normal3;
  double length;
  char buffer[256];

  ReadInputFile();

  AlphaAngle*=M_PI/180.0;
  BetaAngle*=M_PI/180.0;
  GammaAngle*=M_PI/180.0;

  tempd=(cos(AlphaAngle)-cos(GammaAngle)*cos(BetaAngle))/sin(GammaAngle);
  a=A/NrDuplicatesX; b=B/NrDuplicatesY; c=C/NrDuplicatesZ;
  UnitCellBox.ax=a;   UnitCellBox.bx=b*cos(GammaAngle); UnitCellBox.cx=c*cos(BetaAngle);
  UnitCellBox.ay=0.0; UnitCellBox.by=b*sin(GammaAngle); UnitCellBox.cy=c*tempd;
  UnitCellBox.az=0.0; UnitCellBox.bz=0.0;               UnitCellBox.cz=c*sqrt(1.0-SQR(cos(BetaAngle))-SQR(tempd));
  CellBox.ax=A;       CellBox.bx=B*cos(GammaAngle);     CellBox.cx=C*cos(BetaAngle);
  CellBox.ay=0.0;     CellBox.by=B*sin(GammaAngle);     CellBox.cy=C*tempd;
  CellBox.az=0.0;     CellBox.bz=0.0;                   CellBox.cz=C*sqrt(1.0-SQR(cos(BetaAngle))-SQR(tempd));

  Size.x=Size.y=Size.z=0.0;
  shift.x=shift.y=shift.z=0.0;
  edge_vector.x=1.0;
  edge_vector.y=0.0;
  edge_vector.z=0.0;
  pos.x=CellBox.ax*edge_vector.x+CellBox.bx*edge_vector.y+CellBox.cx*edge_vector.z;
  pos.y=CellBox.ay*edge_vector.x+CellBox.by*edge_vector.y+CellBox.cy*edge_vector.z;
  pos.z=CellBox.az*edge_vector.x+CellBox.bz*edge_vector.y+CellBox.cz*edge_vector.z;
  Size.x+=fabs(pos.x);
  Size.y+=fabs(pos.y);
  Size.z+=fabs(pos.z);
  if(pos.x<0.0) shift.x+=pos.x;
  if(pos.y<0.0) shift.y+=pos.y;
  if(pos.z<0.0) shift.z+=pos.z;

  edge_vector.x=0.0;
  edge_vector.y=1.0;
  edge_vector.z=0.0;
  pos.x=CellBox.ax*edge_vector.x+CellBox.bx*edge_vector.y+CellBox.cx*edge_vector.z;
  pos.y=CellBox.ay*edge_vector.x+CellBox.by*edge_vector.y+CellBox.cy*edge_vector.z;
  pos.z=CellBox.az*edge_vector.x+CellBox.bz*edge_vector.y+CellBox.cz*edge_vector.z;
  Size.x+=fabs(pos.x);
  Size.y+=fabs(pos.y);
  Size.z+=fabs(pos.z);
  if(pos.x<0.0) shift.x+=pos.x;
  if(pos.y<0.0) shift.y+=pos.y;
  if(pos.z<0.0) shift.z+=pos.z;

  edge_vector.x=0.0;
  edge_vector.y=0.0;
  edge_vector.z=1.0;
  pos.x=CellBox.ax*edge_vector.x+CellBox.bx*edge_vector.y+CellBox.cx*edge_vector.z;
  pos.y=CellBox.ay*edge_vector.x+CellBox.by*edge_vector.y+CellBox.cy*edge_vector.z;
  pos.z=CellBox.az*edge_vector.x+CellBox.bz*edge_vector.y+CellBox.cz*edge_vector.z;
  Size.x+=fabs(pos.x);
  Size.y+=fabs(pos.y);
  Size.z+=fabs(pos.z);
  if(pos.x<0.0) shift.x+=pos.x;
  if(pos.y<0.0) shift.y+=pos.y;
  if(pos.z<0.0) shift.z+=pos.z;

  max=MAX2(Size.x,MAX2(Size.y,Size.z));

  pos.x=CellBox.ax*CutOutFractionalShift.x+CellBox.bx*CutOutFractionalShift.y+CellBox.cx*CutOutFractionalShift.z;
  pos.y=CellBox.ay*CutOutFractionalShift.x+CellBox.by*CutOutFractionalShift.y+CellBox.cy*CutOutFractionalShift.z;
  pos.z=CellBox.az*CutOutFractionalShift.x+CellBox.bz*CutOutFractionalShift.y+CellBox.cz*CutOutFractionalShift.z;
  center_vector.x=(pos.x-shift.x)*(double)SIZE_X/max;
  center_vector.y=(pos.y-shift.y)*(double)SIZE_Y/max;
  center_vector.z=(pos.z-shift.z)*(double)SIZE_Z/max;

  vtkTransform *transform_from_box_to_diamand_framework=vtkTransform::New();
    transform_from_box_to_diamand_framework->RotateY(CutOutRotationAngle);
  switch(CutOutDirection)
  {
    case X:
      transform_from_box_to_diamand_framework->RotateZ(90);
      break;
    case Y:
      break;
    case Z:
      transform_from_box_to_diamand_framework->RotateX(90);
      break;
  }
  transform_from_box_to_diamand_framework->Translate(-center_vector.x,-center_vector.y,-center_vector.z);

  // create a cylinder with default orientation in the y-direction
  vtkCylinder *CylinderFramework=vtkCylinder::New();
      CylinderFramework->SetRadius(CutOutRadius1*SIZE_X/max);
      CylinderFramework->SetTransform(transform_from_box_to_diamand_framework);

  // create a cylinder with default orientation in the y-direction
  vtkBox *RectangularCylinderFramework=vtkBox::New();
    RectangularCylinderFramework->SetBounds(-CutOutRadius1*SIZE_X/max,CutOutRadius1*SIZE_X/max,-1000,1000,-CutOutRadius2*SIZE_X/max,CutOutRadius2*SIZE_X/max);
    RectangularCylinderFramework->SetTransform(transform_from_box_to_diamand_framework);

  // create a sphere 
  vtkSphere *SphereFramework=vtkSphere::New();
    SphereFramework->SetRadius(CutOutRadius1*SIZE_X/max);
    SphereFramework->SetTransform(transform_from_box_to_diamand_framework);

  // framework-, adsorbate-, and cation-atoms
  pos.x=CellBox.ax*CutOutFractionalShift.x+CellBox.bx*CutOutFractionalShift.y+CellBox.cx*CutOutFractionalShift.z;
  pos.y=CellBox.ay*CutOutFractionalShift.x+CellBox.by*CutOutFractionalShift.y+CellBox.cy*CutOutFractionalShift.z;
  pos.z=CellBox.az*CutOutFractionalShift.x+CellBox.bz*CutOutFractionalShift.y+CellBox.cz*CutOutFractionalShift.z;
  center_vector.x=pos.x;
  center_vector.y=pos.y;
  center_vector.z=pos.z;

  vtkTransform *transform_from_box_to_diamand=vtkTransform::New();
    transform_from_box_to_diamand->RotateY(CutOutRotationAngle);
  switch(CutOutDirection)
  {
    case X:
      transform_from_box_to_diamand->RotateZ(90);
      break;
    case Y:
      break;
    case Z:
      transform_from_box_to_diamand->RotateX(90);
      break;
  }
  transform_from_box_to_diamand->Translate(-center_vector.x,-center_vector.y,-center_vector.z);

  // create a cylinder with default orientation in the y-direction
  vtkCylinder *Cylinder=vtkCylinder::New();
      Cylinder->SetRadius(CutOutRadius1);
      Cylinder->SetTransform(transform_from_box_to_diamand);

  // create a cylinder with default orientation in the y-direction
  vtkBox *RectangularCylinder=vtkBox::New();
    RectangularCylinder->SetBounds(-CutOutRadius1,CutOutRadius1,-1000,1000,-CutOutRadius2,CutOutRadius2);
    RectangularCylinder->SetTransform(transform_from_box_to_diamand);

  // create a sphere 
  vtkSphere *Sphere=vtkSphere::New();
    Sphere->SetRadius(CutOutRadius1);
    Sphere->SetTransform(transform_from_box_to_diamand);

  // first plane normal
  edge_vector.x=0.0;
  edge_vector.y=1.0;
  edge_vector.z=0.0;
  v1.x=CellBox.ax*edge_vector.x+CellBox.bx*edge_vector.y+CellBox.cx*edge_vector.z;
  v1.y=CellBox.ay*edge_vector.x+CellBox.by*edge_vector.y+CellBox.cy*edge_vector.z;
  v1.z=CellBox.az*edge_vector.x+CellBox.bz*edge_vector.y+CellBox.cz*edge_vector.z;

  edge_vector.x=0.0;
  edge_vector.y=0.0;
  edge_vector.z=1.0;
  v2.x=CellBox.ax*edge_vector.x+CellBox.bx*edge_vector.y+CellBox.cx*edge_vector.z;
  v2.y=CellBox.ay*edge_vector.x+CellBox.by*edge_vector.y+CellBox.cy*edge_vector.z;
  v2.z=CellBox.az*edge_vector.x+CellBox.bz*edge_vector.y+CellBox.cz*edge_vector.z;

  plane_normal1.x=v1.y*v2.z-v1.z*v2.y;
  plane_normal1.y=v1.z*v2.x-v1.x*v2.z;
  plane_normal1.z=v1.x*v2.y-v1.y*v2.x;

  length=sqrt(SQR(plane_normal1.x)+SQR(plane_normal1.y)+SQR(plane_normal1.z));
  plane_normal1.x/=length;
  plane_normal1.y/=length;
  plane_normal1.z/=length;

  // second plane normal
  edge_vector.x=0.0;
  edge_vector.y=0.0;
  edge_vector.z=1.0;
  v1.x=CellBox.ax*edge_vector.x+CellBox.bx*edge_vector.y+CellBox.cx*edge_vector.z;
  v1.y=CellBox.ay*edge_vector.x+CellBox.by*edge_vector.y+CellBox.cy*edge_vector.z;
  v1.z=CellBox.az*edge_vector.x+CellBox.bz*edge_vector.y+CellBox.cz*edge_vector.z;

  edge_vector.x=1.0;
  edge_vector.y=0.0;
  edge_vector.z=0.0;
  v2.x=CellBox.ax*edge_vector.x+CellBox.bx*edge_vector.y+CellBox.cx*edge_vector.z;
  v2.y=CellBox.ay*edge_vector.x+CellBox.by*edge_vector.y+CellBox.cy*edge_vector.z;
  v2.z=CellBox.az*edge_vector.x+CellBox.bz*edge_vector.y+CellBox.cz*edge_vector.z;

  plane_normal2.x=v1.y*v2.z-v1.z*v2.y;
  plane_normal2.y=v1.z*v2.x-v1.x*v2.z;
  plane_normal2.z=v1.x*v2.y-v1.y*v2.x;

  length=sqrt(SQR(plane_normal2.x)+SQR(plane_normal2.y)+SQR(plane_normal2.z));
  plane_normal2.x/=length;
  plane_normal2.y/=length;
  plane_normal2.z/=length;

  // third plane normal
  edge_vector.x=1.0;
  edge_vector.y=0.0;
  edge_vector.z=0.0;
  v1.x=CellBox.ax*edge_vector.x+CellBox.bx*edge_vector.y+CellBox.cx*edge_vector.z;
  v1.y=CellBox.ay*edge_vector.x+CellBox.by*edge_vector.y+CellBox.cy*edge_vector.z;
  v1.z=CellBox.az*edge_vector.x+CellBox.bz*edge_vector.y+CellBox.cz*edge_vector.z;
  
  edge_vector.x=0.0;
  edge_vector.y=1.0;
  edge_vector.z=0.0;
  v2.x=CellBox.ax*edge_vector.x+CellBox.bx*edge_vector.y+CellBox.cx*edge_vector.z;
  v2.y=CellBox.ay*edge_vector.x+CellBox.by*edge_vector.y+CellBox.cy*edge_vector.z;
  v2.z=CellBox.az*edge_vector.x+CellBox.bz*edge_vector.y+CellBox.cz*edge_vector.z;

  plane_normal3.x=v1.y*v2.z-v1.z*v2.y;
  plane_normal3.y=v1.z*v2.x-v1.x*v2.z;
  plane_normal3.z=v1.x*v2.y-v1.y*v2.x;

  length=sqrt(SQR(plane_normal3.x)+SQR(plane_normal3.y)+SQR(plane_normal3.z));
  plane_normal3.x/=length;
  plane_normal3.y/=length;
  plane_normal3.z/=length;


  // For the framework, the clipping is done in terms of 150x150x150 times 150/max(length)
  // The clipping is done _after_ the grid has been translated
  vtkPlane *p1=vtkPlane::New();
     p1->SetOrigin(0,0,0);
     p1->SetNormal(plane_normal1.x,plane_normal1.y,plane_normal1.z);
  vtkPlane *p2=vtkPlane::New();
     p2->SetOrigin(VisibilityFractionSurfaceMax.x*CellBox.ax*SIZE_X/max,VisibilityFractionSurfaceMax.x*CellBox.ay*SIZE_Y/max,VisibilityFractionSurfaceMax.x*CellBox.az*SIZE_Z/max);
     p2->SetNormal(-plane_normal1.x,-plane_normal1.y,-plane_normal1.z);
  vtkPlane *p3=vtkPlane::New();
     p3->SetOrigin(0,0,0);
     p3->SetNormal(plane_normal2.x,plane_normal2.y,plane_normal2.z);
  vtkPlane *p4=vtkPlane::New();
     p4->SetOrigin(VisibilityFractionSurfaceMax.y*CellBox.bx*SIZE_X/max,VisibilityFractionSurfaceMax.y*CellBox.by*SIZE_Y/max,VisibilityFractionSurfaceMax.z*CellBox.bz*SIZE_Z/max);
     p4->SetNormal(-plane_normal2.x,-plane_normal2.y,-plane_normal2.z);
  vtkPlane *p5=vtkPlane::New();
     p5->SetOrigin(0,0,0);
     p5->SetNormal(plane_normal3.x,plane_normal3.y,plane_normal3.z);
  vtkPlane *p6=vtkPlane::New();
     p6->SetOrigin(VisibilityFractionSurfaceMax.z*CellBox.cx*SIZE_X/max,VisibilityFractionSurfaceMax.z*CellBox.cy*SIZE_Y/max,VisibilityFractionSurfaceMax.z*CellBox.cz*SIZE_Z/max);
     p6->SetNormal(-plane_normal3.x,-plane_normal3.y,-plane_normal3.z);

  vtkPlaneCollection *BoxClipPlaneCollection=vtkPlaneCollection::New();
   BoxClipPlaneCollection->AddItem(p1);
   BoxClipPlaneCollection->AddItem(p2);
   BoxClipPlaneCollection->AddItem(p3);
   BoxClipPlaneCollection->AddItem(p4);
   BoxClipPlaneCollection->AddItem(p5);
   BoxClipPlaneCollection->AddItem(p6);

  // the isocontour-clipping vtkClipPolyData expects an implicit function: vtkPlanes
  vtkPlanes *BoxClipPlanes=vtkPlanes::New();

  vtkPoints *points=vtkPoints::New();
  vtkDoubleArray *norms=vtkDoubleArray::New();
    norms->SetNumberOfComponents(3);

  points->InsertPoint(0,(-shift.x)*(double)SIZE_X/max,
                        (-shift.y)*(double)SIZE_Y/max,
                        (-shift.z)*(double)SIZE_Z/max);
  norms->InsertTuple3(0,-plane_normal1.x,-plane_normal1.y,-plane_normal1.z);

  points->InsertPoint(1,(CellBox.ax-shift.x)*(double)SIZE_X/max,
                        (CellBox.ay-shift.y)*(double)SIZE_Y/max,
                        (CellBox.az-shift.z)*(double)SIZE_Z/max);
  norms->InsertTuple3(1,plane_normal1.x,plane_normal1.y,plane_normal1.z);

  points->InsertPoint(2,(-shift.x)*(double)SIZE_X/max,
                        (-shift.y)*(double)SIZE_Y/max,
                        (-shift.z)*(double)SIZE_Z/max);
  norms->InsertTuple3(2,-plane_normal2.x,-plane_normal2.y,-plane_normal2.z);

  points->InsertPoint(3,(CellBox.bx-shift.x)*(double)SIZE_X/max,
                        (CellBox.by-shift.y)*(double)SIZE_Y/max,
                        (CellBox.bz-shift.z)*(double)SIZE_Z/max);
  norms->InsertTuple3(3,plane_normal2.x,plane_normal2.y,plane_normal2.z);


  points->InsertPoint(4,(-shift.x)*(double)SIZE_X/max,
                        (-shift.y)*(double)SIZE_Y/max,
                        (-shift.z)*(double)SIZE_Z/max);
  norms->InsertTuple3(4,-plane_normal3.x,-plane_normal3.y,-plane_normal3.z);

  points->InsertPoint(5,(CellBox.cx-shift.x)*(double)SIZE_X/max,
                        (CellBox.cy-shift.y)*(double)SIZE_Y/max,
                        (CellBox.cz-shift.z)*(double)SIZE_Z/max);
  norms->InsertTuple3(5,plane_normal3.x,plane_normal3.y,plane_normal3.z);

  BoxClipPlanes->SetPoints(points);
  BoxClipPlanes->SetNormals(norms);


  // the frameworkbond-clipping vtkClipPolyData expects an implicit function: vtkPlanes
  vtkPlanes *FrameworkBondClipPlanes=vtkPlanes::New();

  vtkPoints *points2=vtkPoints::New();
  vtkDoubleArray *norms2=vtkDoubleArray::New();
    norms2->SetNumberOfComponents(3);

  points2->InsertPoint(0,-epsilon,-epsilon,-epsilon);
  norms2->InsertTuple3(0,-plane_normal1.x,-plane_normal1.y,-plane_normal1.z);

  points2->InsertPoint(1,CellBox.ax+epsilon,CellBox.ay+epsilon,CellBox.az+epsilon);
  norms2->InsertTuple3(1,plane_normal1.x,plane_normal1.y,plane_normal1.z);

  points2->InsertPoint(2,-epsilon,-epsilon,-epsilon);
  norms2->InsertTuple3(2,-plane_normal2.x,-plane_normal2.y,-plane_normal2.z);

  points2->InsertPoint(3,CellBox.bx+epsilon,CellBox.by+epsilon,CellBox.bz+epsilon);
  norms2->InsertTuple3(3,plane_normal2.x,plane_normal2.y,plane_normal2.z);

  points2->InsertPoint(4,-epsilon,-epsilon,-epsilon);
  norms2->InsertTuple3(4,-plane_normal3.x,-plane_normal3.y,-plane_normal3.z);

  points2->InsertPoint(5,CellBox.cx+epsilon,CellBox.cy+epsilon,CellBox.cz+epsilon);
  norms2->InsertTuple3(5,plane_normal3.x,plane_normal3.y,plane_normal3.z);

  FrameworkBondClipPlanes->SetPoints(points2);
  FrameworkBondClipPlanes->SetNormals(norms2);

  // the clipping of adsorbates and cations is done with surface-closing.
  // the vtkClipClosedSurface expects a vtkPlaneCollection
  vtkPlane *q1=vtkPlane::New();
     q1->SetOrigin(0,0,0);
     q1->SetNormal(plane_normal1.x,plane_normal1.y,plane_normal1.z);
  vtkPlane *q2=vtkPlane::New();
     q2->SetOrigin(VisibilityFractionAdsorbateMax.x*CellBox.ax,VisibilityFractionAdsorbateMax.x*CellBox.ay,VisibilityFractionAdsorbateMax.x*CellBox.az);
     q2->SetNormal(-plane_normal1.x,-plane_normal1.y,-plane_normal1.z);
  vtkPlane *q3=vtkPlane::New();
     q3->SetOrigin(0,0,0);
     q3->SetNormal(plane_normal2.x,plane_normal2.y,plane_normal2.z);
  vtkPlane *q4=vtkPlane::New();
     q4->SetOrigin(VisibilityFractionAdsorbateMax.y*CellBox.bx,VisibilityFractionAdsorbateMax.y*CellBox.by,VisibilityFractionAdsorbateMax.y*CellBox.bz);
     q4->SetNormal(-plane_normal2.x,-plane_normal2.y,-plane_normal2.z);
  vtkPlane *q5=vtkPlane::New();
     q5->SetOrigin(0,0,0);
     q5->SetNormal(plane_normal3.x,plane_normal3.y,plane_normal3.z);
  vtkPlane *q6=vtkPlane::New();
     q6->SetOrigin(VisibilityFractionAdsorbateMax.z*CellBox.cx,VisibilityFractionAdsorbateMax.z*CellBox.cy,VisibilityFractionAdsorbateMax.z*CellBox.cz);
     q6->SetNormal(-plane_normal3.x,-plane_normal3.y,-plane_normal3.z);

  vtkPlane *q7=vtkPlane::New();
     q7->SetOrigin(0,0,0);
     q7->SetNormal(plane_normal1.x,plane_normal1.y,plane_normal1.z);
  vtkPlane *q8=vtkPlane::New();
     q8->SetOrigin(VisibilityFractionCationMax.x*CellBox.ax,VisibilityFractionCationMax.x*CellBox.ay,VisibilityFractionCationMax.x*CellBox.az);
     q8->SetNormal(-plane_normal1.x,-plane_normal1.y,-plane_normal1.z);
  vtkPlane *q9=vtkPlane::New();
     q9->SetOrigin(0,0,0);
     q9->SetNormal(plane_normal2.x,plane_normal2.y,plane_normal2.z);
  vtkPlane *q10=vtkPlane::New();
     q10->SetOrigin(VisibilityFractionCationMax.y*CellBox.bx,VisibilityFractionCationMax.y*CellBox.by,VisibilityFractionCationMax.y*CellBox.bz);
     q10->SetNormal(-plane_normal2.x,-plane_normal2.y,-plane_normal2.z);
  vtkPlane *q11=vtkPlane::New();
     q11->SetOrigin(0,0,0);
     q11->SetNormal(plane_normal3.x,plane_normal3.y,plane_normal3.z);
  vtkPlane *q12=vtkPlane::New();
     q12->SetOrigin(VisibilityFractionCationMax.z*CellBox.cx,VisibilityFractionCationMax.z*CellBox.cy,VisibilityFractionCationMax.z*CellBox.cz);
     q12->SetNormal(-plane_normal3.x,-plane_normal3.y,-plane_normal3.z);

  vtkPlaneCollection *PlanesCollectionAdsorbate=vtkPlaneCollection::New();
  switch(ClipAdsorbates)
  {
    case NONE:
      break;
    case X:
      PlanesCollectionAdsorbate->AddItem(q1); PlanesCollectionAdsorbate->AddItem(q2);
      break;
    case Y:
      PlanesCollectionAdsorbate->AddItem(q3); PlanesCollectionAdsorbate->AddItem(q4);
      break;
    case Z:
      PlanesCollectionAdsorbate->AddItem(q5); PlanesCollectionAdsorbate->AddItem(q6);
      break;
    case XY:
      PlanesCollectionAdsorbate->AddItem(q1); PlanesCollectionAdsorbate->AddItem(q2);
      PlanesCollectionAdsorbate->AddItem(q3); PlanesCollectionAdsorbate->AddItem(q4);
      break;
    case XZ:
      PlanesCollectionAdsorbate->AddItem(q1); PlanesCollectionAdsorbate->AddItem(q2);
      PlanesCollectionAdsorbate->AddItem(q5); PlanesCollectionAdsorbate->AddItem(q6);
      break;
    case YZ:
      PlanesCollectionAdsorbate->AddItem(q3); PlanesCollectionAdsorbate->AddItem(q4);
      PlanesCollectionAdsorbate->AddItem(q5); PlanesCollectionAdsorbate->AddItem(q6);
      break;
    case XYZ:
      PlanesCollectionAdsorbate->AddItem(q1); PlanesCollectionAdsorbate->AddItem(q2);
      PlanesCollectionAdsorbate->AddItem(q3); PlanesCollectionAdsorbate->AddItem(q4);
      PlanesCollectionAdsorbate->AddItem(q5); PlanesCollectionAdsorbate->AddItem(q6);
      break;
  }

  vtkPlaneCollection *PlanesCollectionCation=vtkPlaneCollection::New();
  switch(ClipCations)
  {
    case NONE:
      break;
    case X:
      PlanesCollectionCation->AddItem(q7); PlanesCollectionCation->AddItem(q8);
      break;
    case Y:
      PlanesCollectionCation->AddItem(q9); PlanesCollectionCation->AddItem(q10);
      break;
    case Z:
      PlanesCollectionCation->AddItem(q11); PlanesCollectionCation->AddItem(q12);
      break;
    case XY:
      PlanesCollectionCation->AddItem(q7); PlanesCollectionCation->AddItem(q8);
      PlanesCollectionCation->AddItem(q9); PlanesCollectionCation->AddItem(q10);
      break;
    case XZ:
      PlanesCollectionCation->AddItem(q7); PlanesCollectionCation->AddItem(q8);
      PlanesCollectionCation->AddItem(q11); PlanesCollectionCation->AddItem(q12);
      break;
    case YZ:
      PlanesCollectionCation->AddItem(q9); PlanesCollectionCation->AddItem(q10);
      PlanesCollectionCation->AddItem(q11); PlanesCollectionCation->AddItem(q12);
      break;
    case XYZ:
      PlanesCollectionCation->AddItem(q7); PlanesCollectionCation->AddItem(q8);
      PlanesCollectionCation->AddItem(q9); PlanesCollectionCation->AddItem(q10);
      PlanesCollectionCation->AddItem(q11); PlanesCollectionCation->AddItem(q12);
      break;
  }

  // Create the renderer, render window, and interactor
  vtkRenderer *ren1 = vtkRenderer::New();
  vtkRenderWindow *renWin = vtkRenderWindow::New();
    renWin->AddRenderer(ren1);
    renWin->SetSize(600,600);
  vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
    iren->SetRenderWindow(renWin);

    // Create a transfer function mapping scalar value to opacity
    // play around with this tp highlight different things
    // with the current settings the adsorpion sites are clearly visable
    vtkPiecewiseFunction *oTFunSurface = vtkPiecewiseFunction::New();
      oTFunSurface->AddSegment(0, 0.0, 1, 0.0);
      oTFunSurface->AddSegment(1, 0.0, 10000, 0.0);
      oTFunSurface->AddSegment(10000, 0.00, 15000, 0.6);
      oTFunSurface->AddSegment(20000, 0.6, 30000, 0.9);
      oTFunSurface->AddSegment(30000, 0.9, 60000, 1.0);
      oTFunSurface->AddSegment(60000, 1.0, 65000, 1.0);
      oTFunSurface->AddSegment(65000, 0.0, 65536, 0.0);

    vtkPiecewiseFunction *oTFunDensity = vtkPiecewiseFunction::New();
      oTFunDensity->AddSegment(0, 0.0, 1, 0.0);
      oTFunDensity->AddSegment(1, 0.0, 10000, 0.05);
      oTFunDensity->AddSegment(10000, 0.05, 15000, 0.4);
      oTFunDensity->AddSegment(15000, 0.4, 20000, 0.8);
      oTFunDensity->AddSegment(20000, 0.8, 30000, 1.0);
      oTFunDensity->AddSegment(30000, 1.0, 100000, 1.0);


    // Create a transfer function mapping scalar value to color (grey)
    vtkPiecewiseFunction *cTFun = vtkPiecewiseFunction::New();
      cTFun->AddSegment(0, 1.0, 255, 1.0);

    // create a mapping from values to color
    // after many trials, this function looks okay
    vtkColorTransferFunction *cTFunN = vtkColorTransferFunction::New();
      cTFunN->AddRGBPoint(0000.0,0.0,0.0,0.0);
      cTFunN->AddRGBPoint(0001.0,0.5,0.2,0.2);
      cTFunN->AddRGBPoint(5000+1000.0,1.0,0.3,0.3);
      cTFunN->AddRGBPoint(10000+7000.0,1.0,0.6,0.3);
      cTFunN->AddRGBPoint(10000+18000.0,1.0,1.0,0.0);
      cTFunN->AddRGBPoint(10000+19000.0,0.4,1.0,0.6);
      cTFunN->AddRGBPoint(10000+20000.0,0.3,0.7,0.5);
      cTFunN->AddRGBPoint(10000+50000.0,0.3,0.5,1.0);

    // create a mapping from values to color
    // after many trials, this function looks okay
    vtkColorTransferFunction *cTFunN2 = vtkColorTransferFunction::New();
      cTFunN2->AddRGBPoint(0000.0,0.0,0.0,0.0);
      cTFunN2->AddRGBPoint(0001.0,0.5,0.2,0.2);
      cTFunN2->AddRGBPoint(500+2000.0,1.0,0.3,0.3);
      cTFunN2->AddRGBPoint(1000+5000.0,0.1,1.0,0.6);
      cTFunN2->AddRGBPoint(1000+10000.0,0.1,0.7,0.5);
      cTFunN2->AddRGBPoint(1000+25000.0,0.1,0.5,1.0);
      cTFunN2->AddRGBPoint(1000+50000.0,0.0,0.0,1.0);

    vtkPiecewiseFunction *oTFunDensity2 = vtkPiecewiseFunction::New();
      oTFunDensity2->AddSegment(0, 0.0, 1, 0.0);
      oTFunDensity2->AddSegment(1, 0.0, 10000, 0.0);
      oTFunDensity2->AddSegment(10000, 0.0, 15000, 0.4);
      oTFunDensity2->AddSegment(15000, 0.4, 20000, 0.8);
      oTFunDensity2->AddSegment(20000, 0.8, 30000, 1.0);
      oTFunDensity2->AddSegment(30000, 1.0, 100000, 1.0);


  // define a color table to map atoms to colors
  vtkLookupTable *LookupTable = vtkLookupTable::New();
    LookupTable->SetNumberOfColors(33);
    LookupTable->Build();
    LookupTable->SetTableValue( 0,0.0, 0.0, 1.0, 1.0); // blue
    LookupTable->SetTableValue( 1,1.0, 0.0, 0.0, 1.0); // red
    LookupTable->SetTableValue( 2,0.35,0.35,0.35,1.0); // gray
    LookupTable->SetTableValue( 3,1.0, 0.5, 0.0, 1.0); // orange
    LookupTable->SetTableValue( 4,1.0, 1.0, 0.0, 1.0); // yellow
    LookupTable->SetTableValue( 5,0.5, 0.5, 0.2, 1.0); // tan
    LookupTable->SetTableValue( 6,0.6, 0.6, 0.6, 1.0); // silver
    LookupTable->SetTableValue( 7,0.0, 1.0, 0.0, 1.0); // green
    LookupTable->SetTableValue( 8,1.0, 1.0, 1.0, 1.0); // white
    LookupTable->SetTableValue( 9,1.0, 0.6, 0.6, 1.0); // pink
    LookupTable->SetTableValue(10,0.25,0.75,0.75,1.0); // cyan
    LookupTable->SetTableValue(11,0.65,0.0, 0.65,1.0); // purple
    LookupTable->SetTableValue(12,0.5, 0.9, 0.4, 1.0); // lime
    LookupTable->SetTableValue(13,0.9, 0.4, 0.7, 1.0); // mauvre
    LookupTable->SetTableValue(14,0.5, 0.3, 0.0, 1.0); // ochre
    LookupTable->SetTableValue(15,0.5, 0.5, 0.75,1.0); // iceblue
    LookupTable->SetTableValue(16,0.0, 0.0, 0.0, 1.0); // black
    LookupTable->SetTableValue(17,0.88,0.97,0.02,1.0); // yellow2
    LookupTable->SetTableValue(18,0.55,0.9 ,0.02,1.0); // yellow3
    LookupTable->SetTableValue(19,0.0, 0.9 ,0.04,1.0); // green2
    LookupTable->SetTableValue(20,0.0, 0.9 ,0.5 ,1.0); // green3
    LookupTable->SetTableValue(21,0.0, 0.88,1.0 ,1.0); // cyan2
    LookupTable->SetTableValue(22,0.0, 0.76,1.0 ,1.0); // cyan3
    LookupTable->SetTableValue(23,0.02,0.38,0.67,1.0); // blue2
    LookupTable->SetTableValue(24,0.01,0.04,0.93,1.0); // blue3
    LookupTable->SetTableValue(25,0.27,0.0, 0.98,1.0); // violet
    LookupTable->SetTableValue(26,0.45,0.0, 0.9, 1.0); // violet2
    LookupTable->SetTableValue(27,0.9 ,0.0, 0.9, 1.0); // magenta
    LookupTable->SetTableValue(28,1.0 ,0.0, 0.66,1.0); // magenta2
    LookupTable->SetTableValue(29,0.98,0.0, 0.23,1.0); // red2
    LookupTable->SetTableValue(30,0.81,0.0, 0.0 ,1.0); // red3
    LookupTable->SetTableValue(31,0.89,0.35,0.0 ,1.0); // orange2
    LookupTable->SetTableValue(32,0.96,0.72,0.0 ,1.0); // orange2

  vtkPolyData *PolyData=vtkPolyData::New();
  vtkPoints *PointData=vtkPoints::New();
  vtkCellArray *polys = vtkCellArray::New();

  if(Frame==UNIT_CELL)
  {
    // make unitcell-frame
    for(i=0;i<8;i++)
    {
      pos.x=UnitCellBox.ax*x[i][0]+UnitCellBox.bx*x[i][1]+UnitCellBox.cx*x[i][2];
      pos.y=UnitCellBox.ay*x[i][0]+UnitCellBox.by*x[i][1]+UnitCellBox.cy*x[i][2];
      pos.z=UnitCellBox.az*x[i][0]+UnitCellBox.bz*x[i][1]+UnitCellBox.cz*x[i][2];
      PointData->InsertPoint(i,pos.x*SIZE_X/max,pos.y*SIZE_Y/max,pos.z*SIZE_Z/max);
    }

    // generate polyline of the frame
    for(i=0;i<6;i++) polys->InsertNextCell(5,pts[i]);

    PolyData->SetPoints(PointData);
    PolyData->SetLines(polys);

    vtkTubeFilter *TubesFrame=vtkTubeFilter::New();

    #if(VTK_MAJOR_VERSION<6)
     TubesFrame->SetInput(PolyData);
    #else
      TubesFrame->SetInputData(PolyData);
    #endif
      TubesFrame->SetRadius(0.3*FrameThickness);
      TubesFrame->SetNumberOfSides(3*Resolution);
    vtkPolyDataMapper *TubeMapperFrame=vtkPolyDataMapper::New();
      TubeMapperFrame->SetInputConnection(TubesFrame->GetOutputPort());

    index=0;
    for(i=FrameMin[0];i<=FrameMax[0];i++)
      for(j=FrameMin[1];j<=FrameMax[1];j++)
        for(k=FrameMin[2];k<=FrameMax[2];k++)
        {
          pos.x=UnitCellBox.ax*i+UnitCellBox.bx*j+UnitCellBox.cx*k;
          pos.y=UnitCellBox.ay*i+UnitCellBox.by*j+UnitCellBox.cy*k;
          pos.z=UnitCellBox.az*i+UnitCellBox.bz*j+UnitCellBox.cz*k;

          FrameArray[index]=vtkActor::New();
          FrameArray[index]->SetMapper(TubeMapperFrame);
          FrameArray[index]->GetProperty()->SetColor(1,1,1);
          FrameArray[index]->GetProperty()->SetOpacity(FrameOpacity);
          FrameArray[index]->SetPosition(pos.x*SIZE_X/max,pos.y*SIZE_Y/max,pos.z*SIZE_Z/max);
          ren1->AddViewProp(FrameArray[index]);
          index++;
        }
  }
  else if(Frame==FULL_CELL)
  {
    // make unitcell-frame

    // {0,0,0}
    pos.x=UnitCellBox.ax*(FrameMin[0]+x[0][0])+UnitCellBox.bx*(FrameMin[1]+x[0][1])+UnitCellBox.cx*(FrameMin[2]+x[0][2]);
    pos.y=UnitCellBox.ay*(FrameMin[0]+x[0][0])+UnitCellBox.by*(FrameMin[1]+x[0][1])+UnitCellBox.cy*(FrameMin[2]+x[0][2]);
    pos.z=UnitCellBox.az*(FrameMin[0]+x[0][0])+UnitCellBox.bz*(FrameMin[1]+x[0][1])+UnitCellBox.cz*(FrameMin[2]+x[0][2]);
    PointData->InsertPoint(0,pos.x*SIZE_X/max,pos.y*SIZE_Y/max,pos.z*SIZE_Z/max);

    // {1,0,0}
    pos.x=UnitCellBox.ax*(FrameMax[0]+x[1][0])+UnitCellBox.bx*(FrameMin[1]+x[1][1])+UnitCellBox.cx*(FrameMin[2]+x[1][2]);
    pos.y=UnitCellBox.ay*(FrameMax[0]+x[1][0])+UnitCellBox.by*(FrameMin[1]+x[1][1])+UnitCellBox.cy*(FrameMin[2]+x[1][2]);
    pos.z=UnitCellBox.az*(FrameMax[0]+x[1][0])+UnitCellBox.bz*(FrameMin[1]+x[1][1])+UnitCellBox.cz*(FrameMin[2]+x[1][2]);
    PointData->InsertPoint(1,pos.x*SIZE_X/max,pos.y*SIZE_Y/max,pos.z*SIZE_Z/max);

    // {1,1,0}
    pos.x=UnitCellBox.ax*(FrameMax[0]+x[2][0])+UnitCellBox.bx*(FrameMax[1]+x[2][1])+UnitCellBox.cx*(FrameMin[2]+x[2][2]);
    pos.y=UnitCellBox.ay*(FrameMax[0]+x[2][0])+UnitCellBox.by*(FrameMax[1]+x[2][1])+UnitCellBox.cy*(FrameMin[2]+x[2][2]);
    pos.z=UnitCellBox.az*(FrameMax[0]+x[2][0])+UnitCellBox.bz*(FrameMax[1]+x[2][1])+UnitCellBox.cz*(FrameMin[2]+x[2][2]);
    PointData->InsertPoint(2,pos.x*SIZE_X/max,pos.y*SIZE_Y/max,pos.z*SIZE_Z/max);

    // {0,1,0}
    pos.x=UnitCellBox.ax*(FrameMin[0]+x[3][0])+UnitCellBox.bx*(FrameMax[1]+x[3][1])+UnitCellBox.cx*(FrameMin[2]+x[3][2]);
    pos.y=UnitCellBox.ay*(FrameMin[0]+x[3][0])+UnitCellBox.by*(FrameMax[1]+x[3][1])+UnitCellBox.cy*(FrameMin[2]+x[3][2]);
    pos.z=UnitCellBox.az*(FrameMin[0]+x[3][0])+UnitCellBox.bz*(FrameMax[1]+x[3][1])+UnitCellBox.cz*(FrameMin[2]+x[3][2]);
    PointData->InsertPoint(3,pos.x*SIZE_X/max,pos.y*SIZE_Y/max,pos.z*SIZE_Z/max);

    // {0,0,1}
    pos.x=UnitCellBox.ax*(FrameMin[0]+x[4][0])+UnitCellBox.bx*(FrameMin[1]+x[4][1])+UnitCellBox.cx*(FrameMax[2]+x[4][2]);
    pos.y=UnitCellBox.ay*(FrameMin[0]+x[4][0])+UnitCellBox.by*(FrameMin[1]+x[4][1])+UnitCellBox.cy*(FrameMax[2]+x[4][2]);
    pos.z=UnitCellBox.az*(FrameMin[0]+x[4][0])+UnitCellBox.bz*(FrameMin[1]+x[4][1])+UnitCellBox.cz*(FrameMax[2]+x[4][2]);
    PointData->InsertPoint(4,pos.x*SIZE_X/max,pos.y*SIZE_Y/max,pos.z*SIZE_Z/max);

    // {1,0,1}
    pos.x=UnitCellBox.ax*(FrameMax[0]+x[5][0])+UnitCellBox.bx*(FrameMin[1]+x[5][1])+UnitCellBox.cx*(FrameMax[2]+x[5][2]);
    pos.y=UnitCellBox.ay*(FrameMax[0]+x[5][0])+UnitCellBox.by*(FrameMin[1]+x[5][1])+UnitCellBox.cy*(FrameMax[2]+x[5][2]);
    pos.z=UnitCellBox.az*(FrameMax[0]+x[5][0])+UnitCellBox.bz*(FrameMin[1]+x[5][1])+UnitCellBox.cz*(FrameMax[2]+x[5][2]);
    PointData->InsertPoint(5,pos.x*SIZE_X/max,pos.y*SIZE_Y/max,pos.z*SIZE_Z/max);

    // {1,1,1}
    pos.x=UnitCellBox.ax*(FrameMax[0]+x[6][0])+UnitCellBox.bx*(FrameMax[1]+x[6][1])+UnitCellBox.cx*(FrameMax[2]+x[6][2]);
    pos.y=UnitCellBox.ay*(FrameMax[0]+x[6][0])+UnitCellBox.by*(FrameMax[1]+x[6][1])+UnitCellBox.cy*(FrameMax[2]+x[6][2]);
    pos.z=UnitCellBox.az*(FrameMax[0]+x[6][0])+UnitCellBox.bz*(FrameMax[1]+x[6][1])+UnitCellBox.cz*(FrameMax[2]+x[6][2]);
    PointData->InsertPoint(6,pos.x*SIZE_X/max,pos.y*SIZE_Y/max,pos.z*SIZE_Z/max);

    // {0,1,1}
    pos.x=UnitCellBox.ax*(FrameMin[0]+x[7][0])+UnitCellBox.bx*(FrameMax[1]+x[7][1])+UnitCellBox.cx*(FrameMax[2]+x[7][2]);
    pos.y=UnitCellBox.ay*(FrameMin[0]+x[7][0])+UnitCellBox.by*(FrameMax[1]+x[7][1])+UnitCellBox.cy*(FrameMax[2]+x[7][2]);
    pos.z=UnitCellBox.az*(FrameMin[0]+x[7][0])+UnitCellBox.bz*(FrameMax[1]+x[7][1])+UnitCellBox.cz*(FrameMax[2]+x[7][2]);
    PointData->InsertPoint(7,pos.x*SIZE_X/max,pos.y*SIZE_Y/max,pos.z*SIZE_Z/max);

    // generate polyline of the frame
    for(i=0;i<6;i++) polys->InsertNextCell(5,pts[i]);

    PolyData->SetPoints(PointData);
    PolyData->SetLines(polys);

    vtkTubeFilter *TubesFrame=vtkTubeFilter::New();
      #if(VTK_MAJOR_VERSION<6)
        TubesFrame->SetInput(PolyData);
      #else
        TubesFrame->SetInputData(PolyData);
      #endif
      TubesFrame->SetRadius(0.3*FrameThickness);
      TubesFrame->SetNumberOfSides(3*Resolution);
    vtkPolyDataMapper *TubeMapperFrame=vtkPolyDataMapper::New();
      TubeMapperFrame->SetInputConnection(TubesFrame->GetOutputPort());

    FrameArray[0]=vtkActor::New();
    FrameArray[0]->SetMapper(TubeMapperFrame);
    FrameArray[0]->GetProperty()->SetColor(1,1,1);
    FrameArray[0]->GetProperty()->SetOpacity(FrameOpacity);
    FrameArray[0]->SetPosition(0,0,0);
    ren1->AddViewProp(FrameArray[0]);


    vtkTextProperty *multiLineTextPropAxes=vtkTextProperty::New();
      multiLineTextPropAxes->SetFontSize(200);
      multiLineTextPropAxes->SetFontFamilyToCourier();
      multiLineTextPropAxes->BoldOn();
      multiLineTextPropAxes->ItalicOn();
      multiLineTextPropAxes->ShadowOn();
      multiLineTextPropAxes->SetLineSpacing(1.5);
      multiLineTextPropAxes->SetJustificationToLeft();
      multiLineTextPropAxes->SetVerticalJustificationToCentered();
      multiLineTextPropAxes->SetColor(AxesLabelColor.x,AxesLabelColor.y,AxesLabelColor.z);

    vtkCubeAxesActor2D *axes=vtkCubeAxesActor2D::New();
      #if(VTK_MAJOR_VERSION<6)
        axes->SetInput(TubesFrame->GetOutput());
      #else
        axes->SetInputConnection(TubesFrame->GetOutputPort());
      #endif
      axes->SetCamera(ren1->GetActiveCamera());
      axes->SetLabelFormat("%4.3f");
      axes->SetFlyModeToOuterEdges();
      axes->SetRanges(0,1.0,0,1.0,0,1.0);
      axes->SetUseRanges(1);
      axes->SetScaling(1);
      axes->SetFontFactor(AxesFontFactor);
      axes->SetNumberOfLabels(1+NumberOfAxesLabels);
      axes->SetCornerOffset(0.0);
      axes->SetXAxisVisibility(SetXAxisVisibility);
      axes->SetYAxisVisibility(SetYAxisVisibility);
      axes->SetZAxisVisibility(SetZAxisVisibility);
      axes->SetXLabel("X");
      axes->SetYLabel("Y");
      axes->SetZLabel("Z");
      axes->SetAxisLabelTextProperty(multiLineTextPropAxes);
      axes->SetAxisTitleTextProperty(multiLineTextPropAxes);
      axes->GetXAxisActor2D()->SetTickOffset(2);
      axes->GetXAxisActor2D()->SetTickLength(5);
      axes->GetXAxisActor2D()->SetTitleTextProperty(multiLineTextPropAxes);
      axes->GetXAxisActor2D()->SetLabelTextProperty(multiLineTextPropAxes);
      axes->GetXAxisActor2D()->SetLabelFactor(AxesLabelFactor);
      axes->GetXAxisActor2D()->TickVisibilityOn();
      axes->GetYAxisActor2D()->SetTickOffset(2);
      axes->GetYAxisActor2D()->SetTickLength(5);
      axes->GetYAxisActor2D()->SetTitleTextProperty(multiLineTextPropAxes);
      axes->GetYAxisActor2D()->SetLabelTextProperty(multiLineTextPropAxes);
      axes->GetYAxisActor2D()->SetLabelFactor(AxesLabelFactor);
      axes->GetYAxisActor2D()->TickVisibilityOn();
      axes->GetZAxisActor2D()->SetTickOffset(2);
      axes->GetZAxisActor2D()->SetTickLength(5);
      axes->GetZAxisActor2D()->SetTitleTextProperty(multiLineTextPropAxes);
      axes->GetZAxisActor2D()->SetLabelTextProperty(multiLineTextPropAxes);
      axes->GetZAxisActor2D()->SetLabelFactor(AxesLabelFactor);
      axes->GetZAxisActor2D()->TickVisibilityOn();
      axes->GetProperty()->SetColor(AxesColor.x,AxesColor.y,AxesColor.z);
    if(Axes==ON)
      ren1->AddViewProp(axes);
  }

  // render the atoms as spheres
  vtkSphereSource *Atom=vtkSphereSource::New();
    Atom->SetThetaResolution(Resolution);
    Atom->SetPhiResolution(Resolution);
    Atom->SetRadius(2.0*ScaleFactor);

  if(FileExists(FrameworkAtomsFilename)&&(FrameworkAtoms==ON))
  {
    // read in the vtk-files for the framework atoms
    vtkPolyDataReader *ReaderFrameworkAtoms=vtkPolyDataReader::New();
      ReaderFrameworkAtoms->SetFileName(FrameworkAtomsFilename);

    // render the framework atoms as spheres
    vtkGlyph3D *GlyphFramework=vtkGlyph3D::New();
      GlyphFramework->SetInputConnection(ReaderFrameworkAtoms->GetOutputPort());
      GlyphFramework->SetSourceConnection(Atom->GetOutputPort());
      GlyphFramework->SetColorModeToColorByVector();
      GlyphFramework->SetScaleModeToScaleByScalar();

    vtkClipPolyData *CutOutFrameworkAtoms=vtkClipPolyData::New();
      CutOutFrameworkAtoms->SetInputConnection(GlyphFramework->GetOutputPort());
      switch(CutOutType)
      {
        case RECTANGULAR_CYLINDER:
          CutOutFrameworkAtoms->SetClipFunction(RectangularCylinder);
          break;
        case CYLINDER:
          CutOutFrameworkAtoms->SetClipFunction(Cylinder);
          break;
        case SPHERE:
          CutOutFrameworkAtoms->SetClipFunction(Sphere);
          break;
      }
      CutOutFrameworkAtoms->InsideOutOn();

    vtkPolyDataMapper *MapperFrameworkAtoms=vtkPolyDataMapper::New();
      if(FrameworkCutOutShape==ON)
        MapperFrameworkAtoms->SetInputConnection(CutOutFrameworkAtoms->GetOutputPort());
      else
        MapperFrameworkAtoms->SetInputConnection(GlyphFramework->GetOutputPort());
      MapperFrameworkAtoms->SetLookupTable(LookupTable);

    vtkActor *AtomActorFrameworkAtoms=vtkActor::New();
      AtomActorFrameworkAtoms->SetMapper(MapperFrameworkAtoms);
      AtomActorFrameworkAtoms->SetPosition(0,0,0);
      AtomActorFrameworkAtoms->SetScale(SIZE_X/max,SIZE_Y/max,SIZE_Z/max);
      AtomActorFrameworkAtoms->GetProperty()->SetOpacity(FrameworkOpacity);
    ren1->AddActor(AtomActorFrameworkAtoms);
  }
  else
   printf("Skipping rendering of the framework atoms (file '%s' does not exist or 'FrameworkAtoms=OFF')\n",FrameworkAtomsFilename);

  if(FileExists(FrameworkBondsFilename)&&(FrameworkAtoms==ON))
  {
    // read in the vtk-files for the framework atoms
    vtkPolyDataReader *ReaderFrameworkBonds=vtkPolyDataReader::New();
      ReaderFrameworkBonds->SetFileName(FrameworkBondsFilename);

    vtkClipPolyData *ClipSurface=vtkClipPolyData::New();
       ClipSurface->SetClipFunction(FrameworkBondClipPlanes);
       ClipSurface->SetInputConnection(ReaderFrameworkBonds->GetOutputPort());
       ClipSurface->InsideOutOn();

    // render the bonds as white tubes
    vtkTubeFilter *BondsTubeFramework=vtkTubeFilter::New();
      BondsTubeFramework->SetInputConnection(ClipSurface->GetOutputPort());
      BondsTubeFramework->SetRadius(0.125*FrameworkBondThickness);
      BondsTubeFramework->CappingOn();
      BondsTubeFramework->SetNumberOfSides(3*Resolution);

    vtkClipPolyData *CutOutFrameworkBonds=vtkClipPolyData::New();
      CutOutFrameworkBonds->SetInputConnection(BondsTubeFramework->GetOutputPort());
      switch(CutOutType)
      {
        case RECTANGULAR_CYLINDER:
          CutOutFrameworkBonds->SetClipFunction(RectangularCylinder);
          break;
        case CYLINDER:
          CutOutFrameworkBonds->SetClipFunction(Cylinder);
          break;
        case SPHERE:
          CutOutFrameworkBonds->SetClipFunction(Sphere);
          break;
      }
      CutOutFrameworkBonds->InsideOutOn();

    vtkPolyDataMapper *BondsMapperFramework=vtkPolyDataMapper::New();
      if(FrameworkCutOutShape==ON)
        BondsMapperFramework->SetInputConnection(CutOutFrameworkBonds->GetOutputPort());
      else
        BondsMapperFramework->SetInputConnection(BondsTubeFramework->GetOutputPort());
      BondsMapperFramework->ScalarVisibilityOff();

    vtkActor *BondsActorFramework=vtkActor::New();
      BondsActorFramework->SetMapper(BondsMapperFramework);
      BondsActorFramework->GetProperty()->SetColor(0.9,0.9,0.9);
      BondsActorFramework->GetProperty()->SetOpacity(FrameworkOpacity);
      BondsActorFramework->SetPosition(0,0,0);
      BondsActorFramework->SetScale(SIZE_X/max,SIZE_Y/max,SIZE_Z/max);
    ren1->AddActor(BondsActorFramework);
  }
  else
   printf("Skipping rendering of the framework bonds (file '%s' does not exist or 'FrameworkAtoms=OFF')\n",FrameworkBondsFilename);

  if(FileExists(AdsorbateAtomsFilename)&&(Adsorbates==ON))
  {
    // render atom of the adsorbates
    // ======================================================================================================================

    // read in the vtk-files for the adsorbate atoms
    vtkPolyDataReader *AtomReaderAdsorbate=vtkPolyDataReader::New();
      AtomReaderAdsorbate->SetFileName(AdsorbateAtomsFilename);
      AtomReaderAdsorbate->Update();

    // render the atoms as spheres
    vtkGlyph3D *GlyphAdsorbate=vtkGlyph3D::New();
      GlyphAdsorbate->SetInputConnection(AtomReaderAdsorbate->GetOutputPort());
      GlyphAdsorbate->SetSourceConnection(Atom->GetOutputPort());
      GlyphAdsorbate->SetColorModeToColorByVector();
      GlyphAdsorbate->SetScaleModeToScaleByScalar();

    // clip the adsorbates atoms to the unit-cell and polish the cuts with surfaces
    vtkClipClosedSurface *ClipSurfaceAdsorbateAtoms=vtkClipClosedSurface::New();
       ClipSurfaceAdsorbateAtoms->SetClippingPlanes(PlanesCollectionAdsorbate);
       ClipSurfaceAdsorbateAtoms->SetInputConnection(GlyphAdsorbate->GetOutputPort());
       ClipSurfaceAdsorbateAtoms->PassPointDataOn();
       ClipSurfaceAdsorbateAtoms->TriangulationErrorDisplayOn();

    // handle the highlighted atoms
    if(CutOutShape==ON)
    {
      // the highlighted atoms are the ones that are inside the chosen clipping region
      vtkClipPolyData *CutOutHighLightedShapeAdsorbateAtoms=vtkClipPolyData::New();
        CutOutHighLightedShapeAdsorbateAtoms->SetInputConnection(ClipSurfaceAdsorbateAtoms->GetOutputPort());
        switch(CutOutType)
        {
          case RECTANGULAR_CYLINDER:
            CutOutHighLightedShapeAdsorbateAtoms->SetClipFunction(RectangularCylinder);
            break;
          case CYLINDER:
            CutOutHighLightedShapeAdsorbateAtoms->SetClipFunction(Cylinder);
            break;
          case SPHERE:
            CutOutHighLightedShapeAdsorbateAtoms->SetClipFunction(Sphere);
            break;
        }
        CutOutHighLightedShapeAdsorbateAtoms->InsideOutOn();

      vtkPolyDataMapper *MapperHighLightedAdsorbateAtoms=vtkPolyDataMapper::New();
        MapperHighLightedAdsorbateAtoms->SetInputConnection(CutOutHighLightedShapeAdsorbateAtoms->GetOutputPort());
        MapperHighLightedAdsorbateAtoms->SetLookupTable(LookupTable);

      vtkActor *ActorHighLightedAdsorbateAtoms=vtkActor::New();
        ActorHighLightedAdsorbateAtoms->SetMapper(MapperHighLightedAdsorbateAtoms);
        ActorHighLightedAdsorbateAtoms->SetPosition(0,0,0);
        ActorHighLightedAdsorbateAtoms->SetScale(SIZE_X/max,SIZE_Y/max,SIZE_Z/max);
        ActorHighLightedAdsorbateAtoms->GetProperty()->SetOpacity(HighLightedAdsorbateOpacity);
      ren1->AddActor(ActorHighLightedAdsorbateAtoms);
    }

    // handle the lowlighted, regular, atoms
    // the lowlighted atoms are the ones that are outside the chosen clipping region, or all atoms when there is no highlighted volume
    vtkClipPolyData *CutOutLowLightedShapeAdsorbateAtoms=vtkClipPolyData::New();
      CutOutLowLightedShapeAdsorbateAtoms->SetInputConnection(ClipSurfaceAdsorbateAtoms->GetOutputPort());
      switch(CutOutType)
      {
        case RECTANGULAR_CYLINDER:
          CutOutLowLightedShapeAdsorbateAtoms->SetClipFunction(RectangularCylinder);
          break;
        case CYLINDER:
          CutOutLowLightedShapeAdsorbateAtoms->SetClipFunction(Cylinder);
          break;
        case SPHERE:
          CutOutLowLightedShapeAdsorbateAtoms->SetClipFunction(Sphere);
          break;
      }
      CutOutLowLightedShapeAdsorbateAtoms->InsideOutOff();

    vtkPolyDataMapper *MapperLowLightedAdsorbateAtoms=vtkPolyDataMapper::New();
      if(CutOutShape==ON)
        MapperLowLightedAdsorbateAtoms->SetInputConnection(CutOutLowLightedShapeAdsorbateAtoms->GetOutputPort());
      else
        MapperLowLightedAdsorbateAtoms->SetInputConnection(ClipSurfaceAdsorbateAtoms->GetOutputPort());
      MapperLowLightedAdsorbateAtoms->SetLookupTable(LookupTable);

    vtkActor *ActorLowLightedAdsorbateAtoms=vtkActor::New();
      ActorLowLightedAdsorbateAtoms->SetMapper(MapperLowLightedAdsorbateAtoms);
      ActorLowLightedAdsorbateAtoms->SetPosition(0,0,0);
      ActorLowLightedAdsorbateAtoms->SetScale(SIZE_X/max,SIZE_Y/max,SIZE_Z/max);
      ActorLowLightedAdsorbateAtoms->GetProperty()->SetOpacity(AdsorbateOpacity);
    ren1->AddActor(ActorLowLightedAdsorbateAtoms);

    // render bonds of the adsorbates
    // ======================================================================================================================

    // render the bonds as white tubes
    vtkTubeFilter *BondsTubeAdsorbate=vtkTubeFilter::New();
      BondsTubeAdsorbate->SetInputConnection(AtomReaderAdsorbate->GetOutputPort());
      BondsTubeAdsorbate->SetRadius(1.0*ScaleFactor);
      BondsTubeAdsorbate->CappingOn();
      BondsTubeAdsorbate->SetNumberOfSides(Resolution);

    // clip the adsorbates bonds to the unit-cell and polish the cuts with surfaces
    vtkClipClosedSurface *ClipSurfaceAdsorbates=vtkClipClosedSurface::New();
       ClipSurfaceAdsorbates->SetClippingPlanes(PlanesCollectionAdsorbate);
       ClipSurfaceAdsorbates->SetInputConnection(BondsTubeAdsorbate->GetOutputPort());
       ClipSurfaceAdsorbates->TriangulationErrorDisplayOn();

    // handle the highlighted bonds
    if(CutOutShape==ON)
    {
        // the highlighted bonds are the ones that are inside the chosen clipping region
        vtkClipPolyData *CutOutHighLightedShapeAdsorbateBonds=vtkClipPolyData::New();
          CutOutHighLightedShapeAdsorbateBonds->SetInputConnection(ClipSurfaceAdsorbates->GetOutputPort());
          switch(CutOutType)
          {
            case RECTANGULAR_CYLINDER:
              CutOutHighLightedShapeAdsorbateBonds->SetClipFunction(RectangularCylinder);
              break;
            case CYLINDER:
              CutOutHighLightedShapeAdsorbateBonds->SetClipFunction(Cylinder);
              break;
            case SPHERE:
              CutOutHighLightedShapeAdsorbateBonds->SetClipFunction(Sphere);
              break;
          }
          CutOutHighLightedShapeAdsorbateBonds->InsideOutOn();

        vtkPolyDataMapper *BondsMapperAdsorbate=vtkPolyDataMapper::New();
          BondsMapperAdsorbate->SetInputConnection(CutOutHighLightedShapeAdsorbateBonds->GetOutputPort());
          BondsMapperAdsorbate->ScalarVisibilityOff();
        vtkActor *BondsActorAdsorbate=vtkActor::New();
          BondsActorAdsorbate->SetMapper(BondsMapperAdsorbate);
          BondsActorAdsorbate->GetProperty()->SetColor(0.9,0.9,0.9);
          BondsActorAdsorbate->GetProperty()->SetOpacity(HighLightedAdsorbateOpacity);
          BondsActorAdsorbate->SetPosition(0,0,0);
          BondsActorAdsorbate->SetScale(SIZE_X/max,SIZE_Y/max,SIZE_Z/max);
        ren1->AddActor(BondsActorAdsorbate);
    }

    // handle the lowlighted, regular, bonds
    // the lowlighted bonds are the ones that are outside the chosen clipping region, or all bonds when there is no highlighted volume
    vtkClipPolyData *CutOutLowLightedShapeAdsorbateBonds=vtkClipPolyData::New();
      CutOutLowLightedShapeAdsorbateBonds->SetInputConnection(ClipSurfaceAdsorbates->GetOutputPort());
      switch(CutOutType)
      {
        case RECTANGULAR_CYLINDER:
          CutOutLowLightedShapeAdsorbateBonds->SetClipFunction(RectangularCylinder);
          break;
        case CYLINDER:
          CutOutLowLightedShapeAdsorbateBonds->SetClipFunction(Cylinder);
          break;
        case SPHERE:
          CutOutLowLightedShapeAdsorbateBonds->SetClipFunction(Sphere);
          break;
      }
      CutOutLowLightedShapeAdsorbateBonds->InsideOutOff();

    vtkPolyDataMapper *BondsMapperAdsorbates=vtkPolyDataMapper::New();
      if(CutOutShape==ON)
        BondsMapperAdsorbates->SetInputConnection(CutOutLowLightedShapeAdsorbateBonds->GetOutputPort());
      else 
        BondsMapperAdsorbates->SetInputConnection(ClipSurfaceAdsorbates->GetOutputPort());
      BondsMapperAdsorbates->ScalarVisibilityOff();

    vtkActor *BondsActorAdsorbates=vtkActor::New();
      BondsActorAdsorbates->SetMapper(BondsMapperAdsorbates);
      BondsActorAdsorbates->GetProperty()->SetColor(0.9,0.9,0.9);
      BondsActorAdsorbates->GetProperty()->SetOpacity(AdsorbateOpacity);
      BondsActorAdsorbates->SetPosition(0,0,0);
      BondsActorAdsorbates->SetScale(SIZE_X/max,SIZE_Y/max,SIZE_Z/max);
    ren1->AddActor(BondsActorAdsorbates);
  }
  else
   printf("Skipping rendering of the adsorbates (file '%s' does not exist)\n",AdsorbateAtomsFilename);

  if(FileExists(CationAtomsFilename)&&(Cations==ON))
  {
    // render atom of the cations
    // ======================================================================================================================

    // read in the vtk-files for the cations atoms
    vtkPolyDataReader *AtomReaderCations=vtkPolyDataReader::New();
      AtomReaderCations->SetFileName(CationAtomsFilename);
      AtomReaderCations->Update();

    // render the atoms as spheres
    vtkGlyph3D *GlyphCations=vtkGlyph3D::New();
      GlyphCations->SetInputConnection(AtomReaderCations->GetOutputPort());
      GlyphCations->SetSourceConnection(Atom->GetOutputPort());
      GlyphCations->SetColorModeToColorByVector();
      GlyphCations->SetScaleModeToScaleByScalar();

    // clip the cations atoms to the unit-cell and polish the cuts with surfaces
    vtkClipClosedSurface *ClipSurfaceCationAtoms=vtkClipClosedSurface::New();
       ClipSurfaceCationAtoms->SetClippingPlanes(PlanesCollectionCation);
       ClipSurfaceCationAtoms->SetInputConnection(GlyphCations->GetOutputPort());
       ClipSurfaceCationAtoms->PassPointDataOn();
       ClipSurfaceCationAtoms->TriangulationErrorDisplayOn();

    // handle the highlighted atoms
    if(CutOutShape==ON)
    {
      // the highlighted atoms are the ones that are inside the chosen clipping region
      vtkClipPolyData *CutOutHighLightedShapeCationAtoms=vtkClipPolyData::New();
        CutOutHighLightedShapeCationAtoms->SetInputConnection(ClipSurfaceCationAtoms->GetOutputPort());
        switch(CutOutType)
        {
          case RECTANGULAR_CYLINDER:
            CutOutHighLightedShapeCationAtoms->SetClipFunction(RectangularCylinder);
            break;
          case CYLINDER:
            CutOutHighLightedShapeCationAtoms->SetClipFunction(Cylinder);
            break;
          case SPHERE:
            CutOutHighLightedShapeCationAtoms->SetClipFunction(Sphere);
            break;
        }
        CutOutHighLightedShapeCationAtoms->InsideOutOn();

      vtkPolyDataMapper *MapperHighLightedCationAtoms=vtkPolyDataMapper::New();
        MapperHighLightedCationAtoms->SetInputConnection(CutOutHighLightedShapeCationAtoms->GetOutputPort());
        MapperHighLightedCationAtoms->SetLookupTable(LookupTable);

      vtkActor *ActorHighLightedCationAtoms=vtkActor::New();
        ActorHighLightedCationAtoms->SetMapper(MapperHighLightedCationAtoms);
        ActorHighLightedCationAtoms->SetPosition(0,0,0);
        ActorHighLightedCationAtoms->GetProperty()->SetOpacity(HighLightedCationOpacity);
      ren1->AddActor(ActorHighLightedCationAtoms);
    }

    // handle the lowlighted, regular, atoms
    // the lowlighted atoms are the ones that are outside the chosen clipping region, or all atoms when there is no highlighted volume
    vtkClipPolyData *CutOutLowLightedShapeCationAtoms=vtkClipPolyData::New();
      CutOutLowLightedShapeCationAtoms->SetInputConnection(ClipSurfaceCationAtoms->GetOutputPort());
      switch(CutOutType)
      {
        case RECTANGULAR_CYLINDER:
          CutOutLowLightedShapeCationAtoms->SetClipFunction(RectangularCylinder);
          break;
        case CYLINDER:
          CutOutLowLightedShapeCationAtoms->SetClipFunction(Cylinder);
          break;
        case SPHERE:
          CutOutLowLightedShapeCationAtoms->SetClipFunction(Sphere);
          break;
      }
      CutOutLowLightedShapeCationAtoms->InsideOutOff();

    vtkPolyDataMapper *MapperLowLightedCationAtoms=vtkPolyDataMapper::New();
      if(CutOutShape==ON)
        MapperLowLightedCationAtoms->SetInputConnection(CutOutLowLightedShapeCationAtoms->GetOutputPort());
      else
        MapperLowLightedCationAtoms->SetInputConnection(ClipSurfaceCationAtoms->GetOutputPort());
      MapperLowLightedCationAtoms->SetLookupTable(LookupTable);

    vtkActor *ActorLowLightedCationAtoms=vtkActor::New();
      ActorLowLightedCationAtoms->SetMapper(MapperLowLightedCationAtoms);
      ActorLowLightedCationAtoms->SetPosition(0,0,0);
      ActorLowLightedCationAtoms->GetProperty()->SetOpacity(CationOpacity);
    ren1->AddActor(ActorLowLightedCationAtoms);

    // render bonds of the cations
    // ======================================================================================================================

    // render the bonds as white tubes
    vtkTubeFilter *BondsTubeCations=vtkTubeFilter::New();
      BondsTubeCations->SetInputConnection(AtomReaderCations->GetOutputPort());
      BondsTubeCations->SetRadius(1.0*ScaleFactor);
      BondsTubeCations->CappingOn();
      BondsTubeCations->SetNumberOfSides(Resolution);

    // clip the cations bonds to the unit-cell and polish the cuts with surfaces
    vtkClipClosedSurface *ClipSurfaceCations=vtkClipClosedSurface::New();
       ClipSurfaceCations->SetClippingPlanes(PlanesCollectionCation);
       ClipSurfaceCations->SetInputConnection(BondsTubeCations->GetOutputPort());
       ClipSurfaceCations->TriangulationErrorDisplayOn();

    // handle the highlighted bonds
    if(CutOutShape==ON)
    {
        // the highlighted bonds are the ones that are inside the chosen clipping region
        vtkClipPolyData *CutOutHighLightedShapeCationBonds=vtkClipPolyData::New();
          CutOutHighLightedShapeCationBonds->SetInputConnection(ClipSurfaceCations->GetOutputPort());
          switch(CutOutType)
          {
            case RECTANGULAR_CYLINDER:
              CutOutHighLightedShapeCationBonds->SetClipFunction(RectangularCylinder);
              break;
            case CYLINDER:
              CutOutHighLightedShapeCationBonds->SetClipFunction(Cylinder);
              break;
            case SPHERE:
              CutOutHighLightedShapeCationBonds->SetClipFunction(Sphere);
              break;
          }
          CutOutHighLightedShapeCationBonds->InsideOutOn();

        vtkPolyDataMapper *BondsMapperCations=vtkPolyDataMapper::New();
          BondsMapperCations->SetInputConnection(CutOutHighLightedShapeCationBonds->GetOutputPort());
          BondsMapperCations->ScalarVisibilityOff();
        vtkActor *BondsActorCations=vtkActor::New();
          BondsActorCations->SetMapper(BondsMapperCations);
          BondsActorCations->GetProperty()->SetColor(0.9,0.9,0.9);
          BondsActorCations->GetProperty()->SetOpacity(HighLightedCationOpacity);
          BondsActorCations->SetPosition(0,0,0);
        ren1->AddActor(BondsActorCations);
    }

    // handle the lowlighted, regular, bonds
    // the lowlighted bonds are the ones that are outside the chosen clipping region, or all bonds when there is no highlighted volume
    vtkClipPolyData *CutOutLowLightedShapeCationBonds=vtkClipPolyData::New();
      CutOutLowLightedShapeCationBonds->SetInputConnection(ClipSurfaceCations->GetOutputPort());
      switch(CutOutType)
      {
        case RECTANGULAR_CYLINDER:
          CutOutLowLightedShapeCationBonds->SetClipFunction(RectangularCylinder);
          break;
        case CYLINDER:
          CutOutLowLightedShapeCationBonds->SetClipFunction(Cylinder);
          break;
        case SPHERE:
          CutOutLowLightedShapeCationBonds->SetClipFunction(Sphere);
          break;
      }
      CutOutLowLightedShapeCationBonds->InsideOutOff();

    vtkPolyDataMapper *BondsMapperCations=vtkPolyDataMapper::New();
      if(CutOutShape==ON)
        BondsMapperCations->SetInputConnection(CutOutLowLightedShapeCationBonds->GetOutputPort());
      else 
        BondsMapperCations->SetInputConnection(ClipSurfaceCations->GetOutputPort());
      BondsMapperCations->ScalarVisibilityOff();
    vtkActor *BondsActorCations=vtkActor::New();
      BondsActorCations->SetMapper(BondsMapperCations);
      BondsActorCations->GetProperty()->SetColor(0.9,0.9,0.9);
      BondsActorCations->GetProperty()->SetOpacity(CationOpacity);
      BondsActorCations->SetPosition(0,0,0);
    ren1->AddActor(BondsActorCations);
  }
  else
   printf("Skipping rendering of the cations (file '%s' does not exist)\n",CationAtomsFilename);


  if(FileExists(FrameworkSurfaceFilename)&&(FrameworkSurface==ON))
  {
    // Read the data from a vtk file
    // this is the 3d histogram data
    vtkStructuredPointsReader *FrameworkSurfaceReader = vtkStructuredPointsReader::New();
      FrameworkSurfaceReader->SetFileName(FrameworkSurfaceFilename);
      FrameworkSurfaceReader->Update();

    switch(FrameworkRenderingMethod)
    {
      case ISOSURFACE:
        FrameworkSurfaceIsocontour=vtkMarchingCubes::New();
          FrameworkSurfaceIsocontour->SetInputConnection(FrameworkSurfaceReader->GetOutputPort());
          FrameworkSurfaceIsocontour->SetValue(0,50000);
          FrameworkSurfaceIsocontour->ComputeGradientsOn();
          FrameworkSurfaceIsocontour->ComputeScalarsOn();
        FrameworkSurfaceNormals=vtkPolyDataNormals::New();
           FrameworkSurfaceNormals->SetInputConnection(FrameworkSurfaceIsocontour->GetOutputPort());

        CutFullCell=vtkClipDataSet::New();
          CutFullCell->SetInputConnection(FrameworkSurfaceNormals->GetOutputPort());
          CutFullCell->SetClipFunction(BoxClipPlanes);
          CutFullCell->InsideOutOn();

        if(CutOutShape==ON)
        {
          vtkClipDataSet *CutOutHighLightedShape=vtkClipDataSet::New();
            CutOutHighLightedShape->SetInputConnection(CutFullCell->GetOutputPort());
            switch(CutOutType)
            {
              case RECTANGULAR_CYLINDER:
                CutOutHighLightedShape->SetClipFunction(RectangularCylinderFramework);
                break;
              case CYLINDER:
                CutOutHighLightedShape->SetClipFunction(CylinderFramework);
                break;
              case SPHERE:
                CutOutHighLightedShape->SetClipFunction(SphereFramework);
                break;
            }
            CutOutHighLightedShape->InsideOutOn();

          vtkDataSetMapper *FrameworkSurfaceHighLightMapper=vtkDataSetMapper::New();
            if(CutOutShape==ON)
              FrameworkSurfaceHighLightMapper->SetInputConnection(CutOutHighLightedShape->GetOutputPort());
            else
              FrameworkSurfaceHighLightMapper->SetInputConnection(FrameworkSurfaceNormals->GetOutputPort());
            FrameworkSurfaceHighLightMapper->ScalarVisibilityOff();
            FrameworkSurfaceHighLightMapper->SetScalarRange(0,60000);
            FrameworkSurfaceHighLightMapper->ImmediateModeRenderingOn();

          IsoContourArray=vtkActor::New();
          IsoContourArray->SetMapper(FrameworkSurfaceHighLightMapper);
          IsoContourArray->GetProperty()->SetInterpolationToPhong();
          switch(HighLightedIsoSurfaceMaterial)
          {
            case GLASS:
              IsoContourArray->GetProperty()->SetOpacity(HighLightedIsoSurfaceOpacity);
              IsoContourArray->GetProperty()->SetSpecular(0.65);
              IsoContourArray->GetProperty()->SetDiffuse(0.5);
              IsoContourArray->GetProperty()->SetAmbient(0.0);
              IsoContourArray->GetProperty()->SetColor(HighLightedIsoSurfaceColor.x,HighLightedIsoSurfaceColor.y,HighLightedIsoSurfaceColor.z);
              break;
            case METALLIC_PASTEL:
              IsoContourArray->GetProperty()->SetOpacity(HighLightedIsoSurfaceOpacity);
              IsoContourArray->GetProperty()->SetSpecular(0.55);
              IsoContourArray->GetProperty()->SetDiffuse(0.26);
              IsoContourArray->GetProperty()->SetAmbient(0.0);
              IsoContourArray->GetProperty()->SetColor(HighLightedIsoSurfaceColor.x,HighLightedIsoSurfaceColor.y,HighLightedIsoSurfaceColor.z);
              break;
            case TRANSPARENT:
              IsoContourArray->GetProperty()->SetOpacity(HighLightedIsoSurfaceOpacity);
              IsoContourArray->GetProperty()->SetSpecular(0.5);
              IsoContourArray->GetProperty()->SetDiffuse(0.65);
              IsoContourArray->GetProperty()->SetAmbient(0.0);
              IsoContourArray->GetProperty()->SetColor(HighLightedIsoSurfaceColor.x,HighLightedIsoSurfaceColor.y,HighLightedIsoSurfaceColor.z);
              break;
            case BRUSHED_METAL:
              IsoContourArray->GetProperty()->SetOpacity(HighLightedIsoSurfaceOpacity);
              IsoContourArray->GetProperty()->SetSpecular(0.34);
              IsoContourArray->GetProperty()->SetDiffuse(0.39);
              IsoContourArray->GetProperty()->SetAmbient(0.08);
              IsoContourArray->GetProperty()->SetColor(HighLightedIsoSurfaceColor.x,HighLightedIsoSurfaceColor.y,HighLightedIsoSurfaceColor.z);
              break;
            case RASPA:
            default:
              IsoContourArray->GetProperty()->SetOpacity(HighLightedIsoSurfaceOpacity);
              IsoContourArray->GetProperty()->SetSpecular(0.1);
              IsoContourArray->GetProperty()->SetDiffuse(0.5);
              IsoContourArray->GetProperty()->SetAmbient(0.0);
              IsoContourArray->GetProperty()->SetColor(HighLightedIsoSurfaceColor.x,HighLightedIsoSurfaceColor.y,HighLightedIsoSurfaceColor.z);
              break;
          }

          IsoContourArray->SetPosition(shift.x*SIZE_X/max,shift.y*SIZE_X/max,shift.z*SIZE_X/max);
          ren1->AddActor(IsoContourArray);
        }

        CutOutLowLightedShape=vtkClipDataSet::New();
          CutOutLowLightedShape->SetInputConnection(CutFullCell->GetOutputPort());
          switch(CutOutType)
          {
            case RECTANGULAR_CYLINDER:
              CutOutLowLightedShape->SetClipFunction(RectangularCylinderFramework);
              break;
            case CYLINDER:
              CutOutLowLightedShape->SetClipFunction(CylinderFramework);
              break;
            case SPHERE:
              CutOutLowLightedShape->SetClipFunction(SphereFramework);
              break;
          }
          CutOutLowLightedShape->InsideOutOff();

          FrameworkSurfaceLowLightMapper=vtkDataSetMapper::New();
          if(CutOutShape==ON)
            FrameworkSurfaceLowLightMapper->SetInputConnection(CutOutLowLightedShape->GetOutputPort());
          else
            FrameworkSurfaceLowLightMapper->SetInputConnection(CutFullCell->GetOutputPort());
          FrameworkSurfaceLowLightMapper->ScalarVisibilityOff();
          FrameworkSurfaceLowLightMapper->SetScalarRange(0,60000);
          FrameworkSurfaceLowLightMapper->ImmediateModeRenderingOn();

        IsoContourArray=vtkActor::New();
        IsoContourArray->SetMapper(FrameworkSurfaceLowLightMapper);
        IsoContourArray->GetProperty()->SetInterpolationToPhong();
        switch(IsoSurfaceMaterial)
        {
          case GLASS:
            IsoContourArray->GetProperty()->SetOpacity(IsoSurfaceOpacity);
            IsoContourArray->GetProperty()->SetSpecular(0.65);
            IsoContourArray->GetProperty()->SetDiffuse(0.5);
            IsoContourArray->GetProperty()->SetAmbient(0.0);
            break;
          case METALLIC_PASTEL:
            IsoContourArray->GetProperty()->SetOpacity(IsoSurfaceOpacity);
            IsoContourArray->GetProperty()->SetSpecular(0.55);
            IsoContourArray->GetProperty()->SetDiffuse(0.26);
            IsoContourArray->GetProperty()->SetAmbient(0.0);
            break;
          case TRANSPARENT:
            IsoContourArray->GetProperty()->SetOpacity(IsoSurfaceOpacity);
            IsoContourArray->GetProperty()->SetSpecular(0.5);
            IsoContourArray->GetProperty()->SetDiffuse(0.65);
            IsoContourArray->GetProperty()->SetAmbient(0.0);
            break;
          case BRUSHED_METAL:
            IsoContourArray->GetProperty()->SetOpacity(IsoSurfaceOpacity);
            IsoContourArray->GetProperty()->SetSpecular(0.34);
            IsoContourArray->GetProperty()->SetDiffuse(0.39);
            IsoContourArray->GetProperty()->SetAmbient(0.08);
            break;
          case RASPA:
          default:
            IsoContourArray->GetProperty()->SetOpacity(IsoSurfaceOpacity);
            IsoContourArray->GetProperty()->SetSpecular(0.1);
            IsoContourArray->GetProperty()->SetDiffuse(0.5);
            IsoContourArray->GetProperty()->SetAmbient(0.0);
            break;
        }

        IsoContourArray->SetPosition(shift.x*SIZE_X/max,shift.y*SIZE_X/max,shift.z*SIZE_X/max);
        ren1->AddActor(IsoContourArray);
        break;
      case VOLUME_RENDERING:
        // Create a property for the volume and set the transfer functions.
        // Turn shading on and use trilinear interpolation
        vtkVolumeProperty *FrameworkSurfaceProperty = vtkVolumeProperty::New();
          FrameworkSurfaceProperty->SetColor(cTFunN);
          FrameworkSurfaceProperty->SetScalarOpacity(oTFunSurface);
          FrameworkSurfaceProperty->SetInterpolationTypeToLinear();
          FrameworkSurfaceProperty->ShadeOn();
          FrameworkSurfaceProperty->SetSpecular(0.9);

        // Create a ray function - this is a compositing ray function
        vtkVolumeRayCastCompositeFunction *FrameworkSurfacecompositeFunction=vtkVolumeRayCastCompositeFunction::New();
          FrameworkSurfacecompositeFunction->SetCompositeMethodToInterpolateFirst();

        // Create the volume mapper and set the ray function and scalar input
        vtkVolumeRayCastMapper *FrameworkSurfaceMapper = vtkVolumeRayCastMapper::New();
          FrameworkSurfaceMapper->SetInputConnection(FrameworkSurfaceReader->GetOutputPort());
          FrameworkSurfaceMapper->SetVolumeRayCastFunction(FrameworkSurfacecompositeFunction);
          FrameworkSurfaceMapper->SetClippingPlanes(BoxClipPlaneCollection);
          FrameworkSurfaceMapper->SetSampleDistance(SampleDistance);
          FrameworkSurfaceMapper->SetImageSampleDistance(ImageSampleDistance);

        VolumeArray=vtkVolume::New();
          VolumeArray->SetMapper(FrameworkSurfaceMapper);
          VolumeArray->SetProperty(FrameworkSurfaceProperty);
          VolumeArray->SetPosition(shift.x*SIZE_X/max,shift.y*SIZE_X/max,shift.z*SIZE_X/max);
        ren1->AddVolume(VolumeArray);
      break;
    }
  }
  else
   printf("Skipping volume-rendering of the framework surface (file '%s' does not exist)\n",FrameworkSurfaceFilename);

  if(FileExists(DensityFilename)&&(Density==ON))
  {
    // Read the data from a vtk file
    // this is the 3d histogram data
    vtkStructuredPointsReader *DensityReader = vtkStructuredPointsReader::New();
      DensityReader->SetFileName(DensityFilename);
      DensityReader->Update();

    // Create a property for the volume and set the transfer functions.
    // Turn shading on and use trilinear interpolation
    vtkVolumeProperty *DensityProperty = vtkVolumeProperty::New();
      DensityProperty->SetColor(cTFunN2);
      DensityProperty->SetScalarOpacity(oTFunDensity2);
      DensityProperty->SetInterpolationTypeToLinear();
      DensityProperty->ShadeOn();

    // Create a ray function - this is a compositing ray function
    vtkVolumeRayCastCompositeFunction *DensitycompositeFunction=vtkVolumeRayCastCompositeFunction::New();
       DensitycompositeFunction->SetCompositeMethodToInterpolateFirst();

    // Create the volume mapper and set the ray function and scalar input
    // Note: density does not need tobe clipped (zero values will be removed anyway)
    vtkVolumeRayCastMapper *DensityMapper = vtkVolumeRayCastMapper::New();
      DensityMapper->SetInputConnection(DensityReader->GetOutputPort());
      DensityMapper->SetVolumeRayCastFunction(DensitycompositeFunction);

      DensityArray=vtkVolume::New();
      DensityArray->SetMapper(DensityMapper);
      DensityArray->SetProperty(DensityProperty);
      DensityArray->SetPosition(shift.x*SIZE_X/max,shift.y*SIZE_Y/max,shift.z*SIZE_Z/max);
      ren1->AddVolume(DensityArray);
  }
  else
   printf("Skipping volume-rendering of the density (file '%s' does not exist or 'Density=OFF')\n",DensityFilename);

  // handle credits
  if(Credits==ON)
  {
    vtkTextProperty *multiLineTextProp=vtkTextProperty::New();
      multiLineTextProp->SetFontSize(12);
      multiLineTextProp->SetFontFamilyToArial();
      multiLineTextProp->BoldOn();
      multiLineTextProp->SetLineSpacing(1.0);
      multiLineTextProp->SetVerticalJustificationToTop();
      multiLineTextProp->SetColor(0.8,0.8,0.8);
      multiLineTextProp->SetJustificationToLeft();
  
    vtkTextMapper *textMapperL=vtkTextMapper::New();
      textMapperL->SetInput(credits);
      textMapperL->SetTextProperty(multiLineTextProp);
    vtkActor2D *textActorL=vtkActor2D::New();
      textActorL->SetMapper(textMapperL);
      textActorL->GetPositionCoordinate()->SetCoordinateSystemToNormalizedDisplay();
      textActorL->SetPosition(0.005,0.06);
    ren1->AddActor(textActorL);
  }

  // handle title
  if(Title==ON)
  {
    vtkTextProperty *multiLineTextProp2=vtkTextProperty::New();
      multiLineTextProp2->SetFontSize(28);
      multiLineTextProp2->SetFontFamilyToArial();
      multiLineTextProp2->BoldOn();
      multiLineTextProp2->SetLineSpacing(1.0);
      multiLineTextProp2->SetVerticalJustificationToTop();
      multiLineTextProp2->SetColor(0.2,0.2,0.2);
      multiLineTextProp2->SetJustificationToCentered();
      multiLineTextProp2->SetVerticalJustificationToTop();

    vtkTextMapper *textMapperL2=vtkTextMapper::New();
      textMapperL2->SetInput(title);
      textMapperL2->SetTextProperty(multiLineTextProp2);
    vtkActor2D *textActorL2=vtkActor2D::New();
      textActorL2->SetMapper(textMapperL2);
      textActorL2->GetPositionCoordinate()->SetCoordinateSystemToNormalizedDisplay();
      textActorL2->SetPosition(0.5,0.975);
    ren1->AddActor(textActorL2);
  }

  // set viewing angle and reset view-up
  ren1->GetActiveCamera()->Azimuth(Azimuth);
  ren1->GetActiveCamera()->OrthogonalizeViewUp();
  ren1->GetActiveCamera()->Elevation(Elevation);
  ren1->GetActiveCamera()->OrthogonalizeViewUp();
  ren1->GetActiveCamera()->Roll(Roll);
  ren1->GetActiveCamera()->OrthogonalizeViewUp();

  // set additional viewing angles
  ren1->GetActiveCamera()->Azimuth(AdditionalAzimuth);
  ren1->GetActiveCamera()->Elevation(AdditionalElevation);
  ren1->GetActiveCamera()->Roll(AdditionalRoll);

  // set parallel projection and clipping-range of the camera
  ren1->GetActiveCamera()->ParallelProjectionOn();
  ren1->GetActiveCamera()->SetClippingRange(0.001,10000);
  ren1->ResetCamera();
  ren1->GetActiveCamera()->Zoom(ZoomFactor);

  vtkLightKit *light = vtkLightKit::New();
  light->MaintainLuminanceOn();
  light->AddLightsToRenderer(ren1);

  // select a white background color
  ren1->SetBackground(1.0,1.0,1.0);
  ren1->BackingStoreOn();

  // set the size of the picture in pixels
  renWin->SetSize(ImageSizeX,ImageSizeY);

  // set this to 8 for final picture
  renWin->SetAAFrames(AA);
  renWin->Render();

  vtkRenderLargeImage *renderLarge =vtkRenderLargeImage::New();
    renderLarge->SetInput(ren1);
    renderLarge->SetMagnification(magnification);

  vtkPNGWriter *writer_png =vtkPNGWriter::New();
    sprintf(buffer,"%s.png",OutputFileName);
    writer_png->SetFileName(buffer);
    writer_png->SetInputConnection(renderLarge->GetOutputPort());
    writer_png->Write();

  vtkJPEGWriter *writer =vtkJPEGWriter::New();
    sprintf(buffer,"%s.jpg",OutputFileName);
    writer->SetFileName(buffer);
    writer->SetInputConnection(renderLarge->GetOutputPort());
    writer->Write();

  if(Movie==ON)
  {

    for(int i=0;i<=360;i+=1)
    {
      printf("Angle: %d\n",i);
      if(i<10)
        sprintf(buffer,"00%d.jpg",i);
      else if(i<100)
        sprintf(buffer,"0%d.jpg",i);
      else sprintf(buffer,"%d.jpg",i);
      writer->SetFileName(buffer);
      writer->Write();
      ren1->GetActiveCamera()->Azimuth(1.0);
      renWin->Render();
      ren1->Modified();
    }
  }

  // Interact with the data at 3 frames per second
  iren->Initialize();
  iren->Start();

  // Clean up
  ren1->Delete();
  renWin->Delete();
  iren->Delete();
}

// Reads at most 'length-1' characters from a file into a buffer
// The character sequence is zero-terminated
// more characters then the line are discarded
char *ReadLine(char *buffer, size_t length, FILE *file)
{
  char *p;
  size_t last;

  if((p=fgets(buffer,length,file)))
  {
    last=strlen(buffer)-1;

    if(buffer[last]=='\n')
      buffer[last]='\0'; // discard the trailing newline
    else
    {
      fscanf(file,"%*[^\n]");
      (void) fgetc(file); // discard the newline
    }
  }
  return p;
}

void TrimStringInPlace(char *s)
{
  unsigned int i=0,j;

  // Trim spaces and tabs from beginning:
  while(isspace(s[i]))
    i++;
  if(i>0)
  {
    for(j=0;j<strlen(s);j++)
      s[j]=s[j+i];
    s[j]='\0';
  }

  // Trim spaces and tabs from end:
  i=strlen(s)-1;
  while(isspace(s[i]))
    i--;
  if(i<(strlen(s)-1))
    s[i+1]='\0';
}


int RemoveQuotesAroundString(char *string)
{
  char *first,*last;

  if((first=strchr(string,'\'')))
  {
    memmove(string,first+1,strlen(string));
    if((last=strrchr(string,'\''))) *last='\0';
    return ON;
  }

  if((first=strchr(string,'\"')))
  {
    memmove(string,first+1,strlen(string));
    if((last=strrchr(string,'\"'))) *last='\0';
    return ON;
  }
  return OFF;
}

int ReadInputFile(void)
{
  FILE *FilePtr;
  char line[1024];
  char keyword[256],arguments[256],firstargument[256];

  if((FilePtr=fopen("./input","r")))
  {
    printf("reading input\n");
    while(ReadLine(line,16384,FilePtr))
    {
      // extract first word
      strcpy(keyword,"keyword");
      sscanf(line,"%s%[^\n]",keyword,arguments);
      sscanf(arguments,"%s",firstargument);

      if((strcasecmp("BoxLengths",keyword)==0)||(strcasecmp("CellLengths",keyword)==0))
        sscanf(arguments,"%lf %lf %lf",&A,&B,&C);
      if((strcasecmp("BoxAngles",keyword)==0)||(strcasecmp("CellAngles",keyword)==0))
        sscanf(arguments,"%lf %lf %lf",&AlphaAngle,&BetaAngle,&GammaAngle);
      if(strcasecmp("ZoomFactor",keyword)==0) sscanf(arguments,"%lf",&ZoomFactor);
      if(strcasecmp("ScaleFactor",keyword)==0) sscanf(arguments,"%lf",&ScaleFactor);
      if(strcasecmp("VisibilityFractionSurfaceMax",keyword)==0) sscanf(arguments,"%lf %lf %lf",&VisibilityFractionSurfaceMax.x,
         &VisibilityFractionSurfaceMax.y,&VisibilityFractionSurfaceMax.z);
      if(strcasecmp("VisibilityFractionAdsorbateMax",keyword)==0) sscanf(arguments,"%lf %lf %lf",&VisibilityFractionAdsorbateMax.x,
         &VisibilityFractionAdsorbateMax.y,&VisibilityFractionAdsorbateMax.z);
      if(strcasecmp("VisibilityFractionCationMax",keyword)==0) sscanf(arguments,"%lf %lf %lf",&VisibilityFractionCationMax.x,
         &VisibilityFractionCationMax.y,&VisibilityFractionCationMax.z);
      if(strcasecmp("Azimuth",keyword)==0) sscanf(arguments,"%lf",&Azimuth);
      if(strcasecmp("Elevation",keyword)==0) sscanf(arguments,"%lf",&Elevation);
      if(strcasecmp("Roll",keyword)==0) sscanf(arguments,"%lf",&Roll);
      if(strcasecmp("AdditionalAzimuth",keyword)==0) sscanf(arguments,"%lf",&AdditionalAzimuth);
      if(strcasecmp("AdditionalElevation",keyword)==0) sscanf(arguments,"%lf",&AdditionalElevation);
      if(strcasecmp("AdditionalRoll",keyword)==0) sscanf(arguments,"%lf",&AdditionalRoll);
      if(strcasecmp("ImageSize",keyword)==0) sscanf(arguments,"%d %d",&ImageSizeX,&ImageSizeY);
      if(strcasecmp("Magnification",keyword)==0) sscanf(arguments,"%lf",&magnification);
      if(strcasecmp("Resolution",keyword)==0) sscanf(arguments,"%d",&Resolution);
      if(strcasecmp("AA",keyword)==0) sscanf(arguments,"%d",&AA);
      if(strcasecmp("SampleDistance",keyword)==0) sscanf(arguments,"%lf",&SampleDistance);
      if(strcasecmp("ImageSampleDistance",keyword)==0) sscanf(arguments,"%lf",&ImageSampleDistance);
      if(strcasecmp("OutputFileName",keyword)==0) 
      {
        TrimStringInPlace(arguments);
        if(RemoveQuotesAroundString(arguments)==ON)
          strcpy(OutputFileName,arguments);
        else
        {
          TrimStringInPlace(firstargument);
          strcpy(OutputFileName,firstargument);
        }
      }
      if(strcasecmp("TitleString",keyword)==0) 
      {
        TrimStringInPlace(arguments);
        RemoveQuotesAroundString(arguments);
        strcpy(title,arguments);
      }
      if(strcasecmp("Title",keyword)==0) 
      {
        if((strcasecmp("yes",firstargument)==0)||(strcasecmp("on",firstargument)==0)) Title=ON;
        if((strcasecmp("no",firstargument)==0)||(strcasecmp("off",firstargument)==0)) Title=OFF;
      }

      if(strcasecmp("Axes",keyword)==0) 
      {
        if((strcasecmp("yes",firstargument)==0)||(strcasecmp("on",firstargument)==0)) Axes=ON;
        if((strcasecmp("no",firstargument)==0)||(strcasecmp("off",firstargument)==0)) Axes=OFF;
      }
      if(strcasecmp("AxesFontFactor",keyword)==0) sscanf(arguments,"%lf",&AxesFontFactor);
      if(strcasecmp("AxesLabelFactor",keyword)==0) sscanf(arguments,"%lf",&AxesLabelFactor);
      if(strcasecmp("AxesColor",keyword)==0) sscanf(arguments,"%lf %lf %lf",&AxesColor.x,&AxesColor.y,&AxesColor.z);
      if(strcasecmp("AxesLabelColor",keyword)==0) sscanf(arguments,"%lf %lf %lf",&AxesLabelColor.x,&AxesLabelColor.y,&AxesLabelColor.z);
      if(strcasecmp("NumberOfAxesLabels",keyword)==0) sscanf(arguments,"%d",&NumberOfAxesLabels);
      if(strcasecmp("SetXAxisVisibility",keyword)==0) 
      {
        if((strcasecmp("yes",firstargument)==0)||(strcasecmp("on",firstargument)==0)) SetXAxisVisibility=ON;
        if((strcasecmp("no",firstargument)==0)||(strcasecmp("off",firstargument)==0)) SetXAxisVisibility=OFF;
      }
      if(strcasecmp("SetYAxisVisibility",keyword)==0) 
      {
        if((strcasecmp("yes",firstargument)==0)||(strcasecmp("on",firstargument)==0)) SetYAxisVisibility=ON;
        if((strcasecmp("no",firstargument)==0)||(strcasecmp("off",firstargument)==0)) SetYAxisVisibility=OFF;
      }
      if(strcasecmp("SetZAxisVisibility",keyword)==0) 
      {
        if((strcasecmp("yes",firstargument)==0)||(strcasecmp("on",firstargument)==0)) SetZAxisVisibility=ON;
        if((strcasecmp("no",firstargument)==0)||(strcasecmp("off",firstargument)==0)) SetZAxisVisibility=OFF;
      }

      if(strcasecmp("Frame",keyword)==0) 
      {
        if((strcasecmp("NO_FRAME",firstargument)==0)||(strcasecmp("NoFrame",firstargument)==0)) Frame=NO_FRAME;
        if((strcasecmp("FULL_CELL",firstargument)==0)||(strcasecmp("FullCell",firstargument)==0)) Frame=FULL_CELL;
        if((strcasecmp("UNIT_CELL",firstargument)==0)||(strcasecmp("UnitCell",firstargument)==0)) Frame=UNIT_CELL;
      }
      if((strcasecmp("Number_Of_Unit_Cells",keyword)==0)||(strcasecmp("NumberOfUnitCells",keyword)==0))
        sscanf(arguments,"%d %d %d",&NrDuplicatesX,&NrDuplicatesY,&NrDuplicatesZ);
      if((strcasecmp("Frame_Min",keyword)==0)||(strcasecmp("FrameMin",keyword)==0))
        sscanf(arguments,"%d %d %d",&FrameMin[0],&FrameMin[1],&FrameMin[2]);
      if((strcasecmp("Frame_Max",keyword)==0)||(strcasecmp("FrameMax",keyword)==0))
        sscanf(arguments,"%d %d %d",&FrameMax[0],&FrameMax[1],&FrameMax[2]);
      if(strcasecmp("FrameOpacity",keyword)==0) sscanf(arguments,"%lf",&FrameOpacity);
      if(strcasecmp("FrameThickness",keyword)==0) sscanf(arguments,"%lf",&FrameThickness);

      if(strcasecmp("Adsorbates",keyword)==0) 
      {
        if((strcasecmp("yes",firstargument)==0)||(strcasecmp("on",firstargument)==0)) Adsorbates=ON;
        if((strcasecmp("no",firstargument)==0)||(strcasecmp("off",firstargument)==0)) Adsorbates=OFF;
      }
      if(strcasecmp("AdsorbateAtomsFilename",keyword)==0)
      {
        TrimStringInPlace(arguments);
        if(RemoveQuotesAroundString(arguments)==ON)
          strcpy(AdsorbateAtomsFilename,arguments);
        else
        {
          TrimStringInPlace(firstargument);
          strcpy(AdsorbateAtomsFilename,firstargument);
        }
      }
      if(strcasecmp("AdsorbateOpacity",keyword)==0) sscanf(arguments,"%lf",&AdsorbateOpacity);
      if(strcasecmp("HighLightedAdsorbateOpacity",keyword)==0) sscanf(arguments,"%lf",&HighLightedAdsorbateOpacity);
      if(strcasecmp("ClipAdsorbates",keyword)==0)
      {
        if(strcasecmp("XYZ",firstargument)==0) ClipAdsorbates=XYZ;
        if(strcasecmp("XY",firstargument)==0) ClipAdsorbates=XY;
        if(strcasecmp("XZ",firstargument)==0) ClipAdsorbates=XZ;
        if(strcasecmp("YZ",firstargument)==0) ClipAdsorbates=YZ;
        if(strcasecmp("X",firstargument)==0) ClipAdsorbates=X;
        if(strcasecmp("Y",firstargument)==0) ClipAdsorbates=Y;
        if(strcasecmp("Z",firstargument)==0) ClipAdsorbates=Z;
        if((strcasecmp("None",firstargument)==0)||(strcasecmp("off",firstargument)==0)||(strcasecmp("no",firstargument)==0)) ClipAdsorbates=NONE;
      }

      if(strcasecmp("Cations",keyword)==0) 
      {
        if((strcasecmp("yes",firstargument)==0)||(strcasecmp("on",firstargument)==0)) Cations=ON;
        if((strcasecmp("no",firstargument)==0)||(strcasecmp("off",firstargument)==0)) Cations=OFF;
      }
      if(strcasecmp("CationAtomsFilename",keyword)==0)
      {
        TrimStringInPlace(arguments);
        if(RemoveQuotesAroundString(arguments)==ON)
          strcpy(CationAtomsFilename,arguments);
        else
        {
          TrimStringInPlace(firstargument);
          strcpy(CationAtomsFilename,firstargument);
        }
      }
      if(strcasecmp("CationOpacity",keyword)==0) sscanf(arguments,"%lf",&CationOpacity);
      if(strcasecmp("HighLightedCationOpacity",keyword)==0) sscanf(arguments,"%lf",&HighLightedCationOpacity);
      if(strcasecmp("ClipCations",keyword)==0)
      {
        if(strcasecmp("XYZ",firstargument)==0) ClipCations=XYZ;
        if(strcasecmp("XY",firstargument)==0) ClipCations=XY;
        if(strcasecmp("XZ",firstargument)==0) ClipCations=XZ;
        if(strcasecmp("YZ",firstargument)==0) ClipCations=YZ;
        if(strcasecmp("X",firstargument)==0) ClipCations=X;
        if(strcasecmp("Y",firstargument)==0) ClipCations=Y;
        if(strcasecmp("Z",firstargument)==0) ClipCations=Z;
        if((strcasecmp("None",firstargument)==0)||(strcasecmp("off",firstargument)==0)||(strcasecmp("no",firstargument)==0)) ClipCations=NONE;
      }

      if(strcasecmp("FrameworkSurface",keyword)==0) 
      {
        if((strcasecmp("yes",firstargument)==0)||(strcasecmp("on",firstargument)==0)) FrameworkSurface=ON;
        if((strcasecmp("no",firstargument)==0)||(strcasecmp("off",firstargument)==0)) FrameworkSurface=OFF;
      }
      if(strcasecmp("FrameworkSurfaceFilename",keyword)==0)
      {
        TrimStringInPlace(arguments);
        if(RemoveQuotesAroundString(arguments)==ON)
          strcpy(FrameworkSurfaceFilename,arguments);
        else
        {
          TrimStringInPlace(firstargument);
          strcpy(FrameworkSurfaceFilename,firstargument);
        }
      }
      if(strcasecmp("FrameworkRenderingMethod",keyword)==0)
      {
        if((strcasecmp("Iso_Surface",firstargument)==0)||(strcasecmp("IsoSurface",firstargument)==0)) FrameworkRenderingMethod=ISOSURFACE;
        if((strcasecmp("Volume_Rendering",firstargument)==0)||(strcasecmp("VolumeRendering",firstargument)==0)) FrameworkRenderingMethod=VOLUME_RENDERING;
      }
      if(strcasecmp("IsoSurfaceMaterial",keyword)==0)
      {
        if(strcasecmp("GLASS",firstargument)==0) IsoSurfaceMaterial=GLASS;
        if(strcasecmp("TRANSPARENT",firstargument)==0) IsoSurfaceMaterial=TRANSPARENT;
        if(strcasecmp("RASPA",firstargument)==0) IsoSurfaceMaterial=RASPA;
        if((strcasecmp("Brushed_Metal",firstargument)==0)||(strcasecmp("BrushedMetal",firstargument)==0)) IsoSurfaceMaterial=BRUSHED_METAL;
        if((strcasecmp("Metallic_Pastel",firstargument)==0)||(strcasecmp("MetallicPastel",firstargument)==0)) IsoSurfaceMaterial=METALLIC_PASTEL;
      }
      if(strcasecmp("HighLightedIsoSurfaceMaterial",keyword)==0)
      {
        if(strcasecmp("GLASS",firstargument)==0) HighLightedIsoSurfaceMaterial=GLASS;
        if(strcasecmp("TRANSPARENT",firstargument)==0) HighLightedIsoSurfaceMaterial=TRANSPARENT;
        if(strcasecmp("RASPA",firstargument)==0) HighLightedIsoSurfaceMaterial=RASPA;
        if((strcasecmp("Brushed_Metal",firstargument)==0)||(strcasecmp("BrushedMetal",firstargument)==0)) HighLightedIsoSurfaceMaterial=BRUSHED_METAL;
        if((strcasecmp("Metallic_Pastel",firstargument)==0)||(strcasecmp("MetallicPastel",firstargument)==0)) HighLightedIsoSurfaceMaterial=METALLIC_PASTEL;
      }
      if(strcasecmp("IsoSurfaceOpacity",keyword)==0) sscanf(arguments,"%lf",&IsoSurfaceOpacity);
      if(strcasecmp("HighLightedIsoSurfaceOpacity",keyword)==0) sscanf(arguments,"%lf",&HighLightedIsoSurfaceOpacity);
      if(strcasecmp("HighLightedIsoSurfaceColor",keyword)==0) sscanf(arguments,"%lf %lf %lf",&HighLightedIsoSurfaceColor.x,&HighLightedIsoSurfaceColor.y,&HighLightedIsoSurfaceColor.z);
      if(strcasecmp("FrameworkCutOutShape",keyword)==0) 
      {
        if((strcasecmp("yes",firstargument)==0)||(strcasecmp("on",firstargument)==0)) FrameworkCutOutShape=ON;
        if((strcasecmp("no",firstargument)==0)||(strcasecmp("off",firstargument)==0)) FrameworkCutOutShape=OFF;
      }

      if(strcasecmp("FrameworkAtomsFilename",keyword)==0)
      {
        TrimStringInPlace(arguments);
        if(RemoveQuotesAroundString(arguments)==ON)
          strcpy(FrameworkAtomsFilename,arguments);
        else
        {
          TrimStringInPlace(firstargument);
          strcpy(FrameworkAtomsFilename,firstargument);
        }
      }
      if(strcasecmp("FrameworkBondsFilename",keyword)==0)
      {
        TrimStringInPlace(arguments);
        if(RemoveQuotesAroundString(arguments)==ON)
          strcpy(FrameworkBondsFilename,arguments);
        else
        {
          TrimStringInPlace(firstargument);
          strcpy(FrameworkBondsFilename,firstargument);
        }
      }
      if(strcasecmp("FrameworkAtoms",keyword)==0)
      {
        if((strcasecmp("None",firstargument)==0)||(strcasecmp("off",firstargument)==0)||(strcasecmp("no",firstargument)==0)) FrameworkAtoms=OFF;
        if((strcasecmp("All",firstargument)==0)||(strcasecmp("on",firstargument)==0)||(strcasecmp("yes",firstargument)==0)) FrameworkAtoms=ON;
      }
      if(strcasecmp("FrameworkOpacity",keyword)==0) sscanf(arguments,"%lf",&FrameworkOpacity);
      if(strcasecmp("FrameworkBondThickness",keyword)==0) sscanf(arguments,"%lf",&FrameworkBondThickness);

      if(strcasecmp("Density",keyword)==0) 
      {
        if((strcasecmp("yes",firstargument)==0)||(strcasecmp("on",firstargument)==0)) Density=ON;
        if((strcasecmp("no",firstargument)==0)||(strcasecmp("off",firstargument)==0)) Density=OFF;
      }
      if(strcasecmp("DensityFilename",keyword)==0)
      {
        TrimStringInPlace(arguments);
        if(RemoveQuotesAroundString(arguments)==ON)
          strcpy(DensityFilename,arguments);
        else
        {
          TrimStringInPlace(firstargument);
          strcpy(DensityFilename,firstargument);
        }
      }


      if(strcasecmp("Movie",keyword)==0)
      {
        if((strcasecmp("None",firstargument)==0)||(strcasecmp("off",firstargument)==0)||(strcasecmp("no",firstargument)==0)) Movie=OFF;
        if((strcasecmp("All",firstargument)==0)||(strcasecmp("on",firstargument)==0)||(strcasecmp("yes",firstargument)==0)) Movie=ON;
      }
      if(strcasecmp("CutOutShape",keyword)==0)
      {
        if((strcasecmp("None",firstargument)==0)||(strcasecmp("off",firstargument)==0)||(strcasecmp("no",firstargument)==0)) CutOutShape=OFF;
        if((strcasecmp("All",firstargument)==0)||(strcasecmp("on",firstargument)==0)||(strcasecmp("yes",firstargument)==0)) CutOutShape=ON;
      }
      if(strcasecmp("CutOutDirection",keyword)==0)
      {
        if(strcasecmp("X",firstargument)==0) CutOutDirection=X;
        if(strcasecmp("Y",firstargument)==0) CutOutDirection=Y;
        if(strcasecmp("Z",firstargument)==0) CutOutDirection=Z;
      }
      if(strcasecmp("CutOutRotationAngle",keyword)==0) sscanf(arguments,"%lf",&CutOutRotationAngle);
      if(strcasecmp("CutOutRadius1",keyword)==0) sscanf(arguments,"%lf",&CutOutRadius1);
      if(strcasecmp("CutOutRadius2",keyword)==0) sscanf(arguments,"%lf",&CutOutRadius2);
      if(strcasecmp("CutOutFractionalShift",keyword)==0) sscanf(arguments,"%lf %lf %lf",&CutOutFractionalShift.x,&CutOutFractionalShift.y,&CutOutFractionalShift.z);
      if(strcasecmp("CutOutType",keyword)==0)
      {
        if(strcasecmp("Cylinder",firstargument)==0) CutOutType=CYLINDER;
        if((strcasecmp("Rectangular_Cylinder",firstargument)==0||(strcasecmp("RectangularCylinder",firstargument)==0))) CutOutType=RECTANGULAR_CYLINDER;
        if(strcasecmp("Sphere",firstargument)==0) CutOutType=SPHERE;
      }
    }
  }
  return 0;
}
