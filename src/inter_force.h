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

#ifndef INTER_FORCE_H
#define INTER_FORCE_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

void CalculateTotalInterVDWForce(void);
int CalculateTotalInterReplicaVDWForce(void);

int CalculateTotalInterChargeChargeCoulombForce(void);
int CalculateTotalInterReplicaChargeChargeCoulombForce(void);

int CalculateTotalInterChargeBondDipoleCoulombForce(void);
int CalculateTotalInterBondDipoleBondDipoleCoulombForce(void);

int CalculateTotalInterChargeChargeCoulombElectricFieldMC(int New,int excl_ads,int excl_cation);
int CalculateInterElectricFieldFromInducedDipoleMC(int New,int excl_ads,int excl_cation);

int CalculateInterElectricFieldFromInducedDipoles(void);

int CalculateInterChargeInducedDipoleForce(void);
int CalculateInterInducedDipoleInducedDipoleForce(void);

#endif

