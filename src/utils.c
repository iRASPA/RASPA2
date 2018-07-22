/*************************************************************************************************************
    RASPA: a molecular-dynamics, monte-carlo and optimization code for nanoporous materials
    Copyright (C) 2006-2017 David Dubbeldam, Sofia Calero, Thijs Vlugt, Donald E. Ellis, and Randall Q. Snurr.

    D.Dubbeldam@uva.nl            http://molsim.science.uva.nl/
    scaldia@upo.es                http://www.upo.es/raspa/
    t.j.h.vlugt@tudelft.nl        http://homepage.tudelft.nl/v9k6y
    don-ellis@northwestern.edu    http://dvworld.northwestern.edu/
    snurr@northwestern.edu        http://zeolites.cqe.northwestern.edu/

    This file 'utils.c' is part of RASPA-2.0

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
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <ctype.h>
#include <sys/time.h>
#include "molecule.h"
#include "utils.h"
#include "integrate.h"
#include "simulation.h"
#include "framework.h"
#include "framework_energy.h"
#include "framework_force.h"
#include "internal_force.h"
#include "recrossing.h"

#ifdef HAVE_FFTW3
#include <fftw3.h>
#endif

long Factorial(int n)
{
  int c;
  long result = 1;
 
  for (c = 1; c <= n; c++)
    result = result * c;
 
  return result;
}


void BubbleSort(int list[], int n)
{
  long c, d, t;

  for (c = 0 ; c < ( n - 1 ); c++)
  {
    for (d = 0 ; d < n - c - 1; d++)
    {
      if (list[d] > list[d+1])
      {
        /* Swapping */

        t         = list[d];
        list[d]   = list[d+1];
        list[d+1] = t;
      }
    }
  }
}



int isInArrayOfSize(int value,int n,int *array)
{
  int i;

  for(i=0;i<n;i++)
  {
    if(value==array[i]) return TRUE;
  }
  return FALSE;
}

double get_wall_time(void)
{
  struct timeval time;
  if (gettimeofday(&time,NULL)) return 0;
  return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

double get_cpu_time(void)
{
  return (double)clock() / CLOCKS_PER_SEC;
}

REAL Smoothing(REAL theta)
{
  REAL on,off;


  on=170.0*DEG2RAD;
  off=180.0*DEG2RAD;

  if(theta<on) return 1.0;
  return SQR(off-theta)*(off+2.0*theta-3.0*on)/CUBE(off-on);
}

REAL SmoothingDerivative(REAL theta)
{
  REAL on,off;

  on=170.0*DEG2RAD;
  off=180.0*DEG2RAD;

  if(theta<on) return 0.0;
  return 6.0*(off-theta)*(on-theta)/CUBE(off-on);
}

REAL Exp(REAL x)
{
  return (REAL)exp(x);
}

REAL MExp(REAL x)
{
  return (REAL)exp(-x);
}

REAL Exp2(REAL x)
{
  return (REAL)exp(2.0*x);
}

REAL MExp2(REAL x)
{
  return (REAL)exp(-2.0*x);
}

REAL Identity(REAL x)
{
  return (REAL)x;
}

REAL ExpRM(REAL x)
{
  return (REAL)(exp(2.0*x)*exp(-x));
}


// Rotate a point 'p' about the line through 'origin' parallel to 'dir' by an angle 'theta'
VECTOR RotateAboutArbitraryLine(VECTOR p, VECTOR origin, VECTOR dir,REAL theta)
{
  VECTOR vec;
  REAL r;

  // normalize the direction vector
  r=SQR(dir.x)+SQR(dir.y)+SQR(dir.z);
  dir.x/=r;
  dir.y/=r;
  dir.z/=r;

  vec.x=origin.x*(SQR(dir.y)+SQR(dir.z))+
        dir.x*(-origin.y*dir.y-origin.z*dir.z+dir.x*p.x+dir.y*p.y+dir.z*p.z)+
        ((p.x-origin.x)*(SQR(dir.y)+SQR(dir.z))+dir.x*(origin.y*dir.y+origin.z*dir.z-dir.y*p.y-dir.z*p.z))*cos(theta)+
        (origin.y*dir.z-origin.z*dir.y-dir.z*p.y+dir.y*p.z)*sin(theta);

  vec.y=origin.y*(SQR(dir.x)+SQR(dir.z))+
        dir.y*(-origin.x*dir.x-origin.z*dir.z+dir.x*p.x+dir.y*p.y+dir.z*p.z)+
        ((p.y-origin.y)*(SQR(dir.x)+SQR(dir.z))+dir.y*(origin.x*dir.x+origin.z*dir.z-dir.x*p.x-dir.z*p.z))*cos(theta)+
        (-origin.x*dir.z+origin.z*dir.x+dir.z*p.x-dir.x*p.z)*sin(theta);

  vec.z=origin.z*(SQR(dir.x)+SQR(dir.y))+
        dir.z*(-origin.x*dir.x-origin.y*dir.y+dir.x*p.x+dir.y*p.y+dir.z*p.z)+
        ((p.z-origin.z)*(SQR(dir.x)+SQR(dir.y))+dir.z*(origin.x*dir.x+origin.y*dir.y-dir.x*p.x-dir.y*p.y))*cos(theta)+
        (origin.x*dir.y-origin.y*dir.x-dir.y*p.x+dir.x*p.y)*sin(theta);

  return vec;
}

// Rotate a point 'p' about the the bondvector of two beads by an angle 'theta'
VECTOR RotateAboutBondVector(VECTOR p, VECTOR first_bead, VECTOR second_bead, REAL theta)
{
  VECTOR vec,dir;
  REAL r;

  dir.x=second_bead.x-first_bead.x;
  dir.y=second_bead.y-first_bead.y;
  dir.z=second_bead.z-first_bead.z;

  // normalize the direction vector
  r=SQR(dir.x)+SQR(dir.y)+SQR(dir.z);
  dir.x/=r;
  dir.y/=r;
  dir.z/=r;

  vec.x=first_bead.x*(SQR(dir.y)+SQR(dir.z))+
        dir.x*(-first_bead.y*dir.y-first_bead.z*dir.z+dir.x*p.x+dir.y*p.y+dir.z*p.z)+
        ((p.x-first_bead.x)*(SQR(dir.y)+SQR(dir.z))+dir.x*(first_bead.y*dir.y+first_bead.z*dir.z-dir.y*p.y-dir.z*p.z))*cos(theta)+
        (first_bead.y*dir.z-first_bead.z*dir.y-dir.z*p.y+dir.y*p.z)*sin(theta);

  vec.y=first_bead.y*(SQR(dir.x)+SQR(dir.z))+
        dir.y*(-first_bead.x*dir.x-first_bead.z*dir.z+dir.x*p.x+dir.y*p.y+dir.z*p.z)+
        ((p.y-first_bead.y)*(SQR(dir.x)+SQR(dir.z))+dir.y*(first_bead.x*dir.x+first_bead.z*dir.z-dir.x*p.x-dir.z*p.z))*cos(theta)+
        (-first_bead.x*dir.z+first_bead.z*dir.x+dir.z*p.x-dir.x*p.z)*sin(theta);

  vec.z=first_bead.z*(SQR(dir.x)+SQR(dir.y))+
        dir.z*(-first_bead.x*dir.x-first_bead.y*dir.y+dir.x*p.x+dir.y*p.y+dir.z*p.z)+
        ((p.z-first_bead.z)*(SQR(dir.x)+SQR(dir.y))+dir.z*(first_bead.x*dir.x+first_bead.y*dir.y-dir.x*p.x-dir.y*p.y))*cos(theta)+
        (first_bead.x*dir.y-first_bead.y*dir.x-dir.y*p.x+dir.x*p.y)*sin(theta);

  return vec;
}


void ConvertStringToLowercase(char *buffer)
{
  int i;

  for(i=0;i<strlen(buffer);i++)
    buffer[i]=(char)tolower(buffer[i]);
}

void ConvertStringToUppercase(char *buffer)
{
  int i;

  for(i=0;i<strlen(buffer);i++)
    buffer[i]=(char)toupper(buffer[i]);
}


#if defined (__LP64__) || defined (__64BIT__) || defined (_LP64) || (__WORDSIZE == 64)
  /*
     http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/VERSIONS/C-LANG/mt19937-64.c

     A C-program for MT19937-64 (2004/9/29 version).
     Coded by Takuji Nishimura and Makoto Matsumoto.

     This is a 64-bit version of Mersenne Twister pseudorandom number
     generator.

     Before using, initialize the state by using init_genrand64(seed)

     Copyright (C) 2004, Makoto Matsumoto and Takuji Nishimura,
     All rights reserved.

     Redistribution and use in source and binary forms, with or without
     modification, are permitted provided that the following conditions
     are met:

       1. Redistributions of source code must retain the above copyright
          notice, this list of conditions and the following disclaimer.

       2. Redistributions in binary form must reproduce the above copyright
          notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

       3. The names of its contributors may not be used to endorse or promote
          products derived from this software without specific prior written
          permission.

     THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
     "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
     LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
     A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
     CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
     EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
     PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
     PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
     LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
     NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
     SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

     References:
     T. Nishimura, ``Tables of 64-bit Mersenne Twisters''
       ACM Transactions on Modeling and
       Computer Simulation 10. (2000) 348--357.
     M. Matsumoto and T. Nishimura,
       ``Mersenne Twister: a 623-dimensionally equidistributed
         uniform pseudorandom number generator''
       ACM Transactions on Modeling and
       Computer Simulation 8. (Jan. 1998) 3--30.

     Any feedback is very welcome.
     http://www.math.hiroshima-u.ac.jp/~m-mat/MT/emt.html
     email: m-mat @ math.sci.hiroshima-u.ac.jp (remove spaces)
  */

  #define NN 312
  #define MM 156
  #define MATRIX_A 0xB5026F5AA96619E9ULL
  #define UM 0xFFFFFFFF80000000ULL /* Most significant 33 bits */
  #define LM 0x7FFFFFFFULL /* Least significant 31 bits */

  /* The array for the state vector */
  /* mti==NN+1 means mt[NN] is not initialized */
  static unsigned long long mt[NN];
  static unsigned long long mag01[2]={0ULL,MATRIX_A};
  static int mti=NN+1;
  static unsigned long long mt_bak[NN];
  static unsigned long long mag01_bak[2];
  static int mti_bak;

  /* initializes mt[NN] with a seed */
  void InitializeRandomNumberGenerator(unsigned long long seed)
  {
    mt[0]=seed;
    for(mti=1;mti<NN;mti++)
      mt[mti]=(6364136223846793005ULL*(mt[mti-1]^(mt[mti-1]>>62))+mti);
  }

  /* generates a random number on [0, 2^64-1]-interval */
  unsigned long long genrand64_int64(void)
  {
    int i;
    unsigned long long x;

    if(mti>=NN)
    { /* generate NN words at one time */

      /* if init_genrand64() has not been called, */
      /* a default initial seed is used     */
      if(mti==NN+1)
        InitializeRandomNumberGenerator(5489ULL);

      for(i=0;i<NN-MM;i++)
      {
        x=(mt[i]&UM)|(mt[i+1]&LM);
        mt[i]=mt[i+MM]^(x>>1)^mag01[(int)(x&1ULL)];
      }
      for(;i<NN-1;i++)
      {
        x=(mt[i]&UM)|(mt[i+1]&LM);
        mt[i]=mt[i+(MM-NN)]^(x>>1)^mag01[(int)(x&1ULL)];
      }
      x=(mt[NN-1]&UM)|(mt[0]&LM);
      mt[NN-1]=mt[MM-1]^(x>>1)^mag01[(int)(x&1ULL)];

      mti=0;
    }

    x=mt[mti++];

    x^=(x>>29)&0x5555555555555555ULL;
    x^=(x<<17)&0x71D67FFFEDA60000ULL;
    x^=(x<<37)&0xFFF7EEE000000000ULL;
    x^=(x>>43);

    return x;
  }

  /* generates a random number on [0, 2^63-1]-interval */
  long long RandomInteger(void)
  {
    return (long long)(genrand64_int64()>>1);
  }

  /* generates a random number on [0,1]-real-interval */
  double RandomNumber(void)
  {
    return (genrand64_int64()>>11)*(1.0/9007199254740991.0);
  }

  /* generates a random number on [0,1)-real-interval */
  double RandomNumber2(void)
  {
    return (genrand64_int64()>>11)*(1.0/9007199254740992.0);
  }

  /* generates a random number on (0,1)-real-interval */
  double RandomNumber3(void)
  {
    return ((genrand64_int64()>>12)+0.5)*(1.0/4503599627370496.0);
  }

  void StoreRandomNumberStatus(void)
  {
    int i;

    for(i=0;i<NN;i++)
      mt_bak[i]=mt[i];
    mti_bak=mti;
    mag01_bak[0]=mag01[0];
    mag01_bak[1]=mag01[1];
  }

  void RetrieveRandomNumberStatus(void)
  {
    int i;

    for(i=0;i<NN;i++)
      mt[i]=mt_bak[i];
    mti=mti_bak;
    mag01[0]=mag01_bak[0];
    mag01[1]=mag01_bak[1];
  }

  static int versionNumber=1;

  void WriteRestartUtils(FILE *FilePtr)
  {
    REAL Check;

    fwrite(&versionNumber,sizeof(int),1,FilePtr);

    fwrite(mt,NN,sizeof(unsigned long long),FilePtr);
    fwrite(&mti,1,sizeof(int),FilePtr);
    fwrite(mag01,1,sizeof(unsigned long long[2]),FilePtr);

    Check=123456789.0;
    fwrite(&Check,1,sizeof(REAL),FilePtr);
  }

  void ReadRestartUtils(FILE *FilePtr)
  {
    REAL Check;
    int readversionNumber=0;

    fread(&readversionNumber,sizeof(int),1,FilePtr);
    if(readversionNumber > versionNumber)
    {
      fprintf(stderr,"Upgrade to last version of RASPA to read binary restart-file");
      exit(-1);
    }

    fread(mt,NN,sizeof(unsigned long long),FilePtr);
    fread(&mti,1,sizeof(int),FilePtr);
    fread(mag01,1,sizeof(unsigned long long[2]),FilePtr);

    fread(&Check,1,sizeof(REAL),FilePtr);
    if(fabs(Check-123456789.0)>1e-10)
    {
      fprintf(stderr, "Error in binary restart-file (ReadRestartUtils)\n");
      ContinueAfterCrash=FALSE;
    }
  }

  #undef NN
  #undef MM
  #undef MATRIX_A
  #undef UM
  #undef LM

#else

  /*
     A C-program for MT19937, with initialization improved 2002/1/26.
     Coded by Takuji Nishimura and Makoto Matsumoto.

     Before using, initialize the state by using init_genrand(seed)

     Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
     All rights reserved.

     Redistribution and use in source and binary forms, with or without
     modification, are permitted provided that the following conditions
     are met:

       1. Redistributions of source code must retain the above copyright
          notice, this list of conditions and the following disclaimer.

       2. Redistributions in binary form must reproduce the above copyright
          notice, this list of conditions and the following disclaimer in the
          documentation and/or other materials provided with the distribution.

       3. The names of its contributors may not be used to endorse or promote
          products derived from this software without specific prior written
          permission.

     THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
     "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
     LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
     A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
     CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
     EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
     PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
     PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
     LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
     NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
     SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


     Any feedback is very welcome.
     http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
     email: m-mat @ math.sci.hiroshima-u.ac.jp (remove space)
  */

  /* Period parameters */
  #define N 624
  #define M 397
  #define MATRIX_A 0x9908b0dfUL   /* constant vector a */
  #define UPPER_MASK 0x80000000UL /* most significant w-r bits */
  #define LOWER_MASK 0x7fffffffUL /* least significant r bits */

  static unsigned long mt[N]; /* the array for the state vector  */
  static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */
  static unsigned long mag01[2]={0x0UL, MATRIX_A}; /* mag01[x] = x * MATRIX_A  for x=0,1 */
  static unsigned long mt_bak[N];
  static int mti_bak;
  static unsigned long mag01_bak[2];

  /* initializes mt[N] with a seed */
  void InitializeRandomNumberGenerator(unsigned long s)
  {
    mt[0]=s&0xffffffffUL;
    for(mti=1;mti<N;mti++)
    {
      mt[mti]=(1812433253UL*(mt[mti-1]^(mt[mti-1]>>30))+mti);
      /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
      /* In the previous versions, MSBs of the seed affect   */
      /* only MSBs of the array mt[].                        */
      /* 2002/01/09 modified by Makoto Matsumoto             */
      mt[mti]&=0xffffffffUL;
      /* for >32 bit machines */
    }
  }

  /* generates a random number on [0,0xffffffff]-interval */
  unsigned long genrand_int32(void)
  {
    unsigned long y;
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if(mti>=N)
    { /* generate N words at one time */
      int kk;

      if(mti==N+1)   /* if init_genrand() has not been called, */
        InitializeRandomNumberGenerator(5489UL); /* a default initial seed is used */

      for(kk=0;kk<N-M;kk++)
      {
        y=(mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
        mt[kk]=mt[kk+M]^(y>>1)^mag01[y&0x1UL];
      }
      for(;kk<N-1;kk++)
      {
        y=(mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
        mt[kk]=mt[kk+(M-N)]^(y>>1)^mag01[y&0x1UL];
      }
      y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
      mt[N-1]=mt[M-1]^(y>>1)^mag01[y&0x1UL];

      mti=0;
    }

    y = mt[mti++];

    /* Tempering */
    y^=(y>>11);
    y^=(y<<7)&0x9d2c5680UL;
    y^=(y<<15)&0xefc60000UL;
    y^=(y>>18);

    return y;
  }

  /* generates a random number on [0,0x7fffffff]-interval */
  long RandomInteger(void)
  {
    return (long)(genrand_int32()>>1);
  }

  /* generates a random number on [0,1]-real-interval */
  double RandomNumber(void)
  {
    return genrand_int32()*(1.0/4294967295.0);
    /* divided by 2^32-1 */
  }

  /* generates a random number on [0,1)-real-interval */
  double RandomNumber2(void)
  {
    return genrand_int32()*(1.0/4294967296.0);
    /* divided by 2^32 */
  }

  /* generates a random number on (0,1)-real-interval */
  double RandomNumber3(void)
  {
    return (((double)genrand_int32())+0.5)*(1.0/4294967296.0);
    /* divided by 2^32 */
  }

  /* generates a random number on [0,1) with 53-bit resolution*/
  double RandomNumber4(void)
  {
    unsigned long a=genrand_int32()>>5,b=genrand_int32()>>6;
    return(a*67108864.0+b)*(1.0/9007199254740992.0);
  }

  void StoreRandomNumberStatus(void)
  {
    int i;

    for(i=0;i<N;i++)
      mt_bak[i]=mt[i];
    mti_bak=mti;
    mag01_bak[0]=mag01[0];
    mag01_bak[1]=mag01[1];
  }

  void RetrieveRandomNumberStatus(void)
  {
    int i;

    for(i=0;i<N;i++)
      mt[i]=mt_bak[i];
    mti=mti_bak;
    mag01[0]=mag01_bak[0];
    mag01[1]=mag01_bak[1];
  }

  void WriteRestartUtils(FILE *FilePtr)
  {
    fwrite(mt,N,sizeof(unsigned long),FilePtr);
    fwrite(&mti,1,sizeof(int),FilePtr);
    fwrite(mag01,1,sizeof(unsigned long[2]),FilePtr);
  }

  void ReadRestartUtils(FILE *FilePtr)
  {
    fread(mt,N,sizeof(unsigned long),FilePtr);
    fread(&mti,1,sizeof(int),FilePtr);
    fread(mag01,1,sizeof(unsigned long[2]),FilePtr);
  }

  #undef N 624
  #undef M 397
  #undef MATRIX_A
  #undef UPPER_MASK
  #undef LOWER_MASK

#endif

// The polar form of the Box-Muller transformation, See Knuth v2, 3rd ed, p122
REAL RandomGaussianNumber(void)
{
  REAL ran1,ran2,r2;

  do
  {
    ran1=2.0*RandomNumber()-1.0;
    ran2=2.0*RandomNumber()-1.0;
    r2=SQR(ran1)+SQR(ran2);
  }
  while((r2>1.0)||(r2==0.0));
  return ran2*sqrt(-2.0*log(r2)/r2);
}

REAL RandomGaussianNumberMeanVariance(REAL mean,REAL sigma)
{
  REAL ran1,ran2,r2;

  do
  {
    ran1=2.0*RandomNumber()-1.0;
    ran2=2.0*RandomNumber()-1.0;
    r2=SQR(ran1)+SQR(ran2);
  }
  while((r2>1.0)||(r2==0.0));
  return mean+sigma*ran2*sqrt(-2.0*log(r2)/r2);
}

// G. Marsaglia
// "Choosing a Point from the Surface of a Sphere"
// Ann. Math. Stat., Volume 43, page 645-646, 1972

VECTOR RandomNumberOnUnitSphere(void)
{
  REAL ran1,ran2,ranh,ransq;
  VECTOR v;

  do
  {
    ran1=2.0*RandomNumber()-1.0;
    ran2=2.0*RandomNumber()-1.0;
    ransq=SQR(ran1)+SQR(ran2);
  }
  while(ransq>=1.0);

  ranh=2.0*sqrt(1.0-ransq);
  v.x=ran1*ranh;
  v.y=ran2*ranh;
  v.z=1.0-2.0*ransq;
  return v;
}


VECTOR RandomNumberOnCone(VECTOR v,REAL theta)
{
  REAL xprev,yprev,zprev,ran1,cosin;

  cosin=cos(theta);
  ran1=1.0/sqrt(SQR(v.x)+SQR(v.y)+SQR(v.z));
  xprev=v.x*ran1;
  yprev=v.y*ran1;
  zprev=v.z*ran1;

  v=RandomNumberOnUnitSphere();

  ran1=v.x*xprev+v.y*yprev+v.z*zprev;
  v.x=v.x-ran1*xprev;
  v.y=v.y-ran1*yprev;
  v.z=v.z-ran1*zprev;
  ran1=sqrt((1.0-cosin*cosin)/(SQR(v.x)+SQR(v.y)+SQR(v.z)));
  v.x=v.x*ran1+cosin*xprev;
  v.y=v.y*ran1+cosin*yprev;
  v.z=v.z*ran1+cosin*zprev;
  ran1=1.0/sqrt(SQR(v.x)+SQR(v.y)+SQR(v.z));
  v.x=v.x*ran1;
  v.y=v.y*ran1;
  v.z=v.z*ran1;
  return v;
}

// Calculates the matrix from the rotation of source a to image b
// a and b are untouched by this procedure. At the end of this
// procedure, the matrix b is recovered by:
// b[0]=a[0]*rot[0][0]+a[1]*rot[0][1]+a[2]*rot[0][2]
// b[1]=a[0]*rot[1][0]+a[1]*rot[1][1]+a[2]*rot[1][2]
// b[2]=a[0]*rot[2][0]+a[1]*rot[2][1]+a[2]*rot[2][2]
// a=source vector
// b=target vector
// rot=rotation matrix
// note:
// 1) no need to normalize a and b
// 2) beware of the singularity when a==b !!!!
void CalculateRotationMatrix(VECTOR a,VECTOR b,REAL_MATRIX3x3 *rot)
{
  REAL x,y,z,w,s,c,nn;

  nn=(a.x*b.x+a.y*b.y+a.z*b.z)/sqrt((a.x*a.x+a.y*a.y+
      a.z*a.z)*(b.x*b.x+b.y*b.y+b.z*b.z));

  if(nn>0.999999)
  {
    rot->ax=1.0;  rot->bx=0.0;  rot->cx=0.0;
    rot->ay=0.0;  rot->by=1.0;  rot->cy=0.0;
    rot->az=0.0;  rot->bz=0.0;  rot->cz=1.0;
  }
  else if (nn<-0.999999)
  {
    rot->ax=-1.0; rot->bx=0.0;  rot->cx=0.0;
    rot->ay=0.0;  rot->by=-1.0; rot->cy=0.0;
    rot->az=0.0;  rot->bz=0.0;  rot->cz=-1.0;
  }
  else
  {
    x=a.z*b.y-a.y*b.z;
    y=a.x*b.z-a.z*b.x;
    z=a.y*b.x-a.x*b.y;
    s=1.0/sqrt(x*x+y*y+z*z);
    x=x*s;
    y=y*s;
    z=z*s;

    c=(a.x*b.x+a.y*b.y+a.z*b.z)/sqrt((a.x*a.x+a.y*a.y+
       a.z*a.z)*(b.x*b.x+b.y*b.y+b.z*b.z));

    w=1.0-c;
    s=sqrt(1.0-c*c);

    rot->ax=x*x*w+c;    rot->bx=x*y*w+z*s;  rot->cx=x*z*w-y*s;
    rot->ay=x*y*w-z*s;  rot->by=y*y*w+c;    rot->cy=y*z*w+x*s;
    rot->az=x*z*w+y*s;  rot->bz=y*z*w-x*s;  rot->cz=z*z*w+c;
  }
}

VECTOR Perpendicular(VECTOR a,VECTOR b)
{
  VECTOR v,fac;
  int choise;

  fac.x=fabs(b.y*a.z-a.y*b.z);
  fac.y=fabs(b.x*a.z-a.x*b.z);
  fac.z=fabs(b.y*a.x-a.y*b.x);
  choise=1;

  if(fac.y>fac.x)
  {
    choise=2;
    fac.x=fac.y;
  }

  if(fac.z>fac.x) choise=3;

  if(choise==1)
  {
    v.x=1.0;
    v.y=(a.x*b.z-b.x*a.z)/(b.y*a.z-a.y*b.z);
    v.z=(b.x*a.y-a.x*b.y)/(b.y*a.z-a.y*b.z);
  }
  else if (choise==2)
  {
    v.x=(a.y*b.z-b.y*a.z)/(b.x*a.z-a.x*b.z);
    v.y=1.0;
    v.z=(b.y*a.x-a.y*b.x)/(b.x*a.z-a.x*b.z);
  }
  else
  {
    v.x=(b.z*a.y-a.z*b.y)/(b.y*a.x-a.y*b.x);
    v.y=(a.z*b.x-b.z*a.x)/(b.y*a.x-a.y*b.x);
    v.z=1.0;
  }

  fac.x=1.0/sqrt(SQR(v.x)+SQR(v.y)+SQR(v.z));
  v.x*=fac.x;
  v.y*=fac.x;
  v.z*=fac.x;
  return v;
}


// Random-rotation matrix Shoemake's method
void RandomArrayRotationMatrix(VECTOR *Cord,int n)
{
  int i;
  REAL R[3][3];
  POINT p;
  REAL R1,R2;
  REAL X0,Y1,Y2;
  REAL U0,U1,U2,U3;
  REAL COEFI,COEFUU,COEFE;

  X0=RandomNumber();
  Y1=2.0*M_PI*RandomNumber();
  Y2=2.0*M_PI*RandomNumber();
  R1=sqrt(1.0-X0);
  R2=sqrt(X0);
  U0=cos(Y2)*R2;
  U1=sin(Y1)*R1;
  U2=cos(Y1)*R1;
  U3=sin(Y2)*R2;
  COEFI=2.0*U0*U0-1.0;
  COEFUU=2.0;
  COEFE=2.0*U0;
  R[0][0]=COEFI+COEFUU*U1*U1;
  R[1][1]=COEFI+COEFUU*U2*U2;
  R[2][2]=COEFI+COEFUU*U3*U3;
  R[1][2]=COEFUU*U2*U3-COEFE*U1;
  R[2][0]=COEFUU*U3*U1-COEFE*U2;
  R[0][1]=COEFUU*U1*U2-COEFE*U3;
  R[2][1]=COEFUU*U3*U2+COEFE*U1;
  R[0][2]=COEFUU*U1*U3+COEFE*U2;
  R[1][0]=COEFUU*U2*U1+COEFE*U3;

  for(i=0;i<n;i++)
  {
    p.x=Cord[i].x*R[0][0]+Cord[i].y*R[0][1]+Cord[i].z*R[0][2];
    p.y=Cord[i].x*R[1][0]+Cord[i].y*R[1][1]+Cord[i].z*R[1][2];
    p.z=Cord[i].x*R[2][0]+Cord[i].y*R[2][1]+Cord[i].z*R[2][2];
    Cord[i]=p;
  }
}

void RotationAroundXAxis(VECTOR *Cord,int n,REAL theta)
{
  int i;
  REAL w,s,c,rot[3][3];

  c=cos(theta);
  s=sin(theta);

  rot[0][0]=1.0; rot[1][0]=0.0;  rot[2][0]=0.0;
  rot[0][1]=0.0; rot[1][1]=c;    rot[2][1]=-s;
  rot[0][2]=0.0; rot[1][2]=s;    rot[2][2]=c;

  for(i=0;i<n;i++)
  {
    w=Cord[i].x*rot[0][0]+Cord[i].y*rot[0][1]+Cord[i].z*rot[0][2];
    s=Cord[i].x*rot[1][0]+Cord[i].y*rot[1][1]+Cord[i].z*rot[1][2];
    c=Cord[i].x*rot[2][0]+Cord[i].y*rot[2][1]+Cord[i].z*rot[2][2];
    Cord[i].x=w;
    Cord[i].y=s;
    Cord[i].z=c;
  }
}

void RotationAroundYAxis(VECTOR *Cord,int n,REAL theta)
{
  int i;
  REAL w,s,c,rot[3][3];

  c=cos(theta);
  s=sin(theta);

  rot[0][0]=c;   rot[1][0]=0;    rot[2][0]=s;
  rot[0][1]=0;   rot[1][1]=1.0;  rot[2][1]=0;
  rot[0][2]=-s;  rot[1][2]=0;    rot[2][2]=c;

  for(i=0;i<n;i++)
  {
    w=Cord[i].x*rot[0][0]+Cord[i].y*rot[0][1]+Cord[i].z*rot[0][2];
    s=Cord[i].x*rot[1][0]+Cord[i].y*rot[1][1]+Cord[i].z*rot[1][2];
    c=Cord[i].x*rot[2][0]+Cord[i].y*rot[2][1]+Cord[i].z*rot[2][2];
    Cord[i].x=w;
    Cord[i].y=s;
    Cord[i].z=c;
  }
}

void RotationAroundZAxis(VECTOR *Cord,int n,REAL theta)
{
  int i;
  REAL w,s,c,rot[3][3];

  c=cos(theta);
  s=sin(theta);

  rot[0][0]=c;   rot[1][0]=-s;   rot[2][0]=0;
  rot[0][1]=s;   rot[1][1]=c;    rot[2][1]=0;
  rot[0][2]=0;   rot[1][2]=0;    rot[2][2]=1.0;

  for(i=0;i<n;i++)
  {
    w=Cord[i].x*rot[0][0]+Cord[i].y*rot[0][1]+Cord[i].z*rot[0][2];
    s=Cord[i].x*rot[1][0]+Cord[i].y*rot[1][1]+Cord[i].z*rot[1][2];
    c=Cord[i].x*rot[2][0]+Cord[i].y*rot[2][1]+Cord[i].z*rot[2][2];
    Cord[i].x=w;
    Cord[i].y=s;
    Cord[i].z=c;
  }
}

void RotationAroundXYZAxis(VECTOR v,VECTOR *Cord,int n,REAL theta)
{
  int i;
  REAL w,s,c,rot[3][3];

  c=cos(theta);
  w=1.0-c;
  s=sqrt(1.0-c*c);

  if(theta<0.0) s=-s;

  rot[0][0]=(v.x)*(v.x)*w+c;
  rot[0][1]=(v.x)*(v.y)*w+(v.z)*s;
  rot[0][2]=(v.x)*(v.z)*w-(v.y)*s;
  rot[1][0]=(v.x)*(v.y)*w-(v.z)*s;
  rot[1][1]=(v.y)*(v.y)*w+c;
  rot[1][2]=(v.y)*(v.z)*w+(v.x)*s;
  rot[2][0]=(v.x)*(v.z)*w+(v.y)*s;
  rot[2][1]=(v.y)*(v.z)*w-(v.x)*s;
  rot[2][2]=(v.z)*(v.z)*w+c;

  for(i=0;i<n;i++)
  {
    w=Cord[i].x*rot[0][0]+Cord[i].y*rot[0][1]+Cord[i].z*rot[0][2];
    s=Cord[i].x*rot[1][0]+Cord[i].y*rot[1][1]+Cord[i].z*rot[1][2];
    c=Cord[i].x*rot[2][0]+Cord[i].y*rot[2][1]+Cord[i].z*rot[2][2];
    Cord[i].x=w;
    Cord[i].y=s;
    Cord[i].z=c;
  }
}

// by Glenn Murray: 'Rotation About an Arbitrary Axis in 3 Dimensions'
// http://inside.mines.edu/~gmurray/ArbitraryAxisRotation/
// tested with: http://twist-and-shout.appspot.com
void RotationArrayAroundXYZAxisAtPoint(POINT origin,VECTOR v,POINT *Cord,int n,REAL theta)
{
  int i;
  REAL s,c;
  VECTOR p;

  c=cos(theta);
  s=sin(theta);

  for(i=0;i<n;i++)
  {
    p=Cord[i];
    Cord[i].x=(origin.x*(SQR(v.y)+SQR(v.z))-v.x*(origin.y*v.y+origin.z*v.z-v.x*p.x-v.y*p.y-v.z*p.z))*(1.0-c)+p.x*c+(-origin.z*v.y+origin.y*v.z-v.z*p.y+v.y*p.z)*s;
    Cord[i].y=(origin.y*(SQR(v.x)+SQR(v.z))-v.y*(origin.x*v.x+origin.z*v.z-v.x*p.x-v.y*p.y-v.z*p.z))*(1.0-c)+p.y*c+( origin.z*v.x-origin.x*v.z+v.z*p.x-v.x*p.z)*s;
    Cord[i].z=(origin.z*(SQR(v.x)+SQR(v.y))-v.z*(origin.x*v.x+origin.y*v.y-v.x*p.x-v.y*p.y-v.z*p.z))*(1.0-c)+p.z*c+(-origin.y*v.x+origin.x*v.y-v.y*p.x+v.x*p.y)*s;
  }
}

VECTOR RotationAroundXYZAxisAtPoint(POINT origin,VECTOR v,POINT p,REAL theta)
{
  REAL s,c;
  VECTOR Cord;

  c=cos(theta);
  s=sin(theta);

  Cord.x=(origin.x*(SQR(v.y)+SQR(v.z))-v.x*(origin.y*v.y+origin.z*v.z-v.x*p.x-v.y*p.y-v.z*p.z))*(1.0-c)+p.x*c+(-origin.z*v.y+origin.y*v.z-v.z*p.y+v.y*p.z)*s;
  Cord.y=(origin.y*(SQR(v.x)+SQR(v.z))-v.y*(origin.x*v.x+origin.z*v.z-v.x*p.x-v.y*p.y-v.z*p.z))*(1.0-c)+p.y*c+( origin.z*v.x-origin.x*v.z+v.z*p.x-v.x*p.z)*s;
  Cord.z=(origin.z*(SQR(v.x)+SQR(v.y))-v.z*(origin.x*v.x+origin.y*v.y-v.x*p.x-v.y*p.y-v.z*p.z))*(1.0-c)+p.z*c+(-origin.y*v.x+origin.x*v.y-v.y*p.x+v.x*p.y)*s;
  return Cord;
}

VECTOR TransformMapping(REAL_MATRIX3x3 m,VECTOR t)
{
  VECTOR dr;

  //t.x-=BarrierPosition.x;
  //t.y-=BarrierPosition.y;
  //t.z-=BarrierPosition.z;
  dr.x=m.ax*t.x+m.bx*t.y+m.cx*t.z;
  dr.y=m.ay*t.x+m.by*t.y+m.cy*t.z;
  dr.z=m.az*t.x+m.bz*t.y+m.cz*t.z;
  //dr.x+=BarrierPosition.x;
  //dr.y+=BarrierPosition.y;
  //dr.z+=BarrierPosition.z;
  return dr;
}


// angle is counterclockwise
VECTOR RotateVectorAboutX(VECTOR vec,REAL angle)
{
  VECTOR c;

  c.x=vec.x;
  c.y=vec.y*cos(angle)-vec.z*sin(angle);
  c.z=vec.z*cos(angle)+vec.y*sin(angle);
  return c;
}

VECTOR RotateVectorAboutY(VECTOR vec,REAL angle)
{
  VECTOR c;

  c.x=vec.x*cos(angle)-vec.z*sin(angle);
  c.y=vec.y;
  c.z=vec.z*cos(angle)+vec.x*sin(angle);
  return c;
}

VECTOR RotateVectorAboutZ(VECTOR vec,REAL angle)
{
  VECTOR c;

  c.x=vec.x*cos(angle)-vec.y*sin(angle);
  c.y=vec.y*cos(angle)+vec.x*sin(angle);
  c.z=vec.z;
  return c;
}

VECTOR GenerateRandomCylinderPoint(int dir,VECTOR origin,VECTOR rotation,REAL radius)
{
  VECTOR c;

  c.x=c.y=c.z=0.0;
  switch(dir)
  {
    case X_DIR:
      c.x=UnitCellSize[0].x*(RandomNumber()-0.5);
      do
      {
        c.y=radius*2.0*(RandomNumber()-0.5);
        c.z=radius*2.0*(RandomNumber()-0.5);
      }while(sqrt(SQR(c.y)+SQR(c.z))>radius);
      break;
    case Y_DIR:
      c.y=UnitCellSize[0].y*(RandomNumber()-0.5);
      do
      {
        c.x=radius*2.0*(RandomNumber()-0.5);
        c.z=radius*2.0*(RandomNumber()-0.5);
      }while(sqrt(SQR(c.x)+SQR(c.z))>radius);
      break;
    case Z_DIR:
      c.z=UnitCellSize[0].z*(RandomNumber()-0.5);
      do
      {
        c.x=radius*2.0*(RandomNumber()-0.5);
        c.y=radius*2.0*(RandomNumber()-0.5);
      }while(sqrt(SQR(c.x)+SQR(c.y))>radius);
      break;
  }
  c=RotateVectorAboutX(c,rotation.x);
  c=RotateVectorAboutY(c,rotation.y);
  c=RotateVectorAboutZ(c,rotation.z);
  c.x+=origin.x;
  c.y+=origin.y;
  c.z+=origin.z;
  return c;
}

// Ref. K. Shoemake, "Uniform Random Rotations", in D. Kirk, editor, graphic Gems III, pages 124-132, Academic Press, New York, 1992
QUATERNION RandomQuarternion(void)
{
  REAL s;
  REAL sigma1,sigma2;
  REAL theta1,theta2;
  QUATERNION q;

  s=RandomNumber();
  sigma1=sqrt(1.0-s);
  sigma2=sqrt(s);
  theta1=2.0*M_PI*RandomNumber();
  theta2=2.0*M_PI*RandomNumber();

  q.r=sigma2*cos(theta2);
  q.i=sigma1*sin(theta1);
  q.j=sigma1*cos(theta1);
  q.k=sigma2*sin(theta2);
  return q;
}

QUATERNION MultiplyQuarternions(QUATERNION q1,QUATERNION q2)
{
  QUATERNION q;

  q.r=q1.r*q2.r-q1.i*q2.i-q1.j*q2.j-q1.k*q2.k;
  q.i=q1.r*q2.i+q1.i*q2.r+q1.j*q2.k-q1.k*q2.j;
  q.j=q1.r*q2.j-q1.i*q2.k+q1.j*q2.r+q1.k*q2.i;
  q.k=q1.r*q2.k+q1.i*q2.j-q1.j*q2.i+q1.k*q2.r;
  return q;
}

QUATERNION ScaleQuarternion(QUATERNION q1,REAL s)
{
  QUATERNION q;

  q.r=s*q1.r;
  q.i=s*q1.i;
  q.j=s*q1.j;
  q.k=s*q1.k;
  return q;
}

void EulerToQuaternions(QUATERNION *q,VECTOR Ang)
{
  REAL t1,t2,t3;
  t1=0.5*Ang.y;
  t2=0.5*(Ang.x-Ang.z);
  t3=0.5*(Ang.x+Ang.z);
  q->r=cos(t1)*cos(t3);
  q->i=sin(t1)*cos(t2);
  q->j=sin(t1)*sin(t2);
  q->k=cos(t1)*sin(t3);
}


//
// The following prototype declarations need to be added to your code (ie. header file).
//

/*-------------------- Global Function Description Block ----------------------
 *
 *     ***QUINTIC************************************************25.03.98
 *     Solution of a quintic equation by a hybrid method:
 *     first real solution is obtained numerically by the Newton method,
 *     the remaining four roots are obtained analytically by QUARTIC
 *     NO WARRANTY, ALWAYS TEST THIS SUBROUTINE AFTER DOWNLOADING
 *     ******************************************************************
 *     dd(0:4)     (i)  vector containing the polynomial coefficients
 *     sol(1:4)    (o)  results, real part
 *     soli(1:4)   (o)  results, imaginary part
 *     Nsol        (o)  number of real solutions
 *  17-Oct-2004 / Raoul Rausch
 *    Conversion from Fortran to C
 *
 *
 *-----------------------------------------------------------------------------
 */
int quintic(REAL dd[6], REAL sol[5], REAL soli[5], int *Nsol, REAL xstart)
{
  REAL  dd4[5], sol4[4], soli4[4], xnew, xs;//, soli4[4];/*dd[6], sol[5], soli[5],*/
  REAL sum, sum1, eC;
  const REAL eps = 1.e-8;
  int i, Nsol4;

  *Nsol = 0;

  printf("\n Quintic!\n");

  if (dd[5] == 0.0)
  {
    printf("\n ERROR: NOT A QUINTIC EQUATION");
    return 0;
  }

  // Newton iteration of one real root
  xs= xstart;
  xnew = xstart;  //added rr
  do
  {
    xs = xnew;  //added rr
    sum = dd[0];
    for (i=1;i<6;i++)  sum += dd[i]*pow(xs,i);  // Don't know what ** means
    sum1 = dd[1];
    for (i=1;i<5;i++)  sum1 += (REAL)(i+1)*dd[i+1]*pow(xs,i);
    xnew = xs - sum/sum1;
    //if (fabs(xnew-xs) > eps)
    //xs =xnew;
  }while (fabs(xnew-xs) > eps);

  eC = xnew;
  //
  // "eC" is one real root of quintic equation
  // reduce quintic to quartic equation using "eC"
  dd4[4] = dd[5];
  for (i=4;i>0;i--)  dd4[i-1] = dd[i] + eC*dd4[i];

  quartic(dd4, sol4, soli4, &Nsol4);


  sol[0] = eC;
  soli[0] = 0.0;

  for (i=0;i<4;i++)
  {
    sol[i+1] =sol4[i];
    soli[i+1] = soli4[i];
  }
  *Nsol = Nsol4 + 1;

  return 0;
}

/*-------------------- Global Function Description Block ----------------------
 *
 *     ***QUARTIC************************************************25.03.98
 *     Solution of a quartic equation
 *     ref.: J. E. Hacke, Amer. Math. Monthly, Vol. 48, 327-328, (1941)
 *     NO WARRANTY, ALWAYS TEST THIS SUBROUTINE AFTER DOWNLOADING
 *     ******************************************************************
 *     dd(0:4)     (i)  vector containing the polynomial coefficients
 *     sol(1:4)    (o)  results, real part
 *     soli(1:4)   (o)  results, imaginary part
 *     Nsol        (o)  number of real solutions
 *     ==================================================================
 *    17-Oct-2004 / Raoul Rausch
 *    Conversion from Fortran to C
 *
 *
 *-----------------------------------------------------------------------------
 */
 int quartic(REAL dd[5], REAL sol[4], REAL soli[4], int* Nsol)
 {
  REAL AA[4], z[3];
  REAL a, b, c, d, f, p, q, r, zsol, xK2, xL, xK, sqp, sqm;
  int ncube, i;
  *Nsol = 0;

  if (dd[4] == 0.0)
  {
    printf("\n ERROR: NOT A QUARTIC EQUATION");
    return 0;
  }

  a = dd[4];
  b = dd[3];
  c = dd[2];
  d = dd[1];
  f = dd[0];

  p = (-3.0*pow(b,2) + 8.0 *a*c)/(8.0*pow(a,2));
  q = (pow(b,3) - 4.0*a*b*c + 8.0 *d*pow(a,2)) / (8.0*pow(a,3));
  r = (-3.0*pow(b,4) + 16.0 *a*pow(b,2)*c - 64.0 *pow(a,2)*b*d + 256.0 *pow(a,3)*f)/(256.0*pow(a,4));

  // Solve cubic resolvent
  AA[3] = 8.0;
  AA[2] = -4.0*p;
  AA[1] = -8.0*r;
  AA[0] = 4.0*p*r - pow(q,2);

  //printf("\n bcubic %.4e\t%.4e\t%.4e\t%.4e ", AA[0], AA[1], AA[2], AA[3]);
  cubic(AA, z, &ncube);
  //printf("\n acubic %.4e\t%.4e\t%.4e ", z[0], z[1], z[2]);

  zsol = - 1.e99;
  for (i=0;i<ncube;i++)
          zsol = MAX2(zsol, z[i]);  //Not sure C has max fct
  z[0] =zsol;
  xK2 = 2.0*z[0] -p;
  xK = sqrt(xK2);
  xL = q/(2.0*xK);
  sqp = xK2 - 4.0 * (z[0] + xL);
  sqm = xK2 - 4.0 * (z[0] - xL);

  for (i=0;i<4;i++)  soli[i] = 0.0;
  if ( (sqp >= 0.0) && (sqm >= 0.0))
  {
    printf("\n case 1 ");
    sol[0] = 0.5 * (xK + sqrt(sqp));
    sol[1] = 0.5 * (xK - sqrt(sqp));
    sol[2] = 0.5 * (-xK + sqrt(sqm));
    sol[3] = 0.5 * (-xK - sqrt(sqm));
    *Nsol = 4;
  }
  else if ( (sqp >= 0.0) && (sqm < 0.0))
  {
    printf("\n case 2 ");
    sol[0] = 0.5 * (xK + sqrt(sqp));
    sol[1] = 0.5 * (xK - sqrt(sqp));
    sol[2] = -0.5 * xK;
    sol[3] = -0.5 * xK;
    soli[2] =  sqrt(-.25 * sqm);
    soli[3] = -sqrt(-.25 * sqm);
    *Nsol = 2;
  }
  else if ( (sqp < 0.0) && (sqm >= 0.0))
  {
    printf("\n case 3 ");
    sol[0] = 0.5 * (-xK + sqrt(sqm));
    sol[1] = 0.5 * (-xK - sqrt(sqm));
    sol[2] = 0.5 * xK;
    sol[3] = 0.5 * xK;
    soli[2] =  sqrt(-0.25 * sqp);
    soli[3] = -sqrt(-0.25 * sqp);
    *Nsol = 2;
  }
  else if ( (sqp < 0.0) && (sqm < 0.0))
  {
    printf("\n case 4 ");
    sol[0] = -0.5 * xK;
    sol[1] = -0.5 * xK;
    soli[0] =  sqrt(-0.25 * sqm);
    soli[1] = -sqrt(-0.25 * sqm);
    sol[2] = 0.5 * xK;
    sol[3] = 0.5 * xK;
    soli[2] =  sqrt(-0.25 * sqp);
    soli[3] = -sqrt(-0.25 * sqp);
    *Nsol = 0;
  }

  for (i=0;i<4;i++)  sol[i] -= b/(4.0*a);
  return 0;

 }

 /*-------------------- Global Function Description Block ----------------------
  *
  *     ***CUBIC************************************************08.11.1986
  *     Solution of a cubic equation
  *     Equations of lesser degree are solved by the appropriate formulas.
  *     The solutions are arranged in ascending order.
  *     NO WARRANTY, ALWAYS TEST THIS SUBROUTINE AFTER DOWNLOADING
  *     ******************************************************************
  *     A(0:3)      (i)  vector containing the polynomial coefficients
  *     X(1:L)      (o)  results
  *     L           (o)  number of valid solutions (beginning with X(1))
  *     ==================================================================
  *    17-Oct-2004 / Raoul Rausch
  *    Conversion from Fortran to C
  *
  *
  *-----------------------------------------------------------------------------
  */
int cubic(REAL A[4], REAL X[3], int* L)
{
  const REAL PI = 3.1415926535897932;
  const REAL THIRD = 1./3.;
  REAL U[3],W, P, Q, DIS, PHI;
  int i;

  //define cubic root as statement function
  // In C, the function is defined outside of the cubic fct

  // ====determine the degree of the polynomial ====

  if (A[3] != 0.0)
  {
    //cubic problem
    W = A[2]/A[3]*THIRD;
    P = pow((A[1]/A[3]*THIRD - pow(W,2)),3);
    Q = -.5*(2.0*pow(W,3)-(A[1]*W-A[0])/A[3] );
    DIS = pow(Q,2)+P;
    if ( DIS < 0.0 )
    {
      //three real solutions!
      //Confine the argument of ACOS to the interval [-1;1]!
      PHI = acos(MIN2((REAL)1.0,MAX2((REAL)-1.0,Q/sqrt(-P))));
      P=2.0*pow((-P),(5.e-1*THIRD));
      for (i=0;i<3;i++)  U[i] = P*cos((PHI+2.0*((REAL)i)*PI)*THIRD)-W;
      X[0] = MIN2(U[0], MIN2(U[1], U[2]));
      X[1] = MAX2(MIN2(U[0], U[1]),MAX2( MIN2(U[0], U[2]), MIN2(U[1], U[2])));
      X[2] = MAX2(U[0], MAX2(U[1], U[2]));
      *L = 3;
    }
    else
    {
      // only one real solution!
      DIS = sqrt(DIS);
      X[0] = CBRT(Q+DIS)+CBRT(Q-DIS)-W;
      *L=1;
    }
  }
  else if (A[2] != 0.0)
  {
    // quadratic problem
    P = 0.5*A[1]/A[2];
    DIS = pow(P,2)-A[0]/A[2];
    if (DIS > 0.0)
    {
      // 2 real solutions
      X[0] = -P - sqrt(DIS);
      X[1] = -P + sqrt(DIS);
      *L=2;
    }
    else
    {
      // no real solution
      *L=0;
    }
  }
  else if (A[1] != 0.0)
  {
    //linear equation
    X[0] =A[0]/A[1];
    *L=1;
  }
  else
  {
    //no equation
    *L=0;
  }
 /*
  *     ==== perform one step of a newton iteration in order to minimize
  *          round-off errors ====
  */
  for (i=0;i<*L;i++)
  {
    X[i] = X[i] - (A[0]+X[i]*(A[1]+X[i]*(A[2]+X[i]*A[3])))/(A[1]+X[i]*(2.0*A[2]+X[i]*3.0*A[3]));
  //  printf("\n X inside cubic %.15e\n", X[i]);
  }

  return 0;
}

int signR(REAL Z)
{
  int ret;

  if (Z > 0.0)  ret = 1;
  else if (Z < 0.0)  ret = -1;
  else ret =0;

  return ret;
}

REAL CBRT(REAL Z)
{
  REAL ret;
  const REAL THIRD = 1./3.;
  //define cubic root as statement function
  //SIGN has different meanings in both C and Fortran
  // Was unable to use the sign command of C, so wrote my own
  // that why a new variable needs to be introduced that keeps track of the sign of
  // SIGN is supposed to return a 1, -1 or 0 depending on what the sign of the argument is
  ret = (REAL)fabs(pow((REAL)fabs(Z),(REAL)THIRD)) * (REAL)signR((REAL)Z);
  return ret;
}


//#define ALF 1.0e-4
//#define TOLX 1.0e-7
#define ALF 1.0e-8
#define TOLX 1.0e-14
void LineSearch(int np,int nb,int n,REAL *xold,REAL fold,REAL *g,REAL *p,REAL *x,
          REAL *f,REAL stpmax,int *check,REAL (*func)(int,int,REAL []))
{
  int i;
  REAL a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2,slope,sum,temp,test,tmplam;

  fold=0.0;
  fold2=0.0;
  f2=0.0;
  alam2=0.0;

  *check=0;
  for(sum=0.0,i=0;i<n;i++)
    sum+=p[i]*p[i];
  sum=sqrt(sum);
  if(sum>stpmax)
    for(i=0;i<n;i++)
      p[i]*=stpmax/sum;
  for(slope=0.0,i=0;i<n;i++)
    slope+=g[i]*p[i];
  test=0.0;
  for(i=0;i<n;i++)
  {
    temp=fabs(p[i])/MAX2(fabs(xold[i]),(REAL)1.0);
    if(temp>test) test=temp;
  }
  alamin=TOLX/test;
  alam=1.0;
  for(;;)
  {
    for(i=0;i<n;i++)
      x[i]=xold[i]+alam*p[i];
    *f=(*func)(np,nb,x);
    if(alam<alamin)
    {
      for(i=0;i<n;i++)
        x[i]=xold[i];
      *check=1;
      return;
    }
    else if(*f<=fold+ALF*alam*slope) return;
    else
    {
      if(alam==1.0)
        tmplam = -slope/(2.0*(*f-fold-slope));
      else
      {
        rhs1 = *f-fold-alam*slope;
        rhs2=f2-fold2-alam2*slope;
        a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
        b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
        if(a==0.0) tmplam=-slope/(2.0*b);
        else
        {
          disc=b*b-3.0*a*slope;
          if(disc<0.0)
          {
            fprintf(stderr, "Roundoff problem in lnsrch.");
            exit(0);
          }
          else tmplam=(-b+sqrt(disc))/(3.0*a);
        }
        if (tmplam>0.5*alam)
          tmplam=0.5*alam;
      }
    }
    alam2=alam;
    f2 = *f;
    fold2=fold;
    alam=MAX2(tmplam,0.1*alam);
  }
}
#undef ALF
#undef TOLX
#undef NRANSI

#define ITMAX 100
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

REAL brent(int np,int nb,REAL ax, REAL bx, REAL cx, REAL (*f)(int,int,REAL), REAL tol,
  REAL *xmin)
{
  int iter;
  REAL a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  REAL e=0.0;

  d=0.0;
  a=(ax<cx?ax:cx);
  b=(ax>cx?ax:cx);
  x=w=v=bx;
  fw=fv=fx=(*f)(np,nb,x);
  for(iter=0;iter<ITMAX;iter++)
  {
    xm=0.5*(a+b);
    tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
    if(fabs(x-xm)<=(tol2-0.5*(b-a)))
    {
      *xmin=x;
      return fx;
    }
    if(fabs(e)>tol1)
    {
      r=(x-w)*(fx-fv);
      q=(x-v)*(fx-fw);
      p=(x-v)*q-(x-w)*r;
      q=2.0*(q-r);
      if(q>0.0) p=-p;
      q=fabs(q);
      etemp=e;
      e=d;
      if(fabs(p)>=fabs(0.5*q*etemp)||p<=q*(a-x)||p>=q*(b-x))
        d=CGOLD*(e=(x >= xm ? a-x : b-x));
      else
      {
        d=p/q;
        u=x+d;
        if((u-a<tol2)||(b-u<tol2))
          d=SIGN(tol1,xm-x);
      }
    }
    else
    {
      d=CGOLD*(e=(x >= xm ? a-x : b-x));
    }
    u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
    fu=(*f)(np,nb,u);
    if(fu<=fx)
    {
      if (u >= x) a=x; else b=x;
      SHFT(v,w,x,u)
      SHFT(fv,fw,fx,fu)
    }
    else
    {
      if(u<x) a=u;
      else b=u;
      if((fu<=fw)||(w==x))
      {
        v=w;
        w=u;
        fv=fw;
        fw=fu;
      }
      else if((fu<=fv)||(v==x)||(v==w))
      {
        v=u;
        fv=fu;
      }
    }
  }
  fprintf(stderr, "Too many iterations in brent");
  *xmin=x;
  return fx;
}
#undef ITMAX
#undef CGOLD
#undef ZEPS
#undef SHFT

#define ITMAX 100
#define ZEPS 1.0e-10
#define MOV3(a,b,c, d,e,f) (a)=(d);(b)=(e);(c)=(f);

REAL dbrent(int np,int nb,REAL ax, REAL bx, REAL cx, REAL (*f)(int,int,REAL),
  REAL (*df)(int,int,REAL), REAL tol, REAL *xmin)
{
  int iter,ok1,ok2;
  REAL a,b,d,d1,d2,du,dv,dw,dx,e=0.0;
  REAL fu,fv,fw,fx,olde,tol1,tol2,u,u1,u2,v,w,x,xm;

  d=0.0;
  a=(ax < cx ? ax : cx);
  b=(ax > cx ? ax : cx);
  x=w=v=bx;
  fw=fv=fx=(*f)(np,nb,x);
  dw=dv=dx=(*df)(np,nb,x);
  for(iter=0;iter<ITMAX;iter++)
  {
    xm=0.5*(a+b);
    tol1=tol*fabs(x)+ZEPS;
    tol2=2.0*tol1;
    if(fabs(x-xm)<=(tol2-0.5*(b-a)))
    {
      *xmin=x;
      return fx;
    }
    if(fabs(e)>tol1)
    {
      d1=2.0*(b-a);
      d2=d1;
      if(dw!=dx) d1=(w-x)*dx/(dx-dw);
      if(dv!=dx) d2=(v-x)*dx/(dx-dv);
      u1=x+d1;
      u2=x+d2;
      ok1=((a-u1)*(u1-b)>0.0)&&(dx*d1<=0.0);
      ok2 =((a-u2)*(u2-b)>0.0)&&(dx*d2<=0.0);
      olde=e;
      e=d;
      if(ok1||ok2)
      {
        if(ok1&&ok2)
          d=(fabs(d1) < fabs(d2) ? d1 : d2);
        else if (ok1)
          d=d1;
        else
          d=d2;
        if(fabs(d)<=fabs(0.5*olde))
        {
          u=x+d;
          if((u-a<tol2)||(b-u<tol2))
            d=SIGN(tol1,xm-x);
        }
        else
        {
          d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
        }
      }
      else
      {
        d=0.5*(e=(dx>=0.0?a-x:b-x));
      }
    }
    else
    {
      d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
    }
    if (fabs(d) >= tol1)
    {
      u=x+d;
      fu=(*f)(np,nb,u);
    }
    else
    {
      u=x+SIGN(tol1,d);
      fu=(*f)(np,nb,u);
      if(fu>fx)
      {
        *xmin=x;
        return fx;
      }
    }
    du=(*df)(np,nb,u);
    if(fu<=fx)
    {
      if(u>=x) a=x;
      else b=x;
      MOV3(v,fv,dv, w,fw,dw)
      MOV3(w,fw,dw, x,fx,dx)
      MOV3(x,fx,dx, u,fu,du)
    }
    else
    {
      if(u<x) a=u; else b=u;
      if((fu<=fw)||(w==x))
      {
        MOV3(v,fv,dv, w,fw,dw)
        MOV3(w,fw,dw, u,fu,du)
      }
      else if((fu<fv)||(v==x)||(v == w))
      {
        MOV3(v,fv,dv, u,fu,du)
      }
    }
  }
  fprintf(stderr, "Too many iterations in routine dbrent");
  return 0.0;
}
#undef ITMAX
#undef ZEPS
#undef MOV3

#define GOLD 1.618033989
#define GLIMIT 100.0
#define TINY 1.0e-20
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

REAL mnbrak(int np,int nb,REAL *ax, REAL *bx, REAL *cx, REAL *fa, REAL *fb, REAL *fc,
  REAL (*func)(int,int,REAL))
{
  REAL ulim,u,r,q,fu,dum;

  *fa=(*func)(np,nb,*ax);
  *fb=(*func)(np,nb,*bx);
  if(*fb>*fa)
  {
    SHFT(dum,*ax,*bx,dum)
    SHFT(dum,*fb,*fa,dum)
  }
  *cx=(*bx)+GOLD*(*bx-*ax);
  *fc=(*func)(np,nb,*cx);
  while(*fb>*fc)
  {
    r=(*bx-*ax)*(*fb-*fc);
    q=(*bx-*cx)*(*fb-*fa);
    u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/(2.0*SIGN(MAX2(fabs(q-r),(REAL)TINY),q-r));
    ulim=(*bx)+GLIMIT*(*cx-*bx);
    if((*bx-u)*(u-*cx)>0.0)
    {
      fu=(*func)(np,nb,u);
      if(fu<*fc)
      {
        *ax=(*bx);
        *bx=u;
        *fa=(*fb);
        *fb=fu;
        return 0.0;
      }
      else if(fu>*fb)
      {
        *cx=u;
        *fc=fu;
        return 0.0;
      }
      u=(*cx)+GOLD*(*cx-*bx);
      fu=(*func)(np,nb,u);
    }
    else if((*cx-u)*(u-ulim)>0.0)
    {
      fu=(*func)(np,nb,u);
      if(fu<*fc)
      {
        SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx))
        SHFT(*fb,*fc,fu,(*func)(np,nb,u))
      }
    }
    else if((u-ulim)*(ulim-*cx)>=0.0)
    {
      u=ulim;
      fu=(*func)(np,nb,u);
    }
    else
    {
      u=(*cx)+GOLD*(*cx-*bx);
      fu=(*func)(np,nb,u);
    }
    SHFT(*ax,*bx,*cx,u)
    SHFT(*fa,*fb,*fc,fu)
  }
  return 0.0;
}
#undef GOLD
#undef GLIMIT
#undef TINY
#undef SHFT

int ncom;
REAL *pcom,*xicom,(*nrfunc)(int,int,REAL []);
void (*nrdfun)(int,int,REAL [], REAL []);

REAL f1dim(int np,int nb,REAL x)
{
  int j;
  REAL f,*xt;

  xt=(REAL*)calloc(ncom,sizeof(REAL));
  for(j=0;j<ncom;j++)
    xt[j]=pcom[j]+x*xicom[j];
  f=(*nrfunc)(np,nb,xt);
  free(xt);
  return f;
}

REAL df1dim(int np,int nb,REAL x)
{
  int j;
  REAL df1=0.0;
  REAL *xt,*df;

  xt=(REAL*)calloc(ncom,sizeof(REAL));
  df=(REAL*)calloc(ncom,sizeof(REAL));
  for(j=0;j<ncom;j++)
    xt[j]=pcom[j]+x*xicom[j];
  (*nrdfun)(np,nb,xt,df);
  for(j=0;j<ncom;j++)
    df1 += df[j]*xicom[j];
  free(df);
  free(xt);
  return df1;
}


#define TOL 2.0e-4

void dlinmin(int np,int nb,REAL *p, REAL *xi, int n, REAL *fret, REAL (*func)(int,int,REAL []),
  void (*dfunc)(int,int,REAL [], REAL []))
{
  REAL df1dim(int,int,REAL x);
  int j;
  REAL xx,xmin,fx,fb,fa,bx,ax;

  ncom=n;
  pcom=(REAL*)calloc(n,sizeof(REAL));
  xicom=(REAL*)calloc(n,sizeof(REAL));
  nrfunc=func;
  nrdfun=dfunc;
  for(j=0;j<n;j++)
  {
    pcom[j]=p[j];
    xicom[j]=xi[j];
  }
  ax=0.0;
  xx=1.0;
  mnbrak(np,nb,&ax,&xx,&bx,&fa,&fx,&fb,f1dim);
  *fret=dbrent(np,nb,ax,xx,bx,f1dim,df1dim,TOL,&xmin);
  for(j=0;j<n;j++)
  {
    xi[j]*=xmin;
    p[j]+=xi[j];
  }
  free(xicom);
  free(pcom);
}
#undef TOL

#define TOL 2.0e-4
void linmin(int np,int nb,REAL *p, REAL *xi, int n, REAL *fret, REAL (*func)(int,int,REAL []))
{
  int j;
  REAL xx,xmin,fx,fb,fa,bx,ax;

  ncom=n;
  pcom=(REAL*)calloc(n,sizeof(REAL));
  xicom=(REAL*)calloc(n,sizeof(REAL));
  nrfunc=func;
  for(j=0;j<n;j++)
  {
    pcom[j]=p[j];
    xicom[j]=xi[j];
  }
  ax=0.0;
  xx=1.0;
  mnbrak(np,nb,&ax,&xx,&bx,&fa,&fx,&fb,f1dim);
  *fret=brent(np,nb,ax,xx,bx,f1dim,TOL,&xmin);
  for(j=0;j<n;j++)
  {
    xi[j]*=xmin;
    p[j]+=xi[j];
  }
  free(xicom);
  free(pcom);
}
#undef TOL

#define ITMAX 200

void powell(int np,int nb,REAL *p, REAL **xi, int n, REAL ftol, int *iter, REAL *fret,
  REAL (*func)(int,int,REAL []))
{
  int i,ibig,j;
  REAL del,fp,fptt,t,*pt,*ptt,*xit;

  pt=(REAL*)calloc(n,sizeof(REAL));
  ptt=(REAL*)calloc(n,sizeof(REAL));
  xit=(REAL*)calloc(n,sizeof(REAL));
  *fret=(*func)(np,nb,p);
  for(j=0;j<n;j++)
    pt[j]=p[j];
  for(*iter=0;;++(*iter))
  {
    fp=(*fret);
    ibig=0;
    del=0.0;
    for (i=0;i<n;i++)
    {
      for(j=0;j<n;j++)
        xit[j]=xi[j][i];
      fptt=(*fret);
      linmin(np,nb,p,xit,n,fret,func);
      if(fabs(fptt-(*fret))>del)
      {
        del=fabs(fptt-(*fret));
        ibig=i;
      }
    }
    if (2.0*fabs(fp-(*fret)) <= ftol*(fabs(fp)+fabs(*fret)))
    {
      free(xit);
      free(ptt);
      free(pt);
      return;
    }
    if (*iter == ITMAX) fprintf(stderr, "powell exceeding maximum iterations.");
    for(j=0;j<n;j++)
    {
      ptt[j]=2.0*p[j]-pt[j];
      xit[j]=p[j]-pt[j];
      pt[j]=p[j];
    }
    fptt=(*func)(np,nb,ptt);
    if(fptt<fp)
    {
      t=2.0*(fp-2.0*(*fret)+fptt)*SQR(fp-(*fret)-del)-del*SQR(fp-fptt);
      if(t<0.0)
      {
        linmin(np,nb,p,xit,n,fret,func);
        for(j=0;j<n;j++)
        {
          xi[j][ibig]=xi[j][n];
          xi[j][n]=xit[j];
        }
      }
    }
  }
}
#undef ITMAX

void FastFourierTransform(REAL(*fftw_data)[2], unsigned long nn, int isign)
{
  #ifdef HAVE_FFTW3

    if(isign>0)
    {
      fftw_plan plan;
      fftw_complex *result;

      result = (fftw_complex *)fftw_malloc(sizeof (fftw_complex) * nn);
      plan = fftw_plan_dft_1d(nn, fftw_data, fftw_data, FFTW_FORWARD, FFTW_ESTIMATE);
      fftw_execute(plan);
    }

  #else

    unsigned long n,mmax,m,j,istep,i;
    REAL wtemp,wr,wpr,wpi,wi,theta;
    REAL tempr,tempi,temp;
    REAL *data;

    data=(REAL*)fftw_data;


    n=nn << 1;
    j=1;
    for (i=1;i<n;i+=2)
    {
      if(j>i)
      {
        SWAP(data[j-1],data[i-1],temp);
        SWAP(data[j+1-1],data[i+1-1],temp);
      }
      m=n>>1;
      while(m>=2&&j>m)
      {
        j-=m;
        m>>=1;
      }
      j+=m;
    }
    mmax=2;
    while(n>mmax)
    {
      istep=mmax << 1;
      theta=isign*(6.28318530717959/mmax);
      wtemp=sin(0.5*theta);
      wpr = -2.0*wtemp*wtemp;
      wpi=sin(theta);
      wr=1.0;
      wi=0.0;
      for (m=1;m<mmax;m+=2)
      {
        for (i=m;i<=n;i+=istep)
        {
          j=i+mmax;
          tempr=wr*data[j-1]-wi*data[j+1-1];
          tempi=wr*data[j+1-1]+wi*data[j-1];
          data[j-1]=data[i-1]-tempr;
          data[j+1-1]=data[i+1-1]-tempi;
          data[i-1] += tempr;
          data[i+1-1] += tempi;
        }
        wr=(wtemp=wr)*wpr-wi*wpi+wr;
        wi=wi*wpr+wtemp*wpi+wi;
      }
      mmax=istep;
    }
  #endif
}

REAL_MATRIX3x3 ConvertToVoigt2D(REAL_MATRIX9x9 b)
{
  REAL_MATRIX3x3 a;

  a.ax=b.xxxx;
  a.ay=b.xxyy;
  a.az=b.xxxy;

  a.bx=b.yyxx;
  a.by=b.yyyy;
  a.bz=b.yyxy;

  a.cx=b.xyxx;
  a.cy=b.xyyy;
  a.cz=b.xyxy;
  return a;
}

REAL_MATRIX6x6 ConvertToVoigt3D(REAL_MATRIX9x9 b)
{
  REAL_MATRIX6x6 a;

  a.C11=b.xxxx;
  a.C12=b.xxyy;
  a.C13=b.xxzz;
  a.C14=b.xxyz;
  a.C15=b.xxzx;
  a.C16=b.xxxy;

  a.C21=b.yyxx;
  a.C22=b.yyyy;
  a.C23=b.yyzz;
  a.C24=b.yyyz;
  a.C25=b.yyzx;
  a.C26=b.yyxy;

  a.C31=b.zzxx;
  a.C32=b.zzyy;
  a.C33=b.zzzz;
  a.C34=b.zzyz;
  a.C35=b.zzzx;
  a.C36=b.zzxy;

  a.C41=b.yzxx;
  a.C42=b.yzyy;
  a.C43=b.yzzz;
  a.C44=b.yzyz;
  a.C45=b.yzzx;
  a.C46=b.yzxy;

  a.C51=b.zxxx;
  a.C52=b.zxyy;
  a.C53=b.zxzz;
  a.C54=b.zxyz;
  a.C55=b.zxzx;
  a.C56=b.zxxy;

  a.C61=b.xyxx;
  a.C62=b.xyyy;
  a.C63=b.xyzz;
  a.C64=b.xyyz;
  a.C65=b.xyzx;
  a.C66=b.xyxy;
  return a;
}

int CheckEnantioFace(VECTOR posA,VECTOR posB,VECTOR posC,VECTOR posD,VECTOR posE)
{
  REAL angle1,angle2;

  angle1=ReturnDihedralAngle(posA,posB,posC,posD)*RAD2DEG;
  angle2=ReturnDihedralAngle(posA,posB,posC,posE)*RAD2DEG;

  if((angle1>90)&&(angle2<0))
  {
    angle2+=360;
    if(angle1>angle2)
      return ENANTIOFACE_RE;
    else
      return ENANTIOFACE_SI;
  }
  else if ((angle2>90)&&(angle1<0))
  {
    angle1+=360;
    if (angle1>angle2)
      return ENANTIOFACE_RE;
    else
      return ENANTIOFACE_SI;
  }
  else if (angle1>angle2)
  {
    return ENANTIOFACE_RE;
  }
  else
  {
    return ENANTIOFACE_SI;
  }


/*
  ia = itfix(1,ntfix)                   //atom 1182 of framework  (Mof_Osa2)
  ib = itfix(2,ntfix)                   //atom 1094 of framework  (Mof_Mn2)
  ic = itfix(3,ntfix)                   //atom 2206 of framework  (Oxo)
  id1 = itfix(4,ntfix)                  //atom 12 of molecule        (C1_2)
  id2 = itfix(4,ntfix) - 2              //atom 10 of molecule        (C2_2)
  angle1 = geometry(ia,ib,ic,id1)       //calculate dihedral defined by 1182, 1094, 2206, 12 - this is the constrained dihedral in the minimization
  angle2 = geometry(ia,ib,ic,id2)       //calculate dihedral defined by 1182, 1094, 2206, 10 - this is the dihedral using the other C atom in the double bond of the chromene
         if (angle1 .gt. 90 .and. angle2 .lt. 0) then
            angle2 = angle2 + 360
            if (angle1 .gt. angle2) then
               write (iout,10) angle1,angle2
            else
               write (iout,20) angle1,angle2
            end if
         else if (angle2 .gt. 90 .and. angle1 .lt. 0) then
            angle1 = angle1 + 360
            if (angle1 .gt. angle2) then
               write (iout,10) angle1,angle2
            else
               write (iout,20) angle1,angle2
            end if
         else if (angle1 .gt. angle2) then
            write (iout,10) angle1,angle2
         else
            write (iout,20) angle1,angle2
         end if
   10    format (/,' Approaching Enantioface is Re',X,f7.2,X,f7.2)
   20    format (/,' Approaching Enantioface is Si',X,f7.2,X,f7.2)
*/
}

int CheckTypeOfPositionIRMOF(VECTOR pos)
{
  REAL HalfLargeCage = 5.04; //box length *1/4 - half window
  REAL HalfCage_Window = 7.876; //box length *1/4 + half window
  REAL Small_HalfLarge_Window = 17.956; //box length *3/4 - half window
  REAL Small_HalfLarge_2Windows = 20.792; //box length *3/4 + half window
  //REAL Cages_2Windows = 25.8320; //box length

if ((pos.x < HalfLargeCage && pos.y < HalfLargeCage && pos.z < HalfLargeCage)
|| (pos.x < HalfLargeCage && pos.y < HalfLargeCage && pos.z > Small_HalfLarge_2Windows)
|| (pos.x < HalfLargeCage && pos.y > Small_HalfLarge_2Windows && pos.z < HalfLargeCage)
|| (pos.x < HalfLargeCage && pos.y > Small_HalfLarge_2Windows && pos.z > Small_HalfLarge_2Windows)
|| (pos.x > Small_HalfLarge_2Windows && pos.y > Small_HalfLarge_2Windows && pos.z < HalfLargeCage)
|| (pos.x > Small_HalfLarge_2Windows && pos.y < HalfLargeCage && pos.z > Small_HalfLarge_2Windows)
|| (pos.x > Small_HalfLarge_2Windows && pos.y > Small_HalfLarge_2Windows && pos.z > Small_HalfLarge_2Windows)|| (pos.x > Small_HalfLarge_2Windows && pos.y < HalfLargeCage && pos.z < HalfLargeCage)
|| (pos.x < HalfLargeCage && HalfCage_Window < pos.y && pos.y < Small_HalfLarge_Window && HalfCage_Window < pos.z && pos.z < Small_HalfLarge_Window)
|| (pos.x > Small_HalfLarge_2Windows && HalfCage_Window < pos.y && pos.y < Small_HalfLarge_Window && HalfCage_Window < pos.z && pos.z < Small_HalfLarge_Window)|| (HalfCage_Window < pos.x && pos.x < Small_HalfLarge_Window && HalfCage_Window < pos.y && pos.y < Small_HalfLarge_Window && pos.z < HalfLargeCage)
|| (HalfCage_Window < pos.x && pos.x < Small_HalfLarge_Window && HalfCage_Window < pos.y && pos.y < Small_HalfLarge_Window && pos.z > Small_HalfLarge_2Windows)
|| (HalfCage_Window < pos.x && pos.x < Small_HalfLarge_Window && pos.y < HalfLargeCage && HalfCage_Window < pos.z && pos.z < Small_HalfLarge_Window)|| (HalfCage_Window < pos.x && pos.x < Small_HalfLarge_Window && pos.y > Small_HalfLarge_2Windows && HalfCage_Window < pos.z && pos.z < Small_HalfLarge_Window))
return BIG_CAGE;
if ((pos.x >= HalfLargeCage && pos.x <= HalfCage_Window)
|| (pos.x >= Small_HalfLarge_Window && pos.x <= Small_HalfLarge_2Windows)
|| (pos.y >= HalfLargeCage && pos.y <= HalfCage_Window)
|| (pos.y >= Small_HalfLarge_Window && pos.y <= Small_HalfLarge_2Windows)
|| (pos.z >= HalfLargeCage && pos.z <= HalfCage_Window)
|| (pos.z >= Small_HalfLarge_Window && pos.z <= Small_HalfLarge_2Windows))
return WINDOW;

return SMALL_CAGE;
}

void TrimStringInPlace(char *s)
{
    // Trim spaces and tabs from beginning:
    int i=0,j;
    while(isspace(s[i])) {
        i++;
    }
    if(i>0) {
        for(j=0;j<strlen(s);j++) {
            s[j]=s[j+i];
        }
    s[j]='\0';
    }

    // Trim spaces and tabs from end:
    i=strlen(s)-1;
    while(isspace(s[i])) {
        i--;
    }
    if(i<(strlen(s)-1)) {
        s[i+1]='\0';
    }
}

// Another set of functions using the && operator:
void rtrim( char * string, char * trim )
{
    int i;
    for( i = strlen (string) - 1; i >= 0
    && strchr ( trim, string[i] ) != NULL; i-- )
        // replace the string terminator:
        string[i] = '\0';
}

void ltrim( char * string, char * trim )
{
    while ( string[0] != '\0' && strchr ( trim, string[0] ) != NULL )
    {
        memmove( &string[0], &string[1], strlen(string) );
    }
}

void CompressSpacesInString(char *str)
{
        char *dst = str;

        for (; *str; ++str) {
                *dst++ = *str;
                if (isspace(*str)) {
                        do ++str; while (isspace(*str));
                        --str;
                }
        }
        *dst = 0;
}

void StripLeadingAndTrailingQuotesInPlace(char *str)
{
  char *p;

  if(str[0]=='\'')
  {
    memmove(&str[0],&str[1],strlen(str));

    p=strchr(str,'\'');
    if(p) *p='\0';
  }

  if(str[0]=='\"')
  {
    memmove(&str[0],&str[1],strlen(str));

    p=strchr(str,'\"');
    if(p) *p='\0';
  }

}



void StripWhiteSpacesInPlace(char *str)
{
    char *p;

    p = str;
    do
      if(!isspace(*p=*str)) p++;
    while (*str++);
}


void ReplaceCharacterInString(char * string,char searchchar,char replacechar)
{
  // Loop until end of string
  while (*string!='\0')
  {
    if(*string==searchchar)
      *string=replacechar;
    string++;
  }
}

void RemoveWhiteSpacesFromString(char *string)
{
  char *p2;

  p2=string;
  while(*string!='\0')
  {
    if(isspace(*string)) string++;
    else *p2++=*string++;
  }
  *p2='\0';
}

char* StringReverse(char *s)
{
  int i, j;
  char t[256];

  strcpy(t,s);
  for(i = 0 , j = strlen(s) - 1 ; j >= 0 ; i++, j--)
     *(s + i) = *(t + j);
  return s;
}

