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
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <getopt.h>

#include "simulation.h"
#include "run.h"

int main(int argc, char **argv)
{
  int c;
  char *input = NULL, *input_crystal = NULL, *raspa_dir = NULL, *output = NULL;
  bool stream = false;

  // set default for the inputs
  input = strdup("simulation.input");
  input_crystal = strdup("");

  // set default RASPA_DIR
  raspa_dir = getenv("HOME");
  strcat(raspa_dir,"/RASPA/simulations");

  // get the raspa install directory from environement if defined
  if(getenv("RASPA_DIR")&&(strlen(getenv("RASPA_DIR"))>0))
    raspa_dir=getenv("RASPA_DIR");
  // allow multiple versions of raspa to coexist peacefully by using different dirs
  if(getenv("RASPA2_DIR")&&(strlen(getenv("RASPA2_DIR"))>0))
    raspa_dir=getenv("RASPA2_DIR");

  // parse command-line options (":" means the option has an argument)
  while((c=getopt(argc,argv,"a:vhsc:i:d:"))!=-1)
  {
    switch(c)
    {
      case 'a':
        strcpy(FileNameAppend,optarg);
        break;
      case 'h':
        printf("usage: simulate [-hv] [-ifile] [-ddir] [-s [-idata] [-cdata]]\n");
        printf("\t-h help\n");
        printf("\t-v version\n");
        printf("\t-s enables streaming (inputs must be stream, not filepath)\n");
        printf("\t-i the input (either file or stream)\n");
        printf("\t-c if streaming, the crystal structure (as a stream)\n");
        printf("\t-d the raspa directory\n");
        printf("\t-a appends the string to output-files\n");
        return 0;
      case 'i': // set the input-filename
        input = strdup(optarg);
        break;
      case 'd':  // set the raspa-directory
        raspa_dir = strdup(optarg);
        break;
      case 's': // Toggle the streaming boolean
        stream = true;
        break;
      case 'c': // If streaming, pass the input molecule here
        input_crystal = strdup(optarg);
        break;
      case 'v':
        fprintf(stderr, "RASPA 2.0.41 (2021)\n");
        return 0;
      default:
        return 1;
        break;
    }
  }
  output = run(input, input_crystal, raspa_dir, stream);

  // This prints the output, which can be piped into other applications.
  if (stream)
    printf("%s\n", output);
  free(output);

  return 0;
}
