//  ROUTINE COMPUTES INTERMOLECULAR MONOMER/MONOMER
//  TOTAL CORRELATION FUNCTION
#ifndef PCF_CALC_H
#define PCF_CALC_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "nrutility.h"
#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>
using namespace std;
#define pi M_PI
#define MXSZ 500

class PAIRCORRELATION
{
    public:
      PAIRCORRELATION();
      void calculate_hmmr(int,int,double,char*,char*);
      void calculate_hcomr(int,int,double,char*, char*);
    private:
      void file_chck(FILE**, char*);
      long nrst_long(double);
      struct srcv
      {
          int j;
          double x, y, z;
      };
};

#endif // PCF_CALC_H
