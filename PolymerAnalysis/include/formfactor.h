//  ROUTINE COMPUTES INTRAMOLECULAR MONOMER/MONOMER
//  CORRELATION FUNCTION
#ifndef FORMFACTOR_H
#define FORMFACTOR_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "nrutility.h"
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <vector>
using namespace std;
#define MXSZ 500

class FORMFACTOR
{
    public:
      FORMFACTOR();
      void calculate_wmmk(int,int,char*,char*);
      void calculate_wmmk_charge(int,int,vector<int>,char*, char*);
    private:
      void file_chck(FILE**, char*);
      struct srcv
      {
          int j;
          double x, y, z;
      };
};

#endif // FORMFACTOR_H
