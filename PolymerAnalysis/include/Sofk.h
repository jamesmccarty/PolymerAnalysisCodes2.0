#ifndef SOFK_H
#define SOFK_H
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
using namespace std;
#define pi M_PI
#define MXSZ 500

class STRUCTUREFACTOR
{
    public:
      STRUCTUREFACTOR();
      void calculate_Smmk(int,int,double,char*,char*);
      void calculate_Scck(int,int,double,vector<int>,char*, char*);
    private:
      void file_chck(FILE**, char*);
      long nrst_long(double);
      struct srcv
      {
          int j;
          double x, y, z;
      };
};

#endif // SOFK_H
