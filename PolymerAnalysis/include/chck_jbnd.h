//  ROUTINE COMPUTES REJOINED SITE COORDINATES USING PERIODIC BOUNDARY CONDITIONS
#ifndef CHCK_JBOND_H
#define CHCK_JBOND_H
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <cmath>
using namespace std;
#include "nrutility.h"
#define MXTH 2.0
#define MNTH 1.0
#define MXSZ 500

class checkbond
{
    public:
      checkbond();
      void fixbrokenbonds(int,int,double,char*,char*,char*);
    private:
      void file_chck(FILE **, char *);
};

#endif // CHCK_JBOND_H
