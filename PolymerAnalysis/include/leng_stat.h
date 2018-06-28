//  ROUTINE COMPUTES MOLECULAR LENGTH STATISTICS
#ifndef LENG_STAT_H
#define LENG_STAT_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nrutility.h"
#define MXSZ 500
using namespace std;

class lengthstats
{
    public:
      lengthstats();
      void computelengthstats(int,int,char*,char*,char*);
    private:
      void file_chck(FILE **, char *);
};

#endif // LENG_STAT_H
