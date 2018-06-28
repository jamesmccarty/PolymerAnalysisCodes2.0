#ifndef BINAFILE_H
#define BINAFILE_H
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <cstdlib>
using namespace std;
#include "nrutility.h"
#define MXSZ 500

class binafile
{
    public:
      binafile();
      void convertbinary(int,int,char*,char*,char*);
    private:
      void file_chck(FILE**, char*);
      struct srcv
      {
          int j;
          double x, y, z;
      };
};

#endif // BINAFILE_H
