//  ROUTINE COMPUTES STRUCTURE FACTORS
#include "Sofk.h"

STRUCTUREFACTOR::STRUCTUREFACTOR()
{

}
void STRUCTUREFACTOR::calculate_Smmk(int mono, int poly, double L, char base_[], char fram_[])
{
    time_t si, so;
    struct srcv cont;
    char comd_[MXSZ], pfrm_[MXSZ], pcrd_[MXSZ], psok_[MXSZ];
    double **xcor, **ycor, **zcor, *Sk, *Sk_cos, *Sk_sin,
        delk;
    long i, k, m, t, knum, nfrm, lsze;
    double dkx,dky,dkz;
    double Lx,Ly,Lz;
    int nx,ny,nz;
    int nkbound;
    FILE *fram, *cord, *Sofk;

    typedef struct
    {
         double x;
         double y;
         double z;
    } VECTOR;
    VECTOR Wavevector[100000];
    VECTOR k_scale;

    //  INITIALIZE PARAMETERS HERE:
    char cord_[] = "bmon.dat";               //  INPUT: BINARY SITE COORDINATES
    char dirn_[] = "STRUCTURE/";              //  OUTPUT: RESULTS DIRECTORY
    char Sofk_[] = "Smmk.dat";               //  OUTPUT: INTRAMOLECULAR CORRELATION FUNCTION
//    mono = global::Nmono;                                        //  NUMBER OF SITES PER POLYMER
//    poly = global::npolymer;                                       //  NUMBER OF CHAINS
    delk = 1.0;                                      //  INCREMENTAL FACTOR OF WAVE VECTOR
    nkbound = 23;
    Lx=L;Ly=L;Lz=L;
    dkx = 2.0*M_PI/Lx*delk;
    dky = 2.0*M_PI/Ly*delk;
    dkz = 2.0*M_PI/Lz*delk;
    //  OPEN PROGRAM SEQUENCE
    fprintf(stdout, "Program sequence has started.\n");

    //  PREPARE SYSTEM REPOSITORIES
    fprintf(stdout, "System repositories are being prepared.\n");
    sprintf(pfrm_, "%s%s", base_, fram_);
    sprintf(pcrd_, "%s%s", base_, cord_);
    sprintf(psok_, "%s%s", dirn_, Sofk_);
    sprintf(comd_, "mkdir -p %s", dirn_);
    system(comd_);

    //  QUERY STATUS OF INPUT FILES
    fprintf(stdout, "Status of input files is being checked.\n");
    file_chck(&fram, pfrm_);
    file_chck(&cord, pcrd_);
    fprintf(stdout, "Status of input files is satisfactory.\n");

    //  COUNT SIMULATION FRAMES
    fprintf(stdout, "Frames are being counted.\n");
    fram = fopen(pfrm_, "r");
    nfrm = 0;
    while((fscanf(fram, "%li", &i))!=EOF)
    {
        nfrm++;
    }
    fclose(fram);

    //  ALLOCATE MEMORY RESOURCES
    fprintf(stdout, "Memory resources are being allocated.\n");

    knum = 0.0;
    for(nx=0;nx<=nkbound;nx++)
    for(ny=0;ny<=nkbound;ny++)
    for(nz=0;nz<=nkbound;nz++)
    {
        Wavevector[knum].x = nx*dkx;
        Wavevector[knum].y = ny*dky;
	Wavevector[knum].z = nz*dkz;
        knum++;
    }
    printf("number of k's = %ld\n",knum);

    /******************** bubble sort k's ***********************/
    int itop,i1;
    double rxi,ryi,rzi;
    int change;
    itop = knum-1;
    do
    {
    change = false;
    for(i=0;i<itop;i++)
    {
     i1 = i+1;
     if(Wavevector[i].x*Wavevector[i].x + Wavevector[i].y*Wavevector[i].y + Wavevector[i].z*Wavevector[i].z > Wavevector[i1].x*Wavevector[i1].x+Wavevector[i1].y*Wavevector[i1].y+Wavevector[i1].z*Wavevector[i1].z)
     {
       rxi =  Wavevector[i].x;
       ryi =  Wavevector[i].y;
       rzi =  Wavevector[i].z;

       Wavevector[i].x = Wavevector[i1].x;
       Wavevector[i].y = Wavevector[i1].y;
       Wavevector[i].z = Wavevector[i1].z;

       Wavevector[i1].x = rxi;
       Wavevector[i1].y = ryi;
       Wavevector[i1].z = rzi;
       change = true;
     }
    }//end loop i
    itop--;
    }
    while(change==true && itop > 0);

    double kmax;
    int nkmax;
    nkmax = knum-1;
    kmax = sqrt(Wavevector[nkmax].x*Wavevector[nkmax].x+Wavevector[nkmax].y*Wavevector[nkmax].y+Wavevector[nkmax].z*Wavevector[nkmax].z);
    printf("maximum k = %lf\n",kmax);

  /*******************************************/


    xcor = dmatrx(1,poly,1,mono);
    ycor = dmatrx(1,poly,1,mono);
    zcor = dmatrx(1,poly,1,mono);
    Sk = dvectr(0,knum);
    Sk_cos = dvectr(0,knum);
    Sk_sin = dvectr(0,knum);
    for(k=0; k<=knum; k++)
    {
        Sk[k] = 0.0;
        Sk_cos[k] = 0.0;
        Sk_sin[k] = 0.0;
    }

    //  CALCULATE FORM FACTOR
    fprintf(stdout, "Structure factor is being calculated.\n");
    lsze = sizeof(struct srcv);
    cord = fopen(pcrd_, "r");
    for(t=1; t<=nfrm; t++)
    {
        si = time(NULL);
        fprintf(stdout, "CYCLE %5li OF %5li\t", t, nfrm);
        for(m=1; m<=poly; m++)
        {
            for(i=1; i<=mono; i++)
            {
                fread(&cont, lsze, 1, cord);
                xcor[m][i] = cont.x;
                ycor[m][i] = cont.y;
                zcor[m][i] = cont.z;
            }
        }
        for(k=0; k<=knum; k++)
        {
            k_scale.x = Wavevector[k].x*L/L;k_scale.y = Wavevector[k].y*L/L; k_scale.z = Wavevector[k].z*L/L;
            Sk_cos[k] = 0.0;Sk_sin[k] = 0.0;
            for(m=1; m<=poly; m++)
            {
                for(i=1; i<=mono; i++)
                {

                      Sk_cos[k] += cos(k_scale.x*xcor[m][i] + k_scale.y*ycor[m][i] + k_scale.z*zcor[m][i]);
                      Sk_sin[k] += sin(k_scale.x*xcor[m][i] + k_scale.y*ycor[m][i] + k_scale.z*zcor[m][i]);
                }
            }   // end particle loop
        Sk_cos[k] = Sk_cos[k]*Sk_cos[k];
        Sk_sin[k] = Sk_sin[k]*Sk_sin[k];
        Sk[k] += Sk_cos[k] + Sk_sin[k];
        } // end k loop
        so = time(NULL);
        fprintf(stdout, "DELAY: %10li\n", so-si);
    }
    fclose(cord);

    //  NORMALIZE AND PRINT RESULTS FOR FORM FACTOR
    fprintf(stdout, "Results are being printed.\n");
    for(k=0; k<=knum; k++)
    {
        Sk[k] /= (double)(nfrm*poly*mono);
    }
    int kcount;
    double knorm,skmean;
    double kminus,kplus;
    kcount=0; knorm=0.0; skmean=0.0;
    Sofk = fopen(psok_, "w");
    for(k=0; k<knum; k++)
    {
        kcount++;
        knorm += sqrt(Wavevector[k].x*Wavevector[k].x+Wavevector[k].y*Wavevector[k].y+Wavevector[k].z*Wavevector[k].z);
        skmean += Sk[k];
        kplus = (Wavevector[k+1].x*Wavevector[k+1].x+Wavevector[k+1].y*Wavevector[k+1].y+Wavevector[k+1].z*Wavevector[k+1].z);
        kminus = (Wavevector[k].x*Wavevector[k].x+Wavevector[k].y*Wavevector[k].y+Wavevector[k].z*Wavevector[k].z);
        if( fabs(kplus-kminus) > 0.001 || k+1 == knum )
        {
          knorm /= kcount;  skmean /= kcount;
          if(knorm > 2.0*M_PI/L*delk)
          {
          fprintf(Sofk,"% 20.16E\t% 20.16E\n",knorm,skmean);
          }
          kcount=0; knorm=0.0; skmean=0.0;
        }
    }
    fclose(Sofk);

    //  CLOSE PROGRAM SEQUENCE
    fprintf(stdout, "Program sequence has concluded.\n");
}

void STRUCTUREFACTOR::calculate_Scck(int mono, int poly, double L, vector<int> charges, char base_[], char fram_[])
{
    time_t si, so;
    struct srcv cont;
    char comd_[MXSZ], pfrm_[MXSZ], pcrd_[MXSZ], psok_[MXSZ];
    double **xcor, **ycor, **zcor, **type, *Sk, *SkA_cos, *SkA_sin, *SkB_cos, *SkB_sin,
        delk;
    long i, k, m, t, knum, nfrm, lsze;
    double dkx,dky,dkz;
    double Lx,Ly,Lz;
    int nx,ny,nz;
    int nkbound;
    FILE *fram, *cord, *Sofk;

    typedef struct
    {
         double x;
         double y;
         double z;
    } VECTOR;
    VECTOR Wavevector[100000];
    VECTOR k_scale;

      //  INITIALIZE PARAMETERS HERE:
    char cord_[] = "bmon.dat";               //  INPUT: BINARY SITE COORDINATES
    char dirn_[] = "STRUCTURE/";              //  OUTPUT: RESULTS DIRECTORY
    char Sofk_[] = "Scck.dat";               //  OUTPUT: INTRAMOLECULAR CORRELATION FUNCTION
    delk = 1.0;                                      //  INCREMENTAL FACTOR OF WAVE VECTOR
    nkbound = 23;
    int Ntype=0;
    Lx=L;Ly=L;Lz=L;
    dkx = 2.0*M_PI/Lx*delk;
    dky = 2.0*M_PI/Ly*delk;
    dkz = 2.0*M_PI/Lz*delk;
    //  OPEN PROGRAM SEQUENCE
    fprintf(stdout, "Program sequence has started.\n");
    //  PREPARE SYSTEM REPOSITORIES
    fprintf(stdout, "System repositories are being prepared.\n");
    sprintf(pfrm_, "%s%s", base_, fram_);
    sprintf(pcrd_, "%s%s", base_, cord_);
    sprintf(psok_, "%s%s", dirn_, Sofk_);
    sprintf(comd_, "mkdir -p %s", dirn_);
    system(comd_);

    //  QUERY STATUS OF INPUT FILES
    fprintf(stdout, "Status of input files is being checked.\n");
    file_chck(&fram, pfrm_);
    file_chck(&cord, pcrd_);
    fprintf(stdout, "Status of input files is satisfactory.\n");

    //  COUNT SIMULATION FRAMES
    fprintf(stdout, "Frames are being counted.\n");
    fram = fopen(pfrm_, "r");
    nfrm = 0;
    while((fscanf(fram, "%li", &i))!=EOF)
    {
        nfrm++;
    }
    fclose(fram);

    //  ALLOCATE MEMORY RESOURCES
    fprintf(stdout, "Memory resources are being allocated.\n");

    knum = 0.0;
    for(nx=0;nx<=nkbound;nx++)
    for(ny=0;ny<=nkbound;ny++)
    for(nz=0;nz<=nkbound;nz++)
    {
        Wavevector[knum].x = nx*dkx;
        Wavevector[knum].y = ny*dky;
	Wavevector[knum].z = nz*dkz;
        knum++;
    }
    printf("number of k's = %ld\n",knum);

    /******************** bubble sort k's ***********************/
    int itop,i1;
    double rxi,ryi,rzi;
    int change;
    itop = knum-1;
    do
    {
    change = false;
    for(i=0;i<itop;i++)
    {
     i1 = i+1;
     if(Wavevector[i].x*Wavevector[i].x + Wavevector[i].y*Wavevector[i].y + Wavevector[i].z*Wavevector[i].z > Wavevector[i1].x*Wavevector[i1].x+Wavevector[i1].y*Wavevector[i1].y+Wavevector[i1].z*Wavevector[i1].z)
     {
       rxi =  Wavevector[i].x;
       ryi =  Wavevector[i].y;
       rzi =  Wavevector[i].z;

       Wavevector[i].x = Wavevector[i1].x;
       Wavevector[i].y = Wavevector[i1].y;
       Wavevector[i].z = Wavevector[i1].z;

       Wavevector[i1].x = rxi;
       Wavevector[i1].y = ryi;
       Wavevector[i1].z = rzi;
       change = true;
     }
    }//end loop i
    itop--;
    }
    while(change==true && itop > 0);

    double kmax;
    int nkmax;
    nkmax = knum-1;
    kmax = sqrt(Wavevector[nkmax].x*Wavevector[nkmax].x+Wavevector[nkmax].y*Wavevector[nkmax].y+Wavevector[nkmax].z*Wavevector[nkmax].z);
    printf("maximum k = %lf\n",kmax);

  /*******************************************/


    xcor = dmatrx(1,poly,1,mono);
    ycor = dmatrx(1,poly,1,mono);
    zcor = dmatrx(1,poly,1,mono);
    type = dmatrx(1,poly,1,mono);
    Sk = dvectr(0,knum);
    SkA_cos = dvectr(0,knum);
    SkA_sin = dvectr(0,knum);
    SkB_cos = dvectr(0,knum);
    SkB_sin = dvectr(0,knum);
    for(k=0; k<=knum; k++)
    {
        Sk[k] = 0.0;
        SkA_cos[k] = 0.0;
        SkA_sin[k] = 0.0;
        SkB_cos[k] = 0.0;
        SkB_sin[k] = 0.0;
    }

    //  CALCULATE FORM FACTOR
    fprintf(stdout, "Structure factor is being calculated.\n");
    lsze = sizeof(struct srcv);
    cord = fopen(pcrd_, "r");
    for(t=1; t<=nfrm; t++)
    {
        si = time(NULL);
        fprintf(stdout, "CYCLE %5li OF %5li\t", t, nfrm);
        for(m=1; m<=poly; m++)
        {
            for(i=1; i<=mono; i++)
            {
              fread(&cont, lsze, 1, cord);
                xcor[m][i] = cont.x;
                ycor[m][i] = cont.y;
                zcor[m][i] = cont.z;
                type[m][i]= charges[(cont.j-1)];
            }
        }
        for(k=0; k<=knum; k++)
        {
            k_scale.x = Wavevector[k].x*L/L;k_scale.y = Wavevector[k].y*L/L; k_scale.z = Wavevector[k].z*L/L;
            SkA_cos[k] = 0.0;SkA_sin[k] = 0.0;
            SkB_cos[k] = 0.0;SkB_sin[k] = 0.0;
            Ntype=0;
            for(m=1; m<=poly; m++)
            {
                for(i=1; i<=mono; i++)
                {
                      if(type[m][i]==-1)
                      {
                      SkA_cos[k] += cos(k_scale.x*xcor[m][i] + k_scale.y*ycor[m][i] + k_scale.z*zcor[m][i]);
                      SkA_sin[k] += sin(k_scale.x*xcor[m][i] + k_scale.y*ycor[m][i] + k_scale.z*zcor[m][i]);
                      }
                      else if(type[m][i]==1)
                      {
                      SkB_cos[k] += cos(k_scale.x*xcor[m][i] + k_scale.y*ycor[m][i] + k_scale.z*zcor[m][i]);
                      SkB_sin[k] += sin(k_scale.x*xcor[m][i] + k_scale.y*ycor[m][i] + k_scale.z*zcor[m][i]);
                      }
                      else
                      {
                          fprintf(stdout, "Unrecognized monomer charge\n");
                          exit(1);
                      }
                      Ntype++;
                }
            }   // end particle loop
        Sk[k] += SkA_cos[k]*SkA_cos[k] + SkA_sin[k]*SkA_sin[k] + SkB_cos[k]*SkB_cos[k] + SkB_sin[k]*SkB_sin[k] - 2.0*(SkA_cos[k]*SkB_cos[k] + SkA_sin[k]*SkB_sin[k]);
        } // end k loop
        so = time(NULL);
        fprintf(stdout, "DELAY: %10li\n", so-si);
    }
    fclose(cord);

    //  NORMALIZE AND PRINT RESULTS FOR FORM FACTOR
    fprintf(stdout, "Results are being printed.\n");
    for(k=0; k<=knum; k++)
    {
        Sk[k] /= (double)(nfrm*poly*mono);
    }
    int kcount;
    double knorm,skmean;
    double kminus,kplus;
    kcount=0; knorm=0.0; skmean=0.0;
    Sofk = fopen(psok_, "w");
    for(k=0; k<knum; k++)
    {
        kcount++;
        knorm += sqrt(Wavevector[k].x*Wavevector[k].x+Wavevector[k].y*Wavevector[k].y+Wavevector[k].z*Wavevector[k].z);
        skmean += Sk[k];
        kplus = (Wavevector[k+1].x*Wavevector[k+1].x+Wavevector[k+1].y*Wavevector[k+1].y+Wavevector[k+1].z*Wavevector[k+1].z);
        kminus = (Wavevector[k].x*Wavevector[k].x+Wavevector[k].y*Wavevector[k].y+Wavevector[k].z*Wavevector[k].z);
        if( fabs(kplus-kminus) > 0.001 || k+1 == knum )
        {
          knorm /= kcount;  skmean /= kcount;
          fprintf(Sofk,"% 20.16E\t% 20.16E\n",knorm,skmean);
          kcount=0; knorm=0.0; skmean=0.0;
        }
    }
    fclose(Sofk);

    //  CLOSE PROGRAM SEQUENCE
    fprintf(stdout, "Program sequence has concluded.\n");
}

//  FILE EXISTENCE CHECKER
void STRUCTUREFACTOR::file_chck(FILE **fptr, char fstr_[])
{
    if(!(*fptr = fopen(fstr_, "r")))
    {
        fprintf(stdout, "Failed to open file %s\n", fstr_);
        exit(1);
    }
    else
    {
        fclose(*fptr);
    }
}
