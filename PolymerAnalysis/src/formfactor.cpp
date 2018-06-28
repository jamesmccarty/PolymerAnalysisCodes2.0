//  ROUTINE COMPUTES INTRAMOLECULAR CORRELATION FUNCTION
#include "formfactor.h"

FORMFACTOR::FORMFACTOR()
{

}
void FORMFACTOR::calculate_wmmk(int mono, int poly, char base_[], char fram_[])
{
    time_t si, so;
    struct srcv cont;
    char comd_[MXSZ], pfrm_[MXSZ], pcrd_[MXSZ], pwok_[MXSZ];
    double **xcor, **ycor, **zcor, *wkit,
        kmax, kval, delk, xdis, ydis, zdis, dist;
    long i, j, k, m, t, knum, nfrm, lsze;
    FILE *fram, *cord, *wofk;

    //  INITIALIZE PARAMETERS HERE:
    char cord_[] = "bmon.dat";               //  INPUT: BINARY SITE COORDINATES
    char dirn_[] = "STRUCTURE/";              //  OUTPUT: RESULTS DIRECTORY
    char wofk_[] = "wmmk.dat";               //  OUTPUT: INTRAMOLECULAR CORRELATION FUNCTION
    delk = 0.01;                                      //  INCREMENTAL FACTOR OF WAVE VECTOR
    kmax = 12.0;                                      //  MAXIMUM MAGNITUDE OF WAVE VECTOR


    //  OPEN PROGRAM SEQUENCE
    fprintf(stdout, "Program sequence has started.\n");

    //  PREPARE SYSTEM REPOSITORIES
    fprintf(stdout, "System repositories are being prepared.\n");
    sprintf(pfrm_, "%s%s", base_, fram_);
    sprintf(pcrd_, "%s%s", base_, cord_);
    sprintf(pwok_, "%s%s", dirn_, wofk_);
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
    knum = (long)(kmax/delk);
    xcor = dmatrx(1,poly,1,mono);
    ycor = dmatrx(1,poly,1,mono);
    zcor = dmatrx(1,poly,1,mono);
    wkit = dvectr(0,knum);
    for(k=0; k<=knum; k++)
    {
        wkit[k] = 0.0;
    }

    //  CALCULATE FORM FACTOR
    fprintf(stdout, "Form factor is being calculated.\n");
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
            kval = (double)(k*delk);
            for(m=1; m<=poly; m++)
            {
                for(i=1; i<mono; i++)
                {
                    for(j=i+1; j<=mono; j++)
                    {
                        xdis = xcor[m][i] - xcor[m][j];
                        ydis = ycor[m][i] - ycor[m][j];
                        zdis = zcor[m][i] - zcor[m][j];
                        dist = sqrt(xdis*xdis + ydis*ydis + zdis*zdis);
                        if((kval*dist)!=0.0)
                        {
                            wkit[k] = wkit[k] + 2.0*sin(kval*dist)/(kval*dist);
                        }
                        else
                        {
                            wkit[k] = wkit[k] + 2.0;
                        }
                    }
                }
            }
        }
        so = time(NULL);
        fprintf(stdout, "DELAY: %10li\n", so-si);
    }
    fclose(cord);

    //  NORMALIZE AND PRINT RESULTS FOR FORM FACTOR
    fprintf(stdout, "Results are being printed.\n");
    wofk = fopen(pwok_, "w");
    for(k=0; k<=knum; k++)
    {
        kval = (double)(k*delk);
        wkit[k] = wkit[k] + (double)(nfrm*poly*mono);
        fprintf(wofk, "% 20.16E\t% 20.16E\n", kval, wkit[k]/(double)(nfrm*poly*mono));
    }
    fclose(wofk);

    //  CLOSE PROGRAM SEQUENCE
    fprintf(stdout, "Program sequence has concluded.\n");
}

void FORMFACTOR::calculate_wmmk_charge(int mono,int poly, vector<int> charges,char base_[], char fram_[])
{
    time_t si, so;
    struct srcv cont;
    char comd_[MXSZ], pfrm_[MXSZ], pcrd_[MXSZ], pwok_[MXSZ];
    double **xcor, **ycor, **zcor, *wkit,
        kmax, kval, delk, xdis, ydis, zdis, dist;
    long i, j, k, m, t, knum, nfrm, lsze;
    double **seq;
    FILE *fram, *cord, *wofk;

    //  INITIALIZE PARAMETERS HERE:
    char cord_[] = "bmon.dat";               //  INPUT: BINARY SITE COORDINATES
    char dirn_[] = "STRUCTURE/";              //  OUTPUT: RESULTS DIRECTORY
    char wofk_[] = "wcck.dat";               //  OUTPUT: INTRAMOLECULAR CORRELATION FUNCTION
    delk = 0.01;                                      //  INCREMENTAL FACTOR OF WAVE VECTOR
    kmax = 12.0;                                      //  MAXIMUM MAGNITUDE OF WAVE VECTOR
  //  nchargetypes = global::natomtypes;
  //  double charge[nchargetypes] = global::monomercharges;

    //  OPEN PROGRAM SEQUENCE
    fprintf(stdout, "Program sequence has started.\n");

    //  PREPARE SYSTEM REPOSITORIES
    fprintf(stdout, "System repositories are being prepared.\n");
    sprintf(pfrm_, "%s%s", base_, fram_);
    sprintf(pcrd_, "%s%s", base_, cord_);
    sprintf(pwok_, "%s%s", dirn_, wofk_);
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
    knum = (long)(kmax/delk);
    xcor = dmatrx(1,poly,1,mono);
    ycor = dmatrx(1,poly,1,mono);
    zcor = dmatrx(1,poly,1,mono);
    seq = dmatrx(1,poly,1,mono);
    wkit = dvectr(0,knum);
    for(k=0; k<=knum; k++)
    {
        wkit[k] = 0.0;
    }

    //  CALCULATE FORM FACTOR
    fprintf(stdout, "Form factor is being calculated.\n");
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
                seq[m][i] = charges[(cont.j-1)];
            }
        }
        for(k=0; k<=knum; k++)
        {
            kval = (double)(k*delk);
            for(m=1; m<=poly; m++)
            {
                for(i=1; i<mono; i++)
                {
                    for(j=i+1; j<=mono; j++)
                    {
                        xdis = xcor[m][i] - xcor[m][j];
                        ydis = ycor[m][i] - ycor[m][j];
                        zdis = zcor[m][i] - zcor[m][j];
                        dist = sqrt(xdis*xdis + ydis*ydis + zdis*zdis);
                        if((kval*dist)!=0.0)
                        {
                            wkit[k] = wkit[k] + 2.0*seq[m][i]*seq[m][j]*sin(kval*dist)/(kval*dist);
                        }
                        else
                        {
                            wkit[k] = wkit[k] + 0.0;
                        }
                    }
                }
            }
        }
        so = time(NULL);
        fprintf(stdout, "DELAY: %10li\n", so-si);
    }
    fclose(cord);

    //  NORMALIZE AND PRINT RESULTS FOR FORM FACTOR
    fprintf(stdout, "Results are being printed.\n");
    wofk = fopen(pwok_, "w");
    for(k=0; k<=knum; k++)
    {
        kval = (double)(k*delk);
        wkit[k] = wkit[k] + (double)(nfrm*poly*mono);
        fprintf(wofk, "% 20.16E\t% 20.16E\n", kval, wkit[k]/(double)(nfrm*poly*mono));
    }
    fclose(wofk);

    //  CLOSE PROGRAM SEQUENCE
    fprintf(stdout, "Program sequence has concluded.\n");
}

//  FILE EXISTENCE CHECKER
void FORMFACTOR::file_chck(FILE **fptr, char fstr_[])
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
