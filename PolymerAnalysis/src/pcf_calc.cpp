//  ROUTINE COMPUTES INTERMOLECULAR MONOMER/MONOMER
//  TOTAL CORRELATION FUNCTION
#include "pcf_calc.h"

PAIRCORRELATION::PAIRCORRELATION()
{

}
void PAIRCORRELATION::calculate_hmmr(int mono, int poly,double boxs,char base_[], char fram_[])
{
    time_t si, so;
    struct srcv cont;
    char comd_[MXSZ], pfrm_[MXSZ], pcrd_[MXSZ], phor_[MXSZ];
    double **xcor, **ycor, **zcor, *hist, delr, volu,
        xdis, ydis, zdis, dist, invd, invb, rmax, shex, shin, scal;
    long i, j, m, n, t, nbin, nfrm, lsze, ibin;
    FILE *cord, *fram, *hofr;

    //  INITIALIZE ALGORITHM PARAMETERS
    char cord_[] = "bmon.dat";               //  INPUT: BINARY SITE COORDINATES
    char dirn_[] = "STRUCTURE/";              //  OUTPUT: RESULTS DIRECTORY
    char hofr_[] = "hmmr.dat";               //  OUTPUT: TOTAL COORELATION FUNCTION
//    mono = global::Nmono;                                        //  NUMBER OF SITES PER POLYMER
//    poly = global::npolymer;                                       //  NUMBER OF CHAINS
    delr = 0.05;                                       //  HISTOGRAM BIN WIDTH
    scal = 0.5;                                       //  BOX DIMENSION SCALING FACTOR

    //  OPEN PROGRAM SEQUENCE
    fprintf(stdout, "Program sequence has started.\n");

    //  PREPARE SYSTEM REPOSITORIES
    fprintf(stdout, "System repositories are being prepared.\n");
    sprintf(pfrm_, "%s%s", base_, fram_);
    sprintf(pcrd_, "%s%s", base_, cord_);
    sprintf(phor_, "%s%s", dirn_, hofr_);
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
    nbin = (long)(boxs*scal/delr);
    xcor = dmatrx(1, poly, 1, mono);
    ycor = dmatrx(1, poly, 1, mono);
    zcor = dmatrx(1, poly, 1, mono);
    hist = dvectr(0, nbin-1);
    for(i=0; i<nbin; i++)
    {
        hist[i] = 0.0;
    }

    //  COMPUTE HISTOGRAM FOR CORRELATION FUNCTION
    fprintf(stdout, "Histogram is being calculated.\n");
    invd = 1.0/delr;
    invb = 1.0/boxs;
    rmax = (double)(nbin*delr);
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
        for(m=1; m<poly; m++)
        {
            for(n=m+1; n<=poly; n++)
            {
                for(i=1; i<=mono; i++)
                {
                    for(j=1; j<=mono; j++)
                    {
                        xdis = xcor[m][i] - xcor[n][j];
                        ydis = ycor[m][i] - ycor[n][j];
                        zdis = zcor[m][i] - zcor[n][j];
                        xdis = xdis - boxs*(double)nrst_long(xdis*invb);
                        ydis = ydis - boxs*(double)nrst_long(ydis*invb);
                        zdis = zdis - boxs*(double)nrst_long(zdis*invb);
                        dist = sqrt(xdis*xdis + ydis*ydis + zdis*zdis);
                        if(dist<rmax)
                        {
                            ibin = (long)(dist*invd);
                            hist[ibin] = hist[ibin] + 2.0;
                        }
                    }
                }
            }
        }
        so = time(NULL);
        fprintf(stdout, "DELAY: %10li\n", so-si);
    }
    fclose(cord);

    //  NORMALIZE AND PRINT RESULTS FOR CORRELATION FUNCTION
    fprintf(stdout, "Results are being printed.\n");
    hofr = fopen(phor_, "w");
    for(i=0; i<nbin; i++)
    {
        shex = (i+1)*(i+1)*(i+1);
        shin = i*i*i;
        volu = 4.0/3.0*pi*(shex-shin)*delr*delr*delr;
        hist[i] = hist[i]
                    / (double)(mono*poly*mono*poly*nfrm*volu*invb*invb*invb);
        fprintf(hofr, "% 20.16E\t% 20.16E\n", (double)((i+0.5)*delr), hist[i]-1.0);
    }
    fclose(hofr);

    //  CLOSE PROGRAM SEQUENCE
    fprintf(stdout, "Program sequence has concluded.\n");
}

void PAIRCORRELATION::calculate_hcomr(int mono, int poly, double boxs,char base_[], char fram_[])
{
    time_t si, so;
    struct srcv cont;
    char comd_[MXSZ], pfrm_[MXSZ], pcrd_[MXSZ], phor_[MXSZ];
    double *xcor, *ycor, *zcor, *hist, delr, volu,
        xdis, ydis, zdis, dist, invd, invb, rmax, shex, shin, scal;
    long i, m, n, t, nbin, nfrm, lsze, ibin;
    FILE *cord, *fram, *hofr;

    //  INITIALIZE ALGORITHM PARAMETERS
    char cord_[] = "bcom.dat";               //  INPUT: BINARY SITE COORDINATES
    char dirn_[] = "STRUCTURE/";              //  OUTPUT: RESULTS DIRECTORY
    char hofr_[] = "hcomr.dat";               //  OUTPUT: TOTAL COORELATION FUNCTION
    delr = 0.1*sqrt(mono);                            //  HISTOGRAM BIN WIDTH
    scal = 0.5;                                       //  BOX DIMENSION SCALING FACTOR

    //  OPEN PROGRAM SEQUENCE
    fprintf(stdout, "Program sequence has started.\n");

    //  PREPARE SYSTEM REPOSITORIES
    fprintf(stdout, "System repositories are being prepared.\n");
    sprintf(pfrm_, "%s%s", base_, fram_);
    sprintf(pcrd_, "%s%s", base_, cord_);
    sprintf(phor_, "%s%s", dirn_, hofr_);
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
    nbin = (long)(boxs*scal/delr);
    xcor = dvectr(1, poly);
    ycor = dvectr(1, poly);
    zcor = dvectr(1, poly);
    hist = dvectr(0, nbin-1);
    for(i=0; i<nbin; i++)
    {
        hist[i] = 0.0;
    }

    //  COMPUTE HISTOGRAM FOR CORRELATION FUNCTION
    fprintf(stdout, "Histogram is being calculated.\n");
    invd = 1.0/delr;
    invb = 1.0/boxs;
    rmax = (double)(nbin*delr);
    lsze = sizeof(struct srcv);
    cord = fopen(pcrd_, "r");
    for(t=1; t<=nfrm; t++)
    {
        si = time(NULL);
        fprintf(stdout, "CYCLE %5li OF %5li\t", t, nfrm);
        for(m=1; m<=poly; m++)
        {
            fread(&cont, lsze, 1, cord);
            xcor[m] = cont.x;
            ycor[m] = cont.y;
            zcor[m] = cont.z;
        }
        for(m=1; m<poly; m++)
        {
            for(n=m+1; n<=poly; n++)
            {
                xdis = xcor[m] - xcor[n];
                ydis = ycor[m] - ycor[n];
                zdis = zcor[m] - zcor[n];
                xdis = xdis - boxs*(double)nrst_long(xdis*invb);
                ydis = ydis - boxs*(double)nrst_long(ydis*invb);
                zdis = zdis - boxs*(double)nrst_long(zdis*invb);
                dist = sqrt(xdis*xdis + ydis*ydis + zdis*zdis);
                if(dist<rmax)
                {
                    ibin = (long)(dist*invd);
                    hist[ibin] = hist[ibin] + 2.0;
                }
            }
        }
        so = time(NULL);
        fprintf(stdout, "DELAY: %10li\n", so-si);
    }
    fclose(cord);

    //  NORMALIZE AND PRINT RESULTS FOR CORRELATION FUNCTION
    fprintf(stdout, "Results are being printed.\n");
    hofr = fopen(phor_, "w");
    for(i=0; i<nbin; i++)
    {
        shex = (i+1)*(i+1)*(i+1);
        shin = i*i*i;
        volu = 4.0/3.0*pi*(shex-shin)*delr*delr*delr;
        hist[i] = hist[i]
                    / (double)(poly*poly*nfrm*volu*invb*invb*invb);
        fprintf(hofr, "% 20.16E\t% 20.16E\n", (double)((i+0.5)*delr), hist[i]-1.0);
    }
    fclose(hofr);

    //  CLOSE PROGRAM SEQUENCE
    fprintf(stdout, "Program sequence has concluded.\n");
}


//  FILE EXISTENCE CHECKER
void PAIRCORRELATION::file_chck(FILE **fptr, char fstr_[])
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

//  NEAREST INTEGER HANDLER
long PAIRCORRELATION::nrst_long(double argm)
{
    if(argm>0)
    {
        return (long)(argm + 0.5);
    }
    else
    {
        return (long)(argm - 0.5);
    }
}
