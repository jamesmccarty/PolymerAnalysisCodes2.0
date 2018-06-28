//  ROUTINE COMPUTES MOLECULAR LENGTH STATISTICS
#include "leng_stat.h"

lengthstats::lengthstats()
{

}
void lengthstats::computelengthstats(int mono, int poly, char base_[],char cord_[], char fram_[])
{
    char comd_[MXSZ], pcrd_[MXSZ], pfrm_[MXSZ], prbd_[MXSZ], prgc_[MXSZ], pret_[MXSZ];
    double **xste, **yste, **zste, **xmon, **ymon, **zmon,
        *xcom, *ycom, *zcom, xdis, ydis, zdis,
        dist, grgc, grbd, gret;
    long a, i, j, k, nfrm, step;
    FILE *cord, *fram, *frbd, *frgc, *fret;

    //  INITIALIZE ALGORITHM PARAMETERS
    char dirn_[] = "LENGTH_STATS/";            //  OUTPUT: RESULTS DIRECTORY
    char frbd_[] = "monomerbondlength";             //  OUTPUT: MONOMER BOND LENGTH
    char frgc_[] = "radiusofgyration";             //  OUTPUT: MOLECULAR GYRATION RADIUS
    char fret_[] = "endtoenddistance";             //  OUTPUT: POLYMER EXTENSION
  //  mono = global::Nmono;                                        //  NUMBER OF SITES PER POLYMER
  //  poly = global::npolymer;                                       //  NUMBER OF CHAINS

    //  OPEN PROGRAM SEQUENCE
    fprintf(stdout, "Program sequence has started.\n");

    //  PREPARE SYSTEM REPOSITORIES
    fprintf(stdout, "System repositories are being prepared.\n");
    sprintf(pcrd_, "%s%s", base_, cord_);
    sprintf(pfrm_, "%s%s", base_, fram_);
    sprintf(prbd_, "%s%s", dirn_, frbd_);
    sprintf(prgc_, "%s%s", dirn_, frgc_);
    sprintf(pret_, "%s%s", dirn_, fret_);
    sprintf(comd_, "mkdir -p %s", dirn_);
    system(comd_);
    frbd = fopen(prbd_, "w");
    frgc = fopen(prgc_, "w");
    fret = fopen(pret_, "w");

    //  QUERY STATUS OF INPUT FILES
    fprintf(stdout, "Status of input files is being checked.\n");
    file_chck(&fram, pfrm_);
    file_chck(&cord, pcrd_);
    fprintf(stdout, "Status of input files is satisfactory.\n");

    //  COUNT SIMULATION FRAMES
    fprintf(stdout, "Frames are being counted.\n");
    fram = fopen(pfrm_, "r");
    nfrm = 0;
    while((fscanf(fram, "%li", &step))!=EOF)
    {
        nfrm++;
    }
    fclose(fram);

    //  ALLOCATE MEMORY RESOURCES
    fprintf(stdout, "Memory resources are being allocated.\n");
    xste = dmatrx(1,poly,1,mono);
    yste = dmatrx(1,poly,1,mono);
    zste = dmatrx(1,poly,1,mono);
    xmon = dmatrx(1,poly,1,mono);
    ymon = dmatrx(1,poly,1,mono);
    zmon = dmatrx(1,poly,1,mono);
    xcom = dvectr(1,poly);
    ycom = dvectr(1,poly);
    zcom = dvectr(1,poly);

    //  PERFORM CALCULATIONS
    fprintf(stdout, "Calculations are being performed.\n");
    grbd = 0.0;
    grgc = 0.0;
    gret = 0.0;
    fram = fopen(pfrm_, "r");
    cord = fopen(pcrd_, "r");
    for(a=1; a<=nfrm; a++)
    {
        //  READ FRAME CONFIGURATION
        fprintf(stdout, "FRAME %5li of %5li\n", a, nfrm);
        fscanf(fram, "%li", &step);
        for(i=1; i<=poly; i++)
        {
            for(j=1; j<=mono; j++)
            {
                fscanf(cord, "%li %lf %lf %lf", &k, &xdis, &ydis, &zdis);
                xste[i][j] = xdis;
                yste[i][j] = ydis;
                zste[i][j] = zdis;
            }
        }

        //  COMPUTE MOLECULAR CENTER-OF-MASS COORDINATES
        for(i=1; i<=poly; i++)
        {
            xdis = 0.0;
            ydis = 0.0;
            zdis = 0.0;
            for(j=1; j<=mono; j++)
            {
                xdis = xdis + xste[i][j];
                ydis = ydis + yste[i][j];
                zdis = zdis + zste[i][j];
            }
            xcom[i] = xdis/(double)mono;
            ycom[i] = ydis/(double)mono;
            zcom[i] = zdis/(double)mono;
        }

        //  COMPUTE MONOMER COORDINATES
        for(i=1; i<=poly; i++)
        {
            for(j=1; j<=mono; j++)
            {
                xdis = xste[i][j];
                ydis = yste[i][j];
                zdis = zste[i][j];
                xmon[i][j] = xdis;
                ymon[i][j] = ydis;
                zmon[i][j] = zdis;
            }
        }

        //  COMPUTE MONOMER BOND LENGTH
        dist = 0.0;
        for(i=1; i<=poly; i++)
        {
            for(j=1; j<mono; j++)
            {
                xdis = xmon[i][j+1] - xmon[i][j];
                ydis = ymon[i][j+1] - ymon[i][j];
                zdis = zmon[i][j+1] - zmon[i][j];
                dist = dist + xdis*xdis + ydis*ydis + zdis*zdis;
            }
        }
        dist = dist/(double)((mono-1)*poly);
        grbd = grbd + dist;
        fprintf(frbd, "% 10li\t% 20.16E\t% 20.16E\n",
            step, dist, grbd/(double)a);

        //  COMPUTE MOLECULAR RADIUS OF GYRATION
        dist = 0.0;
        for(i=1; i<=poly; i++)
        {
            for(j=1; j<=mono; j++)
            {
                xdis = xste[i][j] - xcom[i];
                ydis = yste[i][j] - ycom[i];
                zdis = zste[i][j] - zcom[i];
                dist = dist + xdis*xdis + ydis*ydis + zdis*zdis;
            }
        }
        dist = dist/(double)(poly*mono);
        grgc = grgc + dist;
        fprintf(frgc, "% 10li\t% 20.16E\t% 20.16E\n",
            step, dist, grgc/(double)a);

        //  COMPUTE MOLECULAR END-TO-END DISTANCE
        dist = 0.0;
        for(i=1; i<=poly; i++)
        {
            xdis = xmon[i][1] - xmon[i][mono];
            ydis = ymon[i][1] - ymon[i][mono];
            zdis = zmon[i][1] - zmon[i][mono];
            dist = dist + xdis*xdis + ydis*ydis + zdis*zdis;
        }
        dist = dist/(double)(poly);
        gret = gret + dist;
        fprintf(fret, "% 10li\t% 20.16E\t% 20.16E\n",
            step, dist, gret/(double)a);
    }
    fclose(fram);
    fclose(cord);
    fclose(frbd);
    fclose(frgc);
    fclose(fret);

    //  CLOSE PROGRAM SEQUENCE
    fprintf(stdout, "Program sequence has concluded.\n");
  //  return(0);
}

//  FILE EXISTENCE CHECKER
void lengthstats::file_chck(FILE **fptr, char fstr_[])
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
