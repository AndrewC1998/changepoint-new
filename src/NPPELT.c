#include <R.h>
#include <Rmath.h>
#include <Rinternals.h> // RK addition
#include <R_ext/RS.h>  // RK addition
#include <R_ext/Lapack.h> // RK addition
#include <R_ext/BLAS.h> // RK addition
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
//#include "costfunctions.c"
#define SWAP(a,b)   { int t; t=a; a=b; b=t; }  // Macro for swapping

static int *checklist;
static double *tmplike;

void FreeNPPELT(error)
  int *error; /* Error code from PELT C function, non-zero => error */
  {
    if(*error==0){
      free((void *)checklist);
      free((void *)tmplike);
    }
  }

  void NPPELT(cost_func, sumstat, n, pen, cptsout, error, minseglen, nquantiles, lastchangelike, lastchangecpts, numchangecpts, MBIC)
  char **cost_func;
  double *sumstat;    //Summary statistic for the time series
  int *n;			// Length of the time series
  double *pen;  // Penalty used to decide if a changepoint is significant
  int *cptsout;    // Vector of identified changepoint locations
  int *error;      //0 by default, nonzero indicates error in code
  int *minseglen; //minimum segment length
  int *nquantiles; //number of quantiles in the empirical distribution approximation (K in EDPELT)
  double *lastchangelike; // stores likelihood up to that time using optimal changepoint locations up to that time
  int *lastchangecpts; // stores last changepoint locations
  int *numchangecpts; //stores the current number of changepoints
  int *MBIC; //// 1 if MBIC penalty, 0 if not
  {

    void (*costfunction)();
    void nonparametric_ed();

    if (strcmp(*cost_func,"nonparametric.ed")==0){
      costfunction = &nonparametric_ed;
    }
    else if (strcmp(*cost_func,"nonparametric.ed.mbic")==0){
      costfunction = &nonparametric_ed;
    }

    int *checklist;
    checklist = (int *)calloc(*n+1,sizeof(int));
    if (checklist==NULL)   {
      *error = 1;
      goto err1;
    }

    int nchecklist;
    double minout;
    double segcost=0;   //cost over specified segment

    double *tmplike;
    tmplike = (double *)calloc(*n+1,sizeof(double));
    if (tmplike==NULL)   {
      *error = 2;
      goto err2;
    }

    int tstar,i,whichout,nchecktmp;


    void min_which();

    lastchangelike[0]= -*pen;
    lastchangecpts[0]=0;
    numchangecpts[0]=0;
    int j;
    int isum;
    double *sumstatout;
    sumstatout = (double *)calloc(*nquantiles,sizeof(double));
     for(j=*minseglen;j<(2*(*minseglen));j++){
       for(isum = 0; isum <*nquantiles; isum++){
         *(sumstatout+isum) = *(sumstat+isum+(*nquantiles*(j))) - *(sumstat+isum+(*nquantiles*(0)));
       }
        costfunction(sumstatout, j, 0, nquantiles, n, &segcost, MBIC);
        lastchangelike[j] = segcost;
     }


    for(j=*minseglen;j<(2*(*minseglen));j++){
      lastchangecpts[j] = 0;
    }
    for(j=*minseglen;j<(2*(*minseglen));j++){
      numchangecpts[j] =1;
    }

    nchecklist=2;
    checklist[0]=0;
    checklist[1]=*minseglen;



    for(tstar=2*(*minseglen);tstar<(*n+1);tstar++){
      for(i=0;i<(nchecklist);i++){
        for(isum = 0; isum <*nquantiles; isum++){
          *(sumstatout+isum) = *(sumstat+isum+(*nquantiles*(tstar))) - *(sumstat+isum+(*nquantiles*(checklist[i])));
        }
        costfunction(sumstatout, tstar, checklist[i], nquantiles, n, &segcost, MBIC);
        tmplike[i] = lastchangelike[checklist[i]] + segcost + *pen;
      }

      min_which(tmplike,&nchecklist,&minout,&whichout); //updates minout and whichout with min and which element
      lastchangelike[tstar]=minout;
      lastchangecpts[tstar]=checklist[whichout];
      numchangecpts[tstar]=numchangecpts[lastchangecpts[tstar]]+1;
      // Update checklist for next iteration, first element is next tau
      nchecktmp=0;
      for(i=0;i<nchecklist;i++){
        if(tmplike[i]<= (lastchangelike[tstar]+*pen)){
          *(checklist+nchecktmp)=checklist[i];
          nchecktmp+=1;
        }
      }
      nchecklist = nchecktmp;
      *(checklist+nchecklist)=tstar-(*minseglen-1); // atleast 1 obs per seg
      nchecklist+=1;
    } // end taustar

    // put final set of changepoints together
    int ncpts=0;
    int last=*n;
    while(last!=0){
      *(cptsout + ncpts) = last;
      last=lastchangecpts[last];
      ncpts+=1;
    }
    err2:  free(checklist);
    err1:  return;
  }
