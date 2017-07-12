// KRUSKAL WALLACE *************************************************************

/* Kruskal-Wallis and normal scores distributions
 approximated with beta. 
 */
double KruskalWallisMaxU(
    int c,
    int n
)
{
  return double(c-1)+1.0/double(n-(c-1));
}

// Probability function for R
DISTS_API void pKruskalWallisR(
    double *Hp,
    int *cp,
    int *np,
    double *Up,
    int *doNormalScorep,
    int *Np,
    double *valuep
)
{
  int N=*Np;
  int i;
  
  for (i=0;i<N;i++)
    valuep[i]=pKruskal_Wallis(Hp[i],cp[i],np[i],Up[i],(bool)doNormalScorep[i]);
  
}


double pKruskal_Wallis(
    double H,  // The Statistic
    int c,   // number of treatments
    int n,   // Total number of observations
    double U,  // Sum (1/ni), where ni is numb obs for each treatment
    bool doNormalScore  // do normal scores
)
{
  double C=c;
  double N=n;
  
  if (H<0.0 || U<=0 || U>KruskalWallisMaxU(c,n))
    return NA_REAL;
  double V=(doNormalScore)?varNormalScores(N,C,U):varKruskal_Wallis(N,C,U);
  double d=((N-C)*(C-1.0)-V)/((N-1.0)*V);
  double f1=(C-1.0)*d;
  double f2=(N-C)*d;
  return pbeta(H/(N-1.0),f1,f2,true,false);
}

DISTS_API void uKruskalWallisR(
    double *Hp,
    int *cp,
    int *np,
    double *Up,
    int *doNormalScorep,
    int *Np,
    double *valuep
)
{
  int N=*Np;
  int i;
  
  for (i=0;i<N;i++)
    valuep[i]=qKruskal_Wallis(Hp[i],cp[i],np[i],Up[i],(bool)doNormalScorep[i]);
  
}


double qKruskal_Wallis(
    double H,  // The Statistic
    int c,   // number of treatments
    int n,   // Total number of observations
    double U,  // Sum (1/ni), where ni is numb obs for each treatment
    bool doNormalScore  // do normal scores
)
{
  
  if (H<0.0 || U<=0 || U>KruskalWallisMaxU(c,n))
    return NA_REAL;
  
  return 1.0-pKruskal_Wallis(H,c,n,U,doNormalScore);
}

DISTS_API void qKruskalWallisR(
    double *pp,
    int *cp,
    int *np,
    double *Up,
    int *doNormalScorep,
    int *Np,
    double *valuep
)
{
  int N=*Np;
  int i;
  
  for (i=0;i<N;i++)
    valuep[i]=xKruskal_Wallis(pp[i],cp[i],np[i],Up[i],(bool)doNormalScorep[i]);
} 

double xKruskal_Wallis(
    double P,
    int c,   // number of treatments
    int n,   // Total number of observations
    double U,  // Sum (1/ni), where ni is numb obs for each treatment
    bool doNormalScore  // do normal scores
)
{
  double C=c;
  double N=n;
  
  if (0>P || P>1 || U<=0 || U>KruskalWallisMaxU(c,n))
    return NA_REAL;
  
  
  double V=(doNormalScore)?varNormalScores(N,C,U):varKruskal_Wallis(N,C,U);
  double d=((N-C)*(C-1.0)-V)/((N-1.0)*V);
  double f1=(C-1.0)*d;
  double f2=(N-C)*d;
  
  return (N-1.0)*qbeta(P,f1,f2,true,false);
}

// Statistics for Kruskal_Wallis and Normal Scores
DISTS_API void sKruskalWallisR(
    int *cp,
    int *np,
    double *Up,
    int *doNormalScorep,
    int *Np,
    double *varp,
    double *modep,
    double *thirdp,
    double *fourthp
)
{
  int N=*Np;
  int i;
  
  for (i=0;i<N;i++) {
    sKruskal_Wallis(cp[i],np[i],Up[i],(bool)doNormalScorep[i],modep+i,thirdp+i,fourthp+i);
    if (Up[i]<=0 || Up[i]>KruskalWallisMaxU(cp[i],np[i])){
      varp[i]=NA_REAL;
    }
    else {
      varp[i]=(doNormalScorep[i])?varNormalScores(np[i],cp[i],Up[i]):varKruskal_Wallis(np[i],cp[i],Up[i]);
    }
  }
} 

void sKruskal_Wallis(
    int c,
    int n,
    double U,
    bool doNormalScore,
    double *mode,
    double *third, // Central moments
    double *fourth
)
{
  int nPoints=128;
  
  if (U<=0) {
    *mode=NA_REAL;
    *third=NA_REAL;
    *fourth=NA_REAL;
  }
  else {
    double minH=xKruskal_Wallis(0.01,c,n,U,doNormalScore);
    double maxH=xKruskal_Wallis(0.99,c,n,U,doNormalScore);
    double md=0.0;
    double maxValue=0.0;
    double m3=0.0;
    double m4=0.0;
    double sum=0.0;
    double mean=(double)(c-1);
    double delta=(maxH-minH)/(double)(nPoints-1);
    double H=minH;
    while (nPoints--) {
      double val=fKruskal_Wallis(H,c,n,U,doNormalScore);
      if (maxValue<val) {
        maxValue=val;
        md=H;
      }
      sum+=val;
      double h=H-mean;
      m3+=val*h*h*h;
      m4+=val*h*h*h*h;
      H+=delta;
    }
    m3/=sum;
    m4/=sum;
    *mode=md;
    *third=m3;
    *fourth=m4;
  }
}


DISTS_API void dKruskalWallisR(
    double *Hp,
    int *cp,
    int *np,
    double *Up,
    int *doNormalScorep,
    int *Np,
    double *valuep
)
{
  int N=*Np;
  int i;
  
  for (i=0;i<N;i++)
    valuep[i]=fKruskal_Wallis(Hp[i],cp[i],np[i],Up[i],(bool)doNormalScorep[i]);
} 

// Use numerical derivitive
double fKruskal_Wallis(
    double H,  // The Statistic
    int c,   // number of treatments
    int n,   // Total number of observations
    double U,  // Sum (1/ni), where ni is numb obs for each treatment
    bool doNormalScore  // do normal scores
)
{
  double delta=0.001;
  return (pKruskal_Wallis(H+delta,c,n,U,doNormalScore)-
          pKruskal_Wallis(H,c,n,U,doNormalScore))/delta;
}


/*
 Kruskal Wallis tau random numbers
 */
DISTS_API void rKruskalWallisR(
    double *randArrayp,
    int *Np,
    int *Mp, // size of parameter arrays
    int *cp,
    int *np,
    double *Up,
    bool *doNormalScorep
)
{
  int N=*Np;
  int M=*Mp;
  int D;
  int j;
  int k;
  int loc;
  int cloc;
  double *tArray;
  
  if (M==1) 
    rKruskal_Wallis(randArrayp,N,*cp,*np,*Up,(bool)*doNormalScorep);
  else { // Allow for random values for each element of nu and lambda
    D=(N/M)+((N%M)?1:0);
    tArray=(double *)S_alloc((long)D,sizeof(double));
    loc=0;
    for (j=0;j<M;j++) {
      rKruskal_Wallis(tArray,D,cp[j],np[j],Up[j],doNormalScorep[j]);
      for (k=0;k<D;k++) {
        cloc=loc+k*M;
        if (cloc<N)
          randArrayp[cloc]=tArray[k];
        else break;
      }
      loc++;
    }
  }
  
}

void rKruskal_Wallis(
    double* randArray,
    int N,  // number of samples
    int c,   // number of treatments
    int n,   // Total number of observations
    double U,  // Sum (1/ni), where ni is numb obs for each treatment
    bool doNormalScore  // do normal scores
)
{
  GetRNGstate();
  
  for (int i=0;i<N;i++){
    randArray[i]=(double)xKruskal_Wallis(unif_rand(),c,n,U,doNormalScore);
  }
  
  PutRNGstate();
}


/* Kruskal-Wallis variance
 Wallace, DL. (1959). Simplified beta-approximations to
 the Kruskal-Wallis H test. JASA 54, 225-230. Equation 6.2 */

double varKruskal_Wallis(
    double N,
    double C,
    double U
)
{
  double f=2.0*(C-1.0);
  f-=0.4*(3.0*C*C-6.0*C+N*(2.0*C*C-6.0*C+1.0))/(N*(N+1.0));
  f-=1.2*U;
  
  return f;
}


/* Normal scores variance. See LU, HT & Smith, PJ. (1979).
 Distribution of the normal scores statistic for nonparametric
 one-way analysis of variance. JASA. 74, 715-772 */
double varNormalScores(
    double N,
    double C,
    double U
)
{
  double alpha=0.375;
  
  double NP=N+1.0;
  double NM=N-1.0;
  double NC=N-C;
  double CM=C-1.0;
  
  double den=1.0-2.0*alpha;
  long n=(long)(0.1+N/2.0);
  
  double e2=0.0;
  double e4=0.0;
  for (long i=1L;i<=n;i++) {
    double e=qnorm(((double)i-alpha)/(N+den),0,1,true,false);
    e*=e;
    e2+=e;
    e4+=e*e;
  }
  
  e2*=4.0*e2;
  double g=(NP*N*NM*NM*2.0*e4-3.0*NM*NM*NM*e2)/(NM*(N-2.0)*(N-3.0)*e2);
  double v=2*CM*NC/NP;
  v-=g*(NP*C*C+2.0*CM*NC-N*NP*U)/(N*NP);
  
  return(v);
}


