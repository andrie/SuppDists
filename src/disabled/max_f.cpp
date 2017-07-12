/****************************************************************************
 Maximum F-Ratio
*/



/*
Lower tail probability by numerical integration using Romberg's rule -- see
Acton, p206
The required probability is
N\int_{0}^{\infty}p(x)[P(Fx)-P(x)]^{N-1}, where p(x) is the chisquare
density and P() its integral.
*/

// The integrand
static double pmaxFRatioIntegrand(
    double x,
    double F,
    int df,
    int N,
    double logC
)
{
  double logV=logC-0.5*x;
  logV+=(0.5*(double)df-1.0)*log(x)+(double)(N-1)*log(pchisq(x*F,df,true,false)-pchisq(x,df,true,false));
  return exp(logV);
}


// Probability function for R
DISTS_API void pmaxFratioR(
    double *Fp,
    int *dfp,
    int *np, // Number of mean squares
    int *Np,
    double *valuep
)
{
  int N=*Np;
  int i;
  
  for (i=0;i<N;i++)
    valuep[i]=pmaxfratio(Fp[i],dfp[i],np[i]);
  
}

double pmaxfratio(
    double F, // Maximum F ratio
    int df,
    int N // Number of mean squares
)
{
  if (F<=0 || df<1 || N<1)
    return NA_REAL;
  
  if (N equals 2) {
    return 1-2*(1.0-pf(F,(double)df,(double)df,true,false));
  }
  
  // Calculate a root of the constant to scale x in pmaxFRatioIntegrand()
  double logC=log((double)(N))-(0.5*(double)df)*LOG2-loggamma(0.5*(double)df);
  
  double upperLimit=qchisq(0.9999,df,true,false);
  double lowLimit=qchisq(0.0001,df,true,false);
  double h=upperLimit-lowLimit; // range of integration
  
  const double MRatioTol=1e-4; // result to this accuracy
  
  const int maxiterate=16;  // Maximum interates allowed -- fcn evaluated at most 2^16 pts
  double A[maxiterate][maxiterate]; // the triangular array of Romberg iterates
  // stop when two diagonal values agree to MRatioTol
  
  double delta=h;   // Initial StepTheKey
  double value=0.0;   // Converged value
  
  A[0][0]=(h/2.0)*(pmaxFRatioIntegrand(lowLimit,F,df,N,logC)+pmaxFRatioIntegrand(upperLimit,F,df,N,logC)); 
  
  double twoPower=1.0;
  int k=0;
  int n=1;
  repeat
    k++;
  delta*=0.5;  // Evaluate pmaxFRatioIntegrand() at half the previous interval
  if (k>1) {
    n*=2;   // Number of new ordinates
  }
  
  twoPower*=2.0;
  double sum=0.0;
  double z=upperLimit-delta; // Start with this value
  int m=n;
  while (m--) {
    double value=pmaxFRatioIntegrand(z,F,df,N,logC);
    
    sum+=value;
    z-=2.0*delta;  // Every other ordinate is a new one
  }
  
  A[0][k]=A[0][k-1]/2.0 + h*sum/twoPower; // Pool old and new ordinate sum
  
  // Fill out the Romberg triangle
  double fourPower=1.0;
  for (int i=1;i<=k;i++) {
    fourPower*=4.0;
    A[i][k-i]=(fourPower*A[i-1][k-i+1]-A[i-1][k-i])/(fourPower-1.0); 
  }
  // Check for convergence
  
  value=A[k][0];
  if (fabs((value - A[k-1][0])/value)<MRatioTol) {
    break;
  }
  until(k>=maxiterate-1);
  
  
  
  return value;
}

/*
Upper tail of maximum F ratio
*/

DISTS_API void umaxFratioR(
    double *Fp,
    int *dfp,
    int *np, // Number of mean squares
    int *Np,
    double *valuep
)
{
  int N=*Np;
  int i;
  
  for (i=0;i<N;i++)
    valuep[i]=qmaxfratio(Fp[i],dfp[i],np[i]);
  
}

double qmaxfratio(
    double F,
    int df,
    int N // Number of mean squares
)
{
  if (F<=0 || df<1 || N<1)
    return NA_REAL;
  
  return 1.0-pmaxfratio(F,df,N);
}


/*
Derivitive of maximum F ratio
*/

static double fmaxFRatioIntegrand(
    double x,
    double F,
    int df,
    int N,
    double logC
)
{
  double logV=logC-0.5*x*(1.0+F);
  logV+=((double)df-1.0)*log(x)+(0.5*(double)df-1.0)*log(F)+
    (double)(N-2)*log(fabs(pchisq(x*F,df,true,false)-pchisq(x,df,true,false))); // fabs because F can be less than 1
  return exp(logV);
}


DISTS_API void dmaxFratioR(
    double *Fp,
    int *dfp,
    int *np, // Number of mean squares
    int *Np,
    double *valuep
)
{
  int N=*Np;
  int i;
  
  for (i=0;i<N;i++)
    valuep[i]=fmaxfratio(Fp[i],dfp[i],np[i]);
  
}

double fmaxfratio(
    double F, // Maximum F ratio
    int dgf,
    int N // Number of mean squares
)
{
  if (F<=0 || dgf<1 || N<1)
    return NA_REAL;
  
  if (N equals 2) {
    return (2.0*df(F,(double)dgf,(double)dgf,false));
  }
  // Calculate a root of the constant to scale x in fmaxFRatioIntegrand()
  double logC=log((double)(N*(N-1)))-(double)dgf*LOG2-2.0*loggamma((double)dgf/2.0);
  
  double upperLimit=qchisq(0.9999,dgf,true,false);
  double lowLimit=qchisq(0.0001,dgf,true,false);
  double h=upperLimit-lowLimit; // range of integration
  
  const double MRatioTol=1e-4; // result to this accuracy
  
  const int maxiterate=16;  // Maximum interates allowed -- fcn evaluated at most 2^16 pts
  double A[maxiterate][maxiterate]; // the triangular array of Romberg iterates
  // stop when two diagonal values agree to MRatioTol
  
  double delta=h;   // Initial StepTheKey
  double value=0.0;   // Converged value
  
  A[0][0]=(h/2.0)*(fmaxFRatioIntegrand(lowLimit,F,dgf,N,logC)+fmaxFRatioIntegrand(upperLimit,F,dgf,N,logC)); 
  
  double twoPower=1.0;
  int k=0;
  int n=1;
  repeat
    k++;
  delta*=0.5;  // Evaluate fmaxFRatioIntegrand() at half the previous interval
  if (k>1) {
    n*=2;   // Number of new ordinates
  }
  
  twoPower*=2.0;
  double sum=0.0;
  double z=upperLimit-delta; // Start with this value
  int m=n;
  while (m--) {
    double value=fmaxFRatioIntegrand(z,F,dgf,N,logC);
    
    sum+=value;
    z-=2.0*delta;  // Every other ordinate is a new one
  }
  
  A[0][k]=A[0][k-1]/2.0 + h*sum/twoPower; // Pool old and new ordinate sum
  
  // Fill out the Romberg triangle
  double fourPower=1.0;
  for (int i=1;i<=k;i++) {
    fourPower*=4.0;
    A[i][k-i]=(fourPower*A[i-1][k-i+1]-A[i-1][k-i])/(fourPower-1.0); 
  }
  // Check for convergence
  
  value=A[k][0];
  if (fabs((value - A[k-1][0])/value)<MRatioTol) {
    break;
  }
  until(k>=maxiterate-1);
  
  return value;
}




/*
Inverse maximum F ratio
Uses Newton iterations
*/

static const int MAXFITER=20;
// Parameters for Johnson fit for df=2,4,8,16,32,64 and K=3,6,9,12
static JohnsonParms parmArray[7][4]={
  {{-1.019,0.510,2.750,0.915,SU},{-1.000,0.512,7.904,2.941,SU},{-0.996,0.512,14.124,5.302,SU},{-0.995,0.512,21.008,7.887,SU}},
  {{-1.542,0.834,1.425,0.616,SU},{-1.403,0.857,2.919,1.519,SU},{-1.372,0.862,4.317,2.235,SU},{-1.358,0.865,5.612,2.858,SU}},
  {{-2.892,1.252,0.927,0.271,SU},{-2.043,1.302,1.625,0.872,SU},{-1.931,1.314,2.184,1.175,SU},{-1.883,1.319,2.649,1.395,SU}},
  {{4.014,1.447,0.862,15.467,SB},{-3.054,1.848,1.127,0.523,SU},{-2.712,1.860,1.439,0.705,SU},{-2.571,1.862,1.681,0.814,SU}},
  {{2.995,1.509,0.896,4.930,SB},{-5.061,2.518,0.884,0.268,SU},{-3.699,2.480,1.126,0.457,SU},{-3.524,2.497,1.261,0.514,SU}},
  {{2.439,1.528,0.927,2.357,SB},{8.346,2.888,0.832,13.668,SB},{-6.664,3.306,0.907,0.209,SU},{-4.496,3.145,1.083,0.348,SU}},
  {{2.326,1.620,0.937,1.515,SB},{8.411,3.231,0.851,7.635,SB},{-7.810,3.947,0.880,0.163,SU},{-4.428,3.497,1.062,0.277,SU}}};

static const double ln2=0.69314718055994172321;

JohnsonParms GetClosestJohnsonParms(
    int df,
    int N
)
{
  int row;
  int col;
  double ddf=(double)df;
  double dN=(double)N;
  dN/=3;
  col=int(floor(0.5+dN))-1;
  col=maxm(0,col);
  col=minm(3,col);
  
  ddf=log(ddf)/ln2;  // log base 2 of ddf
  row=int(floor(0.5+ddf))-1;
  row=maxm(0,row);
  row=minm(6,row);
  return parmArray[row][col];
}

DISTS_API void qmaxFratioR(
    double *pp,
    int *dfp,
    int *np, // Number of mean squares
    int *Np,
    double *valuep
)
{
  int N=*Np;
  int i;
  
  for (i=0;i<N;i++)
    valuep[i]=xmaxfratio(pp[i],dfp[i],np[i]);
  
}

double xmaxfratio(
    double p,
    int dgf,
    int N
)
{
  if (p<0 || p>1 || dgf<1 || N<1)
    return NA_REAL;
  
  // This dist is a folded F
  if (N equals 2) {
    return qf(1.0-0.5*(1.0-p),(double)dgf,(double)dgf,true,false);
  }
  
  // Can actually do a little better, especiall for larger df, but this is safe
  if (dgf>160 || N>24)
    return NA_REAL;
  
  JohnsonParms parms=GetClosestJohnsonParms(dgf,N);
  
  double x=xjohnson(p,parms);
  x=maxm(1.000001,x); // Johnson sometimes comes up with a value less than 1
  
  bool more=true;
  double h;
  double ho=1e6;
  int m=0;
  repeat
    h=(p-pmaxfratio(x,dgf,N))/fmaxfratio(x,dgf,N);
  x+=h;
  more=(fabs(h/x)>TOLNEWTON);
  if (ho<fabs(h)) {
    more=false;
    x-=h;
  }
  ho=fabs(h);
  until(m++>MAXFITER || ! more);
  
  // The max F Ratio must be at least 1
  if (x<1.0) {
    x=1.000001;
  }
  
  return x;
}

void rmaxFratio(
    double *randomArray,
    int N,
    int df,
    int n,
    double *tArray
)
{
  int i;
  int j;
  double maxV;
  double minV;
  
  for (i=0;i<N;i++) {
    if (df<1 || n<1) {
      randomArray[i]=NA_REAL;
    }
    else {
      rdchisq(tArray,n,df);
      maxV=-1.0;
      minV=1e20;
      for (j=0;j<n;j++) {
        maxV=(tArray[j]>maxV)?tArray[j]:maxV;
        minV=(tArray[j]<minV)?tArray[j]:minV;
      }
      randomArray[i]=maxV/minV;
    }
  }
  
}
// Random function for R
DISTS_API void rmaxFratioR(
    int *dfp,
    int *np, // Number of mean squares
    int *Np,
    int *Mp,
    double *valuep
)
{
  int N=*Np;
  int n=0;
  int M=*Mp;
  int D;
  int j;
  int k;
  int loc;
  int cloc;
  double *tArray;
  double *pArray;
  
  if (M==1) {
    pArray=(double *)S_alloc((long)*np,sizeof(double));
    rmaxFratio(valuep,N,*dfp,*np,pArray);
  }
  else { // Allow for random values for each element of nu and lambda
    for (j=0;j<M;j++)
      n=(n>np[j])?n:np[j];
    pArray=(double *)S_alloc((long)n,sizeof(double));
    D=(N/M)+((N%M)?1:0);
    tArray=(double *)S_alloc((long)D,sizeof(double));
    loc=0;
    for (j=0;j<M;j++) {
      rmaxFratio(tArray,D,dfp[j],np[j],pArray);
      for (k=0;k<D;k++) {
        cloc=loc+k*M;
        if (cloc<N)
          valuep[cloc]=tArray[k];
        else break;
      }
      loc++;
    }
  }
  
}





DISTS_API void smaxFratioR(
    int *dfp,
    int *np, // Number of mean squares
    int *Np,
    double *mean,
    double *median,
    double *mode,
    double *variance,
    double *third,
    double *fourth
)
{
  int N=*Np;
  int i;
  
  for (i=0;i<N;i++)
    smaxFratio(dfp[i],np[i],mean+i,median+i,mode+i,variance+i,third+i,fourth+i);
}


static int gdf;
static int gk;
static double gmean;

static double MeanFcn(double x){return x*fmaxfratio(x,gdf,gk);}
static double VarianceFcn(double x){return (x-gmean)*(x-gmean)*fmaxfratio(x,gdf,gk);}
static double ThirdMomentFcn(double x){return (x-gmean)*(x-gmean)*(x-gmean)*fmaxfratio(x,gdf,gk);}
static double FourthMomentFcn(double x){return (x-gmean)*(x-gmean)*(x-gmean)*(x-gmean)*fmaxfratio(x,gdf,gk);}
static double AFunction(double x){return fmaxfratio(x,gdf,gk);}

void smaxFratio(
    int df,
    int N, // Number of mean squares
    double *mean,
    double *median,
    double *mode,
    double *variance,
    double *third,
    double *fourth
)
{
  gdf=df;
  gk=N;
  
  if (df<1 || N<1) {
    *mean=NA_REAL;
    *median=NA_REAL;
    *mode=NA_REAL;
    *variance=NA_REAL;
    *third=NA_REAL;
    *fourth=NA_REAL;
  }
  else {
    double lowX=xmaxfratio(0.01,df,N);
    double highX=xmaxfratio(0.99,df,N);
    
    *mean=FindDistributionStatistic(lowX,highX,MeanFcn);
    gmean=*mean;
    *median=xmaxfratio(0.5,df,N);
    *mode=FindDistributionMode(lowX,highX,AFunction);
    *variance=FindDistributionStatistic(lowX,highX,VarianceFcn);
    *third=FindDistributionStatistic(lowX,highX,ThirdMomentFcn);
    *fourth=FindDistributionStatistic(lowX,highX,FourthMomentFcn);
  }
}

