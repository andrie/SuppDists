// KENDALL's TAU ***************************************************************

/* pkendall -- Kendall's tau probabilities using Edgeworth
expansion. */
static const double CONTCORRK=0.5;

static void fills(
    int *a,
    int k,
    int n
)
{
  int sum=0;
  int i;
  // Sum last n in current array
  for (i=0;i<n;i++) {
    int c=k-i;
    if (0<=c) {
      sum+=a[c];
    }
  }
  // Successively replace a[i] with sum of a[i], a[i-1], ..., a[i-n+1]
  for (i=k;i!=0;i--) {
    long temp=a[i];
    a[i]=sum;
    sum-=temp;
    int c=i-n;
    if (0<=c) {
      sum+=a[c];
    }
  }
  
  // a[0] is always 1
}

/* kendexact -- exact distribution of Kendall's tau. See
Hajek, J. & Sidak, Z. (1967) Theory of rank
tests. Academic Press, NY p140. Should not
be used if numbers will overflow a long --
thus good to n=12 
Note: Although the mean must be zero, the 50th percentile
is not be zero for even n, since M/2 is not integral,
where M=n(n-1)/2, the max value of k.
*/


static double kendexact(
    int N,
    int T,
    bool density // When true, returns the prob at T
)
{
  // a is an array containing n! times the individual probabilities
  // for j=0, ..., k=n(n-1)/2.
  // fills() updates the current row to the one for n.
  
  int *a=(int *)S_alloc((long)(T+1),sizeof(int));
  memset(a,0,(T+1)*sizeof(int));
  
  a[0]=1; // the row for n=1 is 1,0,0,0,0 ...
  int k=1;
  for (int n=2;n<=N;n++) {
    k=(k<=T)?k:T; // No need to calculate beyond T
    fills(a,k,n);
    k+=n; // this will be n(n-1)/2 after n is incremented
  }
  int sum=0;
  if (density) {
    sum=a[T];
  }
  else {
    for (int i=0;i<=T;i++) {
      sum+=a[i];
    }
  }
  
  
  return exp(log((double)sum)-loggamma((double)(N+1)));
}


DISTS_API void pKendallR(
    int *nip,
    double *taup,
    int *Np,
    double *valuep
)
{
  int N=*Np;
  int i;
  
  for (i=0;i<N;i++)
    valuep[i]=pkendall(nip[i],taup[i]);
}


double pkendall(
    int ni,
    double tau
) 
{
  double n=ni;
  double M=(n*(n-1.0))/2.0;
  double mean=M/2.0;
  int k=(int)(.5+(1.0+tau)*mean);
  double t=(double)k;
  
  if (tau>1 || tau <-1 || ni<2)
    return NA_REAL;
  
  if (k<0L) {
    return 0.0;
  }
  if (t>M) {
    return 1.0;
  }
  if (ni<=MAXKENDALEXACT) {
    return kendexact(ni,k,false);
  }
  
  double s=n*(n+1)*(2.0*n+1.0)/6.0;
  double f=((n+1)*3.0*n-1)/5.0;
  double g=(((n*n+2.0)*n-1.0)*3.0*n+1)/7.0;
  double sd=s-n;
  double l4=-1.2*(s*f-n)/(sd*sd);
  double l6=(48.0/7.0)*(s*g-n)/(sd*sd*sd);
  sd=sqrt(sd/12.0);
  double x=(t+CONTCORRK-mean)/sd;
  double z=phi0(x);
  double P=pnorm(x,0,1,true,false);
  P+=((35*l4*l4*phi7(x,z)/56.0+l6*phi5(x,z))/30.0+l4*phi3(x,z))/24.0;
  
  return P;
}

DISTS_API void fourthKendallR(
    int *nip,
    int *Np,
    double *valuep
)
{
  int N=*Np;
  int i;
  
  for (i=0;i<N;i++)
    valuep[i]=fourthkendall(nip[i]);
  
}
// Fourth central moment
double fourthkendall(
    int ni
)
{
  if (ni<4) {
    return NA_REAL;
  }
  else {
    int nPoints=128;
    double m4=0.0;
    double sum=0.0;
    double minTau=xkendall(0.01,ni);
    double maxTau=xkendall(0.99,ni);
    double delta=(maxTau-minTau)/(double)(nPoints-1);
    double tau=minTau;
    while (nPoints--) {
      double val=fkendall(ni,tau);
      m4+=val*tau*tau*tau*tau;
      sum+=val;
      tau+=delta;
    }
    m4/=sum;
    return m4;
  }
}

DISTS_API void dKendallR(
    int *nip,
    double *taup,
    int *Np,
    double *valuep
)
{
  int N=*Np;
  int i;
  
  for (i=0;i<N;i++)
    valuep[i]=fkendall(nip[i],taup[i]);
}

double fkendall(
    int ni,
    double tau
)
{
  double n=ni;
  double M=(n*(n-1.0))/2.0;
  double mean=M/2.0;
  int k=(int)(.5+(1.0+tau)*mean);
  double t=(double)k;
  
  if (tau>1 || tau <-1)
    return NA_REAL;
  
  if (k<0L) {
    return 0.0;
  }
  if (t>M) {
    return 0.0;
  }
  if (ni<=MAXKENDALEXACT) {
    return kendexact(ni,k,true);
  }
  
  double s=n*(n+1)*(2.0*n+1.0)/6.0;
  double f=((n+1)*3.0*n-1)/5.0;
  double g=(((n*n+2.0)*n-1.0)*3.0*n+1)/7.0;
  double sd=s-n;
  double l4=-1.2*(s*f-n)/(sd*sd);
  double l6=(48.0/7.0)*(s*g-n)/(sd*sd*sd);
  sd=sqrt(sd/12.0);
  
  double x=(t+CONTCORRK-mean)/sd;
  double z=phi0(x);
  double P1=pnorm(x,0,1,true,false);
  P1+=((35*l4*l4*phi7(x,z)/56.0+l6*phi5(x,z))/30.0+l4*phi3(x,z))/24.0;
  
  t-=1.0;
  double P0=0;
  if (t>=0.0) {
    x=(t+CONTCORRK-mean)/sd;
    z=phi0(x);
    P0=pnorm(x,0,1,true,false);
    P0+=((35*l4*l4*phi7(x,z)/56.0+l6*phi5(x,z))/30.0+l4*phi3(x,z))/24.0;
  }
  
  return P1-P0;
}

/* qkendall -- Upper tail of Kendall's tau
Retruns Pr(T>=tau|n).*/

DISTS_API void uKendallR(
    int *nip,
    double *taup,
    int *Np,
    double *valuep
)
{
  int N=*Np;
  int i;
  
  for (i=0;i<N;i++)
    valuep[i]=qkendall(nip[i],taup[i]);
}

double qkendall(
    int n,
    double tau
)
{
  if (tau>1 || tau <-1 || n<2)
    return NA_REAL;
  
  return (1.0-pkendall(n,tau));
}

/* xkendall -- returns smallest tau such that pr<=Pr(T<=tau|n). */
DISTS_API void qKendallR(
    int *nip,
    double *pp,
    int *Np,
    double *valuep
)
{
  int N=*Np;
  int i;
  
  for (i=0;i<N;i++)
    valuep[i]=xkendall(pp[i],nip[i]);
}

double xkendall(
    double pr,
    int ni
)
{
  double n=ni;
  double mean=0.25*n*(n-1);
  double sd=n*(2*n+1)*(n+1)/6-n;
  sd=sqrt(sd/12.0);
  long k=(long)(.5+mean+sd*qnorm(pr,0,1,true,false));
  
  double tau=4.0*(double)k/(n*(n-1.0))-1.0;
  bool larger=(pr<=pkendall(ni,tau));
  
  if (0>=pr || pr>=1 || ni<2)
    return NA_REAL;
  
  
  while (larger) {
    if (! k) {
      goto theExit;
    }
    tau=4.0*(double)(--k)/(n*(n-1.0))-1.0;
    larger=(pr<=pkendall(ni,tau));
    if (! larger) {
      ++k;
      goto theExit;
    }
  }
  while (! larger) {
    tau=4.0*(double)(++k)/(n*(n-1.0))-1.0;
    larger=(pr<=pkendall(ni,tau));
    if (larger) {
      goto theExit;
    }
  }
  theExit:
    return 4.0*(double)k/(n*(n-1.0))-1.0; // Tau
}

/*
Kendall tau random numbers
*/

DISTS_API void rKendallR(
    int *nip,
    int *Np,
    int *Mp,
    double *valuep
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
    rkendall(valuep,N,*nip);
  else { // Allow for random values for each element of nu and lambda
    D=(N/M)+((N%M)?1:0);
    tArray=(double *)S_alloc((long)D,sizeof(double));
    loc=0;
    for (j=0;j<M;j++) {
      rkendall(tArray,D,nip[j]);
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

void rkendall(
    double* randArray,
    int N,  // number of samples
    int ni
)
{
  GetRNGstate();
  
  for (int i=0;i<N;i++){
    randArray[i]=(double)xkendall(unif_rand(),ni);
  }
  
  PutRNGstate();
}



