/**************************************************************************************** 
 Hypergeometric distribution 
Given a total of N items, n of which are marked, select a sample of size S, and
find the probability that 

fhypergeometric: that there are x marked items in the sample -- frequency
phypergeometric: that there are x or less in the sample -- cdf.
*/


// Normal approximation to the hypergeometric distribution function due to Peizer
// See Ling, R.F. and Pratt, J.W. (1984) The accuracy of Peizer approximations
// to the hypergeometric distribution, with comparisons to some other 
// approximations. JASA 79-385. 49-60.
double PeizerHypergeometric(
    int x,  // Number of marked items in sample
    int S,  // Sample size
    int n,  // Total number of marked items
    int N  // Total number of items
)
{
  const double oneSix=1.0/6.0;
  
  double dn=(double)n;
  double dm=(double)(N-n);
  double dr=(double)S;
  double ds=(double)(N-S);
  double dN=(double)N;
  double dnp=dn+oneSix;
  double dmp=dm+oneSix;
  double drp=dr+oneSix;
  double dsp=ds+oneSix;
  double dNp=dN-oneSix;
  double A=(double)x+0.5;
  double B=maxm(dn-A,0.5); // prevents B or C from going neg when x=n or x=S
  double C=maxm(dr-A,0.5);
  double D=(dm-dr)+A;
  double Ap=A+oneSix+0.02/(A+0.5)+0.01/(dn+1.0)+0.01/(dr+1.0);
  double Bp=B-oneSix+0.02/(B+0.5)+0.01/(dn+1.0)+0.01/(ds+1.0);
  double Cp=C-oneSix+0.02/(C+0.5)+0.01/(dm+1.0)+0.01/(dr+1.0);
  double Dp=D+oneSix+0.02/(D+0.5)+0.01/(dm+1.0)+0.01/(ds+1.0);
  
  double L=A*log((A*dN)/(dn*dr))+B*log((B*dN)/(dn*ds))+C*log((C*dN)/(dm*dr))+D*log((D*dN)/(dm*ds));
  
  double z=((Ap*Dp-Bp*Cp)/fabs(A*D-B*C))*sqrt(2.0*L*((dm*dn*dr*ds*dNp)/(dmp*dnp*drp*dsp*dN)));
  
  return pnorm(z,0,1,true,false);
  
}

char *hyperNames[]= {
  (char *)"classic",
  (char *)"IAi",
  (char *)"IAii",
  (char *)"IB",
  (char *)"IIA",
  (char *)"IIB",
  (char *)"IIIA",
  (char *)"IIIB",
  (char *)"IV",
  (char *)"no type"
};

// Returns true if the double is an int
bool isint(
    double x
)
{
  return x equals floor(x);
} 


// Finds the type of hypergeometric

hyperType typeHyper(
    double a,     // Sample size
    double m,  // Total number of marked items
    double N   // Total number of items
)
{
  
  hyperType variety;
  
  
  
  if (0.0<a && 0.0<N && 0.0<m && isint(a) && isint(N) && isint(m)) {
    variety=classic;
  }
  
  else
    if (0.0<a && 0.0<N && 0.0<m && isint(m) && m-1.0<a && a<N-(m-1.0)) {
      variety=IAi;
    }
    else
      if (0.0<a && 0.0<N && 0.0<m && isint(a) && a-1.0<m && m<N-(a-1.0)) {
        variety=IAii;
      }
      else
        if (0.0<a && 0.0<N && 0.0<m && ! isint(a) && ! isint(m) && a+m-1.0<N &&
            floor(a) equals floor(m)) {
          variety=IB;  // Specified 1.0<N to avoid problems with small parameters
        }
        else
          if (a<0.0 && N<m+a-1.0 && 0.0<m && isint(m)) { //Kemp&Kemp use b<0 && b!=-1, Ben Bolker mod
            variety=IIA;
          }
          else
            if (a<0.0 && -1.0<N && N<m+a-1.0 && 0.0<m && ! isint(m) && 
                floor(m) equals floor(m+a-1.0-N)) {
              variety=IIB;
            }
            else
              if (0.0<a && N<m-1.0 && m<0.0 && isint(a)) {
                variety=IIIA;
              }
              else
                if (0.0<a && -1.0<N && N<a+m-1.0 && m<0.0 && ! isint(a) && 
                    floor(a) equals floor(a+m-1.0-N)) {
                  variety=IIIB;
                }
                else
                  if (a<0.0 && -1.0<N && m<0.0) { 
                    variety=IV;
                  }
                  else {
                    variety=noType;
                  }
                  
                  return variety;
}

bool checkHyperArgument(
    int k,
    double a,     // Sample size
    double m,  // Total number of marked items
    double N,   // Total number of items
    hyperType variety
)
{
  
  switch (variety) {
  case classic:
    return (maxm(0,(int)(a+m-N))<=k && k<=minm((int)a,(int)m));
    break;
  case IAi:
    return (0<=k && k<=(int)m);
    break;
  case IAii:
    return (0<=k && k<=(int)a);
    break;
  case IB:  // Specified 1.0<N to avoid problems with small parameters
    return (0<=k);
    break;
  case IIA:
    return (0<=k && k<=(int)m);
    break;
  case IIB:
    return (0<=k);
    break;
  case IIIA:
    return (0<=k && k<=(int)a);
    break;
  case IIIB:
    return (0<=k);
    break;
  case IV:
    return (0<=k);
    break;
  case noType:
    break;
  }
  return false;
}


// Reports the type and range to the user
DISTS_API void tghyperR(
    double *ap,     // Sample size
    double *mp,  // Total number of marked items
    double *Np,   // Total number of items
    char **aString
)
{
  double a=*ap;
  double m=*mp;
  double N=*Np;
  
  hyperType variety=typeHyper(a,m,N);
  
  
  switch (variety) {
  case classic:
    snprintf(*aString,127,"type = %s -- %d <= x <= %d",hyperNames[(int)classic],maxm(0,(int)(a+m-N)),minm((int)a,(int)m));
    break;
    
  case IAi:
    snprintf(*aString,127,"type = %s -- 0 <= x <= %d",hyperNames[(int)IAi],(int)m);
    break;
  case IAii:
    snprintf(*aString,127,"type = %s -- 0 <= x <= %d",hyperNames[(int)IAii],(int)a);
    break;
  case IB:  // Specified 1.0<N to avoid problems with small parameters
    snprintf(*aString,127,"type = %s -- x = 0,1,2,...",hyperNames[(int)IB]);
    break;
  case IIA:
    snprintf(*aString,127,"type = %s -- 0 <= x <= %d",hyperNames[(int)IIA],(int)m);
    break;
  case IIB:
    snprintf(*aString,127,"type = %s -- x = 0,1,2,...",hyperNames[(int)IIB]);
    break;
  case IIIA:
    snprintf(*aString,127,"type = %s -- 0 <= x <= %d",hyperNames[(int)IIIA],(int)a);
    break;
  case IIIB:
    snprintf(*aString,127,"type = %s -- x = 0,1,2,...",hyperNames[(int)IIIB]);
    break;
  case IV:
    snprintf(*aString,127,"type = %s -- x = 0,1,2,...",hyperNames[(int)IV]);
    break;
  case noType:
    snprintf(*aString,127,"type = %s",hyperNames[(int)noType]);
  }
  
  
}

DISTS_API void dghyperR(
    int *kp,
    double *ap,
    double *np,
    double *Np,
    int *Mp,
    double *valuep
)
{
  hyperType variety;
  int M=*Mp;
  int i;
  
  for (i=0;i<M;i++) {
    variety=typeHyper(ap[i],np[i],Np[i]);
    if (variety==classic)
      valuep[i]=fhypergeometric(kp[i],(int)ap[i],(int)np[i],(int)Np[i]);
    else if (variety!=noType)
      valuep[i]=fgenhypergeometric(kp[i],ap[i],np[i],Np[i],variety);
    else valuep[i]=NA_REAL;
  }
  
  
}


double fhypergeometric(
    int x,  // Number of marked items in sample
    int S,     // Sample size
    int n,  // Total number of marked items
    int N   // Total number of items
)
{
  double logP=loggamma((double)(n+1))+loggamma((double)(N-n+1))+loggamma((double)(S+1))+
    loggamma((double)(N-S+1));
  logP-=loggamma((double)(x+1))+loggamma((double)(n-x+1))+loggamma((double)(S-x+1))+
    loggamma((double)(N-S-n+x+1))+loggamma((double)(N+1));
  if (logP<-MAXEXP) {
    return 0.0;
  }
  else {
    return exp(logP);
  }
  
}


/*
Essentially, one has a table:

x y | n
? ? | N-n
----------
a N-a| N

(1) x ranges from a=max(0,a+n-N)) to b=min(n,a). If n<a, then the total range
is N-A, for A=a+n-N or n for A=0.
(2) One can interchange n and a if a is smaller, to minimize this range.
(3) If a<N-a, then A=0, otherwise, one can interchange x and y and a and N-a to
get A=0.
Proof: if a<N-a, then A=a+n-N<N-a+n-N=n-a, and n<a hence A=0 by interchange.
*/


DISTS_API void pghyperR(
    int *kp,
    double *ap,
    double *np,
    double *Np,
    int *Mp,
    double *valuep
)
{
  hyperType variety;
  int M=*Mp;
  int i;
  
  for (i=0;i<M;i++) {
    variety=typeHyper(ap[i],np[i],Np[i]);
    if (! checkHyperArgument(kp[i],ap[i],np[i],Np[i],variety))
      valuep[i]=NA_REAL;
    else if (variety==classic)
      valuep[i]=phypergeometric(kp[i],(int)ap[i],(int)np[i],(int)Np[i]);
    else
      valuep[i]=pgenhypergeometric(kp[i],ap[i],np[i],Np[i],variety);
  }
  
  
}

double phypergeometric(
    int x,  // Number of marked items in sample
    int a,     // Sample size
    int n,  // Total number of marked items
    int N   // Total number of items
)
{
  if (x<maxm(0,a-(N-n)) || x>minm(a,n)) 
    return NA_REAL;
  
  // interchange n and a to get the fewest terms to sum
  if (a<n) {
    int k=a;
    a=n;
    n=k;
  }
  
  if (x equals n) {
    return 1.0;
  }
  
  // Switch tails if necessesary to minimize number of terms to sum
  int xmin=maxm(0,n+a-N);
  bool lowerTail=true;
  if (x-xmin>n-x) {      
    x=n-x-1;
    a=N-a;
    xmin=maxm(0,n+a-N);
    lowerTail=false;
  }
  
  int na_N=n+a-N; 
  double logP=loggamma((double)(a+1))+loggamma((double)(N-a+1))+loggamma((double)(n+1))+
    loggamma((double)(N-n+1))-loggamma((double)(N+1))-loggamma((double)(a-xmin+1))-
    loggamma((double)(n-xmin+1))-loggamma(xmin-na_N+1);
  
  if (xmin!=0) {
    logP-=loggamma((double)(xmin+1));
  }
  
  // Use normal approximation if can't do it
  if (! R_FINITE(logP)){
    double p=PeizerHypergeometric(x,a,n,N);
    return lowerTail?p:1.0-p;
  }
  
  double term=1.0;
  double sum=1.0;
  // These are the terms of F[-a,-n;N-n-a+1;a], where F is the Gaussian
  // hypergeometric function -- i.e. coefficients of x^i in the expansion.
  for (int k=xmin;k<x;k++) { 
    term*=((double)(a-k)*(double)(n-k))/((double)(k+1)*(double)(k+1-na_N));
    sum+=term;
  }
  // Use normal aapproximation if can't do it
  if (! R_FINITE(sum)){
    double p=PeizerHypergeometric(x,a,n,N);
    return lowerTail?p:1.0-p;
  }
  
  logP+=log(sum);
  if (logP<-MAXEXP) {
    return lowerTail?0.0:1.0;
  }
  else {
    return lowerTail?exp(logP):1.0-exp(logP);
  }
}


DISTS_API void ughyperR(
    int *kp,
    double *ap,
    double *np,
    double *Np,
    int *Mp,
    double *valuep
)
{
  hyperType variety;
  int M=*Mp;
  int i;
  
  for (i=0;i<M;i++) {
    variety=typeHyper(ap[i],np[i],Np[i]);
    if (! checkHyperArgument(kp[i],ap[i],np[i],Np[i],variety))
      valuep[i]=NA_REAL;
    else if (variety==classic)
      valuep[i]=qhypergeometric(kp[i],(int)ap[i],(int)np[i],(int)Np[i]);
    else
      valuep[i]=qgenhypergeometric(kp[i],ap[i],np[i],Np[i],variety);
  }
  
  
}


double qhypergeometric(
    int x,  // Number of marked items in sample
    int a,     // Sample size
    int n,  // Total number of marked items
    int N   // Total number of items
)
{
  return 1.0-phypergeometric(x,a,n,N);
}

DISTS_API void qghyperR(
    double *pp,
    double *ap,
    double *np,
    double *Np,
    int *Mp,
    double *valuep
)
{
  hyperType variety;
  int M=*Mp;
  int i;
  
  for (i=0;i<M;i++) {
    variety=typeHyper(ap[i],np[i],Np[i]);
    if (variety==classic)
      valuep[i]=xhypergeometric(pp[i],(int)ap[i],(int)np[i],(int)Np[i]);
    else if (variety!=noType)
      valuep[i]=xgenhypergeometric(pp[i],ap[i],np[i],Np[i],variety);
    else valuep[i]=NA_REAL;
  }
  
  
}

// returns smallest x such that p<=Pr(X<=x|a,n,N) 
int xhypergeometric(
    double p,  // cumulative probability
    int a,     // Sample size
    int n,  // Total number of marked items
    int N   // Total number of items
)
{
  double T=qchisq(1.0-p,1,true,false);
  double z=(T*(p*(1.0-p)*(double)(a*(N-a))))/(double)(N-1);
  int x=(int)floor(0.5+p*(double)a+z*z);
  
  int minX=maxm(0,n+a-N); 
  int maxX=minm(n,a);
  x=maxm(x,minX);
  x=minm(x,maxX);
  
  if (0>p || p>1.)
    error("\nProbability must be in the 0 to 1 range");
  
  bool larger=(p<=phypergeometric(x,a,n,N));
  while (larger) {
    if (x equals minX) {
      return x;
    }
    larger=(p<=phypergeometric(--x,a,n,N));
    if (! larger) {
      return ++x;
    }
  }
  while (! larger){
    larger=(p<=phypergeometric(++x,a,n,N));
    if (larger) {
      return x;
    }
  }
  
  return 0;
}


DISTS_API void rghyperR(
    double *ap,
    double *np,
    double *Np,
    int *Mp,
    int *Kp,
    double *valuep
)
{
  hyperType variety;
  
  int M=*Mp;
  int K=*Kp;
  int D;
  int j;
  int k;
  int loc;
  int cloc;
  double *tArray;
  
  if (K==1) {
    variety=typeHyper(*ap,*np,*Np);
    if (variety==classic)
      rhypergeometric(valuep,M,(int)*ap,(int)*np,(int)*Np);
    else if (variety!=noType)
      rgenhypergeometric(valuep,M,*ap,*np,*Np,variety);
    else 
      error("\nParameters are for no recognized type");
  }
  else { // Allow for random values for each element of nu and lambda
    D=(M/K)+((M%K)?1:0);
    tArray=(double *)S_alloc((long)D,sizeof(double));
    loc=0;
    for (j=0;j<K;j++) {
      variety=typeHyper(ap[j],np[j],Np[j]);
      if (variety==classic)
        rhypergeometric(tArray,D,(int)ap[j],(int)np[j],(int)Np[j]);
      else if (variety!=noType)
        rgenhypergeometric(tArray,D,ap[j],np[j],Np[j],variety);
      else 
        error("\nParameters are for no recognized type");
      for (k=0;k<D;k++) {
        cloc=loc+k*K;
        if (cloc<M)
          valuep[cloc]=tArray[k];
        else break;
      }
      loc++;
    }
  }
  
}



/*
Random samples from hypergeometric
*/
void rhypergeometric(
    double* randArray,
    int n,  // number of samples
    int a,     // Sample size
    int m,  // Total number of marked items
    int N   // Total number of items
)
{
  GetRNGstate();
  
  for (int i=0;i<n;i++){ 
    randArray[i]=(double)xhypergeometric(unif_rand(),a,m,N);
  }
  PutRNGstate();
}

DISTS_API void sghyperR(
    double *ap,
    double *mp,
    double *Np,
    int *Mp,
    double *meanp,
    double *medianp,
    double *modep,
    double *variancep,
    double *thirdp,
    double *fourthp
)
{
  int M=*Mp;
  int i;
  hyperType variety;
  
  for (i=0;i<M;i++) {
    variety=typeHyper(ap[i],mp[i],Np[i]);
    sghyper(ap[i],mp[i],Np[i],meanp+i,medianp+i,modep+i,variancep+i,thirdp+i,fourthp+i,variety);
  }
}

void sghyper(
    double a,
    double m,
    double N,
    double *mean,
    double *median,
    double *mode,
    double *variance,
    double *third,
    double *fourth,
    hyperType variety
  
)
{
  bool paramSet=false;
  double n=0;
  double A=0;
  double B=0;
  double T=0;
  
  double m1=0;
  double m2=0;
  double m3=0;
  double m4=0;
  
  
  switch (variety) {
  case IIIB:
    A=minm(m,a);
    n=maxm(m,a);
    B=N-A;
    paramSet=true;
  case IIB:
    if (! paramSet) {
      A=minm(m,a);
      n=maxm(m,a);
      B=N-A;
    }
    T=A+B;
    *mode=(int)n+1;
    *median=(double)xgenhypergeometric(0.5,a,m,N,variety);
    *mean=NA_REAL;
    *variance=NA_REAL;
    *third=NA_REAL;
    *fourth=NA_REAL;
    break;
  case classic:
    n=minm(m,a);
    A=maxm(m,a);
    B=N-maxm(m,a);
    *median=(double)xhypergeometric(0.5,(int)a,(int)m,(int)N);
    paramSet=true;
  case IAi:
    if (! paramSet) {
      n=minm(m,a);
      A=maxm(m,a);
      B=N-A;
      *median=(double)xgenhypergeometric(0.5,a,m,N,variety);
      paramSet=true;
    }
  case IAii:
    if (! paramSet) {
      n=minm(m,a);
      A=maxm(m,a);
      B=N-A;
      *median=(double)xgenhypergeometric(0.5,a,m,N,variety);
      paramSet=true;
    }
  case IIA:
    if (! paramSet) {
      A=minm(m,a);
      n=maxm(m,a);
      B=N-A;
      *median=(double)xgenhypergeometric(0.5,a,m,N,variety);
      paramSet=true;
    }
  case IIIA:
    if (! paramSet) {
      A=minm(m,a);
      n=maxm(m,a);
      B=N-A;
      *median=(double)xgenhypergeometric(0.5,a,m,N,variety);
      paramSet=true;
    }
    
    T=A+B;
    
    if (n<=1.0) {
      *mean=0.0;
    }
    else {
      *mean=m1=(n*A)/T;
    }
    *mode=floor(((n+1.0)*(A+1.0))/(T+2.0));
    if (n<=2.0) {
      *variance=0.0;
    }
    else {
      m2=m1*(B*(T-n))/(T*(T-1.0));
      *variance=m2;
    }
    if (n<=3.0) {
      *third=0.0;
    }
    else {
      m3=m2*((B-A)*(T-2.0*n))/(T*(T-2.0));
      *third=m3;
    }
    if (n<=4.0) {
      *fourth=0.0;
    }
    else {
      m4=(m2/((T-2.0)*(T-3.0)))*(T*(T+1.0-6.0*n)+3.0*A*B*(n-2.0)+
        6.0*n*n+(3.0*A*B*n*(6.0-n))/(T)-(18.0*A*B*n*n)/(T*T));
      *fourth=m4;
    }
    break;
  case IB:
    if (! paramSet) {
      paramSet=true;
    }
  case IV:
    if (! paramSet) {
      paramSet=true;
    }
    
    *median=(double)xgenhypergeometric(0.5,a,m,N,variety);
    
    A=minm(m,a);
    n=maxm(m,a);
    B=N-A;
    T=A+B;
    
    if (T<=0.0) {
      *mean=NA_REAL;
    }
    else {
      *mean=m1=(n*A)/T;
    }
    *mode=floor(((n+1.0)*(A+1.0))/(T+2.0));
    if (T<=1.0) {
      *variance=NA_REAL;
    }
    else {
      m2=m1*(B*(T-n))/(T*(T-1.0));
      *variance=m2;
    }
    if (T<=3.0) {
      *third=NA_REAL;
    }
    else {
      m3=m2*((B-A)*(T-2.0*n))/(T*(T-2.0));
      *third=m3;
    }
    if (T<=4.0) {
      *fourth=NA_REAL;
    }
    else {
      m4=(m2/((T-2.0)*(T-3.0)))*(T*(T+1.0-6.0*n)+3.0*A*B*(n-2.0)+
        6.0*n*n+(3.0*A*B*n*(6.0-n))/(T)-(18.0*A*B*n*n)/(T*T));
      *fourth=m4;
    }
    break;
    
  default:
    break;
  }
  
}


// Generalized hypergeometric

double pgenhypergeometric(
    int x,  
    double a,  
    double n,  
    double N,
    hyperType variety 
)
{
  double logP=0;
  double b=0;
  double temp=0;
  double P=0;
  
  switch (variety) {
  case IAii:
    temp=a;
    a=n;
    n=temp;
    variety =IAi;
  case IAi:
    if (x equals (int)n) {
      return 1.0;
    }
    b=N-a;
    logP=loggamma(b+1.0)+loggamma(N-n+1.0)-loggamma(b-n+1.0)-loggamma(N+1.0);
    break;
  case IIIA:
    temp=a;
    a=n;
    n=temp;
    variety =IIA;
  case IIA:
    if (x equals (int)n) {
      return 1.0;
    }
    b=N-a;
    logP=loggamma(-b+n)+loggamma(-N)-loggamma(-b)-loggamma(-N+n);
    break;
  case IIIB:
    temp=a;
    a=n;
    n=temp;
    variety=IIB;
  case IIB:
    b=N-a;
    // Can't use this because n is not an integer
    //logP=loggamma(-b+n)+loggamma(-N)-loggamma(-b)-loggamma(-N+n);
    P=1.0/GaussianHypergometricFcn(-n,-a,b-n+1.0,1.0);
    break;
  case IB:
    b=N-a;
    logP=loggamma(b+1.0)+loggamma(N-n+1.0)-loggamma(b-n+1.0)-loggamma(N+1.0);
    break;
  case IV:
    b=N-a;
    logP=loggamma(b+1.0)+loggamma(N-n+1.0)-loggamma(b-n+1.0)-loggamma(N+1.0);
    break;
  default:
    break;
  }
  
  double sum=1.0;
  double Tr=1.0;
  double bn=b-n;
  
  for (int i=0;i<x;i++) {
    double r=(double)i;
    double rp=(double)(i+1);
    Tr*=((r-a)*(r-n))/(rp*(bn+rp));
    sum+=Tr;
  }
  
  if (variety equals IIB) {
    P*=sum;
    return minm(P,1.0);  // Occasional numerical error
  }
  else {
    logP+=log(sum);
    
    if (logP<-MAXEXP) {
      return 0.0;
    }
    else {
      return exp(logP);
    }
  }
  
  
}

double qgenhypergeometric(
    int x,  
    double a,  
    double n,  
    double N,
    hyperType variety
)
{
  return 1.0-pgenhypergeometric(x,a,n,N,variety);
}

double fgenhypergeometric(
    int x,  
    double a, 
    double n,  
    double N,
    hyperType variety  
)
{
  double logP=0;
  double b=0;
  double temp=0;
  double P=0;
  
  switch (variety) {
  case IAii:
    temp=a;
    a=n;
    n=temp;
    variety =IAi;
  case IAi:
    b=N-a;
    logP=loggamma(a+1.0)+loggamma(b+1.0)+loggamma(n+1.0)+loggamma(N-n+1.0);
    logP-=loggamma(x+1.0)+loggamma(a-x+1.0)+loggamma(n-x+1.0)+loggamma(b-n+x+1.0)+loggamma(N+1.0);
    break;
  case IIIA:
    temp=a;
    a=n;
    n=temp;
    variety =IIA;
  case IIA:
    b=N-a;
    logP=loggamma(-a+x)+loggamma(-b+n-x)+loggamma(n+1.0)+loggamma(-N);
    logP-=loggamma(x+1.0)+loggamma(-a)+loggamma(n-x+1.0)+loggamma(-b)+loggamma(-N+n);
    break;
  case IIIB:
    temp=a;
    a=n;
    n=temp;
    variety=IIB;
  case IIB:
    b=N-a;
    // (1) -1<b<0. 
    // (2) n is not itegral.
    // (3) must sum terms, cannot use loggamma 
    P=1.0/GaussianHypergometricFcn(-n,-a,b-n+1.0,1.0);
    {
      double Tr=1.0;
      double bn=b-n;
      
      for (int i=0;i<x;i++) {
        double r=(double)i;
        double rp=(double)(i+1);
        Tr*=((r-a)*(r-n))/(rp*(bn+rp));
      }
      P*=Tr;
    }
    
    break;
  case IB:
    b=N-a;
    // Assuming b>0 always
    logP=loggamma(a+1.0)+loggamma(b+1.0)+loggamma(n+1.0)+loggamma(N-n+1.0);
    logP-=loggamma(x+1.0)+loggamma(a-x+1.0)+loggamma(n-x+1.0)+loggamma(b-n+x+1.0)+loggamma(N+1.0);
    break;
  case IV:
    b=N-a;
    logP=loggamma(-a+x)+loggamma(b+1.0)+loggamma(-n+x)+loggamma(N-n+1.0);
    logP-=loggamma(x+1.0)+loggamma(-a)+loggamma(b-n+x+1.0)+loggamma(-n)+loggamma(N+1.0);
    break;
  default:
    break;
  }
  
  
  if (variety equals IIB) {
    return P;
  }
  else {
    if (logP<-MAXEXP) {
      return 0.0;
    }
    else {
      return exp(logP);
    }
  }
}


// returns smallest x such that p<=Pr(X<=x|a,n,N) 
int xgenhypergeometric(
    double p,  // cumulative probability
    double a,  
    double m,  
    double N,
    hyperType variety
)
{
  double b=N-a;
  double n=m;
  
  double m1=(n*a)/N;
  double m2=(m1*(b*(a+b-n)))/(N*(N-1.0));
  
  if (0>p || p>1)
    error("\nProbability must be in the 0 to 1 range");
  
  int x=(int)(0.5+m1+sqrt(m2)*qnorm(p,0,1,true,false));
  x=maxm(0,x);
  
  
  bool larger=(p<=pgenhypergeometric(x,a,m,N,variety));
  while (larger) {
    if (x equals 0) {
      return x;
    }
    larger=(p<=pgenhypergeometric(--x,a,m,N,variety));
    if (! larger) {
      return ++x;
    }
  }
  while (! larger){
    larger=(p<=pgenhypergeometric(++x,a,m,N,variety));
    if (larger) {
      return x;
    }
  }
  
  return 0;
}


/*
Random samples from beta negative binomial
*/
void rgenhypergeometric(
    double* randArray,
    int K,  
    double a,  
    double n,  
    double N,
    hyperType variety
)
{
  GetRNGstate();
  
  for (int i=0;i<K;i++){  
    randArray[i]=(double)xgenhypergeometric(unif_rand(),a,n,N,variety);
  }
  
  PutRNGstate();
}



/*
Integrages a function using Romberg iteration
*/

double Integral(
    double lowX, // Assumed to be less than highX
    double highX,
    double (*function)(double x),
    double Tol
)
{
  
  double h=highX-lowX; // range of integration
  
  /*const double Tol=3e-8;*/ // result to this accuracy -- can't really to better
  
  const int maxiterate=16;  // Maximum interates allowed -- fcn evaluated at most 2^16 pts
  double A[maxiterate][maxiterate]; // the triangular array of Romberg iterates
  // stop when two diagonal values agree to MRatioTol
  
  double delta=h;   // Initial step
  double value=0.0;   // Converged value
  
  A[0][0]=(h/2.0)*((*function)(lowX)+(*function)(highX)); 
  
  double twoPower=1.0;
  int k=0;
  int n=1;
  repeat
    k++;
  delta*=0.5;  // Evaluate (*funciton)() at half the previous interval
  if (k>1) {
    n*=2;   // Number of new evaluations
  }
  
  twoPower*=2.0;
  double sum=0.0;
  double z=highX-delta; // Start with this value
  int m=n;
  while (m--) {
    double value=(*function)(z);
    
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
  if (fabs((value - A[k-1][0])/value)<Tol) {
    break;
  }
  until(k>=maxiterate-1);
  
  
  
  return value;
}


/*
Root of f(root)=0 by Newton iteration
When inverting a distribution, one wants to minimize (p(x)-p)^2, so the Newton
step is h=(p(x)-p)^2/[2(p(x)-p)p'(x)] = (p(x)-p)/[2p'(x)]
Don't use to find the max or min of a function -- the derivative may go to zero
*/

double NewtonRoot(
    double guess,
    bool useLog, // When true, will iterate on z=log(x), keeping x positive
    double (*function)(double x),
    double (*derivative)(double x),
    double TOLN
)
{
  /* const double TOLN=3e-8;*/
  const int MAXITERN=100;
  
  double x=guess;
  double z=(useLog)?log(x):x;
  bool more;
  double h;
  double ho=DBL_MAX;
  double scale=1.0;
  int m=0;
  repeat
    double fcn=(*function)(x);
  double deriv=(*derivative)(x);
  if (useLog) {
    h=scale*0.5*fcn/(x*deriv+DBL_EPSILON*fabs(fcn)); // no small divisors allowed
  }              // This is likely of no use
  else {
    h=scale*0.5*fcn/(deriv+DBL_EPSILON*fabs(fcn));
  }
  
  if (! R_FINITE(h)) {
    error("\nInfinite value in NewtonRoot()");
    return x;
  }
  z-=h;
  more=(fabs(h/z)>TOLN);
  // Iterates increasing
  if (ho<=fabs(h)) {
    scale/=2.0;
    z+=h;   // Restore value
    more=true;
    continue; // retry
  }
  scale=(scale<1.0)?scale*2.0:scale;
  ho=fabs(h);
  
  x=(useLog)?exp(z):z;
  until(m++>MAXITERN || ! more);
  
  if (m>MAXITERN) {
    error("\nIteration limit exceeded in NewtonRoot()");
  }
  
  return x;
}


/*
The Gaussian hypergeometric function, usually denoted as
F[a,b;c;x]
*/
double GaussianHypergometricFcn(
    double a,
    double b,
    double c,
    double x
)
{
  int const MAXITERATES=100;
  
  if (c<0.0 && floor(c) equals c) 
    return NA_REAL;
  
  double sum=0.0;
  double term=1.0;
  int j=1;
  double dj;
  double djm1;
  repeat
    dj=(double)j;
  djm1=dj-1.0;
  sum+=term;
  term*=((a+djm1)*(b+djm1))/(c+djm1)*(x/dj);
  j++;
  until(sum+term equals sum || j>MAXITERATES);
  
  return sum;
}

// The following was moved to here so that ziggurat can use it
// Marsaglia's multiply with cary

static const double RANDCONST=2.32830643654e-10;

ULONG zSeed=362436069, wSeed=521288629;
#define zNew ((zSeed=36969*(zSeed&65535)+(zSeed>>16))<<16)
#define wNew ((wSeed=18000*(wSeed&65535)+(wSeed>>16))&65535)
#define IUNIFORM (zNew+wNew)
#define UNIFORM ((zNew+wNew)*RANDCONST)
#define setseed(A,B) zSeed=A;wSeed=B;
#define getseed(A,B) A=zSeed;B=wSeed;

/* The ziggurat method for RNOR and REXP

Modified very slightly to remove some of the awkwardness of the original C code. In partcular, 
float has been replaced by double, some implicit coercions made explicit, and the 
initialization of jsr moved into zigset to guarantee same sequence from same seed.
REW Mar 01
REW (Mar 2001) cleaned up the code. Made RNOR and REXP inline instead of defines, etc.

Combine the code below with the main program in which you want
normal or exponential variates. Then use of RNOR in any expression
will provide a standard normal variate with mean zero, variance 1,
while use of REXP in any expression will provide an exponential variate
with density exp(-x),x>0.
Before using RNOR or REXP in your main, insert a command such as
zigset(86947731);
with your own choice of seed value>0, rather than 86947731.
(If you do not invoke zigset(...) you will get all zeros for RNOR and REXP.)
For details of the method, see Marsaglia and Tsang, "The ziggurat method
for generating random variables", Journ. Statistical Software.
*/


double nfix(void);
double efix(void);

static bool ziggInitialized=false; // Makes sure zigg is initialized at least once

static ULONG jz,
jsr; // moved initialization
static ULONG jcong;
#define SHR3 (jz=jsr, jsr^=(jsr<<13), jsr^=(jsr>>17), jsr^=(jsr<<5),jz+jsr)
#define CONG (jcong=69069*jcong+1234567)
#define KISS ((IUNIFORM^CONG)+SHR3)
/*Switched to KISS -- see Leong, et.al. (2005) A comment on the implementtion of the Ziggurat Method 
Journal of statistical software 12-7, 1-4*/
/*#define UNI (.5 + (signed)(SHR3) * .2328306e-9) 
#define IUNI SHR3*/
#define UNI (.5 + (signed)(KISS) * .2328306e-9)
#define IUNI KISS

static long hz;

static int  iz;   // was ULONG
static long  kn[128]; // was ULONG
static ULONG ke[256];
static double wn[128],
                fn[128], 
                  we[256],
                    fe[256];

inline double RNOR(void) {
  double xx;
  hz=SHR3; 
  iz=hz&127; 
  xx=(abs(hz)<kn[iz])? double(hz*wn[iz]) : nfix(); // Original had sign conflict
  return xx;           // for abs(hz)<kn[]
}

inline double REXP(void)
{
  double xx;
  jz=SHR3; 
  iz=jz&255; 
  xx=(jz<ke[iz])? double(jz*we[iz]) : efix();
  return xx;
}

#ifdef NOTUSED
// There is noting wrong with these, as modified here, but there is also no advantage.
#define RNOR (hz=SHR3, iz=hz&127, (abs(hz)<kn[iz])? double(hz*wn[iz]) : nfix())
#define REXP (jz=SHR3, iz=jz&255, ( jz <ke[iz])? double(jz*we[iz]) : efix())

#endif

/* nfix() generates variates from the residue when rejection in RNOR occurs. */
double nfix(void) { 
  const double r = 3.442619855899;  /* The starting of the right tail */ 
double x;
double y;

repeat  
  x=hz*wn[iz];
if(iz==0){ /* iz==0, handle the base strip */
repeat 
  x=-log(UNI)/r;    
  y=-log(UNI);   
  solongas(y+y<x*x);
  return (hz>0)? r+x : -r-x;  
}

/* iz>0, handle the wedges of other strips */  
if( fn[iz]+UNI*(fn[iz-1]-fn[iz]) < exp(-.5*x*x) ) 
  return x;

/* start all over */  
hz=SHR3;  
iz=hz&127;  
if(abs((signed)hz)<(signed)kn[iz]) 
  return (hz*wn[iz]); 
forever
  
}



/* efix() generates variates from the residue when rejection in REXP occurs. */

double efix(void)
  
{ 
  double x;
  repeat  
    if(iz==0) 
      return (7.69711-log(UNI));  /* iz==0 */
x=jz*we[iz];    
if( fe[iz]+UNI*(fe[iz-1]-fe[iz]) < exp(-x) ) 
  return (x);
/* Start all over again */  
jz=SHR3;  
iz=(jz&255);  
if(jz<ke[iz]) 
  return (jz*we[iz]);
forever
  
}



/*--------This procedure sets the seed and creates the tables------*/

void zigset(ULONG jsrseed) {  
  const double m1 = 2147483648.0; 
  const double m2 = 4294967296.0;
  double dn=3.442619855899;
  double  tn=dn;
  double  vn=9.91256303526217e-3; 
  double  q; 
  double de=7.697117470131487; 
  double  te=de; 
  double  ve=3.949659822581572e-3;
  int i; 
  
  jsr=123456789; // moved to here to allow seed to control generation
  jsr^=jsrseed;
  setseed(jsrseed,jsrseed) /* Added to support Leong et.al.*/
jcong=jsrseed;
  
  /* Set up tables for RNOR */  
  q=vn/exp(-.5*dn*dn);
  kn[0]=(long)((dn/q)*m1);  
  kn[1]=0;   
  wn[0]=q/m1;  
  wn[127]=dn/m1;
  fn[0]=1.0;  
  fn[127]=exp(-.5*dn*dn); 
  
  for(i=126;i>=1;i--) { 
    dn=sqrt(-2.*log(vn/dn+exp(-.5*dn*dn))); 
    kn[i+1]=(ULONG)((dn/tn)*m1);   
    tn=dn; 
    fn[i]=exp(-.5*dn*dn); 
    wn[i]=dn/m1; 
  }
  
  
  /* Set up tables for REXP */  
  q = ve/exp(-de);
  
  
  ke[0]=(ULONG)((de/q) * m2);  
  ke[1]=0;
  
  we[0]=q/m2;  
  we[255]=de/m2;
  
  fe[0]=1.0;  
  fe[255]=exp(-de);  
  
  for(i=254;i>=1;i--) { 
    de=-log(ve/de+exp(-de)); 
    ke[i+1]= (ULONG)((de/te)*m2);   
    te=de; 
    fe[i]=exp(-de); 
    we[i]=de/m2; 
  }
  
}





/* The following was moved above ziggurat
// Marsaglia's multiply with cary

static const double RANDCONST=2.32830643654e-10;

ULONG zSeed=362436069, wSeed=521288629;
#define zNew ((zSeed=36969*(zSeed&65535)+(zSeed>>16))<<16)
#define wNew ((wSeed=18000*(wSeed&65535)+(wSeed>>16))&65535)
#define IUNIFORM (zNew+wNew)
#define UNIFORM ((zNew+wNew)*RANDCONST)
#define setseed(A,B) zSeed=A;wSeed=B;
#define getseed(A,B) A=zSeed;B=wSeed;
*/

static int nSeed=1020;
static ULONG Q[1020]; // using Q[endQ] to hold variable.
static int endQ=nSeed-1;
static bool QInitialized=false; // Makes sure Q is initialized at least once

#ifdef CANTUSE
// At present time, R allows only a 625 seed array 
static double uval;

DISTS_API int *user_unif_nseed(void)
{
  return &nSeed;
}

DISTS_API int *user_unif_seedloc(void)
{
  return (int *)Q;
}
// R looks for this even when it doesn't use it
DISTS_API double *user_unif_rand(void)
{
  GetRNGstate();
  
  setseed(Q[0],Q[1]);
  uval=UNIFORM;
  getseed(Q[0],Q[1]);
  
  PutRNGstate();
  
  return (double *)&uval;
  
}

void user_unif_init() {}
#endif

void QInit(
    ULONG seed
)
{
  int i;
  ULONG mask = 0xFFFF; 
  setseed(seed&mask,seed>>16);
  Q[endQ]=362436;
  for (i=0; i<endQ; i++)
    Q[i]=IUNIFORM;
  
}


ULONG MWC1019(void){
  //ULONG t;
  int i
  
  ULONG t = 147669672L*Q[i] + Q[endQ]; 
  Q[endQ] = (t >> 16L);
  if(i>0) 
    return(Q[i--] = t);
  i = endQ-1; 
  return(Q[0] = t);
}

DISTS_API void MWC1019R (
    double *randomVector,
    int *Np,
    bool *initializep,
    ULONG *seedp
)
{
  int N=*Np;
  int i;
  ULONG seed=*seedp;
  
  if (*initializep) {
    QInit(seed);
    QInitialized=TRUE; 
  }
  else if (QInitialized==FALSE) { // To insure initialization
    QInit((ULONG)556677);
    QInitialized=TRUE; 
  }
  
  
  for (i=0;i<N;i++)
    randomVector[i]=(double)MWC1019()*RANDCONST;
  
  
  
}

/*
It will provide random 32-bit integers at the rate of 300 million per
second
(on a 850MHz PC).

It requires that you seed Q[0],Q[1],...Q[1018] with 32-bit random integers,

before calling MWC1019( ) in your main. You might use a good RNG such as
KISS
to fill the Q array.

The period of MWC1029 exceeds 10^9824, making it billions and billions ...
and billions times as long as the highly touted longest-period RNG,
the Mersenne twister. It is also several times as fast and takes a few
lines rather than
several pages of code. (This is not to say that the Mersenne twister is
not
a good RNG; it is. I just do not equate complexity of code with
randomness. It is the
complexity of the underlying randomness that counts.)

As for randomness, it passes all tests in The Diehard Battery of Tests of
Randomness
http://stat.fsu.edu/pub/diehard
as well as three new tough tests I have developed with the apparent
property
that a RNG that passes tuftsts.c will pass all the tests in Diehard.

MWC1019 has the property that every possible sequence of 1018 successive
32-bit integers will appear somewhere in the full period, for those
concerned
with the "equi-distribution" in dimensions 2,3,...1016,1017,1018.

I welcome comments on timings or otherwise.

George Marsaglia */

DISTS_API void ziggR(
    double *randomVector,
    int *Np,
    bool *type,
    bool *initilizep,
    ULONG *seedp
)
{
  int N=*Np;
  int i;
  
  if (*initilizep) {
    zigset(*seedp);
    ziggInitialized=true;
  }
  else if (ziggInitialized==false) { // To always insure initialization
    zigset((ULONG)556677);
    ziggInitialized=true;
  }
  
  if (*type==true) {
    for (i=0;i<N;i++) {
      randomVector[i]=RNOR();
    }
  }
  else {
    for (i=0;i<N;i++)
      randomVector[i]=REXP();
  }
  
  
}

DISTS_API void normOrdR(
    double *sp,
    int *np,
    int *n2p
)
{
  
  nscor2(sp,np,n2p);
  
}

// static double correc(int, int);

void nscor2(
    double *s, 
    int *n, 
    int *n2
)
{
  
  /* algorithm as 177.3, applied statistics, v.31, 161-165, 1982.
  
  calculates approximate expected values of normal order statistics.
  claimed accuracy is 0.0001, though usually accurate to 5-6 dec.
  
  *** N.B. This routine was NOT in double precision All constants were f ***
  
  Arguments:
  
  s(n2) = output, the first n2 expected values.
  n  = input, the sample size.
  n2  = input, the number of order statistics required; must
  be <= n/2.
  
  ier removed REW Mar 2001
  
  ier = output, error indicator
  = 0 if no error detected
  = 1 if n <= 1.
  = 2 if n > 2000, in which case the order statistics
  are still calculated, but may be inaccurate.
  = 3 if n2 > n/2 (n.b. this differs from the
  published algorithm which returns an error
  if n2 is not equal to n/2.)
  
  Calls qnorm() [from R] which is an improvement of
  ppnd = applied statistics algorithm 111.
  An alternative is ppnd7 in algorithm AS 241.
  */
  
  /* Initialized data */
  
  const double
  eps[4] = { .419885,.450536, .456936, .468488 },
    dl1[4] = { .112063,.12177, .239299, .215159 },
    dl2[4] = { .080122,.111348,-.211867,-.115049 },
    gam[4] = { .474798,.469051, .208597, .259784 },
    lam[4] = { .282765,.304856,.407708,.414093 };
  const double bb = -.283833;
  const double d = -.106136;
  const double b1 = .5641896;
  
  
  /* Local variables */
  int i, k;
  double e1, e2, ai, an;
  
  /* input parameter checks. */
  
  if (*n2 > *n / 2) {
    error("\nn2>n");
  }
  if (*n <= 1) {
    error("\nn<=1");
  }
  if (*n > 2000) {
    warning("\nValues may be inaccurate because of the size of N");
  }
  
  s[0] = b1;
  if (*n == 2) {
    return;
  }
  
  /* calculate normal tail areas for first 3 order statistics. */
  
  an = (double) (*n);
  k = 3;
  if (*n2 < k)
    k = *n2;
  /* k := min(3, *n2) */
  for (i = 0; i < k; ++i) {
    ai = (double) i+1;
    e1 = (ai - eps[i]) / (an + gam[i]);
    e2 = pow((double) e1, (double) lam[i]);
    s[i] = e1 + e2 * (dl1[i] + e2 * dl2[i]) / an - correc(i+1, *n);
  }
  if (*n2 > k) {
    
    /* calculate normal areas for other cases. */
    
    for (i = 4-1; i < *n2; ++i) {
      ai = (double) i+1;
      e1 = (ai - eps[3]) / (an + gam[3]);
      e2 = pow((double) e1, (double) lam[3] + bb / (ai + d));
      s[i] = e1 + e2 * (dl1[3] + e2 * dl2[3]) / an - correc(i+1, *n);
    }
  }
  /* convert tail areas to normal deviates. */
  
  for (i = 0; i < *n2; ++i)
    s[i] = - qnorm(s[i], 0., 1., 1, 0);
  
  return;
} /* nscor2 */
  
  
  static double correc(int i, int n)
  {
    /* calculates correction for tail area of the i-th largest of n
    order statistics. */
    
    const double
    c1[7] = { 9.5,28.7,1.9,0.,-7.,-6.2,-1.6 },
      c2[7] = { -6195.,-9569.,-6728.,-17614.,-8278.,-3570., 1075. },
      c3[7] = { 93380.,175160.,410400.,2157600.,2.376e6,
                2.065e6,2.065e6 };
    const double mic = 1e-6;
    const double c14 = 1.9e-5;
    
    double an;
    
    if (i * n == 4)  return c14;
    if (i < 1 || i > 7)  return 0;
    if (i != 4 && n > 20) return 0;
    if (i == 4 && n > 40) return 0;
    /* else : */
    an = (double) n;
    an = 1.0 / (an * an);
    i--;
    return((c1[i] + an * (c2[i] + an * c3[i])) * mic);
  } /* correc */
    