#include "wheeler.h"
#include <math.h>
#include <float.h>
//Hornik replaced <new.h> with <new> and then inserted 
// std::set_new_handler(freeStoreException); on line 2283 (March 2008)
#include <new>  
#include <R.h>
#include <Rmath.h>

#include "dists.h"
#include "datatabs.h"

// SuppDists by Robert E. Wheeler, March 2001

bool DllMain(void)					 
{

    return true;
}


static const double LOG10 = 2.3025850929940456840179915;
static const double MAXEXP = LOG10 * DBL_MAX_10_EXP;	// Maximum argument for exp()
static const double TWOPI = 2 * PI;
static const double SQRT2 = 1.414135623730950488;
static const double LNGAMMAHALF = 1.144729885849400174143427/2;	// log of gamma(1/2) = log(sqrt(PI))
static const double LOG2 = 0.6931471805599453094172321;
static const double TOLNEWTON = 3e-8;
static const double LOGSQRT2PI = 0.9189385332046727417803296;
//static const double NA=-1e-12;







/*
	random normal deviates
*/	

void rgauss(
	double* normArray,
	int n,
	double mean,
	double sd
)
{ 
	int i;

	GetRNGstate();
	for (i=0;i<n;i++)
		normArray[i]=rnorm(mean,sd);
	PutRNGstate();


}

// Random chi squre
void	rdchisq(
	double *tArray,
	int n,
	int df
)
{
	int i;
	GetRNGstate();
	for (i=0;i<n;i++)
		tArray[i]=rchisq((double)df);
	PutRNGstate();
}

/*
	Derivitve of chisquared density
*/
double fpchisq(
	double x,
	int df
)
{
	double nu=(double)df/2.0;
	return dchisq(x,df,false)*((nu-1.0)/x-0.5);
}

/******************************************************************************
  	Inverse Gaussian -- Wald
*/


	/* Density function for R */
DISTS_API void dinvGaussR(
	double *xp,
	double *nup,
	double *lambdap,
	int *Np,
	double *valuep
)
{
	int N=*Np;
	int i;

	for (i=0;i<N;i++) {
		valuep[i]=finvGauss(xp[i],nup[i],lambdap[i]);
	}
	
}


   	// Density function	-- assumes x>0
 double finvGauss(
	double x,
	double mu,
	double lambda
)
{
	
	double delta=x/mu-1.0;
	double ratio=lambda/x;
	
	if (x<=0 || mu<=0 || lambda<=0)
		return NA_REAL;

	return sqrt(ratio/(TWOPI*x*x))*exp(-0.5*ratio*delta*delta);
}



	// Probability function for R
DISTS_API void pinvGaussR(
	double *xp,
	double *nup,
	double *lambdap,
	int *Np,
	double *valuep
)
{
	int N=*Np;
	int i;

	for (i=0;i<N;i++)
		valuep[i]=pinvGauss(xp[i],nup[i],lambdap[i]);
	
}



double pinvGauss(
	double x,
	double mu,
	double lambda
)
{

	double a=sqrt(lambda/x);
	double b=x/mu;
	double p1=pnorm(a*(b-1.0),0,1,true,false);
	double p2=pnorm(-a*(b+1.0),0,1,true,false);

	if (x<=0 || mu<=0 || lambda<=0)
		return NA_REAL;

	if (p2 equals 0.0) {
		return p1;
	}
	else {
		double c=2.0*lambda/mu;
		if (c>=MAXEXP) 
			return NA_REAL;

		return p1+exp(c)*p2;
	}
}

	// Upper tail Probability function for R
DISTS_API void uinvGaussR(
	double *xp,
	double *nup,
	double *lambdap,
	int *Np,
	double *valuep
)
{
	int N=*Np;
	int i;

	for (i=0;i<N;i++)
		valuep[i]=qinvGauss(xp[i],nup[i],lambdap[i]);
	
}

 double qinvGauss(
	double x,
	double mu,
	double lambda
)
{

	double a=sqrt(lambda/x);
	double b=x/mu;
	double q=1.0-pnorm(a*(b-1.0),0,1,true,false);
	double p=pnorm(-a*(b+1.0),0,1,true,false);

	if (x<=0 || mu<=0 || lambda<=0)
		return NA_REAL;

	if (p equals 0.0) {
		return q;
	}
	else {
		double c=2.0*lambda/mu;
		if (c>=MAXEXP) 
			return NA_REAL;

		return q-exp(c)*p;
	}
}


 	// Quantile function for R
DISTS_API void qinvGaussR(
	double *pp,
	double *nup,
	double *lambdap,
	int *Np,
	double *valuep
)
{
	int N=*Np;
	int i;

	for (i=0;i<N;i++)
		valuep[i]=xinvGauss(pp[i],nup[i],lambdap[i]);
	
}


 
static double gMu;
static double gLambda;
static double gP;

static double dinvGaussP(double x){return -gP+pinvGauss(x,gMu,gLambda);}
static double finvGaussP(double x){return finvGauss(x,gMu,gLambda);}

/*
	Using Whitmore and Yalovsky for an initial guess. Works well for 
		large t=lambda/mu > 2 perhaps
	Whitmore, G.A. and Yalovsky, M. (1978). A normalizing logarithmic
		transformation for inverse Gaussian random variables, 
		Technometrics 20-2, 207-208
	For small t, with x<0.5 mu, use gamma approx to 1/x -- alpha=1/2 and beta =2/lambda and 1-p
		When x>0.5mu, approx x with gamma for p, and exponentiate -- don't know why this works.
*/
 double xinvGauss(
	double p,
	double mu,
	double lambda
)
{ 	
	gMu=mu;
	gLambda=lambda;
	gP=p;

	if (0>p || p>1 || mu<=0 || lambda<=0)
		return NA_REAL;

//	double mode=-1.5*(mu*mu)/lambda+mu*sqrt(1.0+2.25*(mu*mu)/(lambda*lambda));
	double z;
	if (lambda/mu>2.0) {
		z=(qnorm(p,0,1,true,false)-0.5*sqrt(mu/lambda))/sqrt(lambda/mu);
		z=mu*exp(z);
	}
	else {
		z=lambda/(qgamma(1.0-p,0.5,1.0,true,false)*2.0);
		if (z>mu/2.0) {		// too large for the gamma approx
			z=mu*exp(qgamma(p,0.5,1.0,true,false)*0.1);  // this seems to work for the upper tail ???
		}
	}
	z=NewtonRoot(z,true,dinvGaussP,finvGaussP,3e-8);

	return z;
}



 	// Random function for R
DISTS_API void rinvGaussR(
	double *nup,
	double *lambdap,
	int *Np,
	int *Mp,	// length of nup and lambdap
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
		rinvGauss(valuep,N,*nup,*lambdap);
	else { // Allow for random values for each element of nu and lambda
		D=(N/M)+((N%M)?1:0);
		tArray=(double *)S_alloc((long)D,sizeof(double));
		loc=0;
		for (j=0;j<M;j++) {
			rinvGauss(tArray,D,nup[j],lambdap[j]);
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

   
/*
	random inverse Gaussian values
	Follows Mitchael,J.R., Schucany, W.R. and Haas, R.W. (1976). Generating
	random roots from variates using transformations with multiple roots.
	American Statistician. 30-2. 88-91.
*/	

 void rinvGauss(
	double* normArray,
	int n,
	double mu,
	double lambda
)
{ 
	double b=0.5*mu/lambda;
	double a=mu*b;
	double c=4.0*mu*lambda;
	double d=mu*mu;

	rgauss(normArray,n,0,1);
	GetRNGstate();
	for (int i=0;i<n;i++) {
		if (mu<=0 || lambda<=0) {
			normArray[i]=NA_REAL;
		}
		else {
			double u=unif_rand();
			double v=normArray[i]*normArray[i];	// Chi-square with 1 df
			double x=mu+a*v-b*sqrt(c*v+d*v*v);	// Smallest root
			normArray[i]=(u<(mu/(mu+x)))?x:d/x;	// Pick x with prob mu/(mu+x), else d/x;
			if (normArray[i]<0.0){
				v=x;
			}
		}
	}
	PutRNGstate();
}	


// UTILITIES *******************************************************************


/* 
	Natural logarithm of the gamma function.
	     CACM 291 due to M.C. Pike and I.D. Hill 
	     Accurate to at least 10 places. 
*/
double loggamma(double x)	 {
	const double  T1=1.0/12.0;
	const double  T3=1.0/360.0;
	const double  T5=1.0/1260.0;
	const double  T7=1.0/1680.0;
	const double  T9=1.0/1188.0;

	double f;

	if (x equals 1.0 || x equals 2.0) {
		return 0.0;
	}

	if (x>=7.0) {
		f=0.0;
	}
	else {
		for (f=1.0;x<7.0;x+=1.0) {
			f*=x;
		}
		f=-log(f);
	}

	double z=1.0/(x*x);
	f+=(x-0.5)*log(x)-x+LOGSQRT2PI;
	double t=(T9*z-T7)*z+T5;
	f+=((t*z-T3)*z+T1)/x;

	return (f);
}


/* normdrv -- derivatives of normal */

/* normal density */
double phi0(double x)
{
   double constant=.398942280401433;

   return (constant*exp(-0.5*x*x));
}

/* third derivative z=phi0(x)*/
double phi3(double x,double z)
{
   return (z*x*(3.0-x*x));
}

/* fifth derivative z=phi0(x)*/
double phi5(double x,double z)
{
   double s;

   s=x*x;
   return (-z*x*((s-10.0)*s+15.0));
}

/* seventh derivative z=phi0(x)*/
double phi7(double x,double z)
{
   double s;

   s=x*x;
   return(z*x*(((-s+21.0)*s-105.0)*s+105.0));
}

/* Permute ***************************************************************
|  Randomly pemutes the n integers in a[] using the Fike
|  algorithm.  See Fike, "A permutation generation method"  The Computer
|  Journal, 18-1, Feb 75, 21-22.
*/

void Permute(
	int* a,
	int n
)
{
   int i;
   int j;
   int temp;

	GetRNGstate();
   for (i=1;i<n;i++) {
		j=(int)((double)(1+i)*unif_rand());  
      temp=a[j];
      a[j]=a[i];
      a[i]=temp;
   }
   PutRNGstate();
}






/*
	Derivitive of F density
*/
double   fpfdist(
	double x,
	double nu1,
	double nu2
)
{
	double n1=nu1;
	double n2=nu2;
	double a=0.5*n1-1.0;
	double b=0.5*(n1+n2);
	return df(x,nu1,nu2,false)*(a/x-b*n1/(n1+n2*x));
}







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


 double  pKruskal_Wallis(
	double H,	 // The Statistic
	int c,		 // number of treatments
	int n,		 // Total number of observations
	double U,		// Sum (1/ni), where ni is numb obs for each treatment
	bool doNormalScore	 // do normal scores
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


 double  qKruskal_Wallis(
	double H,	 // The Statistic
	int c,		 // number of treatments
	int n,		 // Total number of observations
	double U,		// Sum (1/ni), where ni is numb obs for each treatment
	bool doNormalScore	 // do normal scores
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

 double  xKruskal_Wallis(
	double P,
	int c,		 // number of treatments
	int n,		 // Total number of observations
	double U,		// Sum (1/ni), where ni is numb obs for each treatment
	bool doNormalScore	 // do normal scores
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
	double *third,	// Central moments
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
 double  fKruskal_Wallis(
	double H,	 // The Statistic
	int c,		 // number of treatments
	int n,		 // Total number of observations
	double U,		// Sum (1/ni), where ni is numb obs for each treatment
	bool doNormalScore	 // do normal scores
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
	int N,	  // number of samples
	int c,		 // number of treatments
	int n,		 // Total number of observations
	double U,		// Sum (1/ni), where ni is numb obs for each treatment
	bool doNormalScore	 // do normal scores
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

/* kendexact -- exact distribution of Kendall's tau.  See
		Hajek, J. & Sidak, Z. (1967) Theory of rank
		tests. Academic Press, NY p140.  Should not
		be used if numbers will overflow a long --
		thus good to n=12 
	Note: Although the mean must be zero, the 50th percentile
		is not be zero for even n, since  M/2 is not integral,
		where M=n(n-1)/2, the max value of k.
*/


static double kendexact(
	int N,
	int T,
	bool density   // When true, returns the prob at T
)
{
	// a is an array containing n! times the individual probabilities
	//	for j=0, ..., k=n(n-1)/2.
	// fills() updates the current row to the one for n.

	int *a=(int *)S_alloc((long)(T+1),sizeof(int));
	memset(a,0,(T+1)*sizeof(int));

	a[0]=1;	// the row for n=1 is 1,0,0,0,0 ...
	int k=1;
	for (int n=2;n<=N;n++) {
		k=(k<=T)?k:T;	// No need to calculate beyond T
		fills(a,k,n);
		k+=n;	// this will be n(n-1)/2 after n is incremented
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


 double  pkendall(
	int ni,
	double  tau
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

double  qkendall(
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

 double  xkendall(
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
	return 4.0*(double)k/(n*(n-1.0))-1.0;  // Tau
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
	int N,	  // number of samples
	int ni
)
{
	GetRNGstate();

	for (int i=0;i<N;i++){
		randArray[i]=(double)xkendall(unif_rand(),ni);
	}

	PutRNGstate();
}



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
	int N	// Number of mean squares
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
	double h=upperLimit-lowLimit;   //  range of integration

	const double MRatioTol=1e-4;  // result to this accuracy

	const int maxiterate=16;	  // Maximum interates allowed -- fcn evaluated at most 2^16 pts
	double A[maxiterate][maxiterate]; // the triangular array of Romberg iterates
									  // stop when two diagonal values agree to MRatioTol
	
	double delta=h;		 // Initial StepTheKey
	double value=0.0;		 // Converged value

	A[0][0]=(h/2.0)*(pmaxFRatioIntegrand(lowLimit,F,df,N,logC)+pmaxFRatioIntegrand(upperLimit,F,df,N,logC)); 

	double twoPower=1.0;
	int k=0;
	int n=1;
	repeat
		k++;
		delta*=0.5;	 // Evaluate pmaxFRatioIntegrand() at half the previous interval
		if (k>1) {
			n*=2;		 // Number of new ordinates
		}

		twoPower*=2.0;
		double sum=0.0;
		double z=upperLimit-delta;	// Start with this value
		int m=n;
		while (m--) {
			double value=pmaxFRatioIntegrand(z,F,df,N,logC);

			sum+=value;
			z-=2.0*delta;	 // Every other ordinate is a new one
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
	int N	// Number of mean squares
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
			(double)(N-2)*log(fabs(pchisq(x*F,df,true,false)-pchisq(x,df,true,false)));  // fabs because F can be less than 1
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
	int N	// Number of mean squares
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
	double h=upperLimit-lowLimit;   //  range of integration

	const double MRatioTol=1e-4;  // result to this accuracy

	const int maxiterate=16;	  // Maximum interates allowed -- fcn evaluated at most 2^16 pts
	double A[maxiterate][maxiterate]; // the triangular array of Romberg iterates
									  // stop when two diagonal values agree to MRatioTol
	
	double delta=h;		 // Initial StepTheKey
	double value=0.0;		 // Converged value

	A[0][0]=(h/2.0)*(fmaxFRatioIntegrand(lowLimit,F,dgf,N,logC)+fmaxFRatioIntegrand(upperLimit,F,dgf,N,logC)); 

	double twoPower=1.0;
	int k=0;
	int n=1;
	repeat
		k++;
		delta*=0.5;	 // Evaluate fmaxFRatioIntegrand() at half the previous interval
		if (k>1) {
			n*=2;		 // Number of new ordinates
		}

		twoPower*=2.0;
		double sum=0.0;
		double z=upperLimit-delta;	// Start with this value
		int m=n;
		while (m--) {
			double value=fmaxFRatioIntegrand(z,F,dgf,N,logC);

			sum+=value;
			z-=2.0*delta;	 // Every other ordinate is a new one
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

	ddf=log(ddf)/ln2;		// log base 2 of ddf
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
	x=maxm(1.000001,x);	// Johnson sometimes comes up with a value less than 1

	bool more=true;
	double h;
	double ho=1e6;
	int m=0;
	repeat
		h=(p-pmaxfratio(x,dgf,N))/fmaxfratio(x,dgf,N);
		x+=h;
		more=(fabs(h/x)>TOLNEWTON);
		if (ho<fabs(h))	{
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

// FRIEDMANS ******************************************************************

	// This global is set when an exact distribution is active.
	// It is reset when either r or n changes.
FriedmanGlobal *FriedmanCurrentGlobal=0L;

/*  Defined in DISTFCNS.H:	
struct FriedmanStrc {
   int* S;		   // Sum of squares for each type
   int	nS;
   double *qdist;  // Cumulative frequencies summed down from the upper tail.
}; 

  This is used in RANK.CPP, hence it is defined in DISTFCNS.H, and
  since it needs FriedmanStrc, that too is defined there
struct FriedmanGlobal {
	int r;
	int n;
	FriedmanStrc *theDist;
};
*/


/* Exact distribution of Friedman test.  
	The output is a FriedmanStrc with the sums of squares in S and 1-F(S) in qdist.
	The Friedman statistic s Xr=12S/nr(r+1) and S=Sum(ri-av)**2 with ri the sum
		of ranks for ith of r treatments.  There are n blocks and av=n(r+1)/2. The 
		max value of S is n^2r(r^2-1)/12. Uses Kendall-Babington Smith algorithm
	    (AMS 1939, 10, 275-287).  See Hajek & Sidak p144 for a formal statement. 

	Spearman's rho is (Xr/(r-1))-1 for n=2, or rho=-1+2S/M, where M=r(r^2-1)/3

	The Kendall-Babington Smith algorithm:
	(1) Let b be a vector of base ranks (these are the centered ranks (ri-av))
		(1-r)/2,(3-r)/2,...,(r-1)/2	    for odd r
		(1-r),(3-r),...,(r-1)			for even r (to keep the values integral)
	(2) Add all permutations of b to b to form the distribution for n=2.
	(3) There will be duplicate types which produce the same S, and differ
		only in their order (the vector a and -a are of the same type).
		Collect these types together and pool (sum) their frequencies to form the types
		frequency.
	(4) Add all permutations of b to each type from n=2 to obtain the distributions for
		n=3, collect them by types and pool the frequencies for each type.
	(5) etc.
	(6) When the desired n is reached, order the results according to S, pooling
		frequencies for types with the same S, and return the results in a FriedmanStrc

  NOTES:
	Hajek & Sidak give a theorm (p144) which shows that the order of a type is not
		important when deriving a new type by adding permutations of b, hence the
		elements may be stored in increasing absolute order, which makes it possible
		to do a preliminary check of two arrays against each other by comparing elements
		sequentially -- if they pass this check, a further check is needed to ensure that
		they are indeed of the same type.  As an example, here are two arrays of different
		type which pass the absolute order test -- they were generated by the algorithm:
			-12,-12,-6,-2, 4, 4,10,14;
			-12,-12, 6, 2,-4,-4,10,14.
		In calculatin the Spearman rho, one only needs to check S, since new types
		will not be generated from the types for n=2.
	The algorithm has been run for a number of different values of r and n -- see
		DoExactFriedman(), and stored as arrays of type FriedmanValues.  A FriedmanStrc
		is created from such arrays in GetFriedmanStatic() to avoid the need to 
		recompute them every time.
	
	Defined in DATATABS.H	  
struct FriedmanValues {
	int S;
	double qdist;
};

*/

static const int MAXTYPE=20000;

struct Ftype {
	int *R;  // Column sums of centered ranks: sorted so |R[i]|<=|R[j]| for i<j 	
	int S;			  // Sum of sqares of R[i]
	double frequency;  	// This is a limiting factor for the algorithm -- must not overflow 
						// hence it is double.  There is the possiblity that small frequencies
						// may not be properly added to larger frequencies, but I have been
						// unable to detect any errors.
};


// Finds next permutation given key following Fike (1975)
//	   The computer Journal 18-1 21-22.
// This generates r! sequences such that 0<=key[i]<=i+1 for
//		i=0...r-1
// PermuteFike maps these into permutations of the original r marks

static bool UpdateTheKey(
	int k,
	int r,
	int *key
)
{
	if (k>=r-1) {
	   return false;
	}
	else 
	if (key[k]) {
		--key[k];
		for (int i=0;i<k;i++) {
			key[i]=i+1;
		}
		return true;
	}
	else {
		return UpdateTheKey(k+1,r,key);
	}
}

// Permutes r-1 values in vec according to key 
// The value in key[i] position is swapped with the value in the i+1th position
// thus key=1,2,3,...,r-1 does nothing, while key=2,3,...,r-1,1, shifts values in
// positions 1 and 2
static void PermuteFike(
	int *vec,
	int *key,
	int r
)
{
	for (int i=0;i<r-1;i++) {
		int temp=vec[i+1];
		vec[i+1]=vec[key[i]];
		vec[key[i]]=temp;
	}
}

	// Comparison function for qsort
	// Sorts by increasing absolute value  
static int AbsIntcmpf(
	const void* ap,
	const void* bp
)
{
	if (abs(*(int*)ap) equals abs(*(int*)bp)) {
		return 0;
	}
	return (abs(*(int*)ap)<abs(*(int*)bp))?-1:1;
}

	// Comparison function for qsort
	// Sorts by increasing  value  
static int IncreasingIntcmpf(
	const void* ap,
	const void* bp
)
{
	if ((*(int*)ap) equals (*(int*)bp)) {
		return 0;
	}
	return ((*(int*)ap)<(*(int*)bp))?-1:1;
}


	// Comparison function for qsort
	// Sorts orderC pairs by increasing  value  
static int IntPairCompare(
	const void* ap,
	const void* bp
)
{
	if (((int*)ap)[1] equals ((int*)bp)[1]) {
		return 0;
	}
	return (((int*)ap)[1]<((int*)bp)[1])?-1:1;
}


	// Checks to see if the Ftype's are the same
	// Assumes that the R arrays are in increasing absolute order
	// The types are the same if the S's are the same and the R arrays are either
	//	the same or negatives of each other
static bool IsSameType(
	int r,
	int n,
	Ftype* aType,
	Ftype* workingType
)
{
	int i;

	if (aType->S!=workingType->S) {
		return false;
	}

		// No need to check further in this case, since only one set of permutations
		// are used -- the second generation types are the final generation
  	if (n equals 2) {
		return true;
	}


	for (i=0;i<r;i++) {
		if (abs(aType->R[i])!=abs(workingType->R[i])) {
			return false;
		}
	}


		// The absolute values in the two arrays are the same, but
		// must check to see that the signs match.  So will sort them
		// in increasing order and compare them.  If this doesn't give
		// agreement, will reverse the signs of one array and try again.
	int *R1=(int *)S_alloc((long)r,sizeof(int));
	int *R2=(int *)S_alloc((long)r,sizeof(int));

	for (i=0;i<r;i++) {
		R1[i]=aType->R[i];
		R2[i]=workingType->R[i];
	}

	qsort((void*)R1,(size_t)r,sizeof(int),IncreasingIntcmpf);
	qsort((void*)R2,(size_t)r,sizeof(int),IncreasingIntcmpf);

	bool sameType=true;
	for (i=0;i<r;i++) {
		if (R1[i]!=R2[i]) {
			sameType=false;
			break;
		}
	}

	if (! sameType) {
		for (i=0;i<r;i++) {
			R2[i]=-R2[i];
		}

		for (i=0;i<r;i++) {
			if (R2[i]!=R1[i]) {
				sameType=false;
				break;
			}
		}

	}


	return sameType;


}

static void DeleteFtype(
	Ftype* aType
)
{
	delete [] aType->R;
	delete aType;
}

static int InsertTypeInList(
	int nWorkingTypes,
	Ftype** workingTypes,
	int r,
	int n,
	Ftype* aType
)
{
	for (int i=0;i<nWorkingTypes;i++) {
		if (IsSameType(r,n,aType,workingTypes[i])) {
				// If an existing type, add frequency and delete aType
			workingTypes[i]->frequency+=aType->frequency;
			DeleteFtype(aType);
			return nWorkingTypes;
		}
	}
		// A new type has been found
	if (nWorkingTypes+1>=MAXTYPE) {
		DeleteFtype(aType);
		for (int i=0;i<nWorkingTypes;i++) {
			DeleteFtype(workingTypes[i]);
		}
		nWorkingTypes=0;
		error("\nInernal error in InsertTypeInList()");   // This will abort, but leave memory dirty
		return nWorkingTypes;
	}
		// Add the aType pointer to the workingTypes list
	workingTypes[nWorkingTypes++]=aType;
	return nWorkingTypes;
}



	// The base contains the centered ranks.
	// The r! permutations of the base are added to the current type array, R, to create
	//	working types.  
	// If the resulting working type is a new type it is added to the workingTypes list.
	// If the resulting working type is not new, it is deleted and the currentType frequency
	//	is added to the existing workingType freqency.

	// NOTE: the frequencies get large and are limiting factors for the algorithm
static int AddPermutations(
	int r,
	int n,
	Ftype* currentType,
	int nWorkingTypes,
	Ftype** workingTypes,
	int* base
)
{
	int i;
		// key is used by PermuteFike() to permute the temp array
	int *key=(int *)S_alloc((long)r,sizeof(int));
	for (i=1;i<r;i++) {
		key[i-1]=i;		// Puts 1...r-1 in key (key[r-1] is not used by UpdateTheKey())
	}


	int *temp=(int *)S_alloc((long)r,sizeof(int));

	repeat 
		for (int i=0;i<r;i++) {
			temp[i]=base[i];	 
		}
		PermuteFike(temp,key,r);	 // Permutes temp according to key

		Ftype* aType=new Ftype;
		aType->R=new int[r];		// This pointer is either added to the workingType list
								// or deleted by InsertTypeInList();

		int S=0;
		for (int j=0;j<r;j++) {
			int R=aType->R[j]=currentType->R[j]+temp[j];
			S+=R*R;
		}
		aType->S=S;
		aType->frequency=currentType->frequency;
			// Keep R in increasing absolute order so IsSameType will work
		qsort((void*)aType->R,(size_t)r,sizeof(int),AbsIntcmpf);
		nWorkingTypes=InsertTypeInList(nWorkingTypes,workingTypes,r,n,aType);
	solongas(UpdateTheKey(0,r,key));	// Cycles key through all possible permutations


	return nWorkingTypes;
}




static FriedmanStrc* MakeFriedmanStrc(
	int nCurrentTypes,
	Ftype** currentTypes,
	double totalFrequency	// Double because totalFrequency is (r!)^(n-1) which gets big.
)
{
	FriedmanStrc* theStrc=new FriedmanStrc;
	int i;
		// Set up orderS for sorting
	int (*orderS)[2];
	orderS=new int[nCurrentTypes][2];
	for (i=0;i<nCurrentTypes;i++) {
		orderS[i][0]=i;		  // Set the index
		orderS[i][1]=currentTypes[i]->S;
	}
		// Sort in increasing S order, carrying index along
	qsort((void*)orderS,(size_t)nCurrentTypes,2*sizeof(int),IntPairCompare);

		// Temp arrays, because after pooling nS may be smaller than nCurrenTypes
	int* ST=new int[nCurrentTypes];
	double* qdistT=new double[nCurrentTypes];

	double sum=0.0;
		// Sum from max S to min S and form qdist
	for (i=nCurrentTypes-1;i>=0;i--) {
		int j=orderS[i][0];	 // Index of ith largest S in currentTypes
		ST[i]=currentTypes[j]->S;
		sum+=(currentTypes[j]->frequency)/totalFrequency;	
		qdistT[i]=sum;
	}
 	delete [] orderS;

		// Pool types with same S
		// Since qdist is cumulative, this means skipping some S
	int j=0;
	int skipped=0;
	double qdist=qdistT[0];
	int S=ST[0];
	for (i=1;i<nCurrentTypes;i++) {
		if (S != ST[i]) {
			ST[j]=S;
			qdistT[j]=qdist;
			S=ST[i];
			qdist=qdistT[i];
			j++;
		}
		else {
			skipped++;
		}
	}
	ST[j]=S;
	qdistT[j]=qdist;
	int nS=theStrc->nS=nCurrentTypes-skipped;

		// Either assign arrays to theStrc or create new ones and copy if nS<nCurrentTypes
	if (nS equals nCurrentTypes) {
		theStrc->S=ST;
		theStrc->qdist=qdistT;
	}
	else {
		theStrc->S=new int[nS];
		theStrc->qdist=new double[nS];
		for (int i=0;i<nS;i++) {
			theStrc->S[i]=ST[i];
			theStrc->qdist[i]=qdistT[i];
		}
		delete [] ST;
		delete [] qdistT;
	}

	return theStrc;
}



	// Creates a FriedmanStrc from data in datatabs.h
	// NOTE: r==2 is always calcualted because it takes little time
	//	The same could be done for small n's

static bool GetFriedmanStatic(
	int r,
	int n,
	FriedmanStrc** theStrc
)
{
	FriedmanValues* theValues=0L;
	bool isInTables=true;

	if (n equals 2) {
		switch (r) {
			case 3:	theValues=(FriedmanValues*)FriedmanData3_2; break;
			case 4: theValues=(FriedmanValues*)FriedmanData4_2; break;
			case 5: theValues=(FriedmanValues*)FriedmanData5_2; break;
			case 6: theValues=(FriedmanValues*)FriedmanData6_2; break;
			case 7: theValues=(FriedmanValues*)FriedmanData7_2; break;
			case 8: theValues=(FriedmanValues*)FriedmanData8_2; break;
			case 9: theValues=(FriedmanValues*)FriedmanData9_2; break;
			case 10: theValues=(FriedmanValues*)FriedmanData10_2; break;
			case 11: theValues=(FriedmanValues*)FriedmanData11_2; break;
				// Tried r==12, but it took 18 hrs and then bombed??
			default:	theValues=0L;
						isInTables=false;
						break;
		}
	}
	else
	if (r equals 3) {
		switch (n) {
			case 3:	theValues=(FriedmanValues*)FriedmanData3_3; break;
			case 4:	theValues=(FriedmanValues*)FriedmanData3_4; break;
			case 5:	theValues=(FriedmanValues*)FriedmanData3_5; break;
			case 6:	theValues=(FriedmanValues*)FriedmanData3_6; break;
			case 7:	theValues=(FriedmanValues*)FriedmanData3_7; break;
			case 8:	theValues=(FriedmanValues*)FriedmanData3_8; break;
			case 9:	theValues=(FriedmanValues*)FriedmanData3_9; break;
			case 10:	theValues=(FriedmanValues*)FriedmanData3_10; break;
			case 11:	theValues=(FriedmanValues*)FriedmanData3_11; break;
			case 12:	theValues=(FriedmanValues*)FriedmanData3_12; break;
			case 13:	theValues=(FriedmanValues*)FriedmanData3_13; break;
			case 14:	theValues=(FriedmanValues*)FriedmanData3_14; break;
			case 15:	theValues=(FriedmanValues*)FriedmanData3_15; break;
			case 16:	theValues=(FriedmanValues*)FriedmanData3_16; break;
			case 17:	theValues=(FriedmanValues*)FriedmanData3_17; break;
			case 18:	theValues=(FriedmanValues*)FriedmanData3_18; break;
			case 19:	theValues=(FriedmanValues*)FriedmanData3_19; break;
			case 20:	theValues=(FriedmanValues*)FriedmanData3_20; break;
			case 21:	theValues=(FriedmanValues*)FriedmanData3_21; break;
			case 22:	theValues=(FriedmanValues*)FriedmanData3_22; break;
			case 23:	theValues=(FriedmanValues*)FriedmanData3_23; break;
			case 24:	theValues=(FriedmanValues*)FriedmanData3_24; break;
			case 25:	theValues=(FriedmanValues*)FriedmanData3_25; break;
			case 26:	theValues=(FriedmanValues*)FriedmanData3_26; break;
			case 27:	theValues=(FriedmanValues*)FriedmanData3_27; break;
			case 28:	theValues=(FriedmanValues*)FriedmanData3_28; break;
			case 29:	theValues=(FriedmanValues*)FriedmanData3_29; break;
			case 30:	theValues=(FriedmanValues*)FriedmanData3_30; break;
			default:	theValues=0L;
						isInTables=false;
						break;
		}
	}
	else 
	if (r equals 4) {
		switch (n) {
			case 3: theValues=(FriedmanValues*)FriedmanData4_3; break;
			case 4: theValues=(FriedmanValues*)FriedmanData4_4; break;
			case 5: theValues=(FriedmanValues*)FriedmanData4_5; break;
			case 6: theValues=(FriedmanValues*)FriedmanData4_6; break;
			case 7: theValues=(FriedmanValues*)FriedmanData4_7; break;
			case 8: theValues=(FriedmanValues*)FriedmanData4_8; break;
			case 9: theValues=(FriedmanValues*)FriedmanData4_9; break;
			case 10: theValues=(FriedmanValues*)FriedmanData4_10; break;
			case 11: theValues=(FriedmanValues*)FriedmanData4_11; break;
			case 12: theValues=(FriedmanValues*)FriedmanData4_12; break;
			case 13: theValues=(FriedmanValues*)FriedmanData4_13; break;
			case 14: theValues=(FriedmanValues*)FriedmanData4_14; break;
			case 15: theValues=(FriedmanValues*)FriedmanData4_15; break;
			default:	theValues=0L;
						isInTables=false;
						break;
		}
	}
	else
	if (r equals 5) {
		switch (n) {
			case 3: theValues=(FriedmanValues*)FriedmanData5_3; break;
			case 4: theValues=(FriedmanValues*)FriedmanData5_4; break;
			case 5: theValues=(FriedmanValues*)FriedmanData5_5; break;
			case 6: theValues=(FriedmanValues*)FriedmanData5_6; break;
			case 7: theValues=(FriedmanValues*)FriedmanData5_7; break;
			case 8: theValues=(FriedmanValues*)FriedmanData5_8; break;
			default:	theValues=0L;
						isInTables=false;
						break;
		}
	}
	else {
		isInTables=false;
	}
	if (isInTables) {
		(*theStrc)=new FriedmanStrc;
		int nS=theValues[0].S;
		(*theStrc)->nS=nS;
		(*theStrc)->S=new int[nS];
		(*theStrc)->qdist=new double[nS];
		for (int i=0;i<nS;i++) {
			(*theStrc)->S[i]=theValues[i+1].S;
			(*theStrc)->qdist[i]=theValues[i+1].qdist;
		}
		return true;
	}
	return false;
}

void freeStoreException() {
	error("\nOut of user-controlled memory");
}

/* 
	Calculates Friedmans Chi exactly
	Starts with a single Ftype type with unpermuted centered ranks and frequency 1
	Adds all possible permutations, collates by type with frequencies cumulated
	 per type.
	Increments n by repeating the above for all possible types, etc.
	This is the Kendall-Babington-Smith algorithm.
*/

FriedmanStrc* FriedmanExact(
	int r,
	int n
)
{
    std::set_new_handler(freeStoreException);

		// if a static array exists, copies it into a FriedmanStrc and returns
	FriedmanStrc* aStrc;
	if (GetFriedmanStatic(r,n,&aStrc)) {
		return aStrc;	
	}

		// Otherwise calculate the FriedmanStrc
	int i;
	Ftype** currentTypes=new Ftype*[MAXTYPE];
	Ftype** workingTypes=new Ftype*[MAXTYPE];

		// The base is an array of unpermuted centered ranks
		//		(1-r)/2,(3-r)/2,...,(r-1)/2	    for odd r
		//		(1-r),(3-r),...,(r-1)			for even r
	int *base=new int[r];
	bool even=(0 equals r%2);
	int step=even?2:1;
	int x=even?1-r:(1-r)/2;
	for (i=0;i<r;i++) {
		base[i]=x;
		x+=step;
	}

		// Set currentTypes[0] to a new Ftype with it's R a copy of the base
	int nCurrentTypes=1;
	currentTypes[0]=new Ftype;
	currentTypes[0]->R=new int[r];
	currentTypes[0]->S=0;
	currentTypes[0]->frequency=1.0;
	for (i=0;i<r;i++) {
		int R=currentTypes[0]->R[i]=base[i];
		currentTypes[0]->S+=R*R;
	}

		// Successively add levels of rankings til n have been added
	int m=n;
	while (--m) {
		int nWorkingTypes=0;
		int j;
			// Add all permutations to each of the current types
			// Allocate the results to elements in the workingTypes list, adding
			//	frequencies or inserting new types as needed.
		for (j=0;j<nCurrentTypes;j++) {
				nWorkingTypes=AddPermutations(r,n,currentTypes[j],nWorkingTypes,
					workingTypes,base);
		}
			// Replace the currentTypes list with the workingTypes list in preparation for
			// adding another level to the rankings
		for (j=0;j<nCurrentTypes;j++) {
			DeleteFtype(currentTypes[j]);
		}
		for (j=0;j<nWorkingTypes;j++) {
			currentTypes[j]=workingTypes[j];
		}
		nCurrentTypes=nWorkingTypes;
	}

	delete [] base;
	delete [] workingTypes;

	double totalFrequency=exp((double)(n-1)*loggamma((double)(r+1))); // (r!)^(n-1)
	FriedmanStrc* theStrc=MakeFriedmanStrc(nCurrentTypes,currentTypes,totalFrequency);

	
	for (i=0;i<nCurrentTypes;i++) {
		DeleteFtype(currentTypes[i]);
	}
	delete [] currentTypes;

	return theStrc;
}


	// Returns true if ok to do exact Friedman calculation
bool DoExactFriedman(
	int r,
	int n,
	bool doRho
)
{
	if (doRho) {  // Spearman's rho -- see limit comment below
		return (1<r && r<12);
	}
	else {
			// These limits can be made larger, but run time becomes excessive
			// Note: if the limits are changed, MAXTYPE may need changing, and 
			// frequency in the Ftype struct's may become too large -- be careful
		switch (r) {
			case 2:	return (n<=100); break;
			case 3: return (n<=30); break;
			case 4: return (n<=15); break;
			case 5: return (n<=8); break; // changed from 10 to 8
				// Can't reasonably do larger cases because of array and time limits
			default: return false; break;
		}
	}
	return false;	
}


/*	This is needed because not all even integers are in S[]
	If lower is true:	
		Finds smallest k such that S[k]>SS -- thus qdist[k] is area above SS and
			1-qdist[k] is area up to and including SS
	If lower is false:
		Finds largest k such that S[k]<=SS -- thus qdist[k] is are above and including
			SS.
*/
static int FriedmanFindS(
	int SS,
	int maxS,
	int *S,
	int nS,
	bool lower
)
{

	double guess=(double)SS/(double)maxS;
	int k=(int)((double)(nS-1)*guess);

	bool larger=(lower)?S[k]>SS:S[k]>=SS;
	while (larger) {
		if (k equals 0) {
			return k;
		}
		else
		if (! lower && SS equals S[k]) {
			return k;
		}
		else {
			larger=S[--k]>SS;
			if (! larger) {
				return (lower)?++k:k;
			}
		}
	}
	while (! larger) {
		if (k equals nS-1) {
			return k;
		}
		else {
			larger=(lower)?S[++k]>SS:S[++k]>=SS;
			if (larger) {
				if (lower || S[k] equals SS) {
					return k;
				}
				else {
					return --k;
				}
			}
		}
	}
	return 0;
}

// deletes S and qdist for the current global
// Also deletes the global struct if deleteAll it true
void ClearFriedmanGlobal(
	bool deleteAll
)
{
	delete [] FriedmanCurrentGlobal->theDist->S;
	delete [] FriedmanCurrentGlobal->theDist->qdist;
	delete [] FriedmanCurrentGlobal->theDist;
	if (deleteAll) {
		delete FriedmanCurrentGlobal;
		FriedmanCurrentGlobal=0L;
	}
}


// Check for existence and feasible of the exact calculation.
// Return true if it exists or it is feasible,
// and then set Q for a given s.
static bool CheckFriedmanExactQ(
	int r,
	int n,
	double s,
	double* Q,
	bool lower,
	bool doRho	// When true s is Spearman's rho
)
{
	if (DoExactFriedman(r,n,doRho)) {
		if (! FriedmanCurrentGlobal || FriedmanCurrentGlobal->r!=r ||
			FriedmanCurrentGlobal->n!=n) {
				if (FriedmanCurrentGlobal) {
					ClearFriedmanGlobal(false);	 // Delete current S and theDist
				}
				else {	// Allocate a new global
					FriedmanCurrentGlobal = new FriedmanGlobal;
				}
				FriedmanCurrentGlobal->theDist=FriedmanExact(r,n);
				FriedmanCurrentGlobal->r=r;
				FriedmanCurrentGlobal->n=n;
		}
		int SS;
		if (doRho) {
			SS=int(0.5+((double)(r*(r*r-1))/6.0)*(1.0+s));
		}
		else {
			SS=(int)(0.5+s*(double)(n*r*(r+1))/12.00);
		}
		if (0 equals r%2) {
			SS*=4;		   // For even r, the ranks were multiplied by 2 to keep them integer.
		}
		FriedmanStrc* theDist=FriedmanCurrentGlobal->theDist;
		int nS=theDist->nS;
		int maxS=theDist->S[nS-1];
		int k=FriedmanFindS(SS,maxS,theDist->S,nS,lower);
		*Q=theDist->qdist[k];
		return true;
	}
	else {
		if (FriedmanCurrentGlobal) {
			ClearFriedmanGlobal(true);	// Delete and zero the current global
		}
		return false;
	}

}


// Check for existence and feasible of the exact calculation.
// Return true if it exists or it is feasible,
// and then set F, the probabilty of exactly X.
static bool CheckFriedmanExactF(
	int r,
	int n,
	double X,
	double* F,
	bool lower,
	bool doRho	// When true s is Spearman's rho
)
{
	if (DoExactFriedman(r,n,doRho)) {
		if (! FriedmanCurrentGlobal || FriedmanCurrentGlobal->r!=r ||
			FriedmanCurrentGlobal->n!=n) {
				if (FriedmanCurrentGlobal) {
					ClearFriedmanGlobal(false);	 // Delete current S and theDist
				}
				else {	// Allocate a new global
					FriedmanCurrentGlobal = new FriedmanGlobal;
				}
				FriedmanCurrentGlobal->theDist=FriedmanExact(r,n);
				FriedmanCurrentGlobal->r=r;
				FriedmanCurrentGlobal->n=n;
		}
		int SS;
		if (doRho) {
			SS=int(0.5+((double)(r*(r*r-1))/6.0)*(1.0+X));
		}
		else {
			SS=(int)(0.5+X*(double)(n*r*(r+1))/12.00);
		}
		if (0 equals r%2) {
			SS*=4;		   // For even r, the ranks were multiplied by 2 to keep them integer.
		}
		FriedmanStrc* theDist=FriedmanCurrentGlobal->theDist;
		int nS=theDist->nS;
		int maxS=theDist->S[nS-1];
		int k=FriedmanFindS(SS,maxS,theDist->S,nS,lower); // Includes prob for x
		*F=theDist->qdist[k];
		if (k<nS-1) {
			*F-=theDist->qdist[k+1];
		}
		return true;
	}
	else {
		if (FriedmanCurrentGlobal) {
			ClearFriedmanGlobal(true);	// Delete and zero the current global
		}
		return false;
	}

}







// Exact or approximation to Friedman's chisquared 
// Continuous approx follows Walsh vol III, p403
// Gives lower tial up to an including X

DISTS_API void pFriedmanR(
	double *Xp,
	int *rp,
	int *np,
	int *Np,
	bool *doRhop,
	double *valuep
)
{
	int N=*Np;
	int i;

	for (i=0;i<N;i++)
		valuep[i]=pfrie(Xp[i],rp[i],np[i],doRhop[i]);
}

 double  pfrie(
	double X,
	int r,
	int n,
	bool doRho	 // When true, X is Spearman's rho
)
{
	double Q;

	if (doRho)
		n=2;
	if (r<3 || n<2)
		return NA_REAL;

	double M=(double)(n*n*r*(r*r-1))/12.0;	// Max value of S
	double S;
	if (doRho) {
		S=(M/2.0)*(1.0+X); // X is rho, so rho to S
	}
	else {
		S=(X*(double)(n*r*(r+1)))/12.0;	// X to S
	}

	if (S>M || S<0)
		return NA_REAL;

	long iS=(long)ceil(S); // Round up to even, because S values must be even
	iS=2*(iS/2);
	S=maxm(1L,iS);


	if (CheckFriedmanExactQ(r,n,X,&Q,true,doRho)) {
		return 1.0-Q;	 // Lower tail including X exactly
	}
	else {
		double W=(S-1.0)/(M+2.0); // Corrected for continuity
		double a=(double)(r-1)-2.0/(double)n;
		double b=a*(double)(n-1);
		return 1.0-pbeta(1.0-W,b/2.0,a/2.0,true,false);
	}
}


DISTS_API void uFriedmanR(
	double *Xp,
	int *rp,
	int *np,
	int *Np,
	bool *doRhop,
	double *valuep
)
{
	int N=*Np;
	int i;

	for (i=0;i<N;i++)
		valuep[i]=qfrie(Xp[i],rp[i],np[i],doRhop[i]);
}

 double  qfrie(
	double X,
	int r,
	int n,
	bool doRho	// When true X is Spearman's rho
)
{
	double Q;

	if (doRho)
		n=2;
	if (r<3 || n<2)
		return NA_REAL;


	double M=(double)(n*n*r*(r*r-1))/12.0;	// Max value of S
	double S;
	if (doRho) {
		S=(M/2.0)*(1.0+X); // X is rho, so rho to S
	}
	else {
		S=(X*(double)(n*r*(r+1)))/12.0;	// X to S
	}

	if (S>M || S<0)
		return NA_REAL;


// Changed to conform to R usage
//		S-=2.0; // Because upper tail, including X, is wanted
	long iS=(long)floor(S);	// Round down to even, because S values must be even
	iS=2*(iS/2);
	S=iS;
	S=maxm(1L,iS);



//	if (CheckFriedmanExactQ(r,n,X,&Q,false,doRho)) {
//		return Q;  // Upper tail area including X exactly
// Changed to conform to R usage
	if (CheckFriedmanExactQ(r,n,X,&Q,true,doRho)) {
		return Q;  // Upper tail area excluding X exactly
	}
	else {
		double W=(S-1.0)/(M+2.0); // Corrected for continuity
		double a=(double)(r-1)-2.0/(double)n;
		double b=a*(double)(n-1);
		return pbeta(1.0-W,b/2.0,a/2.0,true,false);
	}
}

	// Mode of Friedman dist
double modefrie(
	int r,
	int n
)
{
	double maxX=(double)(n*(r-1)); // Max value of X
	double mode=0.0;
	double modeVal=0.0;
	int nPoints=128;
	double delta=maxX/(double)(nPoints-1);
	double X=0.0;
	while (nPoints--){
		double curVal=ffrie(X,r,n,false);
		if (modeVal<curVal) {
			modeVal=curVal;
			mode=X;
		}
		X+=delta;
	}
	return mode;
}


DISTS_API void dFriedmanR(
	double *Xp,
	int *rp,
	int *np,
	int *Np,
	bool *doRhop,
	double *valuep
)
{
	int N=*Np;
	int i;

	for (i=0;i<N;i++)
		valuep[i]=ffrie(Xp[i],rp[i],np[i],doRhop[i]);
}


	// Exact probability or density
 double ffrie(
	double X,
	int r,
	int n,
	bool doRho
)
{
	double F;

	if (doRho)
		n=2;
	if (r<3 || n<2)
		return NA_REAL;

	double M=(double)(n*n*r*(r*r-1))/12.0;	// Max value of S
	double S;
	if (doRho) {
		S=(M/2.0)*(1.0+X); // X is rho, so rho to S
	}
	else {
		S=(X*(double)(n*r*(r+1)))/12.0;	// X to S
	}

	if (S>M || S<0)
		return NA_REAL;


	S-=2.0; // Because upper tail, including X, is wanted
	long iS=(long)floor(S);	// Round down to even, because S values must be even
	iS=2*(iS/2);
	S=iS;
		S=maxm(1L,iS);



	if (CheckFriedmanExactF(r,n,X,&F,false,doRho)) {
		return F;	// Exact probabilty at X
	}
	else {
		double W=(S-1.0)/(M+2.0); // Corrected for continuity
		double a=(double)(r-1)-2.0/(double)n;
		double b=a*(double)(n-1);
			// difference in densities at X
		return pbeta(1.0-W,b/2.0,a/2.0,true,false)-pbeta(1.0-W-2.0/(M+2.0),b/2.0,a/2.0,true,false);
	}
}


DISTS_API void qFriedmanR(
	double *pp,
	int *rp,
	int *np,
	int *Np,
	bool *doRhop,
	double *valuep
)
{
	int N=*Np;
	int i;

	for (i=0;i<N;i++)
		valuep[i]=xfrie(pp[i],rp[i],np[i],doRhop[i]);
}



double  xfrie(
	double P,
	int r,
	int n,
	bool doRho		// When true X is Spearman's rho
)
{
	if (doRho)
		 n=2;
	if (r<3 || n<2)
		return NA_REAL;

		// Initial guess
	double M=(double)(n*n*r*(r*r-1))/12.0;	// Max value of S
	double a=(double)(r-1)-2.0/(double)n;
	double b=a*(double)(n-1);
	double W=1.0-qbeta(1.0-P,b/2.0,a/2.0,true,false); // W is corrected for continuity
	double S=1.0+(M+2.0)*W;
	long iS=(long)ceil(S);


	if (0>P || P>1)
		return NA_REAL;

	iS=2*(iS/2);  // Round up to even, because S values must be even
	S=iS;
	S=maxm(1L,iS);
	double X;
	double step;

	step=12.0/(double)(n*r*(r+1));
	double MX=M*step; // Max vlaue of X
	X=S*step;
	X=maxm(0,X);
	X=minm(MX,X);

		// X must be such that P<=pfrie()
	bool larger=(P<=pfrie(X,r,n,false));
	while (larger) {
		if (X <= 0.0) {
			X=0.0;
			goto theEnd;
		}
		X-=step;
		if (X<0.0) {
			larger=false;
		}
		else
			larger=(P<=pfrie(X,r,n,false));
		if (! larger) {
			X+=step;
			goto theEnd;
		}
	}
	while (! larger) {
		if (X+step >= MX) {
			X=MX;
			goto theEnd;
		}
		X+=step;
		larger=(P<=pfrie(X,r,n,false));
		if (larger) {
			goto theEnd;
		}
	}
	X=0.0;

theEnd:
	if (doRho) {
		return (X/(double)(r-1)-1.0);
	}
	else
		return X;

}


/*
	Friedman random numbers
*/

DISTS_API void rFriedmanR(
	int *rp,
	int *np,
	bool *doRhop,
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
		rfrie(valuep,*Np,*rp,*np,*doRhop);
	else { // Allow for random values for each element of nu and lambda
		D=(N/M)+((N%M)?1:0);
		tArray=(double *)S_alloc((long)D,sizeof(double));
		loc=0;
		for (j=0;j<M;j++) {
			rfrie(tArray,D,rp[j],np[j],doRhop[j]);
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

 void rfrie(
	double* randArray,
	int N,	  // number of samples
	int r,
	int n,
	bool doRho		// When true X is Spearman's rho
)
{
	GetRNGstate();


	for (int i=0;i<N;i++){
		randArray[i]=(double)xfrie(unif_rand(),r,n,doRho);
	}
	PutRNGstate();
}



	// Returns the median of the Friedman dist -- this is never called for Spearman's Rho,
	//	because the median for Rho is 0.
double medianfrie(
	int r,
	int n
)
{
	if (! DoExactFriedman(r,n,false)) {
		return xfrie(0.5,r,n,false);	// continuous approx
	}
	else {
		double high=xfrie(0.5,r,n,false);
		double phigh=pfrie(high,r,n,false);
		double low=high;
		bool same=true;
		double plow;
		double twoS=24.0/(double)(n*r*(r+1)); // Corresponds to 2 on the S scale
		if (0 equals r%2) {
			twoS*=4.0;
		}
			// find the next lower Friedman chi-square
		while (same) {
			low-=twoS;	
			same=(phigh equals (plow=pfrie(low,r,n,false)));
		}
		double alpha=(phigh-0.5)/(phigh-plow);
		return alpha*low+(1.0-alpha)*high;	// interpolate
	}
}
 
DISTS_API void sFriedmanR(
	int *rp,
	int *np,
	bool *rhop,
	int *Np,
	double *meanp,
	double *medianp,
	double *modep,
	double *variancep,
	double *thirdp,
	double *fourthp
)
{
	int N=*Np;
	int i;

	for (i=0;i<N;i++) {
		if (rp[i]<3 || (! rhop[i] && np[i]<2)) {
			meanp[i]=NA_REAL;
			medianp[i]=NA_REAL;
			modep[i]=NA_REAL;
			variancep[i]=NA_REAL;
			thirdp[i]=NA_REAL;
			fourthp[i]=NA_REAL;
		}
		else {
			if (rhop[i]) {
				meanp[i]=0.0;
				medianp[i]=0.0;
				modep[i]=0.0;
				variancep[i]=1.0/(double)(rp[i]-1);
				thirdp[i]=0.0;
				fourthp[i]=(3.0*variancep[i]/(double)(rp[i]-1))*
						((double)(72+rp[i]*(-35+rp[i]*(-38+rp[i]*25)))/(double)(25*rp[i]*(rp[i]*rp[i]-1)));

			}
			else {
				meanp[i]=(double)(rp[i]-1);
				medianp[i]=medianfrie(rp[i],np[i]);
				modep[i]=modefrie(rp[i],np[i]);
				variancep[i]=(double)(2*(rp[i]-1)*(np[i]-1))/(double)np[i];
				thirdp[i]=variancep[i]*(double)(4*(np[i]-2))/(double)np[i];
				fourthp[i]=variancep[i]*((double)(rp[i]-1)/(double)(np[i]*np[i]))*
						((double)(72+rp[i]*(-35+rp[i]*(-38+rp[i]*25)))/(double)(25*rp[i]*(rp[i]*rp[i]-1))+
						(double)(2*(rp[i]-1)*(np[i]-2))+0.5*(double)((rp[i]+3)*(np[i]-2)*(np[i]-3)));
			}
		}
	}

}






		
/********************************************************************************
	Fits a Johnson curve from percentiles
	See Wheeler 1980 Biometrika
	The input is in the struct JohnsonInput
	The output is in the struct JohnsonParms 
*/

 /* 
|	Rotates the 3 element row v into matrix using Givens rotations
|	matrix is 3 by 3 and should be zeroed before the start
|	The weights are on the diagonal
*/

static void	Rotate3(
	double *v,
	double matrix[3][3]
)
{
	double d;
	double dp;
	double c;
	double s;
	double x;
	double r;
	double w;
	int i;
	int j;
	bool skip;

	skip=false;
	w=1.0;
	for (i=0;i<2;i++) {
      if (!skip) {
			if (0.0 equals (x=v[i])) continue;
			d=matrix[i][i];
			dp=d+w*x*x;
			matrix[i][i]=dp;
			c=d/dp;
			s=w*x/dp;
			if (d equals 0.0) skip=true;      /* to avoid 0/0, but d can't be 0 */
			else w*=c;
			for (j=i+1;j<3;j++) {
				matrix[i][j]=s*v[j]+c*(r=matrix[i][j]);
			   v[j]-=x*r;
			}
      }
	}
}



	// Fits the Johnson Su
void JohnsonMomentSu(
	JohnsonParms& parms,
	double mean,
	double sd,
	double sqrtB1,
	double B2
)
{
	const double TOLSU=0.01;
	double B1=sqrtB1*sqrtB1;
	double B2Minus3=B2-3.0;
	
		// First estimate of exp(delta^-2)
	double w=sqrt(sqrt(2.0*B2-2.8*B1-2.0)-1.0);

	double m=0.0;
	double z;
	if (fabs(sqrtB1)>TOLSU) {
		// Johnson Iteration 
		int count=0;
		repeat
			double wPlus1=w+1.0;
			double wMinus1=w-1.0;
			z=wPlus1*B2Minus3;
			double v=w*(6.0+w*(3.0+w));
			double a=8.0*(wMinus1*(3.0+w*(7.0+v))-z);
			double b=16.0*(wMinus1*(6.0+v)-B2Minus3);
			m=(sqrt(a*a-2.0*b*(wMinus1*(3.0+w*(9.0+w*(10.0+v)))-2.0*wPlus1*z))-a)/b;
			z=(4.0*(w+2.0)*m+3.0*wPlus1*wPlus1);
			double h=(2.0*m+wPlus1);
			z=m*wMinus1*(z*z)/(2.0*h*h*h);
			v=w*w;
			w=sqrt(1.0-2.0*(1.5-B2+(B1*(B2-1.5-v*(1.0+0.5*v)))/z));	
			w=sqrt(w-1.0);
		until(fabs(B1-z)<=TOLSU || count++>100);
		if (count>100) {
			error("\nToo many iterations");
			return;
		}
		m/=w;
		m=log(sqrt(m)+sqrt(1.0+m));
		m=(sqrtB1>0.0)?-m:m; 
	}

	parms.delta=sqrt(1.0/log(w));
	double omega=m;
	parms.gamma=omega*parms.delta;

	parms.lambda=sd/(sqrt(0.5*(w-1.0)*(w*cosh(2.0*omega)+1.0)));

	parms.xi=mean+(0.5*sqrt(w)*sinh(omega))*parms.lambda;

	parms.type=SU;
}


  	// Finds first N moments of Johnson SB, puts them in moments[]
bool JohnsonMOM(
	double gamma,
	double delta,
	double* moments
)
{
	const int N=6;
	const int iterLimit=500;
	const double innerTol=1e-8;
	const double outerTol=1e-5;
	const double ln10=2.3025850929940456840179915;
	const double recipSqrtPi=0.5641895835477562869480795;
	const double expA=ln10*DBL_MAX_10_EXP;
	const double expB=-log(DBL_EPSILON);

	bool howExit=true;
	double oldMoments[N];
	memset(oldMoments,0,N*sizeof(double));
	double w=gamma/delta;

	if (w>expA) {
		return false;
	}

	double h=0.75;
	double expW=exp(w)+1.0;
	if (delta<3.0) {
		h=delta/4.0;
	}
	int count=0;
	bool more=false;
		// The outer loop
	repeat
			// Skip this the first time
		if (count) {
			for (int i=0;i<N;i++) {
				oldMoments[i]=moments[i];
			}
			// No convergence, try a smaller h
			h*=0.5;
		 }

		 double t=w;
		 double u=w;
		 double y=h*h;
		 double x=2.0*y;
		 moments[0]=1.0/expW;
		 int i;
		 for (i=1;i<N;i++) {
		 	moments[i]=moments[i-1]/expW;
		 }
		 double v=y;
		 double f=SQRT2*h/delta;

		 int countInner=0;
		 	// Inner loop to evaluate infinite series

		 repeat
			double b[N];
			for (i=0;i<N;i++) {
				b[i]=moments[i];
			}
			u-=f;
			double z=1.0;
			if (u>-expB) {
				z+=exp(u);
			}
			t+=f;
			bool bL=t>expB;
			double s=0;
			if (! bL) {
				s=exp(t)+1.0;
			}
			double p=exp(-v);
			double q=p;

			for (i=0;i<N;i++) {
				double aa=moments[i];
				p/=z;
				double ab=aa;
				aa+=p;
				if (aa equals ab) {
					break;
				}
				if (! bL) {
					q/=s;
					ab=aa;
					aa=aa+q;
					bL=(aa equals ab);
				}
				moments[i]=aa;
  			}
			y+=x;
			v+=y;
			more=false;
			for (i=0;i<N;i++) {
				if (moments[i] equals 0.0) {
					goto errExit;
				}
				if (fabs(moments[i]-b[i])/moments[i] > innerTol) {
					more=true;
				}
			}

		 until(! more || countInner++>iterLimit);
		 if (more) {
		 	goto errExit;
		 }

		 for (i=0;i<N;i++) {
		 	moments[i]*=recipSqrtPi*h;
		 }
		 more=false;
		 for (i=0;i<N;i++) {
		 	if (moments[i] equals 0.0) {
				goto errExit;
			}
			if (fabs(moments[i]-oldMoments[i])/moments[i]>outerTol) {
				more=true;
			}
		 } 

	until(! more || count++>iterLimit);
	if (! more) {
		goto theExit;
	}

errExit:
	howExit=false;
theExit:
	return howExit;
}

	// Fits the Johnson Sb
bool JohnsonMomentSb(
	JohnsonParms& parms,
	double mean,
	double sd,
	double sqrtB1,
	double B2
)
{
	const double TOLSB=0.01;
	const double tTol=1e-2;
	const int iterLimit=50;
	const int N=6;		// Number of moments

	bool howExit=true;
	double absSqrtB1=fabs(sqrtB1);
	double B1=sqrtB1*sqrtB1;
	bool negativeB1=(sqrtB1<0.0);
	
		// Get first estimate of delta
	double delta;
	double x=1.0+0.5*B1;
	double y=absSqrtB1*sqrt(1.0+B1/4.0);
	double w=pow(x+y,1.0/3.0)+pow(x-y,1.0/3.0)-1.0;
	double f=w*w*(3.0+w*(2.0+w))-3.0;
	double e=1.0+B1;
	e=(B2-e)/(f-e);
	if (absSqrtB1<=TOLSB) {
		f=2.0;
	}
	else {
		delta=1.0/sqrt(log(w));
		if (delta<0.64) {
			f=1.25*delta;
		}
		else {
			f=2.0-8.5245/(delta*(delta*(delta-2.163)+11.346));
		}
	}
	f=1.0+e*f;
	if (f<1.8) {
		delta=0.8*(f-1.0);
	}
	else {
		delta=(0.626*f-0.408)*pow((3.0-f),-0.479);
	}

		// Get first estimate of gamma
	double gamma=0.0;
	if (B1>=tTol) {
		if (delta<=1.0) {
			gamma=(0.7466*pow(delta,1.7973)+0.5955)*pow(B1,0.485);
		}
		else {
			if (delta<=2.5) {
				gamma=pow(B1,0.0623*delta+0.4043);
			}
			else {
				gamma=pow(B1,0.0124*delta+0.5291);
			}
			gamma*=(0.9281+delta*(1.0614*delta-0.7077));
		}
	}

	int count=0;
	bool more=false;
	bool errorSet=false;
	double moments[N];
	double h2=0;
	double oldDeltaGamma=100.0;
	double oldDeltaDelta=100.0;
	repeat
			// Get first N moments for latest delta and gamma
		if (JohnsonMOM(gamma,delta,moments)) {
			h2=moments[1]-moments[0]*moments[0];
			if (h2>0.0) {

				double h2a=sqrt(h2)*h2;
				double h2b=h2*h2;
				double h3=moments[2]-moments[0]*(3.0*moments[1]-2.0*moments[0]*moments[0]);
				double rbet=h3/h2a;
				double h4=moments[3]-moments[0]*(4.0*moments[2]-moments[0]*
					(6.0*moments[1]-3.0*moments[0]*moments[0]));
				double bet2=h4/h2b;
				double w=gamma*delta;
				double u=delta*delta;

					// Get derivatives
				double dd[N];
				double deriv[N];

				for (int j=0;j<2;j++) {
					for (int k=0;k<4;k++) {
						double t=(double)k;
						double s;
						if (! j) {
							s=moments[k+1]-moments[k];	
						}
						else {
							s=((w-t)*(moments[k]-moments[k+1])+(1.0+t)*
								(moments[k+1]-moments[k+2]))/u;
						}
						dd[k]=t*s/delta;
					}
					double t=2.0*moments[0]*dd[0];
					double s=moments[0]*dd[1];
					double y=dd[1]-t;
					deriv[j]=(dd[2]-3.0*(s+moments[1]*dd[0]-t*moments[0])-
								1.5*h3*y/h2)/h2a;
					deriv[j+2]=(dd[3]-4.0*(dd[2]*moments[0]+dd[0]*moments[2])+
								6.0*(moments[1]*t+moments[0]*(s-t*moments[0]))-
								2.0*h4*y/h2)/h2b;
				}
				double t=1.0/(deriv[0]*deriv[3]-deriv[1]*deriv[2]);
				double deltaGamma=(deriv[3]*(rbet-absSqrtB1)-deriv[1]*(bet2-B2))*t;
				double deltaDelta=(deriv[0]*(bet2-B2)-deriv[2]*(rbet-absSqrtB1))*t;

					// New estimates of gamma and delta
				gamma-=deltaGamma;
				if (B1 equals 0.0 || gamma<0.0) {
					gamma=0.0;
				}
				delta-=deltaDelta;
				deltaGamma=fabs(deltaGamma);
				deltaDelta=fabs(deltaDelta);
				more=(deltaGamma>tTol || deltaDelta>tTol);
					// error if iterates increase
				errorSet=(deltaGamma>oldDeltaGamma || deltaDelta>oldDeltaDelta);
				oldDeltaGamma=deltaGamma;
				oldDeltaDelta=deltaDelta;
 	  		}
		}
	until(! more || errorSet || count++>iterLimit);
	if (! errorSet && ! more) {
		parms.delta=delta;
		parms.lambda=sd/sqrt(h2);
		if (negativeB1)	{
			gamma=-gamma;
			moments[0]=1.0-moments[0];
		}
		parms.gamma=gamma;
		parms.xi=mean-parms.lambda*moments[0];
		parms.type=SB;

 		goto theExit;
	}

	howExit=false;
theExit:
	return howExit;
}

#ifdef NOTUSED
	// Surprisingly the third and fourth moments are not as accurate as by simple summation
	//	of powers

	// Finds first four moments
	// Uses updating. See Welford, B.P. (1962) Note on a method for calculating corrected
	// sums of squares and products. Technometrics 4-3. 419-420

DISTS_API momentsR(
	double *data,
	int *Np,
	double *meanp,
	double *variancep, // n ! (n-1) divisor used
	double *thirdp,
	double *fourthp
)
{
	int N=*Np;
	int i;
	int n=1;

	double mean=data[0];
	double variance=0;
	double third=0;
	double fourth=0;
	double delta;
	double scale;
	for (i=1;i<N;i++) {
		n++;
		scale=(double)(n-1)/(double)n;
		delta=(data[i]-mean)/n;
		fourth=scale*((((n*(3+(n*(n-3))))*delta*delta+6*variance)*delta-4*third)*delta+fourth);
		third=scale*(((n-2)*delta*delta+3*variance)*delta+third);
		variance=scale*(variance+n*delta*delta);
		mean+=delta;
	}

	*meanp=mean;
	*variancep=variance;
	*thirdp=third;
	*fourthp=fourth;
	
}
#endif

/*
	Fits Johnson curves by moments.
	This is a translation of AS99 by I.D. Hill, et al.
	Throws an exception in case of error
*/

DISTS_API void JohnsonMomentFitR(
	double *meanp,
	double *sdp,
	double *skewp,
	double *kurtp,
	double *gammap,
	double *deltap,
	double *xip,
	double *lambdap,
	int *typep
)
{
	JohnsonParms parms;
	JohnsonMoments moments;

		moments.mean=*meanp;
		moments.sd=*sdp;
		moments.sqrtB1=*skewp;
		moments.B2=3.0+*kurtp;
		parms=JohnsonMomentFit(moments);
		*gammap=parms.gamma;
		*deltap=parms.delta;
		*xip=parms.xi;
		*lambdap=parms.lambda;
		*typep=1+(int)parms.type;
}



 JohnsonParms JohnsonMomentFit(
	JohnsonMoments moments
)
{
	 JohnsonParms parms={0.0,0.0,0.0,0.0,SN};
	double mean=moments.mean;
	double sd=moments.sd;
	double sqrtB1=moments.sqrtB1;
	double B2=moments.B2;

	const double TOLJ=0.1; // changed from 0.01 
	double B1=sqrtB1*sqrtB1;

	if (B2<B1+1.0+TOLJ) {
		error("\nMoment ratio in error");
		return parms;  // Outside the upper limit.
	}

		// Is it a normal ?
	if (fabs(sqrtB1)<=TOLJ && fabs(B2-3.0)<=TOLJ) {
		parms.type=SN;
		parms.gamma=0.0;	// Hill et.al. let gamma delta carry the burden
		parms.delta=1.0;
		parms.lambda=sd;
		parms.xi=mean;
		return parms;
	}
		// Test for position relative to log normal line
		// The log normal line is defined by a pair of equations:
		//	B1=(w-1)(w+2)^2, and B2=w^4+2w^3+3w^2-3.

		// First solve the cubic (w-1)(w+2)^2=B1 then estimate B2 and check the
		//	estimate for closeness to the actual B2

		// This method of solution appears in the statlib copy.  The original
		// publication used a, presumably, less stable method.
	double x=0.5*B1+1.0;
	double y=sqrt(B1+0.25*B1*B1);
	double u=pow(x+y,1.0/3.0);
		// Evaluate w and test B2 equals B2est=w^4+2w^3+3w^2-3
	double w=u+1.0/u-1.0;
	double B2est=w*w*(3.0+w*(2.0+w))-3.0;

	B2=(B2<0.0)?B2est:B2;  // Log fit for neg B2

	double test=B2est-B2;

		 // Is equation close to B2?  If so do a log fit
	if (fabs(test)<TOLJ) {
		parms.type=SL;
		parms.lambda=1.0; 
		parms.delta=1.0/sqrt(log(w));
			// Need the commented out parts only if lambda is not 1
		parms.gamma=0.5*parms.delta*log(w*(w-1.0)/(sd*sd))/*-parms.delta*log(parms.lambda)*/;
		parms.xi=mean-/*parms.lambda* */sd/sqrt(w-1.0);
		return parms;
	}

		 // Su fit
	if (test<=0.0) {
		JohnsonMomentSu(parms,mean,sd,sqrtB1,B2);
		return parms;
	}
		// Sb fit
	else {
		if (JohnsonMomentSb(parms,mean,sd,sqrtB1,B2)) {
			return parms;
		}
		else {
			error("\nCouldn't do an Sb fit");
			return parms;
		}
	}
}


/*
	Fits the Johnson curves from quantiles
*/

DISTS_API void JohnsonFitR(
	double *xnp,
	double *xmp,
	double *x0p,
	double *xkp,
	double *xpp,
	double *gammap,
	double *deltap,
	double *xip,
	double *lambdap,
	int *typep
)
{
	JohnsonParms parms;
	JohnsonInput input;

	input.x0=*x0p;
	input.xk=*xkp;
	input.xm=*xmp;
	input.xn=*xnp;
	input.xp=*xpp;
	parms=JohnsonFit(input);
	*gammap=parms.gamma;
	*deltap=parms.delta;
	*xip=parms.xi;
	*lambdap=parms.lambda;
	*typep=1+(int)parms.type;
}


 JohnsonParms JohnsonFit(
	JohnsonInput input
)
{
	double xn=input.xn;
	double xm=input.xm;
	double x0=input.x0;
	double xk=input.xk;
	double xp=input.xp;
	double zn=1.64485363;

	double t;
	double tu;
	double tb;
	double tbu;
	double delta;
	double gamma=0;
	double xi;
	double lambda;
	double matrix[3][3];
	double array[5][3];
	JohnsonType solution;

	memset(matrix,0,9*sizeof(double));

	const double TOLJN=0.1;

	t=(xn-x0)/(x0-xp);
	tu=(xn-xp)/(xm-xk);
	tb=0.5*(
		((xm-x0)*(xn-xp))/((xn-xm)*(x0-xp))+((xk-x0)*(xp-xn))/((xp-xk)*(x0-xn))
	);

	tbu=tb/tu;
			// Normal solution
	if (fabs(fabs(tbu)-1.0)<TOLJN && fabs(fabs(t-1.0))<TOLJN) {  
		solution=SN;
		delta=1.0;
		gamma=0.0;
	}
	else
		// Log solution
	if (fabs(fabs(tbu)-1.0)<TOLJN) {  
		solution=SL;
		delta=zn/log(t);
		if (! R_FINITE(delta)){
			error("\nInfinite value in SL fit");
		}
	}
	else
		// Bounded solution
	if (tbu>1.0) {  
		solution=SB;
		tb*=0.5;
		double b=tb+sqrt(tb*tb-1.0);
		delta=zn/(2.0*log(b));
		b*=b;
 		if (t>b || t<1.0/b) {
			error("\nBounded solution intermediate values out of range");
		}
		double a=(t-b)/(1-t*b);
		gamma=-delta*log(a);
	}
		// Unbounded solution
	else {   
		solution=SU;
		tu*=0.5;
		double b=tu+sqrt(tu*tu-1.0);
		delta=zn/(2.0*log(b));
		b*=b;
		if (t>b || t<1.0/b) {
			error("\nUnbounded solution intermediate values out of range");
		}
		double a=(1-t*b)/(t-b);
		gamma=-0.5*delta*log(a);
	}

		// Find xi and lambda by least squares
	array[0][1]=zn;
	array[0][2]=xn;
	array[1][1]=zn/2.0;
	array[1][2]=xm;
	array[2][1]=0.0;
	array[2][2]=x0;
	array[3][1]=-zn/2.0;
	array[3][2]=xk;
	array[4][1]=-zn;
	array[4][2]=xp;
	for (int i=0;i<5;i++) {
		array[i][0]=1.0;
		double u=array[i][1];
		if (solution != SN) {
			if (solution equals SL) {
				u=exp(u/delta);
			}
			else {
				u=exp((u-gamma)/delta);
				if (solution equals SB) {
					u=u/(1.0+u);
				}
				else {
					u=(u*u-1.0)/(2.0*u);
				}
			}
		}
		array[i][1]=u;
		Rotate3(array[i],matrix);
	}
	lambda=matrix[1][2];
	xi=matrix[0][2]-lambda*matrix[0][1];

	JohnsonParms output;
	output.gamma=gamma;
	output.delta=delta;
	output.xi=xi;
	output.lambda=lambda;
	output.type=solution;

	return output;
}

/*
	Johnson lower probability
*/

DISTS_API void pJohnsonR(
	double *xp,
	double *gammap,
	double *deltap,
	double *xip,
	double *lambdap,
	int *typep,
	int *Np,
	double *valuep
)
{
	int N=*Np;
	int i;
	JohnsonParms parms;

	for (i=0;i<N;i++) {
		parms.gamma=gammap[i];
		parms.delta=deltap[i];
		parms.xi=xip[i];
		parms.lambda=lambdap[i];
		parms.type=(JohnsonType)(typep[i]-1);
		valuep[i]=pjohnson(xp[i],parms);
	}

}


 double pjohnson(
	double x,
	JohnsonParms parms
)
{
	double u=(x-parms.xi)/parms.lambda;

	switch (parms.type) {
		case SN:
			break;
		case SL: 
				u=log(u);
			break;
		case SU: 
				u+=sqrt(1+u*u);
				u=log(u);	
			break;
		case SB:
				if (u<=0.0 || u>=1.0) {
					error("\nSb values out of range.");
					return 0.0;
				}
				u/=(1-u);
				u=log(u);
			break;
		default:
			error("\nNo type");
			break;
	}

	double z=parms.gamma+parms.delta*u;
	return pnorm(z,0,1,true,false);
}

/*
	Johnson upper probability
*/

DISTS_API void uJohnsonR(
	double *xp,
	double *gammap,
	double *deltap,
	double *xip,
	double *lambdap,
	int *typep,
	int *Np,
	double *valuep
)
{
	int N=*Np;
	int i;
	JohnsonParms parms;

	for (i=0;i<N;i++) {
		parms.gamma=gammap[i];
		parms.delta=deltap[i];
		parms.xi=xip[i];
		parms.lambda=lambdap[i];
		parms.type=(JohnsonType)(typep[i]-1);
		valuep[i]=qjohnson(xp[i],parms);
	}

}

 double qjohnson(
	double x,
	JohnsonParms parms
)
{
	return 1.0-pjohnson(x,parms);
}

/*
	Johnson percentile
*/
DISTS_API void qJohnsonR(
	double *pp,
	double *gammap,
	double *deltap,
	double *xip,
	double *lambdap,
	int *typep,
	int *Np,
	double *valuep
)
{
	int N=*Np;
	int i;
	JohnsonParms parms;

	for (i=0;i<N;i++) {
		parms.gamma=gammap[i];
		parms.delta=deltap[i];
		parms.xi=xip[i];
		parms.lambda=lambdap[i];
		parms.type=(JohnsonType)(typep[i]-1);
		valuep[i]=xjohnson(pp[i],parms);
	}

}

 double xjohnson(
	double p,
	JohnsonParms parms
)
{
	double z=qnorm(p,0,1,true,false);
	double u=(z-parms.gamma)/parms.delta;


	switch (parms.type) {
		case SN:
			break;
		case SL: 
				u=exp(u);
			break;
		case SU: 
				u=exp(u);	
				u=((u*u)-1.0)/(2.0*u);
			break;
		case SB:
				u=exp(u);
				u/=(1.0+u);
			break;
		default:
			break;
	}

	return (parms.xi+parms.lambda*u);

}

/*
	Johnson density
*/

DISTS_API void dJohnsonR(
	double *xp,
	double *gammap,
	double *deltap,
	double *xip,
	double *lambdap,
	int *typep,
	int *Np,
	double *valuep
)
{
	int N=*Np;
	int i;
	JohnsonParms parms;

	for (i=0;i<N;i++) {
		parms.gamma=gammap[i];
		parms.delta=deltap[i];
		parms.xi=xip[i];
		parms.lambda=lambdap[i];
		parms.type=(JohnsonType)(typep[i]-1);
		valuep[i]=fjohnson(xp[i],parms);
	}

}


 double fjohnson(
	double x,
	JohnsonParms parms
)
{
	double u=(x-parms.xi)/parms.lambda;
	double fu=0;
	double ratio=parms.delta/parms.lambda;
	double differential=0;


	switch (parms.type) {
		case SN:
				fu=u;
				differential=ratio;
			break;
		case SL: 
				differential=ratio/u;
				fu=log(u);
			break;
		case SU: 
				fu=u+sqrt(1+u*u);
				differential=ratio/sqrt(1.0+u*u);
				fu=log(fu);	
			break;
		case SB:
				fu=u/(1-u);
				differential=ratio/(u*(1.0-u));
				fu=log(fu);
			break;
		default:
			break;
	}

	double z=parms.gamma+parms.delta*fu;
	return dnorm(z,0,1,false)*differential;
}

/*
	derivitive of Johnson density
*/
double fpjohnson(
	double x,
	JohnsonParms parms
)
{
	double u=(x-parms.xi)/parms.lambda;
	double fu;
	double ratio=parms.delta/parms.lambda;
	double differential=0;

	double z=0;
	double w;

	switch (parms.type) {
		case SN:
				z=parms.gamma+parms.delta*u;
				differential=-ratio*ratio*z;
			break;
		case SL: 
				z=parms.gamma+parms.delta*log(u);
				differential=-(z+1.0/parms.delta)*(ratio/u)*(ratio/u);
			break;
		case SU: 
				fu=u+sqrt(1+u*u);
				z=parms.gamma+parms.delta*log(fu);
				w=1.0/sqrt(1.0+u*u);
				differential=((ratio*w*w)/parms.lambda)*(w/fu-1-z*parms.delta);	
			break;
		case SB:
				fu=u/(1-u);
				z=parms.gamma+parms.delta*log(fu);
				w=1.0/((1.0-u)*(1.0-u));
				differential=((ratio*w)/parms.lambda)*(2.0/fu-(z*parms.delta+1.0)/(u*u));
			break;
		default:
			break;
	}

	return dnorm(z,0,1,false)*differential;
}




	// Translates a normal deviate to a johnson deviate
double xzjohnson(
	double z,
	JohnsonParms parms
)
{
	double u=(z-parms.gamma)/parms.delta;


	switch (parms.type) {
		case SN:
			break;
		case SL: 
				u=exp(u);
			break;
		case SU: 
				u=exp(u);
				u=((u*u)-1.0)/(2.0*u);	
			break;
		case SB:
				u=exp(u);
				u/=(1.0+u);
			break;
		default:
			break;
	}

	return (parms.xi+parms.lambda*u);

}


/*
	Johnson random numbers
*/

DISTS_API void rJohnsonR(
	double *gammap,
	double *deltap,
	double *xip,
	double *lambdap,
	int *typep,
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
	JohnsonParms parms;

	if (M==1) {
		parms.gamma=*gammap;
		parms.delta=*deltap;
		parms.xi=*xip;
		parms.lambda=*lambdap;
		parms.type=(JohnsonType)(*typep-1);
		rjohnson(valuep,N,parms);
	}
	else { // Allow for random values for each element of the parameters
		D=(N/M)+((N%M)?1:0);
		tArray=(double *)S_alloc((long)D,sizeof(double));
		loc=0;
		for (j=0;j<M;j++) {
			parms.gamma=gammap[j];
			parms.delta=deltap[j];
			parms.xi=xip[j];
			parms.lambda=lambdap[j];
			parms.type=(JohnsonType)(typep[j]-1);
			rjohnson(tArray,D,parms);
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

void rjohnson(
	double* johnsonArray,
	int n,
	JohnsonParms parms
)
{

	rgauss(johnsonArray,n,0.0,1.0);
	for (int i=0;i<n;i++) {
		johnsonArray[i]=xzjohnson(johnsonArray[i],parms);
	}
}


double FindDistributionStatistic(
	double lowX,
	double highX,
	double (*function)(double x)
)
{
	return Integral(lowX,highX,function,0.0001);
}

double FindDistributionMode(
	double lowX,
	double highX,
	double (*function)(double x)
)
{
	double mode=0.0;
	double fmode=-1.0;
	int N=128;
	double StepTheKey=(highX-lowX)/(double)(N-1);
	double X=lowX;
	double fvalue;
	for (int i=0;i<N;i++) {
		fvalue=(*function)(X);
		if (fvalue>fmode)  {
			fmode=fvalue;
			mode=X;
		}
		X+=StepTheKey;	
	}

	return mode;
 }

DISTS_API void sJohnsonR(
	double *gammap,
	double *deltap,
	double *xip,
	double *lambdap,
	int *typep,
	int *Np,
	double *meanp,
	double *medianp,
	double *modep,
	double *variancep,
	double *thirdp,
	double *fourthp
)
{
	int N=*Np;
	int i;
	JohnsonParms parms;

	for (i=0;i<N;i++) {
		parms.gamma=gammap[i];
		parms.delta=deltap[i];
		parms.xi=xip[i];
		parms.lambda=lambdap[i];
		parms.type=(JohnsonType)(typep[i]-1);
		sJohnson(parms,meanp+i,medianp+i,modep+i,variancep+i,thirdp+i,fourthp+i);
	}

}

static JohnsonParms gparms;

static double MeanJFcn(double x){return x*fjohnson(x,gparms);}
static double VarianceJFcn(double x){return (x-gmean)*(x-gmean)*fjohnson(x,gparms);}
static double ThirdMomentJFcn(double x){return (x-gmean)*(x-gmean)*(x-gmean)*fjohnson(x,gparms);}
static double FourthMomentJFcn(double x){return (x-gmean)*(x-gmean)*(x-gmean)*(x-gmean)*fjohnson(x,gparms);}
static double AJFunction(double x){return fjohnson(x,gparms);}

void sJohnson(
	JohnsonParms parms,
	double *meanp,
	double *medianp,
	double *modep,
	double *variancep,
	double *thirdp,
	double *fourthp
)
{
	if (fabs(parms.delta)<1e-13) {
		error("\nSorry, can't do it");
		return;
	}

	gparms=parms;



 	double gamma=parms.gamma;
	if (fabs(gamma)<1e-15) {gamma=0.0;}
	double delta=parms.delta;
	double xi=parms.xi;
	if (fabs(xi)<1e-15) {xi=0.0;}
	double lambda=parms.lambda;

	double lowX;
	double highX;
	double mean=0;
	double median=0;
	double mode=0;
	double variance=0;
	double thirdMoment=0;
	double fourthMoment=0;
	double w=exp(1.0/(delta*delta));
	double gd=gamma/delta;
	double egd=exp(-gd);


	switch (parms.type)	{
		case SN:
				mean=xi-gamma*lambda/delta;
				median=mode=mean=mean;
				variance=lambda/delta;
				variance*=variance;
				thirdMoment=0.0;
				fourthMoment=3.0*variance*variance;
			break;
		case SL:
				mean=xi+sqrt(w)*egd*lambda;
				median=xi+egd*lambda;
				mode=xi+lambda*egd/w;
				variance=w*(w-1.0)*egd*egd*lambda*lambda;
				thirdMoment=egd*egd*egd*sqrt(w*w*w)*(w-1.0)*(w-1.0)*(w+2.0);
  				thirdMoment*=lambda*lambda*lambda;
 				fourthMoment=variance*variance*(((w+2.0)*w+3.0)*w*w-3.0);
			break;
		case SU:
				lowX=xjohnson(0.001,parms);
				highX=xjohnson(0.999,parms);
				mean=xi-lambda*sqrt(w)*sinh(gd);
				variance=0.5*lambda*lambda*(w-1.0)*(w*cosh(2.0*gd)+1.0);
				median=xi-lambda*sinh(gd);
				mode=FindDistributionMode(lowX,highX,AJFunction);
				thirdMoment=0.25*sqrt(w)*(w-1.0)*(w-1.0)*(w*(w+2.0)*sinh(3.0*gd)+
					3.0*sinh(gd));
				thirdMoment*=lambda*lambda*lambda;
				thirdMoment=(gamma>=0)?-thirdMoment:thirdMoment;
				fourthMoment=0.125*(w-1.0)*(w-1.0)*(w*w*(((w+2.0)*w+3.0)*w*w-3.0)*cosh(4.0*gd)+
					4.0*w*w*(w+2.0)*cosh(2.0*gd)+3.0*(2.0*w+1.0));
				fourthMoment*=lambda*lambda*lambda*lambda;
			break;
		case SB:
				lowX=xjohnson(0.001,parms);
				highX=xjohnson(0.999,parms);
				mode=FindDistributionMode(lowX,highX,AJFunction);
				mean=FindDistributionStatistic(lowX,highX,MeanJFcn);
				gmean=mean;
				variance=FindDistributionStatistic(lowX,highX,VarianceJFcn);
				thirdMoment=FindDistributionStatistic(lowX,highX,ThirdMomentJFcn);
				fourthMoment=FindDistributionStatistic(lowX,highX,FourthMomentJFcn);
				median=xjohnson(0.5,parms);
			break;
	}
	*meanp=mean;
	*medianp=median;
	*modep=mode;
	*variancep=variance;
	*thirdp=thirdMoment;
	*fourthp=fourthMoment;

}

 /*******************************************************************************
  	Density function of the correlation coefficient
	From eq 6.5, p223 of Johnson and Kotz, Continuous Univariate Distributions,
	 vol 2, 1st edition.  Uses Gaussian hypergeometric function, defined on
	 page 17 (1.104) of Johnson, Kotz, and Kemp, Univariate Discrete Distributions,
	 2nd Ed.  The derivitive is given on page 18.
 */

DISTS_API void dcorrR(
	double *pp,
	double *rhop,
	int *np,
	int *Np,
	double *valuep
)
{
	int N=*Np;
	int i;

	for (i=0;i<N;i++)
		valuep[i]=fcorrelation(pp[i],rhop[i],np[i]);
}

 double fcorrelation(
 	double r,	   // The correlation coefficient
	double rho,	   // Correlation parameter
	int N		   // Sample size
)
{
	const int MAXITERATES=100;

	double n=(double)N;

	if (N<3 || r<-1 || r>1 || rho<-1 || rho>1)
		return NA_REAL;


	if (fabs(r)>=1.0) {
		return 0.0;
	}

	double scale=(n-2.0)/(SQRT2*(n-1.0));
	double logFactor=0.5*(n-1.0)*log(1.0-rho*rho)+0.5*(n-4.0)*log(1.0-r*r);
	logFactor+=(1.5-n)*log(1.0-rho*r)+loggamma(n)-loggamma(n-0.5)-LNGAMMAHALF;

		// Evaluate the hypergometric function F[1/2,1/2;N-1/2,(1+rho*r)/2]	by suming
		//  the series
	double c=n-0.5;
	double x=0.5*(1.0+rho*r);
	double sum=0.0;
	double term=1.0;
	double value;
	int j=1;
	double dj;
	repeat
		dj=(double)j;
		sum+=term;
		value=(double)(2*j-1);
		term*=(0.25*(value*value)/(c+dj-1.0))*(x/dj);
		j++;
	until(sum+term equals sum || j>MAXITERATES);

	return scale*exp(logFactor)*sum;
}



/*
	Distribution function of the correlation coefficient
*/

DISTS_API void pcorrR(
	double *rp,
	double *rhop,
	int *np,
	int *Np,
	double *valuep
)
{
	int N=*Np;
	int i;

	for (i=0;i<N;i++)
		valuep[i]=pcorrelation(rp[i],rhop[i],np[i]);
}


static double grhocorr;
static int gNcorr;

static double fcorrelationP(double r){return fcorrelation(r,grhocorr,gNcorr);}

 double pcorrelation(
	double r,
	double rho,
	int N
)
{
	double P;
	grhocorr=rho;
	gNcorr=N;
	
	if (N<3 || r<-1 || r>1 || rho<-1 || rho>1)
		return NA_REAL;

	P=Integral(-1.0,r,fcorrelationP,3e-8);
	if (P<-0.0001 || P>1.0001)
		return NA_REAL;
	if (P<0) P=0;
	if (P>1) P=1;

	return P;
}

DISTS_API void ucorrR(
	double *rp,
	double *rhop,
	int *np,
	int *Np,
	double *valuep
)
{
	int N=*Np;
	int i;

	for (i=0;i<N;i++)
		valuep[i]=qcorrelation(rp[i],rhop[i],np[i]);
}

 double qcorrelation(
	double r,
	double rho,
	int N
)
{ 	
	if (N<3 || r<-1 || r>1 || rho<-1 || rho>1)
		return NA_REAL;
	return 1.0-pcorrelation(r,rho,N);
}


DISTS_API void qcorrR(
	double *pp,
	double *rhop,
	int *np,
	int *Np,
	double *valuep
)
{
	int N=*Np;
	int i;

	for (i=0;i<N;i++)
		valuep[i]=xcorrelation(pp[i],rhop[i],np[i]);
}


static double gpcorr;

	 // difference between dist fcn and desired p
static double dcorrelationP(double r){return -gpcorr + pcorrelation(r,grhocorr,gNcorr);}

 double xcorrelation(
	double p,
	double rho,
	int N
)
{
	gpcorr=p;
	grhocorr=rho;
	gNcorr=N;
	
	if (N<3 || p<0 || p>1 || rho<-1 || rho>1)
		return NA_REAL;


		// assume arctanh(R) is normal with expectation 0.5*log((1+rho)/(1-rho))
		// and variance 1/(N-3) -- See Johnson and Kotz Vol2, first Ed, 229
	double z=0.5*log((1.0+rho)/(1.0-rho))+qnorm(p,0,1,true,false)/sqrt((double)(N-3));
	z=exp(2.0*z);
	double guess=(z-1.0)/(z+1.0);

	return NewtonRoot(guess,false,dcorrelationP,fcorrelationP,3e-8);
}


/*
	Calculates the correlation coefficient between two n element arrays x[] and y[]
	Following a FORTRAN routine by Allan Miller
*/
static double CorrelationCoefficient(
	double* x,
	double* y,
	int n	
)
{
	double sxx=0.0;
	double syy=0.0;
	double sxy=0.0;
	double xMean=0.0;
	double yMean=0.0;
	double devX;
	double devY;
	double devXu;
	double di;

	for (int i=0;i<n;i++) {
		di=(double)(i+1);
		devX=x[i]-xMean;
		devY=y[i]-yMean;
		xMean+=devX/di;
		yMean+=devY/di;
		sxx+=devX*(devXu=x[i]-xMean);
		syy+=devY*(y[i]-yMean);
		sxy+=devY*devXu;
	}
	return sxy/sqrt(sxx*syy);
} 


/*
	Random correlations with propulation correlation rho
*/
DISTS_API void rcorrR(
	double *rhop,
	int *np,
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
		rcorrelation(valuep,*np,*rhop,N);
	else { // Allow for random values for each element of the arguments
		D=(N/M)+((N%M)?1:0);
		tArray=(double *)S_alloc((long)D,sizeof(double));
		loc=0;
		for (j=0;j<M;j++) {
			rcorrelation(tArray,np[j],rhop[j],D);
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

 void rcorrelation(
	double *randArray,
	long n,
	double rho,
	int N
)
{
	double* x=(double *)S_alloc(n,sizeof(double));
	double* y=(double *)S_alloc(n,sizeof(double));

	if (n<3 || rho<-1 || rho>1) {
		for (int i=0;i<N;i++)
			randArray[i]=NA_REAL;
	}
	else {
		for (int i=0;i<N;i++) {
			rgauss(x,n,0.0,1.0);
				// Construct y with E(y|x)=rho*x, and sd sqrt(1-rho^2)
			rgauss(y,n,0.0,sqrt(1.0-rho*rho));
			for (int j=0;j<n;j++) {
				y[j]+=rho*x[j];
			}

			randArray[i]=CorrelationCoefficient(x,y,n);
		}
	}

}


DISTS_API void scorrR(
	double *rhop,
	int *np,
	int *Np,
	double *meanp,
	double *medianp,
	double *modep,
	double *varp,
	double *thirdp,
	double *fourthp
)
{
	int N=*Np;
	int i;
	double rho;
	int n;
	double m;
	double m2;
	double rho2;
	double rho4;
	double oneMin;
	double oneMin2;
	double oneMin3;
	double oneMin4;

	for (i=0;i<N;i++) {
		rho=rhop[i];
		n=np[i];
		if (n<3 || rho<-1 || rho>1) {
			meanp[i]=NA_REAL;
			medianp[i]=NA_REAL;
			modep[i]=NA_REAL;
			thirdp[i]=NA_REAL;
			fourthp[i]=NA_REAL;
			varp[i]=NA_REAL;
		}
		else {
			grhocorr=rho;
			gNcorr=n;

			m=1.0/(n+6.0);
			m2=m*m;
			rho2=rho*rho;
			rho4=rho2*rho2;
			oneMin=1.0-rho2;
			oneMin2=oneMin*oneMin;
			oneMin3=oneMin2*oneMin;
			oneMin4=oneMin2*oneMin2;

			meanp[i]=rho-0.5*m*rho*oneMin*(1.0+2.25*m*(3.0+rho2)+
				0.375*m2*(121.0+70.0*rho2+25.0*rho4));
			medianp[i]=xcorrelation(0.5,rho,n);
			modep[i]=FindDistributionMode(-1.0,1.0,fcorrelationP);
			thirdp[i]=-m2*rho*oneMin3*(6.0+m*(69.0+88.0*rho2)+
				0.75*m2*(797.0+1691.0*rho2+1560*rho4));
			fourthp[i]=3.0*m2*oneMin4*(1.0+m*(12.0+35.0*rho2)+
				0.25*m2*(436.0+2028*rho2+3025*rho4));
			varp[i]=m*oneMin2*(1.0+0.5*m*(14.0+11.0*rho2)+
				0.5*m2*(98.0+130.0*rho2+75.0*rho4));
		}
	}
}


/**************************************************************************************** 
	Hypergeometric distribution 
	Given a total of N items, n of which are marked, select a sample of size S, and
	 find the probability that 
	 
	fhypergeometric:  that there are x marked items in the sample -- frequency
	phypergeometric:  that there are x or less in the sample -- cdf.
*/


 	// Normal approximation to the hypergeometric distribution function due to Peizer
	// See Ling, R.F. and Pratt, J.W. (1984) The accuracy of Peizer approximations
	//  to the hypergeometric distribution, with comparisons to some other 
	//  approximations. JASA 79-385. 49-60.
double PeizerHypergeometric(
	int x,	  // Number of marked items in sample
	int S,	  // Sample size
	int n,	  // Total number of marked items
	int N	  // Total number of items
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
	double a, 				// Sample size
	double m,      	// Total number of marked items
	double N       		// Total number of items
)
{

	hyperType variety;



	if (0.0<a && 0.0<N && 0.0<m &&  isint(a) && isint(N) && isint(m)) {
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
	if (0.0<a && 0.0<N && 0.0<m &&  ! isint(a) && ! isint(m) && a+m-1.0<N &&
		 		floor(a) equals floor(m)) {
		variety=IB;		// Specified 1.0<N to avoid problems with small parameters
	}
	else
	if (a<0.0 && N<m+a-1.0 && 0.0<m && isint(m))	{ //Kemp&Kemp use b<0 && b!=-1, Ben Bolker mod
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
	double a, 				// Sample size
	double m,      	// Total number of marked items
	double N,       		// Total number of items
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
		case IB:		// Specified 1.0<N to avoid problems with small parameters
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
	double *ap, 				// Sample size
	double *mp,      	// Total number of marked items
	double *Np,       		// Total number of items
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
		case IB:		// Specified 1.0<N to avoid problems with small parameters
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
	int x,		// Number of marked items in sample
	int S, 				// Sample size
	int n,      	// Total number of marked items
	int N       		// Total number of items
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

		x  y | n
		?  ? | N-n
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
	int x,		// Number of marked items in sample
	int a, 				// Sample size
	int n,      	// Total number of marked items
	int N       		// Total number of items
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
		//  hypergeometric function -- i.e. coefficients of x^i in the expansion.
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
	int x,		// Number of marked items in sample
	int a, 				// Sample size
	int n,      	// Total number of marked items
	int N       		// Total number of items
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

	// returns smallest x such that  p<=Pr(X<=x|a,n,N) 
int  xhypergeometric(
	double p,		// cumulative probability
	int a, 				// Sample size
	int n,      	// Total number of marked items
	int N       		// Total number of items
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
	int n,	  // number of samples
	int a, 				// Sample size
	int m,      	// Total number of marked items
	int N       		// Total number of items
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
		return minm(P,1.0);	 // Occasional numerical error
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


	// returns smallest x such that  p<=Pr(X<=x|a,n,N) 
int  xgenhypergeometric(
	double p,		// cumulative probability
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
	double lowX,   // Assumed to be less than highX
	double highX,
	double (*function)(double x),
	double Tol
)
{

	double h=highX-lowX;   //  range of integration

	/*const double Tol=3e-8;*/  // result to this accuracy -- can't really to better

	const int maxiterate=16;	  // Maximum interates allowed -- fcn evaluated at most 2^16 pts
	double A[maxiterate][maxiterate]; // the triangular array of Romberg iterates
									  // stop when two diagonal values agree to MRatioTol
	
	double delta=h;		 // Initial step
	double value=0.0;		 // Converged value

	A[0][0]=(h/2.0)*((*function)(lowX)+(*function)(highX)); 

	double twoPower=1.0;
	int k=0;
	int n=1;
	repeat
		k++;
		delta*=0.5;	 // Evaluate (*funciton)() at half the previous interval
		if (k>1) {
			n*=2;		 // Number of new evaluations
		}

		twoPower*=2.0;
		double sum=0.0;
		double z=highX-delta;	// Start with this value
		int m=n;
		while (m--) {
			double value=(*function)(z);

			sum+=value;
			z-=2.0*delta;	 // Every other ordinate is a new one
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
	bool useLog,	// When true, will iterate on z=log(x), keeping x positive
	double (*function)(double x),
	double (*derivative)(double x),
	double TOLN
)
{
 /*	const double TOLN=3e-8;*/
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
			h=scale*0.5*fcn/(x*deriv+DBL_EPSILON*fabs(fcn));  // no small divisors allowed
		}													  // This is likely of no use
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
		if (ho<=fabs(h))	{
			scale/=2.0;
			z+=h;		 // Restore value
			more=true;
			continue;  // retry
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
#define zNew  ((zSeed=36969*(zSeed&65535)+(zSeed>>16))<<16)
#define wNew  ((wSeed=18000*(wSeed&65535)+(wSeed>>16))&65535)
#define IUNIFORM  (zNew+wNew)
#define UNIFORM   ((zNew+wNew)*RANDCONST)
#define setseed(A,B) zSeed=A;wSeed=B;
#define getseed(A,B) A=zSeed;B=wSeed;

/* The ziggurat method for RNOR and REXP

Modified very slightly to remove some of the awkwardness of the original C code. In partcular, 
float has been replaced by double, some implicit coercions made explicit,  and the 
initialization of jsr moved into zigset to guarantee same sequence from same seed.
REW Mar 01
REW (Mar 2001) cleaned up the code. Made RNOR and REXP inline instead of defines, etc.

Combine the code below with the main program in which you want
normal or exponential variates.   Then use of RNOR in any expression
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

static ULONG	jz,
				jsr; // moved initialization
static ULONG    jcong;
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

static int		iz;		 // was ULONG
static long		kn[128]; // was ULONG
static ULONG	ke[256];
static double	wn[128],
				fn[128], 
				we[256],
				fe[256];

inline double RNOR(void) {
	double xx;
	hz=SHR3; 
	iz=hz&127; 
	xx=(abs(hz)<kn[iz])? double(hz*wn[iz]) : nfix(); // Original had sign conflict
	return xx;										 // for abs(hz)<kn[]
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
#define REXP (jz=SHR3, iz=jz&255, (    jz <ke[iz])? double(jz*we[iz]) : efix())

#endif

/* nfix() generates variates from the residue when rejection in RNOR occurs. */
double nfix(void) {	
	const double r = 3.442619855899; 	/* The starting of the right tail */	
	double  x;
	double  y;

	repeat		
		x=hz*wn[iz];
		if(iz==0){	/* iz==0, handle the base strip */
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
			return (7.69711-log(UNI));		/* iz==0 */
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
	const double	m1 = 2147483648.0; 
	const double	m2 = 4294967296.0;
	double	dn=3.442619855899;
	double 	tn=dn;
	double 	vn=9.91256303526217e-3; 
	double 	q;      
	double	de=7.697117470131487; 
	double 	te=de; 
	double 	ve=3.949659822581572e-3;
	int i;	
		  
	jsr=123456789; // moved to here to allow seed to control generation
	jsr^=jsrseed;
	setseed(jsrseed,jsrseed)  /* Added to support Leong et.al.*/
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
#define zNew  ((zSeed=36969*(zSeed&65535)+(zSeed>>16))<<16)
#define wNew  ((wSeed=18000*(wSeed&65535)+(wSeed>>16))&65535)
#define IUNIFORM  (zNew+wNew)
#define UNIFORM   ((zNew+wNew)*RANDCONST)
#define setseed(A,B) zSeed=A;wSeed=B;
#define getseed(A,B) A=zSeed;B=wSeed;
*/

static int nSeed=1020;
static ULONG Q[1020];   // using Q[endQ] to hold variable.
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
	for (i=0;i<endQ;i++)
		Q[i]=IUNIFORM;

}


ULONG MWC1019(void){
	ULONG t;
	int	i = endQ-1; 

	t = 147669672L*Q[i] + Q[endQ]; 
	Q[endQ] = (t>>32);
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

before calling MWC1019( ) in your main.   You might use a good RNG such as
KISS
to fill the Q array.

The period of MWC1029 exceeds 10^9824, making it  billions and billions ...
and billions times as long as the highly touted longest-period RNG,
the Mersenne twister.   It is also several times as fast and takes a few
lines rather than
 several pages of code.   (This is not to say that the Mersenne twister is
not
a good RNG; it is.  I just do not equate complexity of code  with
randomness.   It is the
complexity of the underlying randomness that counts.)

As for randomness, it passes all tests in The Diehard Battery of Tests of
Randomness
          http://stat.fsu.edu/pub/diehard
as well as three new tough tests I have developed with the apparent
property
that a RNG that passes  tuftsts.c will pass all the tests in Diehard.

MWC1019 has the property that every possible sequence of 1018 successive
32-bit integers will appear somewhere in the full period,  for those
concerned
with the "equi-distribution" in dimensions 2,3,...1016,1017,1018.

I welcome comments  on timings  or otherwise.

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
	else if (ziggInitialized==false) {  // To always insure initialization
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

static double correc(int, int);

void nscor2(
	double *s, 
	int *n, 
	int *n2
)
{

/*     algorithm as 177.3, applied statistics, v.31, 161-165, 1982.

     calculates approximate expected values of normal order statistics.
     claimed accuracy is 0.0001, though usually accurate to 5-6 dec.

 ***  N.B. This routine was NOT in double precision All constants were f ***

     Arguments:

     s(n2)   = output, the first n2 expected values.
     n	     = input, the sample size.
     n2	     = input, the number of order statistics required; must
		      be <= n/2.

	 ier removed REW Mar 2001

     ier     = output, error indicator
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
	dl1[4] = { .112063,.12177,  .239299, .215159 },
	dl2[4] = { .080122,.111348,-.211867,-.115049 },
	gam[4] = { .474798,.469051, .208597, .259784 },
	lam[4] = { .282765,.304856,.407708,.414093 };
    const double bb = -.283833;
    const double d  = -.106136;
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

/*	calculate normal tail areas for first 3 order statistics. */

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

/*	calculate normal areas for other cases. */

	for (i = 4-1; i < *n2; ++i) {
	    ai = (double) i+1;
	    e1 = (ai - eps[3]) / (an + gam[3]);
	    e2 = pow((double) e1, (double) lam[3] + bb / (ai + d));
	    s[i] = e1 + e2 * (dl1[3] + e2 * dl2[3]) / an - correc(i+1, *n);
	}
    }
/*	convert tail areas to normal deviates. */

    for (i = 0; i < *n2; ++i)
	s[i] = - qnorm(s[i], 0., 1., 1, 0);

    return;
} /* nscor2 */


static double correc(int i, int n)
{
/*	calculates correction for tail area of the i-th largest of n
	order statistics. */

    const double
	c1[7] = { 9.5,28.7,1.9,0.,-7.,-6.2,-1.6 },
	c2[7] = { -6195.,-9569.,-6728.,-17614.,-8278.,-3570., 1075. },
	c3[7] = { 93380.,175160.,410400.,2157600.,2.376e6,
		  2.065e6,2.065e6 };
    const double mic = 1e-6;
    const double c14 = 1.9e-5;

    double an;

    if (i * n == 4)		return c14;
    if (i < 1 || i > 7)		return 0;
    if (i != 4 && n > 20)	return 0;
    if (i == 4 && n > 40)	return 0;
    /* else : */
    an = (double) n;
    an = 1. / (an * an);
    i--;
    return((c1[i] + an * (c2[i] + an * c3[i])) * mic);
} /* correc */

