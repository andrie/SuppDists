#define DISTS_API extern "C"

enum hyperType {
	classic,
	IAi,
	IAii,
	IB,
	IIA,
	IIB,
	IIIA,
	IIIB,
	IV,
	noType
};


	/* Prototypes for dist functions */
	// Functions marked with EX may throw exceptions
DISTS_API void ziggR(double *randomVector,int *Np,bool *type,bool *initilizep,ULONG *seedp);
DISTS_API void MWC1019R (double *randomVector,int *Np,bool *initializep,ULONG *seedp);
#ifdef CANTUSE
DISTS_API void user_unif_init(ULONG seed);
DISTS_API int *user_unif_nseed(void);
DISTS_API int *user_unif_seedloc(void);
DISTS_API double *user_unif_rand(void);
#endif

double   gammal(double);

double FindDistributionStatistic(double lowX,double highX,double (*function)(double x));
double FindDistributionMode(double lowX,double highX,double (*function)(double x));



void nscor2(double *s, int *n, int *n2);
DISTS_API void normOrdR(double *sp,int *np,int *n2p);




	// normal
void 	 rgauss(double* normArray,int n,double mean, double sd);

  //	Inverse Gaussian -- Wald
double 	finvGauss(double x,double mu,double lambda);
double 	pinvGauss(double x,double mu,double lambda);
double 	qinvGauss(double x,double mu,double lambda);
double 	xinvGauss(double p,double mu,double lambda);
void 	rinvGauss(double* normArray,int n,double mu,double lambda);

DISTS_API void dinvGaussR(double *xp,double *nup,double *lambdap,int *Np,double *valuep);
DISTS_API void pinvGaussR(double *xp,double *nup,double *lambdap,int *Np,double *valuep);
DISTS_API void uinvGaussR(double *xp,double *nup,double *lambdap,int *Np,double *valuep);
DISTS_API void qinvGaussR(double *xp,double *nup,double *lambdap,int *Np,double *valuep);
DISTS_API void rinvGaussR(double *nup,double *lambdap,int *Np,int *Mp,double *valuep);





 


	// Kruskal Wallace
  double  pKruskal_Wallis(double H,int c,int n,double U,bool doNormalScore);  //EX
  double  qKruskal_Wallis(double H,int c,int n,double U,bool doNormalScore);  //EX
  double  xKruskal_Wallis(double P,int c,int n,double U,bool doNormalScore); //EX
double varKruskal_Wallis(double N,double C,double U);
double varNormalScores(double N,double C,double U);
  void	rKruskal_Wallis(double* randArray,int N,int c,int n,double U,
			bool doNormalScore);					   //EX
  double  fKruskal_Wallis(double H,int c,int n,double U,bool doNormalScore);  //EX
  void	sKruskal_Wallis(int c,int n,double U,bool doNormalScore,double *mode,
			double *third,double *fourth);  // EX
DISTS_API void pKruskalWallisR(double *Hp,int *cp,int *np,double *Up,int *doNormalScorep,int *Np,double *valuep);
DISTS_API void uKruskalWallisR(double *Hp,int *cp,int *np,double *Up,int *doNormalScorep,int *Np,double *valuep);
DISTS_API void qKruskalWallisR(double *Hp,int *cp,int *np,double *Up,int *doNormalScorep,int *Np,double *valuep);
DISTS_API void dKruskalWallisR(double *Hp,int *cp,int *np,double *Up,int *doNormalScorep,int *Np,double *valuep);
DISTS_API void rKruskalWallisR(double *randArrayp,int *Np,int *Mp,int *cp,int *np,double *Up,bool *doNormalScorep);
DISTS_API void sKruskalWallisR(int *cp,int *np,double *Up,int *doNormalScorep,int *Np,double *varp,double *modep,double *thirdp,double *fourthp);

	// Runs test
  double  pruns(int mi,int ni,int r);
  double  qruns(int m,int n,int r);
  int  xruns(double pr,int m,int n);

	// Kendall's Tau

const int MAXKENDALEXACT=12;
  double  pkendall(int ni,double tau);   // EX
  double  qkendall(int n,double tau);	   // EX
  double	xkendall(double pr,int ni);	   // EX
  void	rkendall(double* randArray,int N,int ni);	// EX
  double  fkendall(int ni,double tau);   // EX
double  fourthkendall(int ni);		   // EX

DISTS_API void pKendallR(int *nip,double *taup,int *Np,double *valuep);
DISTS_API void dKendallR(int *nip,double *taup,int *Np,double *valuep);
DISTS_API void qKendallR(int *pp,double *taup,int *Np,double *valuep);
DISTS_API void uKendallR(int *nip,double *taup,int *Np,double *valuep);
DISTS_API void rKendallR(int *nip,int *Np,int *Mp,double *valuep);
DISTS_API void fourthKendallR(int *nip,int *Np,double *valuep);

	// Friedman	This is here, because FriedmanGlobal needs it
struct FriedmanStrc {
   int* S;		  
   int	nS;
   double *qdist;
};

struct FriedmanGlobal {
	int r;
	int n;
	FriedmanStrc *theDist;
};

bool DoExactFriedman(int r,int n,bool doRho);
void ClearFriedmanGlobal(bool deleteAll);

  double  pfrie(double X,int r,int n,bool  doRho);   //EX
  double  qfrie(double X,int r,int n,bool  doRho);   //EX
  double  xfrie(double P,int r,int n,bool  doRho);  //EX
double  medianfrie(int r,int n);				//EX
double  modefrie(int r,int n);					// EX
  void	rfrie(double* randArray,int N,int r,int n,bool doRho);	// EX
  double	ffrie(double X,int r,int n,bool doRho);	   // EX

DISTS_API void pFriedmanR(double *Xp,int *rp,int *np,int *Np,bool *doRhop,double *valuep);
DISTS_API void uFriedmanR(double *Xp,int *rp,int *np,int *Np,bool *doRhop,double *valuep);
DISTS_API void dFriedmanR(double *Xp,int *rp,int *np,int *Np,bool *doRhop,double *valuep);
DISTS_API void qFriedmanR(double *pp,int *rp,int *np,int *Np,bool *doRhop,double *valuep);
DISTS_API void rFriedmanR(int *rp,int *np,bool *doRhop,int *Np,int *Mp,double *valuep);
DISTS_API void sFriedmanR(int *rp,int *np,bool *rhop,int *Np,double *meanp,double *medianp,double *modep,double *variancep,	double *thirdp,double *fourthp);  

	// Maximum F ratios
  double pmaxfratio(double F,int df,int N);
  double qmaxfratio(double F,int df,int N);
  double xmaxfratio(double p,int df,int N);
  double fmaxfratio(double F,int df,int N);
void smaxFratio(int df,int N,double *mean,double *median,double *mode,double *variance,double *third,double *fourth);

DISTS_API void pmaxFratioR(double *Fp,int *dfp,int *np,int *Np,double *valuep);
DISTS_API void umaxFratioR(double *Fp,int *dfp,int *np,int *Np,double *valuep);
DISTS_API void dmaxFratioR(double *Fp,int *dfp,int *np,int *Np,double *valuep);
DISTS_API void qmaxFratioR(double *Fp,int *dfp,int *np,int *Np,double *valuep);
DISTS_API void rmaxFratioR(int *dfp,int *np,int *Np,int *Mp,double *valuep);
DISTS_API void smaxFratioR(int *dfp,int *np,int *Np,double *mean,double *median,double *mode,double *variance,double *third,double *fourth);


	// Johnson curves 

struct JohnsonInput {
	double xn; // distribution value corresponding to zn
	double xm; // distribution value corresponding to zn/2
	double x0; // distribution value corresponding to 0
	double xk; // distribution value corresponding to -zn/2
	double xp; // distribution value corresponding to -zn
};

enum JohnsonType {
	SN,
	SL,
	SU,
	SB
};

struct JohnsonParms {
	double gamma;
	double delta;
	double xi;
	double lambda;
	JohnsonType type;
};

struct JohnsonMoments {
	double mean;
	double sd;
	double sqrtB1;
	double B2;
};

  JohnsonParms JohnsonFit(JohnsonInput input);
  JohnsonParms JohnsonMomentFit(JohnsonMoments moments);

  double pjohnson(double x,JohnsonParms parms);
  double qjohnson(double x,JohnsonParms parms);
  double xjohnson(double p,JohnsonParms parms);
double xzjohnson(double p,JohnsonParms parms);
  double fjohnson(double x,JohnsonParms parms);
double fpjohnson(double x,JohnsonParms parms);
  void   rjohnson(double* johnsonArray,int n,JohnsonParms parms);
void sJohnson(JohnsonParms parms,double *meanp,double *medianp,double *modep,
			  double*variancep,double *thirdp,double *fourthp);

DISTS_API void JohnsonFitR(double *xnp,double *xmp,double *x0p,double *xkp,double *xpp,
	double *gammap,double *deltap,double *xip,double *lambdap,int *typep);
DISTS_API void JohnsonMomentFitR(double *meanp,double *sdp,double *sqrtB1p,double *B2p,
	double *gammap,double *deltap,double *xip,double *lambdap,int *typep);
DISTS_API void pJohnsonR(double *xp,double *gammap,double *deltap,double *xip,
	double *lambdap,int *typep,int *Np,double *valuep);
DISTS_API void uJohnsonR(double *xp,double *gammap,double *deltap,double *xip,
	double *lambdap,int *typep,int *Np,double *valuep);
DISTS_API void dJohnsonR(double *xp,double *gammap,double *deltap,double *xip,
	double *lambdap,int *typep,int *Np,double *valuep);
DISTS_API void qJohnsonR(double *pp,double *gammap,double *deltap,double *xip,
	double *lambdap,int *typep,int *Np,double *valuep);
DISTS_API void rJohnsonR(double *gammap,double *deltap,double *xip,
	double *lambdap,int *typep,int *Np,int *Mp,double *valuep);
DISTS_API void sJohnsonR(double *gammap,double *deltap,double *xip,double *lambdap,
	int *typep,int *Np,double *meanp,double *medianp,double *modep,double *variancep,
	double *thirdp,double *fourthp);
//DISTS_API momentsR(double *data,int *Np,double *meanp,double *variancep,
//	double *thirdp,double *fourthp);

// Correlation coefficient



  double fpcorrelation(double r,double rho,int N);
  double fcorrelation(double r,double rho,int N);
  double pcorrelation(double r,double rho,int N);
  double qcorrelation(double r,double rho,int N);
  double tcorrelation(double r,double rho,int N);
  double xcorrelation(double r,double rho,int N);
  void   rcorrelation(double *randArray,long n,double rho,int N);

DISTS_API void pcorrR(double *rp,double *rhop,int *np,int *Np,double *valuep);
DISTS_API void ucorrR(double *rp,double *rhop,int *np,int *Np,double *valuep);
DISTS_API void qcorrR(double *rp,double *rhop,int *np,int *Np,double *valuep);
DISTS_API void dcorrR(double *rp,double *rhop,int *np,int *Np,double *valuep);
DISTS_API void rcorrR(double *rhop,int *np,int *Np,int *Mp, double *valuep);
DISTS_API void scorrR(double *rhop,int *np,int *Np,double *meanp,double *medianp,double *modep,double *varp,double *thirdp,double *fourthp);

// Hypergeometric

double qhypergeometric(int x,int a,int n,int N);
double phypergeometric(int x,int a,int n,int N);
double fhypergeometric(int x,int a,int n,int N);
int xhypergeometric(double p,int a,int n,int N);
void rhypergeometric(double* randArray,int n,int a,int m,int N);

  double qgenhypergeometric(int x,double a,double m,double N,hyperType variety);
  double pgenhypergeometric(int x,double a,double m,double N,hyperType variety);
  double fgenhypergeometric(int x,double a,double m,double N,hyperType variety);
  int xgenhypergeometric(double p,double a,double m,double N,hyperType variety);
  void rgenhypergeometric(double* randArray,int K,double a,double m,double N,hyperType variety);

void sghyper(double a,double m,double N,double *mean,double *median,double *mode,double *variance,double *third,double *fourth,hyperType variety);

DISTS_API void tghyperR(double *ap,double *mp,double *Np,char **aString);
DISTS_API void pghyperR(int *kp,double *ap,double *np,double *Np,int *Mp,double *valuep);
DISTS_API void ughyperR(int *kp,double *ap,double *np,double *Np,int *Mp,double *valuep);
DISTS_API void qghyperR(double *pp,double *ap,double *np,double *Np,int *Mp,double *valuep);
DISTS_API void dghyperR(int *kp,double *ap,double *np,double *Np,int *Mp,double *valuep);
DISTS_API void rghyperR(double *ap,double *np,double *Np,int *Mp,int *Kp,double *valuep);
DISTS_API void sghyperR(double *ap,double *mp,double *Np,int *Mp,double *meanp,double *medianp,double *modep,double *variancep,double *thirdp,double *fourthp);

// Helper functions
double  Integral(double lowX,double highX,double (*function)(double x));
double  NewtonRoot(double guess,bool useLog,double (*function)(double x),double (*derivative)(double x));
  double  GaussianHypergometricFcn(double a,double b,double c,double x);

double  Integral(double lowX,double highX,double (*function)(double x),double Tol);
double  NewtonRoot(double guess,bool useLog,double (*function)(double x),double (*derivative)(double x),double TOLN);
double  GaussianHypergometricFcn(double a,double b,double c,double x);

 
