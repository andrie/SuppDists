/********************************************************************************
 Fits a Johnson curve from percentiles
See Wheeler 1980 Biometrika
The input is in the struct JohnsonInput
The output is in the struct JohnsonParms 
*/

/* 
| Rotates the 3 element row v into matrix using Givens rotations
| matrix is 3 by 3 and should be zeroed before the start
| The weights are on the diagonal
*/

static void Rotate3(
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
      if (d equals 0.0) skip=true; /* to avoid 0/0, but d can't be 0 */
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
  const int N=6;  // Number of moments
  
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
      if (negativeB1) {
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
// of powers

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
    return parms; // Outside the upper limit.
  }
  
  // Is it a normal ?
  if (fabs(sqrtB1)<=TOLJ && fabs(B2-3.0)<=TOLJ) {
    parms.type=SN;
    parms.gamma=0.0; // Hill et.al. let gamma delta carry the burden
    parms.delta=1.0;
    parms.lambda=sd;
    parms.xi=mean;
    return parms;
  }
  // Test for position relative to log normal line
  // The log normal line is defined by a pair of equations:
  // B1=(w-1)(w+2)^2, and B2=w^4+2w^3+3w^2-3.
  
  // First solve the cubic (w-1)(w+2)^2=B1 then estimate B2 and check the
  // estimate for closeness to the actual B2
  
  // This method of solution appears in the statlib copy. The original
  // publication used a, presumably, less stable method.
  double x=0.5*B1+1.0;
  double y=sqrt(B1+0.25*B1*B1);
  double u=pow(x+y,1.0/3.0);
  // Evaluate w and test B2 equals B2est=w^4+2w^3+3w^2-3
  double w=u+1.0/u-1.0;
  double B2est=w*w*(3.0+w*(2.0+w))-3.0;
  
  B2=(B2<0.0)?B2est:B2; // Log fit for neg B2
  
  double test=B2est-B2;
  
  // Is equation close to B2? If so do a log fit
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
    if (fvalue>fmode) {
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
  
  
  switch (parms.type) {
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
vol 2, 1st edition. Uses Gaussian hypergeometric function, defined on
page 17 (1.104) of Johnson, Kotz, and Kemp, Univariate Discrete Distributions,
2nd Ed. The derivitive is given on page 18.
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
    double r,  // The correlation coefficient
    double rho,  // Correlation parameter
    int N   // Sample size
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
  
  // Evaluate the hypergometric function F[1/2,1/2;N-1/2,(1+rho*r)/2] by suming
  // the series
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


