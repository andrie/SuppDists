#include <math.h>
#include <float.h>
//# Hornik replaced <new.h> with <new> and then inserted
//# std::set_new_handler(freeStoreException); on line 2283 (March 2008)
#include <new>
#include <R.h>
#include <Rmath.h>

#include "dists.h"
#include "datatabs.h"


// FRIEDMANS ******************************************************************

// This global is set when an exact distribution is active.
// It is reset when either r or n changes.
FriedmanGlobal *FriedmanCurrentGlobal=0L;

/* Defined in DISTFCNS.H: 
struct FriedmanStrc {
  int* S;   // Sum of squares for each type
  int nS;
  double *qdist; // Cumulative frequencies summed down from the upper tail.
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
of ranks for ith of r treatments. There are n blocks and av=n(r+1)/2. The 
max value of S is n^2r(r^2-1)/12. Uses Kendall-Babington Smith algorithm
(AMS 1939, 10, 275-287). See Hajek & Sidak p144 for a formal statement. 

Spearman's rho is (Xr/(r-1))-1 for n=2, or rho=-1+2S/M, where M=r(r^2-1)/3

The Kendall-Babington Smith algorithm:
(1) Let b be a vector of base ranks (these are the centered ranks (ri-av))
(1-r)/2,(3-r)/2,...,(r-1)/2  for odd r
(1-r),(3-r),...,(r-1)   for even r (to keep the values integral)
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
they are indeed of the same type. As an example, here are two arrays of different
type which pass the absolute order test -- they were generated by the algorithm:
-12,-12,-6,-2, 4, 4,10,14;
-12,-12, 6, 2,-4,-4,10,14.
In calculatin the Spearman rho, one only needs to check S, since new types
will not be generated from the types for n=2.
The algorithm has been run for a number of different values of r and n -- see
DoExactFriedman(), and stored as arrays of type FriedmanValues. A FriedmanStrc
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
  int *R; // Column sums of centered ranks: sorted so |R[i]|<=|R[j]| for i<j  
  int S;    // Sum of sqares of R[i]
  double frequency;  // This is a limiting factor for the algorithm -- must not overflow 
  // hence it is double. There is the possiblity that small frequencies
  // may not be properly added to larger frequencies, but I have been
  // unable to detect any errors.
};


// Finds next permutation given key following Fike (1975)
//  The computer Journal 18-1 21-22.
// This generates r! sequences such that 0<=key[i]<=i+1 for
//  i=0...r-1
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
// Sorts by increasing value 
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
// Sorts orderC pairs by increasing value 
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
// the same or negatives of each other
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
  // must check to see that the signs match. So will sort them
  // in increasing order and compare them. If this doesn't give
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
    error("\nInernal error in InsertTypeInList()"); // This will abort, but leave memory dirty
    return nWorkingTypes;
  }
  // Add the aType pointer to the workingTypes list
  workingTypes[nWorkingTypes++]=aType;
  return nWorkingTypes;
}



// The base contains the centered ranks.
// The r! permutations of the base are added to the current type array, R, to create
// working types. 
// If the resulting working type is a new type it is added to the workingTypes list.
// If the resulting working type is not new, it is deleted and the currentType frequency
// is added to the existing workingType freqency.

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
    key[i-1]=i;  // Puts 1...r-1 in key (key[r-1] is not used by UpdateTheKey())
  }
  
  
  int *temp=(int *)S_alloc((long)r,sizeof(int));
  
  repeat 
    for (int i=0;i<r;i++) {
      temp[i]=base[i];  
    }
    PermuteFike(temp,key,r);  // Permutes temp according to key
  
  Ftype* aType=new Ftype;
  aType->R=new int[r];  // This pointer is either added to the workingType list
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
  solongas(UpdateTheKey(0,r,key)); // Cycles key through all possible permutations
  
  
  return nWorkingTypes;
}




static FriedmanStrc* MakeFriedmanStrc(
    int nCurrentTypes,
    Ftype** currentTypes,
    double totalFrequency // Double because totalFrequency is (r!)^(n-1) which gets big.
)
{
  FriedmanStrc* theStrc=new FriedmanStrc;
  int i;
  // Set up orderS for sorting
  int (*orderS)[2];
  orderS=new int[nCurrentTypes][2];
  for (i=0;i<nCurrentTypes;i++) {
    orderS[i][0]=i;   // Set the index
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
    int j=orderS[i][0];  // Index of ith largest S in currentTypes
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
// The same could be done for small n's

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
    case 3: theValues=(FriedmanValues*)FriedmanData3_2; break;
    case 4: theValues=(FriedmanValues*)FriedmanData4_2; break;
    case 5: theValues=(FriedmanValues*)FriedmanData5_2; break;
    case 6: theValues=(FriedmanValues*)FriedmanData6_2; break;
    case 7: theValues=(FriedmanValues*)FriedmanData7_2; break;
    case 8: theValues=(FriedmanValues*)FriedmanData8_2; break;
    case 9: theValues=(FriedmanValues*)FriedmanData9_2; break;
    case 10: theValues=(FriedmanValues*)FriedmanData10_2; break;
    case 11: theValues=(FriedmanValues*)FriedmanData11_2; break;
    // Tried r==12, but it took 18 hrs and then bombed??
    default: theValues=0L;
    isInTables=false;
    break;
    }
  }
  else
    if (r equals 3) {
      switch (n) {
      case 3: theValues=(FriedmanValues*)FriedmanData3_3; break;
      case 4: theValues=(FriedmanValues*)FriedmanData3_4; break;
      case 5: theValues=(FriedmanValues*)FriedmanData3_5; break;
      case 6: theValues=(FriedmanValues*)FriedmanData3_6; break;
      case 7: theValues=(FriedmanValues*)FriedmanData3_7; break;
      case 8: theValues=(FriedmanValues*)FriedmanData3_8; break;
      case 9: theValues=(FriedmanValues*)FriedmanData3_9; break;
      case 10: theValues=(FriedmanValues*)FriedmanData3_10; break;
      case 11: theValues=(FriedmanValues*)FriedmanData3_11; break;
      case 12: theValues=(FriedmanValues*)FriedmanData3_12; break;
      case 13: theValues=(FriedmanValues*)FriedmanData3_13; break;
      case 14: theValues=(FriedmanValues*)FriedmanData3_14; break;
      case 15: theValues=(FriedmanValues*)FriedmanData3_15; break;
      case 16: theValues=(FriedmanValues*)FriedmanData3_16; break;
      case 17: theValues=(FriedmanValues*)FriedmanData3_17; break;
      case 18: theValues=(FriedmanValues*)FriedmanData3_18; break;
      case 19: theValues=(FriedmanValues*)FriedmanData3_19; break;
      case 20: theValues=(FriedmanValues*)FriedmanData3_20; break;
      case 21: theValues=(FriedmanValues*)FriedmanData3_21; break;
      case 22: theValues=(FriedmanValues*)FriedmanData3_22; break;
      case 23: theValues=(FriedmanValues*)FriedmanData3_23; break;
      case 24: theValues=(FriedmanValues*)FriedmanData3_24; break;
      case 25: theValues=(FriedmanValues*)FriedmanData3_25; break;
      case 26: theValues=(FriedmanValues*)FriedmanData3_26; break;
      case 27: theValues=(FriedmanValues*)FriedmanData3_27; break;
      case 28: theValues=(FriedmanValues*)FriedmanData3_28; break;
      case 29: theValues=(FriedmanValues*)FriedmanData3_29; break;
      case 30: theValues=(FriedmanValues*)FriedmanData3_30; break;
      default: theValues=0L;
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
        default: theValues=0L;
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
          default: theValues=0L;
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
  //  (1-r)/2,(3-r)/2,...,(r-1)/2  for odd r
  //  (1-r),(3-r),...,(r-1)   for even r
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
    // frequencies or inserting new types as needed.
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
  if (doRho) { // Spearman's rho -- see limit comment below
    return (1<r && r<12);
  }
  else {
    // These limits can be made larger, but run time becomes excessive
    // Note: if the limits are changed, MAXTYPE may need changing, and 
    // frequency in the Ftype struct's may become too large -- be careful
    switch (r) {
    case 2: return (n<=100); break;
    case 3: return (n<=30); break;
    case 4: return (n<=15); break;
    case 5: return (n<=8); break; // changed from 10 to 8
    // Can't reasonably do larger cases because of array and time limits
    default: return false; break;
    }
  }
  return false; 
}


/* This is needed because not all even integers are in S[]
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
    bool doRho // When true s is Spearman's rho
)
{
  if (DoExactFriedman(r,n,doRho)) {
    if (! FriedmanCurrentGlobal || FriedmanCurrentGlobal->r!=r ||
        FriedmanCurrentGlobal->n!=n) {
      if (FriedmanCurrentGlobal) {
        ClearFriedmanGlobal(false);  // Delete current S and theDist
      }
      else { // Allocate a new global
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
      SS*=4;   // For even r, the ranks were multiplied by 2 to keep them integer.
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
      ClearFriedmanGlobal(true); // Delete and zero the current global
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
    bool doRho // When true s is Spearman's rho
)
{
  if (DoExactFriedman(r,n,doRho)) {
    if (! FriedmanCurrentGlobal || FriedmanCurrentGlobal->r!=r ||
        FriedmanCurrentGlobal->n!=n) {
      if (FriedmanCurrentGlobal) {
        ClearFriedmanGlobal(false);  // Delete current S and theDist
      }
      else { // Allocate a new global
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
      SS*=4;   // For even r, the ranks were multiplied by 2 to keep them integer.
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
      ClearFriedmanGlobal(true); // Delete and zero the current global
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

double pfrie(
    double X,
    int r,
    int n,
    bool doRho  // When true, X is Spearman's rho
)
{
  double Q;
  
  if (doRho)
    n=2;
  if (r<3 || n<2)
    return NA_REAL;
  
  double M=(double)(n*n*r*(r*r-1))/12.0; // Max value of S
  double S;
  if (doRho) {
    S=(M/2.0)*(1.0+X); // X is rho, so rho to S
  }
  else {
    S=(X*(double)(n*r*(r+1)))/12.0; // X to S
  }
  
  if (S>M || S<0)
    return NA_REAL;
  
  long iS=(long)ceil(S); // Round up to even, because S values must be even
  iS=2*(iS/2);
  S=maxm(1L,iS);
  
  
  if (CheckFriedmanExactQ(r,n,X,&Q,true,doRho)) {
    return 1.0-Q;  // Lower tail including X exactly
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

double qfrie(
    double X,
    int r,
    int n,
    bool doRho // When true X is Spearman's rho
)
{
  double Q;
  
  if (doRho)
    n=2;
  if (r<3 || n<2)
    return NA_REAL;
  
  
  double M=(double)(n*n*r*(r*r-1))/12.0; // Max value of S
  double S;
  if (doRho) {
    S=(M/2.0)*(1.0+X); // X is rho, so rho to S
  }
  else {
    S=(X*(double)(n*r*(r+1)))/12.0; // X to S
  }
  
  if (S>M || S<0)
    return NA_REAL;
  
  
  // Changed to conform to R usage
  //  S-=2.0; // Because upper tail, including X, is wanted
  long iS=(long)floor(S); // Round down to even, because S values must be even
  iS=2*(iS/2);
  S=iS;
  S=maxm(1L,iS);
  
  
  
  // if (CheckFriedmanExactQ(r,n,X,&Q,false,doRho)) {
  //  return Q; // Upper tail area including X exactly
  // Changed to conform to R usage
  if (CheckFriedmanExactQ(r,n,X,&Q,true,doRho)) {
    return Q; // Upper tail area excluding X exactly
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
  
  double M=(double)(n*n*r*(r*r-1))/12.0; // Max value of S
  double S;
  if (doRho) {
    S=(M/2.0)*(1.0+X); // X is rho, so rho to S
  }
  else {
    S=(X*(double)(n*r*(r+1)))/12.0; // X to S
  }
  
  if (S>M || S<0)
    return NA_REAL;
  
  
  S-=2.0; // Because upper tail, including X, is wanted
  long iS=(long)floor(S); // Round down to even, because S values must be even
  iS=2*(iS/2);
  S=iS;
  S=maxm(1L,iS);
  
  
  
  if (CheckFriedmanExactF(r,n,X,&F,false,doRho)) {
    return F; // Exact probabilty at X
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



double xfrie(
    double P,
    int r,
    int n,
    bool doRho  // When true X is Spearman's rho
)
{
  if (doRho)
    n=2;
  if (r<3 || n<2)
    return NA_REAL;
  
  // Initial guess
  double M=(double)(n*n*r*(r*r-1))/12.0; // Max value of S
  double a=(double)(r-1)-2.0/(double)n;
  double b=a*(double)(n-1);
  double W=1.0-qbeta(1.0-P,b/2.0,a/2.0,true,false); // W is corrected for continuity
  double S=1.0+(M+2.0)*W;
  long iS=(long)ceil(S);
  
  
  if (0>P || P>1)
    return NA_REAL;
  
  iS=2*(iS/2); // Round up to even, because S values must be even
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
    int N,  // number of samples
    int r,
    int n,
    bool doRho  // When true X is Spearman's rho
)
{
  GetRNGstate();
  
  
  for (int i=0;i<N;i++){
    randArray[i]=(double)xfrie(unif_rand(),r,n,doRho);
  }
  PutRNGstate();
}



// Returns the median of the Friedman dist -- this is never called for Spearman's Rho,
// because the median for Rho is 0.
double medianfrie(
    int r,
    int n
)
{
  if (! DoExactFriedman(r,n,false)) {
    return xfrie(0.5,r,n,false); // continuous approx
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
    return alpha*low+(1.0-alpha)*high; // interpolate
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







