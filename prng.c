#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "prng.h"



#define PI 3.141592653589793

int64__t LLMedSelect(int k, int n, int64__t arr[]) {
  int64__t a, temp;

  int i, ir, j, mid, l;

  l=1;
  ir=n;
  for (;;) {
    if (ir <= l+1) {
      if (ir == l+1 && arr[ir] < arr[l]) {
	SWAP(arr[l],arr[ir])
	  }
      return arr[k];
    }
    else
      {
	mid=(l+ir) >> 1;
	SWAP(arr[mid],arr[l+1])
	  if (arr[l] > arr[ir]) {
	    SWAP(arr[l],arr[ir])
	      }
	if (arr[l+1] > arr[ir]) {
	  SWAP(arr[l+1],arr[ir])
	    }
	if (arr[l] > arr[l+1]) {
	  SWAP(arr[l],arr[l+1])
	    }
	i=l+1;
	j=ir;
	a=arr[l+1];
	for (;;) {
	  do i++; while (arr[i] < a);
	  do j--; while (arr[j] > a);
	  if (j < i) break;
	  SWAP(arr[i],arr[j])
	    }
	arr[l+1]=arr[j];
	arr[j]=a;
	if (j >= k) ir=j-1;
	if (j <= k) l=i;
      }
  }

    }


int MedSelect(int k, int n, int arr[]) {
  int a, temp;

  int i, ir, j, mid, l;

  l=1;
  ir=n;
  for (;;) {
    if (ir <= l+1) {
      if (ir == l+1 && arr[ir] < arr[l]) {
	SWAP(arr[l],arr[ir])
	  }
      return arr[k];
    }
    else
      {
	mid=(l+ir) >> 1;
	SWAP(arr[mid],arr[l+1])
	  if (arr[l] > arr[ir]) {
	    SWAP(arr[l],arr[ir])
	      }
	if (arr[l+1] > arr[ir]) {
	  SWAP(arr[l+1],arr[ir])
	    }
	if (arr[l] > arr[l+1]) {
	  SWAP(arr[l],arr[l+1])
	    }
	i=l+1;
	j=ir;
	a=arr[l+1];
	for (;;) {
	  do i++; while (arr[i] < a);
	  do j--; while (arr[j] > a);
	  if (j < i) break;
	  SWAP(arr[i],arr[j])
	    }
	arr[l+1]=arr[j];
	arr[j]=a;
	if (j >= k) ir=j-1;
	if (j <= k) l=i;
      }
  }
}


long hash31(int64__t a, int64__t b, int64__t x)
{

  int64__t result;
  long lresult;

  // return a hash of x using a and b mod (2^31 - 1)
// may need to do another mod afterwards, or drop high bits
// depending on d, number of bad guys
// 2^31 - 1 = 2147483647

  //  result = ((int64__t) a)*((int64__t) x)+((int64__t) b);
  result=(a * x) + b;
  result = ((result >> HL) + result) & MOD;
  lresult=(long) result;

  return(lresult);
}

long fourwise(int64__t a, int64__t b, int64__t c, int64__t d, int64__t x)
{
  int64__t result;
  long lresult;

  // returns values that are 4-wise independent by repeated calls
  // to the pairwise indpendent routine.

  result = hash31(hash31(hash31(a,b,x),c,x),d,x); 
  lresult = (long) result;
  return lresult;
}


/*************************************************************************/
/* First, some pseudo-random number generators sourced from other places */
/*************************************************************************/

// There are *THREE* alternate implementations of PRNGs here.
// One taken from Numerical Recipes in C, the second from www.agner.org
// The third is an internal C random library, srand
// The variable usenric controls which one is used: pick one
// and stick with it, switching between the two will give unpredictable
// results.  This is controlled by the randinit procedure, call it with
// usenric == 1 to use the Numerical Recipes gens
// usenric == 2 to use the agner.org PRNGs or
// usenric == 3 to use the inbuilt C routines

// from the math library:
extern double sqrt(double);

// following definitions needed for the random number generator
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran1(prng_type * prng) {

  // A Random Number Generator that picks a uniform [0,1] random number
  // From Numerical Recipes, page 280
  // Should be called with a NEGATIVE value of idum to initialize
  // subsequent calls should not alter idum

  int j;
  long k;
  float temp;

  if (prng->floatidum <= 0 || !prng->iy) {
    if (-(prng->floatidum) < 1) prng->floatidum=1;
    else prng->floatidum = -(prng->floatidum);
    for (j=NTAB+7;j>=0;j--) {
      k=(prng->floatidum)/IQ;
      prng->floatidum=IA*(prng->floatidum-k*IQ)-IR*k;
      if (prng->floatidum < 0) prng->floatidum+=IM;
      if (j<NTAB) prng->iv[j]=prng->floatidum;
    }
    prng->iy=prng->iv[0];
  }
  k = (prng->floatidum)/IQ;
  prng->floatidum=IA*(prng->floatidum-k*IQ)-IR*k;
  if (prng->floatidum<0) prng->floatidum += IM;
  j = prng->iy/NDIV;
  prng->iy=prng->iv[j];
  prng->iv[j]=prng->floatidum;
  if ((temp=AM*prng->iy) > RNMX) return RNMX;
  else return temp;
}

long ran2(prng_type * prng) {

  // A Random Number Generator that picks a uniform random number
  // from the range of long integers.
  // From Numerical Recipes, page 280
  // Should be called with a NEGATIVE value of idum to initialize
  // subsequent calls should not alter idum
  // This is a hacked version of the above procedure, so proceed with
  // caution.

  int j;
  long k;

  if (prng->intidum <= 0 || !prng->iy) {
    if (-(prng->intidum) < 1) prng->intidum=1;
    else prng->intidum = -(prng->intidum);
    for (j=NTAB+7;j>=0;j--) {
      k=(prng->intidum)/IQ;
      prng->intidum=IA*(prng->intidum-k*IQ)-IR*k;
      if (prng->intidum < 0) prng->intidum+=IM;
      if (j<NTAB) prng->iv[j]=prng->intidum;
    }
    prng->iy=prng->iv[0];
  }
  k = (prng->intidum)/IQ;
  prng->intidum=IA*(prng->intidum-k*IQ)-IR*k;
  if (prng->intidum<0) prng->intidum += IM;
  j = prng->iy/NDIV;
  prng->iy=prng->iv[j];
  prng->iv[j]=prng->intidum;
  return prng->iy;
}

/**********************************************************************/

// Following routines are from www.agner.org

/************************* RANROTB.C ******************** AgF 1999-03-03 *
*  Random Number generator 'RANROT' type B                               *
*                                                                        *
*  This is a lagged-Fibonacci type of random number generator with       *
*  rotation of bits.  The algorithm is:                                  *
*  X[n] = ((X[n-j] rotl r1) + (X[n-k] rotl r2)) modulo 2^b               *
*                                                                        *
*  The last k values of X are stored in a circular buffer named          *
*  randbuffer.                                                           *
*                                                                        *
*  This version works with any integer size: 16, 32, 64 bits etc.        *
*  The integers must be unsigned. The resolution depends on the integer  *
*  size.                                                                 *
*                                                                        *
*  Note that the function RanrotAInit must be called before the first    *
*  call to RanrotA or iRanrotA                                           *
*                                                                        *
*  The theory of the RANROT type of generators is described at           *
*  www.agner.org/random/ranrot.htm                                       *
*                                                                        *
*************************************************************************/

// this should be almost verbatim from the above webpage.
// although it's been hacked with a little bit...

unsigned long rotl (unsigned long x, unsigned long r) {
  return (x << r) | (x >> (sizeof(x)*8-r));}

/* define parameters (R1 and R2 must be smaller than the integer size): */
#define JJ  10
#define R1   5
#define R2   3

/* returns some random bits */
unsigned long ran3(prng_type * prng) {

  unsigned long x;

  /* generate next random number */

  x = prng->randbuffer[prng->r_p1] = rotl(prng->randbuffer[prng->r_p2], R1)
    +  rotl(prng->randbuffer[prng->r_p1], R2);
  /* rotate list pointers */
  if (--prng->r_p1 < 0) prng->r_p1 = KK - 1;
  if (--prng->r_p2 < 0) prng->r_p2 = KK - 1;
  /* conversion to float */
  return x;
}

/* returns a random number between 0 and 1 */
double ran4(prng_type * prng) {

  /* conversion to floating point type */
  return (ran3(prng) * prng->scale);
}

/* this function initializes the random number generator.      */
/* Must be called before the first call to RanrotA or iRanrotA */
void RanrotAInit (prng_type * prng, unsigned long seed) {

  int i;

  /* put semi-random numbers into the buffer */
  for (i=0; i<KK; i++) {
    prng->randbuffer[i] = seed;
    seed = rotl(seed,5) + 97;}

  /* initialize pointers to circular buffer */
  prng->r_p1 = 0;  prng->r_p2 = JJ;

  /* randomize */
  for (i = 0;  i < 300;  i++) ran3(prng);
  prng->scale = ldexp(1.0f, -8.0f * sizeof(unsigned long));
}


/**********************************************************************/
/* These are wrapper procedures for the uniform random number gens    */
/**********************************************************************/

long prng_int(prng_type * prng) {

  // returns a pseudo-random long integer.  Initialise the generator
  // before use!

  long response=0;

  switch (prng->usenric)
    {
    case 1 : response=(ran2(prng)); break;
    case 2 : response=(ran3(prng)); break;
    case 3 : response=(lrand48()); break;
    }
  return response;
}


float prng_float(prng_type * prng) {

  // returns a pseudo-random float in the range [0.0,1.0].
  // Initialise the generator before use!
  float result=0;

  switch (prng->usenric)
    {
    case 1 : result=(ran1(prng)); break;
    case 2 : result=(ran4(prng)); break;
    case 3 : result=(drand48()); break;
    }
  return result;
}

prng_type * prng_Init(long seed, int nric) {

  // Initialise the random number generators.  nric determines
  // which algorithm to use, 1 for Numerical Recipes in C,
  // 0 for the other one.
  prng_type * result;

  result=(prng_type *) calloc(1,sizeof(prng_type));

  result->iy=0;
  result->usenric=nric;
  result->floatidum=-1;
  result->intidum=-1;
  result->iset=0;
  // set a global variable to record which algorithm to use
  switch (nric)
    {
    case 2 :
      RanrotAInit(result,seed);
      break;
    case 1 :
      if (seed>0) {
	// to initialise the NRiC PRNGs, call it with a negative value
	// so make sure it gets a negative value!
	result->floatidum = -(seed);  result->intidum = -(seed);
      } else {
	result->floatidum=seed; result->intidum=seed;
      }
      break;
    case 3 :
      srand48(seed);
      break;
    }

  prng_float(result);
  prng_int(result);
  // call the routines to actually initialise them
  return(result);
}

void prng_Reseed(prng_type * prng, long seed)
{
  switch (prng->usenric)
    {
    case 2 :
      RanrotAInit(prng,seed);
      break;
    case 1 :
      if (seed>0) {
	// to initialise the NRiC PRNGs, call it with a negative value
	// so make sure it gets a negative value!
	prng->floatidum = -(seed);  prng->intidum = -(seed);
      } else {
	prng->floatidum=seed; prng->intidum=seed;
      }
      break;
    case 3 :
      srand48(seed);
      break;
    }
}

void prng_Destroy(prng_type * prng)
{
  free(prng);
}

/**********************************************************************/
/* Next, a load of routines that convert uniform random variables     */
/* from [0,1] to stable distribitions, such as gaussian, levy or      */
/* general                                                            */
/**********************************************************************/

double prng_normal(prng_type * prng) {

  // Pick random values distributed N(0,1) using the Box-Muller transform
  // Taken from numerical recipes in C p289
  // picks two at a time, returns one per call (buffers the other)

  double fac,rsq,v1,v2;
  if (prng->iset == 0) {
    do {
      v1 = 2.0*prng_float(prng)-1.0;
      v2 = 2.0*prng_float(prng)-1.0;
      rsq=v1*v1+v2*v2;
    } while (rsq >= 1.0 || rsq == 0.0);

    fac = sqrt((double) -2.0*log((double)rsq)/rsq);
    prng->gset=v1*fac;
    prng->iset=1;
    return v2*fac;
  }
  else {
    prng->iset = 0;
    return prng->gset;
  }
}

double prng_stabledbn(prng_type * prng, double alpha) {

  // From 'stable distributions', John Nolan, manuscript, p24
  // we set beta = 0 by analogy with the normal and cauchy case
  // identical to the above routine, but returns a double instead
  // of a long double (you'll see this a lot...)

  double theta, W, holder, left, right;

  theta=PI*(prng_float(prng) - 0.5);
  W = -log(prng_float(prng)); // takes natural log

  //  printf("theta %Lf, W = %Lf \n", theta, W);

  // some notes on Nolan's notes:
  // if beta == 0 then c(alpha,beta)=1; theta_0 = 0
  // expression reduces to sin alpha.theta / (cos theta) ^1/alpha
  //  * (cos (theta - alpha theta)/W) ^(1-alpha)/alpha
  left = (sin(alpha*theta)/pow(cos(theta), 1.0/alpha));
  right= pow(cos(theta*(1.0 - alpha))/W, ((1.0-alpha)/alpha));
  holder=left*right;
  return(holder);
}


long double prng_cauchy(prng_type * prng) {

  // return a value from the cauchy distribution
  // using the formula in 'Stable Distributions', p23
  // this is distributed Cauchy(1,0)

  return(tan(PI*(prng_float(prng) - 0.5)));
}


double prng_altstab(prng_type * prng, double p)
{
  double u,v,result;

  u=prng_float(prng);
  v=prng_float(prng);
  result=pow(u,p);
  // result=exp(p*log(u));
  if (v<0.5) result=-result;
  return(result);
}


/*
long double levy() {

  // this would give the levy distribution, except it doesn't get used

  long double z;

  z=gasdev();
  return (1.0/(z*z));
  }

*/

double prng_stable(prng_type * prng, double alpha) {

  // wrapper for the stable distributions above:
  // call the appropriate routine based on the value of alpha given
  // initialising it with the seed in idum

  // randinit must be called before entering this procedure for
  // the first time since it uses the random generators


  if (alpha==2.0)
    return(prng_normal(prng));
  else if (alpha==1.0)
    return(prng_cauchy(prng));
  else if (alpha<0.01)
    return(prng_altstab(prng,-50.0));
  else return (prng_stabledbn(prng,alpha));
}

double zeta(long n, double theta)
{

  // the zeta function, used by the below zipf function
  // (this is not often called from outside this library)
  // ... but have made it public now to speed things up

  int i;
  double ans=0.0;

  for (i=1; i <= n; i++)
    ans += pow(1./(double)i, theta);
  return(ans);
}

double fastzipf(double theta, long n, double zetan, prng_type * prng) {

  // this draws values from the zipf distribution
  // this is mainly useful for test generation purposes
  // n is range, theta is skewness parameter
  // theta = 0 gives uniform dbn,
  // theta > 1 gives highly skewed dbn.
  // original code due to Flip Korn, used with permission

	double alpha;
	double eta;
	double u;
	double uz;
	long double val;

  // randinit must be called before entering this procedure for
  // the first time since it uses the random generators

	alpha = 1. / (1. - theta);
	eta = (1. - pow(2./((double) n), 1. - theta)) 
	  / (1. - zeta(2.,theta)/zetan);

	u = prng_float(prng);
	uz = u * zetan;
	if (uz < 1.) val = 1;
	else if (uz < (1. + pow(0.5, theta))) val = 2;
	else val = 1 + (int64__t)(n * pow(eta*u - eta + 1., alpha));

	return(val);
}

/*
long double zipf(double theta, long n) {

  // this draws values from the zipf distribution
  // this is mainly useful for test generation purposes
  // n is range, theta is skewness parameter
  // theta = 0 gives uniform dbn,
  // theta > 1 gives highly skewed dbn.

	double alpha;
	double zetan;
	double eta;
	double u;
	double uz;
	long double val;

  // randinit must be called before entering this procedure for
  // the first time since it uses the random generators

	alpha = 1. / (1. - theta);
	zetan = zeta(n, theta);
	eta = (1. - pow(2./n, 1. - theta)) / (1. - zeta(2.,theta)/zetan);

	u = randomfl();
	uz = u * zetan;
	if (uz < 1.) val = 1;
	else if (uz < (1. + pow(0.5, theta))) val = 2;
	else val = 1 + (int64__t)(n * pow(eta*u - eta + 1., alpha));

	return(val);
}
*/

/*Tools::Random::Random(uint32_t seed)
 : m_type(Tools::RGT_MERSENNE), m_buffer(0), m_seed(seed)
{
	initMersenne();
}

Tools::Random::Random(uint32_t seed, RandomGeneratorType t)
 : m_type(t), m_buffer(0), m_seed(seed)
{
	switch (m_type)
	{
	default:
		throw Tools::NotSupportedException(
			"Tools::Random::Random: This PRG type is not supported yet."
		);
	case RGT_MERSENNE:
		initMersenne();
		break;
	case RGT_DRAND48:
		initDrand(0x330E);
		break;
	}
}

Tools::Random::Random(uint32_t seed, uint16_t xsubi0)
 : m_type(RGT_DRAND48), m_buffer(0), m_seed(seed)
{
	initDrand(xsubi0);
}

Tools::Random::~Random()
{
	switch (m_type)
	{
	case RGT_MERSENNE:
		delete[] reinterpret_cast<uint32_t*>(m_buffer);
		break;
	case RGT_DRAND48:
		break;
	}
}

void Tools::Random::initMersenne()
{
	m_buffer = new uint32_t[Tools::Random::MERS_N + 1];

	uint32_t& mti = *(reinterpret_cast<uint32_t*>(m_buffer));
	uint32_t* mt = reinterpret_cast<uint32_t*>(m_buffer) + 1;

	mt[0] = m_seed;

	for (mti = 1; mti < MERS_N; mti++)
    {
		mt[mti] = (1812433253UL * (mt[mti - 1] ^ (mt[mti - 1] >> 30)) + mti);
	}

	assert(mti == *(reinterpret_cast<uint32_t*>(m_buffer)));

	m_architecture = Tools::System::getArchitecture();
}

void Tools::Random::initDrand(uint16_t xsubi0)
{
	m_xsubi[0] = xsubi0;
	uint32_t mask = 0xFFFF;
	m_xsubi[1] = static_cast<uint16_t>(m_seed & mask);
	mask = mask << 16;
	m_xsubi[2] = static_cast<uint16_t>((m_seed & mask) >> 16);
		// srand48 needs a 48 bit seed, irrespective
		// of the size of unsigned long and unsigned short.
}

int32_t Tools::Random::nextUniformLong()
{
	if (m_type == RGT_DRAND48)
	{
		return jrand48(m_xsubi);
			// Careful: jrand48 modifies m_xsubi after the call.
	}
	else if (m_type == RGT_MERSENNE)
	{
		// generate 32 random bits
		uint32_t y;
		uint32_t& mti = *(reinterpret_cast<uint32_t*>(m_buffer));
		uint32_t* mt = reinterpret_cast<uint32_t*>(m_buffer) + 1;

		if (mti >= MERS_N)
		{
			// generate MERS_N words at one time
			const uint32_t LOWER_MASK = (1LU << MERS_R) - 1;
				// lower MERS_R bits
			const uint32_t UPPER_MASK = -1L  << MERS_R;
				// upper (32 - MERS_R) bits
			static const uint32_t mag01[2] = {0, MERS_A};
    
			uint32_t kk;
			for (kk = 0; kk < MERS_N - MERS_M; kk++)
			{
				y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
				mt[kk] = mt[kk + MERS_M] ^ (y >> 1) ^ mag01[y & 1];
			}

			for (; kk < MERS_N - 1; kk++)
			{
				y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
				mt[kk] = mt[kk + (MERS_M - MERS_N)] ^ (y >> 1) ^ mag01[y & 1];
			}

			y = (mt[MERS_N - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
			mt[MERS_N - 1] = mt[MERS_M - 1] ^ (y >> 1) ^ mag01[y & 1];
			mti = 0;
		}

		y = mt[mti++];

		// Tempering (May be omitted):
		y ^=  y >> MERS_U;
		y ^= (y << MERS_S) & MERS_B;
		y ^= (y << MERS_T) & MERS_C;
		y ^=  y >> MERS_L;

		return y;
	}
	else
	{
		throw Tools::IllegalStateException(
			"Tools::Random::nextUniformLong: Should never reach here."
		);
	}
}

uint32_t Tools::Random::nextUniformUnsignedLong()
{
	return static_cast<uint32_t>(nextUniformLong());
}

int32_t Tools::Random::nextUniformLong(int32_t low, int32_t high)
{
	return low + static_cast<int32_t>((high - low) * nextUniformDouble());
}

uint32_t Tools::Random::nextUniformUnsignedLong(uint32_t low, uint32_t high)
{
	return low + static_cast<uint32_t>((high - low) * nextUniformDouble());
}

int64__t Tools::Random::nextUniformLongLong()
{
	return static_cast<int64__t>(nextUniformUnsignedLongLong());
}

uint64__t Tools::Random::nextUniformUnsignedLongLong()
{
	uint64__t high = static_cast<uint64__t>(nextUniformUnsignedLong());
	uint64__t low = static_cast<uint64__t>(nextUniformUnsignedLong());
	return (high << 32) | low;
}

int64__t Tools::Random::nextUniformLongLong(int64__t low, int64__t high)
{
	return low + static_cast<int64__t>((high - low) * nextUniformDouble());
}

uint64__t Tools::Random::nextUniformUnsignedLongLong(uint64__t low, uint64__t high)
{
	return low + static_cast<uint64__t>((high - low) * nextUniformDouble());
}

int16_t Tools::Random::nextUniformShort()
{
	return static_cast<int16_t>(nextUniformUnsignedShort());
}

uint16_t Tools::Random::nextUniformUnsignedShort()
{
	return nextUniformUnsignedLong() >> 16;
		// retain the high order bits.
}

double Tools::Random::nextUniformDouble()
{
	if (m_type == RGT_DRAND48)
	{
		return erand48(m_xsubi);
			// Careful: erand48 modifies m_xsubi after the call.
	}
	else if (m_type == RGT_MERSENNE)
	{
		union {double f; uint32_t i[2];} convert;
		uint32_t r = nextUniformUnsignedLong();

		switch (m_architecture)
		{
		case ARCH_LITTLEENDIAN:
			convert.i[0] =  r << 20;
			convert.i[1] = (r >> 12) | 0x3FF00000;
			return convert.f - 1.0;
		case ARCH_BIGENDIAN:
			convert.i[1] =  r << 20;
			convert.i[0] = (r >> 12) | 0x3FF00000;
			return convert.f - 1.0;
		case ARCH_NONIEEE:
		default:
			;
		}

		// This somewhat slower method works for all architectures, including 
		// non-IEEE floating point representation:
		return
			static_cast<double>(r) *
			(1.0 / static_cast<double>(static_cast<uint32_t>(-1L) + 1.0));
	}
	else
	{
		throw Tools::IllegalStateException(
			"Tools::Random::nextUniformDouble: Should never reach here."
		);
	}
}

double Tools::Random::nextUniformDouble(double low, double high)
{
	return (high - low) * nextUniformDouble() + low;
}

// mean 0.0, standard deviation 1.0
double Tools::Random::nextNormalDouble()
{
	static bool haveNextNextGaussian = false;
	static double nextNextGaussian;

	if (haveNextNextGaussian)
	{
		haveNextNextGaussian = false;
		return nextNextGaussian;
	}
	else
	{
		double v1, v2, s;

		do
		{
			v1 = 2 * nextUniformDouble() - 1;   // between -1.0 and 1.0
			v2 = 2 * nextUniformDouble() - 1;   // between -1.0 and 1.0
			s = v1 * v1 + v2 * v2;
		}
		while (s >= 1 || s == 0);

		double multiplier = std::sqrt(-2.0 * std::log(s)/s);
		nextNextGaussian = v2 * multiplier;
		haveNextNextGaussian = true;
		return v1 * multiplier;
	}
}

double Tools::Random::nextNormalDouble(double m, double std)
{
 	return m + (std * nextNormalDouble());
}

int32_t Tools::Random::nextSkewedLong(int32_t low, int32_t high, Level p)
{
	return low + static_cast<int32_t>((high - low) * nextSkewedDouble(p));
}

double Tools::Random::nextSkewedDouble(double low, double high, Level p)
{
	return (high - low) * nextSkewedDouble(p) + low;
}

// WARNING: The inversion method is very slow.
// For discrete distributions use PRGZipf.
double Tools::Random::nextSkewedDouble(Level p)
{
	double HsubV, l;
	uint64__t V = 1000000;

	//HsubV = 0.0;
	//for(uint64__t i = 1; i <= V; i++) HsubV += 1.0 / std::pow(static_cast<double>(i), l);

	switch (p)
	{
	case LVL_VERYLOW:
		HsubV = 1998.54;
		l = 0.5;
		break;
	case LVL_LOW:
		HsubV = 14.3927;
		l = 1.0;
		break;
	case LVL_MEDIUM:
		HsubV = 2.61038;
		l = 1.5;
		break;
	case LVL_HIGH:
		HsubV = 1.64493;
		l = 2.0;
		break;
	case LVL_VERYHIGH:
		HsubV = 1.34149;
		l = 2.5;
		break;
	default:
		throw Tools::IllegalArgumentException(
			"Tools::Random::nextSkewedDouble: Unknown skedeness level."
		);
	}

	double r = nextUniformDouble() * HsubV;
	double sum = 1.0;
	uint64__t i = 1;
	while(sum < r)
	{
		i++;
		sum += 1.0 / std::pow(static_cast<double>(i), l);
	}

	// i follows Zipf distribution and lies between 1 and V
	return (static_cast<double>(i) - 1.0) / (static_cast<double>(V) - 1.0);
}

bool Tools::Random::flipCoin()
{
	if (m_type == RGT_DRAND48)
	{
		if (nextUniformDouble() < 0.5) return true;
		return false;
	}
	else if (m_type == RGT_MERSENNE)
	{
		// the probability of a long being even or odd is 50%
		if ((nextUniformLong() & 1) == 1) return true;
		return false;
	}
	else
	{
		throw Tools::NotSupportedException(
			"Tools::Random::flipCoin: This PRG type is not supported yet."
		);
	}
}

bool Tools::Random::bernulliTrial(double p)
{
	assert(p >= 0.0 && p <= 1.0);

	if (nextUniformDouble() < p) return true;
	return false;	
}

size_t Tools::Random::getSize() const
{
	// the seed is enough.
	return sizeof(uint32_t);
}

uint32_t Tools::Random::getSeed() const
{
	return m_seed;
}

Tools::PRGZipf::PRGZipf(int32_t min, int32_t max, double s, Tools::Random* pRandom)
 : m_min(min), m_max(max), m_s(s), m_pRandom(pRandom)
{
	initLookupTable();
}

Tools::PRGZipf::~PRGZipf()
{
	delete[] m_pLookupTable;
}

void Tools::PRGZipf::initLookupTable()
{
	uint32_t N = m_max - m_min;
	double Hns = 0.0;

	for (uint32_t k = 1; k < N; k++)
		Hns += 1.0 / std::pow(static_cast<double>(k), m_s);

	m_pLookupTable = new double[N];

	double sum = 0.0;
	m_pLookupTable[0] = sum;
	for (uint32_t k = 1; k < N; k++)
	{
		sum += 1.0 / std::pow(static_cast<double>(k), m_s);
		m_pLookupTable[k] = sum / Hns;
	}
}

int32_t Tools::PRGZipf::nextLong()
{
	uint32_t N = m_max - m_min;
	double dart = m_pRandom->nextUniformDouble();
	double* i = std::lower_bound(m_pLookupTable, m_pLookupTable + N, dart);

	assert(i >= m_pLookupTable && i <= m_pLookupTable + N - 1);

	int32_t ret = static_cast<int32_t>(i - m_pLookupTable) + m_min;

	assert(ret >= m_min && ret < m_max);
	return ret;
}

Tools::Architecture Tools::System::getArchitecture()
{
	union {double f; uint32_t i[2];} convert;
	convert.f = 1.0;

	// Note: Old versions of the Gnu g++ compiler may make an error here,
	// compile with the option  -fenum-int-equiv  to fix the problem
	if (convert.i[1] == 0x3FF00000) return ARCH_LITTLEENDIAN;
	else if (convert.i[0] == 0x3FF00000) return ARCH_BIGENDIAN;
	else return ARCH_NONIEEE;
}

Tools::IllegalArgumentException::IllegalArgumentException(std::string s) : m_error(s)
{
}

std::string Tools::IllegalArgumentException::what()
{
	return "IllegalArgumentException: " + m_error;
}

Tools::IllegalStateException::IllegalStateException(std::string s) : m_error(s)
{
}

std::string Tools::IllegalStateException::what()
{
	return "IllegalStateException: " + m_error;
}

Tools::NotSupportedException::NotSupportedException(std::string s) : m_error(s)
{
}

std::string Tools::NotSupportedException::what()
{
	return "NotSupportedException: " + m_error;
}
*/
