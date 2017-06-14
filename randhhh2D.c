/*
Implementation of the 2-dimensional Randomized Hierarchical Heavy Hitters algorithm (RHHH) for IP addresses
-Ran Ben Basat (sran@cs.technion.ac.il) 2017-05-25
-Based on the "Hierarchical Heavy Hitters with the Space Saving Algorithm" paper implementation by Thomas Steinke
	http://people.seas.harvard.edu/~tsteinke/hhh/
*/

#include <stdlib.h>
#include <stdio.h>

//for debugging purposes only
#include <assert.h> 

#include "ulossycount.h"
#define COUNTERSIZE k
#define COUNTERS items
#define COUNT parentg->count

#include "randhhh2D.h"


#ifndef DIMENSION2
#error This is the two-dimensional implementation 
#endif

int min(int a, int b) {return (a <= b ? a : b);}
int max(int a, int b) {return (a >= b ? a : b);}

//The counters
LCU_type * counters[NUM_COUNTERS];
int nrUpdates;
double SSepsval;
#ifdef PROB
double ignoreProb;
double logIgnoreProb;
double minusQuantity;
int nrIgnore;
#endif

double dblmax(double a, double b) {return (a >= b ? a : b);}

double twototheminus(int k) {
	double ans = 1;
	while (k > 0) {ans /= 2; k--;}
	return ans;
}

//initialise
void init(double SSepsilon, double prob) {
	int i;
	SSepsval = SSepsilon;
	srand(SEED);
#ifdef PROB
	ignoreProb = 1 - prob;
	logIgnoreProb = log(ignoreProb);
	minusQuantity = log(RAND_MAX) / logIgnoreProb;
	nrIgnore = log((float)rand()) / logIgnoreProb - minusQuantity;
#endif
	for (i = 0; i < NUM_COUNTERS; i++)
		counters[i] = LCU_Init(SSepsval);
}

//deinitialise
void deinit() {
	int i;
	for (i = 0; i < NUM_COUNTERS; i++)
		LCU_Destroy(counters[i]);
}

//update an input
void update(LCLitem_t item) {
#ifdef  PROB
	if (nrIgnore--) return;
	short i = rand() % NUM_COUNTERS;
	LCU_Update(counters[i], item & masks[i]);
	nrIgnore = log((float)rand()) / logIgnoreProb - minusQuantity;
#else
	short i = rand() % NUM_COUNTERS;
	LCU_Update(counters[i], item & masks[i]);
#endif
}

//we want to sort heavyhitters
int cmpHH(const void * lhs, const void * rhs) {
	if (((const HeavyHitter*) lhs)->item > ((const HeavyHitter*) rhs)->item) return 1;
	if (((const HeavyHitter*) lhs)->item < ((const HeavyHitter*) rhs)->item) return -1;
	if (((const HeavyHitter*) lhs)->mask > ((const HeavyHitter*) rhs)->mask) return 1;
	if (((const HeavyHitter*) lhs)->mask < ((const HeavyHitter*) rhs)->mask) return -1;
	if (((const HeavyHitter*) lhs)->upper != ((const HeavyHitter*) rhs)->upper) return ((const HeavyHitter*) lhs)->upper - ((const HeavyHitter*) rhs)->upper;
	return ((const HeavyHitter*) lhs)->lower - ((const HeavyHitter*) rhs)->lower;
}

typedef struct descendant {
	LCLitem_t id;
	int mask;
} Descendant;

//the two-dimensional output
HeavyHitter * output2(int threshold, int * numhitters) {
	HeavyHitter * output; //This will be the heavy hitters to output
	int n = 0; //the number of items in output
	int outputspace = 5, hpspace = 5;
#ifdef PROB
	float adjustedThreshold = (1-ignoreProb) * threshold / ((float) NUM_COUNTERS);
#else
	float adjustedThreshold = threshold / ((float) NUM_COUNTERS);
#endif	
	int i, j, k, l, m, im;
	LCLitem_t newIP;
	int s;
	Descendant * Hp; //the direct descendants
	int numdescendants;
	int z = 4; //TODO: replace with a function of delta_s

	output = (HeavyHitter*) calloc(sizeof(HeavyHitter), outputspace); //ensure that it is sufficient
	Hp = (Descendant*) calloc(hpspace, sizeof(Descendant));
	
	for (i = 0; i < NUM_COUNTERS; i++) {
		for (j = 0; j < counters[i]->COUNTERSIZE; j++) {
			if (counters[i]->COUNTERS[j].item != LCL_NULLITEM) {
				//Now we just have to check that the counts are sufficient
				//s = amount that needs to be subtracted from the count
				
				if (counters[i]->COUNTERS[j].COUNT < adjustedThreshold) continue; // no point continuing

				numdescendants = 0; // |Hp| = 0
				for (k = 0; k < n; k++) {
					if ((masks[i] & masks[output[k].mask]) == masks[i] && (counters[i]->COUNTERS[j].item & masks[i]) == (output[k].item & masks[i])) {
						//output[k] is a descendant---is it direct?
						//assume it is for now, but check to make sure the assumption was valid for the previous ones
						l = 0;
						for (m = 0; m < numdescendants; m++) {
							Hp[l] = Hp[m];
							if ((masks[i] & masks[Hp[m].mask]) != masks[i] || (Hp[m].id & masks[output[k].mask]) != (output[k].item & masks[output[k].mask])) l++;
						}	
						numdescendants = l;
						Hp[numdescendants].mask = output[k].mask;
						Hp[numdescendants].id = output[k].item;
						numdescendants++;
						while (numdescendants >= hpspace) {hpspace *= 2; Hp = (Descendant*) realloc(Hp, sizeof(Descendant)*hpspace);}
					}
				}

				//now compute s
				s = 0;
				for (k = 0; k < numdescendants; k++) {
					s += LCU_PointEstLow(counters[Hp[k].mask], Hp[k].id);
				}
				for (k = 0; k < numdescendants; k++) {
					for (l = k + 1; l < numdescendants; l++) {

						if ((Hp[k].id & (masks[Hp[k].mask] & masks[Hp[l].mask])) != (Hp[l].id & (masks[Hp[k].mask] & masks[Hp[l].mask]))) continue; //there is no ip common to the subnets k and l

						newIP = (Hp[k].id | Hp[l].id);
						//compute the right mask
						im = 0; while(masks[im] != (masks[Hp[k].mask] | masks[Hp[l].mask])) im++;
						assert(masks[im] == (masks[Hp[k].mask] | masks[Hp[l].mask]));
						s -= LCU_PointEstUpp(counters[im], newIP);
					}
				}

				//now compare to the threshold
#ifdef PROB				
				if ((counters[i]->COUNTERS[j].COUNT - s  + 2*z*sqrt((1-ignoreProb) * 2 * counters[i]->COUNTERS[j].COUNT ) >= adjustedThreshold ) ){
#else
				if ((counters[i]->COUNTERS[j].COUNT - s  + 2*z*sqrt( 2 * counters[i]->COUNTERS[j].COUNT ) >= adjustedThreshold ) ){	
#endif	
					//Add this item to our list of heavy hitters
					while (n >= outputspace) {outputspace *= 2; output = (HeavyHitter *) realloc(output, outputspace * sizeof(HeavyHitter));}
					output[n].item = counters[i]->COUNTERS[j].item;
					output[n].mask = i;
					output[n].upper = counters[i]->COUNTERS[j].COUNT;
					output[n].lower = output[n].upper - counters[i]->COUNTERS[j].delta;
					n++;
				}

			}
		}
	}
	
	free(Hp);

	//now clean up the output
	output = (HeavyHitter*) realloc(output, n * sizeof(HeavyHitter));
	//qsort(output, n, sizeof(HeavyHitter), &cmpHH);

	*numhitters = n;

	return output;
}



