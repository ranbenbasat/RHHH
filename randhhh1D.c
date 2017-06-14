/*
Implementation of the 1-dimensional Randomized Hierarchical Heavy Hitters algorithm (RHHH) for IP addresses
-Ran Ben Basat (sran@cs.technion.ac.il) 2017-05-25
-Based on the "Hierarchical Heavy Hitters with the Space Saving Algorithm" paper implementation by Thomas Steinke -- http://people.seas.harvard.edu/~tsteinke/hhh/
*/

#include <stdlib.h>

#include "ulossycount.h"
#include <math.h>
#define COUNTERSIZE k
#define COUNTERS items
#define COUNT parentg->count


#include "hashtable.h"
#include "randhhh1D.h"


#ifndef VMULT
#define VMULT 1
#endif

#if VMULT>1
#define PROB
#endif

//The counters
LCU_type * counters[NUM_COUNTERS];
double SSepsval;
#ifdef PROB
double ignoreProb;
double logIgnoreProb;
double minusQuantity;
int nrIgnore;
#endif

double max(double a, double b) {return (a >= b ? a : b);}

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
		counters[i] = LCU_Init(max(SSepsilon, twototheminus(leveleps[i])));
}

//deinitialise
void deinit() {
	int i;
	for (i = 0; i < NUM_COUNTERS; i++)
		LCU_Destroy(counters[i]);
}

//update an input
void update(LCLitem_t item, int count) {
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


//the one-dimensional output
HeavyHitter * output(int threshold, int * numhitters, int streamLen) {
	LL ** hashtable; //for transferring counts to parents
	int htsize; //size of hash table;
#ifdef PROB
	float adjustedThreshold = (1-ignoreProb) * threshold / ((float) NUM_COUNTERS);
#else
	float adjustedThreshold = threshold / ((float) NUM_COUNTERS);
#endif
	HeavyHitter * output; //This will be the heavy hitters to output
	int n = 0; //the number of items in output
	int outputspace = 5;
	int i, j;
	LL * tmp;
	int s; //double count
	int z = 4; //TODO: replace with a function of delta_s

	output = (HeavyHitter *) calloc(sizeof(HeavyHitter), outputspace); 

	htsize = (((int)(1.0/SSepsval)) + 1) | 1;
	hashtable = HT_Init(htsize);
	for (i = 0; i < NUM_COUNTERS; i++) {
		for (j = 0; j < counters[i]->COUNTERSIZE; j++) {
			if (counters[i]->COUNTERS[j].item != LCL_NULLITEM) {
				//Now we just have to check that the counts are sufficient
				//first calculate the doublecount s
				tmp = HT_Find(hashtable, ((int64__t) counters[i]->COUNTERS[j].item) << 32 | ((int64__t)i), htsize);
				s = (tmp ? tmp->val : 0); //if there are no doublecounts recoreded assume s=0
				//now erase tmp, we don't need it again
				if (tmp) {
					if (tmp->next) tmp->next->prev = tmp->prev;

					if (tmp->prev) tmp->prev->next = tmp->next; //if prev just update it
					else *(HT_FindEntry(hashtable, ((int64__t)counters[i]->COUNTERS[j].item) << 32 | ((int64__t)i), htsize)) = tmp->next; //else update the start of the list
					free(tmp);
				}

				//now compare to the adjustedThreshold 
#ifdef PROB				
				if (counters[i]->COUNTERS[j].COUNT - s  + 2*z*sqrt((1-ignoreProb) * 2 * counters[i]->COUNTERS[j].COUNT ) >= adjustedThreshold ) {
#else
				if (counters[i]->COUNTERS[j].COUNT - s  + 2*z*sqrt( 2 * counters[i]->COUNTERS[j].COUNT ) >= adjustedThreshold ) {	
#endif	
					while (outputspace <= n) {outputspace *= 2; output = (HeavyHitter *) realloc(output, sizeof(HeavyHitter) * outputspace);}
					//Add this item to our list of heavy hitters
					output[n].item = counters[i]->COUNTERS[j].item;
					output[n].mask = i;
					output[n].upper = counters[i]->COUNTERS[j].COUNT;
					output[n].lower = output[n].upper - counters[i]->COUNTERS[j].delta;
					
					//update s
					s = output[n].lower;
					n++;
				}

				//passs up the doublecount s
				if (s != 0 && i + 1 < NUM_COUNTERS) { //check that there is a positive doublecount and someone to pass up to
					tmp = HT_Insert(hashtable, ((int64__t)(counters[i]->COUNTERS[j].item & masks[i+1])) << 32 | ((int64__t) (i+1)), 0, htsize);
					tmp->val += s;
				}
			}
		}
	}

	HT_Clear(hashtable, htsize);
	
	//now clean up the output
	output = (HeavyHitter *) realloc(output, n * sizeof(HeavyHitter));
	qsort(output, n, sizeof(HeavyHitter), &cmpHH);

	*numhitters = n;

	return output;
}
