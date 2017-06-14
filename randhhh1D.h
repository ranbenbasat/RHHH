/*
Implementation of the 1-dimensional Randomized Hierarchical Heavy Hitters algorithm (RHHH) for IP addresses
-Ran Ben Basat (sran@cs.technion.ac.il) 2017-05-25
-Based on the "Hierarchical Heavy Hitters with the Space Saving Algorithm" paper implementation by Thomas Steinke -- http://people.seas.harvard.edu/~tsteinke/hhh/
*/

#ifndef RANDHHH1D_H
#define RANDHHH1D_H

#include "ulossycount.h"

#ifndef SEED
#define SEED 3421
#endif

//#define FASTRAND 1

#ifndef NUM_MASKS
#define NUM_MASKS 5
#endif
#define NUM_COUNTERS NUM_MASKS

//The masks associated with the counters
//Note that we must ensure that they are in increasing order of generality
extern LCLitem_t masks[NUM_COUNTERS];
extern int leveleps[NUM_COUNTERS];

//initialise
void init(double SSepsilon, double prob);

//deinitialise
void deinit();

void update(LCLitem_t item, int count);


//struct to store a heavy hitter output
typedef struct heavyhitter {
	LCLitem_t item;
	int mask; //The item & mask
	int upper, lower; //Upper and lower count bounds
} HeavyHitter;


//the one-dimensional output
HeavyHitter * output(int threshold, int * numhitters, int streamLen);
#endif

