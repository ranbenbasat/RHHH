/*
Implementation of the 2-dimensional Randomized Hierarchical Heavy Hitters algorithm (RHHH) for IP addresses
-Ran Ben Basat (sran@cs.technion.ac.il) 2017-05-25
-Based on the "Hierarchical Heavy Hitters with the Space Saving Algorithm" paper implementation by Thomas Steinke -- http://people.seas.harvard.edu/~tsteinke/hhh/
*/
#ifndef RANDHHH2D_H
#define RANDHHH2D_H


#define NUM_COUNTERS 25 //number of masks
#define MAX_DEPTH 5 //depth of masks lattice
#define MAX_DESCENDANTS 512 //maximum number of direct descendents of a given ip pair

#include "prng.h"

//make sure we are passing the right #def around
#ifndef DIMENSION2
#error Invalid dimension
#endif

#ifndef SEED
#define SEED 3421
#endif

#if VMULT>1
#define PROB
#endif

//still define things as in lossycount.h
#define LCLitem_t uint64__t

//The masks associated with the counters
//Note that we must ensure that they are in increasing order of generality
extern LCLitem_t masks[NUM_COUNTERS];

//initialise
#ifdef RANDHHH
	void init(double SSepsilon, double prob);
#else
	void init(double epsilon);
#endif

//deinitialise
void deinit();
void update(LCLitem_t item);


//struct to store a heavy hitter output
typedef struct heavyhitter {
	LCLitem_t item; //The item
	int mask; //The mask id
	int upper, lower; //Upper and lower count bounds
} HeavyHitter;

//the two-dimensional output
HeavyHitter * output2(int threshold, int * numhitters);

#endif


