/*
Hashtable implementation for use in the output stage of the HHH algorithm.
-Thomas Steinke (tsteinke@seas.harvard.edu) 2010-10-06
*/

#ifndef HASHTABLE_H
#define HASHTABLE_H

#include "prng.h"

//linked list hash table entry
typedef struct ll LL;
struct ll {
	int64__t item;
	int val;
	LL * next, * prev;
};

LL ** HT_Init(int size); //create a hash table with particular size
void HT_Clear(LL ** hashtable, int size); //destroy. It needs ot know the size
LL * HT_Insert(LL ** hashtable, int64__t item, int val, int size); //if item already exists just return the pointer, otherwise insert and set val
LL * HT_Find(LL ** hashtable, int64__t item, int size); //find the LL corresponding to the item. All updates happen via this, as do removes
LL ** HT_FindEntry(LL ** hashtable, int64__t item, int size); //give the start of the list for the hash corresponding to item
#endif

