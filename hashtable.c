/*
Hash Table implementation. See header.
-Thomas Steinke (tsteinke@seas.harvard.edu) 2010-10-06
*/

#include <stdlib.h>

#include "hashtable.h"


//constants
#define HASHA_VAL 151261303 //hash parameter
#define HASHB_VAL 6722461

LL ** HT_Init(int size) { //create a hash table with particular size
	return (LL **) calloc(size, sizeof(LL *));
}

void HT_Clear(LL ** hashtable, int size) { //destroy
	int i;
	LL * ptr, * tmp;
	for (i = 0; i < size; i++) {
		ptr = hashtable[i];
		hashtable[i] = NULL;
		while (ptr) {
			tmp = ptr->next;
			free(ptr);
			ptr = tmp;
		}
	}
	free(hashtable);
}

LL * HT_Insert(LL ** hashtable, int64__t item, int val, int size) { //if item already exists just return the pointer, otherwise insert and set val
	int hash = (int) hash31(HASHA_VAL, HASHB_VAL, item) % size;

	//ptr is the start of the list
	LL * ptr = hashtable[hash];

	//search thru the list
	while (ptr) {
		if (ptr->item == item) return ptr; //if we have found it, we are done
		ptr = ptr->next; //keep looking
	}
	
	//if we get here then ptr = NULL and we aren't in the list yet so we may prepend it to the list
	ptr = (LL*) calloc(1, sizeof(LL));
	ptr->item = item;
	ptr->val = val;
	ptr->next = hashtable[hash];
	ptr->prev = NULL;
	if (hashtable[hash] != NULL) hashtable[hash]->prev = ptr;
	hashtable[hash] = ptr;
	return ptr;
}

LL * HT_Find(LL ** hashtable, int64__t item, int size) { //find the LL corresponding to the item. All updates happen via this, as do removes
	int hash = (int) hash31(HASHA_VAL, HASHB_VAL, item) % size;

	//ptr is the start of the list
	LL * ptr = hashtable[hash];

	//search thru the list
	while (ptr) {
		if (ptr->item == item) return ptr; //if we have found it, we are done
		ptr = ptr->next; //keep looking
	}
	//can't find
	return NULL;
}

LL ** HT_FindEntry(LL ** hashtable, int64__t item, int size) {//give the start of the list for the hash corresponding to item
	int hash = (int) hash31(HASHA_VAL, HASHB_VAL, item) % size;
	return hashtable + hash;
}

