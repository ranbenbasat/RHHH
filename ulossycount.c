#include <stdlib.h>
#include <stdio.h>
#include "ulossycount.h"
#include "prng.h"


/********************************************************************
Implementation of Frequent algorithm to Find Frequent Items
Based on papers by:
Misra and Gries, 1982
Demaine, Lopez-Ortiz, Munroe, 2002
Karp, Papadimitriou and Shenker, 2003
Implementation by G. Cormode 2002, 2003

Original Code: 2002-11
This version: 2003-10

This work is licensed under the Creative Commons
Attribution-NonCommercial License. To view a copy of this license,
visit http://creativecommons.org/licenses/by-nc/1.0/ or send a letter
to Creative Commons, 559 Nathan Abbott Way, Stanford, California
94305, USA. 
*********************************************************************/

LCU_type * LCU_Init(float fPhi)
{
	int i;
	int k = 1 + (int) 1.0/fPhi;

	LCU_type* result = (LCU_type*) calloc(1,sizeof(LCU_type));

	result->a= (long long) 698124007;
	result->b= (long long) 5125833;
	if (k<1) k=1;
	result->k=k;
	result->n=0;  

	result->tblsz=LCU_HASHMULT*k;  
	result->hashtable=(LCUITEM**) calloc(result->tblsz,sizeof(LCUITEM *));
	result->groups=(LCUGROUP *)calloc(k,sizeof(LCUGROUP));
	result->items=(LCUITEM *) calloc(k,sizeof(LCUITEM));
	result->freegroups=(LCUGROUP **) calloc(k,sizeof(LCUGROUP*));

	for (i=0; i<result->tblsz;i++) 
		result->hashtable[i]=NULL;

	result->root=result->groups;
	result->groups->count=0;
	result->groups->nextg=NULL;
	result->groups->previousg=NULL;

	result->groups->items=result->items;
	for (i=0; i<k;i++)
		result->freegroups[i]=&result->groups[i];
	result->gpt=1; // initialize list of free groups

	for (i=0;i<k;i++) 
	{
		result->items[i].item=LCL_NULLITEM;
		result->items[i].delta=0;
		result->items[i].hash=-1;
		result->items[i].nexti=NULL;
		result->items[i].previousi=NULL;  // initialize values

		result->items[i].parentg=result->groups;
		result->items[i].nexting=&(result->items[i+1]);
		result->items[i].previousing=&(result->items[i-1]); 
		// create doubly linked list
	}
	result->items[0].previousing=&(result->items[k-1]);
	result->items[k-1].nexting=&(result->items[0]);
	// fix start and end of linked list

	return(result);
}  

void LCU_InsertIntoHashtable(LCU_type *lcu, LCUITEM *newi, int i, LCLitem_t newitem){
	newi->nexti=lcu->hashtable[i];
	newi->item=newitem; // overwrite the old item
	newi->hash=i;
	newi->previousi=NULL;
	// insert item into the hashtable
	if (lcu->hashtable[i])
		lcu->hashtable[i]->previousi=newi;
	lcu->hashtable[i]=newi;
}

int LCU_cmp( const void * a, const void * b) {
	LCUITEM * x = (LCUITEM*) a;
	LCUITEM * y = (LCUITEM*) b;
	if (x->parentg->count<y->parentg->count) return -1;
	else if (x->parentg->count>y->parentg->count) return 1;
	else return 0;
}

LCUITEM * LCU_GetNewCounter(LCU_type * lcu) {
	LCUITEM * newi;
	int j;

	newi=lcu->root->items;  // take a counter from the first group
	// but currently it remains in the same group

	newi->nexting->previousing=newi->previousing;
	newi->previousing->nexting=newi->nexting;
	// unhook the new item from the linked list in the hash table	    

	// need to remove this item from the hashtable
	j = newi->hash;
	if (lcu->hashtable[j]==newi)
		lcu->hashtable[j]=newi->nexti;

	if (newi->nexti!=NULL)
		newi->nexti->previousi=newi->previousi;
	if (newi->previousi!=NULL)
		newi->previousi->nexti=newi->nexti;

	return (newi);
}

void LCU_PutInNewGroup(LCU_type * lcu, LCUITEM *newi, LCUGROUP * tmpg){ 
	LCUGROUP * oldgroup;

	oldgroup=newi->parentg;  
	// put item in the tmpg group
	newi->parentg=tmpg;

	if (newi->nexting!=newi) // if the group does not have size 1
	{ // remove the item from its current group
		newi->nexting->previousing=newi->previousing;
		newi->previousing->nexting=newi->nexting;
		oldgroup->items=oldgroup->items->nexting;
	}
	else { // group will be empty
		if (oldgroup->nextg!=NULL) // there is another group
			oldgroup->nextg->previousg=oldgroup->previousg;
		if (lcu->root==oldgroup) // this is the first group
			lcu->root=oldgroup->nextg;
		else
			oldgroup->previousg->nextg=oldgroup->nextg;
		lcu->freegroups[--lcu->gpt]=oldgroup;
		// if we have created an empty group, remove it 
	}	
	newi->nexting=tmpg->items;
	newi->previousing=tmpg->items->previousing;
	newi->previousing->nexting=newi;
	newi->nexting->previousing=newi;
}

void LCU_AddNewGroupAfter(LCU_type * lcu, LCUITEM *newi, LCUGROUP *oldgroup) {
	LCUGROUP *newgroup;

	// remove item from old group...
	newi->nexting->previousing=newi->previousing;
	newi->previousing->nexting=newi->nexting;
	oldgroup->items=newi->nexting;
	//get new group
	newgroup=lcu->freegroups[lcu->gpt++];
	newgroup->count=oldgroup->count+1; // set count to be one more the prev group
	newgroup->items=newi;
	newgroup->previousg=oldgroup;
	newgroup->nextg=oldgroup->nextg;
	oldgroup->nextg=newgroup;
	if (newgroup->nextg!=NULL) // if there is another group
		newgroup->nextg->previousg=newgroup;
	newi->parentg=newgroup;
	newi->nexting=newi;
	newi->previousing=newi;
}

void LCU_IncrementCounter(LCU_type * lcu, LCUITEM *newi)
{
	LCUGROUP *oldgroup;

	oldgroup=newi->parentg;
	if ((oldgroup->nextg!=NULL) && 
		(oldgroup->nextg->count - oldgroup->count==1))
		LCU_PutInNewGroup(lcu, newi,oldgroup->nextg);
	// if the next group exists
	else { // need to create a new group with a differential of one
		if (newi->nexting==newi) // if there is only one item in the group...
			newi->parentg->count++;
		else      
			LCU_AddNewGroupAfter(lcu,newi,oldgroup);
	}
}

void LCU_Update(LCU_type * lcu, LCLitem_t newitem) {
	int h;
	LCUITEM *il;

	lcu->n++;
	h=hash31(lcu->a,lcu->b,newitem) % lcu->tblsz;
	il=lcu->hashtable[h];
	while (il) {
		if (il->item == newitem) 
			break;
		il=il->nexti;
	}
	if (il==NULL) // item is not monitored (not in hashtable) 
	{
		il=LCU_GetNewCounter(lcu);
		/// and put it into the hashtable for the new item 
		il->delta=lcu->root->count;
		// initialize delta with count of first group
		LCU_InsertIntoHashtable(lcu,il,h,newitem);
		// put the new counter into the first group
		// counter is already in first group by defn of how we got it
		LCU_IncrementCounter(lcu, il);
	}
	else 
		LCU_IncrementCounter(lcu, il);
	// if we have an item, we need to increment its counter 
}

int LCU_Size(LCU_type * lcu) {
	return sizeof(LCU_type)+(lcu->tblsz)*sizeof(LCUITEM*) + 
		(lcu->k)*(sizeof(LCUITEM) + sizeof(LCUGROUP) + sizeof(LCUITEM*));
}

void LCU_Destroy(LCU_type * lcu)
{
	free(lcu->freegroups);
	free(lcu->items);
	free(lcu->groups);
	free(lcu->hashtable);
	free (lcu);
} 

int LCU_PointEstUpp(LCU_type * lcu, LCLitem_t newitem) {
	int h;
	LCUITEM *il;

	lcu->n++;
	h=hash31(lcu->a,lcu->b,newitem) % lcu->tblsz;
	il=lcu->hashtable[h];
	while (il) {
		if (il->item == newitem) 
			return il->parentg->count;
		il=il->nexti;
	}
	return (lcu->root != NULL ? lcu->root->count : 0);
}

int thismax(int a, int b) {return (a >= b ? a : b);}

int LCU_PointEstLow(LCU_type * lcu, LCLitem_t newitem) {
	int h;
	LCUITEM *il;

	lcu->n++;
	h=hash31(lcu->a,lcu->b,newitem) % lcu->tblsz;
	il=lcu->hashtable[h];
	while (il) {
		if (il->item == newitem) 
			return thismax(il->parentg->count - il->delta, 0);
		il=il->nexti;
	}
	return 0;
}


