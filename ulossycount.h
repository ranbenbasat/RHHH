// lossycount.h -- header file for Lossy Counting
// see Manku & Motwani, VLDB 2002 for details
// implementation by Graham Cormode, 2002,2003

#ifndef LOSSYCOUNTING_h
#define LOSSYCOUNTING_h

#include "prng.h"

#define LCL_NULLITEM 0x7FFFFFFF
//////////////////////////////////////////////////////
typedef int LCUWT;
//#define LCU_SIZE 101
// if not defined, then it is dynamically allocated based on user parameter
//////////////////////////////////////////////////////

#define LCU_HASHMULT 3
#ifdef LCU_SIZE
#define LCU_TBLSIZE (LCU_HASHMULT*LCU_SIZE)
#endif

#ifdef DIMENSION2
#define LCLitem_t uint64__t
#else
#define LCLitem_t uint32_t
#endif

typedef struct lcu_itemlist LCUITEMLIST;
typedef struct lcu_group LCUGROUP;

typedef struct lcu_item LCUITEM;
//typedef struct lcu_group LCUGROUP;

struct lcu_group 
{
  LCUWT count;
  LCUITEM *items;
  LCUGROUP *previousg, *nextg;
};

struct lcu_item 
{
  LCLitem_t item;
  int hash;
  LCUWT delta;
  LCUGROUP *parentg;
  LCUITEM *previousi, *nexti;
  LCUITEM *nexting, *previousing;
};

typedef struct LCU_type{

  LCUWT n;
  int gpt;
  int k;
  int tblsz;
  long long a,b;
  LCUGROUP * root;
#ifdef LCU_SIZE
  LCUITEM items[LCU_SIZE];
  LCUGROUP groups[LCU_SIZE];
  LCUGROUP *freegroups[LCU_SIZE];
  LCUITEM* hashtable[LCU_TBLSIZE];
#else
  LCUITEM * items;
  LCUGROUP *groups;
  LCUGROUP **freegroups;
  LCUITEM **hashtable;

#endif

} LCU_type;

extern LCU_type * LCU_Init(float fPhi);
extern void LCU_Destroy(LCU_type *);
extern void LCU_Update(LCU_type *, LCLitem_t);
extern int LCU_Size(LCU_type *);
int LCU_PointEstUpp(LCU_type *, LCLitem_t);
int LCU_PointEstLow(LCU_type *, LCLitem_t);
//extern std::map<uint32_t, uint32_t> LCU_Output(LCU_type *,int);

#endif
