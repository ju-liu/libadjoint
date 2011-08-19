#ifndef ADJ_REVOLVE_H
#define ADJ_REVOLVE_H

#ifdef __cplusplus
#include <iostream>
#include "revolve.h"
#endif

typedef struct 
{
  void *ptr;
} CRevolve;

typedef enum 
{ 
  advance, takeshot, restore, firsturn, youturn, terminate, error
} CACTION;

#ifdef __cplusplus
extern "C" 
#endif
CRevolve revolve_create_offline(int st, int sn);

#ifdef __cplusplus
extern "C" 
#endif
CRevolve revolve_create_multistage(int st, int sn, int sn_ram);

#ifdef __cplusplus
extern "C" 
#endif
CRevolve revolve_create_online(int sn);

#ifdef __cplusplus
extern "C" 
#endif
void revolve_destroy(CRevolve r);

#ifdef __cplusplus
extern "C" 
#endif
CACTION revolve(CRevolve r); 

#ifdef __cplusplus
extern "C" 
#endif
int revolve_adjust(CRevolve r, int steps);

#ifdef __cplusplus
extern "C" 
#endif
int revolve_getadvances(CRevolve r);

#ifdef __cplusplus
extern "C" 
#endif
int revolve_getcheck(CRevolve r);

#ifdef __cplusplus
extern "C" 
#endif
int revolve_getcheckram(CRevolve r);

#ifdef __cplusplus
extern "C" 
#endif
int revolve_getcheckrom(CRevolve r);

#ifdef __cplusplus
extern "C" 
#endif
int revolve_getcapo(CRevolve r);

#ifdef __cplusplus
extern "C" 
#endif
int revolve_getfine(CRevolve r);

#ifdef __cplusplus
extern "C" 
#endif
int revolve_getinfo(CRevolve r);

#ifdef __cplusplus
extern "C" 
#endif
int revolve_getoldcapo(CRevolve r);

#ifdef __cplusplus
extern "C" 
#endif
int revolve_getwhere(CRevolve r);

#ifdef __cplusplus
extern "C" 
#endif
int revolve_setinfo(CRevolve r, int inf);

#ifdef __cplusplus
extern "C" 
#endif
const char* revolve_caction_string(CACTION action);

#endif
