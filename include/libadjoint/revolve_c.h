#ifdef __cplusplus
#include <iostream>
#include "revolve.h"
#endif

#ifndef __cplusplus
struct Revolve;
#endif

#ifdef __cplusplus
extern "C" void* revolve_create_offline(int st, int sn);
#else
void* revolve_create_offline(int st, int sn); 
#endif

#ifdef __cplusplus
extern "C" void* revolve_create_multistage(int st, int sn, int sn_ram);
#else
void* revolve_create_multistage(int st, int sn, int sn_ram); 
#endif

#ifdef __cplusplus
extern "C" void* revolve_create_online(int sn);
#else
void* revolve_create_online(int sn);
#endif

#ifdef __cplusplus
extern "C" void revolve_destroy(Revolve *r);
#else
void revolve_destroy(Revolve *r);
#endif

#ifdef __cplusplus
extern "C" int revolve_adjust(Revolve *r, int steps);
#else
int revolve_adjust(struct Revolve *r, int steps); 
#endif

#ifdef __cplusplus
extern "C" int revolve_getadvances(Revolve *r);
#else
int revolve_getadvances(struct Revolve *r); 
#endif

#ifdef __cplusplus
extern "C" int revolve_getcheck(Revolve *r);
#else
int revolve_getcheck(struct Revolve *r); 
#endif

#ifdef __cplusplus
extern "C" int revolve_getcheckram(Revolve *r);
#else
int revolve_getcheckram(struct Revolve *r); 
#endif

#ifdef __cplusplus
extern "C" int revolve_getcheckrom(Revolve *r);
#else
int revolve_getcheckrom(struct Revolve *r); 
#endif

#ifdef __cplusplus
extern "C" int revolve_getcapo(Revolve *r);
#else
int revolve_getcapo(struct Revolve *r); 
#endif

#ifdef __cplusplus
extern "C" int revolve_getfine(Revolve *r);
#else
int revolve_getfine(struct Revolve *r); 
#endif

#ifdef __cplusplus
extern "C" int revolve_getinfo(Revolve *r);
#else
int revolve_getinfo(struct Revolve *r); 
#endif


#ifdef __cplusplus
extern "C" int revolve_getoldcapo(Revolve *r);
#else
int revolve_getoldcapo(struct Revolve *r); 
#endif

#ifdef __cplusplus
extern "C" bool revolve_getwhere(Revolve *r);
#else
bool revolve_getwhere(struct Revolve *r); 
#endif

#ifdef __cplusplus
extern "C" int revolve_setinfo(Revolve *r, int inf);
#else
int revolve_setinfo(struct Revolve *r, int inf); 
#endif
