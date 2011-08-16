#include "../include/libadjoint/revolve_c.h"

extern "C" void* revolve_create_offline(int st, int sn) { return new Revolve(st, sn); }
extern "C" void* revolve_create_multistage(int st, int sn, int sn_ram) { return new Revolve(st, sn, sn_ram); }
extern "C" void* revolve_create_online(int sn) { return new Revolve(sn); }
extern "C" void revolve_destroy(Revolve *r) { delete r; }

extern "C" int revolve_adjust(Revolve* r, int steps) { return r->adjust(steps); }



extern "C" int revolve_getadvances(Revolve *r) { return r->getadvances(); }
extern "C" int revolve_getcheck(Revolve *r) { return r->getcheck(); }
extern "C" int revolve_getcheckram(Revolve *r) { return r->getcheckram(); }
extern "C" int revolve_getcheckrom(Revolve *r) { return r->getcheckrom(); }
extern "C" int revolve_getcapo(Revolve *r)  { return r->getcapo(); }
extern "C" int revolve_getfine(Revolve *r)  { return r->getfine(); }
extern "C" int revolve_getinfo(Revolve *r)  { return r->getinfo(); }
extern "C" int revolve_getoldcapo(Revolve *r) { return r->getoldcapo(); }
extern "C" int revolve_getwhere(Revolve *r) 
{ 
  if (r->getwhere())
    return 1;
  else
    return 0;
}
extern "C" void revolve_set_info(Revolve *r, int inf) { r->set_info(inf); }

