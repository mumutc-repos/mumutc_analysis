#include <iostream>

#include <TMath.h>

#include "HepMC/GenEvent.h"

using namespace std;
using namespace HepMC;

bool Is_q_meson(int pdg_id, int q) {
  //return true if and only if the heavy quark is q.
  return ((TMath::Abs(pdg_id)/100)%100) == q;  
}

bool Is_q_baryon(int pdg_id, int q) {
  //return true if and only if the heaviest quark is q.
  return (TMath::Abs(pdg_id)/1000) == q;  
}

bool Is_q_hadron(int pdg_id, int q) {
  //return true if it is meson or baryon of quark q.
  return Is_q_meson(pdg_id, q) || Is_q_baryon(pdg_id, q);
}

std::set<int> Find_all_objects(int q, GenEvent *evt) {
  std::set<int> intermediate_state;
  std::set<int> final_state;

  for (GenEvent::particle_const_iterator p = evt->particles_begin(); p != evt->particles_end(); ++p) {
    if (Is_q_hadron((*p)->pdg_id(), q)) {
      cout<<"find: "<<q<<" hadron: "<<(*p)->pdg_id()<<endl;
      std::vector<GenParticle*> pvec;
      pvec.push_back((*p));
      while(!pvec.empty()) {
        GenParticle *gp = pvec.back();
        pvec.pop_back();
        if(intermediate_state.find(gp->barcode()) != intermediate_state.end()) {
          continue;
        }
        intermediate_state.insert(gp->barcode());
        if(gp->status() == 1) {
//          cout<<"Processing "<<gp->barcode()<<" PDG: "<<gp->pdg_id()<<endl;
          final_state.insert(gp->barcode());
          continue;
        }
        GenVertex *v = gp->end_vertex();
//        v->print();
        for (GenVertex::particles_out_const_iterator pout = v->particles_out_const_begin(); pout != v->particles_out_const_end(); ++pout) {
          pvec.push_back((*pout)); 
        }
      }
    }
  }
  return final_state;
}
