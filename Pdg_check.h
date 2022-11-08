#ifndef PDG_CHECK_H
#define PDG_CHECK_H

#include "HepMC/GenEvent.h"

using namespace HepMC;

bool Is_q_meson(int pdg_id, int q);

bool Is_q_baryon(int pdg_id, int q);

bool Is_q_hadron(int pdg_id, int q);

std::set<int> Find_all_objects(int q, GenEvent *evt);

#endif
