#ifndef ISOLATION_H
#define ISOLATION_H

#include "TLorentzVector.h"

bool IsIsolated(TLorentzVector pl, TLorentzVector *particles, Int_t n);

#endif
