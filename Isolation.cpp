#include <iostream>

#include "TLorentzVector.h"
#include "TMath.h"

using namespace std;

bool IsIsolated(TLorentzVector pl, TLorentzVector *particles, Int_t n) {
  Float_t Rmax = 0.1; 
  Float_t Ptmin = 0.1;

  Float_t IP=0.0;

  for(Int_t i=0; i<n; i++) {
    if (pl != particles[i]) {
      if (particles[i].Pt()>Ptmin && pl.DeltaR(particles[i])<Rmax) {
         IP = IP + particles[i].Pt();
      }
    }
  }

  IP = IP/pl.Pt();
//  cout<<"IP="<<IP<<endl;
  if (IP<0.1) return 1;
  return 0;
}
