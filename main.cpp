#include <iostream>
#include <sstream>
#include <string>

#include "Isolation.h"
#include "Pdg_check.h"
#include "Btag.h"

#include "TFile.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TMath.h"

#include "HepMC/IO_GenEvent.h"
#include "HepMC/GenEvent.h"

#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"

using namespace std;
using namespace HepMC;
using namespace fastjet;

int main(int argc, char* argv[]) {

if(argc != 5) {
  cout<<"Usage: main [input_file] [output_file] [R] [E]"<<endl;
  cout<<"input_file -- HepMC2 format"<<endl;
  cout<<"output_file -- ROOT format"<<endl;
  cout<<"R -- the R parameter of kt algorithm"<<endl;
  cout<<"E -- the energy cut for jets"<<endl;
  return 0;
}

IO_GenEvent ascii_in(argv[1], std::ios::in);

Int_t nparticles;

Int_t nleptons;
Int_t njets;
Int_t nbjets;
Int_t ncjets;
Int_t nnus;

Float_t ptj[999];
Float_t ptl[999];
Float_t ptnu[999];

Float_t enerj[999];
Float_t enerl[999];
Float_t enernu[999];

Float_t mj[999];

Bool_t btag[999];
Bool_t ctag[999];
Bool_t tflag;

Float_t mjj;
Float_t mj1;
Float_t mj2;
Float_t met, me;

Int_t hb;

TFile *outputFile;
outputFile = TFile::Open(argv[2], "RECREATE");
TTree *tree = new TTree("eetc", "root file");
tree -> Branch("njets", &njets, "njets/I");
tree -> Branch("nbjets", &nbjets, "nbjets/I");
tree -> Branch("ncjets", &ncjets, "ncjets/I");
tree -> Branch("nleptons", &nleptons, "nleptons/I");
tree -> Branch("nnus", &nnus, "nleptons/I");
tree -> Branch("ptj1", &ptj[0], "ptj1/F");
tree -> Branch("ptj2", &ptj[1], "ptj2/F");
tree -> Branch("ptj3", &ptj[2], "ptj3/F");
tree -> Branch("ptl1", &ptl[0], "ptl1/F");
tree -> Branch("ptl2", &ptl[1], "ptl2/F");
tree -> Branch("ptnu1", &ptnu[0], "ptnu1/F");
tree -> Branch("ptnu2", &ptnu[1], "ptnu2/F");
tree -> Branch("enerj1", &enerj[0], "enerj1/F");
tree -> Branch("enerj2", &enerj[1], "enerj2/F");
tree -> Branch("enerj3", &enerj[2], "enerj3/F");
tree -> Branch("enerl1", &enerl[0], "enerl1/F");
tree -> Branch("enerl2", &enerl[1], "enerl2/F");
tree -> Branch("enernu1", &enernu[0], "enernu1/F");
tree -> Branch("enernu2", &enernu[1], "enernu2/F");
tree -> Branch("mjj", &mjj, "mjj/F");
tree -> Branch("mj1", &mj1, "mj1/F");
tree -> Branch("mj2", &mj2, "mj2/F");
tree -> Branch("met", &met, "met/F");
tree -> Branch("me", &me, "me/F");
tree -> Branch("hb", &hb, "hb/I");

int i, j;
int entries;
int nEvents = 0;

stringstream ss1, ss2;
string Rstring;
string Estring;
double R;
double ecut;

double brate = 0.70;
double crate = 0.10;
double qrate = 0.01;

Rstring = argv[3];
ss1<<Rstring;
ss1>>R;

Estring = argv[4];
ss2<<Estring;
ss2>>ecut;

cout<<"input_file: "<<argv[1]<<endl;
cout<<"output_file: "<<argv[2]<<endl;
cout<<"R="<<argv[3]<<endl;
cout<<"Ecut="<<argv[4]<<endl;

GenEvent* evt = ascii_in.read_next_event();

TLorentzVector pvec[999];

TLorentzVector pj[999];
TLorentzVector pl[999];
TLorentzVector pnu[999];
TLorentzVector pnusum;

TLorentzVector plcan;
TLorentzVector temp[3];

while(evt) {//start reading events
  nEvents++;
  vector<PseudoJet> particles;

  nparticles = 0;
  njets = 0;
  nbjets = 0;
  ncjets = 0;
  nleptons = 0;
  nnus = 0;
  mjj = 0.0;
  mj1 = 0.0;
  mj2 = 0.0;
  met = 0.0;
  me = 0.0;
  hb = 0;

  temp[0].SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
  temp[1].SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
  temp[2].SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
  pnusum.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
  plcan.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);

  tflag = false;

  for(i=0; i<999; i++) {
    pvec[i].SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
    pj[i].SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
    pl[i].SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
    pnu[i].SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
    ptj[i] = 0.0;
    ptl[i] = 0.0;
    ptnu[i] = 0.0;
    enerj[i] = 0.0;
    enerl[i] = 0.0;
    enernu[i] = 0.0;
    mj[i] = 0.0;
    btag[i] = false;
    ctag[i] = false;
  }

  //loop all particles
  for(GenEvent::particle_const_iterator p = evt->particles_begin(); p != evt->particles_end(); ++p) {
    if((*p)->status()==1) {
      pvec[nparticles].SetPxPyPzE((*p)->momentum().px(), (*p)->momentum().py(), (*p)->momentum().pz(), (*p)->momentum().e());
      nparticles++;
    }
  }
  //end of looping all particles

  //tag objects originated from b/c hadrons
  std::set<int> bobjects = Find_all_objects(5, evt); 
  std::set<int> cobjects = Find_all_objects(4, evt); 

  class BC4FJ:public PseudoJet::UserInfoBase {
    public: 
      BC4FJ(int barcode): _barcode(barcode) {};
      int _barcode;
  };

  //find isolated lepton and neutrino
  for(GenEvent::particle_const_iterator p = evt->particles_begin(); p != evt->particles_end(); ++p) {
    if((*p)->status()==1) {
      //find isolated leptons
      if((TMath::Abs((*p)->pdg_id())==11 || TMath::Abs((*p)->pdg_id())==13 || TMath::Abs((*p)->pdg_id())==15) && (*p)->momentum().perp()>20.0) { 
        plcan.SetPxPyPzE((*p)->momentum().px(), (*p)->momentum().py(), (*p)->momentum().pz(), (*p)->momentum().e());
        if(IsIsolated(plcan, pvec, nparticles)) {
          pl[nleptons] = plcan;
          ptl[nleptons] = pl[nleptons].Pt();
          enerl[nleptons] = pl[nleptons].E();
          nleptons++;
        }
        else {
          //If a lepton is not isolated, add it to the particles list for running fastjet
          particles.push_back(fastjet::PseudoJet((*p)->momentum().px(), (*p)->momentum().py(), (*p)->momentum().pz(), (*p)->momentum().e()));  
          //record the barcode
          BC4FJ* info = new BC4FJ((*p)->barcode());
          particles.back().set_user_info(info);
        }
      }
      //find neutrino
      else if((TMath::Abs((*p)->pdg_id())==12 || TMath::Abs((*p)->pdg_id())==14 || TMath::Abs((*p)->pdg_id())==16)) { 
        pnu[nnus].SetPxPyPzE((*p)->momentum().px(), (*p)->momentum().py(), (*p)->momentum().pz(), (*p)->momentum().e());
        ptnu[nnus] = pnu[nnus].Pt();
        enernu[nnus] = pnu[nnus].E();
        pnusum = pnusum + pnu[nnus];
        nnus++;
      }
      else {
        //remained particles are added to the list to run fastjet
        particles.push_back(fastjet::PseudoJet((*p)->momentum().px(), (*p)->momentum().py(), (*p)->momentum().pz(), (*p)->momentum().e())); 
        //record the barcode
        BC4FJ* info = new BC4FJ((*p)->barcode());
        particles.back().set_user_info(info);
      }
    }
  }
  //end of find isolated lepton and neutrion

  JetDefinition jet_def(ee_genkt_algorithm, R, 1);
  ClusterSequence cs_ee(particles, jet_def);
  vector<PseudoJet> jets_ee = sorted_by_E(cs_ee.inclusive_jets());

  entries = jets_ee.size();

  //jets
  for(i=0; i<entries; i++) {
    if(jets_ee[i].e()>ecut) {
      vector<PseudoJet> constituents = jets_ee[i].constituents();
      bool is_b_jet = false;
      bool is_c_jet = false;
      for (j=0; j<constituents.size(); j++) {
        assert(constituents[j].has_user_info());
        int barcode = dynamic_cast<const BC4FJ*>(constituents[j].user_info_ptr())->_barcode;
        if (bobjects.find(barcode) != bobjects.end()) {
           is_b_jet = true;
           is_c_jet = false;
           break;
        } 
        if (cobjects.find(barcode) != cobjects.end()) {
           is_b_jet = false;
           is_c_jet = true;
        }
      }
      if (is_b_jet) {
        btag[i] = true; 
        ctag[i] = false;
        nbjets++;
//        cout<<"find b-jet: "<<i<<endl;
      }
      if (is_c_jet) { 
        btag[i] = false;
        ctag[i] = true;
        ncjets++; 
//        cout<<"find c-jet: "<<i<<endl;
      }

      pj[njets].SetPxPyPzE(jets_ee[i].px(), jets_ee[i].py(), jets_ee[i].pz(), jets_ee[i].e());
      ptj[njets] = pj[i].Pt();
      enerj[njets] = pj[i].E();
      mj[njets] = pj[i].M();
      njets++;
    }
  }

  if(mj[1]>mj[0]) {
    temp[0] = pj[0]; pj[0] = pj[1]; pj[1] = temp[0];
    tflag = btag[0]; btag[0] = btag[1]; btag[1] = tflag;
    tflag = ctag[0]; ctag[0] = ctag[1]; ctag[1] = tflag;
  }

  if (btag[0] && Btag(brate)) hb=1; //cout<<"The heaviest jet is b-tagging jet. "<<endl;
  else if (ctag[0] && Btag(crate)) hb=1; //cout<<"The heaviest jet is c-tagging jet. "<<endl;
  else if (Btag(qrate)) hb=1;
  else hb=0;

//  cout<<"hb= "<<hb<<endl;

//  if (btag[1] && Btag(brate)) cout<<"The sub-heaviest jet is b-tagging jet. "<<endl;
//  if (ctag[1] && Btag(crate)) cout<<"The sub-heaviest jet is c-tagging jet. "<<endl;

  ptj[0] = pj[0].Pt();
  ptj[1] = pj[1].Pt();
  mj1 = pj[0].M();
  mj2 = pj[1].M();
  mjj = (pj[0]+pj[1]).M();

  //leptons
  if(nleptons>1) {
    for(i=0; i<nleptons-1; i++) {
      for(j=0; j<nleptons-i; j++) {
        if(ptl[j] < ptl[j+1]) {
          temp[1] = pl[j];
          pl[j] = pl[j+1];
          pl[j+1] = temp[1];
        }
      }
    }
  }

  ptl[0] = pl[0].Pt();
  ptl[1] = pl[1].Pt();
  enerl[0] = pl[0].E();
  enerl[1] = pl[1].E();

  //neutrinos
  if(nnus>1) {
    for(i=0; i<nnus-1; i++) {
      for(j=0; j<nnus-i; j++) {
        if(ptnu[j] < ptnu[j+1]) {
          temp[2] = pnu[j];
          pnu[j] = pnu[j+1];
          pnu[j+1] = temp[2];
        }
      }
    }
  }

  ptnu[0] = pnu[0].Pt();
  ptnu[1] = pnu[1].Pt();
  enernu[0] = pnu[0].E();
  enernu[1] = pnu[1].E();
  met = pnusum.Pt();
  me = pnusum.E();

  tree->Fill();
  delete evt;
  ascii_in >> evt;
}//end of reading events

tree->Write();
outputFile->Close();
return 0;
}//end of main
