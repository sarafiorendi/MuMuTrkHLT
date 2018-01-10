#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TTree.h"
#include "TBranch.h"

#include <iostream>
#include <algorithm>
#include <vector>
#include <map>

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "MyTools/MuMuTrkHLT/src/ntupleTree.h"
#include "TLorentzVector.h"

double muonmass = 0.10565837;
double kaonmass = 0.493677;
double pionmass = 0.139570;

float  b0MassPDG    = 5.27962;

double pt_bins[11]     = { 0, 10, 12, 14, 16, 18, 20, 25, 30, 40, 50} ;
double pt_tk_bins[16]  = { 0, 0.8, 1., 1.2, 1.4, 1.6, 1.8, 2., 2.3, 2.6, 3., 3.5, 4., 5, 7, 10} ;
double ip_tk_bins[6]   = { 0, 0.005, 0.01, 0.015, 0.03, 1.5} ;


bool matchMuon         (double, double, double, std::vector<HLTObjCand>, std::string);
bool matchTrack        (double, double, double, std::vector<HLTObjCand>, std::string);
HLTTkCand matchTk      (double, double, std::vector<HLTTkCand>);
HLTMuCand matchMuon    (double, double, std::vector<HLTMuCand> );

std::string hltDen        = "HLT_DoubleMu4_3_Jpsi_Displaced_v";
std::string hltNum        = "HLT_DoubleMu4_JpsiTrk_Displaced_v";

std::string displacedJpsiFilter1 = "hltDisplacedmumuFilterDoubleMu43Jpsi::HLT";
std::string displacedJpsiFilter2 = "hltDoubleMu43JpsiDisplacedL3Filtered::HLT";
std::string hltTrkFilter         = "hltJpsiTkVertexFilter::HLT";


void readTree_Bp(){

  TChain* tree = new TChain("theNtuples/ntupleTree"); 
  for (int i = 1; i < 10 ; i++){
    tree -> Add(Form("/eos/cms/store/group/phys_bphys/fiorendi/13TeV/data2017/Charmonium/crab_bph_ntuples_passDoubleMu4_3_2016E_23Sep2016/171218_154404/0000/BHltNtuple_%i.root", i));
  }

  TFile* outfile = new TFile("efficiencyHLT_toDel.root","recreate");
  std::cout << "outfile: " << outfile -> GetName() << std::endl;

  float theBmass;
  bool pass;
  TTree* mTree = new TTree("mTree", "mTree");
  mTree -> Branch("bMass" , &theBmass      );
  mTree -> Branch("pass"  , &pass          );
  

  TH1F* b0_den              = new TH1F("b0_den"             ,"b0_den"          ,  200,   4.8,  5.6 );
  TH1F* b0_num              = new TH1F("b0_num"             ,"b0_num"          ,  200,   4.8,  5.6 );
  TH1F* b0_fail             = new TH1F("b0_fail"            ,"b0_fail"         ,  200,   4.8,  5.6 );
 
  TH1F* b0_den_eta          = new TH1F ("b0_den_eta"        , "b0_den_eta"     ,   25,  -2.5,  2.5 );
  TH1F* b0_den_mueta        = new TH1F ("b0_den_mueta"      , "b0_den_mueta"   ,   25,  -2.5,  2.5 );  
  TH2F* b0_den_etaphi       = new TH2F ("b0_den_etaphi"     , "b0_den_etaphi"  ,    6,  -2.4,  2.4, 16, -3.2, 3.2);  
  TH1F* b0_all_eta          = new TH1F ("b0_all_eta"        ,"b0_all_eta"      ,   25,  -2.5,  2.5 );

  TH1F* mumu_mass           = new TH1F("mumu_mass"          ,"mumu_mass"       , 1000,     2,    4 );
  TH1F* trkPt               = new TH1F("trkPt"              ,"trkPt"           ,  100,     0,   10 );

  TEfficiency* eff_eta  	= new TEfficiency( "eff_eta"    , "eff_eta"        ,   25,  -2.5,  2.5 );
  TEfficiency* eff_phi  	= new TEfficiency( "eff_phi"    , "eff_phi"        ,   16,  -3.2,  3.2 );
  TEfficiency* eff_trkPt   	= new TEfficiency( "eff_trkPt"  , "eff_trkPt"      ,   15, pt_tk_bins  );
  TEfficiency* eff_pt   	= new TEfficiency( "eff_pt"     , "eff_pt"         ,   10, pt_bins     );
  TEfficiency* eff_nvtx  	= new TEfficiency( "eff_nvtx"   , "eff_nvtx"       ,   20,     5,   65 );
  TEfficiency* eff_jpsiPt   = new TEfficiency( "eff_jpsiPt" , "eff_jpsiPt"     ,   10, pt_bins     );

  TEfficiency* eff_tkEta  	= new TEfficiency( "eff_tkEta"  , "eff_tkEta"      ,   25,  -2.5,  2.5 );
  TEfficiency* eff_tkPhi 	= new TEfficiency( "eff_tkPhi"  , "eff_tkPhi"      ,   16,  -3.2,  3.2 );
  TEfficiency* eff_tkIP   	= new TEfficiency( "eff_tkIP"   , "eff_tkIP"       ,    5, ip_tk_bins  );
  
  
  
  ntupleEvent* ev      = new ntupleEvent();
  tree -> SetBranchAddress( "event", &ev);

  int nentries = tree->GetEntries();
  std::cout << "Number of entries = " << nentries << std::endl;
    
  //read ntuple tree 
  float themass = 0;
  for (Int_t i=0; i < nentries; i++)
  {
    Int_t IgetEvent   = tree   -> GetEvent(i);
    if (i % 50000 == 0) std::cout << "eventNo: " << i << std::endl;
//     if (ev -> runNumber < 306029)  continue;

    if (!ev-> hltTag.find(hltDen)) continue;

    unsigned int nb0  = ev->bpcands.size();
    for (int ibp = 0; ibp < nb0; ibp++){
    
      BpCand theBp = ev -> bpcands.at(ibp);
      b0_all_eta -> Fill(theBp.BEta);

      pass          = false;
      bool matchTrk = false;

      // offline selections
      if (theBp.BPt <= 8 )                          continue; 
      if (theBp.Mu1Pt  <= 4 || theBp.Mu2Pt <= 4 )   continue; 
      if (theBp.TrkPt <= 1.2 )                      continue; 
      if (theBp.Trkd0Sign <= 2 )                    continue; 
      if (theBp.JpsiTkCosBS <= 0.99)                continue; 
      if (theBp.JpsiTkL / theBp.JpsiTkSigma <= 3.)  continue; 
      if (theBp.JpsiTkCL < 0.1)                     continue; 
      
      // displaced HLT cuts
      if (theBp.JpsiPt < 7.)                    continue; 
      if (theBp.JpsiL / theBp.JpsiSigma < 3.)   continue;
      if (theBp.JpsiCosBS <= 0.9)               continue;
      if (theBp.MuMuCL < 0.1)                   continue;
      if (theBp.MuMuMass > 3.3)                 continue;
      if (theBp.MuMuMass < 2.9)                 continue;
      
    
      bool matchMu1a = matchMuon( theBp.Mu1Eta, theBp.Mu1Phi, theBp.Mu1Pt,ev -> hlt.objects, displacedJpsiFilter1);
      bool matchMu2a = matchMuon( theBp.Mu2Eta, theBp.Mu2Phi, theBp.Mu2Pt,ev -> hlt.objects, displacedJpsiFilter1);
      bool matchMu1b = matchMuon( theBp.Mu1Eta, theBp.Mu1Phi, theBp.Mu1Pt,ev -> hlt.objects, displacedJpsiFilter2);
      bool matchMu2b = matchMuon( theBp.Mu2Eta, theBp.Mu2Phi, theBp.Mu2Pt,ev -> hlt.objects, displacedJpsiFilter2);
      
      if (! (matchMu1a && matchMu2a) ) continue;
      if (! (matchMu1b && matchMu2b) ) continue;
      
      theBmass = theBp.BMass;

      b0_den        -> Fill( theBmass);
      b0_den_eta    -> Fill( theBp.BEta);
      b0_den_mueta  -> Fill( theBp.Mu1Eta);
      b0_den_mueta  -> Fill( theBp.Mu2Eta);
      b0_den_etaphi -> Fill( theBp.BEta, theBp.BPhi); 
      
      mumu_mass -> Fill( theBp.MuMuMass);
      // match offline trk with online
      matchTrk = matchTrack( theBp.TrkEta, theBp.TrkPhi, theBp.TrkPt, ev -> hlt.objects, hltTrkFilter);
      trkPt -> Fill(theBp.TrkPt);
      
      // numerator path is fired and track is matched to the last hlt filter
      pass =  matchTrk && ev-> hltTag.find(hltNum) ;  
      
      if ( pass) b0_num  -> Fill( theBmass);
      if (!pass) b0_fail -> Fill( theBmass);

      // "counting" efficiency
      if (theBmass > 5.16 && theBmass < 5.4) {
        eff_eta     -> Fill(pass, theBp.BEta  );
        eff_phi     -> Fill(pass, theBp.BPhi  );
        eff_pt      -> Fill(pass, theBp.BPt   );
        eff_tkEta   -> Fill(pass, theBp.TrkEta);
        eff_tkPhi   -> Fill(pass, theBp.TrkPhi);
        eff_trkPt   -> Fill(pass, theBp.TrkPt );
        eff_nvtx    -> Fill(pass, ev -> nVtx  );
        eff_tkIP    -> Fill(pass, theBp.TrkIPFromJpsi );
        eff_jpsiPt  -> Fill(pass, theBp.JpsiPt);
      }
    
      mTree    -> Fill();

    }
  }  
  

  outfile        -> cd();
  mTree           -> Write();

  b0_all_eta      -> Write();
  b0_den          -> Write();
  b0_den_eta      -> Write();
  b0_den_mueta    -> Write();
  b0_den_etaphi   -> Write();
  b0_num          -> Write();
  b0_fail         -> Write();
    
  mumu_mass   -> Write();
  trkPt       -> Write();

  eff_eta     -> Write();
  eff_phi     -> Write();
  eff_pt      -> Write();
  eff_tkEta   -> Write();
  eff_tkPhi   -> Write();
  eff_trkPt   -> Write();
  eff_tkIP    -> Write();
  eff_jpsiPt  -> Write();
  eff_nvtx    -> Write();
  
  
  outfile -> Close();
  tree    -> Delete();

}




bool matchTrack(double eta, double phi, double pt, std::vector<HLTObjCand> toc, std::string tagFilterName){

  bool match = false;
  int ntoc = toc.size();

  float minDR = 0.1;
  float theDR = 1;
  for ( std::vector<HLTObjCand>::const_iterator it = toc.begin(); it != toc.end(); ++it ) {
    if ( it->filterTag.compare(tagFilterName) == 0 && fabs(it -> id) == 321) {
      theDR = deltaR(it -> eta, it -> phi, eta, phi);
      if (theDR < minDR){
        minDR = theDR;
        match = true;
      }
    }
  }
  
  return match;
}



bool matchMuon(double eta, double phi, double pt, std::vector<HLTObjCand> toc, std::string tagFilterName){

  bool match = false;
  int ntoc = toc.size();

  float minDR = 0.1;
  float theDR = 1;
  for ( std::vector<HLTObjCand>::const_iterator it = toc.begin(); it != toc.end(); ++it ) {
    if ( it->filterTag.compare(tagFilterName) == 0) {
//       std::cout << it->filterTag << std::endl;
      theDR = deltaR(it -> eta, it -> phi, eta, phi);
      if (theDR < minDR){
        minDR = theDR;
        match = true;
      }
    }
  }
  return match;
}



HLTTkCand matchTk(double eta, double phi, std::vector<HLTTkCand> Tkcands){

  bool match = false;
  int nTk = Tkcands.size();

  float minDR = 0.1;
  float theDR = 100;
  HLTTkCand theTk;
  theTk.pt        = -1000;
  theTk.eta       = -1000;
  theTk.phi       = -1000;
  theTk.charge    = -1000;
  theTk.d0Sign    = -1000;

  for ( std::vector<HLTTkCand>::const_iterator it = Tkcands.begin(); it != Tkcands.end(); ++it ) {
    theDR = deltaR(it -> eta, it -> phi, eta, phi);
    if (theDR < minDR){
      minDR = theDR;
      match = true;
      theTk = *it;
    }
  }
  
  return theTk;
}



HLTMuCand matchMuon(double eta, double phi, std::vector<HLTMuCand> Muoncands){

  bool match = false;
  int nMuon = Muoncands.size();

  float minDR = 0.1;
  float theDR = 100;
  HLTMuCand theMuon;
  theMuon.pt        = -1000;
  theMuon.eta       = -1000;
  theMuon.phi       = -1000;
  theMuon.charge    = -1000;

  for ( std::vector<HLTMuCand>::const_iterator it = Muoncands.begin(); it != Muoncands.end(); ++it ) {
    theDR = deltaR(it -> eta, it -> phi, eta, phi);
    if (theDR < minDR){
      minDR = theDR;
      match = true;
      theMuon = *it;
    }
  }
  
  return theMuon;
}
