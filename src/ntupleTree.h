#ifndef  ntupleTree_h
#define  ntupleTree_h

#include "TROOT.h"
#include "TMath.h"
#include <vector>
#include <string>



class GenParticleCand {
public:
  Int_t   pdgId; 
  Int_t   status; 
  Float_t energy; 
  Float_t pt; 
  Float_t eta; 
  Float_t phi; 
  
  Float_t mumPt ;
  Float_t mumEta;
  Float_t mumPhi;
  Float_t mupPt ;
  Float_t mupEta;
  Float_t mupPhi;
  Float_t piPt  ;
  Float_t piEta ;
  Float_t piPhi ;
  Float_t kPt   ;
  Float_t kEta  ;
  Float_t kPhi  ;

  std::vector<Int_t>  pdgMother; 
  std::vector<Int_t>  pdgRealMother; 

  GenParticleCand(){};
  virtual ~GenParticleCand(){};
  
  ClassDef(GenParticleCand,1)
};




class B0Cand {
public:

  Float_t  Mu1Pt  ;  
  Float_t  Mu2Pt  ;  
  Float_t  Mu1Eta ;  
  Float_t  Mu1Phi ;  
  Float_t  Mu2Eta ;  
  Float_t  Mu2Phi ;  
  Int_t    Mu1Ch  ;  
  Int_t    Mu2Ch  ;  

  Float_t MuMuMass       ;  
  Float_t MuMuCL         ;  
  Float_t JpsiPosition_x ;  
  Float_t JpsiPosition_y ;  
  Float_t JpsiPosition_z ;  
  Float_t JpsiPt         ;  
  Float_t JpsiL          ;  
  Float_t JpsiSigma      ;  
  Float_t JpsiCosBS      ;  

  Float_t TrkPPt         ;  
  Float_t TrkPEta        ;  
  Float_t TrkPPhi        ;  
//   Float_t TrkPHQ         ;  
  Float_t TrkPd0Sign     ;  
  Float_t TrkMPt         ;  
  Float_t TrkMEta        ;  
  Float_t TrkMPhi        ;  
//   Float_t TrkMHQ         ;  
  Float_t TrkMd0Sign     ;  


  Float_t KStarMass        ;  
  Float_t barKStarMass     ;  
  Float_t KStarPt          ;  
  Float_t KStarEta         ;  
  Float_t KStarPhi         ;  
  Float_t B0Mass           ;  
  Float_t barB0Mass        ;  
  Float_t B0StarPt         ;  
  Float_t B0StarEta        ;  
  Float_t B0StarPhi        ;  

  Float_t JpsiTkTkPosition_x        ;  
  Float_t JpsiTkTkPosition_y        ;  
  Float_t JpsiTkTkPosition_z        ;  
  Float_t JpsiTkTkCL                ;  
  Float_t JpsiTkTkL                 ;  
  Float_t JpsiTkTkSigma             ;  
  Float_t JpsiTkTkCosBS             ;  
  
  Float_t cosThetaK                 ;
  Float_t cosThetaKErr              ;
  Float_t cosThetaMup               ;
  Float_t cosThetaMupErr            ;
  Float_t cosThetaMum               ;
  Float_t cosThetaMumErr            ;
  Float_t cosThetaK_Bbar            ;
  Float_t cosThetaKErr_Bbar         ;
  Float_t phiKstMuMuPlane           ;
  Float_t phiKstMuMuPlane_Bbar      ;
  
  
  

  B0Cand(){};
  virtual ~B0Cand(){};

  ClassDef(B0Cand,1)
};


class BpCand {
public:

  Float_t  Mu1Pt  ;  
  Float_t  Mu2Pt  ;  
  Float_t  Mu1Eta ;  
  Float_t  Mu1Phi ;  
  Float_t  Mu2Eta ;  
  Float_t  Mu2Phi ;  
  Int_t    Mu1Ch  ;  
  Int_t    Mu2Ch  ;  

  Float_t MuMuMass       ;  
  Float_t MuMuCL         ;  
  Float_t JpsiPosition_x ;  
  Float_t JpsiPosition_y ;  
  Float_t JpsiPosition_z ;  
  Float_t JpsiPt         ;  
  Float_t JpsiL          ;  
  Float_t JpsiSigma      ;  
  Float_t JpsiCosBS      ;  

  Float_t TrkPt             ;  
  Float_t TrkEta            ;  
  Float_t TrkPhi            ;  
  Float_t Trkd0Sign         ;  
  Int_t   TrkAlgo           ;  
  Float_t TrkIPFromJpsi     ;  
  Float_t TrkIPSignFromJpsi ;  

  Float_t BMass       ;  
  Float_t BPt         ;  
  Float_t BEta        ;  
  Float_t BPhi        ;  

  Float_t JpsiTkPosition_x        ;  
  Float_t JpsiTkPosition_y        ;  
  Float_t JpsiTkPosition_z        ;  
  Float_t JpsiTkCL                ;  
  Float_t JpsiTkL                 ;  
  Float_t JpsiTkSigma             ;  
  Float_t JpsiTkCosBS             ;  

  BpCand(){};
  virtual ~BpCand(){};

  ClassDef(BpCand,1)
};


class L1MuonCand {
public:

  Float_t pt;           
  Float_t eta;          
  Float_t phi;          
  Float_t etaAtVtx;          
  Float_t phiAtVtx;          
  Int_t   charge;      
  Int_t   quality;      
  
  L1MuonCand(){};
  virtual ~L1MuonCand(){};

  ClassDef(L1MuonCand,1)

};


class HLTMuCand {
public:

  Float_t  pt     ;  
  Float_t  eta    ;  
  Float_t  phi    ;  
  Int_t    charge ;  

  HLTMuCand(){};
  virtual ~HLTMuCand(){};

  ClassDef(HLTMuCand,1)

};


class HLTTkCand {
public:

  Float_t  pt     ;  
  Float_t  eta    ;  
  Float_t  phi    ;  
  Int_t    charge ;  
  Float_t  d0Sign ;  

  Float_t  CL      ;  
  Float_t  LSigma  ;  
  Float_t  CosBS   ;  

  HLTTkCand(){};
  virtual ~HLTTkCand(){};

  ClassDef(HLTTkCand,1)

};


class HLTDimuonCand {
public:

  Float_t  DCA     ;  
  Float_t  mu1pt   ;  
  Float_t  mu2pt   ;  
  Float_t  mu1eta  ;  
  Float_t  mu2eta  ;  
  Float_t  mu1phi  ;  
  Float_t  mu2phi  ;  
  Float_t  charge  ;  

  HLTDimuonCand(){};
  virtual ~HLTDimuonCand(){};

  ClassDef(HLTDimuonCand,1)

};



class HLTMuMuVtxCand {
public:

  Float_t  mu1pt   ;  
  Float_t  mu2pt   ;  

  Float_t  CL      ;  
  Float_t  x       ;  
  Float_t  y       ;  
  Float_t  z       ;  
  Float_t  ex      ;  
  Float_t  ey      ;  
  Float_t  ez      ; 
   
  HLTMuMuVtxCand(){};
  virtual ~HLTMuMuVtxCand(){};

  ClassDef(HLTMuMuVtxCand,1)

};


class HLTMuMuTkVtxCand {
public:

  Float_t  mu1pt  ;  
  Float_t  mu2pt  ;  
  Float_t  tkpt   ;  

  Float_t  CL      ;  
  Float_t  normchi2;  

  Float_t  x       ;  
  Float_t  y       ;  
  Float_t  z       ;  
  Float_t  ex      ;  
  Float_t  ey      ;  
  Float_t  ez      ;  

  HLTMuMuTkVtxCand(){};
  virtual ~HLTMuMuTkVtxCand(){};

  ClassDef(HLTMuMuTkVtxCand,1)

};


class HLTObjCand {
public:

  std::string filterTag; // name of filter passed by the object
  Float_t pt;            // pt of the object passing the filter [GeV]
  Float_t eta;           // eta of the object passing the filter
  Float_t phi;           // phi of the object passing the filter
  Int_t   id;            // trigger id of the object passing the filter
  
  HLTObjCand(){};
  virtual ~HLTObjCand(){};

  ClassDef(HLTObjCand,1)

};


class HLTInfo {
public:
  std::vector<std::string>  triggers;  
  std::vector<HLTObjCand>   objects;   

  HLTInfo(){};
  virtual ~HLTInfo(){};
  bool match( const std::string & path ) {
	if (  std::find (  triggers.begin(), triggers.end(), path ) != triggers.end() )  return true;
//     if (! iname.compare("HLT_Mu20_v1") == 0) continue;
	return false;
  }

  bool find( const std::string & path ) {
	for ( std::vector<std::string>::const_iterator it = triggers.begin(); it != triggers.end(); ++it ) {
// 	  std::cout << *it << std::endl;
      if ( it-> compare(path) == 0) return true;
//       if ( it->find ( path ) != std::string::npos ) return true;
	}
	return false;
  }

  ClassDef(HLTInfo,1)

};


class ntupleEvent {
public:

  Int_t   runNumber;             
  Int_t   luminosityBlockNumber; 
  Int_t   eventNumber;           

  unsigned int   prescale_novtx;
  unsigned int   prescale_vtx;
  unsigned int   rej_by_presc_vtx;
  unsigned int   rej_by_presc_novtx;
  Int_t   nVtx;                    
  Int_t   nTrks;   
  Float_t trueNI;   

  Float_t bxId;
  Float_t instLumi; 
  
  Float_t primaryVertex[3];        
  Float_t cov_primaryVertex[3][3]; 


  std::vector <GenParticleCand>     genParticles; 
  std::vector <B0Cand>              b0cands;         
  std::vector <BpCand>              bpcands;         

  std::vector <HLTMuCand>           hlt_mu;      
  std::vector <HLTTkCand>           hlt_tk;      
  std::vector <HLTTkCand>           hlt_glbtk;      
  std::vector <HLTTkCand>           hlt_newtk;      
  std::vector <HLTTkCand>           hlt_pix_tk;      
  std::vector <HLTDimuonCand>       hlt_dimu;
  std::vector <HLTMuMuVtxCand>      hlt_muvtx;      
  std::vector <HLTMuMuTkVtxCand>    hlt_tkvtx;      

  std::vector <L1MuonCand>          L1muons;      
     
  HLTInfo                           hlt;           
  HLTInfo                           hltTag;           

  ntupleEvent(){};
  virtual ~ntupleEvent(){};

  ClassDef(ntupleEvent,1)
};


#endif

