// ---------------------------------------------------------------------
void BHltNtuples::fillB0s  (const edm::Handle<reco::MuonCollection>       & muons ,
                            const edm::Handle<reco::TrackCollection >     & tracks,
                            const edm::Handle<reco::VertexCollection >    & vtxColl,
                            const edm::Event                              & event ,
                            const edm::EventSetup                         & eventSetup)
{

  const float  kStarMassPDG = 0.89166;

  const double thirdTrackMass2  = thirdTrackMass_ *thirdTrackMass_ ; // kaon
  const double fourthTrackMass2 = fourthTrackMass_*fourthTrackMass_; // pion

  //get offline beamspot position
  edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
  event.getByToken(beamspotToken_,recoBeamSpotHandle);
  const reco::BeamSpot& vertexBeamSpot = *recoBeamSpotHandle;

  //get the transient track builder:
  edm::ESHandle<TransientTrackBuilder> theB;
  eventSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);

  //get the b field
  std::string mfName_ = "";
  edm::ESHandle<MagneticField> bFieldHandle;
  eventSetup.get<IdealMagneticFieldRecord>().get("", bFieldHandle);  
  const MagneticField* magField = bFieldHandle.product();
  TSCBLBuilderNoMaterial blsBuilder;

  const reco::VertexCollection pvColl  = *(vtxColl.product())  ;    

  // preselection on tracks
  selTracksDef qualityTracksPlus, qualityTracksMinus;
  for (uint tracksIt =0 ;  tracksIt < tracks->size(); tracksIt++)
  {
    reco::TrackRef checkTrk(tracks,tracksIt) ;                                                
    reco::Track itrk  = tracks->at(tracksIt) ;                                                
    
    if (!checkTrk->quality(reco::TrackBase::highPurity))                continue;
    if (checkTrk->pt() < minPtTrk_ || fabs(checkTrk->eta())  > maxEta_) continue;            
    
    FreeTrajectoryState InitialFTS2 = initialFreeState(itrk, magField);
    TrajectoryStateClosestToBeamLine tscb2( blsBuilder(InitialFTS2, *recoBeamSpotHandle) );
    float trk_d0sig = tscb2.transverseImpactParameter().significance();
    if (trk_d0sig  < mind0Sign_)    continue;

    // check overlap with muon collection
    bool flag = false ;                                                         
    for (reco::MuonCollection::const_iterator mu =muons->begin(); mu!=muons->end(); mu++)                    
    {                                                                
      if (mu->track().get()!= 0 && mu->track().get() == checkTrk.get()) { flag=true; break; }                    
    }                                                                
    if (flag)   continue;            

    if (itrk.charge() == 1)  qualityTracksPlus[tracksIt] = itrk;     
    else                     qualityTracksMinus[tracksIt]= itrk;   
  }


  // Loop on muon collection
  for(std::vector<reco::Muon>::const_iterator mu1=muons->begin(); mu1!=muons->end(); ++mu1) 
  { 
    if( muon::isSoftMuon( (*mu1), pvColl[0] ) && (*mu1).pt() > 0 && fabs( (*mu1).eta() ) < 2.4) 
    {
      // Go on and look for the second muon
      for(std::vector<reco::Muon>::const_iterator mu2=mu1; mu2!=muons->end(); ++mu2) {
        if( mu2 == mu1 ) continue; 
        if( muon::isSoftMuon( (*mu2), pvColl[0]) && (*mu2).pt() > 0 && fabs( (*mu2).eta() ) < 2.4) 
        {
          if(!( mu1->charge() * mu2->charge() < 0 ))         continue; 
//           hists_["mumuMass_all"]->Fill( (mu1->p4() + mu2->p4()).mass() );

          // do jpsi vertex fit
          std::vector<reco::TransientTrack> j_tks;
          j_tks.push_back((*theB).build(mu1->track().get()));
          j_tks.push_back((*theB).build(mu2->track().get()));
          if (j_tks.size()!=2) continue;
       
          KalmanVertexFitter jkvf;
          TransientVertex jtv = jkvf.vertex(j_tks);
          if (!jtv.isValid()) continue;
        
          reco::Vertex jpsivertex = jtv;
          float dimuonCL = 0;
          if( (jpsivertex.chi2()>=0.0) && (jpsivertex.ndof()>0) ) 
            dimuonCL = TMath::Prob(jpsivertex.chi2(), jpsivertex.ndof() );
            
          // calculate three-track transverse momentum
          math::XYZVector jpperp(mu1->px() + mu2->px() ,
                                 mu1->py() + mu2->py() ,
                                 0.);
         
          GlobalPoint jVertex = jtv.position();
          GlobalError jerr    = jtv.positionError();
          
          //calculate decay length  significance w.r.t. the beamspot
          GlobalPoint displacementFromBeamspotJpsi( -1*((vertexBeamSpot.x0() - jVertex.x()) + (jVertex.z() - vertexBeamSpot.z0()) * vertexBeamSpot.dxdz()), 
                                                    -1*((vertexBeamSpot.y0() - jVertex.y()) + (jVertex.z() - vertexBeamSpot.z0()) * vertexBeamSpot.dydz()),
                                                     0);
          reco::Vertex::Point vperpj(displacementFromBeamspotJpsi.x(), displacementFromBeamspotJpsi.y(), 0.);

          float jpsi_ls  = displacementFromBeamspotJpsi.perp() /sqrt(jerr.rerr(displacementFromBeamspotJpsi));
          float jpsi_cos = vperpj.Dot(jpperp)/(vperpj.R()*jpperp.R());

//           hists_["JpsiPt"  ]->Fill( jpperp.R() );
//           hists_["JpsiLS"  ]->Fill( jpsi_ls    );
//           hists_["JpsiCos" ]->Fill( jpsi_cos   );
//           hists_["JpsiCL"  ]->Fill( dimuonCL   );
         
          if (jpperp.R()   < minJpsiPt_    ) continue;
          if (jpsi_ls      < minJpsiLS_    ) continue;
          if (jpsi_cos     < minJpsiCos_   ) continue;
          if (dimuonCL     < minJpsiCL_    ) continue;

          // Loop on track collection - trk 1
          for (selTracksDef::const_iterator tracksIt=qualityTracksPlus.begin(); tracksIt!=qualityTracksPlus.end(); ++tracksIt)
          {
            reco::Track itrk1((*tracksIt).second) ;
            if (overlap(*mu1,itrk1))    continue;
            if (overlap(*mu2,itrk1))    continue;
            
            FreeTrajectoryState InitialFTS = initialFreeState(itrk1, magField);
            TrajectoryStateClosestToBeamLine tscb( blsBuilder(InitialFTS, *recoBeamSpotHandle) );
            float trk1_d0sig = tscb.transverseImpactParameter().significance();

            if (trk1_d0sig < 1.) continue;
//             hists_["trkPt"] -> Fill( itrk1.pt() );
//             hists_["D0sig"] -> Fill( trk1_d0sig );
 
            for (selTracksDef::const_iterator tracksIt2 = qualityTracksMinus.begin(); tracksIt2 != qualityTracksMinus.end(); ++tracksIt2)
            {
              reco::Track itrk2((*tracksIt2).second) ;
              if (itrk2.charge()*itrk1.charge() != -1)            continue;
              if (overlap(*mu1,itrk2))                            continue;
              if (overlap(*mu2,itrk2))                            continue;

              FreeTrajectoryState InitialFTS2 = initialFreeState(itrk2, magField);
              TrajectoryStateClosestToBeamLine tscb2( blsBuilder(InitialFTS2, *recoBeamSpotHandle) );
              float trk2_d0sig = tscb2.transverseImpactParameter().significance();

              if (trk2_d0sig < 1.) continue;

//               hists_["trkPt"] -> Fill(itrk2.pt()   );
//               hists_["D0sig"] -> Fill(trk2_d0sig   );

              reco::Particle::LorentzVector pB, pbarB, p1, p2, p3_k, p4_p, p3_p, p4_k, pKstar, pKstarBar, pJpsi;

              // Combined system
              double e1   = sqrt(mu1->momentum().Mag2()  + MuMass2          );
              double e2   = sqrt(mu2->momentum().Mag2()  + MuMass2          );
              double e3_k = sqrt(itrk1.momentum().Mag2() + thirdTrackMass2  );
              double e3_p = sqrt(itrk1.momentum().Mag2() + fourthTrackMass2 );
              double e4_k = sqrt(itrk2.momentum().Mag2() + thirdTrackMass2  );
              double e4_p = sqrt(itrk2.momentum().Mag2() + fourthTrackMass2 );
            
              p1   = reco::Particle::LorentzVector(mu1->px() , mu1->py() , mu1->pz() , e1  );
              p2   = reco::Particle::LorentzVector(mu2->px() , mu2->py() , mu2->pz() , e2  );
              p3_k = reco::Particle::LorentzVector(itrk1.px(), itrk1.py(), itrk1.pz(), e3_k);
              p3_p = reco::Particle::LorentzVector(itrk1.px(), itrk1.py(), itrk1.pz(), e3_p);
              p4_k = reco::Particle::LorentzVector(itrk2.px(), itrk2.py(), itrk2.pz(), e4_k);
              p4_p = reco::Particle::LorentzVector(itrk2.px(), itrk2.py(), itrk2.pz(), e4_p);
            
              pB        = p1 + p2 + p3_k + p4_p;
              pbarB     = p1 + p2 + p3_p + p4_k;
              pKstar    = p3_k + p4_p;
              pKstarBar = p3_p + p4_k;
              
              pJpsi = p1 + p2;
            
              if (pB.mass() > 7. && pbarB.mass() > 7.) continue;
              if (pB.mass() < 4. && pbarB.mass() < 4.) continue;
              if (fabs(pKstar.mass() - kStarMassPDG) > 0.4 && fabs(pKstarBar.mass() - kStarMassPDG) > 0.4) continue;


              // do the vertex fit
              std::vector<reco::TransientTrack> t_tks;
              t_tks.push_back((*theB).build(mu1->track().get()));
              t_tks.push_back((*theB).build(mu2->track().get()));
              t_tks.push_back((*theB).build(&itrk1));
              t_tks.push_back((*theB).build(&itrk2));
              if (t_tks.size()!=4) continue;
            
              KalmanVertexFitter kvf;
              TransientVertex tv  = kvf.vertex(t_tks);
              reco::Vertex vertex = tv;
              if (!tv.isValid()) continue;
//               hists_["B0InvMass"]->Fill( pB.mass() );

              float JpsiTkTkCL = 0;
              if ((vertex.chi2()>=0.0) && (vertex.ndof()>0) )   
                JpsiTkTkCL = TMath::Prob(vertex.chi2(), vertex.ndof() );
              
              // calculate four-track transverse momentum
              math::XYZVector pperp(mu1->px() + mu2->px() + itrk1.px()+ itrk2.px(),
                                    mu1->py() + mu2->py() + itrk1.py()+ itrk2.py(),
                                    0.);
              // get vertex position and error to calculate the decay length significance
              GlobalPoint secondaryVertex = tv.position();
              GlobalError err             = tv.positionError();
              GlobalPoint displacementFromBeamspot( -1*((vertexBeamSpot.x0() - secondaryVertex.x()) + 
                                                        (secondaryVertex.z() - vertexBeamSpot.z0()) * vertexBeamSpot.dxdz()), 
                                                    -1*((vertexBeamSpot.y0() - secondaryVertex.y()) + 
                                                        (secondaryVertex.z() - vertexBeamSpot.z0()) * vertexBeamSpot.dydz()), 
                                                    0);
              reco::Vertex::Point vperp(displacementFromBeamspot.x(),displacementFromBeamspot.y(),0.);
 
 
              if (displacementFromBeamspot.perp() / sqrt(err.rerr(displacementFromBeamspot)) < 5.) continue;;
 
 
              B0Cand theB0;

 
			  // calculate angular vars for B0
			  {
			  TLorentzVector LoreVecB0;

			  TLorentzVector LoreVecMuMu;
			  LoreVecMuMu.SetXYZM(pJpsi.px(), pJpsi.py(), pJpsi.pz(), pJpsi.mass());
			  TVector3 boostMuMu = LoreVecMuMu.BoostVector();
			  TLorentzVector LoreVecMup;
			  if (mu1->charge() > 0)
				LoreVecMup.SetXYZM(mu1->px(),mu1->py(),mu1->pz(),MuMass);
			  else
				LoreVecMup.SetXYZM(mu2->px(),mu2->py(),mu2->pz(),MuMass);
			  // # Boost mu+ to mumu ref. frame #
			  LoreVecMup.Boost(-boostMuMu);
			  // # Boost B0 to mumu ref. frame #
			  LoreVecB0.SetXYZM(pB.px(), pB.py(), pB.pz(), pB.mass());
			  LoreVecB0.Boost(-boostMuMu);
			  // # Compute angles #
			  double cosThetaMup,  cosThetaMupErr;
			  computeCosAlpha(-LoreVecB0.Px(),-LoreVecB0.Py(),-LoreVecB0.Pz(),
						   LoreVecMup.Px(),LoreVecMup.Py(),LoreVecMup.Pz(),
						   0.0,0.0,0.0,
						   0.0,0.0,0.0,
						   0.0,0.0,0.0,
						   0.0,0.0,0.0,
						   &cosThetaMup,&cosThetaMupErr);

			  TLorentzVector LoreVecKst;
			  LoreVecKst.SetXYZM(pKstar.px(), pKstar.py(), pKstar.pz(),pKstar.mass());
			  TVector3 boostKst = LoreVecKst.BoostVector();
			  TLorentzVector LoreVecK;
			  LoreVecK.SetXYZM(p3_k.px(),p3_k.py(),p3_k.pz(),thirdTrackMass_);
	// 		  LoreVecK.SetXYZM(NTupleIn->kstTrkpPx->at(BestCandIndx),NTupleIn->kstTrkpPy->at(BestCandIndx),NTupleIn->kstTrkpPz->at(BestCandIndx),Utility->kaonMass);


			  // # Boost K+ to K*0 ref. frame #
			  LoreVecK.Boost(-boostKst);
			  // # Boost B0 to K*0 ref. frame #
			  LoreVecB0.SetXYZM(pB.px(), pB.py(), pB.pz(),pB.mass());
			  LoreVecB0.Boost(-boostKst);
			  // # Compute angles #
			  double cosThetaK,  cosThetaKErr, phiKstMuMuPlane;
			  computeCosAlpha(-LoreVecB0.Px(),-LoreVecB0.Py(),-LoreVecB0.Pz(),
						   LoreVecK.Px(),LoreVecK.Py(),LoreVecK.Pz(),
						   0.0,0.0,0.0,
						   0.0,0.0,0.0,
						   0.0,0.0,0.0,
						   0.0,0.0,0.0,
						   &cosThetaK,&cosThetaKErr);

			  // # Angle between [mu+ - mu-] and [K - pi] planes in the B0 ref. frame #
			  LoreVecB0.SetXYZM(pB.px(), pB.py(), pB.pz(),pB.mass());
			  TVector3 boostB0 = LoreVecB0.BoostVector();
			  if (mu1->charge() > 0)
				LoreVecMup.SetXYZM(mu1->px(),mu1->py(),mu1->pz(),MuMass);
			  else
				LoreVecMup.SetXYZM(mu2->px(),mu2->py(),mu2->pz(),MuMass);
			  TLorentzVector LoreVecMum;

			  if (mu1->charge() < 0)
				LoreVecMum.SetXYZM(mu1->px(),mu1->py(),mu1->pz(),MuMass);
			  else
				LoreVecMum.SetXYZM(mu2->px(),mu2->py(),mu2->pz(),MuMass);

			  LoreVecK.SetXYZM(p3_k.px(),p3_k.py(),p3_k.pz(),thirdTrackMass_);
			  TLorentzVector LoreVecPi;
			  LoreVecPi.SetXYZM(p4_p.px(),p4_p.py(),p4_p.pz(),fourthTrackMass_);

			  LoreVecMum.Boost(-boostB0);
			  LoreVecMup.Boost(-boostB0);
			  LoreVecK.Boost(-boostB0);
			  LoreVecPi.Boost(-boostB0);
			  TVector3 MuMuPlane = LoreVecMup.Vect().Cross(LoreVecMum.Vect());
			  TVector3 KstPlane  = LoreVecK.Vect().Cross(LoreVecPi.Vect());
			  if (MuMuPlane.Cross(KstPlane).Dot(-LoreVecB0.Vect()) > 0.0) phiKstMuMuPlane = MuMuPlane.Angle(KstPlane);
			  else                                                        phiKstMuMuPlane = -MuMuPlane.Angle(KstPlane);

			  theB0.cosThetaMup        = cosThetaMup          ;
			  theB0.cosThetaMupErr     = cosThetaMupErr       ;
			  theB0.cosThetaK          = cosThetaK            ;
			  theB0.cosThetaKErr       = cosThetaKErr         ;
			  theB0.phiKstMuMuPlane    = phiKstMuMuPlane      ;

			  }
 

			  // calculate angular vars for B0bar *****************
			  {
			  double cosThetaMum, cosThetaMumErr;
			  double cosThetaK_Bbar, cosThetaKErr_Bbar;
			  double phiKstMuMuPlane_Bbar;

			  TLorentzVector LoreVecB0;

			  TLorentzVector LoreVecMuMu;
			  LoreVecMuMu.SetXYZM(pJpsi.px(), pJpsi.py(), pJpsi.pz(), pJpsi.mass());
			  TVector3 boostMuMu = LoreVecMuMu.BoostVector();
			  TLorentzVector LoreVecMum;
			  if (mu1->charge() < 0)
				LoreVecMum.SetXYZM(mu1->px(),mu1->py(),mu1->pz(),MuMass);
			  else
				LoreVecMum.SetXYZM(mu2->px(),mu2->py(),mu2->pz(),MuMass);

			  // ################################
			  // # Boost mu- to mumu ref. frame #
			  // ################################
			  LoreVecMum.Boost(-boostMuMu);
			  // ###############################
			  // # Boost B0 to mumu ref. frame #
			  // ###############################
			  LoreVecB0.SetXYZM(pbarB.px(), pbarB.py(), pbarB.pz(),pbarB.mass());
			  LoreVecB0.Boost(-boostMuMu);

			  computeCosAlpha(-LoreVecB0.Px(),-LoreVecB0.Py(),-LoreVecB0.Pz(),
						   LoreVecMum.Px(),LoreVecMum.Py(),LoreVecMum.Pz(),
						   0.0,0.0,0.0,
						   0.0,0.0,0.0,
						   0.0,0.0,0.0,
						   0.0,0.0,0.0,
						   &cosThetaMum,&cosThetaMumErr);


			  TLorentzVector LoreVecKst;
			  LoreVecKst.SetXYZM(pKstarBar.px(), pKstarBar.py(), pKstarBar.pz(),pKstarBar.mass());

	// 		  LoreVecKst.SetXYZM(NTupleIn->kstPx->at(BestCandIndx),NTupleIn->kstPy->at(BestCandIndx),NTupleIn->kstPz->at(BestCandIndx),NTupleIn->kstBarMass->at(BestCandIndx));
			  TVector3 boostKst = LoreVecKst.BoostVector();
			  TLorentzVector LoreVecK;
			  LoreVecK.SetXYZM(p4_k.px(), p4_k.py(), p4_k.pz(), thirdTrackMass_);
			  // ##############################
			  // # Boost K- to K*0 ref. frame #
			  // ##############################
			  LoreVecK.Boost(-boostKst);
			  // ##############################
			  // # Boost B0 to K*0 ref. frame #
			  // ##############################
			  LoreVecB0.SetXYZM(pbarB.px(), pbarB.py(), pbarB.pz(),pbarB.mass());
			  LoreVecB0.Boost(-boostKst);
			  // ##################
			  // # Compute angles #
			  // ##################
			  computeCosAlpha(-LoreVecB0.Px(),-LoreVecB0.Py(),-LoreVecB0.Pz(),
						   LoreVecK.Px(),LoreVecK.Py(),LoreVecK.Pz(),
						   0.0,0.0,0.0,
						   0.0,0.0,0.0,
						   0.0,0.0,0.0,
						   0.0,0.0,0.0,
						   &cosThetaK_Bbar,&cosThetaKErr_Bbar);


			  // ######################################################################
			  // # Angle between [mu+ - mu-] and [K - pi] planes in the B0 ref. frame #
			  // ######################################################################
			  LoreVecB0.SetXYZM(pbarB.px(), pbarB.py(), pbarB.pz(),pbarB.mass());
			  TVector3 boostB0 = LoreVecB0.BoostVector();

		  
			  TLorentzVector LoreVecMup;
			  if (mu1->charge() > 0)
				LoreVecMup.SetXYZM(mu1->px(),mu1->py(),mu1->pz(),MuMass);
			  else
				LoreVecMup.SetXYZM(mu2->px(),mu2->py(),mu2->pz(),MuMass);

			  if (mu1->charge() < 0)
				LoreVecMum.SetXYZM(mu1->px(),mu1->py(),mu1->pz(),MuMass);
			  else
				LoreVecMum.SetXYZM(mu2->px(),mu2->py(),mu2->pz(),MuMass);

			  LoreVecK.SetXYZM(p4_k.px(), p4_k.py(), p4_k.pz(), thirdTrackMass_);
			  TLorentzVector LoreVecPi;
			  LoreVecPi.SetXYZM(p3_p.px(), p3_p.py(), p3_p.pz(), fourthTrackMass_);

			  LoreVecMum.Boost(-boostB0);
			  LoreVecMup.Boost(-boostB0);
			  LoreVecK.Boost(-boostB0);
			  LoreVecPi.Boost(-boostB0);
			  TVector3 MuMuPlane = LoreVecMum.Vect().Cross(LoreVecMup.Vect());
			  TVector3 KstPlane  = LoreVecK.Vect().Cross(LoreVecPi.Vect());
			  if (MuMuPlane.Cross(KstPlane).Dot(-LoreVecB0.Vect()) > 0.0) phiKstMuMuPlane_Bbar = MuMuPlane.Angle(KstPlane);
			  else                                                        phiKstMuMuPlane_Bbar = -MuMuPlane.Angle(KstPlane);
			
			  theB0.cosThetaMum          = cosThetaMum            ;
			  theB0.cosThetaMumErr       = cosThetaMumErr         ;
			  theB0.cosThetaK_Bbar       = cosThetaK_Bbar         ;
			  theB0.cosThetaKErr_Bbar    = cosThetaKErr_Bbar      ;
			  theB0.phiKstMuMuPlane_Bbar = phiKstMuMuPlane_Bbar   ;
			
			}
 
 
 
              
              theB0.Mu1Pt    = mu1 -> pt( );
              theB0.Mu1Eta   = mu1 -> eta();
              theB0.Mu1Phi   = mu1 -> phi();
              theB0.Mu1Ch    = mu1 -> charge();
              theB0.Mu2Pt    = mu2 -> pt( );
              theB0.Mu2Eta   = mu2 -> eta();
              theB0.Mu2Phi   = mu2 -> phi();
              theB0.Mu2Ch    = mu2 -> charge();

              theB0.MuMuMass = (mu1->p4() + mu2->p4()).mass();
              theB0.MuMuCL   = dimuonCL;

              theB0.JpsiL          = displacementFromBeamspotJpsi.perp();
              theB0.JpsiSigma      = sqrt(jerr.rerr(displacementFromBeamspotJpsi));
              theB0.JpsiPt         = jpperp.R();
              theB0.JpsiCosBS      = jpsi_cos;
              theB0.JpsiPosition_x = jVertex.x();
              theB0.JpsiPosition_y = jVertex.y();
              theB0.JpsiPosition_z = jVertex.z();
         
              theB0.TrkPPt     = itrk1.pt( ) ;
              theB0.TrkPEta    = itrk1.eta() ;
              theB0.TrkPPhi    = itrk1.phi() ;
              theB0.TrkPd0Sign = trk1_d0sig  ;
//            (itrk1.quality(reco::TrackBase::highPurity)) ? theB0.TrkPHQ = 1 : theB0.TrkPHQ = 0;
              
              theB0.TrkMPt     = itrk2.pt( ) ;
              theB0.TrkMEta    = itrk2.eta() ;
              theB0.TrkMPhi    = itrk2.phi() ;
              theB0.TrkMd0Sign = trk2_d0sig  ;
//            (itrk2.quality(reco::TrackBase::highPurity)) ? theB0.TrkMHQ = 1 : theB0.TrkMHQ = 0;

              theB0.KStarMass          = pKstar.mass()    ;
              theB0.barKStarMass       = pKstarBar.mass() ;
              theB0.KStarPt            = pKstar.pt()      ;
              theB0.KStarEta           = pKstar.eta()     ;
              theB0.KStarPhi           = pKstar.phi()     ;
              theB0.B0Mass             = pB.mass()        ;
              theB0.barB0Mass          = pbarB.mass()     ;
              theB0.B0StarPt           = pB.pt()          ;
              theB0.B0StarEta          = pB.eta()         ;
              theB0.B0StarPhi          = pB.phi()         ;
              
              theB0.JpsiTkTkCL         = JpsiTkTkCL;
              theB0.JpsiTkTkL          = displacementFromBeamspot.perp();
              theB0.JpsiTkTkSigma      = sqrt(err.rerr(displacementFromBeamspot));
              theB0.JpsiTkTkCosBS      = vperp.Dot(pperp)/(vperp.R()*pperp.R());
              theB0.JpsiTkTkPosition_x = secondaryVertex.x();
              theB0.JpsiTkTkPosition_y = secondaryVertex.y();
              theB0.JpsiTkTkPosition_z = secondaryVertex.z();
              
              event_.b0cands.push_back(theB0);
            }
          }
        }
      }
    }
  }
}



