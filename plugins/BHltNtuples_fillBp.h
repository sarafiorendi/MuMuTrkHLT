// ---------------------------------------------------------------------
void BHltNtuples::fillBp  (const edm::Handle<reco::MuonCollection>       & muons ,
                            const edm::Handle<reco::TrackCollection >     & tracks,
                            const edm::Handle<reco::VertexCollection >    & vtxColl,
                            const edm::Event                              & event ,
                            const edm::EventSetup                         & eventSetup)
{

  const double thirdTrackMass2  = thirdTrackMass_ *thirdTrackMass_ ; // kaon

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
  selTracksDef qualityTracks;
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

    qualityTracks[tracksIt] = itrk;     
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
          for (selTracksDef::const_iterator tracksIt=qualityTracks.begin(); tracksIt!=qualityTracks.end(); ++tracksIt)
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
 
              reco::Particle::LorentzVector pB, p1, p2, p3_k, pJpsi;

              // Combined system
              double e1   = sqrt(mu1->momentum().Mag2()  + MuMass2          );
              double e2   = sqrt(mu2->momentum().Mag2()  + MuMass2          );
              double e3_k = sqrt(itrk1.momentum().Mag2() + thirdTrackMass2  );
            
              p1   = reco::Particle::LorentzVector(mu1->px() , mu1->py() , mu1->pz() , e1  );
              p2   = reco::Particle::LorentzVector(mu2->px() , mu2->py() , mu2->pz() , e2  );
              p3_k = reco::Particle::LorentzVector(itrk1.px(), itrk1.py(), itrk1.pz(), e3_k);
            
              pB        = p1 + p2 + p3_k;
              pJpsi = p1 + p2;
            
              if (pB.mass() > 7. || pB.mass() < 4.) continue;

              // do the vertex fit
              std::vector<reco::TransientTrack> t_tks;
              t_tks.push_back((*theB).build(mu1->track().get()));
              t_tks.push_back((*theB).build(mu2->track().get()));
              t_tks.push_back((*theB).build(&itrk1));
              if (t_tks.size()!=3) continue;
            
              KalmanVertexFitter kvf;
              TransientVertex tv  = kvf.vertex(t_tks);
              reco::Vertex vertex = tv;
              if (!tv.isValid()) continue;
//               hists_["B0InvMass"]->Fill( pB.mass() );

              float JpsiTkTkCL = 0;
              if ((vertex.chi2()>=0.0) && (vertex.ndof()>0) )   
                JpsiTkTkCL = TMath::Prob(vertex.chi2(), vertex.ndof() );
              
              // calculate four-track transverse momentum
              math::XYZVector pperp(mu1->px() + mu2->px() + itrk1.px(),
                                    mu1->py() + mu2->py() + itrk1.py(),
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
 
 
              BpCand theBp;

              theBp.Mu1Pt    = mu1 -> pt( );
              theBp.Mu1Eta   = mu1 -> eta();
              theBp.Mu1Phi   = mu1 -> phi();
              theBp.Mu1Ch    = mu1 -> charge();
              theBp.Mu2Pt    = mu2 -> pt( );
              theBp.Mu2Eta   = mu2 -> eta();
              theBp.Mu2Phi   = mu2 -> phi();
              theBp.Mu2Ch    = mu2 -> charge();

              theBp.MuMuMass = (mu1->p4() + mu2->p4()).mass();
              theBp.MuMuCL   = dimuonCL;

              theBp.JpsiL          = displacementFromBeamspotJpsi.perp();
              theBp.JpsiSigma      = sqrt(jerr.rerr(displacementFromBeamspotJpsi));
              theBp.JpsiPt         = jpperp.R();
              theBp.JpsiCosBS      = jpsi_cos;
              theBp.JpsiPosition_x = jVertex.x();
              theBp.JpsiPosition_y = jVertex.y();
              theBp.JpsiPosition_z = jVertex.z();
         
              theBp.TrkPt     = itrk1.pt( ) ;
              theBp.TrkEta    = itrk1.eta() ;
              theBp.TrkPhi    = itrk1.phi() ;
              theBp.Trkd0Sign = trk1_d0sig  ;
//            (itrk1.quality(reco::TrackBase::highPurity)) ? theBp.TrkPHQ = 1 : theBp.TrkPHQ = 0;
              

              theBp.BMass        = pB.mass()        ;
              theBp.BPt           = pB.pt()          ;
              theBp.BEta          = pB.eta()         ;
              theBp.BPhi          = pB.phi()         ;
              
              theBp.JpsiTkCL         = JpsiTkTkCL;
              theBp.JpsiTkL          = displacementFromBeamspot.perp();
              theBp.JpsiTkSigma      = sqrt(err.rerr(displacementFromBeamspot));
              theBp.JpsiTkCosBS      = vperp.Dot(pperp)/(vperp.R()*pperp.R());
              theBp.JpsiTkPosition_x = secondaryVertex.x();
              theBp.JpsiTkPosition_y = secondaryVertex.y();
              theBp.JpsiTkPosition_z = secondaryVertex.z();
              
              event_.bpcands.push_back(theBp);
          }
        }
      }
    }
  }
}



