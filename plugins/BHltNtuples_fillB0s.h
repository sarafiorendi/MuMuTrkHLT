// ---------------------------------------------------------------------
void BHltNtuples::fillB0s  (const edm::Handle<reco::MuonCollection>       & muons ,
                            const edm::Handle<reco::TrackCollection >     & tracks,
                            const edm::Handle<reco::VertexCollection >    & vtxColl,
                            const edm::Event                              & event ,
                            const edm::EventSetup                         & eventSetup)
{

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
  eventSetup.get<IdealMagneticFieldRecord>().get(mfName_, bFieldHandle);  
  const MagneticField* magField = bFieldHandle.product();
  TSCBLBuilderNoMaterial blsBuilder;

  const reco::VertexCollection pvColl  = *(vtxColl.product())  ;    

  // Loop muon collection
  for(std::vector<reco::Muon>::const_iterator mu1=muons->begin(); mu1!=muons->end(); ++mu1) 
  { 
    if( muon::isSoftMuon( (*mu1), pvColl[0] ) && (*mu1).pt() > 0 && fabs( (*mu1).eta() ) < 2.4) 
    {
      // Go on and look for the second tag muon
      for(std::vector<reco::Muon>::const_iterator mu2=mu1; mu2!=muons->end(); ++mu2) {
        if( mu2 == mu1 ) continue; 
        if( muon::isSoftMuon( (*mu2), pvColl[0]) && (*mu2).pt() > 0 && fabs( (*mu2).eta() ) < 2.4) 
        {
          if(!( mu1->charge() * mu2->charge() < 0 ))         continue; 
          hists_["mumuMass_all"]->Fill( (mu1->p4() + mu2->p4()).mass() );

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

          hists_["JpsiPt"  ]->Fill( jpperp.R() );
          hists_["JpsiLS"  ]->Fill( jpsi_ls    );
          hists_["JpsiCos" ]->Fill( jpsi_cos   );
          hists_["JpsiCL"  ]->Fill( dimuonCL   );
         
          if (jpperp.R()   < minJpsiPt_    ) continue;
          if (jpsi_ls      < minJpsiLS_    ) continue;
          if (jpsi_cos     < minJpsiCos_   ) continue;
          if (dimuonCL     < minJpsiCL_    ) continue;

          // Loop on track collection - trk 1
          for (unsigned tracksIt =0 ;  tracksIt < tracks->size(); tracksIt++)
          {
            reco::Track itrk1       = tracks->at(tracksIt) ;                                                
            
            //charge = +1 
            if (itrk1.charge() != 1)    continue;
            if (overlap(*mu1,itrk1))    continue;
            if (overlap(*mu2,itrk1))    continue;
            
            FreeTrajectoryState InitialFTS = initialFreeState(itrk1, magField);
            TrajectoryStateClosestToBeamLine tscb( blsBuilder(InitialFTS, *recoBeamSpotHandle) );
            float trkP_d0sig = tscb.transverseImpactParameter().significance();

            hists_["trkPt"] -> Fill( itrk1.pt() );
            hists_["D0sig"] -> Fill( trkP_d0sig );
 
            // eta and pt cut
            if (fabs(itrk1.eta()) > maxEta_   )                   continue;
            if (itrk1.pt()        < minPtTrk_ )                   continue;
            if (trkP_d0sig        < mind0Sign_)                   continue;
            if (! itrk1.quality(reco::TrackBase::highPurity))     continue;

            for (unsigned tracksIt2 = tracksIt + 1 ;  tracksIt2 < tracks->size(); tracksIt2++)
            {
              reco::Track itrk2       = tracks->at(tracksIt2) ;                                                
              //charge = -1 
              if (itrk2.charge() != -1)                           continue;
              if (itrk2.charge()*itrk1.charge() != -1)            continue;
              if (overlap(*mu1,itrk2))                            continue;
              if (overlap(*mu2,itrk2))                            continue;

              FreeTrajectoryState InitialFTS2 = initialFreeState(itrk2, magField);
              TrajectoryStateClosestToBeamLine tscb2( blsBuilder(InitialFTS2, *recoBeamSpotHandle) );
              float trkM_d0sig = tscb2.transverseImpactParameter().significance();
              hists_["trkPt"] -> Fill(itrk2.pt()   );
              hists_["D0sig"] -> Fill(trkM_d0sig   );
              // eta and pt cut
              if (fabs(itrk2.eta()) > maxEta_   )                 continue;
              if (itrk2.pt()        < minPtTrk_ )                 continue;
              if (trkM_d0sig        < mind0Sign_)                 continue;
              if (! itrk2.quality(reco::TrackBase::highPurity))   continue;

              reco::Particle::LorentzVector pB, pbarB, p1, p2, p3_k, p4_p, p3_p, p4_k, pKstar, pKstarBar;

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
              hists_["B0InvMass"]->Fill( pB.mass() );

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
              

              B0Cand theB0;
              
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
              theB0.TrkPd0Sign = trkP_d0sig  ;
              (itrk1.quality(reco::TrackBase::highPurity)) ? theB0.TrkPHQ = 1 : theB0.TrkPHQ = 0;
              
              theB0.TrkMPt     = itrk2.pt( ) ;
              theB0.TrkMEta    = itrk2.eta() ;
              theB0.TrkMPhi    = itrk2.phi() ;
              theB0.TrkMd0Sign = trkM_d0sig  ;
              (itrk2.quality(reco::TrackBase::highPurity)) ? theB0.TrkMHQ = 1 : theB0.TrkMHQ = 0;

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



