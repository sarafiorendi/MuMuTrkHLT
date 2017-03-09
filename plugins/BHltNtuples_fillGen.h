void BHltNtuples::fillGen(const edm::Handle<reco::GenParticleCollection> & genParticles ,
                          const edm::Event                               & event )
{

  int B0Id     =    511;
  int JpsiId   =    443;
  int kstar    =    313;
  int kplus    =    321;
  int piplus   =    211;
  int muId     =    13;

  for ( size_t i=0; i< genParticles->size(); ++i) 
  { 
    const reco::GenParticle &p = (*genParticles)[i];
    // only save interesting p
    int id = fabs( p.pdgId() );
    if (id != B0Id) continue; 
      
    GenParticleCand theGen;
    theGen.pdgId  = p.pdgId();
    theGen.pt     = p.pt() ;
    theGen.eta    = p.eta();
    theGen.phi    = p.phi();
    theGen.energy = p.energy();
    theGen.status = p.status();
    
    unsigned int n_moms = p.numberOfMothers();
    if (n_moms == 0 ){
      theGen.pdgMother.push_back(0);
      theGen.pdgRealMother.push_back(0);
    }
    else {
      for (unsigned int im=0; im < n_moms; ++im){
        theGen.pdgMother.push_back(p.motherRef(im)->pdgId());
      }
    }
    

	bool boolJpsi  = false;
    bool boolKstar = false;
    bool boolMuM   = false;
    bool boolMuP   = false;
    for ( size_t ides=0; ides < p.numberOfDaughters(); ++ides ) 
    {
  	  const reco::Candidate *des = p.daughter(ides);
	  if( des->pdgId() == JpsiId ) 
	  {
	    boolJpsi = true;
		if (des->numberOfDaughters()!=2)  		continue;  
		for ( size_t imumu=0; imumu < des->numberOfDaughters(); ++imumu ) 
		{
		  const reco::Candidate *mumu = des->daughter(imumu);
		  if ( mumu->pdgId() == -muId )       {
		    theGen.mumPt  = mumu->pt() ;
		    theGen.mumEta = mumu->eta();
		    theGen.mumPhi = mumu->phi();
		  }
		  else if ( mumu->pdgId() == muId )   {
		    theGen.mupPt  = mumu->pt() ;
		    theGen.mupEta = mumu->eta();
		    theGen.mupPhi = mumu->phi();
		  }
		  else continue;		 
		}
	  } 
	  else if( des->pdgId() == -muId ) 
	  {
	    boolMuM = true;
		theGen.mumPt  = des->pt() ;
		theGen.mumEta = des->eta();
		theGen.mumPhi = des->phi();
	  }
	  else if( des->pdgId() == muId ) 
	  {
	    boolMuP = true;
		theGen.mupPt  = des->pt() ;
		theGen.mupEta = des->eta();
		theGen.mupPhi = des->phi();
	  }
	  if( fabs(des->pdgId()) == kstar)  
	  {
		boolKstar = true;
		for ( size_t imumu=0; imumu < des->numberOfDaughters(); ++imumu ) 
		{
		  const reco::Candidate *mumu = des->daughter(imumu);
		  if      ( fabs(mumu->pdgId()) == piplus )   {
		    theGen.piPt  = mumu->pt() ;
		    theGen.piEta = mumu->eta();
		    theGen.piPhi = mumu->phi();
		  }
		  else if ( fabs(mumu->pdgId()) == kplus )   {
		    theGen.kPt  = mumu->pt() ;
		    theGen.kEta = mumu->eta();
		    theGen.kPhi = mumu->phi();
		  }
		  else continue;		 
		}
	  }	      
	} // end for des
 	if (!( (boolJpsi || (boolMuP && boolMuM) ) && boolKstar)) continue;
    event_.genParticles.push_back(theGen);

  }  // end for genParticles

}















