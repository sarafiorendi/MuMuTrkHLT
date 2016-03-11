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
    if (id != B0Id  &&  id != JpsiId  &&  id != kstar  &&  id != kplus  &&  id != piplus  &&  id != muId ) continue; 
      
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
        // if coming from a muon, go back one step ** to be improved **
//         if(n_moms == 1 && fabs(p.motherRef(0)->pdgId()) == muId){
//           for (unsigned int igm = 0; igm < p.motherRef(0)->numberOfMothers(); igm++){
//             theGen.pdgRealMother.push_back(p.motherRef(0)->motherRef(igm)->pdgId());
//           }
//         }
//         else
//           theGen.pdgRealMother.push_back(0);
      }
    }

    event_.genParticles.push_back(theGen);

  }  // end for genParticles

}


