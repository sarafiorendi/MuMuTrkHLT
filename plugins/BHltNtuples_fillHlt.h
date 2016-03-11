void BHltNtuples::fillHlt(const edm::Handle<edm::TriggerResults>   & triggerResults, 
                          const edm::Handle<trigger::TriggerEvent> & triggerEvent  ,
                          const edm::TriggerNames                  & triggerNames  ,
                          const edm::Event                         & event         ,
                          bool                                       isTag         )
{    
   
  for (unsigned int itrig=0; itrig < triggerNames.size(); ++itrig) 
  {
    LogDebug ("triggers") << triggerNames.triggerName(itrig) ;
    if (triggerResults->accept(itrig)) 
    {
      std::string pathName = triggerNames.triggerName(itrig);
      if (isTag) event_.hltTag.triggers.push_back(pathName);
      else       event_.hlt   .triggers.push_back(pathName);
    }
  }
     
     
  const trigger::size_type nFilters(triggerEvent->sizeFilters());
  for (trigger::size_type iFilter=0; iFilter!=nFilters; ++iFilter) 
  {
    std::string filterTag = triggerEvent->filterTag(iFilter).encode();

    trigger::Keys objectKeys = triggerEvent->filterKeys(iFilter);
    const trigger::TriggerObjectCollection& triggerObjects(triggerEvent->getObjects());
    
    for (trigger::size_type iKey=0; iKey<objectKeys.size(); ++iKey) 
    {  
      trigger::size_type objKey = objectKeys.at(iKey);
      const trigger::TriggerObject& triggerObj(triggerObjects[objKey]);
      
      HLTObjCand hltObj;
      
      hltObj.filterTag = filterTag;

      hltObj.pt  = triggerObj.pt();
      hltObj.eta = triggerObj.eta();
      hltObj.phi = triggerObj.phi();
      
      if (isTag)       event_.hltTag.objects.push_back(hltObj);
      else             event_.hlt   .objects.push_back(hltObj);
      
    }       
  }
  
}
