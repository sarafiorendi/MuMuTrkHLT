#include "DataFormats/Math/interface/deltaR.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"


int BHltNtuples::overlap(const reco::Candidate &a, const reco::Track b) 
{
  double eps(1.e-3);
  double dphi = deltaPhi(a.phi(), b.phi());
  dphi *= dphi;
  double deta = a.eta() - b.eta();
  deta *= deta;
  if ((dphi + deta) < eps) {
	return 1;
  }
  return 0;
}


FreeTrajectoryState BHltNtuples::initialFreeState( const reco::Track& tk, const MagneticField* field)
{
  Basic3DVector<float> pos( tk.vertex());
  GlobalPoint gpos( pos);
  Basic3DVector<float> mom( tk.momentum());
  GlobalVector gmom( mom);
  GlobalTrajectoryParameters par( gpos, gmom, tk.charge(), field);
  CurvilinearTrajectoryError err( tk.covariance());
  return FreeTrajectoryState( par, err);
}

std::pair<double,double> BHltNtuples::pionIPBeamSpot(reco::TransientTrack piTT, GlobalPoint BsGp)
{
  std::pair<double,double> measureBS;
  TrajectoryStateClosestToPoint pion_BeamSpot = piTT.trajectoryStateClosestToPoint(BsGp);
  if(pion_BeamSpot.isValid())
  {
	measureBS.first = pion_BeamSpot.perigeeParameters().transverseImpactParameter();
	if(pion_BeamSpot.hasError() && !(pion_BeamSpot.hasError()==0)) 
	{
	  measureBS.second = measureBS.first/pion_BeamSpot.perigeeError().transverseImpactParameterError();
	}
  }
  return measureBS;       
}

std::pair<double,double> pionImpactParameter(reco::TransientTrack piTT, TransientVertex jpsiVtx)
{
  std::pair<double,double> measure;
  std::pair<bool,Measurement1D>  piIP_pair = IPTools::absoluteImpactParameter3D(piTT, jpsiVtx);
  if (piIP_pair.first)
  {
	measure.first  = piIP_pair.second.value();
	measure.second = piIP_pair.second.significance();
  }
  else 
  {
	measure.first  = 0;
	measure.second = 0;
  }  
  return measure;
}

