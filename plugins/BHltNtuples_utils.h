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


void computeCosAlpha (double Vx,
			     double Vy,
			     double Vz,
			     double Wx,
			     double Wy,
			     double Wz,
			     double VxErr2,
			     double VyErr2,
			     double VzErr2,
			     double VxyCov,
			     double VxzCov,
			     double VyzCov,
			     double WxErr2,
			     double WyErr2,
			     double WzErr2,
			     double WxyCov,
			     double WxzCov,
			     double WyzCov,
			     double* cosAlpha,
			     double* cosAlphaErr)
{
  double Vnorm = sqrt(Vx*Vx + Vy*Vy + Vz*Vz);
  double Wnorm = sqrt(Wx*Wx + Wy*Wy + Wz*Wz);
  double VdotW = Vx*Wx + Vy*Wy + Vz*Wz;
  
  if ((Vnorm > 0.) && (Wnorm > 0.))
    {
      *cosAlpha = VdotW / (Vnorm * Wnorm);
      *cosAlphaErr = sqrt( ((Vx*Wnorm - VdotW*Wx) * (Vx*Wnorm - VdotW*Wx) * WxErr2 +
			    (Vy*Wnorm - VdotW*Wy) * (Vy*Wnorm - VdotW*Wy) * WyErr2 +
			    (Vz*Wnorm - VdotW*Wz) * (Vz*Wnorm - VdotW*Wz) * WzErr2 +
			   
			    (Vx*Wnorm - VdotW*Wx) * (Vy*Wnorm - VdotW*Wy) * 2.*WxyCov +
			    (Vx*Wnorm - VdotW*Wx) * (Vz*Wnorm - VdotW*Wz) * 2.*WxzCov +
			    (Vy*Wnorm - VdotW*Wy) * (Vz*Wnorm - VdotW*Wz) * 2.*WyzCov) /
			   (Wnorm*Wnorm*Wnorm*Wnorm) +
			   
			   ((Wx*Vnorm - VdotW*Vx) * (Wx*Vnorm - VdotW*Vx) * VxErr2 +
			    (Wy*Vnorm - VdotW*Vy) * (Wy*Vnorm - VdotW*Vy) * VyErr2 +
			    (Wz*Vnorm - VdotW*Vz) * (Wz*Vnorm - VdotW*Vz) * VzErr2 +
			    
			    (Wx*Vnorm - VdotW*Vx) * (Wy*Vnorm - VdotW*Vy) * 2.*VxyCov +
			    (Wx*Vnorm - VdotW*Vx) * (Wz*Vnorm - VdotW*Vz) * 2.*VxzCov +
			    (Wy*Vnorm - VdotW*Vy) * (Wz*Vnorm - VdotW*Vz) * 2.*VyzCov) /
			   (Vnorm*Vnorm*Vnorm*Vnorm) ) / (Wnorm*Vnorm);
    }
  else
    {
      *cosAlpha = 0.;
      *cosAlphaErr = 0.;
    }
}
