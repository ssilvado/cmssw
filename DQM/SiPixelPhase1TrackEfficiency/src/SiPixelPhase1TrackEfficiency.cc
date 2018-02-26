// -*- C++ -*-
//
// Package:     SiPixelPhase1TrackEfficiency
// Class:       SiPixelPhase1TrackEfficiency
//

// Original Author: Marcel Schneider

#include "DQM/SiPixelPhase1Common/interface/SiPixelPhase1Base.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "TrackingTools/PatternTools/interface/TrajTrackAssociation.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h" 
#include "RecoTracker/MeasurementDet/interface/MeasurementTrackerEvent.h"
#include "TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimatorBase.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"

#include "TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimator.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"

namespace {

class SiPixelPhase1TrackEfficiency final : public SiPixelPhase1Base {
  enum {
    VALID,
    MISSING,
    INACTIVE,
    EFFICIENCY,
    VERTICES
  };

  public:
  explicit SiPixelPhase1TrackEfficiency(const edm::ParameterSet& conf);
  void analyze(const edm::Event&, const edm::EventSetup&) override;

  private:
  edm::EDGetTokenT<edmNew::DetSetVector<SiPixelCluster> > clustersToken_;
  edm::EDGetTokenT<reco::TrackCollection> tracksToken_;
  edm::EDGetTokenT<reco::VertexCollection> vtxToken_; 
  edm::EDGetTokenT<TrajTrackAssociationCollection>  trajTrackCollectionToken_;  
  edm::EDGetTokenT<MeasurementTrackerEvent> tracker_; //new
  bool applyVertexCut_;
  
  //new
 // const std::string propagatorOpposite_;
  //const std::string estimatorName_;
 // edm::ESHandle<Propagator> thePropagatorOpposite;
  ///edm::ESHandle<Chi2MeasurementEstimatorBase> theEstimator;
  const TrackerTopology*                trackerTopology_;
  const Propagator*                     trackerPropagator_;
  const MeasurementEstimator*           chi2MeasurementEstimator_;
};

SiPixelPhase1TrackEfficiency::SiPixelPhase1TrackEfficiency(const edm::ParameterSet& iConfig) :
  SiPixelPhase1Base(iConfig)//,
  //propagatorOpposite_(iConfig.getParameter<std::string>("PropagatorOpposite")),//new
 // estimatorName_(iConfig.getParameter<std::string>("Chi2MeasurementEstimator"))  //new
{
  tracker_ = consumes<MeasurementTrackerEvent>(iConfig.getParameter<edm::InputTag>("tracker")); 
  tracksToken_ = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracks"));
  vtxToken_ = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryvertices"));
  applyVertexCut_=iConfig.getUntrackedParameter<bool>("VertexCut",true);
  trajTrackCollectionToken_ = consumes<TrajTrackAssociationCollection>(iConfig.getParameter<edm::InputTag>("trajectoryInput"));//new
}

void SiPixelPhase1TrackEfficiency::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  if( !checktrigger(iEvent,iSetup,DCS) ) return;

  // get geometry
  edm::ESHandle<TrackerGeometry> tracker;
  iSetup.get<TrackerDigiGeometryRecord>().get(tracker);
  assert(tracker.isValid());

  // get primary vertex
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken( vtxToken_, vertices);

  // TrackerTopology for module informations
  edm::ESHandle<TrackerTopology> trackerTopologyHandle;
  iSetup.get<TrackerTopologyRcd>().get(trackerTopologyHandle);
  trackerTopology_ = trackerTopologyHandle.product();

  // Tracker propagator for propagating tracks to other layers
  edm::ESHandle<Propagator> propagatorHandle;
  iSetup.get<TrackingComponentsRecord>().get("PropagatorWithMaterial", propagatorHandle);
  std::unique_ptr<Propagator> propagatorUniquePtr(propagatorHandle.product() -> clone());
  trackerPropagator_ = propagatorUniquePtr.get();
  const_cast<Propagator*>(trackerPropagator_) -> setPropagationDirection(oppositeToMomentum);

  // Measurement estimator
  edm::ESHandle<Chi2MeasurementEstimatorBase> chi2MeasurementEstimatorHandle;
  iSetup.get<TrackingComponentsRecord>().get("Chi2", chi2MeasurementEstimatorHandle);
  chi2MeasurementEstimator_ = chi2MeasurementEstimatorHandle.product();

  //new
  edm::Handle<TrajTrackAssociationCollection>  trajTrackCollectionHandle;
  iEvent.getByToken(trajTrackCollectionToken_, trajTrackCollectionHandle);
  
  edm::Handle<MeasurementTrackerEvent> trackerMeas;
  iEvent.getByToken(tracker_, trackerMeas);

  //iSetup.get<TrackingComponentsRecord>().get(estimatorName_, theEstimator);
  //iSetup.get<TrackingComponentsRecord>().get(propagatorOpposite_, thePropagatorOpposite);
 

////////////////////////////////////////////////////////////////////////////////////////// 
  for(const auto &pair: *trajTrackCollectionHandle) {
    const edm::Ref<std::vector<Trajectory>> traj = pair.key;
    const reco::TrackRef tracki                   = pair.val;

    //const auto& trajectoryMeasurements = traj -> measurements();

     TrajectoryStateOnSurface tsosPXB2;
        for (const auto &tm : traj->measurements()) {
            if (tm.recHit().get() && tm.recHitR().isValid()) {
                DetId where = tm.recHitR().geographicalId();
                int  source_det_ = where.subdetId();
                int  source_layer_ = trackerTopology_ -> pxbLayer(where);
                int source_det2 = trackerTopology_->layer(where);
                /*if (source_det_ != PixelSubdetector::SubDetector::PixelBarrel ||  source_layer_ != 1) {
                    tsosPXB2 = tm.updatedState().isValid() ? tm.updatedState() : tm.backwardPredictedState();
                    break;
                }*/
                if(source_det_ == PixelSubdetector::SubDetector::PixelBarrel ) {
 			std::cout << "Pixel Barrel " << source_det_ << "  " << source_layer_ << "  " << source_det2 << std::endl; 
		}
                if(source_det_ == PixelSubdetector::SubDetector::PixelEndcap ) {
                        std::cout << "Pixel Endcap " << source_det_ << "  " << source_layer_ << "  " << source_det2 << std::endl;
                }
                if(source_det_ == StripSubdetector::SubDetector::TIB ) {
                        std::cout << "TIB " << source_det_ << "  " << source_layer_ << "  " << source_det2 << std::endl;
                }
                if(source_det_ == StripSubdetector::SubDetector::TOB ) {
                        std::cout << "TOB " << source_det_ << "  " << source_layer_ << "  " << source_det2 << std::endl;
                }
                if(source_det_ == StripSubdetector::SubDetector::TID ) {
                        std::cout << "TID " << source_det_ << "  " << source_layer_ << "  " << source_det2 << std::endl;
                }
                if(source_det_ == StripSubdetector::SubDetector::TEC ) {
                        std::cout << "TEC " << source_det_ << "  " << source_layer_ << "  " << source_det2 << std::endl;
                }
                if (!(source_det_ == PixelSubdetector::SubDetector::PixelBarrel &&  source_layer_ == 1)) {
                      tsosPXB2 = tm.updatedState().isValid() ? tm.updatedState() : tm.backwardPredictedState();
                      //break;
                      std::cout << "starting state on det " << source_det_ << "layer " << source_layer_ << "r = " << tsosPXB2.globalPosition().perp() << "z = " << tsosPXB2.globalPosition().z() << std::endl;
                }
            }
        }
        if (!tsosPXB2.isValid()) std::cout << "WARNING: did not find state for PXB2 Hit" << std::endl;
        if (!tsosPXB2.isValid()) continue; // for now
	
	const GeometricSearchTracker * gst = trackerMeas->geometricSearchTracker();
        //const auto & PXBLs = gst->pixelBarrelLayers();
        //        //for (const auto * PXBLayer : PXBLs) { std::cout << "PXB Layer with radius = " << PXBLayer->specificSurface().radius() << std::endl; }
        const auto *pxbLayer1 = gst->pixelBarrelLayers().front();
        auto compDets = pxbLayer1->compatibleDets(tsosPXB2, *trackerPropagator_, *chi2MeasurementEstimator_);
        for (const auto & detAndState : compDets) {
	             bool isHitValid   = false;
                     bool isHitMissing = false;
                     bool isHitInactive = false;
                     std::cout<<detAndState.second.localPosition().x()<<"   "<<detAndState.second.localPosition().y()<<std::endl;
                     const auto &mdet = trackerMeas->idToDet(detAndState.first->geographicalId());
                     const auto & allhits = mdet.recHits(detAndState.second);
                     std::cout << "Global position debug: x = " << detAndState.second.globalPosition().x() << " y = " <<detAndState.second.globalPosition().y() << " z = " << detAndState.second.globalPosition().z() <<  std::endl;
                     for (const auto & hit : allhits) { 
                    	float distance = std::hypot(std::abs(hit->localPosition().x()-detAndState.second.localPosition().x()), std::abs(hit->localPosition().y()-detAndState.second.localPosition().y()));		        std::cout << "Distance from hit " << distance << std::endl;
		        // Efficiency cut should eta-dependent and maybe with cluster instead of hit
		        if(std::abs(hit->localPosition().x()-detAndState.second.localPosition().x()) < 0.02 && std::abs(hit->localPosition().y()-detAndState.second.localPosition().y()) < 0.03 ){
				std::cout <<  "Hit found!!!" << std::endl;
	                        isHitValid   = hit->getType()==TrackingRecHit::valid;
			}else{
				isHitMissing = true;
			}
		     }
                     if ( (detAndState.first->geographicalId().subdetId() == PixelSubdetector::PixelBarrel && trackerTopology_ -> pxbLayer(detAndState.first->geographicalId()) == 1) ){
    			 if (isHitValid)   {
			        histo[VALID].fill(detAndState.first->geographicalId(), &iEvent);
        			histo[EFFICIENCY].fill(1, detAndState.first->geographicalId(), &iEvent);
      			}
      			if (isHitMissing) {
        			histo[MISSING].fill(detAndState.first->geographicalId(), &iEvent);
        			histo[EFFICIENCY].fill(0, detAndState.first->geographicalId(), &iEvent);
      			}
      			if (isHitInactive)   {
        			histo[INACTIVE].fill(detAndState.first->geographicalId(), &iEvent);
      			}
        	     }
         }

        

}///////////////////////////////////////////////////////////////////////////////////////////////////

  if (!vertices.isValid()) return;

  histo[VERTICES].fill(vertices->size(),DetId(0),&iEvent);

  if (applyVertexCut_ &&  vertices->empty()) return;

  // should be used for weird cuts
  //const auto primaryVertex = vertices->at(0); 

  // get the map
  edm::Handle<reco::TrackCollection> tracks;
  iEvent.getByToken( tracksToken_, tracks);
  if (!tracks.isValid()) return;

  for (auto const & track : *tracks) {
   // std::cout<<track.innerDetId()<<std::endl;

    //this cut is needed to be consisten with residuals calculation
    if (applyVertexCut_ && (track.pt() < 0.75 || std::abs( track.dxy(vertices->at(0).position()) ) > 5*track.dxyError())) continue; 

    bool isBpixtrack = false, isFpixtrack = false;
    int nStripHits = 0;

    // first, look at the full track to see whether it is good
    // auto const & trajParams = track.extra()->trajParams();
    auto hb = track.recHitsBegin();
    for(unsigned int h=0;h<track.recHitsSize();h++){
      
      auto hit = *(hb+h);
      if(!hit->isValid()) continue;

      DetId id = hit->geographicalId();
      uint32_t subdetid = (id.subdetId());

      // count strip hits
      if(subdetid==StripSubdetector::TIB) nStripHits++;
      if(subdetid==StripSubdetector::TOB) nStripHits++;
      if(subdetid==StripSubdetector::TID) nStripHits++;
      if(subdetid==StripSubdetector::TEC) nStripHits++;

      // check that we are in the pixel
      if (subdetid == PixelSubdetector::PixelBarrel) isBpixtrack = true;
      if (subdetid == PixelSubdetector::PixelEndcap) isFpixtrack = true;
    }

    if (!isBpixtrack && !isFpixtrack) continue;

    // then, look at each hit
    for(unsigned int h=0;h<track.recHitsSize();h++){
      auto hit = *(hb+h);

      DetId id = hit->geographicalId();
      uint32_t subdetid = (id.subdetId());
      if (   subdetid != PixelSubdetector::PixelBarrel 
          && subdetid != PixelSubdetector::PixelEndcap) continue;

      bool isHitValid   = hit->getType()==TrackingRecHit::valid;
      bool isHitMissing = hit->getType()==TrackingRecHit::missing;
      bool isHitInactive = hit->getType()==TrackingRecHit::inactive;

      /*
      const SiPixelRecHit* pixhit = dynamic_cast<const SiPixelRecHit*>(hit);
      const PixelGeomDetUnit* geomdetunit = dynamic_cast<const PixelGeomDetUnit*> ( tracker->idToDet(id) );
      const PixelTopology& topol = geomdetunit->specificTopology();

      // this commented part is useful if one wants ROC level maps of hits, however the local position may fall out of a ROC and the ROC maps will look very strange (with no white cross)
      LocalPoint lp;

      if (pixhit) {
        lp = pixhit->localPosition();
      } else {
        lp = trajParams[h].position();
      }

      MeasurementPoint mp = topol.measurementPosition(lp);
      int row = (int) mp.x();
      int col = (int) mp.y();
      */

      if ( !(subdetid == PixelSubdetector::PixelBarrel && trackerTopology_ -> pxbLayer(id) == 1) ){

     if (isHitValid)   {
       	histo[VALID].fill(id, &iEvent);
       	histo[EFFICIENCY].fill(1, id, &iEvent);
      }
      if (isHitMissing) {
        histo[MISSING].fill(id, &iEvent);
        histo[EFFICIENCY].fill(0, id, &iEvent);
      }
      if (isHitInactive)   {
        histo[INACTIVE].fill(id, &iEvent);
      }
	}
    }
  }
  histo[VALID   ].executePerEventHarvesting(&iEvent);
  histo[MISSING ].executePerEventHarvesting(&iEvent);
  histo[INACTIVE].executePerEventHarvesting(&iEvent);
}

} // namespace

DEFINE_FWK_MODULE(SiPixelPhase1TrackEfficiency);
