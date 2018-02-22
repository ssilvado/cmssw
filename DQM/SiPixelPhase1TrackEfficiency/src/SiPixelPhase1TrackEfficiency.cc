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
#include "TrackingTools/TrackRefitter/interface/TrackTransformer.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "RecoTracker/MeasurementDet/interface/MeasurementTrackerEvent.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h" 

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
  TrackTransformer refitter_;
  bool applyVertexCut_;
  
};

SiPixelPhase1TrackEfficiency::SiPixelPhase1TrackEfficiency(const edm::ParameterSet& iConfig) :
 
 SiPixelPhase1Base(iConfig),
 refitter_(iConfig) 
{ 
  tracksToken_ = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracks"));
  vtxToken_ = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryvertices"));
  applyVertexCut_=iConfig.getUntrackedParameter<bool>("VertexCut",true);
   
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
 
  refitter_.setServices(iSetup);

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
    std::vector<Trajectory> traj  = refitter_.transform(track);
    if (traj.size() != 1) continue; 

    TrajectoryStateOnSurface tsosPXB2;
    for (const auto &tm : traj.front().measurements()) {
            if (tm.recHit().get() && tm.recHitR().isValid()) {
                tsosPXB2 = tm.updatedState().isValid() ? tm.updatedState() : tm.backwardPredictedState();
							//DetId where = tm.recHitR().geographicalId();
                //source_det_ = where.subdetId();
                //source_layer_ = theTrkTopo->layer(where);
                //if (source_det_ != PixelSubdetector::SubDetector::PixelBarrel ||  source_layer_ != 1) {
                //    tsosPXB2 = tm.updatedState().isValid() ? tm.updatedState() : tm.backwardPredictedState();
                //    if (debug_) printf("starting state on det %d, layer %d, r = %5.1f, z = %+6.1f\n", 
                 //                       source_det_, source_layer_, tsosPXB2.globalPosition().perp(), tsosPXB2.globalPosition().z());
                  //  break;
                //}
            }
        }

//    const GeometricSearchTracker * gst = tracker->geometricSearchTracker();
        //const auto & PXBLs = gst->pixelBarrelLayers();
        //        //for (const auto * PXBLayer : PXBLs) { std::cout << "PXB Layer with radius = " << PXBLayer->specificSurface().radius() << std::endl; }
  //  const auto *pxbLayer1 = gst->pixelBarrelLayers().front();
   // auto compDets = pxbLayer1->compatibleDets(tsosPXB2, *thePropagatorOpposite, *theEstimator);
    
   // for (const auto & detAndState : compDets) {

//	const auto & tkpos = detAndState.second.localPosition();
    //    const auto & tkerr = detAndState.second.localError().positionError();
  //    }

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
  histo[VALID   ].executePerEventHarvesting(&iEvent);
  histo[MISSING ].executePerEventHarvesting(&iEvent);
  histo[INACTIVE].executePerEventHarvesting(&iEvent);
}

} // namespace

DEFINE_FWK_MODULE(SiPixelPhase1TrackEfficiency);

