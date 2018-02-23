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

  bool applyVertexCut_;
  
};

SiPixelPhase1TrackEfficiency::SiPixelPhase1TrackEfficiency(const edm::ParameterSet& iConfig) :
  SiPixelPhase1Base(iConfig) 
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

    //check all the hits
    std::vector< bool > SubDetCheckerBool = {false, false, false, false, false, false, false};

    for(unsigned int h=0;h<track.recHitsSize();h++){ ///check the location of all the hits
      auto hit = *(hb+h);
      DetId id = hit->geographicalId();

      bool isHitValid   = hit->getType()==TrackingRecHit::valid;
      if (isHitValid){
      if (id() >= 303042564 && id() <= 303087648) SubDetCheckerBool.at(0) = true;  //Layer1               
      if (id() >= 304091140 && id() <= 304201760) SubDetCheckerBool.at(1) = true;  //Layer2                                                                       
      if (id() >= 305139716 && id() <= 305315872) SubDetCheckerBool.at(2) = true;  //Layer3                                                                      
      if (id() >= 306188292 && id() <= 306446368) SubDetCheckerBool.at(3) = true;  //Layer4                                                     
      if ((id() >= 344200196 && id() <= 344426500) || (id() >= 352588804 && id() <= 352815108)) SubDetCheckerBool.at(4) = true;  //Disk1       
      if ((id() >= 344462340 && id() <= 344688644) || (id() >= 352588804 && id() <= 352815108)) SubDetCheckerBool.at(5) = true;  //Disk2                     
      if ((id() >= 344724484 && id() <= 344950788) || (id() >= 353113092 && id() <= 353339396)) SubDetCheckerBool.at(6) = true;  //Disk3                          
      }
    }


    // then, look at each hit
    for(unsigned int h=0;h<track.recHitsSize();h++){
      auto hit = *(hb+h);

      DetId id = hit->geographicalId();
      uint32_t subdetid = (id.subdetId());

      //Check the location of this hit
      std::vector< bool > hitCheckerBool = {false, false, false, false, false, false, false};

      if (id() >= 303042564 && id() <= 303087648) hitCheckerBool.at(0) = true;
      if (id() >= 304091140 && id() <= 304201760) hitCheckerBool.at(1) = true;
      if (id() >= 305139716 && id() <= 305315872) hitCheckerBool.at(2) = true;
      if (id() >= 306188292 && id() <= 306446368) hitCheckerBool.at(3) = true;
      if ((id() >= 344200196 && id() <= 344426500) || (id() >= 352588804 && id() <= 352815108)) hitCheckerBool.at(4) = true;  //Disk1                       
      if ((id() >= 344462340 && id() <= 344688644) || (id() >= 352588804 && id() <= 352815108)) hitCheckerBool.at(5) = true;  //Disk2                        
      if ((id() >= 344724484 && id() <= 344950788) || (id() >= 353113092 && id() <= 353339396)) hitCheckerBool.at(6) = true;  //Disk3                          



      if (   subdetid != PixelSubdetector::PixelBarrel 
          && subdetid != PixelSubdetector::PixelEndcap) continue;

      bool isHitValid   = hit->getType()==TrackingRecHit::valid;
      bool isHitMissing = hit->getType()==TrackingRecHit::missing;
      bool isHitInactive = hit->getType()==TrackingRecHit::inactive;

      // Hp cut
      int TRACK_QUALITY_HIGH_PURITY_BIT = 2;
      int TRACK_QUALITY_HIGH_PURITY_MASK = 1 << TRACK_QUALITY_HIGH_PURITY_BIT;
      if(!((track.qualityMask() & TRACK_QUALITY_HIGH_PURITY_MASK) >> TRACK_QUALITY_HIGH_PURITY_BIT)) 
	{
	  isHitValid = false;

	}

      //Pt cut
      float TRACK_PT_CUT_VAL = 1.0f;
      if(!(TRACK_PT_CUT_VAL < track.pt()))
	{
	  isHitValid = false;

	}

      //Nstrip cut
      int TRACK_NSTRIP_CUT_VAL = 10;
      if(!(TRACK_NSTRIP_CUT_VAL < nStripHits))
	{
	  isHitValid = false;

	}

      //Apply pixel hit cuts
      if (hitCheckerBool.at(0) == true) if (!( //layer1
      					      (SubDetCheckerBool.at(1) == true && SubDetCheckerBool.at(2) == true && SubDetCheckerBool.at(3) == true) ||
      					      (SubDetCheckerBool.at(1) == true && SubDetCheckerBool.at(2) == true && SubDetCheckerBool.at(4) == true) || 
      					      (SubDetCheckerBool.at(1) == true && SubDetCheckerBool.at(4) == true && SubDetCheckerBool.at(5) == true) ||
                                              (SubDetCheckerBool.at(4) == true && SubDetCheckerBool.at(5) == true && SubDetCheckerBool.at(6) == true)))
      					  {
					    isHitValid = false;
      					  }
      if (hitCheckerBool.at(1) == true) if (!( //layer2                                                                                   
                                              (SubDetCheckerBool.at(0) == true && SubDetCheckerBool.at(2) == true && SubDetCheckerBool.at(3) == true) ||
                                              (SubDetCheckerBool.at(0) == true && SubDetCheckerBool.at(2) == true && SubDetCheckerBool.at(4) == true) ||
                                              (SubDetCheckerBool.at(0) == true && SubDetCheckerBool.at(4) == true && SubDetCheckerBool.at(5) == true)))
      					  {
      					    isHitValid = false;
      					  }
      if (hitCheckerBool.at(3) == true) if (!( //layer3                                                                                                            
                                              (SubDetCheckerBool.at(0) == true && SubDetCheckerBool.at(1) == true && SubDetCheckerBool.at(3) == true) ||
                                              (SubDetCheckerBool.at(0) == true && SubDetCheckerBool.at(1) == true && SubDetCheckerBool.at(4) == true)))
                                          {
					    isHitValid = false;
                                          }

      if (hitCheckerBool.at(3) == true) if (!( //layer4                                                                                                           
                                              (SubDetCheckerBool.at(0) == true && SubDetCheckerBool.at(1) == true && SubDetCheckerBool.at(2) == true)))
					  {
					    isHitValid = false;
                                          }

      if (hitCheckerBool.at(4) == true) if (!( //Disk1                                                                                                            
                                              (SubDetCheckerBool.at(0) == true && SubDetCheckerBool.at(1) == true && SubDetCheckerBool.at(2) == true) ||
                                              (SubDetCheckerBool.at(0) == true && SubDetCheckerBool.at(1) == true && SubDetCheckerBool.at(5) == true) ||
                                              (SubDetCheckerBool.at(0) == true && SubDetCheckerBool.at(5) == true && SubDetCheckerBool.at(6) == true)))
                                          {
                                            isHitValid = false;
                                          }
      if (hitCheckerBool.at(5) == true) if (!( //Disk2                                                                                                            
                                              (SubDetCheckerBool.at(0) == true && SubDetCheckerBool.at(1) == true && SubDetCheckerBool.at(1) == true) ||
                                              (SubDetCheckerBool.at(0) == true && SubDetCheckerBool.at(4) == true && SubDetCheckerBool.at(6) == true)))
                                          {
					    isHitValid = false;
                                          }
      if (hitCheckerBool.at(6) == true) if (!( //Disk3                                                                                                    
      					      (SubDetCheckerBool.at(0) == true && SubDetCheckerBool.at(4) == true && SubDetCheckerBool.at(5) == true)))
                                          {
					    isHitValid = false;
                                          }


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

