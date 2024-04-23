// Select events with
// a number of events greater than 0 and at least one hit in each face and maximum 3 in each panel
//

// Mu2e includes.
#include "Offline/MCDataProducts/inc/StrawDigiMC.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art_root_io/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TEnv.h"
#include "TMath.h"
#include "TSystem.h"
#include "Offline/Mu2eUtilities/inc/LsqSums4.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/VirtualDetector.hh"
#include "Offline/MCDataProducts/inc/StrawDigiMC.hh"
#include "Offline/MCDataProducts/inc/StrawGasStep.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art_root_io/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "Offline/RecoDataProducts/inc/StrawDigi.hh"
#include "Offline/Print/inc/StrawDigiPrinter.hh"
#include "Offline/Print/inc/StrawGasStepPrinter.hh"

#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"

#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
// C++ includes
#include <iostream>

using namespace std;

namespace mu2e {

  class CosmicStrawDigiMCFilter : public art::EDFilter {
    public:
      explicit CosmicStrawDigiMCFilter(fhicl::ParameterSet const& pset);
      virtual ~CosmicStrawDigiMCFilter() { }

      bool filter( art::Event& event);

    private:

    unsigned minNDigi_,nHitStations_;
    int diag_, debug_;
    art::InputTag mcDigisTag_;

  };

  CosmicStrawDigiMCFilter::CosmicStrawDigiMCFilter(fhicl::ParameterSet const& pset):
    art::EDFilter{pset},
    debug_(pset.get<int>("debug",0)){
    //mcDigisTag_(pset.get<art::InputTag>("StrawDigiMCCollection","makeSD"))
      //produces<SimParticlePtrCollection>();
  }


  bool CosmicStrawDigiMCFilter::filter(art::Event& evt) {
    bool retval(false);
    auto mcdH = evt.getValidHandle<StrawGasStepCollection>("compressDetStepMCs");
    const std::vector<mu2e::StrawGasStep>* mccol=mcdH.product(); 
    size_t mccolsize=mccol->size();
    GeomHandle<Tracker> handle;

    const Tracker*  _tracker = handle.get();
    std::vector<int> facev;
    std::vector<int> planev;

    for(size_t i=0;i<mccolsize;i++){

      StrawGasStep mcsgs = mccol->at(i);
      const Plane* plane = &_tracker->getPlane(mcsgs.strawId().getPlane());
      const Panel* panel = &plane->getPanel(mcsgs.strawId().getPanel());

      if(mcsgs.strawId().getStation()==0){
	planev.push_back(mcsgs.strawId().getPlane());
	facev.push_back(panel->getStraw(0).id().face());  
      }

    }
    int found_0=0;
    int found_1=0;
    int found_2=0;
    int found_3=0;
    int npan_00=0;
    int npan_01=0;
    int npan_02=0;
    int npan_03=0;
    int npan_04=0;
    int npan_05=0;
    int npan_10=0;
    int npan_11=0;
    int npan_12=0;
    int npan_13=0;
    int npan_14=0;
    int npan_15=0;
    for(size_t i=0;i<mccolsize;i++){
      
      StrawGasStep mcsgs = mccol->at(i);
      uint16_t panel=mcsgs.strawId().getPanel();

      if(mcsgs.strawId().getStation()==0){ 
	if(planev.at(i)==0){

	  if(facev.at(i)==0 and found_0==0){
	    found_0=1;
	 
	  }
	  if(facev.at(i)==1 and found_1==0){
	    found_1=1;
	  
	  }
	  if(panel==0) npan_00++;
	  if(panel==1) npan_01++;
	  if(panel==2) npan_02++;
	  if(panel==3) npan_03++;
	  if(panel==4) npan_04++;
	  if(panel==5) npan_05++;
	}
	if(planev.at(i)==1){

	  if(facev.at(i)==0 and found_2==0){
	    found_2=1;
	 
	  }
	  if(facev.at(i)==1 and found_3==0){
	    found_3=1;
	    
	  }
	  if(panel==0) npan_10++;
	  if(panel==1) npan_11++;
	  if(panel==2) npan_12++;
	  if(panel==3) npan_13++;
	  if(panel==4) npan_14++;
	  if(panel==5) npan_15++;
	}
      }
    }
    if (found_0 == 1 and found_1 == 1 and found_2 ==1 and found_3 ==1  ){
      if(npan_15<=3 and npan_14<=3 and npan_13<=3 and npan_12<=3 and npan_11<=3 and npan_10<=3 and npan_05<=3 and npan_04<=3 and npan_03<=3 and npan_02<=3 and npan_01<=3 and npan_00<=3){
	retval = true;
      }
    }
    return retval;
  }
}


using mu2e::CosmicStrawDigiMCFilter;
DEFINE_ART_MODULE(CosmicStrawDigiMCFilter)
