// Select events with
// a number of events bigger than 0
//

// Mu2e includes.
#include "Offline/MCDataProducts/inc/StrawDigiMC.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
// Framework includes.
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art_root_io/TFileService.h"
// Other includes
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

//#include "Stntuple/base/TStnDataset.hh"
//#include "Stntuple/loop/TStnInputModule.hh"
//#include "Stntuple/loop/TStnAna.hh"
//#include "Stntuple/geom/TDisk.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/VirtualDetector.hh"
//#include "Stntuple/obj/TStnNode.hh"
//#include "Stntuple/obj/TStnHeaderBlock.hh"
//#include "Stntuple/alg/TStntuple.hh"
//#include "Stntuple/geom/TCrvNumerology.hh"
//#include "Stntuple/val/stntuple_val_functions.hh"
// Mu2e includes.
#include "Offline/MCDataProducts/inc/StrawDigiMC.hh"
#include "Offline/MCDataProducts/inc/StrawGasStep.hh"

#include "Offline/MCDataProducts/inc/SimParticle.hh"
// Framework includes.
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art_root_io/TFileService.h"
// Other includes
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
    for(size_t i=0;i<mccolsize;i++){
      
      StrawGasStep mcsgs = mccol->at(i);
    
      if(mcsgs.strawId().getStation()==0){ 
	if(planev.at(i)==0){

	  if(facev.at(i)==0 and found_0==0){
	    found_0=1;
	 
	  }
	  if(facev.at(i)==1 and found_1==0){
	    found_1=1;
	  
	  }
	}
	if(planev.at(i)==1){

	  if(facev.at(i)==0 and found_2==0){
	    found_2=1;
	 
	  }
	  if(facev.at(i)==1 and found_3==0){
	    found_3=1;
	    
	  }
	}
      }
    }
    if (found_0 == 1 and found_1 == 1 and found_2 ==1 and found_3 ==1  ){
      retval = true;
    }
    return retval;
  }
}

using mu2e::CosmicStrawDigiMCFilter;
DEFINE_ART_MODULE(CosmicStrawDigiMCFilter)
