// ======================================================================
//
// ======================================================================
#ifndef __tvst02_mod_TCalStat_hh__
#define __tvst02_mod_TCalStat_hh__

// ROOT includes
#include "TH1F.h"
//#include "TFolder.h"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "fhiclcpp/ParameterSet.h"
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/VirtualDetector.hh"
// #ifndef __CLING__ 
// #include "artdaq-core-mu2e/Overlays/FragmentType.hh"

// typedef artdaq::Fragment::type_t  type_t;

#include "artdaq-core-mu2e/Data/TrackerDataDecoder.hh"
#include "artdaq-core/Data/Fragment.hh"
//  #else 
//  namespace mu2e {
//    class TrackerFragment;
//    class TrackerFragment::TrackerDataPacket;
//  }

// namespace artdaq {
//   class Fragment;
// }
// #endif

// Mu2e includes
#include "Offline/DataProducts/inc/StrawId.hh"
#include "Offline/DataProducts/inc/TrkTypes.hh"

// #include "Stntuple/print/TAnaDump.hh"
#include "Stntuple/mod/mod/THistModule.hh"

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"

#include "Stntuple/loop/TStnModule.hh"

#include "Stntuple/obj/TStnTimeClusterBlock.hh"
#include "Stntuple/obj/TStnHelixBlock.hh"
#include "Stntuple/obj/TStnTrackSeedBlock.hh"
#include "Stntuple/obj/TStnTrackBlock.hh"
#include "Stntuple/obj/TStnClusterBlock.hh"
#include "Stntuple/obj/TCalDataBlock.hh"
#include "Stntuple/obj/TStrawHitBlock.hh"
#include "Stntuple/obj/TGenpBlock.hh"
#include "Stntuple/obj/TSimpBlock.hh"

#include "Stntuple/base/TStnArrayI.hh"

#include "Stntuple/geom/TStnCrystal.hh"
#include "Stntuple/alg/TStnTrackID.hh"
#include "Stntuple/alg/TEmuLogLH.hh"
#include "Offline/RecoDataProducts/inc/StrawDigi.hh"
#include <iostream>
#include <memory>
#include <string>

// ======================================================================
namespace mu2e {

  class TCalStat : public THistModule {

    enum { kNChannels          = 96,
           kMaxNLinks          =  12,
           kMaxNHitsPerChannel = 20,
	   kNErrorBits         = 5,
	   kNBytesErrorBit     = 0x0001,
	   kNWfsErrorBit       = 0x0002,
	   kLinkIDErrorBit     = 0x0004,
	   kChIDErrorBit       = 0x0008,
           kNChHitsErrorBit    = 0x0010
    };

  public:

  

    struct ChannelHist_t {
      TH1F*         nhits;
     

    };
                                        // per-ROC histograms
    struct EventHist_t {
      TH1F*         nhits;
     
    };

    struct RocHist_t {
    
      TH1F*         nhits;
      TH2F*         nh_vs_ch;

      ChannelHist_t channel[kNChannels];
    };

    struct ChannelData_t {
      int      nhits;
    
    };

    struct RocData_t {
     
      int       nhits;
     
      ChannelData_t  channel[kNChannels];
     
    };


//-----------------------------------------------------------------------------
// data part
//-----------------------------------------------------------------------------
                                        // in reality, this is the fragment data, 
                                        // an event can contain multiple fragments

//-----------------------------------------------------------------------------
// talk-to parameters
//-----------------------------------------------------------------------------
   
           // active links - connected ROCs
      art::InputTag    _trkfCollTag;
std::vector<int> _activeLinks;    
 int              _analyzeStation;

//-----------------------------------------------------------------------------
// the rest
//-----------------------------------------------------------------------------
    int              _nActiveLinks;
    int              _referenceChannel[kMaxNLinks][2];
     int              _initialized;   
    const art::Event*  _event;
//-----------------------------------------------------------------------------
// forgetting, for now, about multiple DTC's
//-----------------------------------------------------------------------------
    struct Hist_t {
      EventHist_t   event;
      RocHist_t     roc[kMaxNLinks];
    } _Hist;

    struct EventData_t {
            
      int       nhtot;
    
      RocData_t rdata[kMaxNLinks];

    
    } _event_data;
   

    explicit TCalStat(fhicl::ParameterSet const& pset);
    // explicit TCalStat(const art::EDAnalyzer::Table<Config>& config);
    virtual ~TCalStat() {}
    
    virtual void beginRun(const art::Run& ARun) override;

    virtual void beginJob() override;
    virtual void endJob  () override;
    
    virtual void analyze         (const art::Event& e) override;
    void         analyze_roc(const art::Event& e);

    void         book_channel_histograms(art::TFileDirectory* Dir, int RunNumber, ChannelHist_t* Hist, int Ich);
    void         book_event_histograms  (art::TFileDirectory* Dir, int RunNumber, EventHist_t*   Hist);
    void         book_roc_histograms    (art::TFileDirectory* Dir, int RunNumber, RocHist_t*     Hist, int Link);
    void         book_histograms        (int RunNumber);
  
    void         fill_channel_histograms(ChannelHist_t* Hist, ChannelData_t* Data);
    void         fill_event_histograms  (EventHist_t*   Hist, EventData_t*   Data);
    void         fill_roc_histograms    (RocHist_t*     Hist, RocData_t*     Data);

                                        // returns -1 if in trouble
    int          fill_histograms        ();

    // NWords: number of 2-byte words
   
  private : 
    const Tracker*    _tracker;     // straw tracker
  };
}
#endif
