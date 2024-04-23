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
	   kNPanels            = 6,
	   kNPlanes            = 2
    };

  public:

    struct ChannelHist_t {
      TH1F*         nhits;
    };
                                   
    struct PanelHist_t {
    
      TH1F*         nhits;
      TH1F*         nhitsr;
      TH1F*         nhitsl;
      TH2F*         nh_vs_ch;
      TH2F*         xp_vs_yp;
      TH2F*         x_vs_y;
      TH2F*         xp_vs_yp_nh3;
      TH1F*         x_bias;
      TH1F*         y_bias;
      TH1F*         x;
      TH1F*         y;
      TH2F*         vdrift;
      TH2F*         vdrift_true;
      TH1F*         x_bias_stereo_combo;
      TH1F*         y_bias_stereo_combo;
      TH2F*         x_bias_vs_x;
      TH2F*         x_bias_vs_x_true;

      ChannelHist_t channel[kNChannels];
    };
    

    struct PlaneHist_t {
      TH1F*         nhits;
      TH2F*         x_vs_y;
      PanelHist_t   panel[kNPanels];
     
    };

    struct EventHist_t {
      TH1F*         nhits;
      TH1F*         nhits_nh3;
      TH1F*         nhits_nh3_ex;
      TH2F*         x_vs_y;
      TH2F*         x_vs_y_nh3;
      TH2F*         x_vs_y_nh3_ex;
      TH2F*         x_vs_y_rec_p;
      TH2F*         x_vs_y_rec_pf;
      TH1F*         t0;//distribuzione temporale di ogni hit
      TH1F*         costheta;//distribuzione temporale di ogni hit
      TH2F*         x_vs_y_rec_true;

      PlaneHist_t   plane[kNPlanes];
    };

    struct ChannelData_t {
      int       nhits;

    
    };

    struct PanelData_t {
     
      int            nhits;

     
      ChannelData_t  channel[kNChannels];
     
    };

    struct PlaneData_t {
      int          nhits;
      PanelData_t  panel[kNPanels];
     
    };

    struct Hist_t {
      EventHist_t   event;
      PlaneHist_t   plane[kNPlanes];
    } _Hist;

    struct EventData_t {      
      int         nhits;
      int         nhits_nh3;
      int         nhits_nh3_ex;
      PlaneData_t plane[kNPlanes];    
    } _event_data;

//-----------------------------------------------------------------------------
// talk-to parameters
//-----------------------------------------------------------------------------
   
          
    art::InputTag     _trkfCollTag;
    int               _analyzeStation;
    art::InputTag     _digiCollTag;
    art::InputTag     _recCollTag;
    art::InputTag     _comboCollTag;

    int                _initialized;   

    const art::Event*  _event;


   

    explicit TCalStat(fhicl::ParameterSet const& pset);
    // explicit TCalStat(const art::EDAnalyzer::Table<Config>& config);
    virtual ~TCalStat() {}
    
    virtual void beginRun(const art::Run& ARun) override;

    virtual void beginJob() override;
    virtual void endJob  () override;
    
    virtual void analyze         (const art::Event& e) override;
    void         analyze_panel(const art::Event& e);
    //void         fit_line(double * m, double * sigma_m, double * q, double * sigma_q, double x[], double y[], double dx[], double dy[], int size, double m1);

    void         book_channel_histograms(art::TFileDirectory* Dir, int RunNumber, ChannelHist_t* Hist, int Ich);
    void         book_event_histograms  (art::TFileDirectory* Dir, int RunNumber, EventHist_t*   Hist);
    void         book_panel_histograms    (art::TFileDirectory* Dir, int RunNumber, PanelHist_t*     Hist, int pan);
    void         book_plane_histograms    (art::TFileDirectory* Dir, int RunNumber, PlaneHist_t*     Hist, int pl);

    void         book_histograms        (int RunNumber);
  
    void         fill_channel_histograms(ChannelHist_t* Hist, ChannelData_t* Data);
    void         fill_event_histograms  (EventHist_t*   Hist, EventData_t*   Data);
    void         fill_panel_histograms    (PanelHist_t*     Hist, PanelData_t*     Data);
    void         fill_plane_histograms    (PlaneHist_t*     Hist, PlaneData_t*     Data);
                                        // returns -1 if in trouble
    int          fill_histograms        ();

    // NWords: number of 2-byte words
   
  private : 
    const Tracker*    _tracker;     // straw tracker
  };
}
#endif
