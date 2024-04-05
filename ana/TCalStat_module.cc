#include "TRACE/tracemf.h"
#define TRACE_NAME "TCalStat_module"
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TEnv.h"
#include "TMath.h"
#include "TSystem.h"
#include "Offline/Mu2eUtilities/inc/LsqSums4.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"

#include "Stntuple/base/TStnDataset.hh"
#include "Stntuple/loop/TStnInputModule.hh"
#include "Stntuple/loop/TStnAna.hh"
#include "Stntuple/geom/TDisk.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/VirtualDetector.hh"
#include "Stntuple/obj/TStnNode.hh"
#include "Stntuple/obj/TStnHeaderBlock.hh"
#include "Stntuple/alg/TStntuple.hh"
#include "Stntuple/geom/TCrvNumerology.hh"
#include "Stntuple/val/stntuple_val_functions.hh"
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
#include "tvst02/ana/TCalStat_module.hh"
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

namespace mu2e {

  TCalStat::TCalStat(fhicl::ParameterSet const& PSet) : 
    THistModule      (PSet                      ,"TCalStat"  ) ,
    _trkfCollTag     (PSet.get<art::InputTag>("trkfCollTag")),
   _analyzeStation(PSet.get<int>             ("analyzeStation"))
  {
    _initialized         = 0;
  }
   
 
  void TCalStat::book_channel_histograms(art::TFileDirectory* Dir, int RunNumber, ChannelHist_t* Hist, int I) {

    Hist->nhits   = Dir->make<TH1F>(Form("ch_%02i_nhits",I),Form("run %06i: ch %02i nhits"  ,RunNumber,I), 20, -0.5, 19.5);
    Hist->nhitsr   = Dir->make<TH1F>(Form("ch_%02i_nhitsr",I),Form("run %06i: ch %02i nhitsright"  ,RunNumber,I), 20, -0.5, 19.5);
    Hist->nhitsl   = Dir->make<TH1F>(Form("ch_%02i_nhitsl",I),Form("run %06i: ch %02i nhitsleft"  ,RunNumber,I), 20, -0.5, 19.5);

  }

  void TCalStat::book_panel_histograms(art::TFileDirectory* Dir, int RunNumber, PanelHist_t* Hist, int pan) {
   
    Hist->nhits           = Dir->make<TH1F>(Form("panel_%02i_nhits", pan)   , Form("run %06i, panel%02i: n hits"      ,RunNumber,pan),  300,    0.,   300.);
    Hist->nhitsr          = Dir->make<TH1F>( Form("panel_%02i_nhitsr", pan)   , Form("run %06i, panel%02i: n hitsright"      ,RunNumber,pan),  300,    0.,   300.);
    Hist->nhitsl          = Dir->make<TH1F>(Form("panel_%02i_nhitsl", pan)  , Form("run %06i, panel%02i: n hitsleft"      ,RunNumber,pan),  300,    0.,   300.);
    Hist->nh_vs_ch        = Dir->make<TH2F>(Form("panel_%02i_nh_vs_ch", pan)   , Form("run %06i, panel%02i: nh vs ch"  ,RunNumber,pan),  100,0.,100., 10,0,10);
    Hist->nh_vs_ch_r      = Dir->make<TH2F>( Form("panel_%02i_nh_vs_chr", pan)  , Form("run %06i, panel%02i: nh vs chright"  ,RunNumber,pan),  100,0.,100., 10,0,10);
    Hist->nh_vs_ch_l      = Dir->make<TH2F>( Form("panel_%02i_nh_vs_chl", pan) , Form("run %06i, panel%02i: nh vs chleft"  ,RunNumber,pan),  100,0.,100., 10,0,10);
      for (int i=0; i<kNChannels; i++) {
	art::TFileDirectory chan_dir = Dir->mkdir(Form("ch_%02i",i));
	book_channel_histograms(&chan_dir,RunNumber,&Hist->channel[i],i);
      }
   
    
  }

  void TCalStat::book_plane_histograms(art::TFileDirectory* Dir, int RunNumber, PlaneHist_t* Hist, int pl) {

    Hist->nhits   = Dir->make<TH1F>(Form("plane_%02i_nhits",pl),Form("run %06i: plane %02i nhits"  ,RunNumber,pl), 20, -0.5, 19.5);
     for (int i=0; i<kNPanels; i++) {
	art::TFileDirectory panel_dir = Dir->mkdir(Form("panel_%02i",i));
	book_panel_histograms(&panel_dir,RunNumber,&Hist->panel[i],i);
      }

  }
  void TCalStat::book_event_histograms(art::TFileDirectory* Dir, int RunNumber, EventHist_t* Hist) {
   
    Hist->nhits           = Dir->make<TH1F>("nhits"      , Form("run %06i: nhits total"  ,RunNumber), 1000, 0.,   1000.);
   
  }

  void TCalStat::book_histograms(int RunNumber) {
    art::ServiceHandle<art::TFileService> tfs;
    
    art::TFileDirectory top_dir = tfs->mkdir("trk");

    book_event_histograms(&top_dir,RunNumber,&_Hist.event);

    for (int i=0; i<kNPlanes; i++) {
      art::TFileDirectory plane_dir = top_dir.mkdir(Form("plane_%i",i));
      book_plane_histograms(&plane_dir,RunNumber,&_Hist.plane[i],i);
    }
    
    printf("[mu2e::TCalStat] pointer to the module: 0x%8p\n",(void*) this);
  }


void TCalStat::beginRun(const art::Run& aRun) {
  int rn  = aRun.run();
   if (_initialized != 0) return;
   _initialized = 1;
   book_histograms(rn);
  
}


  void TCalStat::beginJob() {
  }


  void TCalStat::endJob() {
    printf("[mu2e::TCalStat] pointer to the module: 0x%8p\n",(void*) this);
  }


  void TCalStat::fill_channel_histograms(ChannelHist_t* Hist, ChannelData_t* Data) {
    if(Data->nhits>0){
      Hist->nhits->Fill(Data->nhits);
    }
    if(Data->nhitsr>0){
      Hist->nhitsr->Fill   (Data->nhitsr);
    }
    if(Data->nhitsl>0){
      Hist->nhitsl->Fill   (Data->nhitsl);
    }
    
  }


  void TCalStat::fill_panel_histograms(PanelHist_t* Hist, PanelData_t* Rd) {
    if(Rd->nhits>0){
      Hist->nhits->Fill   (Rd->nhits);
    }
    if(Rd->nhitsr>0){
      Hist->nhitsr->Fill   (Rd->nhitsr);
    }
    if(Rd->nhitsl>0){
      Hist->nhitsl->Fill   (Rd->nhitsl);
    }

    for (int ich=0; ich<kNChannels; ich++) {

      ChannelData_t* chd = &Rd->channel[ich];
      ChannelHist_t* hch = &Hist->channel[ich];

      Hist->nh_vs_ch->Fill(ich,chd->nhits);
      Hist->nh_vs_ch_r->Fill(ich,chd->nhitsr);
      Hist->nh_vs_ch_l->Fill(ich,chd->nhitsl);

      fill_channel_histograms(hch,chd);
     
    }
  }
      
  void TCalStat::fill_plane_histograms(PlaneHist_t* Hist, PlaneData_t* Rd) {
    if(Rd->nhits>0){
      Hist->nhits->Fill   (Rd->nhits);
    }
    for (int i=0; i<kNPanels; i++) {
      PanelData_t* chd = &Rd->panel[i];
      PanelHist_t* hch = &Hist->panel[i];
      fill_panel_histograms(hch,chd);
     
    }
  }

  void TCalStat::fill_event_histograms(EventHist_t* Hist, EventData_t* Data) {

  
    Hist->nhits->Fill(Data->nhtot);

  }


int TCalStat::fill_histograms() {

    fill_event_histograms(&_Hist.event,&_event_data);

    for (int i=0; i<kNPlanes; i++) {
      fill_plane_histograms(&_Hist.plane[i],&_event_data.plane[i]);
    }

    return 0;
    }


  void TCalStat::analyze_panel(const art::Event& Evt) {

    return;
  }


void TCalStat::analyze(const art::Event& event) {
  //art::EventNumber_t eventNumber = event.event();

  _event = &event;


  _event_data.nhtot = 0;
  for(int i=0;i<kNPlanes;i++){
    _event_data.plane[i].nhits=0;
    for(int panel=0;panel<kNPanels;panel++){
      _event_data.plane[i].panel[panel].nhits=0;
      _event_data.plane[i].panel[panel].nhitsr=0;
      _event_data.plane[i].panel[panel].nhitsl=0;
    for(int straw=0;straw<kNChannels;straw++){
	_event_data.plane[i].panel[panel].channel[straw].nhits=0;
	_event_data.plane[i].panel[panel].channel[straw].nhitsr=0;
	_event_data.plane[i].panel[panel].channel[straw].nhitsl=0;
      }
    }
  }


  mu2e::GeomHandle<mu2e::Tracker> th;
  _tracker = th.get();
  
  auto handle = event.getValidHandle<StrawGasStepCollection>(_trkfCollTag);
  const std::vector<mu2e::StrawGasStep>* mccol=handle.product(); 

  size_t mccolsize=mccol->size();
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

  if (found_0 == 0 || found_1 == 0 ||found_2 == 0 || found_3 == 0  ){
    return;
  }

  _event_data.nhtot=_event_data.nhtot+mccol->size();


  for (size_t i=0;i<mccolsize;i++) {

    StrawGasStep mcsgs = mccol->at(i);
    uint16_t panel = mcsgs.strawId().getPanel();
    uint16_t plane = mcsgs.strawId().getPlane();
    uint16_t straw = mcsgs.strawId().getStraw();
    double x    = mcsgs.position().x();
    double y    = mcsgs.position().y();
    const Plane* pln = &_tracker->getPlane(mcsgs.strawId().getPlane());
    const Panel* pnl = &pln->getPanel(panel);
    const Straw* strw = &pnl->getStraw(straw);
    double costheta = strw->direction().x();
    double sintheta = strw->direction().y();
    double xp = x*costheta+y*sintheta;

    if(mcsgs.strawId().getStation()==0){ 
      if(xp>strw->getMidPoint().x()){
	_event_data.plane[plane].panel[panel].nhitsr +=1;
	_event_data.plane[plane].panel[panel].channel[straw].nhitsr += 1;
      }else{
	_event_data.plane[plane].panel[panel].nhitsl +=1;
	_event_data.plane[plane].panel[panel].channel[straw].nhitsl += 1;
      }
	_event_data.plane[plane].panel[panel].nhits +=1;
	_event_data.plane[plane].panel[panel].channel[straw].nhits += 1;
	_event_data.plane[plane].nhits +=1;
    }
       std::cout<<"plane nhits right"<<_event_data.plane[plane].panel[panel].nhitsr<<std::endl;
 
  }
  if ( _analyzeStation == 1 and mccolsize>0)  fill_histograms();

  TModule::analyze(event);
}

} // end namespace mu2e

DEFINE_ART_MODULE(mu2e::TCalStat)
