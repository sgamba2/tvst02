///////////////////////////////////////////////////////////////////////////////
// event error codes: 
// 0: OK
// 1: ich>95
// 2: nh >= kMaxNHitsPerChannel
// 3: fragment too long
// 4: wrong active link ID
// 5: nbytes < 2
// 6: ich < 0
// 7: hit->NumADCPackets > 2
///////////////////////////////////////////////////////////////////////////////
#include "TRACE/tracemf.h"
#define TRACE_NAME "TCalStat_module"
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TEnv.h"
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
namespace mu2e {
//-----------------------------------------------------------------------------
// ======================================================================

  TCalStat::TCalStat(fhicl::ParameterSet const& PSet) : 
    THistModule      (PSet                      ,"TCalStat"  ) ,
    _trkfCollTag     (PSet.get<art::InputTag>("trkfCollTag")),
    _activeLinks     (PSet.get<std::vector<int>>("activeLinks"     )),
   _analyzeStation(PSet.get<int>             ("analyzeStation"))
  {

  
    _nActiveLinks        = _activeLinks.size();
    _initialized         = 0;
  }
   
  


//-----------------------------------------------------------------------------
// I : channel number
//-----------------------------------------------------------------------------
  void TCalStat::book_channel_histograms(art::TFileDirectory* Dir, int RunNumber, ChannelHist_t* Hist, int I) {

    Hist->nhits   = Dir->make<TH1F>(Form("ch_%02i_nhits",I),Form("run %06i: ch %02i nhits"  ,RunNumber,I), 20, -0.5, 19.5);

  }

//-----------------------------------------------------------------------------
  void TCalStat::book_roc_histograms(art::TFileDirectory* Dir, int RunNumber, RocHist_t* Hist, int Link) {
   
    Hist->nhits           = Dir->make<TH1F>("nhits"   , Form("run %06i: n hits"      ,RunNumber),  300,    0.,   300.);
    
   
      for (int i=0; i<kNChannels; i++) {
	art::TFileDirectory chan_dir = Dir->mkdir(Form("ch_%02i",i));
	book_channel_histograms(&chan_dir,RunNumber,&Hist->channel[i],i);
      }
   
    
  }

//-----------------------------------------------------------------------------
  void TCalStat::book_event_histograms(art::TFileDirectory* Dir, int RunNumber, EventHist_t* Hist) {
   
    Hist->nhits           = Dir->make<TH1F>("nhits"      , Form("run %06i: nhits total"  ,RunNumber), 1000, 0.,   1000.);
   
  }

//-----------------------------------------------------------------------------
  void TCalStat::book_histograms(int RunNumber) {
    art::ServiceHandle<art::TFileService> tfs;
    
    art::TFileDirectory top_dir = tfs->mkdir("trk");

    book_event_histograms(&top_dir,RunNumber,&_Hist.event);

    for (int i=0; i<_nActiveLinks; i++) {
      int link  = _activeLinks[i];
      art::TFileDirectory roc_dir = top_dir.mkdir(Form("roc_%i",link));
      book_roc_histograms(&roc_dir,RunNumber,&_Hist.roc[link],link);
    }
    
    printf("[mu2e::TCalStat] pointer to the module: 0x%8p\n",(void*) this);
  }

//-----------------------------------------------------------------------------
// for run<=285, wasn't saving the event header and subheader, 64 bytes in total
// or 32 2-byte words
// book histograms only once
//-----------------------------------------------------------------------------
void TCalStat::beginRun(const art::Run& aRun) {
  int rn  = aRun.run();
   if (_initialized != 0) return;
   _initialized = 1;
//-----------------------------------------------------------------------------
// as a last step, book histograms - need to know the number of active links
//-----------------------------------------------------------------------------
   book_histograms(rn);
  
}

//--------------------------------------------------------------------------------
  void TCalStat::beginJob() {
  }

//-----------------------------------------------------------------------------
  void TCalStat::endJob() {
    printf("[mu2e::TCalStat] pointer to the module: 0x%8p\n",(void*) this);
  }

//-----------------------------------------------------------------------------
  void TCalStat::fill_channel_histograms(ChannelHist_t* Hist, ChannelData_t* Data) {
    
      Hist->nhits->Fill(Data->nhits);
    
  }

//-----------------------------------------------------------------------------
  void TCalStat::fill_roc_histograms(RocHist_t* Hist, RocData_t* Rd) {
    
    Hist->nhits->Fill   (Rd->nhits);
    for (int ich=0; ich<kNChannels; ich++) {
      ChannelData_t* chd = &Rd->channel[ich];
      ChannelHist_t* hch = &Hist->channel[ich];
      hch->nhits->Fill(chd->nhits);
    }
      }
      
 

  void TCalStat::fill_event_histograms(EventHist_t* Hist, EventData_t* Data) {

  
    Hist->nhits->Fill(Data->nhtot);

  }

//-----------------------------------------------------------------------------
// fill_roc_histograms also fills the channel histograms
// if in error, only histogram the error code
//-----------------------------------------------------------------------------


int TCalStat::fill_histograms() {

   
   
    fill_event_histograms(&_Hist.event,&_event_data);

    for (int ir=0; ir<_nActiveLinks; ir++) {
      int link = _activeLinks[ir];
      fill_roc_histograms(&_Hist.roc[link],&_event_data.rdata[link]);
    }

    return 0;
    }
//-----------------------------------------------------------------------------
// a fragment may have multiple ROC blocks
//-----------------------------------------------------------------------------
  void TCalStat::analyze_roc(const art::Event& Evt) {

    return;
  }

//--------------------------------------------------------------------------------
// assume that we only have tracker fragment(s)
//-----------------------------------------------------------------------------
void TCalStat::analyze(const art::Event& event) {
  //art::EventNumber_t eventNumber = event.event();

  _event = &event;


  _event_data.nhtot = 0;
  for(int i=0;i<_nActiveLinks;i++){
    _event_data.rdata[i].nhits=0;
    for(int straw=0;straw<kNChannels;straw++){
     _event_data.rdata[i].channel[straw].nhits=0;
    }
  }
mu2e::GeomHandle<mu2e::Tracker> th;
    _tracker = th.get();
    //const Straw* straw = &_tracker->getStraw(sh->strawId());
    //const XYZVectorF* hp = &fStrawHitColl->at(jclosest).pos();
  /*const Straw& straw = _tracker->getStraw(step->strawId());

      plane = straw.id().getPlane();
      panel = straw.id().getPanel();
      layer = straw.id().getLayer();
      ist   = straw.id().getStraw();

      if (i < 0) printf("%5i %5i %5i %5i \n",plane, panel,layer,ist);

      const CLHEP::Hep3Vector& smp = straw.getMidPoint();
      const CLHEP::Hep3Vector& dir = straw.getDirection();

      double x  = step->position().x();
      double y  = step->position().y();
      double dx = x-smp.x();
      double dy = y-smp.y();

      double rho  = dx*smp.unit().y()-dy*smp.unit().x();
      double rho1 = dx*dir.x()+dy*dir.y();

    art::Handle<mu2e::StrawGasStepCollection> sgscH;
    ok = rEvent.getByLabel(_sgsCollTag,sgscH); 
    if (ok) {
      _sgsColl   = sgscH.product();
      _data.nsgs = _sgsColl->size();
    }
    else {
      _sgsColl  = nullptr;
      mf::LogWarning(oname) << " WARNING in " << oname << ":" << __LINE__ 
                            << ": StrawGasStepCollection:" 
                            << _sgsCollTag.encode().data() << " NOT FOUND";
    }
 void GenpHist::fillSgsHistograms(SgsHist_t* Hist, const StrawGasStep* Sgs) {
    Hist->strawId->Fill(Sgs->strawId().asUint16());
    Hist->edep->Fill(Sgs->ionizingEdep());
    Hist->time->Fill(Sgs->time());
  }
  */
  //  auto handle = event.getValidHandle<StrawDigiMCCollection>(_trkfCollTag);
  //auto handle = event.getValidHandle<StrawDigiMCCollection>(_trkfCollTag);
  auto handle = event.getValidHandle<StrawGasStepCollection>(_trkfCollTag);
 
  // CLASSE CHE HA IL FILE mu2e::StrawGasSteps
  //std::vector<mu2e::StrawDigiMC::StrawId> straw;
  //const StrawDigiMCCollection& mccol(*handle);
  const std::vector<mu2e::StrawGasStep>* mccol=handle.product(); //ora dovrei avere la classe
  //if(mccol==nullptr){
  //std::cout<<"No Straw Gas Step found"<<std::endl;
  //  return;
  //}
  // auto const& straws;
  if(mccol->size()>0){
    _event_data.nhtot=_event_data.nhtot+mccol->size();
  }
  for (size_t i=0;i<mccol->size();i++) {
    StrawGasStep mcsgs = mccol->at(i);
    
    
    uint16_t panel = mcsgs.strawId().getPanel();
    uint16_t straw = mcsgs.strawId().getStraw();
    std::cout<<"straw:"<<straw<<" panel: "<<panel<<std:endl;
    _event_data.rdata[panel].nhits += 1;
    _event_data.rdata[panel].channel[straw].nhits += 1;
    std::cout<<"nhits per straw::"<< _event_data.rdata[panel].channel[straw].nhits<<" nhits per panel: "<< _event_data.rdata[panel].nhits  <<std::endl;
    //auto const& sgs=mcdigi.strawGasSteps();
    //const mu2e::StrawDigiMC* sdmc = &mccol->at(i);
    //const StrawGasStep* sgs = sdmc->strawGasSteps();

    //straw.push_back(mccol->strawId()->get());
    // std::cout << straw.at(i) << std::end;
  }
  if ( _analyzeStation == 1 )  fill_histograms();

  TModule::analyze(event);
}

} // end namespace mu2e

DEFINE_ART_MODULE(mu2e::TCalStat)
