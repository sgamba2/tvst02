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
#include "Offline/DataProducts/inc/StrawEnd.hh"
#include "Offline/DataProducts/inc/StrawId.hh"
#include "Offline/DataProducts/inc/StrawIdMask.hh"
#include "Offline/DataProducts/inc/GenVector.hh"
#include "Offline/DataProducts/inc/TrkTypes.hh"
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/RecoDataProducts/inc/StrawHitIndex.hh"
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
#include "Offline/RecoDataProducts/inc/ComboHit.hh"

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
    _trkfCollTag     (PSet.get<art::InputTag>   ("trkfCollTag")),
    _analyzeStation  (PSet.get<int>             ("analyzeStation")),
    _digiCollTag     (PSet.get<art::InputTag>   ("digiCollTag")),
    _recCollTag      (PSet.get<art::InputTag>   ("recCollTag")),
    _comboCollTag    (PSet.get<art::InputTag>   ("comboCollTag"))

  {
    _initialized     = 0;
  }
   
 
  void TCalStat::book_channel_histograms(art::TFileDirectory* Dir, int RunNumber, ChannelHist_t* Hist, int I) {

    Hist->nhits   = Dir->make<TH1F>(Form("ch_%02i_nhits",I),Form("run %06i: ch %02i nhits"  ,RunNumber,I), 30, 0., 30.);
       
  }

  void TCalStat::book_panel_histograms(art::TFileDirectory* Dir, int RunNumber, PanelHist_t* Hist, int pan) {
   
    Hist->nhits           = Dir->make<TH1F>(Form("panel_%02i_nhits", pan)   , Form("run %06i, panel%02i: n hits"      ,RunNumber,pan),  30,    0.,   30.);
    
    Hist->nh_vs_ch        = Dir->make<TH2F>(Form("panel_%02i_nh_vs_ch", pan)   , Form("run %06i, panel%02i: nh vs ch"  ,RunNumber,pan),  100,0.,100., 10,0,10);
   
    Hist->xp_vs_yp        = Dir->make<TH2F>( Form("panel_%02i_xp_vs_yp", pan) , Form("run %06i, panel%02i: xp vs yp"  ,RunNumber,pan),  2000,-1000.,+1000., 1000,0.,1000.);

    Hist->xp_vs_yp_nh3       = Dir->make<TH2F>( Form("panel_%02i_xp_vs_yp_nh3", pan) , Form("run %06i, panel%02i: xp vs yp nhits<=3"  ,RunNumber,pan),  2000,-1000.,+1000., 1000,0.,1000.);

    Hist->x_vs_y        = Dir->make<TH2F>( Form("panel_%02i_x_vs_y", pan) , Form("run %06i, panel%02i: track"  ,RunNumber,pan),  2000,-1000.,+1000., 2000,-1000.,1000.);
    Hist->x_bias   = Dir->make<TH1F>(Form("panel_%02i_x_bias",pan),Form("run %06i: panel %02i xrec-xtrue bias"  ,RunNumber,pan), 600, -300., 300.);
    Hist->y_bias   = Dir->make<TH1F>(Form("panel_%02i_y_bias",pan),Form("run %06i: panel %02i yrec-ytrue bias"  ,RunNumber,pan), 600, -300., 300.);
    Hist->x  = Dir->make<TH1F>(Form("panel_%02i_x",pan),Form("run %06i: panel %02i xrec"  ,RunNumber,pan), 2000, -1000., 1000.);
    Hist->y  = Dir->make<TH1F>(Form("panel_%02i_y",pan),Form("run %06i: panel %02i yrec"  ,RunNumber,pan), 2000, -1000., 1000.);
  Hist->vdrift      = Dir->make<TH2F>( Form("panel_%02i_vdrift", pan) , Form("run %06i, panel%02i: vdrift"  ,RunNumber,pan),  2000,-1000.,+1000., 2000,-100.,100.);
Hist->vdrift_true      = Dir->make<TH2F>( Form("panel_%02i_vdrift_true", pan) , Form("run %06i, panel%02i: vdrift true pos"  ,RunNumber,pan),  2000,-1000.,+1000., 2000,-100.,100.);
 Hist->x_bias_stereo_combo   = Dir->make<TH1F>(Form("panel_%02i_x_bias_stereo_combo",pan),Form("run %06i: panel %02i bias stereo combo"  ,RunNumber,pan), 600, -300., 300.);
    Hist->y_bias_stereo_combo   = Dir->make<TH1F>(Form("panel_%02i_y_bias_stereo_combo",pan),Form("run %06i: panel %02i bias stereo combo"  ,RunNumber,pan), 600, -300., 300.);
    Hist->x_bias_vs_x          = Dir->make<TH2F>( Form("panel_%02i_x_bias_vs_x",pan),  Form("run %06i: panel %02i xbias vs x"  ,RunNumber, pan),  2000,-1000.,+1000., 2000,-1000.,1000.);
Hist->x_bias_vs_x_true          = Dir->make<TH2F>( Form("panel_%02i_x_bias_vs_x_true",pan),  Form("run %06i: panel %02i xbias vs xtrue"  ,RunNumber, pan),  2000,-1000.,+1000., 2000,-1000.,1000.);
    for (int i=0; i<kNChannels; i++) {
      art::TFileDirectory chan_dir = Dir->mkdir(Form("ch_%02i",i));
      book_channel_histograms(&chan_dir,RunNumber,&Hist->channel[i],i);
    }
    
    
  }

  void TCalStat::book_plane_histograms(art::TFileDirectory* Dir, int RunNumber, PlaneHist_t* Hist, int pl) {
    
    Hist->nhits   = Dir->make<TH1F>(Form("plane_%02i_nhits",pl),Form("run %06i: plane %02i nhits"  ,RunNumber,pl), 30, 0., 30.);
    Hist->x_vs_y        = Dir->make<TH2F>( Form("plane_%02i_x_vs_y", pl) , Form("run %06i, plane%02i: track"  ,RunNumber,pl), 2000,-1000.,+1000., 2000,-1000.,1000.);
    for (int i=0; i<kNPanels; i++) {
      art::TFileDirectory panel_dir = Dir->mkdir(Form("panel_%02i",i));
      book_panel_histograms(&panel_dir,RunNumber,&Hist->panel[i],i);
    } 
  }

  void TCalStat::book_event_histograms(art::TFileDirectory* Dir, int RunNumber, EventHist_t* Hist) {
    
    Hist->nhits           = Dir->make<TH1F>("nhits"      , Form("run %06i: nhits total"  ,RunNumber), 1500, 0.,   1500.);
    Hist->x_vs_y          = Dir->make<TH2F>( "x_vs_y" , Form("run %06i, track"  ,RunNumber),  2000,-1000.,+1000., 2000,-1000.,1000.);
    Hist->x_vs_y_nh3      = Dir->make<TH2F>( "x_vs_y_nh3" , Form("run %06i, track nhits per panel<=3"  ,RunNumber),  2000,-1000.,+1000., 2000,-1000.,1000.);
    Hist->x_vs_y_nh3_ex   = Dir->make<TH2F>( "x_vs_y_nh3_ex" , Form("run %06i, track, nhits per panel<=3, extreme removed"  ,RunNumber),  2000,-1000.,+1000., 2000,-1000.,1000.);
    Hist->nhits_nh3_ex    = Dir->make<TH1F>("nhits_nh3_ex"      , Form("run %06i: nhits total,nhits per panel<=3, extreme removed"  ,RunNumber), 30, 0.,   30.);
    Hist->nhits_nh3       = Dir->make<TH1F>("nhits_nh3"      , Form("run %06i: nhits totalnhits per panel<=3"  ,RunNumber), 30, 0.,   30.);
    Hist->x_vs_y_rec_p    = Dir->make<TH2F>( "x_vs_y_rec_p" , Form("run %06i, track, nhits per panel<=3,reco"  ,RunNumber),  2000,-1000.,+1000., 2000,-1000.,1000.);
    Hist->x_vs_y_rec_pf   = Dir->make<TH2F>( "x_vs_y_rec_pf" , Form("run %06i, track, nhits per panel<=3,reco fit"  ,RunNumber),  2000,-1000.,+1000., 2000,-1000.,1000.);
    Hist->t0              = Dir->make<TH1F>("tmean"      , Form("run %06i: (tcal+thv)/2 digi"  ,RunNumber), 100000, 0.,   100000.);
    Hist->costheta              = Dir->make<TH1F>("costheta"      , Form("run %06i: costheta"  ,RunNumber), 2000 , -1.,   1.);
 Hist->x_vs_y_rec_true          = Dir->make<TH2F>( "x_vs_y_rec_true" , Form("run %06i, track rec true"  ,RunNumber),  2000,-1000.,+1000., 2000,-1000.,1000.);

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
      Hist->nhits->Fill(Data->nhits);
    
  }


  void TCalStat::fill_panel_histograms(PanelHist_t* Hist, PanelData_t* Rd) {
    Hist->nhits->Fill   (Rd->nhits);
 
    for (int ich=0; ich<kNChannels; ich++) {

      ChannelData_t* chd = &Rd->channel[ich];
      ChannelHist_t* hch = &Hist->channel[ich];
      Hist->nh_vs_ch->Fill(ich,chd->nhits);
      fill_channel_histograms(hch,chd);
    }
  }
  
  void TCalStat::fill_plane_histograms(PlaneHist_t* Hist, PlaneData_t* Rd) {
    Hist->nhits->Fill   (Rd->nhits);
    for (int i=0; i<kNPanels; i++) {
      PanelData_t* chd = &Rd->panel[i];
      PanelHist_t* hch = &Hist->panel[i];
      fill_panel_histograms(hch,chd);
    }
  }
  
  void TCalStat::fill_event_histograms(EventHist_t* Hist, EventData_t* Data) {
    Hist->nhits->Fill(Data->nhits);
    if(Data->nhits_nh3_ex>0){
      Hist->nhits_nh3_ex->Fill(Data->nhits_nh3_ex);
    }
    if(Data->nhits_nh3>0){
      Hist->nhits_nh3->Fill(Data->nhits_nh3);
    }  
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

  /*  void TCalStat::fit_line(double * m, double * sigma_m, double * q, double * sigma_q, double x[], double y[], double dx[], double dy[], int size, double m1){
    double sx0=0;
    double sx1=0;
    double sx2=0;
    double sxy0=0;
    double sxy1=0;
    for(int i=0;i<size;i++ ){
      double sigma_y=TMath::Sqrt(dy[i]*dy[i]+(m1*dx[i])*(m1*dx[i]));
      sx0=sx0+(1/(sigma_y*sigma_y));
      sx1=sx1+(x[i]/(sigma_y*sigma_y));
      sx2=sx2+((x[i]*x[i])/(sigma_y*sigma_y));
      sxy0=sxy0+(y[i]/(sigma_y*sigma_y));
      sxy1=sxy1+((x[i]*y[i])/(sigma_y*sigma_y));
    }
    double D=sx0*sx2-sx1*sx1;
    *m=(sxy1*sx0-sxy0*sx1)/D;
    *q=(sxy0*sx2-sxy1*sx1)/D;
    *sigma_m=TMath::Sqrt(sx0/D);
    *sigma_q=TMath::Sqrt(sx2/D);
    return;
    }*/

void TCalStat::analyze(const art::Event& event) {

  _event = &event;
  mu2e::GeomHandle<mu2e::Tracker> th;
  _tracker = th.get();
    auto handle_comb = event.getValidHandle<ComboHitCollection>(_comboCollTag);
  const ComboHitCollection& combocol(*handle_comb);
  auto handle = event.getValidHandle<StrawGasStepCollection>(_trkfCollTag);
  auto handle_dig = event.getValidHandle<StrawDigiMCCollection>(_digiCollTag);
 auto handle_rec = event.getValidHandle<ComboHitCollection>(_recCollTag);
  const ComboHitCollection& reccol(*handle_rec);
  const std::vector<mu2e::StrawGasStep>* mccol=handle.product(); 
  const std::vector<mu2e::StrawDigiMC>* digicol=handle_dig.product(); 
  if(digicol->size() > 0){
    for(size_t i =0 ; i<digicol->size(); i++){
      StrawDigiMC mcsgs = digicol->at(i);
      double tmean = (mcsgs.wireEndTime(StrawEnd::cal)+mcsgs.wireEndTime(StrawEnd::hv))/2.;
      _Hist.event.t0->Fill(tmean);
    }
  }
  size_t mccolsize=mccol->size();
  if(fInteractiveMode==1){
    _Hist.event.x_vs_y->Reset();
    _Hist.event.x_vs_y_rec_p->Reset();
  }
  _event_data.nhits = 0;
  _event_data.nhits_nh3=0;
  _event_data.nhits_nh3_ex=0;

  for(int plane=0;plane<kNPlanes;plane++){
    _event_data.plane[plane].nhits=0;
    if(fInteractiveMode==1){
      _Hist.plane[plane].x_vs_y->Reset();
    }
    for(int panel=0;panel<kNPanels;panel++){
      _event_data.plane[plane].panel[panel].nhits=0;
      if(fInteractiveMode==1){
	_Hist.plane[plane].panel[panel].x_vs_y->Reset();
      }
      for(int straw=0;straw<kNChannels;straw++){
	_event_data.plane[plane].panel[panel].channel[straw].nhits=0;
      }
    }
  }
  if (combocol.size() != 4) return;
  if(reccol.size() != 2 ) return;

  std::vector<int> planev;

  for(size_t i=0;i<mccolsize;i++){

    StrawGasStep mcsgs = mccol->at(i);

    if(mcsgs.strawId().getStation()==0){
      planev.push_back(mcsgs.strawId().getPlane());

    }

  }

  int check_true=0;
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
  double yt[12];
  double xt[12];
  //uint16_t pplane[12];
  //uint16_t ppanel[12];
  int nhits_true[12];
  for(int i=0; i<12;i++){
    yt[i]=0;
    xt[i]=0;
    //pplane[i]=999;
    //ppanel[i]=999;
    nhits_true[i]=0;

  }
  for(size_t i=0;i<mccolsize;i++){

    StrawGasStep mcsgs = mccol->at(i);
    uint16_t panel=mcsgs.strawId().getPanel();

    if(mcsgs.strawId().getStation()==0){ 

      if(planev.at(i)==0){
	if(panel==0) npan_00++;
	if(panel==1) npan_01++;
	if(panel==2) npan_02++;
	if(panel==3) npan_03++;
	if(panel==4) npan_04++;
	if(panel==5) npan_05++;
      }

      if(planev.at(i)==1){
	if(panel==0) npan_10++;
	if(panel==1) npan_11++;
	if(panel==2) npan_12++;
	if(panel==3) npan_13++;
	if(panel==4) npan_14++;
	if(panel==5) npan_15++;
      }
    }
  }
  
 
  _event_data.nhits=_event_data.nhits+mccol->size();
  int true_size=0;
  for (size_t i=0;i<mccolsize;i++) {

    StrawGasStep mcsgs = mccol->at(i);
    uint16_t panel = mcsgs.strawId().getPanel();
    uint16_t plane = mcsgs.strawId().getPlane();
    uint16_t straw = mcsgs.strawId().getStraw();
    double x    = mcsgs.position().x();
    double y    = mcsgs.position().y();

    if(mcsgs.strawId().getStation()==0){ 
      
      _Hist.event.x_vs_y->Fill(x,y);
      _Hist.plane[plane].x_vs_y->Fill(x,y);
      _Hist.plane[plane].panel[panel].x_vs_y->Fill(x,y);
    
    
      const Plane* pln = &_tracker->getPlane(plane);
      const Panel* pnl = &pln->getPanel(panel);
      const Straw* strw = &pnl->getStraw(straw);
      double costheta = strw->direction().x();
      double sintheta = strw->direction().y();
      double xp=x*costheta+y*sintheta;
      double yp=-y*costheta+x*sintheta;
   
    
      if(yp<0) yp=yp*(-1.);
  
  
      _Hist.plane[plane].panel[panel].xp_vs_yp->Fill(xp,yp);
      _event_data.plane[plane].panel[panel].nhits +=1;
      _event_data.plane[plane].panel[panel].channel[straw].nhits += 1;
      _event_data.plane[plane].nhits +=1;

    
      if(npan_00<=3 and plane==0 and panel==0){
	_Hist.plane[plane].panel[panel].xp_vs_yp_nh3->Fill(xp,yp);       
      }
      if(npan_01<=3 and plane==0 and panel==1){
	_Hist.plane[plane].panel[panel].xp_vs_yp_nh3->Fill(xp,yp);       
      }
      if(npan_02<=3 and plane==0 and panel==2){
	_Hist.plane[plane].panel[panel].xp_vs_yp_nh3->Fill(xp,yp);       
      }
      if(npan_03<=3 and plane==0 and panel==3){
	_Hist.plane[plane].panel[panel].xp_vs_yp_nh3->Fill(xp,yp);       
      }
      if(npan_04<=3 and plane==0 and panel==4){
	_Hist.plane[plane].panel[panel].xp_vs_yp_nh3->Fill(xp,yp);       
      }
      if(npan_05<=3 and plane==0 and panel==5){
	_Hist.plane[plane].panel[panel].xp_vs_yp_nh3->Fill(xp,yp);       
      }
      if(npan_10<=3 and plane==1 and panel==0){
	_Hist.plane[plane].panel[panel].xp_vs_yp_nh3->Fill(xp,yp);       
      }
      if(npan_11<=3 and plane==1 and panel==1){
	_Hist.plane[plane].panel[panel].xp_vs_yp_nh3->Fill(xp,yp);       
      }
      if(npan_12<=3 and plane==1 and panel==2){
	_Hist.plane[plane].panel[panel].xp_vs_yp_nh3->Fill(xp,yp);       
      }
      if(npan_13<=3 and plane==1 and panel==3){
	_Hist.plane[plane].panel[panel].xp_vs_yp_nh3->Fill(xp,yp);       
      }
      if(npan_14<=3 and plane==1 and panel==4){
	_Hist.plane[plane].panel[panel].xp_vs_yp_nh3->Fill(xp,yp);       
      }
      if(npan_15<=3 and plane==1 and panel==5){
	_Hist.plane[plane].panel[panel].xp_vs_yp_nh3->Fill(xp,yp);       
      }
      if(npan_15<=3 and npan_14<=3 and npan_13<=3 and npan_12<=3 and npan_11<=3 and npan_10<=3 and npan_05<=3 and npan_04<=3 and npan_03<=3 and npan_02<=3 and npan_01<=3 and npan_00<=3 and xp<300 and xp>-300 and yp>300 and yp<500 ){
	_Hist.event.x_vs_y_nh3_ex->Fill(x,y);
	_event_data.nhits_nh3_ex++;
      }
      if(npan_15<=3 and npan_14<=3 and npan_13<=3 and npan_12<=3 and npan_11<=3 and npan_10<=3 and npan_05<=3 and npan_04<=3 and npan_03<=3 and npan_02<=3 and npan_01<=3 and npan_00<=3){
	_Hist.event.x_vs_y_nh3->Fill(x,y);   
	ComboHit ch=combocol[true_size];
	if(combocol[0].strawId().getPanel() != panel and combocol[1].strawId().getPanel() != panel and combocol[2].strawId().getPanel() != panel and combocol[3].strawId().getPanel() != panel and combocol[0].strawId().getPlane()!= plane and combocol[1].strawId().getPlane()!= plane and combocol[2].strawId().getPlane()!= plane and combocol[3].strawId().getPlane()!= plane){
	  check_true=1;
	}
	if((ch.strawId().getPanel()!=panel or ch.strawId().getPlane()!=plane) and check_true!=1){
	  true_size++;
	}
	if(ch.strawId().getPanel()==panel and ch.strawId().getPlane()==plane and check_true!=1){
	  xt[true_size]+=x;
	  yt[true_size]+=y;
	  //ppanel[true_size]=panel;
	  //pplane[true_size]=plane;
	  nhits_true[true_size]++;
	}
      }
    }
  }
  for(int i=0;i<12;i++){
    xt[i] /= nhits_true[i];
    yt[i] /= nhits_true[i];

  }
  if(npan_15<=3 and npan_14<=3 and npan_13<=3 and npan_12<=3 and npan_11<=3 and npan_10<=3 and npan_05<=3 and npan_04<=3 and npan_03<=3 and npan_02<=3 and npan_01<=3 and npan_00<=3  ){
    _event_data.nhits_nh3+=mccolsize;
  }
  
 
  //ora devi raccattare i due punti (x.y.z) e trovare le 2 rette xy e yz
  double x1(0), y1(0), z1(0), x2(0), y2(0), z2(0);
  ComboHit chit1 = reccol[0];
  ComboHit chit2 = reccol[1];
  int plane1=chit1._uvar;
  int plane2=chit2._vvar;

  //posizione vera, condivisa tra 2 panels
  x1=chit1._pos.x();
  y1=chit1._pos.y();
  z1=chit1._pos.z();
  x2=chit2._pos.x();
  y2=chit2._pos.y();
  z2=chit2._pos.z(); //z medio okay
  
  //ricostruzione
  //xy:
  double mxy=(y2-y1)/(x2-x1);
  double qxy=-mxy*x1+y1;
  //yz:
  double myz=(y2-y1)/(z2-z1);
  double qyz=-myz*z1+y1;

  //ora devo trovare i punti di incontro con i pannelli delle combohit


  double z[4];
  double x[4];
  double y[4];
  double costheta[4];
  double sintheta[4];
  double xp[4];
  double yp[4];
  int PL[4];
  int PA[4];
  double thv[4];
  double tcal[4];

  for(int i=0;i<4;i++){
    z[i]=0.;
    x[i]=0.;
    y[i]=0.;
    xp[i]=0.;
    yp[i]=0.;
    costheta[i]=0.;
    sintheta[i]=0.;
    PL[i]=0;
    PA[i]=0;
    thv[i]=0.;
    tcal[i]=0.;
  }


  for(size_t i=0; i<4;i++){
    ComboHit ch = combocol[i];
    z[i]=ch._pos.z();
    //posizione ricostruita
    y[i]=myz*z[i]+qyz;
    x[i]=(myz*z[i]+qyz-qxy)/mxy;
    costheta[i] = ch._udir.x();
    sintheta[i] = ch._udir.y();
    xp[i]=x[i]*costheta[i]+y[i]*sintheta[i];
    yp[i]=-y[i]*costheta[i]+x[i]*sintheta[i];
    PA[i]=ch.strawId().getPanel();
    PL[i]=ch.strawId().getPlane();
    thv[i]=ch._etime[StrawEnd::hv];
    tcal[i]=ch._etime[StrawEnd::cal];
  }
  
  for(int i =0; i<4 ; i++){
    ComboHit ch = combocol[i];
    double xp_t=xt[i]*costheta[i]+yt[i]*sintheta[i];
    double yp_t=-yt[i]*costheta[i]+xt[i]*sintheta[i];
    _Hist.plane[PL[i]].panel[PA[i]].x_bias->Fill(xp[i]-xp_t);
    _Hist.plane[PL[i]].panel[PA[i]].y_bias->Fill(yp[i]-yp_t);
    _Hist.plane[PL[i]].panel[PA[i]].x_bias_vs_x->Fill(xp[i],xp[i]-xp_t);
    _Hist.plane[PL[i]].panel[PA[i]].x_bias_vs_x_true->Fill(xp_t,xp[i]-xp_t);
    _Hist.plane[PL[i]].panel[PA[i]].vdrift_true->Fill(xp_t,thv[i]-tcal[i]);
    _Hist.plane[PL[i]].panel[PA[i]].x->Fill(xp[i]);
    _Hist.plane[PL[i]].panel[PA[i]].y->Fill(yp[i]);
    _Hist.event.x_vs_y_rec_p->Fill(x[i],y[i]);
    _Hist.plane[PL[i]].panel[PA[i]].vdrift->Fill(xp[i],thv[i]-tcal[i]);
    
    if((int)PL[i]==plane1){
      double xp1=x1*costheta[i]+y1*sintheta[i];
      double yp1=-y1*costheta[i]+x1*sintheta[i];
      _Hist.plane[PL[i]].panel[PA[i]].x_bias_stereo_combo->Fill(ch._pos.x()*costheta[i]+ch._pos.y()*sintheta[i]-xp1);
      _Hist.plane[PL[i]].panel[PA[i]].y_bias_stereo_combo->Fill(-ch._pos.y()*costheta[i]+ch._pos.x()*sintheta[i]-yp1);

    }
    if((int)PL[i]==plane2){
      double xp2=x2*costheta[i]+y2*sintheta[i];
      double yp2=-y2*costheta[i]+x2*sintheta[i];
      _Hist.plane[PL[i]].panel[PA[i]].x_bias_stereo_combo->Fill(ch._pos.x()*costheta[i]+ch._pos.y()*sintheta[i]-xp2);
      _Hist.plane[PL[i]].panel[PA[i]].y_bias_stereo_combo->Fill(-ch._pos.y()*costheta[i]+ch._pos.x()*sintheta[i]-yp2);
    }
   

  }
    _Hist.event.costheta->Fill(TMath::ATan(myz));
  if(npan_15<=2 and npan_14<=2 and npan_13<=2 and npan_12<=2 and npan_11<=2 and npan_10<=2 and npan_05<=2 and npan_04<=2 and npan_03<=2 and npan_02<=2 and npan_01<=2 and npan_00<=2 and  fInteractiveMode==1  ){
    std::cout<<"MC points:"<<std::endl;
     for(size_t i =0 ; i<digicol->size(); i++){
      StrawDigiMC mcsgs = digicol->at(i);
      const Plane* pl1 = &_tracker->getPlane(mcsgs.strawId().getPlane());
      const Panel* pa1 = &pl1->getPanel(mcsgs.strawId().getPanel());
      const Straw* strw = &pa1->getStraw(mcsgs.strawId().getStraw());
      std::cout<<"StrawId plane  panel face   dirx   diry  dirz  midx   midy  midz"<<std::endl;
      std::cout<<mcsgs.strawId().getStraw()<<"  "<<mcsgs.strawId().getPlane()<<"  "<<mcsgs.strawId().getPanel()<<" "<<pa1->getStraw(0).id().face()<<"  "<<strw->getDirection()<<"  "<<strw->getMidPoint()<<std::endl;

    }
    std::cout<<"COMBO points:"<<std::endl;
    for(size_t i =0 ; i<combocol.size(); i++){
      ComboHit mcsgs = combocol[i];
      const Plane* pl1 = &_tracker->getPlane(mcsgs.strawId().getPlane());
      const Panel* pa1 = &pl1->getPanel(mcsgs.strawId().getPanel());
      std::cout<<"plane panel face dirx   diry  dirz  midx   midy  midz"<<std::endl;
      std::cout<<mcsgs.strawId().getPlane()<<"  "<<mcsgs.strawId().getPanel()<<" "<<pa1->getStraw(0).id().face()<<"  "<<mcsgs._udir<<"  "<<mcsgs._pos<<std::endl;

    }
    std::cout<<"STEREO points:"<<std::endl;
    for(size_t i =0 ; i<reccol.size(); i++){
      ComboHit mcsgs = reccol[i];
      std::cout<<"plane x  y  z "<<std::endl;
      std::cout<<mcsgs.strawId().getPlane()<<"  "<<mcsgs._udir<<"  "<<mcsgs._pos<<std::endl;

    }
    std::cout<<"STRAIGHT Line:"<<std::endl;
    
    std::cout<<"y= "<<mxy<<"*x+"<<qxy<<std::endl;
    std::cout<<"y= "<<myz<<"*z+"<<qyz<<std::endl;
    std::cout<<"Reconstructed position:"<<std::endl;

    for(size_t i=0; i<4;i++){
      std::cout<<"plane  panel  xrec  yrec  zrec  dirx  diry"<<std::endl;
      std::cout<<PL[i]<<"  "<<PA[i]<<"  "<<x[i]<<"  "<<y[i]<<"  "<<z[i]<<"  "<<costheta[i]<<"  "<<sintheta[i]<<std::endl;
  }


  }




    //_Hist.event.x_vs_y_rec_true->Fill(x1,y1);
    //_Hist.event.x_vs_y_rec_true->Fill(x2,y2);

  //ora devo trovare il punto di intersezione tra i pannelli trovati e la retta ricostruita
  //i punti di incontro sono definiti dalla z del panel
  // posso trovare y sapendo che y=myz * z +qyz
  //so pero che y=mxy x+qxy
  //mxy x +qxy=myz z +qyz
  //x=(myz*z+qyz-qxy)/mxy
  

 /*   
    double xc[1500];
    double yc[1500];
    double midm_x[1500];
    double midm_y[1500];
    double dirm_x[1500];
    double dirm_y[1500];
    double dxc_c[1500];
    double dyc_c[1500];
    double stdmid_y[1500];
    double stdmid_x[1500];
    double dx[1500];
    double dy[1500];
    double xc_c[1500];
    double yc_c[1500];
    int tot_corr=0;
    for(int i=0;i<good_hits;i++){
      xc[i]=0;
      yc[i]=0;
      xc_c[i]=0;
      yc_c[i]=0;
      dxc_c[i]=0;
      dyc_c[i]=0;
      midm_x[i]=0;
      midm_y[i]=0;
      dirm_x[i]=0;
      dirm_y[i]=0;
      dx[i]=0;
      dy[i]=0;
      stdmid_y[i]=0;
      stdmid_x[i]=0;
    }

    //mean value of directions and midpoints in one panel
    int tot=0;
    int k=0;
    while(k<good_hits){
      int pa = PANEL[k];
      int pl= PLANE[k];
      int count=0;
      while(PANEL[k+count]==pa and PLANE[k+count]==pl){
	midm_x[tot]+=mid_x[k+count];
	midm_y[tot]+=mid_y[k+count];
	dirm_y[tot]+=dir_y[k+count];
	dirm_x[tot]+=dir_x[k+count];
	count++;
      }
      midm_x[tot]/=(count);
      midm_y[tot]/=(count);
      dirm_x[tot]/=(count);
      dirm_y[tot]/=(count);
    
      k=k+count;
      tot++;
    }

    //std deviation
    int u=0;
    k=0;
    while(k<good_hits){
      int pa = PANEL[k];
      int pl= PLANE[k];
      int count=0;
      while(PANEL[k+count]==pa and PLANE[k+count]==pl){
	stdmid_x[u]+=TMath::Power(midm_x[u]-mid_x[k+count],2);
	stdmid_y[u]+=TMath::Power(midm_y[u]-mid_y[k+count],2);
	count++;
      }
      if(count>1){
	stdmid_x[u]/=(count*(count-1));
	stdmid_y[u]/=(count*(count-1));
      }
      k=k+count;
      u++;
    }
    //finding intersection points between a straw and the next one	  
    for(int i=0;i<(tot-1);i++ ){ 

      double m1=dirm_y[i]/dirm_x[i]; 
      double q1=-m1*midm_x[i]+midm_y[i];
      double err_mult1=m1*stdmid_x[i];
      double dq1=TMath::Sqrt(TMath::Power(err_mult1,2)+TMath::Power(stdmid_y[i],2));
      double m2=dirm_y[i+1]/dirm_x[i+1];
      double q2=-m2*midm_x[i+1]+midm_y[i+1];
      double err_mult2=m2*stdmid_x[i+1];
      double dq2=TMath::Sqrt(TMath::Power(err_mult2,2)+TMath::Power(stdmid_y[i+1],2));
      xc[i]=(q1-q2)/(m2-m1);
      double err_sum1=TMath::Sqrt(TMath::Power(dq1,2)+TMath::Power(dq2,2));
      dx[i]=(err_sum1)/(m2-m1);
      yc[i]=m2*xc[i]+q2;
      double err_mult3=m2*dx[i];
      dy[i]=TMath::Sqrt(TMath::Power(err_mult3,2)+TMath::Power(dq2,2));
      if(TMath::Sqrt(TMath::Power(xc[i],2)+TMath::Power(yc[i],2))<700){
	xc_c[tot_corr]=xc[i];
	yc_c[tot_corr]=yc[i];
	dyc_c[tot_corr]=dy[i];
	dxc_c[tot_corr]=dx[i];
	tot_corr++;
	_Hist.event.x_vs_y_rec_p->Fill(xc[i],yc[i]);
      }
    }
    //fitting the found points considering errors on x and y (effective error method)
    double m,sigma_m,q,sigma_q;
    double dx1[1500];
    for(int i=0;i<tot_corr;i++){
      dx1[i]=0;
    }
    fit_line(&m,&sigma_m,&q,&sigma_q,xc_c,yc_c,dx1,dyc_c,tot_corr,0);
    for(int i=0;i<5;i++){
      double m1=m;
      fit_line(&m,&sigma_m,&q,&sigma_q,xc_c,yc_c,dxc_c,dyc_c,tot_corr,m1);
    }
    
    //finding intersection points with straws
    for(int i=0;i<good_hits;i++ ){
      double m1=dir_y[i]/dir_x[i]; 
      double q1=-m1*mid_x[i]+mid_y[i];
      double xcf=(q1-q)/(m-m1);
      double dxcf=TMath::Sqrt(TMath::Power(sigma_q/(m-m1),2)+TMath::Power(sigma_m*(q1-q)/TMath::Power((m-m1),2),2));
      double ycf=m*xcf+q;
      double err_mult3=TMath::Sqrt(TMath::Power(sigma_m*xcf,2)+TMath::Power(m*dxcf,2));
      double dycf =TMath::Sqrt(TMath::Power(err_mult3,2)+TMath::Power(sigma_q,2));
      if(TMath::Sqrt(TMath::Power(xcf,2)+TMath::Power(ycf,2))<700){
	_Hist.event.x_vs_y_rec_pf->Fill(xcf,ycf);
	_Hist.plane[PLANE[i]].panel[PANEL[i]].channel[STRAW[i]].x_bias->Fill((xcf*dir_x[i]+ycf*dir_y[i])-xt[i]);
       //sistematics
	_Hist.plane[PLANE[i]].panel[PANEL[i]].channel[STRAW[i]].x->Fill(((xcf)*dir_x[i]+(ycf)*dir_y[i]));

	_Hist.plane[PLANE[i]].panel[PANEL[i]].channel[STRAW[i]].x_bias_p->Fill(((xcf+dxcf)*dir_x[i]+(ycf+dycf)*dir_y[i])-xt[i]);
	_Hist.plane[PLANE[i]].panel[PANEL[i]].channel[STRAW[i]].x_bias_m->Fill(((xcf-dxcf)*dir_x[i]+(ycf-dycf)*dir_y[i])-xt[i]);
      }
    } 
  }
  */
  if ( _analyzeStation == 1 and mccolsize>0)  fill_histograms();
  
  
  TModule::analyze(event);
}
}
  
 // end namespace mu2e

DEFINE_ART_MODULE(mu2e::TCalStat)
