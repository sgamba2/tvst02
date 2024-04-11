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

    Hist->nhits   = Dir->make<TH1F>(Form("ch_%02i_nhits",I),Form("run %06i: ch %02i nhits"  ,RunNumber,I), 30, 0., 30.);
        Hist->x_bias   = Dir->make<TH1F>(Form("ch_%02i_bias",I),Form("run %06i: ch %02i xrec-xtrue bias"  ,RunNumber,I), 600, -300., 300.);


  }

  void TCalStat::book_panel_histograms(art::TFileDirectory* Dir, int RunNumber, PanelHist_t* Hist, int pan) {
   
    Hist->nhits           = Dir->make<TH1F>(Form("panel_%02i_nhits", pan)   , Form("run %06i, panel%02i: n hits"      ,RunNumber,pan),  30,    0.,   30.);
    
    Hist->nh_vs_ch        = Dir->make<TH2F>(Form("panel_%02i_nh_vs_ch", pan)   , Form("run %06i, panel%02i: nh vs ch"  ,RunNumber,pan),  100,0.,100., 10,0,10);
   
    Hist->xp_vs_yp        = Dir->make<TH2F>( Form("panel_%02i_xp_vs_yp", pan) , Form("run %06i, panel%02i: xp vs yp"  ,RunNumber,pan),  2000,-1000.,+1000., 1000,0.,1000.);

    Hist->xp_vs_yp_nh3       = Dir->make<TH2F>( Form("panel_%02i_xp_vs_yp_nh3", pan) , Form("run %06i, panel%02i: xp vs yp nhits<=3"  ,RunNumber,pan),  2000,-1000.,+1000., 1000,0.,1000.);

 Hist->x_vs_y        = Dir->make<TH2F>( Form("panel_%02i_x_vs_y", pan) , Form("run %06i, panel%02i: track"  ,RunNumber,pan),  2000,-1000.,+1000., 2000,-1000.,1000.);
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
   
    Hist->nhits           = Dir->make<TH1F>("nhits"      , Form("run %06i: nhits total"  ,RunNumber), 30, 0.,   30.);
   Hist->x_vs_y        = Dir->make<TH2F>( "x_vs_y" , Form("run %06i, track"  ,RunNumber),  2000,-1000.,+1000., 2000,-1000.,1000.);
   Hist->x_vs_y_nh3        = Dir->make<TH2F>( "x_vs_y_nh3" , Form("run %06i, track nhits per panel<=3"  ,RunNumber),  2000,-1000.,+1000., 2000,-1000.,1000.);
   Hist->x_vs_y_nh3_ex       = Dir->make<TH2F>( "x_vs_y_nh3_ex" , Form("run %06i, track, nhits per panel<=3, extreme removed"  ,RunNumber),  2000,-1000.,+1000., 2000,-1000.,1000.);
   Hist->nhits_nh3_ex           = Dir->make<TH1F>("nhits_nh3_ex"      , Form("run %06i: nhits total,nhits per panel<=3, extreme removed"  ,RunNumber), 30, 0.,   30.);
    Hist->nhits_nh3           = Dir->make<TH1F>("nhits_nh3"      , Form("run %06i: nhits totalnhits per panel<=3"  ,RunNumber), 30, 0.,   30.);
   Hist->x_vs_y_rec_p       = Dir->make<TH2F>( "x_vs_y_rec_p" , Form("run %06i, track, nhits per panel<=3,reco"  ,RunNumber),  2000,-1000.,+1000., 2000,-1000.,1000.);
   Hist->x_vs_y_rec_pf       = Dir->make<TH2F>( "x_vs_y_rec_pf" , Form("run %06i, track, nhits per panel<=3,reco fit"  ,RunNumber),  2000,-1000.,+1000., 2000,-1000.,1000.);


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
    //if(Rd->nhits>0){
      Hist->nhits->Fill   (Rd->nhits);
      //}
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


void TCalStat::analyze(const art::Event& event) {

  _event = &event;
  mu2e::GeomHandle<mu2e::Tracker> th;
  _tracker = th.get();
  
  auto handle = event.getValidHandle<StrawGasStepCollection>(_trkfCollTag);
  const std::vector<mu2e::StrawGasStep>* mccol=handle.product(); 

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
  int found_0=0;
  int found_1=0;
  int found_2=0;
  int found_3=0;

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
	if(facev.at(i)==0) found_0=1;
	if(facev.at(i)==1) found_1=1;
      }
      if(planev.at(i)==1){
	if(panel==0) npan_10++;
	if(panel==1) npan_11++;
	if(panel==2) npan_12++;
	if(panel==3) npan_13++;
	if(panel==4) npan_14++;
	if(panel==5) npan_15++;
	if(facev.at(i)==0) found_2=1;
	if(facev.at(i)==1) found_3=1;
      }
    }
  }

  if(found_0 ==0 or found_1==0 or found_2==0 or found_3==0) printf("non ho filtrato");
  
 
  _event_data.nhits=_event_data.nhits+mccol->size();

  double dir_x[1500] ;
  double  mid_x[1500];
  double  dir_y[1500];
  double  mid_y[1500];
  int good_hits=0;
  uint16_t  PLANE[1500];
  uint16_t  PANEL[1500];
  double xt[1500];
  //double yt[1500];

  uint16_t STRAW[1500];
 
  for(size_t i=0;i<mccolsize;i++){

    dir_x[i]=0;
    dir_y[i]=0;
    mid_x[i]=0;
    mid_y[i]=0;
    PLANE[i]=0;
    PANEL[i]=0;
    xt[i]=0;
    //yt[i]=0;
    STRAW[i]=0;
  }
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
	  _Hist.plane[plane].panel[panel].xp_vs_yp_nh3->Fill(x,y);       
	}
	if(npan_15<=3 and npan_14<=3 and npan_13<=3 and npan_12<=3 and npan_11<=3 and npan_10<=3 and npan_05<=3 and npan_04<=3 and npan_03<=3 and npan_02<=3 and npan_01<=3 and npan_00<=3 and xp<300 and xp>-300 and yp>300 and yp<500 ){
	  _Hist.event.x_vs_y_nh3_ex->Fill(x,y);
	  _event_data.nhits_nh3_ex++;


	}
	if(npan_15<=3 and npan_14<=3 and npan_13<=3 and npan_12<=3 and npan_11<=3 and npan_10<=3 and npan_05<=3 and npan_04<=3 and npan_03<=3 and npan_02<=3 and npan_01<=3 and npan_00<=3){
	  _Hist.event.x_vs_y_nh3->Fill(x,y);    
	  dir_x[good_hits]=strw->direction().x();
	  dir_y[good_hits]= strw->direction().y();
	  mid_x[good_hits]=strw->getMidPoint().x();
	  mid_y[good_hits]=strw->getMidPoint().y();
	  PANEL[good_hits]=panel;
	  PLANE[good_hits]=plane;
	  xt[good_hits]=xp;
	  //yt[good_hits]=yp;
	  STRAW[good_hits]=straw;
	  good_hits++;

	}
    }
  }

  if(npan_15<=3 and npan_14<=3 and npan_13<=3 and npan_12<=3 and npan_11<=3 and npan_10<=3 and npan_05<=3 and npan_04<=3 and npan_03<=3 and npan_02<=3 and npan_01<=3 and npan_00<=3  ){
    _event_data.nhits_nh3+=mccolsize;

  double xc[1500];
  double yc[1500];
  double midm_x[1500];
  double midm_y[1500];
  double dirm_x[1500];
  double dirm_y[1500];
  //double dev_y[1500];
  //double dev_x[1500];
  for(int i=0;i<good_hits;i++){
    xc[i]=0;
    yc[i]=0;
    midm_x[i]=0;
    midm_y[i]=0;
    dirm_x[i]=0;
    dirm_y[i]=0;
  }
  //find the mean
  int tot=0;
  int k=0;
  while(k<good_hits){//1115222344
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

  //compute the points of intersection	ci sara un punto in meno perche prendo l intersezione e i punti sono disallineati perche ho ricostruito le intersezioni non le traccie
  //ora c e da capire come trasportare la traccia ricostruita sulla straw		  
  for(int i=0;i<(tot-1);i++ ){ //good_hits-1
    double m1=dirm_y[i]/dirm_x[i]; //dir_y[i]/dir_x[i];
    double q1=-m1*midm_x[i]+midm_y[i];
    //for(int j=i+1; j<good_hits;j++){
      //if(PANEL[i]!=PANEL[j] or PLANE[i]!=PLANE[j]){
	  double m2=dirm_y[i+1]/dirm_x[i+1];
	  double q2=-m2*midm_x[i+1]+midm_y[i+1];
	  if(m1!=m2){
	    xc[i]=(q1-q2)/(m2-m1);
	    yc[i]=m2*xc[i]+q2;
	    _Hist.event.x_vs_y_rec_p->Fill(xc[i],yc[i]);
	    // break;// perche se no poi trova il punto con ogni straw e diventa linea piatta
	  }
	  //}
	  //}
  }
  //fit funziona, devo implementare errori sugli assi correttamente attraverso la std dev
  double sx0=0;
  double sx1=0;
  double sx2=0;
  double sxy0=0;
  double sxy1=0;
  double sigmay=4;//da ripensarci e da implementare anche errori su x
  //double sigmax=40;
  for(int i=0;i<(tot-1);i++ ){
    sx0=sx0+(1/(sigmay*sigmay));
    sx1=sx1+(xc[i]/(sigmay*sigmay));
    sx2=sx2+((xc[i]*xc[i])/(sigmay*sigmay));
    sxy0=sxy0+(yc[i]/(sigmay*sigmay));
    sxy1=sxy1+((xc[i]*yc[i])/(sigmay*sigmay));


  }
  double D=sx0*sx2-sx1*sx1;
  double m=(sxy1*sx0-sxy0*sx1)/D;
  double q=(sxy0*sx2-sxy1*sx1)/D;
  //la m e la q sono trovate in maniera corretta. ora servirebbe capire come plottarle
  //std::cout<<"m"<<m<<"q"<<q<<std::endl;
  for(int i=0;i<good_hits;i++ ){
    double m1=dir_y[i]/dir_x[i];
    double q1=-m1*mid_x[i]+mid_y[i];//ora riutilizzo le variabili vere delle straw , quello che vorrei capire e' che faccio se una retta non passa per la straw?
    if(m1!=m){
      double xcf=(q1-q)/(m-m1);
      double ycf=m*xcf+q;
       _Hist.event.x_vs_y_rec_pf->Fill(xcf,ycf);
       _Hist.plane[PLANE[i]].panel[PANEL[i]].channel[STRAW[i]].x_bias->Fill((xcf*dir_x[i]+ycf*dir_y[i])-xt[i]);

    }
  }
  
  }
  if ( _analyzeStation == 1 and mccolsize>0)  fill_histograms();


  TModule::analyze(event);
}
}

 // end namespace mu2e

DEFINE_ART_MODULE(mu2e::TCalStat)
