// combination of StrawHits in ComboHits for the calibration, without accounting for the timing and the energy

#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/VirtualDetector.hh"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "TMath.h"

#include <iostream>

namespace mu2e {

  class CombineStrawHitsGeom : public art::EDProducer {

    public:
      struct Config
      {
        using Name    = fhicl::Name;
        using Comment = fhicl::Comment;
        fhicl::Atom<art::InputTag>    CHC     { Name("ComboHitCollection"),    Comment("Input single-straw ComboHit collection") };
	fhicl::Atom<bool>             filter    { Name("filter"),          Comment("filter")};
	fhicl::Atom<bool>             filter2    { Name("filter2"),          Comment("filter2")};
     	fhicl::Atom<bool>             filter3    { Name("filter3"),          Comment("filter3")};

      };

      explicit CombineStrawHitsGeom(const art::EDProducer::Table<Config>& config);
      void produce( art::Event& e);

    private:
    void combine( ComboHitCollection const& chcOrig, ComboHitCollection& chcol);
      void combineHits(const ComboHitCollection& chcOrig, ComboHit& combohit);

   
    art::ProductToken<ComboHitCollection> const _chctoken;
    bool _filter;
    bool _filter2;
    bool _filter3;
  };

  CombineStrawHitsGeom::CombineStrawHitsGeom(const art::EDProducer::Table<Config>& config) :
    EDProducer{config},
    _chctoken{consumes<ComboHitCollection>(config().CHC())},
    _filter(config().filter()),
    _filter2(config().filter2()),
    _filter3(config().filter3())
    {
      produces<ComboHitCollection>();
    }

  void CombineStrawHitsGeom::produce(art::Event& event){   
    //new ComboHit
    auto chcH = event.getValidHandle(_chctoken);
    const ComboHitCollection& chcOrig(*chcH);

    auto chcolNew = std::make_unique<ComboHitCollection>();
    chcolNew->reserve(chcOrig.size());
    chcolNew->setParent(chcH);

   
    ComboHitCollection chcolsort;
    chcolsort.reserve(chcOrig.size());
    chcolsort.setParent(chcOrig.parent());
    std::array<std::vector<uint16_t>,StrawId::_nupanels> panels;
    size_t nsh = chcOrig.size();
    for(uint16_t ish=0;ish<nsh;++ish){
      ComboHit const& ch = chcOrig[ish];
      panels[ch.strawId().uniquePanel()].push_back(ish);
    }
    for (uint16_t ipanel=0;ipanel<StrawId::_nupanels;ipanel++){
      for (uint16_t ish=0;ish<panels[ipanel].size();ish++){
	chcolsort.push_back(chcOrig.at(panels[ipanel][ish]));
      }
    }
    combine(chcolsort, *chcolNew);
    event.put(std::move(chcolNew));
  }


  void CombineStrawHitsGeom::combine( ComboHitCollection const& chcOrig, ComboHitCollection& chcol)  {
    std::vector<bool> isUsed(chcOrig.size(),false);
    for (size_t ich=0;ich<chcOrig.size();++ich) {
      if (isUsed[ich]) continue;
      isUsed[ich] = true;

      const ComboHit& hit1 = chcOrig[ich];
      ComboHit combohit;
      combohit.init(hit1,ich);
      int panel1 = hit1.strawId().uniquePanel();

      for (size_t jch=ich+1;jch<chcOrig.size();++jch) {
        if (isUsed[jch]) continue;
        const ComboHit& hit2 = chcOrig[jch];

        if (hit2.strawId().uniquePanel() != panel1) break;

        bool ok = combohit.addIndex(jch);
        if (!ok){
          std::cout << "CombineStrawHitsGeom past limit" << std::endl;
        }else {
          isUsed[jch]= true;
	}
      }
      combineHits(chcOrig, combohit);

      chcol.push_back(std::move(combohit));
    }
  }


  void CombineStrawHitsGeom::combineHits(const ComboHitCollection& chcOrig, ComboHit& combohit)
  {
    XYZVectorF midpos;
    XYZVectorF dir;
    double errx(0), erry(0);
    double eacc(0),ctacc(0),dtacc(0),twtsum(0),ptacc(0),wacc(0),wacc2(0),wwtsum(0);
    double etacc[2]={0,0};
    combohit._nsh = 0;
    mu2e::GeomHandle<mu2e::Tracker> th;
    const Tracker* _tracker = th.get();
    //nCombo=nStrawHits
    for (size_t ich = 0; ich < combohit.nCombo(); ++ich){
      size_t index = combohit.index(ich);
      if (index > chcOrig.size()){
	std::cout<<"mu2e::CombineStrawHitsGeom: inconsistent index "<<std::endl;
      }
      const ComboHit& ch = chcOrig[index];
      combohit._nsh += ch.nStrawHits();

      uint16_t panel = ch.strawId().getPanel();
      uint16_t plane = ch.strawId().getPlane();
      uint16_t straw = ch.strawId().getStraw();
      const Plane* pln = &_tracker->getPlane(plane);
      const Panel* pnl = &pln->getPanel(panel);
      const Straw* strw = &pnl->getStraw(straw);
      midpos += (XYZVectorF) strw->getMidPoint();//ch.centerPos();
      dir += (XYZVectorF) strw->getDirection();
      eacc += ch.energyDep();
      dtacc += ch.driftTime();
      ptacc += ch.propTime();
      etacc[StrawEnd::cal] += ch.endTime(StrawEnd::cal);
      etacc[StrawEnd::hv] += ch.endTime(StrawEnd::hv);
      double twt = 1.0/(ch.timeVar());
      twtsum += twt;
      ctacc += ch.correctedTime()*twt;
      // weight position values by U (wire direction) error
      double wwt = 1.0/(ch.wireVar());
      wwtsum += wwt;
      wacc  += ch.wireDist()*wwt;
      wacc2 += ch.wireDist()*ch.wireDist()*wwt;
      if(errx<strw->halfLength()){
	errx=strw->halfLength();
      }
    }
   
    if(combohit.nStrawHits() != combohit.nCombo()){
      std::cout<<"mu2e::CombineStrawHitsGeom: inconsistent count "<< std::endl;
    }
    
    double nsh = static_cast<double>(combohit.nStrawHits());
    combohit._edep       = 0; // simple averages
    combohit._etime[StrawEnd::cal] = etacc[StrawEnd::cal]/nsh;
    combohit._etime[StrawEnd::hv] = etacc[StrawEnd::hv]/nsh;
    midpos /= nsh;
    dir    /= nsh;
    errx   *= 2;
    for (size_t ich = 0; ich < combohit.nCombo(); ++ich){
      size_t index = combohit.index(ich);
      const ComboHit& ch = chcOrig[index];
      uint16_t panel = ch.strawId().getPanel();
      uint16_t plane = ch.strawId().getPlane();
      uint16_t straw = ch.strawId().getStraw();
      const Plane* pln = &_tracker->getPlane(plane);
      const Panel* pnl = &pln->getPanel(panel);
      const Straw* strw = &pnl->getStraw(straw);
      erry+=TMath::Power((double) strw->getMidPoint().y()- (double) midpos.y(),2);
    }
    erry=TMath::Sqrt(erry);
    erry /= nsh;
    double costheta = dir.x();
    double sintheta = dir.y();
    double xerrp=TMath::Abs(errx*costheta+erry*sintheta);
    double yerrp=TMath::Abs(-erry*costheta+errx*sintheta);
    combohit._dtime      = dtacc/nsh;
    combohit._ptime      = ptacc/nsh;
    combohit._time       = ctacc/twtsum;
    combohit._timevar    = 1.0/twtsum;
    combohit._wdist      = wacc/wwtsum;
    combohit._edep       =  eacc/nsh; 
    combohit._qual       = 0;

    combohit._pos        = midpos; 
    combohit._udir       = dir; 
    combohit._uvar       = xerrp;
    combohit._vvar       = yerrp; 
    combohit._wvar       = 5.;//mm
  }
}

DEFINE_ART_MODULE(mu2e::CombineStrawHitsGeom)
