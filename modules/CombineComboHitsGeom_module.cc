// combination of StrawHits in ComboHits for the calibration, without accounting for the timing and the energy to get stereo hits

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

  class CombineComboHitsGeom : public art::EDProducer {

    public:
      struct Config
      {
        using Name    = fhicl::Name;
        using Comment = fhicl::Comment;
        fhicl::Atom<art::InputTag>    CHC     { Name("ComboHitCollection"),    Comment("Input single-straw ComboHit collection") };
        fhicl::Atom<int>              minNComboHits    { Name("minNComboHits"),          Comment("Minimum number of combo hits")};
        fhicl::Atom<int>              maxNComboHits    { Name("maxNComboHits"),          Comment("Maximum number of combo hits")};
	fhicl::Atom<int>              maxNStrawHits    { Name("maxNStrawHits"),          Comment("Maximum number of straw hits")};

      };

      explicit CombineComboHitsGeom(const art::EDProducer::Table<Config>& config);
      void produce( art::Event& e);

    private:
    void combine( ComboHitCollection const& chcOrig, ComboHitCollection& chcol);
      void combineHits(ComboHit hit1, ComboHit hit2, ComboHit& combohit);

   
      art::ProductToken<ComboHitCollection> const _chctoken;
      int _maxNComboHits;
      int _minNComboHits;
      int _maxNStrawHits;
  };

  CombineComboHitsGeom::CombineComboHitsGeom(const art::EDProducer::Table<Config>& config) :
    EDProducer{config},
    _chctoken{consumes<ComboHitCollection>(config().CHC())},
    _maxNComboHits(config().maxNComboHits()),
    _minNComboHits(config().minNComboHits()),
    _maxNStrawHits(config().maxNStrawHits())
    {
      produces<ComboHitCollection>();
    }

  void CombineComboHitsGeom::produce(art::Event& event)
  {   
    //old ComboHit
    auto chcH = event.getValidHandle(_chctoken);
    const ComboHitCollection& chcOrig(*chcH);
   
    //new combohit
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


  void CombineComboHitsGeom::combine(ComboHitCollection const& chcOrig, ComboHitCollection& chcol)
  {    
    mu2e::GeomHandle<mu2e::Tracker> th;
    const Tracker* _tracker = th.get();
    
    std::vector<bool> isUsed(chcOrig.size(),false);
    for (size_t ich=0;ich<chcOrig.size();++ich) {

      int check=0;
      const ComboHit& hit1 = chcOrig[ich];
      ComboHit combohit;
      int panel1 = hit1.strawId().getPanel();

      for (size_t jch=ich+1;jch<chcOrig.size();++jch) {
        const ComboHit& hit2 = chcOrig[jch];
        if (hit2.strawId().getPanel() == panel1){
	  std::cout<<"ERROR: 2 combohits with same panel id"<<std::endl;
	}
	int panel2 = hit2.strawId().getPanel();
	int plane2 = hit2.strawId().getPlane();
	int plane1 = hit1.strawId().getPlane();
	const Plane* pl1 = &_tracker->getPlane(plane1);
	const Panel* pa1 = &pl1->getPanel(panel1);
	const Plane* pl2 = &_tracker->getPlane(plane2);
	const Panel* pa2 = &pl2->getPanel(panel2);
	if(plane1 == plane2 and pa1->getStraw(0).id().face()!=pa2->getStraw(0).id().face()){
	  combineHits(hit1, hit2, combohit);
	  check=1;
	  break;
	}

      }
      if(check==1){
	chcol.push_back(std::move(combohit));
      }
    }
  }


  void CombineComboHitsGeom::combineHits(ComboHit hit1, ComboHit hit2, ComboHit& combohit)
  {
    
    double diry1 = hit1._udir.y();
    double dirx1 = hit1._udir.x();
    double midx1 = hit1._pos.x();
    double midy1 = hit1._pos.y();
    double midz1 = hit1._pos.z();
    double diry2 = hit2._udir.y();
    double dirx2 = hit2._udir.x();
    double midx2 = hit2._pos.x();
    double midy2 = hit2._pos.y();
    double midz2 = hit2._pos.z();

    double m1=diry1/dirx1; 
    double q1=-m1*midx1+midy1;

    double m2=diry2/dirx2;
    double q2=-m2*midx2+midy2;

    double x = (q1-q2)/(m2-m1);
    double y = m2*x+q2;

    double z = (midz1+midz2)/2.;

    XYZVectorF midpos(x,y,z);
    XYZVectorF dir(0,0,0);
    combohit._nsh = 1;

    combohit._edep       = 0;
    combohit._etime[StrawEnd::cal] = 0;
    combohit._etime[StrawEnd::hv] = 0;

    combohit._dtime      = 0;
    combohit._ptime      = 0;
    combohit._time       = 0;
    combohit._timevar    = 0;
    combohit._wdist      = 0;
    combohit._edep       = 0; 
    combohit._qual       = 0;

    combohit._pos        = midpos; 
    combohit._udir       = dir; 
    combohit._uvar       = hit1.strawId().getPlane();
    combohit._vvar       = hit2.strawId().getPlane();; 
    combohit._wvar       = 0;//mm
  }
}

DEFINE_ART_MODULE(mu2e::CombineComboHitsGeom)
