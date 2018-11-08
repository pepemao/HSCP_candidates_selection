/*
 numberOfValidTrackerHits() = numberOfValidPixelHits()+numberOfValidStripHits()
*/
// system include files
#include <memory>
#include <fstream>
// user include files
//Tracks
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"

//Muons
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/MuonReco/interface/MuonMETCorrectionData.h"
#include "RecoMuon/MuonIdentification/plugins/MuonTimingProducer.h"
#include "RecoMuon/MuonIdentification/interface/TimeMeasurementSequence.h"
#include "DataFormats/MuonReco/interface/MuonTimeExtra.h"
#include "DataFormats/MuonReco/interface/MuonTimeExtraMap.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"


//Includes for DeDxInfoCollection
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"
#include "DataFormats/TrackReco/interface/TrackDeDxHits.h"
#include "DataFormats/TrackReco/interface/DeDxHit.h"
#include "RecoTracker/DeDx/interface/DeDxTools.h"
#include "TrackingTools/PatternTools/interface/TrajTrackAssociation.h"
#include "AnalysisDataFormats/SUSYBSMObjects/interface/HSCPDeDxInfo.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/ChainEvent.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "AnalysisDataFormats/SUSYBSMObjects/interface/HSCParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "AnalysisDataFormats/SUSYBSMObjects/interface/HSCPIsolation.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2D.h"
#include "DataFormats/SiStripDetId/interface/SiStripSubStructure.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/TrackReco/interface/DeDxHitInfo.h"

#include "Analysis_CommonFunction.h"
#include "Analysis_TOFUtility.h"
#include "TCanvas.h"
#include "TPad.h"


using namespace edm;
using namespace std;

//If use one input file
//TFile f("/afs/cern.ch/work/o/oshkola/MyEDMAnalyzer/CMSSW_8_0_24_patch1/src/Demo/DemoAnalyzer/plugins/Run2015_254232.root");
//TFile f("muonTeV_1.root");

//To use multiple input files
std::vector<std::string> fileNames;

std::string CommonPath = "root://se.cis.gov.pl//cms/store/user/fruboes/HSCP_20170802/Run2016_";
std::string Run2015 = "2015_";
std::string Run2016 = "2016_";
std::string EndFile = ".root";

std::string CommonPathMC = "/eos/cms/store/cmst3/user/querten/15_03_25_HSCP_Run2EDMFiles/";

// class intiialization
class DemoAnalyzer : public edm::EDAnalyzer {
public:
  explicit DemoAnalyzer(const edm::ParameterSet&);
  ~DemoAnalyzer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;


  // ----------member data ---------------------------

  //Initialize output tree
  TTree*  MyTree;
  edm::Service<TFileService> tfs;

  //Labels for usage of MuonTimeExtraMap <- syntax is taken from standard hscp code
  std::string        TOF_Label       = "combined";
  std::string        TOFdt_Label     = "dt";
  std::string        TOFcsc_Label    = "csc";
  /*
    std::string        dEdxS_Label     = "dedxASmi";
    std::string        dEdxS_Legend    = "I_{as}";
    std::string        dEdxM_Label     = "dedxHarm2";
    std::string        dEdxM_Legend    = "I_{h} (MeV/cm)";
  */
}; //end class initialization

//default constructor
DemoAnalyzer::DemoAnalyzer(const edm::ParameterSet& iConfig)
{}
//destructor with no parameters
DemoAnalyzer::~DemoAnalyzer()
{}

void
DemoAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  bool isMC = false;

  unsigned int nruns = 1;

//  std::string Runs2016[nruns] = {"279116","279479", "279588", "279653", "279654", "279656", "279658", "279667", "279681", "279682"};
  std::string Runs2016[nruns] = {"279975"};
  //std::string RunsMC[2] = { "PPStau_13TeV_M1599", "PPStau_13TeV_M494"};
  std::string RunMC_PPStau494 = "PPStau_13TeV_M1599";

  if (!isMC)
  {
    for (unsigned int i = 0; i < nruns; i++)
    {
      std::cout << Runs2016[i] << std::endl;
      fileNames.push_back(CommonPath + Runs2016[i] + EndFile);
    }
  }

  else
  { for (unsigned int i = 1; i < 2; i++)
    {
      fileNames.push_back(CommonPathMC + RunMC_PPStau494 + EndFile);
    }
  }
  // authorize Tree up to 2 Terabytes <- taken from standard hscp code
  TTree::SetMaxTreeSize(1000 * Long64_t(2000000000));
  MyTree = tfs->make<TTree> ("HscpTree", "HscpTree");
 
  const unsigned int MAX_HSCPS =  10000;
  const unsigned int MAX_TRACKS = 10000;
  const unsigned int MAX_MUONS =  10000;

  unsigned int    event_num;
  unsigned int    lumi_sec;
  unsigned int    run_num;
  unsigned int    NHSCPs;
  unsigned int    NTRACKs;
  unsigned int    NMUONs;

  int trkcoll_num_valid_tracker_hits            [MAX_TRACKS];
  int trkcoll_num_lost_tracker_hits             [MAX_TRACKS];
  int trkcoll_num_valid_pixel_hits              [MAX_TRACKS];
  int trkcoll_num_valid_strip_hits              [MAX_TRACKS];
  int trkcoll_num_valid_muon_hits               [MAX_TRACKS];
  int trkcoll_num_lost_muon_hits                [MAX_TRACKS];
  int trkcoll_num_valid_dt_hits                 [MAX_TRACKS];
  int trkcoll_num_valid_csc_hits                [MAX_TRACKS];
  int trkcoll_num_valid_rpc_hits                [MAX_TRACKS];
  int trkcoll_high_purity                       [MAX_TRACKS];
  double trkcoll_pt                             [MAX_TRACKS];
  double trkcoll_pt_err                         [MAX_TRACKS];
  double trkcoll_eta                            [MAX_TRACKS];
  double trkcoll_phi                            [MAX_TRACKS];
  double trkcoll_chi2_ndof                      [MAX_TRACKS];
  double trkcoll_valid_frac                     [MAX_TRACKS];
  double trkcoll_dz                             [MAX_TRACKS];
  double trkcoll_dxy                            [MAX_TRACKS];

  double          mucoll_pt                     [MAX_MUONS];
  double          mucoll_eta                    [MAX_MUONS];
  double          mucoll_phi                    [MAX_MUONS];

  double          trk_valid_frac                [MAX_HSCPS];
  double          hscp_p                        [MAX_HSCPS];
  double          hscp_pt_trkref                [MAX_HSCPS];
  double          hscp_pt_muref                 [MAX_HSCPS];
  double          hscp_pt_comb_mu               [MAX_HSCPS];
  double          hscp_pt_inner_trk             [MAX_HSCPS];
  double          hscp_pterr_trkref             [MAX_HSCPS];
  double          track_eta                     [MAX_HSCPS];
  double          muon_eta                      [MAX_HSCPS];
  double          track_phi                     [MAX_HSCPS];
  double          muon_phi                      [MAX_HSCPS];
  int             num_valid_tracker_hits        [MAX_HSCPS];
  int             num_valid_pixel_hits          [MAX_HSCPS];
  int             num_valid_strip_hits          [MAX_HSCPS];
  int             num_valid_mu_hits             [MAX_HSCPS];
  int             num_valid_dt_hits             [MAX_HSCPS];
  int             num_valid_csc_hits            [MAX_HSCPS];


  MyTree->Branch("event_num"                   , &event_num            , "event_num/I");
  MyTree->Branch("lumi_sec"                    , &lumi_sec             , "lumi_sec/I");
  MyTree->Branch("run_num"                     , &run_num              , "run_num/I");
  MyTree->Branch("NHSCPs"                      , &NHSCPs               , "NHSCPs/I");
  MyTree->Branch("NTRACKs"                     , &NTRACKs              , "NTRACKs/I");
  MyTree->Branch("NMUONs"                     , &NMUONs                , "NMUONs/I");

  MyTree->Branch("trkcoll_num_valid_tracker_hits"               , trkcoll_num_valid_tracker_hits           , "trkcoll_num_valid_tracker_hits[NTRACKs]/I");
  MyTree->Branch("trkcoll_num_lost_tracker_hits"                , trkcoll_num_lost_tracker_hits            , "trkcoll_num_lost_tracker_hits[NTRACKs]/I");
  MyTree->Branch("trkcoll_num_valid_pixel_hits"                 , trkcoll_num_valid_pixel_hits             , "trkcoll_num_valid_pixel_hits[NTRACKs]/I");
  MyTree->Branch("trkcoll_num_valid_strip_hits"                 , trkcoll_num_valid_strip_hits             , "trkcoll_num_valid_strip_hits[NTRACKs]/I");
  MyTree->Branch("trkcoll_num_valid_muon_hits"                  , trkcoll_num_valid_muon_hits              , "trkcoll_num_valid_muon_hits[NTRACKs]/I");
  MyTree->Branch("trkcoll_num_lost_muon_hits"                   , trkcoll_num_lost_muon_hits               , "trkcoll_num_lost_muon_hits[NTRACKs]/I");
  MyTree->Branch("trkcoll_num_valid_dt_hits"                    , trkcoll_num_valid_dt_hits                , "trkcoll_num_valid_dt_hits[NTRACKs]/I");
  MyTree->Branch("trkcoll_num_valid_csc_hits"                   , trkcoll_num_valid_csc_hits               , "trkcoll_num_valid_csc_hits[NTRACKs]/I");
  MyTree->Branch("trkcoll_num_valid_rpc_hits"                   , trkcoll_num_valid_rpc_hits               , "trkcoll_num_valid_rpc_hits[NTRACKs]/I");
  MyTree->Branch("trkcoll_high_purity"                          , trkcoll_high_purity                      , "trkcoll_high_purity[NTRACKs]/I");
  MyTree->Branch("trkcoll_pt"                  , trkcoll_pt            , "trkcoll_pt[NTRACKs]/D");
  MyTree->Branch("trkcoll_pt_err"              , trkcoll_pt_err        , "trkcoll_pt_err[NTRACKs]/D");
  MyTree->Branch("trkcoll_eta"                 , trkcoll_eta           , "trkcoll_eta[NTRACKs]/D");
  MyTree->Branch("trkcoll_phi"                 , trkcoll_phi           , "trkcoll_phi[NTRACKs]/D");
  MyTree->Branch("trkcoll_phi"                 , trkcoll_phi           , "trkcoll_phi[NTRACKs]/D");
  MyTree->Branch("trkcoll_chi2_ndof"           , trkcoll_chi2_ndof     , "trkcoll_chi2_ndof[NTRACKs]/D");
  MyTree->Branch("trkcoll_valid_frac"          , trkcoll_valid_frac    , "trkcoll_valid_frac[NTRACKs]/D");
  MyTree->Branch("trkcoll_dz"                  , trkcoll_dz            , "trkcoll_dz[NTRACKs]/D");
  MyTree->Branch("trkcoll_dxy"                 , trkcoll_dxy           , "trkcoll_dxy[NTRACKs]/D");

  MyTree->Branch("mucoll_pt"                   , mucoll_pt             , "mucoll_pt[NMUONs]/D");
  MyTree->Branch("mucoll_eta"                  , mucoll_eta            , "mucoll_eta[NMUONs]/D");
  MyTree->Branch("mucoll_phi"                  , mucoll_phi            , "mucoll_phi[NMUONs]/D");

  MyTree->Branch("trk_valid_frac"              , trk_valid_frac        , "trk_valid_frac[NHSCPs]/D");
  MyTree->Branch("hscp_p"                      , hscp_p                , "hscp_p[NHSCPs]/D");
  MyTree->Branch("hscp_pt_trkref"              , hscp_pt_trkref        , "hscp_pt_trkref[NHSCPs]/D");
  MyTree->Branch("hscp_pt_muref"               , hscp_pt_muref         , "hscp_pt_muref[NHSCPs]/D");
  MyTree->Branch("hscp_pt_comb_mu"             , hscp_pt_comb_mu       , "hscp_pt_comb_mu[NHSCPs]/D");
  MyTree->Branch("hscp_pt_inner_trk"           , hscp_pt_inner_trk     , "hscp_pt_inner_trk[NHSCPs]/D");
  MyTree->Branch("hscp_pterr_trkref"           , hscp_pterr_trkref     , "hscp_pterr_trkref[NHSCPs]/D");
  MyTree->Branch("track_eta"                   , track_eta             , "track_eta[NHSCPs]/D");
  MyTree->Branch("muon_eta"                    , muon_eta              , "muon_eta[NHSCPs]/D");
  MyTree->Branch("track_phi"                   , track_phi             , "track_phi[NHSCPs]/D");
  MyTree->Branch("muon_phi"                    , muon_phi              , "muon_phi[NHSCPs]/D");
  MyTree->Branch("num_valid_tracker_hits"                       , num_valid_tracker_hits                   , "num_valid_tracker_hits[NHSCPs]/I");
  MyTree->Branch("num_valid_pixel_hits"                         , num_valid_pixel_hits                     , "num_valid_pixel_hits[NHSCPs]/I");
  MyTree->Branch("num_valid_strip_hits"                         , num_valid_strip_hits                     , "num_valid_strip_hits[NHSCPs]/I");
  MyTree->Branch("num_valid_mu_hits"                            , num_valid_mu_hits                        , "num_valid_mu_hits[NHSCPs]/I");
  MyTree->Branch("num_valid_dt_hits"                            , num_valid_dt_hits                        , "num_valid_dt_hits[NHSCPs]/I");
  MyTree->Branch("num_valid_csc_hits"                           , num_valid_csc_hits                       , "num_valid_csc_hits[NHSCPs]/I");

  //Scale Factors for dE/dx measurements
  double dEdxSF [2] = {
    1.00000,   // [0]  unchanged
    1.21836    // [1]  Pixel data to SiStrip data
  };

  TH3F* dEdxTemplates = NULL;
  dedxGainCorrector trackerCorrector;

  //dEdx calibration
  if (!isMC) {
    trackerCorrector.LoadDeDxCalibration("/afs/cern.ch/work/o/oshkola/MyEDMAnalyzer/CMSSW_8_0_24_patch1/src/Demo/DemoAnalyzer/plugins/Data13TeVGains_v2.root");
  } else {
    trackerCorrector.TrackerGains = NULL; //FIXME check gain for MC
  }
  if (!isMC) {
    dEdxTemplates = loadDeDxTemplate("/afs/cern.ch/work/o/oshkola/MyEDMAnalyzer/CMSSW_8_0_24_patch1/src/Demo/DemoAnalyzer/plugins/Data13TeV_Deco_SiStripDeDxMip_3D_Rcd_v2_CCwCI.root", true);
  }
  else
  {
    dEdxTemplates = loadDeDxTemplate("/afs/cern.ch/work/o/oshkola/MyEDMAnalyzer/CMSSW_8_0_24_patch1/src/Demo/DemoAnalyzer/plugins/MC13TeV16_dEdxTemplate.root", true);
  }

  moduleGeom::loadGeometry("/afs/cern.ch/work/o/oshkola/MyEDMAnalyzer/CMSSW_8_0_24_patch1/src/Demo/DemoAnalyzer/plugins/CMS_GeomTree.root");
  muonTimingCalculator tofCalculator;
  tofCalculator.loadTimeOffset("/afs/cern.ch/work/o/oshkola/MyEDMAnalyzer/CMSSW_8_0_24_patch1/src/Demo/DemoAnalyzer/plugins/MuonTimeOffset.txt");

  unsigned int CurrentRun = 0;
 
//For one input file
//fwlite::Event ev(&f);

//For several input files (members of string vector fileNames)
  fwlite::ChainEvent ev(fileNames);
  
 std::cout << "ev size = " << ev.size() << std::endl;

//  ofstream file;
//  file.open("out_info.txt");

// for (Long64_t index = 0; index < ev.size(); index++)
  for (Long64_t index = 0; index < 10000; index++)
  {
    ev.to(index);

    //Event progress counter
    Int_t percent = index * (1.0 / ev.size()) * 100;
    std::cout << "\r" << std::string(percent / 5, '|') << percent << "%";
    std::cout.flush();

    if (CurrentRun != ev.eventAuxiliary().run())
    {
      CurrentRun = ev.eventAuxiliary().run();
      tofCalculator.setRun(CurrentRun);
      trackerCorrector.setRun(CurrentRun);
    }

    event_num = ev.eventAuxiliary().event();
    lumi_sec = ev.eventAuxiliary().luminosityBlock();
    run_num = ev.eventAuxiliary().run();

    //Collection of hscp candidates
    fwlite::Handle<susybsm::HSCParticleCollection> hscpCollH;
    hscpCollH.getByLabel(ev, "HSCParticleProducer");
    if (!hscpCollH.isValid()) continue;
    const susybsm::HSCParticleCollection& hscpColl = *hscpCollH;
    if (!hscpCollH.isValid()) continue;

//Map with information about combined 1/beta
    fwlite::Handle<reco::MuonTimeExtraMap> TOFCollH;
    TOFCollH.getByLabel(ev, "muons", TOF_Label.c_str());
    if (!TOFCollH.isValid()) {printf("Invalid TOF collection\n"); return;}

//Map with information about 1/beta from CSCs
    fwlite::Handle<reco::MuonTimeExtraMap> TOFcscCollH;
    TOFcscCollH.getByLabel(ev, "muons", TOFcsc_Label.c_str());
    if (!TOFcscCollH.isValid()) {printf("Invalid CSC  TOF collection\n"); return;}

//Map with information about 1/beta from DTs
    fwlite::Handle<reco::MuonTimeExtraMap> TOFdtCollH;
    TOFdtCollH.getByLabel(ev, "muons", TOFdt_Label.c_str());
    if (!TOFdtCollH.isValid()) {printf("Invalid DT  TOF collection\n"); return;}

//DeDx collection
    fwlite::Handle<reco::DeDxHitInfoAss> dedxCollH;
    dedxCollH.getByLabel(ev, "dedxHitInfo");
    if (!dedxCollH.isValid()) {printf("Invalid dedx collection\n");  continue;}

//Muon collection
    fwlite::Handle<reco::MuonCollection> MuonCollH;
    MuonCollH.getByLabel(ev, "muons");
    if (!MuonCollH.isValid()) continue;

//Track collecton
    fwlite::Handle<reco::TrackCollection> TrackColH;
    TrackColH.getByLabel(ev, "globalMuons");
    if (!TrackColH.isValid()) { cout << "Track collection is not valid\n"; continue;}
    const reco::TrackCollection& trkcoll = *TrackColH;

//Vertex Collection for dz and ?dxy
    fwlite::Handle< std::vector<reco::Vertex> > vertexCollHandle;
    vertexCollHandle.getByLabel(ev, "offlinePrimaryVertices");
    if (!vertexCollHandle.isValid()) {printf("Vertex Collection NotFound\n"); return;}
    const std::vector<reco::Vertex>& vertexColl = *vertexCollHandle;
    if (vertexColl.size() < 1) {printf("NO VERTEX\n"); return;}

    fwlite::Handle<CSCSegmentCollection> CSCSegmentCollHandle;
    fwlite::Handle<DTRecSegment4DCollection> DTSegmentCollHandle;

    CSCSegmentCollHandle.getByLabel(ev, "cscSegments");
    if (!CSCSegmentCollHandle.isValid()) {printf("CSC Segment Collection not found!\n"); continue;}

    DTSegmentCollHandle.getByLabel(ev, "dt4DSegments");
    if (!DTSegmentCollHandle.isValid()) {printf("DT Segment Collection not found!\n"); continue;}

    fwlite::Handle<susybsm::HSCPIsolationValueMap> IsolationH;
    IsolationH.getByLabel(ev, "HSCPIsolation", "R03"); //New format used for data since 17-07-2015
    if (!IsolationH.isValid())
    {
      IsolationH.getByLabel(ev, "HSCPIsolation03");//Old format used for first 2015B data, Signal and MC Backgrounds
      if (!IsolationH.isValid()) {printf("Invalid IsolationH\n"); continue;}
    }
    const ValueMap<susybsm::HSCPIsolation>& IsolationMap = *IsolationH.product();

    //need to access Track Collection because hitpattern.numberOfMuonHits gives 0 in hscp collection
    //and also i can't extract segment compatibility from hscp collection

    NTRACKs = 0;

    for (size_t i = 0; i < TrackColH->size();  i++)
    {
      const  reco::Track & tr = (*TrackColH)[i];

      if (tr.pt() < 55) continue;
      if (tr.eta() > 2.5) continue;
      if (tr.eta() < -2.5)  continue;

      const reco::HitPattern& hp = tr.hitPattern();

      trkcoll_num_valid_tracker_hits  [NTRACKs] = hp.numberOfValidTrackerHits();
      trkcoll_num_lost_tracker_hits   [NTRACKs] = hp.numberOfLostTrackerHits(reco::HitPattern::HitCategory::TRACK_HITS);
      trkcoll_num_valid_pixel_hits    [NTRACKs] = hp.numberOfValidPixelHits();
      trkcoll_num_valid_strip_hits    [NTRACKs] = hp.numberOfValidStripHits();
      trkcoll_num_valid_muon_hits     [NTRACKs] = hp.numberOfValidMuonHits();
      trkcoll_num_lost_muon_hits      [NTRACKs] = hp.numberOfLostMuonHits();
      trkcoll_num_valid_dt_hits       [NTRACKs] = hp.numberOfValidMuonDTHits();
      trkcoll_num_valid_csc_hits      [NTRACKs] = hp.numberOfValidMuonCSCHits();
      trkcoll_num_valid_rpc_hits      [NTRACKs] = hp.numberOfValidMuonRPCHits();
      trkcoll_high_purity             [NTRACKs] = tr.quality(reco::TrackBase::TrackQuality::highPurity) ? 1 : 0;
      trkcoll_pt                      [NTRACKs] = tr.pt();
      trkcoll_pt_err                  [NTRACKs] = tr.ptError() / tr.pt();
      trkcoll_eta                     [NTRACKs] = tr.eta();
      trkcoll_phi                     [NTRACKs] = tr.phi();
      trkcoll_chi2_ndof               [NTRACKs] = tr.chi2() / tr.ndof();
      trkcoll_valid_frac              [NTRACKs] = tr.validFraction();
      trkcoll_dz                      [NTRACKs] = tr.dz(vertexColl[i].position());
      trkcoll_dxy                     [NTRACKs] = tr.dxy(vertexColl[i].position());

      NTRACKs++;
     }

   
    NMUONs = 0;
   
    vector<double> muon_pt_fromcoll;
    vector<double> segment_compatibility;

    for (size_t i = 0; i < MuonCollH->size();  i++)
    {
      if (!MuonCollH.isValid()) {cout << "muon coll is not valid" << endl; continue;}

      const  reco::Muon & p = (*MuonCollH)[i];
      if (p.pt() < 55) continue;
      if (p.eta() > 2.5) continue;
      if (p.eta() < -2.5)  continue;

      mucoll_pt   [NMUONs] = p.pt();
      mucoll_eta  [NMUONs] = p.eta();
      mucoll_phi  [NMUONs] = p.phi();

      muon_pt_fromcoll.push_back(p.pt());
      segment_compatibility.push_back(muon::segmentCompatibility(p));
      //  cout << "dxy = " << fabs(p.muonBestTrack()->dxy(vertexColl[i].position())) << endl;
      //  cout << "dz = " << fabs(p.muonBestTrack()->dz(vertexColl[i].position())) << endl;
      //  cout << "dxy1 = " << fabs(p.bestTrack()->dxy(vertexColl[i].position())) << endl;
      //  cout << "dz1 = " << fabs(p.bestTrack()->dz(vertexColl[i].position())) << endl;

      NMUONs++;
    }

   
    NHSCPs = 0;

    vector<int> qualind;

    for (unsigned int c = 0; c < hscpColl.size(); c++) {

      susybsm::HSCParticle hscp = hscpColl[c];
      reco::MuonRef muonr  = hscp.muonRef();
      reco::TrackRef track = hscp.trackRef();

      //skip hscp candidates that do not have track
      if (track.isNull()) continue;

      //skip hscp candidates that do not have muon reference
      if (muonr.isNull()) continue;

      if (track->hitPattern().numberOfValidPixelHits() <= 1)         continue;
      if (track->pt() < 55)    continue;
      if  (track->eta() > 2.1) continue;
      if  (track->eta() < -2.1)  continue;
      //if (track->ptError() / track->pt() > 0.25) continue;
      const reco::HitPattern& hp = track->hitPattern();

      //cout << "num of pixel hits = " << hp.numberOfValidPixelHits() << endl;
      // cout << "num of strip hits = " << hp.numberOfValidStripHits() << endl;

      hscp_pt_trkref   [NHSCPs]  = track->pt();
      hscp_pt_muref    [NHSCPs]  = muonr->pt();
      if (muonr->combinedMuon().isNonnull())
      {
        hscp_pt_comb_mu  [NHSCPs]  = muonr->combinedMuon()->pt();
      }
      else hscp_pt_comb_mu  [NHSCPs] = -1;
      if (muonr->innerTrack().isNonnull())
      {
        hscp_pt_inner_trk [NHSCPs] = muonr->innerTrack()->pt();
      }
      else hscp_pt_inner_trk [NHSCPs] = -1;
      hscp_pterr_trkref [NHSCPs] = track->ptError() / track->pt();
      trk_valid_frac   [NHSCPs] = track->validFraction();
      track_eta   [NHSCPs] = track->eta();
      muon_eta    [NHSCPs] = muonr->eta();
      track_phi   [NHSCPs] = track->phi();
      muon_phi    [NHSCPs] = muonr->phi();
      num_valid_pixel_hits [NHSCPs] = hp.numberOfValidPixelHits();
      num_valid_strip_hits [NHSCPs] = hp.numberOfValidStripHits();
      num_valid_tracker_hits [NHSCPs] = hp.numberOfValidTrackerHits();
      num_valid_mu_hits [NHSCPs] = hp.numberOfValidMuonHits();
      num_valid_dt_hits [NHSCPs] = hp.numberOfValidMuonDTHits();
      num_valid_csc_hits [NHSCPs] = hp.numberOfValidMuonCSCHits();

      string muonref, trackref, isstandalone, isglobal, iscombined;
      muonr.isNull() ? muonref = "no" : muonref = "yes";
      track.isNull() ? trackref = "no" : trackref = "yes";
      (muonr.isNonnull() && muonr->combinedMuon().isNonnull()) ? iscombined = "yes" : iscombined = "no";
      (muonr.isNonnull() && muonr->standAloneMuon().isNonnull()) ? isstandalone = "yes" : isstandalone = "no";
      muonr->isGlobalMuon() ? isglobal = "yes" : isglobal = "no";

      // double muon_segment_compatibility = -110;
      // for (unsigned int i = 0; i < muon_pt_fromcoll.size(); i++)
      // {
      //  if (muonr->pt() == muon_pt_fromcoll[i]) { muon_segment_compatibility = segment_compatibility[i]; }

      // }
      // if (muon_segment_compatibility == -100) cout << "SOMETHING IS GOING WRONG WITH SEGMENT COMPATIBILITY \b";

      //define Loose muon
      bool isLoose = muonr->isPFMuon() && (muonr->isGlobalMuon() || muonr->isTrackerMuon());

      //define Medium muon
      bool goodGlobal = muonr->isGlobalMuon() && muonr->globalTrack().isNonnull() && muonr->globalTrack()->normalizedChi2() < 3 &&
                        muonr->combinedQuality().chi2LocalPosition < 12 &&
                        muonr->combinedQuality().trkKink < 20 ;

      bool isMedium = isLoose && muonr->innerTrack()->validFraction() > 0.8;
      //               muon_segment_compatibility > (goodGlobal ? 0.303 : 0.451);

      //define Tight muon
      bool isTight = muonr->isGlobalMuon() &&  muonr->isPFMuon() && muonr->globalTrack().isNonnull() && muonr->globalTrack()->normalizedChi2() < 10. &&
                     muonr->globalTrack()->hitPattern().numberOfValidMuonHits() > 0 && muonr->numberOfMatchedStations() > 1 &&
                     muonr->muonBestTrack().isNonnull() &&
                     fabs(track->dxy(vertexColl[c].position())) < 0.2 &&
                     fabs(track->dz(vertexColl[c].position())) < 0.5 &&
                     // fabs(muonr->muonBestTrack()->dxy(vertexColl[c].position())) < 0.2 &&
                     // fabs(muonr->muonBestTrack()->dz(vertexColl[c].position())) < 0.5 &&

                     //    fabs(muonr->bestTrack()->dxy(vertexColl[c].position())) < 0.2 &&
                     //  fabs(muonr->bestTrack()->dz(vertexColl[c].position())) < 0.5 &&
                     muonr->innerTrack().isNonnull() &&
                     muonr->innerTrack()->hitPattern().numberOfValidPixelHits() > 0 &&
                     muonr->innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5;

      string loose, medium, tight;
      isLoose ? loose = "yes" : loose = "no";
      isMedium ? medium = "yes" : medium = "no";
      isTight ? tight = "yes" : tight = "no";

      /*
      file << ".................Location................." << endl;
      file << " run number = " << run_num << " luminosity block =  " << lumi_sec << " event =  " << event_num  << endl;
      file << "................Kinematics..............." << endl;
      file << " track_pt = " << track->pt()  << " muon_pt = " << muonr->pt()  << " hscp_cand_pt =  " << hscpColl[c].pt() << endl;
      file << " track_eta = " << track->eta() << " muon_eta = " << muonr->eta() << endl;
      file << " track_phi = " << track->phi()  << " muon_phi = " << muonr->phi() << endl;
      for (unsigned int i = 0; i < muon_pt_fromcoll.size(); i++)
      {
        file << "muon_pt from MUON collection = " << muon_pt_fromcoll[i] << endl;
      }
      for (unsigned int i = 0; i < track_pt_fromcoll.size(); i++)
      {
        //if( track->pt() == track_pt_fromcoll[i])
        // {
        file << "track_pt from TRACK collection = " << track_pt_fromcoll[i] << endl;
        file << "track eta from traack collection = " << trk_eta[i] << endl;
        file << "N of valid STRIP hits = " << nstrip[i] << endl;
        file << "N of valid PIXEL hits = " << npixel[i] << endl;
        file << "N of valid DT hits = " << ndt[i] << endl;
        file << "N of valid CSC hits = " << ncsc[i] << endl;
        // }
      }
      file << ".....Check hscp cand-e has ref......" << endl;
      file << " has muon reference? " << muonref << endl;
      file << " has track reference? " << trackref << endl;

      file << "..............Track quality.............." << endl;
      if (track.isNonnull())
      {
        file << " error on track pt = " << track->ptError() / track->pt() << endl;
        file << " track charge = " << track->charge() << endl;
        file << " number of valid track hits = " << track->numberOfValidHits() << endl;
        file << " number of lost track hits = " << track->numberOfLostHits() << endl;
        file << " fraction of valid track hits = " << track->validFraction() << endl;
        file << " chi2/ndof = " << track->normalizedChi2() << endl;
      }
      file << "..............Muon quality.............." << endl;
      file << " is global? " << isglobal << endl;
      file << " is standalone? " << isstandalone << endl;
      file << " is combined? " << iscombined << endl;
      file << " is loose? " << loose << endl;
      file << " is medium? " << medium << endl;
      file << " is tight? " << tight << endl;

      file << ".....Values taken as hscpColl[c].pt()......" << endl;
      if (muonr->combinedMuon().isNonnull())
      {
        file << " combinedMuon()->pt() = " << muonr->combinedMuon()->pt() << endl;
      }
      else file << " Combined muon is null" << endl;

      if (muonr->innerTrack().isNonnull())
      {
        file << " innerTrack->pt() = " << muonr->innerTrack()->pt() << endl;
      }
      else file << "Inner track is null" << endl;
      if (muonr->standAloneMuon().isNonnull())
      {
        file << " standAloneMuon()->pt() = " << muonr->standAloneMuon()->pt() << endl;
      } else file << "standalone muon is null" << endl;
      if (track.isNonnull())
      {
        file << " track->pt() = " << track->pt() << endl;
      } else file << "track is null" << endl;
      file << endl << endl;
      */

      //Ias and Ih DeDx estimators
      const reco::DeDxHitInfo* dedxHits = NULL;
      reco::DeDxHitInfoRef dedxHitsRef = dedxCollH->get(track.key());
      if (!dedxHitsRef.isNull())dedxHits = &(*dedxHitsRef);

      bool  useClusterCleaning = true;

      // Call the function  to compute Ias and Iah
      reco::DeDxData dedxSObjTmp  = computedEdx(dedxHits, dEdxSF, dEdxTemplates, true, useClusterCleaning, false, false, trackerCorrector.TrackerGains, true, true, 99, false, 1, 0.00, NULL);
      reco::DeDxData dedxMObjTmp = computedEdx(dedxHits, dEdxSF, NULL,          true, useClusterCleaning, false      , false, trackerCorrector.TrackerGains, true, true, 99, false, 1, 0.15, NULL);

      reco::DeDxData* dedxSObj  = dedxSObjTmp.numberOfMeasurements() > 0 ? &dedxSObjTmp : NULL;
      reco::DeDxData* dedxMObj  = dedxMObjTmp.numberOfMeasurements() > 0 ? &dedxMObjTmp : NULL;

      //I/B calculations
      const reco::MuonTimeExtra* tof = NULL;
      const reco::MuonTimeExtra* dttof = NULL;
      const reco::MuonTimeExtra* csctof = NULL;

      const reco::MuonTimeExtra* tof_edm = NULL;
      const reco::MuonTimeExtra* dttof_edm = NULL;
      const reco::MuonTimeExtra* csctof_edm = NULL;

      const CSCSegmentCollection& CSCSegmentColl = *CSCSegmentCollHandle;
      const DTRecSegment4DCollection& DTSegmentColl = *DTSegmentCollHandle;

      //Call the functon from TOF_Utility to calculate 1/B
      //  tofCalculator.computeTOF(muonr, CSCSegmentColl, DTSegmentColl, true );
      //tof_edm  = &tofCalculator.combinedTOF;
      // dttof = &tofCalculator.dtTOF;
      // csctof = &tofCalculator.cscTOF;

      //get standard 1/beta (from EDM)
      //tof_edm  = &TOFCollH->get(hscp.muonRef().key());
      // dttof_edm = &TOFdtCollH->get(hscp.muonRef().key());
      // csctof_edm = &TOFcscCollH->get(hscp.muonRef().key());

      //   const reco::HitPattern& hp = track->hitPattern();


      if (isLoose && isMedium && isTight) qualind.push_back(111);
      else if (isLoose && isMedium && !isTight) qualind.push_back(110);
      else if (isLoose && !isMedium && !isTight) qualind.push_back(100);
      else if (!isLoose && !isMedium && !isTight) qualind.push_back(-1);
      else if (isLoose && !isMedium && isTight) qualind.push_back(101);
      else {qualind.push_back(-1);}

      NHSCPs++;

    } //end of loop over hscp candidates;

    MyTree->Fill();
  } //end of the event loop
  //file.close();

#ifdef THIS_IS_AN_EVENT_EXAMPLE
  Handle<ExampleData> pIn;
  iEvent.getByLabel("example", pIn);
#endif

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  ESHandle<SetupData> pSetup;
  iSetup.get<SetupRecord>().get(pSetup);
#endif
}
// ------------ method called once each job just before starting event loop  ------------
void
DemoAnalyzer::beginJob()
{
}
// ------------ method called once each job just after ending the event loop  ------------
void
DemoAnalyzer::endJob()
{
}
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DemoAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(DemoAnalyzer);

