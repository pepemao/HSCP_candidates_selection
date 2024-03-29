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


#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"


#include "Analysis_CommonFunction.h"
#include "Analysis_TOFUtility.h"
#include "TCanvas.h"
#include "TPad.h"

using namespace edm;
using namespace std;

//If use one input file
TFile f("/afs/cern.ch/work/o/oshkola/MYDEMOANALYZER_new/CMSSW_9_4_6_patch1/src/Demo/DemoAnalyzer/GMStau_13TeV_M494.root");
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
  std::string        Trig_label      = "RECO";
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


  bool isMC = true;

  unsigned int nruns = 1;

  std::string Runs2016[nruns] = {"279694"};
  //std::string Runs2016[nruns] = {"279694", "279766", "279931", "279975", "280242", "281707", "282735", "282919", "283059", "283270", "283283", "283478", "283548", "283865", "284037" };

  //std::string Runs2016[nruns] = {"280242" };//std::string RunsMC[2] = { "PPStau_13TeV_M1599", "PPStau_13TeV_M494"};
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

  // authorize the output tree up to 2 Terabytes <- taken from standard hscp code
  TTree::SetMaxTreeSize(1000 * Long64_t(2000000000));
  MyTree = tfs->make<TTree> ("HscpTree", "HscpTree");

  const unsigned int MAX_HSCPS =  10000;
  const unsigned int MAX_TRACKS = 10000;
  const unsigned int MAX_MUONS =  10000;

  unsigned int    event_num;
  unsigned int    lumi_sec;
  unsigned int    run_num;
  int             label;
  unsigned int    NHSCPs;
  unsigned int    NTRACKs;
  unsigned int    NMUONs;


  int    trkcoll_num_valid_tracker_hits  [MAX_TRACKS];
  int    trkcoll_num_lost_tracker_hits   [MAX_TRACKS];
  int    trkcoll_num_valid_pixel_hits    [MAX_TRACKS];
  int    trkcoll_num_valid_strip_hits    [MAX_TRACKS];
  int    trkcoll_num_valid_muon_hits     [MAX_TRACKS];
  int    trkcoll_num_lost_muon_hits      [MAX_TRACKS];
  int    trkcoll_num_valid_dt_hits       [MAX_TRACKS];
  int    trkcoll_num_valid_csc_hits      [MAX_TRACKS];
  int    trkcoll_num_valid_rpc_hits      [MAX_TRACKS];
  int    trkcoll_high_purity             [MAX_TRACKS];
  double trkcoll_pt                      [MAX_TRACKS];
  double trkcoll_pt_err                  [MAX_TRACKS];
  double trkcoll_eta                     [MAX_TRACKS];
  double trkcoll_phi                     [MAX_TRACKS];
  double trkcoll_chi2_ndof               [MAX_TRACKS];
  double trkcoll_valid_frac              [MAX_TRACKS];
  double trkcoll_dz                      [MAX_TRACKS];
  double trkcoll_dxy                     [MAX_TRACKS];
  int    trkcoll_has_hscp                [MAX_TRACKS];

  double mucoll_pt                       [MAX_MUONS];
  double mucoll_eta                      [MAX_MUONS];
  double mucoll_phi                      [MAX_MUONS];


  double hscp_trk_valid_frac             [MAX_HSCPS];
  double hscp_p                          [MAX_HSCPS];
  double hscp_pt_trkref                  [MAX_HSCPS];
  double hscp_pt_muref                   [MAX_HSCPS];
  double hscp_pt_comb_mu                 [MAX_HSCPS];
  double hscp_pt_inner_trk               [MAX_HSCPS];
  double hscp_pterr_trkref               [MAX_HSCPS];
  double hscp_track_eta                  [MAX_HSCPS];
  double hscp_muon_eta                   [MAX_HSCPS];
  double hscp_track_phi                  [MAX_HSCPS];
  double hscp_muon_phi                   [MAX_HSCPS];
  int    hscp_num_valid_tracker_hits     [MAX_HSCPS];
  int    hscp_num_valid_pixel_hits       [MAX_HSCPS];
  int    hscp_num_valid_strip_hits       [MAX_HSCPS];
  int    hscp_num_valid_mu_hits          [MAX_HSCPS];
  int    hscp_num_valid_dt_hits          [MAX_HSCPS];
  int    hscp_num_valid_csc_hits         [MAX_HSCPS];
  double hscp_I_as                       [MAX_HSCPS];
  double hscp_I_h                        [MAX_HSCPS];
  double hscp_ibeta                      [MAX_HSCPS];
  double hscp_ibeta_err                  [MAX_HSCPS];
  double hscp_calo_e_over_trk_p          [MAX_HSCPS];
  double hscp_chi2_ndof                  [MAX_HSCPS];
  double hscp_dxy                        [MAX_HSCPS];
  double hscp_dz                         [MAX_HSCPS];
  double hscp_dz_obs                     [MAX_HSCPS];
  double hscp_dxy_obs                    [MAX_HSCPS];
  double hscp_trk_sum_pt                 [MAX_HSCPS];
  int    hscp_qual_ind                   [MAX_HSCPS];
  double hscp_num_tof_meas               [MAX_HSCPS];
  double hscp_num_dedx_meas              [MAX_HSCPS];

  MyTree->Branch("event_num"                            , &event_num            , "event_num/I");
  MyTree->Branch("lumi_sec"                             , &lumi_sec             , "lumi_sec/I");
  MyTree->Branch("run_num"                              , &run_num              , "run_num/I");
  MyTree->Branch("label"                                , &label                , "label/I");
  MyTree->Branch("NHSCPs"                               , &NHSCPs               , "NHSCPs/I");
  MyTree->Branch("NTRACKs"                              , &NTRACKs              , "NTRACKs/I");
  MyTree->Branch("NMUONs"                               , &NMUONs               , "NMUONs/I");

  MyTree->Branch("trkcoll_num_valid_tracker_hits"     , trkcoll_num_valid_tracker_hits     , "trkcoll_num_valid_tracker_hits[NTRACKs]/I");
  MyTree->Branch("trkcoll_num_lost_tracker_hits"      , trkcoll_num_lost_tracker_hits      , "trkcoll_num_lost_tracker_hits[NTRACKs]/I");
  MyTree->Branch("trkcoll_num_valid_pixel_hits"       , trkcoll_num_valid_pixel_hits       , "trkcoll_num_valid_pixel_hits[NTRACKs]/I");
  MyTree->Branch("trkcoll_num_valid_strip_hits"       , trkcoll_num_valid_strip_hits       , "trkcoll_num_valid_strip_hits[NTRACKs]/I");
  MyTree->Branch("trkcoll_num_valid_muon_hits"        , trkcoll_num_valid_muon_hits        , "trkcoll_num_valid_muon_hits[NTRACKs]/I");
  MyTree->Branch("trkcoll_num_lost_muon_hits"         , trkcoll_num_lost_muon_hits         , "trkcoll_num_lost_muon_hits[NTRACKs]/I");
  MyTree->Branch("trkcoll_num_valid_dt_hits"          , trkcoll_num_valid_dt_hits          , "trkcoll_num_valid_dt_hits[NTRACKs]/I");
  MyTree->Branch("trkcoll_num_valid_csc_hits"         , trkcoll_num_valid_csc_hits         , "trkcoll_num_valid_csc_hits[NTRACKs]/I");
  MyTree->Branch("trkcoll_num_valid_rpc_hits"         , trkcoll_num_valid_rpc_hits         , "trkcoll_num_valid_rpc_hits[NTRACKs]/I");
  MyTree->Branch("trkcoll_high_purity"                , trkcoll_high_purity                , "trkcoll_high_purity[NTRACKs]/I");
  MyTree->Branch("trkcoll_pt"                         , trkcoll_pt                         , "trkcoll_pt[NTRACKs]/D");
  MyTree->Branch("trkcoll_pt_err"                     , trkcoll_pt_err                     , "trkcoll_pt_err[NTRACKs]/D");
  MyTree->Branch("trkcoll_eta"                        , trkcoll_eta                        , "trkcoll_eta[NTRACKs]/D");
  MyTree->Branch("trkcoll_phi"                        , trkcoll_phi                        , "trkcoll_phi[NTRACKs]/D");
  MyTree->Branch("trkcoll_phi"                        , trkcoll_phi                        , "trkcoll_phi[NTRACKs]/D");
  MyTree->Branch("trkcoll_chi2_ndof"                  , trkcoll_chi2_ndof                  , "trkcoll_chi2_ndof[NTRACKs]/D");
  MyTree->Branch("trkcoll_valid_frac"                 , trkcoll_valid_frac                 , "trkcoll_valid_frac[NTRACKs]/D");
  MyTree->Branch("trkcoll_dz"                         , trkcoll_dz                         , "trkcoll_dz[NTRACKs]/D");
  MyTree->Branch("trkcoll_dxy"                        , trkcoll_dxy                        , "trkcoll_dxy[NTRACKs]/D");
  MyTree->Branch("trkcoll_has_hscp"                   , trkcoll_has_hscp                   , "trkcoll_has_hscp[NTRACKs]/I");

  MyTree->Branch("mucoll_pt"                          , mucoll_pt                          , "mucoll_pt[NMUONs]/D");
  MyTree->Branch("mucoll_eta"                         , mucoll_eta                         , "mucoll_eta[NMUONs]/D");
  MyTree->Branch("mucoll_phi"                         , mucoll_phi                         , "mucoll_phi[NMUONs]/D");

  MyTree->Branch("hscp_trk_valid_frac"                , hscp_trk_valid_frac                , "hscp_trk_valid_frac[NHSCPs]/D");
  MyTree->Branch("hscp_p"                             , hscp_p                             , "hscp_p[NHSCPs]/D");
  MyTree->Branch("hscp_pt_trkref"                     , hscp_pt_trkref                     , "hscp_pt_trkref[NHSCPs]/D");
  MyTree->Branch("hscp_pt_muref"                      , hscp_pt_muref                      , "hscp_pt_muref[NHSCPs]/D");
  MyTree->Branch("hscp_pt_comb_mu"                    , hscp_pt_comb_mu                    , "hscp_pt_comb_mu[NHSCPs]/D");
  MyTree->Branch("hscp_pt_inner_trk"                  , hscp_pt_inner_trk                  , "hscp_pt_inner_trk[NHSCPs]/D");
  MyTree->Branch("hscp_pterr_trkref"                  , hscp_pterr_trkref                  , "hscp_pterr_trkref[NHSCPs]/D");
  MyTree->Branch("hscp_track_eta"                     , hscp_track_eta                     , "hscp_track_eta[NHSCPs]/D");
  MyTree->Branch("hscp_muon_eta"                      , hscp_muon_eta                      , "hscp_muon_eta[NHSCPs]/D");
  MyTree->Branch("hscp_track_phi"                     , hscp_track_phi                     , "hscp_track_phi[NHSCPs]/D");
  MyTree->Branch("hscp_muon_phi"                      , hscp_muon_phi                      , "hscp_muon_phi[NHSCPs]/D");
  MyTree->Branch("hscp_num_valid_tracker_hits"        , hscp_num_valid_tracker_hits        , "hscp_num_valid_tracker_hits[NHSCPs]/I");
  MyTree->Branch("hscp_num_valid_pixel_hits"          , hscp_num_valid_pixel_hits          , "hscp_num_valid_pixel_hits[NHSCPs]/I");
  MyTree->Branch("hscp_num_valid_strip_hits"          , hscp_num_valid_strip_hits          , "hscp_num_valid_strip_hits[NHSCPs]/I");
  MyTree->Branch("hscp_num_valid_mu_hits"             , hscp_num_valid_mu_hits             , "hscp_num_valid_mu_hits[NHSCPs]/I");
  MyTree->Branch("hscp_num_valid_dt_hits"             , hscp_num_valid_dt_hits             , "hscp_num_valid_dt_hits[NHSCPs]/I");
  MyTree->Branch("hscp_num_valid_csc_hits"            , hscp_num_valid_csc_hits            , "hscp_num_valid_csc_hits[NHSCPs]/I");
  MyTree->Branch("hscp_I_as"                          , hscp_I_as                          , "hscp_I_as[NHSCPs]/D");
  MyTree->Branch("hscp_I_h"                           , hscp_I_h                           , "hscp_I_h[NHSCPs]/D");
  MyTree->Branch("hscp_ibeta"                         , hscp_ibeta                         , "hscp_ibeta[NHSCPs]/D");
  MyTree->Branch("hscp_ibeta_err"                     , hscp_ibeta_err                     , "hscp_ibeta_err[NHSCPs]/D");
  MyTree->Branch("hscp_calo_e_over_trk_p"             , hscp_calo_e_over_trk_p             , "hscp_calo_e_over_trk_p[NHSCPs]/D");
  MyTree->Branch("hscp_chi2_ndof"                     , hscp_chi2_ndof                     , "hscp_chi2_ndof[NHSCPs]/D");
  MyTree->Branch("hscp_dxy"                           , hscp_dxy                           , "hscp_dxy[NHSCPs]/D");
  MyTree->Branch("hscp_dz"                            , hscp_dz                            , "hscp_dz[NHSCPs]/D");
  MyTree->Branch("hscp_dz_obs"                        , hscp_dz_obs                        , "hscp_dz_obs[NHSCPs]/D");
  MyTree->Branch("hscp_dxy_obs"                       , hscp_dxy_obs                       , "hscp_dxy_obs[NHSCPs]/D");
  MyTree->Branch("hscp_trk_sum_pt"                    , hscp_trk_sum_pt                    , "hscp_trk_sum_pt[NHSCPs]/D");
  MyTree->Branch("hscp_qual_ind"                      , hscp_qual_ind                      , "hscp_qual_ind[NHSCPs]/I");
  MyTree->Branch("hscp_num_tof_meas"                  , hscp_num_tof_meas                  , "hscp_num_tof_meas[NHSCPs]/D");
  MyTree->Branch("hscp_num_dedx_meas"                 , hscp_num_dedx_meas                 , "hscp_num_dedx_meas[NHSCPs]/D");
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

//  For one input file
  fwlite::Event ev(&f);

//For several input files (members of string vector fileNames)

//  fwlite::ChainEvent ev(fileNames);
  std::cout << "ev size = " << ev.size() << std::endl;

//  ofstream file;
//  file.open("out_info.txt");

  for (Long64_t index = 0; index < ev.size(); index++)
// for (Long64_t index = 0; index < 50000; index++)
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
    label = 4;

    /*

        if (run_num == 279694)
        {
          if (lumi_sec == 113 )
          {
            if (event_num != 204899688) continue;
          }
          else if (lumi_sec == 488)
          {
            if (event_num != 887111876) continue;
          }
          else continue;
        }
        else if (run_num == 279766)
        {
          if (lumi_sec != 332) continue;

          if (event_num != 623714738) continue;
        }
        else if (run_num == 279931)
        {
          if (lumi_sec != 240) continue;

          if (event_num != 301499938) continue;
        }
        else if (run_num == 279975)
        {
          if (lumi_sec == 531 )
          {
            if (event_num != 852141622) continue;
          }
          else if (lumi_sec == 692)
          {
            if (event_num != 1.144e+09) continue;
          }
          else continue;
        }
        else if (run_num == 280242)
        {
          if (lumi_sec != 327) continue;

          if (event_num != 592462710) continue;
        }
        else if (run_num == 281707)
        {
          if (lumi_sec != 1039) continue;

          if (event_num != 1.719e+09) continue;
        }
        else if (run_num == 282735)
        {
          if (lumi_sec != 144) continue;

          if (event_num != 293685091) continue;
        }
        else if (run_num == 282919)
        {
          if (lumi_sec != 23) continue;

          if (event_num != 45951248) continue;
        }
        else if (run_num == 283059)
        {
          if (lumi_sec != 370) continue;

          if (event_num != 697932104) continue;
        }
        if (run_num == 283270)
        {
          if (lumi_sec == 1039 )
          {
            if (event_num != 1.713e+09) continue;
          }
          else if (lumi_sec == 477)
          {
            if (event_num != 717611856) continue;
          }
          else continue;
        }
        else if (run_num == 283283)
        {
          if (lumi_sec != 440) continue;

          if (event_num != 893765154) continue;
        }

        else if (run_num == 283478)
        {
          if (lumi_sec != 172) continue;

          if (event_num != 187762323) continue;
        }
        else if (run_num == 283548)
        {
          if (lumi_sec != 174) continue;

          if (event_num != 77304793) continue;
        }
        else if (run_num == 283865)
        {
          if (lumi_sec != 830) continue;

          if (event_num != 1.468e+09) continue;
        }
        else if (run_num == 284037)
        {
          if (lumi_sec != 148) continue;

          if (event_num != 281365725) continue;
        }

    */


    //Collection of hscp candidates
    fwlite::Handle<susybsm::HSCParticleCollection> hscpCollH;
    hscpCollH.getByLabel(ev, "HSCParticleProducer");
    if (!hscpCollH.isValid()) continue;
    const susybsm::HSCParticleCollection& hscpColl = *hscpCollH;


    //Trigger Results
    fwlite::Handle<edm::TriggerResults> trigResults;
    trigResults.getByLabel(ev, "TriggerResults");
    if (!trigResults.isValid()) {cout << "Trig res coll is invalid" << endl; continue;}
    //const edm::TriggerNames& trigNames = *trigResults;
    const edm::TriggerNames& trigNames = ev.triggerNames(*trigResults);
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
    // if (vertexColl.size() < 1) {printf("NO VERTEX\n"); return;}

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


    //const edm::TriggerNames& trigNames = ev.triggerNames(*trigResults);




//   for (unsigned int i = 0, n = trigResults->size(); i < n; ++i)
//   {
    //   std::string nameTriggerAll = trigNames.triggerName(i);
    //    cout << "name = " << nameTriggerAll << endl;



    // std::string nameTrigger= "HLT_IsoMu27_v13";  //1st Single muon trigger
    //       std::string nameTrigger = "HLT_Mu50_v11";      //2nd Single muon trigger



    //    std::string nameTrigger1 = "HLT_Mu18_Mu9_v2";   //Double muon trigger

    //    std::string nameTrigger2 = "HLT_PFMET200_NotCleaned_v5";
//
//              std::string nameTrigger3 = "HLT_CaloMET100_NotCleaned_v3";
//
//              std::string nameTrigger4 = "HLT_MET105_IsoTrk50_v6";


    //     bool passTrigSingle=trigResults->accept(trigNames.triggerIndex(nameTrigger));

    //     bool passTrigDouble=trigResults->accept(trigNames.triggerIndex(nameTrigger1));

    //     bool passTrigPF1=trigResults->accept(trigNames.triggerIndex(nameTrigger2));

//               bool passTrigCalo=trigResults->accept(trigNames.triggerIndex(nameTrigger3));

//              bool passTrigMET=trigResults->accept(trigNames.triggerIndex(nameTrigger4));



    //   if(passTrigSingle) pass_single = true;

    //   if(passTrigDouble) pass_double = true;
//   }









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

      NMUONs++;
    }


    NHSCPs = 0;
    /*

       for (unsigned int v = 0; v < vertexColl.size(); v++)
        {
          cout << "rho = " << vertexColl[v].position().rho() << endl;
      if(vertexColl[i].isFake() || abs(vertexColl[i].z())>24 || vertexColl[i].position().rho()>2 || vertexColl[i].ndof()<=4)continue;

      cout << "v dz = " <<

        }

    */


    for (unsigned int c = 0; c < hscpColl.size(); c++) {

      susybsm::HSCParticle hscp = hscpColl[c];
      reco::MuonRef muonr  = hscp.muonRef();
      reco::TrackRef track = hscp.trackRef();

      //skip hscp candidates that do not have track
      if (track.isNull()) continue;

      //skip hscp candidates that do not have muon reference

      if (muonr.isNull() || !muonr->isStandAloneMuon())continue;

//     if (track->hitPattern().numberOfValidPixelHits() <= 1)         continue;
//      if (track->pt() < 55)    continue;
//      if  (track->eta() > 2.5) continue;
//      if  (track->eta() < -2.5)  continue;
      //if (track->ptError() / track->pt() > 0.25) continue;
      const reco::HitPattern& hp = track->hitPattern();


      susybsm::HSCPIsolation hscpIso = IsolationMap.get((size_t)track.key());
      /*

            if (track->pt() < 55) continue;
            if (track->eta() > 2.1) continue;
            if (track->eta() < -2.1) continue;

            if (abs(track->dz(vertexColl[0].position())) > 0.5) continue;
            if (abs(track->dxy(vertexColl[0].position())) > 0.5) continue;

            if (track->validFraction() < 0.8) continue;
            if ((track->ptError() / track->pt()) > 0.25 ) continue;
            if ((track->chi2() / track->ndof()) > 5) continue;
            if (hp.numberOfValidStripHits()< 6) continue;
            if (hp.numberOfValidPixelHits()< 2) continue;
            if (((hscpIso.Get_ECAL_Energy() + hscpIso.Get_HCAL_Energy()) / track->p()) < 0.3) continue;
            if (hscpIso.Get_TK_SumEt() > 50) continue;

        */


      //  cout << "dz = " << track->dz(vertexColl[0].position()) << endl;



      hscp_p                         [NHSCPs]  = track->p();
      hscp_pt_trkref                 [NHSCPs]  = track->pt();
      hscp_pt_muref                  [NHSCPs]  = muonr->pt();
      if (muonr->combinedMuon().isNonnull())
      {
        hscp_pt_comb_mu              [NHSCPs]  = muonr->combinedMuon()->pt();
      }
      else hscp_pt_comb_mu           [NHSCPs] = -1;
      if (muonr->innerTrack().isNonnull())
      {
        hscp_pt_inner_trk            [NHSCPs] = muonr->innerTrack()->pt();
      }
      else hscp_pt_inner_trk         [NHSCPs] = -1;

      hscp_pterr_trkref              [NHSCPs] = track->ptError() / track->pt();
      hscp_trk_valid_frac            [NHSCPs] = track->validFraction();
      hscp_track_eta                 [NHSCPs] = track->eta();
      hscp_muon_eta                  [NHSCPs] = muonr->eta();
      //  cout << "hscp_eta = " << hscp_track_eta [NHSCPs] << endl;
      hscp_track_phi                 [NHSCPs] = track->phi();
      hscp_muon_phi                  [NHSCPs] = muonr->phi();
      hscp_num_valid_pixel_hits      [NHSCPs] = hp.numberOfValidPixelHits();
      hscp_num_valid_strip_hits      [NHSCPs] = hp.numberOfValidStripHits();
      hscp_num_valid_tracker_hits    [NHSCPs] = hp.numberOfValidTrackerHits();
      hscp_num_valid_mu_hits         [NHSCPs] = hp.numberOfValidMuonHits();
      hscp_num_valid_dt_hits         [NHSCPs] = hp.numberOfValidMuonDTHits();
      hscp_num_valid_csc_hits        [NHSCPs] = hp.numberOfValidMuonCSCHits();
      hscp_chi2_ndof                 [NHSCPs] = track->chi2() / track->ndof();
      hscp_dxy_obs                   [NHSCPs] = track->dxy(vertexColl[0].position());
      hscp_dz_obs                    [NHSCPs] = track->dz(vertexColl[0].position());

     // cout << "default dz = " << track->dz(vertexColl[0].position()) << endl;


      int highestPtGoodVertex = -1;

      double dzmin = 10000;


      for (unsigned int v = 0; v < vertexColl.size(); v++)
      {
        // cout << "rho = " << vertexColl[v].position().rho() << endl;
        if (vertexColl[v].isFake() || abs(vertexColl[v].z()) > 24 || vertexColl[v].position().rho() > 2 || vertexColl[v].ndof() <= 4) { continue;}

//  int highestPtGoodVertex = -1;

        //      double dzmin = 10000;
        if (abs(track->dz (vertexColl[v].position())) < abs(dzmin) ) {
          dzmin = fabs(track->dz (vertexColl[v].position()));
          highestPtGoodVertex = v;

        }
      }
      if (highestPtGoodVertex < 0) {highestPtGoodVertex = 0;}
    //  cout << "needed index = " << highestPtGoodVertex << endl;
    //  cout << "recalc dz = " << track->dz(vertexColl[highestPtGoodVertex].position()) << endl;
      hscp_dxy                       [NHSCPs] = track->dxy(vertexColl[highestPtGoodVertex].position());
      hscp_dz                        [NHSCPs] = track->dz(vertexColl[highestPtGoodVertex].position());

      string muonref, trackref, isstandalone, isglobal, iscombined;
      muonr.isNull() ? muonref = "no" : muonref = "yes";
      track.isNull() ? trackref = "no" : trackref = "yes";
      (muonr.isNonnull() && muonr->combinedMuon().isNonnull()) ? iscombined = "yes" : iscombined = "no";
      (muonr.isNonnull() && muonr->standAloneMuon().isNonnull()) ? isstandalone = "yes" : isstandalone = "no";
      muonr->isGlobalMuon() ? isglobal = "yes" : isglobal = "no";

      double muon_segment_compatibility = -110;
      double muon_dxy_cosm = -110;
      double muon_dz_cosm = - 110;
      for (unsigned int i = 0; i < muon_pt_fromcoll.size(); i++)
      {
        if (muonr->pt() == muon_pt_fromcoll[i])
        {
          muon_segment_compatibility = segment_compatibility[i];
          //   muon_dxy_cosm = muon_dxy[i];
          // muon_dz_cosm = muon_dz[i];
        }

      }
      if (muon_segment_compatibility == -100) cout << "SOMETHING IS GOING WRONG WITH SEGMENT COMPATIBILITY \n";

      //define Loose muon
      bool isLoose = muonr->isPFMuon() && (muonr->isGlobalMuon() || muonr->isTrackerMuon());

      //define Medium muon
      bool goodGlobal = muonr->isGlobalMuon() && muonr->globalTrack().isNonnull() && muonr->globalTrack()->normalizedChi2() < 3 &&
                        muonr->combinedQuality().chi2LocalPosition < 12 &&
                        muonr->combinedQuality().trkKink < 20 ;

      bool isMedium = isLoose && muonr->innerTrack()->validFraction() > 0.8 &&
                      muon_segment_compatibility > (goodGlobal ? 0.303 : 0.451);

      //define Tight muon
      bool isTight = muonr->isGlobalMuon() &&  muonr->isPFMuon() && muonr->globalTrack().isNonnull() && muonr->globalTrack()->normalizedChi2() < 10. &&
                     muonr->globalTrack()->hitPattern().numberOfValidMuonHits() > 0 && muonr->numberOfMatchedStations() > 1 &&
                     muonr->muonBestTrack().isNonnull() &&
                     abs(track->dxy(vertexColl[c].position())) < 0.2 &&
                     abs(track->dz(vertexColl[c].position())) < 0.5 &&
                     muonr->innerTrack().isNonnull() &&
                     muonr->innerTrack()->hitPattern().numberOfValidPixelHits() > 0 &&
                     muonr->innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5;

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
      tofCalculator.computeTOF(muonr, CSCSegmentColl, DTSegmentColl, true );
      tof  = &tofCalculator.combinedTOF;
      // dttof = &tofCalculator.dtTOF;
      // csctof = &tofCalculator.cscTOF;

      //get standard 1/beta (from EDM)
      //tof_edm  = &TOFCollH->get(hscp.muonRef().key());
      // dttof_edm = &TOFdtCollH->get(hscp.muonRef().key());
      // csctof_edm = &TOFcscCollH->get(hscp.muonRef().key());

      //   const reco::HitPattern& hp = track->hitPattern();


      //   cout << "number of TOF measurements = " << tof->nDof() << endl;
      //    cout << "err i/beta = " << tof->inverseBetaErr() << endl;
      //    cout << "Number of dEdx measurements = " << dedxSObjTmp.numberOfMeasurements() << endl;

      hscp_I_as                    [NHSCPs] = dedxSObj->dEdx();
      hscp_I_h                     [NHSCPs] = dedxMObj->dEdx();
      hscp_ibeta                   [NHSCPs] = tof->inverseBeta();
      hscp_ibeta_err               [NHSCPs] = tof->inverseBetaErr();
      hscp_calo_e_over_trk_p       [NHSCPs] = (hscpIso.Get_ECAL_Energy() + hscpIso.Get_HCAL_Energy()) / track->p();
      hscp_trk_sum_pt              [NHSCPs] = hscpIso.Get_TK_SumEt();
      hscp_num_tof_meas            [NHSCPs] = tof->nDof();


      if (isLoose && isMedium && isTight) hscp_qual_ind [NHSCPs] = 111;
      else if (isLoose && isMedium && !isTight) hscp_qual_ind[NHSCPs] = 110;
      else if (isLoose && !isMedium && !isTight) hscp_qual_ind [NHSCPs] = 100;
      else if (isLoose && !isMedium && isTight) hscp_qual_ind [NHSCPs] = 101;
      else if (!isLoose && isMedium && isTight) hscp_qual_ind [NHSCPs] = 1011;
      else if (!isLoose && !isMedium && !isTight) hscp_qual_ind [NHSCPs] = 0;

      NHSCPs++;

    } //end of loop over hscp candidates;


    NTRACKs = 0;

    for (size_t i = 0; i < TrackColH->size();  i++)
    {
      const  reco::Track & tr = (*TrackColH)[i];

      //  if (tr.pt() < 55) continue;
      //  if (tr.eta() > 2.5) continue;
      //  if (tr.eta() < -2.5)  continue;

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


      /*
         double dR = -10;

         if (hscp_track_eta.size() > 0)
         {

           //association between hscp_tracks with tracks from reco::TrackCollection
           for (unsigned int t = 0; t <  hscp_track_eta.size(); t++)
           {
             cout << "trkcoll_eta = " << trkcoll_eta [NTRACKs] << endl;
             cout << "hscp track eta = " << hscp_track_eta[t] << endl;

             cout << "trkcoll_phi = " << trkcoll_phi [NTRACKs] << endl;
             cout << "hscp track phi = " << hscp_track_phi[t] << endl;

             double d_eta = -100;

             if (trkcoll_eta [NTRACKs] * hscp_track_eta[t] >= 0)
             {
               d_eta = abs(hscp_track_eta[t]) - abs(trkcoll_eta [NTRACKs]);
             }
             else if (trkcoll_eta [NTRACKs] * hscp_track_eta[t] <  0)
             {
               d_eta = trkcoll_eta [NTRACKs] + hscp_track_eta[t] ;
             }

             double d_phi = -100;
             if (trkcoll_phi [NTRACKs] * hscp_track_phi[t] >= 0)
             {
               d_phi = abs(hscp_track_phi[t]) - abs(trkcoll_phi [NTRACKs]);
             }
             else if (trkcoll_phi [NTRACKs] * hscp_track_phi[t] <  0)
             {
               d_phi = trkcoll_phi [NTRACKs] + hscp_track_phi[t];
             }

             dR = sqrt(pow(d_eta, 2) + pow(d_phi, 2));


             cout << "DR = " << dR << endl;
           }
         }
         cout << "Dr out of asso loop = " << dR << endl;
      */
      NTRACKs++;
    }

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

