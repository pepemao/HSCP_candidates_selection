
# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('Demo')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
        limit = cms.untracked.int32(-1)
        )
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.options   = cms.untracked.PSet( SkipEvent = cms.untracked.vstring('ProductNotFound'))


process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring(

		#Below in "'"quotes one should put the path to .root file he wants to anayze in the  form of file:/path/file.root;

               # 'file:/afs/cern.ch/cms/Tutorials/TWIKI_DATA/TTJets_8TeV_53X.root'
#               'file:/eos/cms/store/cmst3/user/querten/15_03_25_HSCP_Run2EDMFiles/Run2016_274199.root'
#		'file:/eos/cms/store/cmst3/user/querten/15_03_25_HSCP_Run2EDMFiles/GMStau_13TeV_M247.root'
		'file:/eos/cms/store/cmst3/user/querten/15_03_25_HSCP_Run2EDMFiles/Run2016_274240.root'
#'file:root://se.cis.gov.pl//cms/store/user/fruboes/HSCP_20170802/Run2016_277992.root'
#'file:root://cms-xrd-global.cern.ch//store/mc/RunIISpring15DR74/HSCPppstau_M-871_TuneZ2star_13TeV_pythia6/MINIAODSIM/Asympt25ns_HSCP_customise_MCRUN2_74_V9-v2/60000/F6DAC86B-5617-E511-9A14-0025907DC9DC.root'
			 )
                            )
#process.demo = cms.EDAnalyzer('DemoAnalyzer',

 #                                        minTracks=cms.untracked.uint32(0)

process.demo = cms.EDAnalyzer('DemoAnalyzer'                                  
                              )

process.TFileService = cms.Service("TFileService",
                                       fileName = cms.string('Out_data_2016_279975.root')
                                   )

process.p = cms.Path(process.demo)
