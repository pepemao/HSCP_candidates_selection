
#include "TApplication.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLegend.h"

#define _USE_MATH_DEFINES

#include "stdlib.h"
#include "stdio.h"
#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;

void analyze_tree()
{



	gROOT->Reset();
	gROOT->SetStyle("Plain");

	TChain hscp("demo/HscpTree");

	hscp.Add("Out_data_2016_extrahigh.root");

	hscp.SetBranchStatus("*", 0);

	unsigned int    run_num;
	unsigned int    event_num;
	unsigned int    lumi_sec;

	unsigned int    NHSCPs;
	unsigned int    NTRACKs;
	unsigned int    NMUONs;
	
	double hscp_p                         [1000];
	double hscp_pt_trkref                 [1000];
    double hscp_trk_valid_frac            [1000]; 
    double hscp_pt_muref                  [1000]; 
    double hscp_pt_comb_mu                [1000];
    double hscp_pt_inner_trk              [1000];
    double hscp_pterr_trkref              [1000];
    double hscp_track_eta                 [1000];
    double hscp_muon_eta                  [1000];
    double hscp_track_phi                 [1000]; 
    double hscp_muon_phi                  [1000];
    int    hscp_num_valid_tracker_hits    [1000];
    int    hscp_num_valid_pixel_hits      [1000];
    int    hscp_num_valid_strip_hits      [1000];
    double hscp_I_as                      [1000];
    double hscp_I_h                       [1000];
    double hscp_ibeta                     [1000];
    double hscp_ibeta_err                 [1000];
    double hscp_calo_e_over_trk_p         [1000];
    double hscp_chi2_ndof                 [1000]; 
    double hscp_dxy                       [1000];
    double hscp_dz                        [1000];
    double hscp_trk_sum_pt                [1000];
    int    hscp_qual_ind                  [1000];

    double mucoll_pt                      [1000];
    double mucoll_eta                     [1000];
    double mucoll_phi                     [1000];

    int    trkcoll_num_valid_tracker_hits [1000];
    int    trkcoll_num_valid_pixel_hits   [1000]; 
    int    trkcoll_num_valid_strip_hits   [1000];
    int    trkcoll_num_valid_muon_hits    [1000];
    int    trkcoll_num_valid_dt_hits      [1000];
    int    trkcoll_num_valid_csc_hits     [1000];
    int    trkcoll_num_valid_rpc_hits     [1000];
    double trkcoll_pt                     [1000];
    double trkcoll_pt_err                 [1000];
    double trkcoll_eta                    [1000];
    double trkcoll_phi                    [1000];
    double trkcoll_chi2_ndof              [1000];
    double trkcoll_valid_frac             [1000];
    double trkcoll_dz                     [1000];
    double trkcoll_dxy                    [1000];

	hscp.SetBranchAddress("run_num", &run_num);
	hscp.SetBranchAddress("event_num", &event_num);
	hscp.SetBranchAddress("lumi_sec", &lumi_sec);

	hscp.SetBranchAddress("NHSCPs", &NHSCPs);
	hscp.SetBranchAddress("NTRACKs", &NTRACKs);
	hscp.SetBranchAddress("NMUONs", &NMUONs);

	hscp.SetBranchAddress("hscp_p", hscp_p);
	hscp.SetBranchAddress("hscp_pt_trkref", hscp_pt_trkref);
	hscp.SetBranchAddress("hscp_trk_valid_frac", hscp_trk_valid_frac);
    hscp.SetBranchAddress("hscp_pt_muref", hscp_pt_muref);
    hscp.SetBranchAddress("hscp_pt_comb_mu", hscp_pt_comb_mu);
    hscp.SetBranchAddress("hscp_pt_inner_trk", hscp_pt_inner_trk);
    hscp.SetBranchAddress("hscp_pterr_trkref", hscp_pterr_trkref);
    hscp.SetBranchAddress("hscp_track_eta", hscp_track_eta);
    hscp.SetBranchAddress("hscp_muon_eta", hscp_muon_eta);
    hscp.SetBranchAddress("hscp_track_phi", hscp_track_phi);
    hscp.SetBranchAddress("hscp_muon_phi", hscp_muon_phi);
    hscp.SetBranchAddress("hscp_num_valid_tracker_hits", hscp_num_valid_tracker_hits);
    hscp.SetBranchAddress("hscp_num_valid_pixel_hits", hscp_num_valid_pixel_hits);
    hscp.SetBranchAddress("hscp_num_valid_strip_hits", hscp_num_valid_strip_hits);
    hscp.SetBranchAddress("hscp_I_as", hscp_I_as);
    hscp.SetBranchAddress("hscp_I_h", hscp_I_h); 
    hscp.SetBranchAddress("hscp_ibeta", hscp_ibeta);
    hscp.SetBranchAddress("hscp_ibeta_err", hscp_ibeta_err);
    hscp.SetBranchAddress("hscp_calo_e_over_trk_p", hscp_calo_e_over_trk_p);
    hscp.SetBranchAddress("hscp_chi2_ndof", hscp_chi2_ndof);
    hscp.SetBranchAddress("hscp_dxy", hscp_dxy);
    hscp.SetBranchAddress("hscp_dz", hscp_dz);
    hscp.SetBranchAddress("hscp_trk_sum_pt", hscp_trk_sum_pt);
    hscp.SetBranchAddress("hscp_qual_ind", hscp_qual_ind);

    hscp.SetBranchAddress("mucoll_pt", mucoll_pt);
    hscp.SetBranchAddress("mucoll_eta", mucoll_eta);
    hscp.SetBranchAddress("mucoll_phi", mucoll_phi);

    hscp.SetBranchAddress("trkcoll_num_valid_tracker_hits", trkcoll_num_valid_tracker_hits);
    hscp.SetBranchAddress("trkcoll_num_valid_pixel_hits", trkcoll_num_valid_pixel_hits);
    hscp.SetBranchAddress("trkcoll_num_valid_strip_hits", trkcoll_num_valid_strip_hits);
    hscp.SetBranchAddress("trkcoll_num_valid_muon_hits", trkcoll_num_valid_muon_hits);
    hscp.SetBranchAddress("trkcoll_num_valid_dt_hits", trkcoll_num_valid_dt_hits);
    hscp.SetBranchAddress("trkcoll_num_valid_csc_hits", trkcoll_num_valid_csc_hits);
    hscp.SetBranchAddress("trkcoll_num_valid_rpc_hits", trkcoll_num_valid_rpc_hits);
    hscp.SetBranchAddress("trkcoll_pt", trkcoll_pt);
    hscp.SetBranchAddress("trkcoll_pt_err", trkcoll_pt_err);
    hscp.SetBranchAddress("trkcoll_eta", trkcoll_eta);
    hscp.SetBranchAddress("trkcoll_phi", trkcoll_phi);
    hscp.SetBranchAddress("trkcoll_chi2_ndof", trkcoll_chi2_ndof);
    hscp.SetBranchAddress("trkcoll_valid_frac", trkcoll_valid_frac);
    hscp.SetBranchAddress("trkcoll_dz", trkcoll_dz);
    hscp.SetBranchAddress("trkcoll_dxy", trkcoll_dxy);

	unsigned int NEntries = hscp.GetEntries();
	cout << " n entr = " << NEntries << endl;

    TFile *HistoFile   = new TFile("./hscp_histos.root","RECREATE");

    TH1F * hscp_pt_track = new TH1F("hscp_pt_track"," HSCP track pt " , 50, 50, 40000 );
    TH1F * hscp_pt_comb_mu = new TH1F("hscp_pt_comb_mu"," HSCP comb mu pt " , 50, 50, 40000 );
    //TH1F * hscp_pt_comb_mu = new TH1F("hscp_pt_comb_mu"," HSCP comb mu pt " , 50, 50, 40000 );
    



	for (Int_t i = 0; i < NEntries; i++)
	{
		hscp.GetEntry(i);

		for (int k = 0; k < NHSCPs; k++)
		{

            if(hscp_pt_trkref[k] < 13000 && hscp_pt_comb_mu[k] < 13000) continue;
            
            if(hscp_pt_trkref[k] < 55 ) continue;
            if(hscp_track_eta[k] > 2.1) continue;
            if(hscp_track_eta[k] < -2.1 ) continue;

            if(hscp_dz[k] > 0.5) continue;
            if(hscp_dxy[k] > 0.5) continue;

            if (hscp_trk_valid_frac[k] < 0.8) continue;
            if (hscp_pterr_trkref[k] > 0.25 ) continue;
            if( hscp_chi2_ndof[k] > 5) continue;
            if(hscp_num_valid_strip_hits[k] < 6) continue;
            if(hscp_num_valid_pixel_hits[k] < 2) continue;
            if(hscp_calo_e_over_trk_p[k] < 0.3) continue;
            if(hscp_trk_sum_pt[k] > 50) continue; 
 //           cout<<"err = " <<hscp_pterr_trkref[k]<<endl;

            cout<<"ev nr = " << event_num << endl;
            cout<<"hscp trk pt = " << hscp_pt_trkref[k] << endl;
            cout<<"hscp_pt comb mu = " << hscp_pt_comb_mu[k] << endl;
		    cout<< "I h = " << hscp_I_h[k] << endl;
            cout<< "I as = "<< hscp_I_as[k] << endl;
            cout<< "ibeta = " << hscp_ibeta[k] << endl;    
         
            hscp_pt_track->Fill(hscp_pt_trkref[k]);
            hscp_pt_comb_mu-<Fill(hscp_pt_comb_mu[k]);
		}
		for (int t = 0; t < NTRACKs; t++)
		{
			//cout<<"track_pt = " << trkcoll_pt[t] << endl;
		}
	}
 hscp_p_1->Write();
HistoFile->Close();
}
