
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
#include "stdlib.h"
#include "stdio.h"
#include <iostream>
#include <math.h>

using namespace std;

void analyze_tree()
{
    gROOT->Reset();
    gROOT->SetStyle("Plain");

    TChain hscp("demo/HscpTree");

    const int num_mc_samples = 4;

    const char * MC_samples[num_mc_samples] = {"PP308.root", "PP494.root", "GM308.root", "GM494.root" };

    for (int sample = 0; sample < num_mc_samples; sample++)
    {
        hscp.Add(MC_samples[sample]);
    }

    hscp.SetBranchStatus("*", 0);

    TFile *OutFile   = new TFile("./hscp_histos.root", "RECREATE");

    int scale = 5;
    TCanvas * HSCP_c = new TCanvas("HSCP_c", "Canv", scale * 128, scale * 96);
    HSCP_c->Divide(1, 1);

    unsigned int    run_num;
    unsigned int    event_num;
    unsigned int    lumi_sec;

    int label;
    // label = 1 (PP stau 308 GeV)
    // label = 2 (PP stau 494 GeV)
    // label = 3 (GM stau 308 GeV)
    // label = 4 (PP stau 494 GeV)

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

    hscp.SetBranchAddress("label", &label);

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


    TH1F * hscp_pt_combmu = new TH1F("hscp_pt_combmu", " HSCP comb mu pt " , 50, 50, 40000 );
    TH1F * dr_tracks = new TH1F("dr_tracks", " dR " , 50, -0.5, 0.5 );
    TH1F * dif_eta = new TH1F("dif_eta", " dEta " , 50, -0.01, 0.01 );
    TH1F * dif_phi = new TH1F("dif_hi", " dPhi " , 50, -0.01, 0.01 );
    TH1F * dif_pt = new TH1F("dif_pt", " dPhi " , 50, -4, 4 );
    TH2F *etaeta = new TH2F("etaeta", "", 40, -3, 3, 40, -3, 3);


    vector<double> hscp_eta;
    vector<double> hscp_phi;
    vector<double> hscp_pt;

    double num_init_evts = 25000;

    double num_pp_308_trig = 24453;
    double num_pp_494_trig = 24791;
    double num_gm_308_trig = 24879;
    double num_gm_494_trig = 24938;

    double pp_308_evt_one_q = 0;
    double pp_308_evt_one_q_pt = 0;
    double pp_308_evt_one_good = 0;
    double pp_308_evt_one_good_second_q = 0;
    double pp_308_evt_one_good_second_q_pt = 0;
    double pp_308_evt_at_least_two_good = 0;
    double pp_308_evt_more_than_two_good = 0;

    double pp_494_evt_one_q = 0;
    double pp_494_evt_one_q_pt = 0;
    double pp_494_evt_one_good = 0;
    double pp_494_evt_one_good_second_q = 0;
    double pp_494_evt_one_good_second_q_pt = 0;
    double pp_494_evt_at_least_two_good = 0;
    double pp_494_evt_more_than_two_good = 0;

    double gm_308_evt_one_q = 0;
    double gm_308_evt_one_q_pt = 0;
    double gm_308_evt_one_good = 0;
    double gm_308_evt_one_good_second_q = 0;
    double gm_308_evt_one_good_second_q_pt = 0;
    double gm_308_evt_at_least_two_good = 0;
    double gm_308_evt_more_than_two_good = 0;

    double gm_494_evt_one_q = 0;
    double gm_494_evt_one_q_pt = 0;
    double gm_494_evt_one_good = 0;
    double gm_494_evt_one_good_second_q = 0;
    double gm_494_evt_one_good_second_q_pt = 0;
    double gm_494_evt_at_least_two_good = 0;
    double gm_494_evt_more_than_two_good = 0;

    for (Int_t i = 0; i < NEntries; i++)
    {
        hscp.GetEntry(i);

        cout << " event = " << i << endl;

        bool pass_quality = false;
        bool pass_pt = false;
        bool pass_eta = false;

        int num_pass_q = 0;
        int num_pass_q_pt = 0;
        int num_pass_q_pt_eta = 0;


        if (NHSCPs < 1) continue;

        for (int k = 0; k < NHSCPs; k++)
        {

            pass_quality = abs(hscp_dz[k]) < 0.5 && abs(hscp_dxy[k]) < 0.5 && hscp_trk_valid_frac[k] > 0.8 && hscp_pterr_trkref[k] < 0.25 &&
                           hscp_chi2_ndof[k] < 5 &&  hscp_num_valid_strip_hits[k] >= 6 && hscp_num_valid_pixel_hits[k] >= 2 &&
                           hscp_calo_e_over_trk_p[k] < 0.3 && hscp_trk_sum_pt[k] < 50 ;

            pass_pt = hscp_pt_trkref[k] > 55 ;
            pass_eta = abs(hscp_track_eta[k]) < 2.1;

            if (pass_quality)
            {num_pass_q++;}

            if (pass_quality && pass_pt)
            {num_pass_q_pt++;}

            if (pass_quality && pass_pt && pass_eta)
            {num_pass_q_pt_eta++;}
        }

        bool one_good = pass_quality && pass_pt && pass_eta;

        bool one_good_second_pq = one_good && num_pass_q > 1;

        bool one_good_second_pq_pt = one_good && num_pass_q_pt > 1;

        bool at_least_two_good = one_good && num_pass_q_pt_eta > 1;

        bool more_than_2_good = num_pass_q_pt_eta > 2;


        switch (label)
        {
        case 1:
        {
            if (pass_quality)               {pp_308_evt_one_q++;}
            if (pass_quality && pass_pt)    {pp_308_evt_one_q_pt++;}
            if (one_good)                   {pp_308_evt_one_good++;}
            if (one_good_second_pq)         {pp_308_evt_one_good_second_q++;}
            if (one_good_second_pq_pt)      {pp_308_evt_one_good_second_q_pt++;}
            if (at_least_two_good)          {pp_308_evt_at_least_two_good++;}
            if (more_than_2_good)           {pp_308_evt_more_than_two_good++;}
            break;
        }
        case 2:
        {
            if (pass_quality)               {pp_494_evt_one_q++;}
            if (pass_quality && pass_pt)    {pp_494_evt_one_q_pt++;}
            if (one_good)                   {pp_494_evt_one_good++;}
            if (one_good_second_pq)         {pp_494_evt_one_good_second_q++;}
            if (one_good_second_pq_pt)      {pp_494_evt_one_good_second_q_pt++;}
            if (at_least_two_good)          {pp_494_evt_at_least_two_good++;}
            if (more_than_2_good)           {pp_494_evt_more_than_two_good++;}
            break;
        }
        case 3:
        {
            if (pass_quality)               {gm_308_evt_one_q++;}
            if (pass_quality && pass_pt)    {gm_308_evt_one_q_pt++;}
            if (one_good)                   {gm_308_evt_one_good++;}
            if (one_good_second_pq)         {gm_308_evt_one_good_second_q++;}
            if (one_good_second_pq_pt)      {gm_308_evt_one_good_second_q_pt++;}
            if (at_least_two_good)          {gm_308_evt_at_least_two_good++;}
            if (more_than_2_good)           {gm_308_evt_more_than_two_good++;}
            break;
        }
        case 4:
        {

            if (pass_quality)               {gm_494_evt_one_q++;}
            if (pass_quality && pass_pt)    {gm_494_evt_one_q_pt++;}
            if (one_good)                   {gm_494_evt_one_good++;}
            if (one_good_second_pq)         {gm_494_evt_one_good_second_q++;}
            if (one_good_second_pq_pt )     {gm_494_evt_one_good_second_q_pt++;}
            if (at_least_two_good && label) {gm_494_evt_at_least_two_good++;}
            if (more_than_2_good)           {gm_494_evt_more_than_two_good++;}
            break;
        }
        default:
        {
            cout << "wrong label " << endl;
            break;
        }
        }

        for (int t = 0; t < NTRACKs; t++)
        {

            double d_eta, d_phi;
            if (hscp_track_eta[t]*trkcoll_eta[t] > 0)
            {
                d_eta = abs(hscp_track_eta[t]) - abs(trkcoll_eta[t]);
            }      else {d_eta = hscp_track_eta[t] + trkcoll_eta[t];}
            if (hscp_track_phi[t]*trkcoll_phi[t] > 0)
            {
                d_phi = abs(hscp_track_phi[t]) - abs(trkcoll_phi[t]);
            }
            else {d_phi = hscp_track_phi[t] + trkcoll_phi[t];}


            double Dr = sqrt(pow(d_eta, 2) + pow(d_phi, 2));
            dr_tracks->Fill(Dr);
            dif_eta->Fill(d_eta);
            dif_phi->Fill(d_phi);
            dif_pt->Fill(hscp_pt_trkref[t] - trkcoll_pt[t]);
            etaeta->Fill(hscp_track_eta[t], trkcoll_eta[t]) ;
        }

    } // end of the event loop
    cout << "N(evts) with 1 good  = " << pp_308_evt_one_good << endl;
    cout << "N(evts) with 1 good and 2nd passed quality cuts = " << pp_308_evt_one_good_second_q << endl;
    cout << "N(evts) with 1 good and 2nd passed quality cuts and pt = " << pp_308_evt_one_good_second_q_pt << endl;
    cout << "N(evts) with at least 2 good candidates = " << pp_308_evt_at_least_two_good << endl;
    cout << "N(evts) with more than 2 good candidates = " << pp_308_evt_more_than_two_good << endl;
    cout << "double(num_pp_494_trig / num_init_evts) = " << double(num_pp_494_trig / num_init_evts) << endl;
    cout << "gm_494_evt_one_good  = " << gm_494_evt_one_good << endl;


    const int num_sel_steps = 8;

    const char *Steps_names[num_sel_steps] = {"Raw", "Triggered", "Track quality", "Pt > 55 GeV", "|eta| < 2.1", "1st good && PQ", "1 good & 2nd PQ PT", "A/l 2 good"};

    Double_t eff_PP308[num_sel_steps] = {1.00, num_pp_308_trig / num_init_evts, pp_308_evt_one_q / num_init_evts,
                                         pp_308_evt_one_q_pt / num_init_evts, pp_308_evt_one_good / num_init_evts, pp_308_evt_one_good_second_q / num_init_evts,
                                         pp_308_evt_one_good_second_q_pt / num_init_evts, pp_308_evt_at_least_two_good / num_init_evts
                                        };
    Double_t eff_PP494[num_sel_steps] = {1.00, num_pp_494_trig / num_init_evts, pp_494_evt_one_q / num_init_evts,
                                         pp_494_evt_one_q_pt / num_init_evts, pp_494_evt_one_good / num_init_evts, pp_494_evt_one_good_second_q / num_init_evts,
                                         pp_494_evt_one_good_second_q_pt / num_init_evts, pp_494_evt_at_least_two_good / num_init_evts
                                        };
    Double_t eff_GM308[num_sel_steps] = {1.00, num_gm_308_trig / num_init_evts, gm_308_evt_one_q / num_init_evts,
                                         gm_308_evt_one_q_pt / num_init_evts, gm_308_evt_one_good / num_init_evts, gm_308_evt_one_good_second_q / num_init_evts,
                                         gm_308_evt_one_good_second_q_pt / num_init_evts, gm_308_evt_at_least_two_good / num_init_evts
                                        };
    Double_t eff_GM494[num_sel_steps] = {1.00, num_gm_494_trig / num_init_evts, gm_494_evt_one_q / num_init_evts,
                                         gm_494_evt_one_q_pt / num_init_evts, gm_494_evt_one_good / num_init_evts, gm_494_evt_one_good_second_q / num_init_evts,
                                         gm_494_evt_one_good_second_q_pt / num_init_evts, gm_494_evt_at_least_two_good / num_init_evts
                                        };

    TH1F *pp1 = new TH1F("pp1", "PPStau 308 GeV", num_sel_steps, 0, num_sel_steps);
    pp1->SetLineColor(kBlue);
    pp1->SetStats(0);

    for (int i = 1; i <= num_sel_steps; ++i)
    {
        pp1->SetBinContent(i, eff_PP308[i - 1]);
        pp1->GetXaxis()->SetBinLabel(i, Steps_names[i - 1]);
    }

    TH1F *pp2 = new TH1F("pp2", "PPStau 494 GeV", num_sel_steps, 0, num_sel_steps);
    pp2->SetLineColor(kRed);
    pp2->SetStats(0);

    for (int i = 1; i <= num_sel_steps; ++i)
    {
        pp2->SetBinContent(i, eff_PP494[i - 1]);
        pp2->GetXaxis()->SetBinLabel(i, Steps_names[i - 1]);
    }

    TH1F *gm1 = new TH1F("gm1", "GMStau 308 GeV", num_sel_steps, 0, num_sel_steps);
    gm1->SetLineColor(kBlack);
     gm1->SetStats(0);

    for (int i = 1; i <= num_sel_steps; ++i)
    {
        gm1->SetBinContent(i, eff_GM308[i - 1]);
        gm1->GetXaxis()->SetBinLabel(i, Steps_names[i - 1]);
    }

    TH1F *gm2 = new TH1F("gm2", "GMStau 494 GeV", num_sel_steps, 0, num_sel_steps);
    gm2->SetLineColor(kGreen-3);
    gm2->SetStats(0);

    for (int i = 1; i <= num_sel_steps; ++i)
    {
        gm2->SetBinContent(i, eff_GM494[i - 1]);
        gm2->GetXaxis()->SetBinLabel(i, Steps_names[i - 1]);
    }

    HSCP_c->cd(1);
    pp1->Draw("E0");
    pp2->Draw("E0Same");
    gm1->Draw("E0Same");
    gm2->Draw("E0Same");

    TLegend * legend = new TLegend(.75, .80,  .95, .95);
    legend->SetHeader("Signal samples: ");
    legend->AddEntry(pp1, " PPstau 308 GeV ");
    legend->AddEntry(pp2, " PP stau 494 GeV ");
    legend->AddEntry(gm1, " GMStau 308 GeV ");
    legend->AddEntry(gm2, " GMStau 494 GeV ");
    legend->Draw();
    HSCP_c->Write();


    hscp_pt_combmu->Write();
    dr_tracks->Write();
    dif_eta->Write();
    dif_phi->Write();
    dif_pt->Write();
    etaeta->Write();
    OutFile->Close();
}
