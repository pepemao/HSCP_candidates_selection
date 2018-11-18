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

    typedef unsigned int uint;

    gROOT->Reset();
    gROOT->SetStyle("Plain");

    TChain hscp("demo/HscpTree");

    const int num_mc_samples = 4;

    const char * MC_samples[num_mc_samples] = {"PP308.root", "PP494.root", "GM308.root", "GM494.root" };

    for (int sample = 0; sample < num_mc_samples; sample++)
    {
        hscp.Add(MC_samples[sample]);
    }

    //  hscp.Add("GM494.root");

    hscp.SetBranchStatus("*", 0);

    TFile *OutFile   = new TFile("./hscp_histos.root", "RECREATE");

    int scale = 5;
    TCanvas * HSCP_c = new TCanvas("HSCP_c", "Canv", scale * 128, scale * 96);
    HSCP_c->Divide(1, 1);
    
    uint    run_num;
    uint    event_num;
    uint    lumi_sec;

    int label;
    // label = 1 (PP stau 308 GeV)
    // label = 2 (PP stau 494 GeV)
    // label = 3 (GM stau 308 GeV)
    // label = 4 (PP stau 494 GeV)

    uint    NHSCPs;
    uint    NTRACKs;
    uint    NMUONs;

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

    uint NEntries = hscp.GetEntries();
    cout << " n entr = " << NEntries << endl;

    double Num_init = 25000;
    //PP308, PP494, GM308, GM494
    double Num_pass_trig [4] = {24453, 24791, 24879, 24938};

    double Num_evt [4][6] =  {0};
    // [0] = at least one pass quality cuts
    // [1] = at least one pass quality cuts + pt cut
    // [2] = at least one pass quality cuts + pt cut + eta cut (at least one good cand-te)
    // [3] = one good + second pass quality cuts
    // [4] = one good + second pass quality cuts + pt cut
    // [5] = at least two good cand-s
    //..................
    // [0] [i] = PP 308
    // [1] [i] = PP 494
    // [2] [i] = GM 308
    // [3] [i] = GM 494

    for (uint i = 0; i < NEntries; i++)
    {
        hscp.GetEntry(i);

        bool pass_quality = false;
        bool pass_pt = false;
        bool pass_eta = false;

        vector<int> n_pq;
        vector<int> n_ppt;
        vector<int> n_peta;

        for (uint k = 0; k < NHSCPs; k++)
        {

            pass_quality =     abs(hscp_dz[k]) < 0.5  && abs(hscp_dxy[k]) < 0.5 && hscp_trk_valid_frac[k] > 0.8 && hscp_pterr_trkref[k] < 0.25 &&
                               hscp_chi2_ndof[k] < 5 &&  hscp_num_valid_strip_hits[k] >= 6 && hscp_num_valid_pixel_hits[k] >= 2 &&
                               hscp_calo_e_over_trk_p[k] < 0.3 && hscp_trk_sum_pt[k] < 50 ;

            pass_pt = hscp_pt_trkref[k] > 55;
            pass_eta = abs(hscp_track_eta[k]) < 2.1;

            if (pass_quality)                            {n_pq.push_back(1);}
            if (pass_quality && pass_pt)                 {n_ppt.push_back(1);}
            if (pass_quality && pass_pt && pass_eta)     {n_peta.push_back(1);}
        } // end of the HSCP loop

        switch (label)
        {
        case 1:
        {
            if (n_pq.size() > 0)                                          Num_evt[0][0]++;
            if (n_pq.size() > 0 && n_ppt.size() > 0)                      Num_evt[0][1]++;
            if (n_pq.size() > 0 && n_ppt.size() > 0 && n_peta.size() > 0) Num_evt[0][2]++;
            if (n_pq.size() > 1)                                          Num_evt[0][3]++;
            if (n_pq.size() > 1 && n_ppt.size() > 1)                      Num_evt[0][4]++;
            if (n_pq.size() > 1 && n_ppt.size() > 1 && n_peta.size() > 1) Num_evt[0][5]++;
            break;
        }
        case 2:
        {
            if (n_pq.size() > 0)                                          Num_evt[1][0]++;
            if (n_pq.size() > 0 && n_ppt.size() > 0)                      Num_evt[1][1]++;
            if (n_pq.size() > 0 && n_ppt.size() > 0 && n_peta.size() > 0) Num_evt[1][2]++;
            if (n_pq.size() > 1)                                          Num_evt[1][3]++;
            if (n_pq.size() > 1 && n_ppt.size() > 1)                      Num_evt[1][4]++;
            if (n_pq.size() > 1 && n_ppt.size() > 1 && n_peta.size() > 1) Num_evt[1][5]++;
            break;
        }
        case 3:
        {
            if (n_pq.size() > 0)                                          Num_evt[2][0]++;
            if (n_pq.size() > 0 && n_ppt.size() > 0)                      Num_evt[2][1]++;
            if (n_pq.size() > 0 && n_ppt.size() > 0 && n_peta.size() > 0) Num_evt[2][2]++;
            if (n_pq.size() > 1)                                          Num_evt[2][3]++;
            if (n_pq.size() > 1 && n_ppt.size() > 1)                      Num_evt[2][4]++;
            if (n_pq.size() > 1 && n_ppt.size() > 1 && n_peta.size() > 1) Num_evt[2][5]++;
            break;
        }
        case 4:
        {
            if (n_pq.size() > 0)                                          Num_evt[3][0]++;
            if (n_pq.size() > 0 && n_ppt.size() > 0)                      Num_evt[3][1]++;
            if (n_pq.size() > 0 && n_ppt.size() > 0 && n_peta.size() > 0) Num_evt[3][2]++;
            if (n_pq.size() > 1)                                          Num_evt[3][3]++;
            if (n_pq.size() > 1 && n_ppt.size() > 1)                      Num_evt[3][4]++;
            if (n_pq.size() > 1 && n_ppt.size() > 1 && n_peta.size() > 1) Num_evt[3][5]++;
            break;
        }
        default:
        {
            cout << "wrong label " << endl;
            break;
        }
        }

    } // end of the event loop

    cout << "............PP 308....................\n";
    cout << "1q = " << Num_evt[0][0] << endl;
    cout << "1q_pt = " << Num_evt[0][1] << endl;
    cout << "1q_pt_eta = " << Num_evt[0][2] << endl;
    cout << "1_2_q = " << Num_evt[0][3] << endl;
    cout << "1_2_q_pt = " << Num_evt[0][4] << endl;
    cout << "1_2_q_pt_eta = " << Num_evt[0][5] << endl;
    cout << "............PP 494....................\n";
    cout << "1q = " << Num_evt[1][0] << endl;
    cout << "1q_pt = " << Num_evt[1][1] << endl;
    cout << "1q_pt_eta = " << Num_evt[1][2] << endl;
    cout << "1_2_q = " << Num_evt[1][3] << endl;
    cout << "1_2_q_pt = " << Num_evt[1][4] << endl;
    cout << "1_2_q_pt_eta = " << Num_evt[1][5] << endl;
    cout << "............GM 308....................\n";
    cout << "Num_evt[2][0]  = " <<  Num_init << endl;
    cout << "1q = " << Num_evt[2][0] << endl;
    cout << "1q_pt = " << Num_evt[2][1] << endl;
    cout << "1q_pt_eta = " << Num_evt[2][2] << endl;
    cout << "1_2_q = " << Num_evt[2][3] << endl;
    cout << "1_2_q_pt = " << Num_evt[2][4] << endl;
    cout << "1_2_q_pt_eta = " << Num_evt[2][5] << endl;
    cout << "............GM 494....................\n";
    cout << "1q = " << Num_evt[3][0] << endl;
    cout << "1q_pt = " << Num_evt[3][1] << endl;
    cout << "1q_pt_eta = " << Num_evt[3][2] << endl;
    cout << "1_2_q = " << Num_evt[3][3] << endl;
    cout << "1_2_q_pt = " << Num_evt[3][4] << endl;
    cout << "1_2_q_pt_eta = " << Num_evt[3][5] << endl;

    const int num_sel_steps = 8;
    const char *Steps_names[num_sel_steps] = {"Raw", "Triggered", "1Q", "1Q_Pt", "1Q_Pt_Eta", "1_good_2Q", "1_good_2Q_Pt", "At least 2 good"};

    Double_t eff_PP308[num_sel_steps] = {1.00, Num_pass_trig[0] / Num_init, Num_evt[0][0] / Num_init,
                                         Num_evt[0][1] / Num_init, Num_evt[0][2] / Num_init, Num_evt[0][3] / Num_init,
                                         Num_evt[0][4] / Num_init, Num_evt[0][5] / Num_init};
    Double_t eff_PP494[num_sel_steps] = {1.00, Num_pass_trig[1] / Num_init, Num_evt[1][0] / Num_init,
                                         Num_evt[1][1] / Num_init, Num_evt[1][2] / Num_init, Num_evt[1][3] / Num_init,
                                         Num_evt[1][4] / Num_init, Num_evt[1][5] / Num_init};

    Double_t eff_GM308[num_sel_steps] = {1.00, Num_pass_trig[2] / Num_init, Num_evt[2][0] / Num_init,
                                         Num_evt[2][1] / Num_init, Num_evt[2][2] / Num_init, Num_evt[2][3] / Num_init,
                                         Num_evt[2][4] / Num_init, Num_evt[2][5] / Num_init};

    Double_t eff_GM494[num_sel_steps] = {1.00, Num_pass_trig[3] / Num_init, Num_evt[3][0] / Num_init,
                                         Num_evt[3][1] / Num_init, Num_evt[3][2] / Num_init, Num_evt[3][3] / Num_init,
                                         Num_evt[3][4] / Num_init, Num_evt[3][5] / Num_init};


    //Option 1
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
    gm2->SetLineColor(kGreen - 3);
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
    legend->AddEntry(pp2, " PPstau 494 GeV ");
    legend->AddEntry(gm1, " GMStau 308 GeV ");
    legend->AddEntry(gm2, " GMStau 494 GeV ");
    legend->Draw();
    HSCP_c->Write();

    OutFile->Close();
} // end of main

