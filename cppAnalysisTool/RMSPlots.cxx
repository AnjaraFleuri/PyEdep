#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <vector>
#include <cstring>
#include <fstream>
#include <iostream>
using namespace std;


// Calculate and print RMS/Mean for each histogram
float calculateRMSMean (TH1F* histo) {
    double rms = histo->GetRMS();
    double mean = histo->GetMean();
    double rms_mean_ratio = (rms / mean)*100;
    return rms_mean_ratio ;
}; 

void RMSPlots() {
    string OutputFilename = "RMS.txt";
    ifstream test_file;
    bool file_exists = false;
    test_file.open(OutputFilename.c_str(), std::ifstream::in);
    file_exists = test_file.is_open();
    test_file.close();
    ofstream rmsfile; 
    rmsfile.open(OutputFilename.c_str(), std::ios::trunc);

    if (file_exists){ 
        rmsfile << "#E" << "\t"<< "L1" << "\t" << "L2" <<"\t"<< "Q+L1" << "\t" << "Q+L2 "<<"\t"<< "Q+L1+H" << "\t" << " Q+L2+H" << endl;
    }

    for(int E = 5; E<=50; E+=5){   
        // Open the ROOT file and get the tree
        std::string output1 = "../Inputs/inputFlatLY/outputfile9Feb";
        std::string output2 = "../Inputs/inputRandomLY/output2April";
        std::string filename1 = output1 + std::to_string(E) + "MeV.root";
        std::string filename2 = output2 + std::to_string(E) + "MeV.root";

        TFile *file1 = TFile::Open(filename1.c_str()); 
        TFile *file2 = TFile::Open(filename2.c_str()); 

        if (!file1 || file1->IsZombie()) {
            std::cerr << "Error: Could not open file1!" << std::endl;
            return;
        }
        if (!file2 || file2->IsZombie()) {
            std::cerr << "Error: Could not open file2!" << std::endl;
            return;
        }
        
        TTree *tree1 = (TTree*)file1->Get("Sim");
        if (!tree1) {
            std::cerr << "Error: Could not find tree 'Sim1'!" << std::endl;
            return;
        }
        TTree *tree2 = (TTree*)file2->Get("Sim");
        if (!tree2) {
            std::cerr << "Error: Could not find tree 'Sim2'!" << std::endl;
            return;
        }


        // Declare variables to hold branch data
        float Q_75_old, Q_75_new ;
        float Lsamp1;
        float Lsamp2;
        float L_180;
        std::array<int,8>   N_parList_old, N_parList_new;
        std::array<float,8> Q_depoList_th_75keV_old, Q_depoList_th_75keV_new;

        // Set branch addresses
        tree1->SetBranchAddress("Q_depoTotal_th_75keV", &Q_75_old);
        tree1->SetBranchAddress("L_depoTotal_avg_180PEpMeV", &L_180);
        tree1->SetBranchAddress("N_parList", &N_parList_old);
        tree1->SetBranchAddress("Q_depoList_th_75keV", &Q_depoList_th_75keV_old);

        tree2->SetBranchAddress("Q_depoTotal_th_75keV", &Q_75_new);
        tree2->SetBranchAddress("L_depoTotal_avg_APEX_WP", &Lsamp1);
        tree2->SetBranchAddress("L_depoTotal_avg_APEX_WP2", &Lsamp2);
        tree2->SetBranchAddress("N_parList", &N_parList_new);
        tree2->SetBranchAddress("Q_depoList_th_75keV", &Q_depoList_th_75keV_new);

        TH1F *h_L180 = new TH1F("h_L180", "", 100, 0, E+5);
        TH1F *h_L1 = new TH1F("h_L1", "", 100, 0, E+5);
        TH1F *h_L2 = new TH1F("h_L2", "", 100, 0, E+5);
        TH1F *h_Q_L180 = new TH1F("h_Q_L180", " ", 100, 0, E+3);
        TH1F *h_Q_L1 = new TH1F("h_Q_L1", " ", 100, 0, E+3);
        TH1F *h_Q_L2 = new TH1F("h_Q_L2", " ", 100, 0, E+3);
        TH1F *h_Q_L180_Hadrons = new TH1F("h_Q_L180_Hadrons", "", 100, 0, E+3);
        TH1F *h_Q_L1_Hadrons = new TH1F("h_Q_L1_Hadrons", "", 100, 0, E+3);
        TH1F *h_Q_L2_Hadrons = new TH1F("h_Q_L2_Hadrons", "", 100, 0, E+3);


        // Loop over tree entries
        Long64_t nEntries1 = tree1->GetEntries();
        for (Long64_t i = 0; i < nEntries1; ++i) {
            tree1->GetEntry(i);

            // Compute Q+L for Phase II
            float Q = Q_75_old;
            float L = L_180;
            
            // Compute multiplicity terms
            int N_p = 0; // Number of protons with Q > 75 keV
            int N_n = 0; // Number of neutrons with Q > 75 keV
            int N_alpha = 0; // Number of alphas with Q > 75 keV
            int N_n_captured = 0; // Number of neutrons captured

            if ((N_parList_old)[2]>0) N_n_captured +=N_parList_old[2];
            if ( (Q_depoList_th_75keV_old)[1]>0) N_p +=N_parList_old[1]; // Proton
            else if ( (Q_depoList_th_75keV_old)[2]>0) N_n +=N_parList_old[2]; // Neutron
            else if ( (Q_depoList_th_75keV_old)[6]>0) N_alpha +=N_parList_old[6]; // Alpha

            // Energy in L180, Q+L, Q+L+hadrons  
            float Q_L180 = (Q + L) / 0.96;
            float E_reco_Q_L180_Hadrons = Q_L180 + 7.35 * (N_p + N_alpha) +  7.9 * N_n - 6.1 * N_n_captured;  
            float L180Odep = L/0.42;

            // Fill  histograms
            h_L180->Fill(L180Odep);
            h_Q_L180->Fill(Q_L180);
            h_Q_L180_Hadrons->Fill(E_reco_Q_L180_Hadrons);

        }

        Long64_t nEntries2 = tree2->GetEntries();
        for (Long64_t i = 0; i < nEntries2; ++i) {
            tree2->GetEntry(i);

            // Compute Q+L for Phase II
            float Q = Q_75_new;
            float L1 = Lsamp1;
            float L2 = Lsamp2;
            
            // Compute multiplicity terms
            int N_p = 0; // Number of protons with Q > 75 keV
            int N_n = 0; // Number of neutrons with Q > 75 keV
            int N_alpha = 0; // Number of alphas with Q > 75 keV
            int N_n_captured = 0; // Number of neutrons captured

            if ((N_parList_new)[2]>0) N_n_captured +=N_parList_new[2];
            if ( (Q_depoList_th_75keV_new)[1]>0) N_p +=N_parList_new[1]; // Proton
            else if ( (Q_depoList_th_75keV_new)[2]>0) N_n +=N_parList_new[2]; // Neutron
            else if ( (Q_depoList_th_75keV_new)[6]>0) N_alpha +=N_parList_new[6]; // Alpha

            // Energy in L180, Q+L, Q+L+hadrons  
            float Q_L1 = (Q + Lsamp1) / 0.96;
            float Q_L2 = (Q + Lsamp2) / 0.96;
            float E_reco_Q_L1_Hadrons = Q_L1 + 7.35 * (N_p + N_alpha) +  7.9 * N_n - 6.1 * N_n_captured;  
            float E_reco_Q_L2_Hadrons = Q_L2 + 7.35 * (N_p + N_alpha) +  7.9 * N_n - 6.1 * N_n_captured;  
            float L1dep = Lsamp1/0.42;
            float L2dep = Lsamp2/0.42;

            // Fill Phase II histograms
            h_L1->Fill(L1dep);
            h_L2->Fill(L2dep);
            h_Q_L1->Fill(Q_L1);
            h_Q_L2->Fill(Q_L2);
            h_Q_L1_Hadrons->Fill(E_reco_Q_L1_Hadrons);
            h_Q_L2_Hadrons->Fill(E_reco_Q_L2_Hadrons);

        }


        // Gaussian fits for histograms (focusing on the largest peak)
        auto fitHistogram = [&](TH1F *hist1, TH1F *hist2, TH1F *hist3, std::string name, int color1, int color2, int color3, std::string legendText1, std::string legendText2, std::string legendText3, TLegend *legend, int e) {
                
                int maxBin = hist1->GetMaximumBin();
                double yMax = hist1->GetBinContent(maxBin);
                double xMax = hist1->GetBinCenter(maxBin);
                std::string truEnerg = "True E_{#nu}="+std::to_string(e)+" MeV";
                TLatex *label = new TLatex(E/3, yMax/2.0, truEnerg.c_str());
                label->SetTextColor(kBlue);
                label->SetTextSize(0.04);
                hist1->SetMaximum(yMax+400);
                hist1->SetLineColor(color1);
                hist2->SetLineColor(color2);
                hist3->SetLineColor(color3);
                hist1->SetLineWidth(3);
                hist2->SetLineWidth(3);
                hist3->SetLineWidth(3);
                hist1->SetXTitle("Reco Energy [MeV]");
                hist1->SetYTitle("Events");
                hist1->Draw("HIST");
                hist2->Draw("HIST SAME");
                hist3->Draw("HIST SAME");
                label->Draw();
                legend->AddEntry(hist1, legendText1.c_str(), "l");
                legend->AddEntry(hist2, legendText2.c_str(), "l");
                legend->AddEntry(hist3, legendText3.c_str(), "l");
        };

        gStyle->SetOptStat(0);

        TLegend *legend1 = new TLegend(0.2, 0.8, 0.45, 0.88);
        TLegend *legend2 = new TLegend(0.2, 0.8, 0.45, 0.88);
        TLegend *legend3 = new TLegend(0.2, 0.8, 0.45, 0.88);

        legend1->SetTextSize(0.03);
        legend2->SetTextSize(0.03);
        legend3->SetTextSize(0.03);

        // Remove the border
        legend1->SetBorderSize(0);
        legend2->SetBorderSize(0);
        legend3->SetBorderSize(0);

        // Set the fill style to 0 (transparent)
        legend1->SetFillStyle(0);
        legend2->SetFillStyle(0);
        legend3->SetFillStyle(0);

        TCanvas *c = new TCanvas("c", "", 1400, 400);
        c->Divide(3, 1);

        c->cd(1); 
        fitHistogram(h_L180, h_L1, h_L2, "L", kBlue, kRed, kGreen, "L180/0.42 ", "L/0.42  No Uncertainty", "L/0.42  Pos. Uncertainty", legend1, E);
        legend1->Draw();

        c->cd(2);
        fitHistogram(h_Q_L180, h_Q_L1, h_Q_L2, "Q+L", kBlue, kRed, kGreen, "(Q+L180)/0.96 ", "(Q+L)/0.96  No Uncertainty", "(Q+L)/0.96  Pos. Uncertainty", legend2, E);
        legend2->Draw();

        c->cd(3);
        fitHistogram(h_Q_L180_Hadrons, h_Q_L1_Hadrons, h_Q_L2_Hadrons, "Q+L+multi", kBlue, kRed, kGreen, "Q+L180+#alpha+n+p  ", "Q+L+#alpha+n+p  No Uncertainty", "Q+L+#alpha+n+p  Pos. Uncertainty", legend3, E);
        legend3->Draw();

        c->SetFrameLineWidth(3);

        std::string plotnam1 = "./Plots/result";
        c->SaveAs((plotnam1+E+".png"));  

        float rmsL1 = calculateRMSMean(h_L1);
        float rmsL2 =  calculateRMSMean(h_L2);
        float rmsQL1 = calculateRMSMean(h_Q_L1);
        float rmsQL2 = calculateRMSMean(h_Q_L2);
        float rmsQL1H =  calculateRMSMean(h_Q_L1_Hadrons);
        float rmsQL2H = calculateRMSMean(h_Q_L2_Hadrons);

        rmsfile << setprecision(3);
        rmsfile << E << "\t"<< rmsL1 <<"\t" << rmsL2 <<"\t"<< rmsQL1 << "\t" << rmsQL2 <<"\t"<< rmsQL1H << "\t" << rmsQL2H << endl;


        // Clean up
        delete c; 
        
        delete h_L180 ;
        delete h_L1 ;
        delete h_L2 ;
        delete h_Q_L1;
        delete h_Q_L2;
        delete h_Q_L180;
        delete h_Q_L1_Hadrons;
        delete h_Q_L2_Hadrons;
        delete h_Q_L180_Hadrons;
        file1->Close();
        file2->Close();
        delete file1;
        delete file2;
    }
}
