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


float calculateSigma (TH1F* hist, string name, int energy){
     TF1 *gaus = nullptr;
     if(name=="gaus_L" || name=="gaus_QLH"){
     gaus = new TF1(name.c_str(), "gaus", energy-13.0, energy+5.0);
     }
     else{
     gaus = new TF1(name.c_str(), "gaus", energy-9.0, energy+5.0);
     }
     
     hist->Fit(gaus, "RQ");
     float meangaus = gaus->GetParameter(1);
     float sigmgaus = gaus->GetParameter(2);
     float sigmMean = (sigmgaus/meangaus)*100;

     return sigmMean ;
};

void SigmaPlots() {
    string OutputFilename = "Sigma.txt";
    ifstream test_file;
    bool file_exists = false;
    test_file.open(OutputFilename.c_str(), std::ifstream::in);
    file_exists = test_file.is_open();
    test_file.close();
    ofstream sigmafile; 
    sigmafile.open(OutputFilename.c_str(), std::ios::trunc);

    if (!file_exists){ 
        sigmafile << "#E" << "\t"<< "L" << "\t"<< "Q+L" << "\t" << "Q+L+H" << endl;
    }
 
    for (int E = 5; E <= 50; E += 5) {
        // Construct the input file name
        std::string output = "../Inputs/inputRandomLY/output2April";
        std::string filename = output + std::to_string(E) + "MeV.root";

        // Open the ROOT file
        TFile *file = TFile::Open(filename.c_str());
        if (!file || file->IsZombie()) {
            std::cerr << "Error: Could not open file " << filename << "!" << std::endl;
            continue; }

        // Get the tree
        TTree *tree = (TTree *)file->Get("Sim");
        if (!tree) {
            std::cerr << "Error: Could not find tree 'Sim'!" << std::endl;
            file->Close();
            delete file;
            continue;
        }

        // Declare variables to hold branch data
        float Q_75_new;
        float L;
        std::array<int, 8> N_parList_new;
        std::array<float, 8> Q_depoList_th_75keV_new;

        // Set branch addresses
        tree->SetBranchAddress("Q_depoTotal_th_75keV", &Q_75_new);
        tree->SetBranchAddress("L_depoTotal_avg_APEX_WP", &L);
        tree->SetBranchAddress("N_parList", &N_parList_new);
        tree->SetBranchAddress("Q_depoList_th_75keV", &Q_depoList_th_75keV_new);

        // Histograms
        TH1F *h_Q_75_054_new = new TH1F("h_Q_75_054_new", "", 100, 0, E + 9);
        TH1F *h_L = new TH1F("h_L", "", 100, 0, E + 3);
        TH1F *h_Q_L = new TH1F("h_Q_L", "", 100, 0, E + 9);
        TH1F *h_Q_L_Hadrons = new TH1F("h_Q_L_Hadrons", "", 100, 0, E + 9);

        // Loop over tree entries
        Long64_t nEntries = tree->GetEntries();
        for (Long64_t i = 0; i < nEntries; ++i) {
            tree->GetEntry(i);

            // Compute multiplicity terms
            int N_p = 0;
            int N_n = 0;
            int N_alpha = 0;
            int N_n_captured = 0;

            if ((N_parList_new)[2] > 0) N_n_captured += N_parList_new[2];
            if ((Q_depoList_th_75keV_new)[1] > 0) N_p += N_parList_new[1];
            else if ((Q_depoList_th_75keV_new)[2] > 0) N_n += N_parList_new[2];
            else if ((Q_depoList_th_75keV_new)[6] > 0) N_alpha += N_parList_new[6];

            // Phase II
            float Q_L = (Q_75_new + L) / 0.96;
            float E_reco_Q_L_Hadrons = Q_L +  7.35 * (N_p + N_alpha) +  7.9 * N_n -  6.1 * N_n_captured;
            
            h_L->Fill(L / 0.42);
            h_Q_L->Fill(Q_L);
            h_Q_L_Hadrons->Fill(E_reco_Q_L_Hadrons);
        }

    

        // Gaussian fits for histograms (focusing on the largest peak)
        auto fitHistogram = [&](TH1F *hist, std::string name, int color, std::string legendText, TLegend *legend, int e) {
            if (hist->GetEntries() > 0) {
                
                // Find the bin with the largest content
                int nBins = hist->GetNbinsX();
                int lastPeakBin = -1;
                double maxContent = 0;
                for (int bin = nBins; bin >= 1; --bin) {
                    double content = hist->GetBinContent(bin);
                    if (content > maxContent) {
                        maxContent = content;
                        lastPeakBin = bin;
                    }
                }

                if (lastPeakBin == -1) {
                    std::cerr << "Warning: No peak found in histogram " << hist->GetName() << "!" << std::endl;
                    return;
                }

                // Define fit range around the last peak
                double peakCenter = hist->GetBinCenter(lastPeakBin);
                double fitRange = 0.2 * (hist->GetXaxis()->GetXmax() - hist->GetXaxis()->GetXmin());
                double fitMin = std::max(hist->GetXaxis()->GetXmin(), peakCenter - fitRange);
                double fitMax = std::min(hist->GetXaxis()->GetXmax(), peakCenter + fitRange);
                // Apply the fit
                TF1 *gaus = nullptr;
                if(name=="gaus_L" || name=="gaus_PhaseII_multi"){
                gaus = new TF1(name.c_str(), "gaus", e-13.0, e+5.0);
                }
                else{
                gaus = new TF1(name.c_str(), "gaus", e-9.0, e+5.0);
                }
                
                int maxBin = hist->GetMaximumBin();
                double yMax = hist->GetBinContent(maxBin);
                double xMax = hist->GetBinCenter(maxBin);
                std::string truEnerg = "True E_{#nu}="+std::to_string(e)+" MeV";
                TLatex *label = new TLatex(E/3, yMax/2.0, truEnerg.c_str());
                label->SetTextColor(kBlue);
                label->SetTextSize(0.04);
                hist->SetMaximum(yMax+400);
                hist->Fit(gaus, "RQ");
                hist->SetLineWidth(2);
                hist->SetXTitle("Reco Energy [MeV]");
                hist->SetYTitle("Events");
                TF1 *fit = hist->GetFunction(name.c_str());
                if (fit) {
                    fit->SetLineColor(color);
                    fit->SetLineWidth(3);
                    fit->Draw("SAME");
                    label->Draw();
                    legend->AddEntry(hist, legendText.c_str(), "l");
                    legend->AddEntry(fit, "Gauss Fit", "l");
                }
            } else {
                std::cerr << "Warning: Histogram " << hist->GetName() << " is empty!" << std::endl;
            }
        };

        gStyle->SetOptStat(0);

        TCanvas *c = new TCanvas("c", "", 1400, 400);
        c->Divide(3, 1);

        TLegend *legend1 = new TLegend(0.2, 0.8, 0.45, 0.88);
        TLegend *legend2 = new TLegend(0.2, 0.8, 0.45, 0.88);
        TLegend *legend3 = new TLegend(0.2, 0.8, 0.45, 0.88);
        TLegend *legend4 = new TLegend(0.2, 0.8, 0.45, 0.88);
        legend1->SetTextSize(0.04);
        legend2->SetTextSize(0.04);
        legend3->SetTextSize(0.04);
        legend4->SetTextSize(0.04);

        // Remove the border
        legend1->SetBorderSize(0);
        legend2->SetBorderSize(0);
        legend3->SetBorderSize(0);
        legend4->SetBorderSize(0);

        // Set the fill style to 0 (transparent) 
        legend1->SetFillStyle(0);
        legend2->SetFillStyle(0);
        legend3->SetFillStyle(0);
        legend4->SetFillStyle(0);

        c->cd(1);
        fitHistogram(h_L, "gaus_L", kRed, "L180/0.42 ", legend2, E);
        legend2->Draw();

        c->cd(2);
        fitHistogram(h_Q_L, "gaus_QL", kRed, "(Q75+L180)/0.96 ", legend3, E);
        legend3->Draw(); 

        c->cd(3);
        fitHistogram(h_Q_L_Hadrons, "gaus_QLH", kRed, "Q75+L180+hadrons ", legend4, E);
        legend4->Draw();

        // Save the canvas
        std::string plotnam1 = "./Plots/SimaResult";
        c->SaveAs((plotnam1+E+".png"));

        float sigmaL = calculateSigma(h_L, "gaus_L", E);
        float sigmaQL = calculateSigma(h_Q_L, "gaus_QL", E);
        float sigmaQLH =  calculateSigma(h_Q_L_Hadrons, "gaus_QLH", E);


        sigmafile << setprecision(3);
        sigmafile << E << "\t"<< sigmaL << "\t"<< sigmaQL << "\t" << sigmaQLH << endl;

        

        // Clean up
        delete c;
        delete h_Q_75_054_new;
        delete h_L;
        delete h_Q_L;
        delete h_Q_L_Hadrons;
        file->Close();
        delete file;
    }
}











/*  #include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <vector>
#include <iostream>

#include<sys/types.h>
#include<fcntl.h>
#include <math.h>
#include <stdio.h>
#include <cstring>
#include <fstream>
#include <TFile.h>
#include "TROOT.h"
#include <stdlib.h>

void MakePlots() {
   
    for(int E = 5; E<=50; E+=5){
        std::string output = "outpufileRand";
        // Open the ROOT file and get the tree
        TFile *file2 = TFile::Open(output+E+"MeV.root"); // Replace with your file name
        if (!file2 || file2->IsZombie()) {
            std::cerr << "Error: Could not open file2!" << std::endl;
            return;
        }
        
        TTree *tree2 = (TTree*)file2->Get("Sim");
        if (!tree2) {
            std::cerr << "Error: Could not find tree 'Sim2'!" << std::endl;
            return;
        }

        // Declare variables to hold branch data
        float Q_75_new ;
        float L;
        std::array<int,8>  N_parList_new;
        std::array<float,8> Q_depoList_th_75keV_new;

        // Set branch addresses
        tree2->SetBranchAddress("Q_depoTotal_th_75keV", &Q_75_new);
        tree2->SetBranchAddress("L_depoTotal_avg_APEX_WP", &L);
        tree2->SetBranchAddress("N_parList", &N_parList_new);
        tree2->SetBranchAddress("Q_depoList_th_75keV", &Q_depoList_th_75keV_new);

        // Histograms for Phase II and comparison
        TH1F *h_Q_75_054_new = new TH1F("h_Q_75_054_new", "", 100, 0, E+3);
        TH1F *h_L = new TH1F("h_L", "", 100, 0, E+3);
        TH1F *h_Q_L = new TH1F("h_Q_L", "Phase II FD ", 100, 0, E+3);
        TH1F *h_Q_L_Hadrons = new TH1F("h_Q_L_Hadrons", "", 100, 0, E+3);


        // Loop over tree entries
        Long64_t nEntries2 = tree2->GetEntries();
        for (Long64_t i = 0; i < nEntries2; ++i) {
            tree2->GetEntry(i);

            // Compute Q+L for Phase II
            float Q_75_sum = 0.0;
            float L_sum = 0.0;

            Q_75_sum += Q_75_new;
            L_sum += L;
            
            // Compute multiplicity terms
            int N_p = 0; // Number of protons with Q > 75 keV
            int N_n = 0; // Number of neutrons with Q > 75 keV
            int N_alpha = 0; // Number of alphas with Q > 75 keV
            int N_n_captured = 0; // Number of neutrons captured

            if ((N_parList_new)[2]>0) N_n_captured +=N_parList_new[2];
            if ( (Q_depoList_th_75keV_new)[1]>0) N_p +=N_parList_new[1]; // Proton
            else if ( (Q_depoList_th_75keV_new)[2]>0) N_n +=N_parList_new[2]; // Neutron
            else if ( (Q_depoList_th_75keV_new)[6]>0) N_alpha +=N_parList_new[6]; // Alpha

            // Phase II 
            float Q_L_PhaseII = (Q_75_sum + L_sum) / 0.96;
            float E_reco_Q_L_mult_PhaseII = Q_L_PhaseII + 
                                            7.35 * (N_p + N_alpha) + 
                                            7.9 * N_n - 
                                            6.1 * N_n_captured;  

            float QOver = Q_75_sum/0.54;
            float LOver = L_sum/0.42;
            h_Q_75_054_new->Fill(QOver);
            h_L->Fill(LOver);

            // Fill Phase II histograms
            h_Q_L->Fill(Q_L_PhaseII);
            h_Q_L_Hadrons->Fill(E_reco_Q_L_mult_PhaseII);
        }

        // Calculate and print RMS/Mean for each histogram
        auto calculateRMSMean = [](TH1F* histo, const std::string& name) {
            double rms = histo->GetRMS();
            double mean = histo->GetMean();
            double rms_mean_ratio = rms / mean;
            std::cout << std::setprecision(4);
            std::cout << name << ": RMS = " << rms << ", Mean = " << mean 
                    << ", RMS/Mean = " << rms_mean_ratio << std::endl;
        };

        calculateRMSMean(h_Q_75_054_new, "h_Q_75_054_new");
        calculateRMSMean(h_L, "h_L");
        calculateRMSMean(h_Q_L, "h_Q_L");
        calculateRMSMean(h_Q_L_Hadrons, "h_Q_L_Hadrons");
        
        gStyle->SetOptStat(0);
        int maxBin = h_L->GetMaximumBin();
        double yMax = h_L->GetBinContent(maxBin);

        // Create a canvas and draw histograms
        TCanvas *c1 = new TCanvas("c1", "", 800, 600);
        h_Q_75_054_new->SetLineColor(kOrange+1); h_Q_75_054_new->SetLineWidth(4);
        h_L->SetLineColor(kMagenta); h_L->SetLineWidth(3);
        h_Q_75_054_new->Draw("HIST");
        h_Q_75_054_new->SetMaximum(yMax+100);
        h_Q_75_054_new->SetXTitle("Reco Energy [MeV]");
        h_L->Draw("HIST SAME");

        // Add legend 
        TLegend *legend0 = new TLegend(0.15, 0.7, 0.48, 0.9);
        legend0->AddEntry(h_Q_75_054_new, "Q_75/0.54_new", "l");
        legend0->AddEntry(h_L, "RandomLY/0.42", "l");
        legend0->Draw();

        std::string trueE = "True E_{#nu}=";
        TLatex *label = new TLatex(5, 100, trueE+E+" "+"MeV");
        label->SetTextColor(kBlue);
        label->SetTextSize(0.04);
        label->Draw();



        TCanvas *c2 = new TCanvas("c2", "", 800, 600);
        h_Q_L->SetLineColor(kGreen);
        h_Q_L->SetLineWidth(2);
        h_Q_L_Hadrons->SetLineColor(kOrange);
        h_Q_L_Hadrons->SetLineWidth(5);
        int maxBin1 = h_Q_L->GetMaximumBin();
        double yMax1 = h_Q_L->GetBinContent(maxBin1);

        h_Q_L_Hadrons->Draw("HIST");
        h_Q_L_Hadrons->SetMaximum(yMax1+100);
        h_Q_L_Hadrons->SetXTitle("Reco Energy [MeV]");
        h_Q_L->Draw("HIST SAME");

        // Add legend 
        TLegend *legend1 = new TLegend(0.15, 0.7, 0.48, 0.9);
        legend1->AddEntry(h_Q_L, " Q+RandomLY only", "l");
        legend1->AddEntry(h_Q_L_Hadrons, "Q+RandLY+Multiplicity", "l");
        legend1->Draw();

        TLatex *label1 = new TLatex(5, 100, trueE+E+" "+"MeV");
        label1->SetTextColor(kBlue);
        label1->SetTextSize(0.04);
        label1->Draw();

        c1->SetFrameLineWidth(2);
        c2->SetFrameLineWidth(2);
        //c3->SetFrameLineWidth(2);
        // Save the canvases
        std::string plotnam1 = "./Plots/NEWScenPots/RecoEnergy";
        std::string plotnam2 = "./Plots/NEWScenPots/PhaseII_RecoEnergy";
        c1->SaveAs(plotnam1+E+".png");
        c2->SaveAs(plotnam2+E+".png");
        //c3->SaveAs("./Plots/PhaseI_vs_PhaseII_RecoEnergy.png");

        // Clean up
        delete c1;
        delete c2;
        delete h_Q_75_054_new;
        delete h_L ;
        delete h_Q_L;
        delete h_Q_L_Hadrons;
        file2->Close();
        delete file2;
    }
} */
 