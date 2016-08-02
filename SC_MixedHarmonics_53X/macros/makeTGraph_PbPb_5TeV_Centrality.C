#include "RiceStyle.h"

using namespace std;

double centralityBinCenter[] = {35, 45, 55, 65, 75};

void makeTGraph_PbPb_5TeV_Centrality(){
	
	TFile* file1[10];
	TFile* file2[10];

	for(int i = 0; i < 5; i++){

		file1[i] = new TFile(Form("../rootfiles/SC_MixedHarmonics_PbPb_m3n2_v2_%d.root", i+1));
		file2[i] = new TFile(Form("../rootfiles/SC_MixedHarmonics_PbPb_m4n2_v2_%d.root", i+1));
	}

	TH1D* c2_k_m3n2[5];
	TH1D* c2_m_m3n2[5];

	TH1D* c4_m3n2[5];

	TH1D* c2_k_m4n2[5];
	TH1D* c2_m_m4n2[5];

	TH1D* c4_m4n2[5];

	TGraphErrors* gr1 = new TGraphErrors(5);
	TGraphErrors* gr2 = new TGraphErrors(5);

	for(int i = 0; i < 5; i++){

		c2_k_m3n2[i] = (TH1D*) file1[i]->Get("ana/c2_k");
		c2_m_m3n2[i] = (TH1D*) file1[i]->Get("ana/c2_m");
		c4_m3n2[i] = (TH1D*) file1[i]->Get("ana/c4");

		c2_k_m4n2[i] = (TH1D*) file2[i]->Get("ana/c2_k");
		c2_m_m4n2[i] = (TH1D*) file2[i]->Get("ana/c2_m");
		c4_m4n2[i] = (TH1D*) file2[i]->Get("ana/c4");

		double fourP_m3n2 = c4_m3n2[i]->GetMean();
		double fourP_m3n2_error = c4_m3n2[i]->GetMeanError();

		double fourP_c2_k_m3n2 = c2_k_m3n2[i]->GetMean();
		double fourP_c2_m_m3n2 = c2_m_m3n2[i]->GetMean();

		double total = fourP_m3n2 - fourP_c2_k_m3n2*fourP_c2_m_m3n2;
		double total_error = fourP_m3n2_error;

		gr1->SetPoint(i, centralityBinCenter[i], total);
		gr1->SetPointError(i, 0, total_error);

		double fourP_m4n2 = c4_m4n2[i]->GetMean();
		double fourP_m4n2_error = c4_m4n2[i]->GetMeanError();

		double fourP_c2_k_m4n2 = c2_k_m4n2[i]->GetMean();
		double fourP_c2_m_m4n2 = c2_m_m4n2[i]->GetMean();

		double total = fourP_m4n2 - fourP_c2_k_m4n2*fourP_c2_m_m4n2;
		double total_error = fourP_m4n2_error;

		gr2->SetPoint(i, centralityBinCenter[i], total);
		gr2->SetPointError(i, 0, total_error);

	}


	TFile f1("../dataPoints/PbPb_5TeV_Centrality.root", "RECREATE");
	gr1->Write();
	gr2->Write();

}