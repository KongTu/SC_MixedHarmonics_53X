#include "RiceStyle.h"

using namespace std;

double ntrkBinCenter[] = {131.339, 176.239, 196.654, 231.494, 270.541};

void makeTGraph_pPb_5TeV_Ntrk(){
	
	TFile* file1[10];
	TFile* file2[10];

	for(int i = 0; i < 5; i++){

		file1[i] = new TFile(Form("../rootfiles/SC_MixedHarmonics_pPb_v1_%d.root", i+1));
	}

	TH1D* Ntrk[5];

	TH1D* c2_k_m3n2[5];
	TH1D* c2_m_m3n2[5];

	TH1D* c4_m3n2[5];

	TH1D* c2_k_m4n2[5];
	TH1D* c2_m_m4n2[5];

	TH1D* c4_m4n2[5];

	TGraphErrors* gr1 = new TGraphErrors(5);
	TGraphErrors* gr2 = new TGraphErrors(5);

	for(int i = 0; i < 5; i++){

		Ntrk[i] = (TH1D*) file1[i]->Get("ana_m3n2/Ntrk");
		cout << Ntrk[i]->GetMean() << ", ";

		c2_k_m3n2[i] = (TH1D*) file1[i]->Get("ana_m3n2/c2_k");
		c2_m_m3n2[i] = (TH1D*) file1[i]->Get("ana_m3n2/c2_m");
		c4_m3n2[i] = (TH1D*) file1[i]->Get("ana_m3n2/c4");

		c2_k_m4n2[i] = (TH1D*) file1[i]->Get("ana_m4n2/c2_k");
		c2_m_m4n2[i] = (TH1D*) file1[i]->Get("ana_m4n2/c2_m");
		c4_m4n2[i] = (TH1D*) file1[i]->Get("ana_m4n2/c4");

		double fourP_m3n2 = c4_m3n2[i]->GetMean();
		double fourP_m3n2_error = c4_m3n2[i]->GetMeanError();

		double fourP_c2_k_m3n2 = c2_k_m3n2[i]->GetMean();
		double fourP_c2_m_m3n2 = c2_m_m3n2[i]->GetMean();

		double total = fourP_m3n2 - fourP_c2_k_m3n2*fourP_c2_m_m3n2;
		double total_error = fourP_m3n2_error;

		gr1->SetPoint(i, ntrkBinCenter[i], total);
		gr1->SetPointError(i, 0, total_error);

		double fourP_m4n2 = c4_m4n2[i]->GetMean();
		double fourP_m4n2_error = c4_m4n2[i]->GetMeanError();

		double fourP_c2_k_m4n2 = c2_k_m4n2[i]->GetMean();
		double fourP_c2_m_m4n2 = c2_m_m4n2[i]->GetMean();

		double total = fourP_m4n2 - fourP_c2_k_m4n2*fourP_c2_m_m4n2;
		double total_error = fourP_m4n2_error;

		gr2->SetPoint(i, ntrkBinCenter[i], total);
		gr2->SetPointError(i, 0, total_error);

	}


	TFile f1("../dataPoints/pPb_5TeV_Ntrk.root", "RECREATE");
	gr1->Write();
	gr2->Write();

}