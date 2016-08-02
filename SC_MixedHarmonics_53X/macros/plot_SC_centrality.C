#include "RiceStyle.h"

using namespace std;

void plot_SC_centrality(){

	TFile* file = new TFile("../dataPoints/PbPb_5TeV_Centrality.root");

	TGraphErrors* gr1 = (TGraphErrors*) file->Get("Graph;1");
	TGraphErrors* gr2 = (TGraphErrors*) file->Get("Graph;2");

	TCanvas* c1 = new TCanvas("c1","c1",1,1,750,750);
	gPad->SetLeftMargin(0.20);
	gPad->SetBottomMargin(0.13);
	gPad->SetTopMargin(0.06);
	gPad->SetTicks();

	TGaxis::SetMaxDigits(3);

	TH1D* base1 = makeHist("base1", "", "Centrality", "SC(m,n)", 100,0,100,kBlack);

	base1->GetYaxis()->SetRangeUser(-0.0001,0.0001);
	base1->GetXaxis()->SetRangeUser(25, 85);
	base1->GetXaxis()->SetTitleColor(kBlack);
	
	fixedFontHist1D(base1,1.1,1.25);

	base1->GetYaxis()->SetTitleOffset(1.3);
	base1->GetXaxis()->SetTitleOffset(0.95);
	base1->GetYaxis()->SetTitleSize(base1->GetYaxis()->GetTitleSize()*1.3);
	base1->GetXaxis()->SetTitleSize(base1->GetXaxis()->GetTitleSize()*1.4);
	base1->GetYaxis()->SetLabelSize(base1->GetYaxis()->GetLabelSize()*1.4);
	base1->GetXaxis()->SetLabelSize(base1->GetXaxis()->GetLabelSize()*1.4);


	base1->Draw();


	gr1->SetMarkerStyle(20);
	gr1->SetMarkerSize(1.4);
	gr1->SetMarkerColor(kRed);
	gr1->SetLineColor(kRed);
	gr1->Draw("Psame");

	gr2->SetMarkerStyle(21);
	gr2->SetMarkerSize(1.4);
	gr2->SetMarkerColor(kBlue);
	gr2->SetLineColor(kBlue);
	gr2->Draw("Psame");

}