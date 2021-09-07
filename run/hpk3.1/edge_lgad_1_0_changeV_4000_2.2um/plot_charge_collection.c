#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TStyle.h"
using namespace std;
void temp_plot_charge_collection()
{
	TFile *f_sim = new TFile("./edge_lgad.root");
	TFile *f_sim2 = new TFile("../edge_lgad_1_0_changeV_0000_0um/edge_lgad_${VBIAS}V.root");
	TFile *f_data = new TFile("../../edge_voltage_2019_10_09_12_26_57_HPK-EPI-W2-200-DS-SE5-04.txt.root");
	TCanvas *c = new TCanvas("c", "", 100, 100, 800, 600);
	c->Divide(1,1);
	c->cd(1);

	TTree *chain1 = (TTree *)f_sim->Get("edge");
	TH2F *h1 = new TH2F("TRACS","",4000,-10.5,60.5,1000,-0.01,25);
	chain1->Draw("Qtot*1e3*0.53/50:(z*1e3)>>TRACS","Vbias==${VBIAS}");
		
	TTree *chain2 = (TTree *)f_data->Get("edge");
	TH2F *h2 = new TH2F("data","",4000,-10.5,60.5,1000,-0.01,25);
	chain2->Draw("-Qtot*1e3/50:(((z*1e3)-11953))>>data","Vbias==-${VBIAS}");

	TTree *chain3 = (TTree *)f_sim2->Get("edge");
	TH2F *h3 = new TH2F("TRACS2","",4000,-10.5,60.5,1000,-0.01,25);
	chain3->Draw("Qtot*1e3*0.53/50:(z*1e3)>>TRACS2","Vbias==${VBIAS}");

	h1->GetXaxis()->SetTitle("z or z'[um]");
	h1->GetYaxis()->SetTitle("charge collection[mC]");

	h1->SetMarkerColor(47);
	h1->SetMarkerStyle(42);
	h1->SetMarkerSize(1);
	h1->SetStats(0);
	h1->Draw();
	h2->SetMarkerColor(2);
	h2->SetMarkerStyle(42);
	h2->SetMarkerSize(1);
	h2->Draw("same");
	h3->SetMarkerColor(3);
	h3->SetMarkerStyle(42);
	h3->SetMarkerSize(1);
	h3->Draw("same");

	TLegend *leg1=new TLegend(0.17,0.71,0.4,0.85,NULL,"brNDC");

	leg1->AddEntry(h3,"${VBIAS}V_TRACS(Y-axis scale factor: 0.53)");	
	leg1->AddEntry(h1,"${VBIAS}V_TRACS_doping_plus(Y-axis scale factor: 0.53)");
	leg1->AddEntry(h2,"${VBIAS}V_TCT(z'=z-11953)");

	leg1->SetBorderSize(1);
	leg1->SetTextFont(62);
	leg1->SetTextSize(0.028);
	leg1->SetLineColor(0);
	leg1->SetLineStyle(1);
	leg1->SetLineWidth(3);
	leg1->SetFillColor(0);
	leg1->SetFillStyle(0);

	leg1->Draw();
	c->Print("CC_hpk_lgad_com_TF_bv${VBIAS}.pdf");


}
