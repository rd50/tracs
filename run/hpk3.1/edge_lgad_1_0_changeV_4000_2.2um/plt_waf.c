#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TStyle.h"
using namespace std;
void temp_plt_waf()
{
	TFile *f_sim = new TFile("./edge_lgad_${VBIAS}V.root");
	TFile *f_sim2 = new TFile("../edge_lgad_1_0_changeV_0000_0um/edge_lgad_${VBIAS}V.root");
	TFile *f_data = new TFile("../../edge_voltage_2019_10_09_12_26_57_HPK-EPI-W2-200-DS-SE5-04.txt.root");
	
	TCanvas *c = new TCanvas("c", "", 100, 100, 800, 600);
	c->Divide(1,1);
	c->cd(1);

	TTree *chain1 = (TTree *)f_sim->Get("edge");
	TH2F *h1 = new TH2F("TRACS","",4000,-0.5,10.5,1000,-10,1000);
	chain1->Draw("volt*1e3*0.53:(time-15.35)>>TRACS","z==${Z}");
		
	TTree *chain2 = (TTree *)f_data->Get("edge");
	TH2F *h2 = new TH2F("data","",4000,-0.5,10.5,1000,-10,1000);
	chain2->Draw("-volt*1e3:(time-11.0)>>data","Vbias==-${VBIAS}&&z==${Z_DATA}");

	TTree *chain3 = (TTree *)f_sim2->Get("edge");
	TH2F *h3 = new TH2F("TRACS2","",4000,-0.5,10.5,1000,-10,1000);
	chain3->Draw("volt*1e3*0.53:(time-15.35)>>TRACS2","z==${Z}");

	h1->GetXaxis()->SetTitle("time[ns]");
	h1->GetYaxis()->SetTitle("volt[mv]");

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
	c->Print("waf_hpk_lgad_com_${Z}_bv${VBIAS}_TF.pdf");


}
