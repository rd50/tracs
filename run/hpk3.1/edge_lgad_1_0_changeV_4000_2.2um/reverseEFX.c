#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TStyle.h"
using namespace std;
void reverseEFX()
{
	TFile *f1 = new TFile("./25um/Vbias200V/wf200V");
	TFile *f2 = new TFile("./25um/Vbias120V/wf120V");
	TFile *f3 = new TFile("./25um/Vbias80V/wf80V");
		
	TCanvas *c = new TCanvas("c", "", 100, 100, 800, 600);
	TH1D * horigin0 = (TH1D *)f1->Get("h_d_f_grad_Y");
	TH1D * hcreate0 = new TH1D("200V","",horigin0->GetNbinsX(),horigin0->GetBinLowEdge(1),horigin0->GetBinLowEdge(horigin0->GetNbinsX()));
	TH1D * horigin1 = (TH1D *)f2->Get("h_d_f_grad_Y");
	TH1D * hcreate1 = new TH1D("120V","",horigin1->GetNbinsX(),horigin1->GetBinLowEdge(1),horigin1->GetBinLowEdge(horigin1->GetNbinsX()));	
	TH1D * horigin2 = (TH1D *)f3->Get("h_d_f_grad_Y");
	TH1D * hcreate2 = new TH1D("80V","",horigin2->GetNbinsX(),horigin2->GetBinLowEdge(1),horigin2->GetBinLowEdge(horigin2->GetNbinsX()));	
	for(long i = 1; i <=horigin0->GetNbinsX(); i++)
	{
		hcreate0->SetBinContent(i,horigin0->GetBinContent(horigin0->GetNbinsX()-i+1));
	}
	for(long i = 1; i <=horigin1->GetNbinsX(); i++)
	{
		hcreate1->SetBinContent(i,horigin1->GetBinContent(horigin1->GetNbinsX()-i+1));
	}
	for(long i = 1; i <=horigin2->GetNbinsX(); i++)
	{
		hcreate2->SetBinContent(i,horigin2->GetBinContent(horigin2->GetNbinsX()-i+1));
	}
	hcreate2->SetStats(0);
	hcreate2->SetXTitle("z[um]");
	hcreate2->SetYTitle("Electric field[V/um]");
	hcreate0->SetLineColor(1);
	hcreate1->SetLineColor(2);
	hcreate2->SetLineColor(3);
	hcreate2->Draw();
	hcreate1->Draw("same");
	hcreate0->Draw("same");
	c->BuildLegend();
	c->Print("EF_LGAD.pdf");
}
