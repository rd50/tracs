#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TStyle.h"
using namespace std;
void plt_current()
{
	TFile *f_sim = new TFile("./25um/Vbias200V/current200V_scan0.root");
	TCanvas *c = new TCanvas("c", "", 100, 100, 800, 600);
	c->Divide(1,1);
	c->cd(1);

	TH1D * i_total = (TH1D *)f_sim->Get("i_total");
	TH1D * i_init_elec = (TH1D *)f_sim->Get("i_init_elec");
	TH1D * i_init_hole = (TH1D *)f_sim->Get("i_init_hole");
	TH1D * i_gen_elec = (TH1D *)f_sim->Get("i_gen_elec");
	TH1D * i_gen_hole = (TH1D *)f_sim->Get("i_gen_hole");

	i_total->GetXaxis()->SetTitle("time[s]");
	i_total->GetYaxis()->SetTitle("current[A]");
	i_total->SetTitle("Current distribution after RC Shaping");

	i_total->SetLineColor(1);
	i_total->SetLineWidth(2);
	i_total->SetStats(0);
	i_total->Draw();
	
	i_init_elec->SetLineColor(2);
	i_init_elec->SetLineWidth(2);
	i_init_elec->Draw("same");
	
	i_init_hole->SetLineColor(3);
	i_init_hole->SetLineWidth(2);
	i_init_hole->Draw("same");
	
	i_gen_elec->SetLineColor(4);
	i_gen_elec->SetLineWidth(2);
	i_gen_elec->Draw("same");
	
	i_gen_hole->SetLineColor(6);
	i_gen_hole->SetLineWidth(2);
	i_gen_hole->Draw("same");

	TLegend *leg1=new TLegend(0.5,0.71,0.7,0.88,NULL,"brNDC");

	leg1->AddEntry(i_total,"i_total");
	leg1->AddEntry(i_init_elec,"i_init_elec");
	leg1->AddEntry(i_init_hole,"i_init_hole");
	leg1->AddEntry(i_gen_elec,"i_gen_elec");
	leg1->AddEntry(i_gen_hole,"i_gen_hole");

	leg1->SetBorderSize(1);
	leg1->SetTextFont(80);
	leg1->SetTextSize(0.028);
	leg1->SetLineColor(0);
	leg1->SetLineStyle(1);
	leg1->SetLineWidth(3);
	leg1->SetFillColor(0);
	leg1->SetFillStyle(0);

	leg1->Draw();
	c->Print("current_hpk_lgad_25um_bv200.pdf");


}
