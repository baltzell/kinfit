#include <random> 
#include <vector> 
#include "TLorentzVector.h"
#include "TH1F.h"
#include "TGenPhaseSpace.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "KinFitter.h"

// absolute resolutions (GeV,deg,deg):
const std::vector<double> R = {0.005,0.02,0.02};

TRandom3 rndm3(0);
std::random_device rndm_device;
std::default_random_engine rndm_engine(rndm_device());

auto ptr_to_genvector(default_random_engine e, TLorentzVector *v)
{
    // Something very fishy about this std::random number generator.
    // Only get the expected behavior when it's recreated
    // (reseeded?) every time here, which surely can't be right:
    //std::default_random_engine rndm_engine(rndm_device());
    //e=rndm_engine;

    // smear by relative resolutions:
    //std::normal_distribution<double> varP(1.0,R[0]/v->P());
    //std::normal_distribution<double> varTheta(1.0,R[1]/fabs(v->Theta()));
    //std::normal_distribution<double> varPhi(1.0,R[2]/fabs(v->Phi()));
    //const double p = varP(e) * v->P();
    //const double t = varTheta(e) * v->Theta();
    //const double h = varPhi(e) * v->Phi();

    // smear by absolute resolutions:
    //std::normal_distribution<double> varP(0.0,R[0]);
    //std::normal_distribution<double> varTheta(0.0,R[1]);
    //std::normal_distribution<double> varPhi(0.0,R[2]);
    //const double p = varP(e) + v->P();
    //const double t = varTheta(e) + v->Theta();
    //const double h = varPhi(e) + v->Phi();

    // smear by relative resolutions:
    //const double p = rndm3.Gaus(1.0,R[0]/v->P()) * v->P();
    //const double t = rndm3.Gaus(1.0,R[1]/fabs(v->Theta())) * v->Theta();
    //const double h = rndm3.Gaus(1.0,R[2]/fabs(v->Phi())) * v->Phi();

    // smear by absolute resolutions:
    const double p = rndm3.Gaus(0.0,R[0]) + v->P();
    const double t = rndm3.Gaus(0.0,R[1]) + v->Theta();
    const double h = rndm3.Gaus(0.0,R[2]) + v->Phi();

    const double m = v->M();
    TLorentzVector v_out;
    v_out.SetXYZM( p*sin(t)*cos(h), p*sin(t)*sin(h), p*cos(t), m);
    return v_out;
};

void test()
{
    gStyle->SetOptStat(0);

    const std::vector<double> masses = {0.938, 0.139, 0.139};
    const std::vector<TString> parts = {"p","#pi^{+}","#pi^{-}"};
    const std::vector<TString> kines = {"P","#theta","#phi"};
    const std::vector<TString> units = {"GeV","rad","rad"};

    TLorentzVector target(0.0, 0.0, 0.0, 0.938);
    TLorentzVector beam(0.0, 0.0, 10.6, 10.6);
    TLorentzVector W = beam + target;
    TGenPhaseSpace event;
    event.SetDecay(W, masses.size(), &masses[0]);

    //std::random_device rndm_device;
    //std::default_random_engine rndm_engine(rndm_device());

    auto h_chi = new TH1F("h_chi", ";#chi^{2}/ndf", 501, 0, 10);
    auto h_lik = new TH1F("h_lik", ";Confidence Level", 100, 0, 1);
    auto h_mm_gen = new TH1F("h_mm_gen", ";Missing Mass Squared", 501, -0.3, 0.3);
    auto h_mm_sme = (TH1*)h_mm_gen->Clone("h_mm_sme");
    auto h_mm_fit = (TH1*)h_mm_gen->Clone("h_mm_fit");

    std::vector<TH1*> h_pulls,h_fitres,h_smeres;
    for (int i=0; i<parts.size(); i++) {
        for (int j=0; j<kines.size(); j++) {
            h_pulls.push_back(new TH1F(Form("h_pulls_%d_%d",i,j),
                        Form(";%s_{%s} Pull",parts[i].Data(),kines[j].Data()),201,-10,10));
            h_fitres.push_back(new TH1F(Form("h_fitres_%d_%d",i,j),
                        Form(";%s_{%s} Residual [%s]",parts[i].Data(),kines[j].Data(),units[i].Data()),201,-R[j]*10,R[j]*10));
            h_smeres.push_back(new TH1F(Form("h_smeres_%d_%d",i,j),
                        Form(";%s_{%s} Smearing [%s]",parts[i].Data(),kines[j].Data(),units[i].Data()),201,-R[j]*10,R[j]*10));
        }
    }

    int nevents = 0;
    while (nevents < 10000)
    {
        bool reject = false;

        //std::default_random_engine rndm_engine(rndm_device());

        auto weight = event.Generate();

        std::vector<TLorentzVector> parts_gen,parts_sme,parts_fit;
        for (int ipart=0; ipart<parts.size(); ++ipart) {
            parts_gen.push_back(*(event.GetDecay(ipart)));
            parts_sme.push_back(ptr_to_genvector(rndm_engine, event.GetDecay(ipart)));
            if (parts_sme[ipart].Theta()<8*3.14159/180) {
                reject = true;
                break;
            }
        }

        if (reject) continue;

        nevents++;

        auto kin = new KinFitter({{target,beam},parts_sme},{},{R[0],R[1],R[2],R[0],R[1],R[2],R[0],R[1],R[2]});

        auto missing_gen = target+beam - (parts_gen[0]+parts_gen[1]+parts_gen[2]);
        auto missing_sme = target+beam - (parts_sme[0]+parts_sme[1]+parts_sme[2]);
        auto missing_fit = target+beam - (kin->Ps_y[0]+kin->Ps_y[1]+kin->Ps_y[2]);

        h_mm_gen->Fill(missing_gen.M2(), weight);
        h_mm_sme->Fill(missing_sme.M2(), weight);
        h_mm_fit->Fill(missing_fit.M2(), weight);
        h_chi->Fill(kin->chi2/kin->ndf);
        h_lik->Fill(kin->confLevel);

        for (int ipart=0; ipart<parts.size(); ipart++) {
            for (int jkine=0; jkine<kines.size(); jkine++) {
                h_pulls[ipart*kines.size()+jkine]->Fill(kin->pulls[ipart*kines.size()+jkine]);
            }
            h_fitres[ipart*3+0]->Fill(kin->Ps_y[ipart].Vect().Mag()-parts_sme[ipart].Vect().Mag());
            h_fitres[ipart*3+1]->Fill(kin->Ps_y[ipart].Theta()-parts_sme[ipart].Theta());
            h_fitres[ipart*3+2]->Fill(kin->Ps_y[ipart].Phi()-parts_sme[ipart].Phi());
            h_smeres[ipart*3+0]->Fill(parts_gen[ipart].Vect().Mag()-parts_sme[ipart].Vect().Mag());
            h_smeres[ipart*3+1]->Fill(parts_gen[ipart].Theta()-parts_sme[ipart].Theta());
            h_smeres[ipart*3+2]->Fill(parts_gen[ipart].Phi()-parts_sme[ipart].Phi());
        }
    }

    auto Can_missing = new TCanvas("c","Summary",800,800);
    Can_missing->Divide(1,3);
    Can_missing->cd(1);
    gPad->SetLogy();
    h_mm_gen->Draw("hist");
    h_mm_fit->SetLineColor(kGreen);
    h_mm_fit->SetLineWidth(2);
    h_mm_fit->Draw("hist same");
    h_mm_sme->SetLineColor(kRed);
    h_mm_sme->Draw("hist same");
    auto legend = new TLegend(0.54, 0.87, 0.87, 0.67);
    legend->AddEntry(h_mm_gen, "Gen", "l");
    legend->AddEntry(h_mm_sme, "Smeared", "l");
    legend->AddEntry(h_mm_fit, "Fitted", "l");
    legend->SetFillStyle(0);
    legend->SetLineWidth(0);
    legend->Draw("same ");      
    Can_missing->cd(3);
    h_chi->Draw();
    Can_missing->cd(2);
    gPad->SetLogy();
    h_lik->Draw();

    Can_missing->SaveAs("Can_Missing.pdf");

    TString fitopt="Q";
    auto Can_pulls = new TCanvas("cc","Pulls",900,600);
    Can_pulls->Divide(3,3);
    for (int ipart=0; ipart<parts.size(); ipart++) {
        for (int jkine=0; jkine<kines.size(); jkine++) {
            Can_pulls->cd(ipart*kines.size()+jkine+1);
            h_pulls[ipart*kines.size()+jkine]->Fit("gaus",fitopt);
        }
    }
    Can_pulls->SaveAs("Can_pulls.pdf");     

    auto Can_res = new TCanvas("ccc","Residuals",900,600);
    Can_res->Divide(3,3);
    for (int ipart=0; ipart<parts.size(); ipart++) {
        for (int jkine=0; jkine<kines.size(); jkine++) {
            Can_res->cd(ipart*kines.size()+jkine+1);
            h_fitres[ipart*kines.size()+jkine]->Fit("gaus",fitopt);
        }
    }
    Can_pulls->SaveAs("Can_res.pdf");       

    auto Can_sme = new TCanvas("cccc","Smearing",900,600);
    Can_sme->Divide(3,3);
    for (int ipart=0; ipart<parts.size(); ipart++) {
        for (int jkine=0; jkine<kines.size(); jkine++) {
            Can_sme->cd(ipart*kines.size()+jkine+1);
            h_smeres[ipart*kines.size()+jkine]->Fit("gaus",fitopt);
        }
    }
    Can_sme->SaveAs("Can_sme.pdf");     
}
