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
#include "KinParticle.h"
#include "TRandom3.h"

#include "test.h"

int test_InvMass(const int max_events=10000, const float bg_fraction=0.1)
{
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1111);

    const std::vector<double> masses = {0.150, 0.150};
    const std::vector<TString> parts = {"e", "e"};

    TLorentzVector JPsi;
    JPsi.SetXYZM(0.0, 0.0, 5., 3.);
    TGenPhaseSpace event;
    event.SetDecay(JPsi, masses.size(), &masses[0]);

    auto h_chi = new TH1F("h_chi", ";#chi^{2}/ndf", 101, 0, 10);
    auto h_lik_Signal = new TH1F("h_lik_Signal", ";Confidence Level", 100, 0, 1);
    auto h_lik_BG = new TH1F("h_lik_BG", ";Confidence Level", 100, 0, 1);

    auto h_mm_gen = new TH1F("h_mm_gen", ";Mass [GeV]", 101, 2.5, 3.5);
    auto h_mm_sme = (TH1 *)h_mm_gen->Clone("h_mm_sme");
    auto h_mm_fit = (TH1 *)h_mm_gen->Clone("h_mm_fit");

    std::vector<TH1 *> h_pulls, h_fitres, h_smeres, h_fitgen;
    std::vector<double> resolutions;
    std::vector<int> constraint_idx;
    for (int i = 0; i < parts.size(); i++)
    {
        constraint_idx.push_back(i);
        for (int j = 0; j < KINES.size(); j++)
        {
            resolutions.push_back(RESO[j]);
            h_pulls.push_back(new TH1F(Form("h_pulls_%d_%d", i, j),
                                       Form(";%s_{%s} Pull", parts[i].Data(), KINES[j].Data()), 201, -10, 10));
            h_fitres.push_back(new TH1F(Form("h_fitres_%d_%d", i, j),
                                        Form(";%s_{%s} Residual [%s]", parts[i].Data(), KINES[j].Data(), UNITS[i].Data()), 201, -RESO[j] * 10, RESO[j] * 10));
            h_smeres.push_back(new TH1F(Form("h_smeres_%d_%d", i, j),
                                        Form(";%s_{%s} Smearing [%s]", parts[i].Data(), KINES[j].Data(), UNITS[i].Data()), 201, -RESO[j] * 10, RESO[j] * 10));
            h_fitgen.push_back(new TH1F(Form("h_fitgen%d_%d", i, j),
                                        Form(";%s_{%s} Fit/Gen [%s]", parts[i].Data(), KINES[j].Data(), UNITS[i].Data()), 201, -RESO[j] * 10, RESO[j] * 10));
        }
    }

    bool is_BG = false;
    TLorentzVector BG;

    int nevents = 0;
    while (nevents < max_events)
    {

        if (RNDM3.Uniform(0.0, 1.0) < bg_fraction)
        {
            float mass = RNDM3.Uniform(2.5, 3.5);
            BG.SetXYZM(0.0, 0.0, 5., mass);
            event.SetDecay(BG, masses.size(), &masses[0]);
            is_BG = true;
        }
        else
        {
            event.SetDecay(JPsi, masses.size(), &masses[0]);
            is_BG = false;
        }

        auto weight = event.Generate();

        std::vector<TLorentzVector> parts_gen;
        std::vector<TLorentzVector> parts_sme;

        std::vector<KinParticle> kin_parts_sme;

        for (int ipart = 0; ipart < parts.size(); ++ipart)
        {
            parts_gen.push_back(*(event.GetDecay(ipart)));
            TLorentzVector sme_vector = smear(event.GetDecay(ipart));
            parts_sme.push_back(sme_vector);

            kin_parts_sme.push_back(KinParticle(sme_vector, masses[ipart], RESO));
        }

        nevents++;

        auto kin = new KinFitter({KinParticle(JPsi)}, kin_parts_sme);
        kin->Add_InvMass_Constraint(constraint_idx, 3.0);
        kin->DoFitting(100);

        std::vector<TLorentzVector> parts_fit = kin->GetFitted4Vectors();

        TLorentzVector JPsi_gen(0, 0, 0, 0);
        TLorentzVector JPsi_sme(0, 0, 0, 0);
        TLorentzVector JPsi_fit(0, 0, 0, 0);
        for (int ipart = 0; ipart < parts.size(); ipart++)
        {
            JPsi_gen += parts_gen[ipart];
            JPsi_sme += parts_sme[ipart];
            JPsi_fit += parts_fit[ipart];
        }

        h_mm_gen->Fill(JPsi_gen.M(), weight);
        h_mm_sme->Fill(JPsi_sme.M(), weight);
        h_mm_fit->Fill(JPsi_fit.M(), weight);
        h_chi->Fill(kin->GetChi2() / kin->GetNDF());

         if (is_BG)
            h_lik_BG->Fill(kin->GetConfidenceLevel());
        else
            h_lik_Signal->Fill(kin->GetConfidenceLevel());


        for (int ipart = 0; ipart < parts.size(); ipart++)
        {
            for (int jkine = 0; jkine < KINES.size(); jkine++)
            {
                h_pulls[ipart * KINES.size() + jkine]->Fill(kin->GetPulls()[ipart * KINES.size() + jkine]);
            }
            h_fitres[ipart * 3 + 0]->Fill(parts_fit[ipart].Vect().Mag() - parts_sme[ipart].Vect().Mag());
            h_fitres[ipart * 3 + 1]->Fill(parts_fit[ipart].Theta() - parts_sme[ipart].Theta());
            h_fitres[ipart * 3 + 2]->Fill(parts_fit[ipart].Phi() - parts_sme[ipart].Phi());

            h_smeres[ipart * 3 + 0]->Fill(parts_gen[ipart].Vect().Mag() - parts_sme[ipart].Vect().Mag());
            h_smeres[ipart * 3 + 1]->Fill(parts_gen[ipart].Theta() - parts_sme[ipart].Theta());
            h_smeres[ipart * 3 + 2]->Fill(parts_gen[ipart].Phi() - parts_sme[ipart].Phi());

            h_fitgen[ipart * 3 + 0]->Fill(parts_gen[ipart].Vect().Mag() - parts_fit[ipart].Vect().Mag());
            h_fitgen[ipart * 3 + 1]->Fill(parts_gen[ipart].Theta() - parts_fit[ipart].Theta());
            h_fitgen[ipart * 3 + 2]->Fill(parts_gen[ipart].Phi() - parts_fit[ipart].Phi());
        }
    }

    auto c_missing = new TCanvas("can1", "Summary", 800, 1500);
    c_missing->Divide(1, 3);
    c_missing->cd(1);
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
    c_missing->cd(3);
    h_chi->Draw();
    c_missing->cd(2);
    gPad->SetLogy();
    h_lik_Signal->SetLineColor(kRed);
    h_lik_Signal->Draw();
    h_lik_BG->SetLineColor(kBlue);
    h_lik_BG->Draw("same");
    auto legend1 = new TLegend(0.54, 0.87, 0.87, 0.67);
    legend1->AddEntry(h_lik_Signal, "Signal", "l");
    legend1->AddEntry(h_lik_BG, "BG", "l");
    legend1->SetFillStyle(0);
    legend1->SetLineWidth(0);
    legend1->Draw("same ");
    c_missing->SaveAs("c_Missing_inv.pdf");

    TString fitopt = "Q";
    auto c_pulls = new TCanvas("can2", "Pulls", 900, int(float(1200) * parts.size() / 3));
    c_pulls->Divide(3, parts.size());
    for (int ipart = 0; ipart < parts.size(); ipart++)
    {
        for (int jkine = 0; jkine < KINES.size(); jkine++)
        {
            c_pulls->cd(ipart * KINES.size() + jkine + 1);
            h_pulls[ipart * KINES.size() + jkine]->Fit("gaus", fitopt);
        }
    }
    c_pulls->SaveAs("c_pulls_inv.pdf");

    auto c_res = new TCanvas("can3", "Residuals", 900, 600);
    c_res->Divide(3, 3);
    for (int ipart = 0; ipart < parts.size(); ipart++)
    {
        for (int jkine = 0; jkine < KINES.size(); jkine++)
        {
            c_res->cd(ipart * KINES.size() + jkine + 1);
            h_fitres[ipart * KINES.size() + jkine]->Fit("gaus", fitopt);
        }
    }
    c_pulls->SaveAs("c_res_inv.pdf");

    auto c_sme = new TCanvas("can4", "Smearing", 900, 600);
    c_sme->Divide(3, 3);
    for (int ipart = 0; ipart < parts.size(); ipart++)
    {
        for (int jkine = 0; jkine < KINES.size(); jkine++)
        {
            c_sme->cd(ipart * KINES.size() + jkine + 1);
            h_smeres[ipart * KINES.size() + jkine]->Fit("gaus", fitopt);
            h_fitgen[ipart * KINES.size() + jkine]->SetLineColor(kGreen);
            h_fitgen[ipart * KINES.size() + jkine]->Draw("same");
        }
    }
    c_sme->SaveAs("c_sme_inv.pdf");

    return 0;
}

int main() {
    return test_InvMass();
}

