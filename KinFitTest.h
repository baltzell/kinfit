#ifndef KinFitTest_h
#define KinFitTest_h

#include <vector>
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TRandom3.h"
#include "KinFitter.h"
#include <iostream>

// absolute resolutions (GeV,radians,radians):
const std::vector<double> RESO = {0.15, 0.02, 0.02};
const std::vector<TString> KINES = {"P", "#theta", "#phi"};
const std::vector<TString> UNITS = {"GeV", "rad", "rad"};

TRandom3 RNDM3(0);

TLorentzVector smear(TLorentzVector *v)
{
    // smear by absolute resolutions:
    const double p = RNDM3.Gaus(0.0, RESO[0]) + v->P();
    const double t = RNDM3.Gaus(0.0, RESO[1]) + v->Theta();
    const double h = RNDM3.Gaus(0.0, RESO[2]) + v->Phi();

    TLorentzVector v_out;
    v_out.SetXYZM(p * sin(t) * cos(h), p * sin(t) * sin(h), p * cos(t), v->M());
    return v_out;
};

class KinFitTest
{
public:
    virtual ~KinFitTest() {}

    TString _name;
    TLorentzVector _W;
    std::vector<TString> _parts;
    int _missing;

    TH1 *_h_chi;
    TH1 *_h_lik_Signal;
    TH1 *_h_lik_BG;

    TH1 *_h_mm_gen;
    TH1 *_h_mm_sme;
    TH1 *_h_mm_fit;

    TH1 *_h_inv_gen;
    TH1 *_h_inv_sme;
    TH1 *_h_inv_fit;

    TH1 *_h_E_gen;
    TH1 *_h_E_sme;
    TH1 *_h_E_fit;

    TH1 *_h_Px_gen;
    TH1 *_h_Px_sme;
    TH1 *_h_Px_fit;

    TH1 *_h_Py_gen;
    TH1 *_h_Py_sme;
    TH1 *_h_Py_fit;

    TH1 *_h_Pz_gen;
    TH1 *_h_Pz_sme;
    TH1 *_h_Pz_fit;

    std::vector<TH1 *> _h_pulls;
    std::vector<TH1 *> _h_fitres;
    std::vector<TH1 *> _h_smeres;
    std::vector<TH1 *> _h_fitgen;

    KinFitTest(TString name, std::vector<TString> parts, TLorentzVector W, int missing = 0, float missmass = 0, float invmass = 0)
        : _name(name),
          _parts(parts),
          _W(W),
          _missing(missing)
    {
        _h_chi = new TH1F("h_chi", ";#chi^{2}/ndf", 201, 0, 10);
        _h_lik_Signal = new TH1F("h_lik_Signal", ";Confidence Level", 100, 0, 1);
        _h_lik_BG = new TH1F("h_lik_BG", ";Confidence Level", 100, 0, 1);

        _h_mm_gen = new TH1F("h_mm_gen", ";Missing Mass [GeV^{2}]", 201, missmass - 0.5 * (parts.size() - 1), missmass + 0.5 * (parts.size() - 1));
        _h_mm_sme = (TH1 *)_h_mm_gen->Clone("h_mm_sme");
        _h_mm_fit = (TH1 *)_h_mm_gen->Clone("h_mm_fit");

        _h_inv_gen = new TH1F("h_inv_gen", ";Invariant Mass [GeV]", 201, invmass - 0.5 * (parts.size() - 1), invmass + 0.5 * (parts.size() - 1));
        _h_inv_sme = (TH1 *)_h_inv_gen->Clone("h_inv_sme");
        _h_inv_fit = (TH1 *)_h_inv_gen->Clone("h_inv_sme");

        _h_E_gen = new TH1F("h_E_gen", ";#Delta E [GeV]", 201, -0.5 * (parts.size() - 1), 0.5 * (parts.size() - 1));
        _h_E_sme = (TH1 *)_h_E_gen->Clone("h_E_sme");
        _h_E_fit = (TH1 *)_h_E_gen->Clone("h_E_fit");

        _h_Px_gen = new TH1F("h_Px_gen", ";#Delta Px [GeV]", 201, -0.5 * (parts.size() - 1), 0.5 * (parts.size() - 1));
        _h_Px_sme = (TH1 *)_h_Px_gen->Clone("h_Px_sme");
        _h_Px_fit = (TH1 *)_h_Px_gen->Clone("h_Px_fit");

        _h_Py_gen = new TH1F("h_Py_gen", ";#Delta Py [GeV]", 201, -0.5 * (parts.size() - 1), 0.5 * (parts.size() - 1));
        _h_Py_sme = (TH1 *)_h_Py_gen->Clone("h_Py_sme");
        _h_Py_fit = (TH1 *)_h_Py_gen->Clone("h_Py_fit");

        _h_Pz_gen = new TH1F("h_Pz_gen", ";#Delta Pz [GeV]", 201, -0.5 * (parts.size() - 1), 0.5 * (parts.size() - 1));
        _h_Pz_sme = (TH1 *)_h_Pz_gen->Clone("h_Pz_sme");
        _h_Pz_fit = (TH1 *)_h_Pz_gen->Clone("h_Pz_fit");

        for (int i = 0; i < parts.size() - _missing; i++)
        {
            for (int j = 0; j < KINES.size(); j++)
            {
                _h_pulls.push_back(new TH1F(Form("h_pulls_%d_%d", i, j),
                                            Form(";%s_{%s} Pull", parts[i].Data(), KINES[j].Data()), 201, -10, 10));
                _h_fitres.push_back(new TH1F(Form("h_fitres_%d_%d", i, j),
                                             Form(";%s_{%s} Residual [%s]", parts[i].Data(), KINES[j].Data(), UNITS[i].Data()), 201, -RESO[j] * 10, RESO[j] * 10));
                _h_smeres.push_back(new TH1F(Form("h_smeres_%d_%d", i, j),
                                             Form(";%s_{%s} Smearing [%s]", parts[i].Data(), KINES[j].Data(), UNITS[i].Data()), 201, -RESO[j] * 10, RESO[j] * 10));
                _h_fitgen.push_back(new TH1F(Form("h_fitgen%d_%d", i, j),
                                             Form(";%s_{%s} Fit/Gen [%s]", parts[i].Data(), KINES[j].Data(), UNITS[i].Data()), 201, -RESO[j] * 10, RESO[j] * 10));
            }
        }
    }

    void fill_MissingMass(KinFitter *kin, std::vector<TLorentzVector> parts_gen, std::vector<TLorentzVector> parts_sme, double weight, bool background)
    {
        std::vector<TLorentzVector> parts_fit = kin->GetFitted4Vectors();

        TLorentzVector missing_gen = _W;
        TLorentzVector missing_sme = _W;
        TLorentzVector missing_fit = _W;

        for (int ipart = 0; ipart < _parts.size() - _missing; ipart++)
        {
            missing_gen -= parts_gen[ipart];
            missing_sme -= parts_sme[ipart];
            missing_fit -= parts_fit[ipart];

            for (int jkine = 0; jkine < KINES.size() && kin->GetConfidenceLevel() > 0.0001; jkine++)
            {
                _h_pulls[ipart * KINES.size() + jkine]->Fill(kin->GetPulls()[ipart * KINES.size() + jkine]);
            }
            _h_fitres[ipart * 3 + 0]->Fill(parts_fit[ipart].Vect().Mag() - parts_sme[ipart].Vect().Mag());
            _h_fitres[ipart * 3 + 1]->Fill(parts_fit[ipart].Theta() - parts_sme[ipart].Theta());
            _h_fitres[ipart * 3 + 2]->Fill(parts_fit[ipart].Phi() - parts_sme[ipart].Phi());

            _h_smeres[ipart * 3 + 0]->Fill(parts_gen[ipart].Vect().Mag() - parts_sme[ipart].Vect().Mag());
            _h_smeres[ipart * 3 + 1]->Fill(parts_gen[ipart].Theta() - parts_sme[ipart].Theta());
            _h_smeres[ipart * 3 + 2]->Fill(parts_gen[ipart].Phi() - parts_sme[ipart].Phi());

            _h_fitgen[ipart * 3 + 0]->Fill(parts_gen[ipart].Vect().Mag() - parts_fit[ipart].Vect().Mag());
            _h_fitgen[ipart * 3 + 1]->Fill(parts_gen[ipart].Theta() - parts_fit[ipart].Theta());
            _h_fitgen[ipart * 3 + 2]->Fill(parts_gen[ipart].Phi() - parts_fit[ipart].Phi());
        }

        _h_mm_gen->Fill(missing_gen.M2(), weight);
        _h_mm_sme->Fill(missing_sme.M2(), weight);
        _h_mm_fit->Fill(missing_fit.M2(), weight);

        _h_E_gen->Fill(missing_gen.E(), weight);
        _h_E_sme->Fill(missing_sme.E(), weight);
        _h_E_fit->Fill(missing_fit.E(), weight);

        _h_Px_gen->Fill(missing_gen.Px(), weight);
        _h_Px_sme->Fill(missing_sme.Px(), weight);
        _h_Px_fit->Fill(missing_fit.Px(), weight);

        _h_Py_gen->Fill(missing_gen.Py(), weight);
        _h_Py_sme->Fill(missing_sme.Py(), weight);
        _h_Py_fit->Fill(missing_fit.Py(), weight);

        _h_Pz_gen->Fill(missing_gen.Pz(), weight);
        _h_Pz_sme->Fill(missing_sme.Pz(), weight);
        _h_Pz_fit->Fill(missing_fit.Pz(), weight);

        _h_chi->Fill(kin->GetChi2() / kin->GetNDF());

        if (!background)
            _h_lik_Signal->Fill(kin->GetConfidenceLevel());
        else
            _h_lik_BG->Fill(kin->GetConfidenceLevel());
    }

    void fill_InvariantMass(KinFitter *kin, std::vector<TLorentzVector> parts_gen, std::vector<TLorentzVector> parts_sme, std::vector<int> indices_part, double weight, bool is_background)
    // std::tuple<std::vector<std::vector<double>>, double> fill_InvariantMass(KinFitter *kin, std::vector<TLorentzVector> parts_gen, std::vector<TLorentzVector> parts_sme, std::vector<int> indices_part, double weight, bool is_background)
    {
        int parts_size = _parts.size();
        int kines_size = KINES.size();
        // std::vector<std::vector<double>> pull_vals(parts_size, std::vector<double>(kines_size));

        std::vector<TLorentzVector> parts_fit = kin->GetFitted4Vectors();

        TLorentzVector invariant_gen;
        invariant_gen.SetXYZT(0, 0, 0, 0);
        TLorentzVector invariant_sme;
        invariant_sme.SetXYZT(0, 0, 0, 0);
        TLorentzVector invariant_fit;
        invariant_fit.SetXYZT(0, 0, 0, 0);

        if (kin->HasConverged())
        {
            for (int ipart = 0; ipart < indices_part.size(); ipart++)
            {
                invariant_gen += parts_gen[ipart];
                invariant_sme += parts_sme[ipart];
                invariant_fit += parts_fit[ipart];

                for (int jkine = 0; jkine < KINES.size(); jkine++)
                {
                    if (kin->GetConfidenceLevel() > 0.0001)
                    {
                        _h_pulls[ipart * KINES.size() + jkine]->Fill(kin->GetPulls()[ipart * KINES.size() + jkine]);
                    }
                    // pull_vals[ipart][jkine] = kin->GetPulls()[ipart * KINES.size() + jkine];
                }
                _h_fitres[ipart * 3 + 0]->Fill(parts_fit[ipart].Vect().Mag() - parts_sme[ipart].Vect().Mag());
                _h_fitres[ipart * 3 + 1]->Fill(parts_fit[ipart].Theta() - parts_sme[ipart].Theta());
                _h_fitres[ipart * 3 + 2]->Fill(parts_fit[ipart].Phi() - parts_sme[ipart].Phi());

                _h_smeres[ipart * 3 + 0]->Fill(parts_gen[ipart].Vect().Mag() - parts_sme[ipart].Vect().Mag());
                _h_smeres[ipart * 3 + 1]->Fill(parts_gen[ipart].Theta() - parts_sme[ipart].Theta());
                _h_smeres[ipart * 3 + 2]->Fill(parts_gen[ipart].Phi() - parts_sme[ipart].Phi());

                _h_fitgen[ipart * 3 + 0]->Fill(parts_gen[ipart].Vect().Mag() - parts_fit[ipart].Vect().Mag());
                _h_fitgen[ipart * 3 + 1]->Fill(parts_gen[ipart].Theta() - parts_fit[ipart].Theta());
                _h_fitgen[ipart * 3 + 2]->Fill(parts_gen[ipart].Phi() - parts_fit[ipart].Phi());
            }
            _h_inv_gen->Fill(invariant_gen.M(), weight);
            _h_inv_sme->Fill(invariant_sme.M(), weight);
            _h_inv_fit->Fill(invariant_fit.M(), weight);

            _h_chi->Fill(kin->GetChi2() / kin->GetNDF());

            if (!is_background)
                _h_lik_Signal->Fill(kin->GetConfidenceLevel());
            else
                _h_lik_BG->Fill(kin->GetConfidenceLevel());
        }

        // return std::make_tuple(pull_vals, kin->GetConfidenceLevel());
        return;
    }

    void fill_3C(KinFitter *kin, std::vector<TLorentzVector> parts_gen, std::vector<TLorentzVector> parts_sme, double weight, bool background)
    {
        std::vector<TLorentzVector> parts_fit = kin->GetFitted4Vectors();

        int ipart = 0;
        TLorentzVector missing_gen;
        missing_gen.SetXYZT(0, 0, 0, 0);
        TLorentzVector missing_sme;
        missing_sme.SetXYZT(0, 0, 0, 0);
        TLorentzVector missing_fit;
        missing_fit.SetXYZT(0, 0, 0, 0);

        missing_gen = parts_gen[ipart]-parts_gen[ipart];
        missing_sme = parts_sme[ipart]-parts_gen[ipart];
        missing_fit = parts_fit[ipart]-parts_gen[ipart];

        for (int jkine = 0; jkine < KINES.size() && kin->GetConfidenceLevel() > 0.0001; jkine++)
        {
            _h_pulls[ipart * KINES.size() + jkine]->Fill(kin->GetPulls()[ipart * KINES.size() + jkine]);
        }
        _h_fitres[ipart * 3 + 0]->Fill(parts_fit[ipart].Vect().Mag() - parts_sme[ipart].Vect().Mag());
        _h_fitres[ipart * 3 + 1]->Fill(parts_fit[ipart].Theta() - parts_sme[ipart].Theta());
        _h_fitres[ipart * 3 + 2]->Fill(parts_fit[ipart].Phi() - parts_sme[ipart].Phi());

        _h_smeres[ipart * 3 + 0]->Fill(parts_gen[ipart].Vect().Mag() - parts_sme[ipart].Vect().Mag());
        _h_smeres[ipart * 3 + 1]->Fill(parts_gen[ipart].Theta() - parts_sme[ipart].Theta());
        _h_smeres[ipart * 3 + 2]->Fill(parts_gen[ipart].Phi() - parts_sme[ipart].Phi());

        _h_fitgen[ipart * 3 + 0]->Fill(parts_gen[ipart].Vect().Mag() - parts_fit[ipart].Vect().Mag());
        _h_fitgen[ipart * 3 + 1]->Fill(parts_gen[ipart].Theta() - parts_fit[ipart].Theta());
        _h_fitgen[ipart * 3 + 2]->Fill(parts_gen[ipart].Phi() - parts_fit[ipart].Phi());

        _h_mm_gen->Fill(missing_gen.M2(), weight);
        _h_mm_sme->Fill(missing_sme.M2(), weight);
        _h_mm_fit->Fill(missing_fit.M2(), weight);

        _h_E_gen->Fill(missing_gen.E(), weight);
        _h_E_sme->Fill(missing_sme.E(), weight);
        _h_E_fit->Fill(missing_fit.E(), weight);

        _h_Px_gen->Fill(missing_gen.Px(), weight);
        _h_Px_sme->Fill(missing_sme.Px(), weight);
        _h_Px_fit->Fill(missing_fit.Px(), weight);

        _h_Py_gen->Fill(missing_gen.Py(), weight);
        _h_Py_sme->Fill(missing_sme.Py(), weight);
        _h_Py_fit->Fill(missing_fit.Py(), weight);

        _h_Pz_gen->Fill(missing_gen.Pz(), weight);
        _h_Pz_sme->Fill(missing_sme.Pz(), weight);
        _h_Pz_fit->Fill(missing_fit.Pz(), weight);

        _h_chi->Fill(kin->GetChi2() / kin->GetNDF());

        if (!background)
            _h_lik_Signal->Fill(kin->GetConfidenceLevel());
        else
            _h_lik_BG->Fill(kin->GetConfidenceLevel());
    }

    void plot()
    {
        gStyle->SetOptStat(0);
        gStyle->SetOptFit(1111);

        TCanvas *c_missing = new TCanvas("can_mm", "Summary", 800, 1800);
        c_missing->Divide(1, 3);
        c_missing->cd(1);
        gPad->SetLogy();
        _h_mm_gen->GetYaxis()->SetRangeUser(1., _h_mm_gen->GetMaximum());
        _h_mm_gen->SetFillStyle(1001);
        _h_mm_gen->SetFillColor(kBlue);
        _h_mm_gen->Draw("hist");
        _h_mm_fit->SetLineColor(kGreen);
        _h_mm_fit->SetLineWidth(2);
        _h_mm_fit->Draw("hist same");
        _h_mm_sme->SetLineColor(kRed);
        _h_mm_sme->Draw("hist same");
        auto legend = new TLegend(0.54, 0.87, 0.87, 0.67);
        legend->AddEntry(_h_mm_gen, Form("Gen (%f)", _h_mm_gen->GetEffectiveEntries()), "l");
        legend->AddEntry(_h_mm_sme, Form("Smeared (%f)", _h_mm_sme->GetEffectiveEntries()), "l");
        legend->AddEntry(_h_mm_fit, Form("Fitted (%f)", _h_mm_fit->GetEffectiveEntries()), "l");
        legend->SetFillStyle(0);
        legend->SetLineWidth(0);
        legend->Draw("same ");
        c_missing->cd(3);
        _h_chi->Draw();
        c_missing->cd(2);
        gPad->SetLogy();
        _h_lik_Signal->GetYaxis()->SetRangeUser(1., _h_lik_Signal->GetMaximum());
        _h_lik_Signal->SetLineColor(kRed);
        _h_lik_Signal->Draw();
        _h_lik_BG->SetLineColor(kBlue);
        _h_lik_BG->Draw("same");
        c_missing->Print(Form("%s.pdf(", _name.Data()));

        auto c_invariant = new TCanvas("can_inv", "Summary", 800, 800);
        gPad->SetLogy();
        _h_inv_gen->GetYaxis()->SetRangeUser(1., _h_inv_gen->GetMaximum());
        _h_inv_gen->SetFillStyle(1001);
        _h_inv_gen->SetFillColor(kBlue);
        _h_inv_gen->Draw("hist");
        _h_inv_fit->SetLineColor(kGreen);
        _h_inv_fit->SetLineWidth(2);
        _h_inv_fit->Draw("hist same");
        _h_inv_sme->SetLineColor(kRed);
        _h_inv_sme->Draw("hist same");
        auto legend_inv = new TLegend(0.54, 0.87, 0.87, 0.67);
        legend_inv->AddEntry(_h_inv_gen, Form("Gen (%f)", _h_inv_gen->GetEffectiveEntries()), "l");
        legend_inv->AddEntry(_h_inv_sme, Form("Smeared (%f)", _h_inv_sme->GetEffectiveEntries()), "l");
        legend_inv->AddEntry(_h_inv_fit, Form("Fitted (%f)", _h_inv_fit->GetEffectiveEntries()), "l");
        legend_inv->SetFillStyle(0);
        legend_inv->SetLineWidth(0);
        legend_inv->Draw("same ");
        c_invariant->Print(Form("%s.pdf(", _name.Data()));

        if (_missing == 0)
        {
            auto c_constraint = new TCanvas("c_constraint", "Constraints", 800, 1200);
            c_constraint->Divide(1, 4);
            c_constraint->cd(1);
            gPad->SetLogy();
            _h_E_gen->Draw("hist");
            _h_E_fit->SetLineColor(kGreen);
            _h_E_fit->SetLineWidth(2);
            _h_E_fit->Draw("hist same");
            _h_E_sme->SetLineColor(kRed);
            _h_E_sme->Draw("hist same");
            legend->Draw("same ");
            c_constraint->cd(2);
            gPad->SetLogy();
            _h_Px_gen->Draw("hist");
            _h_Px_fit->SetLineColor(kGreen);
            _h_Px_fit->SetLineWidth(2);
            _h_Px_fit->Draw("hist same");
            _h_Px_sme->SetLineColor(kRed);
            _h_Px_sme->Draw("hist same");
            c_constraint->cd(3);
            gPad->SetLogy();
            _h_Py_gen->Draw("hist");
            _h_Py_fit->SetLineColor(kGreen);
            _h_Py_fit->SetLineWidth(2);
            _h_Py_fit->Draw("hist same");
            _h_Py_sme->SetLineColor(kRed);
            _h_Py_sme->Draw("hist same");
            c_constraint->cd(4);
            gPad->SetLogy();
            _h_Pz_gen->Draw("hist");
            _h_Pz_fit->SetLineColor(kGreen);
            _h_Pz_fit->SetLineWidth(2);
            _h_Pz_fit->Draw("hist same");
            _h_Pz_sme->SetLineColor(kRed);
            _h_Pz_sme->Draw("hist same");
            c_constraint->SaveAs(Form("%s.pdf", _name.Data()));
        }

        TString fitopt = "Q";
        auto c_pulls = new TCanvas("can2", "Pulls", 900, int(float(1200) * (_parts.size() - _missing) / 3));
        c_pulls->Divide(3, _parts.size() - _missing);
        for (int ipart = 0; ipart < _parts.size() - _missing; ipart++)
        {
            for (int jkine = 0; jkine < KINES.size(); jkine++)
            {
                c_pulls->cd(ipart * KINES.size() + jkine + 1);
                _h_pulls[ipart * KINES.size() + jkine]->Fit("gaus", fitopt);
            }
        }
        c_pulls->SaveAs(Form("%s.pdf", _name.Data()));

        auto c_res = new TCanvas("can3", "Residuals", 900, int(float(1200) * (_parts.size() - _missing) / 3));
        c_res->Divide(3, _parts.size() - _missing);
        for (int ipart = 0; ipart < _parts.size() - _missing; ipart++)
        {
            for (int jkine = 0; jkine < KINES.size(); jkine++)
            {
                c_res->cd(ipart * KINES.size() + jkine + 1);
                _h_fitres[ipart * KINES.size() + jkine]->Fit("gaus", fitopt);
            }
        }
        c_res->SaveAs(Form("%s.pdf", _name.Data()));

        auto c_sme = new TCanvas("can4", "Smearing", 900, int(float(1200) * (_parts.size() - _missing) / 3));
        c_sme->Divide(3, _parts.size() - _missing);
        for (int ipart = 0; ipart < _parts.size() - _missing; ipart++)
        {
            for (int jkine = 0; jkine < KINES.size(); jkine++)
            {
                c_sme->cd(ipart * KINES.size() + jkine + 1);
                _h_smeres[ipart * KINES.size() + jkine]->Fit("gaus", fitopt);
                //_h_fitgen[ipart * KINES.size() + jkine]->SetLineColor(kGreen);
                //_h_fitgen[ipart * KINES.size() + jkine]->Draw("same");
            }
        }
        c_sme->SaveAs(Form("%s.pdf)", _name.Data()));
    }
};

#endif
