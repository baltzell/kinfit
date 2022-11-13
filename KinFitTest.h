#ifndef KinFitTest_h
#define KinFitTest_h

#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TRandom3.h"

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
    virtual ~KinFitTest(){}

    TString _name;
    TLorentzVector _W;
    std::vector<TString> _parts;
    int _missing;

    TH1* _h_chi;
    TH1* _h_lik_Signal;
    TH1* _h_lik_BG;

    TH1* _h_mm_gen;
    TH1* _h_mm_sme;
    TH1* _h_mm_fit;

    TH1* _h_E_gen;
    TH1* _h_E_sme;
    TH1* _h_E_fit;

    TH1* _h_Px_gen;
    TH1* _h_Px_sme;
    TH1* _h_Px_fit;

    TH1* _h_Py_gen;
    TH1* _h_Py_sme;
    TH1* _h_Py_fit;

    TH1* _h_Pz_gen;
    TH1* _h_Pz_sme;
    TH1* _h_Pz_fit;
    
    std::vector<TH1*> _h_pulls;
    std::vector<TH1*> _h_fitres;
    std::vector<TH1*> _h_smeres;
    std::vector<TH1*> _h_fitgen;

    KinFitTest(TString name, std::vector<TString> parts, TLorentzVector W, int missing=0)
        :
        _name(name),
        _parts(parts),
        _W(W),
        _missing(missing)
    {
        _h_chi = new TH1F("h_chi", ";#chi^{2}/ndf", 201, 0, 10);
        _h_lik_Signal = new TH1F("h_lik_Signal", ";Confidence Level", 100, 0, 1);
        _h_lik_BG = new TH1F("h_lik_BG", ";Confidence Level", 100, 0, 1);

        _h_mm_gen = new TH1F("h_mm_gen", ";Missing Mass [GeV^{2}]", 201, -0.5 * (parts.size() - 1), 0.5 * (parts.size() - 1));
        _h_mm_sme = (TH1 *)_h_mm_gen->Clone("h_mm_sme");
        _h_mm_fit = (TH1 *)_h_mm_gen->Clone("h_mm_fit");

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

        for (int i = 0; i < parts.size()-_missing; i++)
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

    void fill(KinFitter *kin, std::vector<TLorentzVector> parts_gen, std::vector<TLorentzVector> parts_sme, double weight, bool background)
    {
        std::vector<TLorentzVector> parts_fit = kin->GetFitted4Vectors();
        
        TLorentzVector missing_gen = _W;
        TLorentzVector missing_sme = _W;
        TLorentzVector missing_fit = _W;

        for (int ipart = 0; ipart < _parts.size()-_missing; ipart++)
        {
            missing_gen -= parts_gen[ipart];
            missing_sme -= parts_sme[ipart];
            missing_fit -= parts_fit[ipart];
            
            for (int jkine = 0; jkine < KINES.size(); jkine++)
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
        
        if (background)
            _h_lik_Signal->Fill(kin->GetConfidenceLevel());
        else
            _h_lik_BG->Fill(kin->GetConfidenceLevel());
        

    }


    void plot()
    {
        gStyle->SetOptStat(0);
        gStyle->SetOptFit(1111);

        auto c_missing = new TCanvas("can1", "Summary", 800, 1500);
        c_missing->Divide(1, 3);
        c_missing->cd(1);
        gPad->SetLogy();
        _h_mm_gen->GetYaxis()->SetRangeUser(1., _h_mm_gen->GetMaximum());
        _h_mm_gen->Draw("hist");
        _h_mm_fit->SetLineColor(kGreen);
        _h_mm_fit->SetLineWidth(2);
        _h_mm_fit->Draw("hist same");
        _h_mm_sme->SetLineColor(kRed);
        _h_mm_sme->Draw("hist same");
        auto legend = new TLegend(0.54, 0.87, 0.87, 0.67);
        legend->AddEntry(_h_mm_gen, "Gen", "l");
        legend->AddEntry(_h_mm_sme, "Smeared", "l");
        legend->AddEntry(_h_mm_fit, "Fitted", "l");
        legend->SetFillStyle(0);
        legend->SetLineWidth(0);
        legend->Draw("same ");
        c_missing->cd(3);
        _h_chi->Draw();
        c_missing->cd(2);
        gPad->SetLogy();
        _h_lik_Signal->SetLineColor(kRed);
        _h_lik_Signal->Draw();
        _h_lik_BG->SetLineColor(kBlue);
        _h_lik_BG->Draw("same");
        c_missing->SaveAs(Form("missing_%s.pdf",_name.Data()));

        auto c_constraint = new TCanvas("c_constraint", "c_constraint", 800, 2300);
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
        c_constraint->SaveAs(Form("constraint_%s.pdf",_name.Data()));

        TString fitopt = "Q";
        auto c_pulls = new TCanvas("can2", "Pulls", 900, int(float(1200) * _parts.size() / 3));
        c_pulls->Divide(3, _parts.size());
        for (int ipart = 0; ipart < _parts.size()-_missing; ipart++)
        {
            for (int jkine = 0; jkine < KINES.size(); jkine++)
            {
                c_pulls->cd(ipart * KINES.size() + jkine + 1);
                _h_pulls[ipart * KINES.size() + jkine]->Fit("gaus", fitopt);
            }
        }
        c_pulls->SaveAs(Form("pulls_%s.pdf",_name.Data()));

        auto c_res = new TCanvas("can3", "Residuals", 900, 600);
        c_res->Divide(3, 3);
        for (int ipart = 0; ipart < _parts.size()-_missing; ipart++)
        {
            for (int jkine = 0; jkine < KINES.size(); jkine++)
            {
                c_res->cd(ipart * KINES.size() + jkine + 1);
                _h_fitres[ipart * KINES.size() + jkine]->Fit("gaus", fitopt);
            }
        }
        c_pulls->SaveAs(Form("resolution_%s.pdf",_name.Data()));

        auto c_sme = new TCanvas("can4", "Smearing", 900, 600);
        c_sme->Divide(3, 3);
        for (int ipart = 0; ipart < _parts.size()-_missing; ipart++)
        {
            for (int jkine = 0; jkine < KINES.size(); jkine++)
            {
                c_sme->cd(ipart * KINES.size() + jkine + 1);
                _h_smeres[ipart * KINES.size() + jkine]->Fit("gaus", fitopt);
                _h_fitgen[ipart * KINES.size() + jkine]->SetLineColor(kGreen);
                _h_fitgen[ipart * KINES.size() + jkine]->Draw("same");
            }
        }
        c_sme->SaveAs(Form("smearing_%s.pdf",_name.Data()));

    }

};

#endif

