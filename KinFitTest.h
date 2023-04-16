#ifndef KinFitTest_h
#define KinFitTest_h

#include <vector>
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TRandom3.h"
#include "KinFitter.h"
#include <iostream>


// absolute resolutions (GeV,radians,radians):
//const std::vector<double> RESO = {0.15, 0.02, 0.02};
const std::vector<double> RESO = {0.05, 0.0005, 0.002};
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

    TH1* _h_inv_gen;
    TH1* _h_inv_sme;
    TH1* _h_inv_fit;

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
    std::vector<TH1*> _h_pulls_inv_conf_cut;
    std::vector<TH1*> _h_fitres;
    std::vector<TH1*> _h_smeres;
    std::vector<TH1*> _h_fitgen;
    std::vector<TH2*> _h2_pulls;

    KinFitTest(TString name, std::vector<TString> parts, TLorentzVector W, int missing=0, float missmass=0, float invmass=0)
        :
        _name(name),
        _parts(parts),
        _W(W),
        _missing(missing)
    {
        _h_chi = new TH1F("h_chi", ";#chi^{2}/ndf", 201, 0, 10);
	_h_chi->GetXaxis()->SetTitleSize(0.08);
        _h_lik_Signal = new TH1F("h_lik_Signal", ";Confidence Level", 100, 0, 1);
	_h_lik_Signal->GetXaxis()->SetTitleSize(0.08);
        _h_lik_BG = new TH1F("h_lik_BG", ";Confidence Level", 100, 0, 1);
	_h_lik_BG->GetXaxis()->SetTitleSize(0.08);

        _h_mm_gen = new TH1F("h_mm_gen", ";Missing Mass [GeV^{2}]", 201, missmass - 0.2 * (parts.size() - 1), missmass + 0.2 * (parts.size() - 1));
        _h_mm_gen->GetXaxis()->SetTitleSize(0.08);
	_h_mm_sme = (TH1 *)_h_mm_gen->Clone("h_mm_sme");
        _h_mm_fit = (TH1 *)_h_mm_gen->Clone("h_mm_fit");

        _h_inv_gen = new TH1F("h_inv_gen", ";Invariant Mass [GeV]", 201, invmass - 0.5 * (parts.size() - 1), invmass + 0.5 * (parts.size() - 1));
	//_h_inv_gen->GetXaxis()->SetTitleSize(0.08);
        _h_inv_sme = (TH1 *)_h_inv_gen->Clone("h_inv_sme");
        _h_inv_fit = (TH1 *)_h_inv_gen->Clone("h_inv_sme");

        _h_E_gen = new TH1F("h_E_gen", ";#Delta E [GeV]", 201, -0.5 * (parts.size() - 1), 0.5 * (parts.size() - 1));
	_h_E_gen->GetXaxis()->SetTitleSize(0.08);
        _h_E_gen->GetXaxis()->SetLabelSize(0.07);
	_h_E_gen->GetYaxis()->SetLabelSize(0.07);
	_h_E_sme = (TH1 *)_h_E_gen->Clone("h_E_sme");
        _h_E_fit = (TH1 *)_h_E_gen->Clone("h_E_fit");

        _h_Px_gen = new TH1F("h_Px_gen", ";#Delta Px [GeV]", 201, -0.5 * (parts.size() - 1), 0.5 * (parts.size() - 1));
	_h_Px_gen->GetXaxis()->SetTitleSize(0.08);
	_h_Px_gen->GetXaxis()->SetLabelSize(0.07);
	_h_Px_gen->GetYaxis()->SetLabelSize(0.07);
        _h_Px_sme = (TH1 *)_h_Px_gen->Clone("h_Px_sme");
        _h_Px_fit = (TH1 *)_h_Px_gen->Clone("h_Px_fit");

        _h_Py_gen = new TH1F("h_Py_gen", ";#Delta Py [GeV]", 201, -0.5 * (parts.size() - 1), 0.5 * (parts.size() - 1));
	_h_Py_gen->GetXaxis()->SetTitleSize(0.08);
	_h_Py_gen->GetXaxis()->SetLabelSize(0.07);
	_h_Py_gen->GetYaxis()->SetLabelSize(0.07);
        _h_Py_sme = (TH1 *)_h_Py_gen->Clone("h_Py_sme");
        _h_Py_fit = (TH1 *)_h_Py_gen->Clone("h_Py_fit");

        _h_Pz_gen = new TH1F("h_Pz_gen", ";#Delta Pz [GeV]", 201, -0.5 * (parts.size() - 1), 0.5 * (parts.size() - 1));
	_h_Pz_gen->GetXaxis()->SetTitleSize(0.08);
	_h_Pz_gen->GetXaxis()->SetLabelSize(0.07);
	_h_Pz_gen->GetYaxis()->SetLabelSize(0.07);
        _h_Pz_sme = (TH1 *)_h_Pz_gen->Clone("h_Pz_sme");
        _h_Pz_fit = (TH1 *)_h_Pz_gen->Clone("h_Pz_fit");

        for (int i = 0; i < parts.size()-_missing; i++)
        {
            for (int j = 0; j < KINES.size(); j++)
            {
	      //_h_pulls.push_back(new TH1F(Form("h_pulls_%d_%d", i, j),
	      //            Form("%s;%s Normalized Pull", parts[i].Data(), KINES[j].Data()), 201, -10, 10));
	      TH1* _h_pulls_temp = new TH1F(Form("h_pulls_%d_%d", i, j), 
	             Form("%s;%s Normalized Pull", parts[i].Data(), KINES[j].Data()), 201, -10, 10);
	      _h_pulls_temp->GetXaxis()->SetTitleSize(0.07);
	      _h_pulls.push_back(_h_pulls_temp);
	      TH1* _h_pulls_inv_conf_cut_temp = (TH1 *)_h_pulls_temp->Clone(Form("h_pulls_inv_conf_cut_%d_%d", i, j));
	      _h_pulls_inv_conf_cut.push_back((TH1 *)_h_pulls_temp->Clone(Form("h_pulls_inv_conf_cut_%d_%d", i, j)));
              //_h_fitres.push_back(new TH1F(Form("h_fitres_%d_%d", i, j),
              //              Form("%s;%s Residual [%s]", parts[i].Data(), KINES[j].Data(), UNITS[j].Data()), 201, -RESO[j] * 10, RESO[j] * 10));
	      TH1* _h_fitres_temp = new TH1F(Form("h_fitres_%d_%d", i, j),
					   Form("%s;%s Residual [%s]", parts[i].Data(), KINES[j].Data(), UNITS[j].Data()), 201, -RESO[j] * 10, RESO[j] * 10);
	      _h_fitres_temp->GetXaxis()->SetTitleSize(0.07);
	      _h_fitres.push_back(_h_fitres_temp);
              //_h_smeres.push_back(new TH1F(Form("h_smeres_%d_%d", i, j),
	      //	    Form("%s;%s Smearing [%s]", parts[i].Data(), KINES[j].Data(), UNITS[j].Data()), 201, -RESO[j] * 10, RESO[j] * 10));
	      TH1* _h_smeres_temp = new TH1F(Form("h_smeres_%d_%d", i, j),
					   Form("%s;%s Smearing [%s]", parts[i].Data(), KINES[j].Data(), UNITS[j].Data()), 201, -RESO[j] * 10, RESO[j] * 10);
	      _h_smeres_temp->GetXaxis()->SetTitleSize(0.07);
	      _h_smeres.push_back(_h_smeres_temp);
              _h_fitgen.push_back(new TH1F(Form("h_fitgen%d_%d", i, j),
                            Form("%s;%s Fit/Gen [%s]", parts[i].Data(), KINES[j].Data(), UNITS[j].Data()), 201, -RESO[j] * 10, RESO[j] * 10));
	      if (j == KINES.size() - 1) {
		//_h2_pulls.push_back(new TH2F(Form("h2_pulls_%d_%d_%d", i, j, 0), (Form(";%s_{%s} Pull", parts[i].Data(), KINES[j].Data()), Form(";%s_{%s} Pull", parts[i].Data(), KINES[0].Data())), 201, -10, 10, 201, -10, 10));
		TH2* _h2_pulls_temp = new TH2F(Form("h2_pulls_%d_%d_%d", i, 0, j), Form("%s;%s Normalized Pull;%s Normalized Pull", parts[i].Data(), KINES[0].Data(), KINES[j].Data()), 201, -5, 5, 201, -5, 5);
		_h2_pulls_temp->GetXaxis()->SetTitleSize(0.07);
		_h2_pulls_temp->GetYaxis()->SetTitleSize(0.07);
		_h2_pulls.push_back(_h2_pulls_temp);
	      }
	      else {
		//_h2_pulls.push_back(new TH2F(Form("h2_pulls_%d_%d_%d", i, j, j+1), Form(";%s_{%s} Pull", parts[i].Data(), KINES[j].Data(), Form(";%s_{%s} Pull", parts[i].Data(), KINES[j+1].Data()), 201, -10, 10, 201, -10, 10));
		TH2* _h2_pulls_temp = new TH2F(Form("h2_pulls_%d_%d_%d", i, j+1, j),  Form("%s;%s Normalized Pull;%s Normalized Pull", parts[i].Data(), KINES[j+1].Data(), KINES[j].Data()), 201, -5, 5, 201, -5, 5); 
		_h2_pulls_temp->GetXaxis()->SetTitleSize(0.07);
		_h2_pulls_temp->GetYaxis()->SetTitleSize(0.07);
		_h2_pulls.push_back(_h2_pulls_temp);
	      }
	    }
        }
    }

    void fill_MissingMass(KinFitter *kin, std::vector<TLorentzVector> parts_gen, std::vector<TLorentzVector> parts_sme, double weight, bool background)
    {
        std::vector<TLorentzVector> parts_fit = kin->GetFitted4Vectors();
        
        TLorentzVector missing_gen = _W;
        TLorentzVector missing_sme = _W;
        TLorentzVector missing_fit = _W;
	double parts_fit_phi;
	double parts_rec_phi;
        for (int ipart = 0; ipart < _parts.size()-_missing; ipart++)
        {
            missing_gen -= parts_gen[ipart];
            missing_sme -= parts_sme[ipart];
            missing_fit -= parts_fit[ipart];
	    parts_fit_phi = parts_fit[ipart].Phi();
	    if (parts_fit_phi < -TMath::Pi()/2) {
              parts_fit_phi = TMath::Pi() + (TMath::Abs(TMath::Pi() - TMath::Abs(parts_fit_phi)));
            }
	    parts_rec_phi = parts_sme[ipart].Phi();
	    if (parts_rec_phi < -TMath::Pi()/2) {
	      parts_rec_phi = TMath::Pi() + (TMath::Abs(TMath::Pi() - TMath::Abs(parts_rec_phi)));
	    }
	    //std::cout << "Rec part phi = " << parts_rec_phi << std::endl;
	    //std::cout << "N mc parts = " << parts_gen.size() << std::endl;
	    //std::cout << "Thetas: missing gen: " << missing_gen.Theta() << ", missing rec: " << missing_sme.Theta() << ", missing fit: " << missing_fit.Theta() << std::endl;
	    //std::cout << "mc part px = " << parts_gen[ipart].Px() << ", mc part py = " << parts_gen[ipart].Py() << ", mc part pz = " << parts_gen[ipart].Pz() << std::endl; 
	    
            for (int jkine = 0; jkine < KINES.size() && kin->GetConfidenceLevel()>0.01; jkine++)
            {
                _h_pulls[ipart * KINES.size() + jkine]->Fill(kin->GetPulls()[ipart * KINES.size() + jkine]);

		if (jkine == (KINES.size() - 1)) {
		  _h2_pulls[ipart * KINES.size() + jkine]->Fill(kin->GetPulls()[ipart * KINES.size()], kin->GetPulls()[ipart * KINES.size() + jkine]);
		}
		else {
		  _h2_pulls[ipart * KINES.size() + jkine]->Fill(kin->GetPulls()[ipart * KINES.size() + jkine + 1], kin->GetPulls()[ipart * KINES.size() + jkine]);
		}		
            }
	    for (int jkine = 0; jkine < KINES.size() && kin->GetConfidenceLevel()<=0.01; jkine++)
	      {
		_h_pulls_inv_conf_cut[ipart * KINES.size() + jkine]->Fill(kin->GetPulls()[ipart * KINES.size() + jkine]);
	      }
            _h_fitres[ipart * 3 + 0]->Fill(parts_fit[ipart].Vect().Mag() - parts_sme[ipart].Vect().Mag());
            _h_fitres[ipart * 3 + 1]->Fill(parts_fit[ipart].Theta() - parts_sme[ipart].Theta());
            //_h_fitres[ipart * 3 + 2]->Fill(parts_fit[ipart].Phi() - parts_sme[ipart].Phi());
	    _h_fitres[ipart * 3 + 2]->Fill(parts_fit_phi - parts_rec_phi);

            _h_smeres[ipart * 3 + 0]->Fill(parts_gen[ipart].Vect().Mag() - parts_sme[ipart].Vect().Mag());
            _h_smeres[ipart * 3 + 1]->Fill(parts_gen[ipart].Theta() - parts_sme[ipart].Theta());
            //_h_smeres[ipart * 3 + 2]->Fill(parts_gen[ipart].Phi() - parts_sme[ipart].Phi());
	    _h_smeres[ipart * 3 + 2]->Fill(parts_fit_phi - parts_rec_phi);

            _h_fitgen[ipart * 3 + 0]->Fill(parts_gen[ipart].Vect().Mag() - parts_fit[ipart].Vect().Mag());
            _h_fitgen[ipart * 3 + 1]->Fill(parts_gen[ipart].Theta() - parts_fit[ipart].Theta());
            _h_fitgen[ipart * 3 + 2]->Fill(parts_gen[ipart].Phi() - parts_fit[ipart].Phi());
        }
	//std::cout << "MC MM2 = " << missing_gen.M2() << std::endl;
        //_h_mm_gen->Fill(missing_gen.M2(), weight);
        //_h_mm_sme->Fill(missing_sme.M2(), weight);
        //_h_mm_fit->Fill(missing_fit.M2(), weight);
	_h_mm_gen->Fill(missing_gen.M2());
        _h_mm_sme->Fill(missing_sme.M2());
        _h_mm_fit->Fill(missing_fit.M2());
	//std::cout << "M2: missing gen: " << missing_gen.M2() << ", missing rec: " << missing_sme.M2() << ", missing fit: " << missing_fit.M2() << std::endl;

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
        
        //std::cout<<kin->GetConfidenceLevel()<<"\n";

        if (background)
            _h_lik_Signal->Fill(kin->GetConfidenceLevel());
        else
            _h_lik_BG->Fill(kin->GetConfidenceLevel());
        

    }

    void fill_InvariantMass(KinFitter *kin, std::vector<TLorentzVector> parts_gen, std::vector<TLorentzVector> parts_sme, std::vector<int> indices_part, double weight)
    {
        std::vector<TLorentzVector> parts_fit = kin->GetFitted4Vectors();

        TLorentzVector invariant_gen;
        invariant_gen.SetXYZT(0,0,0,0);
        TLorentzVector invariant_sme;
        invariant_sme.SetXYZT(0,0,0,0);
        TLorentzVector invariant_fit;
        invariant_fit.SetXYZT(0,0,0,0);
        
        for (int ipart = 0; ipart < indices_part.size(); ipart++)
        {
            invariant_gen += parts_gen[ipart];
            invariant_sme += parts_sme[ipart];
            invariant_fit += parts_fit[ipart];
        }
        
        _h_inv_gen->Fill(invariant_gen.M(), weight);
        _h_inv_sme->Fill(invariant_sme.M(), weight);
        _h_inv_fit->Fill(invariant_fit.M(), weight);

    }


    void plot()
    {
      int can_size_x = 800;
      int can_size_y = 800;
        gStyle->SetOptStat(0);
        gStyle->SetOptFit(1111);

        auto c_missing = new TCanvas("can1", "Summary", can_size_x, 1000);
	//c_missing->SetWindowSize(can_size_x, can_size_y);
	//c_missing->SetCanvasSize(can_size_x, can_size_y);
        c_missing->Divide(1, 3);
        c_missing->cd(1);
	gPad->SetBottomMargin(0.2);
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
        legend->AddEntry(_h_mm_gen, "Gen", "l");
        legend->AddEntry(_h_mm_sme, "Rec", "l");
        legend->AddEntry(_h_mm_fit, "Fitted", "l");
        legend->SetFillStyle(0);
        legend->SetLineWidth(0);
        legend->Draw("same ");
        c_missing->cd(3);
	gPad->SetBottomMargin(0.2);
        _h_chi->Draw();
        c_missing->cd(2);
	gPad->SetBottomMargin(0.2);
        gPad->SetLogy();
        //_h_lik_Signal->GetYaxis()->SetRangeUser(1.0, _h_lik_Signal->GetMaximum());
        _h_lik_Signal->SetLineColor(kRed);
        _h_lik_Signal->Draw();
	_h_lik_Signal->GetYaxis()->SetRangeUser(1.0, _h_lik_BG->GetMaximum());
	std::cout << "Conf Level max = " << _h_lik_BG->GetMaximum() << std::endl;
        _h_lik_BG->SetLineColor(kBlue);
        _h_lik_BG->Draw("same");
        c_missing->SaveAs(Form("%s.pdf(",_name.Data()));
	//c_missing->Print(Form("%s.pdf(",_name.Data()));

        //auto c_invariant = new TCanvas("can_inv", "Summary", 800, 800);
	auto c_invariant = new TCanvas("can_inv", "Summary", can_size_x, can_size_y);
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
        legend->Draw("same ");
        c_invariant->SaveAs(Form("%s.pdf",_name.Data()));

        if (_missing==0) {
	  //auto c_constraint = new TCanvas("c_constraint", "Constraints", 800, 1200);
	    auto c_constraint = new TCanvas("c_constraint", "Constraints", can_size_x, can_size_y);
	    c_constraint->Divide(1, 4);
            c_constraint->cd(1);
	    gPad->SetBottomMargin(0.175);
            gPad->SetLogy();
            _h_E_gen->Draw("hist");
            _h_E_fit->SetLineColor(kGreen);
            _h_E_fit->SetLineWidth(2);
            _h_E_fit->Draw("hist same");
            _h_E_sme->SetLineColor(kRed);
            _h_E_sme->Draw("hist same");
            legend->Draw("same ");
            c_constraint->cd(2);
	    gPad->SetBottomMargin(0.175);
            gPad->SetLogy();
            _h_Px_gen->Draw("hist");
            _h_Px_fit->SetLineColor(kGreen);
            _h_Px_fit->SetLineWidth(2);
            _h_Px_fit->Draw("hist same");
            _h_Px_sme->SetLineColor(kRed);
            _h_Px_sme->Draw("hist same");
            c_constraint->cd(3);
	    gPad->SetBottomMargin(0.175);
            gPad->SetLogy();
            _h_Py_gen->Draw("hist");
            _h_Py_fit->SetLineColor(kGreen);
            _h_Py_fit->SetLineWidth(2);
            _h_Py_fit->Draw("hist same");
            _h_Py_sme->SetLineColor(kRed);
            _h_Py_sme->Draw("hist same");
            c_constraint->cd(4);
	    gPad->SetBottomMargin(0.175);
            gPad->SetLogy();
            _h_Pz_gen->Draw("hist");
            _h_Pz_fit->SetLineColor(kGreen);
            _h_Pz_fit->SetLineWidth(2);
            _h_Pz_fit->Draw("hist same");
            _h_Pz_sme->SetLineColor(kRed);
            _h_Pz_sme->Draw("hist same");
            c_constraint->SaveAs(Form("%s.pdf",_name.Data()));
        }

        TString fitopt = "Q";
        //auto c_pulls = new TCanvas("can2", "Pulls", 900, int(float(1200) * (_parts.size()-_missing) / 3));
        auto c_pulls = new TCanvas("can2", "Pulls", can_size_x, can_size_y);
	c_pulls->Divide(3, _parts.size()-_missing);
        for (int ipart = 0; ipart < _parts.size()-_missing; ipart++)
        {
            for (int jkine = 0; jkine < KINES.size(); jkine++)
            {
	      c_pulls->cd(ipart * KINES.size() + jkine + 1);
	      gStyle->SetTitleFontSize(0.08);
	      gPad->SetBottomMargin(0.175);
		_h_pulls[ipart * KINES.size() + jkine]->Fit("gaus", fitopt);
		_h_pulls_inv_conf_cut[ipart * KINES.size() + jkine]->Draw("hist same");
		_h_pulls_inv_conf_cut[ipart * KINES.size() + jkine]->SetLineColor(kGreen);
            }
        }
        c_pulls->SaveAs(Form("%s.pdf",_name.Data()));

        //auto c_res = new TCanvas("can3", "Residuals", 900, int(float(1200) * (_parts.size()-_missing) / 3));
        auto c_res = new TCanvas("can3", "Residuals", can_size_x, can_size_y);
	c_res->Divide(3, _parts.size()-_missing);
        for (int ipart = 0; ipart < _parts.size()-_missing; ipart++)
        {
            for (int jkine = 0; jkine < KINES.size(); jkine++)
            {
                c_res->cd(ipart * KINES.size() + jkine + 1);
		gStyle->SetTitleFontSize(0.08);
		gPad->SetBottomMargin(0.175);
                _h_fitres[ipart * KINES.size() + jkine]->Fit("gaus", fitopt);
            }
        }
        c_res->SaveAs(Form("%s.pdf",_name.Data()));

        //auto c_sme = new TCanvas("can4", "Smearing", 900, int(float(1200) * (_parts.size()-_missing) / 3));
        auto c_sme = new TCanvas("can4", "Smearing", can_size_x, can_size_y);
	c_sme->Divide(3, _parts.size()-_missing);
        for (int ipart = 0; ipart < _parts.size()-_missing; ipart++)
        {
            for (int jkine = 0; jkine < KINES.size(); jkine++)
            {
                c_sme->cd(ipart * KINES.size() + jkine + 1);
		gStyle->SetTitleFontSize(0.08);
		gPad->SetBottomMargin(0.175);
                _h_smeres[ipart * KINES.size() + jkine]->Fit("gaus", fitopt);
		//gStyle->SetStatFontSize(10);
                gStyle->SetStatW(.5);
		//_h_fitgen[ipart * KINES.size() + jkine]->SetLineColor(kGreen);
                //_h_fitgen[ipart * KINES.size() + jkine]->Draw("same");
            }
        }
        c_sme->SaveAs(Form("%s.pdf",_name.Data()));

	//auto c_2D_pulls = new TCanvas("can5", "2D Pulls", 900, int(float(1200) * (_parts.size()-_missing) / 3));
        auto c_2D_pulls = new TCanvas("can5", "2D Pulls", can_size_x, can_size_x);
	c_2D_pulls->Divide(3, _parts.size()-_missing);
        for (int ipart = 0; ipart < _parts.size()-_missing; ipart++)
	  {
            for (int jkine = 0; jkine < KINES.size(); jkine++)
	      {
                c_2D_pulls->cd(ipart * KINES.size() + jkine + 1);
		gStyle->SetTitleFontSize(0.08);
		gPad->SetBottomMargin(0.175);
		gPad->SetLeftMargin(0.15);
                _h2_pulls[ipart * KINES.size() + jkine]->Draw("Colz");
		//_h_pulls[ipart * KINES.size() + jkine]->Fit("gaus", fitopt);
		gStyle->SetLabelSize(0.01, "xyz");
	      }
	  }
	
	c_2D_pulls->SaveAs(Form("%s.pdf)",_name.Data()));
	//c_2D_pulls->SaveAs(Form("2D_pulls.pdf)"));
    }

};

#endif

