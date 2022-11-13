#include <random>
#include <vector>
#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"
#include "TCanvas.h"
#include "KinFitter.h"
#include "KinParticle.h"
#include "KinFitTest.h"

int test_MissingParticle(const int max_events=10000, const float bg_fraction=0.1)
{
    const std::vector<TString> parts = {"#pi^{+}", "#pi^{-}", "p"};
    const std::vector<double> masses = {0.139, 0.139, 0.938};
    const std::vector<double> masses_bg = {0.139, 0.139, 0.938, 0.000};

    TLorentzVector target(0.0, 0.0, 0.0, 0.938);
    TLorentzVector beam(0.0, 0.0, 10.6, 10.6);
    TLorentzVector W = beam + target;
    TGenPhaseSpace event;
    event.SetDecay(W, masses.size(), &masses[0]);
    
    KinFitTest test("MM", parts, W, 1);

    std::vector<double> resolutions;
    std::vector<int> constraint_idx;
    for (int i = 0; i < parts.size() - 1; i++)
    {
        constraint_idx.push_back(i);
        for (int j = 0; j < KINES.size(); j++)
        {
            resolutions.push_back(RESO[j]);
        }
    }

    int nevents = 0;
    while (nevents < max_events)
    {

        const bool is_background = RNDM3.Uniform(0,1) < bg_fraction;

        if (is_background) 
            event.SetDecay(W, masses_bg.size(), &masses_bg[0]);
        else
            event.SetDecay(W, masses.size(), &masses[0]);

        auto weight = event.Generate();

        std::vector<TLorentzVector> parts_gen;
        std::vector<TLorentzVector> parts_sme;
        std::vector<KinParticle> kin_parts_sme;

        for (int ipart = 0; ipart < parts.size() - 1; ++ipart)
        {
            parts_gen.push_back(*(event.GetDecay(ipart)));
            TLorentzVector sme_vector = smear(event.GetDecay(ipart));
            parts_sme.push_back(sme_vector);
            kin_parts_sme.push_back(KinParticle(sme_vector, masses[ipart], RESO));
        }

        nevents++;

        auto kin = new KinFitter({KinParticle(target), KinParticle(beam)}, kin_parts_sme);
        kin->Add_MissingMass_Constraint(constraint_idx, 0.938);
        kin->DoFitting(100);

        test.fill(kin, parts_gen, parts_sme, weight, is_background);
    }

    test.plot();

    return 0;
}

int main() {
    return test_MissingParticle();
}

