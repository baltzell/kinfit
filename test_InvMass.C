#include <random>
#include <vector>
#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"
#include "KinFitter.h"
#include "KinParticle.h"
#include "KinFitTest.h"

int test_InvMass(const int max_events=10000, const float bg_fraction=0.1)
{
    const std::vector<double> masses = {0.150, 0.150};
    const std::vector<TString> parts = {"e", "e"};

    TLorentzVector JPsi;
    JPsi.SetXYZM(0.0, 0.0, 5., 3.);
    TGenPhaseSpace event;
    event.SetDecay(JPsi, masses.size(), &masses[0]);

    KinFitTest test("IM", parts, JPsi, 0);

    std::vector<double> resolutions;
    std::vector<int> constraint_idx;
    for (int i = 0; i < parts.size(); i++)
    {
        constraint_idx.push_back(i);
        for (int j = 0; j < KINES.size(); j++)
        {
            resolutions.push_back(RESO[j]);
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

        test.fill(kin, parts_gen, parts_sme, weight, is_BG);

    }

    test.plot();

    return 0;
}

int main() {
    return test_InvMass();
}

