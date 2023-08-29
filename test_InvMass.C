#include <vector>
#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"
#include "KinFitter.h"
#include "KinParticle.h"
#include "KinFitTest.h"

int test_InvMass(const int max_events=10000, const float bg_fraction=0.1)
{
    const double invmass = 3.0;
    const std::vector<double> masses = {0.150, 0.150, 0.150};
    const std::vector<TString> parts = {"e", "e", "e"};

    TLorentzVector JPsi;
    JPsi.SetXYZM(1.0, 1.0, 5., invmass);
    TGenPhaseSpace event;
    event.SetDecay(JPsi, masses.size(), &masses[0]);

    KinFitTest test("IM", parts, JPsi, 0, 0, invmass);

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

    int nevents = 0;
    while (nevents < max_events)
    {
        const bool is_background = RNDM3.Uniform(0,1) < bg_fraction;

        if (is_background)
        {
            float mass = RNDM3.Uniform(invmass-0.5, invmass+0.5);
            TLorentzVector BG;
            BG.SetXYZM(1.0, 1.0, 5., mass);
            event.SetDecay(BG, masses.size(), &masses[0]);
        }
        else
            event.SetDecay(JPsi, masses.size(), &masses[0]);

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
        kin->Add_InvMass_Constraint(constraint_idx, invmass);
        kin->DoFitting(100);

        //test.fill_MissingMass(kin, parts_gen, parts_sme, weight, is_background);
        //test.fill_InvariantMass(kin, parts_gen, parts_sme, constraint_idx, weight, is_background);
    }

    test.plot();

    return 0;
}

int main() {
    return test_InvMass();
}

