#include <vector>
#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"
#include "KinFitter.h"
#include "KinParticle.h"
#include "KinFitTest.h"

int test_MissAndInv(const int max_events=10000, const float bg_fraction=0.1)
{
    const double invmass = 3.0;
    const double missmass = 0.938;

    const std::vector<TString> parts_intermediate = {"JPsi", "p"};
    const std::vector<TString> parts_final = {"e^{+}", "e^{-}", "p"};
    const std::vector<TString> parts_decay = {"e^{+}", "e^{-}"};

    const std::vector<double> masses_intermediate = {invmass, missmass};
    const std::vector<double> masses_final = {0.01, 0.01, missmass};
    const std::vector<double> masses_decay = {0.01, 0.01};

    TLorentzVector target(0.0, 0.0, 0.0, missmass);
    TLorentzVector beam(0.0, 0.0, 10.6, 10.6);
    TLorentzVector W = beam + target;

    TGenPhaseSpace event_intermediate;
    event_intermediate.SetDecay(W, masses_intermediate.size(), &masses_intermediate[0]);
   
    KinFitTest test("2C", parts_final, W, 1);

    std::vector<double> resolutions;
    std::vector<int> constraint_idx;
    for (int i = 0; i < parts_final.size() - 1; i++)
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
        // generate the intermediate reaction:
        auto weight = event_intermediate.Generate();

        // add on the decay:
        TLorentzVector JPsi = *(event_intermediate.GetDecay(0));
        TGenPhaseSpace decay;
        decay.SetDecay(JPsi, masses_decay.size(), &masses_decay[0]);
        weight *= decay.Generate();

        std::vector<TLorentzVector> parts_gen;
        std::vector<TLorentzVector> parts_sme;
        std::vector<KinParticle> kin_parts_sme;

        for (int ipart = 0; ipart < parts_decay.size(); ++ipart)
        {
            parts_gen.push_back(*(decay.GetDecay(ipart)));
            TLorentzVector sme_vector = smear(decay.GetDecay(ipart));
            parts_sme.push_back(sme_vector);
            kin_parts_sme.push_back(KinParticle(sme_vector, masses_decay[ipart], RESO));
        }

        //parts_gen.push_back(*(event_intermediate.GetDecay(parts_intermediate.size()-1)));
        //TLorentzVector sme_vector = smear(event_intermediate.GetDecay(parts_intermediate.size()-1));
        //parts_sme.push_back(sme_vector);
        //kin_parts_sme.push_back(KinParticle(sme_vector, masses_intermediate[masses_intermediate.size()-1], RESO));

        nevents++;

        auto kin = new KinFitter({KinParticle(target), KinParticle(beam)}, kin_parts_sme);
        kin->Add_MissAndInvMass_Constraint(constraint_idx, invmass, missmass);
        kin->DoFitting(100);

        test.fill(kin, parts_gen, parts_sme, weight, false);
    }

    test.plot();

    return 0;
}

int main() {
    return test_MissAndInv();
}

