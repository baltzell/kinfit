#include <vector>
#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"
#include "KinFitter.h"
#include "KinParticle.h"
#include "KinFitTest.h"

int test_3C(const int max_events = 10000, const float bg_fraction = 0.1)
{
    const std::vector<TString> parts = {"p", "#pi^{+}", "#pi^{-}"};
    const std::vector<double> masses = {0.938, 0.139, 0.139};
    const std::vector<double> masses_bg = {0.938, 0.139, 0.139};

    TLorentzVector target(0.0, 0.0, 0.0, 0.938);
    TLorentzVector beam(0.0, 0.0, 10.6, 10.6);
    TLorentzVector W = beam + target;
    TGenPhaseSpace event;

    KinFitTest test("3C", parts, W);

    std::vector<double> resolutions;
    std::vector<int> constraint_idx;
    constraint_idx.push_back(0);
    for (int i = 0; i < parts.size(); i++)
    {
        for (int j = 0; j < KINES.size(); j++)
        {
            resolutions.push_back(RESO[j]);
        }
    }

    int nevents = 0;
    while (nevents < max_events)
    {
        const bool is_background = RNDM3.Uniform(0, 1) < bg_fraction;

        if (is_background)
            event.SetDecay(W, masses_bg.size(), &masses_bg[0]);
        else
            event.SetDecay(W, masses.size(), &masses[0]);

        auto weight = event.Generate();

        std::vector<TLorentzVector> parts_gen;
        std::vector<KinParticle> parts_gen_input;
        std::vector<TLorentzVector> parts_sme;
        std::vector<KinParticle> kin_parts_sme;

        int ipart = 1;
        parts_gen.push_back(*(event.GetDecay(ipart)));
        TLorentzVector sme_vector = smear(event.GetDecay(ipart));
        parts_sme.push_back(sme_vector);
        parts_gen_input.push_back(KinParticle(*(event.GetDecay(ipart))));
        kin_parts_sme.push_back(KinParticle(sme_vector, masses[ipart], RESO));

        nevents++;

        auto kin = new KinFitter(parts_gen_input, kin_parts_sme);
        kin->Add_3C_Constraint(constraint_idx);
        kin->DoFitting(100);

        test.fill_3C(kin, parts_gen, parts_sme, weight, is_background);

        delete kin;
    }

    test.plot();

    return 0;
}

int main()
{
    return test_3C();
}
