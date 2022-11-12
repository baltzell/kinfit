#ifndef test_h
#define test_h

// absolute resolutions (GeV,radians,radians):
const std::vector<double> RESO = {0.15, 0.02, 0.02};
const std::vector<TString> KINES = {"P", "#theta", "#phi"};
const std::vector<TString> UNITS = {"GeV", "rad", "rad"};

TRandom3 RNDM3(0);

auto smear(TLorentzVector *v)
{
    // smear by absolute resolutions:
    const double p = RNDM3.Gaus(0.0, RESO[0]) + v->P();
    const double t = RNDM3.Gaus(0.0, RESO[1]) + v->Theta();
    const double h = RNDM3.Gaus(0.0, RESO[2]) + v->Phi();

    const double m = v->M();
    TLorentzVector v_out;
    v_out.SetXYZM(p * sin(t) * cos(h), p * sin(t) * sin(h), p * cos(t), m);
    return v_out;
};

#endif
