#ifndef KinParticle_h
#define KinParticle_h

class KinParticle
{

public:
    virtual ~KinParticle() {}
    KinParticle() {}

    KinParticle(TLorentzVector in_vector, double in_mass, TMatrixD in_Cov)
        : _vector(in_vector),
          _mass(in_mass),
          _Cov(in_Cov)
    {
    }

    KinParticle(TLorentzVector in_vector, double in_mass)
        : _vector(in_vector),
          _mass(in_mass)
    {
    }

    KinParticle(TLorentzVector in_vector)
        : _vector(in_vector),
          _mass(in_vector.M())
    {
    }

    KinParticle(TLorentzVector in_vector, double in_mass, std::vector<Double_t> in_sigmas)
        : _vector(in_vector),
          _mass(in_mass)
    {
        TMatrixD _Cov_from_sigmas(in_sigmas.size(), in_sigmas.size());
        for (int ii = 0; ii < in_sigmas.size(); ii++)
        {
            _Cov_from_sigmas[ii][ii] = in_sigmas[ii];
        }
        _Cov_from_sigmas.Sqr();
        SetCovMatrix(_Cov_from_sigmas);
    }

    TLorentzVector GetVector() { return _vector; }
    TMatrixD GetCovMatrix() { return _Cov; }
    double GetMass() { return _mass; }

    void SetVector(TLorentzVector in_vector) { _vector = in_vector; }
    void SetMass(double in_mass) { _mass = in_mass; }
    void SetCovMatrix(TMatrixD in_Cov) { _Cov = in_Cov; }

private:
    TLorentzVector _vector;
    double _mass;
    TMatrixD _Cov = TMatrixD(3, 3);
};

#endif
