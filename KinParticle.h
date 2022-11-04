#ifndef KinParticle_h
#define KinParticle_h

class KinParticle
{

public:
    virtual ~KinParticle() {}
    KinParticle() {}

    KinParticle(TLorentzVector in_vector, TMatrixD in_Cov)
        : _vector(in_vector),
          _Cov(in_Cov)
    {
    }

    KinParticle(TLorentzVector in_vector)
        : _vector(in_vector)
    {
    }

    KinParticle(TLorentzVector in_vector, std::vector<Double_t> in_sigmas)
        : _vector(in_vector)
    {
        TMatrixD _Cov_from_sigmas(in_sigmas.size(), in_sigmas.size());
        for (Int_t ii = 0; ii < in_sigmas.size(); ii++)
        {
            _Cov_from_sigmas[ii][ii] = in_sigmas[ii];
        }
        _Cov_from_sigmas.Sqr();
        SetCovMatrix(_Cov_from_sigmas);
    }

public:
    TLorentzVector GetVector() { return _vector; }
    TMatrixD GetCovMatrix() { return _Cov; }

    void SetVector(TLorentzVector in_vector) { _vector = in_vector; }
    void SetCovMatrix(TMatrixD in_Cov) { _Cov = in_Cov; }

public:
    TLorentzVector _vector;
    TMatrixD _Cov = TMatrixD(3, 3);
};

#endif