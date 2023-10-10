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

    // Constructor with diagonal covariance matrix
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

    // Constructor with covariance matrix constructed using TH3D
    KinParticle(TLorentzVector in_vector, double in_mass, KinCovariance in_Cov, int sector)
        : _vector(in_vector),
          _mass(in_mass)
    {
        _good_cov_matrix = in_Cov.Is_good_to_interpolate(sector, in_vector);
        SetCovMatrix(in_Cov.Interpolate(sector, in_vector));
        SetErrCovMatrix(in_Cov.Interpolate_Error(sector, in_vector));
        SetEntries(in_Cov.GetEntries());
    }

    TLorentzVector GetVector() { return _vector; }
    TMatrixD GetCovMatrix() { return _Cov; }
    TMatrixD GetErrCovMatrix() { return _Err_Cov; }
    double GetMass() { return _mass; }
    double GetEntries() { return _n_entries; }

    void SetVector(TLorentzVector in_vector) { _vector = in_vector; }
    void SetMass(double in_mass) { _mass = in_mass; }
    void SetCovMatrix(TMatrixD in_Cov) { _Cov = in_Cov; }
    void SetErrCovMatrix(TMatrixD in_Err_Cov) { _Err_Cov = in_Err_Cov; }
    void SetEntries(int in_n_entries) { _n_entries = in_n_entries; }

    bool Is_good_cov_matrix() {return  _good_cov_matrix; }

private:
    TLorentzVector _vector;
    double _mass;
    TMatrixD _Cov = TMatrixD(3, 3);
    TMatrixD _Err_Cov = TMatrixD(3, 3);
    int _n_entries;
    bool _good_cov_matrix = true;
};

#endif
