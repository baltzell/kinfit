#ifndef KinConstraint_MissingMass_h
#define KinConstraint_MissingMass_h

#include "KinConstraint.h"

class KinConstraint_MissingMass : public KinConstraint
{
public:
    virtual ~KinConstraint_MissingMass() {}
    KinConstraint_MissingMass() {}

    KinConstraint_MissingMass(std::vector<int> index_in_parts)
    {
        index_Cons_Particles = index_in_parts;
        _nconstraints = 1;
    }

    TVectorD getConstraint(std::vector<TLorentzVector> _p_init_vector, std::vector<TLorentzVector> _p_fin_vector) override
    {
        TVectorD _c;
        _c.ResizeTo(_nconstraints);

        TLorentzVector _p_miss = TLorentzVector(0, 0, 0, 0);

        for (auto idx : index_Cons_Particles)
            _p_miss += _p_fin_vector[idx];
        for (auto p_in : _p_init_vector)
            _p_miss += p_in;

        _c[]=_p_miss.M();
        return _c;
    }

    TMatrixD getDfDx(int idx_part, std::vector<TLorentzVector> _particles) override 
    {
        Double_t theta = _part._vector.Theta();
        Double_t phi = _part._vector.Phi();
        Double_t p = _part._vector.P();
        Double_t E = _part._vector.E();

        Int_t nvars = 3;

        Double_t data[1][3] = {};

        TMatrixD dfdx(_nconstraints, nvars, *data);

        dfdx *= -1.;
        return dfdx;
    }
};

#endif