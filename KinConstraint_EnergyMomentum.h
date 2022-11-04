#ifndef KinConstraint_EnergyMomentum_h
#define KinConstraint_EnergyMomentum_h

#include "KinConstraint.h"

class KinConstraint_EnergyMomentum : public KinConstraint
{
public:
    virtual ~KinConstraint_EnergyMomentum() {}
    KinConstraint_EnergyMomentum() {}

    KinConstraint_EnergyMomentum(std::vector<int> index_in_parts)
    {
        index_Cons_Particles=index_in_parts;
        _nconstraints = 4;
    }

    TVectorD getConstraint( std::vector<TLorentzVector>_p_init_vector, std::vector<TLorentzVector> _p_fin_vector) override
    {
        TVectorD _c;
        _c.ResizeTo(_nconstraints);
       
        TLorentzVector _p_fin = TLorentzVector(0, 0, 0, 0);
        for (auto P : _p_fin_vector)
            _p_fin += P;

        TLorentzVector _p_init = TLorentzVector(0, 0, 0, 0);
        for (auto P : _p_init_vector)
            _p_init += P;

        TLorentzVector _p_diff = _p_init - _p_fin;

        for (Int_t ii = 0; ii < 4; ii++)
        {
            _c[ii] = _p_diff[ii];
        }

        return _c;
    }

    TMatrixD getDfDx(KinParticle _part) override 
    {
        Double_t theta = _part._vector.Theta();
        Double_t phi = _part._vector.Phi();
        Double_t p = _part._vector.P();
        Double_t E = _part._vector.E();

        Int_t nvars = 3;

        Double_t data[4][3] = {{sin(theta) * cos(phi), p * cos(theta) * cos(phi), -p * sin(theta) * sin(phi)},
                               {sin(theta) * sin(phi), p * cos(theta) * sin(phi), p * sin(theta) * cos(phi)},
                               {cos(theta), -p * sin(theta), 0},
                               {p / E, 0, 0}};

       

        TMatrixD dfdx(_nconstraints, nvars, *data);

        dfdx *= -1.;
        return dfdx;
    }

};

#endif