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
        _index_Cons_Particles=index_in_parts;
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

    TMatrixD getDfDx(int idx_part, std::vector<TLorentzVector> _particles) override 
    {


        Double_t theta = _particles[idx_part].Theta();
        Double_t phi = _particles[idx_part].Phi();
        Double_t p = _particles[idx_part].P();
        Double_t E = _particles[idx_part].E();

        Int_t nvars = 3;

        Double_t data[4][3] = {{sin(theta) * cos(phi), p * cos(theta) * cos(phi), -p * sin(theta) * sin(phi)},
                               {sin(theta) * sin(phi), p * cos(theta) * sin(phi), p * sin(theta) * cos(phi)},
                               {cos(theta), -p * sin(theta), 0},
                               {p / E, 0, 0}};

        //cout<<data[0][0]<<" "<<data[0][1]<<" "<<data[0][2]<<endl;

       

        TMatrixD dfdx(_nconstraints, nvars, *data);

        dfdx *= -1.; //overall sign for the constraint to be correct
        return dfdx;
    }

};

#endif