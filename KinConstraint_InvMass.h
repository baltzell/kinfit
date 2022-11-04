#ifndef KinConstraint_InvMass_h
#define KinConstraint_InvMass_h

#include "KinConstraint.h"

class KinConstraint_InvMass : public KinConstraint
{
public:
    virtual ~KinConstraint_InvMass() {}
    KinConstraint_InvMass() {}

    KinConstraint_InvMass(std::vector<int> index_in_parts, double in_inv_mass)
    {
        index_Cons_Particles = index_in_parts;
        _nconstraints = 1;
        _inv_mass = in_inv_mass;
    }

    TVectorD getConstraint( std::vector<TLorentzVector>_p_init_vector, std::vector<TLorentzVector> _p_fin_vector) override
    {
        TVectorD _c;
        _c.ResizeTo(_nconstraints);
       
        TLorentzVector _p_1 = _p_fin_vector[index_Cons_Particles[0]];
        TLorentzVector _p_2 = _p_fin_vector[index_Cons_Particles[1]];

        _c[0]=((_p_1+_p_2).M()-_inv_mass);

        return _c;
    }

    TMatrixD getDfDx(std::vector<KinParticle> _parts) override
    {

        TLorentzVector _p_1 = _parts[index_Cons_Particles[0]];
        TLorentzVector _p_2 = _parts[index_Cons_Particles[1]];

        Double_t theta_1 = _p_1._vector.Theta();
        Double_t phi_1 = _p_1._vector.Phi();
        Double_t p_1 = _p_1._vector.P();
        Double_t E_1 = _p_1._vector.E();

        Double_t theta_2 = _p_2._vector.Theta();
        Double_t phi_2 = _p_2._vector.Phi();
        Double_t p_2 = _p_2._vector.P();
        Double_t E_2 = _p_2._vector.E();



        Double_t data[1][6] = {{ , , , , , }};

        TMatrixD dfdx(_nconstraints, 6, *data);

        dfdx *= -1.;
        return dfdx;
    }
};

#endif