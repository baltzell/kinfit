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

    TMatrixD getDfDx(int idx_part, std::vector<TLorentzVector> _particles) override 
    {
        int other_part = 0;
        if(idx_part==0)
            other_part=1;

        TLorentzVector _p_1 = _particles[idx_part];
        TLorentzVector _p_2 = _particles[other_part]; //Wrong needs to be fixed

        Double_t theta_1 = _p_1.Theta();
        Double_t phi_1 = _p_1.Phi();
        Double_t p_1 = _p_1.P();
        Double_t E_1 = _p_1.E();

        Double_t px_1 = _p_1.Px();
        Double_t py_1 = _p_1.Py();
        Double_t pz_1 = _p_1.Pz();

        Double_t theta_2 = _p_2.Theta();
        Double_t phi_2 = _p_2.Phi();
        Double_t p_2 = _p_2.P();
        Double_t E_2 = _p_2.E();

        Double_t px_2 = _p_2.Px();
        Double_t py_2 = _p_2.Py();
        Double_t pz_2 = _p_2.Pz();



        Double_t data[1][3] = {{2.*(p_1*(E_1+E_2)/E_1) -( (px_1/p_1)*(px_1+px_2) + (py_1/p_1)*(py_1+py_2) + (pz_1/p_1)*(pz_1+pz_2)), 
        -2.*( p_1*cos(phi_1)*cos(theta_1)*(px_1+px_2) + p_1*sin(phi_1)*cos(theta_1)*(py_1+py_2) -1.* p_1*sin(theta_1)*(pz_1+pz_2)), 
        -2.*(-1.*p_1*sin(theta_1)*sin(phi_1)*(px_1+px_2) + p_1*sin(theta_1)*cos(phi_1)*(py_1+py_2))}};

        TMatrixD dfdx(_nconstraints, 3, *data);

        dfdx *= -1.;
 
        return dfdx;
    }
};

#endif