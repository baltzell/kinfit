#ifndef KinConstraint_MissingMass_h
#define KinConstraint_MissingMass_h

#include "KinConstraint.h"

class KinConstraint_MissingMass : public KinConstraint
{
public:

    TLorentzVector _p_miss;
    Double_t _inv_mass;

    virtual ~KinConstraint_MissingMass() {}
    KinConstraint_MissingMass() {}

    KinConstraint_MissingMass(std::vector<int> index_in_parts, Double_t in_inv_mass)
    {
        _index_Cons_Particles = index_in_parts;
        _nconstraints = 1;
        _inv_mass = in_inv_mass;
    }

    TVectorD getConstraint(std::vector<TLorentzVector> init_particles, std::vector<TLorentzVector> fin_particles) override
    {
        TVectorD c;
        c.ResizeTo(_nconstraints);

        for (auto p_in : init_particles)
            _p_miss += p_in;
        for (auto idx : _index_Cons_Particles){
            _p_miss -= fin_particles[idx];
        }

        c[0]=_p_miss.M2()-(_inv_mass*_inv_mass);
        return c;
    }

    TMatrixD getDfDx(int idx_part, std::vector<TLorentzVector> init_particles, std::vector<TLorentzVector> fin_particles) override 
    {
    
        Int_t nvars = 3;

        TLorentzVector current_part = fin_particles[idx_part];

        Double_t theta_c = current_part.Theta();
        Double_t phi_c = current_part.Phi();
        Double_t p_c = current_part.P();
        Double_t E_c = current_part.E();

        Double_t px_c = current_part.Px();
        Double_t py_c = current_part.Py();
        Double_t pz_c = current_part.Pz();
    
        Double_t px_miss = _p_miss.Px();
        Double_t py_miss = _p_miss.Py();
        Double_t pz_miss = _p_miss.Pz();
        Double_t E_miss = _p_miss.E();


        Double_t data[1][3] = {{2.*(p_c*(E_miss)/E_c) -2.*( (px_c/p_c)*(px_miss) + (py_c/p_c)*(py_miss) + (pz_c/p_c)*(pz_miss)),
        -2.*( p_c*cos(phi_c)*cos(theta_c)*(px_miss) + p_c*sin(phi_c)*cos(theta_c)*(py_miss) -1.* p_c*sin(theta_c)*(pz_miss)),
        -2.*(-1.*p_c*sin(theta_c)*sin(phi_c)*(px_miss) + p_c*sin(theta_c)*cos(phi_c)*(py_miss))}};

        TMatrixD dfdx(_nconstraints, nvars, *data);

        dfdx *= -1.; //overall sign to match constraint definition

        return dfdx;
    }
};

#endif