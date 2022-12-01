#ifndef KinConstraint_InvMass_h
#define KinConstraint_InvMass_h

#include "KinConstraint.h"

// WARNING:  This only supports one invariant mass with two daughters.

class KinConstraint_InvMass : public KinConstraint
{
public:

    TLorentzVector _p_inv;
    double _inv_mass;

    KinConstraint_InvMass(std::vector<int> index_in_parts, double in_inv_mass)
    {
        _index_cons_particles = index_in_parts;
        _nconstraints = 1;
        _inv_mass = in_inv_mass;
    }

    TVectorD getConstraint(std::vector<TLorentzVector> init_particles, std::vector<TLorentzVector> fin_particles) override
    {
        _p_inv.SetXYZT(0,0,0,0);

        for (auto i : _index_cons_particles) _p_inv += fin_particles[i];
            
        TVectorD c;
        c.ResizeTo(_nconstraints);
        c[0] = (_p_inv.M2() - (_inv_mass * _inv_mass));

        return c;
    }

    TMatrixD getDfDx(int idx_part, std::vector<TLorentzVector> init_particles, std::vector<TLorentzVector> fin_particles) override
    {
    
        TLorentzVector current_part = fin_particles[idx_part];

        double theta_c = current_part.Theta();
        double phi_c = current_part.Phi();
        double p_c = current_part.P();
        double E_c = current_part.E();

        double px_c = current_part.Px();
        double py_c = current_part.Py();
        double pz_c = current_part.Pz();
    
        double px_inv = _p_inv.Px();
        double py_inv = _p_inv.Py();
        double pz_inv = _p_inv.Pz();
        double E_inv = _p_inv.E();

        double data[1][_nvars] = {
            {
                2*(p_c*(E_inv)/E_c) - 2*( (px_c/p_c)*(px_inv) + (py_c/p_c)*(py_inv) + (pz_c/p_c)*(pz_inv)),
               -2*( p_c*cos(phi_c)*cos(theta_c)*(px_inv) + p_c*sin(phi_c)*cos(theta_c)*(py_inv) - p_c*sin(theta_c)*(pz_inv)),
               -2*(-p_c*sin(theta_c)*sin(phi_c)*(px_inv) + p_c*sin(theta_c)*cos(phi_c)*(py_inv))
            }
        };

        TMatrixD dfdx(_nconstraints, _nvars, *data);

        return dfdx;
    }
};

#endif
