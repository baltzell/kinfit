#ifndef KinConstraint_MissingMass_h
#define KinConstraint_MissingMass_h

#include "KinConstraint.h"

class KinConstraint_MissingMass : public KinConstraint
{
public:

    TLorentzVector _p_miss;
    double _inv_mass;

    KinConstraint_MissingMass(std::vector<int> index_in_parts, double in_inv_mass)
    {
        _index_cons_particles = index_in_parts;
        _nconstraints = 1;
        _inv_mass = in_inv_mass;
    }

    TVectorD getConstraint(std::vector<TLorentzVector> init_particles, std::vector<TLorentzVector> fin_particles) override
    {
        _p_miss.SetXYZT(0,0,0,0);

        for (auto p : init_particles) _p_miss += p;
        for (auto i : _index_cons_particles) _p_miss -= fin_particles[i];

        TVectorD c;
        c.ResizeTo(_nconstraints);
        c[0]=_p_miss.M2()-(_inv_mass*_inv_mass);

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
    
        double px_miss = _p_miss.Px();
        double py_miss = _p_miss.Py();
        double pz_miss = _p_miss.Pz();
        double E_miss = _p_miss.E();

        double data[1][_nvars] = {
            {
                2*(p_c*(E_miss)/E_c) - 2*( (px_c/p_c)*(px_miss) + (py_c/p_c)*(py_miss) + (pz_c/p_c)*(pz_miss)),
               -2*( p_c*cos(phi_c)*cos(theta_c)*(px_miss) + p_c*sin(phi_c)*cos(theta_c)*(py_miss) - p_c*sin(theta_c)*(pz_miss)),
               -2*(-p_c*sin(theta_c)*sin(phi_c)*(px_miss) + p_c*sin(theta_c)*cos(phi_c)*(py_miss))
            }
        };

        TMatrixD dfdx(_nconstraints, _nvars, *data);

        // overall sign to match constraint definition:
        dfdx *= -1.;

        return dfdx;
    }
};

#endif
