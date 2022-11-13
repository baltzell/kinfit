
#ifndef KinConstraint_MissingMass_h
#define KinConstraint_MissingMass_h

#include "KinConstraint.h"

// WARNING:  This only supports one invariant mass with two daughters.
// WARNING:  A naive merging of KC_MissingMass and KC_InvMass ...

class KinConstraint_MissAndInvMass : public KinConstraint
{
public:

    TLorentzVector _p_miss;
    double _inv_mass;
    double _miss_mass;

    KinConstraint_MissingMass(std::vector<int> index_in_parts, double in_inv_mass, double in_miss_mass)
    {
        _index_cons_particles = index_in_parts;
        _nconstraints = 2;
        _inv_mass = in_inv_mass;
        _miss_mass = in_miss_mass;
    }

    TVectorD getConstraint(std::vector<TLorentzVector> init_particles, std::vector<TLorentzVector> fin_particles) override
    {
        // missing mass stuff:
        _p_miss.SetXYZT(0,0,0,0);
        for (auto p : init_particles) _p_miss += p;
        for (auto i : _index_cons_particles) _p_miss -= fin_particles[i];

        // the first two particles are the invariant mass constraint:
        TLorentzVector p_1 = fin_particles[_index_cons_particles[0]];
        TLorentzVector p_2 = fin_particles[_index_cons_particles[1]];

        TVectorD c;
        c.ResizeTo(_nconstraints);

        // the missing mass constraint:
        c[0]=_p_miss.M2()-(_miss_mass*_miss_mass);

        // the invariant mass constraint:
        c[1] = ((p_1 + p_2).M2() - (_inv_mass * _inv_mass));

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

        const int other_part = idx_part==0 ? 1 : 0;;

        TLorentzVector part_1 = fin_particles[idx_part];
        TLorentzVector part_2 = fin_particles[other_part];

        double theta_1 = part_1.Theta();
        double phi_1 = part_1.Phi();
        double p_1 = part_1.P();
        double E_1 = part_1.E();

        double px_1 = part_1.Px();
        double py_1 = part_1.Py();
        double pz_1 = part_1.Pz();

        double theta_2 = part_2.Theta();
        double phi_2 = part_2.Phi();
        double p_2 = part_2.P();
        double E_2 = part_2.E();

        double px_2 = part_2.Px();
        double py_2 = part_2.Py();
        double pz_2 = part_2.Pz();

        double data[2][_nvars] = {
            {
                2*(p_c*(E_miss)/E_c) - 2*( (px_c/p_c)*(px_miss) + (py_c/p_c)*(py_miss) + (pz_c/p_c)*(pz_miss)),
               -2*( p_c*cos(phi_c)*cos(theta_c)*(px_miss) + p_c*sin(phi_c)*cos(theta_c)*(py_miss) - p_c*sin(theta_c)*(pz_miss)),
               -2*(-p_c*sin(theta_c)*sin(phi_c)*(px_miss) + p_c*sin(theta_c)*cos(phi_c)*(py_miss))
            },
            {
                -1*(2*( p_1*(E_1+E_2)/E_1) - 2*( (px_1/p_1)*(px_1+px_2) + (py_1/p_1)*(py_1+py_2) + (pz_1/p_1)*(pz_1+pz_2))),
                -1*(-2*( p_1*cos(phi_1)*cos(theta_1)*(px_1+px_2) + p_1*sin(phi_1)*cos(theta_1)*(py_1+py_2) - p_1*sin(theta_1)*(pz_1+pz_2))),
                -1*(-2*(-p_1*sin(theta_1)*sin(phi_1)*(px_1+px_2) + p_1*sin(theta_1)*cos(phi_1)*(py_1+py_2)))
            }
        };

        TMatrixD dfdx(_nconstraints, _nvars, *data);

        // overall sign to match constraint definition:
        dfdx *= -1.;

        return dfdx;
    }
};

#endif
