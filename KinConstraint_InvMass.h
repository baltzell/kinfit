#ifndef KinConstraint_InvMass_h
#define KinConstraint_InvMass_h

#include "KinConstraint.h"

// WARNING:  This only supports one invariant mass with two daughters.

class KinConstraint_InvMass : public KinConstraint
{
public:

    double _inv_mass;

    KinConstraint_InvMass(std::vector<int> index_in_parts, double in_inv_mass)
    {
        _index_cons_particles = index_in_parts;
        _nconstraints = 1;
        _inv_mass = in_inv_mass;
    }

    TVectorD getConstraint(std::vector<TLorentzVector> init_particles, std::vector<TLorentzVector> fin_particles) override
    {
        TLorentzVector p_1 = fin_particles[_index_cons_particles[0]];
        TLorentzVector p_2 = fin_particles[_index_cons_particles[1]];

        TVectorD c;
        c.ResizeTo(_nconstraints);
        c[0] = ((p_1 + p_2).M2() - (_inv_mass * _inv_mass));

        return c;
    }

    TMatrixD getDfDx(int idx_part, std::vector<TLorentzVector> init_particles, std::vector<TLorentzVector> fin_particles) override
    {
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

        double data[1][_nvars] = {
            {
             2*( p_1*(E_1+E_2)/E_1) - 2*( (px_1/p_1)*(px_1+px_2) + (py_1/p_1)*(py_1+py_2) + (pz_1/p_1)*(pz_1+pz_2)),
            -2*( p_1*cos(phi_1)*cos(theta_1)*(px_1+px_2) + p_1*sin(phi_1)*cos(theta_1)*(py_1+py_2) - p_1*sin(theta_1)*(pz_1+pz_2)),
            -2*(-p_1*sin(theta_1)*sin(phi_1)*(px_1+px_2) + p_1*sin(theta_1)*cos(phi_1)*(py_1+py_2))
            }
        };

        TMatrixD dfdx(_nconstraints, _nvars, *data);

        return dfdx;
    }
};

#endif
