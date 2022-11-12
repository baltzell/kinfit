#ifndef KinConstraint_InvMass_h
#define KinConstraint_InvMass_h

#include "KinConstraint.h"

class KinConstraint_InvMass : public KinConstraint
{
public:
    virtual ~KinConstraint_InvMass() {}
    KinConstraint_InvMass() {}

    Double_t _inv_mass;

    KinConstraint_InvMass(std::vector<int> index_in_parts, Double_t in_inv_mass)
    {
        _index_Cons_Particles = index_in_parts;
        _nconstraints = 1;
        _inv_mass = in_inv_mass;
    }

    TVectorD getConstraint(std::vector<TLorentzVector> init_particles, std::vector<TLorentzVector> fin_particles) override
    {
        TVectorD c;
        c.ResizeTo(_nconstraints);

        TLorentzVector p_1 = fin_particles[_index_Cons_Particles[0]];
        TLorentzVector p_2 = fin_particles[_index_Cons_Particles[1]];

        c[0] = ((p_1 + p_2).M2() - (_inv_mass * _inv_mass));

        return c;
    }

    TMatrixD getDfDx(int idx_part, std::vector<TLorentzVector> init_particles, std::vector<TLorentzVector> fin_particles) override
    {
        int other_part = 0;
        if (idx_part == 0)
            other_part = 1;

        TLorentzVector part_1 = fin_particles[idx_part];
        TLorentzVector part_2 = fin_particles[other_part]; // This only support two particle invariant mass...

        Double_t theta_1 = part_1.Theta();
        Double_t phi_1 = part_1.Phi();
        Double_t p_1 = part_1.P();
        Double_t E_1 = part_1.E();

        Double_t px_1 = part_1.Px();
        Double_t py_1 = part_1.Py();
        Double_t pz_1 = part_1.Pz();

        Double_t theta_2 = part_2.Theta();
        Double_t phi_2 = part_2.Phi();
        Double_t p_2 = part_2.P();
        Double_t E_2 = part_2.E();

        Double_t px_2 = part_2.Px();
        Double_t py_2 = part_2.Py();
        Double_t pz_2 = part_2.Pz();

        Double_t data[1][_nvars] = {{2.*(p_1*(E_1+E_2)/E_1) -2.*( (px_1/p_1)*(px_1+px_2) + (py_1/p_1)*(py_1+py_2) + (pz_1/p_1)*(pz_1+pz_2)),
        -2.*( p_1*cos(phi_1)*cos(theta_1)*(px_1+px_2) + p_1*sin(phi_1)*cos(theta_1)*(py_1+py_2) -1.* p_1*sin(theta_1)*(pz_1+pz_2)),
        -2.*(-1.*p_1*sin(theta_1)*sin(phi_1)*(px_1+px_2) + p_1*sin(theta_1)*cos(phi_1)*(py_1+py_2))}};

        TMatrixD dfdx(_nconstraints, _nvars, *data);

        return dfdx;
    }
};

#endif
