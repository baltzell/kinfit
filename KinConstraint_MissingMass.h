#ifndef KinConstraint_MissingMass_h
#define KinConstraint_MissingMass_h

#include "KinConstraint.h"

class KinConstraint_MissingMass : public KinConstraint
{
public:
    virtual ~KinConstraint_MissingMass() {}
    KinConstraint_MissingMass() {}

    KinConstraint_MissingMass(std::vector<int> index_in_parts, double in_inv_mass)
    {
        _index_Cons_Particles = index_in_parts;
        _nconstraints = 1;
        _inv_mass = in_inv_mass;
    }

    TVectorD getConstraint(std::vector<TLorentzVector> init_particles, std::vector<TLorentzVector> fin_particles) override
    {
        TVectorD c;
        c.ResizeTo(_nconstraints);

        TLorentzVector p_miss = TLorentzVector(0, 0, 0, 0);

        for (auto idx : _index_Cons_Particles)
            p_miss += fin_particles[idx];
        for (auto p_in : init_particles)
            p_miss += p_in;

        c[0]=p_miss.M();
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


        TLorentzVector p_miss = TLorentzVector(0, 0, 0, 0);

        for (auto idx : _index_Cons_Particles)
            p_miss -= fin_particles[idx];
        /*for (auto p_in : particles)
            p_miss += p_in;*/
        
        Double_t px_miss = p_miss.Px();
        Double_t py_miss = p_miss.Py();
        Double_t pz_miss = p_miss.Pz();
        Double_t E_miss = p_miss.E();


        Double_t data[1][3];
        /* Double_t data[1][3] = {{2.*(p_c*(E_miss)/E_1) -2.*( (px_1/p_1)*(px_miss) + (py_1/p_1)*(py_1+py_2) + (pz_1/p_1)*(pz_1+pz_2)),
        -2.*( p_1*cos(phi_1)*cos(theta_1)*(px_1+px_2) + p_1*sin(phi_1)*cos(theta_1)*(py_1+py_2) -1.* p_1*sin(theta_1)*(pz_1+pz_2)),
        -2.*(-1.*p_1*sin(theta_1)*sin(phi_1)*(px_1+px_2) + p_1*sin(theta_1)*cos(phi_1)*(py_1+py_2))}};*/

        TMatrixD dfdx(_nconstraints, nvars, *data);


        return dfdx;
    }
};

#endif