#ifndef KinConstraint_EnergyMomentum_h
#define KinConstraint_EnergyMomentum_h

#include "KinConstraint.h"

class KinConstraint_EnergyMomentum : public KinConstraint
{
public:

    KinConstraint_EnergyMomentum(std::vector<int> index_in_parts)
    {
        _index_cons_particles=index_in_parts;
        _nconstraints = 4;
    }

    TVectorD getConstraint( std::vector<TLorentzVector>init_particles, std::vector<TLorentzVector> fin_particles) override
    {
        TVectorD c;
        c.ResizeTo(_nconstraints);
       
        TLorentzVector p_fin = TLorentzVector(0, 0, 0, 0);
        for (auto P : fin_particles)
            p_fin += P;

        TLorentzVector p_init = TLorentzVector(0, 0, 0, 0);
        for (auto P : init_particles)
            p_init += P;

        TLorentzVector p_diff = p_init - p_fin;

        for (Int_t ii = 0; ii < _nconstraints; ii++)
        {
            c[ii] = p_diff[ii];
        }

        return c;
    }

    TMatrixD getDfDx(int idx_part, std::vector<TLorentzVector> init_particles, std::vector<TLorentzVector> fin_particles) override 
    {


        Double_t theta = fin_particles[idx_part].Theta();
        Double_t phi = fin_particles[idx_part].Phi();
        Double_t p = fin_particles[idx_part].P();
        Double_t E = fin_particles[idx_part].E();

        Double_t data[4][_nvars] = {
            {sin(theta) * cos(phi),  p * cos(theta) * cos(phi), -p * sin(theta) * sin(phi)},
            {sin(theta) * sin(phi),  p * cos(theta) * sin(phi),  p * sin(theta) * cos(phi)},
            {cos(theta),            -p * sin(theta),             0},
            {p / E,                  0,                          0}
        };

        TMatrixD dfdx(_nconstraints, _nvars, *data);

        //overall sign to match constraint definition
        dfdx *= -1.;

        return dfdx;
    }

};

#endif
