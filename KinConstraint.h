#ifndef KinConstraint_h
#define KinConstraint_h

#include "KinParticle.h"
#include <iostream>
using namespace std;

class KinConstraint
{
public:
    std::vector<int> _index_Cons_Particles; //Store only the index of particles used in the constraint
    int _nconstraints;

    virtual ~KinConstraint() {}
    KinConstraint() {}

    // Fully virtual functions to be implememted in each constraint class
    virtual TVectorD getConstraint(std::vector<TLorentzVector> init_particles, std::vector<TLorentzVector> fin_particles) = 0; // Calculate f(x)
    virtual TMatrixD getDfDx(int idx_part, std::vector<TLorentzVector> init_particles, std::vector<TLorentzVector> fin_particles) = 0; // Calculate B for a given particle

    TMatrixD constructBMatrix(std::vector<TLorentzVector> init_particles, std::vector<TLorentzVector> fin_particles) // Concatenate Bs of each fitted particle
    {

        Int_t nparticles = _index_Cons_Particles.size();
        Int_t nvars = 3 * nparticles;

        TMatrixD D_Mat(_nconstraints, nvars);

        int idx_mat = 0;
        for (auto idx_part : _index_Cons_Particles)
        {
            TMatrixD dfdx = this->getDfDx(idx_part, init_particles, fin_particles);
            D_Mat.SetSub(0, 3 * idx_mat, dfdx);

            idx_mat++;
        }

        return D_Mat;
    }
};

#endif