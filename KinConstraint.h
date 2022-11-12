#ifndef KinConstraint_h
#define KinConstraint_h

#include "KinParticle.h"

class KinConstraint
{
public:
    
    // the number of kinematic variables per particle
    static const int _nvars = 3;

    virtual ~KinConstraint() {}

    // Calculate f(x)
    virtual TVectorD getConstraint(std::vector<TLorentzVector> init_particles, std::vector<TLorentzVector> fin_particles) = 0;

    // Calculate B for a given particle
    virtual TMatrixD getDfDx(int idx_part, std::vector<TLorentzVector> init_particles, std::vector<TLorentzVector> fin_particles) = 0;

    // Concatenate Bs of each fitted particle
    TMatrixD constructBMatrix(std::vector<TLorentzVector> init_particles, std::vector<TLorentzVector> fin_particles)
    {
        Int_t nparticles = _index_cons_particles.size();
        Int_t nvars = _nvars * nparticles;

        TMatrixD d_mat(_nconstraints,nvars); 

        int idx_mat = 0;
        for (auto i : _index_cons_particles)
        {
            TMatrixD dfdx = this->getDfDx(i, init_particles, fin_particles);
            d_mat.SetSub(0, _nvars * idx_mat, dfdx);
            idx_mat++;
        }

        return d_mat;
    }

protected:
    // the index of particles used in constraints
    std::vector<int> _index_cons_particles;

    // the number of constraints
    int _nconstraints;
};

#endif
