#ifndef KinConstraint_h
#define KinConstraint_h

#include "KinParticle.h"

class KinConstraint
{
public:
    virtual ~KinConstraint() {}
    KinConstraint() {}


    virtual TVectorD getConstraint(std::vector<TLorentzVector> init_particles, std::vector<TLorentzVector> fin_particles) = 0;
    virtual TMatrixD getDfDx(int idx_part, std::vector<TLorentzVector> init_particles, std::vector<TLorentzVector> fin_particles) = 0;
   

    virtual TMatrixD constructDMatrix(std::vector<TLorentzVector> init_particles, std::vector<TLorentzVector> fin_particles)
    {

        Int_t nparticles = _index_Cons_Particles.size();
        Int_t nvars = 3 * nparticles;

        TMatrixD D_Mat(_nconstraints, nvars);
        
        int idx_mat = 0;
        for ( auto idx_part : _index_Cons_Particles)
        {
            TMatrixD dfdx = this->getDfDx(idx_part, init_particles, fin_particles);
            D_Mat.SetSub(0, 3 * idx_mat, dfdx);
            idx_mat++;
        }
        
        return D_Mat;
    }


public:
    std::vector<int> _index_Cons_Particles;
    int _nconstraints;
    double _inv_mass; // there might be a way to have attributes in the derived class
};

#endif