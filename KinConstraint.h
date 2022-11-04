#ifndef KinConstraint_h
#define KinConstraint_h

#include "KinParticle.h"

class KinConstraint
{
public:
    virtual ~KinConstraint() {}
    KinConstraint() {}


    virtual TVectorD getConstraint(std::vector<TLorentzVector> _p_init, std::vector<TLorentzVector> _p_fin_vector)
    {
        TVectorD _c;
        _c.ResizeTo(_nconstraints);
        return _c;
    }

    virtual TMatrixD getDfDx(int idx_part, std::vector<TLorentzVector> _parts)
    {
        TMatrixD dfdx;
        return dfdx;
    }

    virtual TMatrixD constructDMatrix(std::vector<TLorentzVector> _particles)
    {

        Int_t nparticles = index_Cons_Particles.size();
        Int_t nvars = 3 * nparticles;

        TMatrixD D_Mat(_nconstraints, nvars);
        
        int idx_mat = 0;
        for ( auto idx_part : index_Cons_Particles)
        {
            TMatrixD dfdx = this->getDfDx(idx_part, _particles);
            D_Mat.SetSub(0, 3 * idx_mat, dfdx);
            idx_mat++;
        }
        
        return D_Mat;
    }


public:
    std::vector<int> index_Cons_Particles;
    int _nconstraints;
    double _inv_mass; // there might be a way to have attributes in the derived class
};

#endif