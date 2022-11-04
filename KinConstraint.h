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

    virtual TMatrixD getDfDx(KinParticle _part)
    {
        TMatrixD dfdx;
        return dfdx;
    }

    virtual TMatrixD getDfDx(std::vector<KinParticle> _parts)
    {
        TMatrixD dfdx;
        return dfdx;
    }

    virtual TMatrixD constructDMatrix(std::vector<TLorentzVector> _Particles)
    {

        Int_t nparticles = index_Cons_Particles.size();
        Int_t nvars = 3 * nparticles;

        TMatrixD D_Mat(_nconstraints, nvars);
        TLorentzVector Pcurrent;

        for (Int_t jj = 0; jj < nparticles; ++jj)
        {
            Pcurrent = _Particles[index_Cons_Particles[jj]];
            TMatrixD dfdx = this->getDfDx(KinParticle(Pcurrent));
            D_Mat.SetSub(0, 3 * jj, dfdx);
        }

        return D_Mat;
    }


public:
    std::vector<int> index_Cons_Particles;
    int _nconstraints;
    double _inv_mass;
};

#endif