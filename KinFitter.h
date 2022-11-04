#ifndef KinFitter_h
#define KinFitter_h 1

#include "Rtypes.h"
#include "TObject.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TLorentzVector.h"
#include "TH1D.h"

#include "KinParticle.h"
#include "KinConstraint_EnergyMomentum.h"
#include "KinConstraint_InvMass.h"




class KinFitter
{
    // For now, these methods use the nomenclature in Mike Williams' 2004 analysis note

public:
    double GetConfidenceLevel() { return _confLevel; }
    double GetChi2() { return _chi2; }
    double GetNDF() { return _ndf; }
    TVectorD GetPulls() { return _pulls; }
    std::vector<TLorentzVector> GetFitted4Vectors() { return _Ps_y; }

private:
    Double_t _confLevel;
    TVectorD _pulls;
    Double_t _chi2;
    Int_t _nvars_y;
    Int_t _nconstraints;
    Int_t _ndf;

    int nter = 100;

    TVectorD _eta;
    TVectorD _sigma2_etas;
    TMatrixD _C_eta;
    TVectorD _y;
    TVectorD _sigma2_ys;
    TMatrixD _C_y;
    TVectorD _delta;
    TVectorD _epsilon;

    std::vector<TLorentzVector> _P_inits;
    std::vector<Double_t> _masses_y;   // Probably won't need
    std::vector<TLorentzVector> _Ps_y; // Probably the only things that need to be public

    TMatrixD _A;
    TMatrixD _B;

    TVectorD _ccccc;
    TLorentzVector _P_init;
    TLorentzVector _P_fin;
    TLorentzVector _Pdiff;

    TMatrixD _C_B;
    TVectorD _mu;
    TMatrixD _C_x;

    std::vector<KinConstraint*> _Cons;   

public:
    virtual ~KinFitter() {}

    KinFitter() {}

public:
    void SetInitial(std::vector<KinParticle> in_P_initial)
    {
        for (auto& in_P : in_P_initial)
        {
            _P_inits.push_back(in_P.GetVector());
        }
    }

   void SetFinal(std::vector<KinParticle> in_P_final)
    {
        _nvars_y = 3 * in_P_final.size();
        TMatrixD C_n(_nvars_y, _nvars_y);

        int idx = 0;
        for (auto& in_P : in_P_final)
        {
            _Ps_y.push_back(in_P.GetVector());
            TMatrixD C_in_P = in_P.GetCovMatrix();
            TMatrixDSub(C_n,0+3*idx,2+3*idx,0+3*idx,2+3*idx) = C_in_P; //Create block diagonal matrix with covariances

            idx++;

        }

        

        _ndf = 4;

        _y.ResizeTo(_nvars_y);
        _eta.ResizeTo(_nvars_y);

        _sigma2_ys.ResizeTo(_nvars_y);

        _C_eta.ResizeTo(_nvars_y, _nvars_y);
        _C_eta = C_n;

        TMatrixDDiag diag(_C_eta);
        _sigma2_etas.ResizeTo(_nvars_y);
        _sigma2_etas = diag;

        Init4Vects();
    }

    void Add_EnergyMomentum_Constraint(std::vector<int> index_P_Cons)
    {
            KinConstraint_EnergyMomentum *in_Cons = new KinConstraint_EnergyMomentum(index_P_Cons);
            _Cons.push_back(in_Cons);    
            _nconstraints = 4;        
    }

    void Add_InvMass_Constraint(std::vector<int> index_P_Cons, double in_mass)
    {
            KinConstraint_InvMass *in_Cons = new KinConstraint_InvMass(index_P_Cons, in_mass);
            _Cons.push_back(in_Cons); 
            _nconstraints = 1;           
    }

    // Main, Default Fitter
    void DoFitting(int in_nter=100)
    {
        nter=in_nter;

        TMatrixD *AA = nullptr; // To be address of A
        Int_t iter = 0;
        while (iter < nter)
        { // ... or some other conditions; I ain't your boss.
            SetB();
            SetC();

            ProcessFit(&_B, AA, &_ccccc);
            _Ps_y = get4Vectors(&_y, _masses_y);

            _epsilon += _delta;
            Double_t chi20 = _chi2;
            SetChi2();

            // if( abs(_chi2 - chi20)/chi20 < 0.01 ){
            if (_chi2 - chi20 > 0.)
            {
                if (iter > 0)
                {
                    break;
                }
            }
            
            ++iter;
        }

        PostProcess();
    }

private:
    void Initialize(std::vector<Double_t> sigmas, TMatrixD *C_n = nullptr, TMatrixD *corrMat = nullptr)
    {
        _ndf = abs(_nconstraints);

        _y.ResizeTo(_nvars_y);
        _eta.ResizeTo(_nvars_y);

        _sigma2_ys.ResizeTo(_nvars_y);
        _C_eta.ResizeTo(_nvars_y, _nvars_y);

        // Set up Covariance matrix
        if (C_n == 0)
        {
            // If a covariance matrix is not given, build the covariance matrix from sigmas
            InitErrors(sigmas, corrMat);
        }
        else
        {
            _C_eta = *C_n;
            TMatrixDDiag diag(_C_eta);

            // Get variance from diagonal of covariance matrix
            _sigma2_etas.ResizeTo(_nvars_y);
            _sigma2_etas = diag;
        }
    }

    void InitErrors(std::vector<Double_t> sigmas, TMatrixD *corrMat)
    {
        // Initialize Errors
        TVectorD sigma_etas(sigmas.size());
        for (Int_t ii = 0; ii < sigmas.size(); ++ii)
        {
            sigma_etas[ii] = sigmas[ii];
        }

        _sigma2_etas.ResizeTo(_nvars_y);
        _sigma2_etas = sigma_etas;
        _sigma2_etas.Sqr();

        SetCovMat(corrMat);
    }

    void SetCovMat(TMatrixD *corrMat = nullptr)
    {
        if (corrMat == 0)
        {
            //// Variance is sigma^2
            TMatrixDDiag diag(_C_eta);
            diag = _sigma2_etas;
        }
        else
        {
            // I don't think there's a way to do vec^T * Mat * vec so here's a loop:
            TMatrixD rho = *corrMat;
            TVectorD sigmas = _sigma2_etas;
            sigmas.Sqrt();
            Int_t n = _sigma2_etas.GetNrows();
            for (Int_t ii = 0; ii < n; ++ii)
            {
                for (Int_t jj = 0; jj < n; ++jj)
                {
                    _C_eta[ii][jj] = sigmas[ii] * rho[ii][jj] * sigmas[jj];
                }
            }
        }
    }

    void Init4Vects()
    {
        // Initial 4-Vectors
        for (auto P : _P_inits)
            _P_init += P;
        for (auto P : _Ps_y)
            _masses_y.push_back(P.M());

        _y.ResizeTo(3 * _Ps_y.size());
        _eta.ResizeTo(_y);

        _y = constructTVecFrom4Vecs(_Ps_y);
        _eta = _y;
        _epsilon.ResizeTo(_eta);
    }

    // ---------------------------------------------------------------------------------------------------------------------------------
    // Methods
    // ---------------------------------------------------------------------------------------------------------------------------------

    // ---------------------------------------------------------------------------------------------------------------------------------
    // Setters : Construct needed Vector and Matrices
    // ---------------------------------------------------------------------------------------------------------------------------------

    void SetDelta()
    {
        // Construct delta ( = -_C_eta B^T C_B(  c + A xi ) = -_C_eta B^T * mu )
        _delta.ResizeTo(_y);

        TMatrixD BT = _B;
        BT.T();

        TMatrixD D = _C_eta * BT;

        SetMu();

        _delta = D * _mu;
        _delta *= -1;
    }

    void SetMu()
    {
        // Construct mu ( = C_B(  c + A xi )  )
        _mu.ResizeTo(_ccccc);
        _mu = _ccccc;

        SetC_B();

        _mu *= _C_B;
    }

    void SetC_B()
    {
        // Maybe implement a similarity transformation : B C B^T and B^T C B
        // Construct Matrix C_B ( = inverse(B _C_eta B^T) )
        TMatrixD BT = _B;
        BT.T();

        _C_B.ResizeTo(_B.GetNrows(), _B.GetNrows());
        _C_B = _B * _C_eta * BT;
        _C_B.SetTol(1.e-30); // Set tolerance to be much lower for inverted matrix
        _C_B.Invert();
    }

    void SetC_y()
    {
        // Returns the new covariance matrix from propagation of errors
        // Get Error in fit results: C_y ( = _C_eta - _C_eta B^T C_B B _C_eta + _C_eta B^T C_B A inverse( A^T C_B A ) A^T C_B B _C_eta = _C_eta - _C_eta B^T C_B B _C_eta + _C_eta B^T C_B A C_x A^T C_B B _C_eta)
        _C_y.ResizeTo(_C_eta);

        TMatrixD BT = _B;
        BT.T();

        TMatrixD C_y1 = _C_eta * BT * _C_B;
        TMatrixD C_y2 = C_y1 * _B * _C_eta;

        // Documentation recommended to add matrices this way
        _C_y = _C_eta;
        _C_y -= C_y2;

        if (_A.GetNrows() != 0)
        {
            TMatrixD AT = _A;
            _A.T();
            TMatrixD C_y3 = C_y1 * _A * _C_x * AT * _C_B * _B * _C_eta;
            _C_y += C_y3;
        }
    }

    void SetB(TMatrixD *BB)
    {
        _B.ResizeTo(*BB);
        _B = *BB;
    }

    void SetC(TVectorD *cc)
    {
        _ccccc.ResizeTo(*cc);
        _ccccc = *cc;
    }

    void SetB()
    {
        _B.ResizeTo(_nconstraints, _nvars_y);
        _B=_Cons[0]->constructDMatrix(get4Vectors(&_y, _masses_y));//(&_y, _masses_y);

        //_B = constructDerMatrix(&_y, _masses_y);
    }

    void SetC()
    {
        _ccccc.ResizeTo(_nconstraints);

        _ccccc=_Cons[0]->getConstraint(_P_inits,_Ps_y);

       
    }

    void ProcessFit(TMatrixD *BB, TMatrixD *AA, TVectorD *cc)
    {
        _ccccc.ResizeTo(*cc);
        _B.ResizeTo(*BB);

        _ccccc = *cc; // Just in case these come from outside the class
        _B = *BB;     // Just in case these come from outside the class

        SetDelta();

        _y += _delta;
    }

    void PostProcess()
    {
        // Get change of new fit results (y):  epsilon ( eta - y ) !! Probably won't need
        //_epsilon -= y;

        SetChi2();

        SetConfLevel();

        // Correlations Matrix
        SetC_y();
        _sigma2_ys = TMatrixDDiag(_C_y);

        SetPulls();
    }

    // ---------------------------------------------------------------------------------------------------------------------------------
    // Postfit
    // ---------------------------------------------------------------------------------------------------------------------------------

    void SetChi2()
    {
        // For a kinematic fit, this should follow a chi^2 distribution
        // Construct chi^2 ( = epsilon^T inverse(_C_eta) epsilon = epsilon . epsilon')

        TVectorD eps0 = _epsilon;

        TMatrixD C_nI = _C_eta;
        C_nI.SetTol(1.e-30); // Set tolerance to be much lower for inverted matrix
        C_nI.Invert();

        TVectorD eps1 = C_nI * eps0;

        _chi2 = eps0 * eps1;
    }

    void SetConfLevel()
    {
        // Returns confidence level for a given chi^2 (chi2) and number of degrees of freedom (ndf)
        _confLevel = TMath::Prob(_chi2, _ndf);
    }

    void SetPulls()
    {
        // Returns a vector of pull distributions
        // Get pulls z_i ( = [ eta_i - y_i ] / sqrt( sigma^2_eta_i - sigma^2_y_i )
        _pulls.ResizeTo(_eta);

        TVectorD num = _epsilon;

        //cout<<_sigma2_etas[0]<<" "<<_sigma2_etas[1]<<" "<<_sigma2_etas[2]<<" "<<_sigma2_etas[3]<<" "<<endl;
        //_sigma2_etas.Draw();

        TVectorD denom = _sigma2_etas - _sigma2_ys;
        // TVectorD denom = _sigma2_etas;

        denom.Abs();  // !!! Probably don't need! In case sigmas2_y < sigmas2_n
        denom.Sqrt(); // Element-wise square root!

        _pulls = ElementDiv(num, denom); // Element-wise division!
    }

    // ---------------------------------------------------------------------------------------------------------------------------------
    // These methods are on the bottom because they are very specfic. The rest should be pretty general
    // ---------------------------------------------------------------------------------------------------------------------------------

    std::vector<TLorentzVector> get4Vectors(TVectorD *v, std::vector<Double_t> masses) //Useless function, need to be removed
    {
        std::vector<TLorentzVector> result;
        Int_t nparticles = masses.size();
        for (Int_t ii = 0; ii < nparticles; ++ii)
        {
            Double_t m = masses[ii];
            Double_t p = (*v)[3 * ii];
            Double_t theta = (*v)[3 * ii + 1];
            Double_t phi = (*v)[3 * ii + 2];
            // Double_t E     = v[ 3*ii + 3];
            TVector3 vec3;
            vec3.SetMagThetaPhi(1, theta, phi);
            vec3 *= p;
            result.push_back(TLorentzVector(vec3, sqrt(p * p + m * m)));
        }
        return result;
    }

    TVectorD constructTVecFrom4Vecs(std::vector<TLorentzVector> P4s)
    {
        Int_t nparticles = P4s.size();
        Int_t nvars = 3 * nparticles;
        TVectorD vec(nvars);
        for (Int_t ii = 0; ii < nparticles; ++ii)
        {
            TLorentzVector P = P4s[ii];
            vec[3 * ii + 0] = P.P();
            vec[3 * ii + 1] = P.Theta();
            vec[3 * ii + 2] = P.Phi();
        }
        return vec;
    }


    ClassDef(KinFitter, 1)
};

#endif
