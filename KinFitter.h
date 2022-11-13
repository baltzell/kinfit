#ifndef KinFitter_h
#define KinFitter_h

#include "TVectorD.h"
#include "TMatrixD.h"
#include "TLorentzVector.h"
#include "TH1D.h"

#include "KinParticle.h"
#include "KinConstraint_EnergyMomentum.h"
#include "KinConstraint_InvMass.h"
#include "KinConstraint_MissingMass.h"


/////////////////////////////////////////////////////////////////////////////////////////////////////
//Kinfitter class
//This class implemtent the kinematic fitter describe in https://www.jlab.org/Hall-B/notes/clas_notes03/03-017.pdf
//The variable names match the notation used in the note and are refered to using the equation number of the note.
/////////////////////////////////////////////////////////////////////////////////////////////////////

class KinFitter
{

public:
    double GetConfidenceLevel() { return _confLevel; }
    double GetChi2() { return _chi2; }
    double GetNDF() { return _ndf; }
    TVectorD GetPulls() { return _pulls; }
    std::vector<TLorentzVector> GetFitted4Vectors() { return _Ps_y; }

private:
    double _confLevel; 
    TVectorD _pulls;
    double _chi2;
    int _nvars_y;
    int _nconstraints;
    int _ndf;

    int n_iter = 100;

    TVectorD _eta; //Measured vector Eq.1
    TVectorD _sigma2_etas; //Measured errors for each measured quantities Eq.19
    TMatrixD _C_eta; //Measured covariance matrix
    TMatrixD _C_eta_Inv; //Inverse of the measured covariance matrix, used for the #chi^2 calculation
    TVectorD _y; //Exact vector Eq.1
    TVectorD _sigma2_ys;  //Diagonal of the covariance matrix of the improved measurement Eq.16
    TMatrixD _C_y; //Covariance matrix of the improved measurement Eq.16
    TVectorD _delta; //Deviation from the prefit vector Eq.12
    TVectorD _epsilon; //Deviation from the measured vector Eq.15

    std::vector<TLorentzVector> _P_inits; //Store the 4-Vectors of the initial state particles
    std::vector<double> _masses_y; //Store the masses of the fitted particles
    std::vector<TLorentzVector> _Ps_y; //Store the 4-Vectors of the fitted particles

    TMatrixD _B; //Derivative matrix of the constraints Eq.5
    TMatrixD _previous_B; //Store the previous iteration of the derivative matrix of the constraints Eq.5

    TVectorD _c; //Constraint vector Eq.6

    TMatrixD _C_B; //Useful C_B matrix, defined in Eq.12
    TMatrixD _previous_C_B; //Store the previous iteration of the C_B matrix

    std::vector<KinConstraint *> _Cons; //Vector of contraints

    double alpha = 1.0; //Fit iteration weight

public:
    ~KinFitter() {};

    KinFitter(std::vector<KinParticle> in_P_initial, std::vector<KinParticle> in_P_final)
    {
        for (auto &in_P : in_P_initial)
        {
            _P_inits.push_back(in_P.GetVector());
        }

        // This might need some change once many constraints are set:
        _nvars_y = 3 * in_P_final.size();
        
        TMatrixD C_n(_nvars_y, _nvars_y);

        int idx = 0;
        for (auto &in_P : in_P_final)
        {
            _Ps_y.push_back(in_P.GetVector());
            _masses_y.push_back(in_P.GetMass());
            TMatrixD C_in_P = in_P.GetCovMatrix();

            // Create block diagonal matrix with covariances:
            TMatrixDSub(C_n, 0 + 3 * idx, 2 + 3 * idx, 0 + 3 * idx, 2 + 3 * idx) = C_in_P;

            idx++;
        }

        _y.ResizeTo(_nvars_y);
        _eta.ResizeTo(_nvars_y);
        _delta.ResizeTo(_nvars_y);
        _epsilon.ResizeTo(_nvars_y);
        _pulls.ResizeTo(_nvars_y);
        _sigma2_ys.ResizeTo(_nvars_y);
        _C_eta.ResizeTo(_nvars_y, _nvars_y);
        _C_eta_Inv.ResizeTo(_nvars_y, _nvars_y);
        _sigma2_etas.ResizeTo(_nvars_y);
        _C_y.ResizeTo(_C_eta);

        _C_eta = C_n;

        _C_eta_Inv = _C_eta;
        _C_eta_Inv.SetTol(1.e-30);
        _C_eta_Inv.Invert();

        TMatrixDDiag diag(_C_eta);
        _sigma2_etas = diag;

        _y = constructTVecFrom4Vecs(_Ps_y);
        _eta = _y;
    }
    
    void Add_EnergyMomentum_Constraint(std::vector<int> index_P_Cons)
    {
        KinConstraint_EnergyMomentum *in_Cons = new KinConstraint_EnergyMomentum(index_P_Cons);
        _Cons.push_back(in_Cons);
        _nconstraints = 4;
        _ndf = 4;
    }

    void Add_InvMass_Constraint(std::vector<int> index_P_Cons, double in_mass)
    {
        KinConstraint_InvMass *in_Cons = new KinConstraint_InvMass(index_P_Cons, in_mass);
        _Cons.push_back(in_Cons);
        _nconstraints = 1;
        _ndf = 1;
    }

    void Add_MissingMass_Constraint(std::vector<int> index_P_Cons, double in_mass)
    {
        KinConstraint_MissingMass *in_Cons = new KinConstraint_MissingMass(index_P_Cons, in_mass);
        _Cons.push_back(in_Cons);
        _nconstraints = 1;
        _ndf = 1;
    }

    void DoFitting(int in_nter = 100)
    {
        _chi2 = 10000;

        n_iter = in_nter;
        int iter = 0;
        int iter_inc = 0;
        while (iter < n_iter)
        {

            double chi20 = _chi2;

            /////////////////////////////////////////
            //The math is done here
            /////////////////////////////////////////
            ProcessFit();

            /////////////////////////////////////////
            // Stop conditions of the fit
            /////////////////////////////////////////
            if (_chi2 - chi20 > 0.)
                iter_inc++;

            if (_chi2 - chi20 < 0. && iter > 0)
                iter_inc = 0;

            if (iter_inc > 1)
            {
                if (iter > 0)
                {
                    UndoFit();
                    break;
                }
            }

            if (abs(_chi2 - chi20) / chi20 < 0.001)
            {
                if (iter > 0)
                {
                    break;
                }
            }

            iter++;
        }

        /////////////////////////////////////////
        //Compute the outputs: Pulls, Confidence levels, and fitted vectors
        /////////////////////////////////////////
        PostProcess();
    }

private:
    void ProcessFit()
    {

        _B.ResizeTo(_nconstraints, _nvars_y); // should be moved outside this method when ready
        _previous_B.ResizeTo(_nconstraints, _nvars_y);
        _C_B.ResizeTo(_nconstraints, _nconstraints);
        _previous_C_B.ResizeTo(_nconstraints, _nconstraints);
        _c.ResizeTo(_nconstraints);

        // Store the previous fit parameters
        _previous_C_B = _C_B;
        _previous_B = _B;

        //Compute constraints and derivatives Eq.5 and 6
        _c = _Cons[0]->getConstraint(_P_inits, get4Vectors(&_y, _masses_y)); //As to be calculated before the B matrix
        _B = _Cons[0]->constructBMatrix(_P_inits, get4Vectors(&_y, _masses_y)); // should be moved outside this method when ready

        TMatrixD BT = _B;
        BT.T();

        _C_B = _B * _C_eta * BT;
        _C_B.SetTol(1.e-30); // Set tolerance to be much lower for inverted matrix
        _C_B.Invert();

        _delta = _C_eta * BT * _C_B * _c; //Compute #delta Eq. 12
        _delta *= -1;
        _y += alpha * _delta; //Update y Eq. 14
        _epsilon += _delta; //Update #epsilon Eq. 15

        TVectorD pre_chi2 = _C_eta_Inv * _epsilon;
        _chi2 = _epsilon * pre_chi2;
    }

    void UndoFit() //This method uses stored fit parameters to undo previous fit if the chi_2 increased in the previous iteration
    {
        _y -= alpha * _delta;
        _epsilon -= _delta;

        TVectorD pre_chi2 = _C_eta_Inv * _epsilon;
        _chi2 = _epsilon * pre_chi2;

        _C_B = _previous_C_B;

        _B = _previous_B;
    }

    void PostProcess()
    {
       // Compute Confidence level Eq.17
        _confLevel = TMath::Prob(_chi2, _ndf);

        // Compute pulls
        TMatrixD BT = _B;
        BT.T();
        _C_y = _C_eta;
        _C_y -= _C_eta * BT * _C_B * _B * _C_eta; //Eq.16

        //////////
        //Eq.19
        //////////
        TVectorD num = _epsilon;
        TMatrixDDiag diag(_C_y);
        _sigma2_ys = diag;
        TVectorD denom = _sigma2_etas - _sigma2_ys;
        denom.Sqrt();
        _pulls = ElementDiv(num, denom); 

        // Set final vector - should be changed at some point
        _Ps_y = get4Vectors(&_y, _masses_y);
    }


/////////////
//Some function handling 4Vectors
//We should get ride of them at some point
/////////////
    std::vector<TLorentzVector> get4Vectors(TVectorD *v, std::vector<double> masses)
    {
        std::vector<TLorentzVector> result;
        int nparticles = masses.size();
        for (int ii = 0; ii < nparticles; ++ii)
        {
            double m = masses[ii];
            double p = (*v)[3 * ii];
            double theta = (*v)[3 * ii + 1];
            double phi = (*v)[3 * ii + 2];
            // double E     = v[ 3*ii + 3];
            TVector3 vec3;
            vec3.SetMagThetaPhi(1, theta, phi);
            vec3 *= p;
            result.push_back(TLorentzVector(vec3, sqrt(p * p + m * m)));
        }
        return result;
    }

    TVectorD constructTVecFrom4Vecs(std::vector<TLorentzVector> P4s)
    {
        int nparticles = P4s.size();
        int nvars = 3 * nparticles;
        TVectorD vec(nvars);
        for (int ii = 0; ii < nparticles; ++ii)
        {
            TLorentzVector P = P4s[ii];
            vec[3 * ii + 0] = P.P();
            vec[3 * ii + 1] = P.Theta();
            vec[3 * ii + 2] = P.Phi();
        }
        return vec;
    }
};

#endif
