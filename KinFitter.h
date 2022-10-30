#ifndef KinFitter_h
#define KinFitter_h 1

#include "Rtypes.h"
#include "TObject.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TLorentzVector.h"
#include "TH1D.h"

class KinFitter{
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
        Int_t _nvars_x;
        Int_t _nvars_y;
        Int_t _nconstraints;
        Int_t _nconstraints_extra;
        Int_t _ndf;

        TVectorD _eta;
        TVectorD _sigma2_etas;
        TMatrixD _C_eta;
        TVectorD _y;
        TVectorD _sigma2_ys;
        TMatrixD _C_y;
        TVectorD _x;
        TVectorD _xi;
        TVectorD _delta;
        TVectorD _epsilon;

        std::vector<TLorentzVector> _P_inits;      
        std::vector<Double_t> _masses_x;         // Probably won't need 
        std::vector<Double_t> _masses_y;         // Probably won't need
        std::vector<TLorentzVector> _Ps_y;       //Probably the only things that need to be public 
        std::vector<TLorentzVector> _Ps_x;       //Probably the only things that need to be public 

        TMatrixD _A;
        TMatrixD _B;

        TVectorD _ccccc;
        TLorentzVector _P_init;
        TLorentzVector _P_fin;
        TLorentzVector _Pdiff;

        TMatrixD _C_B;
        TVectorD _mu;
        TMatrixD _C_x;

    public: 
        // ---------------------------------------------------------------------------------------------------------------------------------
        // Ctors and Dtors 
        // ---------------------------------------------------------------------------------------------------------------------------------
        virtual ~KinFitter(){}    

        KinFitter(){}

        // Default Constructor with only 4-Momenta conservation (if nconstraints_extra is 0)
        // Constructor with extra constraints (other than energy and momentum conservation... like an invariant mass)
        // Use this if you have only variances
        KinFitter( std::vector<std::vector<TLorentzVector>> P_InitsFins, std::vector<TLorentzVector> P_Fins_x, std::vector<Double_t> sigmas, Int_t ncons_extra = 0, TMatrixD * C_n = nullptr )
            : 
                _nconstraints_extra( ncons_extra               ),
                _nconstraints(       4 + ncons_extra           ),
                _nvars_x(            3 * P_Fins_x.size()       ), 
                _nvars_y(            3 * P_InitsFins[1].size() ), 
                _Ps_x(               P_Fins_x                  ),
                _Ps_y(               P_InitsFins[1]            ), 
                _P_inits(            P_InitsFins[0]            ) 
        {
            Initialize(sigmas, C_n);
            Init4Vects();  // y and eta defined here

            // Loop until result is satisfactory
            DoFitting(100);
        }

        // Use this constructor if you have custom { x, y, eta, A, B, c, sigmas ... }
        KinFitter( Int_t ncons, std::vector<Double_t> y0s, std::vector<Double_t> x0s, std::vector<Double_t> sigmas, TMatrixD* C_n = nullptr )
            : 
                _nconstraints( ncons      ),
                _nvars_x(      x0s.size() ), 
                _nvars_y(      y0s.size() ) 
        {
            Initialize(sigmas, C_n);

            for( Int_t ii = 0; ii < _nvars_x; ++ii ){ _x[ii]   = x0s[ii]; }
            for( Int_t ii = 0; ii < _nvars_y; ++ii ){ _eta[ii] = y0s[ii]; _y[ii] = y0s[ii];}
            _epsilon.ResizeTo(_eta);
        }

    private:
        void Initialize( std::vector<Double_t> sigmas, TMatrixD *C_n = nullptr, TMatrixD *corrMat = nullptr ){
            _ndf = abs(_nconstraints - _nvars_x); 

            _x.ResizeTo(   _nvars_x );
            _y.ResizeTo(   _nvars_y );
            _eta.ResizeTo( _nvars_y );

            _sigma2_ys.ResizeTo(   _nvars_y );
            _C_eta.ResizeTo( _nvars_y, _nvars_y ); 

            // Set up Covariance matrix
            if( C_n == 0 ){
                // If a covariance matrix is not given, build the covariance matrix from sigmas
                InitErrors(sigmas, corrMat);
            }
            else{
                _C_eta = *C_n;
                TMatrixDDiag diag(_C_eta);

                // Get variance from diagonal of covariance matrix
                _sigma2_etas.ResizeTo( _nvars_y ); 
                _sigma2_etas = diag;
            }
        }

        void InitErrors(std::vector<Double_t> sigmas, TMatrixD* corrMat){
            // Initialize Errors
            TVectorD sigma_etas( sigmas.size() );
            for( Int_t ii = 0; ii < sigmas.size(); ++ii ){ sigma_etas[ii] = sigmas[ii]; }

            _sigma2_etas.ResizeTo( _nvars_y ); 
            _sigma2_etas = sigma_etas;
            _sigma2_etas.Sqr();

            SetCovMat(corrMat);
        }

        void SetCovMat(TMatrixD* corrMat = nullptr){
            if( corrMat == 0 ){
                //// Variance is sigma^2
                TMatrixDDiag diag(_C_eta); 
                diag = _sigma2_etas;
            }
            else{
                // I don't think there's a way to do vec^T * Mat * vec so here's a loop:
                TMatrixD rho = *corrMat;
                TVectorD sigmas = _sigma2_etas;
                sigmas.Sqrt();
                Int_t n = _sigma2_etas.GetNrows();
                for( Int_t ii = 0; ii < n; ++ii ){ 
                    for( Int_t jj = 0; jj < n; ++jj ){ 
                        _C_eta[ii][jj] = sigmas[ii] * rho[ii][jj] * sigmas[jj]; 
                    }
                }
            }
        }

        void Init4Vects(){
            // Initial 4-Vectors
            for( auto P : _P_inits ) _P_init += P;
            for( auto P : _Ps_y ) _masses_y.push_back( P.M() ); 
            for( auto P : _Ps_x ) _masses_x.push_back( P.M() ); 

            _x.ResizeTo( 3 * _Ps_x.size() );
            _y.ResizeTo( 3 * _Ps_y.size() );
            _eta.ResizeTo(_y);

            _x   = constructTVecFrom4Vecs( _Ps_x );
            _y   = constructTVecFrom4Vecs( _Ps_y );
            _eta = _y;
            _epsilon.ResizeTo(_eta);
        }


        // ---------------------------------------------------------------------------------------------------------------------------------
        // Methods 
        // ---------------------------------------------------------------------------------------------------------------------------------

        // Main, Default Fitter
        void DoFitting(Int_t nter){
            // Default Fitting if you want to use default constructed A's, B's, and c's with default loop conditions

            TMatrixD *AA = nullptr; // To be address of A
            Int_t iter = 0;
            while ( iter < nter ) { // ... or some other conditions; I ain't your boss.
                SetB();
                SetC();

                if( _Ps_x.size() > 0 ){
                    SetA();
                    AA = &_A;
                }

                ProcessFit( &_B, AA, &_ccccc);
                _Ps_x  = get4Vectors( &_x, _masses_x );
                _Ps_y  = get4Vectors( &_y, _masses_y );

                _epsilon += _delta;
                Double_t chi20 = _chi2;
                SetChi2();

                //if( abs(_chi2 - chi20)/chi20 < 0.01 ){ 
                if( _chi2 - chi20 > 0. ){ 
                    if( iter > 0 ){
                        break; 
                    }
                }
                ++iter;

            }

            PostProcess();

        }

        // ---------------------------------------------------------------------------------------------------------------------------------
        // Setters : Construct needed Vector and Matrices
        // ---------------------------------------------------------------------------------------------------------------------------------

        void SetDelta(TMatrixD *AA = nullptr ){
            // Construct delta ( = -_C_eta B^T C_B(  c + A xi ) = -_C_eta B^T * mu )
            _delta.ResizeTo(_y);

            TMatrixD BT = _B;
            BT.T();

            TMatrixD D  = _C_eta * BT;

            SetMu(AA);

            _delta =  D * _mu;
            _delta *= -1;
        }

        void SetMu( TMatrixD *AA ){
            // Construct mu ( = C_B(  c + A xi )  )
            _mu.ResizeTo(_ccccc);
            _mu = _ccccc;

            SetC_B();

            if( AA != 0 ){
                SetXi();
                _mu += _A * _xi;
            }

            _mu *= _C_B;
        }

        void SetC_B(){
            // Maybe implement a similarity transformation : B C B^T and B^T C B
            // Construct Matrix C_B ( = inverse(B _C_eta B^T) )
            TMatrixD BT = _B;
            BT.T();

            _C_B.ResizeTo( _B.GetNrows(), _B.GetNrows() );
            _C_B = _B * _C_eta * BT;
            _C_B.SetTol(1.e-30); // Set tolerance to be much lower for inverted matrix
            _C_B.Invert();
        }

        void SetXi(){
            // Construct xi    ( = -inverse(A^T C_B A) A^T C_B c = - C_x A^T C_B c)
            TMatrixD AT = _A;
            AT.T();

            SetC_x();

            TMatrixD C = _C_x * AT * _C_B;      // C_B is defined in getMu

            // BELOW NOT GOOD BECAUSE *= only works when shape is maintained
            //C *= _C_B;
            //C *= AA;
            //C.Invert();
            //C *= AT * _C_B;
            //C *= _C_B;

            _xi.ResizeTo(_x);
            _xi = C*_ccccc;
            _xi *= -1;
        }

        void SetC_x(){
            // Returns the new covariance matrix from propagation of errors
            // Construct C_x    ( = inverse(A^T C_B A) ) 
            TMatrixD AT = _A;
            AT.T();

            _C_x.ResizeTo(_x.GetNrows(), _x.GetNrows());
            _C_x = AT * _C_B * _A;
            _C_x.SetTol(1.e-30); // Set tolerance to be much lower for inverted matrix

            _C_x.Invert();
        }

        void SetC_y( ){
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

            if( _A.GetNrows() != 0 ){
                TMatrixD AT = _A;
                _A.T();
                TMatrixD C_y3 = C_y1 * _A * _C_x * AT * _C_B * _B * _C_eta;    
                _C_y += C_y3;
            }
        }

        void SetB( TMatrixD* BB ){
            _B.ResizeTo( *BB );
            _B = *BB;
        }
        void SetA( TMatrixD* AA ){
            _A.ResizeTo( *AA );
            _A = *AA;
        }
        void SetC( TVectorD* cc ){
            _ccccc.ResizeTo( *cc);
            _ccccc = *cc;
        }

        void SetB(){
            _B.ResizeTo(           _nconstraints, _nvars_y      );
            _B = constructDerMatrix( &_y, _masses_y );
        }

        void SetA(){
            _A.ResizeTo(           _nconstraints, _nvars_x      ); 
            _A = constructDerMatrix( &_x, _masses_x ); 
        }

        void SetC(){
            _ccccc.ResizeTo( _nconstraints );
            // Construct constraint functions f(x,y) == 0
            // c[i] = constraint_i               @ x_0, y_0

            // Conservation of momentum constratints
            SetP_fin();
            _Pdiff = _P_init - _P_fin;  
            for( Int_t ii = 0; ii < std::min(_nconstraints, 4); ++ii){
                _ccccc[ii] = _Pdiff[ii];
            }

            //// Particulars
            //// Extra constraints (invariant mass of 2 photons)
            //if( _nconstraints > 4 ){
            //  Int_t nparticles = _Ps_y.size();
            //  TLorentzVector P_phot1 = _Ps_y[nparticles - 2];
            //  TLorentzVector P_phot2 = _Ps_y[nparticles - 1];
            //  
            //  //Int_t nparticles = _Ps_x.size();
            //  //TLorentzVector P_phot1 = _Ps_x[nparticles - 2];
            //  //TLorentzVector P_phot2 = _Ps_x[nparticles - 1];
            //  
            //  TLorentzVector P_pi0 = P_phot1 + P_phot2;
            //  //TLorentzVector P_pi0 = _Ps_x[0];
            //  
            //  c[4] = (_Pdiff+P_pi0).M2() - (P_pi0).M2() ;
            //  //c[4] = ( P_phot1 + P_phot2 ).M2() - pi0Mass*pi0Mass;
            //  //c[4] = ( _Pdiff + P_phot1 + P_phot2 ).M2() + ( P_phot1 + P_phot2 ).M2() - 2 * pi0Mass*pi0Mass;
            //}
        }

        void SetP_fin(){
            _P_fin = TLorentzVector(0, 0, 0, 0);
            for ( auto P : _Ps_y ) _P_fin += P; 
            for ( auto P : _Ps_x ) _P_fin += P; 
        }

        void ProcessFit( TMatrixD *BB, TMatrixD *AA, TVectorD* cc ){
            _ccccc.ResizeTo(*cc);
            _B.ResizeTo(*BB);

            _ccccc = *cc;              // Just in case these come from outside the class
            _B = *BB;              // Just in case these come from outside the class

            SetDelta( AA );

            _y += _delta; 

            if( AA != 0 ) _x += _xi;     // _xi is defined in "Set Delta" if A is used
        }

        void PostProcess(){
            // Get change of new fit results (y):  epsilon ( eta - y ) !! Probably won't need
            //_epsilon -= y;

            SetChi2();

            SetConfLevel();

            // Correlations Matrix
            SetC_y( );
            _sigma2_ys = TMatrixDDiag(_C_y);

            SetPulls();
        }


        // ---------------------------------------------------------------------------------------------------------------------------------
        // Postfit 
        // ---------------------------------------------------------------------------------------------------------------------------------

        void SetChi2(){
            // For a kinematic fit, this should follow a chi^2 distribution
            // Construct chi^2 ( = epsilon^T inverse(_C_eta) epsilon = epsilon . epsilon')

            TVectorD eps0  = _epsilon;

            TMatrixD C_nI = _C_eta;
            C_nI.SetTol(1.e-30); // Set tolerance to be much lower for inverted matrix
            C_nI.Invert(); 

            TVectorD eps1 = C_nI * eps0; 

            _chi2 = eps0 * eps1;
            //
            ////TVectorD f = _B * _delta + c; 
            ////if( _nvars_x > 0 ){
            ////  f += _A* _xi;
            ////}
            ////_chi2 = _mu *  f ;

            //// Easier evaluation with no additional matrix inversions
            //TMatrixD BT = _B;
            //BT.T();
            //_chi2 = _delta * ( BT * _mu );
            //_chi2 *= -1;
        }

        void SetConfLevel(){
            // Returns confidence level for a given chi^2 (chi2) and number of degrees of freedom (ndf) 
            _confLevel = TMath::Prob( _chi2, _ndf );
        }

        void SetPulls( ){
            // Returns a vector of pull distributions
            // Get pulls z_i ( = [ eta_i - y_i ] / sqrt( sigma^2_eta_i - sigma^2_y_i )
            _pulls.ResizeTo(_eta);

            TVectorD num = _epsilon;

            TVectorD denom = _sigma2_etas - _sigma2_ys;
            //TVectorD denom = _sigma2_etas;

            denom.Abs();                          // !!! Probably don't need! In case sigmas2_y < sigmas2_n
            denom.Sqrt();                         // Element-wise square root!

            _pulls = ElementDiv( num, denom );    // Element-wise division!
        }


        // ---------------------------------------------------------------------------------------------------------------------------------
        // These methods are on the bottom because they are very specfic. The rest should be pretty general
        // ---------------------------------------------------------------------------------------------------------------------------------

        std::vector<TLorentzVector> get4Vectors( TVectorD *v, std::vector<Double_t> masses ){
            std::vector<TLorentzVector> result;
            Int_t nparticles = masses.size();
            for(Int_t ii = 0; ii < nparticles; ++ii){
                Double_t m     = masses[ii];
                Double_t p     = (*v)[ 3* ii];
                Double_t theta = (*v)[ 3* ii + 1];
                Double_t phi   = (*v)[ 3* ii + 2];
                //Double_t E     = v[ 3*ii + 3];
                TVector3 vec3;
                vec3.SetMagThetaPhi(1, theta, phi) ;
                vec3 *= p;
                result.push_back( TLorentzVector( vec3, sqrt( p*p + m*m ) ) );
            }
            return result;
        }


        TVectorD constructTVecFrom4Vecs( std::vector<TLorentzVector> P4s ){
            Int_t nparticles = P4s.size();
            Int_t nvars      = 3 * nparticles;
            TVectorD vec( nvars );
            for( Int_t ii = 0; ii < nparticles; ++ii){
                TLorentzVector P = P4s[ii];
                vec[3*ii+0] = P.P();
                vec[3*ii+1] = P.Theta();
                vec[3*ii+2] = P.Phi();
            }
            return vec;
        }

        TMatrixD constructDerMatrix( TVectorD *vec, std::vector<Double_t> masses = {} ){
            // Constructs the derivative matrices. At the moment, it is identical for A and B

            std::vector<TLorentzVector> fins_vec = get4Vectors( vec, masses );

            Int_t nparticles   = fins_vec.size();
            Int_t nvars        = 3 * nparticles;

            TMatrixD mat( _nconstraints, nvars );
            TLorentzVector Pcurrent;

            // Get Derivative matrix for conservation of 4-momenta 
            for( Int_t jj = 0; jj < nparticles; ++jj){
                Pcurrent = fins_vec[jj];
                TMatrixD dfdx =  getDfDx( &Pcurrent );
                mat.SetSub( 0, 3*jj, dfdx);
            }

            //// If there are more constraints, go back and edit those matrices
            //if( nparticles >= 2 && _nconstraints > 4 ){
            //  
            //  if( masses[nparticles-1] + masses[nparticles -2] < 0.001 && _nconstraints_extra > 0){
            //    // Very specific but particles involved in the extra constraints should come at the end of the list
            //    Int_t iphot1 = nparticles - 2;
            //    Int_t iphot2 = nparticles - 1;
            //    TLorentzVector P_phot1 = fins_vec[iphot1];
            //    TLorentzVector P_phot2 = fins_vec[iphot2];

            //    TMatrixD dfdx =  getDfDx( &P_phot1, &P_phot2 );
            //    mat.SetSub( 0, 3*iphot1, dfdx);

            //    dfdx = getDfDx( &P_phot2, &P_phot1 );
            //    mat.SetSub( 0, 3*iphot2, dfdx);
            //  }
            //  if( nparticles > 3 ){
            //    // Very specific but particles involved in the extra constraints should come at the end of the list
            //    TLorentzVector P1 = fins_vec[0];
            //    TLorentzVector P2 = fins_vec[1];

            //    TMatrixD dfdx =  getDfDx( &P1, &P2 );
            //    mat.SetSub( 0, 0, dfdx);

            //    dfdx = getDfDx( &P2, &P1 );
            //    mat.SetSub( 0, 3, dfdx);
            //  }

            //}

            return mat;
        }

        TMatrixD getDfDx( TLorentzVector *P1, TLorentzVector *P2 = nullptr ){
            // Constraints:
            //   -Conservation of Momentum:
            //       Px = P * sin(theta) * cos(phi);      =>    dPx/dP = sin(theta) * cos(phi);    dPx/dTheta =  P * cos(theta) * cos(phi);   dPx/dphi = -P * sin(theta) * sin(phi);
            //       Py = P * sin(theta) * sin(phi);      =>    dPy/dP = sin(theta) * sin(phi);    dPy/dTheta =  P * cos(theta) * sin(phi);   dPy/dphi =  P * sin(theta) * cos(phi);
            //       Pz = P * cos(theta);                 =>    dPz/dP = cos(theta);               dPz/dTheta = -P * sin(theta);              dPz/dphi =  0; 
            //       E  = sqrt(P^2 + M^2);                =>    dE/dP  = P/sqrt( P^2 + M^2 );      dE/dTheta  =  0;                           dE/dPhi  =  0;
            //   -Invariant pi0 Mass:
            //       Mpi02 = 2 P1 P2 (1 -  [sin(theta1) * cos(phi1) * sin(theta2)*cos(phi2) + sin(theta1) * sin(phi1) * sin(theta2) * sin(phi2) + cos(theta1) * cos(theta2) ] )
            //          =>    Mpi02/dP1     = 2 * P2(1 -  [sin(theta1) * cos(phi1) * sin(theta2)*cos(phi2) + sin(theta1) * sin(phi1) * sin(theta2) * sin(phi2) + cos(theta1) * cos(theta2) ] )
            //          =>    Mpi02/dtheta1 = 2*P1*P2(  cos(theta1) * cos(phi1) * sin(theta2)*cos(phi2) + cos(theta1) * sin(phi1) * sin(theta2) * sin(phi2) - sin(theta1) * cos(theta2) )
            //          =>    Mpi02/dphi1   = 2*P1*P2( -sin(theta1) * sin(phi1) * sin(theta2)*cos(phi2) + sin(theta1) * cos(phi1) * sin(theta2) * sin(phi2)  )

            Double_t theta1 = P1->Theta();
            Double_t phi1   = P1->Phi();
            Double_t p1     = P1->P();
            Double_t E1     = P1->E();

            Int_t nvars = 3; 

            //Double_t zero = 3E-8;
            ////Double_t data[4][3] = { { sin(theta1)*cos(phi1),  p1*cos(theta1)*cos(phi1), -p1*sin(theta1)*sin(phi1) } ,
            ////                        { sin(theta1)*sin(phi1),  p1*cos(theta1)*sin(phi1),  p1*sin(theta1)*cos(phi1) } ,
            ////                        { cos(theta1)          , -p1*sin(theta1)          ,  zero                     } ,
            ////                        { p1/E1                ,  zero                    ,  zero                     } };
            //Double_t data[4][3] = { { sin(theta1)*cos(phi1),  p1*cos(theta1)*cos(phi1), -p1*sin(theta1)*sin(phi1) } ,
            //                        { sin(theta1)*sin(phi1),  p1*cos(theta1)*sin(phi1),  p1*sin(theta1)*cos(phi1) } ,
            //                        { cos(theta1)          , -p1*sin(theta1)          ,  zero                     } ,
            //                        { zero                ,  zero                    ,  zero                     } };

            // this must be for a 4-C fit, all final state particles measured: (?)
            Double_t data[4][3] = { { sin(theta1)*cos(phi1),  p1*cos(theta1)*cos(phi1), -p1*sin(theta1)*sin(phi1) } ,
                { sin(theta1)*sin(phi1),  p1*cos(theta1)*sin(phi1),  p1*sin(theta1)*cos(phi1) } ,
                { cos(theta1)          , -p1*sin(theta1)          ,  0                        } ,
                { p1/E1                ,  0                       ,  0                        } };

            //Double_t data[5][3] = { { sin(theta1)*cos(phi1),  p1*cos(theta1)*cos(phi1), -p1*sin(theta1)*sin(phi1) } ,
            //                        { sin(theta1)*sin(phi1),  p1*cos(theta1)*sin(phi1),  p1*sin(theta1)*cos(phi1) } ,
            //                        { cos(theta1)          , -p1*sin(theta1)          ,  0                        } ,
            //                        { p1/E1                ,  0                       ,  0                        } ,
            //                        { 0                    ,  0                       ,  0                        }  };
            //Double_t data[3][3] = { { sin(theta1)*cos(phi1),  p1*cos(theta1)*cos(phi1), -p1*sin(theta1)*sin(phi1) } ,
            //                        { sin(theta1)*sin(phi1),  p1*cos(theta1)*sin(phi1),  p1*sin(theta1)*cos(phi1) } ,
            //                        { cos(theta1)          , -p1*sin(theta1)          ,  0                        }  };

            //if( P2 != 0 ){
            //  // Fill in bottom row
            //  Double_t theta2 = P2->Theta();
            //  Double_t phi2   = P2->Phi();
            //  Double_t p2     = P2->P();
            //  Double_t E2     = P2->E();
            //  Int_t sgn = -1;  
            //  if( P2->M() + P1->M() < 1.005 ){
            //    sgn = 1;
            //  }
            //  data[4][0] = sgn * 2 * ( p1/E1*E2 - p2*cos( P1->Angle( P2->Vect() ) )); 
            //  data[4][1] = sgn * 2 * p1 * p2 * (   - ( cos(theta1)*cos(phi1) * sin(theta2)*cos(phi2) + cos(theta1)*sin(phi1) * sin(theta2)*sin(phi2) - sin(theta1)*cos(theta2)) ); 
            //  data[4][2] = sgn * 2 * p1 * p2 * (   - (-sin(theta1)*sin(phi1) * sin(theta2)*cos(phi2) + sin(theta1)*cos(phi1) * sin(theta2)*sin(phi2)                          ) ); 
            //}

            TMatrixD dfdx( _nconstraints, nvars, *data );

            dfdx *= -1.;
            return dfdx;
        }

        ClassDef(KinFitter,1)
    };

#endif

