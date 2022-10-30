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
    Double_t confLevel;                     //Probably the only things that need to be public 
    TVectorD pulls;                         //Probably the only things that need to be public 
    Double_t chi2;                          //Probably the only things that need to be public 
    
    Int_t nvars_x;
    Int_t nvars_y;
    
    Int_t nconstraints;
    Int_t nconstraints_extra;
    Int_t ndf;

    TVectorD eta;
    TVectorD sigma2_etas;
    TMatrixD C_eta;
    
    TVectorD y;
    TVectorD sigma2_ys;
    TMatrixD C_y;
    std::vector<std::string> names_y;
    
    TVectorD x;
    TVectorD xi;
    std::vector<std::string> names_x;
    
    TVectorD delta;
    TVectorD epsilon;
    
    std::vector<TLorentzVector> P_inits;      
    std::vector<Double_t> masses_x;         // Probably won't need 
    std::vector<Double_t> masses_y;         // Probably won't need
    std::vector<TLorentzVector> Ps_y;       //Probably the only things that need to be public 
    std::vector<TLorentzVector> Ps_x;       //Probably the only things that need to be public 


    TMatrixD A;
    TMatrixD B;
    
    TVectorD c;
    TLorentzVector P_init;
    TLorentzVector P_fin;
    TLorentzVector Pdiff;

    TMatrixD C_B;
    TVectorD mu;
    TMatrixD C_x;

    TH1D* h_chi2Viter = nullptr; //->
    Int_t niterations = 100; 
    
    Double_t chi2s[100]; 

    

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
        nconstraints_extra( ncons_extra               ),
        
        //nconstraints(       3 + ncons_extra           ),
        nconstraints(       4 + ncons_extra           ),
        nvars_x(            3 * P_Fins_x.size()       ), 
        nvars_y(            3 * P_InitsFins[1].size() ), 
        
        Ps_x(               P_Fins_x                  ),
        Ps_y(               P_InitsFins[1]            ), 
        P_inits(            P_InitsFins[0]            ) 
    {
      Initialize(sigmas, C_n);
      Init4Vects();  // y and eta defined here

      // Loop until result is satisfactory
      DoFitting(100);
    }
    
    
    // Use this constructor if you have custom { x, y, eta, A, B, c, sigmas ... }
    KinFitter( Int_t ncons, std::vector<Double_t> y0s, std::vector<Double_t> x0s, std::vector<Double_t> sigmas, TMatrixD* C_n = nullptr )
      : 
        nconstraints( ncons      ),
        nvars_x(      x0s.size() ), 
        nvars_y(      y0s.size() ) 
    {
      Initialize(sigmas, C_n);

      for( Int_t ii = 0; ii < nvars_x; ii++ ){ x[ii]   = x0s[ii]; }
      for( Int_t ii = 0; ii < nvars_y; ii++ ){ eta[ii] = y0s[ii]; y[ii] = y0s[ii];}
      epsilon.ResizeTo(eta);
    }

    
    
    // ---------------------------------------------------------------------------------------------------------------------------------
    // Initializers  
    // ---------------------------------------------------------------------------------------------------------------------------------
    
    void Initialize( std::vector<Double_t> sigmas, TMatrixD *C_n = nullptr, TMatrixD *corrMat = nullptr ){
      ndf = abs(nconstraints - nvars_x); 

      x.ResizeTo(   nvars_x );
      y.ResizeTo(   nvars_y );
      eta.ResizeTo( nvars_y );
      
      sigma2_ys.ResizeTo(   nvars_y );
      C_eta.ResizeTo( nvars_y, nvars_y ); 
      
      // Set up Covariance matrix
      if( C_n == 0 ){
        // If a covariance matrix is not given, build the covariance matrix from sigmas
        InitErrors(sigmas, corrMat);
      }
      else{
        C_eta = *C_n;
        TMatrixDDiag diag(C_eta);
        
        // Get variance from diagonal of covariance matrix
        sigma2_etas.ResizeTo( nvars_y ); 
        sigma2_etas = diag;
      }
    }

    void InitErrors(std::vector<Double_t> sigmas, TMatrixD* corrMat){
      // Initialize Errors
      TVectorD sigma_etas( sigmas.size() );
      for( Int_t ii = 0; ii < sigmas.size(); ii++ ){ sigma_etas[ii] = sigmas[ii]; }
      
      sigma2_etas.ResizeTo( nvars_y ); 
      sigma2_etas = sigma_etas;
      sigma2_etas.Sqr();
      
      SetCovMat(corrMat);
    }

    void SetCovMat(TMatrixD* corrMat = nullptr){
      if( corrMat == 0 ){
        //// Variance is sigma^2
        TMatrixDDiag diag(C_eta); 
        diag = sigma2_etas;
      }
      else{
        // I don't think there's a way to do vec^T * Mat * vec so here's a loop:
        TMatrixD rho = *corrMat;
        TVectorD sigmas = sigma2_etas;
        sigmas.Sqrt();
        Int_t n = sigma2_etas.GetNrows();
        for( Int_t ii = 0; ii < n; ii++ ){ 
          for( Int_t jj = 0; jj < n; jj++ ){ 
            C_eta[ii][jj] = sigmas[ii] * rho[ii][jj] * sigmas[jj]; 
          }
        }
      }

    }
    
    void Init4Vects(){
      // Initial 4-Vectors

      for( auto P : P_inits ){ P_init += P; }

      for( auto P : Ps_y ){ 
        masses_y.push_back( P.M() ); 
      }
      for( auto P : Ps_x ){ 
        masses_x.push_back( P.M() ); 
      }
      
      x.ResizeTo( 3 * Ps_x.size() );
      y.ResizeTo( 3 * Ps_y.size() );
      eta.ResizeTo(y);

      x   = constructTVecFrom4Vecs( Ps_x );
      y   = constructTVecFrom4Vecs( Ps_y );
      eta = y;
      epsilon.ResizeTo(eta);
    }

    
    // ---------------------------------------------------------------------------------------------------------------------------------
    // Methods 
    // ---------------------------------------------------------------------------------------------------------------------------------

    // Main, Default Fitter
    void DoFitting(Int_t nter){
      // Default Fitting if you want to use default constructed A's, B's, and c's with default loop conditions
      
      TMatrixD *AA = nullptr; // To be address of A
      Int_t iter = 0;
      for( Int_t ii = 0; ii < 100; ii++){ chi2s[ii] = 0.; }
      while( iter < nter ){ // ... or some other conditions; I ain't your boss.
        SetB();
        SetC();
        
        if( Ps_x.size() > 0 ){
          SetA();
          AA = &A;
        }
        
        ProcessFit( &B, AA, &c);
        Ps_x  = get4Vectors( &x, masses_x );
        Ps_y  = get4Vectors( &y, masses_y );
        
        epsilon += delta;
        Double_t chi20 = chi2;
        SetChi2();

        //if( abs(chi2 - chi20)/chi20 < 0.01 ){ 
        chi2s[iter] = chi2;
        if( chi2 - chi20 > 0. ){ 
          if( iter > 0 ){
            niterations = iter; 
            break; 
          }
        }
        iter++;
        
      }

      PostProcess( );

    }
    
    // ---------------------------------------------------------------------------------------------------------------------------------
    // Setters : Construct needed Vector and Matrices
    // ---------------------------------------------------------------------------------------------------------------------------------
    
    void SetDelta(TMatrixD *AA = nullptr ){
      // Construct delta ( = -C_eta B^T C_B(  c + A xi ) = -C_eta B^T * mu )
      delta.ResizeTo(y);
      
      TMatrixD BT = B;
      BT.T();

      TMatrixD D  = C_eta * BT; // C_eta is defined at the ctor level
      
      SetMu(AA);
     
      delta =  D * mu;
      delta *= -1;
    }
    
    void SetMu( TMatrixD *AA ){
      // Construct mu ( = C_B(  c + A xi )  )
      mu.ResizeTo(c);
      mu = c;

      SetC_B();

      if( AA != 0 ){
        SetXi();
        mu += A * xi;
      }

      mu *= C_B;
    }

    void SetC_B(){
      // Maybe implement a similarity transformation : B C B^T and B^T C B
      // Construct Matrix C_B ( = inverse(B C_eta B^T) )
      TMatrixD BT = B;
      BT.T();

      C_B.ResizeTo( B.GetNrows(), B.GetNrows() );
      C_B = B * C_eta * BT;
      C_B.SetTol(1.e-30); // Set tolerance to be much lower for inverted matrix
      C_B.Invert();
    }

    void SetXi(){
      // Construct xi    ( = -inverse(A^T C_B A) A^T C_B c = - C_x A^T C_B c)
      TMatrixD AT = A;
      AT.T();

      SetC_x();

      TMatrixD C   = C_x * AT * C_B;      // C_B is defined in getMu
      //TMatrixD CC = A * C;
      //CC.Print();

      // BELOW NOT GOOD BECAUSE *= only works when shape is maintained
      //C *= C_B;
      //C *= AA;
      //C.Invert();
      //C *= AT * C_B;
      //C *= C_B;

      xi.ResizeTo(x);
      xi = C*c;
      xi *= -1;
    }

    void SetC_x(){
      // Returns the new covariance matrix from propagation of errors
      // Construct C_x    ( = inverse(A^T C_B A) ) 
      TMatrixD AT = A;
      AT.T();

      C_x.ResizeTo(x.GetNrows(), x.GetNrows());
      C_x = AT * C_B * A;
      C_x.SetTol(1.e-30); // Set tolerance to be much lower for inverted matrix

      C_x.Invert();
    }

    void SetC_y( ){
      // Returns the new covariance matrix from propagation of errors
      // Get Error in fit results: C_y ( = C_eta - C_eta B^T C_B B C_eta + C_eta B^T C_B A inverse( A^T C_B A ) A^T C_B B C_eta = C_eta - C_eta B^T C_B B C_eta + C_eta B^T C_B A C_x A^T C_B B C_eta)
      C_y.ResizeTo(C_eta);
      
      TMatrixD BT = B;
      BT.T();

      TMatrixD C_y1 = C_eta * BT * C_B;
      TMatrixD C_y2 = C_y1 * B * C_eta; 

      // Documentation recommended to add matrices this way
      C_y = C_eta;
      C_y -= C_y2;

      if( A.GetNrows() != 0 ){
        TMatrixD AT = A;
        AT.T();

        TMatrixD C_y3 = C_y1 * A * C_x * AT * C_B * B * C_eta;    

        C_y += C_y3;
        //C_y2.Print();
        //C_y3.Print();
      }
      //C_eta.Print();
      //C_y.Print();

    }
    
    void SetB( TMatrixD* BB ){
      B.ResizeTo( *BB );
      B = *BB;
    }
    void SetA( TMatrixD* AA ){
      A.ResizeTo( *AA );
      A = *AA;
    }
    void SetC( TVectorD* cc ){
      c.ResizeTo( *cc);
      c = *cc;
    }
    
    
    void SetB(){
      B.ResizeTo(           nconstraints, nvars_y      );
      B = constructDerMatrix( &y, masses_y );
    }
    
    void SetA(){
      A.ResizeTo(           nconstraints, nvars_x      ); 
      A = constructDerMatrix( &x, masses_x ); 
    }

    void SetC(){
      c.ResizeTo( nconstraints );
      // Construct constraint functions f(x,y) == 0
      // c[i] = constraint_i               @ x_0, y_0

      // Conservation of momentum constratints
      SetP_fin();
      Pdiff = P_init - P_fin;  
      for( Int_t ii = 0; ii < std::min(nconstraints, 4); ii++){
        c[ii] = Pdiff[ii];
      }

      //// Particulars
      //// Extra constraints (invariant mass of 2 photons)
      //if( nconstraints > 4 ){
      //  Int_t nparticles = Ps_y.size();
      //  TLorentzVector P_phot1 = Ps_y[nparticles - 2];
      //  TLorentzVector P_phot2 = Ps_y[nparticles - 1];
      //  
      //  //Int_t nparticles = Ps_x.size();
      //  //TLorentzVector P_phot1 = Ps_x[nparticles - 2];
      //  //TLorentzVector P_phot2 = Ps_x[nparticles - 1];
      //  
      //  TLorentzVector P_pi0 = P_phot1 + P_phot2;
      //  //TLorentzVector P_pi0 = Ps_x[0];
      //  
      //  c[4] = (Pdiff+P_pi0).M2() - (P_pi0).M2() ;
      //  //c[4] = ( P_phot1 + P_phot2 ).M2() - pi0Mass*pi0Mass;
      //  //c[4] = ( Pdiff + P_phot1 + P_phot2 ).M2() + ( P_phot1 + P_phot2 ).M2() - 2 * pi0Mass*pi0Mass;
      //}
    }

    void SetP_fin(){
      P_fin = TLorentzVector(0, 0, 0, 0);

      for( auto P : Ps_y ){ 
        P_fin += P; 
      }
      for( auto P : Ps_x ){ 
        P_fin += P; 
      }
    }
    
    void ProcessFit( TMatrixD *BB, TMatrixD *AA, TVectorD* cc ){
      c.ResizeTo(*cc);
      B.ResizeTo(*BB);
      
      c = *cc;              // Just in case these come from outside the class
      B = *BB;              // Just in case these come from outside the class
      
      SetDelta( AA );
      
      //delta *= 1./10;
      //delta.Print();
      
      y += delta; 

      if( AA != 0 ){
        x    += xi;     // xi is defined in "Set Delta" if A is used
      }
    }

    void PostProcess(){
      // Get change of new fit results (y):  epsilon ( eta - y ) !! Probably won't need
      //epsilon -= y;

      SetChi2();

      SetConfLevel();

      // Correlations Matrix
      SetC_y( );
      sigma2_ys = TMatrixDDiag(C_y);

      SetPulls();
    }


    // ---------------------------------------------------------------------------------------------------------------------------------
    // Postfit 
    // ---------------------------------------------------------------------------------------------------------------------------------

    void SetChi2(){
      // For a kinematic fit, this should follow a chi^2 distribution
      // Construct chi^2 ( = epsilon^T inverse(C_eta) epsilon = epsilon . epsilon')

      TVectorD eps0  = epsilon;
      
      TMatrixD C_nI = C_eta;
      C_nI.SetTol(1.e-30); // Set tolerance to be much lower for inverted matrix
      C_nI.Invert(); 
      
      TVectorD eps1 = C_nI * eps0; 
      
      chi2 = eps0 * eps1;
      //
      ////TVectorD f = B * delta + c; 
      ////if( nvars_x > 0 ){
      ////  f += A* xi;
      ////}
      //////cout << mu * f << endl;
      ////mu.Print();
      ////chi2 = mu *  f ;

      //// Easier evaluation with no additional matrix inversions
      //TMatrixD BT = B;
      //BT.T();
      //chi2 = delta * ( BT * mu );
      //chi2 *= -1;
    }

    void SetConfLevel(){
      // Returns confidence level for a given chi^2 (chi2) and number of degrees of freedom (ndf) 
      // Get confidence levels from chi^2 ( = TMath::Prob( chi^2, ndf ) )
      
      confLevel = TMath::Prob( chi2, ndf );
    }
    
    void SetPulls( ){
      // Returns a vector of pull distributions
      // Get pulls z_i ( = [ eta_i - y_i ] / sqrt( sigma^2_eta_i - sigma^2_y_i )
      pulls.ResizeTo(eta);
      
      TVectorD num = epsilon;
      
      TVectorD denom = sigma2_etas - sigma2_ys;
      //TVectorD denom = sigma2_etas;
      
      denom.Abs();                          // !!! Probably don't need! In case sigmas2_y < sigmas2_n
      denom.Sqrt();                         // Element-wise square root!

      pulls = ElementDiv( num, denom );    // Element-wise division!
    }


    // ---------------------------------------------------------------------------------------------------------------------------------
    // These methods are on the bottom because they are very specfic. The rest should be pretty general
    // ---------------------------------------------------------------------------------------------------------------------------------
    
    std::vector<TLorentzVector> get4Vectors( TVectorD *v, std::vector<Double_t> masses ){
      std::vector<TLorentzVector> result;
      Int_t nparticles = masses.size();
      for(Int_t ii = 0; ii < nparticles; ii++){
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
      Int_t nvars      = 3 * nparticles; // Each particle has 3 variables (p, theta, phi)
      TVectorD vec( nvars );
      for( Int_t ii = 0; ii < nparticles; ii++){
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
      Int_t nvars        = 3 * nparticles; // Each particle has 3 variables (p, theta, phi)
      
      TMatrixD mat( nconstraints, nvars );
      TLorentzVector Pcurrent;
      
      // Get Derivative matrix for conservation of 4-momenta 
      for( Int_t jj = 0; jj < nparticles; jj++){
        //cout << "~~~~~~~~~~~~~ WHATTTTTT !!!!!!~~~~~~~~~~~ " << endl;
        //mat.Print();
        //cout << "~~~~~~~~~~~~~ fhsdfhdsfkJKFHDFJHSDK ~~~~~~~~~~~ " << endl;
        Pcurrent = fins_vec[jj];
        TMatrixD dfdx =  getDfDx( &Pcurrent );
        //dfdx.UnitMatrix();
        mat.SetSub( 0, 3*jj, dfdx);
        //cout << "~~~~~~~~~~~~~ WHATTTTTT ~~~~~~~~~~~ " << endl;
        //mat.Print();
        //cout << "~~~~~~~~~~~~~ fhsdfhdsfkJKFHDFJHSDK ~~~~~~~~~~~ " << endl;
      }

      //// If there are more constraints, go back and edit those matrices
      //if( nparticles >= 2 && nconstraints > 4 ){
      //  
      //  if( masses[nparticles-1] + masses[nparticles -2] < 0.001 && nconstraints_extra > 0){
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

      TMatrixD dfdx( nconstraints, nvars, *data );
      //cout << "======startttt==================================" << endl;
      //dfdx.Print();
      //cout << "=============endddd==========================" << endl;


      dfdx *= -1.;
      return dfdx;
    }


    ClassDef(KinFitter,1)
};

#endif

