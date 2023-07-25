#ifndef KinCovariance_h
#define KinCovariance_h

#include "TH3D.h"
#include "TFile.h"
#include <iostream>


class KinCovariance
{

public:
    virtual ~KinCovariance() {}
    KinCovariance() {}

    KinCovariance(TString in_File_path)
    {
        TFile *in_File = new TFile(in_File_path);
        C_P_hist = *(TH3D*)in_File->Get("C_P_hist");
        C_theta_hist = *(TH3D*)in_File->Get("C_theta_hist");
        C_phi_hist = *(TH3D*)in_File->Get("C_phi_hist");
        C_P_theta_hist = *(TH3D*)in_File->Get("C_P_theta_hist");
        C_P_phi_hist = *(TH3D*)in_File->Get("C_P_phi_hist");
        C_theta_phi_hist = *(TH3D*)in_File->Get("C_theta_phi_hist");
    }

    //Interpolate the covariance matrix to the momenta (in GeV) and angles (in degrees) provided as arguments
    TMatrixD Interpolate(int sector, double p, double theta, double phi)
    {
        TMatrixD Cov_Matrix(3, 3);

        double P_in_vector = p;
        double Theta_in_vector = theta;

        //So far studies are done for sector 2, requiering to do some arythmetics on phi using the sector number
        double Temp_Phi_in_vector = (phi < 0. && sector>1) ? phi+360 : phi ;    
        double Phi_in_vector = (sector==1) ? Temp_Phi_in_vector+60 : Temp_Phi_in_vector-(sector-2)*60.;//fmod(Temp_Phi_in_vector,60.)+60.;

        //std::cout<<"sector "<<sector<<" in phi "<<phi<<" shifted phi "<<Phi_in_vector<<std::endl;

        // Diagonal terms
        // Momentum resolution
        Cov_Matrix[0][0] = C_P_hist.Interpolate(P_in_vector,Theta_in_vector,Phi_in_vector);
        // Theta resolution
        Cov_Matrix[1][1] = C_theta_hist.Interpolate(P_in_vector,Theta_in_vector,Phi_in_vector);
        // Phi resolution
        Cov_Matrix[2][2] = C_phi_hist.Interpolate(P_in_vector,Theta_in_vector,Phi_in_vector);

        // Off Diagonal terms
        // Momentum/Theta covariance
        Cov_Matrix[0][1] = C_P_theta_hist.Interpolate(P_in_vector,Theta_in_vector,Phi_in_vector);
        Cov_Matrix[1][0] = Cov_Matrix[0][1];

        // Momentum/Phi covariance
        Cov_Matrix[0][2] = C_P_phi_hist.Interpolate(P_in_vector,Theta_in_vector,Phi_in_vector);
        Cov_Matrix[2][0] = Cov_Matrix[0][2];

        // Phi/Theta covariance
        Cov_Matrix[1][2] = C_theta_phi_hist.Interpolate(P_in_vector,Theta_in_vector,Phi_in_vector);
        Cov_Matrix[2][1] = Cov_Matrix[1][2];

        //Cov_Matrix.Print();
        return Cov_Matrix;
    }

    TMatrixD Interpolate(int sector, TLorentzVector in_vector)
    {
        return Interpolate(sector, in_vector.P(), in_vector.Theta()*TMath::RadToDeg(), in_vector.Phi()*TMath::RadToDeg());
    }

private:
    TH3D C_P_hist;
    TH3D C_theta_hist;
    TH3D C_phi_hist;
    TH3D C_P_theta_hist;
    TH3D C_P_phi_hist;
    TH3D C_theta_phi_hist;
};

#endif
