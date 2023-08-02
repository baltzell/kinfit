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
        /*C_P_hist = *(TH3D*)in_File->Get("C_P_hist");
        C_theta_hist = *(TH3D*)in_File->Get("C_theta_hist");
        C_phi_hist = *(TH3D*)in_File->Get("C_phi_hist");
        C_P_theta_hist = *(TH3D*)in_File->Get("C_P_theta_hist");
        C_P_phi_hist = *(TH3D*)in_File->Get("C_P_phi_hist");
        C_theta_phi_hist = *(TH3D*)in_File->Get("C_theta_phi_hist");*/
        C_P_hist = *(TH3D *)in_File->Get("C_P_cut");
        C_theta_hist = *(TH3D *)in_File->Get("C_theta_cut");
        C_phi_hist = *(TH3D *)in_File->Get("C_phi_cut");
        C_P_theta_hist = *(TH3D *)in_File->Get("C_P_theta_cut");
        C_P_phi_hist = *(TH3D *)in_File->Get("C_P_phi_cut");
        C_theta_phi_hist = *(TH3D *)in_File->Get("C_theta_phi_cut");
    }

    // Interpolate the covariance matrix to the momenta (in GeV) and angles (in degrees) provided as arguments
    TMatrixD Interpolate(int sector, double p, double theta, double phi)
    {
        TMatrixD Cov_Matrix(3, 3);

        double P_in_vector = p;
        double Theta_in_vector = theta;

        // So far studies are done for sector 2, requiering to do some arythmetics on phi using the sector number
        double Temp_Phi_in_vector = (phi < 0. && sector > 1) ? phi + 360 : phi;
        double Phi_in_vector = (sector == 1) ? Temp_Phi_in_vector + 60 : Temp_Phi_in_vector - (sector - 2) * 60.; // fmod(Temp_Phi_in_vector,60.)+60.;

        // std::cout<<"sector "<<sector<<" in phi "<<phi<<" shifted phi "<<Phi_in_vector<<std::endl;

        // Diagonal terms
        // Momentum resolution
        Cov_Matrix[0][0] = C_P_hist.Interpolate(P_in_vector, Theta_in_vector, Phi_in_vector);
        // Theta resolution
        Cov_Matrix[1][1] = C_theta_hist.Interpolate(P_in_vector, Theta_in_vector, Phi_in_vector);
        // Phi resolution
        Cov_Matrix[2][2] = C_phi_hist.Interpolate(P_in_vector, Theta_in_vector, Phi_in_vector);

        // Off Diagonal terms
        // Momentum/Theta covariance
        Cov_Matrix[0][1] = C_P_theta_hist.Interpolate(P_in_vector, Theta_in_vector, Phi_in_vector);
        Cov_Matrix[1][0] = Cov_Matrix[0][1];

        // Momentum/Phi covariance
        Cov_Matrix[0][2] = C_P_phi_hist.Interpolate(P_in_vector, Theta_in_vector, Phi_in_vector);
        Cov_Matrix[2][0] = Cov_Matrix[0][2];

        // Phi/Theta covariance
        Cov_Matrix[1][2] = C_theta_phi_hist.Interpolate(P_in_vector, Theta_in_vector, Phi_in_vector);
        Cov_Matrix[2][1] = Cov_Matrix[1][2];

        // cout << P_in_vector << " " << Theta_in_vector << " " << Phi_in_vector << endl;
        //  Cov_Matrix.Print();
        return Cov_Matrix;
    }

    TMatrixD Interpolate(int sector, TLorentzVector in_vector)
    {
        return Interpolate(sector, in_vector.P(), in_vector.Theta() * TMath::RadToDeg(), in_vector.Phi() * TMath::RadToDeg());
    }

    bool Is_good_to_interpolate(int sector, TLorentzVector in_vector)
    {
        // very bad implementation
        double phi = in_vector.Phi() * TMath::RadToDeg();
        double Temp_Phi_in_vector = (phi < 0. && sector > 1) ? phi + 360 : phi;
        double Phi_in_vector = (sector == 1) ? Temp_Phi_in_vector + 60 : Temp_Phi_in_vector - (sector - 2) * 60.; // fmod(Temp_Phi_in_vector,60.)+60.;
        return Is_good_to_interpolate(in_vector.P(), in_vector.Theta() * TMath::RadToDeg(), Phi_in_vector);
    }

    bool Is_good_to_interpolate(double p, double theta, double phi)
    {
        // each component of the covariance matrix must be in region where the extrapolation is possible
        return (Is_good_to_interpolate_hist(C_P_hist, p, theta, phi) &&
                Is_good_to_interpolate_hist(C_theta_hist, p, theta, phi) &&
                Is_good_to_interpolate_hist(C_phi_hist, p, theta, phi) &&
                Is_good_to_interpolate_hist(C_P_theta_hist, p, theta, phi) &&
                Is_good_to_interpolate_hist(C_P_phi_hist, p, theta, phi) &&
                Is_good_to_interpolate_hist(C_theta_phi_hist, p, theta, phi));
    }

    bool Is_good_to_interpolate_hist(TH3D hist, double p, double theta, double phi)
    {
        // Ask that interpolation can be done only if each neighboor point have none zero values
        Int_t ubx = hist.GetXaxis()->FindBin(p);
        if (p < hist.GetXaxis()->GetBinCenter(ubx))
            ubx -= 1;
        Int_t obx = ubx + 1;

        Int_t uby = hist.GetYaxis()->FindBin(theta);
        if (theta < hist.GetYaxis()->GetBinCenter(uby))
            uby -= 1;
        Int_t oby = uby + 1;

        Int_t ubz = hist.GetZaxis()->FindBin(phi);
        if (phi < hist.GetZaxis()->GetBinCenter(ubz))
            ubz -= 1;
        Int_t obz = ubz + 1;

        return ((hist.GetBinContent(ubx, uby, ubz) != 0.0) &&
                (hist.GetBinContent(ubx, uby, obz) != 0.0) &&
                (hist.GetBinContent(ubx, oby, ubz) != 0.0) &&
                (hist.GetBinContent(ubx, oby, obz) != 0.0) &&
                (hist.GetBinContent(obx, uby, ubz) != 0.0) &&
                (hist.GetBinContent(obx, uby, obz) != 0.0) &&
                (hist.GetBinContent(obx, oby, ubz) != 0.0) &&
                (hist.GetBinContent(obx, oby, obz) != 0.0));
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
