#ifndef KinCovariance_h
#define KinCovariance_h

#include "TH3D.h"
#include "TFile.h"
#include <iostream>

double Rotate_to_sector2(double phi)
{
    double rotated_phi = (phi < -30.) ? phi + 360 : phi;
    return rotated_phi = (phi > 30.) ? rotated_phi - (60. * ((int)(rotated_phi / 60.) - 1)) : rotated_phi + 60.;
}

class KinCovariance
{

public:
    virtual ~KinCovariance() {}
    KinCovariance() {}

    KinCovariance(TString in_File_path)
    {
        TFile *in_File = new TFile(in_File_path);
        _C_P_hist = *(TH3D *)in_File->Get("C_P_hist");
        _C_theta_hist = *(TH3D *)in_File->Get("C_theta_hist");
        _C_phi_hist = *(TH3D *)in_File->Get("C_phi_hist");
        _C_P_theta_hist = *(TH3D *)in_File->Get("C_P_theta_hist");
        _C_P_phi_hist = *(TH3D *)in_File->Get("C_P_phi_hist");
        _C_theta_phi_hist = *(TH3D *)in_File->Get("C_theta_phi_hist");

        // The covariance matrix histogramms have various name, we must settle on one and stick to it

        /*_C_P_hist = *(TH3D *)in_File->Get("C_P_cut");
        _C_theta_hist = *(TH3D *)in_File->Get("C_theta_cut");
        _C_phi_hist = *(TH3D *)in_File->Get("C_phi_cut");
        _C_P_theta_hist = *(TH3D *)in_File->Get("C_P_theta_cut");
        _C_P_phi_hist = *(TH3D *)in_File->Get("C_P_phi_cut");
        _C_theta_phi_hist = *(TH3D *)in_File->Get("C_theta_phi_cut");*/

        _C_P_err_hist = *(TH3D *)in_File->Get("C_P_err_hist");
        _C_theta_err_hist = *(TH3D *)in_File->Get("C_theta_err_hist");
        _C_phi_err_hist = *(TH3D *)in_File->Get("C_phi_err_hist");
        _C_P_theta_err_hist = *(TH3D *)in_File->Get("C_P_theta_err_hist");
        _C_P_phi_err_hist = *(TH3D *)in_File->Get("C_P_phi_err_hist");
        _C_theta_phi_err_hist = *(TH3D *)in_File->Get("C_theta_phi_err_hist");

        _entries_hist = *(TH3D *)in_File->Get("entries_hist");
    }

    // Interpolate the covariance matrix to the momenta (in GeV) and angles (in degrees) provided as arguments
    TMatrixD Interpolate(int sector, double p, double theta, double phi)
    {
        TMatrixD Cov_Matrix(3, 3);

        double P_in_vector = p;
        double Theta_in_vector = theta;

        // So far studies are done for sector 2, requiering to do some arythmetics on phi using the sector number
        // double Temp_Phi_in_vector = (phi < 0. && sector > 1) ? phi + 360 : phi;
        double Phi_in_vector = Rotate_to_sector2(phi); //(sector == 1) ? Temp_Phi_in_vector + 60 : Temp_Phi_in_vector - (sector - 2) * 60.;

        std::cout << "In cov mat" << std::endl;
        std::cout << "P " << P_in_vector << " Theta " << Theta_in_vector << " Phi " << Phi_in_vector << std::endl;
        // Diagonal terms
        // Momentum resolution
        Cov_Matrix[0][0] = _C_P_hist.Interpolate(P_in_vector, Theta_in_vector, Phi_in_vector);
        // Theta resolution
        Cov_Matrix[1][1] = _C_theta_hist.Interpolate(P_in_vector, Theta_in_vector, Phi_in_vector) * (TMath::DegToRad() * TMath::DegToRad());
        // Phi resolution
        Cov_Matrix[2][2] = _C_phi_hist.Interpolate(P_in_vector, Theta_in_vector, Phi_in_vector) * (TMath::DegToRad() * TMath::DegToRad());

        // Off Diagonal terms
        // Momentum/Theta covariance
        Cov_Matrix[0][1] = _C_P_theta_hist.Interpolate(P_in_vector, Theta_in_vector, Phi_in_vector) * (TMath::DegToRad());
        Cov_Matrix[1][0] = Cov_Matrix[0][1];

        // Momentum/Phi covariance
        Cov_Matrix[0][2] = _C_P_phi_hist.Interpolate(P_in_vector, Theta_in_vector, Phi_in_vector) * (TMath::DegToRad());
        Cov_Matrix[2][0] = Cov_Matrix[0][2];

        // Phi/Theta covariance
        Cov_Matrix[1][2] = _C_theta_phi_hist.Interpolate(P_in_vector, Theta_in_vector, Phi_in_vector) * (TMath::DegToRad() * TMath::DegToRad());
        Cov_Matrix[2][1] = Cov_Matrix[1][2];

        SetEntries(_entries_hist.GetBinContent(_entries_hist.GetBin(_entries_hist.GetXaxis()->FindBin(P_in_vector), _entries_hist.GetYaxis()->FindBin(Theta_in_vector), _entries_hist.GetZaxis()->FindBin(Phi_in_vector))));

        return Cov_Matrix;
    }

    TMatrixD Interpolate(int sector, TLorentzVector in_vector)
    {
        return Interpolate(sector, in_vector.P(), in_vector.Theta() * TMath::RadToDeg(), in_vector.Phi() * TMath::RadToDeg());
    }

    // Methods to get the error of the covariance matrix
    TMatrixD Interpolate_Error(int sector, double p, double theta, double phi)
    {
        // std::cout<<"In err. cov mat"<<std::endl;
        TMatrixD Cov_Err_Matrix(3, 3);
        // return Cov_Err_Matrix;
        double P_in_vector = p;
        double Theta_in_vector = theta;

        // So far studies are done for sector 2, requiering to do some arythmetics on phi using the sector number
        // double Temp_Phi_in_vector = (phi < 0. && sector > 1) ? phi + 360 : phi;
        double Phi_in_vector = Rotate_to_sector2(phi); //(sector == 1) ? Temp_Phi_in_vector + 60 : Temp_Phi_in_vector - (sector - 2) * 60.;

        // std::cout<<"In err. cov mat"<<std::endl;
        // std::cout<<"P "<<P_in_vector<<" Theta "<<Theta_in_vector<<" Phi "<<Phi_in_vector<<std::endl;

        // Diagonal terms
        // Momentum resolution
        Cov_Err_Matrix[0][0] = _C_P_err_hist.Interpolate(P_in_vector, Theta_in_vector, Phi_in_vector);
        // Theta resolution
        Cov_Err_Matrix[1][1] = _C_theta_err_hist.Interpolate(P_in_vector, Theta_in_vector, Phi_in_vector);
        // Phi resolution
        Cov_Err_Matrix[2][2] = _C_phi_err_hist.Interpolate(P_in_vector, Theta_in_vector, Phi_in_vector);

        // Off Diagonal terms
        // Momentum/Theta covariance
        Cov_Err_Matrix[0][1] = _C_P_theta_err_hist.Interpolate(P_in_vector, Theta_in_vector, Phi_in_vector);
        Cov_Err_Matrix[1][0] = Cov_Err_Matrix[0][1];

        // Momentum/Phi covariance
        Cov_Err_Matrix[0][2] = _C_P_phi_err_hist.Interpolate(P_in_vector, Theta_in_vector, Phi_in_vector);
        Cov_Err_Matrix[2][0] = Cov_Err_Matrix[0][2];

        // Phi/Theta covariance
        Cov_Err_Matrix[1][2] = _C_theta_phi_err_hist.Interpolate(P_in_vector, Theta_in_vector, Phi_in_vector);
        Cov_Err_Matrix[2][1] = Cov_Err_Matrix[1][2];

        return Cov_Err_Matrix;
    }

    TMatrixD Interpolate_Error(int sector, TLorentzVector in_vector)
    {
        return Interpolate_Error(sector, in_vector.P(), in_vector.Theta() * TMath::RadToDeg(), in_vector.Phi() * TMath::RadToDeg());
    }

    // Verify that the 8 bins around the point to interpolate are none zero, this is done for each 6 component of the cov. matrix
    bool Is_good_to_interpolate(int sector, TLorentzVector in_vector)
    {
        // very bad implementation
        double phi = in_vector.Phi() * TMath::RadToDeg();
        // double Temp_Phi_in_vector = (phi < 0. && sector > 1) ? phi + 360 : phi;
        // double Phi_in_vector = (sector == 1) ? Temp_Phi_in_vector + 60 : Temp_Phi_in_vector - (sector - 2) * 60.; // fmod(Temp_Phi_in_vector,60.)+60.;
        double Phi_in_vector = Rotate_to_sector2(phi);

        return (Is_good_to_interpolate_hist(_C_P_hist, in_vector.P(), in_vector.Theta() * TMath::RadToDeg(), Phi_in_vector) &&
                Is_good_to_interpolate_hist(_C_theta_hist, in_vector.P(), in_vector.Theta() * TMath::RadToDeg(), Phi_in_vector) &&
                Is_good_to_interpolate_hist(_C_phi_hist, in_vector.P(), in_vector.Theta() * TMath::RadToDeg(), Phi_in_vector) &&
                Is_good_to_interpolate_hist(_C_P_theta_hist, in_vector.P(), in_vector.Theta() * TMath::RadToDeg(), Phi_in_vector) &&
                Is_good_to_interpolate_hist(_C_P_phi_hist, in_vector.P(), in_vector.Theta() * TMath::RadToDeg(), Phi_in_vector) &&
                Is_good_to_interpolate_hist(_C_theta_phi_hist, in_vector.P(), in_vector.Theta() * TMath::RadToDeg(), Phi_in_vector));
    }

    // Verify that the 8 bins around the point to interpolate are none zero, this is done for one specific component of the cov. matrix
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

    int GetEntries()
    {
        return _n_entries;
    }

    void SetEntries(int in_n_entries)
    {
        _n_entries = in_n_entries;
    }

private:
    TH3D _C_P_hist;
    TH3D _C_theta_hist;
    TH3D _C_phi_hist;
    TH3D _C_P_theta_hist;
    TH3D _C_P_phi_hist;
    TH3D _C_theta_phi_hist;

    TH3D _C_P_err_hist;
    TH3D _C_theta_err_hist;
    TH3D _C_phi_err_hist;
    TH3D _C_P_theta_err_hist;
    TH3D _C_P_phi_err_hist;
    TH3D _C_theta_phi_err_hist;

    TH3D _entries_hist;
    int _n_entries;
};

#endif
