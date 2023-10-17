#ifndef KinMomentumCorrections_h
#define KinMomentumCorrections_h

#include "TH3D.h"
#include "TFile.h"
#include <iostream>

class KinMomentumCorrections
{

public:
    virtual ~KinMomentumCorrections() {}
    KinMomentumCorrections() {}

    KinMomentumCorrections(TString in_File_path)
    {
        TFile *in_File = new TFile(in_File_path);
        _P_offset_hist = *(TH3D *)in_File->Get("P_offset_hist");
        _theta_offset_hist = *(TH3D *)in_File->Get("theta_offset_hist");
        _phi_offset_hist = *(TH3D *)in_File->Get("phi_offset_hist");
    }

    // Correct the measured 4-vector by the offset to MC generated 4-vector
    void Correct_4_Vector(int sector, TLorentzVector *in_vector)
    {
        double P_in_vector = in_vector->P();
        double Theta_in_vector = in_vector->Theta() * TMath::RadToDeg();

        // So far studies are done for sector 2, requiering to do some arythmetics on phi using the sector number
        double phi = in_vector->Phi() * TMath::RadToDeg();
        //double Temp_Phi_in_vector = (phi < 0. && sector > 1) ? phi + 360 : phi;
        double Phi_in_vector = Rotate_to_sector2(phi);//(sector == 1) ? Temp_Phi_in_vector + 60 : Temp_Phi_in_vector - (sector - 2) * 60.;

        // 4-vector offset
        double offset_P = _P_offset_hist.Interpolate(P_in_vector, Theta_in_vector, Phi_in_vector);
        double offset_theta = _theta_offset_hist.Interpolate(P_in_vector, Theta_in_vector, Phi_in_vector);
        double offset_phi = _phi_offset_hist.Interpolate(P_in_vector, Theta_in_vector, Phi_in_vector);

		//std::cout<<"Before offsets"<<in_vector->P()<<" "<<in_vector->Theta()*180./3.141592<<" "<<in_vector->Phi()*180./3.141592<<std::endl;
		//std::cout<<"offsets "<<offset_P<<" "<<offset_theta<<" "<<offset_phi<<std::endl;

        in_vector->SetPhi( (phi - offset_phi)*3.141592/180. );
        in_vector->SetTheta( (Theta_in_vector - offset_theta)*3.141592/180. );
        in_vector->SetRho(P_in_vector - offset_P);

		//std::cout<<"after offsets "<<P_in_vector - offset_P<<" "<<Theta_in_vector -offset_theta<<" "<<phi-offset_phi<<std::endl;
		//std::cout<<"after offsets 2 "<<in_vector->P()<<" "<<in_vector->Theta()*180./3.141592<<" "<<in_vector->Phi()*180./3.141592<<std::endl;

        return;
    }

private:
    TH3D _P_offset_hist;
    TH3D _theta_offset_hist;
    TH3D _phi_offset_hist;
};

#endif
