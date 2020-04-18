#include "MiniScatter_ROOTclasses.h"

#include "TH2D.h"
#include "TMath.h"
#include <cmath>

void MSR_plotRZgray_histoScalerVolume(TH2D* rzScaled) {
    double dA, dZ, dV;

    for (Int_t rBinIdx = 1; rBinIdx <= rzScaled->GetYaxis()->GetNbins(); rBinIdx++) {
        dA = TMath::Pi() * (pow(rzScaled->GetYaxis()->GetBinUpEdge(rBinIdx), 2) -
                            pow(rzScaled->GetYaxis()->GetBinLowEdge(rBinIdx),2) ); //[mm^2]
        for (Int_t zBinIdx = 1; zBinIdx <= rzScaled->GetXaxis()->GetNbins(); zBinIdx++) {
            dZ = rzScaled->GetXaxis()->GetBinUpEdge(zBinIdx) - 
                 rzScaled->GetXaxis()->GetBinLowEdge(zBinIdx); //[mm]
            dV = dA*dZ; //[mm^3]

            Int_t binIdx = rzScaled->GetBin(zBinIdx,rBinIdx);
            rzScaled->SetBinContent(binIdx, rzScaled->GetBinContent(binIdx)/dV); //[MeV/mm^3]
        }

    }
}

void MSR_plotZgray_sumR(TH2D* rzHist_Gy, TH1D* centerHist, Int_t yAxis_lastbin, Double_t r0_actual) {
        for (Int_t zIdx = 1; zIdx <= rzHist_Gy->GetNbinsX(); zIdx++ ) {
            double avgDoseZ = 0; //[Gy*mm^2]
            for (Int_t rIdx = 1; rIdx <= yAxis_lastbin; rIdx++) {
                double dA = TMath::Pi() * (pow(rzHist_Gy->GetYaxis()->GetBinUpEdge(rIdx), 2) -
                                           pow(rzHist_Gy->GetYaxis()->GetBinLowEdge(rIdx),2)); // [mm^2]
                avgDoseZ += rzHist_Gy->GetBinContent(zIdx,rIdx)*dA;
            }
            double volTot = TMath::Pi() * pow(r0_actual,2); //[mm^2]
            centerHist->SetBinContent(zIdx,avgDoseZ/volTot);
        }

}