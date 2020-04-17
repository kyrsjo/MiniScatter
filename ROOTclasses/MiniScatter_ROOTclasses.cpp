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
