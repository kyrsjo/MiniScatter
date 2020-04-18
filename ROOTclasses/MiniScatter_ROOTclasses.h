#ifndef MiniScatter_ROOTclasses_h
#define MiniScatter_ROOTclasses_h

#include "TH2D.h"

//Scale a 2D r-z histogram's bincontent by 1/dV,
// as can be used by plotRZgray()
void MSR_plotRZgray_histoScalerVolume(TH2D* rzScaled);

//Integrate a 2D rzHistogram as produced by plotRZgray() to get average dose [Gy] as a function of z-position
void MSR_plotZgray_sumR(TH2D* rzHist_Gy, TH1D* centerHist, Int_t yAxis_lastbin, Double_t r0_actual);

#endif