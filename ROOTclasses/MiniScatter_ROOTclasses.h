#ifndef MiniScatter_ROOTclasses_h
#define MiniScatter_ROOTclasses_h

#include "TH2D.h"

//Scale a 2D r-z histogram's bincontent by 1/dV,
// as can be used by plotRZgray()
void MSR_plotRZgray_histoScalerVolume(TH2D* rzScaled);

#endif