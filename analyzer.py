import ROOT, ROOT.gPad
import sys

names = []
files = []

protonEnergyHistos = []

for arg in sys.argv[1:]:
    print "Reading file '"+arg+"'..."
    names.append(arg)
    files.append(ROOT.TFile(arg,"READ"))
    print files[-1].ls()

    #Analyze
    
    # Protons with dpp/p < 90%
    protonAngle = files[-1].Get("protonAngle")
    protonEnergy = files[-1].Get("protonEnergy")
    
    protonEnergyHistos.append(protonEnergy)

    # c1 = ROOT.TCanvas()
    # protonAngle.Draw()
    # ROOT.gPad.SetLogy()
    # c2 = ROOT.TCanvas()
    # protonEnergy.Draw()
    # ROOT.gPad.SetLogy()
    # raw_input()
    
    print
    print "Integrals for protonEnergy:"
    xaxis = protonEnergy.GetXaxis()
    beamE = 7.0 #[TeV]
    upperE1 = beamE*0.9
    upperE1_bin = xaxis.FindBin(upperE1)
    upperE2 = beamE*0.8
    upperE2_bin = xaxis.FindBin(upperE2)
    lowE = 0.01
    lowE_bin = xaxis.FindBin(lowE)
    zeroE  = 0.0
    zeroE_bin = xaxis.FindBin(zeroE)
    
    print "Int0:"
    print "\t Range: all"
    int0 = protonEnergy.Integral()
    print "\t Integral =" + str(int0)

    print "Int1:"
    print "\t Range: {0}:{1} [TeV] -> {2}:{3} [{4}:{5}]".format(
        zeroE,upperE1,
        xaxis.GetBinLowEdge(zeroE_bin),xaxis.GetBinUpEdge(upperE1_bin),
        zeroE_bin,upperE1_bin
    )
    int1 = protonEnergy.Integral(zeroE_bin,upperE1_bin)
    print "\t int1=" + str(int1)
    
    print "Int2:"
    print "\t Range: {0}:{1} [TeV] -> {2}:{3} [{4}:{5}]".format(
        lowE,upperE1,
        xaxis.GetBinLowEdge(lowE_bin),xaxis.GetBinUpEdge(upperE1_bin),
        lowE_bin,upperE1_bin
    )
    int2 = protonEnergy.Integral(lowE_bin,upperE1_bin)
    print "\t int2=" + str(int2)
    
    print "Int3:"
    print "\t Range: {0}:{1} [TeV] -> {2}:{3} [{4}:{5}]".format(
        lowE,upperE2,
        xaxis.GetBinLowEdge(lowE_bin),xaxis.GetBinUpEdge(upperE2_bin),
        lowE_bin,upperE2_bin
    )
    int3 = protonEnergy.Integral(lowE_bin,upperE2_bin)
    print "\t int3=" + str(int3)
    
    print "Ratio int1/int0 =", int1/int0*100, "%"
    print "Ratio int2/int0 =", int2/int0*100, "%"
    print "Ratio int3/int0 =", int3/int0*100, "%"

print
print

for n,h,idx in zip(names, protonEnergyHistos,xrange(len(names))):
    hNew = None
    if idx == 0:
        hNew = h.DrawNormalized()
    else:
        hNew = h.DrawNormalized("SAME")
    hNew.SetLineColor(idx+1)

    #print idx,n,h
ROOT.gPad.SetLogy()
raw_input()
