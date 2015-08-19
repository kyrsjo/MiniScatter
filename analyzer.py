import ROOT#, ROOT.gPad
import sys

names = []
files = []

protonAngleHistos = []
protonEnergyHistos = []
targetEdep_IELHistos = []

for arg in sys.argv[1:]:
    print "Reading file '"+arg+"'..."
    names.append(arg)
    files.append(ROOT.TFile(arg,"READ"))
    print files[-1].ls()

    #Analyze
    
    # Protons with dpp/p < 90%
    protonAngle = files[-1].Get("protonAngle")
    protonEnergy = files[-1].Get("protonEnergy")
    targetEdep_IEL = files[-1].Get("targetEdep_IEL")

    protonAngleHistos.append(protonAngle)
    protonEnergyHistos.append(protonEnergy)
    targetEdep_IELHistos.append(targetEdep_IEL)
    
    print
    print "Integrals for protonEnergy:"
    xaxis = protonEnergy.GetXaxis()
    beamE = 7.0 #[TeV]
    upperE1 = beamE*0.9
    upperE1_bin = xaxis.FindBin(upperE1)
    upperE2 = beamE*0.8
    upperE2_bin = xaxis.FindBin(upperE2)

    zeroE  = 0.0
    zeroE_bin = xaxis.FindBin(zeroE)
    lowE1 = 0.01
    lowE1_bin = xaxis.FindBin(lowE1)
    lowE2 = 1.0
    lowE2_bin = xaxis.FindBin(lowE2)
    
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
        lowE1,upperE1,
        xaxis.GetBinLowEdge(lowE1_bin),xaxis.GetBinUpEdge(upperE1_bin),
        lowE1_bin,upperE1_bin
    )
    int2 = protonEnergy.Integral(lowE1_bin,upperE1_bin)
    print "\t int2=" + str(int2)
    
    print "Int3:"
    print "\t Range: {0}:{1} [TeV] -> {2}:{3} [{4}:{5}]".format(
        lowE1,upperE2,
        xaxis.GetBinLowEdge(lowE1_bin),xaxis.GetBinUpEdge(upperE2_bin),
        lowE1_bin,upperE2_bin
    )
    int3 = protonEnergy.Integral(lowE1_bin,upperE2_bin)
    print "\t int3=" + str(int3)

    print "Int4:"
    print "\t Range: {0}:{1} [TeV] -> {2}:{3} [{4}:{5}]".format(
        lowE2,upperE2,
        xaxis.GetBinLowEdge(lowE2_bin),xaxis.GetBinUpEdge(upperE2_bin),
        lowE2_bin,upperE2_bin
    )
    int4 = protonEnergy.Integral(lowE2_bin,upperE2_bin)
    print "\t int4=" + str(int4)

    
    print "Ratio int1/int0 =", int1/int0*100, "%"
    print "Ratio int2/int0 =", int2/int0*100, "%"
    print "Ratio int3/int0 =", int3/int0*100, "%"
    print "Ratio int4/int0 =", int4/int0*100, "%"

print
print

def plotHistos(histosList,doLog=True):
    c = ROOT.TCanvas()
    for n,h,idx in zip(names, histosList,xrange(len(names))):
        if h.GetSumOfWeights()==0.0:
            #Empty histogram
            continue
        hNew = None
        if idx == 0:
            hNew = h.DrawNormalized()
        else:
            hNew = h.DrawNormalized("SAME")
        hNew.SetLineColor(idx+1)
    if doLog:
        ROOT.gPad.SetLogy()
    return c;

c1=plotHistos(protonAngleHistos)
c2=plotHistos(protonEnergyHistos)
c3=plotHistos(targetEdep_IELHistos,False)
raw_input()


