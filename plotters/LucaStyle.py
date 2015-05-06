from ROOT import *

# global style for nice plots!
# histoType = 0: histos with errors via sumw2() ==> set markers, not lines!
def SetLucaStyle (histoType):
    #global LS
    LS = TStyle (gStyle) #copy some of the basics of defualt style...
    LS.SetName("LucaStyle")
    LS.SetTitle("Luca Style")

    # pad
    LS.SetOptStat(0)
    LS.SetTickLength(0.02,"X")
    LS.SetTickLength(0.02,"Y")
    LS.SetPadTickY(1)
    LS.SetPadTickX(1)
    LS.SetPadGridY(1);
    LS.SetPadGridX(1);
    #LS.SetPadTopMargin(0.05)
    #LS.SetPadBottomMargin(0.13)
    LS.SetPadLeftMargin(0.16)
    
    LS.SetCanvasDefH(600)
    LS.SetCanvasDefW(600)

    # axis
    LS.SetTitleYOffset(1.4)
    LS.SetTitleXOffset(0.9)
    LS.SetLabelOffset(0.009, "XYZ")
    LS.SetTitleSize(0.050, "XYZ")

    # legend
    LS.SetLegendBorderSize(0)
    LS.SetLegendFont(62)
    #LS.SetLegendFillStyle(0) #NOT IMPLEMENTED in root

    if histoType == 0:
        LS.SetMarkerStyle(8);
        LS.SetMarkerSize(0.8);

    LS.cd() 
    return LS

def SetColors (histo, Mcolor):
    histo.SetLineColor(Mcolor)
    histo.SetMarkerColor(Mcolor)
    return

def SetMarkers (histo, Mtype, Msize):
    histo.SetMarkerSize(Msize)
    histo.SetMarkerStyle(Mstyle)
    return

def SetAllProp (histo, Mcolor, Mtype, Msize):
    SetColors(histo, Mcolor)
    SetMarkers(histo, Mtype, Msize)