import ROOT
import operator
import array


ROOT.gSystem.Load('libRooFit')


class TurnOn:
    def __init__(self, **args):
        self.name        = args.get("Name", "turnon")
        self.legend      = args.get("Legend","Turn-on")
        self.histo       = args.get("Histo", None)
        self.fit         = args.get("Fit", None)
        self.markerColor = args.get("MarkerColor", ROOT.kBlack)
        self.markerStyle = args.get("MarkerStyle", 20)
        self.lineColor   = args.get("LineColor", ROOT.kBlack)
        self.lineStyle   = args.get("LineStyle", 1)
        self.histo.SetName(self.name+"_histo")
        self.fit.SetName(self.name+"_fit")



class TurnOnPlot:
    def __init__(self):
        self.name  = ""
        self.turnons = []
        self.plotDir = "plots/"
        self.xRange = (10, 120)
        self.xTitle = "p_{T}^{offl} [GeV]"
        #self.legendPosition = (0.6,0.2,0.9,0.4)
        self.legendPosition = (0.4,0.2,0.9,0.6)
        self.setPlotStyle()

    def addTurnOn(self, turnon):
        self.turnons.append(turnon)

    def plot(self):
        canvas = ROOT.TCanvas("c_"+self.name, self.name, 800, 800)
        canvas.SetGrid()
        hDummy = ROOT.TH1F("hDummy_"+self.name, self.name, 1, self.xRange[0], self.xRange[1])
        hDummy.SetAxisRange(0, 1.05, "Y")
        hDummy.SetXTitle(self.xTitle)
        hDummy.SetYTitle("Efficiency")
        hDummy.Draw()
        legend = ROOT.TLegend(self.legendPosition[0],self.legendPosition[1],self.legendPosition[2],self.legendPosition[3])
        legend.SetTextFont(42)
        legend.SetFillColor(0)
        legend1 = ROOT.TLegend(0.14, 0.85, 0.80, 1.02)
        legend1.SetBorderSize(0)
        legend1.SetTextFont(62)
        legend1.SetTextSize(0.025)
        legend1.SetLineColor(0)
        legend1.SetLineStyle(1)
        legend1.SetLineWidth(1)
        legend1.SetFillColor(0)
        legend1.SetFillStyle(0)
        legend1.AddEntry("NULL","CMS Simulation: gg #rightarrow H #rightarrow #tau #tau         #sqrt{s}=13 TeV - bx=25ns - PU40","h")
        legend1.AddEntry("NULL","L1 Threshold : 30 GeV","h")
        
        for turnon in self.turnons:
            histo = turnon.histo
            histo.SetMarkerStyle(turnon.markerStyle)
            histo.SetMarkerColor(turnon.markerColor)
            histo.SetLineColor(turnon.markerColor)
            fit = turnon.fit
            fit.SetLineStyle(turnon.lineStyle)
            fit.SetLineColor(turnon.lineColor)
            fit.SetLineWidth(2)
            histo.Draw("p same")
            fit.Draw("l same")
            # legends
            legend.AddEntry(histo, turnon.legend, "pe")
            legend.Draw()
            #if self.name=="turnon_Stage1_Stage2_EB":
        legend1.Draw()
        #print ("DEBUG: " + self.plotDir+"/"+self.name+".eps")
        canvas.Print(self.plotDir+"/"+self.name+".pdf")
        canvas.Print(self.plotDir+"/"+self.name+".png")
        return canvas


    def setPlotStyle(self):
        ROOT.gROOT.SetStyle("Plain")
        ROOT.gStyle.SetOptStat()
        ROOT.gStyle.SetOptFit(0)
        ROOT.gStyle.SetOptTitle(0)
        ROOT.gStyle.SetFrameLineWidth(1)
        ROOT.gStyle.SetPadBottomMargin(0.13)
        ROOT.gStyle.SetPadLeftMargin(0.15)
        ROOT.gStyle.SetPadTopMargin(0.06)
        ROOT.gStyle.SetPadRightMargin(0.05)

        ROOT.gStyle.SetLabelFont(42,"X")
        ROOT.gStyle.SetLabelFont(42,"Y")
        ROOT.gStyle.SetLabelSize(0.04,"X")
        ROOT.gStyle.SetLabelSize(0.04,"Y")
        ROOT.gStyle.SetLabelOffset(0.01,"Y")
        ROOT.gStyle.SetTickLength(0.02,"X")
        ROOT.gStyle.SetTickLength(0.02,"Y")
        ROOT.gStyle.SetLineWidth(1)
        ROOT.gStyle.SetTickLength(0.02 ,"Z")

        ROOT.gStyle.SetTitleSize(0.1)
        ROOT.gStyle.SetTitleFont(42,"X")
        ROOT.gStyle.SetTitleFont(42,"Y")
        ROOT.gStyle.SetTitleSize(0.05,"X")
        ROOT.gStyle.SetTitleSize(0.05,"Y")
        ROOT.gStyle.SetTitleOffset(1.1,"X")
        ROOT.gStyle.SetTitleOffset(1.4,"Y")
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetPalette(1)
        ROOT.gStyle.SetPaintTextFormat("3.2f")
        ROOT.gROOT.ForceStyle()
