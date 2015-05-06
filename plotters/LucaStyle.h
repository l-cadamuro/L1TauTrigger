// global style for nice plots!
// histoType = 0: histos with errors via sumw2() ==> set markers, not lines!
TStyle* SetLucaStyle (int histoType = 0)
{
    TStyle *LS = new TStyle(*gStyle)
    LS->SetName("LucaStyle");
    LS->SetTitle("Luca Style");

    // pad
    LS ->SetOptStat(0);
    LS->SetTickLength(0.02,"X");
    LS->SetTickLength(0.02,"Y");
    LS->SetPadTickY(1);
    LS->SetPadTickX(1);
    LS->SetPadGridY(1);
    LS->SetPadGridX(1);
    //LS->SetPadTopMargin(0.05);
    //LS->SetPadBottomMargin(0.13);
    LS->SetPadLeftMargin(0.16);
    
    // axis
    LS->SetTitleYOffset(1.4);
    LS->SetTitleXOffset(0.9);
    LS->SetLabelOffset(0.009, "XYZ");
    LS->SetTitleSize(0.050, "XYZ")

    LS->SetCanvasDefH(600);
    LS->SetCanvasDefW(600);

    // legend
    LS->SetLegendBorderSize(0);
    //LS->SetLegendFillStyle(0); // NOT AVAILABLE in root, need to set by hand

    // histos
    if (histoType == 0)
    {
        LS->SetMarkerStyle(8);
        LS->SetMarkerSize(0.8);
    }

    LS->cd(); // set as the current style
    return LS;
}

void SetColors (TH1* histo, int color)
{
    histo->SetLineColor(color);
    histo->SetMarkerColor(color);
    return;
}

void SetMarkers (TH1* histo, int type, int size)
{
    histo->SetMarkerSize(size);
    histo->SetMarkerType(type);
    return;
}

void SetAllProp (TH1* histo, int color, int MarType, int MarSize)
{
    SetMarkers(histo, MarType, MarSize);
    SetColors(histo, color);
    return;
}