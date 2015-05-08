// fit the given histo with a gaussian in the range
// mean-rms, mean+rms
// because the fit parameters are not always so good
// also check against fit

void GetCoreGaussianSigma(TString fileName, TString histoName, double lim1, double lim2)
{
    TFile* fIn = new TFile (fileName);
    TH1D* h = (TH1D*)fIn->Get(histoName);

    // get core function
    TF1* fCBFit = h->GetFunction("CBFuncAsymm");
    cout << "From CB fit:    sigma = " << fCBFit->GetParameter(5) << " +/- ";
    if (TMath::IsNaN(fCBFit->GetParError(5) )) cout << "NaN" << endl;
    else cout << fCBFit->GetParError(5) << endl;

    // now fit with gaussian
    double mean = h->GetMean();
    double RMS = h->GetRMS();

    TCanvas* c1 = new TCanvas;
    //TF1* fGaus = new TF1 ("fGaus", "gaus", mean-RMS, mean+0.8*RMS);
    TF1* fGaus = new TF1 ("fGaus", "gaus", lim1, lim2);
    h->Fit("fGaus", "R");
    cout << "In gaussian core fit:   sigma =  " << fGaus->GetParameter(2) << " +/- " << fGaus->GetParError(2) << endl;
    h->Draw();
}