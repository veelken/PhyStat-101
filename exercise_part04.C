
#include <RooAddPdf.h>
#include <RooBreitWigner.h>
#include <RooConstVar.h>
#include <RooDataHist.h>
#include <RooExponential.h>
#include <RooFFTConvPdf.h>
#include <RooFormulaVar.h>
#include <RooGaussian.h>
#include <RooLognormal.h>
#include <RooMinuit.h>
#include <RooPlot.h>
#include <RooProdPdf.h>
#include <RooRealVar.h>

#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1.h>
#include <TLine.h>
#include <TMath.h>
#include <TRandom.h>
#include <TRandom3.h>
#include <TROOT.h> // gROOT

#include <string>

TH1* bookHistogram(const std::string& histogramName, const std::string& histogramTitle, int numBins, double xMin, double xMax)
{
  TH1* histogram = new TH1D(histogramName.data(), histogramTitle.data(), numBins, xMin, xMax);

  // enable computation of bin-errors
  if ( !histogram->GetSumw2N() ) histogram->Sumw2();

  return histogram;
}

double square(double x)
{
  return x*x;
}

void fillHistogram(TH1* histogram, TRandom& rnd, double xMin, double xMax, int numEvents_signal, double mZ, double width, int numEvents_background, double lambda,
		   double eff_muon, double effErr_muon)
{
  //std::cout << "<fillHistogram>:" << std::endl;
  //std::cout << " xMin = " << xMin << std::endl;
  //std::cout << " xMax = " << xMax << std::endl;
  int numEvents = numEvents_signal + numEvents_background;
  assert(eff_muon > 0. && eff_muon <= 1. && effErr_muon >= 0. && effErr_muon <= eff_muon);
  double eff_muon_rnd = rnd.Gaus(eff_muon, effErr_muon);
  if ( eff_muon_rnd < 0.01 ) eff_muon_rnd = 0.01;
  if ( eff_muon_rnd > 1.00 ) eff_muon_rnd = 1.00;
   // signal yield is affected by muon identification efficiency^2, as each signal event contains two muons
  double numEvents_signal_rnd = square(eff_muon_rnd/eff_muon)*numEvents_signal; 
  // assume all muons in background events are genuine muons, so the background yield is affected by muon identification efficiency^2 in the same way as the signal
  double numEvents_background_rnd = square(eff_muon_rnd/eff_muon)*numEvents_background; 
  int numEvents_rnd = rnd.Poisson(numEvents_signal_rnd + numEvents_background_rnd);
  double p_signal = numEvents_signal/((double)numEvents);
  for ( int idxEvent = 0; idxEvent < numEvents_rnd; ++idxEvent ) {
    // keep iterating until an event within the given mass range xMin < x < xMax is obtained
    double x = -1.;
    do {
      // determine if next event is signal or background
      double u = rnd.Rndm();      
      if ( u < p_signal ) { // event is signal
	x = rnd.Gaus(mZ, width);
	//std::cout << "x_gauss = " << x << std::endl;
      } else { // event is background
	assert(lambda > 0.);
	x = rnd.Exp(1./lambda);
	x += xMin;
	//std::cout << "x_exp = " << x << std::endl;
      }
    } while ( !(x > xMin && x < xMax) );
    //std::cout << "filling histogram = " << histogram->GetName() << " with x = " << x << std::endl;
    histogram->Fill(x);
  }
}

void dumpHistogram(const TH1* histogram)
{
  std::cout << "dumping histogram = " << histogram->GetName() << std::endl;
  const TAxis* xAxis = histogram->GetXaxis();
  int numBins = xAxis->GetNbins();
  double sum = 0.;
  double sumErr2 = 0.;
  for ( int idxBin = 1; idxBin <= numBins; ++idxBin ) {
    double binCenter = xAxis->GetBinCenter(idxBin);
    double binContent = histogram->GetBinContent(idxBin);
    double binError = histogram->GetBinError(idxBin);
    std::cout << "bin #" << idxBin << " @ " << binCenter << ": content = " << binContent << " +/- " << binError << std::endl;
    sum += binContent;
    sumErr2 += (binError*binError);
  }
  std::cout << "sum = " << sum << " +/- " << TMath::Sqrt(sumErr2) << std::endl;
}

double normal(double x, double mZ, double sigma)
{
  return (1./(TMath::Sqrt(2.*TMath::Pi())*sigma))*TMath::Gaus(x, mZ, sigma);
}

using namespace RooFit;

void exercise_part04()
{
  // prevent pop-up windows
  gROOT->SetBatch(true);

  // reduce RooFit print-out
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

  // define binning
  // (30 bins from 60 to 120 GeV)
  int numBins = 30;
  double xMin =  60.; // GeV
  double xMax = 120.; // GeV

  // define mass and (intrinsic) width of Z boson
  // (taken from http://pdg.lbl.gov/2018/listings/rpp2018-list-z-boson.pdf )
  double mZ    = 91.2; // GeV
  double width =  4.5; // GeV (includes intrinsic width of Z boson + experimental resolution)

  int numEvents_signal = 1000;
  int numEvents_background = 9000;
  double lambda = 0.1; // parameter of exp(-lambda*x) distribution for background

  double eff_muon = 0.80;
  double effErr_muon = 0.05*eff_muon; // uncertainty on muon identification efficiency

  TRandom3 rnd;
  int numToys = 1000;

  TH1* histogram_numEvents_signal = bookHistogram("histogram_numEvents_signal", "numEvents_signal", 200, 0., 2.*numEvents_signal);
  TH1* histogram_numEventsErr_signal = bookHistogram("histogram_numEventsErr_signal", "numEventsErr_signal", 200, 0., 0.2*numEvents_signal);
  TH1* histogram_pull_signal = bookHistogram("histogram_pull_signal", "pull_signal", 200, -10., +10.);
  TH1* histogram_numEvents_background = bookHistogram("histogram_numEvents_background", "numEvents_background", 200, 0., 2.*numEvents_background);
  TH1* histogram_numEventsErr_background = bookHistogram("histogram_numEventsErr_background", "numEventsErr_background", 200, 0., 0.2*numEvents_background);
  TH1* histogram_pull_background = bookHistogram("histogram_pull_background", "pull_background", 200, -10., +10.);
  TH1* histogram_sf_muon = bookHistogram("histogram_sf_muon", "sf_muon", 200, 0., 2.);

  for ( int idxToy = 0; idxToy < numToys; ++idxToy ) {
    if ( (idxToy % 100) == 0 ) {
      std::cout << "processing toy #" << idxToy << "." << std::endl;
    }

    // generate (pseudo)data
    TH1* histogram_mass = bookHistogram("histogram_mass", "mass", numBins, xMin, xMax);
    fillHistogram(histogram_mass, rnd, xMin, xMax, numEvents_signal, mZ, width, numEvents_background, lambda, eff_muon, effErr_muon);
    //dumpHistogram(histogram_mass);
    RooRealVar observable("observable", "Mass [GeV]", xMin, xMax);
    RooDataHist data("data", "(Pseudo)data", RooArgList(observable), histogram_mass);

    // define fit model
    RooRealVar mean_gauss("mean_gauss", "Gaussian mean", mZ, 0.8*mZ, 1.2*mZ);
    RooRealVar sigma_gauss("sigma_gauss", "Gaussian resolution", width, 0., 2.*width);
    RooGaussian model_S("model_S", "Signal model", observable, mean_gauss, sigma_gauss);

    RooRealVar lambda_exp("lambda_exp", "Background shape", -lambda, -10.*lambda, 0.);
    RooExponential model_B("pdf_B", "Background model", observable, lambda_exp);

    RooRealVar norm_S("norm_S", "Signal yield", numEvents_signal, 0., 2.*histogram_mass->Integral());
    RooRealVar norm_B("norm_B", "Background yield", numEvents_background, 0., 2.*histogram_mass->Integral());
    //-----------------------------------------------------------------------------------------------
    // create constraint term in likelihood function 
    // for systematic uncertainty on muon identification efficiency
    RooRealVar sf_muon("sf_muon", "Data/MC scale-factor for muon identification efficiency", 1., 0., 1./eff_muon);
    RooConstVar const0("const0", "Constant of value 1", 1.); 
    RooConstVar const1("const1", "Uncertainty on data/MC scale-factor for muon identification efficiency", effErr_muon/eff_muon);
    //RooLognormal constraint_muon("constraint_muon", "Constraint on data/MC scale-factor for muon identification efficiency", sf_muon, const0, 1. + const1);
    RooGaussian constraint_muon("constraint_muon", "Constraint on data/MC scale-factor for muon identification efficiency", sf_muon, const0, const1);
    //-----------------------------------------------------------------------------------------------
    RooFormulaVar norm_S_wMuonEff("norm_S_wMuonEff", "Signal yield*muon identification efficiency^2","norm_S*sf_muon*sf_muon", RooArgList(norm_S, sf_muon));
    RooFormulaVar norm_B_wMuonEff("norm_B_wMuonEff", "Background yield*muon identification efficiency^2","norm_B*sf_muon*sf_muon", RooArgList(norm_B, sf_muon));
    RooAddPdf model_SplusB("model_SplusB", "Signal+background model", RooArgList(model_S, model_B), RooArgList(norm_S_wMuonEff, norm_B_wMuonEff));
    RooProdPdf model_SplusB_wMuonEff("model_SplusB_wMuonEff", "Signal+background model with muon identification efficiency uncertainty", model_SplusB, constraint_muon);
    
    // perform fit
    model_SplusB_wMuonEff.fitTo(data, PrintLevel(-1));

    // print signal and background yields determined by fit
    //std::cout << "fit results:" << std::endl;
    //std::cout << " S = " << norm_S.getVal();
    double err_S;
    if ( norm_S.hasAsymError() ) {
      //std::cout << " + " << norm_S.getErrorHi() << " - " << norm_S.getErrorLo();
      if ( norm_S.getVal() > numEvents_signal ) {
	err_S = norm_S.getErrorHi();
      } else {
	err_S = norm_S.getErrorLo();
      }
    } else {
      //std::cout << " +/- " << norm_S.getError();
      err_S = norm_S.getError();
    }
    //std::cout << std::endl;
    double pull_S = (norm_S.getVal() - numEvents_signal)/err_S;
    //std::cout << " B = " << norm_B.getVal();
    double err_B;
    if ( norm_B.hasAsymError() ) {
      //std::cout << " + " << norm_B.getErrorHi() << " - " << norm_B.getErrorLo();
      if ( norm_B.getVal() > numEvents_background ) {
	err_B = norm_B.getErrorHi();
      } else {
	err_B = norm_B.getErrorLo();
      }
    } else {
      //std::cout << " +/- " << norm_B.getError();
      err_B = norm_B.getError();
    }
    //std::cout << std::endl;
    double pull_B = (norm_B.getVal() - numEvents_background)/err_B;

    histogram_numEvents_signal->Fill(norm_S.getVal());
    histogram_numEventsErr_signal->Fill(err_S);
    histogram_pull_signal->Fill(pull_S);
    histogram_numEvents_background->Fill(norm_B.getVal());
    histogram_numEventsErr_background->Fill(err_B);
    histogram_pull_background->Fill(pull_B);
    histogram_sf_muon->Fill(sf_muon.getVal());

    // make control plots (1st toy only)
    if ( idxToy == 0 ) {
      // show (pseudo)data versus fit model,
      // for the values of signal and background yields determined by fit
      RooPlot* frame1 = observable.frame();
      data.plotOn(frame1);
      model_SplusB_wMuonEff.plotOn(frame1);
      model_SplusB_wMuonEff.plotOn(frame1, Components(model_B), LineStyle(ELineStyle::kDashed));
      TCanvas canvas1;
      canvas1.cd();
      frame1->Draw();      
      canvas1.SaveAs("exercise_part04_fit.png");
      
      // show likelihood function
      RooAbsReal* nll = model_SplusB_wMuonEff.createNLL(data, NumCPU(4));
      RooMinuit minuit(*nll);
      minuit.setPrintLevel(-1);
      minuit.migrad();      
      RooPlot* frame2 = norm_S.frame(Bins(10), Range(0.5*numEvents_signal, 1.5*numEvents_signal));
      nll->plotOn(frame2, ShiftToZero());
      RooAbsReal* nll_profiled = nll->createProfile(norm_S);
      nll_profiled->plotOn(frame2, LineColor(kRed));
      frame2->SetMinimum(0.);
      frame2->SetMaximum(5.);
      TCanvas canvas2;
      canvas2.cd();
      TGraph* line_1sigma = new TGraph(2);
      line_1sigma->SetPoint(0, 0.5*numEvents_signal, 1.);
      line_1sigma->SetPoint(1, 1.5*numEvents_signal, 1.);
      line_1sigma->SetLineColor(8);
      line_1sigma->SetLineWidth(1);
      line_1sigma->SetLineStyle(7);
      frame2->addObject(line_1sigma);
      TGraph* line_2sigma = new TGraph(2);
      line_2sigma->SetPoint(0, 0.5*numEvents_signal, 4.);
      line_2sigma->SetPoint(1, 1.5*numEvents_signal, 4.);
      line_2sigma->SetLineColor(8);
      line_2sigma->SetLineWidth(1);
      line_2sigma->SetLineStyle(7);
      frame2->addObject(line_2sigma, "L");
      frame2->Draw();
      canvas2.SaveAs("exercise_part04_nll.png");
      //delete line_1sigma;
      //delete line_2sigma;
    }

    delete histogram_mass;
  }

  // write histograms to ROOT file
  TFile* outputFile = new TFile("exercise_part04.root", "RECREATE");
  outputFile->cd();
  histogram_numEvents_signal->Write();
  histogram_numEventsErr_signal->Write();
  histogram_pull_signal->Write();
  histogram_numEvents_background->Write();
  histogram_numEventsErr_background->Write();
  histogram_pull_background->Write();
  histogram_sf_muon->Write();
  delete outputFile;
}
