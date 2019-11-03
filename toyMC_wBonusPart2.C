
#include <RooAbsPdf.h>
#include <RooAddPdf.h>
#include <RooBreitWigner.h>
#include <RooConstVar.h>
#include <RooDataHist.h>
#include <RooDataSet.h>
#include <RooExponential.h>
#include <RooFFTConvPdf.h>
#include <RooGaussian.h>
#include <RooMinuit.h>
#include <RooPlot.h>
#include <RooRealVar.h>
#include <RooVoigtian.h>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TMath.h>
#include <TRandom.h>
#include <TRandom3.h>
#include <TROOT.h>  // gROOT

#include <string>
#include <assert.h> // assert

// define option to randomly sample from mass distribution for signal
//   1 = use 2-step procedure to sample from Breit-Wigner + Gaussian distribution
//   2 = sample directly from Voigtian distribution
const int mode_signal_fill = 2;

// define option to fit mass distribution for signal
//   1 = use Fourier transform to compute convolution of Breit-Wigner + Gaussian
//   2 = fit Voigtian function
const int mode_signal_fit = 2;

TH1* bookHistogram(const std::string& histogramName, const std::string& histogramTitle, int numBins, double xMin, double xMax)
{
  TH1* histogram = new TH1D(histogramName.data(), histogramTitle.data(), numBins, xMin, xMax);

  // enable computation of bin-errors
  if ( !histogram->GetSumw2N() ) histogram->Sumw2();

  return histogram;
}

void fillHistogram(TH1* histogram, TRandom& rnd, double xMin, double xMax, int numEvents_signal, double mZ, double width, double sigma, int numEvents_background, double lambda)
{
  //std::cout << "<fillHistogram>:" << std::endl;
  //std::cout << " xMin = " << xMin << std::endl;
  //std::cout << " xMax = " << xMax << std::endl;
  int numEvents = numEvents_signal + numEvents_background;
  int numEvents_rnd = rnd.Poisson(numEvents);
  double p_signal = numEvents_signal/((double)numEvents);
  //-----------------------------------------------------------------------------
  // define temporary variables used to randomly sample from Voigtian distribution
  RooRealVar observable("observable", "Mass [GeV]", xMin, xMax);
  RooRealVar mean_bw("mean_bw", "Breit-Wigner mean", mZ, 0.8*mZ, 1.2*mZ);
  RooRealVar width_bw("width_bw", "Breit-Wigner width", width, 0., 2.*width);
  RooRealVar sigma_gauss("sigma_gauss", "Gaussian resolution", sigma, 0., 2.*sigma);
  RooVoigtian model_S("model_S", "Signal model", observable, mean_bw, width_bw, sigma_gauss);
  //-----------------------------------------------------------------------------
  for ( int idxEvent = 0; idxEvent < numEvents_rnd; ++idxEvent ) {
    // keep iterating until an event within the given mass range xMin < x < xMax is obtained
    double x = -1.;
    do {
      // determine if next event is signal or background
      double u = rnd.Rndm();      
      if ( u < p_signal ) { // event is signal
	if ( mode_signal_fill == 1 ) {
	  // draw x from Breit-Wigner distribution
	  x = rnd.BreitWigner(mZ, width);
	  //std::cout << "x_bw = " << x << std::endl;
	  // smear x to account for experimental resolution
	  x = rnd.Gaus(x, sigma);
	  //std::cout << "x_gauss = " << x << std::endl;
	} else if ( mode_signal_fill == 2 ) {
	  // draw x from Voigtian distribution
	  RooDataSet* dataset = model_S.generate(observable, 1);
	  TH1* tmpHistogram = dataset->createHistogram("tmpHistogram", observable, RooFit::Binning(histogram->GetNbinsX(), xMin, xMax));	
	  assert(tmpHistogram->GetEntries() == 1);
	  x = tmpHistogram->GetMean();
	  //std::cout << "x_voigtian = " << tmpHistogram->GetMean() << std::endl;
	  delete dataset;
	  delete tmpHistogram;
	} else {
	  cerr << "Invalid parameter 'mode_signal_fill' = " << mode_signal_fill << " !!" << std::endl;
	  assert(0);
	}
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

void toyMC_wBonusPart2()
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
  double width =  2.5; // GeV 

  // define experimental resolution on Z-boson mass
  double sigma =  2.;  // GeV

  int numEvents_signal = 1000;
  int numEvents_background = 9000;
  double lambda = 0.1; // parameter of exp(-lambda*x) distribution for background

  TRandom3 rnd;
  int numToys = 1000;

  TH1* histogram_numEvents_signal = bookHistogram("histogram_numEvents_signal", "numEvents_signal", 200, 0., 2.*numEvents_signal);
  TH1* histogram_numEventsErr_signal = bookHistogram("histogram_numEventsErr_signal", "numEventsErr_signal", 200, 0., 0.2*numEvents_signal);
  TH1* histogram_pull_signal = bookHistogram("histogram_pull_signal", "pull_signal", 200, -10., +10.);
  TH1* histogram_numEvents_background = bookHistogram("histogram_numEvents_background", "numEvents_background", 200, 0., 2.*numEvents_background);
  TH1* histogram_numEventsErr_background = bookHistogram("histogram_numEventsErr_background", "numEventsErr_background", 200, 0., 0.2*numEvents_background);
  TH1* histogram_pull_background = bookHistogram("histogram_pull_background", "pull_background", 200, -10., +10.);

  for ( int idxToy = 0; idxToy < numToys; ++idxToy ) {
    if ( (idxToy % 100) == 0 ) {
      std::cout << "processing toy #" << idxToy << "." << std::endl;
    }

    // generate (pseudo)data
    TH1* histogram_mass = bookHistogram("histogram_mass", "mass", numBins, xMin, xMax);
    fillHistogram(histogram_mass, rnd, xMin, xMax, numEvents_signal, mZ, width, sigma, numEvents_background, lambda);
    //dumpHistogram(histogram_mass);
    RooRealVar observable("observable", "Mass [GeV]", xMin, xMax);
    RooDataHist data("data", "(Pseudo)data", RooArgList(observable), histogram_mass);

    // define fit model
    RooRealVar mean_bw("mean_bw", "Breit-Wigner mean", mZ, 0.8*mZ, 1.2*mZ);
    RooRealVar width_bw("width_bw", "Breit-Wigner width", width, 0., 2.*width);
    RooBreitWigner pdf_bw("pdf_bw", "Breit-Wigner PDF", observable, mean_bw, width_bw);
    RooConstVar mean_gauss("mean_gauss", "Gaussian mean", 0.);
    RooRealVar sigma_gauss("sigma_gauss", "Gaussian resolution", sigma, 0., 2.*sigma);
    RooGaussian pdf_gauss("pdf_gauss", "Gaussian PDF", observable, mean_gauss, sigma_gauss);
    RooAbsPdf* model_S = 0;
    if ( mode_signal_fit == 1 ) {
      model_S = new RooFFTConvPdf("model_S", "Signal model", observable, pdf_bw, pdf_gauss);
    } else if ( mode_signal_fit == 2 ) {
      model_S = new RooVoigtian("model_S", "Signal model", observable, mean_bw, width_bw, sigma_gauss);
    } else {
      cerr << "Invalid parameter 'mode_signal_fill' = " << mode_signal_fill << " !!" << std::endl;
      assert(0);
    }

    RooRealVar lambda_exp("lambda_exp", "Background shape", -lambda, -10.*lambda, 0.);
    RooExponential model_B("pdf_B", "Background model", observable, lambda_exp);

    RooRealVar norm_S("norm_S", "Signal yield", numEvents_signal, 0., 2.*histogram_mass->Integral());
    RooRealVar norm_B("norm_B", "Background yield", numEvents_background, 0., 2.*histogram_mass->Integral());
    RooAddPdf model_SplusB("model_SplusB", "Signal+background model", RooArgList(*model_S, model_B), RooArgList(norm_S, norm_B));
    
    // perform fit
    model_SplusB.fitTo(data, PrintLevel(-1));

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

    // make control plots (1st toy only)
    if ( idxToy == 0 ) {
      // show (pseudo)data versus fit model,
      // for the values of signal and background yields determined by fit
      RooPlot* frame1 = observable.frame();
      data.plotOn(frame1);
      model_SplusB.plotOn(frame1);
      model_SplusB.plotOn(frame1, Components(model_B), LineStyle(ELineStyle::kDashed));
      TCanvas canvas1;
      canvas1.cd();
      frame1->Draw();      
      canvas1.SaveAs("toyMC_fit.png");
      
      // show likelihood function
      RooAbsReal* nll = model_SplusB.createNLL(data, NumCPU(4));
      RooMinuit minuit(*nll);
      minuit.setPrintLevel(-1);
      minuit.migrad();      
      RooPlot* frame2 = norm_S.frame(Bins(10), Range(0.5*numEvents_signal, 1.5*numEvents_signal));
      nll->plotOn(frame2, ShiftToZero());
      RooAbsReal* nll_profiled = nll->createProfile(norm_S);
      nll_profiled->plotOn(frame2, LineColor(kRed));
      frame2->SetMinimum(0.);
      frame2->SetMaximum(25.);
      TCanvas canvas2;
      canvas2.cd();
      frame2->Draw();
      canvas2.SaveAs("toyMC_nll.png");
    }

    delete histogram_mass;
    delete model_S;
  }

  // write histograms to ROOT file
  TFile* outputFile = new TFile("toyMC.root", "RECREATE");
  outputFile->cd();
  histogram_numEvents_signal->Write();
  histogram_numEventsErr_signal->Write();
  histogram_pull_signal->Write();
  histogram_numEvents_background->Write();
  histogram_numEventsErr_background->Write();
  histogram_pull_background->Write();
  delete outputFile;
}
