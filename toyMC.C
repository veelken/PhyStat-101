
#include <RooAddPdf.h>
#include <RooBreitWigner.h>
#include <RooDataHist.h>
#include <RooExponential.h>
#include <RooFFTConvPdf.h>
#include <RooGaussian.h>
#include <RooMinuit.h>
#include <RooPlot.h>
#include <RooRealVar.h>

#include <TCanvas.h>
#include <TFile.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <TH1.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TPad.h>
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

void fillHistogram(TH1* histogram, TRandom& rnd, double xMin, double xMax, int numEvents_signal, double mZ, double width, int numEvents_background, double lambda)
{
  //std::cout << "<fillHistogram>:" << std::endl;
  //std::cout << " xMin = " << xMin << std::endl;
  //std::cout << " xMax = " << xMax << std::endl;
  int numEvents = numEvents_signal + numEvents_background;
  int numEvents_rnd = rnd.Poisson(numEvents);
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

double sum_fcn(double* x, double* par)
{
  return par[0]*normal(x[0], par[1], par[2]) + par[3]*TMath::Exp(-par[4]*(x[0] - par[5]));
}

void fillHistogram_gauss(TH1* histogram, TRandom& rnd, int numToys, double xMin, double xMax, double mZ, double sigma)
{
  double yMax = normal(mZ, mZ, sigma);
  int idxToy = 0;
  while ( idxToy < numToys ) {
    bool isSelectedToy = false;
    double x;
    while ( !isSelectedToy ) {
      x = xMin + (xMax - xMin)*rnd.Rndm();
      double y = rnd.Rndm()*yMax;
      if ( y < normal(x, mZ, sigma) ) {
	isSelectedToy = true;
      }
    }
    histogram->Fill(x);
    ++idxToy;
  }
}

void fillHistogram_exp(TH1* histogram, TRandom& rnd, int numToys, double xMin, double xMax, double lambda)
{
  double yMax = 1.;
  int idxToy = 0;
  while ( idxToy < numToys ) {
    bool isSelectedToy = false;
    double x;
    while ( !isSelectedToy ) {
      x = xMin + (xMax - xMin)*rnd.Rndm();
      double y = rnd.Rndm()*yMax;
      if ( y < TMath::Exp(-lambda*(x - xMin)) ) {
	isSelectedToy = true;
      }
    }
    histogram->Fill(x);
    ++idxToy;
  }
}

//---------------------------------------------------------------------------------------------------
TH1* compRatioHistogram(const std::string& histogramName_ratio, const TH1* histogram_numerator, const TH1* histogram_denominator)
{
  TH1* histogram_ratio = 0;
  if ( histogram_numerator->GetDimension() == histogram_denominator->GetDimension() &&
       histogram_numerator->GetNbinsX()    == histogram_denominator->GetNbinsX()    ) {
    histogram_ratio = (TH1*)histogram_numerator->Clone(histogramName_ratio.data());
    histogram_ratio->Reset();
    if ( !histogram_ratio->GetSumw2N() ) histogram_ratio->Sumw2();

    int numBins = histogram_denominator->GetNbinsX();
    for ( int idxBin = 0; idxBin <= (numBins + 1); ++idxBin ){
      double binContent_numerator = histogram_numerator->GetBinContent(idxBin);
      double binError_numerator = histogram_numerator->GetBinError(idxBin);

      double binContent_denominator = histogram_denominator->GetBinContent(idxBin);
      double binError_denominator = histogram_denominator->GetBinError(idxBin);

      histogram_ratio->SetBinContent(idxBin, (binContent_numerator - binContent_denominator)/binContent_denominator);
      histogram_ratio->SetBinError(idxBin, binError_numerator/binContent_denominator);
    }
    
    histogram_ratio->SetLineColor(histogram_numerator->GetLineColor());
    histogram_ratio->SetLineWidth(histogram_numerator->GetLineWidth());
    histogram_ratio->SetMarkerColor(histogram_numerator->GetMarkerColor());
    histogram_ratio->SetMarkerStyle(histogram_numerator->GetMarkerStyle());
    histogram_ratio->SetMarkerSize(histogram_numerator->GetMarkerSize());
  }

  return histogram_ratio;
}

TGraphAsymmErrors* convertToGraph(const TH1* histogram, double offsetX = 0.15)
{
  int numPoints = histogram->GetNbinsX();
  TGraphAsymmErrors* graph = new TGraphAsymmErrors(numPoints);
  for ( int idxPoint = 0; idxPoint < numPoints; ++idxPoint ) {
    double idxBin = idxPoint + 1;
    double binCenter = histogram->GetBinCenter(idxBin);
    double binWidth = histogram->GetBinWidth(idxBin);
    assert(binWidth > offsetX);
    double binContent = histogram->GetBinContent(idxBin);
    double binError = histogram->GetBinError(idxBin);
    graph->SetPoint(idxPoint, binCenter + offsetX, binContent);
    graph->SetPointError(idxPoint, 0.5*binWidth + offsetX, 0.5*binWidth - offsetX, binError, binError);
  }
  graph->SetLineColor(histogram->GetLineColor());
  graph->SetLineWidth(histogram->GetLineWidth());
  graph->SetMarkerColor(histogram->GetMarkerColor());
  graph->SetMarkerStyle(histogram->GetMarkerStyle());
  graph->SetMarkerSize(histogram->GetMarkerSize());
  return graph;
}

void showHistograms(double canvasSizeX, double canvasSizeY,
		    TH1* histogram_ref, const std::string& legendEntry_ref,
		    TH1* histogram1, const std::string& legendEntry1,
		    TH1* histogram2, const std::string& legendEntry2,
		    TH1* histogram3, const std::string& legendEntry3,
		    TF1* function, const std::string& legendEntry_function,
		    const std::string& xAxisTitle, double xAxisOffset,
		    bool useLogScale, double yMin, double yMax, const std::string& yAxisTitle, double yAxisOffset,
		    double legendX0, double legendY0, 
		    const std::string& outputFileName)
{
  assert(histogram_ref);

  TCanvas* canvas = new TCanvas("canvas", "canvas", canvasSizeX, canvasSizeY);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2);
  canvas->SetLeftMargin(0.12);
  canvas->SetBottomMargin(0.12);

  TPad* topPad = new TPad("topPad", "topPad", 0.00, 0.35, 1.00, 1.00);
  topPad->SetFillColor(10);
  topPad->SetTopMargin(0.04);
  topPad->SetLeftMargin(0.15);
  topPad->SetBottomMargin(0.03);
  topPad->SetRightMargin(0.05);
  topPad->SetLogy(useLogScale);

  TPad* bottomPad = new TPad("bottomPad", "bottomPad", 0.00, 0.00, 1.00, 0.35);
  bottomPad->SetFillColor(10);
  bottomPad->SetTopMargin(0.02);
  bottomPad->SetLeftMargin(0.15);
  bottomPad->SetBottomMargin(0.24);
  bottomPad->SetRightMargin(0.05);
  bottomPad->SetLogy(false);

  canvas->cd();
  topPad->Draw();
  topPad->cd();

  int colors[6] = { kBlack, kGreen - 6, kBlue - 7,  kMagenta - 7, kCyan - 6, kRed - 6 };
  int markerStyles[6] = { 24, 21, 20, 21, 22, 23 };
  int markerSizes[6] = { 1, 1, 1, 1, 1, 1 };

  TLegend* legend = new TLegend(legendX0, legendY0, legendX0 + 0.26, legendY0 + 0.20, "", "brNDC"); 
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetTextSize(0.045);

  histogram_ref->SetTitle("");
  histogram_ref->SetStats(false);
  histogram_ref->SetMinimum(yMin);
  histogram_ref->SetMaximum(yMax);
  histogram_ref->SetLineColor(colors[0]);
  histogram_ref->SetLineWidth(2);
  histogram_ref->SetMarkerColor(colors[0]);
  histogram_ref->SetMarkerStyle(markerStyles[0]);
  histogram_ref->SetMarkerSize(markerSizes[0]);
  histogram_ref->Draw("hist");
  legend->AddEntry(histogram_ref, legendEntry_ref.data(), "l");

  TAxis* xAxis_top = histogram_ref->GetXaxis();
  xAxis_top->SetTitle(xAxisTitle.data());
  xAxis_top->SetTitleOffset(xAxisOffset);
  xAxis_top->SetLabelColor(10);
  xAxis_top->SetTitleColor(10);

  TAxis* yAxis_top = histogram_ref->GetYaxis();
  yAxis_top->SetTitle(yAxisTitle.data());
  yAxis_top->SetTitleOffset(yAxisOffset);
  yAxis_top->SetTitleSize(0.065);

  TGraphAsymmErrors* graph1 = nullptr;
  if ( histogram1 ) {
    histogram1->SetLineColor(colors[1]);
    histogram1->SetLineWidth(1);
    histogram1->SetMarkerColor(colors[1]);
    histogram1->SetMarkerStyle(markerStyles[1]);
    histogram1->SetMarkerSize(markerSizes[1]);
    //histogram1->Draw("e1psame");
    graph1 = convertToGraph(histogram1, -0.25);
    graph1->Draw("P");
    legend->AddEntry(histogram1, legendEntry1.data(), "p");
  }
  
  TGraphAsymmErrors* graph2 = nullptr;
  if ( histogram2 ) {
    histogram2->SetLineColor(colors[2]);
    histogram2->SetLineWidth(1);
    histogram2->SetMarkerColor(colors[2]);
    histogram2->SetMarkerStyle(markerStyles[2]);
    histogram2->SetMarkerSize(markerSizes[2]);
    //histogram2->Draw("e1psame");
    graph2 = convertToGraph(histogram2, 0.);
    graph2->Draw("P");
    legend->AddEntry(histogram2, legendEntry2.data(), "p");
  }

  TGraphAsymmErrors* graph3 = nullptr;
  if ( histogram3 ) {
    histogram3->SetLineColor(colors[3]);
    histogram3->SetLineWidth(1);
    histogram3->SetMarkerColor(colors[3]);
    histogram3->SetMarkerStyle(markerStyles[3]);
    histogram3->SetMarkerSize(markerSizes[3]);
    //histogram3->Draw("e1psame");
    graph3 = convertToGraph(histogram3, +0.25);
    graph3->Draw("P");
    legend->AddEntry(histogram3, legendEntry3.data(), "p");
  }

  if ( function ) {
    function->SetLineColor(colors[4]);
    function->SetLineWidth(2);
    function->Draw("same");
  }

  legend->Draw();

  canvas->cd();
  bottomPad->Draw();
  bottomPad->cd();

  std::string histogramName_ratio_ref = std::string(histogram_ref->GetName()).append("_div_").append(histogram_ref->GetName());
  TH1* histogram_ratio_ref = compRatioHistogram(histogramName_ratio_ref, histogram_ref, histogram_ref);
  histogram_ratio_ref->SetTitle("");
  histogram_ratio_ref->SetStats(false);
  histogram_ratio_ref->SetMinimum(-0.50);
  histogram_ratio_ref->SetMaximum(+0.50);
  
  TAxis* xAxis_bottom = histogram_ratio_ref->GetXaxis();
  xAxis_bottom->SetTitle(xAxis_top->GetTitle());
  xAxis_bottom->SetLabelColor(1);
  xAxis_bottom->SetTitleColor(1);
  xAxis_bottom->SetTitleOffset(xAxisOffset);
  xAxis_bottom->SetTitleSize(0.08);
  xAxis_bottom->SetLabelOffset(0.02);
  xAxis_bottom->SetLabelSize(0.08);
  xAxis_bottom->SetTickLength(0.055);
      
  TAxis* yAxis_bottom = histogram_ratio_ref->GetYaxis();
  yAxis_bottom->SetTitle("Relative difference");
  yAxis_bottom->SetTitleOffset(0.85);
  yAxis_bottom->SetNdivisions(505);
  yAxis_bottom->CenterTitle();
  yAxis_bottom->SetTitleSize(0.08);
  yAxis_bottom->SetLabelSize(0.08);
  yAxis_bottom->SetTickLength(0.04);  
  
  histogram_ratio_ref->Draw("axis");

  TGraph* graph_line = new TGraph(2);
  graph_line->SetPoint(0, xAxis_top->GetXmin(), 0.);
  graph_line->SetPoint(1, xAxis_top->GetXmax(), 0.);
  graph_line->SetLineColor(8);
  graph_line->SetLineWidth(1);
  graph_line->Draw("L");

  //TH1* histogram_ratio1 = nullptr;
  //TGraphAsymmErrors* graph_ratio1 = nullptr; 
  //if ( histogram1 ) {
  //  std::string histogramName_ratio1 = std::string(histogram1->GetName()).append("_div_").append(histogram_ref->GetName());
  //  histogram_ratio1 = compRatioHistogram(histogramName_ratio1, histogram1, histogram_ref);
  //  //histogram_ratio1->Draw("e1psame");
  //  graph_ratio1 = convertToGraph(histogram_ratio1, -0.25);
  //  graph_ratio1->Draw("P");
  //}
  //
  //TH1* histogram_ratio2 = nullptr;
  //TGraphAsymmErrors* graph_ratio2 = nullptr;
  //if ( histogram2 ) {
  //  std::string histogramName_ratio2 = std::string(histogram2->GetName()).append("_div_").append(histogram_ref->GetName());
  //  histogram_ratio2 = compRatioHistogram(histogramName_ratio2, histogram2, histogram_ref);
  //  //histogram_ratio2->Draw("e1psame");
  //  graph_ratio2 = convertToGraph(histogram_ratio2, 0.);
  //  graph_ratio2->Draw("P");
  //}

  TH1* histogram_ratio3 = nullptr;
  TGraphAsymmErrors* graph_ratio3 = nullptr;
  if ( histogram3 ) {
    std::string histogramName_ratio3 = std::string(histogram3->GetName()).append("_div_").append(histogram_ref->GetName());
    histogram_ratio3 = compRatioHistogram(histogramName_ratio3, histogram3, histogram_ref);
    //histogram_ratio3->Draw("e1psame");
    graph_ratio3 = convertToGraph(histogram_ratio3, +0.25);
    graph_ratio3->Draw("P");
  }

  histogram_ratio_ref->Draw("axissame");
  
  canvas->Update();
  size_t idx = outputFileName.find_last_of('.');
  std::string outputFileName_plot = std::string(outputFileName, 0, idx);
  if ( useLogScale ) outputFileName_plot.append("_log");
  else outputFileName_plot.append("_linear");
  if ( idx != std::string::npos ) canvas->Print(std::string(outputFileName_plot).append(std::string(outputFileName, idx)).data());
  canvas->Print(std::string(outputFileName_plot).append(".png").data());
  //canvas->Print(std::string(outputFileName_plot).append(".pdf").data());
  //canvas->Print(std::string(outputFileName_plot).append(".root").data());
  
  delete legend;
  delete graph1;
  delete graph2;
  delete graph3;
  delete histogram_ratio_ref;
  delete graph_line;
  //delete histogram_ratio1;
  //delete graph_ratio1;
  //delete histogram_ratio2;
  //delete graph_ratio2;
  delete histogram_ratio3;
  delete graph_ratio3;
  delete topPad;
  delete bottomPad;
  delete canvas;  
}
//---------------------------------------------------------------------------------------------------

using namespace RooFit;

void toyMC()
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

  TRandom3 rnd;

  // generate (pseudo)data
  TH1* histogram_mass = bookHistogram("histogram_mass", "mass", numBins, xMin, xMax);
  fillHistogram(histogram_mass, rnd, xMin, xMax, numEvents_signal, mZ, width, numEvents_background, lambda);
  //dumpHistogram(histogram_mass);

  TH1* histogram_gauss = bookHistogram("histogram_mass_exact", "mass, expectation for S", numBins, xMin, xMax);
  fillHistogram_gauss(histogram_gauss, rnd, numEvents_signal, xMin, xMax, mZ, width);

  TH1* histogram_exp = bookHistogram("histogram_mass_exact", "mass, expectation for B", numBins, xMin, xMax);
  fillHistogram_exp(histogram_exp, rnd, numEvents_background, xMin, xMax, lambda);

  TH1* histogram_sum = bookHistogram("histogram_sum", "mass, expectation for S+B", numBins, xMin, xMax);
  histogram_sum->Add(histogram_gauss);
  histogram_sum->Add(histogram_exp);

  TF1* function_sum = new TF1("function_sum", sum_fcn, xMin, xMax, 6);
  function_sum->SetParameter(0, numEvents_signal/(numEvents_signal + numEvents_background)*((xMax - xMin)/numBins));
  function_sum->SetParameter(1, mZ);
  function_sum->SetParameter(2, width);
  function_sum->SetParameter(3, numEvents_background/(numEvents_signal + numEvents_background)*((xMax - xMin)/numBins));
  function_sum->SetParameter(4, lambda);
  function_sum->SetParameter(5, xMin);

  // compare generated (pseudo)data to expected mass distribution for signal + background
  showHistograms(800, 900,
		 histogram_sum, "S+B Expectation",
		 histogram_gauss, "S Expectation", 
		 histogram_exp, "B Expectation", 
		 histogram_mass, "(Pseudo)data", 
		 function_sum, "Gauss + Exp",
		 "Mass [GeV]", 1.30,
		 true, 1.e-6, 1.9e+1, "Toys", 1.15,
		 0.51, 0.74,
		 "toyMC_random_sampling.png");

  RooRealVar observable("observable", "Mass [GeV]", xMin, xMax);
  RooDataHist data("data", "(Pseudo)data", RooArgList(observable), histogram_mass);

  // define fit model
  RooRealVar mean_gauss("mean_gauss", "Gaussian mean", mZ, 0.8*mZ, 1.2*mZ);
  RooRealVar sigma_gauss("sigma_gauss", "Gaussian resolution", width, 0., 2.*width);
  RooGaussian model_S("model_S", "Signal model",observable, mean_gauss, sigma_gauss);

  RooRealVar lambda_exp("lambda_exp", "Background shape", -lambda, -10.*lambda, 0.);
  RooExponential model_B("pdf_B", "Background model", observable, lambda_exp);

  RooRealVar norm_S("norm_S", "Signal yield", numEvents_signal, 0., 2.*histogram_mass->Integral());
  RooRealVar norm_B("norm_B", "Background yield", numEvents_background, 0., 2.*histogram_mass->Integral());
  RooAddPdf model_SplusB("model_SplusB", "Signal+background model", RooArgList(model_S, model_B), RooArgList(norm_S, norm_B));
    
  // perform fit
  model_SplusB.fitTo(data, PrintLevel(-1));

  // print signal and background yields determined by fit
  //std::cout << "fit results:" << std::endl;
  //std::cout << " S = " << norm_S.getVal();
  double pull_S;
  if ( norm_S.hasAsymError() ) {
    //std::cout << " + " << norm_S.getErrorHi() << " - " << norm_S.getErrorLo();
    if ( norm_S.getVal() > numEvents_signal ) {
      pull_S = (norm_S.getVal() - numEvents_signal)/norm_S.getErrorHi();
    } else {
      pull_S = (norm_S.getVal() - numEvents_signal)/norm_S.getErrorLo();
    }
  } else {
    //std::cout << " +/- " << norm_S.getError();
    pull_S = (norm_S.getVal() - numEvents_signal)/norm_S.getError();
  }
  //std::cout << std::endl;
  //std::cout << " B = " << norm_B.getVal();
  double pull_B;
  if ( norm_B.hasAsymError() ) {
    //std::cout << " + " << norm_B.getErrorHi() << " - " << norm_B.getErrorLo();
    if ( norm_B.getVal() > numEvents_background ) {
      pull_B = (norm_B.getVal() - numEvents_background)/norm_B.getErrorHi();
    } else {
      pull_B = (norm_B.getVal() - numEvents_background)/norm_B.getErrorLo();
    }
  } else {
    //std::cout << " +/- " << norm_B.getError();
    pull_B = (norm_B.getVal() - numEvents_background)/norm_B.getError();
  }
  //std::cout << std::endl;
  
  // make control plot of (pseudo)data versus fit model,
  // for the values of signal and background yields determined by fit
  RooPlot* frame1 = observable.frame();
  data.plotOn(frame1);
  model_SplusB.plotOn(frame1);
  model_SplusB.plotOn(frame1, Components(model_B), LineStyle(ELineStyle::kDashed));
  TCanvas canvas1;
  canvas1.cd();
  frame1->Draw();      
  canvas1.SaveAs("toyMC_fit.png");
      
  // make control plot of likelihood function
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
