
#include <TCanvas.h>
#include <TH1.h>
#include <TMath.h>
#include <TRandom.h>
#include <TRandom3.h>
#include <TROOT.h> // gROOT

#include <iostream> // std::cout, std::endl

void showHistogram(double canvasSizeX, double canvasSizeY,
		   TH1* histogram,
		   const std::string& xAxisTitle, double xAxisOffset,
		   bool useLogScale, double yMin, double yMax, const std::string& yAxisTitle, double yAxisOffset,
		   const std::string& outputFileName)
{
  assert(histogram);

  TCanvas* canvas = new TCanvas("canvas", "canvas", canvasSizeX, canvasSizeY);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2);
  canvas->SetLeftMargin(0.14);
  canvas->SetBottomMargin(0.13);
  canvas->SetRightMargin(0.05);
  canvas->SetLogy(useLogScale);

  int colors[6] = { kBlack, kGreen - 6, kBlue - 7,  kMagenta - 7, kCyan - 6, kRed - 6 };
  int markerStyles[6] = { 24, 21, 20, 21, 22, 23 };
  int markerSizes[6] = { 1, 1, 1, 1, 1, 1 };

  histogram->SetTitle("");
  histogram->SetStats(false);
  histogram->SetMinimum(yMin);
  histogram->SetMaximum(yMax);
  histogram->SetLineColor(colors[0]);
  histogram->SetLineWidth(2);
  histogram->SetMarkerColor(colors[0]);
  histogram->SetMarkerStyle(markerStyles[0]);
  histogram->SetMarkerSize(markerSizes[0]);
  histogram->Draw("hist");

  TAxis* xAxis = histogram->GetXaxis();
  xAxis->SetTitle(xAxisTitle.data());
  xAxis->SetTitleOffset(xAxisOffset);
  xAxis->SetTitleSize(0.055);

  TAxis* yAxis = histogram->GetYaxis();
  yAxis->SetTitle(yAxisTitle.data());
  yAxis->SetTitleOffset(yAxisOffset);
  yAxis->SetTitleSize(0.055);
  
  histogram->Draw("axissame");
  
  canvas->Update();
  size_t idx = outputFileName.find_last_of('.');
  std::string outputFileName_plot = std::string(outputFileName, 0, idx);
  if ( useLogScale ) outputFileName_plot.append("_log");
  else outputFileName_plot.append("_linear");
  if ( idx != std::string::npos ) canvas->Print(std::string(outputFileName_plot).append(std::string(outputFileName, idx)).data());
  canvas->Print(std::string(outputFileName_plot).append(".png").data());
  //canvas->Print(std::string(outputFileName_plot).append(".pdf").data());
  //canvas->Print(std::string(outputFileName_plot).append(".root").data());
  
  delete canvas;  
}

void poisson()
{
  // prevent pop-up windows
  gROOT->SetBatch(true);

  // define cross section for SM Higgs boson production via gluon fuion
  // in proton-proton collisions at 13 TeV, computed at 3LO QCD and NLO EW accuracy
  // (taken from https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageAt13TeV )
  double cross_section            = 48.5;    // pb
  double conversion_factor        =  1.e-36; // cm^2 pb-1
  double instantaneous_luminosity =  2.e+34; // cm^-2 s^-1
  double numEvents_expected       = cross_section*conversion_factor*instantaneous_luminosity;
  std::cout << "Average number of Higgs bosons produced per second = " << numEvents_expected << std::endl;

  TH1* histogram = new TH1D("histogram", "Expected number of Higgs bosons", 6, -0.5, 5.5);
  for ( int idxBin = 1; idxBin <= histogram->GetNbinsX(); ++idxBin ) {
    double x = idxBin - 1;
    double y = TMath::PoissonI(x, numEvents_expected);
    histogram->SetBinContent(idxBin, y);
  }
  showHistogram(800, 600, 
		histogram,
		"Number of Higgs bosons", 1.05, 
		true, 1.e-3, 1.e0, "Probability", 1.1, 
		"poisson.png");
}
