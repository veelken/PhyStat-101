
#include <TCanvas.h>
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

double normal_cdf(double x, double mZ, double sigma)
{
  double t = TMath::Sqrt(0.5)*(x - mZ)/sigma;
  return 0.5*TMath::Erf(t);
}

void fillHistogram_Erf(TH1* histogram, double mZ, double sigma)
{
  const TAxis* xAxis = histogram->GetXaxis();
  int numBins = xAxis->GetNbins();
  for ( int idxBin = 1; idxBin <= numBins; ++idxBin ) {
    double binEdge_hi = xAxis->GetBinUpEdge(idxBin);
    double binEdge_lo = xAxis->GetBinLowEdge(idxBin);
    double integral = normal_cdf(binEdge_hi, mZ, sigma) - normal_cdf(binEdge_lo, mZ, sigma);
    histogram->SetBinContent(idxBin, integral);
    histogram->SetBinError(idxBin, 0.);
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

  TLegend* legend = new TLegend(legendX0, legendY0, legendX0 + 0.18, legendY0 + 0.20, "", "brNDC"); 
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
    histogram1->SetLineWidth(2);
    histogram1->SetMarkerColor(colors[1]);
    histogram1->SetMarkerStyle(markerStyles[1]);
    histogram1->SetMarkerSize(markerSizes[1]);
    histogram1->Draw("histsame");
    legend->AddEntry(histogram1, legendEntry1.data(), "l");
    //graph1 = convertToGraph(histogram1, -0.25);
    //graph1->Draw("P");
    //legend->AddEntry(histogram1, legendEntry1.data(), "p");
  }
  
  TGraphAsymmErrors* graph2 = nullptr;
  if ( histogram2 ) {
    histogram2->SetLineColor(colors[2]);
    histogram2->SetLineWidth(2);
    histogram2->SetMarkerColor(colors[2]);
    histogram2->SetMarkerStyle(markerStyles[2]);
    histogram2->SetMarkerSize(markerSizes[2]);
    histogram2->Draw("histsame");
    legend->AddEntry(histogram2, legendEntry2.data(), "l");
    //graph2 = convertToGraph(histogram2, 0.);
    //graph2->Draw("P");
    //legend->AddEntry(histogram2, legendEntry2.data(), "p");
  }

  TGraphAsymmErrors* graph3 = nullptr;
  if ( histogram3 ) {
    histogram3->SetLineColor(colors[3]);
    histogram3->SetLineWidth(2);
    histogram3->SetLineStyle(7);
    histogram3->SetMarkerColor(colors[3]);
    histogram3->SetMarkerStyle(markerStyles[3]);
    histogram3->SetMarkerSize(markerSizes[3]);
    histogram3->Draw("histsame");
    legend->AddEntry(histogram3, legendEntry3.data(), "l");
    //graph3 = convertToGraph(histogram3, +0.25);
    //graph3->Draw("P");
    //legend->AddEntry(histogram3, legendEntry3.data(), "p");
  }

  legend->Draw();

  canvas->cd();
  bottomPad->Draw();
  bottomPad->cd();

  std::string histogramName_ratio_ref = std::string(histogram_ref->GetName()).append("_div_").append(histogram_ref->GetName());
  TH1* histogram_ratio_ref = compRatioHistogram(histogramName_ratio_ref, histogram_ref, histogram_ref);
  histogram_ratio_ref->SetTitle("");
  histogram_ratio_ref->SetStats(false);
  histogram_ratio_ref->SetMinimum(-1.00);
  histogram_ratio_ref->SetMaximum(+1.00);
  
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

  TH1* histogram_ratio1 = nullptr;
  TGraphAsymmErrors* graph_ratio1 = nullptr; 
  if ( histogram1 ) {
    std::string histogramName_ratio1 = std::string(histogram1->GetName()).append("_div_").append(histogram_ref->GetName());
    histogram_ratio1 = compRatioHistogram(histogramName_ratio1, histogram1, histogram_ref);
    histogram_ratio1->Draw("histsame");
    //graph_ratio1 = convertToGraph(histogram_ratio1, -0.25);
    //graph_ratio1->Draw("P");
  }

  TH1* histogram_ratio2 = nullptr;
  TGraphAsymmErrors* graph_ratio2 = nullptr;
  if ( histogram2 ) {
    std::string histogramName_ratio2 = std::string(histogram2->GetName()).append("_div_").append(histogram_ref->GetName());
    histogram_ratio2 = compRatioHistogram(histogramName_ratio2, histogram2, histogram_ref);
    histogram_ratio2->Draw("histsame");
    //graph_ratio2 = convertToGraph(histogram_ratio2, 0.);
    //graph_ratio2->Draw("P");
  }

  TH1* histogram_ratio3 = nullptr;
  TGraphAsymmErrors* graph_ratio3 = nullptr;
  if ( histogram3 ) {
    std::string histogramName_ratio3 = std::string(histogram3->GetName()).append("_div_").append(histogram_ref->GetName());
    histogram_ratio3 = compRatioHistogram(histogramName_ratio3, histogram3, histogram_ref);
    histogram_ratio3->Draw("histsame");
    //graph_ratio3 = convertToGraph(histogram_ratio3, +0.25);
    //graph_ratio3->Draw("P");
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
  delete histogram_ratio1;
  delete graph_ratio1;
  delete histogram_ratio2;
  delete graph_ratio2;
  delete histogram_ratio3;
  delete graph_ratio3;
  delete topPad;
  delete bottomPad;
  delete canvas;  
}
//---------------------------------------------------------------------------------------------------

void interpolateHistogram(TH1* histogram_alpha, const TH1* histogram_plus1p0sigma, const TH1* histogram_nominal, const TH1* histogram_minus1p0sigma, double lambda)
{
  int numBins = histogram_nominal->GetNbinsX();
  assert(histogram_alpha->GetNbinsX() == numBinsX);
  assert(histogram_plus1p0sigma->GetNbinsX() == numBinsX);
  assert(histogram_minus1p0sigma->GetNbinsX() == numBinsX);  
  for ( int idxBin = 1; idxBin <= numBins; ++idxBin ) {
    double binContent_plus1p0sigma = histogram_plus1p0sigma->GetBinContent(idxBin);
    double binError_plus1p0sigma = histogram_plus1p0sigma->GetBinError(idxBin);
    double binContent_nominal = histogram_nominal->GetBinContent(idxBin);
    double binError_nominal = histogram_nominal->GetBinError(idxBin);
    double binContent_minus1p0sigma = histogram_minus1p0sigma->GetBinContent(idxBin);
    double binError_minus1p0sigma = histogram_minus1p0sigma->GetBinError(idxBin);
    double a = 0.;
    if      ( TMath::Abs(lambda) <=  1. ) a = 0.5*lambda*(lambda + 1.);
    else if (            lambda  >  +1. ) a = lambda;
    else if (            lambda  <  -1. ) a = 0.;
    else assert(0);
    double b = 0.;
    if      ( TMath::Abs(lambda) <=  1. ) b = -(lambda - 1.)*(lambda + 1.);
    else                                  b = -(TMath::Abs(lambda) - 1.);
    double c = 0.;
    if      ( TMath::Abs(lambda) <=  1. ) c = 0.5*lambda*(lambda - 1.);
    else if (            lambda  <  -1. ) c = TMath::Abs(lambda);
    else if (            lambda  >  +1. ) c = 0.;
    else assert(0);
    double binContent_lambda = a*binContent_plus1p0sigma + b*binContent_nominal + c*binContent_minus1p0sigma;
    double binError_lambda = a*binError_plus1p0sigma + b*binError_nominal + c*binError_minus1p0sigma;
    histogram_alpha->SetBinContent(idxBin, binContent_lambda);
    histogram_alpha->SetBinError(idxBin, binError_lambda);
  }
}

void verticalTemplateMorphing()
{
  // prevent pop-up windows
  gROOT->SetBatch(true);

  // define binning
  // (30 bins from 60 to 120 GeV)
  int numBins = 30;
  double xMin =  60.; // GeV
  double xMax = 120.; // GeV

  // define mass of Z boson
  // (taken from http://pdg.lbl.gov/2018/listings/rpp2018-list-z-boson.pdf )
  double mZ       = 91.2; // GeV
  double width    =  4.5; // GeV (includes intrinsic width of Z boson + experimental resolution)

  double sigma_mZ = 1.0; // GeV

  TH1* histogram_nominal = bookHistogram("histogram_nominal", "Nominal", numBins, xMin, xMax);
  fillHistogram_Erf(histogram_nominal, mZ, width);

  TH1* histogram_plus1p0sigma = bookHistogram("histogram_plus1p0sigma", "+1 Sigma", numBins, xMin, xMax);
  fillHistogram_Erf(histogram_plus1p0sigma, mZ + sigma_mZ, width);

  TH1* histogram_minus1p0sigma = bookHistogram("histogram_minus1p0sigma", "-1 Sigma", numBins, xMin, xMax);
  fillHistogram_Erf(histogram_minus1p0sigma, mZ - sigma_mZ, width);
  
  TH1* histogram_plus0p5sigma = bookHistogram("histogram_plus0p5sigma", "+1 Sigma", numBins, xMin, xMax);
  interpolateHistogram(histogram_plus0p5sigma, histogram_plus1p0sigma, histogram_nominal, histogram_minus1p0sigma , 0.5);

  showHistograms(800, 900,
		 histogram_nominal,       "Nominal",
		 histogram_plus1p0sigma,  "+1 #sigma", 
		 histogram_minus1p0sigma, "-1 #sigma", 
		 histogram_plus0p5sigma,  "+0.5 #sigma",
		 "Mass [GeV]", 1.30,
		 true, 1.e-6, 1.9e+1, "Events", 1.15,
		 0.76, 0.74,
		 "verticalTemplateMorphing.png");

  delete histogram_nominal;
  delete histogram_plus1p0sigma;
  delete histogram_minus1p0sigma;
  delete histogram_plus0p5sigma;
}
