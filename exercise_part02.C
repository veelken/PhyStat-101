
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
// TMath::Gaus function does not include the 1/sqrt(2*pi)*sigma term,
// which is needed for proper normalization required for a probability density function:
//  integral_{-infinity}^{+infinity} Gaus(x, mZ, sigma) dx = 1,
// so we define the functions with proper normalization here

double normal(double x, double mZ, double sigma)
{
  return (1./(TMath::Sqrt(2.*TMath::Pi())*sigma))*TMath::Gaus(x, mZ, sigma);
}

double normal_fcn(double* x, double* par)
{
  return par[0]*normal(x[0], par[1], par[2]);
}
//---------------------------------------------------------------------------------------------------

void fillHistogram_rejection_sampling(TH1* histogram, TRandom& rnd, int numToys, double xMin, double xMax, double mZ, double sigma)
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
  if ( histogram->Integral() > 0. ) {
    histogram->Scale(1./histogram->Integral());
  }
}

void fillHistogram_inverse_transform_sampling(TH1* histogram, TRandom& rnd, int numToys, double xMin, double xMax, double mZ, double sigma)
{
  // prepare graph of CDF^-1 (inverse of cumulative distribution function)
  std::vector<double> points_x;
  std::vector<double> points_y;
  const double x_step = 1.e-3;
  double sum = 0.;
  for ( double x = xMin; x <= xMax; x += x_step ) {
    points_x.push_back(x);
    points_y.push_back(sum);
    sum += x_step*normal(x, mZ, sigma);
  }
  assert(points_x.size() == points_y.size());
  int numPoints = points_x.size();
  TGraph* graph_inverse_cdf = new TGraph(numPoints);
  for ( int idxPoint = 0; idxPoint < numPoints; ++idxPoint ) {
    // graph of *inverse* of CDF is obtained by swapping x and y
    graph_inverse_cdf->SetPoint(idxPoint, points_y[idxPoint], points_x[idxPoint]);
  }
  for ( int idxToy = 0; idxToy < numToys; ++idxToy ) {
    double y = rnd.Rndm();
    double x = graph_inverse_cdf->Eval(y);
    histogram->Fill(x);
  }
  if ( histogram->Integral() > 0. ) {
    histogram->Scale(1./histogram->Integral());
  }
  delete graph_inverse_cdf;
}

void fillHistogram_rndSum(TH1* histogram, TRandom& rnd, int numToys, double xMin, double xMax, double mZ, double sigma)
{
  const int N = 100;
  double range = TMath::Sqrt(12.*N)*sigma;
  double a = mZ - 0.5*range;
  double b = mZ + 0.5*range;
  for ( int idxToy = 0; idxToy < numToys; ++idxToy ) {
    double sum = 0.;
    for ( int i = 0; i < N; ++i ) {
      sum += a + (b - a)*rnd.Rndm();
    }
    double x = sum/N;
    histogram->Fill(x);
  }
  if ( histogram->Integral() > 0. ) {
    histogram->Scale(1./histogram->Integral());
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

  TH1* histogram_ratio1 = nullptr;
  TGraphAsymmErrors* graph_ratio1 = nullptr; 
  if ( histogram1 ) {
    std::string histogramName_ratio1 = std::string(histogram1->GetName()).append("_div_").append(histogram_ref->GetName());
    histogram_ratio1 = compRatioHistogram(histogramName_ratio1, histogram1, histogram_ref);
    //histogram_ratio1->Draw("e1psame");
    graph_ratio1 = convertToGraph(histogram_ratio1, -0.25);
    graph_ratio1->Draw("P");
  }

  TH1* histogram_ratio2 = nullptr;
  TGraphAsymmErrors* graph_ratio2 = nullptr;
  if ( histogram2 ) {
    std::string histogramName_ratio2 = std::string(histogram2->GetName()).append("_div_").append(histogram_ref->GetName());
    histogram_ratio2 = compRatioHistogram(histogramName_ratio2, histogram2, histogram_ref);
    //histogram_ratio2->Draw("e1psame");
    graph_ratio2 = convertToGraph(histogram_ratio2, 0.);
    graph_ratio2->Draw("P");
  }

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

void exercise_part02()
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
  double mZ    = 91.2; // GeV
  double width =  4.5; // GeV (includes intrinsic width of Z boson + experimental resolution)

  TF1* function_Gaussian = new TF1("function_Gaussian", normal_fcn, xMin, xMax, 2);
  function_Gaussian->SetParameter(0, (xMax - xMin)/numBins);
  function_Gaussian->SetParameter(1, mZ);
  function_Gaussian->SetParameter(2, width);

  TH1* histogram_Erf = bookHistogram("histogram_Erf", "Exact solution", numBins, xMin, xMax);
  fillHistogram_Erf(histogram_Erf, mZ, width);

  TRandom3 rnd;
  int numToys = 100000;

  TH1* histogram_rejection = bookHistogram("histogram_rejection", "Rejection sampling", numBins, xMin, xMax);
  fillHistogram_rejection_sampling(histogram_rejection, rnd, numToys, xMin, xMax, mZ, width);

  TH1* histogram_inverse_transform = bookHistogram("histogram_inverse_transform", "Inverse-transform sampling", numBins, xMin, xMax);
  fillHistogram_inverse_transform_sampling(histogram_inverse_transform, rnd, numToys, xMin, xMax, mZ, width);

  TH1* histogram_rndSum = bookHistogram("histogram_rndSum", "#Sigma u", numBins, xMin, xMax);
  fillHistogram_rndSum(histogram_rndSum, rnd, numToys, xMin, xMax, mZ, width);

  showHistograms(800, 900,
		 histogram_Erf, "Exact solution",
		 histogram_rejection, "Rejection sampling", 
		 histogram_inverse_transform, "Inverse-transform sampling", 
		 histogram_rndSum, "#Sigma u", 
		 function_Gaussian, "Gaussian",
		 "Mass [GeV]", 1.30,
		 true, 1.e-6, 1.9e+1, "Toys", 1.15,
		 0.51, 0.74,
		 "exercise_part02.png");

  delete function_Gaussian;
  delete histogram_Erf;
  delete histogram_rejection;
  delete histogram_inverse_transform;
  delete histogram_rndSum;
}
