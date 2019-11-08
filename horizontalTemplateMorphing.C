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
#include <algorithm> // std::sort

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
		    TH1* histogram1, const std::string& legendEntry1,
		    TH1* histogram2, const std::string& legendEntry2,
		    TH1* histogram3, const std::string& legendEntry3,
		    TH1* histogram4, const std::string& legendEntry4,
		    int colors[], int markerStyles[], int markerSizes[],
		    const std::string& xAxisTitle, double xAxisOffset,
		    bool useLogScale, double yMin, double yMax, const std::string& yAxisTitle, double yAxisOffset,
		    double legendX0, double legendY0, double legendSizeX, double legendSizeY, 
		    const std::string& outputFileName)
{
  assert(histogram_ref);

  TCanvas* canvas = new TCanvas("canvas", "canvas", canvasSizeX, canvasSizeY);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2);
  canvas->SetTopMargin(0.04);
  canvas->SetLeftMargin(0.12);
  canvas->SetBottomMargin(0.15);
  canvas->SetRightMargin(0.04);
  canvas->SetLogy(useLogScale);

  TLegend* legend = new TLegend(legendX0, legendY0, legendX0 + 0.18, legendY0 + 0.20, "", "brNDC"); 
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetTextSize(0.045);

  histogram1->SetTitle("");
  histogram1->SetStats(false);
  histogram1->SetMinimum(yMin);
  histogram1->SetMaximum(yMax);
  histogram1->SetLineColor(colors[0]);
  histogram1->SetLineWidth(2);
  histogram1->SetMarkerColor(colors[0]);
  histogram1->SetMarkerStyle(markerStyles[0]);
  histogram1->SetMarkerSize(markerSizes[0]);
  histogram1->Draw("hist");
  legend->AddEntry(histogram1, legendEntry1.data(), "l");

  TAxis* xAxis = histogram1->GetXaxis();
  xAxis->SetTitle(xAxisTitle.data());
  xAxis->SetTitleOffset(xAxisOffset);
  xAxis->SetTitleSize(0.055);

  TAxis* yAxis = histogram1->GetYaxis();
  yAxis->SetTitle(yAxisTitle.data());
  yAxis->SetTitleOffset(yAxisOffset);
  yAxis->SetTitleSize(0.060);
  yAxis->SetNdivisions(505);

  if ( histogram2 ) {
    histogram2->SetLineColor(colors[1]);
    histogram2->SetLineWidth(2);
    histogram2->SetMarkerColor(colors[1]);
    histogram2->SetMarkerStyle(markerStyles[1]);
    histogram2->SetMarkerSize(markerSizes[1]);
    histogram2->Draw("histsame");
    legend->AddEntry(histogram2, legendEntry2.data(), "l");
  }
  
  if ( histogram3 ) {
    histogram3->SetLineColor(colors[2]);
    histogram3->SetLineWidth(2);
    histogram3->SetMarkerColor(colors[2]);
    histogram3->SetMarkerStyle(markerStyles[2]);
    histogram3->SetMarkerSize(markerSizes[2]);
    histogram3->Draw("histsame");
    legend->AddEntry(histogram3, legendEntry3.data(), "l");
  }

  if ( histogram4 ) {
    histogram4->SetLineColor(colors[3]);
    histogram4->SetLineWidth(2);
    histogram4->SetLineStyle(7);
    histogram4->SetMarkerColor(colors[3]);
    histogram4->SetMarkerStyle(markerStyles[3]);
    histogram4->SetMarkerSize(markerSizes[3]);
    histogram4->Draw("histsame");
    legend->AddEntry(histogram4, legendEntry4.data(), "l");
  }

  histogram1->Draw("axissame");

  legend->Draw();
  
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
  delete canvas;  
}
//---------------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------------
void showGraphs(double canvasSizeX, double canvasSizeY,
		TGraph* graph1, const std::string& legendEntry1,
                TGraph* graph2, const std::string& legendEntry2,
                TGraph* graph3, const std::string& legendEntry3,
                TGraph* graph4, const std::string& legendEntry4,
                int colors[], int markerStyles[], int markerSizes[],
                double xMin, double xMax, const std::string& xAxisTitle, double xAxisOffset,
                bool useLogScale, double yMin, double yMax, const std::string& yAxisTitle, double yAxisOffset,
		double legendX0, double legendY0, double legendSizeX, double legendSizeY, 
                const std::string& outputFileName)
{
  TCanvas* canvas = new TCanvas("canvas", "canvas", canvasSizeX, canvasSizeY);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2);
  canvas->SetTopMargin(0.04);
  canvas->SetLeftMargin(0.12);
  canvas->SetBottomMargin(0.15);
  canvas->SetRightMargin(0.04);
  canvas->SetLogy(useLogScale);

  TLegend* legend = new TLegend(legendX0, legendY0, legendX0 + legendSizeX, legendY0 + legendSizeY, "", "brNDC"); 
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetTextSize(0.045);

  TH1* dummyHistogram = new TH1D("dummyHistogram", "dummyHistogram", 100, xMin, xMax);
  dummyHistogram->SetTitle("");
  dummyHistogram->SetStats(false);
  dummyHistogram->SetMinimum(yMin);
  dummyHistogram->SetMaximum(yMax);

   TAxis* xAxis = dummyHistogram->GetXaxis();
  xAxis->SetTitle(xAxisTitle.data());
  xAxis->SetTitleOffset(xAxisOffset);
  xAxis->SetTitleSize(0.055);

  TAxis* yAxis = dummyHistogram->GetYaxis();
  yAxis->SetTitle(yAxisTitle.data());
  yAxis->SetTitleOffset(yAxisOffset);
  yAxis->SetTitleSize(0.060);

  dummyHistogram->Draw("axis");

  graph1->SetLineColor(colors[0]);
  graph1->SetLineWidth(2);
  graph1->SetMarkerColor(colors[0]);
  graph1->SetMarkerStyle(markerStyles[0]);
  graph1->SetMarkerSize(markerSizes[0]);
  graph1->Draw("L");
  legend->AddEntry(graph1, legendEntry1.data(), "l");

  if ( graph2 ) {
    graph2->SetLineColor(colors[1]);
    graph2->SetLineWidth(2);
    graph2->SetMarkerColor(colors[1]);
    graph2->SetMarkerStyle(markerStyles[1]);
    graph2->SetMarkerSize(markerSizes[1]);
    graph2->Draw("L");
    legend->AddEntry(graph2, legendEntry2.data(), "l");
  }
  
  if ( graph3 ) {
    graph3->SetLineColor(colors[2]);
    graph3->SetLineWidth(2);
    graph3->SetMarkerColor(colors[2]);
    graph3->SetMarkerStyle(markerStyles[2]);
    graph3->SetMarkerSize(markerSizes[2]);
    graph3->Draw("L");
    legend->AddEntry(graph3, legendEntry3.data(), "l");
  }

  if ( graph4 ) {
    graph4->SetLineColor(colors[3]);
    graph4->SetLineWidth(2);
    graph4->SetMarkerColor(colors[3]);
    graph4->SetMarkerStyle(markerStyles[3]);
    graph4->SetMarkerSize(markerSizes[3]);
    graph4->Draw("L");
    legend->AddEntry(graph4, legendEntry4.data(), "l");
  }
  
  dummyHistogram->Draw("axissame");
   
  legend->Draw();

  canvas->Update();
  size_t idx = outputFileName.find_last_of('.');
  std::string outputFileName_plot = std::string(outputFileName, 0, idx);
  if ( useLogScale ) outputFileName_plot.append("_log");
  else outputFileName_plot.append("_linear");
  if ( idx != std::string::npos ) canvas->Print(std::string(outputFileName_plot).append(std::string(outputFileName, idx)).data());
  canvas->Print(std::string(outputFileName_plot).append(".png").data());
  //canvas->Print(std::string(outputFileName_plot).append(".pdf").data());
  //canvas->Print(std::string(outputFileName_plot).append(".root").data());
  
  delete dummyHistogram;
  delete legend;
  delete canvas;  
}
//---------------------------------------------------------------------------------------------------

double square(double x)
{
  return x*x;
}

TGraphAsymmErrors* makeGraphCDF(const TH1* histogram)
{
  // prepare graph of CDF (cumulative distribution function)
  const TAxis* xAxis = histogram->GetXaxis();
  int numBins = xAxis->GetNbins();
  const int numSteps_per_bin = 20;
  TGraphAsymmErrors* graph = new TGraphAsymmErrors(numBins*numSteps_per_bin );
  double sum = 0.;
  double sumErr2 = 0.;
  int idxPoint = 0;
  for ( int idxBin = 1; idxBin <= numBins; ++idxBin ) {
    double binEdgeLo = xAxis->GetBinLowEdge(idxBin);
    double binEdgeHi = xAxis->GetBinUpEdge(idxBin);
    double binContent = histogram->GetBinContent(idxBin);
    double binError = histogram->GetBinError(idxBin);
    const int numSteps = 1;
    for ( int idxStep = 0; idxStep < numSteps_per_bin; ++idxStep ) {
      sum += binContent/numSteps_per_bin;
      sumErr2 += square(binError/numSteps_per_bin);
      double x = binEdgeLo + idxStep*(binEdgeHi - binEdgeLo)/numSteps_per_bin;
      graph->SetPoint(idxPoint, x, sum);
      double binEdgeLo = xAxis->GetBinLowEdge(idxBin);
      double binEdgeHi = xAxis->GetBinUpEdge(idxBin);
      double sumErr = TMath::Sqrt(sumErr2);
      double xErr = 0.5*(binEdgeHi - binEdgeLo)/numSteps_per_bin;
      graph->SetPointError(idxPoint, xErr, xErr, sumErr, sumErr);
      ++idxPoint;
    }
  }
  return graph;
}

TGraphAsymmErrors* makeInverseGraph(const TGraphAsymmErrors* graph)
{
  // prepare graph of inverse function (by swapping x- and y-axis)
  int numPoints = graph->GetN();
  TGraphAsymmErrors* graph_inverse = new TGraphAsymmErrors(numPoints);
  for ( int idxPoint = 0; idxPoint < numPoints; ++idxPoint ) {
    double x, y;
    graph->GetPoint(idxPoint, x, y);
    double xErrLo = graph->GetErrorXlow(idxPoint);
    double xErrHi = graph->GetErrorXhigh(idxPoint);
    double yErrLo = graph->GetErrorYlow(idxPoint);
    double yErrHi = graph->GetErrorYhigh(idxPoint);
    graph_inverse->SetPoint(idxPoint, y, x);
    graph_inverse->SetPointError(idxPoint, yErrLo, yErrHi, xErrLo, xErrHi);
  }
  return graph_inverse;
}

std::vector<TGraph*>
interpolateHistogram(TH1* histogram_interpolated, const TH1* histogram_nominal, const TH1* histogram_endpoint, double a)
{
  assert(a >= 0. && a <= 1.);
  double b = 1. - a;
  TGraphAsymmErrors* graphCDF_nominal = makeGraphCDF(histogram_nominal);
  TGraphAsymmErrors* graphCDFinverse_nominal = makeInverseGraph(graphCDF_nominal);
  TGraphAsymmErrors* graphCDF_endpoint = makeGraphCDF(histogram_endpoint);
  TGraphAsymmErrors* graphCDFinverse_endpoint = makeInverseGraph(graphCDF_endpoint);
  std::vector<double> points_y;
  int numPoints_nominal = graphCDF_nominal->GetN();
  for ( int idxPoint = 0; idxPoint < numPoints_nominal; ++idxPoint ) {
    double x, y;
    graphCDF_nominal->GetPoint(idxPoint, x, y);
    points_y.push_back(y);
  }
  int numPoints_endpoint = graphCDF_endpoint->GetN();
  for ( int idxPoint = 0; idxPoint < numPoints_endpoint; ++idxPoint ) {
    double x, y;
    graphCDF_endpoint->GetPoint(idxPoint, x, y);
    points_y.push_back(y);
  }
  std::sort(points_y.begin(), points_y.end()); // sort y-values in ascending order
  std::vector<double> points_y_woDuplicates;
  double y_last = -1.;
  for ( std::vector<double>::const_iterator y = points_y.begin(); y != points_y.end(); ++y ) {
    const double epsilon = 1.e-3;
    if ( y != points_y.begin() && TMath::Abs((*y) - y_last) < epsilon ) continue;
    points_y_woDuplicates.push_back(*y);
    y_last = (*y);
  }
  int numPoints_interpolated = points_y_woDuplicates.size();
  TGraph* graphCDF_interpolated = new TGraph(numPoints_interpolated);
  for ( int idxPoint = 0; idxPoint < numPoints_interpolated; ++idxPoint ) {
    double y = points_y_woDuplicates[idxPoint];
    double x1 = graphCDFinverse_nominal->Eval(y);
    double x2 = graphCDFinverse_endpoint->Eval(y);
    double x = a*x1 + b*x2;
    graphCDF_interpolated->SetPoint(idxPoint, x, y);
  }
  const TAxis* xAxis = histogram_interpolated->GetXaxis();
  int numBins = xAxis->GetNbins();
  for ( int idxBin = 1; idxBin <= numBins; ++idxBin ) {
    double binEdgeLo = xAxis->GetBinLowEdge(idxBin);
    double binEdgeHi = xAxis->GetBinUpEdge(idxBin);
    double binContent = graphCDF_interpolated->Eval(binEdgeHi) - graphCDF_interpolated->Eval(binEdgeLo);
    assert(binContent >= 0.);
    histogram_interpolated->SetBinContent(idxBin, binContent);
    histogram_interpolated->SetBinError(idxBin, 0.);
  }
  std::vector<TGraph*> graphsCDF;
  graphsCDF.push_back(graphCDF_nominal);
  graphsCDF.push_back(graphCDF_endpoint);
  graphsCDF.push_back(graphCDF_interpolated);  
  return graphsCDF;
}

void horizontalTemplateMorphing()
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

  TH1* histogram_nominal = bookHistogram("histogram_nominal", "Nominal", numBins, xMin, xMax);
  fillHistogram_Erf(histogram_nominal, mZ, width);

  TH1* histogram_endpoint = bookHistogram("histogram_endpoint", "Endpoint", numBins, xMin, xMax);
  fillHistogram_Erf(histogram_endpoint, mZ - 20., width);
  
  TH1* histogram_interpolated = bookHistogram("histogram_interpolated", "Interpolated", numBins, xMin, xMax);
  std::vector<TGraph*> graphsCDF = interpolateHistogram(histogram_interpolated, histogram_nominal, histogram_endpoint, 0.5);
  
  TGraph* graphCDF_nominal = graphsCDF[0];
  TGraph* graphCDF_endpoint = graphsCDF[1];
  TGraph* graphCDF_interpolated = graphsCDF[2];

  int colors[6] = { kBlack, kGreen - 6, kBlue - 7,  kMagenta - 7, kCyan - 6, kRed - 6 };
  int markerStyles[6] = { 24, 21, 20, 21, 22, 23 };
  int markerSizes[6] = { 1, 1, 1, 1, 1, 1 };

  showHistograms(800, 600,
		 histogram_nominal,  "m_{Z} =  91.2 GeV",
		 histogram_endpoint, "m_{Z} =  71.2 GeV",
		 nullptr, "",
		 nullptr, "",
		 colors, markerStyles, markerSizes,
		 "Mass [GeV]", 1.15,
		 false, 0., 0.24, "Events", 1.05,
		 0.68, 0.74, 0.21, 0.12,
		 "horizontalTemplateMorphing_woInterpolatedHistogram.png");

  showGraphs(800, 600,
	     graphCDF_nominal,      "m_{Z} =  91.2 GeV",
	     graphCDF_endpoint,     "m_{Z} =  71.2 GeV",
	     graphCDF_interpolated, "m_{Z} =  81.2 GeV",
	     nullptr, "",
	     colors, markerStyles, markerSizes,
	     xMin, xMax, "Mass [GeV]", 1.15,
	     false, 0., 1.09, "#Sigma Events", 1.05,
	     0.66, 0.24, 0.21, 0.20,
	     "horizontalTemplateMorphing_cdf.png");

  showHistograms(800, 600,
		 histogram_nominal,      "m_{Z} =  91.2 GeV",
		 histogram_endpoint,     "m_{Z} =  71.2 GeV",
		 histogram_interpolated, "m_{Z} =  81.2 GeV",
		 nullptr, "",
		 colors, markerStyles, markerSizes,
		 "Mass [GeV]", 1.15,
		 false, 0., 0.24, "Events", 1.05,
		 0.68, 0.73, 0.21, 0.20,
		 "horizontalTemplateMorphing_wInterpolatedHistogram.png");

  delete histogram_nominal;
  delete histogram_endpoint;
  delete histogram_interpolated;
  delete graphCDF_nominal;
  delete graphCDF_endpoint;
  delete graphCDF_interpolated;
}
