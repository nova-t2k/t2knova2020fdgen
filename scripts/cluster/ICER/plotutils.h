#pragma once

#include <iostream>
#include <limits>
#include <memory>
#include <vector>

#include "TAxis.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TGaxis.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TPad.h"
#include "TSpline.h"

template <typename TH>
double IntegralTH(std::unique_ptr<TH> &h) {
  if (!h) {
    return 0;
  }
  return h->Integral();
}

void StyleAxis(TAxis *a, double scale = 1, double titleoffsetscale = 1,
               double labeloffsetscale = 1) {
  a->SetNdivisions(505);
  a->SetLabelSize(0.05 * scale);
  a->SetLabelOffset(0.01 * labeloffsetscale);
  a->SetTitleSize(0.05 * scale);
  a->SetTitleOffset(1.2 * titleoffsetscale);
}

void StyleAxes(std::unique_ptr<TH1> &a, double scale = 1,
               double titleoffsetscale = 1, double labeloffsetscale = 1) {
  StyleAxis(a->GetYaxis(), scale, titleoffsetscale, labeloffsetscale);
  StyleAxis(a->GetXaxis(), scale, titleoffsetscale, labeloffsetscale);
}

void HideAxis(TAxis *a) {
  a->SetLabelSize(0);
  a->SetLabelOffset(10);
  a->SetTitleSize(0);
  a->SetTitleOffset(10);
}

void StyleAxis(TGaxis *g, TAxis *a) {
  g->SetLabelFont(a->GetLabelFont());
  g->SetLabelSize(a->GetLabelSize());
  g->SetLabelOffset(2 * a->GetLabelOffset());
  g->SetTitleFont(a->GetTitleFont());
  g->SetTitleSize(a->GetTitleSize());
  g->SetTitleOffset(a->GetTitleOffset());
}

void StyleTH1Line(std::unique_ptr<TH1> &h, int color, int width = 1,
                  int style = 1, double alpha = 1) {
  if (!h) {
    return;
  }
  h->SetLineStyle(style);
  h->SetLineColorAlpha(color, alpha);
  h->SetLineWidth(width);
}

void StyleTH1Fill(std::unique_ptr<TH1> &h, int color, int style = 1001,
                  double alpha = 1) {
  if (!h) {
    return;
  }
  h->SetFillStyle(style);
  h->SetFillColorAlpha(color, alpha);
}

double GetMaximumTH1s(
    std::vector<std::reference_wrapper<std::unique_ptr<TH1>>> hs) {
  double max = -std::numeric_limits<double>::max();

  for (auto const &h : hs) {
    if (!h.get()) {
      continue;
    }
    max = std::max(max, h.get()->GetMaximum());
  }

  return max;
}

void DrawTH1s(std::vector<std::reference_wrapper<std::unique_ptr<TH1>>> hs,
              std::string opts = "", bool first = true) {
  for (auto const &h : hs) {
    if (!h.get()) {
      continue;
    }
    if (first) {
      h.get()->GetYaxis()->SetRangeUser(0, GetMaximumTH1s(hs) * 1.05);
    }
    std::string dopt = std::string(!first ? "SAME" : "") + opts;
    h.get()->DrawClone(dopt.c_str());
    first = false;
  }
}

void ScaleTH1s(std::vector<std::reference_wrapper<std::unique_ptr<TH1>>> hs,
               double s, std::string opts = "") {
  for (auto &h : hs) {
    if (!h.get()) {
      continue;
    }
    h.get()->Scale(s, opts.c_str());
  }
}

double GetMaximumTH1s(std::vector<std::unique_ptr<TH1>> &hs) {
  std::vector<std::reference_wrapper<std::unique_ptr<TH1>>> rhs;
  for (auto &h : hs) {
    rhs.push_back(h);
  }
  return GetMaximumTH1s(rhs);
}

void DrawTH1s(std::vector<std::unique_ptr<TH1>> &hs, std::string opts = "",
              bool first = true) {
  std::vector<std::reference_wrapper<std::unique_ptr<TH1>>> rhs;
  for (auto &h : hs) {
    rhs.push_back(h);
  }
  DrawTH1s(rhs, opts, first);
}

void ScaleTH1s(std::vector<std::unique_ptr<TH1>> &hs, double s,
               std::string opts = "") {
  std::vector<std::reference_wrapper<std::unique_ptr<TH1>>> rhs;
  for (auto &h : hs) {
    rhs.push_back(h);
  }
  ScaleTH1s(rhs, s, opts);
}

TCanvas *MakeCanvas(double mt = 0.03, double ml = 0.15, double mr = 0.03,
                    double mb = 0.15, std::string name = "c1") {
  TCanvas *c1 = new TCanvas(name.c_str(), "", 600, 600);
  c1->SetTopMargin(mt);
  c1->SetLeftMargin(ml);
  c1->SetRightMargin(mr);
  c1->SetBottomMargin(mb);
  return c1;
}

TPad *MakeRatioTopPad(double mt = 0.05, double ml = 0.15, double mr = 0.03,
                      double mb = 0.03, double midpoint = 0.4,
                      std::string name = "ptop") {
  TPad *c1 = new TPad(name.c_str(), "", 0, midpoint, 1, 1);
  c1->SetTopMargin(mt);
  c1->SetLeftMargin(ml);
  c1->SetRightMargin(mr);
  c1->SetBottomMargin(mb);
  return c1;
}
TPad *MakeRatioTopPadTopLegend(double mt = 0.3, double ml = 0.15,
                               double mr = 0.03, double mb = 0.03,
                               double midpoint = 0.3,
                               std::string name = "ptop") {
  return MakeRatioTopPad(mt, ml, mr, mb, midpoint, name);
}
TPad *MakeRatioBottomPad(double mt = 0.03, double ml = 0.15, double mr = 0.03,
                         double mb = 0.2, double midpoint = 0.4,
                         std::string name = "pbottom") {
  TPad *c1 = new TPad(name.c_str(), "", 0, 0, 1, midpoint);
  c1->SetTopMargin(mt);
  c1->SetLeftMargin(ml);
  c1->SetRightMargin(mr);
  c1->SetBottomMargin(mb);
  return c1;
}
TPad *MakeRatioBottomPadTopLegend(double mt = 0.03, double ml = 0.15,
                                  double mr = 0.03, double mb = 0.25,
                                  double midpoint = 0.3,
                                  std::string name = "pbottom") {
  return MakeRatioBottomPad(mt, ml, mr, mb, midpoint, name);
}

TCanvas *MakeCanvasTopLegend(std::string name = "c1") {
  return MakeCanvas(0.2, 0.15, 0.03, 0.15, name);
}

TLegend *MakeTopLegend() {
  auto leg = new TLegend(0, 0.81, 0.97, 1);
  leg->SetNColumns(2);
  leg->SetTextSize(0.07);
  leg->SetBorderSize(0);
  leg->SetFillStyle(4001);
  return leg;
}

TGaxis *MakeOtherAxis(double rightmax, double &scale) {
  scale = gPad->GetUymax() / rightmax;

  return new TGaxis(gPad->GetUxmax(), gPad->GetUymin(), gPad->GetUxmax(),
                    gPad->GetUymax(), 0, rightmax, 510, "+L");
}